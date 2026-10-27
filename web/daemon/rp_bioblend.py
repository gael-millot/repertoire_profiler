#!/usr/bin/env python3
"""
Submission daemon of the repertoire_profiler web page (same principle as shaman_bioblend.py of SHAMAN).

The web page (web/app/app.R) writes, for each submission with key KEY:
    <work_dir>/jobs/KEY/sample.zip        fasta files
    <work_dir>/jobs/KEY/metadata.tsv      optional
    <work_dir>/jobs/KEY/params.json       values of the form
    <work_dir>/todo/KEY.json              task: {"key", "mail", "submitted"}

This daemon polls todo/, moves each task to doing/, and in a thread:
    1) creates a Galaxy history and uploads the files,
    2) runs the repertoire_profiler Galaxy tool (galaxy/repertoire_profiler.xml),
    3) waits for the job, downloads the results into jobs/KEY/ and moves the task to done/ (or error/),
    4) sends an e-mail if an address was given.
The task file always holds the current "state" and "message", read by the web page for the "Check my job" tab.
Jobs older than RP_RETENTION_DAYS are removed, locally and in Galaxy.

Configuration (environment variables, never in the code):
    GALAXY_API_KEY      required, key of the Galaxy service account
    GALAXY_URL          default https://galaxy.pasteur.fr
    RP_TOOL_ID          default repertoire_profiler (full Tool Shed id if installed from a Tool Shed)
    RP_MAX_JOBS         max number of jobs running at the same time in Galaxy, default 5
    RP_RETENTION_DAYS   default 15
    RP_PUBLIC_URL       URL of the web page, written in the e-mails
    RP_SMTP_HOST        default smtp.pasteur.fr ("" for no e-mail)
    RP_MAIL_FROM        default repertoire_profiler@pasteur.fr

Usage:
    rp_bioblend.py -w <work_dir> [--once]
"""

import argparse
import datetime
import json
import logging
import os
import re
import shutil
import smtplib
import sys
import threading
import time
from email.message import EmailMessage

import requests
from bioblend import ConnectionError as GalaxyConnectionError
from bioblend.galaxy import GalaxyInstance

KEY_RE = re.compile(r"^[A-Za-z0-9_-]{16,64}$")
POLL_SECONDS = 5
JOB_POLL_SECONDS = 30
RETRY_SECONDS = 300  # when Galaxy is unreachable (e.g., maintenance)

log = logging.getLogger("rp_bioblend")


class Config:
    def __init__(self, work_dir):
        self.work_dir = os.path.abspath(work_dir)
        self.galaxy_url = os.environ.get("GALAXY_URL", "https://galaxy.pasteur.fr")
        self.galaxy_key = os.environ.get("GALAXY_API_KEY", "")
        self.tool_id = os.environ.get("RP_TOOL_ID", "repertoire_profiler")
        self.max_jobs = int(os.environ.get("RP_MAX_JOBS", "5"))
        self.retention_days = int(os.environ.get("RP_RETENTION_DAYS", "15"))
        self.public_url = os.environ.get("RP_PUBLIC_URL", "")
        self.smtp_host = os.environ.get("RP_SMTP_HOST", "smtp.pasteur.fr")
        self.mail_from = os.environ.get("RP_MAIL_FROM", "repertoire_profiler@pasteur.fr")
        if not self.galaxy_key:
            sys.exit("rp_bioblend.py: ERROR: the GALAXY_API_KEY environment variable is not set")

    def path(self, *parts):
        return os.path.join(self.work_dir, *parts)


def now():
    return datetime.datetime.now().isoformat(timespec="seconds")


def read_task(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def write_task(path, task):
    # write then rename, so that the web page never reads a half written file
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(task, f, indent=1)
    os.replace(tmp, path)


class Job(threading.Thread):
    """One submission, from doing/KEY.json to done/ or error/."""

    def __init__(self, cfg, key, slots):
        super().__init__(name="job_" + key, daemon=True)
        self.cfg = cfg
        self.key = key
        self.slots = slots
        self.task_file = cfg.path("doing", key + ".json")
        self.job_dir = cfg.path("jobs", key)
        self.task = read_task(self.task_file)

    def update(self, state, message):
        self.task.update({"state": state, "message": message, "updated": now()})
        write_task(self.task_file, self.task)
        log.info("%s: %s: %s", self.key, state, message)

    def finish(self, folder, state, message):
        self.update(state, message)
        dest = self.cfg.path(folder, self.key + ".json")
        os.replace(self.task_file, dest)
        self.task_file = dest
        self.send_mail(state, message)

    def run(self):
        with self.slots:
            while True:
                try:
                    self.process()
                    return
                except (GalaxyConnectionError, requests.exceptions.RequestException) as e:
                    status = getattr(e, "status_code", None)
                    if status is not None and status < 500:
                        log.exception("%s: Galaxy refused the request", self.key)
                        self.finish("error", "error", "Galaxy refused the job (HTTP {0}). Please contact the administrators.".format(status))
                        return
                    # Galaxy down or under maintenance: keep the task in doing/ and try again later
                    log.warning("%s: Galaxy connection error, retry in %s s: %s", self.key, RETRY_SECONDS, e)
                    self.update("waiting", "The Galaxy server is not reachable for the moment. The job will be submitted again automatically.")
                    time.sleep(RETRY_SECONDS)
                except Exception as e:
                    log.exception("%s: failure", self.key)
                    self.finish("error", "error", "Internal error of the submission system: {0}".format(e))
                    return

    def process(self):
        gi = GalaxyInstance(url=self.cfg.galaxy_url, key=self.cfg.galaxy_key)
        if not self.task.get("job_id"):
            self.submit(gi)
        self.wait(gi)

    def submit(self, gi):
        if not self.task.get("history_id"):
            history = gi.histories.create_history(name="repertoire_profiler_" + self.key)
            self.task["history_id"] = history["id"]
            self.update("uploading", "Uploading the files to Galaxy.")
        hid = self.task["history_id"]
        inputs = {}
        uploads = [("fasta_zip", "sample.zip", "zip"), ("metadata", "metadata.tsv", "tabular"), ("params", "params.json", "json")]
        for name, file_name, file_type in uploads:
            path = os.path.join(self.job_dir, file_name)
            if os.path.isfile(path):
                out = gi.tools.upload_file(path, hid, file_type=file_type, file_name=file_name)
                inputs[name] = {"src": "hda", "id": out["outputs"][0]["id"]}
        for value in inputs.values():
            gi.datasets.wait_for_dataset(value["id"], maxwait=3600)
        run = gi.tools.run_tool(hid, self.cfg.tool_id, inputs)
        self.task["job_id"] = run["jobs"][0]["id"]
        self.task["outputs"] = {o["output_name"]: o["id"] for o in run["outputs"]}
        self.update("running", "The analysis is running on the Institut Pasteur cluster.")

    def wait(self, gi):
        job_id = self.task["job_id"]
        while True:
            job = gi.jobs.show_job(job_id)
            state = job["state"]
            if state == "ok":
                break
            if state in ("error", "deleted", "deleting", "failed", "stopped"):
                details = gi.jobs.show_job(job_id, full_details=True)
                with open(os.path.join(self.job_dir, "error.txt"), "w", encoding="utf-8") as f:
                    f.write("stdout:\n{0}\n\nstderr:\n{1}\n".format(details.get("stdout", ""), details.get("stderr", "")))
                self.finish("error", "error", "The analysis failed. See the error log of the job, and check your files and parameters.")
                return
            if self.task.get("state") != "running":
                self.update("running", "The analysis is running on the Institut Pasteur cluster.")
            time.sleep(JOB_POLL_SECONDS)
        self.update("downloading", "Downloading the results.")
        file_names = {"results": "results.zip", "report": "report.html"}
        for name, dataset_id in self.task["outputs"].items():
            if name in file_names:
                gi.datasets.download_dataset(dataset_id, file_path=os.path.join(self.job_dir, file_names[name]),
                                             use_default_filename=False)
        self.finish("done", "done", "The analysis is finished.")

    def send_mail(self, state, message):
        mail = self.task.get("mail", "")
        if not mail or not self.cfg.smtp_host:
            return
        msg = EmailMessage()
        msg["Subject"] = "repertoire_profiler job {0}: {1}".format(self.key, state)
        msg["From"] = self.cfg.mail_from
        msg["To"] = mail
        text = "{0}\n\nYour job key is: {1}\n".format(message, self.key)
        if self.cfg.public_url:
            text += "Results and status: {0}?key={1}\n".format(self.cfg.public_url, self.key)
        text += "\nResults are kept {0} days.\n".format(self.cfg.retention_days)
        msg.set_content(text)
        try:
            with smtplib.SMTP(self.cfg.smtp_host, 25, timeout=60) as server:
                server.send_message(msg)
        except Exception:
            log.exception("%s: cannot send the e-mail", self.key)


def cleanup(cfg, gi):
    """Remove the jobs (local files and Galaxy history) older than the retention time."""
    limit = time.time() - cfg.retention_days * 86400
    for folder in ("done", "error"):
        for name in os.listdir(cfg.path(folder)):
            if not name.endswith(".json"):
                continue
            path = cfg.path(folder, name)
            if os.path.getmtime(path) > limit:
                continue
            key = name[:-5]
            try:
                task = read_task(path)
                if task.get("history_id"):
                    gi.histories.delete_history(task["history_id"], purge=True)
            except Exception:
                log.exception("%s: cannot delete the Galaxy history", key)
                continue
            shutil.rmtree(cfg.path("jobs", key), ignore_errors=True)
            os.remove(path)
            log.info("%s: removed (retention time)", key)


def start_task(cfg, key, slots, running):
    if key in running and running[key].is_alive():
        return
    job = Job(cfg, key, slots)
    running[key] = job
    job.start()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-w", dest="work_dir", required=True, help="top folder shared with the web page")
    parser.add_argument("--once", action="store_true", help="process the waiting tasks, wait for them and exit (tests)")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(threadName)s %(message)s")
    cfg = Config(args.work_dir)
    for folder in ("todo", "doing", "done", "error", "jobs"):
        os.makedirs(cfg.path(folder), exist_ok=True)

    slots = threading.BoundedSemaphore(cfg.max_jobs)
    running = {}
    # tasks interrupted by a restart of the daemon: resume them
    for name in sorted(os.listdir(cfg.path("doing"))):
        if name.endswith(".json"):
            start_task(cfg, name[:-5], slots, running)

    last_cleanup = 0
    while True:
        for name in sorted(os.listdir(cfg.path("todo"))):
            if not name.endswith(".json"):
                continue
            key = name[:-5]
            src = cfg.path("todo", name)
            if not KEY_RE.match(key) or not os.path.isdir(cfg.path("jobs", key)):
                log.error("invalid task %s, moved to error/", name)
                os.replace(src, cfg.path("error", name))
                continue
            dest = cfg.path("doing", name)
            os.replace(src, dest)
            task = read_task(dest)
            task.update({"state": "queued", "message": "The job is waiting to be sent to Galaxy.", "updated": now()})
            write_task(dest, task)
            start_task(cfg, key, slots, running)
        if time.time() - last_cleanup > 3600:
            try:
                cleanup(cfg, GalaxyInstance(url=cfg.galaxy_url, key=cfg.galaxy_key))
            except Exception:
                log.exception("cleanup failed")
            last_cleanup = time.time()
        if args.once and not any(j.is_alive() for j in running.values()):
            return
        time.sleep(POLL_SECONDS)


if __name__ == "__main__":
    main()
