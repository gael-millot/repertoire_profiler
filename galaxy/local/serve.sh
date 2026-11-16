#!/usr/bin/env bash
#########################################################################
##                                                                     ##
##     serve.sh                                                        ##
##                                                                     ##
##     Local Galaxy web page to run repertoire_profiler (planemo)      ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################

# Starts a local Galaxy with the repertoire_profiler tool (planemo serve), then:
#   1) open http://127.0.0.1:9090 in a web browser (port: RP_PORT variable)
#   2) upload the fasta files (or their zip) and the optional metadata table (Upload button, top left)
#   3) open the repertoire_profiler tool: its form recapitulates nextflow.config, one field per parameter
#   4) click on "Run tool"
#   5) while the job runs (yellow), follow the nextflow processes:
#      - in the Galaxy page: click on the running dataset, then on its i icon (Dataset details): the Tool Standard Output
#        shows the lines of nextflow as they come (one line per started task, refreshed every 3 seconds)
#      - or in a second terminal, the progress table of the nextflow console: python3 galaxy/local/watch.py
#   6) when the job is green, open the "report and result folder" dataset: the eye icon shows the report,
#      the download icon saves the whole result folder as a zip, where the web browser proposes.
# The histories and results are kept between restarts, in the planemo profile (RP_PROFILE variable).
# Stop the server with Ctrl+C.
#
# Usage, from any folder:
#   bash galaxy/local/serve.sh
# Settings: see settings.sh, plus
#   RP_PORT      port of the web server (default: 9090)
#   RP_HOST      address of the web server (default: 127.0.0.1, i.e., this computer only)
#   RP_PROFILE   planemo profile keeping the Galaxy database and datasets (default: repertoire_profiler)
# Requirements: planemo (python3 -m venv ~/.venv && ~/.venv/bin/pip install planemo), git, and internet access.
# The first start installs Galaxy (about 15 minutes), the first run of the tool installs nextflow and apptainer (conda),
# then pulls the container images of the pipeline (about 10 GB, long).

source "$(dirname "${BASH_SOURCE[0]}")/settings.sh"

RP_PORT="${RP_PORT:-9090}"
RP_HOST="${RP_HOST:-127.0.0.1}"
RP_PROFILE="${RP_PROFILE:-repertoire_profiler}"

# profile: sqlite database and datasets kept in ~/.planemo/profiles/<RP_PROFILE>
if ! "${PLANEMO}" profile_list 2>/dev/null | grep -q -x -F "${RP_PROFILE}"; then
    "${PLANEMO}" profile_create "${RP_PROFILE}"
fi

# the HTML report keeps its scripts and styles: the tool is in the sanitization allowlist of this Galaxy
# (on a Galaxy server, the administrators do the same with Admin > Allowlist)
printf 'repertoire_profiler\n' > "${RP_WORK_DIR}/sanitize_allowlist.txt"
export GALAXY_CONFIG_OVERRIDE_SANITIZE_ALLOWLIST_FILE="${RP_WORK_DIR}/sanitize_allowlist.txt"

echo "Galaxy will be at http://${RP_HOST}:${RP_PORT} (profile ${RP_PROFILE})"
exec "${PLANEMO}" serve \
    "${PLANEMO_GALAXY_OPTIONS[@]}" \
    --profile "${RP_PROFILE}" \
    --host "${RP_HOST}" \
    --port "${RP_PORT}" \
    "${RP_TOOL_XML}"
