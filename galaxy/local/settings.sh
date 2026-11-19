#!/usr/bin/env bash
#########################################################################
##                                                                     ##
##     settings.sh                                                     ##
##                                                                     ##
##     Common settings of serve.sh and test.sh (sourced, not run)      ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################

# Settings of the local Galaxy used by planemo to run the repertoire_profiler tool (galaxy/repertoire_profiler.xml).
# All of them can be changed with environment variables, e.g., RP_CPUS=8 bash galaxy/local/serve.sh
#   PLANEMO               planemo executable (default: planemo of the PATH, else ~/.venv/bin/planemo)
#   RP_GALAXY_BRANCH      Galaxy release used by planemo (default: release_26.1)
#   RP_GALAXY_ROOT        folder of this Galaxy, cloned at the first use (default: ~/.planemo/galaxy_<RP_GALAXY_BRANCH>)
#   RP_CONDA_PREFIX       conda used by Galaxy to install the requirements of the tool: nextflow, apptainer, graphviz and python
#                         (default: ~/miniforge3, installed at the first use)
#   RP_APPTAINER_CACHE    folder of the apptainer images of the pipeline (default: ~/apptainer). The first run pulls them (long)
#   RP_APPTAINER_OPTIONS  apptainer run options of the pipeline processes (default: --no-home). The tool adds the bind of
#                         the temporary folder of the job (TMPDIR)
#   RP_CPUS               cores of a repertoire_profiler job (default: all the cores of the computer)
#   RP_MEMORY_MB          memory of a repertoire_profiler job in MB (default: 80% of the memory of the computer)
#   RP_WORK_DIR           folder of the generated job configuration and of the test reports (default: ~/.planemo/repertoire_profiler)
#   RP_TMPDIR             temporary folder of planemo and of the Galaxy jobs, on disk (default: RP_WORK_DIR/tmp)

set -euo pipefail

RP_TOOL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)" # galaxy folder of the repository
RP_TOOL_XML="${RP_TOOL_DIR}/repertoire_profiler.xml"

PLANEMO="${PLANEMO:-$(command -v planemo || echo "${HOME}/.venv/bin/planemo")}"
RP_GALAXY_BRANCH="${RP_GALAXY_BRANCH:-release_26.1}"
RP_GALAXY_ROOT="${RP_GALAXY_ROOT:-${HOME}/.planemo/galaxy_${RP_GALAXY_BRANCH}}"
RP_CONDA_PREFIX="${RP_CONDA_PREFIX:-${HOME}/miniforge3}"
RP_APPTAINER_CACHE="${RP_APPTAINER_CACHE:-${HOME}/apptainer}"
RP_APPTAINER_OPTIONS="${RP_APPTAINER_OPTIONS:---no-home}"
RP_CPUS="${RP_CPUS:-$(nproc)}"
RP_MEMORY_MB="${RP_MEMORY_MB:-$(awk '/^MemTotal:/ {printf "%d", $2 / 1024 * 0.8}' /proc/meminfo)}"
RP_WORK_DIR="${RP_WORK_DIR:-${HOME}/.planemo/repertoire_profiler}"

if [[ ! -x "${PLANEMO}" ]]; then
    echo "ERROR: planemo not found. Install it (python3 -m venv ~/.venv && ~/.venv/bin/pip install planemo) or set the PLANEMO variable." >&2
    exit 1
fi
for i1 in RP_CPUS RP_MEMORY_MB; do
    if [[ ! "${!i1}" =~ ^[1-9][0-9]*$ ]]; then
        echo "ERROR: ${i1} must be a positive integer: ${!i1}" >&2
        exit 1
    fi
done
for i1 in RP_APPTAINER_CACHE RP_APPTAINER_OPTIONS RP_WORK_DIR RP_CONDA_PREFIX; do
    if [[ "${!i1}" == *'"'* || "${!i1}" == *'\'* ]]; then
        echo "ERROR: ${i1} cannot contain a double quote or a backslash: ${!i1}" >&2
        exit 1
    fi
done

mkdir -p "${RP_APPTAINER_CACHE}" "${RP_WORK_DIR}"

# temporary folder of planemo, hence of the Galaxy jobs (work folder of nextflow, builds of the apptainer images):
# on disk, since /tmp can be a small folder in memory (tmpfs, e.g., in WSL)
export TMPDIR="${RP_TMPDIR:-${RP_WORK_DIR}/tmp}"
mkdir -p "${TMPDIR}"

# Galaxy used by planemo: a fixed release, kept between runs (planemo installs its virtual environment and web client at the first use)
if [[ ! -d "${RP_GALAXY_ROOT}" ]]; then
    echo "Cloning Galaxy ${RP_GALAXY_BRANCH} into ${RP_GALAXY_ROOT}"
    # core.autocrlf=false: a global core.autocrlf=true would give CRLF line endings, which break the Galaxy shell scripts
    git -c core.autocrlf=false clone --depth 1 --branch "${RP_GALAXY_BRANCH}" https://github.com/galaxyproject/galaxy.git "${RP_GALAXY_ROOT}"
fi

# Galaxy job configuration: other tools run with 1 core, the repertoire_profiler tool runs with RP_CPUS cores
# (GALAXY_SLOTS) and RP_MEMORY_MB MB (GALAXY_MEMORY_MB), with the apptainer settings read by the tool command.
# tmp_dir: the temporary folder of a job (TMPDIR, TMP and TEMP) is the tmp folder of its job folder, on disk, mounted in the
# containers by the tool command.
# live_tool_output_reporting: the Galaxy page shows the standard output of the running job (one line per nextflow task),
# in the dataset details (i icon of the running dataset), Tool Standard Output. Progress table in a terminal: watch.py
RP_JOB_CONF="${RP_WORK_DIR}/job_conf.yml"
cat > "${RP_JOB_CONF}" <<EOF
# Generated by ${RP_TOOL_DIR}/local/settings.sh: changes are overwritten at the next run
runners:
  local:
    load: galaxy.jobs.runners.local:LocalJobRunner
    workers: 4
execution:
  default: local
  environments:
    local:
      runner: local
      tmp_dir: true
    repertoire_profiler:
      runner: local
      local_slots: ${RP_CPUS}
      tmp_dir: true
      live_tool_output_reporting: true
      env:
        - name: GALAXY_MEMORY_MB
          value: "${RP_MEMORY_MB}"
        - name: RP_APPTAINER_CACHE
          value: "${RP_APPTAINER_CACHE}"
        - name: RP_APPTAINER_OPTIONS
          value: "${RP_APPTAINER_OPTIONS}"
tools:
  - class: local
    environment: local
  - id: repertoire_profiler
    environment: repertoire_profiler
EOF

# options common to planemo serve and planemo test
PLANEMO_GALAXY_OPTIONS=(
    --galaxy_root "${RP_GALAXY_ROOT}"
    --job_config_file "${RP_JOB_CONF}"
    --conda_prefix "${RP_CONDA_PREFIX}"
    --conda_auto_init
    --conda_auto_install
)

echo "planemo:          ${PLANEMO}"
echo "Galaxy:           ${RP_GALAXY_ROOT}"
echo "conda:            ${RP_CONDA_PREFIX}"
echo "apptainer images: ${RP_APPTAINER_CACHE} (options: ${RP_APPTAINER_OPTIONS})"
echo "job resources:    ${RP_CPUS} cores, ${RP_MEMORY_MB} MB"
