#!/usr/bin/env bash
#########################################################################
##                                                                     ##
##     test.sh                                                         ##
##                                                                     ##
##     Lint and test of the repertoire_profiler Galaxy tool (planemo)  ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################

# 1) checks that the form of galaxy/repertoire_profiler.xml still matches nextflow.config (check_tool.py)
# 2) planemo lint of the tool, with the conda requirements
# 3) planemo shed_lint of the Tool Shed repository (.shed.yml at the root of the repository)
# 4) planemo test: runs the tests of the tool in a temporary local Galaxy. Each test runs the whole pipeline
#    (about 15 minutes each, much longer the first time, which pulls the container images).
#    Progress of the nextflow processes of the running test, in a second terminal: python3 galaxy/local/watch.py
#    Report: ${RP_WORK_DIR}/tool_test_output.html
#
# Usage, from any folder:
#   bash galaxy/local/test.sh          # all the steps
#   bash galaxy/local/test.sh lint     # steps 1 to 3 only
# Settings: see settings.sh

MODE="${1:-all}"
source "$(dirname "${BASH_SOURCE[0]}")/settings.sh"
REPO_DIR="$(cd "${RP_TOOL_DIR}/.." && pwd)"

echo -e "\n#### Form of the tool versus nextflow.config\n"
python3 "${RP_TOOL_DIR}/local/check_tool.py" --template "${REPO_DIR}/nextflow.config" --locked "${RP_TOOL_DIR}/locked_params.txt" --tool "${RP_TOOL_XML}"

echo -e "\n#### planemo lint\n"
"${PLANEMO}" lint --fail_level warn --conda_requirements "${RP_TOOL_XML}"

echo -e "\n#### planemo shed_lint\n"
# fail level error: the remote_repository_url warning is expected (planemo expects URLs of tools-iuc like repositories, ending by the local path)
(cd "${REPO_DIR}" && "${PLANEMO}" shed_lint --tools --ensure_metadata --fail_level error .)

if [[ "${MODE}" == "lint" ]]; then
    exit 0
fi

echo -e "\n#### planemo test\n"
cd "${RP_WORK_DIR}"
"${PLANEMO}" test \
    "${PLANEMO_GALAXY_OPTIONS[@]}" \
    --test_output "${RP_WORK_DIR}/tool_test_output.html" \
    --test_output_json "${RP_WORK_DIR}/tool_test_output.json" \
    "${RP_TOOL_XML}"
