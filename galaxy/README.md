# Galaxy tool (local Galaxy with planemo, galaxy.pasteur.fr)

Developer notes. To run the tool on galaxy.pasteur.fr: [HOW TO RUN](../README.md#how-to-run) section of the README.

The *galaxy/repertoire_profiler.xml* Galaxy tool runs the whole pipeline inside a single Galaxy job (nextflow local executor, processes in apptainer containers). Its form, i.e., the web page of the tool in Galaxy, recapitulates the *nextflow.config* file: same sections, one field per parameter, with the name of the parameter in parentheses and its comment as help. The *Run tool* button runs the pipeline. The execution parameters (*Local / Cluster execution* section) are set by Galaxy, and the iTOL upload is disabled.

**Results**: the *report and result folder* dataset of the Galaxy history is the HTML report, which also carries the whole result folder (see [OUTPUT](../README.md#output)): the eye icon shows the report, whose links to the tables, pdf and alignments work in Galaxy, and **the download icon saves the whole result folder as a zip file**, where the web browser proposes (enable "Ask where to save each file" in the browser settings to choose the folder each time). Optionally, the tool also returns the .tsv tables as a collection (for other Galaxy tools) and the result folder as a separate zip dataset. On a Galaxy server, the datasets can also be exported to a storage configured by the administrators (Export datasets).

**Local Galaxy with planemo** ([planemo](https://planemo.readthedocs.io), e.g., in ~/.venv: `python3 -m venv ~/.venv && ~/.venv/bin/pip install planemo`):

<pre>
bash galaxy/local/serve.sh    # local Galaxy with the tool, at http://127.0.0.1:9090 (Ctrl+C to stop)
bash galaxy/local/test.sh     # lint (tool and Tool Shed repository) and tests of the tool (each test runs the whole pipeline)
</pre>

Then, in the web browser: upload the fasta files (or their zip) and the optional metadata table (*Upload* button), open the *repertoire_profiler* tool, set the parameters and click on *Run tool*. The histories are kept between restarts (planemo profile). The first start installs Galaxy (release_26.1, in ~/.planemo), the first run installs nextflow, apptainer, graphviz and python with conda (~/miniforge3), then pulls the container images of the pipeline into ~/apptainer (about 10 GB, long). Settings (cores, memory, folders) are environment variables, described in *galaxy/local/settings.sh*, e.g., `RP_CPUS=8 RP_MEMORY_MB=12000 bash galaxy/local/serve.sh`.

**Tool Shed and galaxy.pasteur.fr**: the *.shed.yml* file at the root of the repository describes the Tool Shed repository: the tool with its test data and the pipeline files (main.nf, nextflow.config, bin, conf and modules folders), in the same layout as this repository. After `bash galaxy/local/test.sh lint`, from the root of the repository (API key of the Tool Shed account):

<pre>
planemo shed_create --shed_target https://toolshed.pasteur.fr --shed_key API_KEY .                 # first upload
planemo shed_update --shed_target https://toolshed.pasteur.fr --shed_key API_KEY --check_diff .    # next versions
</pre>

The Galaxy administrators then install the repository from the Tool Shed and set the RP_APPTAINER_CACHE job environment variable (folder of the apptainer images, shared by the jobs), as described at the top of *galaxy/repertoire_profiler.xml*. If the computing nodes have no internet access, the images are pulled beforehand into this folder with `bash galaxy/pull_images.sh <folder>`.

- *galaxy/repertoire_profiler.xml*: the Galaxy tool. Inputs: fasta files or their zip, optional metadata table, the parameters of *nextflow.config*, and an optional JSON parameter file (whose values replace those of the form, e.g., to reuse a set of parameters).
- *galaxy/fill_config.py*: writes the *nextflow.config* of a job from the form (parameters are in env {} blocks and cannot be changed with --param on the command line). User values are strictly checked (no quote, backslash or $), so that no Groovy code can be injected into the config.
- *galaxy/locked_params.txt*: parameters that users cannot set (execution, paths, keys).
- *galaxy/pull_images.sh*: pulls the container images of *nextflow.config* into a folder, with the file names used by nextflow.
- *galaxy/local/*: planemo scripts (not in the Tool Shed repository). *check_tool.py*, run by *test.sh*, checks that the form still matches *nextflow.config*: run it after any change of the parameters of *nextflow.config*.
