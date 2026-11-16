# Public web page (developer notes)

The pipeline can also be run from a public web page, without any account. The design is the same as [SHAMAN](https://github.com/aghozlane/shaman):

<pre>
browser -> web page (web/app/app.R, Shiny)  -- writes -->  todo/KEY.json + jobs/KEY/{sample.zip, metadata.tsv, params.json}
                                                                  |
           web/daemon/rp_bioblend.py (same container, Galaxy service account, BioBlend)
                                                                  |
           galaxy.pasteur.fr: repertoire_profiler tool (galaxy/repertoire_profiler.xml)  -- Slurm -->  Maestro
</pre>

- The values of the web form are sent to the Galaxy tool ([galaxy/README.md](../galaxy/README.md)) as its JSON parameter file.
- *web/app/app.R*: the web page. The form is built from *nextflow.config* at start-up, one field per parameter, with the comment of the parameter as help. Users get a key to follow the job and download the results ("Check my job" tab, or link in the e-mail).
- *web/daemon/rp_bioblend.py*: sends each job to Galaxy, waits, downloads the results, sends an e-mail and removes old jobs. Settings are environment variables (see the top of the file). The Galaxy API key must never be written in the code.
- *web/Dockerfile*: image of the web page and the daemon:

```
docker build -f web/Dockerfile -t repertoire_profiler_web .
docker run -p 3838:3838 -e GALAXY_API_KEY=... -e RP_PUBLIC_URL=https://... -v rp_work:/srv/rp_work repertoire_profiler_web
```
