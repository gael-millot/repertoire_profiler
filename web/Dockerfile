# Web page of repertoire_profiler: Shiny form (web/app/app.R) + submission daemon to Galaxy (web/daemon/rp_bioblend.py).
# Build from the root of the repository:
#   docker build -f web/Dockerfile -t repertoire_profiler_web .
# Run:
#   docker run -p 3838:3838 -e GALAXY_API_KEY=... -e RP_PUBLIC_URL=https://... -v rp_work:/srv/rp_work repertoire_profiler_web
FROM rocker/shiny:4.4.2

RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-pip \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 install --no-cache-dir --break-system-packages bioblend==1.9.0 \
    && rm -rf /srv/shiny-server/*

# only the files used by the web page and the daemon (the pipeline itself runs in Galaxy)
COPY nextflow.config /opt/repertoire_profiler/nextflow.config
COPY galaxy/fill_config.py galaxy/locked_params.txt /opt/repertoire_profiler/galaxy/
COPY web/daemon/rp_bioblend.py /opt/repertoire_profiler/web/daemon/
COPY web/app/app.R /srv/shiny-server/app.R
COPY web/start.sh /usr/bin/start.sh

RUN chmod +x /usr/bin/start.sh \
    && mkdir -p /srv/rp_work && chown -R shiny:shiny /srv/rp_work /srv/shiny-server

ENV RP_WORK_DIR=/srv/rp_work \
    RP_REPO_DIR=/opt/repertoire_profiler \
    RP_PYTHON=python3 \
    RP_MAX_UPLOAD_MB=50 \
    RP_RETENTION_DAYS=15

VOLUME /srv/rp_work
EXPOSE 3838
CMD ["/usr/bin/start.sh"]
