FROM dynverse/dynwrap_latest:latest

ARG GITHUB_PAT

RUN R -e 'devtools::install_cran(c("dplyr", "purrr"))'

RUN R -e 'devtools::install_github("tianlt/dbcti")'

COPY definition.yml run.R dbcti.R /code/

ENTRYPOINT ["/code/run.R"]