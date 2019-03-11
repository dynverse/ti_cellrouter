FROM dynverse/dynwrapr:v0.1.0

RUN R -e 'remotes::install_cran(c("reshape", "reshape2", "pheatmap", "tsne", "igraph", "ggplot2", "mclust", "Rtsne", "cccd", "irlba"))'

RUN git clone --depth 1 https://github.com/edroaldo/cellrouter.git && find cellrouter -type f | grep -v "^cellrouter/CellRouter" | xargs rm

RUN apt-get update && apt-get install -y default-jre

COPY run.R example.R definition.yml /code/

ENTRYPOINT ["/code/run.R"]
