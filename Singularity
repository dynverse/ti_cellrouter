#######################################################################################
## DO NOT EDIT THIS FILE! This file was automatically generated from the dockerfile. ##
## Run dynwrap:::.container_dockerfile_to_singularityrecipe() to update this file.   ##
#######################################################################################

Bootstrap: shub

From: dynverse/dynwrap:bioc

%labels
    version 0.1.2

%post
    chmod -R a+r /code
    chmod a+x /code
    R -e 'devtools::install_cran(c("reshape", "reshape2", "pheatmap", "tsne", "igraph", "ggplot2", "mclust", "grid", "Rtsne", "cccd", "irlba"))'
    git clone https://github.com/edroaldo/cellrouter.git && find cellrouter -type f | grep -v "^cellrouter/CellRouter" | xargs rm
    apt-get update && apt-get install -y default-jre

%files
    . /code

%runscript
    exec Rscript /code/run.R

