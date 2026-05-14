FROM ghcr.io/bioconductor/bioconductor_docker:RELEASE_3_21

COPY . /opt/HPCell
RUN R -q -e "desc <- read.dcf('/opt/HPCell/DESCRIPTION'); pkgs <- unique(trimws(unlist(strsplit(paste(desc[1, c('Depends', 'Imports')], collapse = ','), ',')))); pkgs <- sub('\\\\s*\\\\(.*\\\\)', '', pkgs); pkgs <- pkgs[pkgs != '' & pkgs != 'R']; repos <- c(CRAN = 'https://cloud.r-project.org', BioCsoft = 'https://bioconductor.org/packages/3.21/bioc', BioCann = 'https://bioconductor.org/packages/3.21/data/annotation', BioCexp = 'https://bioconductor.org/packages/3.21/data/experiment', BioCworkflows = 'https://bioconductor.org/packages/3.21/workflows'); install.packages(pkgs, repos = repos)"
RUN R CMD INSTALL /opt/HPCell

WORKDIR /work
CMD ["R"]
