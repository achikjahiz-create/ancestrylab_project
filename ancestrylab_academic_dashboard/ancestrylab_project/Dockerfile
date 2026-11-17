FROM rocker/shiny:4.4.2
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev libgsl-dev gfortran make g++ \
    libfontconfig1-dev libfreetype6-dev libpng-dev libtiff-dev libjpeg-dev libharfbuzz-dev libfribidi-dev python3 python3-requests zip \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /srv/ancestrylab
COPY . /srv/ancestrylab
RUN Rscript scripts/install.R
EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('app', host='0.0.0.0', port=as.integer(Sys.getenv('PORT', '3838')))"]
