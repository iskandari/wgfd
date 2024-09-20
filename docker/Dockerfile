# Start from the rocker/r-ver:latest image
FROM rocker/geospatial:latest

LABEL maintainer="alexander.tedeschi@cornell.edu" \
      version="1.0" \
      description="Docker image for vol2birdR"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libproj-dev \
    gsl-bin \
    libgsl-dev \
    pkg-config \
    libbz2-dev \
    && rm -rf /var/lib/apt/lists/*  # Clean up cache to reduce image size

RUN \
  apt-get update && apt-get install --no-install-recommends -y \
  python3 python3-pip python3-setuptools && \
  rm -rf /var/lib/apt/lists/* && \
  pip3 install boto boto3 pytz astral


# Install the vol2birdR R package and MistNet
RUN R -e "install.packages('vol2birdR', dependencies=TRUE, repos='https://cloud.r-project.org/');" \
    && R -e "library(vol2birdR); install_mistnet(); mistnet_exists(); install_mistnet_model()"

    
# Install additional R packages
RUN R -e "install.packages(c('terra', 'raster', 'logger', 'lubridate', 'ggplot2', 'dplyr', \
                             'mgcv', 'units', 'gridExtra', 'pracma', 'parallel', \
                             'sf'), dependencies=TRUE, repos='https://cloud.r-project.org/')"

                             # Install the 'remotes' package to enable installation from GitHub
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org/')"

# Install 'bioRad' from GitHub
RUN R -e "remotes::install_github('adokter/bioRad')"

RUN R -e "install.packages('microbenchmark')"

RUN mkdir -p /opt/occult
RUN mkdir -p /opt/gis
RUN mkdir -p /opt/output
COPY occult /opt/occult
COPY gis  /opt/gis
COPY wgfd.R /opt/wgfd.R
COPY closest_records_KRIW.csv /opt/closest_records_KRIW.csv
COPY closest_records_KCYS.csv /opt/closest_records_KCYS.csv

WORKDIR /opt

ENTRYPOINT ["Rscript", "/opt/wgfd.R"]