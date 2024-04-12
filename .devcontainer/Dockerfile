# R 4.0.0
FROM drclab/drupyter:5.7.2 as drubase
#--------------------------------
RUN apt install -y docker.io
#--------------------------------
USER drclab
#----------
CMD ["R"]

## set MRAN install date to 20th April 2020

# RUN echo "options(repos = c(CRAN='https://mran.microsoft.com/snapshot/2020-04-20'), download.file.method = 'libcurl')" >> /etc/R/Rprofile.site \
#   && apt-get update \
#   && apt-get install -y apt-utils \
#   && apt-get install -y libv8-dev libcurl4 libcurl4-openssl-dev libxml2-dev libxml2 \
#   && install2.r --ncpus 4 --deps NA StanHeaders rstan methods abind mvtnorm plyr functional reshape2 ggplot2

