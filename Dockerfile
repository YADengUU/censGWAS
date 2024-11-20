FROM rocker/r-base:4.2.2
MAINTAINER Yaqi Deng <yaqi.deng@igp.uu.se>

RUN apt update && apt-get install -y\
    libcurl4-openssl-dev \
    libxml2-dev
RUN R -e "install.packages('remotes')"
RUN R -e "install.packages('broom')"
RUN R -e "install.packages('censReg')"
RUN R -e "install.packages('RNOmni')"

RUN R -e "remotes::install_github('YADengUU/censGWAS')"

RUN plink2 --version || echo "PLINK not found, please install PLINK on your system."

ENTRYPOINT ["R"]
CMD ["--help"]
