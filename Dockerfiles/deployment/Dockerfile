FROM continuumio/miniconda3:latest

COPY ./ ./

RUN until apt-get update; \
    do \
      echo "== apt-get update failed, retrying =="; \
      sleep 5; \
    done \
    && until apt-get install -y procps wget gzip; \
    do \
      echo "== apt-get install failed, retrying =="; \
      sleep 5; \
    done \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda config --add channels conda-forge && \
		conda config --add channels bioconda && \
		conda config --add channels default

RUN conda install -y -c bioconda python=3.8 pblat pandas screed && \
conda clean -a && \
python setup.py install
