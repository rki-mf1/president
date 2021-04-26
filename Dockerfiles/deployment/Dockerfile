FROM continuumio/miniconda3:4.8.2
COPY ./ ./
RUN apt update && apt install -y procps wget gzip && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN conda config --add channels conda-forge && \
		conda config --add channels bioconda && \
		conda config --add channels default
RUN conda install -y -c bioconda python=3.8 pblat diamond pandas screed && \
conda clean -a && \
python setup.py install
