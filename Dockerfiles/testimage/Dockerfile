FROM continuumio/miniconda3

RUN pip install flake8 \
pytest \
pytest-cov \
pytest-flake8 \
pytest-pydocstyle \
pandas \
screed \
lxml \
mock \
--quiet && \
conda install -c bioconda pblat diamond -y
