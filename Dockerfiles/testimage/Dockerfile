FROM continuumio/miniconda3:latest

RUN pip install \
'flake8<=3.8.4' \
pytest \
pytest-cov \
pytest-flake8 \
pytest-pydocstyle \
pandas \
screed \
--quiet && \
conda install -c bioconda pblat -y
