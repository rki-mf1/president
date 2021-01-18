FROM continuumio/miniconda3

RUN pip install flake8 \
pytest \
pytest-cov \
pytest-flake8 \
pytest-pydocstyle \
pytest-benchmark \
pandas \
screed \
--quiet && \
conda install -c bioconda pblat -y
