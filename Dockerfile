FROM rocker/r-ubuntu

RUN  apt-get update
RUN  apt-get install -y --no-install-recommends \
    software-properties-common \
    dirmngr \
    wget \
	build-essential \
    python \
	python-setuptools \
	python-dev \
	python-pip

# Clean up
RUN apt-get autoremove -y && autoclean -y