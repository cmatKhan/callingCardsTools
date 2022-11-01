FROM rocker/r-ubuntu

RUN  apt-get update
RUN  apt-get install -y --no-install-recommends \
    software-properties-common \
    dirmngr \
    wget \
	build-essential \
    python3.9 \
	python3-pip

# Clean up
RUN apt-get autoremove -y

RUN pip install --upgrade pip

RUN wget https://github.com/cmatKhan/callingCardsTools/archive/refs/tags/v0.0.7.tar.gz

RUN pip install v0.0.7.tar.gz

RUN rm v0.0.7.tar.gz