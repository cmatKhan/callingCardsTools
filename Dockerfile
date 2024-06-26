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

#RUN wget https://github.com/cmatKhan/callingCardsTools/archive/refs/tags/v0.0.10.tar.gz

#RUN pip install v0.0.10.tar.gz

#RUN rm v0.0.10.tar.gz

#COPY dist/callingCardsTools-0.1.0.tar.gz .

RUN pip install callingcardstools==1.1.0