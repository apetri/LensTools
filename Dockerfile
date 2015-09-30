FROM python:2
MAINTAINER Andrea Petri <apetri@phys.columbia.edu>

#Install liblapack-dev
RUN apt-get update
RUN apt-get -y install liblapack-dev

#Install scipy with apt-get
RUN apt-get -y install python-scipy

#Clone the LensTools repository 
RUN git clone https://github.com/apetri/LensTools

#Install requirements and LensTools
RUN cd /LensTools ; git checkout docker ; pip install -r requirements.txt ; python setup.py install

#Point PYTHONPATH so that the python interpreter finds the scipy installation
ENV PYTHONPATH=/usr/lib/python2.7/dist-packages:$PYTHONPATH

#Nice colors
ENV PS1="\t \[$(tput bold)\]\[$(tput setaf 3)\][\[$(tput setaf 5)\]\u\[$(tput setaf 3)\]@\[$(tput setaf 5)\]\h \[$(tput setaf 2)\]\W\[$(tput setaf 3)\]]\\$ \[$(tput sgr0)\]"