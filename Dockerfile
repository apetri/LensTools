FROM python:2
MAINTAINER Andrea Petri <apetri@phys.columbia.edu>

#Install liblapack-dev
RUN apt-get update
RUN apt-get -y install liblapack-dev

#Install numpy,scipy,matplotlib
RUN apt-get -y install python-numpy
RUN apt-get -y install python-scipy
RUN apt-get -y install python-matplotlib

#Clone the LensTools repository 
RUN git clone https://github.com/apetri/LensTools

#Install requirements and LensTools
RUN cd /LensTools ; git checkout docker ; pip install -r requirements.txt ; python setup.py install
