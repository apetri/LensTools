FROM python:2
MAINTAINER Andrea Petri <apetri@phys.columbia.edu>

#Install liblapack-dev
RUN apt-get update
RUN apt-get -y install liblapack-dev

#Install requirements
RUN pip install -r requirements.txt

#Install lenstools
RUN python setup.py install
