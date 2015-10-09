# Use phusion/baseimage as base image. To make your builds reproducible, make
# sure you lock down to a specific version, not to `latest`!
# See https://github.com/phusion/baseimage-docker/blob/master/Changelog.md for
# a list of version numbers.

FROM phusion/baseimage
MAINTAINER Andrea Petri <apetri@phys.columbia.edu>

#Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

#This computer
ENV THIS="dock-lenstools"

#Refresh package cache
RUN apt-get update

########################################################
######Dependencies to install with apt-get##############
########################################################

#OpenMPI
RUN apt-get -y install libopenmpi-dev openmpi-bin

#Python, with headers and pip
RUN apt-get -y install python-dev python-pip 

#LAPACK and scipy
RUN apt-get -y install liblapack-dev python-scipy

#Clone the LensTools repository 
#RUN git clone https://github.com/apetri/LensTools

#Install requirements and LensTools
#RUN cd /LensTools ; git checkout docker ; pip install -r requirements.txt ; python setup.py install

#Point PYTHONPATH so that the python interpreter finds the scipy installation
#ENV PYTHONPATH=/usr/lib/python2.7/dist-packages:$PYTHONPATH

#Nice colors
#ENV PS1="\t \[$(tput bold)\]\[$(tput setaf 3)\][\[$(tput setaf 5)\]\u\[$(tput setaf 3)\]@\[$(tput setaf 5)\]\h \[$(tput setaf 2)\]\W\[$(tput setaf 3)\]]\\$ \[$(tput sgr0)\]"

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*