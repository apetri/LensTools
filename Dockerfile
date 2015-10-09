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

#pkg-config
RUN apt-get -y install pkg-config

#git
RUN apt-get -y install git

#gfortran and OpenMPI
RUN apt-get -y install gfortran
RUN apt-get -y install libopenmpi-dev openmpi-bin

#Python, with headers and pip
RUN apt-get -y install python-dev python-pip 

#LAPACK and scipy
RUN apt-get -y install liblapack-dev python-scipy

#Clone the LensTools repository, install it along with its requirements 
RUN git clone https://github.com/apetri/LensTools
RUN cd /LensTools ; git checkout docker-ubuntu ; pip install -r requirements.txt ; python setup.py install

#Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*