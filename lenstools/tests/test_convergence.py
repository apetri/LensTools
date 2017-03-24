import os

from .. import ConvergenceMap

from .. import dataExtern

import numpy as np
from astropy.units import deg,rad

import matplotlib.pyplot as plt


test_map = ConvergenceMap.load(os.path.join(dataExtern(),"conv.fit"))

#Set bin edges
l_edges = np.arange(200.0,50000.0,200.0)
thresholds_mf = np.arange(-2.0,2.0,0.2)
thresholds_pk = np.arange(-1.0,5.0,0.2)

def test_visualize():

	assert test_map.data.dtype == np.float

	test_map.setAngularUnits(deg)
	test_map.visualize()
	test_map.savefig("map.png")
	test_map.setAngularUnits(deg)

def test_save():
	test_map.save("conv_save.fits")

def test_power():

	#Compute
	l,Pl = test_map.powerSpectrum(l_edges)
	assert type(l)==np.ndarray
	assert type(Pl)==np.ndarray

	#Visualize
	fig,ax = plt.subplots()
	ax.plot(l,l*(l+1)*Pl/(2.0*np.pi))
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$l$")
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$")

	plt.savefig("power_spectrum.png")
	plt.clf()

def test_cross():

	#Load
	conv1 = ConvergenceMap.load(os.path.join(dataExtern(),"conv1.fit"))
	conv2 = ConvergenceMap.load(os.path.join(dataExtern(),"conv2.fit"))

	#Cross
	l,Pl = conv1.cross(conv2,l_edges=l_edges)

	#Visualize
	fig,ax = plt.subplots()
	ax.plot(l,np.abs(l*(l+1)*Pl/(2.0*np.pi)))
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$l$")
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$")

	plt.savefig("cross_spectrum.png")
	plt.clf()

def test_bispectrum():

	#Measure bispectrum
	l,be = test_map.bispectrum(l_edges,configuration="equilateral",scale=lambda x:(x**4)/(2.*np.pi)**2)
	l,bf = test_map.bispectrum(l_edges,configuration="folded",ratio=0.5,scale=lambda x:(x**4)/(2.*np.pi)**2)

	#Plot
	fig,ax = plt.subplots()
	ax.plot(l,np.abs(be),label=r"$\ell_1=\ell_2=\ell_3=\ell$")
	ax.plot(l,np.abs(bf),label=r"$\ell_1=\ell,\ell_2=\ell_3=\ell/2$")
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$\ell$")
	ax.set_ylabel(r"$(\ell_1\ell_2\ell_3)^{4/3}B(\ell_1,\ell_2,\ell_3)/(2\pi)^2$")
	ax.legend()

	plt.savefig("kappa_bispectrum.png")
	plt.clf()


def test_pdf():

	#Compute
	v,p = test_map.pdf(thresholds_mf,norm=True)

	#Visualize
	fig,ax = plt.subplots()
	ax.plot(v,p)
	ax.set_xlabel(r"$\nu=\kappa/\sigma$")
	ax.set_ylabel(r"$P(\nu)$")

	plt.savefig("pdf.png")
	plt.clf()



def test_minkowski():

	#Compute
	nu,V0,V1,V2 = test_map.minkowskiFunctionals(thresholds_mf,norm=True)

	#Assert computation went OK
	assert hasattr(test_map,"gradient_x")
	assert hasattr(test_map,"gradient_y")
	assert hasattr(test_map,"hessian_xx")
	assert hasattr(test_map,"hessian_yy")
	assert hasattr(test_map,"hessian_xy")

	#Visualize
	fig,ax = plt.subplots(1,3,figsize=(24,8))
	ax[0].plot(nu,V0)
	ax[1].plot(nu,V1)
	ax[2].plot(nu,V2)

	ax[0].set_xlabel(r"$\nu=\kappa/\sigma$")
	ax[0].set_ylabel(r"$V_0(\nu)$")

	ax[1].set_xlabel(r"$\nu=\kappa/\sigma$")
	ax[1].set_ylabel(r"$V_1(\nu)$")

	ax[2].set_xlabel(r"$\nu=\kappa/\sigma$")
	ax[2].set_ylabel(r"$V_2(\nu)$")

	fig.tight_layout()

	plt.savefig("minkowski.png")
	plt.clf()

def test_peaks():

	#Compute
	nu,pk = test_map.peakCount(thresholds_pk,norm=True)

	#Check if computation went OK
	assert type(nu)==np.ndarray
	assert type(pk)==np.ndarray

	#Visualize
	fig,ax = plt.subplots()
	ax.plot(nu,pk)
	ax.set_xlabel(r"$\nu=\kappa/\sigma$")
	ax.set_ylabel(r"$dN/d\nu$")

	plt.savefig("peaks.png")


def test_peak_locations():

	#Thresholds for high peaks
	high_thresholds = np.arange(0.3,0.6,0.01)

	#Find the peak locations
	values,locations = test_map.locatePeaks(high_thresholds)

	#Visualize the result
	fig,ax = plt.subplots(1,2,figsize=(16,8))
	test_map.visualize(fig=fig,ax=ax[0],colorbar=True)
	test_map.visualize(fig=fig,ax=ax[1])
	ax[1].scatter(locations[:,0].value,locations[:,1].value,color="black")
	ax[1].set_xlim(0.0,test_map.side_angle.value)
	ax[1].set_ylim(0.0,test_map.side_angle.value)

	#And save it
	fig.tight_layout()
	fig.savefig("peak_locations.png")


def test_getValues():

	b = np.linspace(0.0,test_map.side_angle.value,test_map.data.shape[0])
	xx,yy = np.meshgrid(b,b) * deg

	new_values = test_map.getValues(xx,yy)
	assert (new_values==test_map.data)[:-1,:-1].all()

def test_gradient_partial():

	b = np.linspace(0.0,test_map.side_angle.value,test_map.data.shape[0])
	xx,yy = np.meshgrid(b,b) * deg

	gx,gy = test_map.gradient()
	gxp,gyp = test_map.gradient(xx,yy)

	assert (gx==gxp)[:-1,:-1].all()
	assert (gy==gyp)[:-1,:-1].all()


def test_hessian_partial():

	b = np.linspace(0.0,test_map.side_angle.value,test_map.data.shape[0])
	xx,yy = np.meshgrid(b,b) * deg

	hxx,hyy,hxy = test_map.hessian()
	hxxp,hyyp,hxyp = test_map.hessian(xx,yy)

	assert (hxx==hxxp)[:-1,:-1].all()
	assert (hyy==hyyp)[:-1,:-1].all()
	assert (hxy==hxyp)[:-1,:-1].all()


def test_cut():

	b = np.array([0.0,test_map.side_angle.value/2,0.0,test_map.side_angle.value/2]) * deg
	cut_map = test_map.cutRegion(b)
	cut_map.visualize()
	cut_map.savefig("map_cut.png")


def test_translate():

	b = np.array([0.5,test_map.side_angle.value+0.5,0.0,test_map.side_angle.value]) * deg
	translated_map = test_map.cutRegion(b)
	translated_map.visualize()
	translated_map.savefig("map_translated.png")




