c	program main
	integer function f77main(thI, tomegamI,
     .     tomegavI, tomegakI, tomegaQI,
     .     twQI, twQpI, z1, z2, D1, D2,
     .     tzmaxact, tzminact, tiwmodeI)
	implicit none
c   JMK:	include 'snlens.inc'
	include 'parameters.inc'
	double precision tzmaxact,tzminact,tomegamI,
     .		tomegavI,tomegaQI,twQI,twQpI,thI,tomegakI
	integer tiwmodeI
c   JMK:
	double precision zmaxact,zminact,omegamI,omegavI,omegaQI,
     .		wQI,wQpI,hI,omegakI 
	integer iwmodeI
	common /cfinsert/ zmaxact,zminact,omegamI,omegavI,omegaQI,
     .		wQI,wQpI,hI,omegakI,iwmodeI

	double precision D1,D2,Tgrowthfac,z1,z2
	
c JMK:
c	double precision omegamI,omegavI,omegaQI,wQI,wQpI,hI

	double precision omegampass,omegavpass,omegaQpass,
     .		wQpass,wQPpass,hpass,omegakpass
      	common /setparm5/ omegampass,omegavpass,omegaQpass,
     .		wQpass,wQPpass,hpass,omegakpass



	zmaxact=tzmaxact
	zminact=tzminact
	omegamI=tomegamI
	omegavI=tomegavI
	omegaQI=tomegaQI
	wQI=twQI
	wQpI=twQpI
	hI=thI
	omegakI=tomegakI
	iwmodeI=tiwmodeI
	

	hpass=hI
	omegampass=omegamI
	omegavpass=omegavI
	omegakpass=omegakI
	omegaQpass=omegaQI
	wQpass=wQI
	wQPpass=wQpI


c	write(*,*) 'input z'
c	read(*,*) z1
c	read(*,*) z2
	
	call growthini(1.0d0)

Cf2py intent(out) D1	
	D1=Tgrowthfac(z1,omegamI,omegavI,omegakI,omegaQI)
Cf2py intent(out) D2
	D2=Tgrowthfac(z2,omegamI,omegavI,omegakI,omegaQI)

c	write(*,*) D1,D2
c	write(*,*) (D1/D2),(D1/D2)**2

	f77main = 0
	end function f77main

c	stop
c	end 


	