      subroutine growthini(ak)
c Note that ak here is in unit of 1/Mpc.
c This is based on growthq.f that I copied from Marilena/ISW/Code on 11/19/07.
c See comments in growfk.f.
c The equation I want to solve is Dodelson eq. 7.73, except that I modify 
c the third term: (- 1.5 Omega_m H_0^2 / (a^5 H^2)) * f * delta
c The GR case is f=1.
c
c      program growthq
c this program computes the growth factor for a quintessence
c universe
      implicit none
c JMK:      include 'snlens.inc'
	  include 'parameters.inc'
c JMK:
	  double precision zmaxact,zminact,omegamI,omegavI,omegaQI,
     .		wQI,wQpI,hI,omegakI 
	  integer iwmodeI
	  common /cfinsert/ zmaxact,zminact,omegamI,omegavI,omegaQI,
     .		wQI,wQpI,hI,omegakI,iwmodeI

      integer i
      double precision zmin,zmax,amax,amin,da,aa
      double precision g0,z
      integer neq,ind,ier
      parameter (neq=1)
      double precision y(neq),c(24),work(neq,9)
      double precision tol
      external deriv
      double precision astart,aend
      double precision g(na)
      double precision fg,zaa
      common /growthpass/ zaa(na),fg(na)

      integer igrowthinimode
      parameter (igrowthinimode=2)
c igrowthinimode=1 means use an old approximation that I think I got from a paper by Wang, Steinhardt and company.
c igrowthinimode=2 means integrating the exact (small wavelength limit) equation (Dodelson eq. 7.73).
      integer neqp
      parameter (neqp=2)
      double precision yp(neqp),cp(24),workp(neqp,9)
      external derivDodelson

      double precision ak,akpass
      common /akpasspass/ akpass

      integer igrid
      parameter (igrid=2)
c igrid=1 means using a uniform grid in a.
c igrid=2 means using a uniform grid in z.
c Should use igrid=2 if zminact is close to -1.
      double precision dz

      akpass=ak

c      zmin=0.0d0
c      zmax=zsI

      zmin=zminact
c need to set zmin=0.0, because in universal.f, 
c growthfac evaluated at z=0.0 is actually required.

      zmax=zmaxact

c      zmax=1.1d0
c      na=1000

      if (igrid.eq.1) then
      	amax=1.0d0/(1.0d0+zmin)
        amin=1.0d0/(1.0d0+zmax)
        da=(amax-amin)/(na-1)
      elseif (igrid.eq.2) then
	dz=(zmax-zmin)/(na-1)
      endif

c      open(unit=1,file='growthq.dat',form='formatted')
c      rewind 1

      do i=1,na
         if (igrid.eq.1) then
         	aa=amin+da*(i-1)
	 elseif (igrid.eq.2) then
		aa=zmax-dz*(i-1)
		aa=1.0d0/(1.0d0+aa)
	 endif

         tol=1.0d-12

c It seems when wQpI is positive and large, I need to set astart to be larger
c in order for dverk not to have problems. Same if wQI is significantly larger than -1.
c So, I will just set astart=0.01d0 which seems to avoid problems for most wQI and wQpI 
c of interests. And I have also checked that in cases where both small and larger astart's
c work, using either astart doesn't really change the resulting growth factor very much.
c         astart=0.0001d0
         astart=0.01d0
         aend=aa

         ier=0
         ind=1

	 if (igrowthinimode.eq.1) then

           y(1)=astart

           call dverk(neq,deriv,astart,y,aend,tol,ind,c
     .      ,neq,work,ier)

           if (ind.lt.0.or.ier.gt.0) then
              write(*,*) 'dverk error, ind, ier=',ind,ier
           end if

	   g(i)=y(1)

           if (i.eq.na) then
              g0=y(1)
           endif

	 elseif (igrowthinimode.eq.2) then

	   yp(1)=astart
	   yp(2)=1.0d0	
c yp(1)=delta, yp(2)=d(delta)/da

	   call dverk(neqp,derivDodelson,astart,yp,aend,tol,ind,cp
     .      ,neqp,workp,ier)

           if (ind.lt.0.or.ier.gt.0) then
              write(*,*) 'dverk error, ind, ier=',ind,ier
           end if

	   g(i)=yp(1)

           if (i.eq.na) then
              g0=yp(1)
           endif

	 endif
         
      enddo

      do i=1,na
         fg(i)=g(na-i+1)
      enddo

      do i=1,na
	 if (igrid.eq.1) then
         	aa=amax-da*(i-1)
         	z=1.0d0/aa-1.0d0
	 elseif (igrid.eq.2) then
		z=zmin+dz*(i-1)
		aa=1.0d0/(1.0d0+z)
	 endif
         g(i)=fg(i)
         fg(i)=g(i)/aa
         zaa(i)=z
c         write(1,*) z,g(i)/aa,g(i)/g0
c second column is what PD called g(omega)
c and third column is the growth factor normalized
c to 1 at z=0
      enddo

c      close(1)
      
      return
c      stop
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine deriv(n,x,y,dydx)
      implicit none
      integer n
      double precision x,y(n),dydx(n),growthfacQ

      dydx(1)=growthfacQ(x)*y(1)/x

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine derivDodelson(n,x,y,dydx)
      	implicit none
      	integer n
      	double precision x,y(n),dydx(n)
	double precision Dodfac1,Dodfac2

	dydx(1)=y(2)
	dydx(2)=Dodfac1(x)*y(1)-Dodfac2(x)*y(2)

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function Dodfac1(a)
	implicit none
c JMK:	include 'snlens.inc'
	include 'parameters.inc'
c JMK:
	double precision zmaxact,zminact,omegamI,omegavI,omegaQI,
     .		wQI,wQpI,hI,omegakI
	integer iwmodeI
	common /cfinsert/ zmaxact,zminact,omegamI,omegavI,omegaQI,
     .		wQI,wQpI,hI,omegakI,iwmodeI
	
	double precision a,Dodfac1
	double precision omegampass,omegavpass,omegaQpass,
     .		wQpass,wQPpass,hpass,omegakpass
      	common /setparm5/ omegampass,omegavpass,omegaQpass,
     .		wQpass,wQPpass,hpass,omegakpass
c This gives (3 Omegam H0^2)/(2 a^5 H^2)
	double precision z,hzz
	integer iwmode
	double precision omegam,omegav,omegaQ,
     .		wQ,wQp,h,omegak

	double precision akpass
	common /akpasspass/ akpass
	double precision ak,fmodfact,fmod

	ak=akpass
	fmodfact=fmod(ak,a)

	h=hpass
	omegam=omegampass
	omegav=omegavpass
	omegak=omegakpass
	omegaQ=omegaQpass
	wQ=wQpass
	wQp=wQPpass

	iwmode=iwmodeI

	z=1.0d0/a - 1

	call HzzQ(z,hzz,h,
     .            omegam,omegav,omegak,
     .            omegaQ,wQ,wQp,iwmode)

	Dodfac1=1.5d0*omegam/a**5/(hzz/100.0d0/h)**2
     .		*fmodfact

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function Dodfac2(a)
	implicit none
c JMK:	include 'snlens.inc'
	include 'parameters.inc'
c JMK:
	double precision zmaxact,zminact,omegamI,omegavI,omegaQI,
     .		wQI,wQpI,hI,omegakI
	integer iwmodeI
	common /cfinsert/ zmaxact,zminact,omegamI,omegavI,omegaQI,
     .		wQI,wQpI,hI,omegakI,iwmodeI

	double precision a,Dodfac2
	double precision omegampass,omegavpass,omegaQpass,
     .		wQpass,wQPpass,hpass,omegakpass
      	common /setparm5/ omegampass,omegavpass,omegaQpass,
     .		wQpass,wQPpass,hpass,omegakpass
c This gives (1/H)*(dH/da) + 3/a
	double precision z,hzz,dhda
	integer iwmode
	double precision omegam,omegav,omegaQ,
     .		wQ,wQp,h,omegak

	h=hpass
	omegam=omegampass
	omegav=omegavpass
	omegak=omegakpass
	omegaQ=omegaQpass
	wQ=wQpass
	wQp=wQPpass

	iwmode=iwmodeI

	z=1.0d0/a - 1

	call HzzQ(z,hzz,h,
     .            omegam,omegav,omegak,
     .            omegaQ,wQ,wQp,iwmode)

	call dHdaQ(z,dhda,h,
     .            omegam,omegav,omegak,
     .            omegaQ,wQ,wQp,iwmode)

	Dodfac2=dhda/hzz + 3.0d0/a

	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function growthfacQ(a)
      implicit none
      double precision alpha,a,growthfacQ,omegamAA
      double precision omegampass,omegavpass,omegaQpass,
     .		wQpass,wQPpass,hpass,omegakpass
      common /setparm5/ omegampass,omegavpass,omegaQpass,
     .		wQpass,wQPpass,hpass,omegakpass

      if (abs(wQPpass).gt.1.0d-12) then
	write(*,*) 'this code does not work for nonzero wp'
	stop
      endif

      alpha=3.0d0/(5.0d0-wQpass/(1.0d0-wQpass)) + 3.0d0/125.0d0 *
     .      (1.0d0 - wQpass)*(1.0d0 - 3.0d0*wQpass/2.0d0)/
     .      (1.0d0 - 6.0d0*wQpass/5.0d0)**3 * (1.0d0 - omegampass)
      growthfacQ=omegamAA(a)**alpha
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function omegamAA(a)
      implicit none
      double precision a,omegamAA
      double precision omegampass,omegavpass,omegaQpass,
     .		wQpass,wQPpass,hpass,omegakpass
      common /setparm5/ omegampass,omegavpass,omegaQpass,
     .		wQpass,wQPpass,hpass,omegakpass

      omegamAA=omegampass/(omegampass+omegakpass*a+omegavpass*a**3
     .         +omegaqpass/a**(3.0d0*wQpass))

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function fgQQ(z)
      implicit none
c JMK:     include 'snlens.inc'
	  include 'parameters.inc'
c JMK:
      double precision zmaxact,zminact,omegamI,omegavI,omegaQI,
     .		wQI,wQpI,hI,omegakI
      integer iwmodeI
	  common /cfinsert/ zmaxact,zminact,omegamI,omegavI,omegaQI,
     .		wQI,wQpI,hI,omegakI,iwmodeI

      double precision fgQQ,z
      double precision fg,zaa
      common /growthpass/ zaa(na),fg(na)
      integer iok,i
      double precision zatol,wL,wR
      
      zatol=1.0d-10

      iok=0

      do i=1,na-1
         if (z.ge.zaa(i).and.z.le.zaa(i+1)) then
            if ((zaa(i+1)-zaa(i)).gt.zatol) then
               wL=(z-zaa(i))/(zaa(i+1)-zaa(i))
               wR=(zaa(i+1)-z)/(zaa(i+1)-zaa(i))
            else
               wL=0.0
               wR=1.0
            endif
            fgQQ=fg(i)*wR+fg(i+1)*wL
            iok=1
            goto 10
         endif
      enddo 
      
      if (abs(z-zaa(na)).lt.1.0d-5) then
         fgQQ=fg(na)
         iok=1
         goto 10
      endif

      if (iok.eq.0) then
         write(*,*) 'z outside range of power table'
         write(*,*) 'z asked for',z
         write(*,*) 'largest table z',zaa(na)
         write(*,*) 'smallest table z',zaa(1)
         if (z.lt.zaa(1)) then
            fgQQ=fg(1)
         elseif (z.gt.zaa(na)) then
            fgQQ=fg(na)
         endif
      endif
      
 10   continue

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function fmod(ak,a)
      implicit none
c JMK:     include 'snlens.inc'
	  include 'parameters.inc'
c JMK:
      double precision zmaxact,zminact,omegamI,omegavI,omegaQI,
     .		wQI,wQpI,hI,omegakI
      integer iwmodeI
      common /cfinsert/ zmaxact,zminact,omegamI,omegavI,omegaQI,
     .		wQI,wQpI,hI,omegakI,iwmodeI

      double precision a,ak,fmod
      integer imode
      parameter (imode=5)
c imode determines the choice of the function f(ak,a)
c For imode=1: basically White and Kochanek
      double precision alpha
      parameter (alpha=1.0d0)
      double precision amass2
      parameter (amass2=(0.05d0*0.7d0)**2)
c      parameter (amass2=(0.0001d0*0.7d0)**2)
c According to White and Kochanek, alpha=1 is allowed as long as amass < 0.1 h/Mpc.
c According to Verde et al., alpha_Verde = - alpha_White, and
c              lambda_Verde = 1/amass.
c	       From Verde et al. Fig 2, looks like lambda_Verde > 20 is OK,
c	       which means amass < 1/20 ~ 0.05 h/Mpc.
c	       Let's use that!
c 
c For imode=2: basically an expansion in (aH/k) like in Wagoner et al.
      double precision Afac,Bfac
      parameter (Afac=-3.0d0)
      parameter (Bfac=0.0d0)
      double precision ratio,hzz,z
c For imode=3, basically same as imode=1, except make amass propto a.
c For imode=4, basically same as imode=2, except make Afac propto a, Bfac propto a^2.
c For imode=5, same as GR.

      double precision omegampass,omegavpass,omegaQpass,
     .		wQpass,wQPpass,hpass,omegakpass
      common /setparm5/ omegampass,omegavpass,omegaQpass,
     .		wQpass,wQPpass,hpass,omegakpass
      integer iwmode
      double precision omegam,omegav,omegaQ,
     .		wQ,wQp,h,omegak

      if (imode.eq.1) then
          fmod=(1.0d0-alpha) + alpha*ak**2/
     .			(ak**2 + amass2)
      elseif (imode.eq.2) then
	  h=hpass
      	  omegam=omegampass
          omegav=omegavpass
          omegak=omegakpass
          omegaQ=omegaQpass
          wQ=wQpass
          wQp=wQPpass
          iwmode=iwmodeI

	  z=1.0d0/a - 1.0d0
          call HzzQ(z,hzz,h,
     .            omegam,omegav,omegak,
     .            omegaQ,wQ,wQp,iwmode)
          ratio=a*hzz/ak/cspeed
	  fmod=1.0d0 + Afac*ratio + Bfac*ratio**2
      elseif (imode.eq.3) then
	  fmod=(1.0d0-alpha) + alpha*ak**2/
     .			(ak**2 + amass2*a**2)
      elseif (imode.eq.4) then
	  h=hpass
      	  omegam=omegampass
          omegav=omegavpass
          omegak=omegakpass
          omegaQ=omegaQpass
          wQ=wQpass
          wQp=wQPpass
          iwmode=iwmodeI

	  z=1.0d0/a - 1.0d0
          call HzzQ(z,hzz,h,
     .            omegam,omegav,omegak,
     .            omegaQ,wQ,wQp,iwmode)
          ratio=a*hzz/ak/cspeed
	  fmod=1.0d0 + Afac*ratio*a + Bfac*ratio**2*a*a
      elseif (imode.eq.5) then
	  fmod=1.0d0
      endif

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
