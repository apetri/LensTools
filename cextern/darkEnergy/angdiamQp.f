c This is based on angdiamQ.f, which assumed wQ=constant.
c Here, I want to allow for a time-varying wQ.
c What does that mean?
c It means that p = w rho.
c Recall that d(rho * a^3) = - p d(a^3)
c i.e. d(rho a^3) = - w*rho*(3 a^2 da)
c i.e. 3*rho*a^2*da + a^3*drho = -3*w*rho*a^2*da
c i.e. drho/rho = (-3*w*a^2*da - 3*a^2*da)/a^3
c               = -3*(w+1)*da/a
c
c if iwmode=1, w=constant, we have dln(rho) = dln(a^{-3(1+w)})
c              i.e. rho ~ a^-3(1+w)
c
c if iwmode=2, w=wQ+wQp*z, we have w=wQ+wQp*(1/a - 1)
c	       w=(wQ-wQp) + wQp/a
c 	       therefore, dln(rho) = -3*[(1+wQ-wQp) + wQp/a]*da/a
c			           = -3*(1+wQ-wQp)dlna - 3*wQp*da/a^2
c				   = dln(a^-3*(1+wQ-wQp)) + d(3*wQp/a)
c			  dln(rho) = d(ln(a^-3*(1+wQ-wQp)) + 3*wQp/a)
c			  ln(rho) = ln(a^-3*(1+wQ-wQp)) + 3*wQp/a + constant
c			  rho ~ a^-3*(1+wQ-wQp) * exp(3*wQp/a)
c                             ~ (1+z)^3(1+wQ-wQp) * exp(3*wQp*z)   [the factor of exp(3*wQp) is just an overall constant]
c 		the above agrees with expression given in Kim's 0304509.
c		it also has the nice property that rho(z=0) = 1.
c
c if iwmode=3, w=wQ+wQp*(1-a) = wQ+wQp - wQp*a
c	       therefore, dln(rho) = -3*[(1+wQ+wQp) - wQp*a]*da/a
c	       i.e. dln(rho) = d(lna^-3(1+wQ+wQp)) + d(3*wQp*a)
c		    ln(rho) = ln(a^-3(1+wQ+wQp)) + 3*wQp*a
c		    rho ~ a^-3(1+wQ+wQp) * exp(3*wQp*a)
c		        ~ (1+z)^3(1+wQ+wQp) * exp(3*wQp/(1+z)) * exp(-3*wQp)  [extra factor of exp(-3*wQp) is there to
c									       make rho(z=0) = 1]
c			~ (1+z)^3(1+wQ+wQp) * exp(-3*wQp*z/(1+z))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine HzzQ(z,hz,h,omegam,omegav,omegak,
     .               omegaQ,wQ,wQp,iwmode)
c gives Hubble constant at z
      implicit none
      double precision z,hz,h,omegam,omegav,omegak,h00
      double precision omegaQ,wQ

      double precision wQp
      integer iwmode

      h00=100.0d0*h

      if (abs(omegam+omegak+omegav+omegaQ-1).gt.1.0d-5) 
     .     then
         write(*,*) 'Omegas do not sum to one'
         write(*,*) 'error in HzzQ'
         write(*,*) omegam+omegak+omegav+omegaQ
         stop
      endif

      if (iwmode.eq.1) then
         hz=h00*(omegam*(1.0d0+z)**3 + omegak*(1.0d0+z)**2
     .        + omegav + omegaQ*(1.0d0+z)**(3.0d0*(1.0d0+wQ))
     .        )**0.5d0
      elseif (iwmode.eq.2) then
	hz=h00*(omegam*(1.0d0+z)**3 + omegak*(1.0d0+z)**2
     .        + omegav + omegaQ*( 
     .	      (1.0d0+z)**(3.0d0*(1.0d0+wQ-wQp))*
     .	      exp(3.0d0*wQp*z))
     .        )**0.5d0
      elseif (iwmode.eq.3) then
	hz=h00*(omegam*(1.0d0+z)**3 + omegak*(1.0d0+z)**2
     .        + omegav + omegaQ*( 
     .	      (1.0d0+z)**(3.0d0*(1.0d0+wQ+wQp))*
     .	      exp(-3.0d0*wQp*z/(1.0d0+z)))
     .        )**0.5d0
      endif

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chiRQ(z,chi,h,omegam,omegav,omegak,
     .              omegaQ,wQ,wQp,iwmode)
c integrates chi directly
c see comments in rrchiQ
c
c this gives the radial comoving distance to an object at z
c where z is input in rrchi

c chi in unit of Mpc, if want h^-1 Mpc, then
c take chi and multiply by h.

      implicit none
      double precision z,chi,h,omegam,omegav,omegak
      double precision omegaQ,wQ
      double precision h00,c00,chm00

      double precision wQp
      integer iwmode

      double precision omegampass,omegavpass,
     .       omegakpass,omegaQpass,wQpass,wQppass
      common /derivsQpass/ omegampass,omegavpass,
     .       omegakpass,omegaQpass,wQpass,wQppass
      integer iwmodepass
      common /derivsQpassint/ iwmodepass

      integer neq,ind,ier
      parameter (neq=1)
      double precision y(neq),c(24),work(neq,9)
      double precision tol
      external derivsQ
      double precision xstart,xend


      if (abs(omegam+omegak+omegav+omegaQ-1).gt.1.0d-5) 
     .     then
         write(*,*) 'Omegas do not sum to one'
         write(*,*) 'error in chiRQ'
         stop
      endif

      omegampass=omegam
      omegakpass=omegak
      omegavpass=omegav
      omegaQpass=omegaQ
      wQpass=wQ
      wQppass=wQp
      iwmodepass=iwmode

      h00=100.0d0*h
      c00=2.99792458d5

      chm00=c00/h00

      tol=1.0d-12

      xstart=0.0d0
      xend=z
      ier=0
      ind=1
      y(1)=0.0
      chi=0.0
      call dverk(neq,derivsQ,xstart,y,xend,tol,ind,c
     .      ,neq,work,ier)
      if (ind.lt.0.or.ier.gt.0) then
         write(*,*) 'dverk error, ind, ier=',ind,ier
      end if

      chi=y(1)*chm00

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rrchiQ(z,rchi,h,omegam,omegav,omegak,
     .                 omegaQ,wQ,wQp,iwmode)
c this is an alternative version of rrchi in angdiam.f
c here, I am integrating directly
c
c this gives the comoving coordinate distance rchi for a given z
c such that at angle dtheta, the object at z has a comoving
c size of dtheta * rchi

c Note that rchi is not the radial comoving distance to the
c object, the two are equal only for flat universes
c
c assumes a=1 today
c rchi in unit of Mpc, if want h^-1 Mpc, then
c take rchi and multiply by h.
      

c the strategy is to do this in two steps:
c let ds^2 = -dt^2 + a^2 (t) (d chi^2 + rchi^2 d theta^2)
c where 
c rchi = chi  for flat universe 
c   (sqrt{omegak}/(c H0^-1))^-1 sinh[(sqrt{omegak}/(c H0^-1)) chi] 
c   for omegak > 0
c   (sqrt{-omegak}/(c H0^-1))^-1 sin[(sqrt{-omegak}/(c H0^-1)) chi] 
c   for omegak < 0
c
c Then compute chi first using
c int dchi = c H0^-1 \int_a(z)^1 (da / [a^2 (H/H0)]) 
c  or = c H0^-1 \int_0^z dz / (H/H0) 
c and then evaluate rchi using the above.

      implicit none
      double precision z,rchi,h,omegam,omegav,omegak,
     .                 omegaQ,wQ,wQp
      integer iwmode
      double precision chi,h00,c00,chm00,kappa,kappa2

      if (abs(omegam+omegak+omegav+omegaQ-1).gt.1.0d-5) 
     .     then
         write(*,*) 'Omegas do not sum to one'
         write(*,*) 'error in rrchiQ'
         stop
      endif

      call chiRQ(z,chi,h,omegam,omegav,omegak,omegaQ,wQ,
     .	wQp,iwmode)

      h00=100.0d0*h
      c00=2.99792458d5

      chm00=c00/h00

      kappa=-omegak/chm00/chm00
      
      
      if (omegak.eq.0) then
c flat universe
         rchi=chi
      elseif (omegak.gt.0.0) then
c open universe
         kappa2=(-kappa)**0.5
         rchi=sinh(kappa2*chi)/kappa2
      elseif (omegak.lt.0.0) then
         kappa2=(kappa)**0.5
         rchi=sin(kappa2*chi)/kappa2
      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine derivsQ(n,x,y,dydx)
      implicit none
      integer n
      double precision x,y(n),dydx(n)

      double precision omegampass,omegavpass,
     .       omegakpass,omegaQpass,wQpass,wQppass
      common /derivsQpass/ omegampass,omegavpass,
     .       omegakpass,omegaQpass,wQpass,wQppass
      integer iwmodepass
      common /derivsQpassint/ iwmodepass

      double precision hz,h

      h=1.0d0

      call HzzQ(x,hz,h,omegampass,omegavpass,omegakpass
     .        ,omegaQpass,wQpass,wQppass,iwmodepass)

      hz=hz/100.0d0

      dydx(1)=1.0d0/hz

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rrchikappa(chi,h,omegak,rchi)
c gives rchi for an input chi
      implicit none
      double precision chi,h,omegak,rchi
      double precision h00,c00,chm00,kappa,kappa2

      h00=100.0d0*h
      c00=2.99792458d5

      chm00=c00/h00

      kappa=-omegak/chm00/chm00


      if (omegak.eq.0) then
c flat universe
         rchi=chi
      elseif (omegak.gt.0.0) then
c open universe
         kappa2=(-kappa)**0.5
         rchi=sinh(kappa2*chi)/kappa2
      elseif (omegak.lt.0.0) then
         kappa2=(kappa)**0.5
         rchi=sin(kappa2*chi)/kappa2
      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dchidwQcompute(z,dchidwQ,
     .     h,omegam,omegav,omegak,
     .              omegaQ,wQ,wQp,iwmode)

      implicit none
      double precision z,dchidwQ
      double precision h,omegam,omegav,omegak
      double precision omegaQ,wQ,wQp
      integer iwmode
      double precision h00,c00,chm00

      double precision omegampass,omegavpass,
     .       omegakpass,omegaQpass,wQpass,wQppass
      common /derivsdchidwQpass/ omegampass,omegavpass,
     .       omegakpass,omegaQpass,wQpass,wQppass
      integer iwmodepass
      common /derivsdchidwQpassint/ iwmodepass

      integer neq,ind,ier
      parameter (neq=1)
      double precision y(neq),c(24),work(neq,9)
      double precision tol
      external derivsdchidwQ
      double precision xstart,xend

      omegampass=omegam
      omegakpass=omegak
      omegavpass=omegav
      omegaQpass=omegaQ
      wQpass=wQ
      wQppass=wQp
      iwmodepass=iwmode

      h00=100.0d0*h
      c00=2.99792458d5

      chm00=c00/h00

      tol=1.0d-12

      xstart=0.0d0
      xend=z
      ier=0
      ind=1
      y(1)=0.0
      dchidwQ=0.0
      call dverk(neq,derivsdchidwQ,
     .     xstart,y,xend,tol,ind,c
     .      ,neq,work,ier)
      if (ind.lt.0.or.ier.gt.0) then
         write(*,*) 'dverk error, ind, ier=',ind,ier
      end if

      dchidwQ=y(1)*(-chm00)/2.0

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine derivsdchidwQ(n,x,y,dydx)
      implicit none
      integer n
      double precision x,y(n),dydx(n)

      double precision omegampass,omegavpass,
     .       omegakpass,omegaQpass,wQpass,wQppass
      common /derivsdchidwQpass/ omegampass,omegavpass,
     .       omegakpass,omegaQpass,wQpass,wQppass
      integer iwmodepass
      common /derivsdchidwQpassint/ iwmodepass	

      double precision hz,h

      h=1.0d0

      call HzzQ(x,hz,h,omegampass,omegavpass,omegakpass
     .        ,omegaQpass,wQpass,wQppass,iwmodepass)

      hz=hz/100.0d0

      if (iwmodepass.eq.1) then
         dydx(1)=1.0d0/hz**3 * omegaQpass * 
     .     (1.0+x)**(3.0*(1.0+wQpass)) * 3.0 *
     .     log(1.0 + x)
      elseif (iwmodepass.eq.2) then
         dydx(1)=1.0d0/hz**3 * omegaQpass * 
     .	   (1.0+x)**(3.0*(1.0+wQpass-wQppass))*
     .	   exp(3.0*wQppass*x) * 3.0 *
     .	   log(1.0 + x)
      elseif (iwmodepass.eq.3) then
         dydx(1)=1.0d0/hz**3 * omegaQpass * 
     .	   (1.0+x)**(3.0*(1.0+wQpass+wQppass))*
     .	   exp(-3.0*wQppass*x/(1.0+x)) * 3.0 *
     .	   log(1.0 + x)
      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dchidwQpcompute(z,dchidwQp,
     .     h,omegam,omegav,omegak,
     .              omegaQ,wQ,wQp,iwmode)

      implicit none
      double precision z,dchidwQp
      double precision h,omegam,omegav,omegak
      double precision omegaQ,wQ,wQp
      integer iwmode
      double precision h00,c00,chm00

      double precision omegampass,omegavpass,
     .       omegakpass,omegaQpass,wQpass,wQppass
      common /derivsdchidwQppass/ omegampass,omegavpass,
     .       omegakpass,omegaQpass,wQpass,wQppass
      integer iwmodepass
      common /derivsdchidwQppassint/ iwmodepass

      integer neq,ind,ier
      parameter (neq=1)
      double precision y(neq),c(24),work(neq,9)
      double precision tol
      external derivsdchidwQp
      double precision xstart,xend

      omegampass=omegam
      omegakpass=omegak
      omegavpass=omegav
      omegaQpass=omegaQ
      wQpass=wQ
      wQppass=wQp
      iwmodepass=iwmode

      h00=100.0d0*h
      c00=2.99792458d5

      chm00=c00/h00

      tol=1.0d-12

      xstart=0.0d0
      xend=z
      ier=0
      ind=1
      y(1)=0.0
      dchidwQp=0.0
      call dverk(neq,derivsdchidwQp,
     .     xstart,y,xend,tol,ind,c
     .      ,neq,work,ier)
      if (ind.lt.0.or.ier.gt.0) then
         write(*,*) 'dverk error, ind, ier=',ind,ier
      end if

      dchidwQp=y(1)*(-chm00)/2.0

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine derivsdchidwQp(n,x,y,dydx)
      implicit none
      integer n
      double precision x,y(n),dydx(n)

      double precision omegampass,omegavpass,
     .       omegakpass,omegaQpass,wQpass,wQppass
      common /derivsdchidwQppass/ omegampass,omegavpass,
     .       omegakpass,omegaQpass,wQpass,wQppass
      integer iwmodepass
      common /derivsdchidwQppassint/ iwmodepass	

      double precision hz,h

      h=1.0d0

      call HzzQ(x,hz,h,omegampass,omegavpass,omegakpass
     .        ,omegaQpass,wQpass,wQppass,iwmodepass)

      hz=hz/100.0d0

      if (iwmodepass.eq.1) then
         dydx(1)=0.0d0
      elseif (iwmodepass.eq.2) then
         dydx(1)=1.0d0/hz**3 * omegaQpass * 
     .	   (1.0+x)**(3.0*(1.0+wQpass-wQppass))*
     .	   exp(3.0*wQppass*x) * ( -3.0 *
     .	   log(1.0 + x) + 3.0*x)
      elseif (iwmodepass.eq.3) then
         dydx(1)=1.0d0/hz**3 * omegaQpass * 
     .	   (1.0+x)**(3.0*(1.0+wQpass+wQppass))*
     .	   exp(-3.0*wQppass*x/(1.0+x)) * (3.0 *
     .	   log(1.0 + x) - 3.0*x/(1.0 + x))
      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dchidomegaQcompute(z,dchidomegaQ,
     .     h,omegam,omegav,omegak,
     .              omegaQ,wQ,wQp,iwmode)

      implicit none
      double precision z,dchidomegaQ
      double precision h,omegam,omegav,omegak
      double precision omegaQ,wQ,wQp
      integer iwmode
      double precision h00,c00,chm00

      double precision omegampass,omegavpass,
     .       omegakpass,omegaQpass,wQpass,wQppass
      common /derivsdchidomegaQpass/ omegampass,omegavpass,
     .       omegakpass,omegaQpass,wQpass,wQppass
      integer iwmodepass
      common /derivsdchidomegaQpassint/ iwmodepass

      integer neq,ind,ier
      parameter (neq=1)
      double precision y(neq),c(24),work(neq,9)
      double precision tol
      external derivsdchidomegaQ
      double precision xstart,xend

      omegampass=omegam
      omegakpass=omegak
      omegavpass=omegav
      omegaQpass=omegaQ
      wQpass=wQ
      wQppass=wQp
      iwmodepass=iwmode

      h00=100.0d0*h
      c00=2.99792458d5

      chm00=c00/h00

      tol=1.0d-12

      xstart=0.0d0
      xend=z
      ier=0
      ind=1
      y(1)=0.0
      dchidomegaQ=0.0
      call dverk(neq,derivsdchidomegaQ,
     .     xstart,y,xend,tol,ind,c
     .      ,neq,work,ier)
      if (ind.lt.0.or.ier.gt.0) then
         write(*,*) 'dverk error, ind, ier=',ind,ier
      end if

      dchidomegaQ=y(1)*(-chm00)/2.0

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine derivsdchidomegaQ(n,x,y,dydx)
      implicit none
      integer n
      double precision x,y(n),dydx(n)

      double precision omegampass,omegavpass,
     .       omegakpass,omegaQpass,wQpass,wQppass
      common /derivsdchidomegaQpass/ omegampass,omegavpass,
     .       omegakpass,omegaQpass,wQpass,wQppass
      integer iwmodepass
      common /derivsdchidomegaQpassint/ iwmodepass

      double precision hz,h

      h=1.0d0

      call HzzQ(x,hz,h,omegampass,omegavpass,omegakpass
     .        ,omegaQpass,wQpass,wQppass,iwmodepass)

      hz=hz/100.0d0

      if (iwmodepass.eq.1) then
         dydx(1)=1.0d0/hz**3 * 
     .     ((1.0+x)**(3.0*(1.0+wQpass)) - (1.0+x)**3)
      elseif (iwmodepass.eq.2) then
	 dydx(1)=1.0d0/hz**3 * 
     .     ((1.0+x)**(3.0*(1.0+wQpass-wQppass))*
     .	   exp(3.0*wQppass*x) - (1.0+x)**3)
      elseif (iwmodepass.eq.3) then
         dydx(1)=1.0d0/hz**3 * 
     .     ((1.0+x)**(3.0*(1.0+wQpass+wQppass))*
     .	   exp(-3.0*wQppass*x/(1.0+x)) - (1.0+x)**3)
      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine dHdaQ(z,dhda,h,omegam,omegav,omegak,
     .               omegaQ,wQ,wQp,iwmode)
c gives dH(a)/da
	implicit none
      	double precision z,dhda,h,omegam,omegav,omegak,h00
      	double precision omegaQ,wQ

      	double precision wQp
      	integer iwmode
	double precision factor

	h00=100.0d0*h

	if (abs(omegam+omegak+omegav+omegaQ-1).gt.1.0d-5) 
     .     	then
         	write(*,*) 'Omegas do not sum to one'
         	write(*,*) 'error in dHdaQ'
         	write(*,*) omegam+omegak+omegav+omegaQ
         	stop
      	endif

	if (iwmode.eq.1) then
         	factor=(omegam*(1.0d0+z)**3 + omegak*(1.0d0+z)**2
     .        		+ omegav + omegaQ*(1.0d0+z)**(3.0d0*(1.0d0+wQ))
     .        		)**0.5d0
		dhda= -(1.0d0+z)**2 * h00 /2.0d0/factor*
     .			(3.0d0*omegam*(1.0d0+z)**2 + 
     .			 2.0d0*omegak*(1.0d0+z) +
     .			 3.0d0*(1.0d0+wQ)*omegaQ*(1.0d0+z)**(2.0d0+3.0d0*wQ))
	elseif (iwmode.eq.2) then
		factor=(omegam*(1.0d0+z)**3 + omegak*(1.0d0+z)**2
     .        		+ omegav + omegaQ*( 
     .	      		(1.0d0+z)**(3.0d0*(1.0d0+wQ-wQp))*
     .	      		exp(3.0d0*wQp*z))
     .        		)**0.5d0
		dhda= -(1.0d0+z)**2 * h00 /2.0d0/factor*
     .			(3.0d0*omegam*(1.0d0+z)**2 + 
     .			 2.0d0*omegak*(1.0d0+z) +
     .			 omegaQ*(
     .			 3.0d0*(1.0d0+wQ-wQp)*(1.0d0+z)**(2.0d0+3.0d0*(wQ-wQp))*
     .			 exp(3.0d0*wQp*z) + 
     .			 (1.0d0+z)**(3.0d0*(1.0d0+wQ-wQp))*3.0d0*wQp*exp(3.0d0*wQp*z)))
	elseif (iwmode.eq.3) then
		factor=(omegam*(1.0d0+z)**3 + omegak*(1.0d0+z)**2
     .        		+ omegav + omegaQ*( 
     .	      		(1.0d0+z)**(3.0d0*(1.0d0+wQ+wQp))*
     .	      		exp(-3.0d0*wQp*z/(1.0d0+z)))
     .        		)**0.5d0
		dhda= -(1.0d0+z)**2 * h00 /2.0d0/factor*
     .			(3.0d0*omegam*(1.0d0+z)**2 + 
     .			 2.0d0*omegak*(1.0d0+z) +
     .			 omegaQ*(
     .			 3.0d0*(1.0d0+wQ+wQp)*(1.0d0+z)**(2.0d0+3.0d0*(wQ+wQp))*
     .			 exp(-3.0d0*wQp*z/(1.0d0+z)) +
     .			 (1.0d0+z)**(3.0d0*(1.0d0+wQ+wQp))*
     .			 (-3.0d0*wQp/(1.0d0+z) + 3.0d0*wQp*z/(1.0d0+z)**2)*
     .			 exp(-3.0d0*wQp*z/(1.0d0+z))))
	endif

	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
