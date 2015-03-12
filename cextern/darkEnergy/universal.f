      function phiPD(x,aneff,omegam,
     .             omegav,omegak,omegaQ)
      implicit none
      double precision x,aneff,omegam,omegav,omegak,omegaQ
      double precision APD,BPD,alphaPD,betaPD,VPD
      double precision g,g3,fac,phiPD,growthfac

      APD = 0.482d0*(1.0d0+aneff/3.0d0)**(-0.947)
      BPD = 0.226d0*(1.0d0+aneff/3.0d0)**(-1.778)
      alphaPD = 3.310d0*(1.0d0+aneff/3.0d0)**(-0.244)
      betaPD = 0.862d0*(1.0d0+aneff/3.0d0)**(-0.287)
      VPD = 11.55d0*(1.0d0+aneff/3.0d0)**(-0.423)
  
      g=growthfac(0.0d0,omegam,omegav,omegak,omegaQ)
      g3=g**3

      fac=1.0d0 + BPD*betaPD*x + (APD*x)**(alphaPD*betaPD)
      fac=fac/(1.0d0 + 
     .      ( (APD*x)**alphaPD*g3/(VPD*x**0.5) )**betaPD)
      phiPD = x*(fac**(1.0d0/betaPD))

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function growthfac(z,omegam,omegav,omegak,omegaQ)
c this gives an approximate growth factor (from Carroll et al. 
c ann. rev. of A. and A. 1992)
c put z to 0 to get what people call g (omega)
c
c note that this growthfactor is not normalized to 1 at z=0 !
      implicit none
      double precision z,omegam,omegav,omegak,omegaQ
      double precision growthfac
      double precision omegamAA,omegavAA,AA,AA3,fgQQ

c      if (omegaQ.ne.0.0) then
c         write(*,*) 'can only allow lambda not Q'
c         stop
c      endif
      if (abs(omegam+omegav+omegak+omegaQ-1.0d0).gt.1.0d-5) then
c      if ((omegam+omegav+omegak+omegaQ).ne.1.0) then
         write(*,*) 'omegam omegav omegak do not sum to 1'
         stop
      endif

      if (omegaQ.eq.0.0) then

      AA=1.0d0/(1.0d0+z)
      AA3=AA**3

      omegamAA=omegam/(AA+omegam*(1.0d0-AA)+omegav*(AA3-AA))
      omegavAA=AA3*omegav/(AA+omegam*(1.0d0-AA)+omegav*(AA3-AA))
      
      growthfac=2.5d0/(1.0d0+z)*omegamAA/
     .    (omegamAA**(4.0d0/7.0d0)-omegavAA+
     .    (1.0d0+omegamAA/2.0d0)*(1.0d0+omegavAA/70.0d0))

      else

         growthfac=fgQQ(z)/(1.0d0+z)

      endif
      

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function Tgrowthfac(z,omegam,omegav,omegak,omegaQ)
      implicit none
      double precision z,omegam,omegav,omegaQ,Tgrowthfac
      double precision growthfac,omegak

      Tgrowthfac=growthfac(z,omegam,omegav,omegak,omegaQ)
     .          /growthfac(0.0d0,omegam,omegav,omegak,omegaQ)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function phiJMW(x)
      implicit none
      double precision phiJMW,x
      double precision top,bottom

      top=1.0+0.6*x+x*x-0.2*x**3-1.5*(x**3.5)+x**4
      bottom=1.0+0.0037*x**3
      
      phiJMW=x*(top/bottom)**0.5
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
