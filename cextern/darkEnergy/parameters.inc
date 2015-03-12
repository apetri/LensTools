c  This is slim version of snlens.inc, including only fundamental parameters that never change:
c  i.e. number of interpolation points (na), speed of light (cspeed), and max and min z over which range is calculated (zmaxact, zminact). 
c  The only parameter that really needs to be declared like this is "na", because it gives the size of alloctated arrays.

c	double precision zmaxact
c	parameter (zmaxact=65.0d0)
c	double precision zminact
c	parameter (zminact=0.0d0)
c Lam:	parameter (zminact=-0.99d0)


c   integer iwmodeI
c   parameter (iwmodeI=3)

        integer na
        parameter (na=5000)	

	double precision cspeed
	parameter (cspeed=2.99792458d5)
