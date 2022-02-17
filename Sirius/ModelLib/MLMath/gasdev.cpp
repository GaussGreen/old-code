#include "stdafx.h"
#include <math.h>


#include "ran2.h"

double gasdev(long *idum)
{
	double ran2(long *idum);
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if (*idum < 0) iset=0;
	if  (iset == 0) {
		do {
			v1=2.0*ran2(idum)-1.0;
			v2=2.0*ran2(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}	

	

double gasdev1(long *idum,double lo,double lup)
{	
	double ran2(long *idum);
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;
	double rand1,rand2;
	
	if (*idum < 0) iset=0;
	if  (iset == 0) {
		do {

			rand1 = (lup-lo)*ran2(idum)+lo;
//			rand2 = (lup-lo)*(1.0-ran2(&idum)+0.5)+lo;
			rand2 = ran2(idum);

			v1=2.0*rand1-1.0;
			v2=2.0*rand2-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);


		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}	





