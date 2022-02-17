/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file hypergeometric.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include <cmath>
#include <complex>

#include "gpclosedforms/gaussian_integrals.h"
#include "gpbase/numericconstant.h"

#include <glob/expt.h>   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_EPS 1.0e-12
#define ARM_CF_MAXIT 1000



complex<double>  Hypergeometric1F1_serie ( complex<double>  a, complex<double>  b, complex<double>  z)
{
	
	const int NbMax = 85;

	complex<double>  an(a);
	complex<double>  bn(b);
	complex<double>	 n(1.,0.);
	complex<double>  sum_totale(1,0.);
	complex<double>  sum_partielle(1,0.);
	complex<double>  un(1,0);

	int ctr = 0;
	do 
	{
		ctr++;
		sum_partielle*=z*an/(bn*n);
		sum_totale+=sum_partielle;
		an+=un;
		bn+=un;
		n+=un;
		if ( ctr> NbMax) break;

	}
	while (abs(sum_partielle)>ARM_CF_EPS);
	return sum_totale;
}
/*
complex<double>  Hypergeometric1F1_ContinuousFraction ( complex<double>  a, complex<double>  b, complex<double>  z)
{

	const int nb(50);
	vector<complex<double> > an(nb);

	for( int i=0 ; i< nb; ++i)	
		an[i] = ( (i+1.0) + a)* z /((i+2.0)*((i+1.0)+ b ));

	complex<double>	tmp = an[nb-1];
	for ( i=1; i<nb; ++i){
		tmp *= -1.0;
		tmp += 1.0 + ( nb-1-i>0?1.0:0.0) * an[nb-1-i];
		tmp /= an[nb-1-i];
	}

	return 1.0+tmp;
}
*/
complex<double>  Hypergeometric1F1 ( complex<double>  a, complex<double>  b, complex<double>  z)		// Kummer Function
{
	return Hypergeometric1F1_serie(a,b,z);
}

complex<double>  Hypergeometric2F1_serie ( complex<double>  a, complex<double>  b, complex<double>  c, complex<double>  z)
{
	complex<double>  an(a);
	complex<double>  bn(b);
	complex<double>  cn(c);
	complex<double>	 n(1.,0.);
	complex<double>  sum_totale(1,0.);
	complex<double>  sum_partielle(1,0.);
	complex<double>  un(1,0);
	do 
	{
		sum_partielle*=z*an*bn/(cn*n);
		sum_totale+=sum_partielle;
		an+=un;
		bn+=un;
		cn+=un;
		n+=un;
	}
	while (abs(sum_partielle)>ARM_CF_EPS);
	return sum_totale;
}

complex<double>  Hypergeometric2F1 ( complex<double>  a, complex<double>  b, complex<double>  c, complex<double>  z)
{
	complex<double>  un(1,0);
	if (abs(z)>=1)
		return  pow(un-z,-a)*Hypergeometric2F1_serie(a,b,c,z/(z-un));
	else
		return Hypergeometric2F1_serie(a,b,c,z);
	
}

complex<double>  Hypergeometric2F0_serie ( complex<double>  a, complex<double>  b, complex<double>  z)
{
	complex<double>  an(a);
	complex<double>  bn(b);
	complex<double>	 n(1.,0.);
	complex<double>  sum_totale(1,0.);
	complex<double>  sum_partielle(1,0.);
	complex<double>  un(1,0);
	do 
	{
		sum_partielle*=z*an*bn/n;
		sum_totale+=sum_partielle;
		an+=un;
		bn+=un;
		n+=un;
	}
	while (abs(sum_partielle)>ARM_CF_EPS);
	return sum_totale;
}

double  Hypergeometric2F1_serie ( double  a, double  b, double  c, double  z)
{
	double  an=a;
	double  bn=b;
	double  cn=c;
	double	n=1.;
	double  sum_totale=1;
	double  sum_partielle=1;
	double  un=1.;
	do 
	{
		sum_partielle*=z*an*bn/(cn*n);
		sum_totale+=sum_partielle;
		an+=un;
		bn+=un;
		cn+=un;
		n+=un;
	}
	while (fabs(sum_partielle)>ARM_CF_EPS);
	return sum_totale;
}

double  Hypergeometric2F1 ( double  a, double  b, double  c, double  z)
{
	if(z<-0.5)
		return pow(1.-z,-a)*Hypergeometric2F1_serie(a,c-b,c,z/(z-1.));
	else
		return Hypergeometric2F1_serie(a,b,c,z);
}

long double   Hypergeometric2F0_serie ( long double   a, long double   b, long double   z)
{
	long double  an(a);
	long double  bn(b);
	long double	 n(1.);
	long double  sum_totale(1.);
	long double  sum_partielle(1.);
	long double  un(1.);
	do 
	{
		sum_partielle*=z*an*bn/n;
		sum_totale+=sum_partielle;
		an+=un;
		bn+=un;
		n+=un;
	}
	while (abs(sum_partielle)>ARM_CF_EPS);
	return sum_totale;
}

complex<double>  Hypergeometric2F0 ( complex<double>  a, complex<double>  b, complex<double>  z)
{
	return Hypergeometric2F0_serie(a,b,z);
}

long double  Hypergeometric2F0 ( long double  a, long double  b, long double  z)
{
	return Hypergeometric2F0_serie(a,b,z);
}

complex<double>  HypergeometricU ( complex<double>  a, complex<double>  b, complex<double>  z)
{
	complex<double>  un(1,0);
	return Hypergeometric2F0(a,un+a-b,-un/z)*pow(z,-a);
}

///  |z|>1 
long double   HypergeometricU ( long double   a, long double   b, long double   z)
{
	long double  un(1.);
	return Hypergeometric2F0(a,un+a-b,-un/z)*pow(z,-a);
}

/***************************AppellF1 hypergeometric*******************/

complex<double>  HypergeometricAppellF1_serie ( complex<double>  a, complex<double>  b1,complex<double>  b2, complex<double>  c, complex<double>  x,complex<double> y, int Nb)
{
	if ((abs(x)>=1)||(abs(y)>=1) )	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"HypergeometricAppellF1_serie: x or y  has an absolute value >=1" );
	}
	int i,j;
	complex<double> sumtotale,currentterm,firstterm,amnfirst,cmnfirst,nx,my,amn,cmn,b1n,b2m;
	complex<double>  un(1,0);
	complex<double>  zero(0,0);
	currentterm=un;
	nx=un;
	my=un;
	b2m=b2;
	firstterm=un;
	b1n=b1;
	sumtotale=zero;
	amnfirst=a;
	cmnfirst=c;
	for(i=0;i<Nb;i++)
	{
		currentterm=firstterm;
		b2m=b2;
		amn=amnfirst;
		cmn=cmnfirst;
		my=un;
		for(j=0;j<Nb;j++)
		{
			sumtotale+=currentterm;
			currentterm*=amn*b2m*y/cmn/my;
			amn+=un;
			b2m+=un;
			cmn+=un;
			my+=un;
		}
		firstterm*=x*amnfirst*b1n/cmnfirst/nx;
		amnfirst+=un;
		cmnfirst+=un;
		b1n+=un;
		nx+=un;
	}
	return sumtotale;

}

complex<double>  HypergeometricAppellF1 ( complex<double>  a, complex<double>  b1,complex<double>  b2, complex<double>  c, complex<double>  z1,complex<double> z2, int Nb)
{
	complex<double> un(1.,0);
	double abs1=abs(z1);
	double abs2=abs(z1);
	if((abs1<1)&&(abs2<1)) return  HypergeometricAppellF1_serie (a,b1,b2,c,z1,z2,Nb);
	if((abs1<1)&&(abs2>1)) return  pow(un-z1,c-a-b1)*pow(un-z2,-b2)*HypergeometricAppellF1_serie (c-a,c-b1-b2,b2,c,z1,(z2-z1)/(z2-un),Nb);
	if((abs1>1)&&(abs2<1)) return  pow(un-z1,-b1)*pow(un-z2,c-a-b2)*HypergeometricAppellF1_serie (c-a,c-b1-b2,b2,c,(z1-z2)/(z1-un),z2,Nb);
	else
	{
		if((abs((z2-z1)/(z2-un))<1)&&(abs(z2/(z2-un))<1)) return pow(un-z2,-a)*HypergeometricAppellF1_serie (a,b1,c-b1-b2,c,(z2-z1)/(z2-un),z2/(z2-un),Nb);
		if((abs((z1-z2)/(z1-un))<1)&&(abs(z1/(z1-un))<1))return pow(un-z1,-a)*HypergeometricAppellF1_serie (a,b1,c-b1-b2,c,z1/(z1-un),(z1-z2)/(z1-un),Nb);
		
		if((abs(z2/(z2-un))<1)&&(abs(z1/(z1-un))<1)) return pow(un-z1,-b1)*pow(un-z2,-b2)*HypergeometricAppellF1_serie (c-a,b1,b2,c,z1/(z1-un),z2/(z2-un),Nb);
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"HypergeometricAppellF1: x or y  has an absolute value >=1 and I do not know how to do" );
	}
}





CC_END_NAMESPACE()


#undef ARM_CF_EPS
#undef ARM_CF_MAXIT

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/