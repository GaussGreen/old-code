/***********************************************************
*
*	kutility.h
*
************************************************************/


#ifndef _BAS_H_
#define	_BAS_H_
#include "kplatform.h"
#include <iostream>
#include <valarray>
#include <vector>  

//using namespace std;


const	int		K_NO_INTERP =	0;
const	int		K_LINEAR_INTERP = 1;
const	int		K_STAIR_INTERP = 2;
const	int		K_FLAT_FWD_INTERP = 3;
const	int		K_COUPON_B30	= 4;
const	int		K_COUPON_ACT	= 5;

const	int		K_BOND_STUB	= 0;
const	int		K_SIMPLE_STUB	= 1;
const	int		K_NO_STUB		= 2;

const	int		K_SUCCESS = 1;   
const	int		K_FAILURE	= 0;
//tree trimming parameter
const	int		K_NUM_SIGMA = 5;
const	double	K_PI = 3.141592653;
const	double	K_JUMPCOEFF = 1.5;


template	<class T>
std::ostream & operator<<(std::ostream &out, const std::vector<T> &m);

template <class T>
T	Max(const T&a, const T&b){return (a > b ? a:b );}
template <class T>
T	Min(const T&a, const T&b){return (a < b ? a:b );}

template	<class T>
std::ostream & operator<<(std::ostream &out, const std::vector<T> &m)
{
	out<<"[ ";
	for(int i=0; i<m.size(); i++)
		out<<m[i]<<" ";
	out<<"]";
	return(out);
}
//return C(k, N)p^k(1-p)^(N-k)
double	BiNomial(int k, int N, double p);
double	Normal(double);
double	NormalInv(double);
double	BiNormal(double a, double b, double rho);
double	BlackOption(double	fwd, 
					double	strike, 
					double	t, 
					double	vol, 
					char	coP);
#endif



