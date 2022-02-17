/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2005
 *
 *  solve the tridiagonal matrix problem
 *
 *	\file tridiagonalsolve.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date dec 2005
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include <cmath>

#include <vector>
#include "gpbase/numericconstant.h"
#include "gpbase/gplinalgconvert.h"
#include "gpbase/gpvector.h"

#include "expt.h"   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)
/*
a :diagonale inferieure
b: diagonale
c: diagonale superieure
r: vecteur constant
u: resultat
gam: vecteur tampon
  */

void tridiagonalsolve(const vector<double>& a, const vector<double>& b, const vector<double>& c,const vector<double>& r, vector<double>& u,vector<double>& gam)
{
	int n=b.size();
	int i;
	if (b[0] == 0.0) throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"b[1] == 0 in tridiagonalsolve" );;
	gam[0]=b[0];
	u[0]=r[0]/gam[0];
	for (i=1;i<n;i++) {
		gam[i]=b[i]-a[i-1]*c[i-1]/gam[i-1];
		if (gam[i] == 0.0) 
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"gam[i] == 0 in tridiagonalsolve" );
		}
		u[i]=(r[i]-a[i-1]*u[i-1])/gam[i];
	}
	for (i=n-2;i>=0;i--)
		u[i] -= c[i]*u[i+1]/gam[i];

}




CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/