
#ifdef unix
#include <sys/types.h>
#include <unistd.h>
#endif

#include "firsttoinc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "util.h"
#include "expt.h"
#include "volspline.h"

#include "volflat.h"

//#include "swap.h"
#include "gpbase/gplinalgconvert.h"

using ARM::To_ARM_Vector;
using ARM::To_ARM_GP_Vector;
using ARM::To_ARM_GP_Vector;
using ARM::To_pARM_Vector;
using ARM::To_pARM_GP_Vector;
using ARM::ARM_GP_Vector;

#include "gpclosedforms\nonparametric_spline.h"

using ARM::cubicspline_precompute;
using ARM::cubicspline_inter;

void ARM_VolSplineInterpol::Init(void)
{
    SetName(ARM_VOL_SPLINE_INTERPOL);

    SetStrikes(NULL);

	itsDer2			= NULL;
	itsVolsbyRaw	= NULL;

	SetInterpType(K_SPLINE);
}


ARM_VolSplineInterpol::ARM_VolSplineInterpol(void)
{
	Init();
}

ARM_VolSplineInterpol::ARM_VolSplineInterpol(const ARM_Date& asOf, ARM_Vector* yearTerms, 
											 ARM_Vector* strikes, ARM_Matrix* volatilities,
											 int strikeType,
											 int volType,
											 ARM_Currency* ccy) :
ARM_VolLInterpol(asOf,yearTerms,strikes,volatilities,strikeType,volType,ccy), itsDer2(NULL), itsVolsbyRaw(NULL)
{
	SetInterpType(K_SPLINE);

	buildSpline();
}

ARM_VolSplineInterpol::ARM_VolSplineInterpol(const ARM_Date& asOf, ARM_Vector* yearTerms, 
											 ARM_Vector* volatilities, 
											 int strikeType,
											 int volType,
											 ARM_Currency* ccy) : 
ARM_VolLInterpol(asOf,yearTerms,volatilities,strikeType,volType,ccy)
{
	Init();
}

ARM_VolSplineInterpol::ARM_VolSplineInterpol(ARM_VolFlat* volFlat, 
											 ARM_Vector* expiries,
											 ARM_Vector* undTenors) :
ARM_VolLInterpol(volFlat,expiries,undTenors)
{
	Init();
}

ARM_VolSplineInterpol::ARM_VolSplineInterpol(const ARM_VolSplineInterpol& rhs) : 
ARM_VolLInterpol(rhs)
{
	Init();

	BitwiseCopy(&rhs);
}

ARM_VolSplineInterpol::~ARM_VolSplineInterpol()
{
	if(itsDer2)
		delete [] itsDer2;

	if(itsVolsbyRaw)
		delete [] itsVolsbyRaw;
}

ARM_VolSplineInterpol& ARM_VolSplineInterpol::operator = (const ARM_VolSplineInterpol& rhs)
{
    (*this).ARM_VolCurve::operator = (rhs);

    BitwiseCopy(&rhs);

    return(*this);
}

void ARM_VolSplineInterpol::BitwiseCopy(const ARM_Object * srcVolspline)
{
	ARM_VolLInterpol::BitwiseCopy(srcVolspline);

	ARM_VolSplineInterpol * Volspline = (ARM_VolSplineInterpol*)(srcVolspline);

	if(itsDer2) 
		delete [] itsDer2;

	if(itsVolsbyRaw)
		delete [] itsVolsbyRaw;

	if(Volspline->itsDer2)
	{
		int i, size = (int)GetExpiryTerms()->size();
		itsDer2 = new ARM_Vector [size];
		for(i = 0; i < size; i++)
			itsDer2[i] = (*(ARM_Vector*)(Volspline->itsDer2[i].Clone()));
	}

	if(Volspline->itsVolsbyRaw)
	{
		int i, size = (int)GetExpiryTerms()->size();
		itsVolsbyRaw = new ARM_Vector [size];
		for(i = 0; i < size; i++)
			itsVolsbyRaw[i] = (*(ARM_Vector*)(Volspline->itsVolsbyRaw[i].Clone()));
	}
}

void ARM_VolSplineInterpol::Copy(const ARM_Object* vollint)
{
    ARM_VolLInterpol::Copy(vollint);

    BitwiseCopy(vollint);
}

ARM_Object* ARM_VolSplineInterpol::Clone(void)
{
    ARM_VolSplineInterpol* theClone = new ARM_VolSplineInterpol();

    theClone->Copy(this);

    return(theClone);
}

void ARM_VolSplineInterpol::buildSpline()
{
	if(itsDer2)
	{
		delete [] itsDer2;
		itsDer2 = NULL;

		delete [] itsVolsbyRaw;
		itsVolsbyRaw = NULL;
	}

	if(GetStrikes())
	{
		ARM_Vector * strikes = GetStrikes();
		ARM_Matrix * vols = GetVolatilities();

		int i, size = (int)GetExpiryTerms()->size();
		int k, ksize = strikes->size();

		itsDer2 = new ARM_Vector [size];
		itsVolsbyRaw = new ARM_Vector [size];

		for(i = 0; i < size; i++)
		{
			itsVolsbyRaw[i].Resize(ksize);
			itsDer2[i].Resize(ksize);

			for(k = 0; k < ksize; k++) itsVolsbyRaw[i][k] = (*vols).Elt(i,k);

			cubicspline_precompute(strikes->begin(), itsVolsbyRaw[i].begin(), ksize, 1e30,1e30, itsDer2[i].begin());
		}
	}
}

double ARM_VolSplineInterpol::VolatilityFunction(double m1, double m2)
{
    return(VolatilityFunction(m1, m2, 0.0));
}

double ARM_VolSplineInterpol::VolatilityFunction(double m1, double K, double m2)
{
	if(itsDer2 == NULL) return ARM_VolLInterpol::VolatilityFunction(m1,K,m2);

	ARM_Vector* strikes = GetStrikes();

    ARM_Matrix* vols = GetVolatilities();

    ARM_Vector* terms = GetExpiryTerms();

    if (( K == 0.0 ) && (!( m2 == 0.0 ))) // Here we are AT THE MONEY
    {
       K = m2;
    }
    else
    {
       m2 = 0.0; // m2 : is not relevant
    }

	int i, tsize = terms->size();

	if(m1 <= (*terms)[0])
		return cubicspline_inter(strikes->begin(), itsVolsbyRaw[0].begin(), itsDer2[0].begin(), itsDer2[0].size(), K);

	if(m1 >= (*terms)[tsize-1])
		return cubicspline_inter(strikes->begin(), itsVolsbyRaw[tsize-1].begin(), itsDer2[tsize-1].begin(), itsDer2[tsize-1].size(), K);

	i = 0;
	while((*terms)[i] < m1) i++;

	double vol1 = cubicspline_inter(strikes->begin(), itsVolsbyRaw[i-1].begin(), itsDer2[i-1].begin(), itsDer2[i-1].size(), K);
	double vol2 = cubicspline_inter(strikes->begin(), itsVolsbyRaw[i].begin(), itsDer2[i].begin(), itsDer2[i].size(), K);

	return vol1 + (m1 - (*terms)[i-1])*(vol2-vol1)/((*terms)[i]-(*terms)[i-1]);
}


