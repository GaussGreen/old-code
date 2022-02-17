/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file sabr_calibration.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include <glob/firsttoinc.h>
#include <glob/expt.h>   // for the exceptions

#include <cmath>
#include <complex>


#include "gpclosedforms/qgm_skewcalibration.h"
#include "gpclosedforms/optimization1.h"

#include "nag.h"
#include "nage04.h"

#include "gpbase/port.h"
#include "gpbase/removenagwarning.h"

using namespace std;

CC_BEGIN_NAMESPACE(ARM)


QGM_ParameterSet::QGM_ParameterSet(Optimization_Result_Set* r)
{
	_X = (ARM_GP_Vector*) r->OptimalParamSet->Clone();
	_objective=r->OptimalObjective;
}

QGM_ParameterSet* QGM_CalibrateSkew(ARM_GP_Vector* tagetVect,
	  ARM_GP_Vector* weights,
      ARM_GP_Vector* precisions,
	  ARM_GP_Vector* InitVector,
	  ARM_GP_Vector* LBoundVector,
	  ARM_GP_Vector* UBoundVector,
	  int algorithm)
{
	int sizeVector = tagetVect->size();

	class objectiveFuntion : public Optimization_ObjectiveFuntion
	{
    private:
        ARM_GP_Vector itsPresicions;
        ARM_GP_Vector itsWieghts;
        ARM_GP_Vector itsTagetVect;
	  
	public:
		void NAG_CALL operator() (Integer m, Integer n, 
			double x[], /// input
			double f[],	/// output (f(x))
			double fjac[],  /// output  (Df(x,i))
			Integer tdfjac, Nag_Comm *comm)
			
		{
			int i;
            double sumWeight = itsWieghts.sum();
			for (int l =0; l<m*n; ++l)
				fjac[l] = 0;
			for(i=0;i<m;i++)
			{
				f[i] = 0;
                double precision = itsPresicions[i];
                double sqrtWeight = sqrt(itsWieghts[i]/sumWeight);
                double a = sqrtWeight/precision;
				for(int j=i; j<m; ++j)
				{               
					f[i]+=x[j];
					fjac[m*i+j] = 1.0/(m-i)*a;
				}
				f[i]=f[i]/double((m-i))*a;
				double target = itsTagetVect[i];
                f[i]-= target*a;
			}
		}
		objectiveFuntion(ARM_GP_Vector& precisions, ARM_GP_Vector& weights,ARM_GP_Vector& tagetVect)
        :itsPresicions(precisions),
        itsWieghts(weights),
        itsTagetVect(tagetVect)
        {}	
	};
	objectiveFuntion func(*precisions,*weights,*tagetVect);

	Optimization_Result_Set* result= OptimizeWithDerivatives(&ARM_GP_Vector(sizeVector,0.0),
		&ARM_GP_Vector(sizeVector,0.0),
		&ARM_GP_Vector(sizeVector,1.0),
		&func,
		InitVector,
		LBoundVector,
		UBoundVector,
		algorithm,
		false,
		//"C:\\temp\\QGMnag.txt",
		"",
		0.0001,  // Tolerance
		50);	// max_iter

	QGM_ParameterSet* setptr= new QGM_ParameterSet(result);
	delete result;
	return setptr;


}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
