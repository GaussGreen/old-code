/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file smoother.cpp
 *
 *  \brief
 *
 *	\author  J-M Prié
 *	\version 1.0
 *	\date february 2005
 */

/// to remove identified warning
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/env.h"
#include "gpbase/gpvector.h"

#include "gpnummethods/smoother.h"

/// gpnumlib
#include "gpnumlib/solver.h"
#include "gpnumlib/newtonraphson.h"


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_SmootherLinear
///	Routine: Compute
///	Returns: void
///	Action : Compute a correction assuming a linear exercise
///          function around the boundary located between
///          exerFct[exerIdx] and exerFct[exerIdx+1]
////////////////////////////////////////////////////
double ARM_SmootherLinear::Compute(const ARM_GP_Vector& exerFct, int exerIdx, double coef, ARM_GP_Vector& smooth) const
{
    double f0=exerFct[exerIdx];
    double f1=exerFct[exerIdx+1];

    /// Space step is cst on a slice then we can say x0=0, x1=1 and exerFct(z)=0
    double z=-f0/(f1-f0);
    double fmid = 0.5*(f0+f1);

    /// Compute the average value of the positive exercise fct between
    /// middle of the interval and the boundary
    double exerFctAvge = 0.5*fmid*(0.5-z);

    if(z < 0.5)
        smooth[exerIdx] += coef*exerFctAvge;
    else
        smooth[exerIdx+1] -= coef*exerFctAvge;

    return z;
}


////////////////////////////////////////////////////
///	Class  : ARM_SmootherQuadratic
///	Routine: Compute
///	Returns: void
///	Action : Compute a correction assuming a quadratic exercise
///          function around the unique boundary located between
///          exerFct[exerIdx] and exerFct[exerIdx+1]
////////////////////////////////////////////////////
double ARM_SmootherQuadratic::Compute(const ARM_GP_Vector& exerFct, int exerIdx, double coef, ARM_GP_Vector& smooth) const
{
    if(exerIdx > 0 && exerIdx < exerFct.size()-1)
    {
        /// Find a, b & c such that exerFct(x)=a.x.x + 2.bp.x + c
        /// Space step is cst on a slice then we can say x0=0, x1=1, x2=2 and exerFct(z)=0
        double f0=exerFct[exerIdx];
        double f1=exerFct[exerIdx+1];
        double f2=exerFct[exerIdx+2];
        double df20=0.5*(f2-f0);
        double df10=f1-f0;
        double a=df20-df10;
        double bp=df10-0.5*df20;
        double c=f0;
        double dp=sqrt(bp*bp-a*c);
        double z0=(-bp+dp)/a;
        double z1=(-bp-dp)/a;
        double z;
        if(0.0<=z0 && z0<=1.0)
            z=z0;
        else if(0.0<=z1 && z1<=1.0)
            z=z1;
        else
            /// No way then back to linear smoothing
            return ARM_SmootherLinear::Compute(exerFct,exerIdx,coef,smooth);

        /// Compute the average value of the positive exercise fct between
        /// middle of the interval and the boundary
        double z2=z*z;
        double z3=z2*z;

        double exerFctAvge=a/3*(0.125-z3) + bp*(0.25-z2) + c*(0.5-z);
      
        if(z < 0.5)
            smooth[exerIdx] += coef*exerFctAvge;
        else
            smooth[exerIdx+1] -= coef*exerFctAvge;

        return z;
    }
    else
        return ARM_SmootherLinear::Compute(exerFct,exerIdx,coef,smooth);
}


////////////////////////////////////////////////////
///	Class  : ARM_SmootherCubic
///	Routine: Compute
///	Returns: void
///	Action : Compute a correction assuming a cubic exercise
///          function around the boundary located between
///          exerFct[exerIdx] and exerFct[exerIdx+1]
////////////////////////////////////////////////////
double ARM_SmootherCubic::Compute(const ARM_GP_Vector& exerFct, int exerIdx, double coef, ARM_GP_Vector& smooth) const
{
    size_t nbStates=exerFct.size();
    if(exerIdx > 1 && exerIdx < nbStates-1)
    {
        double f_1=exerFct[exerIdx-1];
        double f0=exerFct[exerIdx];
        double f1=exerFct[exerIdx+1];
        double f2=exerFct[exerIdx+2];
        double df_10=f_1-f0;
        double df10=f1-f0;
        double df20=0.5*(f2-f0);
        double b=0.5*(df_10+df10);
        double a=(df20-df10-b)/3;
        double c=df10-a-b;
        double d=f0;

        /// Solve exerFct(x)=0 with an initial guess = solution of linear approximation
        ARM_SmootherCubic::ExerciseFunction fctToSlove(a,b,c,d);
        T_NewtonRaphsonSolverNoThrow<ARM_SmootherCubic::ExerciseFunction> solver(fctToSlove);

        /// To initialize departure point
		solver.setInitialGuess(-f0/(f1-f0));
		double z = solver.Solve();
        if(z<0.0 || z>1.0)
            /// No way then back to linear smoothing
            return ARM_SmootherLinear::Compute(exerFct,exerIdx,coef,smooth);

        /// Compute the average value of the positive exercise fct between
        /// middle of the interval and the boundary
        double z2=z*z;
        double z3=z2*z;
        double z4=z3*z;

        double exerFctAvge=0.25*a*(0.0625-z4) + b/3*(0.125-z3)+ 0.5*c*(0.25-z2) + d*(0.5-z);
      
        if(z < 0.5)
            smooth[exerIdx] += coef*exerFctAvge;
        else
            smooth[exerIdx+1] -= coef*exerFctAvge;

        return z;
    }
    else if(exerIdx > 0 && exerIdx < nbStates-1)
        return ARM_SmootherQuadratic::Compute(exerFct,exerIdx,coef,smooth);
    else
        return ARM_SmootherLinear::Compute(exerFct,exerIdx,coef,smooth);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
