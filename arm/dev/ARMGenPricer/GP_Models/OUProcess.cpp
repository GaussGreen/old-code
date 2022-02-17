/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file OUProcess.cpp
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date July 2004
 */


#include "gpmodels/ouprocess.h"
#include "gpbase/typedef.h"
#include "gpbase/curve.h"

#include <math.h>

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Struct : ARM_OUProcess
///	Routine: Drift
///	Returns: value of the drift
///	Action : Relative drift of the process from a to
///          b>=a : X(b) = Drift(a,b)*X(a) + Variance(a,b)
////////////////////////////////////////////////////
double ARM_OUProcess::Drift(double a,double b,double mrs)
{
    if(fabs(mrs)>K_NEW_DOUBLE_TOL)
        return exp(-mrs*(b-a)/K_YEAR_LEN);
    else
        return 1.0;
}


////////////////////////////////////////////////////
///	Class  : ARM_OUProcess
///	Routine: Variance
///	Returns: value of the variance
///	Action : Compute the variance of the process
//           from a to b>=a : X(b) = Drift(a,b)*X(a) + Variance(a,b)
////////////////////////////////////////////////////
double ARM_OUProcess::Variance(double a,double b,const ARM_Curve& sigma,double mrs)
{
    if(fabs(mrs)<=K_NEW_DOUBLE_TOL)
        mrs=(mrs>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

    double scale = 2*mrs;

    return IntegrateScaledSigma(a,b,sigma,scale) * exp(-scale*b/K_YEAR_LEN);
}


////////////////////////////////////////////////////
///	Class  : ARM_OUProcess
///	Routine: BasicVariance
///	Returns: value of the integral
///	Action : Integrate sigma(t)^2*exp(scale*t)
///          between [a,b]
///          Sigma curve is assumed to be stepwise
///          right constant
////////////////////////////////////////////////////
double ARM_OUProcess::IntegrateScaledSigma(double a, double b, const ARM_Curve& sigma, double scale)
{

    const ARM_GP_Vector& sigmaValues = sigma.GetOrdinates(); // Ui
    const ARM_GP_Vector& sigmaTimes = sigma.GetAbscisses(); // Sigma(Ui)

    /// Find Un(a)-1 < a <= Un(a) and Un(b)-1< b <= Un(b)
    int i,nbU = sigmaTimes.size();
    for(i=0;i<nbU;++i)
        if(a <= sigmaTimes[i] + K_NEW_DOUBLE_TOL) break;
    int na=i;
    for(;i<nbU;++i)
        if(b <= sigmaTimes[i] + K_NEW_DOUBLE_TOL) break;
    int nb=i;

    double U,lastU,lastExpScale,newExpScale,sig;
    double timeScale=scale/K_YEAR_LEN;
    bool notTiny=(fabs(scale)>0.1*K_NEW_DOUBLE_TOL);

    double value=0.0;

    lastExpScale=exp(timeScale*a);
    lastU=a;

    /// Between a and Un(b)-1 then b
    for(i=na;i<=nb;++i)
    {
        U=(i<nb ? sigmaTimes[i] : b);
        sig=sigmaValues[i<nbU ? i : nbU-1];
        newExpScale=exp(timeScale*U);
        value += sig*sig*(notTiny ? newExpScale-lastExpScale : U-lastU);
        lastExpScale=newExpScale;
        lastU=U;
    }

    return value/(notTiny ? scale : K_YEAR_LEN);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

