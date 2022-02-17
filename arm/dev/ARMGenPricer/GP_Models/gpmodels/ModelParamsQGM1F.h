/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ModelParamsQGM1F.h,v $
 * Revision 1.1  2004/07/30 14:45:48  jmprie
 * Initial revision
 *
 *
 */



/*! \file ModelParamsQGM1F.h
 *
 *  \brief 
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date July 2004
 */


#ifndef _INGPMODELS_MODELPARAMSQGM1F_H
#define _INGPMODELS_MODELPARAMSQGM1F_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpinfra/modelparams.h"

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_ModelParamsQGM1F
// \brief Class for model parameters of 1F Quadratic Gaussian Model
//-----------------------------------------------------------------------------
class ARM_ModelParamsQGM1F : public ARM_ModelParams 
{
public:
	ARM_ModelParamsQGM1F( const ARM_ModelParamsQGM1F& rhs );
	ARM_ModelParamsQGM1F( const ARM_ModelParamVector& params=ARM_ModelParamVector() );
	virtual ~ARM_ModelParamsQGM1F();
    ARM_ModelParamsQGM1F& operator = (const ARM_ModelParamsQGM1F& rhs);

	/// How many factors?
    virtual size_t FactorCount() const { return 1; };

    /// Drift from a to b of the state variable
    double StateLocalDrift(double a,double b) const;

    /// Variance in [a,b] of the state variable
    double StateLocalVariance(double a,double b) const;

    /// Variance in [a,b] of Zc(.,T1)/Zc(.,T2)
    double FwdZcLocalVariance(double a,double b,double T1,double T2) const;

    /// Covariance in [a,b] of Zc(.,T1)/Zc(.,U1) and Zc(.,T2)/Zc(.,U2)
    double FwdZcLocalCovariance(double a,double b,double T1,double U1,double T2,double U2,ARM_GP_Vector& vars) const;

    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0);
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model ,int factorNb=0) {};

	// Standard ARM object support
	virtual ARM_Object* Clone() const;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

