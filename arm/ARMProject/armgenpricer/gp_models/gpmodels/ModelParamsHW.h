/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ModelParamsHW.h,v $
 * Revision 1.1  2003/10/13 07:51:48  jmprie
 * Initial revision
 *
 *
 */



/*! \file ModelParamsHW.h
 *
 *  \brief 
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPMODELS_MODELPARAMSHW_H
#define _INGPMODELS_MODELPARAMSHW_H

#include "gpbase/port.h"
#include "gpinfra/ModelParams.h"

CC_BEGIN_NAMESPACE( ARM )


//-----------------------------------------------------------------------------
// \class ARM_ModelParamsHW
// \brief Interface class for model parameters of Hull & White models
//-----------------------------------------------------------------------------
class ARM_ModelParamsHW : public ARM_ModelParams 
{
public:

    static const double MrsMinValue;
    static const double VOL_LIMIT;


	ARM_ModelParamsHW( const ARM_ModelParamsHW& rhs );
	ARM_ModelParamsHW( const ARM_ModelParamVector& params=ARM_ModelParamVector() );
	virtual ~ARM_ModelParamsHW();
    ARM_ModelParamsHW& operator = (const ARM_ModelParamsHW& rhs);

	/// Discretization Dates of the model params
	virtual ARM_GP_VectorPtr ModelParamsTimeSteps() const = 0;

	/// How many factors?
	virtual size_t FactorCount() const = 0;

    /// Variance in [a,b] of Zc(.,T1)/Zc(.,T2)
    virtual double FwdZcLocalVariance(double a,double b,double T1,double T2) const = 0;

    /// Covariance in [a,b] of Zc(.,T1)/Zc(.,U1) and Zc(.,T2)/Zc(.,U2)
    virtual double FwdZcLocalCovariance(double a,double b,double T1,double U1,double T2,double U2,std::vector<double>& vars) const = 0;

    /// Covariance in [0,t] between the state variable and Zc(.,T)
    virtual double StateZcCovariance(double t,double T) const = 0;

    /// Variance spread of Zc(.,T1) in [0,t1] and and Zc(.,T2) in [0,t2]
    virtual double ZcVarianceSpread(double t1,double t2,double T1,double T2) const = 0;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

