/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file SFRM.h
 *
 *  \brief TargetFuncHW.h
 *	\author  E Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPMODELS_TARGETFUNCHW_H
#define _INGPMODELS_TARGETFUNCHW_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
#include "gpbase/env.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparamtype.h"
#include "gpinfra/typedef.h"
#include "gpnumlib/numfunction.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_PricingModel;

///////////////////////////////////////////////////////////////
/// \class VarDiffFunc
/// \brief
///  Function to solve using a solver monodim.
///////////////////////////////////////////////////////////////

class VarDiffFunc : public UnaryFuncWithDerivative<double,double>
{
private:    
    ARM_PricingModel* itsPricingModel;
    ARM_ModelParams*   itsParams;
    ARM_DbleToDbleFunctor* itsDerivative;
    ARM_ModelParamType::ParamNb 
                        itsParamType;
    DbleBinaryFunctor* itsBinaryFunc;

    size_t itsIndex;
    double itsTime;
    double itsTarget;

public:

    VarDiffFunc(ARM_ModelParams* params,
        ARM_PricingModel* model,
        ARM_ModelParamType::ParamNb paramType = ARM_ModelParamType::Volatility,
        DbleBinaryFunctor* Func = new ARM_BinFuncMinus );

    virtual ~VarDiffFunc();

    inline const ARM_PricingModel* GetPricingModel() const {return itsPricingModel;}
    void SetParamType(ARM_ModelParamType::ParamNb paramType ) {itsParamType = paramType;}
    void SetParams( size_t index, double target, double time ) {itsIndex = index, itsTarget = target, itsTime = time ;};
    double operator () ( double x ) const ;
    inline virtual ARM_DbleToDbleFunctor* Derivative() const { return itsDerivative; }
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
