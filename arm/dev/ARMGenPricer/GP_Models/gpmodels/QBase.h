/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file QBase.h
 *
 *  \brief base class for the q model
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date July 2004
 */


#ifndef _INGPMODELS_QBASE_H
#define _INGPMODELS_QBASE_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpinfra/pricingmodelir.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_ModelParamsQBase;

//-----------------------------------------------------------------------------
// \class ARM_QModelBase
// \brief
//  base Class for the Q Model
//-----------------------------------------------------------------------------

class ARM_QModelBase : public ARM_PricingModel
{
private:
    ARM_GP_Vector* itsPrecomputedFwds;
    void CopyNoCleanUp(const ARM_QModelBase& rhs);
    void CleanUp();

public:
	ARM_QModelBase(const ARM_ZeroCurvePtr& zc, const ARM_ModelParams& params);
	ARM_QModelBase(const ARM_QModelBase& rhs);
    ARM_QModelBase& operator = (const ARM_QModelBase& rhs);
	virtual ~ARM_QModelBase();

	///Calibration purpose
	virtual void ValidateCalibMethod(ARM_CalibMethod& calibMethod);

	/// accessors
	ARM_GP_Vector* GetPrecomputedFwds() const { return itsPrecomputedFwds; }

	/// common functions
	virtual double ComputeFwdAtTime( double evalTime ) const = 0;
	virtual double MappingFunction( double x, double x0, double q0 ) const = 0;
	virtual double MappingFunctionDerivative( double x, double r0, double q0 ) const;

	/// numerical method bool flags + validation
    virtual bool SupportBackwardInduction() const {	return true; }
    virtual bool SupportForwardInduction()  const {	return true; }
	virtual bool SupportAnalyticMarginal()  const {	return true;}
    virtual bool NeedArrowDebreuPrices() const { return true; }

    /// Default initialisation of the model
	virtual ARM_PricingStatesPtr Init( const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos );
	virtual void PostInit();
	virtual void TreeStatesToModelStates( ARM_PricingStatesPtr& states, int timeIndex ) const;
	virtual void MCModelStatesFromToNextTime( ARM_PricingStatesPtr& states,int timeIndex ) const;
    virtual bool NeedsToCholeskyDecomposeFactors() const {return false;}
	virtual bool IsMeanRevertingCompatible() const { return true;}

    /// Calibration purpose
    virtual void Re_InitialiseCalibParams( ARM_ModelFitter& modelFitter ){};
    virtual void PreProcessing( ARM_ModelFitter& modelFitter ) {};
    virtual void PostProcessing( const ARM_ModelFitter& modelFitter ) {};
    virtual void AdviseCurrentCalibSecIndex( size_t index,ARM_ModelFitter& modelFitter){};
    virtual void AdviseCurrentCalib( ARM_ModelFitter& modelFitter ){};

	// This function tells to MC to return or not an integrated process 
	virtual ARM_BoolVector NeedMCIntegProcess() const;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
