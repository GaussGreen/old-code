/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Smiled_Fx.h
 *
 *  \brief 
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date June 2006
 */


#ifndef _INGPMODELS_SMILED_FX_H
#define _INGPMODELS_SMILED_FX_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "gpbase/vectormanip.h"

#include "gpinfra/pricingmodelir.h"

#include "gpmodels/EqFxBase.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/ModelParamsSmiledFRM.h"
#include "gpmodels/ProcessBuilderSmiledFRM.h"

#include "gpcalib/numerical.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_ModelParamsSmiled_Fx : public ARM_ModelParams_Fx, public ARM_ModelParamsSmiled 
{
public:
	ARM_ModelParamsSmiled_Fx(
		const ARM_ModelParamVector& params, 
		const ARM_ZeroCurvePtr& domCurve, 
		const ARM_ZeroCurvePtr& forCurve, 
		double spot,
		size_t dim, 
		CorrelType correlType, 
		double ReCorrel)
	:
	ARM_ModelParams_Fx(
		domCurve,
		forCurve,
		spot),
	ARM_ModelParamsSmiled(
		params, 
		dim, 
		correlType, 
		false,
		ReCorrel)
	{
	}

	ARM_ModelParamsSmiled_Fx(const ARM_ModelParamsSmiled_Fx& rhs) :
	ARM_ModelParams_Fx(rhs),
	ARM_ModelParamsSmiled(rhs)
	{
	}

	ASSIGN_OPERATOR(ARM_ModelParamsSmiled_Fx)
	virtual ~ARM_ModelParamsSmiled_Fx()
	{
	}

	virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {}
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};

	/// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_ModelParamsSmiled_Fx(*this); };
    virtual string toString(const string& indent="",const string& nextIndent="") const {
		string str;
		str += ARM_ModelParams_Fx::toString(indent, nextIndent); 
		str += ARM_ModelParamsSmiled::toString(indent, nextIndent);
		return str; 
	};
};

class ARM_2IRFXModel;

//-----------------------------------------------------------------------------
// \class ARM_Smiled_Fx
// \brief
//  In this class we define a multi factor smiled FX model:
//  S(t,Ti) = Fi(X_i(t))
// 
//  with dX_i(t)/X_i(t) = sigma_i(t)*dW_i(t)
//  and d<W_i(t),W_i(t)> = ro_ij*dt
//
//  Each FX forward is a function of gausssian process.
//  And the n gaussian process are simulated with p factors
//
//  This model will be used then in a 1RFX multi asset model
//  It could also be used in a standalone FX model (the domestic and foreign 
/// interest rates are assumed to be deterministic).
//-----------------------------------------------------------------------------

class ARM_SmiledModel_Fx :  public ARM_EqFxBase
{
private:
	ARM_PricingModelPtr itsIRModel;

	ARM_NumericalModelFitter*				itsNumericalModelFitter;	//assoc
	vector< ARM_ProcessBuilder* >			itsProcess;
	mutable ARM_VectorVector				itsEigenValues;

	bool									itsCalibrationStatus;
	bool									itsIsCalibrated;

	ARM_GP_Vector							itsResetDates;
	ARM_GP_Vector							itsSettlementDates;

	ARM_GP_Vector							itsTotalVar;

	size_t									itsTimeStepsNb;
	size_t									itsGridSize;
	double									itsStdDevNb;

	bool									itsSkipPDE;
	bool									itsWithRescalling;

	ARM_2IRFXModelPtr						its2IRFXModel;
public:

	inline double	getTerminalTime() const								{ return itsResetDates.Elt( itsResetDates.size() -1 ); }

	inline const ARM_GP_Vector& getResetTimes()	const					{ return itsResetDates; }
	inline void	 setResetTimes( const ARM_GP_Vector& resetTimes )		{ itsResetDates = resetTimes;	}

	inline const ARM_GP_Vector& getSettlementTimes()						{ return itsSettlementDates; }
	inline void setSettlementTimes( const ARM_GP_Vector& settlementTimes )	{ itsSettlementDates = settlementTimes; }

	inline bool		DoesResetDateExist( double date ) const				{ return ExistsInVector( itsResetDates, date );	}

	ARM_SmiledModel_Fx(const ARM_ZeroCurvePtr& zc,
		ARM_ModelParamsSmiled_Fx* modelParam,
		size_t timeStepNb,
		size_t grideSize,
		double stdDevNb,
		bool skipPDE,
		bool withRescalling,
		const ARM_2IRFXModel* Model2IRFX=NULL);
	
	ARM_SmiledModel_Fx( const ARM_SmiledModel_Fx& rhs );

	ASSIGN_OPERATOR(ARM_SmiledModel_Fx)
	virtual ~ARM_SmiledModel_Fx();

	void SetModel2IRFX(ARM_2IRFXModel* Model2IRFX) { its2IRFXModel = CreateClonedPtr(Model2IRFX);}

	/// Forward function
	virtual ARM_VectorPtr Forward(
		const string& modelName, 
        double evalTime,
	    double expiryTime,
	    double settlementTime,
	    double payTime,
        const ARM_PricingStatesPtr& states) const;

	// Call function (vectorial strike version)
	virtual ARM_VectorPtr CallVectorial(
		const string& modelName,
        double evalTime,
	    double expiryTime,
	    double settlementTime,
		const ARM_GP_Vector& strikePerState,
		int callPut,
	    double payTime,
        const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const;

	/// Default Digital call function (vectorial strike version)
	virtual ARM_VectorPtr DigitalVectorial(
		const string& modelName,
        double evalTime,
	    double expiryTime,
	    double settlementTime,
		const ARM_GP_Vector& strikePerState,
		double notional,
		int callPut,
	    double payTime,
        const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const;

	virtual void ModelStateLocalVariances( 
		const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void	NumMethodStateLocalVariances(	const ARM_GP_Vector& timeSteps,
												ARM_MatrixVector& localVariances ) const;		

	virtual void	NumMethodStateGlobalVariances( const ARM_GP_Vector& timeSteps,
												ARM_MatrixVector& globalVariances ) const;

	virtual void IntegratedLocalDrifts(
		const ARM_GP_Vector& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	virtual bool	SupportBackwardInduction() const {	return true;}
    virtual bool	SupportForwardInduction()  const {	return true;}
	virtual bool	SupportAnalyticMarginal()  const {	return true;}

	const ARM_GP_Vector& GetTotalVar() const { return itsTotalVar; }

	// numerical calibration
	virtual void	setNumericalModelFitter(ARM_NumericalModelFitter*);
	virtual void	setCalibrationStatus( bool calibrationStatus )		{ itsCalibrationStatus = calibrationStatus;};

	/// function for the monte carlo
	virtual ARM_PricingStatesPtr	FirstPricingStates( size_t bucketSize ) const;
	virtual size_t ModelStatesSize() const;
	virtual ARM_PricingStatesPtr Init( const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos );
	virtual ARM_PricingStatesPtr	Induct(ARM_PricingStatesPtr& states,double toTime);
	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;
	virtual ARM_GP_Vector* ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos );

	/// function for the generic tree
	virtual void VolatilitiesAndCorrelations( const ARM_GP_Vector& timeSteps, ARM_GP_MatrixPtr& vols, 
											ARM_GP_MatrixPtr& d1Vols, ARM_GP_MatrixPtr& correls, 
											bool linearVol = true) const;
	void EulerLocalDrifts( const ARM_GP_Vector& timeSteps, ARM_GP_MatrixPtr& relativeDrifts, 
						ARM_GP_MatrixPtr& absoluteDrifts) const;

	/// ARM Support
	virtual ARM_Object* Clone() const { return new ARM_SmiledModel_Fx( *this ); }
	virtual string ExportShortName() const { return "LSMFX";}
	virtual	string			toString(const string& indent="",const string& nextIndent="") const;

	/// To keep links with its IRModel
	virtual void UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel );
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
