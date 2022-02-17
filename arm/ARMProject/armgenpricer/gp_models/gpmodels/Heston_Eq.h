/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Heston_Eq.h
 *
 *  \brief 
 *
 *	\author  E. Ezzine
 *	\version 1.0
 *	\date May 2006
 */


#ifndef _INGPMODELS_HESTON_EQ_H
#define _INGPMODELS_HESTON_EQ_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"

#include "gpmodels/EqFxBase.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/Heston_ModelParams.h"

/// forward declaration
class ARM_ZeroCurve;

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_ModelParamsHeston_Eq
// \brief Class for model parameters of Heston Eq
//-----------------------------------------------------------------------------

class ARM_ModelParamsHeston_Eq : public ARM_ModelParams_Eq, public ARM_Heston_ModelParams
{
private:
	void ValidateModelParams() const;

public:
	ARM_ModelParamsHeston_Eq( const ARM_ModelParamVector& params=ARM_ModelParamVector(), ARM_ZeroCurvePtr domCurve=ARM_ZeroCurvePtr(NULL), double spot=100.0 );
	ARM_ModelParamsHeston_Eq( const ARM_ModelParamsHeston_Eq& rhs );
	ASSIGN_OPERATOR(ARM_ModelParamsHeston_Eq)
	virtual ~ARM_ModelParamsHeston_Eq();
    

	/// How many factors?
    virtual size_t FactorCount() const { return 2; }
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {}
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};

	/// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_ModelParamsHeston_Eq(*this); };
    virtual string toString(const string& indent="",const string& nextIndent="") const;
};


class ARM_HestonModel_Eq :  public ARM_EqFxBase
{
public:

	ARM_HestonModel_Eq(const ARM_ZeroCurvePtr& zc, 
		ARM_ModelParamsHeston_Eq* modelParam);

	ARM_HestonModel_Eq( const ARM_HestonModel_Eq& rhs )
	:	ARM_EqFxBase(rhs) {}
	 
	ASSIGN_OPERATOR(ARM_HestonModel_Eq)
	virtual ~ARM_HestonModel_Eq(){}

	/// How many factors?
    virtual size_t FactorCount() const { return 2; }

	/// convention support : calendar + gap support
	virtual string GetSettlementCalendar(const string& modelName="") const;
	virtual double GetSettlementGap(const string& modelName="") const;
	virtual int GetType() const ;

	/// ARM Support
	virtual ARM_Object* Clone() const { return new ARM_HestonModel_Eq( *this ); }
	virtual string ExportShortName() const { return "LHEQM";}

	virtual ARM_PricingStatesPtr Init( const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos );
	virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;

	virtual bool SupportBackwardInduction() const {	return false; }
    virtual bool SupportForwardInduction()  const {	return true; }
	virtual bool SupportAnalyticMarginal()  const {	return false;}
    

	virtual ARM_VectorPtr DiscountFactor( 
		const string& curveName,
        double evalTime, 
		const std::vector<double>&  maturityTime,
        const ARM_PricingStatesPtr& states) const;

	/// Forward function
	virtual ARM_VectorPtr Forward(	const string& modelName, 
									double evalTime,
									double expiryTime,
									double settlementTime,
									double payTime,
									const ARM_PricingStatesPtr& states) const;
	/*
	/// Call function (vectorial strike version)
	virtual ARM_VectorPtr CallVectorial(
									const string& modelName,
									double evalTime,
									double expiryTime,
									double settlementTime,
									const std::vector<double>& strikePerState,
									int callPut,
									double payTime,
									const ARM_PricingStatesPtr& states,
									ARM_PricingContext* context) const;
	*/
	//NDC
	/// Call function (vectorial strike version)
	virtual ARM_VectorPtr CallVectorial(
									const string& modelName,
									double evalTime,
									const std::vector<double>& expiryTime,
									const std::vector<double>& settlementTime,
									const std::vector<double>& strikePerState,
									int callPut,
									const std::vector<double>& payTime,
									const ARM_PricingStatesPtr& states,
									ARM_PricingContext* context) const;

	/// Greek function (vectorial strike version)
	virtual ARM_VectorPtr GreekVectorial(
									const string& modelName,
									double evalTime,
									const std::vector<double>& expiryTime,
									const std::vector<double>& settlementTime,
									const std::vector<double>& strikesPerState,
									int callPut,
									const std::vector<double>& payTime,
									const string& greekType,
									const ARM_PricingStatesPtr& states,
									ARM_PricingContext* context=NULL) const;

	virtual double ComputeFwdAtTime( double evalTime  ) const;

	/// Monte-Carlo
	virtual void ModelStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps);

	virtual void NumMethodStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps );

	virtual void ModelStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances ) const;			

    virtual void NumMethodStateLocalGlobalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances, ARM_MatrixVector& variances ) const;	

	virtual void NumMethodStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances ) const;		

	virtual void NumMethodStateGlobalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& globalVariances ) const;		

	virtual bool NeedsToCholeskyDecomposeFactors( ) const { return false; }					

	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;

	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const;

	virtual double VarianceToTime(double var, double minTime=0.0, double maxTime=5*K_YEAR_LEN) const;
};



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
