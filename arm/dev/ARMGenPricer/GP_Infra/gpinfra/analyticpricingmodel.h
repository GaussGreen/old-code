/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: AnalyticPricingModel.h,v $
 * Revision 1.1  2004/04/16 17:33:13  jmprie
 * Initial revision
 *
 *
 *
 */


/*! \file AnalyticPricingModel.h
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date April 2004
 */


#ifndef _INGPINFRA_ANALYTICPRICINGMODEL_H
#define _INGPINFRA_ANALYTICPRICINGMODEL_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
#include "gpbase/env.h"

#include "pricingmodel.h"
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM ) /// macro for namespace ... define namespace only if supported


///////////////////////////////////////////////////////
/// \class ARM_AnalyticPricingModel
/// \brief
/// This abstract class is the standard
/// interface for a pricing model only usable with
/// analytical formulas (induction impossible) 
///////////////////////////////////////////////////////
class ARM_AnalyticPricingModel : public ARM_PricingModel
{
public:
	ARM_AnalyticPricingModel(const ARM_ZeroCurvePtr& zc=ARM_ZeroCurvePtr(NULL), const ARM_ModelParams* params=NULL); 
	ARM_AnalyticPricingModel(const ARM_AnalyticPricingModel& rhs);
	virtual ~ARM_AnalyticPricingModel();
    ARM_AnalyticPricingModel& operator = (const ARM_AnalyticPricingModel& rhs);

    /// Non relevant pure virtual fcts...
    ///...for calibration
    virtual void Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter);
    virtual void PreProcessing(ARM_ModelFitter& modelFitter);
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter);
    virtual void AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter);
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter);
    virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb=0 );
    ///...for pricing process   
	virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);
	virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;
	virtual ARM_GP_Vector* ComputeModelTimes( const ARM_TimeInfoPtrVector& timeInfos );
// FIXMEFRED: mig.vc8 (24/05/2007 10:42:18):cast
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const {return static_cast<ARM_VectorPtr>(NULL);}
	virtual bool SupportBackwardInduction() const;
	virtual bool SupportForwardInduction() const;
	virtual bool SupportAnalyticMarginal() const;
	virtual void SetNumMethod(const ARM_NumMethodPtr& numMethodPtr );

    ///...for pricing computation
	virtual void IntegratedLocalDrifts(const ARM_GP_Vector& timeSteps,ARM_GP_MatrixPtr& relativeDrifts,ARM_GP_MatrixPtr& absoluteDrifts) const;
	virtual void ModelStateLocalVariances( const ARM_GP_Vector& timeSteps,	ARM_MatrixVector& localVariances ) const;
	virtual void NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps,	ARM_MatrixVector& localVariances ) const;
	virtual bool NeedsToCholeskyDecomposeFactors( ) const;
    virtual double VarianceToTime(double var,double minTime=0.0,double maxTime=5*K_YEAR_LEN) const;
	
	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;
	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states,int timeIndex) const;

    /// Standard ARM object support
	virtual string toString(const string& indent="",const string& nextIndent="") const { return indent + string("ARM_AnalyticPricingModel : abstract class !"); }
};



CC_END_NAMESPACE()

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
