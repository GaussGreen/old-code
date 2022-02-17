/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file HW2F.h
 *
 *  \brief 
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPMODELS_HW2F_H
#define _INGPMODELS_HW2F_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"

#include "gpbase/port.h"

#include "HW.H"
#include "gpinfra/typedef.h"


#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )


//-----------------------------------------------------------------------------
// \class HW2FVarianceFunction
// \brief
//  Variance function to find the time to reach a given variance
//-----------------------------------------------------------------------------

class HW2FVarianceFunction
{
private:
    const ARM_ModelParams* itsModelParams;

public:
    HW2FVarianceFunction(const ARM_ModelParams* modelParams);
    virtual ~HW2FVarianceFunction();
    double operator () ( double x ) const;
};


class csecurity;

//-----------------------------------------------------------------------------
// \class ARM_HullWhite2F
// \brief
//  2 factors Hull & White pricing model for closed form,
//  backward and forward diffusion abilities
//-----------------------------------------------------------------------------

class ARM_HullWhite2F : public ARM_HullWhite
{
private:
	// for fast calibration purpose
	vector<ARM_ModelParamType::ParamNb>		itsCalibParam;
	vector<csecurity*>						itsPF;
	vector<size_t>							itsSize;

private :
    ARM_VectorPtr LogNorBondVanillaSwaption(
		const string& curveName,
		double evalTime,
		double swapResetTime,
		const ARM_GP_Vector& fixNotional,
		const ARM_GP_Vector& floatNotional,
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& floatResetTimes,
		const ARM_GP_Vector& floatStartTimes,
		const ARM_GP_Vector& floatEndTimes,
		const ARM_GP_Vector& floatIntTerms,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
        const ARM_PricingStatesPtr& states,
		bool isConstantNotional,
		bool isConstantSpread,
		bool isConstantStrike) const;

    ARM_VectorPtr IntegrationVanillaSwaption(
		const string& curveName,
		double evalTime,
		double swapResetTime,
		const ARM_GP_Vector& fixNotional,
		const ARM_GP_Vector& floatNotional,
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& floatResetTimes,
		const ARM_GP_Vector& floatStartTimes,
		const ARM_GP_Vector& floatEndTimes,
		const ARM_GP_Vector& floatIntTerms,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
        const ARM_PricingStatesPtr& states,
		bool isConstantNotional,
		bool isConstantSpread,
		bool isConstantStrike) const;

	ARM_VectorPtr IntegrationVanillaSwaptionNew(
		const string& curveName,
		double evalTime,
		double swapResetTime,
		const ARM_GP_Vector& fixNotional,
		const ARM_GP_Vector& floatNotional,
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& floatResetTimes,
		const ARM_GP_Vector& floatStartTimes,
		const ARM_GP_Vector& floatEndTimes,
		const ARM_GP_Vector& floatIntTerms,	
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
        const ARM_PricingStatesPtr& states,
		bool isConstantNotional,
		bool isConstantSpread,
		bool isConstantStrike) const;

public:
	ARM_HullWhite2F( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params=NULL, const ARM_BoolVector& soFormulaFlags=ARM_BoolVector(2,true) );
	ARM_HullWhite2F(const ARM_HullWhite2F& rhs);
	virtual ~ARM_HullWhite2F();

    ARM_HullWhite2F& operator = (const ARM_HullWhite2F& rhs);

	virtual ARM_VectorPtr DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr VanillaSwaption(
		const string& curveName,
		double evalTime,
		double swapResetTime,
		const ARM_GP_Vector& fixNotional,
		const ARM_GP_Vector& floatNotional,
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& floatResetTimes,
		const ARM_GP_Vector& floatStartTimes,
		const ARM_GP_Vector& floatEndTimes,
		const ARM_GP_Vector& floatIntTerms,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
		const ARM_PricingStatesPtr& states,
		bool isConstantNotional = true,
		bool isConstantSpread = true,
		bool isConstantStrike = true) const;

	virtual ARM_VectorPtr  VanillaSpreadOptionLet(
		const string& curveName,
		double evalTime,
		int callPut,
		double startTime,
		double endTime,
		double resetTime,
		double payTime,
		double payPeriod,
		double notional,
		double coeffLong,
		double coeffShort,
		const ARM_GP_Vector& strikes,
		double swapLongFloatStartTime,
		double swapLongFloatEndTime,
		const ARM_GP_Vector& swapLongFixPayTimes,
		const ARM_GP_Vector& swapLongFixPayPeriods,
		double swapShortFloatStartTime,
		double swapShortFloatEndTime,
		const ARM_GP_Vector& swapShortFixPayTimes,
		const ARM_GP_Vector& swapShortFixPayPeriods,
		const ARM_PricingStatesPtr& states) const;


	/// Variable Notional swaptions can be priced in HW2F
	virtual bool ClosedFormulaSwaptionFlag(
		bool isConstantNominal,
		bool isConstantSpread,
		bool isConstantstrike) const { return isConstantSpread&&isConstantstrike;}

private:
	// for variable notio swaption (numerical integral)
	// called by method VanillaSwaption
	virtual ARM_VectorPtr VariableNotionalSwaption(
		const string& curveName,
		double evalTime,
		double swapResetTime,
		const ARM_GP_Vector& fixNotional,
		const ARM_GP_Vector& floatNotional,
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& floatResetTimes,
		const ARM_GP_Vector& floatStartTimes,
		const ARM_GP_Vector& floatEndTimes,
		const ARM_GP_Vector& floatIntTerms,	
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
		const ARM_PricingStatesPtr& states,
		bool isConstantNotional = true ,
		bool isConstantSpread = true ,
		bool isConstantStrike = true ) const;

public:
    // Give local drifts and variances w.r.t. a given schedule
    virtual void IntegratedLocalDrifts(
		const ARM_GP_Vector& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

    virtual void NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void NumMethodStateGlobalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& globalVariances ) const;

	virtual void ModelStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

    virtual double VarianceToTime(double var,double minTime,double maxTime) const;
	
	ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;

	/// function for the generic tree (2nd generation)
	virtual void VolatilitiesAndCorrelations( const ARM_GP_Vector& timeSteps, 
		ARM_GP_MatrixPtr& vols,
		ARM_GP_MatrixPtr& d1Vols,
		ARM_GP_MatrixPtr& correls,
		bool linearVol = true) const;

	/// function for the generic tree (2nd generation)
    virtual void EulerLocalDrifts(const ARM_GP_Vector& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	virtual bool NeedsToCholeskyDecomposeFactors( ) const { return true; }

	virtual ARM_BoolVector NeedMCIntegProcess() const { return ARM_BoolVector(2, true); };

	/// -------------- methods not implemented 
	/// -------------- but forced to be redefined for safety!
     /// Calibration purpose
    virtual void Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter){};
    virtual void PreProcessing(ARM_ModelFitter& modelFitter){};
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter){};
    virtual void AdviseCurrentCalibSecIndex(size_t index, ARM_ModelFitter& modelFitter){};
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter){};
	virtual void ValidateCalibMethod(ARM_CalibMethod& calibMethod);
	virtual bool ValidateModelParams(const ARM_ModelParams& params) const;

	/// method to advise the break point times
    virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb=0 );

	/////////// FAST CALIBRATION
	inline ARM_ModelParamType::ParamNb	GetParamType(size_t n) const					{return itsCalibParam[n];};
	inline void							StoreParamType(ARM_ModelParamType::ParamNb pt)	{itsCalibParam.push_back(pt);};
	inline size_t						GetSize(size_t n) const							{return itsSize[n];};
	inline void							StoreSize(size_t k)								{itsSize.push_back(k);}
	inline csecurity*					GetSec(size_t n,size_t index) const				{return itsPF[n*itsSize[0]+index];};
	inline void							StoreSec(csecurity* sec)						{itsPF.push_back(sec);};
	inline size_t						GetNbStored() const								{return itsSize.size();};

	void								FreeFastCalib();

	inline bool							IsPresentCalibParam(ARM_ModelParamType::ParamNb pt)
	{
		for (size_t i=0;i<GetNbStored();i++)
			if (GetParamType(i)== pt)
				return true;
		return false;
	}
	inline bool							IsMissingCalibParam(ARM_ModelParamType::ParamNb pt)
	{
		for (size_t i=0;i<GetNbStored();i++)
			if (GetParamType(i)== pt)
				return false;
		return true;
	}

    /// Standard ARM object support
	virtual ARM_Object* Clone() const;
	virtual string ExportShortName() const { return "LHW2M";}
	virtual string toString(const string& indent="",const string& nextIndent="") const;
};

////////////////// FOR FAST CALIBRATION
struct ARM_VanillaSwaptionArg;
struct ARM_VanillaSpreadOptionArg;

class csecurity
{
protected:
	ARM_PricingModel* m_pmodel;
	double m_reset;
	double m_price;
	double m_A;
	double m_B;
	double m_C;
	double m_D;
	int	   m_callPut;
public:
	csecurity(ARM_PricingModel* pmodel, double price){m_pmodel = pmodel;m_price = price;};
	virtual ~csecurity(){};

	void setMarketPrice(double mktprice){m_price = mktprice;};
	ARM_GP_Vector* getRateWeight(double tStart,double tEnd,const ARM_GP_Vector& payTimes,const ARM_GP_Vector& periods,double mrs1,double mrs2);
	ARM_GP_Vector* getChangeNumWeight(double tPay,double tStart,double tEnd,const ARM_GP_Vector& payTimes,const ARM_GP_Vector& periods,double mrs1,double mrs2);
	double GetA(){ return m_A;}
	double GetB(){ return m_B;}
	double GetC(){ return m_C;}
	double GetD(){ return m_D;}
	double GetResetTime(){ return m_reset;}
	virtual void BuildCoeffs(const ARM_GP_Vector& lv, const ARM_GP_Vector& lv0, const ARM_GP_Vector& lv1) = 0;

	/// VNS extension
	ARM_GP_Vector* getRateWeight(double t_start,const ARM_GP_Vector& floatPayTimes,const ARM_GP_Vector& floatNotionals,
		const ARM_GP_Vector& fixPayTimes,const ARM_GP_Vector& fixPayNotionals,const ARM_GP_Vector& fixPeriods,
		double mrs1,double mrs2);
};

class cswaption:public csecurity
{
private:
	double m_fwd;
	double m_level;
	double m_strike;
	double m_w1;
	double m_w2;
public:
	cswaption(ARM_PricingModel* model, double price, ARM_VanillaSwaptionArg* sw, double mrs1, double mrs2);
	virtual ~cswaption(){};
	virtual void BuildCoeffs(const ARM_GP_Vector& lv, const ARM_GP_Vector& lv0, const ARM_GP_Vector& lv1);
};

class cspreadoption:public csecurity
{
private:
	double m_fwdLong;
	double m_fwdShort;
	double m_levLong;
	double m_levShort;

	double m_level;
	double m_strike;

	double m_wLong1;
	double m_wLong2;
	double m_wShort1;
	double m_wShort2;

	double m_cLong1;
	double m_cLong2;
	double m_cShort1;
	double m_cShort2;
public:
	cspreadoption(ARM_PricingModel* pmodel, double price, ARM_VanillaSpreadOptionArg* sw, double mrs1, double mrs2,ARM_Currency* ccy);
	virtual ~cspreadoption(){};
	virtual void BuildCoeffs(const ARM_GP_Vector& lv, const ARM_GP_Vector& lv0, const ARM_GP_Vector& lv1);
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
