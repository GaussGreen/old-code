/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file QBaseIR.h
 *
 *  \brief base class for the q model
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date July 2004
 */


#ifndef _INGPMODELS_QBASEIR_H
#define _INGPMODELS_QBASEIR_H

#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "QBase.h"
#include "gpinfra/pricingmodelir.h"

#include <map>
CC_USING_NS(std,map)

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_QModelBaseIR
// \brief Interest rate class for the Q Model
//  base Class for the Q Model
//-----------------------------------------------------------------------------

class ARM_QModelBaseIR: public ARM_QModelBase,
						public ARM_PricingFunctionIR
{

protected:

    /// Nested class that gives target & external stochastic payoff component
    /// to allow calibration of zero-coupon bond approximated formula
    class ARM_QDfTarget
    {
    private :
        double itsTargetPayoff;
        ARM_GP_VectorPtr itsPayoffFactor;

    public:
        ARM_QDfTarget(double targetPayoff,const ARM_GP_VectorPtr& payoffFactor)
            : itsTargetPayoff(targetPayoff), itsPayoffFactor(payoffFactor) {}

        inline const double GetTargetPayoff() const { return itsTargetPayoff; }
        inline void SetTargetPayoff(double targetPayoff) { itsTargetPayoff=targetPayoff; }
        inline const ARM_GP_VectorPtr GetPayoffFactor() const { return itsPayoffFactor; }
        inline void SetPayoffFactor(const ARM_GP_VectorPtr& payoffFactor) { itsPayoffFactor=payoffFactor; }
    };

private:

    /// To indicate is Q Model is an H&W like (except for swaption
    /// where true H&W closed from is not used because of possible use of basis curve)
	bool itsDegenerateInHW;


    typedef map< double, ARM_QDfTarget > ARM_QDfTargetMap;
    typedef ARM_QDfTargetMap::iterator  ARM_QDfTargetMapIter;
    typedef ARM_QDfTargetMap::const_iterator  ARM_QDfTargetMapConstIter;

    /// A map to manage efficiency & simultaneously DFs
    ARM_QDfTargetMap itsQDfTargetMap;


    void CopyNoCleanUp(const ARM_QModelBaseIR& rhs);

public:

	ARM_QModelBaseIR(const ARM_ZeroCurvePtr& zc, const ARM_ModelParams& params, bool degenerateInHW = false );
	ARM_QModelBaseIR(const ARM_QModelBaseIR& rhs);
    ARM_QModelBaseIR& operator = (const ARM_QModelBaseIR& rhs);
	virtual ~ARM_QModelBaseIR();

	inline bool IsDegenerateInHW() const { return itsDegenerateInHW; }

    inline const ARM_QDfTarget* GetTargetPayoff(double T) const;
    inline void SetTargetPayoff(double T, double targetPayoff, const ARM_GP_VectorPtr& payoffFactor);
    void ResetTargetPayoffs() { itsQDfTargetMap.clear(); }

	/// common functions
	virtual double ComputeFwdAtTime( double evalTime ) const;
    double DecayFactor( double a, double b, double T, const ARM_ModelParams* modelParams ) const;
    double SuperDecayFactor( double a, double b, double T, const ARM_ModelParams* modelParams ) const;
    double PartialIntegratedDrift( int i, double mrs, double a, double b, double T, const ARM_ModelParams* modelParams ) const;
    double IntegratedDrift( double t,double T, const ARM_ModelParams* modelParams ) const;
    double VolZc(double t, double T) const;

	virtual ARM_GP_MatrixPtr MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const;

    /// Pricing of forward Libor
	virtual ARM_VectorPtr Libor( 
		const string& curveName,
        double evalTime,
		double fwdStartTime,
        double fwdEndTime,
		double period,
        double fwdResetTime,    // for convexity adjustment...
        double payTime,         //... in derived classes
        const ARM_PricingStatesPtr& states) const;
    
	/// function to avoid computing twice the fixLegAnnuity
	/// if you want to keep the annuity value, make sure you have cloned it before!
	virtual ARM_VectorPtr SwapRateInPlaceWithComputedAnnuity(
		const string& curveName, 
		double evalTime,
		double floatStartTime, 
		double floatEndTime, 
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& fwdStartTimes,
        const ARM_GP_Vector& fwdEndTimes,
        const ARM_GP_Vector& fwdPayPeriods,
		const ARM_GP_Vector& floatPayTimes,
        const ARM_GP_Vector& floatPayPeriods,
        const ARM_GP_Vector& margin,
        bool isDbleNotional,
		const ARM_VectorPtr& PreComputedAnnuity,
		const ARM_PricingStatesPtr& states) const;


	/// Defauft annuity provided but may be redefined
	virtual ARM_VectorPtr Annuity(
		const string& curveName, 
        double evalTime,
		const ARM_GP_Vector& fixPayTimes,
        const ARM_GP_Vector& fixPayPeriods,
        const ARM_PricingStatesPtr& states) const;

	/// Defauft Swap
	virtual ARM_VectorPtr NPVSwap(
		const string& curveName, 
        double evalTime,
		double floatStartTime,
        double floatEndTime, 
		const ARM_GP_Vector& fixPayTimes,
        const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& fwdStartTimes, 
		const ARM_GP_Vector& fwdEndTimes, 
		const ARM_GP_Vector& fwdPayPeriods, 
		const ARM_GP_Vector& floatPayTimes, 
		const ARM_GP_Vector& floatPayPeriods, 
		const ARM_GP_Vector& margin,
        bool isDbleNotional,
		const ARM_GP_Vector& FixNotional,
		const ARM_GP_Vector& FloatNotional,
		const ARM_GP_Matrix& strikesPerState,
        int payRec,
        const ARM_PricingStatesPtr& states) const;


	virtual ARM_GP_MatrixPtr NPVSwapLeg(
		const string& curveName, 
		double evalTime,
		const ARM_GP_Vector& fwdStartTimes, 
		const ARM_GP_Vector& fwdEndTimes, 
		const ARM_GP_Vector& fwdPayPeriods, 
		const ARM_GP_Vector& PayTimes, 
		const ARM_GP_Vector& PayPeriods, 
		const ARM_GP_Vector& margin, 
		const ARM_GP_Vector& notional, 
		const ARM_PricingStatesPtr& states) const ;

	virtual ARM_GP_MatrixPtr NPVFixLeg(
		const string& curveName, 
		double evalTime,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& FixNotional,
		const ARM_GP_Matrix& strikesPerState,
		int   payRec,
		const ARM_PricingStatesPtr& states) const;


    virtual double ImpliedVol(const ARM_VanillaArg& arg) const;

	virtual double MappingFunction( double x, double x0, double q0 ) const;
    virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb =0);
	virtual int GetType() const;

	/// global constant
	static const double RatePeriod;
	static const double RateUpperBound;
	static const double Rate0LowerBound;
};

inline const ARM_QModelBaseIR::ARM_QDfTarget* ARM_QModelBaseIR::GetTargetPayoff(double T) const
{
    ARM_QDfTargetMapConstIter found = itsQDfTargetMap.find(T);
    if(found != itsQDfTargetMap.end())
        return &found->second;
    else
        return NULL;
}

inline void ARM_QModelBaseIR::SetTargetPayoff(double T, double targetPayoff, const ARM_GP_VectorPtr& payoffFactor)
{
    ARM_QDfTargetMapIter found = itsQDfTargetMap.find(T);
    if(found != itsQDfTargetMap.end())
    {
        /// Just replace it
        found->second.SetTargetPayoff(targetPayoff);
        found->second.SetPayoffFactor(payoffFactor); /// ptr may be shared many times !
    }
    else
    {
        ;
        std::pair< double,ARM_QDfTarget > value(T,ARM_QDfTarget(targetPayoff,payoffFactor));
        std::pair< ARM_QDfTargetMapIter,bool > result = itsQDfTargetMap.insert(value);
        if(!result.second)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": can't insert data for DF calibration");
    }
}

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
