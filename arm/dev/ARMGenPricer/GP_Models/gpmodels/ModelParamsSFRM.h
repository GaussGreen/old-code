/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelparamsSFRM.h
 *
 *  \brief general class for model params in SFRM!
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */


#ifndef _INGPMODELS_MODELPARAMSSFRM_H
#define _INGPMODELS_MODELPARAMSSFRM_H

#include "gpbase/port.h"

#include "gpinfra/modelparams.h"
#include "gpinfra/curvemodelparam.h"

#include "gpcalib/typedef.h"

#include <vector>
CC_USING_NS(std,vector)


/// forward declaration
class ARM_IRIndex;

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_SFRM;

struct ARM_VanillaSwaptionArg;
struct ARM_VanillaCapDigitalArg;
struct ARM_VanillaArg;
struct MeanRevFunc;


//-----------------------------------------------------------------------------
// \class ARM_ModelParamsSFRM
// \brief Interface class for model parameters of the SFRM model
//-----------------------------------------------------------------------------
class ARM_ModelParamsSFRM : public ARM_ModelParams 
{
/// standard private section
private:
	void UseShiftModelParam(const ARM_PricingModel& model);
	/// validation of the model params
	void ValidateModelParams() const;
	void UpdateRealizedCorrel();					
													
	
	size_t itsFactorsNb;								/// nb of factor for the multi_dimension SFRM Model
	ARM_IRIndex* itsIRIndex;							/// SFRM ir index (cloned by the constructor)
	ARM_VectorPtrVector itsMuCoeff;					    /// mu coeff for the vol Swap/FRA relationship
    ARM_VectorPtr itsOneMu;	                            /// Local vector to store itsMuCoeff[index] for each GenSecurity
    ARM_GP_Vector* itsFwdValues;                            /// Local vector to store Fwds values.
	ARM_ModelParamType::ParamNb itsCorrelType;
	void UpdateCorrelation( ARM_ModelParamType::ParamNb correlType );

/// protected for easy access by derived classes!
protected:
	/// for position in ModelParams
    inline ARM_Curve* GetVolCurve() const { return ((ARM_CurveModelParam& )GetModelParam(ARM_ModelParamType::Volatility)).GetCurve(); }

public:
	ARM_ModelParamsSFRM( const ARM_ModelParamVector& params, ARM_IRIndex* index, size_t factorsNb = 1);
	ARM_ModelParamsSFRM( const ARM_ModelParamsSFRM& rhs);
	ARM_ModelParamsSFRM& operator=( const ARM_ModelParamsSFRM& rhs );
	virtual ~ARM_ModelParamsSFRM();

    /// Initial function to cumpte FwdsValues 
    /// to swap from beta to shift or vice versa
    void PreComputeFwds(const ARM_SFRM& model,const ARM_Portfolio& shiftConvPort);
    void ConvertToShiftParam(const ARM_SFRM& model,const ARM_Portfolio& shiftConvPort);
    void ConvertToBetaParam(const ARM_SFRM& model,const ARM_Portfolio& shiftConvPort);

	/// pricing function
    ARM_GP_Vector* InstantaneousVolatility(double t, double T) const;  
    double IntegratedVol(double s, double t, double T)const ;
	double IntegratedVariance(double s, double t, double T) const;
	double IntegratedCovariance(double s, double t, double T1, double T2) const;
	double IntegratedCovariance(double s, double t, double T1, double T2, double T3) const;
    double IntegratedCorrelation( double s, double t, double T1, double T2) const;
    double InstantaneousCorrelation(double T1, double T2) const;
	ARM_VectorPtr IntegratedVolPerFactor(double s, double t, double T);
	ARM_VectorPtr IntegratedVariancePerFactor(double s, double t, double T);
	ARM_VectorPtr IntegratedTreeStatesVariancePerFactor(double s, double t);
    ARM_GP_Vector* InterpolateCorrelation(double T) const;
    ARM_VectorPtr FwdsVarCorrelCoeffs(double T_1, double T_2);

    virtual double VarianceToTime(double var,double minTime,double maxTime) const = 0;

	double ShiftValue(double t);
    double BetaValue(double t);

    /// function to handle the mean reversion
    /// currently supported is constant mean reversion
	double IntegrateSquaredExpMRV(double time1, double time2 ) const;
    double ExpFunction( double time ) const;

	ARM_VectorPtr ComputeVolSwapvolFRA( const ARM_VanillaSwaptionArg& arg, 
		const ARM_PricingModel& model, bool isConstantNotional = true, bool useFixFrequency = false);

	ARM_VectorPtr ComputeVolSensiVolFra( const ARM_VanillaSwaptionArg& arg, 
		const ARM_PricingModel& model );

	ARM_VectorPtr ComputeVolSwapvolFRAUnderSwapProba( const ARM_VanillaSwaptionArg& arg, 
		const ARM_PricingModel& model, const ARM_VanillaSwaptionArg& argProba );

	ARM_VectorPtr IntegratedSwaptionVolVec(double s, double t, const ARM_VanillaSwaptionArg& arg,
        const ARM_VectorPtr Mu);
	void DumpSwaptionVolsAndTimes( 
		const ARM_VanillaSwaptionArg& arg,
		ARM_VectorPtr times,
		ARM_VectorPtr volatilities);

	ARM_VectorPtr IntegratedCapVolVec(double s, double t, const ARM_VanillaCapDigitalArg& arg);
	double LocalVolatity(double s, double t, const ARM_VanillaArg& arg, const ARM_VectorPtr Mu);
	double AverageShift(const ARM_GP_Vector& times );

	/// Discretization Dates of the model params
	virtual ARM_GP_VectorPtr ModelParamsTimeSteps() const;

	/// function that are different in DIAG and ROW model
    virtual double VolatilityFunction(double t, double T) const = 0;
	virtual double IntegratedLocalVariance(double s, double t) const= 0;
	virtual double MaturityTerm(double T) const =0;
	virtual double MaturityTermSquared(double T) const =0;
	virtual double VolatilitySpotSquared(double s) const =0;

    /// initialise the news parameters to calibrate (not pure virtual to allow the use of 
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {}
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 );
    void PostProcessing();
	
	/// accessors
	inline const ARM_IRIndex* GetIRIndex() const                    { return itsIRIndex; }
	inline ARM_IRIndex* GetIRIndex()                                { return itsIRIndex; }
	virtual size_t FactorCount() const                              { return itsFactorsNb; }
    inline const ARM_VectorPtr GetOneMu() const                     { return itsOneMu;}
    void SetOneMu(ARM_VectorPtr OneMu)                              { itsOneMu = OneMu;}
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

