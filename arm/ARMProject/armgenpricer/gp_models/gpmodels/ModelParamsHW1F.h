/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ModelParamsHW1F.h,v $
 * Revision 1.1  2003/10/13 07:51:48  jmprie
 * Initial revision
 *
 *
 */



/*! \file ModelParamsHW1F.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou JM Prie
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPMODELS_MODELPARAMSHW1F_H
#define _INGPMODELS_MODELPARAMSHW1F_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "modelparamshw.h"
#include "gpbase/port.h"
#include "gpcalib/typedef.h"			/// for calibParamVector
#include "gpinfra/modelparamtype.h"		/// for ARM_ModelParamType


CC_BEGIN_NAMESPACE( ARM )

class ARM_CurveModelParam;

//-----------------------------------------------------------------------------
// \class ARM_ModelParamsHW1F
// \brief Interface class for model parameters of the Hull & White 1F model
//-----------------------------------------------------------------------------
class ARM_ModelParamsHW1F : public ARM_ModelParamsHW 
{
private:
	int itsVolatilityType;

public:
	ARM_ModelParamsHW1F( const ARM_ModelParamsHW1F& rhs );
	ARM_ModelParamsHW1F( const ARM_ModelParamVector& params=ARM_ModelParamVector(), int volatilityType = ARM_ModelParamType::Volatility );
	virtual ~ARM_ModelParamsHW1F();
    ARM_ModelParamsHW1F& operator = (const ARM_ModelParamsHW1F& rhs);

	/// Dates at which modelparams values change
	virtual ARM_GP_VectorPtr ModelParamsTimeSteps() const;

	/// How many factors?
	virtual size_t FactorCount() const { return 1; }

    /// Coefficient of the state variable (for Zc closed form formula)
    virtual double BetatT(double t,double T) const = 0;

    /// Drift from a to b of the state variable
    virtual double StateLocalDrift(double a,double b) const = 0;

    /// Variance in [a,b] of the state variable
    virtual double StateLocalVariance(double a,double b,double c) const = 0;

    /// Variance in [a,b] of Zc(.,T1)/Zc(.,T2)
    virtual double FwdZcLocalVariance(double a,double b,double T1,double T2) const = 0;

    /// Covariance in [a,b] of Zc(.,T1)/Zc(.,U1) and Zc(.,T2)/Zc(.,U2)
    virtual double FwdZcLocalCovariance(double a,double b,double T1,double U1,double T2,double U2,std::vector<double>& vars) const = 0;

    /// Covariance in [0,t] between the state variable and Zc(.,T)
    virtual double StateZcCovariance(double t,double T) const = 0;

    /// Variance spread of Zc(.,T1) in [0,t1] and and Zc(.,T2) in [0,t2]
    virtual double ZcVarianceSpread(double t1,double t2,double T1,double T2) const = 0;

// FIXMEFRED: mig.vc8 (22/05/2007 18:09:59):cast
	inline int GetVolatilityType() const { return itsVolatilityType; }
	inline void SetVolatilityType(int volatilityType) { itsVolatilityType = volatilityType; }

	bool IsLn() const;

	/// Static functions for integration of two Hull & White 1F standard model
    /// i.e. for cst MRS mean-reverting and perfectly correlated models (brownian correlation = 1)
	static double HW1FStateCovariance( const ARM_ModelParamsHW1F* lhs, const ARM_ModelParamsHW1F* rhs, double a, double b, double T );
    static double HW1FZcCovariance( const ARM_ModelParamsHW1F* lhs, const ARM_ModelParamsHW1F* rhs, double a, double b, double T1, double T2=-1 );
    static double HW1FStateZcCovariance( const ARM_ModelParamsHW1F* lhs, const ARM_ModelParamsHW1F* rhs, double a, double b, double T, double U );
    static double HW1FEqFxZcCovariance( const ARM_CurveModelParam& eqFxVol, const ARM_ModelParamsHW1F* rhs, double a, double b, double T );
    static double HW1FEqFxStateCovariance( const ARM_CurveModelParam& eqFxVol, const ARM_ModelParamsHW1F* rhs, double a, double b, double T );
	static double HW1FStateCovarianceWithVec(
		const std::vector<double>& sigmaTimes1,
		const std::vector<double>& sigmaValues1,
		const std::vector<double>& timesTimes2,
		const std::vector<double>& sigmaValues2,
		double MRSValue1,
		double MRSValue2,
		double a,
		double b,
		double T);

	static double BetatT(const ARM_ModelParam& mrsParam,double t,double T);
	static double DerivBetatT(const  ARM_ModelParam& mrsParam, double t,double T);
};


//-----------------------------------------------------------------------------
// \class ARM_ModelParamsHW1FStd
// \brief Model parameters of the Hull & White 1F standard model
//  (variable sigma and constant MRS)
//-----------------------------------------------------------------------------
class ARM_ModelParamsHW1FStd : public ARM_ModelParamsHW1F
{
private:
    double Variance(double a,double b,double scale) const;
	double VarianceWithZeroMeanReversion(double a,double b) const;

public:
	ARM_ModelParamsHW1FStd( const ARM_ModelParamsHW1FStd& rhs );
	ARM_ModelParamsHW1FStd( const ARM_ModelParamVector& params=ARM_ModelParamVector(), int volatilityType = ARM_ModelParamType::Volatility );
	virtual ~ARM_ModelParamsHW1FStd();
    ARM_ModelParamsHW1FStd& operator = (const ARM_ModelParamsHW1FStd& rhs);

    virtual double BetatT(double t,double T) const
	{ return ARM_ModelParamsHW1F::BetatT( GetModelParam(ARM_ModelParamType::MeanReversion),t,T); }

    virtual double StateLocalDrift(double a,double b) const;
    virtual double StateLocalVariance(double a,double b,double c) const;
    virtual double FwdZcLocalVariance(double a,double b,double T1,double T2) const;
    virtual double FwdZcLocalCovariance(double a,double b,double T1,double U1,double T2,double U2,std::vector<double>& vars) const;
    virtual double StateZcCovariance(double t,double T) const;
    virtual double ZcVarianceSpread(double t1,double t2,double T1,double T2) const;

	// Standard ARM object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const;

	/// -------------- methods not implemented 
	/// -------------- but forced to be redefined for safety!
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0){};
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ){};
};

//-----------------------------------------------------------------------------
// \class ARM_ModelParamsHW1FExt
// \brief Model parameters of the Hull & White 1F extended model
//  (variable sigma and MRS)
//-----------------------------------------------------------------------------
class ARM_ModelParamsHW1FExt : public ARM_ModelParamsHW1F
{
private:
    double Lambda(double a,double b) const;
    double ScaledBetatT(double a,double b,double scale,double& expLambda) const;
    double Variance(double a,double b) const;
    std::vector<double>& VecBetatT(double t,double T) const;
    double IntegLambdaBeta(double a,double b,double T,std::vector<double>& betatTs) const;

public:
	ARM_ModelParamsHW1FExt( const ARM_ModelParamsHW1FExt& rhs );
	ARM_ModelParamsHW1FExt( const ARM_ModelParamVector& params=ARM_ModelParamVector(), int volatilityType = ARM_ModelParamType::Volatility );
	virtual ~ARM_ModelParamsHW1FExt();
    ARM_ModelParamsHW1FExt& operator = (const ARM_ModelParamsHW1FExt& rhs);

    double GetLambda( double a, double b) const {return Lambda(a,b);}
    double GetScaledBetatT(double a,double b,double scale,double& expLambda) const{ 
                                    return ScaledBetatT(a,b,scale,expLambda);}

    virtual double BetatT(double t,double T) const;
    virtual double StateLocalDrift(double a,double b) const;
    virtual double StateLocalVariance(double a,double b,double c) const;
    virtual double FwdZcLocalVariance(double a,double b,double T1,double T2) const;
    virtual double FwdZcLocalCovariance(double a,double b,double T1,double U1,double T2,double U2,std::vector<double>& vars) const;
    virtual double StateZcCovariance(double t,double T) const;
    virtual double ZcVarianceSpread(double t1,double t2,double T1,double T2) const;

   	// Standard ARM object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const;

	/// -------------- methods not implemented 
	/// -------------- but forced to be redefined for safety!
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0){};
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ){};
};

//-----------------------------------------------------------------------------
// \class ARM_CalibParamsHW1FExt
// \brief Model parameters of the Hull & White 1F extended model
//  (variable sigma and MRS)
//-----------------------------------------------------------------------------
class ARM_CalibParamsHW1FExt : public ARM_ModelParamsHW1FExt
{
private:
    double ScaledBetatTToMRS(double a,double b,double scale,double& expLambda) const;
    double VarianceToVolatility(double a,double b) const;

public:
	ARM_CalibParamsHW1FExt( const ARM_CalibParamsHW1FExt& rhs );
	ARM_CalibParamsHW1FExt( const ARM_ModelParamVector& params=ARM_ModelParamVector(), int volatilityType = ARM_ModelParamType::Volatility );
	virtual ~ARM_CalibParamsHW1FExt();
    ARM_CalibParamsHW1FExt& operator = (const ARM_CalibParamsHW1FExt& rhs);
    double ExpLambdaIntegral(double a,double b,
                                std::vector<double>& times,
                                std::vector<double>& values) const;
    double SquaredIntegral(double a,double b,std::vector<double>& times, 
                           std::vector<double>& values) const;

    virtual double FwdZcLocalVariance(double a,double b,double T1,double T2) const;

    /// initialise the news parameters to calibrate
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0);
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 );
    void VarianceToSigma(ARM_ModelParamVector& CalibParamVector,ARM_PricingModel* model);
    virtual void AdviseCurrentCalibSecIndex(size_t index){};

	// Standard ARM object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const { return "ARM_CalibParamsHW1FExt"; }
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

