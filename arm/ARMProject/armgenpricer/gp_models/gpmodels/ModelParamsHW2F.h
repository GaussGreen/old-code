/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ModelParamsHW2F.h,v $
 * Revision 1.1  2003/10/21 07:51:48  jmprie
 * Initial revision
 *
 *
 */



/*! \file ModelParamsHW2F.h
 *
 *  \brief 
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPMODELS_MODELPARAMSHW2F_H
#define _INGPMODELS_MODELPARAMSHW2F_H
/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"

#include "ModelParamsHW.h"
#include "gpbase/port.h"
#include "gpbase/gpmatrixtriangular.h"

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_ModelParamsHW2F
// \brief
//  Class for model parameters of the Hull & White 2F model
//  Mean reversion speeds and correlation are assumed constants
//  sigma1(t) = sigma(t)
//  sigma2(t) = VolatilityRatio * sigma(t)
//  MRS2 = MRS1 + MeanReversionSpread
//-----------------------------------------------------------------------------
class ARM_ModelParamsHW2F : public ARM_ModelParamsHW 
{
public:
	ARM_ModelParamsHW2F( const ARM_ModelParamsHW2F& rhs );
	ARM_ModelParamsHW2F( const ARM_ModelParamVector& params=ARM_ModelParamVector() );
	virtual ~ARM_ModelParamsHW2F();
    ARM_ModelParamsHW2F& operator = (const ARM_ModelParamsHW2F& rhs);

	/// Dates at which modelparams values change
	virtual ARM_GP_VectorPtr ModelParamsTimeSteps() const;

	/// How many factors?
	virtual size_t FactorCount() const { return 2; }

    /// Coefficient of the state variable (for Zc closed form formula)
    virtual std::vector<double>& BetatT(double t,double T) const = 0;

    /// Drift from a to b of the state variable
    virtual std::vector<double>& StateLocalDrift(double a,double b) const = 0;

	virtual ARM_GP_MatrixPtr StateLocalDriftVec( const std::vector<double>& timeSteps ) const = 0;

    /// Variance in [a,b] of the state variable
    virtual ARM_GP_TriangularMatrix* StateLocalVariance(double a,double b) const = 0;

    /// Variance in [a,b] of Zc(.,T1)/Zc(.,T2)
    virtual double FwdZcLocalVariance(double a,double b,double T1,double T2) const = 0;

    /// Covariance in [a,b] of Zc(.,T1)/Zc(.,U1) and Zc(.,T2)/Zc(.,U2)
    virtual double FwdZcLocalCovariance(double a,double b,double T1,double U1,double T2,double U2,std::vector<double>& vars) const = 0;

    /// Covariance in [0,t] between the state variable and Zc(.,T)
    virtual double StateZcCovariance(double t,double T) const = 0;

    /// Variance spread of Zc(.,T1) in [0,t1] and and Zc(.,T2) in [0,t2]
    virtual double ZcVarianceSpread(double t1,double t2,double T1,double T2) const = 0;

	/// Static function (MRSs are assumed constant !)
	static std::vector<double>& BetatT(const ARM_ModelParam& mrsParam,const ARM_ModelParam& mrsSpreadParam,double t,double T);


};

//-----------------------------------------------------------------------------
// \class ARM_ModelParamsHW2FStd
// \brief
//  Class for model parameters of the Hull & White 2F model
//  Mean reversion speeds and correlation are assumed constants
//  sigma1(t) = sigma(t)
//  sigma2(t) = VolatilityRatio * sigma(t)
//  MRS2 = MRS1 + MeanReversionSpread
//-----------------------------------------------------------------------------
class ARM_ModelParamsHW2FStd : public ARM_ModelParamsHW2F 
{
private:
    std::vector<double>& Variance(double a,double b,std::vector<double>& scale) const;

public:
	ARM_ModelParamsHW2FStd( const ARM_ModelParamsHW2FStd& rhs );
	ARM_ModelParamsHW2FStd( const ARM_ModelParamVector& params=ARM_ModelParamVector() );
	virtual ~ARM_ModelParamsHW2FStd();
    ARM_ModelParamsHW2FStd& operator = (const ARM_ModelParamsHW2FStd& rhs);

   /// Coefficient of the state variable (for Zc closed form formula)
    virtual std::vector<double>& BetatT(double t,double T) const
	{ return ARM_ModelParamsHW2F::BetatT( GetModelParam(ARM_ModelParamType::MeanReversion),GetModelParam(ARM_ModelParamType::MeanReversionSpread),t,T); }

    /// Drift from a to b of the state variable
    virtual std::vector<double>& StateLocalDrift(double a,double b) const;

	virtual ARM_GP_MatrixPtr StateLocalDriftVec( const std::vector<double>& timeSteps ) const;

    /// Variance in [a,b] of the state variable
    virtual ARM_GP_TriangularMatrix* StateLocalVariance(double a,double b) const;

    /// Variance in [a,b] of Zc(.,T1)/Zc(.,T2)
    virtual double FwdZcLocalVariance(double a,double b,double T1,double T2) const;

    /// Covariance in [a,b] of Zc(.,T1)/Zc(.,U1) and Zc(.,T2)/Zc(.,U2)
    virtual double FwdZcLocalCovariance(double a,double b,double T1,double U1,double T2,double U2,std::vector<double>& vars) const;

    /// Covariance in [0,t] between the state variable and Zc(.,T)
    virtual double StateZcCovariance(double t,double T) const;

    /// Variance spread of Zc(.,T1) in [0,t1] and and Zc(.,T2) in [0,t2]
    virtual double ZcVarianceSpread(double t1,double t2,double T1,double T2) const;


	// Standard ARM object support
	virtual ARM_Object* Clone() const;
	virtual string toString() const;

	/// -------------- methods not implemented 
	/// -------------- but forced to be redefined for safety!
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0){};
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ){};
};

//-----------------------------------------------------------------------------
// \class ARM_ModelParamsHW2FExt
// \brief
//  Class for model parameters of the Hull & White 2F model
//  Mean reversion speeds is assumed constant
//  sigma1(t) = sigma(t)
//  sigma2(t) = VolatilityRatio(t) * sigma(t)
//  MRS2 = MRS1 + MeanReversionSpread
//-----------------------------------------------------------------------------
class ARM_ModelParamsHW2FExt : public ARM_ModelParamsHW2F 
{
public:

	ARM_ModelParamsHW2FExt( const ARM_ModelParamsHW2FExt& rhs );
	ARM_ModelParamsHW2FExt( const ARM_ModelParamVector& params=ARM_ModelParamVector() );
	virtual ~ARM_ModelParamsHW2FExt();
    ARM_ModelParamsHW2FExt& operator = (const ARM_ModelParamsHW2FExt& rhs);

   /// Coefficient of the state variable (for Zc closed form formula)
    virtual std::vector<double>& BetatT(double t,double T) const
	{ return ARM_ModelParamsHW2F::BetatT( GetModelParam(ARM_ModelParamType::MeanReversion),GetModelParam(ARM_ModelParamType::MeanReversionSpread),t,T); }

    /// Drift from a to b of the state variable
    virtual std::vector<double>& StateLocalDrift(double a,double b) const;

	virtual ARM_GP_MatrixPtr StateLocalDriftVec( const std::vector<double>& timeSteps ) const;

    /// Variance in [a,b] of the state variable
    virtual ARM_GP_TriangularMatrix* StateLocalVariance(double a,double b) const;

    /// Variance in [a,b] of Zc(.,T1)/Zc(.,T2)
    virtual double FwdZcLocalVariance(double a,double b,double T1,double T2) const;

    /// Covariance in [a,b] of Zc(.,T1)/Zc(.,U1) and Zc(.,T2)/Zc(.,U2)
    virtual double FwdZcLocalCovariance(double a,double b,double T1,double U1,double T2,double U2,std::vector<double>& vars) const;

    /// Covariance in [0,t] between the state variable and Zc(.,T)
    virtual double StateZcCovariance(double t,double T) const;

    /// Variance spread of Zc(.,T1) in [0,t1] and and Zc(.,T2) in [0,t2]
    virtual double ZcVarianceSpread(double t1,double t2,double T1,double T2) const;


	// Standard ARM object support
	virtual ARM_Object* Clone() const;
	virtual string toString() const;

	/// -------------- methods not implemented 
	/// -------------- but forced to be redefined for safety!
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0){};
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ){};

	//fast calibration purpose
	std::vector<double>& Variance(double a,double b,std::vector<double>& scale) const;
};




CC_END_NAMESPACE()





#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

