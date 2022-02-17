/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Mixture_Fx.h
 *
 *  \brief 
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date 17 August 2006
 */


#ifndef _INGPMODELS_MIXTURE_FX_H
#define _INGPMODELS_MIXTURE_FX_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"

#include "gpinfra/modelparams.h"
#include "gpinfra/pricingmodelir.h"

#include "gpmodels/EqFxBase.h"
#include "gpmodels/ModelParams_EqFxBase.h"


/// forward declaration
class ARM_ZeroCurve;

CC_BEGIN_NAMESPACE( ARM )


class ARM_ParamsMixture_Fx : public ARM_RootObject
{
public:
	ARM_ParamsMixture_Fx(
		const std::vector<double>& lags,
		const std::vector<double>& volATM,
		const std::vector<double>& decVol,
		const std::vector<double>& shift,
		const std::vector<double>& lambda,
		const string& interpolName);
		
	ARM_ParamsMixture_Fx(const ARM_ParamsMixture_Fx& rhs)
		:
	itsVolATM(rhs.itsVolATM),
	itsDecVol(rhs.itsDecVol),
	itsShift(rhs.itsShift),
	itsLambda(rhs.itsLambda)
	{
	}

	ASSIGN_OPERATOR(ARM_ParamsMixture_Fx)
	virtual ~ARM_ParamsMixture_Fx() {};

	virtual ARM_Object* Clone() const { return new ARM_ParamsMixture_Fx( *this ); }
	virtual string ExportShortName() const { return "LPMFX";}
	virtual string toString(const string& indent, const string& nextIndent) const;

	const std::vector<double>& GetLags() const { return itsLags; }
	const std::vector<double>& GetVolATM() const { return itsVolATM; }
	const std::vector<double>& GetDecVol() const { return itsDecVol; }
	const std::vector<double>& GetShift() const { return itsShift; }
	const std::vector<double>& GetLambda() const { return itsLambda; }
	const string& GetInterpolName() const { return itsInterpolName; }

private:
	std::vector<double> itsLags;
	std::vector<double> itsVolATM;
	std::vector<double> itsDecVol;
	std::vector<double> itsShift;
	std::vector<double> itsLambda;
	string itsInterpolName;
};


class ARM_MixtureModel_Fx :  public ARM_EqFxBase
{
private:
	/// sigma is sigma2(Richard's paper)
	CC_IS_MUTABLE bool itsIsSigmaInputed;

	void Validate() const;
public:
	typedef ARM_ModelParams_Fx_T<ARM_ModelParams> ARM_ModelParamsMixture_Fx;

	ARM_MixtureModel_Fx(
		const ARM_ZeroCurvePtr& zc,
		ARM_ModelParamsMixture_Fx* modelParam);

	ARM_MixtureModel_Fx(
		const ARM_ZeroCurvePtr& domZc,
		const ARM_ZeroCurvePtr& forZc,
		double spot,
		ARM_ParamsMixture_Fx* mixtureParams);
	
	ARM_MixtureModel_Fx( const ARM_MixtureModel_Fx& rhs )
	:	ARM_EqFxBase(rhs),
		itsIsSigmaInputed(rhs.itsIsSigmaInputed)
	{}
	ASSIGN_OPERATOR(ARM_MixtureModel_Fx)
	virtual ~ARM_MixtureModel_Fx(){};

	//  Update the density functor at expiryTime
	virtual void UpdateDensityFunctor(double fwd, double expiryTime, double tenor);

	// Call function (vectorial strike version)
	virtual ARM_VectorPtr CallVectorial(
		const string& modelName,
        double evalTime,
	    double expiryTime,
	    double settlementTime,
		const std::vector<double>& strikePerState,
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
		const std::vector<double>& strikePerState,
		double notional,
		int callPut,
	    double payTime,
		ARM_DigitType digitType,
		double epsilon,
        const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const;

	virtual ARM_PricingStatesPtr Init( const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos );

	virtual void ModelStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	 virtual double ImpliedVol(const ARM_VanillaArg& arg) const;

	/// function for the monte carlo
	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;

	/// ARM Support
	virtual ARM_Object* Clone() const { return new ARM_MixtureModel_Fx( *this ); }
	virtual string ExportShortName() const { return "LFXMX";}

	static double CallMixturePrice(double fwd, double strike, int callPut, double stdDev1, double stdDev2, double alpha, double lambda);
	static double DigitalMixturePrice(double fwd, double strike, int callPut, double stdDev1, double stdDev2, double alpha, double lambda);
	static double CalibVol2(double fwd, double strike, int callPut, double stdDev1, double stdDevATM, double alpha, double lambda);
	static void CalibMixture(
		double fwd, 
		double T,
		int callPut, 
		const vector<double>& K, 
		const vector<double>& vols,
		const vector<double>& decVol,
		const vector<double>& alpha,
		const vector<double>& lambda,
		vector<double>& outParams);
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
