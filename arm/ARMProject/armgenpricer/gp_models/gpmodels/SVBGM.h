/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 */

#ifndef _INGPMODELS_MODELSVBGM_H
#define _INGPMODELS_MODELSVBGM_H

#include "gpinfra/modelparams.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparamsvec.h"
#include "gpbase/curve.h"
#include "typedef.h"

#include "gpmodels/modelparamssvbgm.h"

#include "gpinfra/pricingmodelir.h"
#include "gpcalib/numerical.h"

#include "gpclosedforms/sabrimpliedvol.h"
#include "gpclosedforms/smile_sabr.h"
#include "nag.h"
#include "nage04.h"
#include "gpclosedforms/optimization1.h"
#include "gpnumlib/levmarq.h"

CC_BEGIN_NAMESPACE( ARM )

#define DefaultAlpha 999.

class ARM_ShiftedSABR
{
protected:
	
	double		itsShift;
	int			itsShiftSign;
	double		itsATMVol;
	double		itsAlpha;
	double		itsNu;
	double		itsRho;

	double		itsResetTime;
	double		itsFwd;
	double		itsStrike;
	int			itsSens;

	bool		itsIsInit;

public:
	ARM_ShiftedSABR()
	{
		itsIsInit	= false;
		itsAlpha	= DefaultAlpha;
	}

	virtual ~ARM_ShiftedSABR()
	{
	}

public:

	double		GetShift() const		{return itsShift;};
	double		GetAlpha() const		{return itsAlpha;};
	double		GetNu() const			{return itsNu;};
	double		GetRho() const			{return itsRho;};

	double		GetFwd() const			{return itsFwd;};
	double		GetResetTime() const	{return itsResetTime;};
	double		GetStrike() const		{return itsStrike;};
	int			GetSens() const			{return itsSens;};

	void		build(double resetTime, double forward, double strike, int sens, double atmvol, double alpha,
					double rho, double nu, double shift);

	void		calibre(double resetTime, double forward, double atmvol, std::vector<double> * mktVols, std::vector<double> * strikes,
					double initshift, double initrho, double initnu);

	double		price();

	double		price(double resetTime, double forward, double strike, int sens, double atmvol, double alpha, 
					double rho, double nu, double shift);

	void		SetShift(double shift);
	void		SetAlpha(double alpha);
	void		SetRho(double rho);
	void		SetNu(double nu);

	void		SetStrike(double strike)	{itsStrike = strike;};
	void		SetSens(int sens)			{itsSens = sens;};

private:

	void		computeAlpha(double rho, double nu);
	double		sabrvol(double fwd, double stk, double rho, double nu);

	class ObjectiveFunction : 
		public Optimization_ObjectiveFuntion,
		public ARM_LEVMARQFunc 
	{
	private:
		ARM_ShiftedSABR *	_this;
		bool				itsCalibShift;
		bool				itsCalibRho;
		bool				itsCalibNu;
		double				itsEPSforDer;
		std::vector<double> *		itsStrikes;
		std::vector<double> *		itsMktVols;

	public:
		ObjectiveFunction(ARM_ShiftedSABR * shsabr, bool calibshift, bool calibrho, bool calibnu, double eps, std::vector<double> * strikes, std::vector<double> * mktvols)
		{
			_this			= shsabr;
			itsCalibShift	= calibshift;
			itsCalibRho		= calibrho;
			itsCalibNu		= calibnu;
			itsEPSforDer	= eps;
			itsStrikes		= strikes;
			itsMktVols		= mktvols;
		}

	public:
		void NAG_CALL operator ()(Integer m, Integer n, double x[], double f[], double fjac[],
				Integer tdfjac, Nag_Comm * comm)
		{
			int k = 0, i = 0;

			if(itsCalibShift)	_this->SetShift(x[k++]);
			if(itsCalibRho)		_this->SetRho(x[k++]);
			if(itsCalibNu)		_this->SetNu(x[k]);

			ComputeImpliedVols(f);

			k = 0;

			// derivée par rapport au shift

			if(itsCalibShift)
			{
				_this->SetShift(x[k] + itsEPSforDer);
				ComputeImpliedVols(fjac,i,n);
				Computefjac(m,f,fjac,i,n);
				_this->SetShift(x[k]);
				_this->SetAlpha(DefaultAlpha);
				k++; i++;
			}

			if(itsCalibRho)
			{
				_this->SetRho(x[k] + itsEPSforDer);
				ComputeImpliedVols(fjac,i,n);
				Computefjac(m,f,fjac,i,n);
				_this->SetRho(x[k]);
				_this->SetAlpha(DefaultAlpha);
				k++; i++;
			}

			if(itsCalibNu)
			{
				_this->SetNu(x[k] + itsEPSforDer);
				ComputeImpliedVols(fjac,i,n);
				Computefjac(m,f,fjac,i,n);
				_this->SetNu(x[k]);
			}

			_this->SetAlpha(DefaultAlpha);
		}

		void operator()(double p[], double hx[], int m, int n, void * adata = NULL) const
		{
			int k = 0, i = 0;

			if(itsCalibShift)	_this->SetShift(p[k++]);
			if(itsCalibRho)		_this->SetRho(p[k++]);
			if(itsCalibNu)		_this->SetNu(p[k]);
			
			int size = itsStrikes->size();

			for(k = 0; k < size; k++)
			{
				_this->SetStrike((*itsStrikes)[k]);
				_this->SetSens(_this->GetStrike() > _this->GetFwd() ? 1 : -1);

				double price		= _this->price();
				double vol			= VanillaImpliedVol_BS(_this->GetFwd(), _this->GetStrike(), _this->GetResetTime(), price, _this->GetSens());
				
				hx[k] = vol - (*itsMktVols)[k];
			}	
		}

	private:

		void ComputeImpliedVols(double f[], int k0 = 0, int kdecal = 1)
		{
			int k, size = itsStrikes->size();

			for(k = 0; k < size; k++)
			{
				_this->SetStrike((*itsStrikes)[k]);
				_this->SetSens(_this->GetStrike() > _this->GetFwd() ? 1 : -1);

				double price		= _this->price();
				double vol			= VanillaImpliedVol_BS(_this->GetFwd(), _this->GetStrike(), _this->GetResetTime(), price, _this->GetSens());
				
				f[k*kdecal + k0] = vol;
			}
		}

		void Computefjac(int m, double f[], double fjac[], int k0, int kdecal)
		{
			for(int k = 0; k < m; k++)
			{
				fjac[k*kdecal+k0] = (fjac[k*kdecal+k0] - f[k]) / itsEPSforDer;
			}
		}
	};
};

inline void ARM_ShiftedSABR::build(double resetTime, double forward, double strike, int sens, double atmvol, 
								   double alpha, double rho, double nu, double shift)
{
	itsResetTime	= resetTime;
	itsFwd			= forward;
	itsStrike		= strike;
	itsSens			= sens;
	itsATMVol		= atmvol;
	itsAlpha		= alpha;
	itsRho			= rho;
	itsNu			= nu;
	itsShiftSign	= shift < 0. ? -1 : 1;
	itsShift		= fabs(shift) < 0.001 ? 0.001 * itsShiftSign : CC_Min<double>(CC_Max<double>(shift,-5.),1.);

	itsIsInit		= true;
}

inline void ARM_ShiftedSABR::SetShift(double shift)
{
	itsShiftSign	= shift < 0. ? -1 : 1;
	itsShift		= fabs(shift) < 0.001 ? 0.001 * itsShiftSign : CC_Min<double>(CC_Max<double>(shift,-5.),1.);
}

inline void ARM_ShiftedSABR::SetRho(double rho)
{
	itsRho	= rho > 0.999 ? 0.999 : rho < -0.999 ? -0.999 : rho;
}

inline void ARM_ShiftedSABR::SetNu(double nu)
{
	itsNu	= nu < 1e-8 ? 1e-8 : nu;
}

inline void ARM_ShiftedSABR::SetAlpha(double alpha)
{
	itsAlpha	= alpha < 1e-8 ? 1e-8 : alpha;
}

inline double ARM_ShiftedSABR::price(double resetTime, double forward, double strike, int sens, double atmvol,
									 double alpha, double rho, double nu, double shift)
{
	build(resetTime, forward, strike, sens, atmvol, alpha, rho, nu, shift);

	return price();
}

inline double ARM_ShiftedSABR::price()
{
	if(itsIsInit == false) return 0.;

	double fwd	= itsFwd / fabs(itsShift);
	double stk	= itsShift < 0. ? fwd * (1. + fabs(itsShift)) - itsStrike : itsStrike - (itsShift - 1.) * fwd;
	double vol	= sabrvol(fwd, stk, itsRho, itsNu);

	return BS(fwd, stk, itsResetTime, vol, itsSens * itsShiftSign);
}

inline void ARM_ShiftedSABR::computeAlpha(double rho, double nu)
{
	double a, b, c, d, atmvol = itsATMVol * fabs(itsShift);
	
	a	= 0.25 * itsResetTime * rho * nu;

	b	= 1. + itsResetTime * (2. - 3. * rho * rho) * nu * nu / 24.;

	c	= - atmvol;

	d	= b * b - 4. * a * c;

	if(d < 0. || fabs(a) < 1e-12)
	{
		itsAlpha = atmvol;
	}
	else
	{
		double rac1 = (- b + sqrt(d)) / (2. * a);
		double rac2 = (- b - sqrt(d)) / (2. * a);

		itsAlpha = (rac1 > 0. ? rac1 : rac2 > 0. ? rac2 : atmvol);
	}
}

inline double ARM_ShiftedSABR::sabrvol(double fwd, double strike, double rho, double nu)
{
	if(fabs(itsAlpha - DefaultAlpha) < K_DOUBLE_TOL) computeAlpha(rho, nu);

	if(fabs(fwd - strike) < 0.0009)
	{
		double zeta = 2. * (nu / itsAlpha) * (fwd - strike) / (fwd + strike);
		double fac	= 1. + (rho * nu * itsAlpha /4. + (2. - 3. * rho * rho) * nu * nu / 24.) * itsResetTime;

		return (1. - zeta / 12. * (6. * rho + zeta * (-2. + 3. * rho * rho))) * itsAlpha * fac;
	}
	else
	{
		double fav		=	(fwd + strike)/2.;

		double Cf		=	fav;
		double dCf		=	1.;
		double d2Cf		=	0.;

		// rapport C'(fav)/C(fav)
		double g1		=	1. / fav;
		// rapport C"(fav)/C(fav)
		double g2		=	0.;
		// intégrale de (1/C) entre f et K
		double ln		=	log(fwd / strike);

		double ksi		=	(nu / itsAlpha) * (fwd - strike) / fav;
		double xksi		=	log((sqrt(1.- 2. * rho * ksi + ksi * ksi) + ksi - rho) / (1.- rho));

		double coeff;

		coeff	= (2. * g2 - g1 * g1 + 1. / (fav * fav)) * itsAlpha * itsAlpha * Cf * Cf / 24.;
		coeff	+= rho * nu * itsAlpha * g1 * Cf / 4.;
		coeff	+= (2.- 3. * rho * rho) * nu * nu / 24.;

		coeff	= 1.+ coeff * itsResetTime;

		return itsAlpha * log(fwd / strike) * ksi * coeff / (ln * xksi);		
	}
}

class ARM_SVBGM : public ARM_PricingModelIR 
{
protected:
	
	bool							itsCalibrationStatus;
	bool							itsCalibratedStatus;
	
	// les échéanciers
	std::vector<double>					itsResetTimes;
	std::vector<double>					itsStartTimes;
	std::vector<double>					itsEndTimes;
	std::vector<double>					itsDelta;
	std::vector<double>					itsFwdRate;
	std::vector<double>					itsFromShiftedToRate;
	
	bool							itsAllowInterpol;	// interpolation des zc

	ARM_GP_Matrix					itsWeight;

	bool							itsComputeDriftAdj;
	
	ARM_VectorVector				itsEigenValues;
	ARM_MatrixVector				itsModelLocalRealVar;

	ARM_VanillaSecDensityPtrVector	itsCalibSecDensities;

	bool							itsProxyStatus;
	int								itsNbEffectiveReset;

public:
	ARM_SVBGM( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params = NULL, bool AllowInterpol = false, bool ComputeDrift = true, bool Proxy = false);
	
	ARM_SVBGM( const ARM_SVBGM& rhs);
	
	virtual ~ARM_SVBGM();

	virtual ARM_Object*				Clone() const	{ return new ARM_SVBGM(*this);};
	
	ARM_SVBGM&						operator = (const ARM_SVBGM& rhs);

public:
	// méthodes virtuelles pures
	virtual bool					ValidateModelParams(const ARM_ModelParams& params) const;
	virtual void					ValidateCalibMethod(ARM_CalibMethod& calibMethod);
	virtual void					SetNumeraire(const ARM_NumerairePtr& numerairePtr);

	
    virtual void					PreProcessing(ARM_ModelFitter& modelFitter);
    virtual void					PostProcessing(const ARM_ModelFitter& modelFitter);
    virtual void					AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter);
    virtual void					AdviseCurrentCalib(ARM_ModelFitter& modelFitter) {};
	virtual void					Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter) {};
	
    
	virtual bool					SupportBackwardInduction() const {	return false;}
    virtual bool					SupportForwardInduction()  const {	return true;}
	virtual bool					SupportAnalyticMarginal()  const {	return false;}

	virtual void					PostInit();
	
	virtual void					setCalibrationStatus(bool calibrationStatus);

	virtual void					setCalibratedStatus(bool calibratedStatus);

	virtual ARM_PricingStatesPtr	FirstPricingStates( size_t bucketSize ) const;

	virtual ARM_PricingStatesPtr	Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);

	virtual ARM_PricingStatesPtr	Induct(ARM_PricingStatesPtr& states,double toTime);

	virtual std::vector<double>&			ComputeModelTimes(const ARM_TimeInfoPtrVector& timeInfos);
	
	virtual ARM_VectorPtr			ComputeNumeraireTimes(const ARM_TimeInfoPtrVector& timeInfos) const;

	virtual void					AdviseBreakPointTimes(const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb = 0);

	//////////////////////////////
	//
	//		pour passer les échéanciers !
	//
	//////////////////////////////

	virtual void					setNumericalModelFitter(ARM_NumericalModelFitter*);

	//////////////////////////////
	//
	//		pour le monte carlo
	//
	//////////////////////////////

	virtual void					ModelStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps);

	virtual void					NumMethodStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps );

	virtual void					ModelStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances ) const;

	virtual void					NumMethodStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances ) const;		

	virtual void					NumMethodStateGlobalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& globalVariances ) const;		

	virtual bool					NeedsToCholeskyDecomposeFactors( ) const { return false; }					

	virtual void					MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;

	// Pas utile !!
	virtual void					TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const;
	virtual double					VarianceToTime(double var, double minTime=0.0, double maxTime=5*K_YEAR_LEN) const;


	///////////////////////
	//
	//		les prix
	//
	///////////////////////

	virtual ARM_VectorPtr			DiscountFactor(	const string& curveName, double evalTime, double maturityTime, const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr			Libor(const string& curveName,
										double evalTime,
										double fwdStartTime,
										double fwdEndTime,
										double period,
										double resetTime,
										double payTime,
										const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr			SwapRate(const string& curveName, 
										double evalTime,
										double floatStartTime,
										double floatEndTime, 
										const std::vector<double>& fixPayTimes,
										const std::vector<double>& fixPayPeriods,
										const std::vector<double>& fwdStartTimes,
										const std::vector<double>& fwdEndTimes,
										const std::vector<double>& fwdPayPeriods, 
										const std::vector<double>& floatPayTimes,
										const std::vector<double>& floatPayPeriods,
										const std::vector<double>& margin,
										bool isDbleNotional,
										const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr			VanillaCaplet(const string& curveName, 
										double evalTime,
										double payTime, 
										double period,
										double payNotional,
										double fwdResetTime, 
										double fwdStartTime,
										double fwdEndTime,
										double fwdPeriod,
										const std::vector<double>& strikesPerState,
										int capFloor,
										const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr			VanillaSwaption(const string& curveName,
										double evalTime,
										double swapResetTime,
										const std::vector<double>& fixNotional,
										const std::vector<double>& floatNotional,
										double floatStartTime,
										double floatEndTime,
										const std::vector<double>& floatResetTimes,
										const std::vector<double>& floatStartTimes,
										const std::vector<double>& floatEndTimes,
										const std::vector<double>& floatIntTerms,
										const std::vector<double>& fixPayTimes,
										const std::vector<double>& fixPayPeriods,
										const ARM_GP_Matrix& strikesPerState,
										int callPut,
										const ARM_PricingStatesPtr& states,
										bool isConstantNotional = true,
										bool isConstantSpread = true,
										bool isConstantStrike = true) const;

	virtual ARM_VectorPtr			VanillaSpreadOptionLet(const string& curveName,
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
										const std::vector<double>& strikes,
										double swapLongFloatStartTime,
										double swapLongFloatEndTime,
										const std::vector<double>& swapLongFixPayTimes,
										const std::vector<double>& swapLongFixPayPeriods,
										double swapShortFloatStartTime,
										double swapShortFloatEndTime,
										const std::vector<double>& swapShortFixPayTimes,
										const std::vector<double>& swapShortFixPayPeriods,
										const ARM_PricingStatesPtr& states) const;
	
	/// Standard ARM object support
	virtual string					toString(const string& indent="",const string& nextIndent="") const;

protected:
	
	double			getTerminalTime() const;

	bool			DoesResetTimeExist(double time) const;

	void			computeWeights();

	bool			IsOnSamePath(size_t i, size_t j) const;

	ARM_VectorPtr	DiscountFactorNoInterpol(const string& curveName,
											double evalTime, 
											double maturityTime,
											const ARM_PricingStatesPtr& states) const;

	ARM_VectorPtr	ForwardDiscountFactor(const string& curveName,
											double evalTime, 
											double startTime,
											double endTime,
											const ARM_PricingStatesPtr& states) const;

	ARM_VectorPtr	ForwardDiscountFactorFromIdx(const string& curveName,
											double evalTime,
											size_t IdxFrom,
											size_t IdxTo,
											size_t modelNb,
											const ARM_PricingStatesPtr& states) const;

	void			kthLocalCalibration(int k, double& shift, double& alpha, double& rho, double& nu);

	void			getkthStrikesAndVols(int k, std::vector<double>& strikes, std::vector<double>& vols, double& atmvol);
};

inline void ARM_SVBGM::PostInit()
{
}

inline void ARM_SVBGM::setCalibrationStatus(bool calibrationStatus)
{
	itsCalibrationStatus = calibrationStatus;
}

inline void ARM_SVBGM::setCalibratedStatus(bool calibratedStatus)
{
	itsCalibratedStatus = calibratedStatus;
}

inline void ARM_SVBGM::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   
	ARM_THROW( ERR_INVALID_ARGUMENT,"ARM_SVBGM::TreeStatesToModelStates :  not implemented ARM_SVBGM Model!");
}

inline double ARM_SVBGM::VarianceToTime(double var,double minTime,double maxTime) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT,"ARM_SVBGM::VarianceToTime not implemented ARM_SVBGM Model!");

}

inline double ARM_SVBGM::getTerminalTime() const
{
	return itsEndTimes.Elt( itsEndTimes.size() -1 );
}

inline bool ARM_SVBGM::DoesResetTimeExist(double time) const
{
	return ExistsInVector( itsResetTimes, time );
}

inline bool ARM_SVBGM::IsOnSamePath(size_t i, size_t j) const
{
	return itsWeight(i,j) == 1.;
}

inline int IdxFromValue(const std::vector<double>& dateVector, double date, double tol)
{
	if(fabs(tol) < K_DOUBLE_TOL)
	{
		if(ExistsInVector(dateVector, date))
			return IdxFromValue(dateVector, date);
		else
			return -1;
	}
	else
	{
		int k, size = dateVector.size();

		for(k = 0; k < size; k++)
		{
			if(fabs(date - dateVector[k]) < tol) return k;
		}
	}

	return -1;
}

CC_END_NAMESPACE()

#endif
