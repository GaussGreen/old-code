/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file MSV1F.h
 *
 *  \brief 
 *
 *	\author  A. Triki
 *	\version 1.0
 *	\date November 2005
 */


#ifndef _INGPMODELS_HWSV1F_H
#define _INGPMODELS_HWSV1F_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"

#include "gpmodels/hwsv.h"

#include "gpnumlib/ran2.h"

#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

//// forward declaration
class ARM_HullWhite;

class ARM_RiccatiHWSV1F : public ARM_RiccatiHWSV
{
public:

	/// Nested structure to solve time dependent Riccati system
	struct ARM_SystemData
	{
		/// Model dependent csts (dependency to vol, VoV, & correl schedules + Im(u) vector)
		/// CstB = -exp(-Mrs*T0)/Lambda
		double itsTime;
		double itsYf;
		double itsMrs;				/// Lambda (rate MRS)
		double itsKappaTheta;		/// kappa (MRS of variance) * theta (long term variance)

		double itsModA;				/// -0.5*Nu²

		double itsModB1;			/// Kappa + Rho(t)*Nu(t)*Epsilon(t)
		double itsModB2;			/// Rho(t)*Nu(t)*Epsilon(t)*CstB
		double itsModB2x;			/// -Rho(t)*Nu(t)*Epsilon(t)*F0 (Heston formula)
		ARM_GP_Vector itsModB2y;	/// -Rho(t)*Nu(t)*Epsilon(t)*F0*Im(u) (Heston formula)

		double itsModC;				/// 0.5*F0²*Epsilon(t)²
		ARM_GP_Vector itsModCx;		/// 0.5*F0²*Epsilon(t)²*Im(u)²
		ARM_GP_Vector itsModCy;		/// 0.5*F0²*Epsilon(t)²*Im(u) (Heston formula)
	}; 

	/// Nested structure to solve Riccati equation with stepwise constant coefficients
	/// Variance mean reversion is linked by kappa(t) = lambda - rho(t)*nu(t)*epsilon(t)/lambda
	struct ARM_FunctionData
	{	
		double              itsVolOfVol_t;
		double				itsNu_2;
		double              itsVol_t;
		double              itsRho_t;
		double              itsKappa_t;
		double              itsTheta_t;
		double              itsDeltaX;
		double              itsDeltaY;
		double              itsRoot_S2_X; // s2
		double              itsRoot_S2_Y;

		/////// Temporary Variable To avoid recalculation when updating Ux and Uy
		double				itsTempB;
		double				itsTempC;

		/////// Coefficients of the Riccati Equation
		///////  A X² + B X + C
		double				itsA;		// -0.5*Nu²
		double				itsBx;
		double				itsBy;
		double				itsCx;
		double				itsCy;

		double				itsFactX;	// 0.5*Nu²/delta
		double				itsFactY;

		double				itsFactVolX;
		double				itsFactVolY;

		double				itsKappaA;			// 2*kappa/Nu²
		double				itsKappaBetaVol;	// kappa*Integ{ti->ti+1,vol(s)*exp(-lambda*(Te-s))}

		double				itsBeta;
		double				itsBetaVolInst;
		double				itsBetaVol;
		double				itsBetaVolDeltaX;
		double				itsBetaVolDeltaY;
		double				itsExpBetaVolDeltaX;
		double				itsExpBetaVolDeltaY;


		double				itsCstX;				//// piece-wise function to ensure continuity
		double				itsCstY;
		////// Add Functions to updates some specific data

	}; 

protected :
	const vector< ARM_SystemData >& itsSystemDatas;

public:

	ARM_RiccatiHWSV1F::ARM_RiccatiHWSV1F(const vector< ARM_SystemData >& systemDatas,size_t solverType=ARM_ODEFunc::RK5Adaptative,bool isStdFormula=false);
	ARM_RiccatiHWSV1F::ARM_RiccatiHWSV1F(const ARM_RiccatiHWSV1F& rhs);
	ARM_RiccatiHWSV1F& ARM_RiccatiHWSV1F::operator=(const ARM_RiccatiHWSV1F& rhs);
	virtual ~ARM_RiccatiHWSV1F() {}

    virtual string toString(const string& indent="", const string& nextIndent="") const {return string("ARM_RiccatiHWSV1F");}
    virtual ARM_Object* Clone() const {return new ARM_RiccatiHWSV1F(*this);}

	virtual size_t GetSystemSize() const {return itsSystemDatas.size();}
	virtual double GetYf(size_t i) const {return itsSystemDatas[i].itsYf;}
	virtual size_t GetNbPhis() const {return itsSystemDatas[0].itsModB2y.size();}

	virtual int ComputeNbSteps(const ARM_GP_Vector& solverParams);

	virtual void derivs(double t, ARM_GP_Vector* yt, ARM_GP_Vector* dyt) const;
};

struct ARM_RiccatiHWSV1F_SO : public ARM_RiccatiHWSV1F
{
	/// Nested structure to solve time dependent Riccati system
	struct ARM_SystemData
	{
		/// quadratic coefficients
		double itsModB;				/// -nu²(t)

		/// unit coefficients (Heston formula)
		double itsModD;				/// -mu0²*vol(t)²
		ARM_GP_Vector itsModDx;		/// { modD * Im(u) }
	}; 

	const vector< ARM_SystemData >& itsSystemDatasSO;

	ARM_RiccatiHWSV1F_SO::ARM_RiccatiHWSV1F_SO(const vector< ARM_RiccatiHWSV1F::ARM_SystemData >& systemDatas,
		const vector< ARM_SystemData >& systemDatasSO,bool isStdFormula=false)
		: ARM_RiccatiHWSV1F(systemDatas,ARM_ODEFunc::RK5Adaptative,isStdFormula), itsSystemDatasSO(systemDatasSO) {}
	virtual ~ARM_RiccatiHWSV1F_SO() {}

	virtual void derivs(double t, ARM_GP_Vector* yt, ARM_GP_Vector* dyt) const;
};

//-----------------------------------------------------------------------------
// \class ARM_HWSV1F
// \brief
//  1 factor Hull & White Stochastic  pricing model for closed form,
//  backward and forward diffusion abilities
//-----------------------------------------------------------------------------

class ARM_HWSV1F : public ARM_HWSV
{
protected:

	/// Nested structure to solve NR equation in SLN moment matching
	struct SLNMatch;
	struct SLNMatchD1
	{
		private:
			SLNMatch *itsF;
		public:
			SLNMatchD1(SLNMatch *F) : itsF(F) {};
			double operator () ( double x ) const
			{
				/// Exact derivative may be a little bit time consumming...
				double h=fabs(x)*0.001;
				if(h<1.0e-6)
					h=1.0e-6;
				return ((*itsF)(x+h)-(*itsF)(x-h))/(2.0*h);
			}
	};
	struct SLNMatch
	{
		private:
			double itsESt;
			double itsESt2;
			double itsESt3;
			double itsEST;
			double itsEST2;
			double itsEST3;
			mutable double itsVar;
			SLNMatchD1 *itsDerivative;
		public:
			SLNMatch(double ESt,double ESt2,double ESt3,double EST,double EST2,double EST3)
				: itsESt(ESt),itsESt2(ESt2),itsESt3(ESt3),itsEST(EST),itsEST2(EST2),itsEST3(EST3)
			{itsDerivative = new SLNMatchD1(this);}
			~SLNMatch() { delete itsDerivative;}
			double GetVar() const { return itsVar; }
			double operator () ( double x ) const
			{
				double x2 = x*x;
				double x3 = x2*x;
				double EStx = itsESt+x, ESTx = itsEST+x;
				double EStx2 = itsESt2+2*itsESt*x+x2,ESTx2 = itsEST2+2*itsEST*x+x2;
				double EStx3 = itsESt3+3*itsESt2*x+3*itsESt*x2+x3,ESTx3 = itsEST3+3*itsEST2*x+3*itsEST*x2+x3;
				double A = log(ESTx/EStx);
				double B = log(ESTx2/EStx2);
				double C = log(ESTx3/EStx3);
				itsVar = B-2*A;
				return (C - 3*(B-A))/pow(itsVar,1.5);

			}
			SLNMatchD1* Derivative() const {return itsDerivative;}
	};

	CC_IS_MUTABLE vector< ARM_RiccatiHWSV1F::ARM_SystemData >		itsSystemDatas;
	CC_IS_MUTABLE vector< ARM_RiccatiHWSV1F_SO::ARM_SystemData >	itsSystemDatasSO;
	
	/// Function parameters for all sampling intervals
    CC_IS_MUTABLE vector< ARM_RiccatiHWSV1F::ARM_FunctionData > itsFunctionDatas;

	/// To maintain continuity of the angular part of Gamma(u=(1.0,t)) & Gamma(u=(0,t))
	CC_IS_MUTABLE ARM_GP_Vector itsAngularShifts1;
	CC_IS_MUTABLE ARM_GP_Vector itsPrevValues1;
	CC_IS_MUTABLE ARM_GP_Vector itsAngularShifts2;
	CC_IS_MUTABLE ARM_GP_Vector itsPrevValues2;

	CC_IS_MUTABLE ARM_GP_Matrix itsRealizedVar;
	CC_IS_MUTABLE size_t		itsVarIdx;

	CC_IS_MUTABLE vector<ARM_PricingStatesPtr> itsRealizedStates;

	/// Collect & save model parameters for function definition
	void ComputeRiccatiSchedule(double T, ARM_GP_Vector& schedule, double Tref, bool isSOFormula=false, double Tstart=0.0) const;

	/// Analytical Riccati version with sampled mrs exponentials
	void InitAnalyticalData(double evalTime, double F0, double Te, double expT0,
							bool isSOFormula,
							ARM_GP_Vector& schedule=ARM_GP_Vector(0), double mu0=0.0) const;

	ARM_VectorPtr ComputeAnalyticalOptionPrice(	double evalTime,
												const ARM_VectorPtr& F0,
												double T0,
												double Te,
												const ARM_GP_Vector& newStrikes,
												int	RecPay,
												double payNotional,
												const ARM_PricingStatesPtr& states) const;

	/// Analytical Riccati version with kappa linked to rho, nu and vol
    virtual void InitFunctionData(	double F0,
									double Te,
									double expTe) const;

	void UpdateUDependentData(		double F0,
									double Ux,
									double Uy,
									double expT0) const;

	virtual void GenerateFunction (		double expTe,
										double& PsiX,
										double& PsiY,
										ARM_GP_Vector& prevValues,
										ARM_GP_Vector& angularShifts) const;

	ARM_VectorPtr ComputeOptionPrice(	double F0,
										double T0,
										double Te,
										double newStrike,
										int	RecPay,
										double payNotional,
										const ARM_PricingStatesPtr& states) const;

	/// Numerical version and Riccati ODE system solved by Runge-Kutta
	void InitSystemData(double evalTime,double F0, double mrs, double Te, double expT0, size_t maxSize,
						bool isStdFormula, bool isSOFormula=false, double Mu0=0.0) const;
	void UpdateUDependentSystemData(const ARM_GP_Vector& ImU, bool isStdFormula, bool isSOFormula=false) const;

	ARM_VectorPtr ComputeRungeKuttaOptionPrice(	double evalTime,
												const ARM_VectorPtr& F0,
												double T0,
												double Te,
												const ARM_GP_Vector& newStrikes,
												int	RecPay,
												double payNotional,
												const ARM_PricingStatesPtr& states) const;

	ARM_VectorPtr ComputeSpreadOptionPrice(	double evalTime,
											const ARM_GP_Vector& F0,
											const ARM_GP_Vector& Mu0,
											double mrs,
											double Tp, double Te,
											const ARM_GP_Vector& newStrikes,
											int	callPut,
											double payNotional,
											double payPeriod,
											const ARM_PricingStatesPtr& states,
											const ARM_IntVector& statusITM) const;

	ARM_VectorPtr ComputeCapletF0(	double evalTime,
									double fwdStartTime, 
									double fwdEndTime,
									double fwdPeriod,
									const ARM_GP_Vector& strikes,
									const ARM_PricingStatesPtr& states,
									ARM_GP_Vector& newStrikes) const;

	ARM_VectorPtr ComputeSwaptionF0(	double evalTime,
										double startTime, 
										double endTime,
										const ARM_GP_Vector& fixPayTimes,
										const ARM_GP_Vector& fixPayPeriods,
										int callPut,
										const ARM_GP_Matrix& strikes,
										bool isConstantNotional,
										const ARM_GP_Vector& fixNotional,
										const ARM_GP_Vector& floatNotional,
										size_t refNotionalIdx,
										const ARM_PricingStatesPtr& states,
										ARM_GP_Vector& newStrikes) const;

	ARM_VectorPtr ComputeSwapRateF0(	double evalTime,
										double startTime, 
										double payTime, 
										const ARM_GP_Vector& fixPayTimes,
										const ARM_GP_Vector& fixPayPeriods,
										double mrs,
										const ARM_PricingStatesPtr& states,
										ARM_VectorPtr& swapRate) const;

	void UpdateStdHWModel(double volFactor=1.0) const;

public:

	ARM_HWSV1F(){}

	//ARM_HWSV1F(const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params=NULL,
	//	int solverType=ARM_ODEFunc::RK5Adaptative,const ARM_GP_Vector& solverParams=ARM_GP_Vector(0),
	//	int forumalType=HWSVNumericals::Lewis,const ARM_GP_Vector& formulaParams=ARM_GP_Vector(0),
	//	double maxDecay=0.0);
	ARM_HWSV1F(const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params=NULL,
		int solverType=ARM_ODEFunc::RK5Adaptative,const ARM_GP_Vector& solverParams=ARM_GP_Vector(0),
		int forumalType=HWSVNumericals::Lewis,const ARM_GP_Vector& formulaParams=ARM_GP_Vector(0),
		int formulaTypeSO=HWSVNumericals::Heston,const ARM_GP_Vector& formulaParamsSO=ARM_GP_Vector(0),
		double maxDecay=0.0,double maxDecaySO=0.0);

	ARM_HWSV1F(const ARM_HWSV1F& rhs);
    ARM_HWSV1F& operator = (const ARM_HWSV1F& rhs);
	virtual ~ARM_HWSV1F() {}

	/// DF
	virtual ARM_VectorPtr DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const;

	/// only for variable notional swaptions (numerical integration)
	/// if std swaption, this method will call the ARM_HullWhite method
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


	virtual ARM_VectorPtr VanillaCaplet(
		const string& curveName, 
		double evalTime,
		double payTime, 
		double period,
        double payNotional,
		double fwdResetTime, 
		double fwdStartTime,
        double fwdEndTime,
		double fwdPeriod,
		const ARM_GP_Vector& strikesPerState,
        int capFloor,
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr MaxRate(		
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
        double firstReset,
		double firstStart,
		const ARM_VectorPtr& firstRate,
		int MaxOrMin,
		int ResetFreq,
		const ARM_VectorPtr& strikes,
		int CapOrFloor,
		double RhoMinMax,
		bool IsAccrued,
		double MinAccrued,
		double MaxAccrued,
        const ARM_PricingStatesPtr& states) const;

ARM_VectorPtr  VanillaSpreadOptionLet(
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


	double GetMrs() const;

    // Give local drifts and variances w.r.t. a given schedule
	virtual void NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual size_t ModelStatesSize() const { return 3;} // X, V and Phi

	ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;

	ARM_VectorPtr LocalDiscounts(size_t timeIdx, double dt, const ARM_PricingStatesPtr& states) const;

     /// Calibration purpose
	virtual bool ValidateModelParams(const ARM_ModelParams& params) const;
    virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb=0 );

	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;


	//////////////////////////////////////////////////////////////////////////////////////////////
    /// Standard ARM object support
	//////////////////////////////////////////////////////////////////////////////////////////////
	virtual ARM_Object* Clone() const;
	virtual string ExportShortName() const { return "LHWS1";}
	virtual string toString(const string& indent="",const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
