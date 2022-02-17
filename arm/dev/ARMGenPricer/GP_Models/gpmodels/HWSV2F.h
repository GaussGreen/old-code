/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file MSV2F.h
 *
 *  \brief 
 *
 *	\author  J-M Prié
 *	\version 1.0
 *	\date September 2006
 */


#ifndef _INGPMODELS_HWSV2F_H
#define _INGPMODELS_HWSV2F_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"

#include "gpnumlib/odefunctions.h"
#include "gpmodels/enummodel.h"
#include "gpmodels/hwsv.h"

#include "gpinfra/pricingmodelir.h"
#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

//// forward declaration
class ARM_HullWhite;

class ARM_RiccatiHWSV2F : public ARM_RiccatiHWSV
{
public:
	/// Nested structure to solve time dependent Riccati system
	struct ARM_SystemData
	{
		/// cstB1 = -exp(-mrs1*T0)/mrs1, cstB2 = -exp(-mrs2*T0)/mrs2
		double itsTime;
		double itsYf;
		double itsMrs1;				/// 1st rate factor MRS
		double itsMrs2;				/// 2nd rate factor MRS
		double itsKappaTheta;		/// kappa (MRS of variance) * theta (long term variance)

		/// quadratic coefficients
		double itsModA;				/// -0.5*nu²(t)
		double itsModB;				/// -nu²(t)

		/// linear coefficients (Heston formula)
		double itsModB1;			/// kappa + nu(t)*(rho13(t)*vol1(t)/mrs1 + rho23(t)*vol2(t)/mrs2)
		double itsModB21;			/// cstB1*nu(t)*rho13(t)*vol1(t)
		double itsModB22;			/// cstB2*nu(t)*rho23(t)*vol2(t)
		double itsModB21x;			/// -nu(t)*rho13(t)*vol1(t)*F10
		double itsModB22x;			/// -nu(t)*rho23(t)*vol2(t)*F20
		ARM_GP_Matrix itsModB2y;	/// {modB21x * Im(u),
									///  modB22x * Im(u)}

		/// unit coefficients (Heston formula)
		double itsModC11;			/// 0.5*F10²*vol1(t)²
		double itsModC12;			/// F10*F20*vol1(t)*vol2(t)*rho12(t)
		double itsModC22;			/// 0.5*F20²*vol2(t)²
		ARM_GP_Matrix itsModC;		/// {{modC11 * Im(u)²,modC12 * Im(u)²,modC22 * Im(u)²},
									///  {modC11 * Im(u), modC12 * Im(u), modC22 * Im(u) }}

		double itsModD11;			/// -mu110²*vol1(t)²
		double itsModD12;			/// -(mu120+mu210)*vol1(t)*vol2(t)*rho12(t)
		double itsModD22;			/// -mu220²*vol2(t)²
		ARM_GP_Matrix itsModD;		/// { modD11 * Im(u), modD12 * Im(u),modD22 * Im(u) }
	}; 

protected :
	const vector< ARM_SystemData >& itsSystemDatas;

public:

	ARM_RiccatiHWSV2F::ARM_RiccatiHWSV2F(const vector< ARM_SystemData >& systemDatas,bool isStdFormula=false);
	ARM_RiccatiHWSV2F::ARM_RiccatiHWSV2F(const ARM_RiccatiHWSV2F& rhs);
	ARM_RiccatiHWSV2F& ARM_RiccatiHWSV2F::operator=(const ARM_RiccatiHWSV2F& rhs);
	virtual ~ARM_RiccatiHWSV2F() {}

    virtual string toString(const string& indent="", const string& nextIndent="") const {return string("ARM_RiccatiHWSV2F");}
    virtual ARM_Object* Clone() const {return new ARM_RiccatiHWSV2F(*this);}

	virtual size_t GetSystemSize() const {return itsSystemDatas.size();}
	virtual double GetYf(size_t i) const {return itsSystemDatas[i].itsYf;}
	virtual size_t GetNbPhis() const {return itsSystemDatas[0].itsModB2y.cols();}

	virtual void derivs(double t, ARM_GP_Vector* yt, ARM_GP_Vector* dyt) const;
};

struct ARM_RiccatiHWSV2F_SO : public ARM_RiccatiHWSV2F
{
	/// Nested structure to solve time dependent Riccati system
	struct ARM_SystemData
	{
		/// quadratic coefficients
		double itsModB;				/// -nu²(t)

		/// unit coefficients (Heston formula)
		double itsModD11;			/// -mu110²*vol1(t)²
		double itsModD12;			/// -(mu120+mu210)*vol1(t)*vol2(t)*rho12(t)
		double itsModD22;			/// -mu220²*vol2(t)²
		ARM_GP_Matrix itsModD;		/// { modD11 * Im(u), modD12 * Im(u),modD22 * Im(u) }
	}; 

	const vector< ARM_SystemData >& itsSystemDatasSO;

	ARM_RiccatiHWSV2F_SO::ARM_RiccatiHWSV2F_SO(const vector< ARM_RiccatiHWSV2F::ARM_SystemData >& systemDatas,
		const vector< ARM_SystemData >& systemDatasSO,bool isStdFormula=false)
		: ARM_RiccatiHWSV2F(systemDatas,isStdFormula), itsSystemDatasSO(systemDatasSO) {}
	virtual ~ARM_RiccatiHWSV2F_SO() {}

	virtual void derivs(double t, ARM_GP_Vector* yt, ARM_GP_Vector* dyt) const;
};




//-----------------------------------------------------------------------------
// \class ARM_HWSV2F
// \brief
//  1 factor Hull & White Stochastic  pricing model for closed form,
//  backward and forward diffusion abilities
//-----------------------------------------------------------------------------

class ARM_HWSV2F : public ARM_HWSV
{
private:
	CC_IS_MUTABLE vector< ARM_RiccatiHWSV2F::ARM_SystemData >		itsSystemDatas;
	CC_IS_MUTABLE vector< ARM_RiccatiHWSV2F_SO::ARM_SystemData >	itsSystemDatasSO;

	/// Collect & save model parameters for function definition
	void ComputeRiccatiSchedule(double T, ARM_GP_Vector& schedule, double Tref, bool isSOFormula=false, double Tstart=0.0) const;

	void InitAnalyticalData(double evalTime,const ARM_VectorPtr& F0, double Te, double exp1T0, double exp2T0,
							bool isSOFormula, ARM_GP_Vector& schedule) const;

	ARM_VectorPtr ComputeAnalyticalOptionPrice(	double evalTime,
												const ARM_VectorPtr& F0,
												double mrs1, double mrs2,
												double T0, double Te,
												const ARM_GP_Vector& newStrikes,
												int	RecPay,
												double payNotional,
												const ARM_PricingStatesPtr& states) const;

	void InitSystemData(double evalTime,const ARM_VectorPtr& F0, double mrs1, double mrs2,
						double Te, double exp1Tn, double exp2Tn, size_t maxSize,
						bool isStdFormula=true, bool isSOFormula=false,
						const ARM_GP_MatrixPtr& Mu0=ARM_GP_MatrixPtr(NULL)) const;

	void UpdateUDependentSystemData(const ARM_GP_Vector& ImU, bool isStdFormula=true, bool isSOFormula=false) const;

	ARM_VectorPtr ComputeRungeKuttaOptionPrice(	double evalTime,
												const ARM_VectorPtr& F0,
												double mrs1, double mrs2,
												double T0, double Te,
												const ARM_GP_Vector& newStrikes,
												int	RecPay,
												double payNotional,
												const ARM_PricingStatesPtr& states) const;

	ARM_VectorPtr ComputeRungeKuttaSpreadOptionPrice(	double evalTime,
														const ARM_GP_Vector& F0,
														const ARM_GP_Vector& Mu0,
														double mrs1, double mrs2,
														double Tp, double Te,
														const ARM_GP_Vector& newStrikes,
														int	callPut,
														double payNotional,
														double payPeriod,
														const ARM_PricingStatesPtr& states,
														const ARM_IntVector& statusITM) const;
	ARM_VectorPtr ComputeCapletF0(double evalTime,
								double fwdStartTime, 
								double fwdEndTime,
								double fwdPeriod,
								const ARM_GP_Vector& strikes,
								double mrs1,double mrs2,
								const ARM_PricingStatesPtr& states,
								ARM_GP_Vector& newStrikes) const;

	ARM_VectorPtr ComputeSwaptionF0( double evalTime,
								double startTime, 
								double endTime,
								const ARM_GP_Vector& fixPayTimes,
								const ARM_GP_Vector& fixPayPeriods,
								int callPut,
								const ARM_GP_Matrix& strikes,
								bool isConstantNotional,
								const ARM_GP_Vector& fixNotional,
								const ARM_GP_Vector& floatNotional,
								int refNotionalIdx,
								double mrs1, double mrs2,
								const ARM_PricingStatesPtr& states,
								ARM_GP_Vector& newStrikes) const;

	ARM_VectorPtr ComputeSwapRateF0(double evalTime,
								double startTime, 
								double payTime, 
								const ARM_GP_Vector& fixPayTimes,
								const ARM_GP_Vector& fixPayPeriods,
								double mrs1, double mrs2,
								const ARM_PricingStatesPtr& states,
								ARM_VectorPtr& swapRate) const;

	void UpdateStdHWModel(double volFactor=1.0) const;

	void GetMrs(double& mrs1, double& mrs2) const;

public:

	ARM_HWSV2F(const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params=NULL,
		const ARM_GP_Vector& solverParams=ARM_GP_Vector(0),
		int formulaType=HWSVNumericals::Lewis,const ARM_GP_Vector& formulaParams=ARM_GP_Vector(0),
		int formulaTypeSO=HWSVNumericals::Heston,const ARM_GP_Vector& formulaParamsSO=ARM_GP_Vector(0),
		double maxDecay=0.0,double maxDecaySO=0.0);

	ARM_HWSV2F(const ARM_HWSV2F& rhs);
    ARM_HWSV2F& operator = (const ARM_HWSV2F& rhs);
	virtual ~ARM_HWSV2F() {}

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

	virtual ARM_VectorPtr VanillaSpreadOptionLet(const string& curveName,
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

	

    // Give local drifts and variances w.r.t. a given schedule
	virtual void NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual size_t ModelStatesSize() const { return 6;} // X1, X2, V, Phi11, Phi12 and Phi22

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
	virtual string ExportShortName() const { return "LHWS2";}
	virtual string toString(const string& indent="",const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
