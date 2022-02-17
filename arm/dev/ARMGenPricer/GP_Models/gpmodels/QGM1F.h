/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file QGM1F.h
 *
 *  \brief 
 *
 *	\author  JM Prie, Amine Triki
 *	\version 1.0
 *	\date July 2004
 */


#ifndef _INGPMODELS_QGM1F_H
#define _INGPMODELS_QGM1F_H


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/env.h"
#include "gpinfra/pricingmodelir.h"
#include "gpbase/gpmatrix.h"

#include <vector>
CC_USING_NS(std,vector)

#include <map>
CC_USING_NS(std,map)

CC_BEGIN_NAMESPACE( ARM )

class ARM_ModelParamsQGM1F;

//-----------------------------------------------------------------------------
// \class ARM_QGM1F
// \brief
//  Class for 1F Quadratic Gaussian Model
//-----------------------------------------------------------------------------

class ARM_QGM1F : public ARM_PricingModelIR
{
private :
    void CopyNoCleanUp(const ARM_QGM1F& rhs);

    /// Internal structures for function description
    struct ARM_FunctData
    {
        double              itsTime;
        double              itsMrs;
        double              its2Vol2;
        double              itsDelta;
        int                 itsNbRoots;
        double              itsRootInf;
        double              itsRootSup;
    };

    struct ARM_FunctType
    {
        ARM_FunctType(const ARM_FunctData* data=NULL) : itsData(data), itsFlag(0) {}

        double GetBCoef(double tinf, double tsup, bool isStorable=false);
        double GetB(double tinf, double tsup, bool isStorable=false);
        double GetCACoef(double tinf, double tsup, bool isStorable=false);
        double GetVar(double tinf, double tsup, bool isStorable=false);

        virtual void InitCst(double t, double a) = 0;

        virtual double ComputeA(double t) const = 0;
        virtual double ComputeBCoef(double tinf, double tsup) const = 0;
        virtual double ComputeB(double tinf, double tsup) const = 0;
        virtual double ComputeCACoef(double tinf, double tsup) const = 0;
        virtual double ComputeVar(double tinf, double tsup) const = 0;

        const ARM_FunctData*  itsData;

        double          itsCst;
        short           itsFlag;
        double          itsBCoef;
        double          itsB;
        double          itsCACoef;
        double          itsVar;
    };

    struct ARM_Funct0 : public ARM_FunctType
    {
        /// Discriminant < 0 : no solution
        ARM_Funct0(const ARM_FunctData* data=NULL) : ARM_FunctType(data) {}

        virtual void InitCst(double t, double a);

        virtual double ComputeA(double t) const;
        virtual double ComputeBCoef(double tinf, double tsup) const;
        virtual double ComputeB(double tinf, double tsup) const;
        virtual double ComputeCACoef(double tinf, double tsup) const;
        virtual double ComputeVar(double tinf, double tsup) const;
    };

    struct ARM_Funct1 : public ARM_FunctType
    {
        /// Discriminant = 0 : one solution
        ARM_Funct1(const ARM_FunctData* data=NULL) : ARM_FunctType(data) {}

        virtual void InitCst(double t, double a);

        virtual double ComputeA(double t) const;
        virtual double ComputeBCoef(double tinf, double tsup) const;
        virtual double ComputeB(double tinf, double tsup) const;
        virtual double ComputeCACoef(double tinf, double tsup) const;
        virtual double ComputeVar(double tinf, double tsup) const;
    };

    struct ARM_Funct2 : public ARM_FunctType
    {
        /// Discriminant > 0 : two solutions
        ARM_Funct2(const ARM_FunctData* data=NULL) : ARM_FunctType(data) {}
    };

    struct ARM_Funct2Out : public ARM_Funct2
    {
        /// Solution outside roots
        ARM_Funct2Out(const ARM_FunctData* data=NULL) : ARM_Funct2(data) {}

        virtual void InitCst(double t, double a);

        virtual double ComputeA(double t) const;
        virtual double ComputeBCoef(double tinf, double tsup) const;
        virtual double ComputeB(double tinf, double tsup) const;
        virtual double ComputeCACoef(double tinf, double tsup) const;
        virtual double ComputeVar(double tinf, double tsup) const;
    };

    struct ARM_Funct2In : public ARM_Funct2
    {
        /// Solution inside roots
        ARM_Funct2In(const ARM_FunctData* data=NULL) : ARM_Funct2(data) {}

        virtual void InitCst(double t, double a);

        virtual double ComputeA(double t) const;
        virtual double ComputeBCoef(double tinf, double tsup) const;
        virtual double ComputeB(double tinf, double tsup) const;
        virtual double ComputeCACoef(double tinf, double tsup) const;
        virtual double ComputeVar(double tinf, double tsup) const;
    };

    struct ARM_Funct2Sup : public ARM_Funct2
    {
        /// Solution equal to the greater root
        ARM_Funct2Sup(const ARM_FunctData* data=NULL) : ARM_Funct2(data) {}

        virtual void InitCst(double t, double a);

        virtual double ComputeA(double t) const;
        virtual double ComputeBCoef(double tinf, double tsup) const;
        virtual double ComputeB(double tinf, double tsup) const;
        virtual double ComputeCACoef(double tinf, double tsup) const;
        virtual double ComputeVar(double tinf, double tsup) const;
    };

    struct ARM_Funct2Inf : public ARM_Funct2
    {
        /// Solution equal to the lower root
        ARM_Funct2Inf(const ARM_FunctData* data=NULL) : ARM_Funct2(data) {}

        virtual void InitCst(double t, double a);

        virtual double ComputeA(double t) const;
        virtual double ComputeBCoef(double tinf, double tsup) const;
        virtual double ComputeB(double tinf, double tsup) const;
        virtual double ComputeCACoef(double tinf, double tsup) const;
        virtual double ComputeVar(double tinf, double tsup) const;
    };

    typedef ARM_CountedPtr< ARM_FunctType >             ARM_FunctTypePtr;
    typedef map< double,vector< ARM_FunctTypePtr > >    ARM_FunctionsMap;
    typedef ARM_FunctionsMap::iterator                  ARM_FunctionsMapIter;

    /// Function parameters for all sampling intervals
    vector< ARM_FunctData > itsFunctionDatas;

    /// For each maturity, function description over all sampling intervals
    /// Mutable because GenerateFunction() will update this map
    CC_IS_MUTABLE ARM_FunctionsMap itsFunctions;

    /// Precomputed surface of the reconstruction formula
    ARM_GP_Vector      itsTime;        // t
    ARM_GP_Vector      itsMaturity;    // T
    ARM_GP_Matrix       itsA;           // a(t,T)
    ARM_GP_Matrix       itsB;           // b(t,T)
    ARM_GP_Matrix       itsC;           // c(t,T)

    /// Collect & save model parameters for function definition
    void InitFunctionDatas();

    /// Build the function solution of the Riccati equation
    ARM_FunctionsMapIter GenerateFunction(double maturity) const;

    double IntegrateB2(double tinf, double tsup,
        double time,vector< ARM_FunctTypePtr >& timeFunction,
        double maturity,vector< ARM_FunctTypePtr >& maturityFunction) const;

    double IntegrateBBCoef(double tinf, double tsup,int tIdx,
        double numeraireTime,vector< ARM_FunctTypePtr >& numeraireTimeFunction) const;

    /// Computation of A,B & C for the Zc formula :
    /// Zc(t,T) = Zc(0,T)/Zc(0,t) * exp{-A(t,T).X(t)^2 - B(t,T).X(t) - C(t,T)}
    double ComputeA(double time,const vector< ARM_FunctTypePtr >& function) const;
    double ComputeB(double time,double maturity,vector< ARM_FunctTypePtr >& function) const;
    double ComputeC(double time,vector< ARM_FunctTypePtr >& timeFunction,
                    double maturity,vector< ARM_FunctTypePtr >& maturityFunction) const;

    /// Computation of mean and variance of X(maturity) knowing states at time
    ARM_VectorPtr ConditionalMean(double time,int timeIdx,double maturity,int maturityIdx,
        double numeraireTime,vector< ARM_FunctTypePtr >& numeraireTimeFunction,
	    const ARM_PricingStatesPtr& states,double absoluteDrift) const;

    double Variance(double time,int timeIdx,double maturity,int maturityIdx,
        double numeraireTime,vector< ARM_FunctTypePtr >& numeraireTimeFunction) const;

    double SwaptionRecPayoff(double x,
                             const ARM_GP_Vector& AtT0Tp,
                             const ARM_GP_Vector& BtT0Tp,
                             const ARM_GP_Vector& CtT0Tp,
                             const ARM_GP_Vector& zcCoef,
                             const ARM_GP_Vector& zc0T0Tp,
                             bool isSynchroEnd) const;

    ARM_GP_VectorPtr LocateExerciseSignedZeros(double firstRoot,double endRoot,
                                                 const ARM_GP_Vector& ATTpT0,
                                                 const ARM_GP_Vector& BTTpT0,
                                                 const ARM_GP_Vector& CTTpT0,
                                                 const ARM_GP_Vector& zc0TpT0,
                                                 const ARM_GP_Vector& zcCoef,
                                                 bool isSynchroEnd, int payRec,
                                                 double Amax,double Bmax,
                                                 double minStep,double maxMove) const;

	ARM_GP_VectorPtr LocateExerciseZeros(double rminmin,double rmaxmax,
                                          const ARM_GP_Vector& AtT0Tp,
                                          const ARM_GP_Vector& BtT0Tp,
                                          const ARM_GP_Vector& CtT0Tp,
                                          const ARM_GP_Vector& zc0T0Tp,
                                          const ARM_GP_Vector& zcCoef,
                                          bool isSynchroEnd, int payRec,
                                          double Amax,double Bmax,
                                          bool& isFirstExer) const;

public:
	ARM_QGM1F(const ARM_ZeroCurvePtr& zc, const ARM_ModelParamsQGM1F& params);
	ARM_QGM1F(const ARM_QGM1F& rhs);
	virtual ~ARM_QGM1F();

    ARM_QGM1F& operator = (const ARM_QGM1F& rhs);

	void InitABC(const ARM_GP_Vector& times, const vector< ARM_GP_Vector >& maturities);
	
    /// Computation of A(t,T), B(t,T) and C(t,T)
    double A(double t, double T) const;
    double B(double t, double T) const;
    double C(double t, double T) const;

    ARM_VectorPtr XLaw(double time, double maturity, double numeraireTime,
                       const ARM_PricingStatesPtr& states,double absoluteDrift=0) const;

    /// Implied vol fcts for testing purpose
    double CapletImpliedVol(
		    double payTime, 
		    double period,
            double payNotional,
		    double fwdResetTime, 
		    double fwdStartTime,
            double fwdEndTime,
		    double fwdPeriod,
		    double strike,
            int capFloor,
		    double price) const;

    double SwaptionImpliedVol(
		    double swapResetTime,
            double swapNotional,
		    double floatStartTime,
		    double floatEndTime,
		    const ARM_GP_Vector& fixPayTimes,
		    const ARM_GP_Vector& fixPayPeriods,
		    double strike,
            int callPut,
		    double price) const;

    /// Reconstruction formula
	virtual ARM_VectorPtr DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr Libor( 
		const string& curveName,
        double evalTime,
		double fwdStartTime,
        double fwdEndTime,
		double period,
        double resetTime,
        double payTime,
        const ARM_PricingStatesPtr& states) const;

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
		bool isConstantNotional = true ,
		bool isConstantSpread = true ,
		bool isConstantStrike = true ) const;
    
    virtual ARM_VectorPtr VanillaDigital(
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

	/// Variable Notional swaptions can be priced in QGM1F
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
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
		const ARM_PricingStatesPtr& states,
		bool isConstantNotional = true ,
		bool isConstantSpread = true ,
		bool isConstantStrike = true ) const;

public:
    /// Calibration purpose
    virtual void Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter){};
    virtual void PreProcessing(ARM_ModelFitter& modelFitter);
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter);
    virtual void AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter);
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter);
    virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb=0 );
	virtual void ValidateCalibMethod(ARM_CalibMethod& calibMethod);
    virtual bool ValidateModelParams(const ARM_ModelParams& params) const;

    /// Default initialisation of the model
	virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);
    virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;

    virtual ARM_GP_Vector* ComputeModelTimes( const ARM_TimeInfoPtrVector& timeInfos );
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const {return static_cast<ARM_VectorPtr>(NULL);}
    virtual bool SupportBackwardInduction() const {	return true;}
    virtual bool SupportForwardInduction()  const {	return true;}
	virtual bool SupportAnalyticMarginal()  const {	return true;}

    // Give local drifts and variances w.r.t. a given schedule
	virtual double ComputeLocalDrift(double t, double u, int lastIndex, const vector< ARM_FunctTypePtr >& function) const;
	virtual double ComputeLocalVariance(double t, double u, int lastIndex, const vector< ARM_FunctTypePtr >& function) const;
	virtual double ComputeDeterministicDrift(double evalTime, double maturity,double T, vector< ARM_FunctTypePtr >& function) const;

    virtual void IntegratedLocalDrifts(
		const ARM_GP_Vector& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	virtual void NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void ModelStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

    virtual void NumMethodStateGlobalVariances(
        const ARM_GP_Vector& timeSteps,
        ARM_MatrixVector& variances) const;

    virtual bool NeedsToCholeskyDecomposeFactors( ) const {return false;}

	virtual ARM_BoolVector NeedMCIntegProcess() const { return ARM_BoolVector(1, true); };

    virtual double VarianceToTime(double var,double minTime=0.0,double maxTime=5*K_YEAR_LEN) const {return 0;}

	virtual ARM_VectorPtr  VanillaSpreadOptionLet(const string& curveName,
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

	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;

	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const;

    /// Standard ARM object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LQM1F";}
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
