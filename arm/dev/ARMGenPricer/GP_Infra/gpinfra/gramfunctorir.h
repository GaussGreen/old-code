/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gramfunctorir.h
 *
 *  \brief gram functor are functor for the grammar of the generic pricer
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */



#ifndef _INGPINFRA_GRAMFUNCTORIR_H
#define _INGPINFRA_GRAMFUNCTORIR_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "gramfunctorbase.h"

#include "gpinfra/irrate.h"


CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_PricingFunctionDF;
class ARM_PricingFunctionIR;
class ARM_MultiAssetsModel;

/// Discount Factor
struct ARM_GP_DF : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_DF( *this ) ); }

private:
	void GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	string                  itsCurveName;
	double                  itsEvalTime;
	double                  itsMaturityTime;
	ARM_PricingModel*       itsModel;

	/// name for error writting
    static string itsFuncName;
};

/// Libor Rate
struct ARM_GP_Libor : public ARM_GramFctor
{
	/// default constructor to allow its creation!
	ARM_GP_Libor() : ARM_GramFctor() {}

	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_Libor( *this ) ); }

private:
	void GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	string                  itsCurveName;
	double                  itsEvalTime;
	double                  itsFwdStartTime;
	double                  itsFwdEndTime;
	double                  itsResetTime;
	double                  itsPayTime;
	double                  itsPeriod;
	ARM_PricingFunctionIR*  itsModelIR;

	/// name for error writting
    static string itsFuncName;
};

/// Annuity key word
struct ARM_GP_Annuity : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_Annuity( *this ) ); }

private:
	void GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	string                  itsCurveName;
	double                  itsEvalTime;
	ARM_GP_Vector           itsFixPayTimes;
	ARM_GP_Vector           itsFixPayPeriods;
	ARM_PricingModelIR*		itsModelIR;
	ARM_GP_VectorPtr        itsVariableNotional;
	bool					itsIsVariableNotional;
	/// name for error writting
    static string itsFuncName;
};

/// Swap Rate
struct ARM_GP_SwapRate : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_SwapRate( *this ) ); }

private:
	/// arg is given as ARM_GramFctorArgVector& and not const ARM_GramFctorArgVector& because it modifies the arg!
	void SetDefaults( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// name for error writting
    static string itsFuncName;

protected:
	/// arg is given as ARM_GramFctorArgVector& and not const ARM_GramFctorArgVector& because it modifies the arg!
	void GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	
    /// members to store
	string                  itsCurveName;
	double                  itsEvalTime;
	double                  itsFloatStartTime;
	double                  itsFloatEndTime;
	ARM_GP_Vector           itsFwdResetTimes;    /// For the interpolation of the margin
	ARM_GP_Vector           itsFixResetTimes;	/// For the interpolation of strikes
	ARM_GP_Vector           itsFixPayTimes;
	ARM_GP_Vector           itsFixPayPeriods;
	ARM_GP_Vector           itsFwdStartTimes;
	ARM_GP_Vector           itsFwdEndTimes;
	ARM_GP_Vector           itsFwdPayPeriods;
	ARM_GP_Vector           itsFloatPayTimes;
	ARM_GP_Vector           itsFloatPayPeriods;
    ARM_GP_VectorPtr        itsMarginVector;
    bool                    itsDbleNotional;
	ARM_PricingFunctionIR*  itsModelIR;
	size_t					itsFloatDateStripOffset;
	size_t					itsFixDateStripOffset;
	
};

/// Spread
struct ARM_GP_Spread : public ARM_GramFctor
{
	/// default constructor to allow its creation!
	ARM_GP_Spread() : ARM_GramFctor() {}

	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_Spread( *this ) ); }

private:
	void GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	string                  itsCurveName;
	double                  itsEvalTime;

	double					itsCoeff1;
	double                  itsFloatStartTime1;
	double                  itsFloatEndTime1;
	ARM_GP_Vector           itsFwdResetTimes1;		/// For the interpolation of the margin 1
	ARM_GP_Vector           itsFixResetTimes1;		/// For the interpolation of strikes 1
	ARM_GP_Vector           itsFixPayTimes1;
	ARM_GP_Vector           itsFixPayPeriods1;
	ARM_GP_Vector           itsFwdStartTimes1;
	ARM_GP_Vector           itsFwdEndTimes1;
	ARM_GP_Vector           itsFwdPayPeriods1;
	ARM_GP_Vector           itsFloatPayTimes1;
	ARM_GP_Vector           itsFloatPayPeriods1;
	ARM_GP_VectorPtr        itsMarginVector1;

	double					itsCoeff2;
	double                  itsFloatStartTime2;
	double                  itsFloatEndTime2;
	ARM_GP_Vector           itsFwdResetTimes2;		/// For the interpolation of the margin 2
	ARM_GP_Vector           itsFixResetTimes2;		/// For the interpolation of strikes 2
	ARM_GP_Vector           itsFixPayTimes2;
	ARM_GP_Vector           itsFixPayPeriods2;
	ARM_GP_Vector           itsFwdStartTimes2;
	ARM_GP_Vector           itsFwdEndTimes2;
	ARM_GP_Vector           itsFwdPayPeriods2;
	ARM_GP_Vector           itsFloatPayTimes2;
	ARM_GP_Vector           itsFloatPayPeriods2;
    ARM_GP_VectorPtr        itsMarginVector2;

    ARM_PricingFunctionIR*  itsModelIR;

	/// name for error writting
    static string itsFuncName;
};

/// Swap operator
struct ARM_GP_Swap : public ARM_GP_SwapRate
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_Swap( *this ) ); }

private:
	/// arg is given as ARM_GramFctorArgVector& and not const ARM_GramFctorArgVector& because it modifies the arg!
	void GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	/// arg is given as ARM_GramFctorArgVector& and not const ARM_GramFctorArgVector& because it modifies the arg!
	void SetDefaults( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	ARM_GP_VectorPtr    itsFixNotionalVector;
	ARM_GP_VectorPtr    itsFloatNotionalVector;
	double              itsStrikeDouble;
	ARM_GP_MatrixPtr    itsStrikeMatrix;
	int                 itsPayRec; 
    
	/// name for error writting
    static string itsFuncName;

};

/// Swap operator
struct ARM_GP_BasisSwap : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_BasisSwap( *this ) ); }

private:

	/// arg is given as ARM_GramFctorArgVector& and not const ARM_GramFctorArgVector& because it modifies the arg!
	void SetDefaults( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );
	/// arg is given as ARM_GramFctorArgVector& and not const ARM_GramFctorArgVector& because it modifies the arg!
	void GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	/// members to store

	string                  itsDomCurveName;
	string                  itsForCurveName;
	string                  itsFxCurveName;
	double                  itsEvalTime;
	double                  itsStartTime;
	double                  itsEndTime;	

	ARM_GP_Vector           itsDomResetTimes;	    
	ARM_GP_Vector           itsDomFwdStartTimes;
	ARM_GP_Vector           itsDomFwdEndTimes;
	ARM_GP_Vector           itsDomFlowStartTimes;			
	ARM_GP_Vector           itsDomFlowEndTimes;	
	ARM_GP_Vector           itsDomFwdPeriods;	
	ARM_GP_Vector           itsDomPayTimes;
	ARM_GP_Vector           itsDomPayPeriods;
	ARM_GP_VectorPtr        itsDomMarginVector;
	ARM_GP_VectorPtr        itsDomNotionalVector;
	bool                    itsIsDomFlottant;

	ARM_GP_Vector           itsForResetTimes;       
	ARM_GP_Vector           itsForFwdStartTimes;
	ARM_GP_Vector           itsForFwdEndTimes;
	ARM_GP_Vector           itsForFlowStartTimes;   		
	ARM_GP_Vector           itsForFlowEndTimes;	    
	ARM_GP_Vector           itsForFwdPeriods;	
	ARM_GP_Vector           itsForPayTimes;
	ARM_GP_Vector           itsForPayPeriods;
	ARM_GP_VectorPtr        itsForMarginVector;
	ARM_GP_VectorPtr        itsForNotionalVector;
	bool                    itsIsForFlottant;
	string                  itsExNotionalType;
	ARM_GP_Vector           itsFxResetTimes;
	ARM_GP_Vector           itsFxSettlTimes;  
	ARM_GP_Vector           itsFxPayTimes;

	ARM_GP_MatrixPtr		itsDomStrikeMatrix;
	ARM_GP_MatrixPtr		itsForStrikeMatrix;
	int						itsPayRec; 
	ARM_MultiAssetsModel*   itsModelIR;
	
/// name for error writting
    static string itsFuncName;

/// number of arguments
	static long itsNbArg;

};

/// Euro Cap Floor /Digital Cap Floor-let operator
struct ARM_GP_CapDigital_Common :  public ARM_GramFctor
{
	typedef ARM_VectorPtr (ARM_PricingFunctionIR::*ARM_CapNDigital_PricingFunc)(
		const string& curveName,
		double evalTime,
		double payTime,			
		double period,
        double payNotional,
		double fwdResetTime,
		double fwdStartTime,
        double fwdEndTime,
        double fwdPeriod,
        const ARM_GP_Vector& strike,
        int capFloor,
        const ARM_PricingStatesPtr& states) const;    

	/// constructor ( no copy constructor, destructor and assignment as this is a very light object)
	/// that should not be copied over! One created once!
	ARM_GP_CapDigital_Common( const ARM_CapNDigital_PricingFunc& pricingFunc, const string& name );
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const{	return ARM_GramFctorPtr( new ARM_GP_CapDigital_Common( *this ) ); 	}

private:
	void GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	string                      itsCurveName;
	double                      itsEvalTime;
	double                      itsFwdStartTime;
	double                      itsFwdEndTime;
	double                      itsFwdPeriod;
	double                      itsResetTime;
	double                      itsPayTime;
	double                      itsNotional;
	double                      itsPeriod;
	double                      itsStrikeDouble;
	ARM_VectorPtr               itsStrikeVector;
	int                         itsCapFloor;
	ARM_PricingFunctionIR*      itsModelIR;

	ARM_CapNDigital_PricingFunc itsPricingFunc;

	/// name for error writting
    string itsFuncName;
};


/// Swaption operator
struct ARM_GP_Swaption :  public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_Swaption( *this ) ); }
	
private:
	void GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	string                  itsCurveName;
	double                  itsEvalTime;
	double                  itsSwapResetTime;
	double                  itsFloatStartTime;
	double                  itsFloatEndTime;
	bool					itsIsConstantNotional;
	ARM_GP_Vector			itsFloatResetTimes;
	ARM_GP_Vector			itsFloatStartTimes;
	ARM_GP_Vector			itsFloatPayTimes;
	ARM_GP_Vector			itsFloatIntTerms;
	ARM_GP_VectorPtr		itsFixNotionalVector;
	ARM_GP_VectorPtr		itsFloatNotionalVector;
	ARM_GP_MatrixPtr        itsStrikeMatrix;
	int                     itsCallPut; // on the swap rate
	ARM_GP_Vector           itsFixPayTimes;
	ARM_GP_Vector           itsFixPayPeriods;
	ARM_PricingFunctionIR*  itsModelIR;

	/// name for error writting
    static string itsFuncName;
};

/// Swaption operator
struct ARM_GP_ImpliedVol :  public ARM_GramFctor
{
	/// constructor ( no copy constructor, destructor and assignment as this is a very light object)
	/// that should not be copied over! One created once!
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const{	return ARM_GramFctorPtr( new ARM_GP_ImpliedVol( *this ) ); 	}

private:
	void GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	string                      itsCurveName;
	double                      itsEvalTime;
	double                      itsFwdStartTime;
	double                      itsFwdEndTime;
	double                      itsFwdPeriod;
	double                      itsResetTime;
	double                      itsPayTime;
	double                      itsNotional;
	double                      itsPeriod;
	double                      itsStrikeDouble;
	ARM_VectorPtr               itsStrikeVector;
	int                         itsCapFloor;
	ARM_PricingFunctionIR*      itsModelIR;

	/// name for error writting
    static string itsFuncName;
};


struct ARM_GP_Cap : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_Cap( *this ) ); }

private:
	void GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	string                  itsCurveName;
	double                  itsEvalTime;
	ARM_GP_Vector           itsPayTimes;
	ARM_GP_Vector           itsPeriods;
	double                  itsNotional;
	ARM_GP_Vector           itsFwdResetTimes;
	ARM_GP_Vector           itsFwdStartTimes;
	ARM_GP_Vector           itsFwdEndTimes;
	ARM_GP_Vector           itsFwdPeriods;
	double                  itsStrikeDouble;
	ARM_VectorPtr           itsStrikeVector;
	int                     itsCapFloor;    
	ARM_PricingFunctionIR*  itsModelIR;

	/// name for error writting
    static string itsFuncName;
};

/// Price To Yield operator
struct ARM_GP_PTYAndYTPFctor  : public ARM_GramFctor
{
	typedef double (*ARM_BondAnalytic_PricingFunc)(double, double, double, int, double, double, double );

	ARM_GP_PTYAndYTPFctor( const ARM_BondAnalytic_PricingFunc& pricingFunc, const string& funcName );
	
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );

// FIXMEFRED: mig.vc8 (24/05/2007 10:42:56):cast
	virtual ARM_GramFctorPtr Clone() const { return static_cast<ARM_GramFctorPtr>(new ARM_GP_PTYAndYTPFctor( *this )); }

private:
	ARM_PricingFunctionIR* itsModelIR;

	/// function to price the functor
	ARM_BondAnalytic_PricingFunc itsPricingFunc;

	/// name for error writting
	string itsFuncName;
};

/// SumOption operator

struct ARM_GP_SumOption : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_SumOption( *this ) ); }

private:
	void GrabInputs(ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults(ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// members to store
	string                  itsCurveName;
	double                  itsEvalTime;
	ARM_GP_Vector			itsCoeffs;
	ARM_GP_Vector           itsFwdResetTimes;
	ARM_GP_Vector           itsFwdStartTimes;
	ARM_GP_Vector           itsFwdEndTimes;
	double					itsPayTime;
	ARM_GP_Vector           itsFwdPeriods;
	double                  itsStrikeDouble;
	ARM_VectorPtr           itsStrikeVector;
	int                     itsCapFloor;    
	ARM_GP_VectorPtr		itsCoeffVector;
	ARM_PricingFunctionIR*  itsModelIR;

	/// name for error writting
    static string itsFuncName;
};


/// Spread Option Operator

struct ARM_GP_SpreadOption : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_SpreadOption( *this ) ); }
	virtual ~ARM_GP_SpreadOption ();

private:
	void GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );


	/// members to store
	string                  itsCurveName;
	double                  itsEvalTime;
	int						itsCallPut;
	/*---- to be removed ----*
	double					itsStartTime;
	double					itsEndTime;
	 *-----------------------*/
	ARM_GP_Vector           itsResetTimes;
	ARM_GP_Vector           itsPayTimes;
	ARM_GP_Vector           itsPayPeriods;
	ARM_GP_VectorPtr		itsNotionals;// could be a ARM_GP_VectorPtr
	ARM_GP_VectorPtr        itsCoeffLong;
	ARM_GP_VectorPtr        itsCoeffShort;
	ARM_GP_MatrixPtr		itsStrikes;
	ARM_GP_Vector			itsSwapLongFloatStartTime;
	ARM_GP_Vector			itsSwapLongFloatEndTime;
	ARM_VectorVector		itsSwapLongFixPayTimes;
	ARM_VectorVector		itsSwapLongFixPayPeriods;
	ARM_GP_Vector			itsSwapShortFloatStartTime;
	ARM_GP_Vector			itsSwapShortFloatEndTime;
	ARM_VectorVector		itsSwapShortFixPayTimes;
	ARM_VectorVector		itsSwapShortFixPayPeriods;
	double					itsLeveragePrev;
	ARM_PricingFunctionIR*  itsModelIR;

	/// name for error writting
    static string itsFuncName;
	
};


/// Corridor operator

struct ARM_GP_Corridor : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( 
		ARM_GramFctorArgVector& arg, 
		ARM_PricingModel* mod,
		double evalDate, 
		const ARM_PricingStatesPtr& states, 
		vector<ARM_ExpNodePtr>& nodes );

	virtual ARM_GramFctorArg EvalCorridor( 
		ARM_GramFctorArgVector& arg, 
		ARM_PricingModel* mod,
		double evalDate, 
		const ARM_PricingStatesPtr& states, 
		vector<ARM_ExpNodePtr>& nodes );

	virtual ARM_GramFctorArg EvalCMSSpreadCorridor( 
		ARM_GramFctorArgVector& arg, 
		ARM_PricingModel* mod,
		double evalDate, 
		const ARM_PricingStatesPtr& states, 
		vector<ARM_ExpNodePtr>& nodes );
	
    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_Corridor( *this ) ); }
	virtual ~ARM_GP_Corridor ();

private:
	typedef ARM_GramFctorArg (ARM_GP_Corridor::*CorridorFunc)(
		ARM_GramFctorArgVector& arg, 
		ARM_PricingModel* mod,
		double evalDate, 
		const ARM_PricingStatesPtr& states, 
		vector<ARM_ExpNodePtr>& nodes );

	void GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	void SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );


	/// members to store
	string						itsCurveName;
	double						itsEvalTime;
	int							itsRcvPay;
	int							itsPayIndexDayCount;
	ARM_GP_Vector				itsStartTimes;
	ARM_GP_Vector				itsEndTimes;
	ARM_GP_Vector				itsPayTimes;
	int							itsPayIndexType;
	ARM_GP_Vector				itsPayIndexResetTimes;
	ARM_GP_Vector				itsPayIndexStartTimes;
	//ARM_GP_Vector				itsPayIndexEndTimes;
	ARM_GP_Vector				itsPayIndexTerms;
	//ARM_GP_Vector				itsPayPeriods;
	ARM_GP_Vector				itsPayFixIndexValue;
	ARM_VectorVector			itsFixingTimes;
	int							itsFirstIndexType;
	ARM_VectorVector			itsFirstIndexStartTimes;
	ARM_VectorVector			itsFirstIndexEndTimes;
	ARM_VectorVector			itsFirstIndexTerms;
	ARM_SwapRatePtrVectorVector	itsFirstIndexSwapRates;
	int							itsSecondIndexType;
	ARM_VectorVector			itsSecondIndexStartTimes;
	ARM_VectorVector			itsSecondIndexEndTimes;
	ARM_VectorVector			itsSecondIndexTerms;
	ARM_SwapRatePtrVectorVector	itsSecondIndexSwapRates;
	ARM_VectorVector			itsFixingWeights;
	ARM_GP_Vector				itsFixValues;
	ARM_GP_Vector				itsPayIndexMults;
	ARM_GP_Vector				itsSpreads;
	ARM_GP_Vector				itsNotionals;
	ARM_VectorVector			itsCoeffs1;
	ARM_VectorVector			itsCoeffs2;
	ARM_VectorVector			itsBarriersDown;
	ARM_VectorVector			itsBarriersUp;

	int							itsThirdIndexType;
	ARM_SwapRatePtrVectorVector	itsThirdIndexSwapRates;
	ARM_VectorVector			itsBarriersDown3;
	ARM_VectorVector			itsBarriersUp3;


	ARM_PricingFunctionIR*		itsModelIR;

	/// case of a CMS spread corridor with variable payments (CMS or LIBOR)
	ARM_SwapRatePtrVector		itsPayRates; /// vector of payed rates (CMS or LIBOR)

	CorridorFunc			itsCorridorFunc;

	/// name for error writting
    static string itsFuncName;
	
};

struct ARM_GP_MaxRate : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_MaxRate( *this ) ); }

private:
	/// arg is given as ARM_GramFctorArgVector& and not const ARM_GramFctorArgVector& because it modifies the arg!
	void SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// name for error writting
    static string itsFuncName;

protected:
	/// arg is given as ARM_GramFctorArgVector& and not const ARM_GramFctorArgVector& because it modifies the arg!
	void GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	
    /// members to store
	string					itsCurveName;
	double					itsEvalTime;
	double                  itsFloatStartTime;
	double                  itsFloatEndTime;
	ARM_GP_Vector           itsFwdResetTimes;    /// For the interpolation of the margin
	ARM_GP_Vector           itsFixResetTimes;	/// For the interpolation of strikes
	ARM_GP_Vector           itsFixPayTimes;
	ARM_GP_Vector           itsFixPayPeriods;
	ARM_GP_Vector           itsFwdStartTimes;
	ARM_GP_Vector           itsFwdEndTimes;
	ARM_GP_Vector           itsFwdPayPeriods;
	ARM_GP_Vector           itsFloatPayTimes;
	ARM_GP_Vector           itsFloatPayPeriods;
    ARM_GP_VectorPtr        itsMarginVector;
    bool                    itsDbleNotional;

	double					itsFirstResetTime;
	double					itsFirstStartTime;
	ARM_VectorPtr			itsFirstRate;
	int						itsMaxOrMin;
	int						itsCapOrFloor;
	int						itsResetFreq;

	ARM_VectorPtr			itsStrikes;	/// To support deterministic or state dependent strikes

	double					itsRhoMinMax;

	bool					itsIsAccrued;
	double					itsMinAccrued;
	double					itsMaxAccrued;
	
	ARM_PricingFunctionIR*  itsModelIR;
	size_t					itsFloatDateStripOffset;
	size_t					itsFixDateStripOffset;
	
};

struct ARM_GP_DoubleDigital : public ARM_GramFctor
{
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes );
    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );
	virtual ARM_GramFctorPtr Clone() const { return ARM_GramFctorPtr( new ARM_GP_DoubleDigital( *this ) ); }

private:
	/// arg is given as ARM_GramFctorArgVector& and not const ARM_GramFctorArgVector& because it modifies the arg!
	void SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, vector<ARM_ExpNodePtr>& nodes );

	/// name for error writting
    static string itsFuncName;

	/// arg is given as ARM_GramFctorArgVector& and not const ARM_GramFctorArgVector& because it modifies the arg!
	void GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, double evalDate, vector<ARM_ExpNodePtr>& nodes );
	
    /// members to store
	string					itsModelName;
	double					itsEvalTime;

	ARM_VectorPtr			itsFirstRate;
	ARM_VectorPtr			itsFirstStrikeDown;
	ARM_VectorPtr			itsFirstStrikeUp;
	double					itsFirstStrikeSpread;

	ARM_VectorPtr			itsSecondRate;
	ARM_VectorPtr			itsSecondStrikeDown;
	ARM_VectorPtr			itsSecondStrikeUp;
	double					itsSecondStrikeSpread;

	ARM_PricingFunctionIR*  itsModelIR;
	
};

inline ARM_GP_Vector* SubVector(size_t offset,const ARM_GP_Vector& v, size_t len = -1)
{
	ARM_GP_Vector* res;
	size_t size = (len>0?len:v.size()-offset);
	if (size>0)
	{
		res=new ARM_GP_Vector(size);
		for (size_t i = 0;i < size;i++)
		{
			(*res)[i] = v[i+offset];
		}
		return res;
	}
	else
	{
		return NULL;
	}
}

CC_END_NAMESPACE()

#endif


/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

