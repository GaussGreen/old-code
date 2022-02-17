/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file vanillaswaption.cpp
 *
 *  \brief vanilla swaption
 *	\author  E.M Ezzine E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/ostringstream.h"
#include <iomanip> /// for setprecision()
CC_USING_NS(std,setw)
CC_USING_NS(std,setprecision)

#include "gpcalib/vanillaswaption.h"

/// gpbase
#include "gpbase/gpvector.h"
#include "gpbase/curve.h"
#include "gpbase/checkarg.h"
#include "gpbase/argconvdefault.h"

/// gpinfra
#include "gpinfra/pricingmodelir.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/argconvdefault.h"
#include "gpinfra/zccrvfunctor.h"

// gpclosedform
#include "gpclosedforms/vanilla_normal.h"


CC_BEGIN_NAMESPACE( ARM )

const int NBROWS = 2;
const int NBCOLUMNS = 5;

////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArg
///	Routine: default constructor, copy constructor,
///          assigment, destructor, clone
///	Returns: 
///	Action :
////////////////////////////////////////////////////
void ARM_VanillaSwaptionArg::CopyNoCleanUp(const ARM_VanillaSwaptionArg& rhs)
{
	itsStartTime        =   rhs.itsStartTime;
    itsEndTime          =   rhs.itsEndTime;
    itsFixFrequency     =   rhs.itsFixFrequency;
	itsFloatFrequency   =   rhs.itsFloatFrequency;
	itsFixDayCount		=   rhs.itsFixDayCount;
	itsFloatDayCount	=   rhs.itsFloatDayCount;
	itsResetTime		=	rhs.itsResetTime;
	itsIsConstantNotional = rhs.itsIsConstantNotional;
	itsIsConstantSpread = rhs.itsIsConstantSpread;
	itsIsConstantStrike = rhs.itsIsConstantStrike;
	itsAsofDate			= rhs.itsAsofDate;
	itsGenSecString		= rhs.itsGenSecString;
	itsAtTheMoneyFlag	= rhs.itsAtTheMoneyFlag;

#if defined(__GP_STRICT_VALIDATION)
	ThrowErrorOnNullObject( "FixPayTimes",		rhs.itsFixPayTimes );
	ThrowErrorOnNullObject( "FixPayPeriods",	rhs.itsFixPayPeriods );
	ThrowErrorOnNullObject( "FloatResetTimes",	rhs.itsFloatResetTimes );
	ThrowErrorOnNullObject( "FloatStartTimes",	rhs.itsFloatStartTimes );
	ThrowErrorOnNullObject( "FloatEndTimes",	rhs.itsFloatEndTimes );
	ThrowErrorOnNullObject( "FloatIntTerms",	rhs.itsFloatIntTerms );
    ThrowErrorOnNullObject( "Strikes",	        rhs.itsStrikes );
	ThrowErrorOnNullObject( "FixNotional",		rhs.itsFixNominal );
    ThrowErrorOnNullObject( "FloatNotional",	rhs.itsFloatNominal );
	ThrowErrorOnNullObject( "FloatNotional",	rhs.itsFundingSpread );

#endif


	/// for fast access no test of NULL pointor in release!	
	itsFixPayTimes      =  static_cast<ARM_GP_Vector*>( rhs.itsFixPayTimes->Clone());
	itsFixPayPeriods    =  static_cast<ARM_GP_Vector*>( rhs.itsFixPayPeriods->Clone());
	itsFloatResetTimes  =  static_cast<ARM_GP_Vector*>( rhs.itsFloatResetTimes->Clone());
	itsFloatStartTimes  =  static_cast<ARM_GP_Vector*>( rhs.itsFloatStartTimes->Clone());
	itsFloatEndTimes    =  static_cast<ARM_GP_Vector*>( rhs.itsFloatEndTimes->Clone());
	itsFloatIntTerms    =  static_cast<ARM_GP_Vector*>( rhs.itsFloatIntTerms->Clone());
    itsStrikes          =  static_cast<ARM_GP_Vector*>( rhs.itsStrikes->Clone());
	itsFixNominal       =  static_cast<ARM_GP_Vector*>( rhs.itsFixNominal->Clone());
	itsFloatNominal     =  static_cast<ARM_GP_Vector*>( rhs.itsFloatNominal->Clone());
	itsFundingSpread	=  static_cast<ARM_GP_Vector*>(rhs.itsFundingSpread->Clone());



#if defined(__GP_STRICT_VALIDATION)
    //CC_NS(ARM_Check,CheckSameArgSize)(*itsStrikes, *itsFixPayTimes, "itsStrikes", "FixPayTimes",__LINE__,__FILE__);
#endif
}


ARM_VanillaSwaptionArg::ARM_VanillaSwaptionArg(const ARM_VanillaSwaptionArg& arg)
:	ARM_VanillaArg(arg),
	itsFixPayTimes( NULL ),
	itsFixPayPeriods( NULL ),  
	itsFloatResetTimes( NULL ),
	itsFloatStartTimes( NULL ),
	itsFloatEndTimes( NULL ),
	itsFloatIntTerms( NULL ),
    itsStrikes(NULL),
	itsFloatNominal(NULL),
	itsFixNominal(NULL),
	itsFundingSpread (NULL),
    itsFixFrequency(K_SEMIANNUAL),
	itsFloatFrequency(K_SEMIANNUAL),
	itsFixDayCount(K30_360),
	itsFloatDayCount(KACTUAL_360),
	itsGenSecString(NULL)
{
    CopyNoCleanUp(arg);
}

ARM_VanillaSwaptionArg& ARM_VanillaSwaptionArg::operator=(const ARM_VanillaSwaptionArg& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaArg::operator=(rhs);
		CleanUp();
        CopyNoCleanUp(rhs);
	}
	return *this;
}

ARM_VanillaSwaptionArg::~ARM_VanillaSwaptionArg()
{
	CleanUp();
};

ARM_Object* ARM_VanillaSwaptionArg::Clone() const
{
	return new ARM_VanillaSwaptionArg(*this);
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArg
///	Routine: CleanUp
///	Returns: 
///	Action : Delete the various data members
////////////////////////////////////////////////////
void ARM_VanillaSwaptionArg::CleanUp()
{
	delete itsFixPayTimes;
	delete itsFixPayPeriods;
	delete itsFloatResetTimes;
	delete itsFloatStartTimes;
	delete itsFloatEndTimes;
	delete itsFloatIntTerms; 
    delete itsStrikes;
	delete itsFixNominal;
	delete itsFloatNominal;
	delete itsFundingSpread;
}
////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArg
///	Routine: ImpliedVol
///	Returns: Exeption
///	Action :  No Implied Volatility 
////////////////////////////////////////////////////
double ARM_VanillaSwaptionArg::ImpliedVol(ARM_PricingModel* model) const
{
    /// create a dumState to avoid dummy price..
	double impliedVol;
	ARM_PricingFunctionIR* IRModel = dynamic_cast<ARM_PricingFunctionIR*>(model);
	
    double totalPrice = 0.0;
    if(IRModel)
    {
		impliedVol = IRModel->ImpliedVol(*this);
    }
    else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Model is not an interest rate model: not derived from ARM_PricingFunctionIR. So cannot price caplet, please advise");

    return impliedVol;
}
	
////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArg
///	Routine: Price
///	Returns: 
///	Action : price a swaption with a model checking that it
///				is derived from an interest rate model
////////////////////////////////////////////////////
double ARM_VanillaSwaptionArg::Price(ARM_PricingModel* model) const
{
	/// force to use a real dumStates that is not null!
	ARM_PricingFunctionIR* IRModel = dynamic_cast<ARM_PricingFunctionIR*>(model);
    ARM_VectorPtr Price;
	double price;
    ARM_PricingStatesPtr dumStates( new ARM_PricingStates(1,1,0) );    
    if(IRModel)
	{
		if (GetAtTheMoneyFlag())
		{
			ARM_PricingModelIR* newModel = dynamic_cast<ARM_PricingModelIR*>(model);
			ARM_VectorPtr annuity  = newModel->Annuity(GetCurveName(),GetEvalTime(),*itsFixPayTimes,*itsFixPayPeriods,dumStates);
#if defined(__GP_STRICT_VALIDATION)
			if ((*annuity)[0] < K_NEW_DOUBLE_TOL)
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					"The Swap Annuity is Null!");
#endif
			ARM_VectorPtr ZCStart	= newModel->GetDiscountFunctor()->DiscountFactor(GetCurveName(),GetEvalTime(),itsStartTime,dumStates);			
			ARM_VectorPtr ZCEnd	= newModel->GetDiscountFunctor()->DiscountFactor(GetCurveName(),GetEvalTime(),itsEndTime,dumStates);			
			double floatLeg = (*ZCStart)[0] - (*ZCEnd)[0];
			double newStrike = floatLeg/(*annuity)[0];
			int sizeVector = itsStrikes->size();
			for(int l=0; l<itsStrikes->size();++l)
				(*itsStrikes)[l] = newStrike;
			itsAtTheMoneyFlag = false;
		}
		bool closedFormulaFlag = IRModel->ClosedFormulaSwaptionFlag(itsIsConstantNotional,itsIsConstantSpread,itsIsConstantStrike);
		if(closedFormulaFlag)
		{
			Price = IRModel->VanillaSwaptionScalar(
			GetCurveName(),
			GetEvalTime(),
			itsResetTime,
			*itsFixNominal,
			*itsFloatNominal,
			itsStartTime,
			itsEndTime,
			*itsFloatResetTimes,
			*itsFloatStartTimes,
			*itsFloatEndTimes,
			*itsFloatIntTerms,
			*itsFixPayTimes,
			*itsFixPayPeriods,
			*itsStrikes,
			GetCallPut(),
			dumStates,
			itsIsConstantNotional,
			itsIsConstantSpread,
			itsIsConstantStrike);
			price = (*Price)[0];
		}
		else
		{
			if (model->GetNumMethod()==ARM_NumMethodPtr(NULL))
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					"The Model needs a numerical Method to Price this Swaption");

			ARM_GenSecurityPtr GenSecPtr = VanillaSwaptionToGenSec();
			(*itsGenSecString) = GenSecPtr->toString();
			ARM_GenPricer* genPricer = new ARM_GenPricer( &(*GenSecPtr),model );
			price = genPricer->Price();
			delete genPricer;
		}
	}
    else
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Model is not an interest rate model: not derived from ARM_PricingFunctionIR. So cannot price swaption, please advise");
  
    return price;
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArg
///	Routine: Price
///	Returns: 
///	Action : Compute the vega of the swaption
////////////////////////////////////////////////////

double ARM_VanillaSwaptionArg::Vega(ARM_PricingModel* model) const
{
	ARM_PricingModelIR* newModel = dynamic_cast<ARM_PricingModelIR*>(model);
	ARM_ZeroCurvePtr zcCurve = newModel->GetZeroCurve();
	double asOfDate = zcCurve->GetAsOfDate().GetJulian();

	double annuity = 0.0, floatLeg = 0.0;

	for (int i = 0; i < itsFloatStartTimes->size(); ++i)
	{
		floatLeg += (*itsFloatNominal)[i]*(zcCurve->DiscountPrice((*itsFloatStartTimes)[i]/K_YEAR_LEN));
		floatLeg += -(*itsFloatNominal)[i]*(zcCurve->DiscountPrice((*itsFloatEndTimes)[i]/K_YEAR_LEN));
	}

	for (i = 0; i < itsFixPayTimes->size(); ++i)
	{
		annuity += (*itsFixNominal)[i]*(*itsFixPayPeriods)[i]*(zcCurve->DiscountPrice((*itsFixPayTimes)[i]/K_YEAR_LEN));
	}

	double swapRate = floatLeg/annuity;
	double strike = (*itsStrikes)(0);

	double Expiry = itsResetTime / K_YEAR_LEN;

	double price = Price(model)/annuity;

	double vol = VanillaImpliedVol_N(swapRate, price, strike, Expiry, GetCallPut());

	return VegaVanillaOption_N(swapRate,vol,strike, Expiry, GetCallPut());
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArg
///	Routine: VanillaSwaptionToGenSec
///	Returns: 
///	Action : Returns a Generic Security Swaption for variable Nominal, Strikes or Spread Options
////////////////////////////////////////////////////
ARM_GenSecurityPtr ARM_VanillaSwaptionArg::VanillaSwaptionToGenSec() const
{
	vector<string> text( NBCOLUMNS*NBROWS );
	vector<ARM_GP_VALUE_TYPE> format(NBCOLUMNS*NBROWS, ARM_MISSING_TYPE );
	
	string PayRec = (GetCallPut()==K_CALL)?"P":"R";
	string fixFreq=ARM_ArgConvReverse_MatFrequency.GetString(GetFixFrequency());
	string floatFreq=ARM_ArgConvReverse_MatFrequency.GetString(GetFloatFrequency());
	string fixDayCount =ARM_ArgConvReverse_LgNameDayCount.GetString(GetFixDayCount());
	string floatDayCount =ARM_ArgConvReverse_LgNameDayCount.GetString(GetFloatDayCount());


	/// Names for th constant Manager
	vector<string> names(3,"");
	names[0] = "Notional";
	names[1] = "Spread";
	names[2] = "Strike";


	CC_NS(std,vector)<string>::iterator textIter			= text.begin();
	CC_NS(std,vector)<ARM_GP_VALUE_TYPE>::iterator formatIter	= format.begin();


	size_t descSize = NBCOLUMNS;
    vector< string > colNamesVec(descSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(descSize, ARM_STRING); 
    for(size_t i=0;i<descSize; ++i)
        colNamesVec[i] = SwaptionColNamesTable[i];
    ARM_RowInfo rowInfoNames(colNamesVec,colTypeVec);
	
    vector< string > rowDescVec(descSize);
    vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 
	
	CC_Ostringstream ResetDateString;
	double Reset = GetResetTime() + GetAsofDate();
	ResetDateString << setw(16) << CC_NS(std,fixed) << Reset; 
	rowDescVec[ResetDate] = ResetDateString.str();
    rowTypeVec[ResetDate] = ARM_DATE_TYPE;


	CC_Ostringstream StartDateString;
	double Start = GetStartTime() + GetAsofDate();
    StartDateString<<setw(16)<<CC_NS(std,fixed)<<Start;
	rowDescVec[StartDate] = StartDateString.str();
    rowTypeVec[StartDate] = ARM_DATE_TYPE;
	
	CC_Ostringstream PayDateString;
	double End = GetEndTime() + GetAsofDate();
    PayDateString<<setw(16)<<CC_NS(std,fixed)<<End;
	rowDescVec[PayDate] = PayDateString.str();
    rowTypeVec[PayDate] = ARM_DATE_TYPE;

	CC_Ostringstream SwapDesc;
	SwapDesc <<"Swap("<<GetCurveName()<<","<<SwaptionColNamesTable[StartDate]<<"[i],"<<SwaptionColNamesTable[PayDate]<<"[i],"<<names[2]<<","<<PayRec<<","<<fixFreq<<","<<fixDayCount<<","<<floatFreq<<","<<floatDayCount<<","<<names[1]<<","<<names[0]<<")";
	rowDescVec[Swap] = SwapDesc.str();
    rowTypeVec[Swap] = ARM_STRING;
	
	
	CC_Ostringstream OptionDesc;
	OptionDesc<<"MAX("<<SwaptionColNamesTable[Swap]<<"[i],0)";
	rowDescVec[Option] = OptionDesc.str();
	rowTypeVec[Option] = ARM_STRING;

	ARM_RowInfo RowInfoContents = ARM_RowInfo(rowDescVec,rowTypeVec);


	CC_NS( std,copy)( rowInfoNames.first.begin(), rowInfoNames.first.end(), textIter );
	CC_NS( std,copy)( rowInfoNames.second.begin(), rowInfoNames.second.end(), formatIter );
	textIter   += NBCOLUMNS;
	formatIter += NBCOLUMNS;

	CC_NS( std,copy)( RowInfoContents.first.begin(), RowInfoContents.first.end(), textIter );
	CC_NS( std,copy)( RowInfoContents.second.begin(), RowInfoContents.second.end(), formatIter );
	textIter   += NBCOLUMNS;
	formatIter += NBCOLUMNS;

	/// Object Manager
	vector<ARM_GramFctorArg> objVector;
	ARM_Object* object=NULL;

	/// Notional, Margin and Strikes are Vectors to avoid a second interpolaton
	if(!itsIsConstantNotional)
	{
		ARM_GP_CurvePtr notionalCurve = ARM_GP_CurvePtr (new ARM_Curve(*GetFixPayTimes(),*GetFixNotional(),new ARM_StepUpRightOpenCstExtrapolDble, false));
		//ARM_VectorPtr notionalVector((ARM_VectorPtr(static_cast<ARM_GP_Vector*> (GetFixNotional()->Clone()))));
		objVector.push_back(ARM_GramFctorArg(notionalCurve));
	}
	else
		objVector.push_back(ARM_GramFctorArg(GetNotional()));

	if(!itsIsConstantSpread)
	{
		ARM_VectorPtr spreadVector((ARM_VectorPtr(static_cast<ARM_GP_Vector*> (GetFundingSpread()->Clone()))));
		objVector.push_back(ARM_GramFctorArg(spreadVector));

	}
	else
		objVector.push_back(ARM_GramFctorArg(0.0));

	if(!itsIsConstantStrike)
	{
		ARM_VectorPtr strikesVector((ARM_GP_Vector*) GetStrikes()->Clone());
		objVector.push_back(ARM_GramFctorArg(strikesVector));
	}
	else
		objVector.push_back(ARM_GramFctorArg((*GetStrikes())[0]));	

	ARM_CstManagerPtr cstManager(new ARM_CstManager(names,objVector));

	ARM_DealDescriptionPtr pDealDescription( new ARM_DealDescription( text, format, NBROWS,NBCOLUMNS) );		
	return ARM_GenSecurityPtr( new ARM_GenSecurity( pDealDescription, GetCurveName(),cstManager)  );
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArg
///	Routine: VanillaSwaptionToGenSec
///	Returns: 
///	Action : Returns Column Names Table for the Generic Security
////////////////////////////////////////////////////
const string ARM_VanillaSwaptionArg::SwaptionColNamesTable [] =
{
    "ResetDate",
    "StartDate",
    "PayDate",
    "Swap",
    "Option"
};


////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArg
///	Routine: toString
///	Returns: string
///	Action : stringify the object to give details about it
////////////////////////////////////////////////////
string ARM_VanillaSwaptionArg::toString(const string& indent, const string& nextIndent) const
{ 
	return "ARM_VanillaSwaptionArg"; 
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

