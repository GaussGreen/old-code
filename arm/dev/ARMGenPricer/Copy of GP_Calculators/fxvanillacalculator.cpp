/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file FXVanillaCalculator.cpp
 *
 *  \brief file for the fxvanilla calculator 
 *	\author  K.Belkheir & E.Ezzine
 *	\version 1.0
 *	\date January 2007
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/fxvanillacalculator.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/autocleaner.h"
#include "gpbase/ostringstream.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/datestrip.h"
#include "gpbase/singleton.h"
#include "gpbase/globalconstant.h"
#include "gpbase/cloneutilityfunc.h"


#include "gpcalib/calibmethod.h"

#include "gpnummethods/cfmethod.h"

/// gpinfra
#include "gpinfra/mktdatamanagerrep.h"
#include "gpinfra/argconvdefault.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/cstmanager.h"
#include "gpinfra/discretisationscheme.h"

/// gpmodels
#include "gpmodels/HybridIRFX.h"
#include "gpmodels/fxname.h"


/// kernel
#include <inst/forex.h>

/// STL
#include <iomanip> /// for setprecision()
#include <list>
CC_USING_NS(std,list)

CC_BEGIN_NAMESPACE( ARM )

/// Default MDM key names
const string YC_BASIS_KEY_NAME              = "YC_BASIS_";
const string FXMODEL_KEY_NAME               = "FXMOD_";
const string IRMODEL_KEY_NAME				= "IRMOD_";
const string CORREL_KEY_NAME                = "CORREL_";


/////////////////////////////////////////////////////////////////
///	Class  : ARM_FXVanillaCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_FXVanillaCalculator::ARM_FXVanillaCalculator( const ARM_FXVanillaCalculator& rhs )
:	
ARM_GenCalculator( rhs ),
	itsDateStrip(CreateClonedPtr(&*rhs.itsDateStrip)),     
	itsStartDate(rhs.itsStartDate),          
	itsEndDate(rhs.itsEndDate),   
	itsDayCount(rhs.itsDayCount),        
	itsFrequency(rhs.itsFrequency),   
	itsFX1name(rhs.itsFX1name),        
	itsFX2name(rhs.itsFX2name),        
	itsPayCcy(rhs.itsPayCcy),  
	itsExpiryGap(rhs.itsExpiryGap),
	itsSetlmentGap(rhs.itsSetlmentGap),
	itsPaymentGap(rhs.itsPaymentGap),
	itsResetCalendar(rhs.itsResetCalendar),
	itsPayCalendar(rhs.itsPayCalendar),
	itsCallPut(rhs.itsCallPut),
	itsCallPut2(rhs.itsCallPut2),
	itsNominalCv(rhs.itsNominalCv),
	itsStrikeCv(rhs.itsStrikeCv),
	itsAlphaCv(rhs.itsAlphaCv),
	itsBetaCv(rhs.itsBetaCv),
	itsStrike2Cv(rhs.itsStrike2Cv),
	itsLeverageCv(rhs.itsLeverageCv),
	itsVanillaType(rhs.itsVanillaType),
	itsMinMax(rhs.itsMinMax),
	itsDigitType(rhs.itsDigitType),
	itsEpsilon(rhs.itsEpsilon),
	itsKmin1Cv(rhs.itsKmin1Cv),
	itsKmax1Cv(rhs.itsKmax1Cv),
	itsKmin2Cv(rhs.itsKmin2Cv),
	itsKmax2Cv(rhs.itsKmax2Cv),
	itsNleft1(rhs.itsNleft1),
	itsNcenter1(rhs.itsNcenter1),
	itsNright1(rhs.itsNright1),
	itsNleft2(rhs.itsNleft2),
	itsNcenter2(rhs.itsNcenter2),
	itsNright2(rhs.itsNright2),
	itsIntermediatePrices(CreateClone(rhs.itsIntermediatePrices)),
	itsHasBeenComputed(rhs.itsHasBeenComputed)
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_FXVanillaCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////
ARM_FXVanillaCalculator::ARM_FXVanillaCalculator(const ARM_Date& asOfDate,
		const ARM_DateStrip& dateStrip,
		const ARM_Date& startDate,
		const ARM_Date& endDate,
		int expiryGap,
		int setlmentGap,
		int paymentGap,
		int frequency,
		int dayCount, 
		string resetCal,
		string payCal,
		const string fx1Name,
		const string fx2Name,
		const ARM_Currency& payCcy,
		const ARM_Curve& nominal,
		int callPut,
		const ARM_Curve& strike,
		const ARM_Curve& alpha,
		const ARM_Curve& beta,
		const ARM_Curve& strike2,
		const ARM_Curve& leverage,
		int callPut2,
		ARM_BasketType minMax,
		ARM_DigitType digitType,
		double epsilon,
		ARM_VanillaType vanillaType)
:	
	ARM_GenCalculator(asOfDate),
	itsDateStrip( CreateClonedPtr(const_cast<ARM_DateStrip*>(&dateStrip)) ),
	itsStartDate(startDate),          
	itsEndDate(endDate),                      
	itsExpiryGap(expiryGap),
	itsSetlmentGap(setlmentGap),
	itsPaymentGap(paymentGap),
	itsFrequency(frequency),  
	itsDayCount(dayCount),
	itsFX1name(fx1Name),
	itsFX2name(fx2Name),
	itsPayCcy(payCcy),
	itsResetCalendar(resetCal),
	itsPayCalendar(payCal),
	itsNominalCv(nominal),
	itsCallPut(callPut),
	itsStrikeCv(strike),
	itsAlphaCv(alpha),
	itsBetaCv(beta),
	itsStrike2Cv(strike2),
	itsLeverageCv(leverage),
	itsCallPut2(callPut2),
	itsMinMax(minMax),
	itsDigitType(digitType),
	itsEpsilon(epsilon),
	itsVanillaType(vanillaType),
	itsHasBeenComputed(false),
	itsIntermediatePrices(NULL)
{
	DatesStructure();

	// Set Keys
	ARM_FXName  fxName1(fx1Name.substr(fx1Name.size()-6));    
	ARM_FXName  fxName2(fx2Name.substr(fx2Name.size()-6));
	string Fx1Fx2Name("(" + fxName1.GetMktName() + "," + fxName2.GetMktName() + ")" );
	string payCcyName( payCcy.GetCcyName() );
    ARM_StringVector mdmKeys(NbKeys);

	mdmKeys[YcBasisPayKey]    = YC_BASIS_KEY_NAME + payCcyName;
    mdmKeys[Fx1ModelKey]      = FXMODEL_KEY_NAME  + fxName1.GetMktName();
	mdmKeys[Fx2ModelKey]      = FXMODEL_KEY_NAME  + fxName2.GetMktName();
	mdmKeys[IrConvAdjModelKey]= IRMODEL_KEY_NAME  + payCcyName;
    mdmKeys[CorrelMatrixKey]  = CORREL_KEY_NAME   + Fx1Fx2Name;

    SetKeys(mdmKeys);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_FXVanillaCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds calculator from FxSpreadStripOption
/////////////////////////////////////////////////////////////////
ARM_FXVanillaCalculator::ARM_FXVanillaCalculator(const ARM_FxSpreadStripOption& FxStrip,
		ARM_BasketType minMax,
		ARM_DigitType digitType,
		ARM_VanillaType vanillaType)
:	
	ARM_GenCalculator(FxStrip.GetAsOfDate()),
	itsDateStrip( CreateClonedPtr(FxStrip.GetSchedule()) ),
	itsStartDate(FxStrip.GetStartDate()),          
	itsEndDate(FxStrip.GetExpiryDate()),                      
	itsExpiryGap(FxStrip.GetExpiryGap()),
	itsSetlmentGap(FxStrip.GetSetlmentGap()),
	itsPaymentGap(FxStrip.GetPaymentGap()),
	itsFrequency(FxStrip.GetFrequency()),  
	itsDayCount(FxStrip.GetDayCount()),
	itsFX1name(FxStrip.GetFx1Name()),
	itsFX2name(FxStrip.GetFx2Name()),
	itsPayCcy(FxStrip.GetPaymentCcy()),
	itsResetCalendar(string(FxStrip.GetSchedule()->GetResetCalendar())),
	itsPayCalendar(string(FxStrip.GetSchedule()->GetPayCalendar())),
	itsNominalCv(FxStrip.GetNotional()),
	itsCallPut(FxStrip.GetOptionType()),
	itsStrikeCv(FxStrip.GetStrikes()),
	itsAlphaCv(FxStrip.GetAlpha()),
	itsBetaCv(FxStrip.GetBeta()),
	itsLeverageCv(FxStrip.GetLeverage()), 
	itsMinMax(minMax),
	itsDigitType(digitType),
	itsEpsilon(FxStrip.GetEpsilon()),
	itsVanillaType(vanillaType),
	itsHasBeenComputed(false),
	itsIntermediatePrices(NULL)
{
	// default
	itsCallPut2		= itsCallPut;
	itsStrike2Cv	= itsStrikeCv;

	DatesStructure();

	// Set Keys
	ARM_FXName  fxName1(itsFX1name.substr(itsFX1name.size()-6));    
	ARM_FXName  fxName2(itsFX2name.substr(itsFX2name.size()-6));
	string Fx1Fx2Name("(" + fxName1.GetMktName() + "," + fxName2.GetMktName() + ")" );
	string payCcyName( itsPayCcy.GetCcyName() );
    ARM_StringVector mdmKeys(NbKeys);

	mdmKeys[YcBasisPayKey]    = YC_BASIS_KEY_NAME + payCcyName;
    mdmKeys[Fx1ModelKey]      = FXMODEL_KEY_NAME  + fxName1.GetMktName();
	mdmKeys[Fx2ModelKey]      = FXMODEL_KEY_NAME  + fxName2.GetMktName();
	mdmKeys[IrConvAdjModelKey]= IRMODEL_KEY_NAME  + payCcyName;
    mdmKeys[CorrelMatrixKey]  = CORREL_KEY_NAME   + Fx1Fx2Name;

    SetKeys(mdmKeys);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_FXVanillaCalculator
///	Routine: init
///	Returns: void
///	Action : initialize the fxvanillacalculator calculator
/////////////////////////////////////////////////////////////////
void ARM_FXVanillaCalculator::Init(	const ARM_MarketData_ManagerRep& mktDataManager,
				const ARM_Curve& kmin1Cv,
				const ARM_Curve& kmax1Cv,
				const ARM_Curve& kmin2Cv,
				const ARM_Curve& kmax2Cv,
				int nLeft1,
				int nCenter1,
				int nRight1,
				int nLeft2,
				int nCenter2,
				int nRight2)
{
	/// Register input objects but no more than the internal number of access keys
    size_t nbToReg = GetKeys().size();
    for(size_t i(0); i<nbToReg; ++i)
	{
		if(!mktDataManager.TestIfKeyMissing(GetKeys()[i]))
			GetMktDataManager()->RegisterData(GetKeys()[i],mktDataManager.GetData(GetKeys()[i]));
	}
	GetMktDataManager()->SetDetailMode(mktDataManager.GetDetailMode());

	/// Check market datas
    CheckMktDataAndTimeIt();

	/// Create a HybridIRFX model
    CreateAndSetModelAndTimeIt();

	/// Set the NumMethod
	itsKmin1Cv = kmin1Cv;
	itsKmax1Cv = kmax1Cv;
	itsKmin2Cv = kmin2Cv;
	itsKmax2Cv = kmax2Cv;

	itsNleft1   = nLeft1;
	itsNcenter1 = nCenter1;
	itsNright1  = nRight1;
	itsNleft2   = nLeft2;
	itsNcenter2 = nCenter2;
	itsNright2  = nRight2;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_FXVanillaCalculator
///	Routine: CheckMktData
///	Returns: void
///	Action : check if FXVanilla market datas are consistent
/////////////////////////////////////////////////////////////////
void ARM_FXVanillaCalculator::CheckMktData()
{
    /// Market datas checking
	ARM_EqFxBase* modelfx = dynamic_cast< ARM_EqFxBase* >( GetMktDataManager()->GetData(GetKeys()[Fx1ModelKey]) );
    if(!modelfx)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : FX model for key=" + GetKeys()[Fx1ModelKey] + " is expected in the Market Data Manager");

	modelfx = dynamic_cast< ARM_EqFxBase* >( GetMktDataManager()->GetData(GetKeys()[Fx2ModelKey]) );
    if(!modelfx)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : FX model for key=" + GetKeys()[Fx2ModelKey] + " is expected in the Market Data Manager");

	ARM_PricingModelIR* modelir = dynamic_cast< ARM_PricingModelIR* >( GetMktDataManager()->GetData(GetKeys()[IrConvAdjModelKey]) );
    if(!modelir)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : IR model for key=" + GetKeys()[IrConvAdjModelKey] + " is expected in the Market Data Manager");

	ARM_Curve* correlCurve = dynamic_cast< ARM_Curve* >( GetMktDataManager()->GetData(GetKeys()[CorrelMatrixKey]) );
	if(!correlCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Correl matrix curve for key=" + GetKeys()[CorrelMatrixKey] +" is expected in the Market Data Manager");
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_FXVanillaCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : Create and Set the model with its numMethod
/////////////////////////////////////////////////////////////////
void ARM_FXVanillaCalculator::CreateAndSetModel()
{
	int size = 3;
	ARM_StringVector names(size);
	ARM_StringVectorVector depends(size);
	vector< ARM_PricingModelPtr > models(size);
	//mktdata
	ARM_EqFxBase* modelfx1 = dynamic_cast< ARM_EqFxBase* >( GetMktDataManager()->GetData(GetKeys()[Fx1ModelKey]) );	
	ARM_EqFxBase* modelfx2 = dynamic_cast< ARM_EqFxBase* >( GetMktDataManager()->GetData(GetKeys()[Fx2ModelKey]) );
	ARM_ZeroCurve*  payCurve = dynamic_cast< ARM_ZeroCurve* >( GetMktDataManager()->GetData(GetKeys()[YcBasisPayKey]) );
	ARM_PricingModelIR*  modelir = dynamic_cast< ARM_PricingModelIR* >( GetMktDataManager()->GetData(GetKeys()[IrConvAdjModelKey]) );

	models[0] = ARM_PricingModelPtr((ARM_EqFxBase*)modelfx1->Clone());
	models[1] = ARM_PricingModelPtr((ARM_EqFxBase*)modelfx2->Clone());
	models[2] = ARM_PricingModelPtr((ARM_PricingModelIR*)modelir->Clone());

	ARM_FXName  fxName1(itsFX1name.substr(itsFX1name.size()-6));    
	ARM_FXName  fxName2(itsFX2name.substr(itsFX2name.size()-6));
	names[0] = fxName1.GetMktName();
	names[1] = fxName2.GetMktName();
	names[2] = itsPayCcy.GetCcyName();
	ARM_Curve* correlFX1FX2Curve = dynamic_cast< ARM_Curve* >( GetMktDataManager()->GetData(GetKeys()[CorrelMatrixKey]) );

	size_t sizeCv = correlFX1FX2Curve->size();
	ARM_GP_Vector ordinates(sizeCv,1.0);
	ARM_Curve diagonal(correlFX1FX2Curve->GetAbscisses(),ordinates,new ARM_LinInterpCstExtrapolDble);
	vector<ARM_Curve*>  curveVector(9,&diagonal);//(IR,FX1,FX2)x(IR,FX1,FX2)
	curveVector[5]= correlFX1FX2Curve;
	curveVector[7]= correlFX1FX2Curve;
	ARM_CurveMatrix curveMatrix(curveVector,3,3);
	
	/// Create a modelnamemap
	ARM_ModelNameMap modelMap( names, models, depends );

	/// Finally, We create and set the brand HybridIRFX !
    ARM_PricingModelPtr model = ARM_PricingModelPtr( new ARM_HybridIRFX( modelMap,curveMatrix) );

	model->SetNumMethod( ARM_NumMethodPtr( new ARM_CFMethod()));
	model->SetPayModelName(string(itsPayCcy.GetCcyName()) );
	model->SetZeroCurve(ARM_ZeroCurvePtr((ARM_ZeroCurve*)payCurve->Clone() ));

	/// Set the model
	SetPricingModel(model);
}

////////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCCalculator
///	Routine: ComputePricingData
///	Returns: 
///	Action : pricing function called from the addin GetPricingData()
////////////////////////////////////////////////////////////////////
void ARM_FXVanillaCalculator::ComputePricingData() const
{
	if (!itsHasBeenComputed)
		const_cast<ARM_FXVanillaCalculator*>(this)->PriceAndTimeIt();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_FXVanillaCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create theDatestrip 
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_FXVanillaCalculator::DatesStructure() const 
{
	int fwdRule	=	K_MOD_FOLLOWING;	// for forward dates
	int intRule	=	K_ADJUSTED;			// for interest dates
	int resetTiming	=	K_ARREARS;
	int payTiming	=	K_ARREARS;
	int stubRule	=	K_SHORTSTART;


	ARM_INDEX_TYPE indexType = ((ARM_Currency*)GetCurrencyUnit())->GetVanillaIndexType();
	char* DefaultresetCalendar = GetCurrencyUnit()->GetResetCalName(indexType);
	CC_NS(std,auto_ptr)<char> holdresetCalendar(DefaultresetCalendar);
	char* DefaultpayCalendar  = GetCurrencyUnit()->GetPayCalName(indexType);
	CC_NS(std,auto_ptr)<char> holdpayCalendar(DefaultpayCalendar);

	const char* resetCalendar	= itsResetCalendar == "" ? DefaultresetCalendar : itsResetCalendar.c_str();
	const char* payCalendar		= itsPayCalendar == "" ? DefaultpayCalendar: itsPayCalendar.c_str();

	if (itsDateStrip == ARM_DateStripPtr(NULL))
	{
		ARM_DateStrip Sched(itsStartDate,itsEndDate,itsFrequency,itsDayCount,resetCalendar,fwdRule,intRule,
			stubRule,-fabs(itsExpiryGap),itsFrequency,GETDEFAULTVALUE, payCalendar,payTiming,payTiming);
		const_cast< ARM_FXVanillaCalculator* >(this)->itsDateStrip = ARM_DateStripPtr( new ARM_DateStrip(Sched));
	}
	ARM_DateStripVector SchedVect(1);
    SchedVect[0] = &*itsDateStrip;
    ARM_DateStripCombiner EventSchedule(SchedVect,"ResetDate");
	
	return EventSchedule;
};
/////////////////////////////////////////////////////////////////
///	Class  : ARM_FXVanillaCalculator
///	Routine: Price
///	Returns: a double
///	Action : price the FXVanilla option
/////////////////////////////////////////////////////////////////
double ARM_FXVanillaCalculator::Price()
{

	ARM_Curve* correlCv = dynamic_cast< ARM_Curve* >( GetMktDataManager()->GetData(GetKeys()[CorrelMatrixKey]) );
	double asofdate = GetMktDataManager()->GetAsOfDate().GetJulian();
	double evalTime = 0.0;
	
	// loop for the pricing of the strip
	ARM_GP_Vector* resetDates = itsDateStrip->GetResetDates() ;
	double cumprice = 0;
	//intermediatePrices
	size_t size = resetDates->size();
	itsIntermediatePrices= new ARM_Vector(size, 0.0);
	
	//model
	ARM_HybridIRFX* model = dynamic_cast <ARM_HybridIRFX*> (&*GetPricingModel());
	ARM_CFMethod* cfnumMethod = dynamic_cast<ARM_CFMethod*>(&*model->GetNumMethod());
	ARM_GP_Matrix integParameters(5,2);
	integParameters(2,0)=itsNcenter1;
	integParameters(3,0)=itsNleft1;
	integParameters(4,0)=itsNright1;
	integParameters(2,1)=itsNcenter2;
	integParameters(3,1)=itsNleft2;
	integParameters(4,1)=itsNright2;

	//loop
	for ( size_t i(0); i<size; i++)
	{
		double expiryTime  = (*resetDates)[i] - asofdate;
		double strike   = itsStrikeCv.Interpolate(expiryTime);
		ARM_GP_Vector strikePerState(1,strike);
		double alpha = itsAlphaCv.Interpolate(expiryTime);
		double beta = itsBetaCv.Interpolate(expiryTime);
		ARM_Date expiryDate((*resetDates)[i]);
		double settlementTime = expiryDate.GapBusinessDay(itsSetlmentGap, const_cast <char*> (itsResetCalendar.c_str()) ).GetJulian() - asofdate;
		double payTime =expiryDate.GapBusinessDay(itsPaymentGap, const_cast <char*> (itsResetCalendar.c_str()) ).GetJulian() - asofdate;
		double nominal = itsNominalCv.Interpolate(expiryTime);
		ARM_GP_Vector* interestTerms = itsDateStrip->GetInterestTerms() ;
		double intTerm  = (*interestTerms)[i];
		double leverage = itsLeverageCv.Interpolate(expiryTime);
		double coef = intTerm*nominal*leverage;

		//NumMethod setting
		double kmin1 = itsKmin1Cv.Interpolate(expiryTime);
		double kmax1 = itsKmax1Cv.Interpolate(expiryTime);
		double kmin2 = itsKmin2Cv.Interpolate(expiryTime);
		double kmax2 = itsKmax2Cv.Interpolate(expiryTime);

		integParameters(0,0)=kmin1;
		integParameters(1,0)=kmax1;
		integParameters(0,1)=kmin2;
		integParameters(1,1)=kmax2;

		cfnumMethod->SetCFParameters(integParameters);

		if(itsVanillaType == ARM_FXVanillaType::spread)
		{
			ARM_GP_VectorPtr price = model->CallSpreadVectorial(
			itsFX1name,
			itsFX2name,
			evalTime,
			expiryTime,
			settlementTime,
			settlementTime,
			payTime,
			strikePerState,
			alpha,
			beta,	
			ARM_PricingStatesPtr(NULL) );

			(*itsIntermediatePrices)[i] = coef*(*price)[0];//storage for ComputeAll
			cumprice += (*itsIntermediatePrices)[i];
		}
		if(itsVanillaType == ARM_FXVanillaType::quotient)
		{
			ARM_GP_VectorPtr price = model->CallQuotientVectorial(
			itsFX1name,
			itsFX2name,
			evalTime,
			expiryTime,
			settlementTime,
			settlementTime,
			payTime,
			strikePerState,
			ARM_PricingStatesPtr(NULL) );

			(*itsIntermediatePrices)[i] = coef*(*price)[0];//storage for ComputeAll
			cumprice += (*itsIntermediatePrices)[i];
		}
		else if(itsVanillaType == ARM_FXVanillaType::vanilla)//vanilla in the final version
		{
			ARM_GP_VectorPtr price = model-> CallVectorial(itsFX1name,
				evalTime,
				expiryTime,
				settlementTime,
				strikePerState,
				itsCallPut,
				payTime,
				ARM_PricingStatesPtr(NULL) );

			(*itsIntermediatePrices)[i] = coef*(*price)[0];//storage for ComputeAll
			cumprice += (*itsIntermediatePrices)[i];
		}
		else if(itsVanillaType == ARM_FXVanillaType::digit)
		{
			ARM_GP_VectorPtr price = model->DigitalVectorial(itsFX1name,
				evalTime,
				expiryTime,
				settlementTime,
				strikePerState,
				1.0,
				itsCallPut,
				payTime,
				itsDigitType,
				itsEpsilon,
				ARM_PricingStatesPtr(NULL) );

			(*itsIntermediatePrices)[i] = coef*(*price)[0];//storage for ComputeAll
			cumprice += (*itsIntermediatePrices)[i];
		}
		else if(itsVanillaType == ARM_FXVanillaType::digitspread)
		{
			double signed_eps = itsCallPut*itsEpsilon;
			double norm = 1.0/itsEpsilon;
			ARM_GP_Vector strikeDownPerState, strikeUpPerState;

			switch (itsDigitType)
			{
			case ARM_FXDigitType::backward:
				{
					strikeDownPerState = strikePerState - signed_eps;
					strikeUpPerState = strikePerState;
					break;
				}
			case ARM_FXDigitType::forward:
				{
					strikeDownPerState = strikePerState;
					strikeUpPerState = strikePerState + signed_eps;
					break;
				}
			case ARM_FXDigitType::centred:
				{
					norm = 1/(2*itsEpsilon);
					strikeDownPerState = strikePerState - signed_eps;
					strikeUpPerState = strikePerState + signed_eps;
					break;
				}
			}
			ARM_GP_VectorPtr callDown = model->CallSpreadVectorial(itsFX1name,
				itsFX2name,
				evalTime,
				expiryTime,
				settlementTime,
				settlementTime,
				payTime,
				strikeDownPerState,
				alpha,
				beta,	
				ARM_PricingStatesPtr(NULL) );
			ARM_GP_VectorPtr callUp = model->CallSpreadVectorial(itsFX1name,
				itsFX2name,
				evalTime,
				expiryTime,
				settlementTime,
				settlementTime,
				payTime,
				strikeUpPerState,
				alpha,
				beta,	
				ARM_PricingStatesPtr(NULL) );
			 (*itsIntermediatePrices)[i] = coef*norm*( (*callDown)[0] - (*callUp)[0] );//storage for ComputeAll
			 cumprice += (*itsIntermediatePrices)[i];
		}
		else if(itsVanillaType == ARM_FXVanillaType::perf)
		{
			double barrier  = itsStrike2Cv.Interpolate(expiryTime);
			ARM_GP_Vector invbarrierPerState(1,1.0/barrier);
			ARM_GP_VectorPtr call = model-> CallVectorial(itsFX1name,
				evalTime,
				expiryTime,
				settlementTime,
				invbarrierPerState,
				itsCallPut,
				payTime,
				ARM_PricingStatesPtr(NULL) );
			ARM_GP_VectorPtr digit = model-> DigitalVectorial(itsFX1name,
				evalTime,
				expiryTime,
				settlementTime,
				invbarrierPerState,
				1.0,
				itsCallPut,
				payTime,
				itsDigitType,
				itsEpsilon,
				ARM_PricingStatesPtr(NULL) );

			(*itsIntermediatePrices)[i] = coef*( strike*(*call)[0] + itsCallPut*(strike/barrier-1)*(*digit)[0] );//storage for ComputeAll
			cumprice += (*itsIntermediatePrices)[i];
		}
	}
	itsHasBeenComputed = true;
	return cumprice;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_FXVanillaCalculator
///	Routine: ComputeAll
///	Returns: an ARM_Vector
///	Action : return the intermediate prices
/////////////////////////////////////////////////////////////////
ARM_Vector* ARM_FXVanillaCalculator::ComputeAll()
{	
	if(!itsHasBeenComputed)
		double x = Price();
		
	return (ARM_Vector*)itsIntermediatePrices->Clone();
}
////////////////////////////////////////////////////
///	Class   : ARM_FXVanillaCalculator
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_FXVanillaCalculator::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream prdcData;

    /// GenCalculator general datas dump
	prdcData << "\n\nCommon Calculator Part\n" << ARM_GenCalculator::toString(indent,nextIndent) << "\n\n";

    return prdcData.str();
}


ARM_DateStripCombiner ARM_FXVanillaCalculator::CustomDatesStructure(const ARM_DateStripVector& dateStrips ) const
{
	return ARM_DateStripCombiner();
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

