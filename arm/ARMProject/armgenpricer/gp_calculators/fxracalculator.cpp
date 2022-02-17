/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file FXRACalculator.cpp
 *
 *  \brief file for the fxvanilla calculator 
 *	\author  K.Belkheir & E.Ezzine
 *	\version 1.0
 *	\date January 2007
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/FXRACalculator.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/autocleaner.h"
#include "gpbase/ostringstream.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/datestrip.h"
#include "gpbase/singleton.h"
#include "gpbase/globalconstant.h"
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/stringconvert.h"

/// gpcalib
#include "gpcalib/calibmethod.h"
/// gpnummethods
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
//#include <inst/forex.h>
#include <util/fromto.h>

/// STL
#include <iomanip> /// for setprecision()
#include <list>
CC_USING_NS(std,list)

CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////////////////////////
///	Class  : ARM_FXRACalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_FXRACalculator::ARM_FXRACalculator( const ARM_FXRACalculator& rhs )
:	
ARM_FXVanillaCalculator( rhs ),
	itsIntRule(rhs.itsIntRule),
	itsStubType(rhs.itsStubType),
	itsResetTiming(rhs.itsResetTiming),
	itsFixingfrequency(rhs.itsFixingfrequency),
	itsPayIdx(rhs.itsPayIdx),
	itsPayIdxSpread(rhs.itsPayIdxSpread),
	itsPayIdxIT(rhs.itsPayIdxIT),
	itsIrIdx(rhs.itsIrIdx),
	itsIrIdxIT(rhs.itsIrIdxIT),
	itsFxDownBarrierCv(rhs.itsFxDownBarrierCv),
	itsFxUpBarrierCv(rhs.itsFxUpBarrierCv),
	itsIrDownBarrierCv(rhs.itsIrDownBarrierCv),
	itsIrUpBarrierCv(rhs.itsIrUpBarrierCv),
	itsFixingTimes(rhs.itsFixingTimes),
	itsIRindexResetTimes(rhs.itsIRindexResetTimes),
	itsIRindexStartTimes(rhs.itsIRindexStartTimes),
	itsIRindexEndTimes(rhs.itsIRindexStartTimes),
	itsIRindexTerms(rhs.itsIRindexTerms),
	itsAllFixingPrices(CreateClone(rhs.itsAllFixingPrices))
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_FXRACalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////
ARM_FXRACalculator::ARM_FXRACalculator(const ARM_Date& asOfDate,
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
		ARM_VanillaType vanillaType,
		int	intRule,
		int stubType,
		int resetTiming,
		int fixingfrequency,
		string payIdx,
		double payIdxSpread,
		string payIdxIT,
		string irIdx,
		string irIdxIT,
		const ARM_Curve& fxDownBarrierCv,
		const ARM_Curve& fxUpBarrierCv,
		const ARM_Curve& irDownBarrierCv,
		const ARM_Curve& irUpBarrierCv)
:	
	ARM_FXVanillaCalculator(asOfDate,
		dateStrip,
		startDate,
		endDate,
		payCcy,
		expiryGap,
		setlmentGap,
		frequency,
		dayCount, 
		resetCal,
		payCal,
		fx1Name,
		fx2Name,
		leverage,
		nominal,
		strike,
		callPut,	
		vanillaType,
		alpha,
		beta,
		digitType,
		epsilon,
		strike2,
		callPut2,
		minMax),
		itsIntRule(intRule),
		itsStubType(stubType),
		itsResetTiming(resetTiming),
		itsFixingfrequency(fixingfrequency),
		itsPayIdx(payIdx),
		itsPayIdxSpread(payIdxSpread),
		itsPayIdxIT(payIdxIT),
		itsIrIdx(irIdx),
		itsIrIdxIT(irIdxIT),
		itsFxDownBarrierCv(fxDownBarrierCv),
		itsFxUpBarrierCv(fxUpBarrierCv),
		itsIrDownBarrierCv(irDownBarrierCv),
		itsIrUpBarrierCv(irUpBarrierCv),
		itsAllFixingPrices(NULL)
{
		ARM_FXVanillaCalculator::DatesStructure();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_FXRACalculator
///	Routine: init
///	Returns: void
///	Action : initialize the FXRACalculator calculator
/////////////////////////////////////////////////////////////////
void ARM_FXRACalculator::Init(	const ARM_MarketData_ManagerRep& mktDataManager,
				int nbPoints1,
				int nbPoints2)
{

	//Temporally because the strikes will disappear someday
	ARM_Curve kmin1Cv = ARM_FlatCurve(1e-4);
	ARM_Curve kmin2Cv = ARM_FlatCurve(1e-4);
	ARM_Curve kmax1Cv = ARM_FlatCurve(1e+4);
	ARM_Curve kmax2Cv = ARM_FlatCurve(1e+4);
	int nLeft1 = 0;
	int nCenter1 = nbPoints1;
	int nRight1 = 0;
	int nLeft2 = 0;
	int nCenter2 = nbPoints2;
	int nRight2 = 0;

	ARM_FXVanillaCalculator::Init(  mktDataManager,
									kmin1Cv,
									kmax1Cv,
									kmin2Cv,
									kmax2Cv,
									nLeft1,
									nCenter1,
									nRight1,
									nLeft2,
									nCenter2,
									nRight2);
	/// Check market datas
    CheckMktDataAndTimeIt();

	/// Create a HybridIRFX model
    CreateAndSetModelAndTimeIt();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_FXRACalculator
///	Routine: CheckMktData
///	Returns: void
///	Action : check if FXVanilla market datas are consistent
/////////////////////////////////////////////////////////////////
void ARM_FXRACalculator::CheckMktData()
{
    /// Market datas checking
	ARM_EqFxBase* modelfx = dynamic_cast< ARM_EqFxBase* >( GetMktDataManager()->GetData(GetKeys()[Fx1ModelKey]) );
    if(!modelfx)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : FX model for key=" + GetKeys()[Fx1ModelKey] + " is expected in the Market Data Manager");

	modelfx = dynamic_cast< ARM_EqFxBase* >( GetMktDataManager()->GetData(GetKeys()[Fx2ModelKey]) );
    if(!modelfx)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : FX model for key=" + GetKeys()[Fx2ModelKey] + " is expected in the Market Data Manager");

	ARM_PricingModelIR* modelir = dynamic_cast< ARM_PricingModelIR* >( GetMktDataManager()->GetData(GetKeys()[IrModelKey]) );
    if(!modelir)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : IR model for key=" + GetKeys()[IrModelKey] + " is expected in the Market Data Manager");

	ARM_Curve* correlCurve = dynamic_cast< ARM_Curve* >( GetMktDataManager()->GetData(GetKeys()[CorrelMatrixKey]) );
	if(!correlCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Correl matrix curve for key=" + GetKeys()[CorrelMatrixKey] +" is expected in the Market Data Manager");
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_FXRACalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : Create and Set the model with its numMethod
/////////////////////////////////////////////////////////////////
void ARM_FXRACalculator::CreateAndSetModel()
{
	ARM_FXVanillaCalculator::CreateAndSetModel();
}

////////////////////////////////////////////////////////////////////
///	Class  : ARM_FXRACalculator
///	Routine: ComputePricingData
///	Returns: 
///	Action : pricing function called from the addin GetPricingData()
////////////////////////////////////////////////////////////////////
void ARM_FXRACalculator::ComputePricingData() const
{
	if (!itsHasBeenComputed)
		const_cast<ARM_FXRACalculator*>(this)->PriceAndTimeIt();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_FXRACalculator
///	Routine: FixingStructure
///	Returns: nothing
///	Action : create the fixings data 
/////////////////////////////////////////////////////////////////
void ARM_FXRACalculator::FixingsStructure(ARM_Date startDate, ARM_Date endDate, char* ccyName, ARM_IRIndex irIndex, string irIndexTerm, ARM_PricingModel* model)  
{
	/// Creation of the intermediate datestrip for fixingDates
	ARM_INDEX_TYPE defaultIndex = LIBOR3M;//GetDefaultIndexFromCurrency( ccyName );
	char* payCalendar	= itsPayCcy.GetPayCalName(defaultIndex);
	char* resetCalendar = itsPayCcy.GetResetCalName(defaultIndex);
	int freq		    = itsFixingfrequency;
    int dayCount        = itsDayCount;
	double resetGap=0;
	double paymentGap=0;
    ARM_DateStripPtr sched(NULL);
	sched = ARM_DateStripPtr( new ARM_DateStrip(startDate, endDate, freq, dayCount, resetCalendar,
		K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTEND, resetGap, freq, paymentGap, payCalendar,K_ARREARS) );
	
	/// handle fixingTimes, IRindexStartTimes, IRindexEndTimes and IRindexTerms
	size_t i,nbFixing = sched->GetResetDates()->size();
    itsFixingTimes = std::vector<double>(nbFixing);
	itsIRindexResetTimes = std::vector<double>(nbFixing);
	itsIRindexStartTimes = std::vector<double>(nbFixing);
	itsIRindexEndTimes = std::vector<double>(nbFixing);
	itsIRindexTerms = std::vector<double>(nbFixing);

	for(i=0;i<nbFixing;++i)
	{
		//Fixing Times
		itsFixingTimes[i] = model->GetTimeFromDate((*(sched->GetResetDates()))[i]);
		//IRindexStartTimes
		ARM_Date fwdStart((*(sched->GetResetDates()))[i]);//I use the reset date?
		itsIRindexStartTimes[i] = itsFixingTimes[i];
		//IRindexEndTimes
		ARM_Date fwdEnd((*(sched->GetResetDates()))[i]);//I use the reset date?
		//fwdEnd.AddPeriod(irIndexType);
		fwdEnd.AddPeriod(irIndexTerm);
		fwdEnd.GoodBusinessDay(irIndex.GetFwdRule() * irIndex.GetIntRule(), resetCalendar); 
		itsIRindexEndTimes[i] = model->GetTimeFromDate(fwdEnd);
		//IRindexResetTimes
		itsIRindexResetTimes[i] = itsIRindexEndTimes[i];//for the moment ir index reset in arrears and at EndTimes
		//IRindexTerms
		itsIRindexTerms[i] = CountYearsWithoutException( dayCount, fwdStart, fwdEnd );//CountBusinessDays(fwdStart, fwdEnd, payCalendar);
	}

	/// delete char* for memory leak
	delete payCalendar;
	delete resetCalendar;
};

/////////////////////////////////////////////////////////////////
///	Class  : ARM_FXRACalculator
///	Routine: Price
///	Returns: a double
///	Action : price the FXVanilla option
/////////////////////////////////////////////////////////////////
double ARM_FXRACalculator::Price()
{
	ARM_Curve* correlCv = dynamic_cast< ARM_Curve* >( GetMktDataManager()->GetData(GetKeys()[CorrelMatrixKey]) );
	double asofdate = GetMktDataManager()->GetAsOfDate().GetJulian();
	double evalTime = 0.0;

	//dates
	std::vector<double>* startDates = itsDateStrip->GetFwdRateStartDates();
	std::vector<double>* endDates = itsDateStrip->GetFwdRateEndDates();
	std::vector<double>* payDates = itsDateStrip->GetPaymentDates();
	std::vector<double>* interestTerms = itsDateStrip->GetInterestTerms() ;

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

	//Vector issued from Curves
	std::vector<double> fxDownBarriers, fxUpBarriers, irDownBarriers, irUpBarriers, notionals;
	
	//PAY and IR index creation
	char* ccyName			= itsPayCcy.GetCcyName();
	//PAY
	double payIndexTerm		= StringMaturityToYearTerm(itsPayIdxIT);
	//string payIndexTermStr	= ConvertYearTermToStringMatu(payIndexTerm);
	int payIndexType		= ARM_ArgConv_IndexClass.GetNumber(itsPayIdx);
	//IR
	double irIndexTerm		= StringMaturityToYearTerm(itsIrIdxIT);
	string irIndexTermStr	= ConvertYearTermToStringMatu(irIndexTerm);
	int irIndexType_temp		= ARM_ArgConv_IndexClass.GetNumber(itsIrIdx);
	ARM_INDEX_TYPE irIndexType;
	if (irIndexType_temp == K_LIBOR)
		irIndexType = (ARM_INDEX_TYPE)FromIndexAndTermToIndexType(irIndexTermStr, ccyName);
	else
		irIndexType = EURIBOR1Y;
	//ARM_IRIndex irIndex(irIndexType);

	//intermediatePrices
	size_t size = startDates->size();
	itsIntermediatePrices= new ARM_Vector(size, 0.0);
	
	itsAllFixingPrices= new ARM_Vector();
	
	// loop for the pricing of the strip
	double cumprice = 0;
	for ( size_t i(0); i<size; i++)
	{
		//Fixing info
		ARM_Date startDate = (*startDates)[i];
		ARM_Date endDate = (*endDates)[i];
		ARM_Date payDate = (*payDates)[i];
		FixingsStructure(startDate, endDate, ccyName, ccyName/*irIndex*/, irIndexTermStr, model);

		//from dates to time
		double startTime = (*startDates)[i] - asofdate;
		double endTime   = (*endDates)[i] - asofdate;
		double payTime   = (*payDates)[i] - asofdate;

		//from curve to vector
		double lag = startTime;
		double temp = itsFxDownBarrierCv.Interpolate(lag);
		fxDownBarriers.push_back(temp);
		temp = itsFxUpBarrierCv.Interpolate(temp);
		fxUpBarriers.push_back(temp);
		temp = itsIrDownBarrierCv.Interpolate(lag);
		irDownBarriers.push_back(temp);
		temp = itsIrUpBarrierCv.Interpolate(lag);
		irUpBarriers.push_back(temp);
		temp = itsLeverageCv.Interpolate(lag);
		temp = temp*itsNominalCv.Interpolate(lag);
		notionals.push_back(temp);
	
		//NumMethod setting
		double kmin1 = itsKmin1Cv.Interpolate(startTime);
		double kmax1 = itsKmax1Cv.Interpolate(startTime);
		double kmin2 = itsKmin2Cv.Interpolate(startTime);
		double kmax2 = itsKmax2Cv.Interpolate(startTime);
		integParameters(0,0)=kmin1;
		integParameters(1,0)=kmax1;
		integParameters(0,1)=kmin2;
		integParameters(1,1)=kmax2;
		cfnumMethod->SetCFParameters(integParameters);
	
		//price
		ARM_GP_VectorPtr price = model->RangeAccrualVectorial(
			ccyName,
			evalTime,
			startTime,
			endTime,
			payTime,
			itsFixingTimes,
			payIndexType, 
			payIndexTerm,
			itsFX1name,
			irIndexType, 
			itsIRindexResetTimes,
			itsIRindexStartTimes,
			itsIRindexEndTimes,
			itsIRindexTerms,
			fxDownBarriers,
			fxUpBarriers,
			irDownBarriers,
			irUpBarriers,
			notionals,
			ARM_PricingStatesPtr(NULL),
			itsAllFixingPrices);

			(*itsIntermediatePrices)[i] = (*price)[0] * (*interestTerms)[i];//storage for ComputeAll
			cumprice += (*itsIntermediatePrices)[i];
	}
	itsHasBeenComputed = true;
	return cumprice;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_FXRACalculator
///	Routine: ComputeAll
///	Returns: an ARM_Vector
///	Action : return the intermediate prices
/////////////////////////////////////////////////////////////////
ARM_Vector* ARM_FXRACalculator::ComputeAll()
{	
	if(!itsHasBeenComputed)
		double x = Price();
		
	return (ARM_Vector*)itsIntermediatePrices->Clone();
}
////////////////////////////////////////////////////
///	Class   : ARM_FXRACalculator
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_FXRACalculator::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream prdcData;

    /// GenCalculator general datas dump
	prdcData << "\n\nCommon Calculator Part\n" << ARM_GenCalculator::toString(indent,nextIndent) << "\n\n";

    return prdcData.str();
}

/*
ARM_DateStripCombiner ARM_FXRACalculator::CustomDatesStructure(const ARM_DateStripVector& dateStrips ) const
{
	return ARM_DateStripCombiner();
}
*/

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

