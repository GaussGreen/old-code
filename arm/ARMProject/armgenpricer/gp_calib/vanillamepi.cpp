/*!
 *
 * Copyright (c) CDC IXIS CM June 2005 Paris
 *
 *	\file vanillamepi.cpp
 *
 *  \brief vanilla mepi
 *
 *	\author  A SCHAULY
 *	\version 1.0
 *	\date June 2005
 */

#include "gpbase/removeidentifiedwarning.h"

/// gpcalib
#include "gpcalib/vanillamepi.h"
#include "gpcalib/vanillaargnumeric.h"

#include "gpcalib/typedef.h"

/// gpnummethods
#include "gpinfra/gensecurity.h"
#include "gpinfra/pricingmodel.h"

#include "gpbase/port.h"
#include "gpinfra/typedef.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/pricingmodelir.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/genpricer.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/gpvector.h"
#include "gpbase/checkarg.h"
#include "gpbase/singleton.h"
#include "gpbase/datestrip.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/countedptr.h"
#include "gpbase/autocleaner.h"

/// gpnummethods
#include "gpnummethods/treemethod.h"
#include "gpnummethods/treebase.h"
#include "gpnummethods/treefactory.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/reconnector.h"
#include "gpnummethods/smoother.h"

CC_BEGIN_NAMESPACE( ARM )

const int TREE_NBSTEPS_PER_YEAR	=30;
const double STD_DEV_RATIO		=5.0;
const double MIN_STD_DEV		=0.001;


const string ARM_VanillaMepi::VanillaMepiColNamesTable [] =
{
		"EventDate",
		"NextDate",
		"MaturityDate",
		"Underlying",
		"CapiFactor",
		"PaidCapiFactor",
		"Cash",
		"Portfolio",
		"AvgPortfolio",
		"FixedProtection",
		"ActifRisque",
		"CF"
};

const string ARM_VanillaMepi::VanillaMepiConstantsTable [] = 
{
	"RiskFactor",
	"Strike",
	"MaxBorrow",
	"StartingPortfolio",
	"StartingCash",
	"MinInvested",
	"LeverageCost",
	"CashSpread",
	"Fees",
	"AlreadyAsianed",
	"AsianDatesNb"
};

////////////////////////////////////////////////////
///	Struct : ARM_VanillaMepi
///	Routine: default constructor, copy constructor,
///          assigment, destructor, clone
///	Returns: 
///	Action :
////////////////////////////////////////////////////
void ARM_VanillaMepi::CopyNoCleanUp(const ARM_VanillaMepi& rhs)
{
	itsEquityModelName=rhs.itsEquityModelName;
	itsStartDate=rhs.itsStartDate;
	itsEndDate=rhs.itsEndDate;
	itsStrike=rhs.itsStrike;
	itsMaxBorrow=rhs.itsMaxBorrow;
	itsProtectionCurveStart=rhs.itsProtectionCurveStart;
	itsProtectionCurveEnd=rhs.itsProtectionCurveEnd;
	itsStartingPortfolio=rhs.itsStartingPortfolio;
	itsStartingCash=rhs.itsStartingCash;
	itsMinInvested=rhs.itsMinInvested;
	itsLeverageCost=rhs.itsLeverageCost;
	itsCashSpread=rhs.itsCashSpread;
	itsFees=rhs.itsFees;
	itsBorrowingRate=rhs.itsBorrowingRate;
	itsRiskFactor=rhs.itsRiskFactor;
	itsInitialProtection=rhs.itsInitialProtection;
	itsFinalProtection=rhs.itsFinalProtection;
	itsAlreadyAsianed=rhs.itsAlreadyAsianed;
	itsAsianingPeriodNb=rhs.itsAsianingPeriodNb;
	itsMaturity=rhs.itsMaturity;
	itsProtectionStep=rhs.itsProtectionStep;
	itsEventDates=rhs.itsEventDates;


}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaMepi
///	Routine: constructor
///          
///	Returns: 
///	Action :
////////////////////////////////////////////////////


ARM_VanillaMepi::ARM_VanillaMepi(const string& curveName, const string& equityModelName,	double startDate,	double endDate, long resetFreq, 
		double riskFactor, double strike, double maxBorrow, double protectionCurveStart, double protectionCurveEnd, double startingPortfolio, double startingCash, 
		double minInvested, double leverageCost, double cashSpread, double fees, double alreadyAsianed,	long AsianingPeriodNb) : 
	ARM_VanillaArgNumeric(curveName, startDate, 0,endDate), 
		itsStrike(strike), itsFees(fees), itsLeverageCost(leverageCost), itsStartingCash(startingCash), 
		itsEquityModelName( equityModelName ),itsStartingPortfolio(startingPortfolio), itsCashSpread(cashSpread), 
		itsMinInvested(minInvested), itsStartDate(startDate), itsEndDate(endDate), itsProtectionCurveStart(protectionCurveStart), 
		itsProtectionCurveEnd(protectionCurveEnd), itsRiskFactor(riskFactor), itsMaxBorrow(maxBorrow), itsAsianingPeriodNb(AsianingPeriodNb), itsAlreadyAsianed(alreadyAsianed)
{
	ARM_DateStrip dateStrip( itsStartDate, itsEndDate, resetFreq, 3, curveName.c_str(), K_MOD_FOLLOWING, K_UNADJUSTED, K_SHORTSTART, 0, resetFreq, 0, curveName.c_str(), K_ARREARS, K_ARREARS, false );

	itsEventDates = ARM_GP_VectorPtr(new ARM_GP_Vector(*dateStrip.GetResetDates()));//static_cast<ARM_GP_Vector*> ( dateStrip.GetResetDates()->Clone() ));

	if( fabs( itsEventDates->Elt(0) - itsStartDate ) > 5 )
		itsEventDates->insert( itsEventDates->begin(), itsStartDate );

	itsEventDates->insert( itsEventDates->begin(), 1 ); 

	itsMaturity = endDate - startDate;
	itsProtectionStep = ( itsProtectionCurveEnd - itsProtectionCurveStart ) / itsMaturity;
}


ARM_VanillaMepi& ARM_VanillaMepi::operator=(const ARM_VanillaMepi& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaArg::operator=(rhs);
	 	CleanUp();
		CopyNoCleanUp(rhs);
	}
	return *this;
}

ARM_VanillaMepi::ARM_VanillaMepi(const ARM_VanillaMepi& rhs) : ARM_VanillaArgNumeric( rhs )
{
	CopyNoCleanUp(rhs);
}


void ARM_VanillaMepi::CleanUp()
{
	
}



ARM_VanillaMepi::~ARM_VanillaMepi()
{
	CleanUp();
}

ARM_Object* ARM_VanillaMepi::Clone() const
{
	return new ARM_VanillaMepi(*this);
}


double ARM_VanillaMepi::ComputeDelta(ARM_PricingModel* model, double bump ) const
{
	const ARM_CstManagerPtr cstManager( getConstantManager() );
	ARM_CstManager::iterator pfiter = cstManager->find( VanillaMepiConstantsTable[StartingPortfolio] );
	double portfolioValue = (*pfiter).second->GetDouble();
	(*pfiter).second->SetDouble( portfolioValue + bump );
	double Delta = Price( model );
	(*pfiter).second->SetDouble( portfolioValue );
	Delta -= Price( model );
	return Delta;
}


ARM_RowInfo ARM_VanillaMepi::NumArgColumnNames() const
{
	size_t colNamesSize = sizeof(VanillaMepiColNamesTable)/sizeof(VanillaMepiColNamesTable[0]);
    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0;i<colNamesSize; ++i)
        colNamesVec[i] = VanillaMepiColNamesTable[i];

    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
}



ARM_RowInfo ARM_VanillaMepi::NumArgMiddleRows( size_t eventIdx, const ARM_GP_VectorPtr& eventDates) const
{
    size_t descSize = sizeof(VanillaMepiColNamesTable)/sizeof(VanillaMepiColNamesTable[0]);
	size_t nbEvents = eventDates->size();

    vector< string > rowDescVec(descSize);
    vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 

	string curveName = GetCurveName();
	double firstDate = eventDates->Elt(1);
	
	// EventDate
    double eventDate = eventDates->Elt(eventIdx);
    CC_Ostringstream eventDateDesc;
    eventDateDesc << CC_NS(std,fixed) << eventDate;
    rowDescVec[EventDate] = eventDateDesc.str();
    rowTypeVec[EventDate] = ARM_DATE_TYPE;

	// NextDate
    double nextEventDate = ( eventIdx < nbEvents-1 ) ? eventDates->Elt(eventIdx+1) : eventDates->Elt(eventIdx);
    CC_Ostringstream nextDateDesc;
    nextDateDesc << CC_NS(std,fixed) << nextEventDate;
    rowDescVec[NextDate] = nextDateDesc.str();
    rowTypeVec[NextDate] = ARM_DATE_TYPE;

	// MaturityDate
    double MatDate = eventDates->Elt(nbEvents-1);
    CC_Ostringstream MaturityDateDesc;
    MaturityDateDesc << CC_NS(std,fixed) << MatDate;
    rowDescVec[MaturityDate] = MaturityDateDesc.str();
    rowTypeVec[MaturityDate] = ARM_DATE_TYPE;

	// Underlying
	CC_Ostringstream UnderlyingDesc;
	UnderlyingDesc << "Spot(" << itsEquityModelName <<")";
	rowDescVec[Underlying] = UnderlyingDesc.str();
	rowTypeVec[Underlying] = ARM_STRING;

	// Capifactor
	CC_Ostringstream CapiFactorDesc;
	CapiFactorDesc << "DF(" << curveName << "," << VanillaMepiColNamesTable[NextDate] << "[i])";
	rowDescVec[CapiFactor] = CapiFactorDesc.str();
	rowTypeVec[CapiFactor] = ARM_STRING;

	// PaidCapiFactor
	CC_Ostringstream PaidCapiFactorDesc;
	PaidCapiFactorDesc << "IF(" << VanillaMepiColNamesTable[Cash] << "[i]<0," << VanillaMepiColNamesTable[CapiFactor] << "[i]/(1+" 
		<< VanillaMepiConstantsTable[LeverageCost];
	PaidCapiFactorDesc << "*" << VanillaMepiColNamesTable[CapiFactor] << "[i]/12)," << VanillaMepiColNamesTable[CapiFactor] 
		<< "[i]/(1-" << VanillaMepiConstantsTable[CashSpread];
	PaidCapiFactorDesc << "*" << "CapiFactor" << "[i]/12))";
	rowDescVec[PaidCapiFactor] = PaidCapiFactorDesc.str();
	rowTypeVec[PaidCapiFactor] = ARM_STRING;

	// Cash
	CC_Ostringstream CashDesc;
	(eventIdx == 1) ? (CashDesc << VanillaMepiConstantsTable[StartingCash] ) : (CashDesc << "(" << VanillaMepiColNamesTable[Portfolio] 
		<< "[i-1]-" << VanillaMepiColNamesTable[ActifRisque] << "[i-1])/" << VanillaMepiColNamesTable[PaidCapiFactor] << "[i-1]" );
	rowDescVec[Cash] = CashDesc.str();
	rowTypeVec[Cash] = ARM_STRING;

    // Portfolio
	CC_Ostringstream PortfolioDesc;
	(eventIdx == 1) ? ( PortfolioDesc << VanillaMepiConstantsTable[StartingPortfolio] ) : ( PortfolioDesc << VanillaMepiColNamesTable[ActifRisque] 
		<< "[i-1]*" << VanillaMepiColNamesTable[Underlying] << "[i]/" << VanillaMepiColNamesTable[Underlying] << "[i-1]+" 
		<< VanillaMepiColNamesTable[Cash] << "[i]-Fees/12"); 
	rowDescVec[Portfolio] = PortfolioDesc.str();
	rowTypeVec[Portfolio] = ARM_STRING;

	// AvgPortfolio
	// The value associated with AlreadyAsianed should be equal to 0 , when not in the asianing period
	CC_Ostringstream AvgPortfolioDesc;
	if(eventIdx == nbEvents-1)
	{
		size_t asianingNb = MIN( itsAsianingPeriodNb , nbEvents-1 );
		AvgPortfolioDesc << "(" << VanillaMepiConstantsTable[AlreadyAsianed] << "+"<< VanillaMepiColNamesTable[Portfolio] << "[i]";
		for( long i=1 ; i<asianingNb;++i)
			AvgPortfolioDesc << "+" << VanillaMepiColNamesTable[Portfolio] << "[i-" << i << "]";

		AvgPortfolioDesc << ")/" << VanillaMepiConstantsTable[AsianDatesNb];

		rowTypeVec[AvgPortfolio] = ARM_STRING;
	}
	else
	{
		AvgPortfolioDesc << "0";
		rowTypeVec[AvgPortfolio] = ARM_DOUBLE;
	}
	rowDescVec[AvgPortfolio] = AvgPortfolioDesc.str();

	/// FixedProtection
	CC_Ostringstream FixedProtectionDesc;
	double protection = itsProtectionCurveStart + ( eventDates->Elt(eventIdx) - itsStartDate ) * itsProtectionStep;
	FixedProtectionDesc << CC_NS(std,fixed) << protection;
	rowDescVec[FixedProtection] = FixedProtectionDesc.str();
	rowTypeVec[FixedProtection] = ARM_DOUBLE;
	
	/// ActifRisque
	CC_Ostringstream ActifRisqueDesc;
	ActifRisqueDesc << "Min(Max(" <<VanillaMepiConstantsTable[MinInvested] << "*" << VanillaMepiColNamesTable[Portfolio] 
		<< "[i]," << VanillaMepiConstantsTable[RiskFactor] << "*(" << VanillaMepiColNamesTable[Portfolio] << "[i]-" 
		<< VanillaMepiColNamesTable[FixedProtection] << "[i])),";
	ActifRisqueDesc << VanillaMepiColNamesTable[Portfolio] << "[i]+" <<VanillaMepiConstantsTable[MaxBorrow] << ")";
	rowDescVec[ActifRisque] = ActifRisqueDesc.str();
	rowTypeVec[ActifRisque] = ARM_STRING;
	
	/// CF
	CC_Ostringstream CFDesc;
	if( eventIdx == nbEvents-1 ) 
	{
		CFDesc << "Max(" << VanillaMepiColNamesTable[AvgPortfolio] << "[i]-" << VanillaMepiConstantsTable[Strike] << ",0)";
		rowTypeVec[CF] = ARM_STRING;
	}
	else
	{
		CFDesc << "0";
		rowTypeVec[CF] = ARM_DOUBLE;
	}
	rowDescVec[CF] = CFDesc.str();

	return ARM_RowInfo(rowDescVec,rowTypeVec);
}



ARM_GP_VectorPtr ARM_VanillaMepi::DatesStructure() const
{
	return ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>( itsEventDates->Clone() ) );
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaArgNumeric
///	Routine: securityCstManager
///	Returns: 
///	Action : gives the constant manager of the security
////////////////////////////////////////////////////
const ARM_CstManagerPtr ARM_VanillaMepi::securityCstManager() const
{
	size_t i;
	size_t cstManagerSize = sizeof(VanillaMepiConstantsTable)/sizeof(VanillaMepiConstantsTable[0]);
	vector<string> names(cstManagerSize);
	vector<double> values(cstManagerSize);

	for( i=0 ; i < cstManagerSize ; ++i )
		names[i] = VanillaMepiConstantsTable[i];

	values[RiskFactor] = itsRiskFactor;
	values[Strike] = itsStrike;
	values[MaxBorrow] = itsMaxBorrow;
	values[StartingPortfolio] = itsStartingPortfolio;
	values[StartingCash] = itsStartingCash;
	values[MinInvested] = itsMinInvested;
	values[LeverageCost] = itsLeverageCost;
	values[CashSpread] = itsCashSpread;
	values[Fees] = itsFees;
	values[AlreadyAsianed] = itsAlreadyAsianed;
	values[AsianDatesNb] = (double) itsAsianingPeriodNb;

// FIXMEFRED: mig.vc8 (22/05/2007 18:11:21): cast
	return static_cast<ARM_CstManagerPtr>(new ARM_CstManager(names,values));
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaArgNumeric
///	Routine: AnalyticPrice
///	Returns: 
///	Action : analytic price of spreadoption with a model checking 
////////////////////////////////////////////////////
double ARM_VanillaMepi::AnalyticPrice(ARM_PricingModel* model) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : No Analytic Price for MEPI " );

}
////////////////////////////////////////////////////
///	Struct : ARM_VanillaMepi
///	Routine: toString
///	Returns: string
///	Action : stringify the object to give details about it
////////////////////////////////////////////////////
string ARM_VanillaMepi::toString(const string& indent, const string& nextIndent) const
{ 
	return GetGenSecurity()->toString();
}

double ARM_VanillaMepi::Price(ARM_PricingModel* model ) const
{
	ARM_GenPricer genPricer( &*GetGenSecurity(), model );
	return genPricer.Price();
}



////////////////////////////////////////////////////////////////////////////////////////////


///                                       VanillaMepiDelta


////////////////////////////////////////////////////////////////////////////////////////////


const string ARM_VanillaMepiDelta::VanillaMepiColNamesTable [] =
{
		"EventDate",
		"NextDate",
		"MaturityDate",
		"Underlying",
		"CapiFactor",
		"FixedProtection",
		"PaidCapiFactor1",
		"Cash1",
		"Portfolio1",
		"AvgPortfolio1",
		"ActifRisque1",
		"PaidCapiFactor2",
		"Cash2",
		"Portfolio2",
		"AvgPortfolio2",
		"ActifRisque2",
		"CF"
};

const string ARM_VanillaMepiDelta::VanillaMepiConstantsTable [] = 
{
	"RiskFactor",
	"Strike",
	"MaxBorrow",
	"StartingPortfolio",
	"StartingCash",
	"MinInvested",
	"LeverageCost",
	"CashSpread",
	"Fees",
	"AlreadyAsianed",
	"DeltaShift",
	"AsianDatesNb"
};

////////////////////////////////////////////////////
///	Struct : ARM_VanillaMepiDelta
///	Routine: default constructor, copy constructor,
///          assigment, destructor, clone
///	Returns: 
///	Action :
////////////////////////////////////////////////////
void ARM_VanillaMepiDelta::CopyNoCleanUp(const ARM_VanillaMepiDelta& rhs)
{
	itsEquityModelName=rhs.itsEquityModelName;
	itsStartDate=rhs.itsStartDate;
	itsEndDate=rhs.itsEndDate;
	itsStrike=rhs.itsStrike;
	itsMaxBorrow=rhs.itsMaxBorrow;
	itsProtectionCurveStart=rhs.itsProtectionCurveStart;
	itsProtectionCurveEnd=rhs.itsProtectionCurveEnd;
	itsStartingPortfolio=rhs.itsStartingPortfolio;
	itsStartingCash=rhs.itsStartingCash;
	itsMinInvested=rhs.itsMinInvested;
	itsLeverageCost=rhs.itsLeverageCost;
	itsCashSpread=rhs.itsCashSpread;
	itsFees=rhs.itsFees;
	itsBorrowingRate=rhs.itsBorrowingRate;
	itsRiskFactor=rhs.itsRiskFactor;
	itsInitialProtection=rhs.itsInitialProtection;
	itsFinalProtection=rhs.itsFinalProtection;
	itsAlreadyAsianed=rhs.itsAlreadyAsianed;
	itsAsianingPeriodNb=rhs.itsAsianingPeriodNb;
	itsMaturity=rhs.itsMaturity;
	itsProtectionStep=rhs.itsProtectionStep;
	itsEventDates=rhs.itsEventDates,itsEventDates;
	itsDeltaShift=rhs.itsDeltaShift;


}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaMepiDelta
///	Routine: constructor
///          
///	Returns: 
///	Action :
////////////////////////////////////////////////////


ARM_VanillaMepiDelta::ARM_VanillaMepiDelta(const string& curveName, const string& equityModelName,	double startDate,	double endDate, long resetFreq, 
		double riskFactor, double strike, double maxBorrow, double protectionCurveStart, double protectionCurveEnd, double startingPortfolio, double startingCash, 
		double minInvested, double leverageCost, double cashSpread, double fees, double alreadyAsianed,long AsianingPeriodNb,double deltaShift) : 
	ARM_VanillaArgNumeric(curveName, startDate, 0,endDate), 
		itsStrike(strike), itsFees(fees), itsLeverageCost(leverageCost), itsStartingCash(startingCash), 
		itsEquityModelName( equityModelName ),itsStartingPortfolio(startingPortfolio), itsCashSpread(cashSpread), 
		itsMinInvested(minInvested), itsStartDate(startDate), itsEndDate(endDate), itsProtectionCurveStart(protectionCurveStart), 
		itsProtectionCurveEnd(protectionCurveEnd), itsRiskFactor(riskFactor), itsMaxBorrow(maxBorrow), 
		itsAsianingPeriodNb(AsianingPeriodNb), itsAlreadyAsianed(alreadyAsianed),itsDeltaShift(deltaShift)
{
	ARM_DateStrip dateStrip( itsStartDate, itsEndDate, resetFreq, 3, curveName.c_str(), K_MOD_FOLLOWING, K_UNADJUSTED, K_SHORTSTART, 0, resetFreq, 0, curveName.c_str(), K_ARREARS, K_ARREARS, false );

	itsEventDates = ARM_GP_VectorPtr(new ARM_GP_Vector(*dateStrip.GetResetDates()));//ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*> ( dateStrip.GetResetDates()->Clone() ));

	if( fabs( itsEventDates->Elt(0) - itsStartDate ) > 5 )
		itsEventDates->insert( itsEventDates->begin(), itsStartDate );

	itsEventDates->insert( itsEventDates->begin(), 1 ); 

	itsMaturity = endDate - startDate;
	itsProtectionStep = ( itsProtectionCurveEnd - itsProtectionCurveStart ) / itsMaturity;
}


ARM_VanillaMepiDelta& ARM_VanillaMepiDelta::operator=(const ARM_VanillaMepiDelta& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaArg::operator=(rhs);
	 	CleanUp();
		CopyNoCleanUp(rhs);
	}
	return *this;
}

ARM_VanillaMepiDelta::ARM_VanillaMepiDelta(const ARM_VanillaMepiDelta& rhs) : ARM_VanillaArgNumeric( rhs )
{
	CopyNoCleanUp(rhs);
}


void ARM_VanillaMepiDelta::CleanUp()
{
	
}



ARM_VanillaMepiDelta::~ARM_VanillaMepiDelta()
{
	CleanUp();
}

ARM_Object* ARM_VanillaMepiDelta::Clone() const
{
	return new ARM_VanillaMepiDelta(*this);
}


double ARM_VanillaMepiDelta::ComputeDelta(ARM_PricingModel* model, double bump ) const
{
	const ARM_CstManagerPtr cstManager( getConstantManager() );
	ARM_CstManager::iterator pfiter = cstManager->find( VanillaMepiConstantsTable[StartingPortfolio] );
	double portfolioValue = (*pfiter).second->GetDouble();
	(*pfiter).second->SetDouble( portfolioValue + bump );
	double Delta = Price( model );
	(*pfiter).second->SetDouble( portfolioValue );
	Delta -= Price( model );
	return Delta;
}


ARM_RowInfo ARM_VanillaMepiDelta::NumArgColumnNames() const
{
	size_t colNamesSize = sizeof(VanillaMepiColNamesTable)/sizeof(VanillaMepiColNamesTable[0]);
    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0;i<colNamesSize; ++i)
        colNamesVec[i] = VanillaMepiColNamesTable[i];

    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
}



ARM_RowInfo ARM_VanillaMepiDelta::NumArgMiddleRows( size_t eventIdx, const ARM_GP_VectorPtr& eventDates) const
{
    size_t descSize = sizeof(VanillaMepiColNamesTable)/sizeof(VanillaMepiColNamesTable[0]);
	size_t nbEvents = eventDates->size();

    vector< string > rowDescVec(descSize);
    vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 

	string curveName = GetCurveName();
	double firstDate = eventDates->Elt(1);
	
	// EventDate
    double eventDate = eventDates->Elt(eventIdx);
    CC_Ostringstream eventDateDesc;
    eventDateDesc << CC_NS(std,fixed) << eventDate;
    rowDescVec[EventDate] = eventDateDesc.str();
    rowTypeVec[EventDate] = ARM_DATE_TYPE;

	// NextDate
    double nextEventDate = ( eventIdx < nbEvents-1 ) ? eventDates->Elt(eventIdx+1) : eventDates->Elt(eventIdx);
    CC_Ostringstream nextDateDesc;
    nextDateDesc << CC_NS(std,fixed) << nextEventDate;
    rowDescVec[NextDate] = nextDateDesc.str();
    rowTypeVec[NextDate] = ARM_DATE_TYPE;

	// MaturityDate
    double MatDate = eventDates->Elt(nbEvents-1);
    CC_Ostringstream MaturityDateDesc;
    MaturityDateDesc << CC_NS(std,fixed) << MatDate;
    rowDescVec[MaturityDate] = MaturityDateDesc.str();
    rowTypeVec[MaturityDate] = ARM_DATE_TYPE;

	// Underlying
	CC_Ostringstream UnderlyingDesc;
	UnderlyingDesc << "Spot(" << itsEquityModelName <<")";
	rowDescVec[Underlying] = UnderlyingDesc.str();
	rowTypeVec[Underlying] = ARM_STRING;

	// Capifactor
	CC_Ostringstream CapiFactorDesc;
	CapiFactorDesc << "DF(" << curveName << "," << VanillaMepiColNamesTable[NextDate] << "[i])";
	rowDescVec[CapiFactor] = CapiFactorDesc.str();
	rowTypeVec[CapiFactor] = ARM_STRING;

	/// FixedProtection
	CC_Ostringstream FixedProtectionDesc;
	double protection = itsProtectionCurveStart + ( eventDates->Elt(eventIdx) - itsStartDate ) * itsProtectionStep;
	FixedProtectionDesc << CC_NS(std,fixed) << protection;
	rowDescVec[FixedProtection] = FixedProtectionDesc.str();
	rowTypeVec[FixedProtection] = ARM_DOUBLE;


	// PaidCapiFactor1
	CC_Ostringstream PaidCapiFactorDesc1;
	PaidCapiFactorDesc1 << "IF(" << VanillaMepiColNamesTable[Cash1] << "[i]<0," << VanillaMepiColNamesTable[CapiFactor] << "[i]/(1+" << VanillaMepiConstantsTable[LeverageCost];
	PaidCapiFactorDesc1 << "*" << VanillaMepiColNamesTable[CapiFactor] << "[i]/12)," << VanillaMepiColNamesTable[CapiFactor] << "[i]/(1-" << VanillaMepiConstantsTable[CashSpread];
	PaidCapiFactorDesc1 << "*" << "CapiFactor" << "[i]/12))";
	rowDescVec[PaidCapiFactor1] = PaidCapiFactorDesc1.str();
	rowTypeVec[PaidCapiFactor1] = ARM_STRING;

	// Cash1
	CC_Ostringstream CashDesc1;
	(eventIdx == 1) ? (CashDesc1 << VanillaMepiConstantsTable[StartingCash] ) : (CashDesc1 << "(" << VanillaMepiColNamesTable[Portfolio1] << "[i-1]-" << VanillaMepiColNamesTable[ActifRisque1] << "[i-1])/" << VanillaMepiColNamesTable[PaidCapiFactor1] << "[i-1]" );
	rowDescVec[Cash1] = CashDesc1.str();
	rowTypeVec[Cash1] = ARM_STRING;

    // Portfolio1
	CC_Ostringstream PortfolioDesc1;
	(eventIdx == 1) ? ( PortfolioDesc1 << VanillaMepiConstantsTable[StartingPortfolio] << "+" << VanillaMepiConstantsTable[DeltaShift] << "/2" ) : ( PortfolioDesc1 << VanillaMepiColNamesTable[ActifRisque1] << "[i-1]*" << VanillaMepiColNamesTable[Underlying] << "[i]/" << VanillaMepiColNamesTable[Underlying] << "[i-1]+" << VanillaMepiColNamesTable[Cash1] << "[i]-Fees/12"); 
	rowDescVec[Portfolio1] = PortfolioDesc1.str();
	rowTypeVec[Portfolio1] = ARM_STRING;

	// AvgPortfolio1
	// The value associated with AlreadyAsianed should be equal to 0 , when not in the asianing period
	CC_Ostringstream AvgPortfolioDesc1;
	if(eventIdx == nbEvents-1)
	{
		size_t asianingNb = MIN( itsAsianingPeriodNb , nbEvents-1 );
		AvgPortfolioDesc1 << "(" << VanillaMepiConstantsTable[AlreadyAsianed] << "+"<< VanillaMepiColNamesTable[Portfolio1] << "[i]";
		for( long i=1 ; i<asianingNb;++i)
			AvgPortfolioDesc1 << "+" << VanillaMepiColNamesTable[Portfolio1] << "[i-" << i << "]";

		AvgPortfolioDesc1 << ")/" << VanillaMepiConstantsTable[AsianDatesNb];

		rowTypeVec[AvgPortfolio1] = ARM_STRING;
	}
	else
	{
		AvgPortfolioDesc1 << "0";
		rowTypeVec[AvgPortfolio1] = ARM_DOUBLE;
	}
	rowDescVec[AvgPortfolio1] = AvgPortfolioDesc1.str();

	
	/// ActifRisque1
	CC_Ostringstream ActifRisqueDesc1;
	ActifRisqueDesc1 << "Min(Max(" <<VanillaMepiConstantsTable[MinInvested] << "*" << VanillaMepiColNamesTable[Portfolio1] << "[i]," << VanillaMepiConstantsTable[RiskFactor] << "*(" << VanillaMepiColNamesTable[Portfolio1] << "[i]-" << VanillaMepiColNamesTable[FixedProtection] << "[i])),";
	ActifRisqueDesc1 << VanillaMepiColNamesTable[Portfolio1] << "[i]+" <<VanillaMepiConstantsTable[MaxBorrow] << ")";
	rowDescVec[ActifRisque1] = ActifRisqueDesc1.str();
	rowTypeVec[ActifRisque1] = ARM_STRING;

	// PaidCapiFactor2
	CC_Ostringstream PaidCapiFactorDesc2;
	PaidCapiFactorDesc2 << "IF(" << VanillaMepiColNamesTable[Cash2] << "[i]<0," << VanillaMepiColNamesTable[CapiFactor] << "[i]/(1+" << VanillaMepiConstantsTable[LeverageCost];
	PaidCapiFactorDesc2 << "*" << VanillaMepiColNamesTable[CapiFactor] << "[i]/12)," << VanillaMepiColNamesTable[CapiFactor] << "[i]/(1-" << VanillaMepiConstantsTable[CashSpread];
	PaidCapiFactorDesc2 << "*" << "CapiFactor" << "[i]/12))";
	rowDescVec[PaidCapiFactor2] = PaidCapiFactorDesc2.str();
	rowTypeVec[PaidCapiFactor2] = ARM_STRING;

	// Cash2
	CC_Ostringstream CashDesc2;
	(eventIdx == 1) ? (CashDesc2 << VanillaMepiConstantsTable[StartingCash] ) : (CashDesc2 << "(" << VanillaMepiColNamesTable[Portfolio2] << "[i-1]-" << VanillaMepiColNamesTable[ActifRisque2] << "[i-1])/" << VanillaMepiColNamesTable[PaidCapiFactor2] << "[i-1]" );
	rowDescVec[Cash2] = CashDesc2.str();
	rowTypeVec[Cash2] = ARM_STRING;

    // Portfolio2
	CC_Ostringstream PortfolioDesc2;
	(eventIdx == 1) ? ( PortfolioDesc2 << VanillaMepiConstantsTable[StartingPortfolio] << "-" << VanillaMepiConstantsTable[DeltaShift] << "/2" ) : ( PortfolioDesc2 << VanillaMepiColNamesTable[ActifRisque2] << "[i-1]*" << VanillaMepiColNamesTable[Underlying] << "[i]/" << VanillaMepiColNamesTable[Underlying] << "[i-1]+" << VanillaMepiColNamesTable[Cash2] << "[i]-Fees/12"); 
	rowDescVec[Portfolio2] = PortfolioDesc2.str();
	rowTypeVec[Portfolio2] = ARM_STRING;

	// AvgPortfolio2
	// The value associated with AlreadyAsianed should be equal to 0 , when not in the asianing period
	CC_Ostringstream AvgPortfolioDesc2;
	if(eventIdx == nbEvents-1)
	{
		size_t asianingNb = MIN( itsAsianingPeriodNb , nbEvents-1 );
		AvgPortfolioDesc2 << "(" << VanillaMepiConstantsTable[AlreadyAsianed] << "+"<< VanillaMepiColNamesTable[Portfolio2] << "[i]";
		for( long i=1 ; i<asianingNb;++i)
			AvgPortfolioDesc2 << "+" << VanillaMepiColNamesTable[Portfolio2] << "[i-" << i << "]";

		AvgPortfolioDesc2 << ")/" << VanillaMepiConstantsTable[AsianDatesNb];

		rowTypeVec[AvgPortfolio2] = ARM_STRING;
	}
	else
	{
		AvgPortfolioDesc2 << "0";
		rowTypeVec[AvgPortfolio2] = ARM_DOUBLE;
	}
	rowDescVec[AvgPortfolio2] = AvgPortfolioDesc2.str();

	
	/// ActifRisque2
	CC_Ostringstream ActifRisqueDesc2;
	ActifRisqueDesc2 << "Min(Max(" <<VanillaMepiConstantsTable[MinInvested] << "*" << VanillaMepiColNamesTable[Portfolio2] << "[i]," << VanillaMepiConstantsTable[RiskFactor] << "*(" << VanillaMepiColNamesTable[Portfolio2] << "[i]-" << VanillaMepiColNamesTable[FixedProtection] << "[i])),";
	ActifRisqueDesc2 << VanillaMepiColNamesTable[Portfolio2] << "[i]+" <<VanillaMepiConstantsTable[MaxBorrow] << ")";
	rowDescVec[ActifRisque2] = ActifRisqueDesc2.str();
	rowTypeVec[ActifRisque2] = ARM_STRING;


	// CF
	CC_Ostringstream CFDesc;
	if( eventIdx == nbEvents-1 ) 
	{
		CFDesc << "(Max(" << VanillaMepiColNamesTable[AvgPortfolio1] << "[i]-" << VanillaMepiConstantsTable[Strike] << ",0)-Max(" << VanillaMepiColNamesTable[AvgPortfolio2] << "[i]-" << VanillaMepiConstantsTable[Strike] << ",0) )/" << VanillaMepiConstantsTable[DeltaShift];
		rowTypeVec[CF] = ARM_STRING;
	}
	else
	{
		CFDesc << "0";
		rowTypeVec[CF] = ARM_DOUBLE;
	}
	rowDescVec[CF] = CFDesc.str();

	return ARM_RowInfo(rowDescVec,rowTypeVec);
}



ARM_GP_VectorPtr ARM_VanillaMepiDelta::DatesStructure() const
{
	return ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>( itsEventDates->Clone() ) );
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaArgNumeric
///	Routine: securityCstManager
///	Returns: 
///	Action : gives the constant manager of the security
////////////////////////////////////////////////////
const ARM_CstManagerPtr ARM_VanillaMepiDelta::securityCstManager() const
{
	size_t i;
	size_t cstManagerSize = sizeof(VanillaMepiConstantsTable)/sizeof(VanillaMepiConstantsTable[0]);
	vector<string> names(cstManagerSize);
	vector<double> values(cstManagerSize);

	for( i=0 ; i < cstManagerSize ; ++i )
		names[i] = VanillaMepiConstantsTable[i];

	values[RiskFactor] = itsRiskFactor;
	values[Strike] = itsStrike;
	values[MaxBorrow] = itsMaxBorrow;
	values[StartingPortfolio] = itsStartingPortfolio;
	values[StartingCash] = itsStartingCash;
	values[MinInvested] = itsMinInvested;
	values[LeverageCost] = itsLeverageCost;
	values[CashSpread] = itsCashSpread;
	values[Fees] = itsFees;
	values[AlreadyAsianed] = itsAlreadyAsianed;
	values[AsianDatesNb] = (double) itsAsianingPeriodNb;
	values[DeltaShift] = itsDeltaShift;

// FIXMEFRED: mig.vc8 (22/05/2007 18:11:41):cast
	return static_cast<ARM_CstManagerPtr>(new ARM_CstManager(names,values));
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaArgNumeric
///	Routine: AnalyticPrice
///	Returns: 
///	Action : analytic price of spreadoption with a model checking 
////////////////////////////////////////////////////
double ARM_VanillaMepiDelta::AnalyticPrice(ARM_PricingModel* model) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : No Analytic Price for MEPI " );

}
////////////////////////////////////////////////////
///	Struct : ARM_VanillaMepiDelta
///	Routine: toString
///	Returns: string
///	Action : stringify the object to give details about it
////////////////////////////////////////////////////
string ARM_VanillaMepiDelta::toString(const string& indent, const string& nextIndent) const
{ 
	return GetGenSecurity()->toString();
}

double ARM_VanillaMepiDelta::Price(ARM_PricingModel* model ) const
{
	ARM_GenPricer genPricer( &*GetGenSecurity(), model );
	return genPricer.Price();
}







CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


