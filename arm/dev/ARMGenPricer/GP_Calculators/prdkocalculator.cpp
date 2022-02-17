/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file prdccalculator.cpp
 *
 *  \brief file for the Power Reverse Dual Currencies Calculator
 *	\author  E.Ezzine
 *	\version 1.0
 *	\date February 2006
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/prdkocalculator.h"
#include "gpcalculators/argconvdefault.h"

/// gpbase
#include "gpbase/datestripcombiner.h"
#include "gpbase/globalconstant.h"


/// gpinfra
#include "gpinfra/dealdescription.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/cstmanager.h"
#include "gpinfra/pricingadviser.h"
#include "gpcalib/calibmethod.h"

/// Kernel
#include <inst/option.h>
#include <inst/portfolio.h>
#include <inst/forex.h>
#include <mod/xbsfx.h>

/// STL
#include <iomanip> /// for setprecision()
#include <list>
CC_USING_NS(std,list)


CC_BEGIN_NAMESPACE( ARM )

const double FX_DEFAULT_WEIGHT      = 1.0;
const double FX_DEFAULT_PRICE       = 1.0e+30;


const string ARM_PRDKOCalculator::PRDKOColNamesTable [] =
{
	"Barrier",
    "FwdFX",
    "PRDKO",
    "PRDKOption"
};
/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDKOCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_PRDKOCalculator::ARM_PRDKOCalculator( const ARM_PRDKOCalculator& rhs )
:	ARM_PRDCalculator( rhs ),
	
	itsvBarrier(rhs.itsvBarrier),
	itsNbCallable(rhs.itsNbCallable),
	itsvCpnIsCallable(rhs.itsvCpnIsCallable)
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDKOCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////
ARM_PRDKOCalculator::ARM_PRDKOCalculator(const ARM_Date& asofDate,
		const ARM_Date& startDate,
		const ARM_Date& fixEndDate,
		const ARM_Date& endDate,
		const ARM_Currency& domCcy,
		const ARM_Currency& forCcy,
		const ARM_Currency& fundCcy,
		int cpnDayCount,
		int cpnFreq,
		int FxResetGap,
		int stubRule,
		int resetTiming,
		const string& cpnResetCal,
		const string& cpnPayCal,
		const ARM_Curve& cpnnominal,
		const ARM_Curve& domesticCpn,
		const ARM_Curve& foreignCpn,
		const ARM_Curve& minCpn,
		const ARM_Curve& maxCpn,
		const ARM_Curve& initialFx,
		const ARM_Curve& barrier,
		const ARM_Curve& noticeType,
		int fundFreq,
		int fundDayCount,
		int compFreq,
		int compType,
		const ARM_Curve& fundnominal,
		const ARM_Curve& fundSpread,
		int exerciseFreq,
		int noticeGap,
		int payRec,
		size_t nbNCall,
		const ARM_Curve& fees,
		int redemptionGap ,
		double redemptionStrike,
		const ARM_RedemptionType& redemptionType,		
		const ARM_StringVector& columnsToPrice,
		bool fxLocalModelFlag,
		bool basisIRCalibFlag,
		const ARM_Curve& fundlevrage)
:	
	ARM_PRDCalculator( asofDate,
		startDate,
		fixEndDate,
		endDate,
		domCcy,
		forCcy,
		fundCcy,
		cpnDayCount,
		cpnFreq,
		FxResetGap,
		stubRule,
		resetTiming,
		cpnResetCal,
		cpnPayCal,
		cpnnominal,
		domesticCpn,
		foreignCpn,
		minCpn,
		maxCpn,
		initialFx,
		fundFreq,
		fundDayCount,
		compFreq,
		compType,
		fundnominal,
		fundSpread,
		exerciseFreq,
		noticeGap,
		payRec,
		nbNCall,
		fees,
		redemptionGap ,
		redemptionStrike,
		redemptionType,
		columnsToPrice,
		fxLocalModelFlag,
		basisIRCalibFlag,
		fundlevrage),
	itsvBarrier(0),	
	itsNbCallable(0),
	itsvCpnIsCallable(0)

{
	/// convert barrier from curve to vector
	ComputeProductVectorsFromCurves(barrier,noticeType);

	/// Build an internal cst manager to save
    /// Fx options date strip + strikes, nominals & leverage profiles
	 ARM_CstManagerPtr fxDataManager(NULL);
     fxDataManager = ComputeCstManager();

    /// Create the Generic Security paid in domestic=coupon currency
    CreateAndSetDealDescriptionAndTimeIt(GetKeys()[YcBasisDomKey],itsColumnsToPrice,fxDataManager);

    /// To boost time, disable intermediatePayoffs & snapshots computation
    GetGenSecurity()->SetOtherPayoffsFlag(false);
	
}
	/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: ComputeProductVectorsFromCurves
///	Returns: void
///	Action : compute vector attributes from ARM_Curve 
/////////////////////////////////////////////////////////////////
void ARM_PRDKOCalculator::ComputeProductVectorsFromCurves(const ARM_Curve& barrier,
														  const ARM_Curve& noticeType)
{
	size_t i;
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();	
	
	/// Convert ARM_Curves into relevant vectors
	/// this will be much easier for computations
	ARM_GP_Vector* cpnResetDates = itsStructDateStrip->GetResetDates() ;
	double lag;
	itsCpnSize = cpnResetDates->size();
	itsvCpnIsCallable = ARM_BoolVector(itsCpnSize,false);
	for (i= 0; i<itsCpnSize; i++)
	{
		lag = (*cpnResetDates)[i]-asOfDate;
		itsvBarrier.push_back(barrier.Interpolate(lag) );
		ARM_PRDNoticeType type = (ARM_PRDNoticeType)(int)noticeType.Interpolate(lag);
		if(type == ARM_MixPRDNoticeType::Call)
			itsvCpnIsCallable[i] = true;
	}

	itsvCpnIsCallable = ARM_BoolVector(itsExerSize,false);
	for (i= 0; i<itsExerSize; i++)
	{
		lag = (*cpnResetDates)[itsvCpnIndex[i]]-asOfDate;
		ARM_PRDNoticeType type = (ARM_PRDNoticeType)(int)noticeType.Interpolate(lag);
		if(type == ARM_MixPRDNoticeType::Call)
			itsvCpnIsCallable[i] = true;
	}


}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CheckData
///	Returns: void
///	Action : check if PRDC datas are consistent
/////////////////////////////////////////////////////////////////
void ARM_PRDKOCalculator::CheckData()
{
    /// Check if input columns to price are present in the deal description
    size_t colNamesSize = sizeof(PRDKOColNamesTable)/sizeof(PRDKOColNamesTable[0]);
    for(size_t i=0;i<itsColumnsToPrice.size();++i)
    {
        for(size_t j=0;j<colNamesSize;++j)
            if(itsColumnsToPrice[i] == PRDKOColNamesTable[j])
                break;
        if(j==itsColumnsToPrice.size())
			ARM_PRDCalculator::CheckData();
    }
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDKOCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_PRDKOCalculator::ColumnNames() const
{
	ARM_RowInfo  rowInfo = ARM_PRDCalculator::ColumnNames();

    size_t colNamesSize = sizeof(PRDKOColNamesTable)/sizeof(PRDKOColNamesTable[0]);
    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0;i<colNamesSize; ++i){
        rowInfo.first.push_back(PRDKOColNamesTable[i]);
		rowInfo.second.push_back(ARM_STRING);
	}

    return rowInfo;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDKOCalculator
///	Routine: InitPriceableColumns
///	Returns: void
///	Action : initialise to 0 column to be able to sum
/////////////////////////////////////////////////////////////////
void ARM_PRDKOCalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
	string zeroValue("0");
    rowDescVec[PRDKOption] = zeroValue;
    rowTypeVec[PRDKOption] = ARM_DOUBLE;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDKOCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_PRDKOCalculator::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{
	ARM_RowInfo  rowInfo = ARM_PRDCalculator::MiddleRows(eventIdx,datesStructure );

	size_t eventSize = datesStructure.GetDateStrip(0)->GetResetDates()->size();
	size_t descSize = sizeof(PRDKOColNamesTable)/sizeof(PRDKOColNamesTable[0]);

	vector< string > rowDescVec(descSize);
	vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 

	/// Set default 0 value for each column to be able to sum it
	InitPriceableColumns(rowDescVec,rowTypeVec);

	/// Get the model names for domestic and forex models
	string basisDomModelName = GetKeys()[YcBasisDomKey];
	string fxModelName = GetKeys()[ForexKey];

	string nextExerIdx("[i+1]");
	bool isLastEvent = (eventIdx+1>=eventSize);
	
	/// Barrier
	CC_Ostringstream barrierDesc;
	double barrier = itsvBarrier[itsvCpnIndex[eventIdx]];
	barrierDesc << CC_NS(std,fixed) << barrier;
	rowDescVec[Barrier] = barrierDesc.str();
	rowTypeVec[Barrier] = ARM_DOUBLE;

	/// FwdFx
	CC_Ostringstream fwdfxDesc;
	fwdfxDesc << "FWD("<< fxModelName << "," << PRDColNamesTable[EventDate] << "[i])";
	rowDescVec[FwdFX] = fwdfxDesc.str();
	rowTypeVec[FwdFX] = ARM_STRING;

	string zeroValue("0.0");
	string FxModelName(fxModelName);
	CC_Ostringstream prdkoDesc;
	if(itsvCpnIsCallable[itsvCpnIndex[eventIdx]])
	{
		if(itsvIsExerDate[eventIdx]){
			prdkoDesc << "Exercise(" << zeroValue << "," << PRDColNamesTable[PRDSwap] << "[i],";
			isLastEvent ? prdkoDesc << zeroValue << ")" : prdkoDesc << "" << PRDKOColNamesTable[PRDKO] << "[i+1])" ;
		}
		else{
			isLastEvent ? prdkoDesc << zeroValue  : prdkoDesc << "PV(" << PRDKOColNamesTable[PRDKO] << "[i+1])" ;
		}
			rowDescVec[PRDKO] = prdkoDesc.str();
			rowTypeVec[PRDKO] = ARM_STRING;

	}
	else {
		if(itsvIsExerDate[eventIdx]){
			prdkoDesc << "if("<< PRDKOColNamesTable[FwdFX] << "[i]>=" << PRDKOColNamesTable[Barrier] << "[i],";
			prdkoDesc << PRDColNamesTable[PRDSwap] << "[i],";
			isLastEvent ? prdkoDesc << zeroValue << ")" : prdkoDesc << "PV(" << PRDKOColNamesTable[PRDKO] << "[i+1]))" ;
		}
		else{
			isLastEvent ? prdkoDesc << zeroValue  :  prdkoDesc << "PV(" << PRDKOColNamesTable[PRDKO] << "[i+1])" ;
		}
		rowDescVec[PRDKO] = prdkoDesc.str();
		rowTypeVec[PRDKO] = ARM_STRING;
	}

		bool isFirstEvent = (eventIdx==GetNbNoCall());
		CC_Ostringstream prdcOptionDesc;
		isFirstEvent ? prdcOptionDesc << PRDKOColNamesTable[PRDKO] << "[i]" : prdcOptionDesc << zeroValue ;
		rowDescVec[PRDKOption] = prdcOptionDesc.str();
		rowTypeVec[PRDKOption] = ARM_STRING;

	for(size_t i=0;i<descSize; ++i){
		rowInfo.first.push_back(rowDescVec[i]);
		rowInfo.second.push_back(rowTypeVec[i]);
	}

	return rowInfo;
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDKOCalculator
///	Routine: GetStrikeAtATSFx
///	Returns: double
///	Action : compute market target prices of the FX option
/////////////////////////////////////////////////////////////////
void ARM_PRDKOCalculator::UpdateFxOption(size_t eventIdx,ARM_Option*& fxOption)
{
	/// strike = moneyness * fwdFx
	if( !itsvCpnIsCallable[itsvCpnIndex[eventIdx]]){
		double barrier = itsvBarrier[itsvCpnIndex[eventIdx]];
		fxOption->SetStrike(barrier);
	}
}
/*/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDKOCalculator
///	Routine: CreateFxOption
///	Returns: a portfolio
///	Action : create the list of ATM FX options
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_PRDKOCalculator::CreateFxOption()
{
    double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

    //ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );

	/// To support new FX GP models !
    ARM_DFBSModel* fxBSModel = dynamic_cast< ARM_DFBSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );
    
    ARM_Option* fxOption;
    list < ARM_Security* > fxOptionList;

    const ARM_DealDescription& dealDesc = GetGenSecurity()->GetDealDescription();
    size_t i,nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names
    size_t nbFutureExer=itsvIsExerDate.size();
    if(nbFutureExer+1 != nbEvents)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : mismatch between exercise flags & notice dates in the deal description");

    /// Restore forex
	ARM_Forex* forex = static_cast< ARM_Forex* >( GetMktDataManager()->GetData(GetKeys()[ForexKey]) );
	
	ARM_GP_Vector* exerResetDates = itsExerciseDateStrip->GetResetDates();
    if(itsCalibType == ARM_PRCSCalibTypes::ATMCalib)
    {
        /// Portfolio = market standard FX options (ATM)

		/// Compute standard vol expiries
		ARM_BSModel* bsModel = fxBSModel->CreateFxBSModel();
		ARM_VolCurve* fxBSVol = bsModel->GetVolatility();
		size_t nbExp = fxBSVol->GetExpiryTerms()->GetSize();
		ARM_GP_Vector fxStdExp(nbExp);
		for(i=0;i<nbExp;++i)
			fxStdExp[i] = asOfDate + K_YEAR_LEN * (*(fxBSVol->GetExpiryTerms()))[i];
		delete bsModel;

        /// Check if the 1st notice date may replace the 1st standard fx option expiry
        ARM_Date firstNoticeDate;
        size_t fxExpIdx=0;
        for(size_t exerIdx=0; exerIdx<itsExerSize; ++exerIdx)
        {
            if(itsvIsExerDate[exerIdx] && !itsvCpnIsFixed [exerIdx])
            {
                firstNoticeDate = ARM_Date( (*itsExerciseDateStrip->GetResetDates())[exerIdx]);
                if(firstNoticeDate.GetJulian() >= fxStdExp[fxExpIdx]-K_NEW_DOUBLE_TOL)
                {
                    /// Replace the 1st standard expiry by the 1st notice
                    fxOption    = new ARM_Option(forex,firstNoticeDate,K_MARKET_RATE,K_CALL,K_EUROPEAN);
                    fxOptionList.push_back(static_cast< ARM_Security* >(fxOption));
                    ++fxExpIdx;
                    break;
                }
            }
        }

        /// Create a fx option for each following standard expiry greater than the 1st notice
        /// but lower or equal to 30y
        double maxFxOptionExpiry = asOfDate + 11111; /// 11111/365 is > 30y!!
        for(;fxExpIdx < nbExp && fxStdExp[fxExpIdx] < maxFxOptionExpiry;++fxExpIdx)
        {
            if(fxStdExp[fxExpIdx] > firstNoticeDate.GetJulian() + K_NEW_DOUBLE_TOL)
            {
                fxOption    = new ARM_Option(forex,ARM_Date(fxStdExp[fxExpIdx]),K_MARKET_RATE,K_CALL,K_EUROPEAN);
                fxOptionList.push_back(static_cast< ARM_Security* >(fxOption));
            }
        }
    }
    else
    {
        /// Portfolio = FX options of the PRD leg, calibration strike depends on calibType
        ARM_Date fxResetDate,fxPayDate;
        double fxStrike;
		bool latestIsCallable = false;
        size_t strikeIdx=0,nbStrikes = itsCalibDatas.size();
        for(size_t exerIdx=0; exerIdx<itsExerSize; ++exerIdx)
        {
            if(itsvIsExerDate[exerIdx])
            {
				if(itsvCpnIsCallable[exerIdx] && !itsvCpnIsFixed [itsvCpnIndex[exerIdx]]){
					/// Here is a Fx option for a floored coupon
					fxResetDate = ARM_Date((*(itsForexDateStrip->GetResetDates()))[itsvCpnIndex[exerIdx]]);
					fxPayDate   =ARM_Date((*(itsStructDateStrip->GetPaymentDates()))[itsvCpnIndex[exerIdx]]);
					latestIsCallable = true;
				}
				else if(itsvCpnIsCallable[exerIdx] && itsvCpnIsFixed [itsvCpnIndex[exerIdx]]){
					continue;
				}
				else if(latestIsCallable){
					latestIsCallable = false;
					continue;
				}
				else{
					fxResetDate = ARM_Date((*(itsExerciseDateStrip->GetResetDates()))[exerIdx]);
				}

                if(itsCalibType == ARM_PRCSCalibTypes::ATSFxEquivCalib){
                    fxStrike = ComputeFxEquivStrike(exerIdx+1);
                }
                else if(itsCalibType == ARM_PRCSCalibTypes::ATSFxProfileCalib){
                    fxStrike = itsCalibDatas[strikeIdx];
                    if(strikeIdx+1 < nbStrikes)
                        ++strikeIdx;
                }
                else{
					double lowfxStrike    = itsvLowStrikeCpn[itsvCpnIndex[exerIdx]];
					double upfxStrike     = itsvHighStrikeCpn[itsvCpnIndex[exerIdx]];
					
					/// Here is a Fx option for a capped coupon
                    if(itsvCpnIsCapped[itsvCpnIndex[exerIdx]])
                        fxStrike = 0.5*( lowfxStrike + upfxStrike);
                }

                fxOption = new ARM_Option(forex,fxResetDate,fxStrike,K_CALL,K_EUROPEAN);
				if(itsvCpnIsCallable[exerIdx])
					fxOption->SetPayDate(fxPayDate);

                fxOptionList.push_back(static_cast< ARM_Security* >(fxOption));
            }
        }
    }


    ARM_StdPortfolio* fxPort = new ARM_StdPortfolio(fxOptionList);
    for(i=0;i<fxPort->size();++i)
    {
        fxPort->SetWeight(FX_DEFAULT_WEIGHT,i);
        fxPort->SetPrice((i+1)*FX_DEFAULT_PRICE,i);
    }

    return ARM_StdPortfolioPtr(fxPort);
}

*/
////////////////////////////////////////////////////
///	Class   : ARM_PRDCCalculator
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_PRDKOCalculator::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream prdcData;


    /// PRDC Calculator specific datas dump
    prdcData << indent <<"\n\n =======> POWER REVERSE DUAL CURRENCIES CALCULATOR <====== \n";
	prdcData << indent <<"\t - Knouk out/in\n\n";

    prdcData << indent << "Columns to Price :\n";
    for(size_t i=0;i<itsColumnsToPrice.size();++i)
    {
            prdcData << indent << itsColumnsToPrice[i] << "\n";
    }
/// GenCalculator general datas dump
	prdcData << ARM_PRDCalculator::toString(indent,nextIndent) << "\n\n";
    return prdcData.str();
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

