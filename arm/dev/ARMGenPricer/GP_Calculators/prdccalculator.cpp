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
#include "gpcalculators/prdccalculator.h"
#include "gpcalculators/basisconverter.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/autocleaner.h"
#include "gpbase/ostringstream.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/datestrip.h"
#include "gpbase/singleton.h"
#include "gpbase/utilityport.h"  
#include "gpbase/curveconvert.h"  
#include "gpbase/gplinalgconvert.h"
#include "gpbase/datestripconvert.h"
#include "gpbase/globalconstant.h"

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
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/cstmanager.h"
#include "gpinfra/discretisationscheme.h"

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/modelfitter.h"
#include "gpcalib/VanillaSwaption.h"
#include "gpcalib/VanillaFxOption.h"
#include "gpcalib/VanillaIrFxSwaption.h"
#include "gpcalib/KernelToGP.h"
#include "gpcalib/stripper.h"


/// gpmodels
#include "gpmodels/2irfxModel.h"
#include "gpmodels/q1f.h"
#include "gpmodels/q1f_fx.h"
#include "gpmodels/modelparamsq1f.h"
#include "gpmodels/eqfx_modelfactory.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/forwardmarginbasis.h"
#include "gpmodels/local_sln_model.h"
#include "gpmodels/local_sln_modelparams.h"
#include "gpmodels/MarketIRModel.h"
#include "gpmodels/MarketHybridModel.h"

/// gpnummethods
#include "gpnummethods/treebase.h"
#include "gpnummethods/treefactory.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/reconnector.h"
#include "gpnummethods/smoother.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"

/// kernel
#include <crv/volflat.h>
#include <inst/powrev.h>
#include <mod/xbsfx.h>
#include <inst/forex.h>
#include <inst/swaption.h>
#include <inst/option.h>
#include <inst/portfolio.h>
#include <inst/fixleg.h>
#include <inst/swapleg.h>
#include <inst/swap.h>
#include <util/fromto.h>
#include <inst/optionportfolio.h>


/// STL
#include <iomanip> /// for setprecision()
#include <list>
CC_USING_NS(std,list)


CC_BEGIN_NAMESPACE( ARM )

const string ARM_PRDCCalculator::PRDCColNamesTable [] =
{
    "PRDCBermuda",
    "PRDCOption"
};



/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_PRDCCalculator::ARM_PRDCCalculator( const ARM_PRDCCalculator& rhs )
:	ARM_PRDCalculator( rhs )
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////
ARM_PRDCCalculator::ARM_PRDCCalculator(const ARM_Date& asofDate,
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
	ARM_PRDCalculator(asofDate,
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
		fundlevrage)
{

	
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
///	Class  : ARM_PRDCCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////
ARM_PRDCCalculator::ARM_PRDCCalculator(
                ARM_PowerReverse* powRev,
                ARM_DFBSModel* model,
                const ARM_ObjectVector& otherMktDatas,
                const ARM_GP_Vector& schedulerDatas,
                const ARM_GP_Vector& truncatorDatas,
                const ARM_StringVector& columnsToPrice,
                bool markovianDriftSamplerFlag,
                bool fxLocalModelFlag,
                ARM_PRDCCalibType calibType,
                const ARM_GP_Vector& calibDatas,
                bool basisIRCalibFlag,
				ARM_BasisType basisType)
:	
	ARM_PRDCalculator(powRev,
                 model,
                 otherMktDatas,
                 schedulerDatas,
                 truncatorDatas,
                 columnsToPrice,
                 markovianDriftSamplerFlag,
                fxLocalModelFlag,
                calibType,
                calibDatas,
                basisIRCalibFlag,
				basisType)
{
	 ARM_CstManagerPtr fxDataManager(NULL);
     fxDataManager = fxDataManager=ComputeCstManager();

    /// Create the Generic Security paid in domestic=coupon currency
    CreateAndSetDealDescriptionAndTimeIt(GetKeys()[YcBasisDomKey],itsColumnsToPrice,fxDataManager);

    /// To boost time, disable intermediatePayoffs & snapshots computation
    GetGenSecurity()->SetOtherPayoffsFlag(false);

    /// Create calibration sets
    CreateAndSetCalibrationAndTimeIt();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_PRDCCalculator::ColumnNames() const
{
	ARM_RowInfo  rowInfo = ARM_PRDCalculator::ColumnNames();

    size_t colNamesSize = sizeof(PRDCColNamesTable)/sizeof(PRDCColNamesTable[0]);
    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0;i<colNamesSize; ++i){
        rowInfo.first.push_back(PRDCColNamesTable[i]);
		rowInfo.second.push_back(ARM_STRING);
	}

    return rowInfo;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CheckData
///	Returns: void
///	Action : check if PRDC datas are consistent
/////////////////////////////////////////////////////////////////
void ARM_PRDCCalculator::CheckData()
{
    /// Check if input columns to price are present in the deal description
    size_t colNamesSize = sizeof(PRDCColNamesTable)/sizeof(PRDCColNamesTable[0]);
    for(size_t i=0;i<itsColumnsToPrice.size();++i)
    {
        for(size_t j=0;j<colNamesSize;++j)
            if(itsColumnsToPrice[i] == PRDCColNamesTable[j])
                break;
        if(j==itsColumnsToPrice.size())
			ARM_PRDCalculator::CheckData();
    }
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCCalculator
///	Routine: InitPriceableColumns
///	Returns: void
///	Action : initialise to 0 column to be able to sum
/////////////////////////////////////////////////////////////////
void ARM_PRDCCalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
    string zeroValue("0");
    rowDescVec[PRDCOption] = zeroValue;
    rowTypeVec[PRDCOption] = ARM_DOUBLE;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_PRDCCalculator::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{
	ARM_RowInfo  rowInfo = ARM_PRDCalculator::MiddleRows(eventIdx,datesStructure );

	size_t eventSize = datesStructure.GetDateStrip(0)->GetResetDates()->size();
	size_t descSize = sizeof(PRDCColNamesTable)/sizeof(PRDCColNamesTable[0]);

	vector< string > rowDescVec(descSize);
	vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 

	/// Set default 0 value for each column to be able to sum it
	InitPriceableColumns(rowDescVec,rowTypeVec);

	/// Get the model names for domestic and forex models
	string basisDomModelName = GetKeys()[YcBasisDomKey];
	string fxModelName = GetKeys()[ForexKey];

	string nextExerIdx("[i+1]");
	bool isLastEvent = (eventIdx+1>=eventSize);


	/// PRDC option description
	CC_Ostringstream prdcBermudaDesc;
	if(itsvIsExerDate[eventIdx])
	{
		/// Actual exercise at this event date
		prdcBermudaDesc << "EXERCISE(0," << PRDColNamesTable[PRDSwap] << "[i],";
		if(isLastEvent)
			/// No residual option
			prdcBermudaDesc << "0)";
		else
			/// Only one description to get the bermuda price
			prdcBermudaDesc << PRDCColNamesTable[PRDCBermuda] << nextExerIdx << ")";
	}
	else
	{
		/// No exercise allowed at this event date
		if(isLastEvent)
			/// No residual option
			prdcBermudaDesc << "0";
		else
			/// Just actualise residual option
			prdcBermudaDesc << "PV(" << PRDCColNamesTable[PRDCBermuda] << nextExerIdx << ")";
	}

	rowDescVec[PRDCBermuda] = prdcBermudaDesc.str();
	rowTypeVec[PRDCBermuda] = ARM_STRING;

	/// The single PRDC option a first line to get the price
	bool isFirstEvent = (eventIdx==itsNbNoCall);
	if(isFirstEvent)
	{
		CC_Ostringstream prdcOptionDesc;
		prdcOptionDesc << PRDCColNamesTable[PRDCBermuda] << "[i]";
		rowDescVec[PRDCOption] = prdcOptionDesc.str();
		rowTypeVec[PRDCOption] = ARM_STRING;
	}

	for(size_t i=0;i<descSize; ++i){
		rowInfo.first.push_back(rowDescVec[i]);
		rowInfo.second.push_back(rowTypeVec[i]);
	}

	return rowInfo;
}

////////////////////////////////////////////////////
///	Class   : ARM_PRDCCalculator
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_PRDCCalculator::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream prdcData;


    /// PRDC Calculator specific datas dump
    prdcData << indent <<"\n\n =======> POWER REVERSE DUAL CURRENCIES CALCULATOR <====== \n";
	prdcData << indent <<"\t - Callable\n\n";

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

