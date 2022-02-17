/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gencalculator.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou, J-M Prié
 *	\version 1.0
 *	\date March 2004
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/gencalculator.h"
#include "gpcalculators/basisconverter.h"

/// gpbase
#include "gpbase/datestripcombiner.h"
#include "gpbase/datestrip.h"
#include "gpbase/singleton.h"
#include "gpbase/valuetype.h"
#include <glob/expt.h>
#include "gpbase/timer.h"
#include "gpbase/warning.h"

/// gpinfra
#include "gpinfra/mktdatamanagerrep.h"
#include "gpinfra/mktdatamanager.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/modelnrefcall.h"
#include "gpinfra/gramnode.h"
#include "gpinfra/pricingmodel.h"

/// gpcalib
#include "gpcalib/calibmethod.h"


CC_BEGIN_NAMESPACE( ARM )



////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: CopyNoCleanUp
///	Returns: void
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_GenCalculator::CopyNoCleanUp(const ARM_GenCalculator& rhs)
{
    itsAccountingPricingFlag = rhs.itsAccountingPricingFlag;    
	itsLiborFxFixings = rhs.itsLiborFxFixings ? (ARM_FixingSched *) rhs.itsLiborFxFixings->Clone() : NULL;
    itsCurrency				= rhs.itsCurrency? (ARM_Currency*) rhs.itsCurrency->Clone(): NULL;
    itsFundingCcy           = rhs.itsFundingCcy;
    itsBasisCcy             = rhs.itsBasisCcy;
    itsDomesticCcy			= rhs.itsDomesticCcy;
    itsForeignCcy			= rhs.itsForeignCcy;

    itsGenSecurity          = ARM_GenSecurityPtr((rhs.itsGenSecurity != ARM_GenSecurityPtr(NULL)) ? (ARM_GenSecurity*) rhs.itsGenSecurity->Clone() : NULL);
    itsCalibMethod	        = ARM_CalibMethodPtr((rhs.itsCalibMethod != ARM_CalibMethodPtr(NULL)) ? (ARM_CalibMethod*) rhs.itsCalibMethod->Clone() : NULL);
    itsPricingModel	        = ARM_PricingModelPtr((rhs.itsPricingModel != ARM_PricingModelPtr(NULL)) ? (ARM_PricingModel*) rhs.itsPricingModel->Clone() : NULL);

    /// Clone to maintain keys consistency
    itsMktDataManager       =   rhs.itsMktDataManager != ARM_MarketData_ManagerRepPtr(NULL)
                                ? ARM_MarketData_ManagerRepPtr( (ARM_MarketData_ManagerRep*) const_cast<ARM_GenCalculator&>(rhs).itsMktDataManager->Clone() )
                                : rhs.itsMktDataManager;
    itsMDMKeys              = rhs.itsMDMKeys;
	itsPricingData			= rhs.itsPricingData;
	itsFullTerminated		= rhs.itsFullTerminated;
	itsUpdate				= rhs.itsUpdate;
	itsPorS                 = rhs.itsPorS;
	itsWarning				= rhs.itsWarning? static_cast<ARM_Warning*>( rhs.itsWarning->Clone() ) : NULL;
	itsDegeneratedCalculator = rhs.itsDegeneratedCalculator;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: default constructor
///	Returns: void
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_GenCalculator::ARM_GenCalculator(const ARM_Date& asOfDate)
:	ARM_RootObject(),
    itsAccountingPricingFlag(ARM_DISC_PRICING_METH), 
    itsLiborFxFixings(NULL),
	itsGenSecurity(NULL),
	itsCalibMethod(NULL),
	itsPricingModel(NULL),
	itsMktDataManager(NULL),
	itsCurrency( new ARM_Currency ),
	itsFundingCcy(),
	itsDomesticCcy(),
	itsForeignCcy(),
	itsPricingData(),
	itsFullTerminated(false),
	itsUpdate(false),
	itsPorS(1),
	itsWarning(NULL),
	itsDegeneratedCalculator(No)
{
    itsMktDataManager = ARM_MarketData_ManagerRepPtr(new ARM_MarketData_ManagerRep(asOfDate));
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_GenCalculator::ARM_GenCalculator( const ARM_MarketData_ManagerRep& mktDataManager )
:	ARM_RootObject(),
    itsAccountingPricingFlag(ARM_DISC_PRICING_METH),
    itsLiborFxFixings(NULL),
	itsGenSecurity(NULL),
	itsCalibMethod(NULL),
	itsPricingModel(NULL),
	itsMktDataManager(NULL),
	itsCurrency( new ARM_Currency ),
	itsFullTerminated(false),
	itsUpdate(false),
	itsPorS(1),
	itsWarning(NULL),
	itsDegeneratedCalculator(No)
{
    itsMktDataManager = ARM_MarketData_ManagerRepPtr( (ARM_MarketData_ManagerRep*) const_cast<ARM_MarketData_ManagerRep&>(mktDataManager).Clone() );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_GenCalculator::ARM_GenCalculator( const ARM_GenCalculator& rhs )
:	ARM_RootObject( rhs )
{
    CopyNoCleanUp(rhs);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: assignment operator
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_GenCalculator& ARM_GenCalculator::operator=( const ARM_GenCalculator& rhs )
{
	if( this != & rhs )
	{
		ARM_RootObject::operator=(rhs);
		delete itsCurrency;
		delete itsWarning;
		itsWarning = 0;
        CopyNoCleanUp(rhs);
	}
	return *this;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: Destructor
///	Returns: void
///	Action : destroys the object
/////////////////////////////////////////////////////////////////
ARM_GenCalculator::~ARM_GenCalculator()
{
    delete itsLiborFxFixings;
	itsLiborFxFixings = NULL;
	delete itsCurrency;
	delete itsWarning;
	itsWarning = 0;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: ValidateGenSecurity
///	Returns: void
///	Action : Check if the security is able to replace the built-in one
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::ValidateGenSecurity(const ARM_GenSecurityPtr& genSecurity)
{
    if(itsGenSecurity != ARM_GenSecurityPtr(NULL))
    {
	    ARM_RowInfo inColumnNames=ColumnNames();

	    ARM_StringVector::iterator inTextIter			            = inColumnNames.first.begin();
	    ARM_StringVector::iterator inTextIterEnd		            = inColumnNames.first.end();
	    CC_NS(std,vector)<ARM_GP_VALUE_TYPE>::iterator inFormatIter	= inColumnNames.second.begin();

        ARM_DealDescriptionPtr outColumnNames=genSecurity->GetDealDescription().GetRow(0);

        size_t nbCols;
        ARM_StringVector outText = outColumnNames->GetText();
        CC_NS(std,vector)<ARM_GP_VALUE_TYPE> outFormat = outColumnNames->GetFormat();
        if((nbCols=inColumnNames.first.size()) == outText.size())
        {
	        ARM_StringVector::iterator outTextIter			            = outText.begin();
	        CC_NS(std,vector)<ARM_GP_VALUE_TYPE>::iterator outFormatIter	= outFormat.begin();
            for(;inTextIter!=inTextIterEnd;++inTextIter,++inFormatIter,++outTextIter,++outFormatIter)
                if(*inTextIter != *outTextIter || *inFormatIter != *outFormatIter)
                    break;
        }

        if(inTextIter != inTextIterEnd)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : can't set the given security because its format is not incompatible" );
    }
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: SetGenSecurity
///	Returns: void
///	Action : affectation of a Generic Security
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::SetGenSecurity( const ARM_GenSecurityPtr& genSecurity )
{
    /// Check the input GS consistency
    ValidateGenSecurity(genSecurity);

    itsGenSecurity = genSecurity;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: SetCalibMethod
///	Returns: void
///	Action : affectation of a Calib Methods
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::SetCalibMethod( const ARM_CalibMethodPtr& calibMethod )
{
    itsCalibMethod = calibMethod;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: SetPricingModel
///	Returns: void
///	Action : affectation of a Pricing Model with validation
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::SetPricingModel( const ARM_PricingModelPtr& pricingModel )
{
    itsPricingModel = pricingModel;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: SetMktDataManager
///	Returns: void
///	Action : affectation of a MDM
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::SetMktDataManager(const ARM_MarketData_ManagerRepPtr& mktDataManager)
{
    itsMktDataManager = mktDataManager;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: ComputeDomesticBasis
///	Returns: void
///	Action : Convert le basis from funding Ccy to domestic Ccy
/////////////////////////////////////////////////////////////////
ARM_GP_Vector ARM_GenCalculator::ComputeDomesticBasis() const 
{
	ARM_BasisConverter basisConveter(GetDomesticCcy(),
					GetFundingCcy(),
					GetRefDateStrip(),
					GetOriginalFundDateStrip(),
					GetRefFundDateStrip(),
					GetDomesticZeroCurve(),	
					GetForeignZeroCurve(),
					GetDomesticDiscountZeroCurve(),
					GetForeignDiscountZeroCurve(),
					*GetForex(),
					GetvCpnNominal(),
					GetvFundNominal(),
					GetvFundSpread());

	return  basisConveter.ComputeDomMargin();
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: FillRowInfo
///	Returns: void
///	Action : writes the row info into the textIter and formatIter
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::FillRowInfo( const ARM_RowInfo& rowInfo,
	CC_NS(std,vector)<string>::iterator& textIter, 
	CC_NS(std,vector)<ARM_GP_VALUE_TYPE>::iterator& formatIter )
{
	CC_NS( std,copy)( rowInfo.first.begin(), rowInfo.first.end(), textIter );
	CC_NS( std,copy)( rowInfo.second.begin(), rowInfo.second.end(), formatIter );

	size_t colsNb = rowInfo.first.size();

#if defined(__GP_STRICT_VALIDATION)
	if(rowInfo.first.size() != rowInfo.second.size() )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : text and format of different size" );
#endif

	textIter   += colsNb;
	formatIter += colsNb;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: CreateAndSetDealDescription
///	Returns: void
///	Action : default implementation for the deal description
///			based on the column names, first row, middlerows, lastrows
///			and datesStructure
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::CreateAndSetDealDescription(const string& payModelName, const ARM_StringVector& PricedColumns, ARM_CstManagerPtr cstManager, bool fixBoundary, bool otherPayoffs)
{
	ARM_RowInfo columnNames					= ColumnNames();
	ARM_DateStripCombiner datesStructure	= DatesStructure();
	ARM_VectorPtr eventDates				= datesStructure.GetMergeData();
	
	size_t colsNb	= columnNames.first.size();
	size_t rowsNb	= eventDates->size();

	/// fill first row with column names
	vector<string> text( colsNb * (rowsNb+1) );
	vector<ARM_GP_VALUE_TYPE> format( colsNb * (rowsNb+1), ARM_MISSING_TYPE );

	CC_NS(std,vector)<string>::iterator textIter			= text.begin();
	CC_NS(std,vector)<ARM_GP_VALUE_TYPE>::iterator formatIter	= format.begin();
	FillRowInfo( columnNames, textIter, formatIter );

    /// First skip all event dates before asOfDate
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
    size_t i=0;
	while(i < rowsNb && (*eventDates)[i] < asOfDate)
        ++i;
    if(i>=rowsNb)
    {
		CC_Ostringstream os;
		os << ARM_USERNAME << " : No event date posterior to the MDM asOfDate (" << GetMktDataManager()->GetAsOfDate().toString() << ")";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
    }


	/// Then fills rows one by one
	for(; i<rowsNb; ++i )
		FillRowInfo( MiddleRows( i, datesStructure ) , textIter, formatIter );


	ARM_DealDescriptionPtr pDealDescription( new ARM_DealDescription( text, format, rowsNb+1, colsNb, PricedColumns ) );

	itsGenSecurity = ARM_GenSecurityPtr( new ARM_GenSecurity( pDealDescription, payModelName, cstManager, !fixBoundary, otherPayoffs)  );
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: CreateAndSetDealDescription
///	Returns: void
///	Action : default implementation for the deal description
///			based on the column names, first row, middlerows, lastrows
///			and datesStructure
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::CreateAndSetCustomDealDescription(const ARM_DateStripVector& dateStrips, const string& payModelName, const ARM_StringVector& PricedColumns, ARM_CstManagerPtr cstManager, bool fixBoundary, bool otherPayoffs)
{
	ARM_RowInfo columnNames					= ColumnNames();
	ARM_DateStripCombiner datesStructure	= CustomDatesStructure(dateStrips);
	ARM_VectorPtr eventDates				= datesStructure.GetMergeData();
	
	size_t colsNb	= columnNames.first.size();
	size_t rowsNb	= eventDates->size();

	/// fill first row with column names
	vector<string> text( colsNb * (rowsNb+1) );
	vector<ARM_GP_VALUE_TYPE> format( colsNb * (rowsNb+1), ARM_MISSING_TYPE );

	CC_NS(std,vector)<string>::iterator textIter			= text.begin();
	CC_NS(std,vector)<ARM_GP_VALUE_TYPE>::iterator formatIter	= format.begin();
	FillRowInfo( columnNames, textIter, formatIter );

    /// First skip all event dates before asOfDate
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
    size_t i=0;
	while(i < rowsNb && (*eventDates)[i] < asOfDate)
        ++i;
    if(i>=rowsNb)
    {
		CC_Ostringstream os;
		os << ARM_USERNAME << " : No event date posterior to the MDM asOfDate (" << GetMktDataManager()->GetAsOfDate().toString() << ")";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
    }


	/// Then fills rows one by one
	for(; i<rowsNb; ++i )
		FillRowInfo( MiddleRows( i, datesStructure ) , textIter, formatIter );


	ARM_DealDescriptionPtr pDealDescription( new ARM_DealDescription( text, format, rowsNb+1, colsNb, PricedColumns ) );

	itsGenSecurity = ARM_GenSecurityPtr( new ARM_GenSecurity( pDealDescription, payModelName, cstManager, !fixBoundary, otherPayoffs)  );
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: Update
///	Returns: void
///	Action : update internal datas depending on the Market Data Manager
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::Update()
{
    /// Check data consistency before...
    CheckMktDataAndTimeIt();

    ///... update model datas and...
    UpdateModelAndTimeIt();

    ///... calibration datas (portfolio target prices)
    UpdateCalibrationAndTimeIt();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: Initialize
///	Returns: nothing
///	Action : Initialize and complete MktDataManager with Market Data
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::Initialize(ARM_MarketData_ManagerRep* mktDataManager)
{
	if(itsMktDataManager->GetAsOfDate() != mktDataManager->GetAsOfDate())
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : inconsistency between Market Data Managers, please, check ,as of date" );
	/// Register input objects but no more than the internal number of access keys
    size_t nbToReg = itsMDMKeys.size();
    for(int i=0; i<nbToReg; ++i)
	{
		if(!mktDataManager->TestIfKeyMissing(itsMDMKeys[i]))
			itsMktDataManager->RegisterData(itsMDMKeys[i],mktDataManager->GetData(itsMDMKeys[i]));
	}

	itsMktDataManager->SetDetailMode(mktDataManager->GetDetailMode());

	/// Check input datas before...
    CheckMktDataAndTimeIt();

    ///... creation of internal objects depending of the MarketDataManager
    CreateAndSetModelAndTimeIt();

    CreateAndSetCalibrationAndTimeIt();

    itsUpdate = true;	
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: Init
///	Returns: void
///	Action : initialisation of market data dependant objects : model &
///          calibration datas (Summit interface)
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::Init(const vector< ARM_Object* >& marketDatas, bool toCalibrate)
{
    size_t i;
    /// Save market datas in the MarketDataManager
    ARM_Date asOfDate;
    for(i=0;i<marketDatas.size();++i)
    {
        if(dynamic_cast< ARM_ZeroCurve* >(marketDatas[i]))
        {
            asOfDate = static_cast< ARM_ZeroCurve* >(marketDatas[i])->GetAsOfDate();
            break;
        }
    }

    if(itsMktDataManager == ARM_MarketData_ManagerRepPtr(NULL))
        itsMktDataManager = ARM_MarketData_ManagerRepPtr(new ARM_MarketData_ManagerRep(asOfDate));
    else if(itsMktDataManager->GetAsOfDate() != asOfDate)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : inconsistency between Market Data Manager and given zero curve as of dates" );

    /// Register input objects but no more than the internal number of access keys
    size_t nbToReg = (marketDatas.size() < itsMDMKeys.size() ? marketDatas.size() : itsMDMKeys.size());
    for(i=0;i<nbToReg;++i)
        itsMktDataManager->RegisterData(itsMDMKeys[i],marketDatas[i]);

    /// Check input datas before...
    CheckMktDataAndTimeIt();
    ///... creation of internal objects depending of the MarketDataManager
	if(toCalibrate)
	{
		CreateAndSetModelAndTimeIt();
		CreateAndSetCalibrationAndTimeIt();
	}

    itsUpdate = true;	
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: Update
///	Returns: void
///	Action : update market data dependant objects : model &
///          calibration datas (Summit interface)
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::Update(const vector< ARM_Object* >& marketDatas)
{
	if (itsUpdate)
	{
		/// Replace market datas in the MarketDataManager but no more than the internal number of access keys
		size_t nbToReg = (marketDatas.size() < itsMDMKeys.size() ? marketDatas.size() : itsMDMKeys.size());
		for(size_t i=0;i<nbToReg;++i)
			itsMktDataManager->RegisterData(itsMDMKeys[i],marketDatas[i]);

		/// Update the internal objects depending of the MarketDataManager
		Update();
	}
	else
		Init(marketDatas);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: Update
///	Returns: void
///	Action : update market data dependant objects : model &
///          calibration datas
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::Update(const ARM_MarketData_ManagerRep *mktDataMgr)
{
	if (!mktDataMgr)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : invalid market data manager, unable to update calculator !" );

	if(itsMktDataManager->GetAsOfDate() != mktDataMgr->GetAsOfDate())
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : inconsistency between Market Data Managers, please, check ,as of date" );

	/// Register input objects but no more than the internal number of access keys
    size_t nbKeys = itsMDMKeys.size();
    for(int i=0; i<nbKeys; ++i)
	{
		if (!mktDataMgr->TestIfKeyMissing(itsMDMKeys[i]))
			itsMktDataManager->RegisterData(itsMDMKeys[i],mktDataMgr->GetData(itsMDMKeys[i]));
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : invalid market data manager, inconsistent keys between market data manager !" );
	}

	if (itsUpdate)
		/// Update the internal objects depending of the MarketDataManager
		Update();
	else
	{
	    CheckMktDataAndTimeIt();
		CreateAndSetModelAndTimeIt();
		CreateAndSetCalibrationAndTimeIt();

		itsUpdate = true;	
	}
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: UpDateWithSubMktDataManager
///	Returns: void
///	Action : update market data dependant objects : model &
///          calibration datas (Summit interface)
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::UpDateWithSubMktDataManager(const ARM_MarketData_ManagerRep& mktDataManager)
{
    /// Replace market datas in the MarketDataManager but no more than the internal number of access keys
	ARM_MarketData_ManagerImp::iterator pos = mktDataManager.GetMktDataManager()->begin();
	for(; pos != mktDataManager.GetMktDataManager()->end(); ++pos)
	{
		if(!itsMktDataManager->TestIfKeyMissing((*pos).first))
			itsMktDataManager->RegisterData((*pos).first,&*(*pos).second);
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : you 're not authorized to update Mktmanager with a new key(" +(*pos).first+ "), please advice !" );
	}
		
	/// Update the internal objects depending of the MarketDataManager
    Update();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: UpDateOnlyMktDataManager
///	Returns: void
///	Action : update the market data dependant objects
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::UpDateMktDataManagerOnly(const ARM_MarketData_ManagerRep& mktDataManager)
{
    /// Replace market datas in the MarketDataManager but no more than the internal number of access keys
	ARM_MarketData_ManagerImp::iterator pos = mktDataManager.GetMktDataManager()->begin();
	for(; pos != mktDataManager.GetMktDataManager()->end(); ++pos)
	{
		if(!itsMktDataManager->TestIfKeyMissing((*pos).first))
			itsMktDataManager->RegisterData((*pos).first,&*(*pos).second);
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : you 're not authorized to update Mktmanager with a new key(" +(*pos).first+ "), please advice !" );
	}	
	
	CheckMktDataAndTimeIt();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: InitializeMktDataManagerOnly
///	Returns: void
///	Action : update the market data dependant objects
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::InitializeMktDataManagerOnly(const ARM_MarketData_ManagerRep& mktDataManager)
{
    /// Replace market datas in the MarketDataManager but no more than the internal number of access keys
	size_t i, size = itsMDMKeys.size();
	for(i=0; i<size; ++i)
		if(!mktDataManager.TestIfKeyMissing(itsMDMKeys[i]))
			itsMktDataManager->RegisterData(itsMDMKeys[i],mktDataManager.GetData(itsMDMKeys[i]));		
	
	CheckMktDataAndTimeIt();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: toString
///	Returns: string
///	Action : stringify the object
/////////////////////////////////////////////////////////////////

string ARM_GenCalculator::toString(const string& indent,const string& nextIndent) const
{
    CC_Ostringstream os;
	if( itsWarning )
	{
		os << " #    #    ##    #####   #    #     #    #    #   #### \n";
		os << " #    #   #  #   #    #  ##   #     #    ##   #  #    #\n";
		os << " #    #  #    #  #    #  # #  #     #    # #  #  #	  \n";
		os << " # ## #  ######  #####   #  # #     #    #  # #  #  ###\n";
		os << " ##  ##  #    #  #   #   #   ##     #    #   ##  #    #\n";
		os << " #    #  #    #  #    #  #    #     #    #    #   #### \n";
		os << itsWarning->toString(indent, nextIndent );
	}


	os << "\n\n =======> GENERIC CALCULATOR <====== \n";
	
    if(itsGenSecurity != ARM_GenSecurityPtr(NULL))
         os << indent << itsGenSecurity->toString();

    if(itsCalibMethod != ARM_CalibMethodPtr(NULL))
    {
        os << "\n" ;
        string tmpIndent = indent+ "|\t";
        os << "  GLOBAL CALIB METHOD OBJECT \n" ;
        os << itsCalibMethod->toString(tmpIndent);
        os << tmpIndent << "\n" ;
        os << "  ======> END OF INFOS CALIB METHOD <================== \n\n";
    }

	/// part on the pricing model
    if(itsPricingModel != ARM_PricingModelPtr(NULL))
        os << indent << itsPricingModel->toString(indent);

	/// part on the mkt data manager
    if(itsMktDataManager != ARM_MarketData_ManagerRepPtr(NULL))
        os << indent << itsMktDataManager->toString(indent);

	/// pricing data part
	ComputePricingData();
	os << itsPricingData.toString(indent,nextIndent);

    return os.str();
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: SetCurrencyUnit
///	Returns: void
///	Action : set the currency unit of the generic security
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::SetCurrencyUnit( ARM_Currency* ccy )
{	
	if( ccy )
	{
		delete itsCurrency;
		itsCurrency = (ARM_Currency*) ccy->Clone();
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: CallFuncAndTimeIt
///	Returns: void
///	Action : call func And Time It
/////////////////////////////////////////////////////////////////

void ARM_GenCalculator::CallFuncAndTimeIt( const GenCalcVoidToVoidFunc& func, const string& funcName )
{
	ARM_Timer timer;
	timer.ClockStartTime();
	(this->*func)();
	timer.ClockEndTime();
	itsPricingData[ "Duration_" + funcName ] = timer.GetDuration();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: CallFuncAndTimeIt
///	Returns: void
///	Action : call func And Time It
/////////////////////////////////////////////////////////////////

double ARM_GenCalculator::CallFuncAndTimeIt( const GenCalcDbleToVoidFunc& func, const string& funcName )
{
	ARM_Timer timer;
	timer.ClockStartTime();
	double result = (this->*func)();
	timer.ClockEndTime();
	itsPricingData[ "Duration_" + funcName ] = timer.GetDuration();
	return result;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: CreateAndSetCalibrationAndTimeIt
///	Returns: void
///	Action : call CreateAndSetCalibration And Time It
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::CreateAndSetCalibrationAndTimeIt()
{	
	CallFuncAndTimeIt( &ARM_GenCalculator::CreateAndSetCalibration, "CreateAndSetCalibration" ); 
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: CreateAndSetModelAndTimeIt
///	Returns: void
///	Action : call CreateAndSetModel And Time It
/////////////////////////////////////////////////////////////////
void ARM_GenCalculator::CreateAndSetModelAndTimeIt()
{	CallFuncAndTimeIt( &ARM_GenCalculator::CreateAndSetModel, "CreateAndSetModel" ); }

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: UpdateModelAndTimeIt
///	Returns: void
///	Action : call UpdateModel And Time It
/////////////////////////////////////////////////////////////////

void ARM_GenCalculator::UpdateModelAndTimeIt()
{	CallFuncAndTimeIt( &ARM_GenCalculator::UpdateModel, "UpdateModel" ); }

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: CalibrateAndTimeIt
///	Returns: void
///	Action : call Calibrate And Time It
/////////////////////////////////////////////////////////////////

void ARM_GenCalculator::CalibrateAndTimeIt()
{
	CallFuncAndTimeIt( &ARM_GenCalculator::Calibrate, "Calibrate" ); 

	delete itsWarning;
	itsWarning = 0;
	if( itsCalibMethod != ARM_CalibMethodPtr(NULL) )
		if( itsCalibMethod->GetWarning() )
			itsWarning = static_cast<ARM_Warning*>( itsCalibMethod->GetWarning()->Clone() );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: CheckDataAndTimeIt
///	Returns: void
///	Action : call CheckData And Time It
/////////////////////////////////////////////////////////////////

void ARM_GenCalculator::CheckDataAndTimeIt()
{	CallFuncAndTimeIt( &ARM_GenCalculator::CheckData, "CheckData" ); }


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: CheckMktDataAndTimeIt
///	Returns: void
///	Action : call CheckMktData And Time It
/////////////////////////////////////////////////////////////////

void ARM_GenCalculator::CheckMktDataAndTimeIt()
{	CallFuncAndTimeIt( &ARM_GenCalculator::CheckMktData, "CheckMktData" ); }

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: PriceAndTimeIt
///	Returns: void
///	Action : call Price And Time It
/////////////////////////////////////////////////////////////////
double ARM_GenCalculator::PriceAndTimeIt()
{	return CallFuncAndTimeIt( &ARM_GenCalculator::Price, "Price" ); }

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: UpdateCalibrationAndTimeIt
///	Returns: void
///	Action : call UpdateCalibration And Time It
/////////////////////////////////////////////////////////////////

void ARM_GenCalculator::UpdateCalibrationAndTimeIt(bool isUpdateStrike)
{
	ARM_Timer timer;
	timer.ClockStartTime();
	UpdateCalibration(isUpdateStrike);
	timer.ClockEndTime();
	itsPricingData[ "Duration_UpdateCalib" ] = timer.GetDuration();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: CreateAndSetDealDescription
///	Returns: void
///	Action : call CreateAndSetDealDescription And Time It
/////////////////////////////////////////////////////////////////

void ARM_GenCalculator::CreateAndSetDealDescriptionAndTimeIt(const string& payModelName, const ARM_StringVector& PricedColumns, ARM_CstManagerPtr cstManager, bool fixBoundary, bool otherPayoffs)
{
	ARM_Timer timer;
	timer.ClockStartTime();
	CreateAndSetDealDescription(payModelName, PricedColumns, cstManager,fixBoundary,otherPayoffs);
	timer.ClockEndTime();
	itsPricingData[ "Duration_CreateAndSetDealDes" ] = timer.GetDuration();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: CreateAndSetDealDescription
///	Returns: void
///	Action : call CreateAndSetDealDescription And Time It
/////////////////////////////////////////////////////////////////

void ARM_GenCalculator::CreateAndSetCustomDealDescriptionAndTimeIt(const ARM_DateStripVector& dateStrips, const string& payModelName, const ARM_StringVector& PricedColumns, ARM_CstManagerPtr cstManager, bool fixBoundary, bool otherPayoffs)
{
	ARM_Timer timer;
	timer.ClockStartTime();
	CreateAndSetCustomDealDescription(dateStrips, payModelName, PricedColumns, cstManager,fixBoundary,otherPayoffs);
	timer.ClockEndTime();
	itsPricingData[ "Duration_CreateAndSetCustDealDes" ] = timer.GetDuration();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: GetSubGenSecurity
///	Returns: ARM_GenSecurity*
///	Action : ability to get a sub security
/////////////////////////////////////////////////////////////////

ARM_GenSecurity* ARM_GenCalculator::GetSubGenSecurity(int prodIdx) const 
{
	/// Build the sub-generic security
	size_t startRowIdx = 1; /// start with the first actual row (=> skip column names)
	size_t endRowIdx = itsGenSecurity->GetDealDescription().GetRowsNb() - 1;
	ARM_DealDescriptionPtr subDealDesc = itsGenSecurity->GetDealDescription().GetSubDescription(startRowIdx,endRowIdx,prodIdx+1);
	ARM_GenSecurity* genSec = new ARM_GenSecurity(subDealDesc,GetGenSecurity()->GetPayModelName());
	return genSec;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

