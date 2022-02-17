/*!
 *
 * Copyright (c) CDC IXIS CM February 2005 Paris
 *
 *	\file vanillaargnumeric.cpp
 *
 *  \brief vanilla spreadoption
 *
 *	\author  Y KHLIF
 *	\version 1.0
 *	\date February 2005
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpcalib/vanillaargnumeric.h"
#include "gpcalib/vanillaarg.h"

#include "gpbase/datestripcombiner.h"
#include "gpbase/singleton.h"
#include "gpbase/autocleaner.h"
#include "gpbase/countedptr.h"
#include "gpbase/port.h"

/// gpnummethods
#include "gpnummethods/mcmethod.h"

/// gpnumlib
#include "gpnumlib/randomgenfactory.h"
#include "gpnumlib/compositegen.h"
#include "gpnumlib/antitheticgen.h"

/// gpnummethods
#include "gpnummethods/treemethod.h"
#include "gpnummethods/treebase.h"
#include "gpnummethods/treefactory.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/reconnector.h"
#include "gpnummethods/smoother.h"
#include "gpnummethods/normalcentredsampler.h"
#include "gpnummethods/scheduler.h"

#include "gpinfra/genpricer.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/pricingadviser.h"


CC_BEGIN_NAMESPACE( ARM )

const double	MC_NB_INTER_STEPS		= 1;
const int		MC_SIMULATION			=500;

const int		TREE_NBSTEPS_PER_YEAR	=2;
const double	STD_DEV_RATIO			=3.0;
const double	MIN_STD_DEV				=0.001;


////////////////////////////////////////////////////
///	Struct : ARM_VanillaArgNumeric
///	Routine: default constructor, copy constructor,
///          assigment, destructor, clone
///	Returns: 
///	Action :
////////////////////////////////////////////////////
void ARM_VanillaArgNumeric::CopyNoCleanUp(const ARM_VanillaArgNumeric& rhs)
{
	itsGenSecurity          = ARM_GenSecurityPtr((rhs.itsGenSecurity != ARM_GenSecurityPtr(NULL)) ? (ARM_GenSecurity*) rhs.itsGenSecurity->Clone() : NULL);
	itsNumMethod	        = ARM_NumMethodPtr((rhs.itsNumMethod != ARM_NumMethodPtr(NULL)) ? (ARM_NumMethod*) rhs.itsNumMethod->Clone() : NULL);
	itsCstManager			= ARM_CstManagerPtr((rhs.itsCstManager != ARM_CstManagerPtr(NULL)) ? (ARM_CstManager*) rhs.itsCstManager->Clone() : NULL);

	itsGenSecurity->SetCstManager(itsCstManager);
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaArgNumeric
///	Routine: constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_VanillaArgNumeric::ARM_VanillaArgNumeric(const string& curveName,
					double evalTime,
					int CallPut,
					double expiryTime)
:	ARM_VanillaArg( curveName, evalTime, CallPut,expiryTime ),	
	isGenSecurityCreated(false),
	itsGenSecurity(NULL)
{	
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaArgNumeric
///	Routine: Copy constructor
///	Returns: 
///	Action : does nothing but necessary for linking
////////////////////////////////////////////////////
ARM_VanillaArgNumeric::ARM_VanillaArgNumeric(const ARM_VanillaArgNumeric& rhs)
:	ARM_VanillaArg(rhs),
	isGenSecurityCreated(rhs.isGenSecurityCreated)
{
	CopyNoCleanUp(rhs);
} 

////////////////////////////////////////////////////
///	Struct : ARM_VanillaArgNumeric
///	Routine: operator Equal
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
	
ARM_VanillaArgNumeric& ARM_VanillaArgNumeric::operator=(const ARM_VanillaArgNumeric& rhs)
{
	if( this != &rhs )
	{
		ARM_VanillaArg::operator =( rhs);
		isGenSecurityCreated = rhs.isGenSecurityCreated;
		CopyNoCleanUp(rhs);
	}
	return *this;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: ValidateGenSecurity
///	Returns: void
///	Action : Check if the security is able to replace the built-in one
/////////////////////////////////////////////////////////////////
void ARM_VanillaArgNumeric::ValidateGenSecurity(const ARM_GenSecurityPtr& genSecurity)
{
    if(itsGenSecurity != ARM_GenSecurityPtr(NULL))
    {
	    ARM_RowInfo inColumnNames=NumArgColumnNames();

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
void ARM_VanillaArgNumeric::SetGenSecurity( const ARM_GenSecurityPtr& genSecurity )
{
    /// Check the input GS consistency
    ValidateGenSecurity(genSecurity);

    itsGenSecurity = genSecurity;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: SetGenSecurity
///	Returns: void
///	Action : affectation of a Generic Security
/////////////////////////////////////////////////////////////////
void ARM_VanillaArgNumeric::SetNumMethod( const ARM_NumMethodPtr& nummethod )
{
    itsNumMethod = nummethod;
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaArgNumeric
///	Routine: CreateAndSetNumMethod
///	Returns: 
///	Action : price a swaption with a model checking that it
///				is derived from an interest rate model
////////////////////////////////////////////////////
void ARM_VanillaArgNumeric::CreateAndSetNumMethod(const ARM_Date& asOfDate)
{
	if(GetItsPricingDirection() ==  ARM_NumMethod::GP_FWDLOOKING)
	{
		// MRGK5 Random Generator
		ARM_RandomGeneratorPtr  pBaseRandomGen( ARM_RandGenFactory.Instance()->CreateRandGen( 
					ARM_RandGenFactoryImp::MRGK5,
					ARM_RandGenFactoryImp::UnknownTransformAlgo ) );

		/// antithetic box muller
		ARM_RandomGeneratorPtr normRandGen( ARM_RandGenFactory.Instance()->CreateRandGen( 
			ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
			ARM_RandGenFactoryImp::BoxMuller,
			pBaseRandomGen ) );


		//// antithetic variates!
		ARM_RandomGeneratorPtr numRandGen = ARM_RandomGeneratorPtr( ARM_RandGenFactory.Instance()->CreateRandGen( 
		ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
		ARM_RandGenFactoryImp::AntitheticOne,
		normRandGen ) );

		ARM_TimeStepPerYearScheduler scheduler(0);
		ARM_NormalCentredSamplerND sampler(&scheduler);

		ARM_MCMethod mcMethod(MC_SIMULATION,numRandGen, &sampler);
		((ARM_VanillaArgNumeric*)this)->SetNumMethod(ARM_NumMethodPtr(new ARM_MCMethod(mcMethod)));
	}
	else if (GetItsPricingDirection() ==  ARM_NumMethod::GP_BCKWDLOOKING)
	{
		ARM_GP_VectorPtr eventDates             = DatesStructure();
		size_t eventSize	= (*eventDates).size();

		double lastEventTime = (*eventDates)[eventSize-1]-asOfDate.GetJulian();
		int nbSteps=static_cast<int>(floor(TREE_NBSTEPS_PER_YEAR*lastEventTime/K_YEAR_LEN));

		//int nbSteps = TREE_NBSTEPS_PER_YEAR;
		int initialStep =static_cast<int> (nbSteps/10.0);
		int schedulerType=ARM_SchedulerBase::ConstantVariance;
		std::vector<double> schedulerDatas(3);
		schedulerDatas[0] = nbSteps;
		schedulerDatas[1] = initialStep;
		schedulerDatas[2] = MIN_STD_DEV;
		int samplerType=ARM_SamplerBase::NormalCentred;
		std::vector<double> samplerDatas;
		int truncatorType=ARM_TruncatorBase::StandardDeviation;
		std::vector<double> truncatorDatas(1,STD_DEV_RATIO);
		int reconnectorType=ARM_ReconnectorBase::Mean;
		int smootherType=ARM_SmootherBase::DoNothing;
		ARM_TreeBase* tree = ARM_TreeFactory.Instance()->CreateTreeND(2,schedulerType,schedulerDatas,
			samplerType,samplerDatas,truncatorType,truncatorDatas,false,reconnectorType,smootherType);
		((ARM_VanillaArgNumeric*)this)->SetNumMethod(ARM_NumMethodPtr(tree));
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenCalculator
///	Routine: FillRowInfo
///	Returns: void
///	Action : writes the row info into the textIter and formatIter
/////////////////////////////////////////////////////////////////
void ARM_VanillaArgNumeric::FillRowInfo( const ARM_RowInfo& rowInfo,
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
///	Routine: CreateAndSetGenSec
///	Returns: void
///	Action : default implementation for the deal description
///			based on the column names, first row, middlerows, lastrows
///			and datesStructure
/////////////////////////////////////////////////////////////////
void ARM_VanillaArgNumeric::CreateAndSetGenSec(const string& curvename, const ARM_Date& asOfDate)
{
	if(!isGenSecurityCreated)
	{
		ARM_RowInfo columnNames					= NumArgColumnNames();
		ARM_GP_VectorPtr eventDates             = DatesStructure();
					
		
		size_t colsNb	= columnNames.first.size();
		size_t rowsNb	= eventDates->size();

		/// fill first row with column names
		vector<string> text( colsNb * (rowsNb+1) );
		vector<ARM_GP_VALUE_TYPE> format( colsNb * (rowsNb+1), ARM_MISSING_TYPE );

		CC_NS(std,vector)<string>::iterator textIter			= text.begin();
		CC_NS(std,vector)<ARM_GP_VALUE_TYPE>::iterator formatIter	= format.begin();
		FillRowInfo( columnNames, textIter, formatIter );

		/// First skip all event dates before asOfDate
		//double asOfDate = 100;//GetMktDataManager()->GetAsOfDate().GetJulian();
		size_t i=0;
		while(i < rowsNb && (*eventDates)[i] < asOfDate.GetJulian())
			++i;
		if(i>=rowsNb)
		{
			CC_Ostringstream os;
			os << ARM_USERNAME << " : No event date in calibration is posterior to the asOfDate ";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}


		/// Then fills rows one by one
		for(; i<rowsNb; ++i )
			FillRowInfo( NumArgMiddleRows( i, eventDates ) , textIter, formatIter );


		ARM_DealDescriptionPtr pDealDescription( new ARM_DealDescription( text, format, rowsNb+1, colsNb ) );
		const ARM_CstManagerPtr cstManager = securityCstManager();

		setConstantManager( cstManager );
		itsGenSecurity = ARM_GenSecurityPtr( new ARM_GenSecurity( pDealDescription, curvename, cstManager)  );
		isGenSecurityCreated = true;
	}
}
////////////////////////////////////////////////////
///	Struct : ARM_VanillaArgNumeric
///	Routine: ImpliedVol
///	Returns: Exeption
///	Action :  No Implied Volatility 
////////////////////////////////////////////////////
double ARM_VanillaArgNumeric::ImpliedVol(ARM_PricingModel* model) const
{
    CC_Ostringstream os;
	os << ARM_USERNAME << " : No formula is valid to calculate Implied Volatilty";
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );

}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaArgNumeric
///	Routine: Price
///	Returns: 
///	Action : price an equity option with a model checking that it
///				is derived from an equity model
////////////////////////////////////////////////////
double ARM_VanillaArgNumeric::Price(ARM_PricingModel* model) const
{
	double price=0.0;
	if(model)
	{
		string curveName = GetCurveName();
		ARM_Date asofDate = model->GetAsOfDate();
		const_cast<ARM_VanillaArgNumeric* const>(this)->CreateAndSetGenSec(curveName, asofDate);
		ARM_GenSecurityPtr genSec = GetGenSecurity();
		//ARM_NumMethodPtr numMethod = GetNumMethod();
	
		if( itsNumMethod==ARM_NumMethodPtr(NULL))
		{
			const_cast<ARM_VanillaArgNumeric* const>(this)->CreateAndSetNumMethod(asofDate);				
		}
			

		if( (model->GetNumeraire()) == ARM_NumerairePtr(NULL))
		{
			ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
			model->SetNumeraire(numeraire);
		}			
		model->SetNumMethod(itsNumMethod);
					
		ARM_GenPricer* genPricer = new ARM_GenPricer( &*genSec,model);

		ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );

		price = genPricer->Price();		
	}
    
    else
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"Model is not a model: not derived from ");
  
    return price;
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaEqOption
///	Routine: toString
///	Returns: string
///	Action : stringify the object to give details about it
////////////////////////////////////////////////////
string ARM_VanillaArgNumeric::toString(const string& indent, const string& nextIndent) const
{ 
	return "ARM_VanillaArgNumeric"; 
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

