
#include "gpbase/removeidentifiedwarning.h"
#include "gpinfra/gramfunctoreq.h"

/// gpbase
#include "gpbase/checkinputs.h"
#include "gpbase/env.h"
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/datestrip.h"
#include "gpbase/singleton.h"
#include "gpbase/stringconvert.h"

/// gpinfra
#include "gpinfra/argconvdefault.h"
#include "gpinfra/gramfunctorargcheck.h"
#include "gpinfra/gramfunctorconv.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingfunctionequity.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/pricingstates.h"

///gpmodel
#include "gpmodels/HybridIRFX.h"
#include "gpmodels/ModelParams_EqFXBase.h"
#include "gpmodels/fxname.h"
//
////gpcalculators
#include "gpcalculators/forexvanilla.h"
#include "gpcalculators/fxvanillafactory.h"

/// kernel
//#include <inst/irindex.h>
#include <util/fromto.h>

CC_BEGIN_NAMESPACE( ARM )

#define ARM_CF_EPS 1.0e-13

////////////////////////////////////////////////////
///	Class  : ARM_GramFctor
///	Member : itsFuncName is more for debugging output
////////////////////////////////////////////////////

/// the corresponding static string
string ARM_GP_SpotFctor::itsFuncName			= "Spot";
string ARM_GP_FwdFctor::itsFuncName				= "Forward";
string ARM_GP_CallFctor::itsFuncName			= "Call";
string ARM_GP_GreekFctor::itsFuncName			= "Greek";
string ARM_GP_EqDigitalFctor::itsFuncName		= "Equity Digital";
string ARM_GP_CallStripFctor::itsFuncName		= "Call Strip";
string ARM_GP_CallSpreadFctor::itsFuncName	    = "CallSpread";
string ARM_GP_RangeAccrualFctor::itsFuncName	= "RangeAccrual";


////////////////////////////////////////////////////
////////// Spot Factor functor /////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_GP_SpotFctor
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_SpotFctor::SetDefaults( double evalDate,  ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 2, ARM_GP_SpotFctor::itsFuncName  );
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,		ARM_GP_SpotFctor::itsFuncName  );
	GPAF_CheckArgType( arg[1], GFAT_DATE_OR_DOUBLE_TYPE,	ARM_GP_SpotFctor::itsFuncName );

	itsModelName			= arg[0].GetString();
	itsModelEquity			= dynamic_cast< ARM_PricingFunctionEquity* >(mod);
    if(!itsModelEquity)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Equity or Fx compatible model required for SPOT keyword");
	ARM_Date evalDateDate(evalDate);
	SetSettlementDate( nodes, evalDateDate, itsModelEquity->GetSettlementCalendar(itsModelName), itsModelEquity->GetSettlementGap(itsModelName), arg, 1, "ARM_GP_SpotFctor::SetDefaults" );
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_SpotFctor
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_SpotFctor::GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		SetDefaults(evalDate, arg, mod, nodes );

		/// grab the inputs after the setDefaults
		itsEvalTime			    = GetAndCheckEvalTime( evalDate, mod, ARM_GP_SpotFctor::itsFuncName );
		ARM_Date settlementDate = arg[1].GetDate();
		itsSettlementTime		= mod->GetTimeFromDate(settlementDate);

		/// do not do it again!
		SetAlreadyComputed(true);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_SpotFctor
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the equity spot function and return 
///				the corresponding values
////////////////////////////////////////////////////
ARM_GramFctorArg ARM_GP_SpotFctor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	
	ARM_VectorPtr spot =  itsModelEquity->Forward( itsModelName, itsEvalTime, itsEvalTime, itsEvalTime, itsEvalTime, states );
	ARM_GramFctorArg res(spot);//ARM_GP_VectorPtr(new ARM_GP_Vector(*spot)));
	return res;
}



////////////////////////////////////////////////////
///	Class  : ARM_GP_SpotFctor
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_SpotFctor::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );
	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( 1, itsEvalTime )), ARM_AdditionalTimeInfoPtr(NULL) );
}


////////////////////////////////////////////////////
///////////  Forward value        //////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_GP_FwdFctor
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_FwdFctor::SetDefaults( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 4, ARM_GP_FwdFctor::itsFuncName  );
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,		ARM_GP_FwdFctor::itsFuncName  );
	GPAF_CheckArgType( arg[1], GFAT_DATE_TYPE,			ARM_GP_FwdFctor::itsFuncName  );
	GPAF_CheckArgType( arg[2], GFAT_DATE_OR_DOUBLE_TYPE,ARM_GP_FwdFctor::itsFuncName  );
	GPAF_CheckArgType( arg[3], GFAT_DATE_OR_DOUBLE_TYPE,ARM_GP_FwdFctor::itsFuncName  );

	itsModelName			= arg[0].GetString();
	itsModelEquity			= dynamic_cast< ARM_PricingFunctionEquity* >(mod);
    if(!itsModelEquity)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Equity or Fx compatible model required for SPOT keyword");
	ARM_Date expiryDate   = arg[1].GetDate();
	SetSettlementDate( nodes, expiryDate, itsModelEquity->GetSettlementCalendar(itsModelName), itsModelEquity->GetSettlementGap(itsModelName), arg, 2, "ARM_GP_FwdFctor::SetDefaults" );
	SetSettlementDate( nodes, expiryDate, itsModelEquity->GetSettlementCalendar(itsModelName), itsModelEquity->GetSettlementGap(itsModelName), arg, 3, "ARM_GP_FwdFctor::SetDefaults" );
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_FwdFctor
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_FwdFctor::GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		SetDefaults( arg, mod, nodes );

		/// grab inputs
		itsEvalTime			    = GetAndCheckEvalTime( evalDate, mod, ARM_GP_FwdFctor::itsFuncName );
		ARM_Date expiryDate	    = arg[1].GetDate();
		itsExpiryTime			= mod->GetTimeFromDate(expiryDate);
		ARM_Date settlementDate	= arg[2].GetDate();
		itsSettlementTime	    = mod->GetTimeFromDate(settlementDate);
		ARM_Date payDate	    = arg[3].GetDate();
		itsPayTime	            = mod->GetTimeFromDate(payDate);

		/// some validation
		CheckNbSmaller( itsEvalTime, itsExpiryTime,	"EvalTime", "ExpiryTime",  ARM_GP_FwdFctor::itsFuncName );
		CheckNbSmaller( itsExpiryTime, itsSettlementTime,	"ExpiryTime", "SettlementTime",  ARM_GP_FwdFctor::itsFuncName );
		CheckNbSmaller( itsExpiryTime, itsPayTime,	"ExpiryTime", "PayTime",  ARM_GP_FwdFctor::itsFuncName );

		/// do not do it again!
		SetAlreadyComputed(true);
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_GP_FwdFctor
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the equity forward function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_FwdFctor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );
	ARM_VectorPtr forward =  itsModelEquity->Forward( itsModelName, itsEvalTime, itsExpiryTime, itsSettlementTime, itsPayTime, states );

	/// return result
	return ARM_GramFctorArg(forward);
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_FwdFctor
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_FwdFctor::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );

	/// return result
	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( 1, itsExpiryTime )),ARM_AdditionalTimeInfoPtr( new ARM_AdditionalTimeInfo() ) );
}


////////////////////////////////////////////////////
///////////		 Base Equity      //////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_GP_EquityOptionAndGreeksFctor
///	Routine: ARM_GP_EquityOptionAndGreeksFctor
///	Returns: 
///	Action : Default constructor
////////////////////////////////////////////////////
ARM_GP_EquityOptionAndGreeksFctor::ARM_GP_EquityOptionAndGreeksFctor()
: ARM_GramFctor(), itsModelEquity(NULL), itsFXname(), itsEvalTime(0.), 
	itsExpiryTimeVector(), itsSettlementTimeVector(), itsStrikeVector(), itsPayTimeVector(),
	itsCallPut(1)
{
	// break point space
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_EquityOptionAndGreeksFctor
///	Routine: ~ARM_GP_EquityOptionAndGreeksFctor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_GP_EquityOptionAndGreeksFctor::~ARM_GP_EquityOptionAndGreeksFctor()
{
	itsModelEquity = NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_EquityOptionAndGreeksFctor
///	Routine: ARM_GP_EquityOptionAndGreeksFctor
///	Returns: 
///	Action : copy constructor
////////////////////////////////////////////////////
ARM_GP_EquityOptionAndGreeksFctor::ARM_GP_EquityOptionAndGreeksFctor(const ARM_GP_EquityOptionAndGreeksFctor& from)
: ARM_GramFctor(from), itsModelEquity(from.itsModelEquity), itsFXname(from.itsFXname), itsEvalTime(from.itsEvalTime), 
	itsExpiryTimeVector(from.itsExpiryTimeVector), itsSettlementTimeVector(from.itsSettlementTimeVector), 
	itsStrikeVector(from.itsStrikeVector), itsPayTimeVector(from.itsPayTimeVector),
	itsCallPut(from.itsCallPut)
{
	// break point space
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_EquityOptionAndGreeksFctor
///	Routine: operator=
///	Returns: 
///	Action : Assignement operator
////////////////////////////////////////////////////
ARM_GP_EquityOptionAndGreeksFctor& ARM_GP_EquityOptionAndGreeksFctor::operator=(const ARM_GP_EquityOptionAndGreeksFctor& from)
{
	if ( &from != this )
	{
		ARM_GramFctor::operator=(from);

		itsModelEquity			= from.itsModelEquity;
		itsFXname				= from.itsFXname;
		itsEvalTime				= from.itsEvalTime;
		itsExpiryTimeVector		= from.itsExpiryTimeVector;
		itsSettlementTimeVector	= from.itsSettlementTimeVector;
		itsStrikeVector			= from.itsStrikeVector;
		itsPayTimeVector		= from.itsPayTimeVector;
		itsCallPut				= from.itsCallPut;
	}

  return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_EquityOptionAndGreeksFctor
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////
void ARM_GP_EquityOptionAndGreeksFctor::SetDefaults( 
										ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
										vector< ARM_ExpNodePtr >& nodes 
										)
{
	this->checkArgument(arg);

	itsFXname		= arg[0].GetString();
	itsModelEquity	= dynamic_cast< ARM_PricingFunctionEquity* >(mod);
    if(!itsModelEquity)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Equity or Fx compatible model required for SPOT keyword");


	double defaultGap;
	string calendar;

	ARM_EqFxBase* modelFX = dynamic_cast< ARM_EqFxBase* >(mod);
	if(modelFX)//Pure FX
	{
		string marketFxName = ARM_FXName(itsFXname).GetMktName();

		defaultGap	= modelFX->GetSettlementGap(marketFxName);
		calendar	= modelFX->GetSettlementCalendar(marketFxName);
	}
	else {
		ARM_HybridIRFX* hybridModel	= dynamic_cast< ARM_HybridIRFX* >(mod);
		if(hybridModel)//For the pure quanto case
		{
			string marketFxName = ARM_FXName(itsFXname).GetMktName();

			defaultGap	= hybridModel->GetSettlementGap(marketFxName);
			calendar	= hybridModel->GetSettlementCalendar(marketFxName);
		}
		else {
			defaultGap	= itsModelEquity->GetSettlementGap(itsFXname);
			calendar	= itsModelEquity->GetSettlementCalendar(itsFXname);
		}
	}

	if ( (arg[1].GetType() == GFAT_DATE_TYPE) || (arg[1].GetType() == GFAT_DOUBLE_TYPE) )
	{
		ARM_Date expiryDate   = ( (arg[1].GetType() == GFAT_DATE_TYPE) ? arg[1].GetDate() : ARM_Date(arg[1].GetDouble()) );
		
		//Settlement and Payment Dates
		SetSettlementDate( nodes, expiryDate, calendar, defaultGap, arg, 4, "ARM_GP_EquityOptionAndGreeksFctor::SetDefaults" );
		SetSettlementDate( nodes, expiryDate, calendar, defaultGap, arg, 5, "ARM_GP_EquityOptionAndGreeksFctor::SetDefaults" );
	}
	else if ( arg[1].GetType() ==  GFAT_VECTOR_TYPE ) 
	{
		ARM_VectorPtr expiryDates = arg[1].GetVector();
		size_t sz = expiryDates->size();
		if ( sz > 0 ) {
			if ( (arg[4].GetType() ==  GFAT_DOUBLE_TYPE) || (arg[4].GetType() ==  GFAT_DATE_TYPE) )
			{
				ARM_VectorPtr settlements(NULL);

				if (arg[4].GetType() ==  GFAT_DOUBLE_TYPE) 
				{
					settlements = ARM_VectorPtr(new std::vector<double>( (*expiryDates) ) );
					double value = arg[4].GetDouble();
					for ( size_t i(0); i < sz; ++i) {
						(*settlements)[i] = SetSettlementDate((*settlements)[i],calendar,defaultGap,value);
					}
				}
				else {
					double value = arg[4].GetDate().GetJulian();
					settlements = ARM_VectorPtr(new std::vector<double>(sz,value) );
				}
				
				arg[4].SetVector(settlements);
			}

			if ( (arg[5].GetType() ==  GFAT_DOUBLE_TYPE) || (arg[5].GetType() ==  GFAT_DATE_TYPE) )
			{
				ARM_VectorPtr payments(NULL);

				if (arg[5].GetType() ==  GFAT_DOUBLE_TYPE) 
				{
					payments = ARM_VectorPtr(new std::vector<double>( (*expiryDates) ) );
					double value = arg[5].GetDouble();
					for ( size_t i(0); i < sz; ++i) {
						(*payments)[i] = SetSettlementDate((*payments)[i],calendar,defaultGap,value);
					}
				}
				else {
					double value = arg[5].GetDate().GetJulian();
					payments = ARM_VectorPtr(new std::vector<double>(sz,value) );
				}

				arg[5].SetVector(payments);
			}
		}

	}
	else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "expiryDate should be double vector or date" );

}

////////////////////////////////////////////////////
///	Class  : ARM_GP_EquityOptionAndGreeksFctor
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_EquityOptionAndGreeksFctor::GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes, size_t statesSize  )
{

	/*
	std::string fileName = "c:/Debug/StatAverageVector_operator";
#if defined( _DEBUG )
	fileName += ".dbg";
#else
	fileName += ".release";
#endif
	fileName += ".txt";

	std::ofstream os;
	os.open(fileName.c_str(),ios::out | ios::app );

	if (!os) {
		std::string msg = "Failed to open file : " + fileName;
		os.close();
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg.c_str() );
	}
	os << "statesSize = " << statesSize << std::endl;
	*/

	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		SetDefaults( arg, mod, nodes );

		/// handle vectorial strike
		if( arg[2].GetType() ==  GFAT_DOUBLE_TYPE )
		{
			//os << std::endl << "Strike ( arg[2].GetType() ==  GFAT_DOUBLE_TYPE )" << std::endl;
			double strikeDouble	= arg[2].GetDouble();
			//os << "strikeDouble	= " << strikeDouble << std::endl;
			itsStrikeVector		= ARM_VectorPtr(new std::vector<double>(statesSize,strikeDouble));
			//os << "itsStrikeVector size = " << itsStrikeVector->size() << std::endl;

		}
		else if( arg[2].GetType() ==  GFAT_VECTOR_TYPE )
		{
			//os << std::endl << "Strike ( arg[2].GetType() ==  GFAT_VECTOR_TYPE )" << std::endl;
			itsStrikeVector	= arg[2].GetVector();
			/*
			if ( !itsStrikeVector.IsNull() )
				os << "itsStrikeVector size = " << itsStrikeVector->size() << std::endl;
			*/
		}
		else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "strike should be vector or double" );

		/// grab inputs
		/// get model
		itsEvalTime			    = GetAndCheckEvalTime( evalDate, mod, string("ARM_GP_EquityOptionAndGreeksFctor") );

		string callPutString    = arg[3].GetString();
		itsCallPut			    = ARM_ArgConv_CallPut.GetNumber( callPutString );

		if ( (arg[1].GetType() ==  GFAT_DATE_TYPE) || (arg[1].GetType() ==  GFAT_DOUBLE_TYPE) )
		{
			//os << "(arg[1].GetType() ==  GFAT_DATE_TYPE) || (arg[1].GetType() ==  GFAT_DOUBLE_TYPE)" << std::endl;
			ARM_Date expiryDate     = ( (arg[1].GetType() == GFAT_DATE_TYPE) ? arg[1].GetDate() : ARM_Date(arg[1].GetDouble()) );//arg[1].GetDate();
			double expiryTime	= mod->GetTimeFromDate(expiryDate);
			ARM_Date settlementDate = ( (arg[4].GetType() == GFAT_DATE_TYPE) ? arg[4].GetDate() : ARM_Date(arg[4].GetDouble()) );//arg[4].GetDate();
			double settlementTime	= mod->GetTimeFromDate(settlementDate);
			ARM_Date payDate        = ( (arg[5].GetType() == GFAT_DATE_TYPE) ? arg[5].GetDate() : ARM_Date(arg[5].GetDouble()) );//arg[5].GetDate();
			double payTime		= mod->GetTimeFromDate(payDate);

			/*
			os << "Expiry     = (" << expiryDate.toString() << ", " << expiryDate.GetJulian() << ", " << expiryTime << ")" << std::endl;
			os << "Settlement = (" << settlementDate.toString() << ", " << settlementDate.GetJulian() << ", " << settlementTime << ")" << std::endl;
			os << "payDate    = (" << payDate.toString() << ", " << payDate.GetJulian() << ", " << payTime << ")" << std::endl;
			*/

			/// Validation
			CheckNbSmaller( itsEvalTime, expiryTime,	"EvalTime", "ExpiryTime",  string("ARM_GP_EquityOptionAndGreeksFctor") );
			CheckNbSmaller( expiryTime, settlementTime,	"ExpiryTime", "SettlementTime",  string("ARM_GP_EquityOptionAndGreeksFctor") );
			CheckNbSmaller( expiryTime, payTime,	"ExpiryTime", "PayTime",  string("ARM_GP_EquityOptionAndGreeksFctor") );

			itsExpiryTimeVector		= ARM_VectorPtr(new std::vector<double>(statesSize ,expiryTime));
			itsSettlementTimeVector	= ARM_VectorPtr(new std::vector<double>(statesSize ,settlementTime));
			itsPayTimeVector		= ARM_VectorPtr(new std::vector<double>(statesSize ,payTime));
			/*
			os << "Expiry size = " << itsExpiryTimeVector->size() << std::endl;
			os << "Settlement size = " << itsSettlementTimeVector->size() << std::endl;
			os << "PayTime size    = " << itsPayTimeVector->size() << std::endl;
			*/

		}
		else if ( arg[1].GetType() ==  GFAT_VECTOR_TYPE ) 
		{
			//os << "( arg[1].GetType() ==  GFAT_VECTOR_TYPE )" << std::endl;
			itsExpiryTimeVector		= arg[1].GetVector();
			if ( itsExpiryTimeVector.IsNull() )
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "expiry vector ptr is NULL" );
			size_t i, sz = itsExpiryTimeVector->size();
			//os << "Expiry size = " << sz << std::endl;

			if (sz>0) {
				if ( sz != statesSize)
					Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "something wrong in ARM_GP_EquityOptionAndGreeksFctor::GrabInputs(...)" );

				itsSettlementTimeVector	= arg[4].GetVector();
				itsPayTimeVector		= arg[5].GetVector();
				/*
				if ( !itsSettlementTimeVector.IsNull() )
					os << "Settlement size = " << itsSettlementTimeVector->size() << std::endl;
				if ( !itsPayTimeVector.IsNull() )
					os << "PayTime size    = " << itsPayTimeVector->size() << std::endl;
				*/

				// Convert vector of julian in vector of year fraction
				//os << "[" << std::endl;
				double time;
				for ( i = 0; i < sz; ++i) {
					
					ARM_Date dummyDate = ARM_Date( (*itsExpiryTimeVector)[i] );
					time = mod->GetTimeFromDate(dummyDate);
					(*itsExpiryTimeVector)[i] = time;

					//os << i << "Expiry = (" << dummyDate.toString() << ", " << dummyDate.GetJulian() << ", " << time << ")" << std::endl;

					dummyDate = ARM_Date( (*itsSettlementTimeVector)[i] );
					time = mod->GetTimeFromDate(dummyDate);
					(*itsSettlementTimeVector)[i] = time;

					//os << i << "Settlement = (" << dummyDate.toString() << ", " << dummyDate.GetJulian() << ", " << time << ")" << std::endl;

					dummyDate = ARM_Date( (*itsPayTimeVector)[i] );
					time = mod->GetTimeFromDate(dummyDate);
					(*itsPayTimeVector)[i] = time;

					//os << i << "PayTime = (" << dummyDate.toString() << ", " <<  dummyDate.GetJulian() << ", " << time << ")" << std::endl;

				}
				//os << "]" << std::endl;

				// Check sizes
				if ( (itsSettlementTimeVector->size() != sz) || (itsPayTimeVector->size()!= sz) ) throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "expiryTimeVector, itsSettlementTimeVector and itsPayTimeVector must have the same size" );

				// Validation
				for (i = 0; i< sz; ++i) {
					CheckNbSmaller( itsEvalTime, (*itsExpiryTimeVector)[i],	"EvalTime", "ExpiryTime",  string("ARM_GP_EquityOptionAndGreeksFctor") );
					CheckNbSmaller( (*itsExpiryTimeVector)[i], (*itsSettlementTimeVector)[i],	"ExpiryTime", "SettlementTime",  string("ARM_GP_EquityOptionAndGreeksFctor") );
					CheckNbSmaller( (*itsExpiryTimeVector)[i], (*itsPayTimeVector)[i],	"ExpiryTime", "PayTime",  string("ARM_GP_EquityOptionAndGreeksFctor") );
				}
			}
			//else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "expiry vector is empty" );
		}
		else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "expiry should be vector or date" );

		GrabMoreInputs(arg);

		SetAlreadyComputed(true);
/*
		os << std::endl;
		os.close();
*/
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_EquityOptionAndGreeksFctor
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_EquityOptionAndGreeksFctor::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );
	/// return result
	// NDC : CHECK THAT
	return ARM_NodeInfo( itsExpiryTimeVector,ARM_AdditionalTimeInfoPtr(NULL) );
	//NDC : the following line leads to something wrong, i don't know why
	//return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( 1, (*itsExpiryTimeVector)[0])),ARM_AdditionalTimeInfoPtr(NULL) );
}



////////////////////////////////////////////////////
///////////		 Call Option      //////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_GP_CallFctor
///	Routine: checkArgument
///	Returns: void
///	Action : check arguments number and type.
////////////////////////////////////////////////////
void ARM_GP_CallFctor::checkArgument( ARM_GramFctorArgVector& arg )
{
	/// checking of the size and type
 	GPAF_CheckArgSize( arg, 6, ARM_GP_CallFctor::itsFuncName);
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,			ARM_GP_CallFctor::itsFuncName);
	GPAF_CheckArgType( arg[1], GFAT_DATE_OR_VECTOR_TYPE,	ARM_GP_CallFctor::itsFuncName);
	GPAF_CheckArgType( arg[2], GFAT_VECTOR_TYPE,			ARM_GP_CallFctor::itsFuncName);
	GPAF_CheckArgType( arg[3], GFAT_STRING_TYPE,			ARM_GP_CallFctor::itsFuncName);
	GPAF_CheckArgType( arg[4], GFAT_DATE_OR_VECTOR_TYPE,	ARM_GP_CallFctor::itsFuncName );
	GPAF_CheckArgType( arg[5], GFAT_DATE_OR_VECTOR_TYPE,	ARM_GP_CallFctor::itsFuncName );
}



////////////////////////////////////////////////////
///	Class  : ARM_GP_CallFctor
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the call function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_CallFctor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	size_t statesSize = (false == states.IsNull())? MAX(states->size(),1) : 1;
	ARM_VectorPtr callValue( new std::vector<double>(statesSize,0.0) );

	GrabInputs( arg, mod, evalDate, nodes, statesSize);

	callValue =  itsModelEquity->CallVectorial( itsFXname, itsEvalTime, *itsExpiryTimeVector, *itsSettlementTimeVector, *itsStrikeVector, itsCallPut, *itsPayTimeVector, states );

	if (callValue.IsNull())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "callValue vector ptr is NULL" );

	/*
	std::string fileName = "c:/Debug/StatAverageVector_operator";
#if defined( _DEBUG )
	fileName += ".dbg";
#else
	fileName += ".release";
#endif
	fileName += ".txt";

	std::ofstream os;
	os.open(fileName.c_str(),ios::out | ios::app );

	if (!os) {
		std::string msg = "Failed to open file : " + fileName;
		os.close();
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg.c_str() );
	}
	

	os << "statesSize = " << statesSize << ",  callValue size = " << callValue->size() << std::endl;

	os << "[ ";
	for (size_t i(0); i < callValue->size(); ++i)
	{
		os << (*callValue)[i] << ", ";

	}
	os << "]" << std::endl;

	os.close();
	*/

	/// return result
	return ARM_GramFctorArg(callValue);
}



////////////////////////////////////////////////////
///////////		 Call Greek      //////////////////
////////////////////////////////////////////////////

void ARM_GP_GreekFctor::GrabMoreInputs( ARM_GramFctorArgVector& arg )
{
	greekType_ = arg[6].GetString();
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_GreekFctor
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////
void ARM_GP_GreekFctor::checkArgument( ARM_GramFctorArgVector& arg )
{
	/// checking of the size and type
 	GPAF_CheckArgSize( arg, 7, ARM_GP_GreekFctor::itsFuncName);
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,			ARM_GP_GreekFctor::itsFuncName);
	GPAF_CheckArgType( arg[1], GFAT_DATE_OR_VECTOR_TYPE,	ARM_GP_GreekFctor::itsFuncName);
	GPAF_CheckArgType( arg[2], GFAT_VECTOR_TYPE,			ARM_GP_GreekFctor::itsFuncName);
	GPAF_CheckArgType( arg[3], GFAT_STRING_TYPE,			ARM_GP_GreekFctor::itsFuncName);
	GPAF_CheckArgType( arg[4], GFAT_DATE_OR_VECTOR_TYPE,	ARM_GP_GreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[5], GFAT_DATE_OR_VECTOR_TYPE,	ARM_GP_GreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[6], GFAT_STRING_TYPE,			ARM_GP_GreekFctor::itsFuncName );
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_GreekFctor
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the call function and return 
///				the corresponding values
////////////////////////////////////////////////////
ARM_GramFctorArg ARM_GP_GreekFctor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	size_t statesSize = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;
	ARM_VectorPtr greekValue( new std::vector<double>(statesSize,0.0) );

	GrabInputs( arg, mod, evalDate, nodes, statesSize);

	greekValue =  itsModelEquity->GreekVectorial( itsFXname, itsEvalTime, *itsExpiryTimeVector, *itsSettlementTimeVector, *itsStrikeVector, itsCallPut, *itsPayTimeVector, greekType_, states );

	if (greekValue.IsNull())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "greekValue vector ptr is NULL" );

	/// return result
	return ARM_GramFctorArg(greekValue);
}



////////////////////////////////////////////////////
///////////		 CallSpread Option      //////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_GP_CallSpreadFctor
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_CallSpreadFctor::SetDefaults( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 9, ARM_GP_CallSpreadFctor::itsFuncName);
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,			ARM_GP_CallSpreadFctor::itsFuncName );//modelName1
	GPAF_CheckArgType( arg[1], GFAT_STRING_TYPE,			ARM_GP_CallSpreadFctor::itsFuncName );//modelName2
	GPAF_CheckArgType( arg[2], GFAT_DATE_TYPE,				ARM_GP_CallSpreadFctor::itsFuncName );//ResetDate
	GPAF_CheckArgType( arg[3], GFAT_VECTOR_TYPE,			ARM_GP_CallSpreadFctor::itsFuncName );//Strike
	GPAF_CheckArgType( arg[4], GFAT_DOUBLE_TYPE,			ARM_GP_CallSpreadFctor::itsFuncName );//alpha
	GPAF_CheckArgType( arg[5], GFAT_DOUBLE_TYPE,			ARM_GP_CallSpreadFctor::itsFuncName );//beta
	GPAF_CheckArgType( arg[6], GFAT_DATE_OR_DOUBLE_TYPE,	ARM_GP_CallSpreadFctor::itsFuncName );//Settlement gap1
	GPAF_CheckArgType( arg[7], GFAT_DATE_OR_DOUBLE_TYPE,	ARM_GP_CallSpreadFctor::itsFuncName );//Settlement gap2
	GPAF_CheckArgType( arg[8], GFAT_DATE_OR_DOUBLE_TYPE,	ARM_GP_CallSpreadFctor::itsFuncName );//Pay gap

	itsFX1Name			= arg[0].GetString();
	itsFX2Name			= arg[1].GetString();

	//Spread de FX for the moment
	itsHybridModel		= dynamic_cast< ARM_HybridIRFX* >(mod);
    if(!itsHybridModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only Spread of FX for the moment");

	ARM_FXName  fx1Name(itsFX1Name.substr(itsFX1Name.size()-6));    
	ARM_FXName  fx2Name(itsFX1Name.substr(itsFX2Name.size()-6));

	ARM_Date expiryDate   = arg[2].GetDate();
	SetSettlementDate( nodes, expiryDate, itsHybridModel->GetSettlementCalendar(fx1Name.GetMktName()), itsHybridModel->GetSettlementGap(fx1Name.GetMktName()), arg, 6, "ARM_GP_CallSpreadFctor::SetDefaults" );
	SetSettlementDate( nodes, expiryDate, itsHybridModel->GetSettlementCalendar(fx2Name.GetMktName()), itsHybridModel->GetSettlementGap(fx2Name.GetMktName()), arg, 7, "ARM_GP_CallSpreadFctor::SetDefaults" );
	//payment date (by convention equals to the settlemet date of the first FX for the moment)
	SetSettlementDate( nodes, expiryDate, itsHybridModel->GetSettlementCalendar(fx1Name.GetMktName()), itsHybridModel->GetSettlementGap(fx1Name.GetMktName()), arg, 8, "ARM_GP_CallSpreadFctor::SetDefaults" );
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_CallSpreadFctor
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_CallSpreadFctor::GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes  )
{
	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		SetDefaults( arg, mod, nodes );

		/// handle vectorial strike
		if( arg[3].GetType() ==  GFAT_DOUBLE_TYPE )
		{
			itsStrikeVector	= ARM_VectorPtr( new std::vector<double>( 1,arg[3].GetDouble() ) );
		}
		else if( arg[3].GetType() ==  GFAT_VECTOR_TYPE )
		{
			itsStrikeVector	= arg[3].GetVector();
		}
		else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "strike should be vector or double" );
		/// handle vectorial alpha
		if( arg[4].GetType() ==  GFAT_DOUBLE_TYPE )
		{
			itsAlphaDouble = arg[4].GetDouble();
		}
		else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "alpha should be a double" );
		/// handle vectorial beta
		if( arg[5].GetType() ==  GFAT_DOUBLE_TYPE )
		{
			itsBetaDouble   = arg[5].GetDouble();
		}
		else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "beta should be a double" );

		/// grab inputs
		/// get model
		itsEvalTime			    = GetAndCheckEvalTime( evalDate, mod, ARM_GP_CallSpreadFctor::itsFuncName );
		ARM_Date expiryDate     = arg[2].GetDate();
		itsExpiryTime			= mod->GetTimeFromDate(expiryDate);
		ARM_Date settlementDate1 = arg[6].GetDate();
		itsSettlementTime1		= mod->GetTimeFromDate(settlementDate1);
		ARM_Date settlementDate2 = arg[7].GetDate();
		itsSettlementTime2		= mod->GetTimeFromDate(settlementDate2);
		ARM_Date payDate        = arg[8].GetDate();
		itsPayTime		        = mod->GetTimeFromDate(payDate);

		/// Validation
		CheckNbSmaller( itsEvalTime, itsExpiryTime,	"EvalTime", "ExpiryTime",  ARM_GP_CallSpreadFctor::itsFuncName );
		CheckNbSmaller( itsExpiryTime, itsSettlementTime1,	"ExpiryTime", "SettlementTime",  ARM_GP_CallSpreadFctor::itsFuncName );
		CheckNbSmaller( itsExpiryTime, itsSettlementTime2,	"ExpiryTime", "SettlementTime",  ARM_GP_CallSpreadFctor::itsFuncName );
		CheckNbSmaller( itsExpiryTime, itsPayTime,	"ExpiryTime", "PayTime",  ARM_GP_CallSpreadFctor::itsFuncName );

		SetAlreadyComputed(true);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_CallSpreadFctor
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the call spread function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_CallSpreadFctor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );
	/// init of call vectorial
	size_t payoffSize = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;
	ARM_VectorPtr callspreadValue;
	callspreadValue =  itsHybridModel->CallSpreadVectorial( itsFX1Name,
		itsFX2Name, 
		itsEvalTime, 
		itsExpiryTime, 
		itsSettlementTime1,
		itsSettlementTime2,
		itsPayTime,
		*itsStrikeVector, 
		itsAlphaDouble, 
		itsBetaDouble, 
		states );
	/// return result
	return ARM_GramFctorArg(callspreadValue);
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_CallSpreadFctor
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_CallSpreadFctor::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );

	/// return result
	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( 1, itsExpiryTime)),ARM_AdditionalTimeInfoPtr(NULL) );
}



////////////////////////////////////////////////////
///////////		Equity Digital	  //////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_GP_EqDigitalFctor
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_EqDigitalFctor::SetDefaults( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 9, ARM_GP_EqDigitalFctor::itsFuncName);
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,			ARM_GP_EqDigitalFctor::itsFuncName);
	GPAF_CheckArgType( arg[1], GFAT_DATE_TYPE,				ARM_GP_EqDigitalFctor::itsFuncName);
	GPAF_CheckArgType( arg[2], GFAT_VECTOR_TYPE,			ARM_GP_EqDigitalFctor::itsFuncName);
	GPAF_CheckArgType( arg[3], GFAT_STRING_TYPE,			ARM_GP_EqDigitalFctor::itsFuncName);
	GPAF_CheckArgType( arg[4], GFAT_DOUBLE_TYPE,			ARM_GP_EqDigitalFctor::itsFuncName);
	GPAF_CheckArgType( arg[5], GFAT_DATE_OR_DOUBLE_TYPE,	ARM_GP_EqDigitalFctor::itsFuncName );
	GPAF_CheckArgType( arg[6], GFAT_DATE_OR_DOUBLE_TYPE,	ARM_GP_EqDigitalFctor::itsFuncName );
	GPAF_CheckArgType( arg[7], GFAT_STRING_TYPE,			ARM_GP_EqDigitalFctor::itsFuncName);
	GPAF_CheckArgType( arg[8], GFAT_DOUBLE_TYPE,			ARM_GP_EqDigitalFctor::itsFuncName);

	itsModelName			= arg[0].GetString();
	itsModelEquity			= dynamic_cast< ARM_PricingFunctionEquity* >(mod);
    if(!itsModelEquity)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Equity or Fx compatible model required for SPOT keyword");

	ARM_Date expiryDate   = arg[1].GetDate();
	SetSettlementDate( nodes, expiryDate, itsModelEquity->GetSettlementCalendar(itsModelName), itsModelEquity->GetSettlementGap(itsModelName), arg, 5, "ARM_GP_EqDigitalFctor::SetDefaults" );
	SetSettlementDate( nodes, expiryDate, itsModelEquity->GetSettlementCalendar(itsModelName), itsModelEquity->GetSettlementGap(itsModelName), arg, 6, "ARM_GP_EqDigitalFctor::SetDefaults" );
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_EqDigitalFctor
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_EqDigitalFctor::GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		SetDefaults( arg, mod, nodes );

		/// handle vectorial strike
		if( arg[2].GetType() ==  GFAT_DOUBLE_TYPE )
		{
			itsStrikeDouble = arg[2].GetDouble();
			itsStrikeVector	= ARM_VectorPtr(NULL);
		}
		else if( arg[2].GetType() ==  GFAT_VECTOR_TYPE )
		{
			const double UNASSIGNED = -1111;
			itsStrikeDouble =  UNASSIGNED;
			itsStrikeVector	= arg[2].GetVector();
		}
		else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "strike should be vector or double" );
		/// handle epsilon
		if( arg[8].GetType() ==  GFAT_DOUBLE_TYPE )
		{
			itsEpsilon = arg[8].GetDouble();
		}

		/// grab inputs
		itsEvalTime			    = GetAndCheckEvalTime( evalDate, mod, ARM_GP_EqDigitalFctor::itsFuncName );
		string callPutString    = arg[3].GetString();
		itsNotional				= arg[4].GetDouble();
		itsCallPut			    = ARM_ArgConv_CallPut.GetNumber( callPutString );
		ARM_Date expiryDate     = arg[1].GetDate();
		itsExpiryTime			= mod->GetTimeFromDate(expiryDate);
		ARM_Date settlementDate = arg[5].GetDate();
		itsSettlementTime	    = mod->GetTimeFromDate(settlementDate);
		ARM_Date payDate        = arg[6].GetDate();
		itsPayTime	            = mod->GetTimeFromDate(payDate);
		itsDigitType			= (ARM_DigitType)ARM_ArgConv_DigitType.GetNumber( arg[7].GetString());

		/// Validation
		CheckNbSmaller( itsEvalTime, itsExpiryTime,	"EvalTime", "ExpiryTime",  ARM_GP_EqDigitalFctor::itsFuncName );
		CheckNbSmaller( itsExpiryTime, itsSettlementTime,	"ExpiryTime", "SettlementTime",  ARM_GP_EqDigitalFctor::itsFuncName );
		CheckNbSmaller( itsExpiryTime, itsPayTime,	"ExpiryTime", "PayTime",  ARM_GP_EqDigitalFctor::itsFuncName );

		SetAlreadyComputed(true);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_EqDigitalFctor
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the Equity digital function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_EqDigitalFctor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	
	ARM_VectorPtr digitalValue;
	
	if( itsStrikeVector	 == ARM_VectorPtr(NULL) )
	{
		digitalValue =  itsModelEquity->DigitalScalar( itsModelName, itsEvalTime, itsExpiryTime, itsSettlementTime, itsStrikeDouble, itsNotional, itsCallPut, itsPayTime, itsDigitType, itsEpsilon, states );
	}
	else
	{
		digitalValue =  itsModelEquity->DigitalVectorial( itsModelName, itsEvalTime, itsExpiryTime, itsSettlementTime, *itsStrikeVector, itsNotional, itsCallPut, itsPayTime, itsDigitType, itsEpsilon, states );
	}

	/// return result
	return ARM_GramFctorArg(digitalValue);

}


////////////////////////////////////////////////////
///	Class  : ARM_GP_EqDigitalFctor
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_EqDigitalFctor::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );

	/// return result
	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( 1, itsExpiryTime)), ARM_AdditionalTimeInfoPtr(NULL) );
}



////////////////////////////////////////////////////
/////////	Call Option Strip     //////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_GP_CallStripFctor
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_CallStripFctor::SetDefaults( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 14, ARM_GP_CallStripFctor::itsFuncName);
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,			ARM_GP_CallStripFctor::itsFuncName);
	GPAF_CheckArgType( arg[1], GFAT_DATE_TYPE,				ARM_GP_CallStripFctor::itsFuncName);
	GPAF_CheckArgType( arg[2], GFAT_DATEORMATU_TYPE,        ARM_GP_CallStripFctor::itsFuncName );
	GPAF_CheckArgType( arg[4], GFAT_STRING_TYPE,			ARM_GP_CallStripFctor::itsFuncName);
	GPAF_CheckArgType( arg[5], GFAT_STRING_TYPE,			ARM_GP_CallStripFctor::itsFuncName);
	GPAF_CheckArgType( arg[11],GFAT_DATESTRIP_TYPE,		    ARM_GP_CallStripFctor::itsFuncName);

	itsModelName			= arg[0].GetString();
	itsModelEquity			= dynamic_cast< ARM_PricingFunctionEquity* >(mod);
    if(!itsModelEquity)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Equity or Fx compatible model required for SPOT keyword");

	SetEndDatePayCal(nodes,arg,2,mod,itsModelName);

	SetResetDaysGap(nodes,arg,6,mod,itsModelName );
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_CallStripFctor
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_CallStripFctor::GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes  )
{
	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		SetDefaults( arg, mod, nodes );

        size_t dateStripOffset = static_cast<size_t>(arg[12].GetDouble());
		int dateStripWidth = static_cast<int>(arg[13].GetDouble());

        ARM_DateStripPtr sched(NULL);
        if(arg[11].GetType() == GFAT_DATESTRIP_TYPE)
            sched = arg[11].GetDateStrip();

        double resetGap = GetResetDaysGap(arg[6],mod,itsModelName);
        double paymentGap = arg[7].GetDouble();
        if(paymentGap < resetGap)
            paymentGap = resetGap;


        if(sched == ARM_DateStripPtr(NULL))
        {
			/*ARM_Currency* ccy	= mod->GetCurrency( itsModelName );
		    char* ccyName		= ccy->GetCcyName();
		    ARM_INDEX_TYPE defaultIndex = GetDefaultIndexFromCurrency( ccyName );
		    char* payCalendar	= ccy->GetPayCalName(defaultIndex);
		    char* resetCalendar = ccy->GetResetCalName(defaultIndex);

            /// Build the option strip schedule
		    ARM_Date startDate  = arg[1].GetDate();
		    ARM_Date endDate    = arg[2].GetDate();
		    int freq		    = GetFrequency(arg[5], mod, itsModelName, ARM_GP_CallStripFctor::itsFuncName, GetFixedFrequencyFtor(ccy));
            int dayCount        = KNOBASE;
		    sched = ARM_DateStripPtr( new ARM_DateStrip(startDate, endDate, freq, dayCount, resetCalendar,
			    K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, resetGap, freq, paymentGap, payCalendar) );

	        /// delete char* for memory leak
	        delete payCalendar;
	        delete resetCalendar;*/
        }

        size_t i,nbOptions = sched->GetResetDates()->size();
        if(dateStripOffset >= nbOptions)
            dateStripOffset = nbOptions-1;
        int nbActOptions=nbOptions-dateStripOffset;
		if (dateStripWidth != -1)
			nbActOptions = dateStripWidth < nbActOptions ? dateStripWidth : nbActOptions;


        /// Fill arguments from the dateStrip, at least expiries are input
		itsExpiryTimes.resize(nbActOptions);
		for(i=0;i<nbActOptions;++i)
			itsExpiryTimes[i] = mod->GetTimeFromDate((*(sched->GetResetDates()))[i+dateStripOffset]);

        int setllementGap   = -resetGap;
	    ARM_PricingFunctionEquity* itsModelEquity = dynamic_cast< ARM_PricingFunctionEquity* >(mod);
        const string settlementCalendar(itsModelEquity->GetSettlementCalendar(itsModelName));

		itsSettlementTimes.resize(nbActOptions);
        if(sched->GetFlowStartDates())
        {            
		    for(i=0;i<nbActOptions;++i)
			    itsSettlementTimes[i] = mod->GetTimeFromDate((*(sched->GetFlowStartDates()))[i+dateStripOffset]);
        }
        else
        {
            /// Build settlement dates w.r.t. market conventions
            for(i=0;i<nbActOptions;++i)
		    {
                ARM_Date settlementDate((*(sched->GetResetDates()))[i+dateStripOffset]);
                settlementDate.GapBusinessDay( setllementGap, const_cast<char*>(settlementCalendar.c_str()) ) ;
                itsSettlementTimes[i] = mod->GetTimeFromDate(settlementDate.GetJulian());
            }
        }

		itsPayTimes.resize(nbActOptions);
        if(sched->GetPaymentDates())
        {
		    for(i=0;i<nbActOptions;++i)
			    itsPayTimes[i] = mod->GetTimeFromDate((*(sched->GetPaymentDates()))[i+dateStripOffset]);
        }
        else
        {
            /// Build payment dates w.r.t.input payment gap between settle & pay
		    for(i=0;i<nbActOptions;++i)
            {
                if(paymentGap != 0.0)
                {
                    ARM_Date paymentDate((*(sched->GetResetDates()))[i+dateStripOffset]);
                    paymentDate.GapBusinessDay( setllementGap, const_cast<char*>(settlementCalendar.c_str()) ) ;
                    paymentDate.GapBusinessDay( paymentGap, const_cast<char*>(settlementCalendar.c_str()) ) ;
                    itsPayTimes[i] = mod->GetTimeFromDate(paymentDate.GetJulian());
                }
                else
                    itsPayTimes[i] = itsSettlementTimes[i];
            }
        }


        /// Fill strikes and leverages (= initial leverages * nominals)
        itsStrikes  = GetProfileValues(arg[3],dateStripOffset,itsExpiryTimes,"strike");
        itsNominals = GetProfileValues(arg[8],dateStripOffset,itsExpiryTimes,"nominal");

        ARM_VectorPtr leverages = GetProfileValues(arg[10],dateStripOffset,itsExpiryTimes,"leverage");
        for(i=0;i<nbActOptions;++i)
            (*itsNominals)[i] *= (*leverages)[i];


		/// Fill evalTime & C/P
		itsEvalTime			    = GetAndCheckEvalTime( evalDate, mod, ARM_GP_CallStripFctor::itsFuncName );
		string callPutString    = arg[4].GetString();
		itsCallPut			    = ARM_ArgConv_CallPut.GetNumber( callPutString );


		SetAlreadyComputed(true);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_CallStripFctor
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the call function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_CallStripFctor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	
	ARM_VectorPtr callValue;

    /// Only scalar (deterministic strike profile) version is used at the moment
	callValue =  itsModelEquity->CallStripScalar( itsModelName, itsEvalTime, itsExpiryTimes, itsSettlementTimes, itsStrikes->GetValues(), itsCallPut, itsPayTimes, itsNominals->GetValues(), states );

	/// return result
	return ARM_GramFctorArg(callValue);
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_CallStripFctor
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_CallStripFctor::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );

	/// return result
	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( itsExpiryTimes )) ,ARM_AdditionalTimeInfoPtr(NULL) );
}


////////////////////////////////////////////////////
///////////		 RangeAccrual Option      //////////////////
////////////////////////////////////////////////////

// Parameters
// 0  Curve						string
// 1  StartDate					date
// 2  EndDate					date
// 3  PaymentGap				double
// 4  PaymentDayCount			string
// 5  PAYindexType				string
// 6  PAYindexTerm				string
// 7  PAYindexTiming			string
// 8  PAYindexSpread			double
// 9  FixingFrequency			string
// 10 FXccy						string
// 11 FXdownBarrier				double
// 12 FXupBarrier				double
// 13 Notional					double
// 14 IRindexType				string
// 15 IRindexFixingTerm			string
// 16 IRdownBarrier				double
// 17 IRupBarrier				double

//
// Cash flow is :
// (PayIndex + Spread) * Notional *	PeriodRation * Indic{FXBarrierDown <= FX <= FXBarrierUp} * Indic{IRBarrierDown <= IEIndex <= IRBarrierUp}
//

////////////////////////////////////////////////////
///	Class  : ARM_GP_RangeAccrualFctor
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_RangeAccrualFctor::SetDefaults( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 18, ARM_GP_RangeAccrualFctor::itsFuncName);
	GPAF_CheckArgType( arg[0],  GFAT_STRING_TYPE,			ARM_GP_RangeAccrualFctor::itsFuncName );//Curve name
	GPAF_CheckArgType( arg[1],  GFAT_DATE_TYPE,				ARM_GP_RangeAccrualFctor::itsFuncName );//StartDate
	GPAF_CheckArgType( arg[2],  GFAT_DATE_TYPE,				ARM_GP_RangeAccrualFctor::itsFuncName );//EndDate
	GPAF_CheckArgType( arg[3],  GFAT_DATE_OR_DOUBLE_TYPE,	ARM_GP_RangeAccrualFctor::itsFuncName );//PaymentGap
	GPAF_CheckArgType( arg[4],  GFAT_STRING_TYPE,			ARM_GP_RangeAccrualFctor::itsFuncName );//PaymentDayCount
	GPAF_CheckArgType( arg[5],  GFAT_STRING_TYPE,			ARM_GP_RangeAccrualFctor::itsFuncName );//PayIndexType
	GPAF_CheckArgType( arg[6],  GFAT_STRING_TYPE,			ARM_GP_RangeAccrualFctor::itsFuncName );//PayIndexTerm
	GPAF_CheckArgType( arg[7],  GFAT_STRING_TYPE,			ARM_GP_RangeAccrualFctor::itsFuncName );//PayIndexTiming
	GPAF_CheckArgType( arg[8],  GFAT_DOUBLE_TYPE,			ARM_GP_RangeAccrualFctor::itsFuncName );//PayIndexSpread
	GPAF_CheckArgType( arg[9],  GFAT_STRING_TYPE,			ARM_GP_RangeAccrualFctor::itsFuncName );//FixingFrequency
	GPAF_CheckArgType( arg[10], GFAT_STRING_TYPE,			ARM_GP_RangeAccrualFctor::itsFuncName );//FXccy
	GPAF_CheckArgType( arg[11], GFAT_DOUBLE_TYPE,			ARM_GP_RangeAccrualFctor::itsFuncName );//FXDownBarrier
	GPAF_CheckArgType( arg[12], GFAT_DOUBLE_TYPE,			ARM_GP_RangeAccrualFctor::itsFuncName );//FXUpBarrier
	GPAF_CheckArgType( arg[13], GFAT_DOUBLE_TYPE,			ARM_GP_RangeAccrualFctor::itsFuncName );//Notional
	GPAF_CheckArgType( arg[14], GFAT_STRING_TYPE,			ARM_GP_RangeAccrualFctor::itsFuncName );//IRindexType
	GPAF_CheckArgType( arg[15], GFAT_STRING_TYPE,			ARM_GP_RangeAccrualFctor::itsFuncName );//IRindexFixingTerm
	GPAF_CheckArgType( arg[16], GFAT_DOUBLE_TYPE,			ARM_GP_RangeAccrualFctor::itsFuncName );//IRDownBarrier
	GPAF_CheckArgType( arg[17], GFAT_DOUBLE_TYPE,			ARM_GP_RangeAccrualFctor::itsFuncName );//IRUpBarrier

	// curve name
	itsCurveName		= arg[0].GetString();
	// fxName
	itsFXname			= arg[10].GetString();

	//Validation of model
	itsModelEquity			= dynamic_cast< ARM_PricingFunctionEquity* >(mod);
    if(!itsModelEquity)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Equity or Fx compatible model required for RangeAccrual");
	
	//payment date (by convention equals to the settlement date of the first FX for the moment)
	ARM_FXName  fxname(itsFXname.substr(itsFXname.size()-6));   
	ARM_Date endDate		= arg[2].GetDate();
	SetSettlementDate( nodes, endDate, itsModelEquity->GetSettlementCalendar(fxname.GetMktName()), itsModelEquity->GetSettlementGap(fxname.GetMktName()), arg, 3, "ARM_GP_RangeAccrualFctor::SetDefaults" );
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_RangeAccrualFctor
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_RangeAccrualFctor::GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes  ){}
//{
//	/// if it is not already computed
//	if(!GetAlreadyComputed())
//	{
//		SetDefaults( arg, mod, nodes );
//
//// PAY INDEX
//
//		itsPAYindexType				= ARM_ArgConv_IndexClass.GetNumber(arg[5].GetString());
//		string payIndexTerm;
//
//		ARM_Currency* ccy			= mod->GetCurrency( itsCurveName );
//		char* ccyName				= ccy->GetCcyName();
//
//		itsPAYindexTerm				= GetTerm(arg[6], mod, itsCurveName, ARM_GP_RangeAccrualFctor::itsFuncName, GetFloatFrequencyFtor(ccy));
//		int payIndexTiming			= ARM_ArgConv_Timing.GetNumber(arg[7].GetString());
//
//// IR INDEX
//
//		itsIRindexType				= ARM_ArgConv_IndexClass.GetNumber(arg[14].GetString());
//		double fixIndexTerm			= GetTerm(arg[15], mod, itsCurveName, ARM_GP_RangeAccrualFctor::itsFuncName, GetFloatFrequencyFtor(ccy));
//		string fixIndexTermStr		= ConvertYearTermToStringMatu(fixIndexTerm);
//
//		ARM_INDEX_TYPE fixIndexType;
//		if (itsIRindexType == K_LIBOR)
//			fixIndexType = (ARM_INDEX_TYPE)FromIndexAndTermToIndexType(fixIndexTermStr, ccyName);
//		else
//			fixIndexType = EURIBOR1Y;
//		ARM_IRIndex fixIndex(fixIndexType);
//
//		int payDayCount				= GetDayCount( arg[4], mod, itsCurveName, GetFloatDayCountFtor(ccy));
//		
//		/// handle PAYindexSpread
//		if( arg[8].GetType() ==  GFAT_DOUBLE_TYPE )
//			itsPAYindexSpread = arg[8].GetDouble();
//		else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "alpha should be a double" );
//		/// handle vectorial FXdownBarrier
//		if( arg[11].GetType() ==  GFAT_DOUBLE_TYPE )
//			itsFXDownBarrierVector	= ARM_VectorPtr( new std::vector<double>( 1,arg[11].GetDouble() ) );
//		else if( arg[11].GetType() ==  GFAT_VECTOR_TYPE )
//			itsFXDownBarrierVector	= arg[11].GetVector();
//		else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "FXdownBarrier should be vector or double" );
//		/// handle vectorial FXupBarrier
//		if( arg[12].GetType() ==  GFAT_DOUBLE_TYPE )
//			itsFXUpBarrierVector	= ARM_VectorPtr( new std::vector<double>( 1,arg[12].GetDouble() ) );
//		else if( arg[12].GetType() ==  GFAT_VECTOR_TYPE )
//			itsFXUpBarrierVector	= arg[12].GetVector();
//		else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "FXupBarrier should be vector or double" );
//		/// handle vectorial IRdownBarrier
//		if( arg[16].GetType() ==  GFAT_DOUBLE_TYPE )
//			itsIRDownBarrierVector	= ARM_VectorPtr( new std::vector<double>( 1,arg[16].GetDouble() ) );
//		else if( arg[16].GetType() ==  GFAT_VECTOR_TYPE )
//			itsIRDownBarrierVector	= arg[16].GetVector();
//		else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "IRdownBarrier should be vector or double" );
//		/// handle vectorial IRupBarrier
//		if( arg[17].GetType() ==  GFAT_DOUBLE_TYPE )
//			itsIRUpBarrierVector	= ARM_VectorPtr( new std::vector<double>( 1,arg[17].GetDouble() ) );
//		else if( arg[17].GetType() ==  GFAT_VECTOR_TYPE )
//			itsIRUpBarrierVector	= arg[17].GetVector();
//		else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "IRupBarrier should be vector or double" );
//		/// handle notional
//		if( arg[13].GetType() ==  GFAT_DOUBLE_TYPE )
//			itsNotionalVector	= ARM_VectorPtr( new std::vector<double>( 1,arg[13].GetDouble() ) );
//		else if( arg[13].GetType() ==  GFAT_VECTOR_TYPE )
//			itsNotionalVector	= arg[13].GetVector();
//		else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "notional should be vector or double" );
//
//		/// handle dates
//		itsEvalTime			    = GetAndCheckEvalTime( evalDate, mod, ARM_GP_RangeAccrualFctor::itsFuncName );
//		ARM_Date startDate		= arg[1].GetDate();
//		itsStartTime			= mod->GetTimeFromDate(startDate);
//		ARM_Date endDate		= arg[2].GetDate();
//		itsEndTime				= mod->GetTimeFromDate(endDate);
//		ARM_Date payDate        = arg[3].GetDate();
//		itsPAYtime	            = mod->GetTimeFromDate(payDate);
//
//        /// Creation of the intermediate datestrip for fixingDates
//        ARM_DateStripPtr sched(NULL);
//		ARM_INDEX_TYPE defaultIndex = GetDefaultIndexFromCurrency( ccyName );
//		char* payCalendar	= ccy->GetPayCalName(defaultIndex);
//		char* resetCalendar = ccy->GetResetCalName(defaultIndex);
//		int freq		    = GetFrequency(arg[9], mod, itsFXname, ARM_GP_RangeAccrualFctor::itsFuncName, GetFixedFrequencyFtor(ccy));
//        int dayCount        = payDayCount;
//		double resetGap=0;
//		double paymentGap=0;
//		sched = ARM_DateStripPtr( new ARM_DateStrip(startDate, endDate, freq, dayCount, resetCalendar,
//			K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTEND, resetGap, freq, paymentGap, payCalendar,K_ARREARS) );
//		
//		/// handle fixingTimes, IRindexStartTimes, IRindexEndTimes and IRindexTerms
//		size_t i,nbFixing = sched->GetResetDates()->size();
//		itsFixingTimes.resize(nbFixing);
//		itsIRindexResetTimes.resize(nbFixing);
//		itsIRindexStartTimes.resize(nbFixing);
//		itsIRindexEndTimes.resize(nbFixing);
//		itsIRindexTerms.resize(nbFixing);
//
//		for(i=0;i<nbFixing;++i)
//		{
//			//Fixing Times
//			itsFixingTimes[i] = mod->GetTimeFromDate((*(sched->GetResetDates()))[i]);
//			//IRindexStartTimes
//			ARM_Date fwdStart((*(sched->GetResetDates()))[i]);//I use the reset date?
//			itsIRindexStartTimes[i] = itsFixingTimes[i];
//			//IRindexEndTimes
//			ARM_Date fwdEnd((*(sched->GetResetDates()))[i]);//I use the reset date?
//		    //fwdEnd.AddPeriod(fixIndexType);
//			fwdEnd.AddPeriod(fixIndexTermStr);
//		    fwdEnd.GoodBusinessDay(fixIndex.GetFwdRule() * fixIndex.GetIntRule(), resetCalendar); 
//		    itsIRindexEndTimes[i] = mod->GetTimeFromDate(fwdEnd);
//			//IRindexResetTimes
//			itsIRindexResetTimes[i] = itsIRindexEndTimes[i];//for the moment ir index reset in arrears and at EndTimes
//			//IRindexTerms
//			itsIRindexTerms[i] = CountYearsWithoutException( dayCount, fwdStart, fwdEnd );//CountBusinessDays(fwdStart, fwdEnd, payCalendar);
//		}
//
//		/// delete char* for memory leak
//	    delete payCalendar;
//	    delete resetCalendar;
//
//
//		/// Validation
//		CheckNbSmaller( itsEvalTime, itsStartTime, "EvalTime", "StartTime",  ARM_GP_RangeAccrualFctor::itsFuncName );
//		CheckNbSmaller( itsStartTime, itsEndTime, "StartTime", "EndTime",  ARM_GP_RangeAccrualFctor::itsFuncName );
//		CheckNbSmaller( itsEndTime, itsPAYtime,	"EndTime", "PayTime",  ARM_GP_RangeAccrualFctor::itsFuncName );
//
//		SetAlreadyComputed(true);
//	}
//}


////////////////////////////////////////////////////
///	Class  : ARM_GP_RangeAccrualFctor
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the range accrual function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_RangeAccrualFctor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );
	/// init of call vectorial
	size_t payoffSize = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;
	ARM_VectorPtr rangeAccrualValue;
	rangeAccrualValue =  itsModelEquity->RangeAccrualVectorial( 
		itsCurveName,
		itsEvalTime,
		itsStartTime,
		itsEndTime,
		itsPAYtime,
		itsFixingTimes,
		itsPAYindexType, 
		itsPAYindexTerm,
		itsFXname,
		itsIRindexType,
		itsIRindexResetTimes,
        itsIRindexStartTimes,
        itsIRindexEndTimes,
        itsIRindexTerms,
		*itsFXDownBarrierVector, 
		*itsFXUpBarrierVector, 
		*itsIRDownBarrierVector, 
		*itsIRUpBarrierVector, 
		*itsNotionalVector, 
		states );
	/// return result
	return ARM_GramFctorArg(rangeAccrualValue);
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_RangeAccrualFctor
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_RangeAccrualFctor::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );
	/// return result
	return ARM_NodeInfo( ARM_GP_VectorPtr( new ARM_GP_Vector(itsFixingTimes) ),ARM_AdditionalTimeInfoPtr(NULL) );
}


#define ARM_CF_EPS 1.0e-13

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

