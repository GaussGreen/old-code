#pragma warning(disable : 4786)

#include "firstToBeIncluded.h"
#include "arm_local_xxxproject.h"

/// ARM Kernel
#include <inst/security.h>
#include <gpinfra/gramfunctorargdict.h>
#include <gpinfra\mktdatamanagerrep.h>

#include <gpbase/gpvector.h>
#include <gpbase/countedptr.h>
#include <gpbase/cloneutilityfunc.h>

#include <GP_Help\gphelp\crmcookies.h>

#include <xxxproject/mktdatas.h>
#include <xxxproject/vanillapricer.h>
#include <xxxproject/pricerbuilder.h>
#include <xxxproject/scenario.h>
#include <xxxproject/argconvdefault.h>
#include <xxxproject/hedge.h>

#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_local_wrapper.h"

#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_persistent.h>



#include <string>
CC_USING_NS(std,string)

#include <set>
CC_USING_NS(std,set)

/// Objects are in namespace, hence the using directive!
using ARM::ARM_MktData;
using ARM::ARM_YC_Subject;
using ARM::ARM_BS_Subject;
using ARM::ARM_CAP_Subject;
using ARM::ARM_OSW_Subject;
using ARM::ARM_FX_Subject;
using ARM::ARM_VFX_Subject;
using ARM::ARM_MIX_Subject;
using ARM::ARM_COR_Subject;
using ARM::ARM_DAT_Subject;
using ARM::ARM_SO_Subject;
using ARM::ARM_INF_Subject;
using ARM::ARM_IBS_Subject;

using ARM::ARM_GP_T_Vector;
using ARM::ARM_VanillaPricer;
using ARM::ARM_PricerBuilder;
using ARM::ARM_Pricer;
using ARM::ARM_Scenario;
using ARM::ARM_0D_Scenario;
using ARM::ARM_1D_Scenario;
using ARM::ARM_2D_Scenario;
using ARM::ARM_ND_Scenario;
using ARM::ARM_Hedge;
using ARM::CreateClone;

using ARM::ARM_GramFctorArgDict;
using ARM::ARM_MarketData_ManagerRep;
using ARM::ARM_CRMCookies;
//using ARM::ARM_GramFctorArg;

extern long ARMLOCAL_MktDatas_Create(
	const double &		  asOfDate,
	const vector<string>& keys,
	const vector<long>& mktDatas,
	ARM_result& result,
	long objId)
		{
		/// input checks
		if( !GlobalPersistanceOk( result ) )
			return ARM_KO;

		/// used in the MACRO ARM_RESULT
		CCString msg ("");
		ARM_MktData* theMktData= NULL;

		try
		{
			if (keys.size() != mktDatas.size())
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_INPUT,
							"The keys and the market data should have the same size.");
			}

			vector<ARM_Object*> mktDatasObj(mktDatas.size());

			for (int i = 0; i < keys.size(); ++i)
				mktDatasObj[i] = LOCAL_PERSISTENT_OBJECTS->GetObject(mktDatas[i]);

			char * startDate = new char[11];
			Local_XLDATE2ARMDATE(asOfDate,startDate);

			theMktData = new ARM_MktData((ARM_Date) startDate, keys,mktDatasObj);

			// assign object
			if( !assignObject( theMktData, result, objId ) ){
				return ARM_KO; }
			else{
				return ARM_OK; }
		}
	
		catch(Exception& x)
		{
			delete theMktData;
			x.DebugPrint();
			ARM_RESULT();
		}
	}


extern long ARMLOCAL_ARM_XXX_Price (long secId, long mktDt, ARM_result& result)
{
	ARM_Object* sec			= NULL;
	const ARM_MktData* mkt		= NULL;

	double price		= 0.0;
	CCString msg ("");

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 ){
		result.setMsg("ARM_ERR: Pb with accessing objects");
		return(ARM_KO);
	}

	try
	{
		sec =	dynamic_cast<ARM_Object*>(LOCAL_PERSISTENT_OBJECTS->GetObject(secId));
		mkt =	dynamic_cast<ARM_MktData* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDt));
	
		string						tmpCurrency("EUR");

		if (!sec)
		{
			result.setMsg("ARM_ERR: Sec Id is not a security.");
			return(ARM_KO);
		}

		if (!mkt)
		{
			result.setMsg("ARM_ERR: MktDatas Id is not a MktDatas.");
			return(ARM_KO);
		}



		ARM_Pricer* pricer = ARM_PricerBuilder::BuildPricer(sec);
		pricer->SetMkt( (ARM_MktData * ) mkt ) ;

		double price = pricer->Price( ); 

		// save the result!
		result.setDouble(price);
		return ARM_OK;
	}

	// first catch ARM type Exception
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg("ARM_ERR: unrecognized failure in ARMLOCAL_ARM_Price");
        return(ARM_KO);
	}
}

extern long ARMLOCAL_ARM_XXX_Scenario(	const double&		shift,
										const string&		currency,
										const string&		typeScenario,
										const string&		subTypeScenario,
										const string&		stressOrder,
										const long&			relatif,
										const long&			cumulInv,
										const long&			perturbative,
										ARM_result&			result,
										long				objId)	{

		if( !GlobalPersistanceOk( result ) )		return ARM_KO;
		CCString msg ("");
		int volType;


		ARM_Scenario* theScenario = NULL;

		try
		{
			if		(	typeScenario		== string("DELTA_YC")	)
				theScenario	= new  ARM_1D_Scenario< ARM_YC_Subject >(shift, currency,  relatif,  cumulInv, perturbative );

			else if (	typeScenario		==	string("DELTA_BS")	)
				theScenario	=new ARM_1D_Scenario< ARM_BS_Subject >	(shift, currency,  relatif,  cumulInv, perturbative );
			
			else if (	typeScenario		==	string("DELTA_FX")	)
				theScenario	= new  ARM_0D_Scenario< ARM_FX_Subject >(shift, currency,  relatif,  cumulInv, perturbative );

			else if (	typeScenario		==	string("CEGA")		) 
				theScenario = new  ARM_2D_Scenario< ARM_COR_Subject >(shift, currency,  relatif,  cumulInv, perturbative );

			
			else if (	typeScenario		==	string("VEGA_CAP")	) {
				if	(	subTypeScenario		==	string("")	)
					ARMTHROW(ERR_INVALID_ARGUMENT,"The SubType Scenario has to be specified.");
				
				volType	=	ARM::ARM_ArgConv_MktVolType.GetNumber(subTypeScenario);
				switch( volType ){
					case ATM:
						theScenario	=	new  ARM_2D_Scenario<ARM_CAP_Subject<ATM> >(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					case RHO:
						theScenario	=	new  ARM_2D_Scenario<ARM_CAP_Subject<RHO> >(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					case NU:	
						theScenario	=	new  ARM_2D_Scenario<ARM_CAP_Subject<NU> >(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					case BETA:
						theScenario	=	new  ARM_2D_Scenario<ARM_CAP_Subject<BETA> >(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					}
				}
			else if (	typeScenario		==	string("VEGA_OSW")	) {
				if	(	subTypeScenario		==	string("")	)
					ARMTHROW(ERR_INVALID_ARGUMENT,"The SubType Scenario has to be specified.");
				
				volType	=	ARM::ARM_ArgConv_MktVolType.GetNumber(subTypeScenario);
				switch( volType ){
					case ATM:
						theScenario	=	new  ARM_2D_Scenario<ARM_OSW_Subject<ATM> >(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					case RHO:
						theScenario	=	new  ARM_2D_Scenario<ARM_OSW_Subject<RHO> >(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					case NU:	
						theScenario	=	new  ARM_2D_Scenario<ARM_OSW_Subject<NU> >(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					case BETA:
						theScenario	=	new  ARM_2D_Scenario<ARM_OSW_Subject<BETA> >(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					}			
				}
			else if (	typeScenario		==	string("VEGA_VFX")	) {
				if	(	subTypeScenario		==	string("")	)
					ARMTHROW(ERR_INVALID_ARGUMENT,"Implemented but commented.");
				
				volType	=	ARM::ARM_ArgConv_MktVolType.GetNumber(subTypeScenario);
				switch( volType ){
					case PIV:
						theScenario	=	new  ARM_2D_Scenario<ARM_VFX_Subject<PIV> >(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					case RR:
						theScenario	=	new  ARM_2D_Scenario<ARM_VFX_Subject<RR> >(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					case STR:	
						theScenario	=	new  ARM_2D_Scenario<ARM_VFX_Subject<STR> >(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					}		
				}
			else if (	typeScenario	==	string("VEGA_MIX")	){
				if	(	subTypeScenario	==	string("")	)
					ARMTHROW(ERR_INVALID_ARGUMENT,"The SubType Scenario has to be specified.");

				volType	=	ARM::ARM_ArgConv_MktVolType.GetNumber(subTypeScenario);
				switch( volType ){
					case VOL:
						theScenario	=	new  ARM_2D_Scenario<ARM_MIX_Subject<VOL> >(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					case SMILE:
						theScenario	=	new  ARM_2D_Scenario<ARM_MIX_Subject<SMILE> >(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					case SHIFT:	
						theScenario	=	new  ARM_2D_Scenario<ARM_MIX_Subject<SHIFT> >(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					case Q:	
						theScenario	=	new  ARM_2D_Scenario<ARM_MIX_Subject<Q	>	>(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					}
				}	
			else if (	typeScenario	==	string("VEGA_SO")	){
				if	(	subTypeScenario	==	string("")	)
					ARMTHROW(ERR_INVALID_ARGUMENT,"The SubType Scenario has to be specified.");

				volType	=	ARM::ARM_ArgConv_MktVolType.GetNumber(subTypeScenario);
				switch( volType ){
					case ADJ:
						theScenario	=	new  ARM_2D_Scenario<ARM_SO_Subject<ADJ> >(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					case ATM:
						theScenario	=	new  ARM_2D_Scenario<ARM_SO_Subject<ATM> >(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					}			
				}
			else if (	typeScenario	==	string("THETA")	) {
				char * shiftDate = new char[11];
				Local_XLDATE2ARMDATE(shift,shiftDate);
				theScenario= new  ARM_0D_Scenario< ARM_DAT_Subject >(((ARM_Date) shiftDate).GetJulian(), currency,  relatif,  cumulInv, perturbative );
			}
			
			else if (	typeScenario	==	string("DELTA_INF")	){
				theScenario	= new  ARM_1D_Scenario< ARM_INF_Subject >(shift, currency,  relatif,  cumulInv, perturbative );
			
			}
			
			else if (	typeScenario	==	string("VEGA_INF")	){
				if	(	subTypeScenario	==	string("")	)
					ARMTHROW(ERR_INVALID_ARGUMENT,"The SubType Scenario has to be specified.");

				volType	=	ARM::ARM_ArgConv_MktVolType.GetNumber(subTypeScenario);
				switch( volType ){
					case CPI:
						theScenario	=	new  ARM_2D_Scenario<ARM_IBS_Subject<CPI> >(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					case YOY:
						theScenario	=	new  ARM_2D_Scenario<ARM_IBS_Subject<YOY> >(shift, currency,  relatif,  cumulInv, perturbative );
						break;
					}			
				}		

			theScenario->InitScenario(stressOrder);

			if( !assignObject( theScenario, result, objId ) )
				return ARM_KO; 
			else
				return ARM_OK; 
			
		}

		catch(Exception& x)
		{
			delete theScenario;
			x.DebugPrint();
			ARM_RESULT();
		}
}

extern long ARMLOCAL_Hedge_Create (		const long&			secId, 
										const long&			scenId,
										const long&			mktDt,
 										ARM_result&			result, 
										long				objId)
{
	ARM_Object*			sec			=	NULL;
	ARM_Scenario*		sce			=	NULL;
	ARM_MktData*		mkt			=	NULL;
	ARM_Hedge*			theHedge	=	NULL;

	CCString msg ("");

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 ){
		result.setMsg("ARM_ERR: Pb with accessing objects");
		return(ARM_KO);
	}

	try
	{
		sec =	dynamic_cast<ARM_Object*>(LOCAL_PERSISTENT_OBJECTS->GetObject(secId));
		if (!sec)
		{
			result.setMsg("ARM_ERR: MktDatas Id is not a MktDatas.");
			return(ARM_KO);
		}

		sce =	dynamic_cast<ARM_Scenario*>(LOCAL_PERSISTENT_OBJECTS->GetObject(scenId));
		if (!sce)
		{
			result.setMsg("ARM_ERR: MktDatas Id is not a MktDatas.");
			return(ARM_KO);
		}

		mkt =	dynamic_cast<ARM_MktData*>(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDt));
		if (!mkt)
		{
			result.setMsg("ARM_ERR: MktDatas Id is not a MktDatas.");
			return(ARM_KO);
		}

		theHedge	=	new ARM_Hedge(sec, mkt);
		theHedge	->	ComputeHedge(sce);

		// assign object
		if( !assignObject( theHedge, result, objId ) )
			return ARM_KO; 
		else
			return ARM_OK; 

  	}

	// first catch ARM type Exception
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg("ARM_ERR: unrecognized failure in ARMLOCAL_ARM_Price");
        return(ARM_KO);
	}
}


extern long ARMLOCAL_Scenari_Compose(	const long &		scenariId1,
										const long &		scenariId2,
										ARM_result&			result,
										long				objId){

	ARM_Scenario*		scenId1				=	NULL;
	ARM_Scenario*		scenId2				=	NULL;
	ARM_ND_Scenario*	theScenariCompose	=	NULL;

	CCString msg ("");

	if ( CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0 ){
		result.setMsg("ARM_ERR: Pb with accessing objects");
		return(ARM_KO);
	}

	try
	{
		scenId1 =	dynamic_cast<ARM_Scenario*>(LOCAL_PERSISTENT_OBJECTS->GetObject( scenariId1 ));
		if (!scenId1 )
		{
			result.setMsg("ARM_ERR: Scenario Id is not a Scenario.");
			return(ARM_KO);
		}

		scenId2 =	dynamic_cast<ARM_Scenario*>(LOCAL_PERSISTENT_OBJECTS->GetObject( scenariId2 ));
		if (!scenId2 )
		{
			result.setMsg("ARM_ERR: Scenario Id is not a Scenario.");
			return(ARM_KO);
		}

		vector<ARM_Scenario*> Scen(2);
		Scen[0] = scenId1;
		Scen[1] = scenId2;


		if ( theScenariCompose ) {	delete theScenariCompose; theScenariCompose=NULL; }
		theScenariCompose	=	new ARM_ND_Scenario( Scen );

		if( !assignObject( theScenariCompose, result, objId ) )
			return ARM_KO; 
		else
			return ARM_OK; 

  	}

	catch(Exception& x)
	{
		delete scenId1;
		delete scenId2;

		x.DebugPrint();
		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg("ARM_ERR: unrecognized failure in ARMLOCAL_Scenari_Create");
        return(ARM_KO);
	}

}

long ARMLOCAL_Hedge_GetData(	long				hedgeId,
								const string&		key,
								ARM_GramFctorArg&	argResult,
								ARM_result&			result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_Hedge*		theHedge	=	NULL;

		theHedge =	dynamic_cast<ARM_Hedge *>(LOCAL_PERSISTENT_OBJECTS->GetObject(hedgeId));
		if (!theHedge)
		{
			result.setMsg("ARM_XXX_ERR: MktDatas Id is not a MktDatas.");
			return(ARM_KO);
		}

		argResult = theHedge->GetDictionnary().GetData(key);

		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARM_GP_XXXMktDataFromMktDataMgerFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_MktData* theMktData= NULL;
	ARM_MarketData_ManagerRep* mktDataManager= NULL;
	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, " Model Params bumping" );

		//Model Pricing
		long objetId = genericParams->GetParamValue("MktDataMgerId").GetObjectId();
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &mktDataManager, objetId,"Marhet Data Manager", result ) ) return ARM_KO;

		vector<string >	keys;
		vector<ARM_Object* > mktObjects;
		ARM_Object* object =	NULL;

		ARM_MarketData_ManagerRep::const_iterator iter = mktDataManager->begin();
		for (iter = mktDataManager->begin(); iter != mktDataManager->end() ; ++iter )
		{
			keys.push_back(*iter);
			object = mktDataManager->GetData(*iter);
			mktObjects.push_back(object);
		}

		theMktData = new ARM_MktData(ARM_Date(), keys,mktObjects);
		// assign object
		return assignObject( theMktData, result, objId ) ? ARM_OK : ARM_KO; 
	}

	catch(Exception& x)
	{
		delete	theMktData;
		theMktData = NULL;

		x.DebugPrint();
		ARM_RESULT();
	}
}


long XXX_MktData_ApplyScenarioFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_MktData* theNewMktData= NULL;
	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "XXX_MktData_ApplyScenario" );

		// Mkt Data
		long mktDataId = genericParams->GetParamValue("MktDataId").GetObjectId();
		ARM_MktData* mktData =	dynamic_cast<ARM_MktData*>(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataId));
		if (!mktData)
		{
			result.setMsg("ARM_ERR: MktDataId is not a MktData.");
			return(ARM_KO);
		}

		long scenarioId = genericParams->GetParamValue("ScenarioId").GetObjectId();
		ARM_Scenario* scenario =	dynamic_cast<ARM_Scenario*>(LOCAL_PERSISTENT_OBJECTS->GetObject(scenarioId));
		if (!mktData)
		{
			result.setMsg("ARM_ERR: MktDataId is not a MktData.");
			return(ARM_KO);
		}
		
		int nbShift = genericParams->GetParamValue("NbShift").GetDouble();

		theNewMktData = mktData->CreateCopy();

		ARM_Hedge theHedge(NULL, theNewMktData);

		if (nbShift < 0)
			nbShift = scenario->GetNbShift();

		theHedge.ApplyScenario(scenario, nbShift);
	
		// assign object
		return assignObject( theNewMktData, result, objId ) ? ARM_OK : ARM_KO; 
	}

	catch(Exception& x)
	{
		delete	theNewMktData;
		theNewMktData = NULL;

		x.DebugPrint();
		ARM_RESULT();
	}
}
