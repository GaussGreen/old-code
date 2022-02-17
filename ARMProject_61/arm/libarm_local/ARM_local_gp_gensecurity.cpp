/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_gensecurity.cpp,v $
 * Revision 1.1  2003/13/09 15:08:43  ebenhamou
 * Initial version
 *
 */


/*! \file ARM_local_gp_gensecurity.cpp
 *
 *  \brief file for the generic security of the generic pricer local addins functions
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2003
 */


#include "firstToBeIncluded.h"

///////////////////////////////////////////////
/// WARNING include these headers FIRST
/// because it uses the new style headers
/// CCxl uses cstyle headers and leads to conflict
/// if defined first
///////////////////////////////////////////////


#include <GP_Base\gpbase\datestrip.h>

#include <GP_Infra\gpinfra\typedef.h>
#include <GP_Infra\gpinfra\gensecurity.h>
#include <GP_Infra\gpinfra\dealdescription.h>
#include <GP_Infra\gpinfra\modelnrefcall.h>
#include <GP_Infra\gpinfra\gramnode.h>
#include <GP_Infra\gpinfra\gramfunctorarg.h>
#include <GP_Infra\gpinfra\pricingadviser.h>
#include <GP_Infra\gpinfra\cstmanager.h>
#include <GP_Infra\gpinfra\gensecmanipulator.h>
#include <GP_Infra\gpinfra\pricingmodel.h>
#include <GP_Infra\gpinfra\genpricer.h>
#include <GP_Infra\gpinfra\dealdescriptionpreprocessor.h>


#include <GP_Help\gphelp\gramhelp.h>
#include <GP_Help\gphelp\crmcookies.h>

#include "ARM_local_gp_gensecurity.h"
#include "ARM_local_class.h"
#include "ARM_local_persistent.h"
#include "ARM_result.h"
#include "ARM_local_glob.h"
#include "ARM_local_wrapper.h"

#include "gpbase/typedef.h"




/// to debug release only crash
#include "gpbase/eventviewerfwd.h"

//// using the namespace directive to access ARM object!
using ARM::ARM_DealDescription;
using ARM::ARM_GenSecurity;
using ARM::ARM_DealDescriptionPtr;
using ARM::ARM_GramHelper;
using ARM::ARM_CRMCookies;
using ARM::ARM_CstManager;
using ARM::ARM_CstManagerPtr;
using ARM::ARM_GenSecManipulator;
using ARM::ARM_GenPricer;
using ARM::ARM_PricingModel;
using ARM::ARM_DescriptionPreprocessor;
using ARM::ARM_GramFctorArg;
using ARM::ARM_VectorPtr;
using ARM::ARM_GP_MatrixPtr;
using ARM::ARM_GP_CurvePtr;
using ARM::ARM_GP_Vector;
using ARM::ARM_GP_Matrix;
using ARM::ARM_Curve;
using ARM::ARM_DateStrip;
using ARM::ARM_DateStripPtr;


/// to debug release only crash
using ARM::ARM_EventViewerImp;
using ARM::ARM_TheEventViewer;


////////////////////////////////////////////
//// Function to create a deal description
////////////////////////////////////////////
extern long ARMLOCAL_DealDes_Create(
	const VECTOR<string>&	vecString,
	const VECTOR<ARM_GP_VALUE_TYPE>& vecTypes,
	long rowsNb,
	long colsNb,
    const string& ,	/// not used!
	long ,			/// not used!
	bool ,			/// not used!
	bool ,			/// not used!
	const VECTOR<string>&	pricedColumns,
	bool ,			/// not used!
	ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_DealDescription* deal = NULL;

	try
	{
		deal = new ARM_DealDescription( vecString, vecTypes, (size_t ) rowsNb, (size_t ) colsNb );

		/// assign object
		if( !assignObject( deal, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete deal;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to create a GramFunctorArg 
////////////////////////////////////////////
////////////////////////////////////////////
extern long ARMLOCAL_ObjManager_Create(
	const VECTOR<string>&	cstNames,
	const VECTOR<long>&    objIds,
	ARM_result&	result, 
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_CstManager* cstManager= NULL;
	vector<ARM_GramFctorArg> objVector;
	ARM_Object* object=NULL;

	for (int i=0;i<objIds.size();i++)
	{
		object = LOCAL_PERSISTENT_OBJECTS->GetObject(objIds[i]);
		ARM_GP_Vector* tmpVector = dynamic_cast< ARM_GP_Vector* >(object);
		ARM_GP_Matrix* tmpMatrix = dynamic_cast< ARM_GP_Matrix* >(object);
		ARM_Curve* tmpCurve = dynamic_cast< ARM_Curve* >(object);
		if (tmpVector )
			objVector.push_back(ARM_GramFctorArg(ARM_VectorPtr(static_cast<ARM_GP_Vector*> (tmpVector->Clone()))));
		else if (tmpMatrix )
			objVector.push_back(ARM_GramFctorArg(ARM_GP_MatrixPtr(static_cast<ARM_GP_Matrix*> (tmpMatrix->Clone()))));
		else if(tmpCurve )
			objVector.push_back(ARM_GramFctorArg(ARM_GP_CurvePtr(static_cast<ARM_Curve*> (tmpCurve->Clone()))));
        else if( dynamic_cast< ARM_DateStrip* >(object) )
			objVector.push_back(ARM_GramFctorArg(ARM_DateStripPtr(static_cast<ARM_DateStrip*> (object->Clone()))));
	}
	try
	{
		cstManager = new ARM_CstManager( cstNames, objVector );
		/// assign object
		if( !assignObject( cstManager, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete cstManager;
		x.DebugPrint();
		ARM_RESULT();
	}

}



////////////////////////////////////////////
//// Function to create a cst manager
////////////////////////////////////////////
extern long ARMLOCAL_CstManager_Create(
	const VECTOR<string>&	cstNames,
	const VECTOR<double>&	values,
	ARM_result& result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_CstManager* cstManager= NULL;

	try
	{
		cstManager = new ARM_CstManager( cstNames, values );

		/// assign object
		if( !assignObject( cstManager, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete cstManager;
		x.DebugPrint();
		ARM_RESULT();
	}
}

extern long ARMLOCAL_GenSec_GetCstManager(
	    const long& GenSecId,
	    ARM_result& result, 
	    long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");


	try
	{
		ARM_GenSecurity* genSec = NULL;
		if( !GetObjectFromId( &genSec, GenSecId, ARM_GENSECURITY ) )
		{
			result.setMsg ("ARM_ERR: generic Security is not of a good type");
			return ARM_KO;
		}

        if(genSec->GetCstManager() != ARM_CstManagerPtr(NULL))
        {
            ARM_CstManager* cstManager= new ARM_CstManager(*(genSec->GetCstManager()));
		    /// assign object
		    if( !assignObject( cstManager, result, objId ) ){
			    return ARM_KO; }
		    else{
			    return ARM_OK; }
        }
        else
        {
			result.setMsg ("ARM_ERR: no cst manager found for this generic Security");
			return ARM_KO;
        }
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}



////////////////////////////////////////////
//// Function to create a generic security
////////////////////////////////////////////
extern long ARMLOCAL_GenSec_Create(
	const VECTOR<string>&	vecString,
	const VECTOR<ARM_GP_VALUE_TYPE>& vecTypes,
	long rowsNb,
	long colsNb,
    const string& payCurveName,
	long cstManagerId,
	bool ExercBoundaryResetFlag,
	bool otherPayoffsFlag,
	const VECTOR<string>&	pricedColumns,
	bool ivFlag,
	ARM_result&	result,
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenSecurity* genSec = NULL;

	try
	{


/// comment this out to activate tracer
/// #define SHOW_EVENT_IN_ARMLOCAL_GENSEC_CREATE

/// to debug release only crash
#if defined(SHOW_EVENT_IN_ARMLOCAL_GENSEC_CREATE)
		if( ARM_EventViewerImp::DebugIsOn )
		{
			CC_Ostringstream os;
			os << "\ntrying to build a generic security : \n";
			time_t aClock;
			time( &aClock);	/// Get time in seconds
			os << "time = " << asctime( localtime( &aClock ) ) << "\n";

            char fOutName[200];

            ARM_GetTmpAbsFile("DebugReleasePb", fOutName);

			ARM_TheEventViewer.Instance()->WriteMessageToFile(fOutName, os.str() );
		}
#endif

		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Generic Security" );

#if defined(SHOW_EVENT_IN_ARMLOCAL_GENSEC_CREATE)
		if( ARM_EventViewerImp::DebugIsOn )
		{
			string msg( "\nRegistering cookies" );

            char fOutName[200];

            ARM_GetTmpAbsFile("DebugReleasePb", fOutName);

			ARM_TheEventViewer.Instance()->WriteMessageToFile(fOutName, msg );
		}
#endif

		ARM_DealDescriptionPtr deal( new ARM_DealDescription( 
			vecString, vecTypes, (size_t ) rowsNb, (size_t ) colsNb, pricedColumns ) );

#if defined(SHOW_EVENT_IN_ARMLOCAL_GENSEC_CREATE)
		if( ARM_EventViewerImp::DebugIsOn )
		{
			string msg( "\nCreated deal description" );

            char fOutName[200];

            ARM_GetTmpAbsFile("DebugReleasePb", fOutName);

			ARM_TheEventViewer.Instance()->WriteMessageToFile(fOutName, msg );
		}
#endif

		ARM_CstManager* cstManager = NULL;
		if( !GetObjectFromIdWithDynamicCastCheckwNull( &cstManager, cstManagerId  ) )
		{
			result.setMsg ("ARM_ERR: cst manager is not of a good type");
			return ARM_KO;
		};

		ARM_CstManagerPtr cstManagerPtr = ARM_CstManagerPtr( NULL );
		if( cstManager )
			cstManagerPtr = ARM_CstManagerPtr( (ARM_CstManager* ) ( cstManager->Clone() ) );
		genSec = new ARM_GenSecurity( deal , payCurveName, cstManagerPtr, ExercBoundaryResetFlag, otherPayoffsFlag, ivFlag );

#if defined(SHOW_EVENT_IN_ARMLOCAL_GENSEC_CREATE)
		if( ARM_EventViewerImp::DebugIsOn )
		{
			string msg( "\nCreated generic security" );

            char fOutName[200];

            ARM_GetTmpAbsFile("DebugReleasePb", fOutName);

            ARM_TheEventViewer.Instance()->WriteMessageToFile(fOutName, msg );
		}
#endif

		/// assign object
		if( !assignObject( genSec, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete genSec;
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete genSec;
		result.setMsg ("ARM_ERR: unrecognized failure in creating the generic security");
		return ARM_KO;
	}

}

////////////////////////////////////////////
//// Function to create a new generic security
//// from an old one
////////////////////////////////////////////
extern long ARMLOCAL_GenSec_ChangeAmericanIntoTrigger_Common(
	const long& GenSecId, const long& PricingModelId, 
	ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	ARM_GenSecurity* genSec = NULL;
	ARM_GenSecurity* newGenSec = NULL;
	ARM_PricingModel* pricingModel = NULL;
	ARM_GenPricer* genpricer = NULL;

	CCString msg ("");

	try
	{
		ARM_GenSecManipulator genSecManipulator;
		
		if( !GetObjectFromId( &genSec, GenSecId, ARM_GENSECURITY ) )
		{
			result.setMsg ("ARM_ERR: generic Security is not of a good type");
			return ARM_KO;
		}

		if( !GetObjectFromId( &pricingModel, PricingModelId, ARM_PRICINGMODEL ) )
		{
			result.setMsg ("ARM_ERR: pricing model is not of a good type");
			return ARM_KO;
		}

		ARM_DealDescription *newDealDesc = new ARM_DealDescription( genSec->GetDealDescription() );
		ARM_DealDescriptionPtr newDealDescPtr( newDealDesc );

		newGenSec = new ARM_GenSecurity( newDealDescPtr, genSec->GetPayModelName(), genSec->GetCstManager() , 
			false );

		/// creates locally a genpricer
		/// and uses it to price; so that exercise boundary is computed. 
		genpricer = new ARM_GenPricer( newGenSec,pricingModel );

		/// price!
		genpricer->Price();
		delete genpricer;

		genSecManipulator.ChangeAmericanIntoTrigger( *newGenSec, newDealDescPtr );

		/// assign object
		if( !assignObject( newGenSec, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
	catch (...)
	{
		delete genSec;
		result.setMsg ("ARM_ERR: unrecognized failure in creating the generic security");
		return ARM_KO;
	}
	return ARM_KO;
}

////////////////////////////////////////////
//// Function to create a new generic security
//// from an old one
////////////////////////////////////////////
extern long ARMLOCAL_GenSec_ChangeMAXPVIntoExercise_Common(
	const long& GenSecId, 
	ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	ARM_GenSecurity* genSec = NULL;
	ARM_GenSecurity* newGenSec = NULL;

	CCString msg ("");

	try
	{
		ARM_GenSecManipulator genSecManipulator;
		
		if( !GetObjectFromId( &genSec, GenSecId, ARM_GENSECURITY ) )
		{
			result.setMsg ("ARM_ERR: generic Security is not of a good type");
			return ARM_KO;
		}

		ARM_DealDescriptionPtr newDealDescPtr = ARM_DescriptionPreprocessor::ChangeMaxPVPreprocessing( genSec->GetDealDescription() );

		newGenSec = new ARM_GenSecurity( newDealDescPtr, genSec->GetPayModelName(), genSec->GetCstManager() , 
			false );

		/// assign object
		if( !assignObject( newGenSec, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
	catch (...)
	{
		delete genSec;
		result.setMsg ("ARM_ERR: unrecognized failure in creating the generic security");
		return ARM_KO;
	}
	return ARM_KO;
}



////////////////////////////////////////////
//// Function to set parse tree on and off
////////////////////////////////////////////
extern long ARMLOCAL_GenSec_SetParseTreeFlag(
	long GenSecId,
	bool ParseTreeFlag,
	ARM_result&	result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");


	try
	{
		ARM_GenSecurity* genSec = NULL;
		if( !GetObjectFromId( &genSec, GenSecId, ARM_GENSECURITY ) )
		{
			result.setMsg ("ARM_ERR: generic Security is not of a good type");
			return ARM_KO;
		};

		genSec->SetParseTreePrint( ParseTreeFlag );
		string txt( "ParseTree :" );
		txt += ParseTreeFlag? "On" : "Off";
		result.setString(txt.c_str());
		
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to create a helper object
////////////////////////////////////////////
extern long ARMLOCAL_GramHelper(
	const string& FuncName,
	ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GramHelper* helper= NULL;

	try
	{
		helper = new ARM_GramHelper( FuncName );

		/// assign object
		if( !assignObject( helper, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete helper;
		x.DebugPrint();
		ARM_RESULT();
	}
}



////////////////////////////////////////////
//// Function to set parse tree on and off
////////////////////////////////////////////
extern long ARMLOCAL_GenSec_SetCircularRefFlag(
	bool CircularRefFlag,
	ARM_result&	result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_GenSecurity::SetLookForCircularRef(CircularRefFlag);
		string txt( "Circular Ref :" );
		txt += CircularRefFlag? "On" : "Off";
		result.setString(txt.c_str());
		
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to get the deal description table
////////////////////////////////////////////
extern long ARMLOCAL_GenSec_GetDealDesTable(
	    long GenSecId,
		VECTOR< string >& dealDesText,
		VECTOR< ARM_GP_VALUE_TYPE >& dealDesFormat,
        long& nbRows,
        long& nbCols,
		ARM_result&	result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");


	try
	{
		ARM_GenSecurity* genSec = NULL;
		if( !GetObjectFromId( &genSec, GenSecId, ARM_GENSECURITY ) )
		{
			result.setMsg ("ARM_ERR: generic Security is not of a good type");
			return ARM_KO;
		}

        const ARM_DealDescription dealDesc(genSec->GetDealDescription());

        nbRows = (long) dealDesc.GetRowsNb();
        nbCols = (long) dealDesc.GetColsNb();

        dealDesText = dealDesc.GetText();

        dealDesFormat = dealDesc.GetFormat();

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to extract a sub deal description in the gen security
////////////////////////////////////////////

extern long ARMLOCAL_GenSec_ExtractSubDealDes(
	const long& GenSecId,
	const string& ColName,
	const vector<string>& cols,
	const string& payCurveName,
	ARM_result&	result,
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenSecurity* newGenSec = NULL;

	try
	{
		ARM_GenSecurity* genSec = NULL;
		if( !GetObjectFromId( &genSec, GenSecId, ARM_GENSECURITY ) )
		{
			result.setMsg ("ARM_ERR: generic Security is not of a good type");
			return ARM_KO;
		}

        const ARM_DealDescription& dealDesc = genSec->GetDealDescription();
		size_t colIndex = dealDesc.GetColIndex(ColName);
		ARM_DealDescriptionPtr subDealDesc = dealDesc.GetSubDescription(1,dealDesc.GetRowsNb()-1,colIndex+1,cols);
		string payCurveNameUsed  = genSec->GetPayModelName() != ""? genSec->GetPayModelName() :  payCurveName;
		newGenSec = new ARM_GenSecurity(subDealDesc,payCurveNameUsed ,genSec->GetCstManager());

		/// assign object
		if( !assignObject( newGenSec, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}
