/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_mktdata.cpp,v $
 * Revision 1.1  2004/03/10 15:08:43  ebenhamou
 * Initial version
 *
 */

/*! \file ARM_local_gp_mktdata.cpp
 *
 *  \brief file for the mkt data of the generic pricer local addins functions
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2003
 */


#include "firstToBeIncluded.h"
#include "ARM_local_gp_mktdata.h"
#include "ARM_local_wrapper.h"
#include <GP_Infra\gpinfra\mktdatamanagerrep.h>
#include <GP_Infra\gpinfra\mktdatamanager.h>


#include <crv\zerocurv.h>
#include <crv\volint.h>
#include <mod\model.h>
#include <glob\dates.h>


//// using the namespace directive to access ARM object!
using ARM::ARM_MarketData_ManagerRep;
using ARM::ARM_MarketData_ManagerImp;

////////////////////////////////////////////
//// Function to create a mkt data manager
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_Create(
	const ARM_Date& asOf, 
    ARM_result&	result, 
    long			objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_MarketData_ManagerRep* mktDataManager= NULL;

	try
	{
        mktDataManager = new ARM_MarketData_ManagerRep(asOf);

		/// assign object
		if( !assignObject( mktDataManager, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete mktDataManager;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to register a zero curve
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_RegisterZCCurve(
	long mktDataManagerId,
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& source, 
	long zcCurveId,
    ARM_result&	result, 
    long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	ARM_MarketData_ManagerRep* mktDataManager = NULL;

	try
	{
		ARM_MarketData_ManagerRep* previousMktDataManager = NULL;
		if( !GetObjectFromId( &previousMktDataManager, mktDataManagerId, ARM_MKTDATAMANAGER ) )
		{
			result.setMsg ("ARM_ERR: mkt data manager is not of a good type");
			return ARM_KO;
		};

		ARM_ZeroCurve* zcCurve = NULL;
		if( !GetObjectFromId( &zcCurve, zcCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: zero curve is not of a good type");
			return ARM_KO;
		};

		/// clone for safety
		mktDataManager = (ARM_MarketData_ManagerRep*) previousMktDataManager->Clone();
		mktDataManager->RegisterZCCurve( indexName, ccy, cvName, asOf, source, zcCurve );

		/// assign object
		if( !assignObject( mktDataManager, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete mktDataManager;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to get a zero coupon curve
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_GetZCCurve(
	long mktDataManagerId,
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& source, 
	ARM_result&	result, 
    long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_MarketData_ManagerRep* mktDataManager = NULL;
		if( !GetObjectFromId( &mktDataManager, mktDataManagerId, ARM_MKTDATAMANAGER ) )
		{
			result.setMsg ("ARM_ERR: mkt data manager is not of a good type");
			return ARM_KO;
		};

		/// clone the object
		ARM_ZeroCurve* zcCurve = mktDataManager->GetZCCurveAndClone( indexName, ccy, cvName, asOf, source );

		/// assign object (we use GetName to avoid strict typing!)
		if( !assignObject( zcCurve, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure in getting a zc curve from mkt data manager");
		return ARM_KO;
	}
}






////////////////////////////////////////////
//// Function to get a zero coupon curve
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_GetZCCurveKey(
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& source, 
	ARM_result&	result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		string mssg = ARM_MarketData_ManagerImp::CoinZCObjKey( indexName, ccy, cvName, asOf, source );
		result.setString( mssg.c_str() );
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}



////////////////////////////////////////////
//// Function to register a VolCurve
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_RegisterVolCurve( 
	long mktDataManagerId,
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType, 
	const string& volType, 
	const string& source, 
	long volCurveId,
    ARM_result&	result, 
    long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	/// create empty mktDataManager
	ARM_MarketData_ManagerRep* mktDataManager = NULL;

	try
	{
		ARM_MarketData_ManagerRep* previousMktDataManager = NULL;
		if( !GetObjectFromId( &previousMktDataManager, mktDataManagerId, ARM_MKTDATAMANAGER ) )
		{
			result.setMsg ("ARM_ERR: mkt data manager is not of a good type");
			return ARM_KO;
		};

		ARM_VolCurve* volCurve = NULL;
		if( !GetObjectFromId( &volCurve, volCurveId, ARM_VOL_CURVE ) )
		{
			result.setMsg ("ARM_ERR: vol curve is not of a good type");
			return ARM_KO;
		};

		/// clone for safety
		mktDataManager = (ARM_MarketData_ManagerRep*) previousMktDataManager->Clone();
		mktDataManager->RegisterVolCurve( indexName, ccy, cvName, asOf, volMktType, volType, source, volCurve );

		/// assign object
		if( !assignObject( mktDataManager, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete mktDataManager;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to get a vol curve
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_GetVolCurve( 
	long mktDataManagerId,
	const string& indexName,
	const string& ccy,
	const string& cvName,
	const ARM_Date& asOf,
	const string& volMktType,
	const string& volType,
	const string& source,
	ARM_result&	result, 
    long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_MarketData_ManagerRep* mktDataManager = NULL;
		if( !GetObjectFromId( &mktDataManager, mktDataManagerId, ARM_MKTDATAMANAGER ) )
		{
			result.setMsg ("ARM_ERR: mkt data manager is not of a good type");
			return ARM_KO;
		};

		/// clone the object
		ARM_VolCurve* volCurve = mktDataManager->GetVolCurveAndClone( indexName, ccy, cvName, asOf, volMktType, volType, source );

		/// assign object (we use GetName to avoid strict typing!)
		if( !assignObject( volCurve, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure in getting a vol curve curve from mkt data manager");
		return ARM_KO;
	}

}



////////////////////////////////////////////
//// Function to get a vol curve
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_GetVolCurveKey(
	const string& indexName,
	const string& ccy,
	const string& cvName,
	const ARM_Date& asOf,
	const string& volMktType,
	const string& volType,
	const string& source,
	ARM_result&	result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		string mssg = ARM_MarketData_ManagerImp::CoinVolObjKey( indexName, ccy, cvName, asOf, volMktType, volType, source );
		result.setString( mssg.c_str() );
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to register a Vol Mkt Model
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_RegisterVolMktModel( 
	long mktDataManagerId,
	const string& indexName,
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf, 
	const string& volMktType,
	const string& source, 
	long volMktModelId,
    ARM_result&	result, 
    long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	/// create empty mktDataManager
	ARM_MarketData_ManagerRep* mktDataManager = NULL;

	try
	{
		ARM_MarketData_ManagerRep* previousMktDataManager = NULL;
		if( !GetObjectFromId( &previousMktDataManager, mktDataManagerId, ARM_MKTDATAMANAGER ) )
		{
			result.setMsg ("ARM_ERR: mkt data manager is not of a good type");
			return ARM_KO;
		};

		ARM_Model* volMktModel = NULL;
		if( !GetObjectFromId( &volMktModel, volMktModelId, ARM_MODEL ) )
		{
			result.setMsg ("ARM_ERR: vol mkt model is not of a good type");
			return ARM_KO;
		};

		/// clone for safety
		mktDataManager = (ARM_MarketData_ManagerRep*) previousMktDataManager->Clone();
		mktDataManager->RegisterVolMktModel( indexName, ccy, cvName, asOf, volMktType, source, volMktModel );

		/// assign object
		if( !assignObject( mktDataManager, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}

	catch(Exception& x)
	{
		delete mktDataManager;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to get a vol mkt model
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_GetVolMktModel(
	long mktDataManagerId,
	const string& indexName,
	const string& ccy,
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType,
	const string& source,
	ARM_result&	result, 
    long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_MarketData_ManagerRep* mktDataManager = NULL;
		if( !GetObjectFromId( &mktDataManager, mktDataManagerId, ARM_MKTDATAMANAGER ) )
		{
			result.setMsg ("ARM_ERR: mkt data manager is not of a good type");
			return ARM_KO;
		};

		/// clone the object
		ARM_Model* mktModel = mktDataManager->GetVolMktModelAndClone( indexName, ccy, cvName, asOf, volMktType, source );

		/// assign object
		if( !assignObject( mktModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}



////////////////////////////////////////////
//// Function to get a vol curve
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_GetVolMktModelKey(
	const string& indexName,
	const string& ccy,
	const string& cvName,
	const ARM_Date& asOf,
	const string& volMktType,
	const string& source,
	ARM_result&	result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		string mssg = ARM_MarketData_ManagerImp::CoinVolMktModelObjKey( indexName, ccy, cvName, asOf, volMktType, source );
		result.setString( mssg.c_str() );
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}




////////////////////////////////////////////
//// Function to register data
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_RegisterData( 
	long mktDataManagerId,
	const string& objectKey, 
	long objectId,
    bool isFill,
    ARM_result&	result, 
    long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	/// create empty mktDataManager
	ARM_MarketData_ManagerRep* mktDataManager = NULL;

	try
	{
		ARM_MarketData_ManagerRep* previousMktDataManager = NULL;
        string updateMsg("Update:");
		if( !GetObjectFromId( &previousMktDataManager, mktDataManagerId, ARM_MKTDATAMANAGER ) )
		{
			result.setMsg ("ARM_ERR: mkt data manager is not of a good type");
			return ARM_KO;
		};

		/// not type validation as it can be anything!
		ARM_Object* object = (ARM_Object*) LOCAL_PERSISTENT_OBJECTS->GetObject(objectId);

        if(isFill)
            mktDataManager = previousMktDataManager;
        else
            /// clone for safety
            mktDataManager = (ARM_MarketData_ManagerRep*) previousMktDataManager->Clone();

		mktDataManager->RegisterData(objectKey,object);
        updateMsg += " Done!!";
        if(isFill)
        {
            result.setString(updateMsg.c_str());
			return ARM_OK;
        }

		/// assign object
		if( !assignObject( mktDataManager, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete mktDataManager;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to get any object 
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_GetData(
	long mktDataManagerId,
	const string& objectKey,
	ARM_result&	result, 
    long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_MarketData_ManagerRep* mktDataManager = NULL;
		if( !GetObjectFromId( &mktDataManager, mktDataManagerId, ARM_MKTDATAMANAGER ) )
		{
			result.setMsg ("ARM_ERR: mkt data manager is not of a good type");
			return ARM_KO;
		};

		/// clone the object
		ARM_Object* object = mktDataManager->GetDataAndClone(objectKey );

		/// assign object (we use GetName to avoid strict typing!)
		if( !assignObject( object, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		result.setMsg ("ARM_ERR: unrecognized failure in getting data from mkt data mger");
		return ARM_KO;
	}
}


////////////////////////////////////////////
//// Function to set detail mode on and off
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_SetDetailMode(
	long mktDataManagerId,
	bool detailModeFlag,
	ARM_result&	result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_MarketData_ManagerRep* mktDataManager = NULL;
		if( !GetObjectFromId( &mktDataManager, mktDataManagerId, ARM_MKTDATAMANAGER ) )
		{
			result.setMsg ("ARM_ERR: mkt data manager is not of a good type");
			return ARM_KO;
		};

		string txt( "Detail Mode" );
		txt += mktDataManager->SetDetailMode(detailModeFlag)? ": On" : " Off";
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
//// Function to change the date on a mkt data manager
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_ChangeDate(
	long mktDataManagerId,
	const ARM_Date& newDate,
	ARM_result&	result, 
    long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_MarketData_ManagerRep* mktDataManager = NULL;

	try
	{
		ARM_MarketData_ManagerRep* previousMktDataManager = NULL;
		if( !GetObjectFromId( &previousMktDataManager, mktDataManagerId, ARM_MKTDATAMANAGER ) )
		{
			result.setMsg ("ARM_ERR: mkt data manager is not of a good type");
			return ARM_KO;
		};

		/// clone for safety
		mktDataManager = (ARM_MarketData_ManagerRep*) previousMktDataManager->Clone();
		mktDataManager->SetAsOfDate( newDate );

		/// assign object
		if( !assignObject( mktDataManager, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete mktDataManager;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to reset all data
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_ResetAllData(
	long mktDataManagerId,
	ARM_result&	result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_MarketData_ManagerRep* mktDataManager = NULL;
		if( !GetObjectFromId( &mktDataManager, mktDataManagerId, ARM_MKTDATAMANAGER ) )
		{
			result.setMsg ("ARM_ERR: mkt data manager is not of a good type");
			return ARM_KO;
		};

		mktDataManager->ResetAllData();
		result.setString( "Reseted All Data" );
		
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}



////////////////////////////////////////////
//// Function to reset only the current mkt data manager rep
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_ResetMyData(
	long mktDataManagerId,
    ARM_result&	result, 
    long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	ARM_MarketData_ManagerRep* mktDataManager = NULL;

	try
	{
		ARM_MarketData_ManagerRep* previousMktDataManager = NULL;
		if( !GetObjectFromId( &previousMktDataManager, mktDataManagerId, ARM_MKTDATAMANAGER ) )
		{
			result.setMsg ("ARM_ERR: mkt data manager is not of a good type");
			return ARM_KO;
		};

		/// clone for safety
		mktDataManager = (ARM_MarketData_ManagerRep*) previousMktDataManager->Clone();
		mktDataManager->ResetMyData();

		/// assign object
		if( !assignObject( mktDataManager, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete mktDataManager;
		x.DebugPrint();
		ARM_RESULT();
	}
}

