/*! \file ARM_local_wrapper.h
 *
 *  \brief file to factorise the interface
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date June 2003
 */

#ifndef ARM_LOCAL_WRAPPER_H
#define ARM_LOCAL_WRAPPER_H

#include "firstToBeIncluded.h"
#include "ARM_result.h"
#include <glob\armglob.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <GP_Base\gpbase\rootobject.h>

using ARM::ARM_RootObject;

/*!
 * function to check the global persistance
 */
bool GlobalPersistanceOk( ARM_result& result);


/*!
 * General function for assigning an object
 * in a function wrap to excel
 * Beaware that this does not do the clean up
 * of pointor... hence if by any chance
 * you assign pointor, you will need to clean 
 * the object
 */


template< class T >
	bool assignObject( T* newObj, 
		ARM_result& result, 
		long objId )
{
	/// get if necessary the short name from the newObject
	ARM_RootObject* rootObj = dynamic_cast<ARM_RootObject*>( newObj);
	if( rootObj )
		result.setShortName( rootObj->GetExportShortName().c_str() );

	if( (objId == ARM_NULL_OBJECT_ID ) )
	{
		CREATE_GLOBAL_OBJECT();
		long newObjId = LOCAL_PERSISTENT_OBJECTS->SetPersistent( (ARM_Object*) newObj  );
		
		if (newObjId == RET_KO)
		{
			delete newObj;
			result.setMsg ("ARM_ERR: Pb with inserting object");				
			return false;
		}
		result.setLong(newObjId );
	}
	else
	{
		T* oldObj= (T*) LOCAL_PERSISTENT_OBJECTS->GetARMObject(objId);

		if (typeid(oldObj) == typeid(newObj))
		{
			delete oldObj;
			
			/// use GetARMObject to avoid confusion with SDK GetObject (basically inappropriate name)
			/// thanks to the inline code should be as quick as GetObject!
			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*) newObj, objId);

			result.setLong(objId );

		}
		// If the class of the two objects are different we create a brand new id.
		// The previous objects wil be destroyed it the fillXLL (FreeCurCellContent)
		else
		{
			long newObjId = LOCAL_PERSISTENT_OBJECTS->SetPersistent( (ARM_Object*) newObj  );
			result.setLong(newObjId);
		}
	}
	return true;
}



/*!
 * function to check vector size
 */
template <class T, class U> bool sameVectorSize(
	  const VECTOR<T>& X,
	  const CCString& XName,
	  const VECTOR<U>& Y,
	  const CCString& YName,	
	  ARM_result& result
										  )
{
	if ( X.size() != Y.size() )
	{
		result.setMsg( CCString( "ARM_ERR: " ) + CCString( XName ) + CCString( " size is " ) 
			+ CCString( X.size() ) + CCString( "vs " ) 
			+ YName + CCString( " size " ) + CCString( Y.size() ) );
		return false;
	}
	else return true;
}



/*!
 * version that get the object from its Id
 * with no test specific to ARM_NULL_OBJECT 
 */
template <class T>
bool GetObjectFromId( T** myObj, long ObjId, ARM_CLASS_NAME className )
{
	*myObj =  (T *) LOCAL_PERSISTENT_OBJECTS->GetObject(ObjId);

	/// validation test to check that it is of the good type:
	/// it can be either a derived object from the root class
	/// or of the type class name... for object with no derivation
	if( *myObj )
	{
		if ( ( LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK( *myObj, className)
			||	LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK( *myObj, className ) ) == false)
		{
			*myObj = NULL;
			return false;
		}
		else
			return true;
	}
	else
		return false;

}


template <class T>
	bool GetObjectFromIdWithDynamicCastCheck( T** myObj, long ObjId)
{
	ARM_Object* previousObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(ObjId);
	*myObj = dynamic_cast<T*>(previousObj);
	return *myObj != NULL;
}

template <class T>
	bool GetObjectFromIdWithDynamicCastCheckandMsge( T** myObj, long ObjId, const CCString& msg,ARM_result& result)
{
	ARM_Object* previousObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(ObjId);
	*myObj = dynamic_cast<T*>(previousObj);
	bool isObjNull =  *myObj == NULL;
	if(isObjNull)
		result.setMsg (CCString("ARM_ERR:") + CCString( msg ) + CCString(" is not of a good type"));

	return isObjNull;
}

template <class T>
	bool GetObjectFromIdWithDynamicCastCheckwNullandMsge( T** myObj, long ObjId, const CCString& msg,ARM_result& result)
{
	if( ObjId == ARM_NULL_OBJECT )
		return false;
	else
	{
		delete *myObj;
		*myObj = NULL;
		return GetObjectFromIdWithDynamicCastCheckandMsge( myObj, ObjId ,msg,result);
	}
}


template <class T>
	bool GetObjectFromIdWithDynamicCastCheckwNull( T** myObj, long ObjId)
{
	if( ObjId == ARM_NULL_OBJECT )
	{
		myObj = NULL;
		return true;
	}
	else
		return GetObjectFromIdWithDynamicCastCheck( myObj, ObjId );
}




/*!
 * version that test for ARM_NULL_OBJECT 
 */
template <class T>
	bool GetObjectFromIdwNull( T** myObj, long ObjId, ARM_CLASS_NAME className )
{
	//// this version tests whether the objId is an ARM_NULL_OBJECT
	/// if this is the case, then it retun a NULL pointor and carries on
	if( ObjId == ARM_NULL_OBJECT )
	{
		myObj = NULL;
		return true;
	}
	else
		return GetObjectFromId( myObj, ObjId, className );
}



#endif /* ARM_LOCAL_WRAPPER_H */