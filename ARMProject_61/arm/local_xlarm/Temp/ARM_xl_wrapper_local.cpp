#include <ARM\libarm_local\firstToBeIncluded.h>
#include "ARM_xl_wrapper_local.h"
#include "ARM_local_interglob.h"
#include "ARM\local_xlarm\ARM_local_interface.h"

/// FIX FIX to debug release only crash
#include "gpbase/eventviewerfwd.h"

/// FIX FIX to debug release only crash
using ARM::ARM_EventViewerImp;
using ARM::ARM_TheEventViewer;

/*! 
 * This function is factorise most of the
 * excel interface code
 * since any addins just at the end of the day
 * call a functor of the type ARMResultLong2LongFunc
 * with a specific context
 *
 * In order to use this function
 * one has to define properly the functor
 */
void fillXL_Result( const CCString& curClass,
				   ARMResultLong2LongFunc& CreateObjFctor, 
				   ARM_result& C_result,
				   XLOPER& XL_result,
				   bool PersistentInXL,
				   bool GetNameFromResult )
{
	/// to make the infrastructure very robust we put a try catch!

#if !defined(_DEBUG)
	try
	{
#endif		
		/// comment this out to activate tracing
		/// #define SHOW_EVENT_IN_FILLXL_RESULT
		
		CCString stringId;

		const long CreateBrandNew = ARM_NULL_OBJECT_ID;
		long objId;
		long retCode;
		
		/// case of a simple call
		/// no need to know if persistent or not
		if( !PersistentInXL )
		{
			
			retCode = CreateObjFctor(  C_result, CreateBrandNew );
			CCString className = GetNameFromResult? C_result.getShortName() : curClass;

			if(retCode == ARM_OK)
			{
				objId = C_result.getLong ();
				stringId = LocalMakeObjectId (objId, className);
			}
		}
		/// case where we manage persistence
		else
		{
            stringId = GetLastCurCellEnvValue();

			/// first case, the cell was empty here
			if (!stringId)
			{
				retCode = CreateObjFctor(  C_result, CreateBrandNew );
				CCString className = GetNameFromResult? C_result.getShortName() : curClass;
				
				if ( retCode == ARM_OK)
				{
					objId = C_result.getLong ();
					LocalSetCurCellEnvValue (className, objId); 
					stringId = LocalMakeObjectId (objId, className);
				}
			}
			
			/// ok, cell already used
			else
			{
				CCString prevClass = LocalGetStringObjectClass (stringId);
				objId	= LocalGetNumObjectId (stringId);
				retCode = CreateObjFctor( C_result, objId );
				CCString className = GetNameFromResult? C_result.getShortName() : curClass;
				
				/// is it the same type as the previous object!
				if ( className == prevClass )
				{
					if ( retCode == ARM_OK )
					{ 
						LocalSetCurCellEnvValue (className, objId);
						stringId = LocalMakeObjectId (objId, className);
					}
				}
				else 
				{
					// ARM_Object is already freed in CreateObjFctor() !!
					FreeCurCellContentWithoutFreeingARMObject();
					
                    if ( retCode == ARM_OK )
					{
						objId = C_result.getLong ();
						LocalSetCurCellEnvValue (className, objId);
						stringId = LocalMakeObjectId (objId, className);
					}
				}
			}
		}


		/// return the XLOPER object
		if ( retCode == ARM_OK )
		{	
			FreeCurCellErr ();
		
            XL_result.xltype = xltypeStr;
			XL_result.val.str = XL_StrC2StrPascal (stringId);
			XL_result.xltype |= xlbitDLLFree;
		}
		else
		{
			ARM_ARG_ERR();
		}
#if !defined(_DEBUG)
	}
	/// catch ARM_Exception
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_ARG_ERR();
	}

	/// catch the rest
	catch( ... )
	{
		string mssg( "ARM_ERR: unrecognized failure in fillXL_Result with class " );
		char* classChar = curClass;
		mssg += classChar;
		delete classChar;
		C_result.setMsg( mssg.c_str() );
		
        ARM_ARG_ERR();
	}
#endif
}


/// version that gets the short name directly from the C_result
void fillXL_Result_withName(ARMResultLong2LongFunc& CreateObjFctor, 
							ARM_result& C_result,
							XLOPER& XL_result,
							bool PersistentInXL )
{
	fillXL_Result(
		"NOT USED",
		CreateObjFctor, 
		C_result,
		XL_result,
		PersistentInXL,
		true );
}
