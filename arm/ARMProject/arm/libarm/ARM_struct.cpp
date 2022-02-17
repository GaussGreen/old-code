#include "ARM_interglob.h"
#include "ARM_struct.h"

long ARM_STRUCTURE (VECTOR<long> instsId, VECTOR<double>& coeffs,
			        ARM_result& result, long objId)
{
	long idx = 0;
	long ARM_REQUEST_ID;
	long ARM_REQUEST_NBPAR = 3;
	CCString stringObjectId;

	/*--- parameters checking ---*/
	if(instsId.size () != coeffs.size ())
	{
		result.setMsg ("ARM_ERR: Assets and weights arrays must have same size");
		return ARM_KO;
	}
			
	try
	{
		ARM_CORBA_init ();

		if(CORBA_OBJECT_INTERFACE)
		{
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
			ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;
			
			stringObjectId = GetLastCurCellEnvValue ();

			if (objId != -1)
			{
				ARM_REQUEST_ID = RPC_SETSTRUCTURE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
				REQ_SetLong (reqIn, objId, 0);
				idx = 1;
			}
			else
			{
				ARM_REQUEST_ID = RPC_CREATESTRUCTURE;
				reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
				idx = 0;
			}

			REQ_SetLong (reqIn, instsId.size (), idx);
			REQ_SetLongVector (reqIn, instsId, ++idx);
			REQ_SetDoubleVector (reqIn, coeffs, ++idx);
			
			/*--- CORBA Server Call ---*/

 	  		CORBA_OBJECT_INTERFACE->Send(*reqIn, reqOut);

			result.set (reqOut);

			/*--- requests freeing ---*/
			REQ_Delete (reqOut);
			REQ_Delete (reqIn);

			ARM_RESULT();
		}
	}
	catch(const CORBA::Exception &e)
	{
		CORBA_ERR();
	}

	return ARM_KO;
}

// EOF %M%