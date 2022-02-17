/*----------------------------------------------------------------------------*/
/*                                                                            */
/* Set of functions setting ARM_CORBA_PARAM parameters (for client)           */
/*                                                                            */
/*----------------------------------------------------------------------------*/


#define ARM_CORBA_PARAM_clnt_cpp
#include "ARM_CORBA_PARAM_clnt.h"


void PARAM_GetLongVector(const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, VECTOR<long>& vecVals)
{
	vecVals.clear ();
 
	for (CORBA::Long i = 0; i < param->val.longArray().length(); i++)
	{
		vecVals.push_back ((param->val.longArray())[i]);
	}
}


void PARAM_SetLongVector(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, const VECTOR<long>& vecVals)
{
	param->IsParam = ARM_CorbaRequestModule::OBJ_PARAM;

	param->nbValues = vecVals.size ();

	long nbVal = vecVals.size ();
	long* arVals = new long [nbVal];
	for(CORBA::Long i = 0; i < nbVal; i++)
	{
		arVals[i] = vecVals [i];
	} 

	ARM_CorbaRequestModule::ARM_CORBA_LONG_ARRAY_TYPE arrayLong(nbVal, nbVal, arVals, 0);

	param->val.longArray(arrayLong);

	delete arVals;
}


void PARAM_SetDoubleVector(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, const VECTOR<double>& vecVals)
{
	param->IsParam = ARM_CorbaRequestModule::OBJ_PARAM;

	param->nbValues = vecVals.size ();

	long nbVal = vecVals.size ();
	CORBA::Double* arVals = new CORBA::Double [nbVal];
	for(CORBA::Long i = 0; i < nbVal; i++)
	{
		arVals[i] = vecVals [i];
	}

	ARM_CorbaRequestModule::ARM_CORBA_DOUBLE_ARRAY_TYPE arrayDouble(nbVal, nbVal, arVals, 0);

	param->val.doubleArray(arrayDouble);

	delete arVals;
}


void PARAM_GetDoubleVector(const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, VECTOR<double>& vecVals)
{
	vecVals.clear ();
 
	for (CORBA::Long i = 0; i < param->val.doubleArray().length(); i++)
	{
		vecVals.push_back ((param->val.doubleArray())[i]);
	}
}


void PARAM_GetStringVector(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, VECTOR<CCString>& vecVals)
{
	vecVals.clear ();
 
	for (CORBA::Long i = 0; i < param->val.stringArray().length(); i++)
	{
		vecVals.push_back (CCString ((const char*)(param->val.stringArray())[i]));
	}
}


void PARAM_SetStringVector(ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, const VECTOR<CCString>& vecVals)
{
	param->IsParam = ARM_CorbaRequestModule::OBJ_PARAM;

	param->nbValues = vecVals.size ();

	long nbVal = vecVals.size ();
	char** arVals = new char* [nbVal];
	for(CORBA::Long i = 0; i < nbVal; i++)
	{
		arVals[i] = (char*)vecVals [i];
	} 

	ARM_CorbaRequestModule::ARM_CORBA_STRING_ARRAY_TYPE arrayString(nbVal, nbVal, arVals, 0);

	param->val.stringArray(arrayString);

	for(i = 0; i < nbVal; i++)
	{
		if(arVals[i])
		{
			free (arVals[i]);
			arVals[i] = NULL;
		}
	}

	delete arVals;
}


void PARAM_Print_message (const ARM_CorbaRequestModule::ARM_CORBA_PARAM* param, long idx)
{
	MSG_printf_message (MSG_INFO, "====> ARM_CORBA_PARAM : %ld", idx);
 
	MSG_printf_message (MSG_INFO, "nbValues = %d", PARAM_GetNbValues (param));
 
	CORBA::ULong i = 0;
 
	switch(PARAM_GetValueType (param))
    	{
		case ARM_CorbaRequestModule::ARM_CORBA_DOUBLE:
		{
			MSG_printf_message (MSG_INFO, "val = %lf", param->val.doubleVal ());
		};
		break;
 
		case ARM_CorbaRequestModule::ARM_CORBA_LONG:
		{
			MSG_printf_message (MSG_INFO, "val = %ld", param->val.longVal ());
		};
        	break;

		case ARM_CorbaRequestModule::ARM_CORBA_STRING:
		{
			MSG_printf_message (MSG_INFO, "val = %s", param->val.stringVal ());
		};
		break;
 
		case ARM_CorbaRequestModule::ARM_CORBA_LONG_ARRAY:
		{
			for(i = 0; i < param->val.longArray().length(); i++)
			{
				MSG_printf_message (MSG_INFO, "val[%d] = %ld", i, param->val.longArray ()[i]);
			}
		};
		break;
 
		case ARM_CorbaRequestModule::ARM_CORBA_DOUBLE_ARRAY:
		{
			for (i = 0; i < param->val.doubleArray().length (); i++)
			{
				MSG_printf_message (MSG_INFO, "val[%d] = %lf", i, param->val.doubleArray ()[i]);
			}
		};
		break;

		case ARM_CorbaRequestModule::ARM_CORBA_STRING_ARRAY:
		{
			for(i = 0; i < param->val.stringArray().length(); i++)
			{
				MSG_printf_message (MSG_INFO, "val[%d] = %s", i, (const char*)(param->val.stringArray ()[i]));
			}
		};
		break;
 
		default:
		break;
	}
 
	MSG_printf_message (MSG_INFO, "<==== ARM_CORBA_PARAM");
}


// EOF %M%
