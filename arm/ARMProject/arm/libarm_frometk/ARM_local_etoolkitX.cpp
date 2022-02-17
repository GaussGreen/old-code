#pragma warning(disable : 4786)
#pragma warning(disable : 4541)

#include <ARM\libarm_local\firstToBeIncluded.h>

#include <comdef.h>

#include <ARM\libarm_frometk\eToolkitX.h>
#include <ARM\libarm_frometk\ARM_local_etoolkitX.h>

#include <libCCatl\CCatl.h>


void ARM_Etoolkit::Connect()
{
	HRESULT hr;

    Deconnect();

	_variant_t applicationName=(_bstr_t)"Arm_Local";
	_variant_t programMode=(_bstr_t)"STANDARD";

	_variant_t vDatabaseContext=(_bstr_t)itsDataBaseContext.c_str();
	_variant_t vUserName=(_bstr_t)itsUserName.c_str();
	_variant_t vPassWord=(_bstr_t)itsPassword.c_str();
	
	_putenv(itsItConfigDomainsDir.c_str());
	_putenv(itsItDomainName.c_str());
	
	try
	{
		hr = CoInitialize(NULL); 

        SUCCEEDED(hr) ? 0 : throw hr;
		
        hr = CoCreateInstance(CLSID_eToolkit, NULL, CLSCTX_LOCAL_SERVER     , 
			   IID_IeToolkit, (void**)&itsEToolkitptr);
		
        SUCCEEDED(hr) ? 0 : throw hr;

		itsEToolkitptr->Connect(&vUserName, &vPassWord, &applicationName, 
                                &programMode,&vDatabaseContext);	
	}

	catch (...)
	{
        // TMP: May be it's better to throw an exception
    }
}



void ARM_Etoolkit::Deconnect()
{
	try
	{
		if (itsEToolkitptr!=NULL)
		{
			itsEToolkitptr->Disconnect();
			itsEToolkitptr->Release();
			itsEToolkitptr=NULL;
		}
	}

    catch(...)
	{
		itsEToolkitptr=NULL;
	}

	itsEToolkitptr=NULL;
}


void ARM_Etoolkit::Shutdown()
{
	try
	{
		if (itsEToolkitptr!=NULL)
		{
			itsEToolkitptr->ShutDown();
//			itsEToolkitptr->Release();
			itsEToolkitptr=NULL;
		}
	}

    catch(...)
	{
		itsEToolkitptr=NULL;
	}

	itsEToolkitptr=NULL;
}



void ARM_Etoolkit::Execute(CCString command, CCString xmlRequest, CCString & xmlResponse, CCString & messageList)
{
	static char* strnothing = "";

	VARIANT var_xmlResponse;
	VARIANT var_messageList;
	VARIANT var_command;
	VARIANT var_xmlRequest;

   	VariantInit(&var_xmlResponse);
   	VariantInit(&var_messageList);
   	VariantInit(&var_command);
   	VariantInit(&var_xmlRequest);
	
	CCString2VARIANT (command, &var_command);
	CCString2VARIANT (xmlRequest, &var_xmlRequest);
	CCString2VARIANT (strnothing, &var_xmlResponse); //Les variant doivent etre
	CCString2VARIANT (strnothing, &var_messageList);//initialisés sinon pas d'appel à l'activeX !!?

	if (itsEToolkitptr != NULL)
		itsEToolkitptr->Execute(&var_command,&var_xmlRequest,&var_xmlResponse,&var_messageList);

	VARIANT2CCString (var_xmlResponse, xmlResponse);
//	VARIANT2CCString (var_messageList, messageList);

   	VariantClear(&var_xmlResponse);
   	VariantClear(&var_messageList);
   	VariantClear(&var_command);
   	VariantClear(&var_xmlRequest);
}




void ARM_Etoolkit::GetCommsetName(CCString varName, CCString varName2, CCString varAsOf, CCString varType, CCString varCvName, CCString& response, CCString& messageList)
{
    static char* strnothing = "";

	VARIANT var_xmlResponse;
	VARIANT var_messageList;

	VARIANT var_name;
	VARIANT var_name2;
	VARIANT var_type;
	VARIANT var_cvname;
	VARIANT var_asof;

   	VariantInit(&var_xmlResponse);
   	VariantInit(&var_messageList);
   	VariantInit(&var_name);
   	VariantInit(&var_name2);
   	VariantInit(&var_type);
   	VariantInit(&var_cvname);
   	VariantInit(&var_asof);

	CCString2VARIANT (varName, &var_name);
	CCString2VARIANT (varName2, &var_name2);
	CCString2VARIANT (varType, &var_type);
	CCString2VARIANT (varCvName, &var_cvname);
	CCString2VARIANT (varAsOf, &var_asof);
	
	CCString2VARIANT (strnothing, &var_xmlResponse); //Les variant doivent etre
	CCString2VARIANT (strnothing, &var_messageList);//initialisés sinon pas d'appel à l'activeX !!?

	if (itsEToolkitptr != NULL)
		itsEToolkitptr->GetCommsetName(&var_name,&var_name2,&var_asof,&var_type,&var_cvname,&var_xmlResponse,&var_messageList);

	VARIANT2CCString (var_xmlResponse, response);
	VARIANT2CCString (var_messageList, messageList);

   	VariantClear(&var_xmlResponse);
   	VariantClear(&var_messageList);
   	VariantClear(&var_name);
   	VariantClear(&var_name2);
   	VariantClear(&var_type);
   	VariantClear(&var_cvname);
   	VariantClear(&var_asof);
}


void ARM_Etoolkit::GetRefRate(CCString source, CCString ccy, CCString index, CCString tenor, CCString& response, CCString& messageList)
{
    static char* strnothing = "";

	VARIANT var_xmlResponse;
	VARIANT var_messageList;

	VARIANT var_source;
	VARIANT var_ccy;
	VARIANT var_index;
	VARIANT var_tenor;

   	VariantInit(&var_xmlResponse);
   	VariantInit(&var_messageList);
   	VariantInit(&var_source);
   	VariantInit(&var_ccy);
   	VariantInit(&var_index);
   	VariantInit(&var_tenor);

	CCString2VARIANT (source, &var_source);
	CCString2VARIANT (ccy, &var_ccy);
	CCString2VARIANT (index, &var_index);
	CCString2VARIANT (tenor, &var_tenor);
	
	CCString2VARIANT (strnothing, &var_xmlResponse); //Les variant doivent etre
	CCString2VARIANT (strnothing, &var_messageList);//initialisés sinon pas d'appel à l'activeX !!?

	if (itsEToolkitptr != NULL)
		itsEToolkitptr->getREFRATE(&var_source,&var_ccy,&var_index,&var_tenor,&var_messageList, &var_xmlResponse);

	VARIANT2CCString (var_xmlResponse, response);
//	VARIANT2CCString (var_messageList, messageList);

   	VariantClear(&var_xmlResponse);
   	VariantClear(&var_messageList);
   	VariantClear(&var_source);
   	VariantClear(&var_ccy);
   	VariantClear(&var_index);
   	VariantClear(&var_tenor);
}
