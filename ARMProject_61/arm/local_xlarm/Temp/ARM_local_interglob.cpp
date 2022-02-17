#include <ARM\libarm_local\firstToBeIncluded.h>
#include <Winsock2.h>
#include <RC4\RC4.h>

#include "ARM_local_interglob.h"
#include "XL_local_xlarm_common.h"
#include "ARM_local_license.h"
#include <ARM\libarm_local\ARM_local_init.h>
#include <ARM\libarm_local\ARM_local_persistent.h>
#include <ARM\libarm_local\ARM_local_glob.h>

#include <ARM\libarm_frometk\ARM_local_etoolkit.h>

#include <ARM\libarm_local\arm_local_class.h>
//#include <ARM\libarm\ARM_result.h>
#include <glob\expt.h>
#include <gpbase\env.h>

#include <libCCtools++\CCString.h>
#include <libCCxll\CCxll.h>

#include <time.h>
#include <cstdio>


#include "ICMKernel\glob\icm_enums.h"

#define MATURITY_SIZE		10

// pour NAG et la compil Multithreaded DLL
//FIXMEFRED: problem to compile
//#ifndef _PURIFY
//#define _PURIFY 
//	extern "C" 
//	{
//		int __mb_cur_max;
//		unsigned short* _pctype;
//	}
//#endif

ARM_OBJECTS_LIST_t* ARM_OBJECTS_LIST;
ARM_ERRORS_LIST_t* ARM_ERRORS_LIST;

ARMLOCAL_Init* armlocal_init = NULL;


//	--------------------------------------------------------------------------------------------------

VECTOR<CCString> GetIpAddress(CCString hostname) 
{ 
	HOSTENT*			lpHost=NULL; 
	struct				sockaddr_in   dest; 
	VECTOR<CCString>	csAdrIp;
	CCString			str;
	CCString			csTemp = hostname; 
	WORD				wVersionRequested; 
	WSADATA				wsaData; 
	int					err; 

	wVersionRequested = MAKEWORD( 2, 2 ); //WinSock DLL supports 2.2 
	err = WSAStartup( wVersionRequested, &wsaData ); 
	if ( err != 0 ) 
	{ 
		csAdrIp.push_back("ERR_INIT");	
		return csAdrIp; 
	} 

	lpHost = gethostbyname((const char*)hostname); 
	if (lpHost != NULL) 
	{ 
		for(int i=0; lpHost->h_addr_list[i] != NULL ;i++) 
		{
			memcpy(&(dest.sin_addr), lpHost->h_addr_list[i],lpHost->h_length); 
			str = inet_ntoa(dest.sin_addr); 
			csAdrIp.push_back (str);
		}
	}
	return csAdrIp;
} 


CCString PcName() 
{ 
	char	szUserName[100]; 
	DWORD	dwSize; 
	CCString ccUserName;

	dwSize = sizeof(szUserName); 
	GetComputerName(szUserName,&dwSize);
	ccUserName = szUserName;
	return ccUserName; 
} 


VECTOR<CCString> GetPcIp() 
{
	CCString pcName;
	CCString ccIpAddress;

	pcName = PcName(); 
	return GetIpAddress(pcName); 
} 


int GetLicense()
{
	FILE* Fp;
	if ((Fp = fopen("C:/Program Files/ARM/LicenseForLocal.txt","r")) == NULL)
		return 0;

	char myCrypted [128];
	
	int i = fread(myCrypted,sizeof(char),128,Fp);
	fclose(Fp);

	rc4_key rKey;
	unsigned char myKey [9] = "33export";

	prepare_key((unsigned char*) myKey,strlen((const char*) myKey),&rKey);
	
	rc4((unsigned char*)myCrypted,i,&rKey);

	CCString tmpLastDate(11,myCrypted+19);
	CCString myLastDate(10,myCrypted+19);

	const char* cLastDate = tmpLastDate;
	long sumChar (0);

	for (i=0;i<10;i++)
		sumChar+=(unsigned long)cLastDate[i];

	if (cLastDate[10] != (sumChar % 97))
		return 0;

	myLastDate.Replace(".",'/');
	myLastDate.Replace(".",'/');

	VECTOR<CCString> vectStr;
	myLastDate.Parser('.',vectStr);

	for(i=0;i<vectStr.size();i++)
	{
		if (atoi(vectStr[i]) == 0)
			return 0;
	}

	double myDate = Local_ARMDATE2XLDATE(myLastDate);

	char* cDateNow = DAT_gmt_to_fstring(DAT_now,"%d.%m.%Y");
	
	double myDateNow = Local_ARMDATE2XLDATE(cDateNow);

	if (myDateNow > myDate)
		return 0;
	else
		return 1;
}


int isAuthorized()
{

	int val = 0;

	CCString domain = getenv ("USERDOMAIN");
	domain.toUpper ();

	VECTOR<CCString> ccIpAddress = GetPcIp();
	
	VECTOR<CCString> vIpAddress;

	if ((domain == "MARCHES") || (domain == "CM") || (domain == "CDC_USA_NT") || (domain == "CMNANET") || (domain == "CHILD") || (domain == "CIB"))
	{
		int real_size = ccIpAddress.size ();

		for (int i=0; i<real_size; i++)
		{
			ccIpAddress[i].Parser('.',vIpAddress);

			// Test sur les Adresses IP à Paris Austerlitz (10.14)
			if (
				(atoi((const char*) vIpAddress[0]) == 10) &&
				(atoi((const char*) vIpAddress[1]) >= 14) &&
				(atoi((const char*) vIpAddress[1]) <= 15) 
			   )
				 val = 1;

			// Test sur les Adresses IP Natexis (10.248.86.32) 
			if (
				(atoi((const char*) vIpAddress[0]) == 10) &&
				(atoi((const char*) vIpAddress[1]) == 248)
			   )
				 val = 1;
			// Test sur les Adresses IP à Paris (10.9.16->10.9.24)
			else if (
				(atoi((const char*) vIpAddress[0]) == 10) &&
				(atoi((const char*) vIpAddress[1]) == 9) &&
				(atoi((const char*) vIpAddress[2]) >= 16) &&
				(atoi((const char*) vIpAddress[2]) <= 24)
			   )
				 val = 1;

			// Test sur les Adresses IP à Fkt (10.9.103)
			else if (
				(atoi((const char*) vIpAddress[0]) == 10) &&
				(atoi((const char*) vIpAddress[1]) == 9) &&
				(atoi((const char*) vIpAddress[2]) == 103)
			   )
				 val = 1;

			// Test sur les adresses IP du réseau Back (10.9.64->10.9.67)
			else if (
					 (atoi((const char*) vIpAddress[0]) == 10) &&
					 (atoi((const char*) vIpAddress[1]) == 9) &&
					 (atoi((const char*) vIpAddress[2]) >= 64) &&
					 (atoi((const char*) vIpAddress[2]) <= 67)
					)
				val = 1;

			// Test sur les adresses IP du réseau CNCE (10.9.39)
			else if (
					 (atoi((const char*) vIpAddress[0]) == 10) &&
					 (atoi((const char*) vIpAddress[1]) == 9) &&
					 (atoi((const char*) vIpAddress[2]) == 39)
					)
				val = 1;

			// Test sur les adresses IP de Londres (10.9.30->10.9.31)
			else if (
					 (atoi((const char*) vIpAddress[0]) == 10) &&
					 (atoi((const char*) vIpAddress[1]) == 9) &&
					 (atoi((const char*) vIpAddress[2]) >= 30) &&
					 (atoi((const char*) vIpAddress[2]) <= 31)
					)
				val = 1;

			// Tests sur les adresses IP de Tokyo (158.156.162)
			else if (
					 (atoi((const char*) vIpAddress[0]) == 158) &&
					 (atoi((const char*) vIpAddress[1]) == 156) &&
					 (atoi((const char*) vIpAddress[2]) == 162)
					)
				val = 1;

			// Tests sur les adresses IP du site de BackUp Arcueil
			else if (
					 (atoi((const char*) vIpAddress[0]) == 10) &&
					 (atoi((const char*) vIpAddress[1]) == 9) &&
					 (atoi((const char*) vIpAddress[2]) == 44) 
					)
				val = 1;

			// Tests sur les adresses IP du site de BackUp
			else if (
					 (atoi((const char*) vIpAddress[0]) == 10) &&
					 (atoi((const char*) vIpAddress[1]) == 9) &&
					 (atoi((const char*) vIpAddress[2]) == 48) 
					)
				val = 1;

			// Tests sur les adresses IP du site de HongKong
			else if (
					 (atoi((const char*) vIpAddress[0]) == 10) &&
					 (atoi((const char*) vIpAddress[1]) == 9) &&
					 (atoi((const char*) vIpAddress[2]) == 200) 
					)
				val = 1;

			// Tests sur les adresses IP du site de Tokyo
			else if (
					 (atoi((const char*) vIpAddress[0]) == 10) &&
					 (atoi((const char*) vIpAddress[1]) == 9) &&
					 (atoi((const char*) vIpAddress[2]) == 205) 
					)
				val = 1;

			// Tests sur les nouvelles adresses IP du site de Tokyo
			else if (
					 (atoi((const char*) vIpAddress[0]) == 10) &&
					 (atoi((const char*) vIpAddress[1]) == 19) &&
					 (atoi((const char*) vIpAddress[2]) >= 1) &&
					 (atoi((const char*) vIpAddress[2]) <= 4)
					)
				val = 1;

			// Tests sur les adresses IP du site de Milan
			else if (
					 (atoi((const char*) vIpAddress[0]) == 10) &&
					 (atoi((const char*) vIpAddress[1]) == 9) &&
					 (atoi((const char*) vIpAddress[2]) == 194) 
					)
				val = 1;

			// Tests sur les nouvelles adresses IP du site de Milan
			else if (
					 (atoi((const char*) vIpAddress[0]) == 10) &&
					 (atoi((const char*) vIpAddress[1]) == 19) &&
					 (atoi((const char*) vIpAddress[2]) >= 160) &&
					 (atoi((const char*) vIpAddress[2]) <= 167)
					)
				val = 1;

			// Tests sur les adresses IP du site de Madrid
			else if (
					 (atoi((const char*) vIpAddress[0]) == 10) &&
					 (atoi((const char*) vIpAddress[1]) == 9) &&
					 (atoi((const char*) vIpAddress[2]) == 198) 
					)
				val = 1;

			// Tests sur les nouvelles adresses IP du site de Madrid
			else if (
					 (atoi((const char*) vIpAddress[0]) == 10) &&
					 (atoi((const char*) vIpAddress[1]) == 19) &&
					 (atoi((const char*) vIpAddress[2]) >= 168) &&
					 (atoi((const char*) vIpAddress[2]) <= 175)
					)
				val = 1;

			// Tests sur les adresses IP du site de NY1
			else if (
					 (atoi((const char*) vIpAddress[0]) == 10) &&
					 (atoi((const char*) vIpAddress[1]) == 93) &&
					 (atoi((const char*) vIpAddress[2]) >= 80) &&
					 (atoi((const char*) vIpAddress[2]) <= 83)
					)
				val = 1;
		
			// Tests sur les adresses IP du site de NY2
			else if (
					 (atoi((const char*) vIpAddress[0]) == 10) &&
					 (atoi((const char*) vIpAddress[1]) == 96) &&
					 (atoi((const char*) vIpAddress[2]) >= 80) &&
					 (atoi((const char*) vIpAddress[2]) <= 83)
					)
				val = 1;
		}

		// Test sur le fichier de license
		if (val == 0)
			val = GetLicense();
	}
	else
		val = GetLicense();

	return val;
}


void LOCALARM_IniFileRead(void)
{
	armlocal_init = new ARMLOCAL_Init(XLLOCALARM_INIFILE_NAME);

#ifndef _MT
    MSG_printf_message(MSG_INFO, "\n======> Begin: [ARMLOCAL_IniFileRead ()] \n");
#endif

    try
	{
		armlocal_init->Init ();
	}

	catch(CInitReaderException*)
	{
		armlocal_init->default_currency    = DEFAULT_CURRENCY;
		armlocal_init->data_folder         = DATA_FOLDER;

#ifndef _MT
		MSG_printf_message(MSG_INFO, " EXCEPTION : ARM_IniFileRead (): catch() %s \n", armlocal_init->GetFileName ());
	
        MSG_printf_message(MSG_INFO, "             [ARM_IniFileRead ()]: catch() %s", (const char *) (*pEx));

        MSG_printf_message(MSG_INFO, "             [ARM_IniFileRead ()]: catch() %s\n", armlocal_init->GetFileName ());
#endif
	}

#ifndef _MT
    MSG_printf_message(MSG_INFO, "\n<====== End: [ARM_IniFileRead ()]\n");
#endif

	ARM_result C_result;

	long retCode = ARMLOCAL_ARM_SetDefaultCurrency ((const char*)(armlocal_init->default_currency.c_str()), C_result);
}


void FreeCell (const CCString& stringId)
{
	LOCAL_PERSISTENT_OBJECTS->FreeObject(LocalGetNumObjectId (stringId));
}


void FreeCellObject (const CCString& stringId)
{
	if(!stringId)
	{
		return;
	}
	else
	{
		FreeCell (stringId);
	}
}


#ifndef VBARM_LOCAL_EXPORTS


long LOCALARM_PersistentListsInit ()
{
	CREATE_GLOBAL_OBJECT();

	if(ARM_OBJECTS_LIST)
	{
		delete ARM_OBJECTS_LIST;
		ARM_OBJECTS_LIST = NULL;
	}

	ARM_OBJECTS_LIST = new ARM_OBJECTS_LIST_t ();
	
	MSG_printf_message (MSG_INFO, "Init: ObjectsList created");

	if(ARM_ERRORS_LIST)
	{
		delete ARM_ERRORS_LIST;
		ARM_ERRORS_LIST = NULL;
	}
	
	ARM_ERRORS_LIST = new ARM_ERRORS_LIST_t ();
	
	MSG_printf_message (MSG_INFO, "Init: ErrorsList created");

	return ARM_OK;
}



long LOCALARM_PersistentListsClear ()
{
	if(ARM_OBJECTS_LIST)
	{
		ARM_OBJECTS_LIST_t::iterator i = ARM_OBJECTS_LIST->begin ();
		ARM_OBJECTS_LIST_t::iterator e = ARM_OBJECTS_LIST->end ();
	
		while(i != e)
		{
            const char* toto = (*i).first;
            free((void*)toto);

			ARM_message* objectId = (*i).second;
			if(objectId)
			{
				delete objectId;
			}
			++i;
		}

		ARM_OBJECTS_LIST->clear ();
	}
	
	MSG_printf_message (MSG_INFO, "Init: ObjectsList cleared");

	if(ARM_ERRORS_LIST)
	{
		ARM_ERRORS_LIST_t::iterator i = ARM_ERRORS_LIST->begin ();
		ARM_ERRORS_LIST_t::iterator e = ARM_ERRORS_LIST->end ();
	
		while(i != e)
		{
            const char* toto = (*i).first;
            free((void*)toto);

			ARM_message* objectId = (*i).second;
			if(objectId)
			{
				delete objectId;
			}
			++i;
		}

		ARM_ERRORS_LIST->clear ();
	}
	
	MSG_printf_message (MSG_INFO, "Init: ErrorsList cleared");

	return ARM_OK;
}


long LOCALARM_PersistentListsDelete ()
{

	LOCALARM_PersistentListsClear ();
		
	if (LOCAL_PERSISTENT_OBJECTS)
	{
		delete LOCAL_PERSISTENT_OBJECTS;
		LOCAL_PERSISTENT_OBJECTS = NULL;
	}

	if(ARM_OBJECTS_LIST)
	{
		delete ARM_OBJECTS_LIST;
		ARM_OBJECTS_LIST = NULL;
	}
	
	MSG_printf_message (MSG_INFO, "Init: ObjectsList deleted");

	if(ARM_ERRORS_LIST)
	{
		delete ARM_ERRORS_LIST;
		ARM_ERRORS_LIST = NULL;
	}
	
	MSG_printf_message (MSG_INFO, "Init: ErrorsList deleted");

	return ARM_OK;
}


long LOCALARM_DeconnexionEToolkit ()
{
	return deconnection_etoolkit();
}




CCString GetCurStringCellCoordinates ()
{
	CCString caller = XL_getCaller ();

	return (caller);
}



void SetCurCellErrValue (const CCString& errValue)
{
	CCString caller = GetCurStringCellCoordinates ();

	if (strcmp((const char*)caller,"caller unknown") != 0)
	{
		if(ARM_ERRORS_LIST)
		{
			if(ARM_ERRORS_LIST->count ((const char*)caller) > 0)
			{
				ARM_ERRORS_LIST_t::iterator i = ARM_ERRORS_LIST->find ((const char*)caller);
				ARM_message* err = (*i).second;

				if(err)
				{
					delete err;
					err = NULL;
				}
				(*ARM_ERRORS_LIST)[(const char*)caller] = BuildErrorMessage (errValue,true);
			}
			else
			{
				int size	= ARM_ERRORS_LIST->size();
				(*ARM_ERRORS_LIST)[(const char*)strdup ((const char*)caller)] = BuildErrorMessage (errValue,true);
			}
		}
	}
	else
	{
		ARM_message* err = BuildErrorMessage (errValue,true);
		delete err;
	}
}


void UnsetCurCellValEnv ()
{
	CCString caller = GetCurStringCellCoordinates ();

	if (strcmp((const char*)caller,"caller unknown") != 0)
	{
		ARM_OBJECTS_LIST_t::iterator i = ARM_OBJECTS_LIST->find ((const char*)caller);
		ARM_message* err = (*i).second;
		ARM_OBJECTS_LIST->erase ((const char*)caller);
		if(err)
		{
			delete err;
		}
	}
}




CCString GetLastCurCellEnvValue ()
{
	CCString caller = GetCurStringCellCoordinates ();
	return (GetEnvVar (caller));
}


long FreeCurCellContent ()
{
	CCString envVal = GetLastCurCellEnvValue ();

	if(!envVal)
	{
		return ARM_OK;
	}
	else
	{
		UnsetCurCellValEnv ();
		FreeCellObject (envVal);
	}

	return ARM_OK;
}


long FreeCurCellContentWithoutFreeingARMObject()
{
	CCString envVal = GetLastCurCellEnvValue ();

	if(!envVal)
	{
		return ARM_OK;
	}
	else
	{
		UnsetCurCellValEnv ();
	}

	return ARM_OK;
}


void FreeCurCellErr ()
{
}


void PrintErrorsList ()
{
	ARM_ERRORS_LIST_t::iterator i = ARM_ERRORS_LIST->begin ();
	ARM_ERRORS_LIST_t::iterator e = ARM_ERRORS_LIST->end ();
	
	while(i != e)
	{
		MSG_printf_message (MSG_TRACE, "ARM_ERRORS_LIST[%s] = %s", (*i).first, (const char*)(*i).second->getMsg ());
		++i;
	}
}



CCString GetErrValue (const CCString& cellCoord)
{
	CCString lastCurCellErrValue;
	
	if(ARM_ERRORS_LIST)
	{
		PrintErrorsList ();
		if(ARM_ERRORS_LIST->count ((const char*)cellCoord) > 0)
		{
			ARM_message* error = (*ARM_ERRORS_LIST)[(const char*)cellCoord];
			if(error)
				return error->getMsg();
		}
	}

	lastCurCellErrValue.Set ("NO CORRESPONDING ERROR MESSAGE FOUND !!!");

	return (lastCurCellErrValue);
}


void  LocalSetCurCellEnvValue (const CCString& curClass, long objId)
{
	CCString caller = GetCurStringCellCoordinates ();

	if (strcmp((const char*)caller,"caller unknown") != 0)
	{

		CCString stringObjectId = LocalMakeObjectId (objId, curClass);

		if(ARM_OBJECTS_LIST)
		{
			if(ARM_OBJECTS_LIST->count ((const char*)caller) > 0)
			{
				ARM_OBJECTS_LIST_t::iterator i = ARM_OBJECTS_LIST->find ((const char*)caller);
				ARM_message* err = (*i).second;
				//ARM_OBJECTS_LIST->erase ((const char*)caller);
				if(err)
				{
					delete err;
				}
				(*ARM_OBJECTS_LIST)[(const char*)caller] = BuildErrorMessage (stringObjectId,false);	
			}
			else
			{
				(*ARM_OBJECTS_LIST)[(const char*)strdup ((const char*)caller)] = BuildErrorMessage (stringObjectId,false);
			}
		}
	}
}




#endif



/// Ne sert pas qu'à construire les messages d'erreur
ARM_message* BuildErrorMessage (const CCString& errValue, bool storeInLog)
{
	ARM_message* msg = new ARM_message (errValue);

	if(!msg)
	{
		return NULL;
	}
	else
	{
		/// get the message
		std::string msgString = CCSTringToSTLString(msg->getMsg());

		if (storeInLog)
		{
			/// and store it in the log file
			char username[50];
			char filename[254];

			DWORD nbChar = sizeof(username);

			GetUserName(username, &nbChar);

            char fOutName[200];

            ARM_GetTmpAbsFile("LogAllErrors_", fOutName);

			strcpy(filename, fOutName);
			strcat(filename, username);
			strcat(filename, ".txt");

			FILE* errorLog = fopen(filename, "a" );

			fprintf(errorLog , "%s\n\n", msgString.c_str() );
			
            fclose(errorLog);
		}
	}

	return(msg);
}



void UnsetCurCellValEnv (const CCString & caller)
{
	//CCString caller = GetCurStringCellCoordinates ();

	ARM_OBJECTS_LIST_t::iterator i = ARM_OBJECTS_LIST->find ((const char*)caller);
	ARM_message* err = (*i).second;
	ARM_OBJECTS_LIST->erase ((const char*)caller);
	if(err)
	{
		delete err;
	}
}


void LocalSetCurCellEnvValue (const CCString& curClass, long objId,const CCString& caller)
{
	//CCString caller = GetCurStringCellCoordinates ();

	CCString stringObjectId = LocalMakeObjectId (objId, curClass);

	if(ARM_OBJECTS_LIST)
	{
		if(ARM_OBJECTS_LIST->count ((const char*)caller) > 0)
		{
			ARM_OBJECTS_LIST_t::iterator i = ARM_OBJECTS_LIST->find ((const char*)caller);
			ARM_message* err = (*i).second;
			//ARM_OBJECTS_LIST->erase ((const char*)caller);
			if(err)
			{
				delete err;
			}
			(*ARM_OBJECTS_LIST)[(const char*)caller] = BuildErrorMessage (stringObjectId,false);	
		}
		else
		{
			(*ARM_OBJECTS_LIST)[(const char*)strdup ((const char*)caller)] = BuildErrorMessage (stringObjectId,false);
		}
	}
}


CCString GetEnvVar (const CCString& cellCoord)
{
	CCString lastCurCellValue;

	if(ARM_OBJECTS_LIST)
	{
		if(ARM_OBJECTS_LIST->count ((const char*)cellCoord) > 0)
		{
			ARM_message* error = (*ARM_OBJECTS_LIST)[(const char*)cellCoord];
			if(error)
			{
				lastCurCellValue = error->getMsg ();
			}
		}
	}

	return (lastCurCellValue);
}



long FreeCurCellContent (const CCString & envVal)
{
	//CCString envVal = GetLastCurCellEnvValue ();

	if(!envVal)
	{
		return ARM_OK;
	}
	else
	{
		UnsetCurCellValEnv (envVal);
		FreeCellObject (envVal);
	}

	return ARM_OK;
}


CCString LocalGetStringObjectClass (const CCString& stringObjectId)
{
	char buf[6];

	strncpy (buf, (const char*)stringObjectId, 5);
	buf[5] = '\0';

	CCString stringObjectClass (buf);

	return (stringObjectClass);
}


CCString LocalMakeObjectId (long objId, const CCString& objClass)
{
	// TMP: objId a tester
	char buf[50];
	char sHour[4], sMin[4], sSec[4];

    struct tm *newtime;
    time_t long_time;

    time( &long_time );                /* Get time as long integer. */
    newtime = localtime( &long_time ); /* Convert to local time. */
	
	if (newtime->tm_hour < 10)
		sprintf (sHour, "0%1d", newtime->tm_hour);
	else
		sprintf (sHour, "%2d", newtime->tm_hour);
	
	if (newtime->tm_min < 10)
		sprintf (sMin, "0%1d", newtime->tm_min);
	else
		sprintf (sMin, "%2d", newtime->tm_min);
	
	if (newtime->tm_sec < 10)
		sprintf (sSec, "0%1d", newtime->tm_sec);
	else
		sprintf (sSec, "%2d", newtime->tm_sec);

	sprintf (buf, "%ld %s:%s:%s", objId,sHour,sMin,sSec);

	CCString stringObjectId = objClass + "_" + buf;

	return (stringObjectId);
}


//////////////////////////////////////////////
/// Version to handle reading of previous class
//////////////////////////////////////////////

CCString LocalGetStringObjectClassAndError(const CCString& stringObjectId)
{
	if(		!stringObjectId
		|| (stringObjectId == "ARM_ERR")
		|| (stringObjectId == "ARM_FATAL_ERR")
		|| (stringObjectId == "OBSOLETE")
		|| (stringObjectId == "ERROR")
		)
	{
		return CCString(" ");
	}

	CCString curClass = LocalGetStringObjectClass (stringObjectId);
	return curClass;
}



//	if the string is empty (means: the cell was previously empty)
//	then returns -1 
long LocalGetNumObjectId (const CCString& stringObjectId)
{
	/// we suppose that any ARM Object start by L
	/// and is followed by 4 capital letters or numbers characters followed by '_'
	/// if this is not the case, return an error

	if ( (const char*)stringObjectId == NULL )
	{
		return ARM_KO;
    }

	if ( stringObjectId[0] != 'L' )
	{
		return ARM_KO;
    }
	if (stringObjectId.GetLen() <6) 
		return ARM_KO ;
	for ( int i = 1; i < 5; ++i)
	{
		if ((( stringObjectId[i] < 'A' ) || ( stringObjectId[i] > 'Z' ))
			&&  
			(( stringObjectId[i] < '0' ) || ( stringObjectId[i] > '9' ))
		   )
		{
		   return ARM_KO;
		}
	}
	
	if ( stringObjectId[5] != '_' )
	{
		return ARM_KO;
    }

	char buf[50];
	const char* aString = (const char *) stringObjectId;
	strncpy (buf, &(aString[6]),49);buf[49]=0;

	if ( strlen (buf) == 0 )
	{
	   return ARM_KO;
	}

	char strId[50];

	// API documentation: equivalent to sscanf(buf, "%s", strId);
	// "%[a-z] : read all characters between a and z , then stop
	// "%[^ ] : read all characters that are not "space", then end.. 
	sscanf(buf, "%[^ ]", strId);

	long ret = atol(strId);
	if (ret==0 && strcmp(strId,"0")!=0) return ARM_KO; 
	return ret ;
}


long ARM_ConvPriceYield (const CCString& aPriceYield, ARM_result& result)
{
	CCString tmp = aPriceYield;
	tmp.toUpper ();

	if(tmp == "PRICE")
	{
		return K_PRICE;
	}
	if(tmp == "P")
	{
		return K_PRICE;
	}
	if(tmp == "0")
	{
		return K_PRICE;
	}
	if(tmp == "YIELD")
	{
		return K_YIELD;
	}
	if(tmp == "Y")
	{
		return K_YIELD;
	}
	if(tmp == "1")
	{
		return K_YIELD;
	}
    
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Parameter - Valid are Yield, Price, Y, P, 0, 1");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvIrIndName (const CCString& aIrIndName, ARM_result& result)
{
	CCString tmp = aIrIndName;
	tmp.toUpper ();

    if ( tmp == "FIXED" )
	{
	   return K_FIXED;
	}
	if(tmp == "PIBOR1M")
	{
		return K_PIBOR1M;
	}
	if(tmp == "PIBOR3M")
	{
		return K_PIBOR3M;
	}
	if(tmp == "PIBOR6M")
	{
		return K_PIBOR6M;
	}
	if(tmp == "PIBOR1Y")
	{
		return K_PIBOR1Y;
	}
	if(tmp == "LIBOR1M")
	{
		return K_LIBOR1M;
	}
	if(tmp == "LIBOR3M")
	{
		return K_LIBOR3M;
	}
	if(tmp == "LIBOR6M")
	{
		return K_LIBOR6M;
	}
	if(tmp == "LIBOR12M")
	{
		return K_LIBOR1Y;
	}
	if(tmp == "LIBOR1Y")
	{
		return K_LIBOR1Y;
	}
	if(tmp == "EURIBOR1M")
	{
		return K_EURIBOR1M;
	}
	if(tmp == "EURIBOR3M")
	{
		return K_EURIBOR3M;
	}
	if(tmp == "EURIBOR6M")
	{
		return K_EURIBOR6M;
	}
	if(tmp == "EURIBOR1Y")
	{
		return K_EURIBOR1Y;
	}
	if(tmp == "EURIBOR12M")
	{
		return K_EURIBOR1Y;
	}
	if(tmp == "LIVRETA")
	{
		return K_LIVRET_A;
	}
	if(tmp == "EUR3M")
	{
		return K_EUR3M;
	}
	if(tmp == "EUR12")
	{
		return K_EUR12;
	}
	if(tmp == "EUR1M")
	{
		return K_EUR1M;
	}
	if(tmp == "EURIB")
	{
		return K_EURIBOR6M;
	}
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Type Name - Valid are Fixed, Euribor1m, Euribor3m, Euribor6m, Euribor12m, Euribor1y, Libor1m, Libor3m, Libor6m, Libor12m, Libor1y, Eur1m, Eurib");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvIrIndNameToFreq (const CCString& aIrIndName, ARM_result& result)
{
	CCString tmp = aIrIndName;
	tmp.toUpper ();

	if(tmp == "PIBOR1M")
	{
		return K_MONTHLY;
	}
	if(tmp == "PIBOR3M")
	{
		return K_QUARTERLY;
	}
	if(tmp == "PIBOR6M")
	{
		return K_SEMIANNUAL;
	}
	if(tmp == "PIBOR12M")
	{
		return K_ANNUAL;
	}
	if(tmp == "PIBOR1Y")
	{
		return K_ANNUAL;
	}
	if(tmp == "LIBOR1M")
	{
		return K_MONTHLY;
	}
	if(tmp == "LIBOR3M")
	{
		return K_QUARTERLY;
	}
	if(tmp == "LIBOR6M")
	{
		return K_SEMIANNUAL;
	}
	if(tmp == "LIBOR12M")
	{
		return K_ANNUAL;
	}
	if(tmp == "LIBOR1Y")
	{
		return K_ANNUAL;
	}
	if(tmp == "EURIBOR1M")
	{
		return K_MONTHLY;
	}
	if(tmp == "EURIBOR3M")
	{
		return K_QUARTERLY;
	}
	if(tmp == "EURIBOR6M")
	{
		return K_SEMIANNUAL;
	}
	if(tmp == "EURIBOR12M")
	{
		return K_ANNUAL;
	}
	if(tmp == "EURIBOR1Y")
	{
		return K_ANNUAL;
	}
	if(tmp == "EURIB")
	{
		return K_SEMIANNUAL;
	}
	if(tmp == "EUR1M")
	{
		return K_MONTHLY;
	}
	
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Type Name - Valid are Euribor1m, Euribor3m, Euribor6m, Euribor12m, Euribor1y, Libor1m, Libor3m, Libor6m, Libor12m, Libor1y, Eur1m, Eurib");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvTmIxName (const CCString& aTmIxName, ARM_result& result)
{ 
	CCString tmp = aTmIxName;
	tmp.toUpper ();

    if(tmp == "T4M")
	{
         return 26;
	}
    if(tmp == "TAM")
	{
         return 28;
	}
    if(tmp == "TMP")
	{
         return 30;
	}
    if(tmp == "T4M_FIXED")
	{
         return 27;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Type Name - Valid are T4M, TAM, TMP, T4M_FIXED");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvYieldOrVol (const CCString& yieldOrVol, ARM_result& result)
{ 
	CCString tmp = yieldOrVol;
	tmp.toUpper ();

    if(tmp == "Y")
	{
         return 0;
	}
    if(tmp == "V")
	{
         return 1;
	}
    
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Type Name - Valid are Y, V");

	return ARM_DEFAULT_ERR;
}

long ARM_ConvInfSwapTypeCommon(const CCString& input, ARM_result& result, bool allowGeneric = false )
{
	CCString tmp = input;
	tmp.toUpper();

    if(tmp == "ZC")
         return K_ZEROCOUPON_LEG;
    if(tmp == "YTY")
         return K_YEARTOYEAR_LEG;
    if(tmp == "OAT")
         return K_OATTYPE_LEG;

    if( allowGeneric  )
		if(tmp == "GENERIC")
			return K_GENERICINF_LEG;
    
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Type Name - Valid are ZC, YTY, OAT, GENERIC");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvInfSwapType(const CCString& input, ARM_result& result )
{
	return ARM_ConvInfSwapTypeCommon( input, result, false );
}


long ARM_ConvCalcMod (const CCString& calcMod, ARM_result& result)
{ 
	CCString tmp = calcMod;
	tmp.toUpper ();

    if(tmp == "NOR")
	{
         return 0;
	}
    if(tmp == "LOGNOR")
	{
         return 1;
	}
	if(tmp == "SQR")
	{
         return 2;
	}
    
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Type Name - Valid are NOR, LOGNOR, SQR");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvFrequency (const CCString& aFrequency, ARM_result& result)
{ 
	CCString tmp = aFrequency;
	tmp.toUpper ();

	if(tmp == "-1")
	{
		return K_DEF_FREQ;
	}
    if( (tmp == "A") || (tmp == "Y") )
	{
        return K_ANNUAL;
	}
    if(tmp == "S")
	{
        return K_SEMIANNUAL;
	}
    if(tmp == "Q")
	{
        return K_QUARTERLY;
	}
    if(tmp == "B")
	{
        return K_BIMONTHLY;
	}
    if(tmp == "M")
	{
        return K_MONTHLY;
	}
    if(tmp == "W")
	{
		return K_WEEKLY;
	}
    if(tmp == "D")
	{
		return K_DAILY;
	}
    if(tmp == "Z")
	{
        return K_ZEROCOUPON;
	}
    if(tmp == "1")
	{
        return K_ANNUAL;
	}
    if(tmp == "2")
	{
		return K_SEMIANNUAL;
	}
    if(tmp == "4")
	{
        return K_QUARTERLY;
	}
    if(tmp == "6")
	{
        return K_BIMONTHLY;
	}
    if(tmp == "12")
	{
		return K_MONTHLY;
	}
    if(tmp == "52")
	{
        return K_WEEKLY;
	}
    if(tmp == "365")
	{
        return K_DAILY;
	}
    if(tmp == "0")
	{
         return K_ZEROCOUPON;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Type Name - Valid are A, S, Q, B(imonth), M, D, Z");

	return ARM_DEFAULT_ERR;
}




 



long ARM_ConvCMIndName (const CCString& aIrType, ARM_result& result)
{
	CCString tmp = aIrType;
	tmp.toUpper ();

   
   
    if(tmp == "CMS1")
	{
        return K_CMS1;
	}
    if(tmp == "CMS2")
	{
        return K_CMS2;
	}
    if(tmp == "CMS3")
	{
        return K_CMS3;
	}
    if(tmp == "CMS4")
	{
        return K_CMS4;
	}
    if(tmp == "CMS5")
	{
        return K_CMS5;
	}
    if(tmp == "CMS6")
	{
        return K_CMS6;
	}   
    if(tmp == "CMS7")
	{
        return K_CMS7;
	}   
    if(tmp == "CMS8")
	{
        return K_CMS8;
	}   
    if(tmp == "CMS9")
	{
        return K_CMS9;
	}   
 
    if(tmp == "CMS10")
	{
        return K_CMS10;
	}
    if(tmp == "CMS11")
	{
        return K_CMS11;
	}
    if(tmp == "CMS12")
	{
        return K_CMS12;
	}
    if(tmp == "CMS13")
	{
        return K_CMS13;
	}
    if(tmp == "CMS14")
	{
        return K_CMS14;
	}
    if(tmp == "CMS15")
	{
        return K_CMS15;
	}
    if(tmp == "CMS16")
	{
        return K_CMS16;
	}
    if(tmp == "CMS17")
	{
        return K_CMS17;
	}
    if(tmp == "CMS18")
	{
        return K_CMS18;
	}
    if(tmp == "CMS19")
	{
        return K_CMS19;
	}
    if(tmp == "CMS20")
	{
        return K_CMS20;
	}
    if(tmp == "CMS21")
	{
        return K_CMS21;
	}
    if(tmp == "CMS22")
	{
        return K_CMS22;
	}
    if(tmp == "CMS23")
	{
        return K_CMS23;
	}
    if(tmp == "CMS24")
	{
        return K_CMS24;
	}
    if(tmp == "CMS25")
	{
        return K_CMS25;
	}
    if(tmp == "CMS26")
	{
        return K_CMS26;
	}
    if(tmp == "CMS27")
	{
        return K_CMS27;
	}
    if(tmp == "CMS28")
	{
        return K_CMS28;
	}
    if(tmp == "CMS29")
	{
        return K_CMS29;
	}
    if(tmp == "CMS30")
	{
        return K_CMS30;
	}
   
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Type index - Valid are CMS1, CMS2, ... CMS10 etc ...");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvCompMeth (const CCString& aCompMeth, ARM_result& result)
{
	CCString tmp = aCompMeth;
	tmp.toUpper ();

	if(tmp == "C")
	{
        return K_COMP_CONT;
	}
    if(tmp == "P")
	{
        return K_COMP_PROP;
	}
    if(tmp == "A")
	{
        return K_COMP_ANNUAL;
	}
    if(tmp == "S")
	{
        return K_COMP_SEMIANNUAL;
	}
    if(tmp == "Q")
	{
        return K_COMP_QUARTERLY;
	}
    if(tmp == "B")
	{
        return K_COMP_BIMONTHLY;
	}
    if(tmp == "M")
	{
        return K_COMP_MONTHLY;
	}
    if(tmp == "D360")
	{
        return K_COMP_DAILY_360;
	}
    if(tmp == "D365")
	{
        return K_COMP_DAILY_365;
	}
    if(tmp == "0")
	{
        return K_COMP_CONT;
	}
    if(tmp == "-1")
	{
        return K_COMP_PROP;
	}
    if(tmp == "1")
	{
        return K_COMP_ANNUAL;
	}
    if(tmp == "2")
	{
        return K_COMP_SEMIANNUAL;
	}
    if(tmp == "4")
	{
        return K_COMP_QUARTERLY;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Compounding Method - Valid are: 'C'ontinuous, 'P'roportional (Lineaire), 'A', 'S', 'Q', 'B', 'M', 'D360', 'D365'");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvForwardYieldMethod (const CCString& aForwardYieldMeth, ARM_result& result)
{
	CCString tmp = aForwardYieldMeth;
	tmp.toUpper ();

	if(tmp == "MON")
	{
        return -1;
	}
	if(tmp == "-1")
	{
		return -1;
	}
    if(tmp == "CONT")
	{
        return 0;
	}
	if(tmp[0] == '0')
	{
		return 0;
	}
	if(tmp == "ACTU")
	{
        return 1;
	}
	if(tmp[0] == '1')
	{
		return 1;
	}
	
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Forward Yield Method - Valid are: 'MON', '-1', 'CONT', '0', 'ACTU', '1'");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvCallOrPut (const CCString& aCorP, ARM_result& result)
{
	CCString tmp = aCorP;
	tmp.toUpper ();

    if(tmp == "CALL")
	{
        return K_CALL;
	}
    if(tmp == "PUT")
	{
        return K_PUT;
	}
    if(tmp == "C")
	{
        return K_CALL;
	}
    if(tmp == "P")
    {
		return K_PUT;
	}
    if(tmp == "1")
	{
        return K_CALL;
	}
    if(tmp == "-1")
	{
        return K_PUT;
	}
 
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Option Type - Valid are: (Call, C, 1) Or (Put, P, -1)");

	return ARM_DEFAULT_ERR;
}

long ARM_ConvKO_Or_NoKO (const CCString& param, ARM_result& result)
{
	CCString tmp = param;
	//tmp.toUpper ();

    if(tmp == "KO")
	{
        return K_KO;
	}
    if(tmp == "No_KO_ACC")
	{
        return K_No_KO_ACC;
	}
	if(tmp == "No_KO_No_ACC")
	{
        return K_No_KO_No_ACC;
	}
    if(tmp == "1")
	{
        return K_KO;
	}
    if(tmp == "0")
    {
		return K_No_KO_ACC;
	}
    if(tmp == "-1")
	{
        return K_No_KO_No_ACC;
	}
 
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Option Type - Valid are: (KO, 1) Or (No_KO_ACC, 0) or (No_KO_No_ACC, -1)");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvKOStyle (const CCString& param, ARM_result& result)
{
	qKoStyle ret; 
	ICM_EnumsCnv::cnv(CCSTringToSTLString(param),ret); 
	return ret ; 
}

long ARM_ConvAccStyle (const CCString& param, ARM_result& result)
{
	qAccelerationStyle ret; 
	ICM_EnumsCnv::cnv(CCSTringToSTLString(param),ret); 
	return ret ; 
}


long ARM_ConvParam (const CCString& param, ARM_result& result)
{
	CCString tmp = param;
	tmp.toUpper ();

    if(tmp == "PRICE")
	{
        return 0;
	}
    if(tmp == "DELTA")
	{
        return 1;
	}
    
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Param - Valid are: PRICE, DELTA");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvMktType (const CCString& aMkt, ARM_result& result)
{
	CCString tmp = aMkt;
	tmp.toUpper ();

    if( (tmp == "F") || (tmp == "FUT") || (tmp[0] == 'F') )
	{
        return K_FUT;
	}
    if( (tmp == "M") || (tmp == "MM") || (tmp[0] == 'M') )
	{
        return K_MM;
	}
    if( (tmp == "S") || (tmp == "SWAP") || (tmp[0] == 'S') )
	{
        return K_SWAP;
	}
    if( (tmp == "B") || (tmp == "BOND")  || (tmp[0] == 'B'))
	{
        return K_BOND;
	}
    if( (tmp == "0") || (tmp[0] == '0') )
	{
        return K_FUT;
	}
    if( (tmp == "1") || (tmp[0] == '1') ) 
	{
        return K_MM;
	}
    if( (tmp == "2") || (tmp[0] == '2') ) 
	{
        return K_SWAP;
	}
    if( (tmp == "3") || (tmp[0] == '3') ) 
	{
        return K_BOND;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Market - Valid are: (FUT,0) - (MM,1) - (SWAP,2) - (BOND,3)");

	return ARM_DEFAULT_ERR;
}


/** 
long ARM_Conv_Sectorial_Correlation (const CCString& aIntMeth, ARM_result& result)
{
	CCString tmp = aIntMeth;
	tmp.toUpper ();



    if ((tmp == "FULL") || (tmp == "DIFF_INTER_DIFF_INTRA"))
	{
		return TFCT_FULL;
	}

    if (tmp == "SAME_INTER_DIFF_INTRA")
	{
        return TFCT_SAME_INTER_DIFF_INTRA;
	}

    if (tmp == "SAME_INTER_SAME_INTRA")
	{
        return TFCT_SAME_INTER_SAME_INTRA;
	}

    if( (tmp == "0") || (tmp[0] == '0') )
	{
        return TFCT_FULL;
	}
    if( (tmp == "1") || (tmp[0] == '1') ) 
	{
        return TFCT_SAME_INTER_DIFF_INTRA;
	}
    if( (tmp == "2") || (tmp[0] == '2') ) 
	{
        return TFCT_SAME_INTER_SAME_INTRA;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Sectorial Correlation - Valid are: (DIFF_INTER_DIFF_INTRA, 0) (SAME_INTER_DIFF_INTRA, 1) (SAME_INTER_SAME_INTRA, 2");

	return ARM_DEFAULT_ERR;
}
**/ 

long ARM_ConvInterpMethod (const CCString& Type, ARM_result& result)
{
	CCString tmp = Type;
	tmp.toUpper ();

	if (tmp == "FIX_RIGHT_MATU")
	{
		return K_ICM_STEPUP_RIGHT_MATURITY;
	}

	if (tmp == "FIX_LEFT_MATU")
	{
		return K_ICM_STEPUP_LEFT_MATURITY;
	}

    if (( tmp[0] == 'C' ) || ( tmp[0] == '0' ))
	{
		return K_CONTINUOUS;
	}

    if (( tmp[0] == 'L' ) || ( tmp[0] == '1' ))
	{
        return K_LINEAR;
	} 

    if (( tmp[0] == 'S' ) || ( tmp[0] == '2' ))
	{
        return K_SPLINE;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Interpol Method - Valid are: (CONT,0) (LIN,1) (SPL, 2");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvWhatIsInterp (const CCString& aIntMeth, ARM_result& result)
{
	CCString tmp = aIntMeth;
	tmp.toUpper ();

	if(tmp==CCString("DEFAULT"))
	{
		return FXINTERP_STRIKE;
	}

    if(tmp[0] == 'S')
	{
		return FXINTERP_STRIKE;
	}
    if(tmp[0] == 'D')
	{
        return FXINTERP_DELTA;
	}
    if(tmp[0] == 'E')
	{
        return FXINTERP_DELTA_SMOOTH;
	}
    if(tmp[0] == '0')
	{
        return FXINTERP_STRIKE;
	}
    if(tmp[0] == '1')
	{
        return FXINTERP_DELTA;
	}
    if(tmp[0] == '2')
	{
        return FXINTERP_DELTA_SMOOTH;
	}
	
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid what is interpolated - Valid are: (S,0) - (D,1) - (E,2)");

	return ARM_DEFAULT_ERR;
}
/*MICHAEL*/


long ARM_ConvCPInterpMethod( const CCString& InfMethod, ARM_result& result, long Flag )
{
	CCString tmp = InfMethod;
	tmp.toUpper();

	switch( Flag )
	{ 
		/// all flags allowed for daily interpolation
		case K_DAILY:
		{
			if( tmp == "STEPWISESTART" || tmp == "START" )
					return K_CPISTEPWISESTART;
			if( tmp == "STEPWISEMID" || tmp == "MID" )
				return K_CPISTEPWISEMIDDLE;
			if( tmp == "STEPWISEEND" || tmp == "END")
				return K_CPISTEPWISEEND;
			if( tmp == "STEPWISE")
				return K_CPISTEPWISE;
			if( tmp == "CPILINEAR")
				return K_CPILINEAR;
			if( tmp == "ZCLINEAR")
				return K_ZCLINEAR;
			if( tmp == "ZCCTFWD")
				return K_ZCCTFWD;
			if( tmp == "-1")
				return -1;

			// other cases are errors
			result.setRetCode (ARM_DEFAULT_ERR);
			result.setMsg ("ARM_ERR: Invalid Daily Interpol Method - Valid are: \nSTEPWISESTART, STEPWISEMID, STEPWISEEND, STEPWISE, CPILINEAR, ZCLINEAR, ZCCTFWD");
			return ARM_DEFAULT_ERR;
		}

		/// all flags allowed for monthly interpolation
		case K_MONTHLY:
		{
			if( tmp == "CPILINEAR" )
				return K_CPILINEAR;
			if( tmp == "ZCLINEAR")
				return K_ZCLINEAR;
			if( tmp == "ZCCTFWD")
				return K_ZCCTFWD;
			if( tmp == "-1")
				return -1;

			// other cases are errors
			result.setRetCode (ARM_DEFAULT_ERR);
			result.setMsg ("ARM_ERR: Invalid Daily Interpol Method - Valid are: \nCPILINEAR, ZCLINEAR, ZCCTFWD");
			return ARM_DEFAULT_ERR;
		}

		// other cases are errors
		default:
		{
			result.setRetCode (ARM_DEFAULT_ERR);
			result.setMsg ("ARM_ERR: Invalid Interpol Method - Valid are: \nMONTHLY or DAILY");
			return ARM_DEFAULT_ERR;
		}
	}
}


// String to long conversion for monthly interpolation for CPI
long ARM_ConvCPIMonthlyInterpMethod( const CCString& InfMethod, ARM_result& result)
{
	return ARM_ConvCPInterpMethod( InfMethod, result, K_MONTHLY );
}

// String to long conversion for daily interpolation for CPI
long ARM_ConvCPIDailyInterpMethod( const CCString& InfMethod, ARM_result& result)
{
	return ARM_ConvCPInterpMethod( InfMethod, result, K_DAILY );
}



long ARM_ConvCPIExtrapolMethod( const CCString& InfMethod, ARM_result& result)
{
	CCString tmp = InfMethod;
	tmp.toUpper();

	if( tmp  == "LASTTWO")
		return K_LASTTWO;
	if( tmp  == "MIDDLE")
		return K_MIDDLE;
	if( tmp  == "FIRSTTWO")
		return K_FIRSTTWO;
	
	if( tmp  == "-1")
		return -1;
	/// other cases are error
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Extrapol Method - Valid are: \nLASTTWO, MIDDLE, FIRSTTWO");
	return ARM_DEFAULT_ERR;
}


long ARM_Strike( const CCString& InfMethod, ARM_result& result)
{
	return ARM_ConvCPInterpMethod( InfMethod, result, K_DAILY );
}



long ARM_ConvCapOrFloor (const CCString& aCorF, ARM_result& result)
{ 
	CCString tmp = aCorF;
	tmp.toUpper ();


    if ( tmp == "CAP" )
	{
       return(K_CAP);
	}

    if ( tmp == "FLOOR" )
	{
       return(K_FLOOR);
	}

    if ( tmp == "CAPFLOOR" )
	{
       return(K_CAPFLOOR);
	}

    if ( tmp == "C" )
	{
       return(K_CAP);
	}

    if ( tmp == "F" )
	{
       return(K_FLOOR);
	}

    if ( tmp == "CF" )
	{
       return(K_CAPFLOOR);
	}

    if ( tmp == "1" )
	{
       return(K_CAP);
	}

    if ( tmp == "-1" )
	{
       return(K_FLOOR);
}

    if ( tmp == "0" )
	{
       return(K_CAPFLOOR);
	}

	result.setRetCode(ARM_DEFAULT_ERR);
	result.setMsg("ARM_ERR: Invalid CapOrFloor - Valid are : Cap, C, 1 or Floor, F, -1, CAPFLOOR, CF, 0");

	return(ARM_DEFAULT_ERR);
}



long ARM_ConvExerciseType(const CCString& aOptType, ARM_result& result)
{
	CCString tmp = aOptType;
	tmp.toUpper ();

    if(tmp == "E")
	{
		return K_EUROPEAN;
	}
    if(tmp == "A")
	{
        return K_AMERICAN;
	}
    if(tmp == "B")
	{
        return K_BERMUDAN;
	}
 
	
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Exercise Type - Valid are: A{merican}, E{uropean}, B{ermudan}");

	return ARM_DEFAULT_ERR;
}

long ARM_ConvRecOrPay (const CCString& aRorP, ARM_result& result)
{
	CCString tmp = aRorP;
	tmp.toUpper ();

    if(tmp == "R")
	{
        return K_RCV;
	}
    if(tmp == "P")
	{
		return K_PAY;
	}
    if(tmp == "REC")
	{
		return K_RCV;
	}
    if(tmp == "PAY")
	{
		return K_PAY;
	}
    if(tmp == "1")
	{
		return K_RCV;
	}
    if(tmp == "-1")
	{
		return K_PAY;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid RcvOrPay - Valid are REC, R, r, 1 or PAY, P, p, -1");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvLongOrShort (const CCString& aLorP, ARM_result& result)
{
	CCString tmp = aLorP;
	tmp.toUpper ();

    if(tmp == "L")
	{
        return K_RCV;
	}
    if(tmp == "S")
	{
		return K_PAY;
	}
    if(tmp == "LONG")
	{
		return K_RCV;
	}
    if(tmp == "SHORT")
	{
		return K_PAY;
	}
    if(tmp == "1")
	{
		return K_RCV;
	}
    if(tmp == "-1")
	{
		return K_PAY;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid LongOrShort - Valid are LONG, L, l, 1 or SHORT, S, s, -1");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvFwdRule (const CCString& aFwdRule, ARM_result& result)
{
	CCString tmp = aFwdRule;
	tmp.toUpper ();

    if(tmp == "F")
	{
        return K_FOLLOWING;
	}
    if(tmp == "MF")
	{
        return K_MOD_FOLLOWING;
	}
    if(tmp == "P")
	{
        return K_PREVIOUS;
	}
    if(tmp == "MP")
	{
        return K_MOD_PREVIOUS;
	}
    if(tmp == "1")
	{
        return K_FOLLOWING;
	}
    if(tmp == "2")
	{
        return K_MOD_FOLLOWING;
	}
    if(tmp == "-1")
	{
        return K_PREVIOUS;
	}
    if(tmp == "-2")
	{
        return K_MOD_PREVIOUS;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Fwd Rule - Valid are: F(ollowing) (1), M(odified)F(ollowing) (2),  P(revious) (-1), M((odified)P(revious) (-2) ");

	return ARM_DEFAULT_ERR;
}

long ARM_ConvRule (const CCString& rule, ARM_result& result)
{
	CCString tmp = rule;
	tmp.toUpper ();

    if(tmp == "F")
	{
        return 1;
	}
	if(tmp == "1")
	{
        return 1;
	}
    if(tmp == "P")
	{
        return -1;
	}
	if(tmp == "-1")
	{
        return -1;
	}
    if(tmp == "MF")
	{
        return 2;
	}
	if(tmp == "2")
	{
        return 2;
	}
    if(tmp == "MP")
	{
        return -2;
	}
	if(tmp == "-2")
	{
        return -2;
	}
    if(tmp == "NONE")
	{
        return 0;
	}
    
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Rule - Valid are: F, MF, P, MP or NONE");

	return ARM_DEFAULT_ERR;
}
 


 
long ARM_ConvObjectClass (const CCString& idUnderlying, ARM_result& result)
{
	CCString tmp = idUnderlying;
	tmp.toUpper ();

	if(tmp == "BOND")
	{
		return 1;
	}
	if(tmp == "PORT")
	{
		return 1;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Object Class");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvSvtyParam (const CCString& aParam, ARM_result& result)
{
	CCString tmp = aParam;
	tmp.toUpper ();

	if(tmp == "DELTA")
	{
		return K_DELTA;
	}
	if(tmp == "GAMMA")
	{
		return K_GAMMA;
	}
	if(tmp == "RHO")
	{
		return K_RHO;
	}
	if(tmp == "THETA")
	{
		return K_THETA;
	}
	if(tmp == "VEGA")
	{
		return K_VEGA;
	}
	if(tmp == "KAPPA")
	{
		return K_KAPPA;
	}

    if (( tmp == "CORREL") || ( tmp == "CORR" ))
    {
       return(K_CORREL);
    }

	if(tmp[0] == '0')
	{
		return K_DELTA;
	}
	if(tmp[0] == '1')
	{
		return K_GAMMA;
	}
	if(tmp[0] == '2')
	{
		return K_RHO;
	}
	if(tmp[0] == '3')
	{
		return K_THETA;
	}
	if(tmp[0] == '4')
	{
		return K_VEGA;
	}
	if(tmp[0] == '5')
	{
		return K_KAPPA;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvCurvSvtyParam (const CCString& aParam, ARM_result& result)
{
	CCString tmp = aParam;
	tmp.toUpper ();

	if ((tmp == "TAUX") || (tmp[0] == 'Y') || (tmp[0] == 'y')) 
	{
		return 10;
	}
	if((tmp == "VOL")|| (tmp[0] == 'V') || (tmp[0] == 'v'))
	{
		return 11;
	}
	long buf;
	if(sscanf (tmp, "%ld", &buf) == 1)
	{
		return (buf + 10);
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter");

	return ARM_DEFAULT_ERR;
}

long ARM_ConvUpDownDouble(const CCString& aUpDownDouble, ARM_result& result)
{
	CCString tmp = aUpDownDouble;
	tmp.toUpper ();

	if(tmp == "UP")
	{
		return K_UP;
	}
	if(tmp == "DOWN")
	{
		return K_DOWN;
	}
	if(tmp == "DOUBLE")
	{
		return K_DOUBLE;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are Up, Down, Double, 1, -1, 2");

	return ARM_DEFAULT_ERR;
}

long ARM_ConvInOut (const CCString& aInOut, ARM_result& result)
{
	CCString tmp = aInOut;

    tmp.toUpper ();



	if ( tmp == "IN" )
	{
	   return(K_IN);
	}

	if ( tmp == "OUT" )
	{
	   return(K_OUT);
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are In, 1, Out, -1");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvStdBoost(const CCString& aInOut, ARM_result& result)
{
	CCString tmp = aInOut;

	tmp.toUpper ();



	if ( tmp == "GLOB" )
	{
	   return(K_IN);
	}
    else if ( tmp == "BOOST" )
	{
	   return(K_OUT);
	}
    else if ( tmp == "STD" )
	{
	   return(K_STD);
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are In, 1, Out, -1");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvVolType (const CCString& aVolType, ARM_result& result)
{
	CCString tmp = aVolType;
	tmp.toUpper ();

	if ( tmp == "ATM" )
	{
	   return(K_ATMF_VOL);
	}

	if ( tmp[0] == 'A' )
	{
	   return(K_ATMF_VOL);
	}

	if ( tmp == "SMILE" )
	{
	   return(K_SMILE_VOL);
	}

	if ( tmp[0] == 'S' )
	{
	   return(K_SMILE_VOL);
	}

    if ( tmp == "FXSPLI" )
	{
	   return(K_FX_VOL_SP_INTERP);
	}

    if ( tmp[0] == 'F' )
	{
	   return(K_FX_VOL_SP_INTERP);
	}

	result.setRetCode(ARM_DEFAULT_ERR);
	result.setMsg("Invalid Parameter - Valid are ATM, Smile, FXSPLI, A, S, F");

	return(ARM_DEFAULT_ERR);
}


long ARM_ConvYesOrNo (const CCString& YesOrNo, ARM_result& result)
{
	CCString tmp = YesOrNo;
	tmp.toUpper ();

	if ( (tmp == "YES") || (tmp == "1") || (tmp[0] == 'Y') )
	{
		return K_YES;
	}

	if ( (tmp == "NO") || (tmp == "0") || (tmp[0] == 'N') )
	{
		return K_NO;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are Yes, No, Y, N, 0 or 1");

	return ARM_DEFAULT_ERR;
}
bool ARM_ConvYesOrNo(const std::string& s)
{
	if (!s.empty()) 
	{
		if ( s=="Y" || s=="1" || s=="YES" ) return true; 
		if ( s=="N" || s=="0" || s=="NO" ) return false; 
	}
	ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_ConvYesOrNo("<<s<<"): Valid are Yes, No, Y, N, 0 or 1"); 
}

long ARM_ConvCorrPfType (const CCString& corrPfType, ARM_result& result)
{
	CCString tmp = corrPfType;
	tmp.toUpper ();

	if ( (tmp == "0") || (tmp == "SIGMA") )
	{
		return 0;
	}

	if ( (tmp == "1") || (tmp == "MR") )
	{
		return 1;
	}

	if ( (tmp == "2") || (tmp == "BETA") )
	{
		return 2;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are SIGMA, 0, MR, 1, BETA or 2");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvSmiledModelFlag (const CCString& SmiledModelType, ARM_result& result)
{
	CCString tmp = SmiledModelType;
	tmp.toUpper ();

	if(tmp == "LD")
	{
		return K_LD;
	}
	if(tmp == "SABR_G")
	{
		return K_SABR_GEO;
	}
	if(tmp == "SABR_A")
	{
		return K_SABR_ARITH;
	}
	if(tmp == "SABR_IMPLNVOL")
	{
		return K_SABR_IMPLNVOL;
	}
	if(tmp == "SABR_WEIGHT")
	{
		return K_SABR_WEIGHT;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are LD, SABR_G, SABR_A, SABR_IMPLNVOL, SABR_WEIGHT");

	return ARM_DEFAULT_ERR;
}




long ARM_GetDimType( const CCString& aRule, ARM_result& result)
{
	CCString tmp = aRule;
	tmp.toUpper ();

	if ( tmp == "TENOR" )
		return K_TENOR;
	if ( tmp == "TIMETOSTART" )
		return K_TIMETOSTART;
	if ( tmp == "TIMETOEXPIRY" )
		return K_TIMETOEXPIRY;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Dim Type - Valid are Tenor and TimeToStart TimeToExpiry case insentitive");

	return ARM_DEFAULT_ERR;
}


extern long ARM_NotionalType(const CCString& aRule, ARM_result& result)
{
	CCString tmp = aRule;
	tmp.toUpper ();


    if ( tmp == "NXNONE" ||tmp == "NONE" )
	{
       return(K_NX_NONE);
	}
    
    if ( tmp == "NXEND" )
	{
        return(K_NX_END);
	}

    if ( tmp == "NXINFEND" )
	{
        return(K_NX_INFEND);
	}

    if ( tmp == "NXASINFEND" )
	{
        return(K_NX_ASINFEND);
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Notional Exchange Rule - Valid are  NXNONE, NXEND, NXASINFEND or NXINFEND");

	return ARM_DEFAULT_ERR;
}


long IsCall (const CCString& aString)
{
	CCString tmp = aString;
	tmp.toUpper ();

	if(tmp[0] == 'C')
	{
		return (1);
	}

	return (0);
}



long IsPut (const CCString& aString)
{
	CCString tmp = aString;
	tmp.toUpper ();

	if(tmp[0] == 'P')
	{
		return (1);
	}

	return (0);
}



long IsCap (const CCString& aString)
{
	CCString tmp = aString;
	tmp.toUpper ();

	if(tmp[0] == 'C')
	{
		return (1);
	}

	return (0);
}



long IsFloor (const CCString& aString)
{
	CCString tmp = aString;
	tmp.toUpper ();

	if(tmp[0] == 'F')
	{
		return (1);
	}

	return (0);
}



long IsReceive (const CCString& aString)
{
	CCString tmp = aString;
	tmp.toUpper ();

	if(tmp[0] == 'R')
	{
		return (1);
	}

	return (0);
}



long IsPay (const CCString& aString)
{
	CCString tmp = aString;
	tmp.toUpper ();

	if(tmp[0] == 'P')
	{
		return (1);
	}

	return (0);
}



long IsBermudan (const CCString& aString)
{
	CCString tmp = aString;
	tmp.toUpper ();

	if(tmp[0] == 'B')
	{
		return (1);
	}

	return (0);
}



long StrikeCode (const CCString& aString, ARM_result& result)
{
	CCString tmp = aString;
	tmp.toUpper ();

    if((!tmp) || (tmp[0] == 'P') || (tmp[0] == '1') || (tmp[0] == 'S'))
	{
		return (1);
	}
	if((tmp[0] == 'Y') || (tmp[0] == 'R') || (tmp[0] == '0'))
	{
		return (0);
	}
    if((tmp[0] == 'M') || (tmp[0] == 'T') || (tmp[0] == '3'))
	{
		return (3);
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are ATM, Smile, A, S, 0, 1");

	return ARM_DEFAULT_ERR;
}


long ARM_StrikeCode (const CCString& aString, ARM_result& result)
{
	return StrikeCode( aString, result );
}



long ARM_ConvIrIndIdToFreqId (long IrIndId, ARM_result& result)
{
	if(IrIndId == K_PIBOR1M)
	{
		return K_MONTHLY;
	}
	if(IrIndId == K_PIBOR3M)
	{
		return K_QUARTERLY;
	}
	if(IrIndId == K_PIBOR6M)
	{
		return K_SEMIANNUAL;
	}
	if(IrIndId == K_PIBOR1Y)
	{
		return K_ANNUAL;
	}
	if(IrIndId == K_LIBOR1M)
	{
		return K_MONTHLY;
	}
	if(IrIndId == K_LIBOR3M)
	{
		return K_QUARTERLY;
	}
	if(IrIndId == K_LIBOR6M)
	{
		return K_SEMIANNUAL;
	}
	if(IrIndId == K_LIBOR1Y)
	{
		return K_ANNUAL;
	}
	if(IrIndId == K_EURIBOR1M)
	{
		return K_MONTHLY;
	}
	if(IrIndId == K_EURIBOR3M)
	{
		return K_QUARTERLY;
	}
	if(IrIndId == K_EURIBOR6M)
	{
		return K_SEMIANNUAL;
	}
	if(IrIndId == K_EURIBOR1Y)
	{
		return K_ANNUAL;
	}
	
	return ARM_DEFAULT_ERR;
}





long ARM_ConvAutoMode (const CCString& aAutoMode, ARM_result& result)
{
	CCString tmp = aAutoMode;
	tmp.toUpper ();

	if(tmp == "SWOPT")
	{
		return INTERFACE_AUTO_MODE::AUTO_MODE_IDX;
	}
	if(tmp == "IRG")
	{
		return INTERFACE_AUTO_MODE::AUTO_MODE_IRG;
	}
	if(tmp == "SWOPT_IRG")
	{
		return INTERFACE_AUTO_MODE::AUTO_MODE_SWOPT_IRG;
	}
	if(tmp == "IRG_SWOPT")
	{
		return INTERFACE_AUTO_MODE::AUTO_MODE_IRG_SWOPT;
	}
	
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are SWOPT, IRG, SWOPT_IRG or IRG_SWOPT");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvAutoMode2 (const CCString& aAutoMode, ARM_result& result)
{
	CCString tmp = aAutoMode;
	tmp.toUpper ();

	if(tmp == "SWOPT")
	{
		return INTERFACE_AUTO_MODE::AUTO_MODE_IDX;
	}
	if(tmp == "IRG")
	{
		return INTERFACE_AUTO_MODE::AUTO_MODE_IRG;
	}
	if(tmp == "SWOPT_IRG")
	{
		return INTERFACE_AUTO_MODE::AUTO_MODE_SWOPT_IRG;
	}
	if(tmp == "IRG_SWOPT")
	{
		return INTERFACE_AUTO_MODE::AUTO_MODE_IRG_SWOPT;
	}
	if(tmp == "DUAL")
	{
		return INTERFACE_AUTO_MODE::AUTO_MODE_DUAL;
	}
	if(tmp == "MEAN")
	{
		return INTERFACE_AUTO_MODE::AUTO_MODE_AUTOMEAN;
	}
	if(tmp == "DIAG")
	{
		return INTERFACE_AUTO_MODE::AUTO_MODE_DIAG;
	}
	if(tmp == "CAPDIAG")
	{
		return INTERFACE_AUTO_MODE::AUTO_MODE_CAPDIAG;
	}
	if(tmp == "SWOPT_IRG_DUAL_CAPDIAG")
	{
		return (INTERFACE_AUTO_MODE::AUTO_MODE_SWOPT_IRG
			    | INTERFACE_AUTO_MODE::AUTO_MODE_DUAL 
				| INTERFACE_AUTO_MODE::AUTO_MODE_CAPDIAG);
	}
	if(tmp == "SWOPT_IRG_DUAL")
	{
		return (INTERFACE_AUTO_MODE::AUTO_MODE_SWOPT_IRG
			    | INTERFACE_AUTO_MODE::AUTO_MODE_DUAL);
	}
	if(tmp == "SWOPT_IRG_DUAL_DDIAG")
	{
		return (INTERFACE_AUTO_MODE::AUTO_MODE_SWOPT_IRG
			    | INTERFACE_AUTO_MODE::AUTO_MODE_DUAL
				| INTERFACE_AUTO_MODE::AUTO_MODE_DIAG
				| INTERFACE_AUTO_MODE::AUTO_MODE_CAPDIAG);
	}
	if(tmp == "SWOPT_IRG_AUTOMEAN")
	{
		return (INTERFACE_AUTO_MODE::AUTO_MODE_SWOPT_IRG
			    | INTERFACE_AUTO_MODE::AUTO_MODE_AUTOMEAN);
	}
	if(tmp == "SWOPT_IRG_AUTOMEAN_BASIS")
	{
		return (INTERFACE_AUTO_MODE::AUTO_MODE_SWOPT_IRG
			    | INTERFACE_AUTO_MODE::AUTO_MODE_AUTOMEAN
				| INTERFACE_AUTO_MODE::AUTO_MODE_BASISDISC);
	}
	if(tmp == "SWOPT_IRG_AUTOMEAN_DIAG")
	{
		return (INTERFACE_AUTO_MODE::AUTO_MODE_SWOPT_IRG
			    | INTERFACE_AUTO_MODE::AUTO_MODE_AUTOMEAN
				| INTERFACE_AUTO_MODE::AUTO_MODE_DIAG);
	}
	if(tmp == "SWOPT_IRG_AUTOMEAN_DIAG_BASIS")
	{
		return (INTERFACE_AUTO_MODE::AUTO_MODE_SWOPT_IRG
			    | INTERFACE_AUTO_MODE::AUTO_MODE_AUTOMEAN
				| INTERFACE_AUTO_MODE::AUTO_MODE_DIAG
				| INTERFACE_AUTO_MODE::AUTO_MODE_BASISDISC);
	}
	if(tmp == "SWOPT_BASIS")
	{
		return (INTERFACE_AUTO_MODE::AUTO_MODE_IDX
				| INTERFACE_AUTO_MODE::AUTO_MODE_BASISDISC);
	}
	if(tmp == "SWOPT_BASIS_DIAG")
	{
		return (INTERFACE_AUTO_MODE::AUTO_MODE_IDX
				| INTERFACE_AUTO_MODE::AUTO_MODE_BASISDISC
				| INTERFACE_AUTO_MODE::AUTO_MODE_DIAG);
	}
	if(tmp == "SWOPT_IRG_DUAL_CAPDIAG_BASIS")
	{
		return (INTERFACE_AUTO_MODE::AUTO_MODE_SWOPT_IRG
				| INTERFACE_AUTO_MODE::AUTO_MODE_DUAL
				| INTERFACE_AUTO_MODE::AUTO_MODE_CAPDIAG
				| INTERFACE_AUTO_MODE::AUTO_MODE_BASISDISC);
	}
	if(tmp == "SWOPT_IRG_DUAL_BASIS")
	{
		return (INTERFACE_AUTO_MODE::AUTO_MODE_SWOPT_IRG
				| INTERFACE_AUTO_MODE::AUTO_MODE_DUAL
				| INTERFACE_AUTO_MODE::AUTO_MODE_BASISDISC);
	}
	if(tmp == "SWOPT_IRG_DUAL_DDIAG_BASIS")
	{
		return (INTERFACE_AUTO_MODE::AUTO_MODE_SWOPT_IRG
				| INTERFACE_AUTO_MODE::AUTO_MODE_DUAL
				| INTERFACE_AUTO_MODE::AUTO_MODE_DIAG
				| INTERFACE_AUTO_MODE::AUTO_MODE_CAPDIAG
				| INTERFACE_AUTO_MODE::AUTO_MODE_BASISDISC);
	}
	
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are SWOPT, IRG, SWOPT_IRG, IRG_SWOPT, DUAL, MEAN, DIAG, SWOPT_IRG_DUAL_DIAG, SWOPT_IRG_AUTOMEAN_DIAG or SWOPT_IRG_DUAL_DIAG_BASIS");

	return ARM_DEFAULT_ERR;
}

long ARM_ConvFineMode (const CCString& aFineMode, ARM_result& result)
{
	CCString tmp = aFineMode;
	tmp.toUpper ();


    if (tmp == "L_NONE")
        return K_L_NONE;
    else if (tmp == "L_DEFAULT")
        return K_L_DEFAULT;
    else if (tmp == "L_QUICK")
        return K_L_QUICK;
    else if (tmp == "L_2")
        return K_L_2;
    else if (tmp == "L_4")
        return K_L_4;
    else if (tmp == "L_8")
        return K_L_8;
    else if (tmp == "L_12")
        return K_L_12;
    else if (tmp == "L_16")
        return K_L_16;
    else if (tmp == "L_24")
        return K_L_24;
    else if (tmp == "L_7DAY")
        return K_L_7DAY;
    else if (tmp == "L_5DAY")
        return K_L_5DAY;
    else if (tmp == "L_3DAY")
        return K_L_3DAY;
    else if (tmp == "L_2DAY")
        return K_L_2DAY;
    else if (tmp == "L_1DAY")
        return K_L_1DAY; 
    else if (tmp == "L_A")
           return K_L_A; 
    else if (tmp == "L_B")
           return K_L_B; 
    else if (tmp == "L_C")
           return K_L_C; 
    else 
    {
        result.setRetCode (ARM_DEFAULT_ERR);
        result.setMsg ("Invalid Parameter - Valid are L_QUICK, L_DEFAULT, L_2, L_4, L_8, L_12, L_16, L_24");

        return ARM_DEFAULT_ERR;
    }
}


long ARM_ConvMcMethod (const CCString& aMcMethod, ARM_result& result)
{
	CCString tmp = aMcMethod;
	tmp.toUpper ();


    if ((tmp == "SCRAMBLE") || (tmp == "SC"))
        return K_MC_SCRAMBLE;
    else if ((tmp == "FAURE") || (tmp == "F"))
        return K_MC_FAURE;
    else if ((tmp == "SIMPLE") || (tmp == "SI"))
        return K_MC_SIMPLE;
    else 
    {
        result.setRetCode (ARM_DEFAULT_ERR);
        result.setMsg ("Invalid Parameter - Valid are SCRAMBLE, FAURE, SIMPLE, SC, F, SI");

        return ARM_DEFAULT_ERR;
    }
}

long ARM_ConvShapeType (const CCString& aShapeType, ARM_result& result)
{
	CCString tmp = aShapeType;
	tmp.toUpper ();


    if ((tmp == "ROW") || (tmp == "R"))
        return K_ROW;
    else if ((tmp == "DIAG") || (tmp == "D"))
        return K_DIAG;
    else
    {
        result.setRetCode (ARM_DEFAULT_ERR);
        result.setMsg ("Invalid Parameter - Valid are ROW, DIAG, R, D");

        return ARM_DEFAULT_ERR;
    }
}

long ARM_ConvCalculationMethod (const CCString& aCalcMeth, ARM_result& result)
{
	CCString tmp = aCalcMeth;
	tmp.toUpper ();

    if(tmp == "L")
	{
		return K_LININTERPOL_REF;
	}
    if(tmp == "S")
	{
        return K_DISCRETE_REF;
	}
    if(tmp == "LIN")
	{
        return K_LININTERPOL_REF;
	}
    if(tmp == "STEP")
	{
        return K_DISCRETE_REF;
	}
    if(tmp == "STEPUP_RIGHT")
	{
        return K_STEPUP_RIGHT;
	}
    if(tmp == "STEPUP_LEFT")
	{
        return K_STEPUP_LEFT;
	}
    if(tmp == "0")
	{
        return K_LININTERPOL_REF;
	}
    if(tmp == "1")
	{
        return K_DISCRETE_REF;
	}
    if(tmp == "LIN_CST_EXTRA")
	{
        return K_LIN_INTER_CST_EXTRA_REF;
	}
    if(tmp == "PERF_DISC_REF")
	{
        return K_PERFECT_DISCRETE_REF;
	}	
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Interpol Method - Valid are: (LIN,0) - (STEP,1) - STEPUP_RIGHT - STEPUP_LEFT, LIN_CST_EXTRA, PERF_DISC_REF");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvSmileNoSmile (CCString* aSmileNotSmile, ARM_result& result)
{
	aSmileNotSmile->toUpper ();


    if ((*aSmileNotSmile == "S") || (*aSmileNotSmile == "SMILE") ||
		(*aSmileNotSmile == "N") || (*aSmileNotSmile == "NOSMILE") ||
		(*aSmileNotSmile == "A") || (*aSmileNotSmile == "ATMANDSMILE"))
		return 1;
    else
    {
        result.setRetCode (ARM_DEFAULT_ERR);
        result.setMsg ("Invalid Parameter - Valid are SMILE, NOSMILE, ATMANDSMILE, S, N, A");

        return ARM_DEFAULT_ERR;

    }
}

long ARM_ConvAmortMethod (const CCString& aAmortMethod, ARM_result& result)
{
	CCString tmp = aAmortMethod;
	tmp.toUpper ();

	if ( (tmp == "FIXED") || (tmp == "F") )
		return K_AMORT_FIXED;

	if ( (tmp == "FIXEDEND") || (tmp == "FE") )
		return K_AMORT_FIXEDEND;

	if ( (tmp == "PERCENTAGE") || (tmp == "P") )
		return K_AMORT_PERCENTAGE;

	if ( (tmp == "ANNUITY") || (tmp == "A") )
		return K_AMORT_ANNUITY;

	if ( (tmp == "ANNUITYREDUCED") || (tmp == "AR") )
		return K_AMORT_ANNUITY_REDUCED;

	if ( (tmp == "MORTGAGE") || (tmp == "M") )
		return K_AMORT_MORTGAGE;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are FIXED, FIXEDEND, PERCENTAGE, ANNUITY, ANNUITYREDUCED, MORTGAGE, F, FE, P, A, AR or M");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvModelType (const CCString& aModelType, ARM_result& result)
{
	CCString tmp = aModelType;
	tmp.toUpper ();

	if ( (tmp == "LOG") || (tmp == "L") )
		return K_LOG;

	if ( (tmp == "NOR") || (tmp == "N") )
		return K_NOR;

	if ( (tmp == "2LOG") || (tmp == "2") )
		return K_2LOG;

	if ( (tmp == "2SABR") || (tmp == "B") )
		return K_2SABR;

	if ( (tmp == "GCOP") || (tmp == "G") )
		return K_GCOP;

	if (tmp == "2LOG_EXTENDED")
		return K_2LOG_EXTENDED;

	if (tmp == "2SABR_CALIB_DEGEN")
		return K_2SABR_CALIB_DEGEN;

	if (tmp == "2SABR_CALIB_NODEGEN")
		return K_2SABR_CALIB_NODEGEN;

	if (tmp == "2SABR_NOCALIB_DEGEN")
		return K_2SABR_NOCALIB_DEGEN;

	if (tmp == "2SABR_NOCALIB_NODEGEN")
		return K_2SABR_NOCALIB_NODEGEN;

	if (tmp == "GCOP_CALIB_DEGEN")
		return K_GCOP_CALIB_DEGEN;

	if (tmp == "GCOP_CALIB_NODEGEN")
		return K_GCOP_CALIB_NODEGEN;

	if (tmp == "GCOP_NOCALIB_DEGEN")
		return K_GCOP_NOCALIB_DEGEN;

	if (tmp == "GCOP_NOCALIB_NODEGEN")
		return K_GCOP_NOCALIB_NODEGEN;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are 2LOG or 2, LOG or L, NOR or N, 2SABR or B, GCOP or G, 2LOG_EXTENDED, 2SABR_[NO]CALIB_[NO]DEGEN, GCOP_[NO]CALIB_[NO]DEGEN");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvVolType2 (const CCString& aVolType, ARM_result& result)
{
	CCString tmp = aVolType;
	tmp.toUpper ();

	if ( (tmp == "INPUT") || (tmp == "I") )
		return K_INPUTED;

	if ( (tmp == "COMPUTE") || (tmp == "C") )
		return K_COMPUTED;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are INPUT, I, COMPUTE or C");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvDatePowerReverseDataType(const CCString& aValuesType, ARM_result& result)
{
	CCString tmp = aValuesType;
	tmp.toUpper ();

	if ( (tmp == "CANCELLATIONDATE") || (tmp == "CD") )
		return K_CANCEL_DATES;

	if ( (tmp == "NOTIFICATIONDATE") || (tmp == "ND") )
		return K_NOTIF_DATES;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are CANCELLATIONDATE(CD), NOTIFICATIONDATE(ND)");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvDateStripDataType(const CCString& aValuesType, ARM_result& result)
{
	CCString tmp = aValuesType;
	tmp.toUpper ();

	if ( (tmp == "STARTDATE") || (tmp == "SD") )
		return K_START_DATES;

	if ( (tmp == "ENDDATE") || (tmp == "ED") )
		return K_END_DATES;

	if ( (tmp == "FWDSTARTDATE") || (tmp == "FSD") )
		return K_FWD_START_DATES;

	if ( (tmp == "FWDENDDATE") || (tmp == "FED") )
		return K_FWD_END_DATES;

	if ( (tmp == "RESETDATE") || (tmp == "RD") )
		return K_RESET_DATES;

	if ( (tmp == "PAYDATE") || (tmp == "PD") )
		return K_PAY_DATES;

	if ( (tmp == "INTERESTDAYS") || (tmp == "ID") )
		return K_INT_DAYS;

	if ( (tmp == "INTERESTTERMS") || (tmp == "IT") )
		return K_INT_TERMS;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are STARTDATE(SD), ENDDATE(ED), FWDSTARTDATE(FSD), FWDENDDATE(FED), RESETDATE(RD), PAYDATE(PD), INTERESTTERMS(IT) or INTERESTDAYS(ID)");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvTypeValues (const CCString& aValuesType, ARM_result& result)
{
	CCString tmp = aValuesType;
	tmp.toUpper ();

	if ( (tmp == "FWDRATE") || (tmp == "FR") )
		return K_FWD_RATES;

	if ( (tmp == "RAWFWDRATE") || (tmp == "RF") )
		return K_RAW_FWD_RATES;

	if ( (tmp == "INTERESTDAYS") || (tmp == "ID") )
		return K_INT_DAYS;

	if ( (tmp == "INTERESTTERMS") || (tmp == "IT") )
		return K_INT_TERMS;

	if ( (tmp == "FLOWVALUE") || (tmp == "FV") )
		return K_FLOW_VALUES;

    if ( (tmp == "FLOWVALUEPV") || (tmp == "FPV") )
		return K_FLOWS_PV_VALUES;

	if ( (tmp == "NOTIONALVALUE") || (tmp == "NV") )
		return K_NOTIONAL_VALUES;

	if ( (tmp == "AMORTVALUE") || (tmp == "AV") )
		return K_AMORT_VALUES;

    if ( (tmp == "VOLFWD") || ( tmp == "VF") )
		return K_VOL_FWD;

    if ( (tmp == "VOLCAP") || ( tmp == "VLC" ) )
		return K_VOL_CAP;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are FWDRATE(FR), RAWFWDRATE(RF), INTERESTDAYS(ID), INTERESTTERMS(IT), FLOWVALUE(FV), NOTIONALVALUE(NV), AMORTVALUE(AV), VOLFWD(VF), FLOWVALUEPV(FPV), VOLCAP(VLC)");

	return ARM_DEFAULT_ERR;
}

long ARM_ConvTypeDates (const CCString& aDatesType, ARM_result& result)
{
	CCString tmp = aDatesType;
	tmp.toUpper ();

	if ( (tmp == "STARTDATE") || (tmp == "SD") )
		return K_START_DATES;

	if ( (tmp == "ENDDATE") || (tmp == "ED") )
		return K_END_DATES;

	if ( (tmp == "RESETDATE") || (tmp == "RD") )
		return K_RESET_DATES;

	if ( (tmp == "PAYDATE") || (tmp == "PD") )
		return K_PAY_DATES;

	if ( (tmp == "FWDSTARTDATE") || (tmp == "FSD") )
		return K_FWD_START_DATES;

	if ( (tmp == "FWDENDDATE") || (tmp == "FED") )
		return K_FWD_END_DATES;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are STARTDATE(SD), ENDDATE(ED), RESETDATE(RD), PAYDATE(PD), FWDSTARTDATE(FSD), FWDENDDATE(FED)");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvSigRhoNu (const CCString& param, ARM_result& result)
{
	CCString tmp = param;
	tmp.toUpper ();

	if ( (tmp == "SIGMA") || (tmp == "S") )
		return ARM_SIGMA;

	if ( (tmp == "RHO") || (tmp == "R") )
		return ARM_RHO;

	if ( (tmp == "NU") || (tmp == "N") )
		return ARM_NU;

	if ( (tmp == "BETA") || (tmp == "B") )
		return ARM_BETA;
	

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are SIGMA, RHO, BETA or NU");

	return ARM_DEFAULT_ERR;
}

long ARM_ConvPricerType (const CCString& param, ARM_result& result)
{
	CCString tmp = param;
	tmp.toUpper ();

    if (tmp == "FNMC")
		return 1;

	if (tmp == "RNMC")
		return 2;

	if (tmp == "LSMC")
		return 7;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are RNMC or LSMC");

	return ARM_DEFAULT_ERR;
}

long ARM_ConvBootstrap (const CCString& param, ARM_result& result)
{
	CCString tmp = param;
	tmp.toUpper ();

	if (tmp == "B")
		return 0;

	if (tmp == "G")
		return 1;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are G or B");

	return ARM_DEFAULT_ERR;
}


extern CCString ARM_ConvTypeDeal(const CCString& typeDeal, ARM_result& result)
{
	CCString tmp = typeDeal;
	tmp.toUpper ();


    if (tmp == "SWAP")
	{
		return(LOCAL_SWAP_CLASS);
	}
    if (tmp == "SWOPT")
	{
		return(LOCAL_SWAPTION_CLASS);
	}
    if ((tmp == "FXOPT") || (tmp == "FXOPTSTRIP"))
	{
		return(LOCAL_OPTION_CLASS);
	}
    if ((tmp == "PRCS") || (tmp == "PRCS2"))
	{
		return(LOCAL_POWER_REVERSE_CLASS);
	}
    if (tmp == "CDS")
	{
		return(LOCAL_CDS_CLASS);
	}
    if (tmp == "FTD")
	{
		return(LOCAL_FTD_CLASS);
	}
    if (tmp == "NTD")
	{
		return(LOCAL_NTHTD_CLASS);
	}
    if (tmp == "CDO")
	{
		return(LOCAL_MEZ_CLASS);
	}
    if (tmp == "CDO2")
	{
		return(LOCAL_MEZ_CLASS);
	}
	if (tmp == "ACCRUALOPTION")
	{
		return(LOCAL_FLEX_SWAPTION_CLASS);
	}
    if (tmp == "CRF")
	{
		return(LOCAL_GC_CRF_CLASS);
	}
	if (tmp == "CRA")
	{
		return (LOCAL_OPTION_PORTFOLIO_CLASS);
	}
	if (tmp == "IRG")
	{
		return (LOCAL_CAPFLOOR_CLASS);
	}
	if (tmp == "SPDOPT") 
	{
		return (LOCAL_SPREAD_OPTION_CLASS);
	}
	if ( (tmp == "EXOTIC") || (tmp == "MEMORYSO") || (tmp == "MEMORYCAP") )
	{
		return (LOCAL_PF_CLASS);
	}
	if (tmp == "MATURITYCAP") 
	{
		return (LOCAL_GC_MATURITYCAP_CLASS);
	}
 	if (tmp == "RFTARN") 
 	{
 		return (LOCAL_GC_TARN_CLASS);
 	}
 	if (tmp == "BERM") 
 	{
 		return (LOCAL_GC_BERMUDASWAPTION_CLASS);
 	}
 	if (tmp == "ALMCAPTION") 
 	{
 		return (LOCAL_GC_CAPTION_CLASS);
 	}
	if (tmp == "CDSOPT") 
	{
		return (LOCAL_OPTION_CLASS);
	}
	if( (tmp == "EXOTIC.RA") || (tmp == "SWAP.RA") )
	{
		return (LOCAL_PF_CLASS);
	}
	if (tmp == "GLOBALCAP")
	{
		return (LOCAL_GLOBALCAP_CLASS);
	}
	if (tmp == "FXSTRIP")
	{
		return (LOCAL_FXSTRIP_CLASS);
	}
	if (tmp == "CFXSTRIP")
	{
		return (LOCAL_GC_FXVANILLA_CLASS);
	}
	if (tmp == "RNG_DOUBLE")
	{
		return (LOCAL_PF_CLASS);
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Deal type - Valid are SWAP, SWOPT, FXOPT, PRCS, ACCRUALOPTION, CRF, CDS, FTD, NTD, MEZ, CRA, IRG, SPDOPT, EXOTIC, CDSOPT, MATURITYCAP, RFTARN, ALMCAPTION, CDO2, FXOPTSTRIP, RNG_DOUBLE");

	return ARM_DEFAULT_ERR;
}


long ARM_Convhedge(const CCString& typeHedge, ARM_result& result)
{
	CCString tmp = typeHedge;
	tmp.toUpper ();

	if (tmp == "IRDELTA")
		return 1;

	if (tmp == "IRVEGA")
		return 2;

	if (tmp == "FXDELTA")
		return 3;

	if (tmp == "FXVEGA")
		return 4;

	if (tmp == "MEANREV")
		return 5;

	if (tmp == "CORREL")
		return 6;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are IRDELTA, IRVEGA, FXDELTA, FXVEGA, MEANREV or CORREL");

	return ARM_DEFAULT_ERR;
}



long LocalExtractCurveFromFileMO (const CCString& curveFileName, std::vector<CCString>& matu, std::vector<double>& rate,long adjOrNotId)
{
	FILE* curveFile = fopen ((const char*)curveFileName, "r");

	if(curveFile == NULL)
	{
		return ARM_KO;
	}
	
	char sEch[MATURITY_SIZE];
	char buffer[MATURITY_SIZE];
	double val;
	double adj;

	matu.clear ();
	rate.clear ();

	int status = 0;

	while(status != EOF)
	{
		status = fscanf(curveFile, "%s",buffer);
		status = fscanf(curveFile, "%s",sEch);
		status = fscanf(curveFile, "%lf",&val);
		status = fscanf(curveFile, "%lf/n",&adj);

		matu.push_back (CCString (sEch));
		if (adjOrNotId == K_YES)
			rate.push_back (val+adj/100.);
		else
			rate.push_back (val);
	}

	fclose (curveFile);

	return ARM_OK;
}


long LocalExtractCurveFromFileBS (const CCString& curveFileName, std::vector<CCString>& matu, std::vector<double>& rate)
{
	FILE* curveFile = fopen ((const char*)curveFileName, "r");

	if(curveFile == NULL)
	{
		return ARM_KO;
	}
	
	char sEch[MATURITY_SIZE];
	char buffer[MATURITY_SIZE];
	double val;

	matu.clear ();
	rate.clear ();

	int status = 0;

	while(status != EOF)
	{
		status = fscanf(curveFile, "%s",buffer);
		status = fscanf(curveFile, "%s",sEch);
		status = fscanf(curveFile, "%s",buffer);
		status = fscanf(curveFile, "%lf/n",&val);

		matu.push_back (CCString (sEch));
		rate.push_back (val);
	}

	fclose (curveFile);

	return ARM_OK;
}


long LocalExtractVolFromFile (const CCString& FileName, std::vector<CCString>& yearTerm, std::vector<CCString>& tenor, std::vector<double>& vol)
{
	FILE* Fp = fopen ((const char*)FileName, "r");

	if(Fp == NULL)
	{
		return ARM_KO;
	}

	int indiceI;
	int indiceJ;
	int i, j;
	int compteur(0);

	char sEch[50];
	char buffer[50];
	double val;

	int rc = fscanf(Fp, "%s",buffer);

	while (rc != EOF)
	{
		rc = fscanf(Fp, "%s",sEch);

		indiceI = -1;
		for (i=0;i<tenor.size();i++)
		{
			if (strcmp((const char*)tenor[i],sEch) == 0)
			{
				indiceI = i;
				i = tenor.size();
			}
		}
		if (indiceI == -1)
		{
			indiceI = tenor.size();
			tenor.push_back(sEch);
		}

		rc = fscanf(Fp, "%s",sEch);
		indiceJ = -1;
		for (j=0;j<yearTerm.size();j++)
		{
			if (strcmp((const char*)yearTerm[j],sEch) == 0)
			{
				indiceJ = j;
				j = yearTerm.size();
			}
		}
		if (indiceJ == -1)
		{
			indiceJ = yearTerm.size();
			yearTerm.push_back(sEch);
		}

		rc = fscanf(Fp, "%lf",&val);

		vol.push_back(val);

		compteur++;
		rc = fscanf(Fp, "%s",buffer);
	}

	fclose(Fp);

	return ARM_OK;
}



long LocalExtractSmileFromFile (const CCString& FileName, std::vector<CCString>& yearTerm, std::vector<double>& strike, std::vector<double>& vol)
{
	FILE* Fp = fopen ((const char*)FileName, "r");

	if(Fp == NULL)
	{
		return ARM_KO;
	}
	int indiceI;
	int indiceJ;
	int i, j;
	int compteur(0);

	char sEch[50];
	char buffer[50];
	double val;

	int rc = fscanf(Fp, "%s",buffer);
	while (rc != EOF)
	{
		rc = fscanf(Fp, "%s",sEch);

		indiceI = -1;
		for (i=0;i<yearTerm.size();i++)
		{
			if (strcmp((const char*)yearTerm[i],sEch) == 0)
			{
				indiceI = i;
				i = yearTerm.size();
			}
		}
		if (indiceI == -1)
		{
			indiceI = yearTerm.size();
			yearTerm.push_back(sEch);
		}

		rc = fscanf(Fp, "%s",sEch);
		val = atof(sEch);
		indiceJ = -1;
		for (j=0;j<strike.size();j++)
		{
			if (strike[j] == val)
			{
				indiceJ = j;
				j = strike.size();
			}
		}
		if (indiceJ == -1)
		{
			indiceJ = strike.size();
			strike.push_back(val);
		}

		rc = fscanf(Fp, "%lf",&val);

		vol.push_back(val);

		compteur++;
		rc = fscanf(Fp, "%s",buffer);
	}

	fclose(Fp);

	return ARM_OK;
}



long ExtractVectorDoubleFromFile (const CCString& FileName, std::vector<double>& vect, long& vecSize)
{
    CCString clientViewFileName = CCString(VIEW_FILE_CLIENT_LOCATION)+CCString(VIEW_FILE_PREFIX)+FileName;

	FILE* Fp = fopen ((const char*)clientViewFileName, "r");

	if(Fp == NULL)
	{
		return ARM_KO;
	}

	double val;
	int tmpSize;

	vect.clear ();

	int status = 0;

	status = fscanf(Fp, "%d\n",&tmpSize);

	while(status != EOF)
	{
		status = fscanf(Fp, "%lf\n",&val);

		/// we should push value only if we are not at the 
		/// end of the file!
		if(status != EOF )
			vect.push_back (val);
	}

	vecSize = tmpSize;

	fclose (Fp);

	DeleteFile (clientViewFileName);

	return ARM_OK;
}


long ExtractVectorDateFromFile (const CCString& FileName, std::vector<CCString>& vect, long& vecSize)
{
    CCString clientViewFileName = CCString(VIEW_FILE_CLIENT_LOCATION)+CCString(VIEW_FILE_PREFIX)+FileName;

	FILE* Fp = fopen ((const char*)clientViewFileName, "r");

	if(Fp == NULL)
	{
		return ARM_KO;
	}

	char myDate[MATURITY_SIZE];
	int tmpSize;

	vect.clear ();

	int status = 0;


	status = fscanf(Fp, "%d\n",&tmpSize);

	while(status != EOF)
	{
		status = fscanf(Fp, "%s\n",myDate);

		/// we should push value only if we are not at the 
		/// end of the file!
		if(status!= EOF )
			vect.push_back (myDate);
	}

	fclose (Fp);

	vecSize = tmpSize;

	DeleteFile (clientViewFileName);

	return ARM_OK;
}


long LocalExtractInfoFromFilePRCS(const CCString& FileName, std::vector<CCString>& listeString, std::vector<double>* listeDouble, int& nbCol)
{
    CCString clientViewFileName = CCString(VIEW_FILE_CLIENT_LOCATION)+CCString(VIEW_FILE_PREFIX)+FileName;

	FILE* Fp = fopen ((const char*)clientViewFileName, "r");

	if(Fp == NULL)
	{
		return ARM_KO;
	}

	int nbLine;

	int status = 0;

	char buffer[50];

	status = fscanf(Fp, "%s",buffer);
	status = fscanf(Fp, "%s",buffer);
	status = fscanf(Fp, "%s",buffer);
	status = fscanf(Fp, "%s\n",buffer);
	status = fscanf(Fp, "NbLin=%d\n",&nbLine);
	status = fscanf(Fp, "NbCol=%d\n",&nbCol);

	double dVal;

	while(status != EOF)
	{
		status = fscanf(Fp, "%s",buffer);
		listeString.push_back(CCString (buffer));

		status = fscanf(Fp, "%lf",&dVal);
		listeDouble[0].push_back(dVal);

		for (int i=2;i<=6;i++)
		{
			status = fscanf(Fp, "%s",buffer);
			status = fscanf(Fp, "%s",buffer);
			listeDouble[i-1].push_back(Local_ARMDATE2XLDATE(CCString (buffer)));
		}

		for (i=7;i<=nbCol-2;i++)
		{
			status = fscanf(Fp, "%lf",&dVal);
			listeDouble[i-1].push_back(dVal);
		}

		status = fscanf(Fp, "%lf/n",&dVal);
		listeDouble[nbCol-2].push_back(dVal);
	}

	fclose (Fp);

	DeleteFile (clientViewFileName);

	return ARM_OK;
}

// Fonctions de conversion pour le crédit


extern long ARM_ConvStringToBool (const CCString& aBoolean, ARM_result& result)
{
	CCString tmp = aBoolean;
	tmp.toUpper ();

	if (tmp == "TRUE")
		return (long)true;

	if (tmp == "FALSE")
		return (long)false;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are TRUE, FALSE");

	return ARM_DEFAULT_ERR;
}

extern long ARM_ConvCBOptionType (const CCString& aOptionType, ARM_result& result)
{
	CCString tmp = aOptionType;
	tmp.toUpper ();

	if (tmp == "CONVERSION")
		return 1;

	if (tmp == "ISSUERCALL")
		return 2;

	if (tmp == "BEARERPUT")
		return 3;

	if (tmp == "NONE")
		return 0;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Option Type : Invalid Parameter - Valid are CONVERSION, ISSUERCALL, BEARERPUT, NONE");

	return ARM_DEFAULT_ERR;
}

extern long ARM_ConvCBOptionStrikeType (const CCString& aStrikeType, ARM_result& result)
{
	CCString tmp = aStrikeType;
	tmp.toUpper ();

	if (tmp == "NBSHARES")
		return 2;

	if (tmp == "BONDPRICE")
		return 3;

	if (tmp == "BONDYIELD")
		return 4;

	if (tmp == "NONE")
		return 0;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Strike Type : Invalid Parameter - Valid are STOCKPRICE, NBSHARES, BONDPRICE, BONDYIELD, NONE");

	return ARM_DEFAULT_ERR;
}

extern long ARM_ConvCBOptionAccruedOnEx (const CCString& aAccrued, ARM_result& result)
{
	CCString tmp = aAccrued;
	tmp.toUpper ();

	if (tmp == "ACC")
		return 1;

	if (tmp == "CURR")
		return 2;

	if (tmp == "NONE")
		return 0;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Accrued On Exercise : Invalid Parameter - Valid are ACC, CURR, NONE");

	return ARM_DEFAULT_ERR;
}

extern long ARM_ConvCBOptionBarrierType (const CCString& aBarrierType, ARM_result& result)
{
	CCString tmp = aBarrierType;
	tmp.toUpper ();

	if (tmp == "UPANDIN")
		return 1;

	if (tmp == "DOWNANDIN")
		return 2;

	if (tmp == "NONE")
		return 0;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Barrier Type : Invalid Parameter - Valid are UPANDIN, DOWNANDIN, NONE");

	return ARM_DEFAULT_ERR;
}


extern long ARM_ConvCBOptionBarrierStrikeType (const CCString& aStrikeType, ARM_result& result)
{
	CCString tmp = aStrikeType;
	tmp.toUpper ();

	if (tmp == "STOCKPRICE")
		return 1;

	if (tmp == "CONVOVERACCVAL")
		return 5;

	if (tmp == "NONE")
		return 0;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Barrier strike type : Invalid Parameter - Valid are STOCKPRICE, CONVOVERACCVAL, NONE");

	return ARM_DEFAULT_ERR;
}
 

extern long ARM_ConvPaysOnDef (const CCString& aRorP, ARM_result& result)
{
	CCString tmp = aRorP;
	tmp.toUpper ();

    if(tmp == "FULL")
	{
        return 0;
	}
    if(tmp == "NOS")
	{
		return 1;
	}
    if(tmp == "REC")
	{
		return 2;
	}
    if(tmp == "LOSS")
	{
		return 3;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Pays on default - Valid are FULL, NOS, REC, LOSS");

	return ARM_DEFAULT_ERR;
}

extern long ARM_CheckDate (const double& date, ARM_result& result)
{

    if ((date >= 25569 /* 01/01/1970 */) && (date <= 73051 /* 01/01/2100 */))
	{
        return 1;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Date - must be between 01/01/1970 & 01/01/2100 ");

	return ARM_DEFAULT_ERR;
}
 
extern long ICM_ConvCptType (const CCString& option, ARM_result& result)
{
	CCString tmp = option;
	tmp.toUpper ();

      if(tmp == "NPV"){return qCMPPRICE;}
    if(tmp == "FEELEG")	{return qCMPFEELEGPV;}
    if(tmp == "DEFLEG"){return qCMPDEFLEGPV;}
    if(tmp == "ACCRUED_PV"){return qCMPACCRUED ;} // qCMPACCRUEDPV
    if(tmp == "ACCRUED"){return qCMPACCRUED;}
	if (tmp=="PREMIUM") return qCMPPREMIUM ;
	if (tmp=="FWDSPREAD") return qCMPFWDSPREAD; 
	if (tmp=="DURATION") return qCMPDURATION;
	if (tmp=="CORREL_UP") return qCMPCORRELUP;		// for tranches
	if (tmp=="CORREL_DOWN") return qCMPCORRELDOWN;
	if (tmp=="AVGCORRDEF") return qCMPAVGCORRDEF;	// for CDO2
	if (tmp=="EL") return qCMPEL;	// for CDO
	if (tmp=="FEES") return qCMPFEES;	// for CDO
	if (tmp=="SPREAD") return qCMPSPREAD;
	if (tmp=="FLATCORR") return qCMPFLATCORR ;
	if(tmp == "NONE"){return 999;}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid - Valid are: NPV,SPREAD, FEES, FEELEG, DEFLEG,ACCRUED,ACCRUED_PV,PREMIUM,FWDSPREAD,DURATION,CORREL_DOWN,CORREL_UP, AVGCORRDEF, EL");

	return ARM_DEFAULT_ERR;
}

extern long ICM_ConvGreekType (const CCString& option, ARM_result& result)
{
	CCString tmp = option;
	tmp.toUpper ();

    if(tmp == "DELTA"){	return ICM_GREEK_DELTA_TYPE;}
    if(tmp == "GAMMA"){	return ICM_GREEK_GAMMA_TYPE;}
    if(tmp == "VEGA"){ return ICM_GREEK_VEGA_TYPE;}
    if(tmp == "THETA"){	return ICM_GREEK_THETA_TYPE;}
    if(tmp == "RHO"){return ICM_GREEK_RHO_TYPE;}
    if(tmp == "SLN_SIGMA"){return ICM_SLN_SIGMA;}
    if(tmp == "SLN_M"){return ICM_SLN_M;}
    if(tmp == "SLN_DELTA"){return ICM_SLN_DELTA;}
    if(tmp == "NONE"){return 999;}
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Greek Type - Valid are: DELTA, GAMMA,VEGA,THETA, RHO");

	return ARM_DEFAULT_ERR;
}


extern long ARM_ConvTypeGenDates (const CCString& aDatesType, ARM_result& result)
{
	CCString tmp = aDatesType;
	tmp.toUpper ();

	if  (tmp == "ACC_START_DATE") 
		return 0;

	if  (tmp == "ACC_END_DATE") 
		return 1;

	if  (tmp == "PAY_DATE") 
		return 2;

	if  (tmp == "OBS_START_DATE") 
		return 3;

	if  (tmp == "OBS_END_DATE") 
		return 4;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are ACC_START_DATE, ACC_END_DATE, PAY_DATE, OBS_START_DATE, OBS_END_DATE");

	return ARM_DEFAULT_ERR;
}

qCDS_ADJ ARM_ConvAdjCalCDS (const CCString& aDatesType, ARM_result& result)
{
	CCString tmp = aDatesType;
	tmp.toUpper ();
	std::string liste;
	bool ok; 
	qCDS_ADJ ret; 
	ICM_EnumsCnv::cnv(tmp.c_str(),ret,ok,liste); 
	if (ok) return ret ; 
	ICMRESULT(result,"Invalid Parameters - Valid are " << liste); 
	result.setRetCode (ARM_DEFAULT_ERR);
	return (qCDS_ADJ)ARM_DEFAULT_ERR;
}



extern CCString ARM_ConvTypeModel(const CCString& typeModel, ARM_result& result)
{
	CCString tmp = typeModel;
	tmp.toUpper ();


    if (tmp == "MMC")
	{
       return(LOCAL_MULTICURVESMODEL_CLASS);
	}
    if (tmp == "CDSOPT")
	{
       return(LOCAL_DEFPROBMODEL_CLASS);
	}
    
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Model type - Valid are  MMC");

	return ARM_DEFAULT_ERR;
}

extern long ARM_ConvTypePayLag (const CCString& aDatesType,ARM_result& result)
{
	CCString tmp = aDatesType;
	tmp.toUpper ();

	if  (tmp == "NXT_PAY_DATE") 
		return 0;

	if  (tmp == "LAST_PAY_DATE") 
		return 1;

	if  (tmp == "SPEC_PAY_DATE") 
		return 2;

	if  (tmp == "NUMERIC_LAG") 
		return 3;

	if  (tmp == "DEFAULT_LAG") 
		return 4;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are NXT_PAY_DATE, LAST_PAY_DATE, SPEC_PAY_DATE, NUMERIC_LAG, DEFAULT_LAG");

	return ARM_DEFAULT_ERR;
}

long ARM_ConvCreditPricerType (const CCString& PricerType, ARM_result& result)
{
	CCString tmp = PricerType;
	tmp.toUpper ();

	if  (tmp == "OLD_PRICER_DEFLEG") 
		return ICM_PRICERDEFAULTLEG;

	if  (tmp == "OLD_PRICER_DEFCDS") 
		return ICM_PRICERDEFAULTTWOLEGS;

	// if  (tmp == "OLD_TREEPRICER") 
	// 	return ICM_PRICERTREE;

	if  (tmp == "MC_CDO2") 
		return ICM_MC_PRICER_CDO2;

	if  (tmp == "OLD_DIFFPRICER") 
		return ICM_DIFFUSIONPRICER;

	if  (tmp == "COPULA") 
		return ICM_PRICERHOMOGENEOUS;

	if  (tmp == "CLD_CDS") 
		return ICM_PRICER_CDS;

	if  (tmp == "CLD_CLN") 
		return ICM_GENERIC_NAME1;
	
	if  (tmp == "CDSOPTION") 
		return ICM_PRICER_CDSOPTION;
	
	if  (tmp == "CDS_INDEX" || tmp == "INDEX") 
		return ICM_PRICER_CDS_INDEX;

	if (tmp == "ANALYTIC_CDO2") 
		return ICM_PRICER_ANALYTIC_CDO2;

	if (tmp == "ANALYTIC_CDO2_STRIKE") 
		return ICM_PRICER_ANALYTIC_CDO2_STRIKE;

	if  (tmp == "INDEXOPTION") 
		return ICM_PRICER_INDEXOPTION;

	if (tmp == "COPULA_SMILE") 
		return ICM_PRICER_HOMOGENEOUS_SMILE;

	if (tmp == "COPULA_MC") 
		return ICM_PRICER_MONTE_CARLO;
	
	if (tmp == "CAPFLOORCMCDS") 
		return ICM_PRICER_CAPFLOORCMCDS;

	if (tmp == "BASKET") 
		return ICM_DIFFUSIONPRICER;

	if (tmp == "CDS_BIN_TREE") 
		return ICM_BINOMIAL_TREE_PRICER_CDS;

// 	if (tmp == "CDS_TRI_TREE") 
// 		return ICM_TRINOMIAL_TREE_PRICER_CDS;

	if (tmp=="SPREADOPTION_BIN_TREE")
		return ICM_BINOMIAL_TREE_PRICER_SPREADOPTION ;

	if (tmp == "CPPI") 
		return ICM_PRICER_CPPI;

	if (tmp == "CORRIDOR") 
		return ICM_PRICER_CORRIDOR;

//	if  (tmp == "INDEXOPTION_SMILE") 
//		return ICM_PRICER_INDEXOPTION_SMILE;

//	if  (tmp == "INDEXOPTION_SABR") 
//		return ICM_PRICER_INDEXOPTION_SABR;

	if  (tmp == "GENERIC_MC_CDO") 
		return ICM_PRICER_GENERIC_MC_CDO;

	if  (tmp == "GENERIC_MC_CDO2") 
		return ICM_PRICER_GENERIC_MC_CDO2;

	if (tmp == "COPULA_RF") 
		return ICM_PRICER_HOMOGENEOUS_SMILE_RF;

	if (tmp == "IR_MODEL")
		return ICM_PRICER_IR_YCMOD;

	if (tmp == "COPULA_SM_FWD")
		return ICM_PRICER_VANILLA_FWD;

	if (tmp == "COPULA_SMILE_COLLAT_FWD")
		return ICM_PRICER_HOMOGENEOUS_SMILE_COLLAT_FWD;

	if (tmp == "CUSTOMIZED_CREDIT_PRICING")
		return ICM_CREDIT_MANAGER;


//	if  (tmp == "INDEX_GEN_BLACK_SCHOLES") 
//		return ICM_PRICER_INDEX_CAPGEN_BLACKSCHOLES;

//	if  (tmp == "INDEX_GEN_BLACK_SCHOLES_STRIKES") 
//		return ICM_PRICER_INDEX_CAPGEN_BLACKSCHOLES_STRIKES;

//	if  (tmp == "INDEX_GEN_SMILE_LNS") 
//		return ICM_PRICER_INDEX_CAPGEN_SMILE;

	if (tmp == "IMPLIED_LOSS_TREE")
		return ICM_IMPLIED_LOSS_TREE_PRICER;

	if (tmp == "LSS_GAP_OPTION")
		return ICM_PRICER_LSS_GAP_OPTION;

	if (tmp == "CPDO") 
		return ICM_PRICER_CPDO;
	if (tmp=="MC_RESTRIKABLE_CDO")
		return ICM_PRICER_TRANCHE_RESTRIKABLE ;

	if (tmp == "BASE") 
		return ICM_BASEPRICER ;

	if (tmp == "MC_CDO") 
		return ICM_PRICER_MC_CDO;

	if (tmp == "HW_CDO_TREE") 
		return ICM_PRICER_TREE_HW_CDO;
	
	if (tmp == "UNKNOWN") 
		return ARM_OBJECT;


	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are CLD_CDS,COPULA,MC_CDO2,CDSOPTION,INDEXOPTION,CAPFLOORCMCDS,CDS_BIN_TREE,CDS_TRI_TREE,GENERIC_MC_CDO,GENERIC_MC_CDO2,COPULA_RF");

	return ARM_DEFAULT_ERR;
}
long ARM_ConvCreditPricerType (const std::string&name)
{
	ARM_result result; 
	long value = ARM_ConvCreditPricerType (name.c_str(), result); 
	if (value==ARM_DEFAULT_ERR) ICMTHROW(ERR_INVALID_ARGUMENT,CCSTringToSTLString(result.getMsg())); 
	return value; 
}
extern long ARM_ConvCreditMezzType (const CCString& MezzType, ARM_result& result)
{
	CCString tmp = MezzType;
	tmp.toUpper ();

	if  (tmp == "RUNNING") 
		return 0;

	if  (tmp == "FUNDED") 
		return 1;

	if  (tmp == "INFINE") 
		return 2;

	if  (tmp == "ZC") 
		return 3;

	if  (tmp == "CMTRANCHE") 
		return 4;
/*
	if  (tmp == "CMAT") 
		return 3;

	if  (tmp == "ZC") 
		return 4;

	if  (tmp == "CMTRANCHE") 
		return 5;
*/	
	result.setRetCode (ARM_DEFAULT_ERR);
//	result.setMsg ("Invalid Parameter - Valid are RUNNING, FUNDED, INFINE, CMAT, ZC, CMTRANCHE");
	result.setMsg ("Invalid Parameter - Valid are RUNNING, FUNDED, INFINE, ZC, CMTRANCHE");

	return ARM_DEFAULT_ERR;
}


extern long ARM_ConvLegType (const CCString& Type, ARM_result& result)
{
	CCString tmp = Type;
	tmp.toUpper ();

	if  (tmp == "NONE") {return qNone_Leg;}
	if  (tmp == "RUNNING") {return qRunning_Leg;}
	if  (tmp == "FUNDED") {return qFunded_Leg;}
	if  (tmp == "INFINE") {return qInfine_Leg;}
	if  (tmp == "ZC") {return qZeroCoupon_Leg;}
	if  (tmp == "CM") {return qConstant_Maturity_Leg;}
	if  (tmp == "RECOVERY") {return qStandart_Recovery_Leg;}
	if  (tmp == "FWD_IR") {return qForward_IRrates_Leg;}
	if  (tmp == "CMCDS") {return qCMCDS_Leg;}
	if  (tmp == "SWAPLEG") {return qSwapLeg;}
	if  (tmp == "INFLEG") {return qInflation_Leg;}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are NONE,RUNNING,FUNDED,INFINE,ZC,CM,RECOVERY,FWD_IR,CMCDS");

	return ARM_DEFAULT_ERR;
}

extern long ARM_ConvCalibType (const CCString& Type, ARM_result& result)
{
	CCString tmp = Type;
	tmp.toUpper ();

	if  (tmp == "BASE") {return qCAL_BASE_CORRELATION;}
	if  (tmp == "BASETS") {return qCAL_BASE_CORRELATION_TS;}
	if  (tmp == "BASETSR") {return qCAL_BASE_CORRELATION_TSR;}
	if  (tmp == "BASE_LINEAR") {return qCAL_BASE_CORRELATION_LINEAR_MATURITY;}
	if  (tmp == "PWC") {return qCAL_PWC_CORREL;}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are BASE,BASETS,BASETSR,BASE_LINEAR,PWC");

	return ARM_DEFAULT_ERR;
}

extern long ARM_ConvCalibMeth (const CCString& Type, ARM_result& result)
{
	CCString tmp = Type;
	tmp.toUpper ();

	if  (tmp == "DICHO") {return qDICHOTOMY;}
	if  (tmp == "NEWTON") {return qNEWTON;}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are DICHO,NEWTON");

	return ARM_DEFAULT_ERR;
}

// End Fonctions de convertions pour le credit -----------------------------------


long ARM_ConvSummitManual (const CCString& aSummitManual, ARM_result& result)
{
	CCString tmp = aSummitManual;

    tmp.toUpper ();

	if ( tmp == "SUMMIT" )
	{
	   return(K_SUMMIT);
	}
	else if ( tmp == "MANUAL" )
	{
	   return(K_MANUAL);
	}
	else if( tmp == "FASTER" )
	{
		/// "scotch"
		return K_SUMMIT+2;
	}
	else if( tmp == "NEWMODE" )
	{
		/// "scotch"
		return K_SUMMIT+3;
	}

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are SUMMIT, 0, MANUAL, 1, FASTER, 2, NEWMODE,3");

	return ARM_DEFAULT_ERR;
}


long ARM_ConvReplicMode( const CCString& ReplicMode, ARM_result& result)
{
	CCString tmp = ReplicMode;
	tmp.toUpper ();

	if (tmp == "CONST_STEP")
		return K_CONST_STEP;

	if (tmp == "CONST_PREC")
		return K_CONST_PREC;

	if (tmp == "ALTN_PREC")
		return K_ALTN_PREC;

	if (tmp == "NUMINT_GL")
		return K_NUMINT_GL;

	if (tmp == "NUMINT_GB")
		return K_NUMINT_GB;

	if (tmp == "SWAP_SETTLE")
		return K_SWAP_SETTLE;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are CONST_STEP, CONST_PREC, ALTN_PREC, NUMINT_GL, NUMINT_GB, SWAP_SETTLE");

	return ARM_DEFAULT_ERR;
}

long ARM_ConvStopMode( const CCString& StopMode, ARM_result& result)
{
	CCString tmp = StopMode;
	tmp.toUpper ();

	if (tmp == "PRICE")
		return K_STOP_PRICE;

	if (tmp == "PRICE_WEIGHT")
		return K_STOP_PRICE_WEIGHT;


	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are PRICE, PRICE_WEIGHT");

	return ARM_DEFAULT_ERR;
}



long ARM_ConvSeasonalityMode(const CCString& seasonMode, ARM_result& result)
{
	CCString tmp = seasonMode;
	tmp.toUpper ();

	if (tmp == "PLUS")
		return K_SEASONADJ_ADDITIVE;
	if (tmp == "MULT")
		return K_SEASONADJ_MULTIPLICATIVE;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are PLUS, MULT");
	return ARM_DEFAULT_ERR;
}


extern long ARM_ConvCopulaPricerType (const CCString& PricerType, ARM_result& result)
{
	CCString tmp = PricerType;
	tmp.toUpper ();

	if  (tmp == "GAUSSIAN") 
		return 1;	// qGAUSSIAN;

	if  (tmp == "STUDENT") 
		return 2;	// qSTUDENT;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are GAUSSIAN, STUDENT");

	return ARM_DEFAULT_ERR;
}

extern long ARM_ConvIndexMethod (const CCString& Method, ARM_result& result)
{
	CCString tmp = Method;
	tmp.toUpper ();

    if(tmp == "AVG")
	{
        return 0 ;
	}
    if(tmp == "HOMOTHETIE")
	{
        return 1 ;
	}
     
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid Method Type - Valid are: AVG Or HOMOTHETIE");

	return ARM_DEFAULT_ERR;
}
extern long ARM_ConvMaturityCapMode (const CCString& MaturityCapMode, ARM_result& result)
{
	CCString tmp = MaturityCapMode;
	tmp.toUpper ();

	if  (tmp == "TRI") 
		return K_TRI;

	if  (tmp == "INFINE") 
		return K_INFINE;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are TRI, INFINE");

	return ARM_DEFAULT_ERR;
}

extern long ARM_ConvMaturityCapCalibrationMode (const CCString& MaturityCapCalibrationMode, ARM_result& result)
{
	CCString tmp = MaturityCapCalibrationMode;
	tmp.toUpper ();

	if  (tmp == "ATM") 
		return K_MCAP_ATM;

	if  (tmp == "FLAT") 
		return K_MCAP_FLAT;

	if  (tmp == "EX_BOUNDARY") 
		return K_MCAP_EX_BOUNDARY;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are ATM, FIRST, EX_BOUNDARY");

	return ARM_DEFAULT_ERR;
}

extern long ARM_ConvSummitFormulae(const CCString& param, ARM_result& result)
{
	CCString tmp = param;
	tmp.toUpper();
	if (tmp == "SUMREP")
		return K_CONV_ADJ_SUMREP;
	if (tmp == "EXP")
		return K_CONV_ADJ_EXP;
	if (tmp == "SUMANA")
		return K_CONV_ADJ_SUMANA;
	if (tmp == "SUMEXP")
		return K_CONV_ADJ_SUMEXP;
	if (tmp == "LIBORORG")
		return K_CONV_ADJ_LIBORORG;
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are SUMREP, EXP, SUMANA, SUMEXP, LIBORORG");

	return ARM_DEFAULT_ERR;
}

//	------------------------------------------------------------------------------------------------
long ARM_ConvUnderlyingMatuStyle  (const CCString& param, ARM_result& result)
{
	qUnderlying_Maturity_Style ret; 
	ICM_EnumsCnv::cnv(CCSTringToSTLString(param),ret); 
	return ret ; 
}
	
extern long ARM_ConvTypeUpOrLow (const CCString& UpOrLow,ARM_result& result)
{
	CCString tmp = UpOrLow;
	tmp.toUpper ();

	if  (tmp == "UP") 
		return 0;

	if  (tmp == "LOW") 
		return 1;

	if  (tmp == "DOWN") 
		return 1;

	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("Invalid Parameter - Valid are UP and LOW, DOWN");

	return ARM_DEFAULT_ERR;
}

long ARM_ConvINFDigitalPayoffType (const CCString& payOffType, ARM_result& result)
{ 
	CCString tmp = payOffType;
	tmp.toUpper ();

    if(tmp == "SIMPLE")
	{
        return K_SIMPLE;
	}
    if(tmp == "PRODUCT")
	{
        return K_PRODUCT;
	}
    
	result.setRetCode (ARM_DEFAULT_ERR);
	result.setMsg ("ARM_ERR: Invalid PayoffType - Valid are : SIMPLE or PRODUCT");

	return ARM_DEFAULT_ERR;
}
