#include "ARM_CORBA_INIT.h"

#include "ARM_interglob.h"
#include "ARM_glob.h"

#include <iosfwd>
#include <fstream>

///using std::ostream;

#define ARM_PB_CORBA_LOG_FILENAME	"c:\\temp\\xlarmcorba.log"





ARM_CorbaRequestModule::ARM_CORBA_REQUEST_CALL* CORBA_OBJECT_INTERFACE;

static int ARM_CORBA_initialized = 0;
static int ARM_VIRTUAL_DISCONNECT = 0;



void ARM_PB_TRACE(CORBA::ORB_var& orb)
{
	ofstream file(ARM_PB_CORBA_LOG_FILENAME); 


    streambuf* sbuf = cout.rdbuf(); 
   
   

	cout.rdbuf()->setbuf((char *) file.rdbuf(), 10000); 
	
	// ...

//  Not in ORBIXE2A orb->setDiagnostics(255);

	//...

    // cout.rdbuf(); 
}



long ARM_CORBA_init(void)
{
	if ( ARM_CORBA_initialized == 1 )
	{
	   return(ARM_OK);
	}

	int argc = 0;

	CORBA::ORB_var orb = CORBA::ORB_init(argc, NULL, "Orbix");

	MSG_printf_message(MSG_INFO, "\n\n===========> ARM_init: ORB_init OK <============\n\n");

	CORBA::Object_var obj = orb->resolve_initial_references ("NameService");
	CosNaming::NamingContext_var nc = CosNaming::NamingContext::_narrow (obj);

	MSG_printf_message(MSG_INFO, "\n\n ---------------> ARM_init: name service called \n\n");

	try
	{
//		if (nc)
		{
			CosNaming::Name appName;
			appName.length (1);
			appName[0].id = CORBA::string_dup ("ARM_CORBA_REQUEST_CALL");
			appName[0].kind = CORBA::string_dup ("");
			CORBA::Object_var cObj = nc->resolve (appName);

			CORBA_OBJECT_INTERFACE = ARM_CorbaRequestModule::ARM_CORBA_REQUEST_CALL::_narrow (cObj);

			ARM_CORBA_initialized = 1;

			MSG_printf_message (MSG_INFO, "\n\n *---> ARM_init:  CORBA_OBJECT_INTERFACE created \n\n");

			return(ARM_OK);
		}
	}

	catch(CORBA::SystemException& e)
	{
		MSG_printf_message (MSG_INFO, "====>???? Corba: unexpected system exception: %s", (const char*)&e);

		return(ARM_KO);
	}

	return(ARM_KO);
}



long ARM_CORBA_exit()
{
	ARM_ExitArm();

	MSG_printf_message (MSG_INFO, "\n====> Arm server killed");

	if(CORBA_OBJECT_INTERFACE)
	{
		delete CORBA_OBJECT_INTERFACE;
	}

	MSG_printf_message (MSG_INFO, "\n====> CORBA_OBJECT_INTERFACE deleted");

	return(ARM_OK);
}



void ARM_SetVirtualDisconnect ()
{
	ARM_VIRTUAL_DISCONNECT = 1;
}

void ARM_UnsetVirtualDisconnect ()
{
	ARM_VIRTUAL_DISCONNECT = 0;
}

int ARM_IsVirtualDisconnect ()
{
	return(ARM_VIRTUAL_DISCONNECT == 1);
}

int ARM_IsCORBA_initialized ()
{
	return (ARM_CORBA_initialized == 1);
}
// EOF %M%
