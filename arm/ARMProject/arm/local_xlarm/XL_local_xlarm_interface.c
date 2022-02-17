#include "XL_local_xlarm_interface.h"

#include "XL_local_xlarm_init.h"

#include <libCCmessage\CCmessage.h>
#include "ARM_local_help.h"

#ifdef DEBUG
char* XLLOCALARM_revision = "%I% debug";
#else
char* XLLOCALARM_revision = "%I%";
#endif  // DEBUG

static FILE* fic;


LPSTR g_rgWorksheetFuncs[][g_rgWorksheetFuncsCols] =
{
#include "XL_local_xlarm_interface_YieldCurve.h"
#include "XL_local_xlarm_interface_Utilities.h"
#include "XL_local_xlarm_interface_Models.h"
#include "XL_local_xlarm_interface_Securities.h"
#include "XL_local_xlarm_interface_Structures.h"
#include "XL_local_xlarm_interface_Credit.h"
#include "XL_local_xlarm_interface_GenPricer.h"
#include "XL_local_xlarm_interface_GenCalib.h"
#include "XL_local_xlarm_interface_Calculators.h"
#include "XL_local_xlarm_interface_Inflation.h"
#include "XL_local_xlarm_interface_Option.h"
#include "XL_local_xlarm_interface_ClosedForms.h"
#include "XL_local_xlarm_interface_XXXProject.h"
#include "../MercureARM_XLL/XL_local_xlarm_interface_Mercure.h"
};



int g_rgWorksheetFuncsRows = sizeof( g_rgWorksheetFuncs ) / sizeof( g_rgWorksheetFuncs[0] );

LPSTR g_rgCommandFuncs[][g_rgCommandFuncsCols] =
{
		{       " Local_DisplayErrorMessage",
                " A",
                " Local_DisplayErrorMessage",
                " ",
                " 2",
                XLLOCALARM_UTIL_GROUP,
                " "
        },
		{       " Local_FreeAllObjects",
                " A",
                " Local_FreeAllObjects",
                " ",
                " 2",
                XLLOCALARM_UTIL_GROUP,
                " "
        },
		{		" Local_View",
                " A",
                " Local_View",
                " ",
                " 2",
                XLLOCALARM_UTIL_GROUP,
                " "
        },
		{		" Local_View_XML",
                " A",
                " Local_View_XML",
                " ",
                " 2",
                XLLOCALARM_MERCURE_GROUP,
                " "
        },
		{       " XLLOCALARM_Exit",
                " A",
                " XLLOCALARM_Exit",
                " ",
                " 2",
                XLLOCALARM_UTIL_GROUP,
                " "
        },
		{       " Local_ARM_Help",
                " A",
                " Local_ARM_Help",
                " ",
                " 2",
                XLLOCALARM_UTIL_GROUP,
                " "
        },
		{       " Local_ARM_ProdConnect",
                " A",
                " Local_ARM_ProdConnect",
                " ",
                " 2",
                XLLOCALARM_UTIL_GROUP,
                " "
        },
		{       " Local_ARM_RepliConnect",
                " A",
                " Local_ARM_RepliConnect",
                " ",
                " 2",
                XLLOCALARM_UTIL_GROUP,
                " "
        },
		{       " Local_ARM_InfocConnect",
                " A",
                " Local_ARM_InfocConnect",
                " ",
                " 2",
                XLLOCALARM_UTIL_GROUP,
                " "
        },
		{       " Local_ARM_CalypsoDevConnect",
                " A",
                " Local_ARM_CalypsoDevConnect",
                " ",
                " 2",
                XLLOCALARM_UTIL_GROUP,
                " "
        },
		{       " Local_ARM_CalypsoProdConnect",
                " A",
                " Local_ARM_CalypsoProdConnect",
                " ",
                " 2",
                XLLOCALARM_UTIL_GROUP,
                " "
        },
		{       " Local_ARM_CalypsoRecConnect",
                " A",
                " Local_ARM_CalypsoRecConnect",
                " ",
                " 2",
                XLLOCALARM_UTIL_GROUP,
                " "
        },
		{       " Local_ARM_RecConnect",
                " A",
                " Local_ARM_RecConnect",
                " ",
                " 2",
                XLLOCALARM_UTIL_GROUP,
                " "
        },
		{       " Local_ARM_ShutDownETK",
                " A",
                " Local_ARM_ShutDownETK",
                " ",
                " 2",
                XLLOCALARM_UTIL_GROUP,
                " "
        },
		{       " Local_ARM_SwitchToETK",
                " A",
                " Local_ARM_SwitchToETK",
                " ",
                " 2",
                XLLOCALARM_UTIL_GROUP,
                " "
        },
		{       " Local_ARM_SwitchToETK_WithFallBack",
                " A",
                " Local_ARM_SwitchToETK_WithFallBack",
                " ",
                " 2",
                XLLOCALARM_UTIL_GROUP,
                " "
        },
		{       " Local_ARM_SwitchToFLATFILE",
                " A",
                " Local_ARM_SwitchToFLATFILE",
                " ",
                " 2",
                XLLOCALARM_UTIL_GROUP,
                " "
        },
		{       " Local_ARM_deconnection_etoolkit",
                " A",
                " Local_ARM_deconnection_etoolkit",
                " ",
                " 2",
                XLLOCALARM_UTIL_GROUP,
                " "
        },
		{       " Local_ARM_SwitchToWSETK",
                " A",
                " Local_ARM_SwitchToWSETK",
                " ",
                " 2",
                XLLOCALARM_UTIL_GROUP,
                " "
        },
};

int g_rgCommandFuncsRows = sizeof( g_rgCommandFuncs ) / sizeof( g_rgCommandFuncs[0] );

LPSTR g_rgMenu[][g_rgMenuCols] =
{
        {		" &Arm Local",
                " ",
                " ",
                " ARM Excel client",
                " "
		},
		{		" Error &Message",
                " Local_DisplayErrorMessage",
                " ",
                " Display error message",
                " " 
		},
        {		" &View",
                " Local_View",
                " ",
                " View", 
                " "
		},
        {		" View &XML",
                " Local_View_XML",
                " ",
                " View XML output", 
                " "
		},
		{		" Switch to ETK",
                " Local_ARM_SwitchToETK",
                " ",
                " Retrieves MarketData with ETK",
                " " 
		},
		{		" Switch to Flat File",
                " Local_ARM_SwitchToFLATFILE",
                " ",
                " Retrieves MarketData with ETK",
                " " 
		},
		{		" Connexion Prod",
                " Local_ARM_ProdConnect",
                " ",
                " Connexion à la base Summit de prod",
                " " 
		},
		{		" Connexion Replication",
                " Local_ARM_RepliConnect",
                " ",
                " Connexion à la base Summit de Rplication",
                " " 
		},
		{		" Connexion Infocentre",
                " Local_ARM_InfocConnect",
                " ",
                " Connexion à la base Summit Infocentre",
                " " 
		},
		{		" Connexion Recette",
                " Local_ARM_RecConnect",
                " ",
                " Connexion à la base Summit de recette",
                " " 
		},
		{		" Déconnexion ETK",
                " Local_ARM_deconnection_etoolkit",
                " ",
                " Déconnexion de la base active",
                " " 
		},
		{		" Calypso: Connexion dev ",
                " Local_ARM_CalypsoDevConnect",
                " ",
                " Connexion à la base Calypso de dev ",
                " " 
		},		
		{		" Calypso: Connexion rec ",
                " Local_ARM_CalypsoRecConnect",
                " ",
                " Connexion à la base Calypso de recette ",
                " " 
		},		
		{		" Calypso: Connexion prod ",
                " Local_ARM_CalypsoProdConnect",
                " ",
                " Connexion à la base Calypso de production ",
                " " 
		},		
		{		" &Clean ARM environment",
                " Local_FreeAllObjects",
                " ",
                " Clean ARM environment",
                " " 
		},
/*		{		" &Exit from ARM",
                " XLLOCALARM_Exit",
                " ",
                " Exit from ARM Excel",
                " " 
		},
*/		{		" &Help",
                " Local_ARM_Help",
                " ",
                " ARM Help",
                " " 
		},
};


/*---- Calculate the number of items -----*/

int g_rgMenuRows = sizeof(g_rgMenu)/sizeof(g_rgMenu[0]);


int XLLOCALARM_xlAutoClose (void)
{
	XLLOCALARM_PersistentListsDelete ();
	XLOCALLARM_deconnexioneToolkit ();
	XLLOCALARM_EndCRMTracing ();
	XLLOCALARM_ReleaseGrandPrix ();

    MSG_printf_message (MSG_INFO, "Fin de xlarm local version %s", "%I%");

    if (fic)
    {
       fclose (fic);
    }

	return 1;
}


int XLLOCALARM_xlAutoOpen(void)
{
	char* user = getenv ("USERNAME");
	char* debug_mode = getenv ("ARM_DEBUG_MODE");

	MSG_set_stdlog (stderr);
	
    if ((fic = fopen (XLLOCALARM_LOG_FILENAME, "w")))
	{
	   setbuf (fic, NULL);
	
       MSG_set_stdlog (fic);
	}

	if (debug_mode)
	{
	   MSG_enable_debug ();
		
       MSG_enable_verbose ();
	}
	
    MSG_set_stderr_quiet ();

	MSG_printf_message (MSG_INFO, "Bienvenue %s dans xlarm local version %s", user ? user : DEFAULT_NULL_STRING, "%I%");
		
	XLLOCALARM_PersistentListsInit ();

	XLLOCALARM_InitFileRead();

	XLLOCALARM_InitCalendar();

	XLLOCALARM_InitEnvVariables();

    XLLOCALARM_InitGrandPrix();

	XLLOCALARM_xlFirstCalculate (); 

	return(1);
}



__declspec(dllexport) int WINAPI XLLOCALARM_xlFirstCalculate(void)
{
	XLOPER xResult;

	Excel4 (xlcCalculateNow, &xResult, 0);
	
	return 0;
}



__declspec(dllexport) int WINAPI XLLOCALARM_xlCalculateNow (void)
{
	XLOPER xResult;

	Excel4 (xlcCalculateDocument, &xResult, 0);	

	Excel4 (xlcCalculateDocument, &xResult, 0);

	return 0;
}


/*------------------------------------------------------------------------------------*/
/*---- End Of File ----*/ 
