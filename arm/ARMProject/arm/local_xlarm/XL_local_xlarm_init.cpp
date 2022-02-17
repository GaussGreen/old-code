#pragma warning(disable : 4786)

#include "XL_local_api.h"
#include "XL_local_xlarm_init.h"
#include <ARM\libarm_local\ARM_local_persistent.h>
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_gp_genericaddin.h>
#include <ARM\libarm_frometk\ARM_local_etoolkit.h>

#include "ARM_local_interglob.h"

#include <glob\calend.h>
#include <glob\dates.h>

#include <GP_Help\gphelp\crmcookies.h>
#include <GP_Base\gpbase\env.h>
using ARM::ARM_CRMCookies;

#include <GP_Infra\gpinfra\gramfunction.h>
#include <GP_Help\gphelp\splashwindow.h>

#include "XL_local_xlarm_resource.h"

ARM_Date TMP_DATE;

void XLLOCALARM_InitFileRead ()
{
	LOCALARM_IniFileRead ();
}


void XLLOCALARM_InitEnvVariables ()
{
	ARM::LOCALARM_InitEnvVariables ();
}

/// function that is responsible to initialise the various
/// global variables of Grand Prix such as the grammar function
/// table
void XLLOCALARM_InitGrandPrix ()
{
	ARM_GenericAddinsDesc::CreateTheAddingTable();
	ARM::ARM_GramFunction::InitTheFuncTable ();
	ARM::ARM_SplashWindow::SetResourceHandle(IDR_SPLASHWIN); 
	ARM::ARM_SplashWindow::SetModuleHandle( g_hInst);
}

long XLLOCALARM_PersistentListsInit ()
{
	return LOCALARM_PersistentListsInit ();
}

long XLOCALLARM_PersistentListsClear ()
{
	return LOCALARM_PersistentListsClear ();
}

long XLOCALLARM_deconnexioneToolkit ()
{
	return kill_etoolkit();
}

long XLLOCALARM_PersistentListsDelete ()
{
	return LOCALARM_PersistentListsDelete ();
}

/// to end customer relationship management tracing
/// currently only done for the generic pricer
void XLLOCALARM_EndCRMTracing()
{
	ARM_CRMCookies.Instance()->EndTracing();
}


void XLLOCALARM_ReleaseGrandPrix()
{
	ARM_GenericAddinsDesc::ReleaseTheAddingTable();
	ARM::ARM_GramFunction::ReleaseTheFuncTable();
}

void XLLOCALARM_InitCalendar ()
{
    ARM_Calendar* cal;

    int rc;
    char calenDarFileName[200];


    rc = ARM_GetCalendarFile(calenDarFileName);

    if ( rc == RET_OK )
    {
       cal = new ARM_Calendar(calenDarFileName);

       TMP_DATE.SetCalendar(cal);
    }
}

void XLLOCALARM_FreeAllObjects ()
{
	ARM_result result;

	ARMLOCAL_FreeAllObjects (result);
}


// EOF %M%