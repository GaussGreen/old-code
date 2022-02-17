/* -------------------------------------------------------------------------
 
   File: %M% 
   Path: %P% 
   Description: Constantes communes a tous les fichiers 
   Created: 99/12/02
   Author: Charles-Emmanuel MUSY
   Modified: %E% %U% 
   Last maintained by: Charles-Emmanuel MUSY
   Revision: %I% 
 
   -------------------------------------------------------------------------
 
   Note:
 
   ------------------------------------------------------------------------- */

#ifndef XL_LOCAL_XLARM_COMMON
#define XL_LOCAL_XLARM_COMMON

#define XLLOCALARM_TITLE_STR				" Arm Local"
#define MERCURE_ARM_XLL_TITLE_STR			"MercureARM"
#define MERCURE_ARM_XLL_NAME				"MercureARM_XLL.xll"

#define XLLOCALARM_APP_NAME					" local_xlarm"
#define XLLOCALARM_LOG_FILENAME				"c:\\temp\\xlarm_local.log"
#define XLLOCALARM_XLL_TOPIC				" Arm Local"
#define XLLOCALARM_XLL_NAME					"local_xlarm.xll"
#define XLLOCALARM_XLL_DESCRIPTION			" Local Arm objects interface"

#define XLLOCALARM_CURVE_GROUP				" Arm Yield and Vol Curves"
#define XLLOCALARM_SEC_GROUP				" Arm Securities"
#define XLLOCALARM_MODEL_GROUP				" Arm Models"
#define XLLOCALARM_GENCALIB_GROUP           " Arm GenCalibration"
#define XLLOCALARM_GENMODELS_GROUP          " Arm GenModels"
#define XLLOCALARM_CLOSEDFORME_GROUP        " Arm ClosedForms"
#define XLLOCALARM_GENCALCULATOR_GROUP      " Arm GenCalculators" 
#define XLLOCALARM_GENNUMMETHODS_GROUP		" Arm GenNumMethods" 
#define XLLOCALARM_GLOBAL_GROUP				" Arm Global" 
#define XLLOCALARM_MKTMANAGER_GROUP			" Arm MktManager" 
#define XLLOCALARM_GETFROMID_GROUP			" Arm GetFromId" 
#define XLLOCALARM_SETTOID_GROUP			" Arm SetToId" 
#define XLLOCALARM_IRFUT_GROUP				" Arm IrFut" /* Unused XLLOCALARM_GROUP*/
#define XLLOCALARM_REV_GROUP				" Arm Reverse" /* Unused XLLOCALARM_GROUP*/
#define XLLOCALARM_REFVAL_GROUP				" Arm RefVal" /* Unused XLLOCALARM_GROUP*/
#define XLLOCALARM_XSTYLE_GROUP				" Arm XStyle" /* Unused XLLOCALARM_GROUP*/
#define XLLOCALARM_OPTION_GROUP				" Arm Option, CapFloor, Swaption"
#define XLLOCALARM_UTIL_GROUP				" Arm Utilities"
#define XLLOCALARM_STRUCT_GROUP				" Arm Structure and PF"
#define XLLOCALARM_CREDIT_GROUP				" Arm Credit"
#define XLLOCALARM_MERCURE_GROUP			" Mercure"

#define EDITOR								"notepad "
#define XML_EDITOR							"iexplore "

#define VIEW_FILE_PREFIX				"VF"
#define VIEW_FILE_SERVER_LOCATION		"/soft/ARM/CORBA/tmp/"
//#define VIEW_FILE_SERVER_LOCATION		"/products/ARM4/CORBA/tmp/" (for Tokyo!)
#define VIEW_FILE_CLIENT_LOCATION		"c:\\temp\\"

#define ZC_FILE_PREFIX					"ICVF"
#define ZC_FILE_SERVER_LOCATION			"/soft/ARM/CORBA/tmp/"
//#define ZC_FILE_SERVER_LOCATION			"/products/ARM4/CORBA/tmp/" (for Tokyo!)
#define ZC_FILE_CLIENT_LOCATION			"c:\\temp\\"

#define HR_FILE_PREFIX					"RHIST_"
#define HR_FILE_SERVER_LOCATION			"/soft/ARM/CORBA/tmp/"
//#define HR_FILE_SERVER_LOCATION			"/products/ARM4/CORBA/tmp/" (for Tokyo!)
#define HR_FILE_CLIENT_LOCATION			"c:\\temp\\"

#define ZC_SUMMIT_FILE_LOCATION			"\\ZC\\"
#define YLD_SUMMIT_FILE_LOCATION		"\\YLD\\"
#define FX_SUMMIT_FILE_LOCATION			"\\FX\\"
#define VOL_SUMMIT_FILE_LOCATION		"\\VOL\\"
#define SMILE_SUMMIT_FILE_LOCATION		"\\SMILE\\"
#define FXVOL_SUMMIT_FILE_LOCATION		"\\FXVOL\\"
#define FXSMILE_SUMMIT_FILE_LOCATION	"\\FXSMILE\\"
#define CORREL_SUMMIT_FILE_LOCATION		"\\INDEXCOR\\"
#define SUMMIT_FILE_LOCATION_HISTO		"histo\\"

#define ID_VIEW_CLOSE                   1002
#define ID_VIEW_TEXT                    1003

#endif	/* XL_LOCAL_XLARM_COMMON */

/* EOF %M% */ 
