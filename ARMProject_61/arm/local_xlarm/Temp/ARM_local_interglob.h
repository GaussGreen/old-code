#ifndef _ARM_LOCAL_INTERGLOB_H
#define _ARM_LOCAL_INTERGLOB_H


#include <ARM\libarm_local\ARM_local_persistent.h>
//#include <ARM\libarm_local\ARM_local_init.h>

#define ARM_DEFAULT_ERR			-990  

#ifdef STL_WIN32
#include <map>
#include <vector>
#define MAP		std::map
#else
#include <map.h>
#include <vector.h>
#define MAP		map
#endif	// STL_WIN32

/* mimic the class AUTO_MODE from modes.h   */
struct INTERFACE_AUTO_MODE
{
    enum _type
    {
        AUTO_MODE_NONE      = 1,
        AUTO_MODE_IRG       = AUTO_MODE_NONE << 1,          // calibrate on swaptions
        AUTO_MODE_IDX       = AUTO_MODE_IRG << 1,           // calibrate on caps
        AUTO_MODE_SWOPT_IRG = AUTO_MODE_IDX <<1,            // calibrate on swaptions, optimize on caps
        AUTO_MODE_IRG_SWOPT = AUTO_MODE_SWOPT_IRG << 1,     // calibrate on caps, optimize on swaptions
        AUTO_MODE_DUAL      = AUTO_MODE_IRG_SWOPT << 1,     // calibrate swaptions and captions, 2 vols curves
        AUTO_MODE_AUTOMEAN  = AUTO_MODE_DUAL << 1,          // optimize the mean reversion rate
        AUTO_MODE_BASISDISC = AUTO_MODE_AUTOMEAN << 1,		// basis discounting
		AUTO_MODE_CAPDIAG   = AUTO_MODE_BASISDISC << 1,
		AUTO_MODE_DIAG		= AUTO_MODE_CAPDIAG << 1
    };
};




//extern ARMLOCAL_Init* armlocal_init;

#define XLLOCALARM_INIFILE_NAME		"C:\\Program Files\\Arm\\local_xlarm.ini"
#define DEFAULT_CURRENCY			"EUR"
#define DATA_FOLDER			        "P:\\Samba\\arm-import"


#define ARM_OK					0
#define ARM_KO					-1





#define ARM_NULL_OBJECT_ID		-1

/// The default size for error message is 1000
// 
//	JLA: we don't want handle message size limitation.
//
// #define ARM_RESULT()					result.setRetCode(ARM_KO); \
// 										char* var = new char[1000]; \
// 										x.GetErrorMessage(var); \
// 										msg = var; \
// 										result.setMsg(msg); \
// 										long retCodeTmp = result.getRetCode (); \
// 										delete var;  \
// 										return (retCodeTmp); 
// 

#define ARM_RESULT()			{ result.setRetCode(ARM_KO); \
								std::string var = x.GetErrorString(); \
								msg = var.c_str(); \
								result.setMsg(msg); \
								long retCodeTmp = result.getRetCode ();  \
								return (retCodeTmp); }

//	This macro allows to fill an ARM_Result with a stream syntax
//	usage : 
//	ICMRESULT( result, "Some message, i="<<i<<" and string="<<myArg" ) ;
#ifndef	ICMRESULT
	#define ICMRESULT(argICMRESULT,msgICMRESULT) \
	{ \
		std::stringstream _sstr; _sstr<<msgICMRESULT<<std::endl; \
		argICMRESULT.setMsgString(_sstr.str()); \
	}
#endif	// ICMRESULT


//extern ARM_message* BuildErrorMessage (const string& errValue, bool storeInLog = false);

extern string GetCurStringCellCoordinates ();
extern string GetLastCurCellEnvValue ();
extern string GetLastCurActiveCellEnvValue ();

extern string LocalMakeObjectId (long objId, const string& objClass);
std::string	LocalGetObjectLabel(long objId,const std::string& objClass); 
extern string LocalGetStringObjectClass (const string& stringObjectId);
extern string LocalGetStringObjectClassAndError(const string& stringObjectId);
extern long LocalGetNumObjectId (const string& stringObjectId);

void LocalSetCurCellEnvValue (const string& curClass, long objId);
void LocalSetCurCellEnvValue (const string& curClass, long objId,const string& caller);

extern void SetCurCellErrValue (const string& errValue);
extern void FreeCurCellErr ();
extern long FreeCurCellContent ();
extern long FreeCurCellContentWithoutFreeingARMObject();
extern long FreeCurCellContent (const string & envVal);
extern string GetEnvVar (const string& cellCoord);
extern string GetErrValue (const string& cellCoord);
extern void PrintErrorsList ();

extern long LOCALARM_PersistentListsInit ();
extern long LOCALARM_PersistentListsDelete ();
extern long LOCALARM_PersistentListsClear ();
extern void LOCALARM_IniFileRead(void);
extern long LOCALARM_DeconnexionEToolkit ();



extern long IsCall (const string& aString);
extern long IsPut (const string& aString);
extern long IsCap (const string& aString);
extern long IsFloor (const string& aString);
extern long IsReceive (const string& aString);
extern long IsPay (const string& aString);
extern long IsBermudan (const string& aString);


#endif	// _ARM_LOCAL_INTERGLOB_H

// EOF %M%
