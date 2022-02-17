#ifndef ARM_LOCAL_PARSEXML_COMMON_H
#define ARM_LOCAL_PARSEXML_COMMON_H


#include <ARM\libarm_local\firstToBeIncluded.h>

#include <libCCTools++/CCString.h>
#include <ARM\libarm_local\ARM_local_glob.h>

#include <ARM\local_xlarm\ARM_local_interglob.h>
/// headers to remove definition of va_start and va_end as this is redefined later on!
/// handle this with care after sorting out why this is so!!
#include <ARM\libarm_local\undef_va_vars.h>

#include <libCCdate\CCdate.h>
#include <util\fromto.h>

#define SUMMIT_DEFAULT_CURVE "MO"

class ICM_Credit_Index;
class ICM_ProportionsInfo;
class ARM_Currency;

void Local_XLDATE2ARMDATE (double xldate, char* myArmDate);
double XML_doubleNodeTreating(void* listItem, const CCString& nodeName);
char* ExtractCorr(char* valeur);
int SummitStub2ARMStub(char* stub);
void Local_XML_TRACE(CCString fichier, CCString message);
int Summit_InterpMethod(const char* interpol);
wchar_t * constchar2wchar(const char* input);
void CvtStrProp(char* input,char**& labels, ARM_Vector*& prop,ARM_Vector*& correl_low,ARM_Vector*& correl_up,ICM_Credit_Index**& vIndex);
// void CvtStrPropInfo(char* input,ICM_ProportionsInfo*& info);
bool ExistXMLNode (void* node, string nodeName);
CCString GetIndexNameFromARMIndex(ARM_INDEX_TYPE index);

void DeduceRefDateAndStub(ARM_Date StartDate,
						  ARM_Date EndDate,
						  ARM_Date StartStubDate,
						  ARM_Date EndStubDate,
						  const int& payfreq,
						  const int & RefDay,
						  const std::string&  ccy,
						  char* ReferenceDate,
						  int& stub);

class SUMMIT_NAMES
{
	private:
	static vector<string> itsListNames ; 

	public :

    static string STARTDATE;
	static string ENDDATE;
	static string FIXINGDATE;
	static string INTDAYS;
	static string RATE;
	static string FXRATE;
	static string FORWARD;
	static string SPREAD;
	static string DECOMPRATE;
	static string NOTIONAL;
	static string PAYDATE;
	static string INTERIMINTEREST;
	static string FLOWS;
	static string STRIKE;
	static string TYPE;
	static string DAYS;
	static string ZERORATE;
	static string DISCFACTOR;
	static string PV;
	static string AICRATE;

	static string TENOR_FWD;
	static string PARTICIPATION_RATE;
	static string FLOOR;
	static string CAP;
	
	static string EVENT_DATE;
	static string EVENT_ADATE;

	static vector<string>& list_names(void);
	static int nbParameters();
};


#endif 

