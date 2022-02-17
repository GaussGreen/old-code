
#include <ARM\libarm_frometk\arm_local_parsexml_common.h>
#include <ICMKernel\inst\icm_credit_index.h>
#include "ICMKernel\inst\icm_ftd.h"
#include <ICMKernel\util\icm_utils.h>

#ifndef XML_DEFINE
#define XML_DEFINE
#import <msxml3.dll> raw_interfaces_only
using namespace MSXML2;
#endif

#include <stdio.h>

string SUMMIT_NAMES::STARTDATE="StartDate";
string SUMMIT_NAMES::ENDDATE="EndDate";
string SUMMIT_NAMES::FIXINGDATE="FixingDate";
string SUMMIT_NAMES::INTDAYS="IntDays";
string SUMMIT_NAMES::RATE="Rate";
string SUMMIT_NAMES::FXRATE="FxRate";
string SUMMIT_NAMES::FORWARD="Forward";
string SUMMIT_NAMES::SPREAD="Spread";
string SUMMIT_NAMES::DECOMPRATE="DecompRate";
string SUMMIT_NAMES::NOTIONAL="Notional";
string SUMMIT_NAMES::PAYDATE="PayDate";
string SUMMIT_NAMES::INTERIMINTEREST="InterimInterest";
string SUMMIT_NAMES::FLOWS="Flows";
string SUMMIT_NAMES::STRIKE="Strike";
string SUMMIT_NAMES::TYPE="Type";
string SUMMIT_NAMES::DAYS="Days";
string SUMMIT_NAMES::ZERORATE="ZeroRate";
string SUMMIT_NAMES::DISCFACTOR="DiscFactor";
string SUMMIT_NAMES::PV="PV";
string SUMMIT_NAMES::AICRATE="AICRate";

string SUMMIT_NAMES::TENOR_FWD="Amount1";
string SUMMIT_NAMES::PARTICIPATION_RATE="Rate1";
string SUMMIT_NAMES::FLOOR="Rate2";
string SUMMIT_NAMES::CAP="Rate3";

string SUMMIT_NAMES::EVENT_DATE="Date";
string SUMMIT_NAMES::EVENT_ADATE="ADate";

vector<string> SUMMIT_NAMES::itsListNames ; 

vector<string>&
SUMMIT_NAMES::list_names() 
{
	if (itsListNames.empty())
	{
		itsListNames.push_back(STARTDATE); 
		itsListNames.push_back(ENDDATE); 
		itsListNames.push_back(FIXINGDATE); 
		itsListNames.push_back(INTDAYS); 
		itsListNames.push_back(RATE); 
		itsListNames.push_back(FXRATE); 
		itsListNames.push_back(FORWARD); 
		itsListNames.push_back(SPREAD); 
		itsListNames.push_back(DECOMPRATE); 
		itsListNames.push_back(NOTIONAL); 
		itsListNames.push_back(PAYDATE); 
		itsListNames.push_back(INTERIMINTEREST); 
		itsListNames.push_back(FLOWS); 
		itsListNames.push_back(STRIKE); 
		itsListNames.push_back(TYPE); 
		itsListNames.push_back(DAYS); 
		itsListNames.push_back(ZERORATE); 
		itsListNames.push_back(DISCFACTOR); 
		itsListNames.push_back(PV); 
		itsListNames.push_back(AICRATE); 
		itsListNames.push_back(TENOR_FWD); 
		itsListNames.push_back(PARTICIPATION_RATE); 
		itsListNames.push_back(FLOOR); 
		itsListNames.push_back(CAP); 
		itsListNames.push_back(EVENT_DATE); 
		itsListNames.push_back(EVENT_ADATE); 
	}
	return itsListNames; 
}

int
SUMMIT_NAMES::nbParameters()
{
	return list_names().size(); 
}

double XML_doubleNodeTreating(void* listItem, const CCString& nodeName)
{
	MSXML2::IXMLDOMNode* theNode;
	MSXML2::IXMLDOMNode* listItem2 = (MSXML2::IXMLDOMNode*)listItem;

	double res = 0.0;

	listItem2->selectSingleNode(_bstr_t(nodeName), &theNode);
	if (theNode!=NULL)
	{
		BSTR resultat = NULL;
		theNode->get_text(&resultat);

		_bstr_t ff(resultat,false);
		char * ff1=(char *)ff;

		res = atof((const char*) ff1);

		theNode->Release();
		theNode=NULL;
		if (resultat) SysFreeString(resultat);
	}

	return res;
}

int SummitStub2ARMStub(char* stub)
{
	if (strcmp(stub,"SS") == NULL)
		return K_SHORTSTART;
	else if (strcmp(stub,"LS") == NULL)
		return K_LONGSTART;
	else if (strcmp(stub,"SE") == NULL)
		return K_SHORTEND;
	else 
		return K_LONGEND;
}

char* ExtractCorr(char* valeur)
{
	char* tmp = new char[255];
	memset(tmp,'\0',sizeof(char)*255);
	strcpy(tmp,valeur+sizeof(char)*9 /* size ISSRCORR */ );
	int pos = strcspn( tmp, "/" );
	tmp[pos] = '\0';

	return (tmp);
}

void Local_XLDATE2ARMDATE (double xldate, char* myArmDate)
{
	int y, m, d;
	char tmpDay[4], tmpMonth[4], tmpYear[5];
	long long_xldate = (long)xldate;

	DAT_ssdate_to_struct ((double)long_xldate, &y, &m, &d);

	if(d < 10)
	{
		sprintf (tmpDay, "0%1d.", d);
	}
	else
	{
		sprintf (tmpDay, "%2d.", d);
	}
	
	if(m < 10)
	{
		sprintf (tmpMonth, "0%1d.", m);
	}
	else
	{
		sprintf (tmpMonth, "%2d.", m);
	}

	sprintf (tmpYear, "%4d", y);

	sprintf(myArmDate,"%s%s%s",tmpDay,tmpMonth,tmpYear);
}

void Local_XLDATE2ARMDATE (double xldate, ARM_Date& armDate)
{
	int y, m, d;
	// jla: ???  double->long->double ?
	long long_xldate = (long)xldate; 
	DAT_ssdate_to_struct ((double)long_xldate, &y, &m, &d);
	armDate.ChangeDate(d,m,y); 
}

void Local_XML_TRACE(CCString fichier, CCString message)
{
	#ifdef _DEBUG

	FILE *stream;
	char* tmp = fichier.c_str();

	char buf[100000];
	memset(buf,'\0',100000*sizeof(char));

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

	char* tmp2 = message.c_str();

	sprintf (buf, "%s:%s:%s:%s\n",sHour,sMin,sSec,tmp2);

	CCString Fic = (CCString) buf;

    stream  = fopen( tmp, "a+" );

	if (stream==NULL)
	    stream  = fopen( tmp, "w+" );

   int numwrite = fwrite( buf, sizeof( char ), Fic.GetLen(), stream );

   delete[] tmp;
   delete[] tmp2;

   fclose(stream);

   #endif
}


int Summit_InterpMethod(const char* interpol)
{
	if (strcmp(interpol,"CONTINUOUS") == NULL)
		return K_CONTINUOUS;
	else if (strcmp(interpol,"LINEAR") == NULL)
		return K_LINEAR;
	else if (strcmp(interpol,"SQRT") == NULL)
		return K_SQRT;
	else if (strcmp(interpol,"SLOPELIN") == NULL)
		return K_SLOPELIN;
	else if (strcmp(interpol,"SLOPESQRT") == NULL)
		return K_SLOPESQRT;
	else if (strcmp(interpol,"CPILINEAR") == NULL)
		return K_CPILINEAR;
	else if (strcmp(interpol,"CPISTEPWISE") == NULL)
		return K_CPISTEPWISE;
	else if (strcmp(interpol,"ZCLINEAR") == NULL)
		return K_ZCLINEAR;
	else if (strcmp(interpol,"ZCCTFWD") == NULL)
		return K_ZCCTFWD;
	else if (strcmp(interpol,"CPISTEPWISESTART") == NULL)
		return K_CPISTEPWISESTART;
	else if (strcmp(interpol,"CPISTEPWISEEND") == NULL)
		return K_CPISTEPWISEEND;
	else if (strcmp(interpol,"K_CPISTEPWISEMIDDLE") == NULL)
		return K_CPISTEPWISEMIDDLE;
	else 
		return K_LINEAR;
}



wchar_t * constchar2wchar(const char* input)
{
	wchar_t * xmlWCharText = NULL;


	if ( input == NULL )  
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					   "constchar2wchar : MB_ERR_INVALID_CHARS");

	long size = strlen(input);
	long size_max = 50*1024*1024;	//(environ 50 Mo)
	int output = 0;

	if ((size<3)||(size>size_max))  throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
									"constchar2wchar : MB_ERR_INVALID_CHARS");
	
	if ( input[size-3] == '>' )
	   size = size-2;
	else if ( input[size-2] == '>' )
		size = size-1;

	try
	{
        xmlWCharText = new wchar_t[size + 1];

	    if ( xmlWCharText == NULL ) 
           throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
									    "constchar2wchar : Memory Overflow");

        output = MultiByteToWideChar(CP_ACP, 0, (LPCSTR)input, -1, xmlWCharText, size);
	}

	catch(...)
	{
	    if (output == ERROR_INSUFFICIENT_BUFFER)
		    throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				     "constchar2wchar : ERROR_INSUFFICIENT_BUFFER");
	    else if (output == ERROR_INVALID_FLAGS)
		    throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				     "constchar2wchar : ERROR_INVALID_FLAGS");
	    else if (output == ERROR_INVALID_PARAMETER)
		    throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				     "constchar2wchar : ERROR_INVALID_PARAMETER");
	    else if (output == ERROR_NO_UNICODE_TRANSLATION)
		    throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				     "constchar2wchar : ERROR_NO_UNICODE_TRANSLATION (invalid characters )");
	    else throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				     "constchar2wchar : unknown exception");

	}

    xmlWCharText[size] = '\0';

	return(xmlWCharText);
}



void CvtStrProp(char* input,char**& labels, ARM_Vector*& prop,
                ARM_Vector*& correl_low,
                ARM_Vector*& correl_up,
                ICM_Credit_Index**& vIndex)
{
	if (strlen(input) == 0)
	{
		labels = NULL;
		prop = NULL;
		return;
	}

	int maxsize = strlen(input);	
	int nbseparator = 0;
	int i = 0,j=0;
	int nbelt = 0;
	int il=0;

	for (i=0;i<maxsize;i++)
	{ if (input[i] == ';') nbseparator++;}

	nbelt = nbseparator + 1;

	prop = new ARM_Vector(nbelt,0.);
	correl_low = new ARM_Vector(nbelt,CREDIT_DEFAULT_VALUE);
	correl_up = new ARM_Vector(nbelt,CREDIT_DEFAULT_VALUE);
	vIndex = new ICM_Credit_Index*[nbelt];

	labels = new char*[nbelt];
	int indice = 0;

	int size_Index = 125;
	// char** NomIndex = new char*[size_Index];
	std::vector<std::string> NomIndex(size_Index); 

	CCString chaine = (CCString) input;
	vector<CCString> list1;

	chaine.Parser(';',list1);	

	//for (il=0; il<size_Index; il++) NomIndex[il] = new char[100];


	for (i=0;i<list1.size();i++)
	{
		vector<CCString> list2;
		list1[i].ParserWithEmpty(':',list2);	
		labels[i] = new char[255];

		CCString chaine = list2[1] + (CCString) "_" + list2[0] + (CCString) "_" +  list2[3];

		strcpy(labels[i],(const char*)chaine);
		prop->InitElt(i,atof((const char*)list2[5])/100.);

		if (!(list2[7] == "") && !(list2[8] == ""))
		{
			correl_low->InitElt(i,atof((const char*)list2[7])/100.);
			correl_up->InitElt(i,atof((const char*)list2[8])/100.);
			//delete[] vIndex; vIndex=NULL;return;
		}

		vIndex[i] = NULL;		

		for (int il=0; il<size_Index; il++)
			NomIndex[il]= CCSTringToSTLString(list2[6]);
		ARM_Vector YT(1); YT.Elt(0) = 5.0;
		ARM_Vector spread(1); spread.Elt(0) = 0.;
		vIndex[i] = new ICM_Credit_Index(KACTUAL_360, 
									 K_QUARTERLY, 
									 K_QUARTERLY, 
									 YT,
									 ARM_DEFAULT_COUNTRY,
									 NomIndex[0],
									 NomIndex,
									 qAVERAGE,spread,NULL,0,0,K_ARREARS,0,K_ARREARS,0,qCredit_Adjust20,5,2);

	}

	//FreePointerTabChar(NomIndex,size_Index);

	if (vIndex[0] == NULL) {delete[] vIndex; vIndex=NULL;}

	bool flgForcedStrike = false;

	for (i=0;i<nbelt;i++)
	{
		if ((correl_low->Elt(i)!=CREDIT_DEFAULT_VALUE) || 
			(correl_up->Elt(i)!=CREDIT_DEFAULT_VALUE))	
		{
			flgForcedStrike = true;
			break;
		}
	}

	if (flgForcedStrike) 
	{
		if (vIndex)
		{
		for (i=0;i<list1.size();i++)
		{
			if (vIndex[i])
				delete vIndex[i];
			vIndex[i] = NULL;
		}
		delete[] vIndex; vIndex=NULL;
		}
	}
} 

/**

void CvtStrPropInfo(char* input,ICM_ProportionsInfo*& info)
{
	if (strlen(input) == 0) return;

	vector<string> Currency;
	vector<string> IndexName;
	vector<string> BaseCorrelName;
	vector<string> Vtype;
	vector<string> Terms;
	vector<double> Proportions;
	vector<string> IndexDefProb;
	vector<double> StrikeLow;
	vector<double> StrikeUp;

	info = new ICM_ProportionsInfo();

	int maxsize = strlen(input);	
	int nbseparator = 0;
	int i = 0,j=0;
	int nbelt = 0;
	int il=0;

	for (i=0;i<maxsize;i++)
	{ if (input[i] == ';') nbseparator++;}

	nbelt = nbseparator + 1;

	int indice = 0;

	CCString chaine = (CCString) input;
	vector<CCString> list1;

	chaine.Parser(';',list1);	

	IndexName.resize(list1.size());
	BaseCorrelName.resize(list1.size());
	Proportions.resize(list1.size());
	StrikeLow.resize(list1.size());
	StrikeUp.resize(list1.size());
	IndexDefProb.resize(list1.size());
	Currency.resize(list1.size());
	Vtype.resize(list1.size());
	Terms.resize(list1.size());

	for (i=0;i<list1.size();i++)
	{
		vector<CCString> list2;
		list1[i].ParserWithEmpty(':',list2);	

		CCString chaine = list2[1] + (CCString) "_" + list2[0] + (CCString) "_" +  list2[3];

		IndexName[i] = (string)((const char*)chaine);
		BaseCorrelName[i] = (string)((const char*)list2[1]);
		Proportions[i] = atof((const char*)list2[5])/100.;
		Currency[i] = (string)((const char*)list2[0]);
		IndexDefProb[i] = (string)((const char*)list2[6]);
		Vtype[i] = (string)((const char*)list2[2]);
		Terms[i] = (string)((const char*)list2[3]);

		if (!(list2[7] == "") && !(list2[8] == ""))
		{
			StrikeLow[i] = atof((const char*)list2[7])/100.;
			StrikeUp[i]=atof((const char*)list2[8])/100.;
		}
		else
		{
			StrikeLow[i] = CREDIT_DEFAULT_VALUE;
			StrikeUp[i] = CREDIT_DEFAULT_VALUE;
		}

	}

	info->SetIndexName(IndexName);
	info->SetBaseCorrelName(BaseCorrelName);
	info->SetProportions(Proportions);
	info->SetStrikeLow(StrikeLow);
	info->SetStrikeUp(StrikeUp);
	info->SetIndexDefProb(IndexDefProb);
	info->SetCurrency(Currency);
	info->SetVtype(Vtype);
	info->SetTerms(Terms);

}
**/ 

bool ExistXMLNode (void* nod, string nodeName)
{
	MSXML2::IXMLDOMNode* node = (MSXML2::IXMLDOMNode*)nod;
	MSXML2::IXMLDOMNode* tmpNode = NULL;
	BSTR resultat = NULL;

	try
	{
		node->selectSingleNode((_bstr_t)nodeName.c_str(), &tmpNode);
		if (tmpNode != NULL)
			return true;
		else
			return false;
	}
	catch(...)
	{
		string msg="Pb in getting node " + nodeName;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 (char*) msg.c_str());
	}

	return false;
}

// 
//	

void DeduceRefDateAndStub(ARM_Date StartDate,
						  ARM_Date EndDate,
						  ARM_Date StartStubDate,
						  ARM_Date EndStubDate,
						  const int& payFreq,
						  const int & RefDay,
						  const std::string& ccy,
						  char* ReferenceDate,
						  int& stub)
{
	ARM_Date dateToday ; // (19,12,2005);

	if (StartStubDate != dateToday)
		StartStubDate.JulianToStrDate(ReferenceDate);
	else if (EndStubDate != dateToday)
	{
		ARM_Date tmpDate (StartDate);
		tmpDate.AddPeriodMult(payFreq,1,ccy);
		tmpDate.SetDay(EndStubDate.GetDay());
		tmpDate.JulianToStrDate(ReferenceDate);
	}
	else 
		strcpy(ReferenceDate,"NULL");

	if (EndStubDate == dateToday)
	{
		if (EndDate.GetDaysInMonth(EndDate.GetMonth(),EndDate.GetYear())>=RefDay)
			EndDate.SetDay(RefDay);
		else
			EndDate.SetDay(EndDate.GetDaysInMonth(EndDate.GetMonth(),EndDate.GetYear()));
	}

	if (strcmp((const char*) ReferenceDate, "NULL") != 0)
	{
		ARM_Date tmpDate (EndStubDate);
		tmpDate.AddPeriodMult(payFreq,1,ccy);
		if (tmpDate > EndDate)
			stub = K_SHORTEND;
		else
			stub = K_LONGEND;

		ARM_Date firstRefDate (ReferenceDate);
		ARM_Date refDate (ReferenceDate);
		if (refDate.GetDay() != RefDay)
		{
			if (refDate.GetDaysInMonth(refDate.GetMonth(),refDate.GetYear())>=RefDay)
				refDate.ChangeDate(RefDay,refDate.GetMonth(),refDate.GetYear());
			else
				refDate.ChangeDate(refDate.GetDaysInMonth(refDate.GetMonth(),refDate.GetYear()),refDate.GetMonth(),refDate.GetYear());

			if (fabs(refDate.GetJulian()-firstRefDate.GetJulian()) > 10)
			{
				if ( refDate.GetMonth() == MAXMONTH )
				{
					refDate.ChangeDate(refDate.GetDay(),1,refDate.GetYear()+1);
				}
				else
				{
					refDate.ChangeDate(refDate.GetDay(),refDate.GetMonth()+1,refDate.GetYear());
				}
			}
			refDate.JulianToStrDate(ReferenceDate);
		}
	}
}


void DeduceRefDateAndStub2(const ARM_Date& StartDate,
						  const ARM_Date& EndDate,
						  const ARM_Date*StartStubDate,
						  const ARM_Date*EndStubDate,
						  const int& payFreq,
						  const int& RefDay,
						  const ARM_Currency& ccy,
						  // output 
						  bool &isRefDate,
						  ARM_Date& refDate,
						  int& stub)
{
	/**
	isRefDate=false; 
	
	 
	if (StartStubDate )  
	{ 
		refDate=StartStubDate; isRefDate=true; 
	}
	else if (EndStubDate )
	{	
		refDate=StartDate; 
		refDate.AddPeriodMult(payFreq,1,ccy.GetCcyName());
		refDate.SetDay(EndStubDate.GetDay());
		isRefDate=true; 
	}

	//	Si pas de Stub End Date : on fait en sorte que la EndDate 
	//	tombe sur le jour de roll , ou bien sur le dernier jour du mois. 
	ARM_Date tempEndDate; 
	if (EndStubDate == 0)
	{
		if (EndDate.GetDaysInMonth(EndDate.GetMonth(),EndDate.GetYear())>=RefDay)
			EndDate.SetDay(RefDay);
		else
			EndDate.SetDay(EndDate.GetDaysInMonth(EndDate.GetMonth(),EndDate.GetYear()));
	}

	if (isRefDate) 
	{
		stub = K_LONGEND; 

		ARM_Date firstRefDate (ReferenceDate);
		ARM_Date refDate (ReferenceDate);
		if (refDate.GetDay() != RefDay)
		{
			if (refDate.GetDaysInMonth(refDate.GetMonth(),refDate.GetYear())>=RefDay)
				refDate.ChangeDate(RefDay,refDate.GetMonth(),refDate.GetYear());
			else
				refDate.ChangeDate(refDate.GetDaysInMonth(refDate.GetMonth(),refDate.GetYear()),refDate.GetMonth(),refDate.GetYear());

			if (fabs(refDate.GetJulian()-firstRefDate.GetJulian()) > 10)
			{
				if ( refDate.GetMonth() == MAXMONTH )
				{
					refDate.ChangeDate(refDate.GetDay(),1,refDate.GetYear()+1);
				}
				else
				{
					refDate.ChangeDate(refDate.GetDay(),refDate.GetMonth()+1,refDate.GetYear());
				}
			}
			refDate.JulianToStrDate(ReferenceDate);
		}
	}
	**/ 
}


CCString GetIndexNameFromARMIndex(ARM_INDEX_TYPE index)
{
	CCString INDX = "EURIB";

    switch(index)
    {
        case K_LIBOR3M :
		case K_LIBOR6M :
		case K_PIBOR3M :
        {
            INDX = "LIBOR";
        };
        break;
        case K_EURIBOR6M :
        {
            INDX = "EURIB";
        };
        break;
        default :
        {
            return(INDX);
        };
        break;
    }

    return(INDX);
}