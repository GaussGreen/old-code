/*-----------------------------------------------------------------------------*

     Global E-ToolKit header

 *-----------------------------------------------------------------------------*/
#ifndef ARM_LOCAL_ETOOLKITX_H
#define ARM_LOCAL_ETOOLKITX_H

#include <libCCTools++/CCString.h>
#include <string>

using namespace std;
struct IeToolkit;


class ARM_Etoolkit;

// TMP const char* strnothing = "";

extern ARM_Etoolkit* eToolkitPtr;


#ifndef DEFAULT_USER
	
    #define DEFAULT_USER	"NeSertARien"
	#define DEFAULT_PWD		"NeSertARien"

	#define SUMMIT_PROD_CONNEXION_CONTEXT		"msbdb_sumotc:MX_SV_TPSMTOTC"
	#define SUMMIT_PROD_CONNEXION_USERNAME		DEFAULT_USER
	#define SUMMIT_PROD_CONNEXION_PASSWD		DEFAULT_PWD
	#define SUMMIT_PROD_IT_CONFIG_DOMAINSDIR	"IT_CONFIG_DOMAINS_DIR=\\\\cm.net\\ShareParis\\Soft\\Orbix\\orbix621\\config\\Prod\\MX_DO_FRMSCUPIDON62\\etc\\domains"
	#define SUMMIT_PROD_IT_DOMAIN_NAME			"IT_DOMAIN_NAME=MX_DO_FRMSCUPIDON62"

	#define SUMMIT_REPLI_CONNEXION_CONTEXT		"msbdb_sumotc:MX_SV_SESMTOTC"
	#define SUMMIT_REPLI_CONNEXION_USERNAME		DEFAULT_USER
	#define SUMMIT_REPLI_CONNEXION_PASSWD		DEFAULT_PWD
	#define SUMMIT_REPLI_IT_CONFIG_DOMAINSDIR	"IT_CONFIG_DOMAINS_DIR=\\\\cm.net\\ShareParis\\Soft\\Orbix\\orbix621\\config\\Prod\\MX_DO_FRMSCUPIDON62\\etc\\domains"
	#define SUMMIT_REPLI_IT_DOMAIN_NAME			"IT_DOMAIN_NAME=MX_DO_FRMSCUPIDON62"

	#define SUMMIT_REC342_CONNEXION_CONTEXT		"msudb_summit_v5:MU_SV_SUMFR"
	#define SUMMIT_REC342_CONNEXION_USERNAME	DEFAULT_USER
	#define SUMMIT_REC342_CONNEXION_PASSWD		DEFAULT_PWD
	#define SUMMIT_REC342_IT_CONFIG_DOMAINSDIR	"IT_CONFIG_DOMAINS_DIR=\\\\cm.net\\ShareParis\\Soft\\Orbix\\orbix621\\config\\Rec\\MU_DO_FRMSLAMENTIN62\\etc\\domains"
	#define SUMMIT_REC342_IT_DOMAIN_NAME		"IT_DOMAIN_NAME=MU_DO_FRMSLAMENTIN62"

	#define SUMMIT_INFOC_CONNEXION_CONTEXT		"msbdb_sumotc:MX_SV_ICSMTOTC"
	#define SUMMIT_INFOC_CONNEXION_USERNAME		DEFAULT_USER
	#define SUMMIT_INFOC_CONNEXION_PASSWD		DEFAULT_PWD
	#define SUMMIT_INFOC_IT_CONFIG_DOMAINSDIR	"IT_CONFIG_DOMAINS_DIR=\\\\cm.net\\ShareParis\\Soft\\Orbix\\orbix621\\config\\Prod\\MX_DO_FRMSCUPIDON62\\etc\\domains"
	#define SUMMIT_INFOC_IT_DOMAIN_NAME			"IT_DOMAIN_NAME=MX_DO_FRMSCUPIDON62"

	#define ETK_MKTDATA_RETRIEVAL_IS_ON  1
	#define ETK_MKTDATA_RETRIEVAL_IS_OFF 0


#endif



class ARM_Etoolkit
{
public:

	ARM_Etoolkit(void)
	{
		itsEToolkitptr = NULL;
		itsUserName = "";
		itsPassword = "";
		itsDataBaseContext = "";
		itsItConfigPath = "";
		itsItDaemonPort = "";
		itsItNamesServerHost = "";
		itsItConfigDomainsDir = "";
		itsItDomainName = "";

		Init();
	}
	
   ~ARM_Etoolkit()
	{
		Deconnect();

/*		if (itsUserName)
			delete itsUserName;
		itsUserName = NULL;

		if (itsPassword)
			delete itsPassword;
		itsPassword = NULL;

		if (itsDataBaseContext)
			delete itsDataBaseContext;
		itsDataBaseContext = NULL;

		if (itsItConfigPath)
			delete itsItConfigPath;
		itsItConfigPath = NULL;

		if (itsItDaemonPort)
			delete itsItDaemonPort;
		itsItDaemonPort = NULL;

		if (itsItNamesServerHost)
			delete itsItNamesServerHost;
		itsItNamesServerHost = NULL;

		if (itsItConfigDomainsDir)
			delete itsItConfigDomainsDir;
		itsItConfigDomainsDir = NULL;

		if (itsItDomainName)
			delete itsItDomainName;
		itsItDomainName = NULL;
*/	}

	void Init(const CCString& userName = ( (const CCString&)SUMMIT_PROD_CONNEXION_USERNAME),
			  const CCString& passWord = ( (const CCString&)SUMMIT_PROD_CONNEXION_PASSWD),
			  const CCString& databasecontext = ( (const CCString&)SUMMIT_PROD_CONNEXION_CONTEXT),
			  const CCString& itConfigDomainsDir = ( (const CCString&)SUMMIT_PROD_IT_CONFIG_DOMAINSDIR),
			  const CCString& itDomainName = ( (const CCString&)SUMMIT_PROD_IT_DOMAIN_NAME)
			  )
	{
		itsEToolkitptr = NULL;

/*		if (itsUserName)
		   delete itsUserName;
*/		itsUserName = (const char *) userName;

/*		if (itsPassword)
			delete itsPassword;
*/		itsPassword = (const char *) passWord;

/*		if (itsDataBaseContext)
			delete itsDataBaseContext;
*/		itsDataBaseContext = (const char *) databasecontext;

/*		if (itsItConfigPath)
			delete itsItConfigPath;
*/		itsItConfigPath = "";

/*		if (itsItDaemonPort)
			delete itsItDaemonPort;
*/		itsItDaemonPort = "";

/*		if (itsItNamesServerHost)
			delete itsItNamesServerHost;
*/		itsItNamesServerHost = "";

/*		if (itsItConfigDomainsDir)
			delete itsItConfigDomainsDir;
*/		itsItConfigDomainsDir = (const char *) itConfigDomainsDir;

/*		if (itsItDomainName)
			delete itsItDomainName;
*/		itsItDomainName = (const char *) itDomainName;
	}

	void Connect(const CCString& userName,
				 const CCString& passWord,
				 const CCString& databasecontext,
				 const CCString& itConfigDomainsDir,
				 const CCString& itDomainName)
	{
		Init(userName, passWord, databasecontext, itConfigDomainsDir, itDomainName);

		Connect();
	}


	void Connect();

	int isConnecte()
	{
		if ( itsEToolkitptr != NULL )
		   return 1;
		else
		   return 0;
	}

	void Deconnect();

	void Shutdown();

	void Execute(CCString command, CCString xmlRequest, CCString & xmlResponse, CCString & messageList);

	void GetCommsetName(CCString varName, CCString varName2, CCString varAsOf, CCString varType, CCString varCvName, CCString& response, CCString& messageList);

	void GetRefRate(CCString source, CCString ccy, CCString index, CCString tenor, CCString& response, CCString& messageList);

private:

	IeToolkit* itsEToolkitptr;

	string itsUserName;
	string itsPassword;
	string itsDataBaseContext;
	string itsItConfigPath;
	string itsItDaemonPort;
	string itsItNamesServerHost;

	string itsItConfigDomainsDir;
	string itsItDomainName;
};


#endif