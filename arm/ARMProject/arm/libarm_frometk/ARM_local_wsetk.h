#ifndef ARM_LOCAL_WSETK_H
#define ARM_LOCAL_WSETK_H

#pragma warning(disable : 4002)

#ifdef _DEBUG
// probleme dans atlbase.h
#define va_alist ""
#endif

#include <atlbase.h>

#include <libCCTools++/CCString.h>

#include <string>
// using namespace std;

// #import <msxml3.dll> raw_interfaces_only
// using namespace MSXML2;


// #import "C:\Program Files\Common Files\MSSoap\Binaries\mssoap30.dll" named_guids \
//   exclude("IStream", "IErrorInfo", "ISequentialStream", "_LARGE_INTEGER", \
//   "_ULARGE_INTEGER", "tagSTATSTG", "_FILETIME")

// using namespace MSSOAPLib30;

namespace MSSOAPLib30
{
	struct ISoapClient; 
}; 
class ARM_EtkSoapClient;

extern ARM_EtkSoapClient* eToolkitSoapPtr;

class ARM_EtkSoapClient
{
private:
	CComPtr<MSSOAPLib30::ISoapClient>	itsSoapClient;
	const TCHAR*			itsURL;
	const TCHAR*			itsBase;
	TCHAR					itsUserName[50];

public:
	ARM_EtkSoapClient();
	ARM_EtkSoapClient(const CCString& pProd, const CCString& pBase, const CCString& pUserName);
	~ARM_EtkSoapClient();

	void Execute(CCString command, CCString xmlRequest, CCString & xmlResponse, CCString & messageList);
	std::vector<std::string> GetTradeList(CCString xmlRequest);
	void GetCommsetName(CCString name, CCString name2, CCString asOf, CCString type, CCString cvname, CCString & xmlResponse, CCString & messageList);
	void GetRefRate(CCString source, CCString ccy, CCString index, CCString tenor, CCString & xmlResponse, CCString & messageList);

};


#endif // ARM_LOCAL_WSETK_H