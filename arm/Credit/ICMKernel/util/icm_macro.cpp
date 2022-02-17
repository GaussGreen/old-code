//
//
//

#include <winsock2.h>
#include <windows.h>	// for ::OutputDebugString
#include <winbase.h>	// for ::GetModuleHandleEx
#include <lmcons.h.>
#include <memory>
#include "ICMKernel\util\icm_macro.h"

#define stringer( x ) "" #x ""


#if _MSC_VER >= 1300    // for VC 7.0
  // from ATL 7.0 sources
  #ifndef _delayimp_h
  extern "C" IMAGE_DOS_HEADER __ImageBase;
  #endif
#endif

HMODULE GetCurrentModule()
{
#if _MSC_VER < 1300    // earlier than .NET compiler (VC 6.0)

  // Here's a trick that will get you the handle of the module
  // you're running in without any a-priori knowledge:
  // http://www.dotnet247.com/247reference/msgs/13/65259.aspx

  MEMORY_BASIC_INFORMATION mbi;
  static int dummy;
  VirtualQuery( &dummy, &mbi, sizeof(mbi) );

  return reinterpret_cast<HMODULE>(mbi.AllocationBase);

#else    // VC 7.0

  // from ATL 7.0 sources

  return reinterpret_cast<HMODULE>(&__ImageBase);
#endif
}



//	-------------------------------------------------------
//
//		This is singleton socket , allow sending to 
//		localhost & remote host
// 
class ICM_Socket 
{
public:
private:
	bool send(SOCKET&s,const std::string&addr,int port,const std::string& msg)
	{
		if (s==INVALID_SOCKET) s = socket(AF_INET,SOCK_STREAM,0); 
		if (s!=INVALID_SOCKET) 
		{
			// int tempo=1; 
			// setsockopt(s,IPPROTO_TCP,TCP_NODELAY,(char *)&tempo,sizeof(tempo));
			SOCKADDR_IN addrIn ;
			addrIn.sin_family=AF_INET;
			addrIn.sin_addr.s_addr=inet_addr(addr.c_str());   
			addrIn.sin_port=htons(port) ;
			int err = connect(s,(struct sockaddr*)&addrIn,sizeof(addrIn)); 
			if ( (err==0) || (WSAGetLastError()==WSAEISCONN) )
				err=::send(s,msg.c_str(),msg.size()+1,0) ; 
			if (err==SOCKET_ERROR) 
			{
				int i=WSAGetLastError(); 
				shutdown(s,2); 
				closesocket(s); 
				s=INVALID_SOCKET ;
				return false; 
			}
			return true; 
		}
		return false; 
	}
public:
	void localsend(const std::string&msg ) 
	{
		if (localWait) { localWait--; return; }
		if (!send(localSock,localHost,localPort,msg)) 
				localWait=100; 
	}
	void remotesend(const std::string&msg )
	{
		if (remoteWait) { remoteWait--; return; }
		if (!send(remoteSock,remoteHost,remotePort,msg)) 
				remoteWait=10; 
	}
public:
	static ICM_Socket& get()
	{
		if (itsInstance.get()==0) itsInstance=std::auto_ptr<ICM_Socket>(new ICM_Socket); 
		return *itsInstance.get(); 
	}
private:
	ICM_Socket() : localSock(INVALID_SOCKET), remoteSock(INVALID_SOCKET) 
	{ 
		WSADATA init ; 
		WSAStartup(MAKEWORD(2,2),&init); 
		localSock=socket(AF_INET,SOCK_STREAM,0); 
		remoteSock=socket(AF_INET,SOCK_STREAM,0); 
		localHost = "127.0.0.1" ;
		remoteHost= "10.14.21.34" ;
		localWait=0; 
		remoteWait=0; 
		localPort=1972; 
		remotePort=2903; 
	} 
	~ICM_Socket() 
	{
		if (remoteSock!=INVALID_SOCKET) { shutdown(remoteSock,2); closesocket(remoteSock);}
		if (localSock!=INVALID_SOCKET)  { shutdown(localSock,2); closesocket(localSock); }
		WSACleanup(); 
	}
private: 
	WSADATA init; 
	SOCKET localSock,remoteSock; 
	std::string localHost,remoteHost; 
	long localWait,remoteWait; 
	long localPort,remotePort; 
private:
	static std::auto_ptr<ICM_Socket> itsInstance; 
	friend class std::auto_ptr<ICM_Socket> ;
} ; 
std::auto_ptr<ICM_Socket> ICM_Socket::itsInstance; 
//	---------------------------------------------------------------------------
void 
IcmLogger::logDebugger(qIcmLoggerLevel level,const std::string& file,long line,const std::string&msg)
{
	std::stringstream sstr1; sstr1<<file<<" ("<<line<<") :"<<msg ; 
	std::string msg1 = sstr1.str() ;
	::OutputDebugString(msg1.c_str()); 

	unsigned long l=UNLEN ; 
	std::string user ; 
	user.resize(UNLEN +1 ); 
	::GetUserName(&(*user.begin()) , &l ) ;
	l =MAX_COMPUTERNAME_LENGTH  ;
	std::string computer ;
	computer.resize(MAX_COMPUTERNAME_LENGTH +1); 
	::GetComputerName(  &(*computer.begin()) , &l  ) ;
	long pid= GetCurrentProcessId(); 
	std::stringstream sstr2 ; 
	sstr2<<msg1<<":icm:"<<level<<":" <<user<<":"<<computer<<":"<<pid ;
	ICM_Socket::get().localsend(sstr2.str()); 
	if (level!=SILENT) 
		ICM_Socket::get().remotesend(sstr2.str()); 
}

//	---------------------------------------------------------------------------
void 
IcmLogger::logStream(std::ostream&o,const std::string&  file,long line,const std::string& msg)
{
	o<<file<<" ("<<line<<") :"<<msg ;
	std::stringstream sstr; sstr<<file<<" ("<<line<<") :"<<msg; 
	::OutputDebugString(sstr.str().c_str()); 
}
//	---------------------------------------------------------------------------
class nullbuf : public std::streambuf
{
public:
int sync() { return 0; }
int overflow( int ch ) { return 0; }
};
//	-- global 
nullbuf nbuf; 
std::ostream nout(&nbuf); 

//	---------------------------------------------------------------------------
class PreciseChronoImpl
{
protected:
	LARGE_INTEGER m_liStart;
	LARGE_INTEGER m_liStop;
	LONGLONG m_llFrequency;
public:
	PreciseChronoImpl(void)
	{
		LARGE_INTEGER liFrequency;
		QueryPerformanceFrequency(&liFrequency);
		m_llFrequency = liFrequency.QuadPart;
	}
	void Start(void)
	{
		Sleep(0);
		QueryPerformanceCounter(&m_liStart);
	}
	void Stop(void) 
	{
		QueryPerformanceCounter(&m_liStop);
	}
	double GetDurationInSeconds(void) const { return (m_liStop.QuadPart-m_liStart.QuadPart)/ m_llFrequency;}
	double GetDurationInTicks(void) const { return (m_liStop.QuadPart-m_liStart.QuadPart) ; }
	double GetDurationInMilliseconds(void) const { return (m_liStop.QuadPart-m_liStart.QuadPart)*1000.0 / m_llFrequency;}
	double GetDurationInMicroseconds(void) const { return (m_liStop.QuadPart-m_liStart.QuadPart)*1000000.0 / m_llFrequency;}
};

//	---------------------------------------------------------------------------
PreciseChrono::PreciseChrono(void) { impl=new PreciseChronoImpl(); }
PreciseChrono::~PreciseChrono(void) { if(impl) delete impl; impl=0 ;} 
void PreciseChrono::Start(void) {impl->Start(); }
void PreciseChrono::Stop(void) {impl->Stop(); }
double PreciseChrono::GetDurationInSeconds(void) const { return impl->GetDurationInSeconds(); }
double PreciseChrono::GetDurationInTicks(void) const { return impl->GetDurationInTicks(); }
double PreciseChrono::GetDurationInMilliseconds(void) const { return impl->GetDurationInMilliseconds(); }
double PreciseChrono::GetDurationInMicroseconds(void) const { return impl->GetDurationInMicroseconds(); }
//	---------------------------------------------------------------------------
QuantifyerInfo::QuantifyerInfo(const std::string&arg,const std::string& _file,long _line) 
	: file(_file),line(_line),callCount(0),avgTime(0),maxTime(0),minTime(0),elapsedTime(0),msg(arg) {} 
QuantifyerInfo::~QuantifyerInfo() 
{
	std::stringstream sstr; 
	sstr<<"Quantifyer:"<<msg<<":elapsedTime:"<<elapsedTime<<":call#:"<<callCount<<":avgTime:"<<avgTime<<":maxTime#:"<<maxTime<<":minTime#:"<<minTime<<std::endl;
	IcmLogger::logDebugger(SILENT,file,line,sstr.str()); 
}
//	---------------------------------------------------------------------------
Quantifyer::Quantifyer(QuantifyerInfo&	ref) : quantifyerInfo(ref) { chrono.Start(); } 
Quantifyer::~Quantifyer() 
{
	chrono.Stop() ; 
	quantifyerInfo.callCount++; 
	double duration = chrono.GetDurationInMilliseconds(); 
	if (quantifyerInfo.callCount==1) 
	{ 
		quantifyerInfo.maxTime=quantifyerInfo.minTime=quantifyerInfo.avgTime=duration; 
	}
	else 
	{
		if (duration>quantifyerInfo.maxTime) quantifyerInfo.maxTime=duration; 
		if (duration<quantifyerInfo.minTime) quantifyerInfo.minTime=duration; 
		quantifyerInfo.avgTime = (quantifyerInfo.callCount-1.)/quantifyerInfo.callCount*quantifyerInfo.avgTime
			+ duration/quantifyerInfo.callCount ;
	}
	quantifyerInfo.elapsedTime+=duration ;
} 
//	---------------------------------------------------------------------------
static VersionInfoSingleton versionInfoSingleton; 
//	---------------------------------------------------------------------------
std::string 
VersionInfoSingleton::version() 
{
	return "2.2.3" ;
}
//	---------------------------------------------------------------------------
std::string 
VersionInfoSingleton::longversion() 
{
	std::string ver("ICMKernel-"); 
	ver += version() ;
#ifndef icmBUILD
	ver+="[build:free]";
#else
	ver+="[build:" ;
	// char s[10]; sprintf(s,"%s",\" icmUSERNAME \" ); 
	ver+=  icmUSERNAME  ;
	ver+="@"; 
	ver+= icmCOMPUTERNAME  ;
	ver+="]" ;
#endif
#ifdef	_DEBUG
	ver +="[Debug]"; 
#else
	ver +="[Release]"; 
#endif	
	ver +=" Compiled " ; 
	ver +=  __DATE__ ; 
	ver +=" " ; 
	ver +=  __TIME__  ; 
	return ver ;
}
//	---------------------------------------------------------------------------
VersionInfoSingleton::VersionInfoSingleton()
{
	std::string file; 
	file.resize(_MAX_PATH +1); 
	::GetModuleFileName( GetCurrentModule() , &(*file.begin()) , file.size() ) ;
	ICMMSG(INFO,"Running "<<longversion()<<", bin="<<file); 
}
//	---------------------------------------------------------------------------
VersionInfoSingleton::~VersionInfoSingleton()
{
	// ICMMSG(INFO,"Exiting "<<longversion()); 
}

