#include <winsock2.h>
//#include <windows.h>	// for ::OutputDebugString
//#include <winbase.h>	// for ::GetModuleHandleEx
#include <lmcons.h.>
//#include <memory>
//#include <iostream>
//#include <fstream>

#include "tech_macro.h"


using namespace std;


//	-------------------------------------------------------
//
//		This is singleton socket , allow sending to 
//		localhost & remote host
// 
class _Socket 
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
	void remotesend(const std::string&msg )
	{
		if (remoteWait) { remoteWait--; return; }
		if (!send(remoteSock,remoteHost,remotePort,msg)) 
				remoteWait=0; 
	}
public:
	static _Socket& get()
	{
		if (itsInstance.get()==0) 
			itsInstance=std::auto_ptr<_Socket>(new _Socket); 
		return *itsInstance.get(); 
	}
	void setRemoteHost(std::string ip) { remoteHost= ip;}
private:
	_Socket() : remoteSock(INVALID_SOCKET) 
	{ 
		WSADATA init ; 
		WSAStartup(MAKEWORD(2,2),&init); 
		remoteSock=socket(AF_INET,SOCK_STREAM,0); 
		remoteWait=0; 
		remotePort=2904;	//+1
	} 
	~_Socket() 
	{
		if (remoteSock!=INVALID_SOCKET) { shutdown(remoteSock,2); closesocket(remoteSock);}
		WSACleanup(); 
	}
private: 
	WSADATA init; 
	SOCKET remoteSock; 
	std::string remoteHost; 
	long remoteWait; 
	long remotePort; 
private:
	static std::auto_ptr<_Socket> itsInstance; 
	friend class std::auto_ptr<_Socket> ;
} ; 
std::auto_ptr<_Socket> _Socket::itsInstance; 


///////////////////////////////////////////////////////////////////////////////////
// Impl Logger

Logger::Logger()
{
	isEnabled_ = false;
	useBuffer_ = true;
	nbLogs_ = 0;
	bufferSize_= 2;
}


void Logger::setConfig(const std::string& listenerIp, bool useBuffer, int bufferSize)
{
	useBuffer_ = useBuffer;
	nbLogs_ = 0;
	bufferSize_ = bufferSize;
	_Socket::get().setRemoteHost(listenerIp);
}


Logger& Logger::instance()
{
	// Meyers singleton	(see "Modern C++ Design", Andrei Alexandrescu)
	static Logger instance_;
	return instance_;
}

void Logger::addLog(const std::string& msg)
{
	if(isEnabled_)
	{
		map<string, int>::iterator logElem = logs_.find(msg);
		if(logElem != logs_.end())	{
			logElem->second += 1;
		}else{
			logs_[msg] = 1;
		}
		nbLogs_ += 1;
		if(nbLogs_ >= bufferSize_ && useBuffer_)	//ajouter compteur
			sendLogs();
	}
}

void Logger::sendLogs()
{
	if(isEnabled_)
	{
		unsigned long l=UNLEN ; 
		std::string	user ; 
		user.resize(UNLEN +1 ); 
		::GetUserName(&(*user.begin()) , &l ) ;
		
		l = MAX_COMPUTERNAME_LENGTH  ;
		std::string computer ;
		computer.resize(MAX_COMPUTERNAME_LENGTH +1); 
		::GetComputerName(&(*computer.begin()) , &l  ) ;

		long pid= GetCurrentProcessId(); 

		map<string, int>::const_iterator mapIter;
		std::stringstream ssSend;
		for(mapIter = logs_.begin(); mapIter != logs_.end(); ++mapIter)
		{
			ssSend << mapIter->first << "|" << mapIter->second << "|";
		}
		ssSend << ":arm:" << user << ":" << computer << ":" << pid ; 

		_Socket::get().remotesend(ssSend.str()); 
		// clear the logs
		logs_.clear();
		nbLogs_ = 0;
	}
}

