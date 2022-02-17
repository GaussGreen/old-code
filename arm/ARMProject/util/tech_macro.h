//	NATIXIS
//
//

#ifndef TECH_MACRO_H_
#define TECH_MACRO_H_

#include <expt.h>
#include <crtdbg.h>
#include <sstream>
#include <list>
#include <map>

/////////////////////////////////////////////////
// Singleton Logger

class Logger
{ 
public:	
	// singleton
	static Logger& instance();
	
	void enable() {isEnabled_=true;}
	void disable() {isEnabled_=false;}
	void setConfig(const std::string& listenerIp, bool useBuffer, int bufferSize=0);
	void addLog(const std::string& msg);
	void sendLogs();
private:
	std::map<std::string, int> logs_;
	bool isEnabled_;
	bool useBuffer_;
	int nbLogs_;
	int bufferSize_;

	Logger();
	Logger(const Logger&) {}
public:
	virtual ~Logger() 
	{
		sendLogs();
	}	
};

#ifndef	ADD_LOG
#define ADD_LOG(msgLOG)							\
{												\
	std::stringstream _sstr;					\
	_sstr<<msgLOG<<std::endl;					\
	Logger::instance().addLog(_sstr.str());		\
} 
#endif	//LOG

#ifndef	SEND_LOGS
#define SEND_LOGS								\
{												\
	Logger::instance().sendLogs();				\
} 
#endif	//LOG

#endif	//TECH_MACRO_H_

