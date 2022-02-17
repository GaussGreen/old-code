#include "esl_log.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdarg.h>

#include <DRDiagnosticMsgCB.h> // DR_Error hookup


#include <stdexcept>
#include <iostream>

Log eslLog;

/* ****************************** Stringifier *********************************/

std::string Stringifier::get(const std::ostream& os) {
	const Stringifier::Stream *s = dynamic_cast<const Stringifier::Stream*>(&os);
	if (s) return s->str();
	return "<Problem with Stringifier>";
}

/* ********************************** Log *************************************/

char const* Log::levelStrings[] =
{ " EMERG",
  " ALERT",
  "  CRIT",
  " ERROR",
  "  WARN",
  "NOTICE",
  "  INFO",
  " DEBUG",
  "  NONE"
};

//-----------------------------------------------------------------------------
// Log system
//-----------------------------------------------------------------------------

extern "C" {

extern DRDiagnosticMsgCB EslLogCallbackFnp;

#ifdef WIN32
#define vsnprintf _vsnprintf
#endif

static int defaultStreamCallback(char const* fmt, va_list args)
{
    if (eslLog.active())
    {
        char buf[1024];
        int cnt = vsnprintf(buf, sizeof(buf), fmt, args);
        if (cnt > 0) eslLog.write(Log::error, buf, cnt);
    }
    return 0;
}

}


Log::Log() {
    EslLogCallbackFnp = defaultStreamCallback;
	for (unsigned i=0; i<=size; ++i) devices[i]=0;
}

Log::~Log() {
    EslLogCallbackFnp = 0;
}

void Log::write(Log::Level l, const char* msg, size_t msgLen) {
	Device **cur = devices;
	while (*cur) (*cur++)->write(l, msg, msgLen);
}

void Log::insert(Log::Device *l) {
	for (unsigned i=0; i<size; ++i) {
		if (devices[i]==0) {
			devices[i]=l;
			return;
		}
	}
	throw EslEx("Cannot insert device into Log: too many devices");
}

void Log::remove(Log::Device *l) {
	bool found=false;
	for (unsigned i=0; i<size; ++i) {
		if (devices[i]==l) { found=true; }
		if (found) { devices[i]=devices[i+1];}
	}
}

void Log::clear(void) {
	for (unsigned i=0; i<size; ++i) devices[i]=0;
}

//-----------------------------------------------------------------------------
// Log message
//-----------------------------------------------------------------------------

/* ******************************** FileLogger ********************************/

void FileLogger::openFile(bool append) {
	m_fp = ::fopen(m_fName.c_str(), (append ? "ab" : "wb") );
	if (!m_fp) {
        std::cerr<<"Cannot create log file \""<<m_fName<<"\""<<std::endl;
		return;
	}
	if (!append) { m_currentSize=0; return; }
	
//	struct stat st;
	m_currentSize=0; //!!!!
/*	if (::fstat(::fileno(m_fp), &st) == 0) {
		m_currentSize = st.st_size;
	} else {
		cerr<<"Cannot access to log file \""<<m_fName<<"\""<<endl;
		fclose(m_fp); 
		m_fp=0;
		return;
	}*/
}

void FileLogger::closeFile(void) {
	if (!m_fp) return;
	if (::fclose(m_fp)!=0) {
        std::cerr<<"Cannot close log file \""<<m_fName<<"\""<<std::endl;
	}
	m_fp = 0;
}

FileLogger::FileLogger(Log& log, char const* fName, bool append, size_t maxSize)
: Device(log), m_fName(fName), m_maxSize(maxSize), m_currentSize(0) {
	openFile(append);
}

void FileLogger::write(Log::Level l, const char* buf, size_t len) {
	if (::fwrite(buf, len, 1, m_fp)!=1 
	|| ::fwrite("\n", 1, 1, m_fp)!=1 
	|| ::fflush(m_fp)!=0) {
        std::cerr<<"Cannot write to log file"<<std::endl;
	}	
	m_currentSize += len;
	if (m_maxSize==0 || m_currentSize < m_maxSize) return;

	closeFile();
    std::string oldFName = m_fName + ".bak";
	::rename(m_fName.c_str(), oldFName.c_str());
	openFile(true);
}


