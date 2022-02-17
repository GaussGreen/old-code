#ifndef ESL_LOG_DOT_H
#define ESL_LOG_DOT_H

#include <string>
#include <sstream>
#include <exception>
#include <iostream>

//--------------------------------- Stringifier --------------------------------

class Stringifier {
public:
	struct Stream : public std::stringstream { ///< tmp stream to build the msg
		std::ostream& os(void) { return *this; } ///< force type conversion
	};
	static std::string get(const std::ostream& os);
};

/// St2Ststands for STream TO STring
#define St2St(x) Stringifier::get(Stringifier::Stream().os() << x)

//-------------------------------- ESL Exception -------------------------------

class EslEx : public std::exception { ///< ESL Exception
	std::string msg;
public:
	virtual const char *what(void) const throw () { return msg.c_str(); }
	EslEx(void) {}
	EslEx(const char * msg) : msg(msg) {}
	EslEx(const std::string& msg) : msg(msg) {}
	EslEx(const std::string& msg, const std::exception& e) { 
		this->msg = msg + "\n" + e.what();
	}
	virtual ~EslEx(void) throw () {}
};

#define EslExS(stream) EslEx(St2St(stream))


//-----------------------------------------------------------------------------
// ESL Exception
//-----------------------------------------------------------------------------

/** Simple exception class that uses output stream semantics to
 *  build the message. Currently does not support user defined
 *  manipulators with non-default signature.
 */
class EslException : public std::exception
{
public:
	EslException(){}
	EslException(char const* msg) { m_msg << msg; }
	EslException(std::string const& msg) { m_msg << msg; }
	EslException(std::string const& msg, const std::exception& e) {
        m_msg << msg << "\n" << e.what();
    }
    EslException(const EslException& e) : std::exception(), m_msg(e.m_msg.str()) {}

    const EslException& operator=(const EslException& e) {
        if (this != &e)
            m_msg.str(e.m_msg.str());
        return *this;
    }

    template <typename T>
    EslException& operator<< (T const& t) {
        m_msg << t;
        return *this;
    }

    EslException& operator<< (std::ostream& (*fp)(std::ostream& s)) {
        fp(m_msg);
        return *this;
    }
    EslException& operator<< (std::ios& (*fp)(std::ios& s)) {
        fp(m_msg);
        return *this;
    }
    EslException& operator<< (std::ios_base& (*fp)(std::ios_base& s)) {
        fp(m_msg);
        return *this;
    }

	virtual const char* what() const throw () { 
        EslException& me = const_cast<EslException&>(*this);
        me.m_what = m_msg.str(); 
        return m_what.c_str(); 
    }

	virtual ~EslException() throw () {}

private:
    std::ostringstream m_msg;
    std::string m_what;
};
 



//-----------------------------------------------------------------------------
// Log system
//-----------------------------------------------------------------------------

/// Message logging stream - a collection of logging devices
class Log 
{
public:

    /// Log levels based on unix syslog
    enum Level
    { 
		emerg=0, alert, crit, error, warn, notice, info, debug, none 
	};

	/// A simple message output device - pure interface plus registration
	struct Device
    {
        /** Construct the device and attach it to the given stream.
         */
        Device(Log& log) : m_log(log) {
            log.insert(this);
        }

        /** Remove the device from the stream before destroying it
         */
        virtual ~Device() {
            m_log.remove(this);
        }
        
        /** Write the the message to the output. The level argument is 
         *  provided for possible additional filtering/formating at the
         *  device level.
         */
		virtual void write(Level l, const char* msg, size_t msgLen) = 0;

    private:

        Log& m_log;
	};

    /** Log message class that uses output stream semantics to
     *  build the message. Currently does not support user defined
     *  manipulators with non-default signature.
     */
    class Msg
    {
    public:
	    Msg(Level level, Log& log) : m_level(level), m_log(log) {}

        template <typename T>
        Msg& operator<< (T const& t) {
            if (m_level <= m_log.getLevel())
                m_msg << t;
            return *this;
        }

        Msg& operator<< (std::ostream& (*fp)(std::ostream& s)) {
            if (m_level <= m_log.getLevel())
                fp(m_msg);
            return *this;
        }

        Msg& operator<< (std::ios& (*fp)(std::ios& s)) {
            if (m_level <= m_log.getLevel())
                fp(m_msg);
            return *this;
        }

        Msg& operator<< (std::ios_base& (*fp)(std::ios_base& s)) {
            if (m_level <= m_log.getLevel())
                fp(m_msg);
            return *this;
        }

        /// Output accumulated message in dtor
	    virtual ~Msg() {
            std::string s = m_msg.str();
		    m_log.write(m_level, s.c_str(), s.size()); 
        }

    private:

        Level               m_level;
        Log&                m_log;
        std::ostringstream  m_msg;
    };


    /// Create empy logger object
	Log();

    /// Delete logger object
	virtual ~Log();

	/** Send the message to all registered devices 
     *  TODO - make it nothrow
     */
	void write(Level l, const char* msg, size_t msgLen);

    /** Return true if the log stream has devices attached (can actually print)
     */
    bool active() const {
        return devices[0] != 0;
    }

    /** Get streams current log level
     */
    Level getLevel() const {
        return threshold;
    }

    /** Set stream's log level
     */
    void setLevel(Level l) {
        threshold = l;
    }

    /** Add logging device into the list
     *  @param l        logging device
     *  NOTE - should be private but breaks VCC
     */
	void insert(Device *l);

    /** Remove logging device from the list
     *  @param l        logging device
     *  NOTE - should be private but breaks VCC
     */
	void remove(Device *l);

    /** Remove all logging devices from the list */
	void clear(void);

    static char const* levelStrings[];  ///< human readable strings for log levels

private:

    enum {size=5};                      ///< VCC bug

    Level           threshold;          ///< do not display messages above threshold
    Device*         devices[size+1];    ///< list of devices to send msg to
};


extern Log eslLog; ///< define one Log as a global variable


/// Generic message logging macro
#define LOG_MSG(level, log) \
    if (level > log.getLevel()) ; \
    else Log::Msg(level, log) << Log::levelStrings[level] << " - "

#define LOG_RAW_MSG(level, log) \
    if (level > log.getLevel()) ; \
    else Log::Msg(level, log)

/// Default logger message logging macro
#define ESL_LOG(level)  LOG_MSG(level, eslLog)
#define ESL_OUT(level)  LOG_RAW_MSG(level, eslLog)


/*
#define ESL_LOG(level, stream) \
	if (Log::level <= eslLog.threshold) eslLog.write( \
	St2St(Log::levelStrings[Log::level] << " - " << stream << "\n"));

#define GEN_LOG(logObject, level, stream) \
	if (Log::level <= eslLog.threshold) (logObject).write( \
	St2St(Log::levelStrings[Log::level] << " - " << stream << "\n"));
*/

//-----------------------------------------------------------------------------
// Concrete log devices
//-----------------------------------------------------------------------------

/// FileLogger writes log messages to the named file.
class FileLogger : public Log::Device {
public:
    /** @param name     name of the log file
        @param append   if false new log is started
        */
	FileLogger(Log& log, char const* fname, bool append=true, size_t max_size=0);
    virtual ~FileLogger(void) { closeFile(); }

    virtual void write(Log::Level l, const char* buff, size_t len);

private:    	
    FILE*        m_fp;
    std::string  m_fName;
    size_t       m_maxSize;
    size_t       m_currentSize;
    
	void openFile(bool append);
	void closeFile(void);
};


/// StreamLogger writes log messages to the output stream
class StreamLogger : public Log::Device { 
public:
    /** @param os   output stream */
	StreamLogger(Log& log, std::ostream& os) : Device(log), m_os(os) {} 

	/* inline for performance */
    virtual void  write(Log::Level, const char* buf, size_t len) { 
        m_os.write(buf, len); 
    }

private:
    std::ostream& m_os;
};

/// StringLogger accumulates log messages into a string
class StringLogger : public Log::Device { 
public:
    /** @param os   output stream */
	StringLogger(Log& log) : Device(log) {} 

    virtual void  write(Log::Level l, const char* buf, size_t len) { 
        m_messages.append(buf, len); 
    }

    /// get accumulated logs 
    char const* str() {
        return m_messages.c_str();
    }

	std::string pop(void) {
		std::string copy = m_messages;
    	m_messages.erase();
        return copy;
    }

private:
    std::string m_messages;
};

#endif
