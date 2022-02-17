///////////////////////////////////////////////////////////////////////////////
// Name:        ito33/log.h
// Purpose:     logging functions and macros
// Author:      Vadim Zeitlin
// Created:     2004-10-22
// RCS-ID:      $Id: log.h,v 1.17 2006/08/11 22:51:48 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
///////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/log.h
    @brief Declares logging functions and macros.

    These macros just hide the real logging framework behind them. Currently we
    use rLog but it could change later.
 */

#ifndef _ITO33_LOG_H_
#define _ITO33_LOG_H_

#ifndef ITO33_NO_LOG

#include "ito33/common.h"

#ifdef _WIN32
  // include windows.h before rlog does it unsafely:
  #include "ito33/win32/winwrap.h"
#endif // _WIN32

#include <rlog/rlog.h>
#include <rlog/RLogChannel.h>
#include <rlog/StdioNode.h>

#include "ito33/string.h"
#include "ito33/debug.h"

namespace ito33
{

/**
    Contains classes and functions for message logging.
 */
namespace log
{

/**
    Contains information about the location of the log message in the source
    code.

    This is a struct with at least the following fields
      - fileName
      - lineNum
      - functionName (may be empty with older compilers)
 */
typedef ::rlog::PublishLoc Location;

/**
    Log messages level or severity.

    This enum has at least the following elements:
      - Log_Critical
      - Log_Error
      - Log_Warning
      - Log_Info
      - Log_Debug
 */
typedef ::rlog::LogLevel Level;


/**
    Common base class for all log sinks.

    Log sinks are the classes which receive log output and send it somewhere.
 */
class Sink : public ::rlog::RLogNode
{
public:
  Sink() { }

  /**
      Subscribes to all log messages.
   */
  void SubscribeAll()
  {
    addPublisher(rlog::GetGlobalChannel(""));
  }

protected:
  /**
      This function must be overridden in the derived class to actually log the
      message.

      Notice that this function may be called in a context of another thread.

      @param why the level of the message (error, info, debug, ...)
      @param who name, or kind, of this message (gui, rtd, odbc, ...)
      @param what the message itself
      @param where describes where was it generated from
      @param when the time of the message generation
   */
  virtual void OnLog(Level why,
                     const char *who,
                     const char *what,
                     const Location& where,
                     time_t when) = 0;

  // implement RLogNode::publish() in terms of our OnLog()
  virtual void publish(const rlog::RLogData& data)
  {
    const rlog::PublishLoc * const pub = data.publisher;
    CHECK_VOID( pub, "log message without associated publisher?" );

    const rlog::RLogChannel * const ch = pub->channel;
    CHECK_VOID( ch, "log message without associated channel?" );

    OnLog(ch->logLevel(), ch->name().c_str(), data.msg, *pub, data.time);
  }

  NO_COPY_CLASS(Sink);
};

/**
    This sink implementation outputs the messages to stderr.
 */
class StdioSink : public Sink
{
public:
  StdioSink(int fd = 2 /* stderr */) : m_sinkReal(fd) { }

protected:
  // we override publish() directly so OnLog() shouldn't be ever invoked
  virtual void
  OnLog(Level, const char *, const char *, const Location&, time_t)
  {
    FAIL( "shouldn't be called" );
  }

  virtual void publish(const rlog::RLogData& data)
  {
    // we can't call protected method directly but here we really, really need
    // to so use a cast to allow it
    (static_cast<StdioSink *>(static_cast<rlog::RLogNode *>(&m_sinkReal)))->
      publish(data);
  }

  rlog::StdioNode m_sinkReal;
};

#ifdef _WIN32

/**
    This class output the messages to Windows debug output.
 */
class Win32DebugSink : public Sink
{
public:
  Win32DebugSink() { }

protected:
  virtual void OnLog(Level /* why */,
                     const char *who,
                     const char *what,
                     const Location& where,
                     time_t /* when */)
  {
    // no need for time stamping because any debug output viewer does it anyhow

    std::string loc;
    if ( where.functionName )
      loc = where.functionName;
    if ( loc.empty() )
      loc = String::Printf("at %s:%d", where.fileName, where.lineNum);
    else
      loc = "in " + loc + "()";

    ::OutputDebugStringA
    (
      String::Printf("[%s] %s: %s\r\n", loc.c_str(), who, what).c_str()
    );
  }

  NO_COPY_CLASS(Win32DebugSink);
};

#endif // _WIN32

/**
    This typedef stands for the sink class to use for debug messages.

    It uses stderr under Unix and debug output under Windows.
 */
#ifdef _WIN32
  typedef Win32DebugSink DebugSink;
#else // !_WIN32
  typedef StdioSink DebugSink;
#endif // _WIN32/!_WIN32

} // namespace log

} // namespace ito33

/**
    Defines a logging category.

    The category defined by this macro can be used with ITO33_TRACE and
    ITO33_TRACE_CATEGORY below. This macro can be used only once for the given
    category name, see ITO33_DECLARE_LOG_CATEGORY() if you need to use the same
    log category in multiple files.

    @param name the name of the category object in C++
    @param category the unique string identifying this category
 */
#define ITO33_DEFINE_LOG_CATEGORY(name, category)                             \
  rlog::RLogChannel *name = RLOG_CHANNEL(category)


/**
    Declares a logging category defined elsewhere.

    This macro can be used to reuse the logging category defined using
    ITO33_DEFINE_LOG_CATEGORY() elsewhere.

    @param name the name of the category object in C++
 */
#define ITO33_DECLARE_LOG_CATEGORY(name) extern rlog::RLogChannel *name

/**
    Generic tracing macro.

    The first parameter must be a category object, the second one -- a
    printf()-like format string and the subsequent ones the parameters to be
    inserted into it.
 */
#define ITO33_TRACE rLog

/**
    A shortcut for tracing messages of specific category.

    This is exactly equivalent to ITO33_TRACE with the first parameter set to
    the value of category. The stuff to trace should follow inside parentheses:
    @code
        ITO33_TRACE_CATEGORY(CatObj)("format string", ...);
    @endcode
 */
#define ITO33_TRACE_CATEGORY(category) \
  _rMessage ( LOGID, category, RLOG_COMPONENT)


#ifdef _MSC_VER
  #ifdef NDEBUG
    #pragma comment(lib, "rlog.lib")
  #else
    #pragma comment(lib, "rlogd.lib")
  #endif
#endif

#else // ITO33_NO_LOG

#define ITO33_DEFINE_LOG_CATEGORY(name, category) struct ito33_log_Dummy_##name
#define ITO33_DECLARE_LOG_CATEGORY(name) struct ito33_log_Dummy_##name
#define ITO33_TRACE
#define ITO33_TRACE_CATEGORY(category)

#endif // !ITO33_NO_LOG/ITO33_NO_LOG

#endif // _ITO33_LOG_H_
