/////////////////////////////////////////////////////////////////////////////
// Name:        rtd.cpp
// Purpose:     implementation of the test RTD server
// Author:      Vadim Zeitlin
// Created:     2004-09-20
// RCS-ID:      $Id: rtd.cpp,v 1.6 2005/12/16 16:14:05 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/common.h"
#include "ito33/debug.h"
#include "ito33/thread.h"

#include "ito33/XL/rtd.h"

#include "rtd_h.h"

using namespace ito33;
using namespace ito33::XL;

#ifdef _MSC_VER
  #pragma warning(disable:4355)
#endif

class TestTopic;

class UpdateTimer : public Runnable
{
public:
  UpdateTimer(TestTopic *topic) { m_topic = topic; m_stop = false; }

  virtual void Run();
  void Stop() { m_stop = true; }

private:
  TestTopic *m_topic;
  bool m_stop;
};

class TestTopic : public XL::RTDTopic
{
public:
  TestTopic(XL::RTDServer *server, XL::RTDTopicId id, double value)
    : XL::RTDTopic(id),
      m_server(server),
      m_value(value),
      m_thread(m_timer = new UpdateTimer(this))
  {
  }

  virtual ~TestTopic()
  {
    m_timer->Stop();

    m_thread.Join();
  }

  virtual std::string GetValue() const
  {
    return String::Printf("%.2f", m_value);
  }

  void DoUpdate()
  {
    if ( rand() > RAND_MAX / 2 )
      m_value += 0.5;
    else
      m_value -= 0.5;

    m_server->Update(this);
  }

private:
  XL::RTDServer *m_server;
  UpdateTimer *m_timer;
  JoinableThread m_thread;
  double m_value;
};

class ErrorTopic : public XL::RTDTopic
{
public:
  ErrorTopic(XL::RTDTopicId idTopic, const std::string& msg)
    : XL::RTDTopic(idTopic),
      m_msg(msg)
  {
  }

  virtual std::string GetValue() const { return m_msg; }
  virtual bool HasInitialValue() const { return true; }

private:
  const std::string m_msg;

  NO_COPY_CLASS(ErrorTopic);
};


void UpdateTimer::Run()
{
  while ( !m_stop )
  {
    m_topic->DoUpdate();

    Thread::Sleep(static_cast<unsigned long>(((7.*rand())/RAND_MAX)*1000000));
  }
}


class TestServerImpl : public COM::ImplementCoClassFor<XL::RTDServer,
                                                       TestServerImpl>
{
public:
protected:
    virtual RTDTopicPtr
        OnConnect(const XL::RTDTopicStrings& strings, XL::RTDTopicId id)
  {
    if ( strings.size() == 2 )
    {
      if ( strings[0] == L"spot" )
      {
        double value;
        if ( String::ToDouble(std::string(String::WC2MB(strings[1])), &value) )
          return RTDTopicPtr(new TestTopic(this, id, value));
      }
      else
      {
        return RTDTopicPtr
               (
                 new ErrorTopic(id, String::Printf("ERR: Unknown topic \"%ls\"",
                                                   strings[0].c_str()))
               );
      }
    }
    else
    {
      return RTDTopicPtr
             (
               new ErrorTopic(id,
                              String::Printf("ERR: Wrong number of arguments %d",
                                             strings.size()))
             );
    }

    return RTDTopicPtr();
  }
};

DEFINE_COM_COCLASS_FULL(CLSID_Server, "ITO33RTDTest.Server", 1, TestServerImpl);

