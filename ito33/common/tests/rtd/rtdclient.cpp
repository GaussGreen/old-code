/////////////////////////////////////////////////////////////////////////////
// Name:        tests/rtd/rtdclient.cpp
// Purpose:     test RTD client emulating the way Excel uses RTD server
// Author:      Vadim Zeitlin
// Created:     2006-05-19
// RCS-ID:      $Id: rtdclient.cpp,v 1.4 2006/05/20 00:22:55 zeitlin Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#define INITGUID

#include "ito33/com/init.h"
#include "ito33/com/safearray.h"
#include "ito33/com/dispatch.h"
#include "ito33/com/unknown_impl.h"

#include "ito33/XL/rtdiface.h"

#include <string>
#include <vector>

#include <stdio.h>
#include <conio.h>

#import "rtd.tlb"

DEFINE_GUID(IID_IRTDUpdateEvent,0xA43788C1,0xD91B,0x11D3,0x8F,0x39,0x00,0xC0,0x4F,0x36,0x51,0xB8);

using namespace ito33;

// ----------------------------------------------------------------------------
// constants
// ----------------------------------------------------------------------------

const long TOPIC_ID = 17;

// ----------------------------------------------------------------------------
// class used to pass information from update thread
// ----------------------------------------------------------------------------

class ThreadData
{
public:
  ThreadData() { m_updated = false; }

  // set the updated flag, will be reset by next read
  void SetUpdated()
  {
    Lock<CriticalSection> lock(m_cs);
    m_updated = true;
  }

  // check if we were updated, resets the flag
  bool IsUpdated()
  {
    Lock<CriticalSection> lock(m_cs);
    if ( !m_updated )
      return false;

    m_updated = false;
    return true;
  }

private:
  // updated flag
  bool m_updated;

  // protects m_updated
  CriticalSection m_cs;


  NO_COPY_CLASS(ThreadData);
};

// ----------------------------------------------------------------------------
// our RTD callback class
// ----------------------------------------------------------------------------

class ClientRTDUpdate;
DEFINE_COM_IFACE_FULL(RTDUpdateEvent, ClientRTDUpdate, IDispatch);

class ClientRTDUpdate : public COM::ImplementCustomUnknown
                                    <
                                        IRTDUpdateEvent,
                                        TYPE_LIST_1(IDispatch),
                                        COM::DisposalPolicy::DoNothing
                                    >
{
public:
  ClientRTDUpdate(ThreadData& data)
    : m_interval(60),
      m_data(data)
  {
  }

  // NB: this can be called from another thread context!
  virtual HRESULT STDMETHODCALLTYPE UpdateNotify()
  {
    puts("Update notification received.");

    m_data.SetUpdated();

    return S_OK;
  }

  virtual HRESULT STDMETHODCALLTYPE get_HeartbeatInterval(long *plRetVal)
  {
    *plRetVal = m_interval;

    return S_OK;
  }

  virtual HRESULT STDMETHODCALLTYPE put_HeartbeatInterval(long plRetVal)
  {
    m_interval = plRetVal;

    return S_OK;
  }

  virtual HRESULT STDMETHODCALLTYPE Disconnect()
  {
    puts("IRTDUpdateEvent::Disconnect()");

    return S_OK;
  }

private:
  long m_interval;
  ThreadData& m_data;

  NO_COPY_CLASS(ClientRTDUpdate);
};

// ============================================================================
// implementation
// ============================================================================

int main(int argc, char **argv)
{
  try
  {
    // parse command line
    // ------------------

    if ( argc < 3 )
    {
      fprintf(stderr, "Usage: %s server topic [arguments...]\n", argv[0]);
      return 1;
    }

    const std::string server = argv[1];
    const std::vector<std::string> topics(argv + 2, argv + argc);


    // connect to the server
    // ---------------------

    COM::Initialize initCOM;

    ITO33RTDTest::IRtdServerPtr srv(server.c_str());

    ThreadData data;
    ClientRTDUpdate cb(data);
    long rc = srv->ServerStart(static_cast<ITO33RTDTest::IRTDUpdateEvent *>(
                                static_cast<IUnknown *>(&cb)));
    if ( !rc )
    {
      fprintf(stderr, "RTD server startup failed with error code %ld\n", rc);
      return 2;
    }

    puts("RTD server startup successful");

    const size_t numTopics = topics.size();
    COM::SafeArray<VARIANT> aStrings(numTopics);
    COM::SafeArrayAccessor<VARIANT> strings(aStrings);
    for ( size_t n = 0; n < numTopics; ++n )
    {
      strings[n].vt = VT_BSTR;
      strings[n].bstrVal = COM::BasicString(topics[n].c_str()).Detach();
    }

    VARIANT_BOOL gotValue;
    _variant_t value = srv->ConnectData(TOPIC_ID, aStrings, &gotValue);

    printf("Connected to RTD server, ");

    if ( gotValue == VARIANT_TRUE )
    {
      printf("initial topic value is \"%s\"\n",
             (const char *)String::WC2MB(value.bstrVal));
    }
    else
    {
      puts("no initial value");
    }


    // monitor data changes
    // --------------------

    puts("Press a key to exit...");
    for ( ;; )
    {
      if ( srv->Heartbeat() != 1 )
      {
        fprintf(stderr, "RTD server doesn't respond to hearbeat any more\n");
        break;
      }

      ::Sleep(200);

      if ( _kbhit() )
      {
        _getch();
        break;
      }

      if ( data.IsUpdated() )
      {
        long count = 1;
        COM::SafeArray<VARIANT> aData(srv->RefreshData(&count));

        ASSERT_MSG( aData.GetDim() == 2, "unexpected data array format" );
        ASSERT_MSG( aData.GetCount(0) == 1, "unexpected number of topics" );
        ASSERT_MSG( aData.GetCount(1) == 2, "unexpected number of columns" );

        VARIANT *data = aData.Lock();
        ASSERT_MSG( data[0].vt == VT_I4 &&
                      data[0].lVal == TOPIC_ID, "unexpected topic" );
        ASSERT_MSG( data[1].vt == VT_R8, "unexpected topic type" );
        printf("New value:\t%g\n", data[1].dblVal);
        aData.Unlock();
      }
    }

    // shutdown
    // --------

    srv->DisconnectData(TOPIC_ID);
    puts("Disconnected from RTD server");

    printf("Shutting down RTD server...");
    srv->ServerTerminate();
    puts(" done.");

    return 0;
  }
  catch ( const _com_error& ce )
  {
    printf("COM exception: %#08x", ce.Error());

    IErrorInfo *pErrorInfo;
    if ( SUCCEEDED(::GetErrorInfo(0, &pErrorInfo)) && pErrorInfo )
    {
      BSTR bstr;
      pErrorInfo->GetDescription(&bstr);
      if ( *bstr )
        printf(":\n\n%s", (char *)_bstr_t(bstr));

      pErrorInfo->Release();
    }

    putchar('\n');
  }
  catch ( ito33::Exception& e )
  {
    fprintf(stderr, "ito33::Exception:\n%s\n", e.GetFullMessage().c_str());
  }
  catch ( std::exception& e )
  {
    fprintf(stderr, "std::exception: %s\n", e.what());
  }
  catch ( ... )
  {
    fprintf(stderr, "Unexpected exception.\n");
  }

  return -1;
}
