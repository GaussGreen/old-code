///////////////////////////////////////////////////////////////////////////////
// Name:        ito33/XL/rtd.h
// Purpose:     implementation of RTDServer
// Author:      Vadim Zeitlin
// Created:     2004-09-20
// RCS-ID:      $Id: rtd.cpp,v 1.23 2006/05/14 13:03:49 zeitlin Exp $
// Copyright:   (c) 2004 Trilemma LLP
///////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#define INITGUID

#include "ito33/beforestd.h"
#include <algorithm>
#include "ito33/afterstd.h"

#include "ito33/exception.h"
#include "ito33/log.h"
#include "ito33/date.h"

#include "ito33/XL/rtd.h"

#include "ito33/com/bstr.h"
#include "ito33/com/server.h"
#include "ito33/com/safearray.h"

using namespace ito33;
using namespace ito33::XL;

ITO33_DEFINE_LOG_CATEGORY(LogRTD, "rtd");

#define RTD_TRACE     ITO33_TRACE_CATEGORY(LogRTD)

// ----------------------------------------------------------------------------
// some globals
// ----------------------------------------------------------------------------

DEFINE_GUID(IID_IRTDUpdateEvent,0xA43788C1,0xD91B,0x11D3,0x8F,0x39,0x00,0xC0,0x4F,0x36,0x51,0xB8);
DEFINE_GUID(IID_IRtdServer,0xEC0E6191,0xDB51,0x11D3,0x8F,0x3E,0x00,0xC0,0x4F,0x36,0x51,0xB8);

// ----------------------------------------------------------------------------
// shared pointers
// ----------------------------------------------------------------------------

namespace ito33
{
  ITO33_IMPLEMENT_SHAREDPTR(XL::RTDTopic);
}

// ============================================================================
// RTDTopic implementation
// ============================================================================

RTDTopic::~RTDTopic()
{
  // nothing to do here
}

// ============================================================================
// RTDServer implementation
// ============================================================================

// ----------------------------------------------------------------------------
// dtor
// ----------------------------------------------------------------------------

RTDServer::~RTDServer()
{
  ASSERT_MSG( m_topics.empty(), "memory leak: RTDServer::m_topics not freed" );
}

// ----------------------------------------------------------------------------
// other public methods
// ----------------------------------------------------------------------------

void RTDServer::Update(RTDTopic *topic)
{
  CHECK_VOID( topic, "RTDTopic::Update(): NULL topic" );

  Lock<CriticalSection> lock(m_csUpdate);

  // add the topic to the set of updates unless it is already there
  m_updates.insert(topic);

  // also let Excel know about it if we can
  if ( m_pRTDUpdate )
    m_pRTDUpdate->UpdateNotify();
}

// ----------------------------------------------------------------------------
// helpers
// ----------------------------------------------------------------------------

void RTDServer::ClearAllTopics()
{
  m_topics.clear();
}

// ----------------------------------------------------------------------------
// IRTDServer methods implementation
// ----------------------------------------------------------------------------

STDMETHODIMP
RTDServer::ServerStart(IRTDUpdateEvent *callbackObject, long *status)
{
  CHECK( callbackObject && status, E_POINTER,
            "RTDServer::ServerStart(): NULL pointer" );

  // for each call to ServerStart(), Excel calls ServerTerminate(), even if
  // ServerStart() failed, so we have to lock the server (i.e. prevent the DLL
  // from being unloaded by incrementing the active object number) in any case
  // as ServerTerminate() doesn't know whether we succeeded or failed here and
  // calls DecNumberOfActiveObjects() unconditionally
  COM::IncNumberOfActiveObjects();

  HRESULT hr = S_OK;

  try
  {
    m_pRTDUpdate.Assign(callbackObject);

    if ( !OnStart() )
      hr = E_FAIL;
  }
  catch ( const Exception& e )
  {
    hr = SetErrorInfo(e);
  }
  catch ( ... )
  {
    hr = E_UNEXPECTED;
  }

  // any return value < 1 indicates an error here
  if ( FAILED(hr) )
  {
    RTD_TRACE("startup failed");

    *status = 0;
  }
  else // server started
  {
    RTD_TRACE("startup ok");

    *status = 1;
  }

  return hr;
}

STDMETHODIMP
RTDServer::ConnectData(long topicId,
                       SAFEARRAY **ppStrings,
                       VARIANT_BOOL *getNewValues,
                       VARIANT *rtdData)
{
  CHECK( ppStrings && *ppStrings && getNewValues && rtdData, E_POINTER,
            "RTDServer::ConnectData(): NULL pointer" );

  try
  {
    // just a safety check: we're expecting a list of strings
    ASSERT_MSG( (*ppStrings)->cDims == 1,
                  "RTDServer::ConnectData(): expect 1D array, got ...?" );

    // get all components of the topic name in a vector
    COM::SafeArray<VARIANT> aStrings(*ppStrings);
    COM::SafeArrayAccessor<VARIANT> strings(aStrings);

    const size_t count = aStrings.GetCount();
    RTDTopicStrings topicComponents(count);
    for ( size_t n = 0; n < count; ++n )
    {
      topicComponents[n] = strings[n].bstrVal;
    }

    // TODO: what to do if we've already got this topic?

    // create a topic object for this topic string
    RTDTopicPtr topic(OnConnect(topicComponents, topicId));
    if ( !topic )
    {
      RTD_TRACE("no such topic %ld", topicId);

      return E_FAIL;
    }

    // remember this topic with its id
    m_topics[topicId] = topic;

    // do we already have any new results for this topic?
    *getNewValues = topic->HasInitialValue() ? VARIANT_TRUE : VARIANT_FALSE;
    if ( *getNewValues )
    {
      rtdData->vt = VT_BSTR;
      rtdData->bstrVal = COM::BasicString(topic->GetValue().c_str()).Detach();
    }

    RTD_TRACE("topic %ld connected", topicId);

    return S_OK;
  }
  catch ( ... )
  {
    RTD_TRACE("connecting to topic %ld failed", topicId);

    return E_UNEXPECTED;
  }
}

STDMETHODIMP RTDServer::RefreshData(long *topicCount, SAFEARRAY **data)
{
  CHECK( topicCount && data, E_POINTER,
            "RTDServer::RefreshData(): NULL pointers" );

  Updates updates;
  size_t count;

  // block inside which m_csUpdate is locked
  {
    Lock<CriticalSection> lock(m_csUpdate);

    count = m_updates.size();
    if ( !count )
    {
      // this does happen sometimes, no idea why XL decides to call us when we
      // didn't ask it to but it does (this happens when switching from/to XL
      // to/from another application usually)
      RTD_TRACE("RTDServer::RefreshData() called but no updates");

      *topicCount = 0;
      return S_OK;
    }

    // there are 2 reasons for doing this here:
    //  1. even if an error happens, we should still flush m_updates, so ensure
    //     it is going to be empty
    //  2. we must not call GetValue() below while keeping m_csUpdate lock
    //     because this can easily result in a deadlock if the user topic
    //     implementation uses (another) critical section to protect something
    //     in its GetValue() and code which calls our Update(): as the order of
    //     acquisition of the locks is not necessarily the same, we could well
    //     deadlock in GetValue() call if we Update() is waiting for m_csUpdate,
    //     so we must release this lock a.s.a.p.
    m_updates.swap(updates);
  }


  try
  {
    // fill the 2D data safearray with our pairs id / value
    const size_t cols = 2;
    const size_t rows = count;

    HRESULT hr = S_OK;

    // initialize the bounds structure
    SAFEARRAYBOUND bounds[2];
    bounds[0].lLbound =
    bounds[1].lLbound = 0;
    bounds[0].cElements = cols;
    bounds[1].cElements = static_cast<ULONG>(rows);

    // and create the SAFEARRAY.
    *data = ::SafeArrayCreate(VT_VARIANT, SIZEOF(bounds), bounds);
    if ( !*data )
      throw COM_EXCEPTION("SafeArrayCreate", E_UNEXPECTED);

    // lock the array and copy the data.
    VARIANT *arrayData;
    hr = ::SafeArrayAccessData(*data, reinterpret_cast<void **>(&arrayData));
    if ( FAILED(hr) )
      throw COM_EXCEPTION("SafeArrayAccessData", hr);

    // now copy data for all updates topics
    *topicCount = 0;
    for ( Updates::const_iterator i = updates.begin(),
                                end = updates.end();
          i != end;
          ++i, ++*topicCount )
    {
      // the id
      VARIANT& currentId = *arrayData++;
      currentId.vt = VT_I4;
      currentId.lVal = (*i)->GetId();

      // and the value
      VARIANT& currentVal = *arrayData++;

      std::string value;
      try
      {
        value = (*i)->GetValue();
      }
      catch ( ito33::Exception& e )
      {
        value = e.GetErrorMessage();
      }

      // try to guess the correct type of the value: it's important to return
      // numbers as such to Excel as otherwise it doesn't display them properly
      // and operations on them such as comparisons don't work neither
      if ( String::ToLong(value, &currentVal.lVal) )
      {
        currentVal.vt = VT_I4;
      }
      else if ( String::ToDouble(value, &currentVal.dblVal) )
      {
        currentVal.vt = VT_R8;
      }
      else // not a numeric type
      {
        // try to interpret it as date using both local and standard ISO formats
        Date dt;
        const char *p = dt.Parse(value.c_str());
        if ( !p )
          p = dt.Parse(value.c_str(), "%Y-%m-%d");
        if ( p && !*p ) // check whether we parsed the entire string
        {
          currentVal.vt = VT_DATE;
          currentVal.date = dt.GetOleDate();
        }
        else // if all else fails, interpret it as a string
        {
          currentVal.vt = VT_BSTR;
          currentVal.bstrVal = COM::BasicString(value.c_str()).Detach();
        }
      }
    }

    ASSERT_MSG( *topicCount == (LONG)count,
                    "RefreshData(): not all topics copied?" );

    hr = ::SafeArrayUnaccessData(*data);
    if ( FAILED(hr) )
      throw COM_EXCEPTION("SafeArrayUnaccessData", hr);

    RTD_TRACE("%ld topics fetched", *topicCount);

    return S_OK;
  }
  catch ( ... )
  {
    RTD_TRACE("refreshing data failed");

    return E_UNEXPECTED;
  }
}

STDMETHODIMP RTDServer::Heartbeat(long *rc)
{
  CHECK( rc, E_POINTER, "RTDServer::Heartbeat(): NULL output pointer" );

  RTD_TRACE("toc-toc");

  // return 1 to indicate that the server is running
  *rc = 1;

  return S_OK;
}

STDMETHODIMP RTDServer::DisconnectData(long topicId)
{
  Topics::iterator i = m_topics.find(topicId);
  if ( i != m_topics.end() )
  {
    RTDTopicPtr topic(i->second);

    // remove any pending updates for this topic
    {
      Lock<CriticalSection> lock(m_csUpdate);
      m_updates.erase(topic.get());
    }

    // delete the topic and remove the pointer to it from map
    m_topics.erase(i);
  }
  //else: it is possible to get DisconnectData() for a topic we don't have if
  //      we returned E_FAIL for it from ConnectData()

  RTD_TRACE("topic %ld disconnected", topicId);

  return S_OK;
}

STDMETHODIMP RTDServer::ServerTerminate()
{
  try
  {
    // whether we succeed or fail, we can be unloaded now, so do this first
    COM::DecNumberOfActiveObjects();

    // forget the update pointer, can't do updates any more
    m_pRTDUpdate.Release();

    // destroy all the topics (may be left if DisconnectData() wasn't called)
    ClearAllTopics();

    RTD_TRACE("server shut down ok");

    // notify the derived class about server shut down
    OnTerminate();

    return S_OK;
  }
  catch ( ... )
  {
    RTD_TRACE("server shut down failed");

    return E_UNEXPECTED;
  }
}

