///////////////////////////////////////////////////////////////////////////////
// Name:        ito33/XL/rtd.h
// Purpose:     declaration of RTDServer and related classes
// Author:      Vadim Zeitlin
// Created:     2004-09-20
// RCS-ID:      $Id: rtd.h,v 1.6 2005/06/08 11:44:53 zeitlin Exp $
// Copyright:   (c) 2004 Trilemma LLP
///////////////////////////////////////////////////////////////////////////////

/**
    @name ito33/XL/rtd.h
    @brief Declares a base class for Excel real-time data servers.
 */

#ifndef _ITO33_XL_RTD_H_
#define _ITO33_XL_RTD_H_

#include "ito33/beforestd.h"
#include <map>
#include <set>
#include <vector>
#include "ito33/afterstd.h"
#include "ito33/string.h"

#include "ito33/sharedptr.h"
#include "ito33/thread.h"

#include "ito33/com/traits.h"
#include "ito33/com/dispatch.h"
#include "ito33/com/errorinfo.h"
#include "ito33/com/ptr.h"

#include "ito33/XL/rtdiface.h"

DEFINE_COM_TRAITS(IRtdServer, IUnknown);

namespace ito33
{

/**
    This namespace contains everything related to working with Microsoft Excel.
 */
namespace XL
{

/// Type of topic ids used by Excel.
typedef long RTDTopicId;

/// A list of strings used by the user to choose the topic.
typedef std::vector<std::wstring> RTDTopicStrings;


/**
    A topic provided by RTD server.

    The OnConnect() method of RTDServer returns a new RTDTopic object for each
    topic supported by the server.
 */
class RTDTopic
{
public:
  /**
      Topic constructor associates it with the given id.

      @param id topic id provided by Excel and stored in this object.
   */
  RTDTopic(RTDTopicId id) : m_id(id) { }

  /// Virtual dtor
  virtual ~RTDTopic();

  /// Gets the Excel id of this topic.
  RTDTopicId GetId() const { return m_id; }

  /**
      Return the current value of this topic.
   */
  virtual std::string GetValue() const = 0;

  /**
      Do we have an initial value to display?

      Return true from here if this topic should be updated immediately after
      creation, false if the RTD server should wait until Update() is called
      for it before refreshing it the first time.
   */
  virtual bool HasInitialValue() const { return false; }

private:
  RTDTopicId m_id;

  // it doesn't make sense to copy topics, there is always a single topic for
  // the given id, you should manipulate pointers to them, not the objects
  // themselves
  NO_COPY_CLASS(RTDTopic);
};

typedef SharedPtr<RTDTopic, Policy::MTSafe> RTDTopicPtr;

/**
    Base class for Excel RTD servers.

    Derive a new class from this one and call its Update() method to let Excel
    know that new data is available.

    To define a new RTD server class you must specialize ImplementCoClassOn for
    RTDServer, i.e. write
      @code
        class MyServer : public COM::ImplementCoClassOn<XL::RTDServer>
        {
          ...
        };
      @endcode

    Also, don't forget to use the DEFINE_COM_COCLASS_FULL macro for MyServer!
 */
class RTDServer : public COM::ImplementCustomUnknown
                         <
                            IRtdServer,
                            TYPE_LIST_2(IDispatch, ISupportErrorInfo)
                         >
{
public:
  /**
      Constructor.

      We shouldn't do anything non trivial here, wait until OnStart() is called
      for any real initialization.
   */
  RTDServer() { }

  /**
      Virtual destructor.
   */
  virtual ~RTDServer();


  /**
      Notify Excel that the value of the given topic has changed.

      This method is thread safe in the sense that it can be called from
      another thread without any locking, all calls to Update() are serialized.
   */
  void Update(RTDTopic *topic);


protected:
  /**
      @name RTD server callbacks.

      These methods may be defined in the derived class to do something
      whenever Excel calls the corresponding method of IRtdServer interface.
   */
  //@{

  /**
      Called when the first RTD topic is requested by Excel.

      @return true if ok, false to fail server initialization.
   */
  virtual bool OnStart() { return true; }

  /**
      Called when Excel doesn't need us any longer.
   */
  virtual void OnTerminate() { }

  /**
      Called when Excel tries to subscribe to a new topic.

      @param strings a non empty array of strings specified by the user
      @param id the id for the new topic
      @return a new topic allocated with "new" or NULL if we don't support it
   */
  virtual RTDTopicPtr OnConnect(const RTDTopicStrings& topic,
                                RTDTopicId id) = 0;

  /**
      Called when Excel is not interested in the given topic any longer.

      Note that this method is called if the Excel cell containing the RTD()
      call for this topic is cleared, for example, but is not called if the
      user simply closes the workbook or Excel, so you should @b not rely on it
      being called.

      The topic object is going to be deleted after this function returns, so
      if you just want to do some cleanup, it is better, in the light of above,
      to do it in RTDTopic destructor as it will be also called when the server
      is shutting down.
   */
  virtual void OnDisconnect(RTDTopic * /* topic */) { }

  //@}

private:
  typedef std::map< RTDTopicId, RTDTopicPtr > Topics;
  typedef std::set<RTDTopic *> Updates;


  // implement IRtdServer methods
  STDMETHODIMP ServerStart(IRTDUpdateEvent *callbackObject, long *status);
  STDMETHODIMP ConnectData(long topicId,
                           SAFEARRAY **strings,
                           VARIANT_BOOL *getNewValues,
                           VARIANT *rtdData);
  STDMETHODIMP RefreshData(long *topicCount, SAFEARRAY **data);
  STDMETHODIMP DisconnectData(long topicId);
  STDMETHODIMP Heartbeat(long *rc);
  STDMETHODIMP ServerTerminate();

  // clear the topics map, i.e. destroy all topics and empty the map
  void ClearAllTopics();


  // m_topics contains all currently used topics
  Topics m_topics;

  // m_updates contains all topics which have been updated since last
  // RefreshData() call
  Updates m_updates;

  // the callback object used by Update()
  COM::Ptr<IRTDUpdateEvent> m_pRTDUpdate;

  // critical section protects access to m_updates from Update() and
  // RefreshData()
  CriticalSection m_csUpdate;
};

} // namespace XL

}// namespace ito33

#endif // _ITO33_XL_RTD_H_
