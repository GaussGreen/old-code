/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/coclass.h
// Purpose:     COM coclasses support stuff
// Author:      Vadim Zeitlin
// Created:     29.01.03
// RCS-ID:      $Id: coclass.h,v 1.17 2006/04/14 01:50:35 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/com/coclass.h
    @brief Support for coclasses.

    COM coclasses are the objects which can be created directly. They are
    identified by their CLSID or, sometimes, by their progid. The global
    DllGetClassObject() must know about all coclasses to be able to create the
    objects of the corresponding type and this is where CoClassInfo is used.
 */

#ifndef _ITO33_COM_COCLASS_H_
#define _ITO33_COM_COCLASS_H_

#include "ito33/com/server.h"
#include "ito33/com/classfactory_impl.h"

// disable warning about
//
//  "typedef-name 'identifier1' used as synonym for class-name 'identifier2'"
//
// which VC++ gives when weu se BaseUnknown below and which doesn't make sense
// at all here
#ifdef _MSC_VER
  #pragma warning(disable:4097)
#endif

namespace ito33
{

namespace COM
{

/**
    All information we need about a single coclass.

    CoClassInfo is used to get information about a coclass and, most
    importantly, to create the objects of the coclass' type. To make it
    possible to find the coclass from the given CLSID, all coclasses are kept
    in a (single) linked list and we provide the Find() method to find a
    coclass by either prog or clsid.

    The objects of this class are usually global static and are created using
    DEFINE_COM_COCLASS() macro and they are never destroyed until the program
    end.
 */
class CoClassInfo
{
public:
  /// the type of the function which may be used to create a new COM objcet
  typedef IUnknown *(*Create)();

  /**
      Creates an object containing information about the coclass.

      This ctor shouldn't be used directly, use DEFINE_COM_COCLASS() macro
      instead.
   */
  CoClassInfo(const CLSID& clsid,
        const char *progid,
        int version,
        Create create)
    : m_clsid(clsid), m_progid(progid), m_version(version), m_create(create)
  {
    // insert us in the head of the linked list (it is simpler than
    // inserting in the tail and the order doesn't matter here anyhow)
    m_next = ms_first;
    ms_first = this;
  }

  /**
      CLSID of this coclass: its unique identifier.

      Each coclass is given an UUID in the .idl file, this is the same UUID
      which should be used here (CLSID is just a synonym for UUID)
   */
  const CLSID& GetClsid() const { return m_clsid; }

  /**
      progid is a short huma readable string identifying the coclass

      It is usually of the form vendor.application or application.category,
      e.g. Word.Document. It is used to identify a class, like the CLSID,
      but with less precision because the conflicts are possible when using
      progid. OTOH, it is human readable (unlike CLSID). Nevertheless, it
      shouldn't appear directly in the UI, in particular this string is never
      translated.

      Also, the progid is used as the name of one of the registry keys
      associated with the COM server and so there are several restrictions
      associated with it:
      - it must have no more than 39 characters
      - it can't contain any punctuation except for the periods
      - but it must contain at least one period
      - it must start with a letter (i.e. not a digit nor period)
      - it must be as unique as possible
   */
  const char *GetProgid() const { return m_progid; }

  /**
      version is used together with progid to form the full progid.

      What we call progid here is, in fact, the "version indepedent progid"
      in COM terms and the full progid is the combination of that progid with
      the version, e.g. the full progid may be Word.Document.8 if the version
      is 8.
   */
  int GetVersion() const { return m_version; }

  /**
      Create a new object implementing this coclass interface(s).

      This function is used by DllGetClassObject() and our class factory
      implementation to create the object of this coclass.
   */
  IUnknown *CreateInstance() const { return (*m_create)(); }


  /**
      Return the coclass corresponding to the given CLSID.

      If there is a coclass whose GetClsid() is equal to the given CLSID,
      return it, otherwise returns @c NULL.
   */
  static const CoClassInfo *Find(const CLSID& clsid);

  /**
      Get the pointer to the next node in the linked list.

      The CoClassInfo structures are meant to be arranged in a linked list
      and, in fact, the ctor always adds the new objects to the list so it
      will always be like this. The next pointer either points to the
      information about the next coclass in the chain or is @c NULL.
    */
  const CoClassInfo *GetNext() const { return m_next; }

  /**
      Get the pointer to the first node in the linked list.

      This method is used to start iterating over all elements of the linked
      list and GetNext() is used to continue it.
   */
  static const CoClassInfo *GetFirst() { return ms_first; }


private:
  // we cannot have -- and don't need --- neither assignment operator nor the
  // copy ctor
  CoClassInfo(const CoClassInfo&);
  CoClassInfo& operator=(const CoClassInfo&);

  // the CLSID
  const CLSID& m_clsid;

  // version independent progid
  const char *m_progid;

  // the version part of progid
  int m_version;

  // function to create the object of this coclass.
  Create m_create;

  // the next object in the linked list
  CoClassInfo *m_next;

  // the hea of the linked list (may be NULL)
  static CoClassInfo *ms_first;
};

/**
    Macro which must be used for all the existing coclasses.

    If you forget to use this macro, the objects of the corresponding coclass
    won't be created during run-time because DllGetClassObject() wouldn't know
    about them (and also because DllRegisterServer() would never register them
    in the first place)

    This macro must be used at the file level.

    @param clsid the coclass uuid (@c CLSID_XXX)
    @param progid the version independent progid
    @param version the version of the full progid
    @param classname the name of the class which implements this coclass
 */
#define DEFINE_COM_COCLASS_FULL(clsid, progid, version, classname)          \
  struct CreatorFor ## classname                                          \
  {                                                                       \
    static IUnknown *Create() { return new classname; }                 \
  };                                                                      \
                                      \
  static ::ito33::COM::CoClassInfo CoClassInfoFor ## classname              \
           (                                                    \
              clsid,                                          \
              progid,                                         \
              version,                                        \
              CreatorFor ## classname::Create                 \
           )

/**
    Easier to use and shorter form of DEFINE_COM_COCLASS_FULL().

    Instead of taking many different arguments this macro takes only one which
    is the name of the coclass. However this only works if the standard naming
    conventions are used, i.e. the class implementing the COM coclass has the
    same name with "Impl" suffix and so on.

    It also allows to give only the second part of progid and constructs the
    full progid using COM_APPNAME.

    @param name the name of the coclass
    @param version the COM version of the coclass
 */
#define DEFINE_COM_COCLASS(name, version)                                   \
  DEFINE_COM_COCLASS_FULL(CLSID_ ## name,                                 \
              COM_APPNAME "." #name,                          \
              version,                                        \
              name ## Impl)

#ifndef COM_APPNAME
  /**
      This macro should be defined in the user code to contain the app name.

      you should define COM_APPNAME to contain the application name before
      using DEFINE_COM_COCLASS()
   */
  #define COM_APPNAME ""
#endif

/**
    This is a more general form of ImplementCoClass which may be used if the
    custom IUnknown implementation we want to create objects of is already
    defined elsewhere.

    This class adds support for IClassFactory in its QueryInterface() which is
    required for the coclasses because the COM clients always start by asking
    for a class factory object.

    Note that this class allows the implementation class to support some extra
    interfaces in the same manner as ImplementCustomUnknown does.

    The template parameters is the base ImplementCustomUnknown-derived class.
 */
template <class Unknown, class ImplClass>
class ImplementCoClassFor : public Unknown
{
public:
  /// provide a more readable name for the base ImplementUnknown
  typedef Unknown BaseUnknown;

  /// Remove this object from the "active" list artificially in ctor
  ImplementCoClassFor()
  {
    // the existence of this object shouldn't prevent the server from shutting
    // down
    DecNumberOfActiveObjects();
  }

  /// Readd this object to the "active" list
  virtual ~ImplementCoClassFor()
  {
    // if we don't do this, the number of active objects could become negative
    // when it is decremented in the base class as we had artificially
    // decreased it in the ctor
    IncNumberOfActiveObjects();
  }

  /**
      @name Special reference counting for coclasses.

      We subvert the normal reference counting in the case of class factories
      which is needed to make the out of process servers work. The problem here
      is that we create the class factories ourselves on startup (to register
      them using CoRegisterClassObject()) and if we didn't do anything, these
      objects would always stay alive and the server would never terminate.
      OTOH, if we didn't count these objects as alive at all, the server could
      terminate even if an (external) client still had a lock on one of our
      class factories.

      Logically, the solution is to count a class factory in an EXE server
      active only when its ref count is > 1 and not when it is > 0. And this is
      what we implement here.
   */
  //@{

  STDMETHODIMP_(ULONG) AddRef()
  {
    ULONG cRef = BaseUnknown::AddRef();
    if ( cRef == 2 )
    {
      Server *server = Server::Get();
      if ( server )
        server->Lock();
    }

    return cRef;
  }

  STDMETHODIMP_(ULONG) Release()
  {
    ULONG cRef = BaseUnknown::Release();
    if ( cRef == 1 )
    {
      // undo the Lock() above: this Unlock() may trigger the server shut down
      Server *server = Server::Get();
      if ( server )
        server->Unlock();
    }

    return cRef;
  }

  //@}

  /// override QueryInterface() to know about IClassFactory
  STDMETHODIMP QueryInterface(REFIID riid, void **ppObj);

private:
  NO_COPY_CLASS(ImplementCoClassFor);
};

/**
    The base which should be used for the coclass classes instead of
    ImplementUnknown.

    The template parameters are:
        - Iface the (abstract) interface implemented by the class
        - ImplClass the (concrete) class implementing this interface
        - ExtraIfaces type list of extra interfaces supported by the class
 */
template <class Iface,
          class ImplClass,
          class ExtraIfaces = NoExtraIfaces>
class ImplementCoClass : public ImplementCoClassFor
                                <
                                  ImplementCustomUnknown
                                  <
                                    Iface,
                                    ExtraIfaces,
                                    DisposalPolicy::Default<Iface>
                                  >,
                                  ImplClass
                                >
{
};

// ----------------------------------------------------------------------------
// implementation from now on
// ----------------------------------------------------------------------------

template <class Unknown, class ImplClass>
STDMETHODIMP
ImplementCoClassFor<Unknown, ImplClass>::QueryInterface(REFIID riid,
                                                        void **ppObj)
{
  if ( riid == IID_IClassFactory )
  {
    // TODO: we could cache the pointer here, should we?
    *ppObj = new ImplementClassFactory<ImplClass>;

    return S_OK;
  }

  return BaseUnknown::QueryInterface(riid, ppObj);
}

} // namespace COM

} // namespace ito33

#ifdef _MSC_VER
  #pragma warning(default:4097)
#endif

#endif // _ITO33_COM_COCLASS_H_

