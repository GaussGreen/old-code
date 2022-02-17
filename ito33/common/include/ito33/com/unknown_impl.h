/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/unknown_impl.h
// Purpose:     declaration and definition of COM::ImplementUnknown<> class
// Author:      Vadim Zeitlin
// Created:     23.12.02
// RCS-ID:      $Id: unknown_impl.h,v 1.26 2004/10/05 09:13:36 pedro Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/com/unknown_impl.h
    @brief Defines COM::ImplementUnknown<> template class.

    ImplementUnknown<> provides an as simple as possible imlpementation of the
    IUnknown interface methods. It should be used as the base class for the
    real implementation class.

    This implementation supports different deletion polices (the pointer is
    simply deleted when the ref count reaches 0 by default because this is the
    most common case but this can be changed using the appropriate template
    parameters) and classes supporting multiple interfaces via yet another
    template parameter which is a typelist containing the list of "extra"
    interfaces (there is always the main one from which the implementation
    class derives). Extra interfaces support may require carrying some extra
    data specific to tose interfaces as well and to allow this we provide
    ImplementInterface template class which can be specialized for the
    interfaces of interest to inject arbitrary data into the implementation
    class.
 */

#ifndef _ITO33_COM_UNKNOWNIMPL_H_
#define _ITO33_COM_UNKNOWNIMPL_H_

#include "ito33/debug.h"

#include "ito33/thread.h"
#include "ito33/typelist.h"

#include "ito33/com/traits.h"

namespace ito33
{

namespace COM
{

/**
    This functions returns the number of active COM objects.

    When this count reaches 0, the DLL can be safely unloaded. This function
    doesn't make much sense for an EXE server and shouldn't be called at all in
    that case -- COM::Server::OnNoMoreObjects() callback is used in that case
    instead.

    @return the count of the COM objects still alive.
 */
extern unsigned long GetNumberOfActiveObjects();

// these 2 functions increment/decrement the number of active objects and are
// not meant to be called by anybody but ImplementUnknown<> methods
extern void IncNumberOfActiveObjects();
extern void DecNumberOfActiveObjects();

/**
    Variable which may be set to force a debug break on AddRef/Release done
    on the specified object.

    When debugging memory leaks (caused by too many AddRef()s) or crashes
    (caused by too many Release()s) it is very useful to track all changes
    of the reference count of the given object and to do it you can simply
    set the value of this variable to its address -- then you'll break into
    the debugger each time the ref count is changed.
 */
#ifndef NDEBUG
extern IUnknown *g_objectDebugRef;      // defined in com/utils.cpp
#endif // Debug

// ----------------------------------------------------------------------------
// ImplementUnknown helpers
// ----------------------------------------------------------------------------

/**
    This is just a placeholder used by default for ImplementUnknown second
    template parameter.

    Unfortunately VC++ doesn't support partial template specialization so we
    can't avoid not defining this class at all and specializing
    ImplementUnknown for NoExtraIfaces to not use QueryInterface() on it at
    all. But we can hope that the optimizer can throw away the empty class.
 */
typedef Type::Null NoExtraIfaces;

/**
    This class is used to provide support for the extra interfaces in
    ImplementCustomUnknown.

    Due to VC++ bug we must provide the default implementation for this
    template (normally we would have just forward declared it so that any
    attempt to use it for the interface it wasn't specialized for would
    immediately result in an error message) but this default implementation is
    unusable, you @b must specialize it for any interface you use with
    ImplementCustomUnknown!

    Also, again due to VC++ limitations (lack of support for partial
    specialization), your specialization must have the same form as below i.e.
    the real class must be a template taking one parameter nested inside
    ImplementInterface itself.

    The template parameter of this class is the interface we implement.
 */
template <class Iface>
struct ImplementInterface
{
  /**
      The real implementation of ImplementInterface.

      As explained above, to work around VC++ bugs you must define your real
      class as a template taking the base class as parameter.
   */
  template <class Base>
  struct RealImpl
  {
  };
};

// ----------------------------------------------------------------------------
// ImplementUnknown and ImplementCustomUnknown classes
// ----------------------------------------------------------------------------

/**
    ImplementUnknown provides an easy way to implement the IUnknown methods in
    a derived class.

    To use ImplementUnknown<> specialization for the given interface, the
    COM::Traits<> specialization for it must be defined. This is typically done
    using the DEFINE_COM_TRAITS() macro.

    Note that this class currently only supports single inheritance and that
    all object of this class must be allocated on the heap using operator new
    as they will be deleted when their ref count reaches 0. If this is
    unacceptable, use ImplementCustomUnknown<> with a non default disposal
    policy.

    Template parameter Iface is the interface this class derives from (and
    which should be implemented by the class deriving from this one)
 */
template <class Iface>
class ImplementUnknown : public Iface
{
public:
  /// make the interface available outside of this template
  typedef Iface Interface;

  /// typedef for the base class of this interface
  typedef typename COM::Traits<Iface>::Base BaseInterface;

  /// ctor sets the reference count to 1
  ImplementUnknown() : m_nCount(1) { IncNumberOfActiveObjects(); }

  /// must have virtual dtor, as any base class
  virtual ~ImplementUnknown() { }

  /// @name IUnknown methods
  //@{

  /// used to verify if this interface supports the given IID
  STDMETHODIMP QueryInterface(REFIID riid, void **ppObj);

  /// increments the reference count
  STDMETHODIMP_(ULONG) AddRef();

  /// decrements the reference count and deletes the object if it reaches 0
  STDMETHODIMP_(ULONG) Release();

  //@}

  /**
      Returns true if the interface supports the given IID.

      This template function works by walking the inheritance tree until it
      reaches IUnknown for which is has a specialization.
   */
  static bool IsSupportedInterface(REFIID riid);

protected:
  /**
      This function is called when our ref count reaches 0.

      Having it here is inefficient because it means that potentially
      unnecessary code is used (the policy class in ImplementCustomUnknown
      possibly duplicates it) and we also pay for the (very often
      unnecessary) virtual function call overhead. However it's not that bad
      as the objects are not deleted that often and, anyhow, there is no way
      to do it otherwise without compiler support for partial template
      specialization which VC++ lacks.
   */
  virtual void DoDispose();

private:
  /// the reference count, the object is destroyed when it reaches zero
  unsigned long m_nCount;
};

/**
    Specialization of ImplementInterface for IUnknown.

    This specialization is trivial but allows us to work around yet another
    internal compiler error in VC++ and also shows how to write specializations
    of ImplementInterface.
 */
template <>
struct ImplementInterface<IUnknown>
{
  /**
      Real implementation (which must be a nested class for VC++).

      Notice how we cheat here: we know that IUnknown is always the last
      interface in the list and so Base must be just a simple interface (and
      not some complicated class inheriting from the interface via many
      intermediaries) and so we can inherit from ImplementUnknown<> knowing
      that it, in turn, inherits from Interface which is by now the same as
      Base.
   */
  template <class Base>
  struct RealImpl : public ImplementUnknown<Base>
  {
  };
};


/**
    This namespace contains the disposal policies defined for ImplementUnknown.

    Using a namespace makes it possible to make it more clear that we are
    referring to the policy classes when using them as template arguments
    without having to give unpronounceable names to the policy classes.
 */
namespace DisposalPolicy
{

/**
    This is the default disposal policy for ImplementCustomUnknown and simply
    deletes the object using operator delete.

    This only works with the objects allocated on the heap!

    Note that as ImplementUnknown doesn't derive from any base class with a
    virtual dtor (IUnknown doesn't have a dtor) we must pass the type of the
    object itself to this policy so that the correct destructor is invoked by
    operator delete.
 */
template <class T>
struct Delete
{
  /// simplest possible Dispose() implementation: delete the object
  static void Dispose(T *object) { delete object; }
};

/**
    This disposal policy doesn't do anything at all.

    It can be used with the objects which will be deleted anyhow later on and
    which are not allocated on the heap making usage of Delete policy
    impossible.
 */
struct DoNothing
{
  /// this Dispose() implementation doesn't do anything at all (by design)
  static void Dispose(void * /* object */ ) { }
};

/// symbolic name for the default disposal policy
template <class Iface>
struct Default : public Delete< ImplementUnknown<Iface> >
{
};

} // namespace DisposalPolicy


/**
    This namespace contains helper templates needed by ImplementCustomUnknown.
 */
namespace Impl
{
  /**
      Helper class needed by VC++.

      When an inner template class directly inherits from another template
      class VC++ dies with an internal compiler error and we need an extra
      level of indirection to keep it happy.
   */
  template <class T> struct DeriveFrom : T { };

  /**
      DeriveFromAll helper.

      What you see here are theu usual ugly work arounds for the lack of
      partial template specialization in VC++. It hurts the eye but it does
      work -- and simple, elegant, standard-conforming code doesn't, so ...
   */

  /// first, we need a forward declarations (otherwise VC++ dies)
  template <typename TL>
  struct DeriveFromAllHelper;

  /**
      next we provide a generic implementation (before the specialization or
      VC++ dies)
   */
  template <typename TL>
  struct DeriveFromAllHelper
  {
    /**
        Recursively-defined class inheriting from the implementation for
        ths first typelist element and itself for the remaining elements.
     */
    template <class Root>
    struct Inner : ImplementInterface<typename TL::Head>::
            template RealImpl<typename
                DeriveFromAllHelper<typename TL::Tail>::
                  template Inner<Root>
            >
    {
    };
  };

  /**
      Specialization stopping the recursion.
   */
  template <>
  struct DeriveFromAllHelper<Type::Null>
  {
    /**
        The last class in the hierarchy is trivial.

        Notice that here Root and Interface are actually the same thing.
     */
    template <class Root>
    struct Inner : public ImplementInterface<IUnknown>::template RealImpl<Root>
    {
    };
  };


  /**
      DeriveFromAll may be used to generate a linear hierarchy of classes.

      Instantiating DeriveFromAll for the given type list containing the
      interfaces and the base class will result in a class deriving from all
      ImplementInterface instantiations for each interface and, at the root
      of the hierarchy, from the base class.
   */
  template <typename TL, class Root> 
  struct DeriveFromAll : public DeriveFromAllHelper<TL>::template Inner<Root>
  {
  };
} // namespace Impl

/**
    More flexible IUnknown implementation.

    NB: this comment is out of date, will fix a.s.a.p.

    This class provides slightly more flexibility for the situations when the
    stock IUnknown implementation is not suitable, for example if you don't
    want to systematically delete the object when its ref count reaches 0. This
    class allows you to specify a policy, i.e. a class which has a Dispose()
    method (see the standard Delete and DoNothing policies) as the second
    template parameter instead of always deleting the object.

    The last template parameter may be used to tackle on support for some
    utility interfaces to the implementation class. For example, currently it
    is used to implement IErrorInfo-related stuff (see ito33/com/errorinfo.h)

    Note that it is unfortunately necessary to have a separate class for these
    features as we wouldn't be able to specialized ImplementUnknown<> for
    IUnknown, as we have to for IsSupportedInterface() implementation, if
    ImplementUnknown had more than one template parameter as the partial
    template specialization is not supported by VC++.

    Template parameters:
        - Iface is the interface this class derives from (and which should be
          implemented by the class deriving from this one)
        - Extra is an optional typelist such that ImplementInterface
          is specialized for each its element; by default it is NoExtraIfaces
          and so no extra interfaces are supported
        - DisposalPolicy is the policy class determining what to do when our
          reference count reaches 0 and is Delete by default meaning that
          the object simply deletes itself (and so this only works with the
          objects allocated on the heap!)
 */
template <class Iface,
          class Extra = NoExtraIfaces,
          class DispPolicy = DisposalPolicy::Delete< ImplementUnknown<Iface> > >
class ImplementCustomUnknown : public Impl::DeriveFromAll
                                      <
                                        Extra,
                                        Iface
                                      >
{
public:
  /// provide typedef for the extra class so that it can be used outside
  typedef Extra ExtraIfaces;

  /// more readable synonym for our base class
  typedef ImplementUnknown<Iface> BaseUnknown;

protected:
  /// override to not always delete the object
  virtual void DoDispose() { DispPolicy::Dispose(this); }
};

// ============================================================================
// implementation from now on
// ============================================================================

// ----------------------------------------------------------------------------
// ImplementUnknown
// ----------------------------------------------------------------------------

// this method uses compile-time recursion, it stops thanks to IUnknown
// specialization below
template <class Iface>
inline bool
ImplementUnknown<Iface>::IsSupportedInterface(REFIID riid)
{
  return riid == COM::Traits<Interface>::Uuid() ||
        ImplementUnknown<BaseInterface>::IsSupportedInterface(riid);
}

template <class Iface>
STDMETHODIMP
ImplementUnknown<Iface>::QueryInterface(REFIID riid, void **ppObj)
{
  ASSERT_MSG( ppObj, "ppObj can't be NULL in IUnknown::QueryInterface" );

  if ( IsSupportedInterface(riid) )
  {
    AddRef();

    *ppObj = static_cast<void *>(this);
    return S_OK;
  }

  *ppObj = NULL;

  return E_NOINTERFACE;
}

template <class Iface>
inline STDMETHODIMP_(ULONG) ImplementUnknown<Iface>::AddRef()
{
#ifndef NDEBUG
  if ( this == g_objectDebugRef )
    Trap();
#endif // Debug

  return ++m_nCount;
}

template <class Iface>
inline STDMETHODIMP_(ULONG) ImplementUnknown<Iface>::Release()
{
#ifndef NDEBUG
  if ( this == g_objectDebugRef )
    Trap();
#endif // Debug

  if ( !--m_nCount )
  {
    // NB: order of calls below is important for ImplementCoClass to work!

    // kill the object first
    DoDispose();

    // now that the object is really dead, don't count it as alive any more
    DecNumberOfActiveObjects();

    return 0;
  }

  return m_nCount;
}

template <class Iface>
inline void ImplementUnknown<Iface>::DoDispose()
{
  // simply destroy the object -- it better be allocated on the heap
  delete this;
}

// ----------------------------------------------------------------------------
// speciailization of ImplementUnknown for IUnknown
// ----------------------------------------------------------------------------

// this is where the recursion ends
template <>
inline bool
ImplementUnknown<IUnknown>::IsSupportedInterface(REFIID riid)
{
  return IsEqualIID(riid, IID_IUnknown) != 0;
}

} // namespace COM

} // namespace ito33

#endif // _ITO33_COM_UNKNOWNIMPL_H_

