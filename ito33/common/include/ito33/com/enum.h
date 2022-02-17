/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/enum.h
// Purpose:     COM enumerators
// Author:      Vadim Zeitlin
// Created:     25.03.03
// RCS-ID:      $Id: enum.h,v 1.11 2006/08/19 22:19:05 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/com/enum.h
    @brief COM enumerators support.

    COM uses IEnumXXX classes to iterate over the collections of objects. We
    provide here an IEnumVARIANT (the only enumerator used by OLE Automation
    implementation which can be used to iterate over any standard C++
    container.

    Note that we could easily be more general and provide templates for
    implementing IEnumXXX iterating on values of any type but for now such
    flexibility doesn't seem to be needed. When/if it becomes necessary, we'd
    just need to add more template parameters (and as they could have the
    default values we wouldn't even break compatibility with the old code) to
    the classes below.
 */

#ifndef _ITO33_COM_ENUM_H_
#define _ITO33_COM_ENUM_H_

#include "ito33/sharedptr.h"

#include "ito33/com/unknown_impl.h"
#include "ito33/com/variant.h"
#include "ito33/com/ptr.h"

DEFINE_COM_TRAITS(IEnumVARIANT, IUnknown);

namespace ito33
{

namespace COM
{

namespace Private
{
    /**
        Traits for the class being enumerated.

        In the case of primitive types we can simply return the object being
        enumerated as itself (as is done in the specialization of this class
        below), but when enumerating the objects we should create interface
        pointers wrapping them 

        Note that we suppose that all object classes in collections are OLE
        Automation compatible and so implement IDispatch (this does make sense
        as collections are OLE Automation notion anyhow).
     */
    template <class T> struct EnumTraits
    {
        static T GetValue(T value) { return value; }
    };

#if 0
    /**
        Specialization of EnumTraits for object pointers.
     */
    template <class T> struct EnumTraits< COM::Ptr<T> >
    {
        static IDispatch *GetValue(const COM::Ptr<T>& ptr)
        {
            // it is important to return VT_DISPATCH and not just VT_UNKNOWN as
            // otherwise VB complains about "object required"
            return static_cast<IDispatch *>(ptr.Get());
        }
    };
#endif

    /**
        Specialization of EnumTraits for objects.

        If the container stores directly the objects, we have to create a copy
        of the original object in this case which means that modifying
        collection elements is not going to work. For this you need to have a
        container of shared pointers.
     */
    template <class T> struct EnumTraits<T *>
    {
        typedef T Iface;
        typedef typename COM::Traits<Iface>::Impl Impl;
        typedef typename COM::Traits<Iface>::Class Class;

        static IDispatch *GetValue(const Class& obj)
        {
            return new Impl(shared_ptr<Class>(new Class(obj)));
        }

        static IDispatch *GetValue(const shared_ptr<Class>& ptr)
        {
            return new Impl(ptr);
        }
    };
} // namespace Private


/**
    StdEnum class implements IEnumVARIANT using a standard C++ container.

    This is a template class which can be used with any standard container.
 */
template <class T, class U = typename T::value_type>
class StdEnum : public ImplementUnknown<IEnumVARIANT>
{
public:
    /// ContainerType is the type of the container we iterate on
    typedef T ContainerType;

    /// COMElementType is the type of the enumerated elements on COM side
    typedef U COMElementType;

    /// ElementType is the type of the elements we enumerate
    typedef typename ContainerType::value_type ElementType;

    /// The traits of the enumerated type
    typedef typename Private::EnumTraits<COMElementType> Traits;


    /**
        Creates an enumerator for the elements in the given container.

        We don't copy the container but just keep a reference to it, so the
        containers life time should be greater than that of the enumerator
        object and it also shouldn't change while we exist.

        @param container contains the items to enumerate
     */
    StdEnum(const ContainerType& container);

    /**
        @name IEnumVARIANT methods

        See the MSDN docs for further details about these methods.
     */
    //@{

    /**
        Returns next batch of elements.

        @param celt the number of elements to try to get, must be positive
        @param rgelt the array where the elements are to be put
        @param pceltFetched pointer to where the number of elements returned
                            is to be put, may only be @c NULL if celt is 1
        @return S_OK if the number of elements returned is celt, S_FALSE if it
                is less than celt or E_XXX if an error occured
     */
  STDMETHODIMP Next(ULONG celt, VARIANT *rgelt, ULONG *pceltFetched);

    /**
        Skip the given number of elements.

        @param celt the number of elements by which the current position should
                    be advanced
        @return S_OK if ok or S_FALSE if there are not enough elements
     */
  STDMETHODIMP Skip(ULONG celt);

    /**
        Rewinds the enumerators back to the beginning.

        @return always S_OK in this implementation
     */
  STDMETHODIMP Reset(void);

    /**
        Creates and returns an exact copy of this enumerator.

        @param ppEnum the out pointer to returned enumerator
        @return S_OK if ok, E_XXX if an error occurs
     */
  STDMETHODIMP Clone(IEnumVARIANT **ppEnum);

    //@}

private:
    // common part of ctor and Reset()
    void Init();


    // the type for the iterators for our container
    typedef typename ContainerType::const_iterator IteratorType;

    // the items we enumerate
    const ContainerType& m_container;

    // the current position
    IteratorType m_current;


    // enumerators can't be copied, use Clone() instead
    StdEnum(const StdEnum&);
    StdEnum& operator=(const StdEnum&);
};

// ----------------------------------------------------------------------------
// helper functions
// ----------------------------------------------------------------------------

/**
    Creates a standard enumerator for the given container.

    The advantage of using this function is that it can automatically deduce
    the template parameter from its argument type.
 */
template <class U, class T>
inline
StdEnum<T, U> *CreateStdEnum(const T& container)
{
    return new StdEnum<T, U>(container);
}

// ----------------------------------------------------------------------------
// inline functions implementation
// ----------------------------------------------------------------------------

template <class T, class U>
inline
void StdEnum<T, U>::Init()
{
    m_current = m_container.begin();
}

template <class T, class U>
inline
StdEnum<T, U>::StdEnum(const typename StdEnum::ContainerType& container)
          : m_container(container)
{
    Init();
}

template <class T, class U>
STDMETHODIMP StdEnum<T, U>::Next(ULONG celt, VARIANT *rgelt, ULONG *pceltFetched)
{
    if ( !rgelt || (!pceltFetched && celt > 1) )
        return E_POINTER;

    if ( pceltFetched )
        *pceltFetched = 0;

    for ( const IteratorType end = m_container.end(); celt--; ++m_current )
    {
        if ( m_current == end )
        {
            // no more elements to return
            return S_FALSE;
        }

        CopyToVariant(*rgelt++, Traits::GetValue(*m_current));

        if ( pceltFetched )
            ++(*pceltFetched);
    }

    return S_OK;
}

template <class T, class U>
STDMETHODIMP StdEnum<T, U>::Skip(ULONG celt)
{
    // we might have a template function specialized for random access
    // iterators...
    for ( const IteratorType end = m_container.end(); celt--; ++m_current )
    {
        if ( m_current == end )
        {
            // no more elements to skip
            return S_FALSE;
        }
    }

    return S_OK;
}

template <class T, class U>
STDMETHODIMP StdEnum<T, U>::Reset(void)
{
    Init();

    return S_OK;
}

template <class T, class U>
STDMETHODIMP StdEnum<T, U>::Clone(IEnumVARIANT **ppEnum)
{
    if ( !ppEnum )
        return E_POINTER;

    try
    {
        StdEnum<T, U> *stdEnum = new StdEnum<T, U>(m_container);
        stdEnum->m_current = m_current;
        *ppEnum = stdEnum;

        return S_OK;
    }
    catch ( std::bad_alloc& )
    {
        return E_OUTOFMEMORY;
    }
    catch ( ... )
    {
        return E_UNEXPECTED;
    }
}

} // namespace COM

} // namespace ito33

#endif // _ITO33_COM_ENUM_H_

