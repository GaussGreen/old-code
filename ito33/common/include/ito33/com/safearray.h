/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/safearray.h
// Purpose:     COM::SafeArray template class encapsulates OLE SAFEARRAY
// Author:      Vadim Zeitlin
// Created:     13.01.03
// RCS-ID:      $Id: safearray.h,v 1.19 2004/10/05 09:13:36 pedro Exp $
// Copyright:   (c) 2002,2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/com/safearray.h
    @brief A safer and more convenient to use wrappear for SAFEARRAY.

    SAFEARRAY is the struct used in OLE Automation to pass arrays of the
    elements of any types and arbitrary number of dimensions thru COM. The
    SafeArray<> template class provides a convenient and safer wrapper for
    this abomination.
 */

#ifndef _ITO33_COM_SAFEARRAY_H_
#define _ITO33_COM_SAFEARRAY_H_

#include "ito33/debug.h"
#include "ito33/error.h"

#include "ito33/array.h"
#include "ito33/gettext.h"
#include "ito33/com/exception.h"
#include "ito33/com/variant.h"

namespace ito33
{

namespace COM
{

/**
    SafeArray is a convenient wrapper for SAFEARRAY.

    Note that this class doesn't follow the usual ownership pattern: it never
    takes the ownership of SAFEARRAY, even if it creates a new one itself. I.e.
    this class will never call SafeArrayDestroy() itself. Such behaviour may be
    unusual but useful for the COM servers because they are never going to need
    the SAFEARRAYs themselves: they either create a new one and return it to
    the client (which will destroy it) or fill in the existing one passed to
    them from the client (and which must not be destroyed in the server). For
    the more usual ownership pattern, use SafeArrayPtr.

    This class doesn't provide any way to access the array contents, you should
    use SafeArrayAccessor for this.

    The template parameter T must be of the type supported by SAFEARRAY and, if
    SafeArray is used for working with an existing SAFEARRAY object, even of
    exactly the same type. If this is not the case, an exception will be thrown
    from the ctor.
 */
template <typename T>
class SafeArray
{
public:
  /// the type of the contained elements
  typedef T ElementType;

  /**
      @name Constructors

      Note that this class has no (user-defined) destructor and hence no
      virtual destructor neither and so shouldn't be used polymorphically.
   */
  //@{

  /**
      This constructor creates a new one dimensional array.

      Note that for this function to work, VarType<T> specialization must be
      defined (this is done in variant.h for all standard primitive types).

      The array lower bound will be 0 and upper bound will be count - 1.

      @param count the (fixed) number of array elements
   */
  SafeArray(size_t count);

  /**
      Associates this object with an existing SAFEARRAY.

      Note that we do @b not take ownership of the SAFEARRAY and will not
      destroy it.

      The ctor will throw an exception if pSafeArray is @c NULL or the type
      of its elements is different from ElementType.

      @param pSafeArray pointer to OLE SAFEARRAY, must not be @c NULL
   */
  SafeArray(SAFEARRAY *pSafeArray);

  //@}

  /// @name Operations
  //@{

  /**
      Locks the array for faster access to its elements.

      This shouldn't be called directly, use SafeArrayAccessor instead
      which prevents you from forgetting to call Unlock() (which may lead to
      very hard to debug problems later).

      Note that while the array is locked, it can't be modified.

      Finally, also note that the returned pointer must not be deleted nor
      freed.

      @return the pointer to array data which may be used as a C array
   */
  T *Lock();

  /**
      Unlocks the array allowing modifying it again.

      Unlock() must be called exactly the same number of times as Lock(). Use
      SafeArrayAccessor to ensure that this happens automatically.
   */
  void Unlock();

  //@}

  /// @name Accessors
  //@{

  /**
      Returns the number of the arrays dimensions.

      Note that the safe arrays may be one or multi dimensional and have
      different (and different from 0) bounds for each direction.
   */
  size_t GetDim() const;

  /**
      Returns the number of elements in the array for the given dimension.

      The dimensions are always counted from 0 which corresponds to the first
      dimension.
   */
  size_t GetCount(size_t dim = 0) const;

  /**
      Returns the pointer to the real SAFEARRAY associated with this object.

      This function should only be used to return a SAFEARRAY from some
      function, SAFEARRAY must not be modified directly via this pointer.

      @return a pointer to SAFEARRAY, never @c NULL
   */
  SAFEARRAY *Get() const;

  /**
      Implicit conversion to SAFEARRAY.

      As all implicit conversions, this one is dangerous but its convenience
      outweighs the dangers because it makes it possible to pass SafeArray
      directly to all functions taking a SAFEARRAY which saves us many, many
      calls to Get().

      @sa Get()
   */
  operator SAFEARRAY **();

  //@}

private:
  // the OLE SAFEARRAY
  SAFEARRAY *m_pSafeArray;

  // SafeArray can't be copied
  SafeArray(const SafeArray&);
  SafeArray& operator=(const SafeArray&);
};


/**
    SafeArrayPtr is a simple SAFEARRAY smart pointer.

    Unlike SafeArray itself, it will destroy the array given to it in the ctor.
 */
class SafeArrayPtr
{
public:
  /**
      Constructor takes ownership of an existing SAFEARRAY.

      @param pSA the SAFEARRAY which will be destroyed in our dtor
                 (must not be @c NULL)
   */
  SafeArrayPtr(SAFEARRAY *pSA) : m_pSA(pSA) { }

  /// Destructor destroys the safe array
  ~SafeArrayPtr()
  {
    HRESULT hr = ::SafeArrayDestroy(m_pSA);
    if ( FAILED(hr) )
    {
#ifndef _NDEBUG
      // we can't throw from a dtor so don't let the exception propagate but
      // at least still log it
      std::string
        msg = COM_EXCEPTION("SafeArrayDestroy", hr).GetErrorMessage();
      ::OutputDebugString((msg + "\r\n").c_str());
#endif // _NDEBUG
    }
  }

  /// Implicit conversion to SAFEARRAY **
  operator SAFEARRAY **() { return &m_pSA; }

private:
  SAFEARRAY *m_pSA;

  SafeArrayPtr(const SafeArrayPtr&);
  SafeArrayPtr& operator=(const SafeArrayPtr&);
};


/**
    SafeArrayAccessor provides an efficient way to access SAFEARRAY.

    To access SAFEARRAY data efficiently, it must be locked while this is done.
    This class locks the provided SAFEARRAY in its ctor and unlocks it in the
    dtor ensuring that no "lock leaks" occur.

    Note that SafeArrayAccessor only works with one dimensional arrays and will
    throw an exception if the array has more dimensions.

    @todo: data access for multi dimensional arrays
 */
template <typename T>
class SafeArrayAccessor
{
public:
  /// the type of the elements of the array we access
  typedef T ElementType;

  /**
      Ctor calls SafeArrayAccessData locking the array for data access.

      While the array is locked, its data cannot be modified, so don't leave
      SafeArrayAccessor objects around. Normal usage pattern is to create
      such objects on the stack and let them go out of scope as soon as
      possible.

      An exception is thrown if an error occurs.

      @param sa the array to access
   */
  SafeArrayAccessor(SafeArray<T>& sa);

  /**
      Ctor allowing to access a raw SAFEARRAY.

      Having this ctor allows to directly access SAFEARRAY without having to
      create a SafeArray wrapper object just for this. In fact, what this
      method does is to create such wrapper internally so it's not more
      efficient than doing it in the user code and then using the ctor above,
      but it definitely is simpler to use.

      An exception will be thrown if an error occurs.

      @param pSafeArray pointer to raw SAFEARRAY, can't be @c NULL
   */
  SafeArrayAccessor(SAFEARRAY *pSafeArray);

  /**
      Dtor calls SafeArrayUnaccessData() unlocking the array.

      Any pointers to the array data elements which could have been stored
      are invalidated when the SafeArray object is destroyed (so it's a good
      idea to not keep such pointers at all).

      Also note that the destructor is not virtual and hence this class is
      not supposed to be used polymorphically.
   */
  ~SafeArrayAccessor();

  /**
      @name Miscellaneous accessors

      For convenience, we provide some wrappers for SafeArray<> functions.
   */
  //@{

  /**
      Returns the number of the elements in the array.

      Note that this only works for one dimensional arrays.
   */
  size_t GetCount() const;

  //@}

  /**
      @name Data access

      Array elements can be accessed by index, as with a normal C array. No
      index checks are done in release build.
   */
  //@{

  /**
      Overloaded operator[] allows for an array-like access.

      Note that this only works for one dimensional arrays.
   */
  T& operator[](size_t n);

  /**
      Overloaded operator[] allows for an array-like access.

      Note that this only works for one dimensional arrays.
   */
  const T& operator[](size_t n) const;

  //@}

private:
  // init m_values -- m_array must be already initialized
  void InitValues();


  // the array we're working with
  SafeArray<T> *m_array;

  // the pointer to the values of the array
  T *m_values;

  // if true, we delete m_array when we're done with it
  bool m_ownsArray;


  // semantics of copying SafeArrayAccessors is not very clear so, in doubt,
  // disable it
  SafeArrayAccessor(const SafeArrayAccessor&);
  SafeArrayAccessor& operator=(const SafeArrayAccessor&);
};

/**
    Converts a C array to a SAFEARRAY.

    This function may throw if a COM error occurs.

    @param size the number of elements in the array data
    @param data the array containing size elements
 */
template <typename T>
inline SAFEARRAY *ConvertToSafeArray(size_t size, const T *data)
{
  SafeArray<T> array(size);
  SafeArrayAccessor<T> arrayData(array);
  for ( size_t n = 0; n < size; ++n )
  {
    arrayData[n] = data[n];
  }

  return array.Get();
}

/**
    Converts a CountedArray to SAFEARRAY.

    This function may throw if a COM error occurs.

    @param array CountedArray object
    @return new SAFEARRAY containing all array elements
 */
template <typename T>
inline SAFEARRAY *ConvertToSafeArray(const CountedArray<T>& array)
{
  return ConvertToSafeArray(array.GetCount(), array.Get());
}

/**
    Converts a SAFEARRAY to our CountedArray.

    May only throw std::bad_alloc.

    @param sa accessor for the safe array (used to infer type automatically)
    @return new CountedArray object containing the safe array data
 */
template <typename T>
static CountedArray<T> ConvertFromSafeArray(COM::SafeArrayAccessor<T>& sa)
{
  const size_t count = sa.GetCount();

  CountedArray<T> a(count);
  for ( size_t n = 0; n < count; ++n )
  {
    a[n] = sa[n];
  }

  return a;
}

/**
    Creates a SAFEARRAY from CountedArray and checks for errors.

    @param ppSA pointer to the SAFEARRAY to be created
    @param array the data to be copied to ppSA
    @return the COM error code or S_OK if ok
 */
template <typename T>
inline HRESULT CreateSafeArray(SAFEARRAY **ppSA, const CountedArray<T>& array)
{
  try
  {
    *ppSA = ConvertToSafeArray(array);

    return S_OK;
  }
  catch ( COM::Exception& e )
  {
    return e.GetHresult();
  }
  catch ( ... )
  {
    return E_UNEXPECTED;
  }
}

// ============================================================================
// implementation only from now on
// ============================================================================

// ----------------------------------------------------------------------------
// SafeArray
// ----------------------------------------------------------------------------

template <typename T>
inline
SafeArray<T>::SafeArray(SAFEARRAY *pSafeArray)
{
  // check that the array is not NULL and has valid type
  VARTYPE vt;
  if ( !pSafeArray ||
        FAILED(SafeArrayGetVartype(pSafeArray, &vt)) ||
            !IsCompatible<T>::With(vt) )
  {
    // this is not a COM exception, really, but it makes sense to emulate
    // one here
    throw COM_EXCEPTION("SafeArrayGetVartype", E_INVALIDARG);
  }

  m_pSafeArray = pSafeArray;
}

template <typename T>
inline
SafeArray<T>::SafeArray(size_t count)
{
  // we must use the "Ex" version because otherwise the array is created
  // without FADF_HAVEVARTYPE flag and so creating SafeArray from
  // m_pSafeArray pointer later would fail (because SafeArrayGetVartype()
  // wouldn't work)
  m_pSafeArray = ::SafeArrayCreateVectorEx
                   (
                      VarType<T>::Value,
                      0,
                      static_cast<ULONG>(count),
                      NULL
                   );
  if ( !m_pSafeArray )
  {
    throw COM_EXCEPTION("SafeArrayCreateVectorEx", E_UNEXPECTED);
  }
}

template <typename T>
T *SafeArray<T>::Lock()
{
  T *data;
  HRESULT hr = ::SafeArrayAccessData(m_pSafeArray,
                    reinterpret_cast<void **>(&data));

  if ( FAILED(hr) )
    throw COM_EXCEPTION("SafeArrayAccessData", hr);

  return data;
}

template <typename T>
void SafeArray<T>::Unlock()
{
  HRESULT hr = ::SafeArrayUnaccessData(m_pSafeArray);

  if ( FAILED(hr) )
    throw COM_EXCEPTION("SafeArrayUnaccessData", hr);
}

template <typename T>
inline
SAFEARRAY *SafeArray<T>::Get() const
{
  return m_pSafeArray;
}

template <typename T>
inline
SafeArray<T>::operator SAFEARRAY **()
{
  return &m_pSafeArray;
}

template <typename T>
inline
size_t SafeArray<T>::GetDim() const
{
  // this should be the same but faster than ::SafeArrayGetDim(m_pSafeArray)
  return m_pSafeArray->cDims;
}

template <typename T>
inline
size_t SafeArray<T>::GetCount(size_t dim) const
{
  ASSERT_MSG( dim < GetDim(), "invalid dimension in SafeArray::GetCount" );

  return m_pSafeArray->rgsabound[dim].cElements;
}

// ----------------------------------------------------------------------------
// SafeArrayAccessor
// ----------------------------------------------------------------------------

template <typename T>
inline
void SafeArrayAccessor<T>::InitValues()
{
  m_values = m_array->Lock();
}

template <typename T>
inline
SafeArrayAccessor<T>::SafeArrayAccessor(SafeArray<T>& sa)
{
  m_array = &sa;

  m_ownsArray = false;

  InitValues();
}

template <typename T>
inline
SafeArrayAccessor<T>::SafeArrayAccessor(SAFEARRAY *pSafeArray)
{
  m_array = new SafeArray<T>(pSafeArray);

  if ( m_array->GetDim() != 1 )
  {
    delete m_array;

    throw ito33::EXCEPTION_MSG
       (
        ITO33_BAD_PARAM,
        TRANS("Array must be one dimensional.")
       );
  }

  m_ownsArray = true;

  InitValues();
}

template <typename T>
inline
SafeArrayAccessor<T>::~SafeArrayAccessor()
{
  m_array->Unlock();

  if ( m_ownsArray )
    delete m_array;
}

template <typename T>
inline
size_t SafeArrayAccessor<T>::GetCount() const
{
  return m_array->GetCount();
}

template <typename T>
inline
T& SafeArrayAccessor<T>::operator[](size_t n)
{
  ASSERT_MSG( n < m_array->GetCount(), "invalid SafeArray index" );

  return m_values[n];
}

template <typename T>
inline
const T& SafeArrayAccessor<T>::operator[](size_t n) const
{
  ASSERT_MSG( n < m_array->GetCount(), "invalid SafeArray index" );

  return m_values[n];
}


} // namespace COM

} // namespace ito33

#endif // _ITO33_COM_SAFEARRAY_H_


