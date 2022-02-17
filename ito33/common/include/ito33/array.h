/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/array.h
// Purpose:     An auto_ptr<>-like classes for the arrays
// Author:      Vadim Zeitlin
// Created:     24.01.03
// RCS-ID:      $Id: array.h,v 1.22 2005/02/21 14:29:06 nabil Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/array.h
    @brief  Declares the template Array class.

    The Array class should be used for the arrays which never change size as in
    this case it is (slightly) more efficient than vector. CountedArray is a
    minor variation on the same theme -- the only difference is that it also
    stores the number of elements in the array.
 */

#ifndef _ITO33_ARRAY_H_
#define _ITO33_ARRAY_H_

#include "ito33/beforestd.h"
#include "ito33/afterstd.h"

#include "ito33/debug.h"

// disable the intel compiler remark about non virtual dtor for base class
#ifdef __INTEL_COMPILER  
#pragma warning(disable:444)
#endif

namespace ito33
{

/**
    Simple template container class behaving like a C array.

    The differences between Array<> and a C array are that Array also stores
    its size and that it automatically frees its memory in dtor. So, in fact,
    this class more or less emulates VLAs from the latest C standard.

    This array cannot be extended so its size must be specified either during
    its creation or by assigning another array to it later. Its copy semantics
    are the same as for auto_ptr<>, in other words it does destructive copy.
 */
template <typename T>
class Array
{
public:
  /// the type of the array elements
  typedef T ElementType;

  /// @name Constructors, assignment operator and destructor.
  //@{

  /**
      Default ctor.

      The array is created empty.
   */
  Array();

  /**
      Ctor takes ownership of the pointer, ie we will free it later.

      The pointer must be either @c NULL or allocated with @c new[].

      @param data pointer to the array elements
   */
  Array(T *data);

  /**
      Ctor allocates memory for count elements of the given type.

      The array elements have undefined value, the caller is responsible for
      initializing them.

      @param count number of elements in the array
   */
  explicit Array(size_t count);

  /**
      Copy ctor takes ownership of the other objects data.

      As with std::auto_ptr<>, the other object is reset to the initial,
      uninitialized, state.
   */
  Array(const Array& other);

  /**
      Copy ctor takes ownership of the other objects data.

      As with std::auto_ptr<>, the other object is reset to the initial,
      uninitialized, state.
   */
  Array& operator=(const Array& other);

  /**
      Dtor frees the memory taken by the array elements using @c delete[]
   */
  virtual ~Array();

  //@}


  /// @name Accessors
  //@{

  /// returns the n-th element of the array (readonly)
  const T& operator[](size_t n) const;

  /// returns the n-th element of the array
  T& operator[](size_t n);

  /// get the pointer to data
  T *Get();

  /// get the const pointer to data
  const T *Get() const;

  //@}


  /// @name Standard functions
  //@{

  /// releases ownership of the held pointer and returns it
  T* Release();

  /// Swap the underlying data pointer with another array
  virtual void Swap(Array& other);

  //@}


protected:
  /// allocates space for count elements (may throw std::bad_alloc)
  T *Alloc(size_t count);

  /// delete m_array, doesn't reset it to NULL
  void Free();

  /// the array data (not private because also accessed by CountedArray)
  T *m_array;
};

/**
    CountedArray is the same as a simple Array except that it also stores the
    number of its elements.

    In debug mode, the indices are checked in operator[] but in release mode no
    checks are done, so don't count on it.
 */
template <typename T>
class CountedArray : public Array<T>
{
public:
  /**
      Ctor takes ownership of the pointer, ie we will free it later.

      The pointer must be either @c NULL or allocated with @c new[].

      @param count the number of elements in the data array
      @param data pointer to the array elements
   */
  CountedArray(size_t count, T *data);

  /**
      Ctor creates an uninitialized array of the given size.

      The array elements have undefined value, the caller is responsible for
      initializing them.

      @param count the number of elements in the array to create
   */
  explicit CountedArray(size_t count = 0);


  /**
      Copy the data from C array to this object.

      Previous contents of the array is lost.
   */
  void Copy(size_t count, T *pData);


  /// Swap the underlying data pointer with another array
  virtual void Swap(CountedArray& other);


  /// @name Accessors
  //@{

  /// returns the number of elements in the array
  size_t GetCount() const;

  /// returns the n-th element of the array (readonly)
  const T& operator[](size_t n) const;

  /// returns the n-th element of the array
  T& operator[](size_t n);

  //@}

private:
  // the number of elements in the array
  size_t m_count;

  // bring non-dependent stuff used by our methods in this scope to fix
  // compilation with compilers properly implementing 2 phase lookup (such as
  // g++ 3.4)
  using Array<T>::Alloc;
  using Array<T>::Free;

  using Array<T>::m_array;
};

// ============================================================================
// inline methods implementation
// ============================================================================

// ----------------------------------------------------------------------------
// Array
// ----------------------------------------------------------------------------

template <typename T>
inline
T *Array<T>::Alloc(size_t count)
{
  return new T[count];
}

template <typename T>
inline
void Array<T>::Free()
{
  delete [] m_array;
}

template <typename T>
inline
Array<T>::Array()
{
  m_array = NULL;
}

template <typename T>
inline
Array<T>::Array(T *data)
{
  m_array = data;
}

template <typename T>
inline
Array<T>::Array(size_t count)
{
  m_array = count ? Alloc(count) : NULL;
}

template <typename T>
inline
Array<T>::Array(const Array<T>& other)
{
  m_array = other.m_array;

  const_cast<Array<T>&>(other).m_array = NULL;
}

template <typename T>
inline
Array<T>& Array<T>::operator=(const Array<T>& other)
{
  Free();

  m_array = other.m_array;

  const_cast<Array<T>&>(other).m_array = NULL;

  return *this;
}

template <typename T>
inline
T& Array<T>::operator[](size_t n)
{
  return m_array[n];
}

template <typename T>
inline
const T& Array<T>::operator[](size_t n) const
{
  return m_array[n];
}

template <typename T>
inline
T *Array<T>::Get()
{
  return m_array;
}

template <typename T>
inline
const T *Array<T>::Get() const
{
  return m_array;
}

template<typename T>
inline
T* Array<T>::Release()
{
  T* p = m_array;

  m_array = 0;

  return p;
}

template <typename T>
inline
void Array<T>::Swap(Array& other)
{
  T* pTmp = other.m_array;
  other.m_array = m_array;
  m_array = pTmp;
}

template <typename T>
inline
Array<T>::~Array()
{
  Free();
}



// ----------------------------------------------------------------------------
// CountedArray
// ----------------------------------------------------------------------------

template <typename T>
inline
CountedArray<T>::CountedArray(size_t count, T *data)
       : Array<T>(data)
{
  m_count = count;
}

template <typename T>
inline
CountedArray<T>::CountedArray(size_t count)
       : Array<T>(count)
{
  m_count = count;
}

template <typename T>
inline
void CountedArray<T>::Copy(size_t count, T *pData)
{
  Free();

  // don't complain if the pointer is NULL if we don't have to use it anyhow
  m_count = count;
  if ( !count )
  {
    m_array = NULL;
  }
  else // do copy data
  {
    if ( !pData )
    {
      FAIL( "NULL data pointer in CountedArray<>::Copy()" );

      m_array = NULL;
    }
    else
    {
      m_array = Alloc(count);
      memcpy(m_array, pData, sizeof(T)*count);
    }
  }
}

template <typename T>
inline
size_t CountedArray<T>::GetCount() const
{
  return m_count;
}

template <typename T>
inline
void CountedArray<T>::Swap(CountedArray& other)
{
  T* pTmp = other.m_array;
  other.m_array = m_array;
  m_array = pTmp;

  size_t nCountTmp = other.m_count;
  other.m_count = m_count;
  m_count = nCountTmp;
}

template <typename T>
inline
T& CountedArray<T>::operator[](size_t n)
{
  ASSERT_MSG( n < m_count, "index out of range in CountedArray<>" );

  return m_array[n];
}

template <typename T>
inline
const T& CountedArray<T>::operator[](size_t n) const
{
  ASSERT_MSG( n < m_count, "index out of range in CountedArray<>" );

  return m_array[n];
}

// ============================================================================
// Specializations of generic functions
// ============================================================================

// ----------------------------------------------------------------------------
// Array
// ----------------------------------------------------------------------------

template <typename T>
void swap(Array<T>& a, Array<T>& b)
{
  a.Swap(b);
}

// ----------------------------------------------------------------------------
// CountedArray
// ----------------------------------------------------------------------------

template <typename T>
void swap(CountedArray<T>& a, CountedArray<T>& b)
{
  a.Swap(b);
}

// restore the intel compiler remark about non virtual dtor for base class
#ifdef __INTEL_COMPILER  
#pragma warning(default:444)
#endif

} // namespace ito33

#endif // _ITO33_ARRAY_H_

