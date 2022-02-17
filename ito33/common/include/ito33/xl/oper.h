/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/XL/oper.h
// Purpose:     XL::Oper class wrapping XLOPER struct from Excel SDK
// Author:      Vadim Zeitlin
// Created:     2006-04-14
// RCS-ID:      $Id: oper.h,v 1.20 2006/08/21 15:05:12 willy Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/XL/oper.h
    @brief Convenient wrapper class around XLOPER struct from Excel SDK.

    This header provides classes and functions wrapping Excel C API. XL::Oper
    replaces XLOPER struct and XL::Call() is a safe version of Excel4()
    function.

    MT-safety note: currently these classes are not MT-safe so all Excel
    interaction must be performed by a single thread.
 */

#ifndef _ITO33_XL_OPER_H_
#define _ITO33_XL_OPER_H_

#include "ito33/date.h"

#include <set>
#include <boost/preprocessor/repetition.hpp>

#include "ito33/win32/winwrap.h"    // must be included before xlcall.h

#define bool boolean
#include <xlcall.h>
#undef bool

#include "ito33/XL/exception.h"

/**
    Maximal number of columns in a table supported by XL::Oper::AsTable().
 */
#ifndef XL_MAX_TABLE_COLUMNS
  #define XL_MAX_TABLE_COLUMNS 10
#endif

namespace ito33
{

namespace XL
{

#ifdef DOXYGEN

/**
    Type-traits template allows to work with different types in the same way.
 */
template <typename T>
struct Traits
{
  /**
      The Excel type corresponding to this type.
   */
  enum { Type = xltypeNil };

  /**
      The Excel signature for the values of this type.
   */
  enum { Sig = 'Z' };

  /**
      Extract the value of this type from the given XLOPER (which is supposed
      to be of correct type).
   */
  static T Extract(const XLOPER *xlo);

  /**
      Returns true if this XLOPER value is equal to the given one.

      Uses simple comparison for the primitive types and deep comparison (e.g.
      strcmp) for the other ones.
   */
  static bool IsSameAs(const XLOPER *xlo, T value);
};

#else // DOXYGEN

template <typename T> struct Traits;

// Excel doesn't make any distinction between const and non-const pointers
template <typename T>
struct Traits<const T *> : Traits<T *>
{
};

template <typename T>
struct ShallowCompare
{
  static bool IsSameAs(const XLOPER *xlo, T value)
  {
    return Traits<T>::Extract(xlo) == value;
  }
};

template <>
struct Traits<bool> : ShallowCompare<bool>
{
  enum { Type = xltypeBool };
  enum { Sig = 'A' };

  static bool Extract(const XLOPER *xlo) { return xlo->val.boolean != 0; }
};

struct NumTraits : ShallowCompare<double>
{
  enum { Type = xltypeNum };
  enum { Sig = 'B' };
};

template <>
struct Traits<double> : NumTraits
{
  static double Extract(const XLOPER *xlo) { return xlo->val.num; }
};

// dates are represented as "serials", i.e. doubles, in Excel
template <>
struct Traits<Date> : NumTraits
{
  static Date Extract(const XLOPER *xlo)
  {
    return Date(static_cast<unsigned long>(xlo->val.num));
  }
};

struct StrTraits
{
  enum { Type = xltypeStr };
  enum { Sig = 'C' }; // 'F' is a synonym

  static bool IsSameAs(const XLOPER *xlo, const char *value)
  {
    // first byte is the length, skip it
    return xlo->val.str && strncmp(xlo->val.str + 1, value, *xlo->val.str) == 0;
  }
};

template <>
struct Traits<char *> : StrTraits
{
};

template <>
struct Traits<std::string> : StrTraits
{
  static std::string Extract(const XLOPER *xlo)
  {
    return xlo->val.str ? std::string(xlo->val.str + 1,
                                      static_cast<unsigned char>(*xlo->val.str))
                        : std::string();
  }

};

template <>
struct Traits<unsigned char *> : StrTraits
{
  enum { Sig = 'D' }; // 'G' is a synonym
};

struct IntTraits : ShallowCompare<int>
{
  enum { Type = xltypeInt };

  static int Extract(const XLOPER *xlo) { return xlo->val.w; }
};

template <>
struct Traits<unsigned short> : IntTraits
{
  enum { Sig = 'H' };
};

template <>
struct Traits<short> : IntTraits
{
  enum { Sig = 'I' };
};

template <>
struct Traits<int> : IntTraits
{
  enum { Sig = 'J' };
};

template <>
struct Traits<long> : Traits<int>
{
  enum { Sig = 'J' };
};

// notice the difference between using the parameters of type "XLOPER *" and
// "XL::Oper *" in the function signature: the former is just the general
// XLOPER which can contain references to the other cells while the latter is
// always a dereferenced value (what some XL documentation calls struct OPER)
// and this is generally what you need, unless you really have to know whether
// the function argument comes directly from the formula or from another cell
class Oper;

template <>
struct Traits<Oper *>
{
  enum { Sig = 'P' };   // corresponds to OPER, i.e. dereferenced XLOPER
};

template <>
struct Traits<XLOPER *>
{
  enum { Sig = 'R' };
};

#endif // DOXYGEN/!DOXYGEN

/**
    Allows to pass Oper objects to Excel4().

    This class is mostly used internally but is also useful if you have to use
    XL::Call() directly.
 */
class ByRef
{
public:
  /**
      Associates this object with the given XL::Oper.
   */
  ByRef(Oper& oper) : m_oper(oper)
  {
  }

  /**
      Destructor marks the object as being filled by Excel.

      This ensures that we don't try to free memory allocated by Excel and is
      the reason for this class existence (instead of just allowing casting
      Oper to non-const XLOPER).
   */
  ~ByRef();

  /**
      Implicit conversion to non-const XLOPER.

      This conversion allows to implicitly pass ByRef() objects to XL::Call()
      while the absence of such conversion in XL::Oper itself ensures that we
      can't pass it to Call() directly.
   */
  operator XLOPER *() const;

private:
  Oper& m_oper;

  NO_COPY_CLASS(ByRef);
};

/**
      Define a convenient Ref class corresponds to XLMREF limited
      to 1 range only.
      
      This class must be the parameter of Oper constructor if one wants to 
      create an Oper of xltypeRef type.
*/

class Ref
{
public:

  /**
      Default Constructor.
      
      No reference is set.    
   */
  Ref():m_SheetId(0),m_FirstRow(0),m_LastRow(0),m_FirstCol(0),m_LastCol(0),hasRefSet( false)
  {
  }

  /**
    Constructor from a reference to a range of cells

    Have default argument
  */
  Ref(unsigned rwfirst,unsigned colfirst,
      unsigned nbrows = 1, unsigned nbcols = 1,
      unsigned long idSheet =0)
  {
    m_SheetId = idSheet;
    InitRange(rwfirst,colfirst,nbrows,nbcols);    
  }

  //--------------------------------------------------------------------------
  // Getter functions
  //--------------------------------------------------------------------------
  
  /**
      @return return the sheet id of the reference
  */
  unsigned long GetSheetId() const {return m_SheetId;}
  
  /**
      It must be called to know if GetFirstRow(),GetLastCol(),GetLastRow() and GetfirstColl()
      methods return garbage or a real range address.
            
      @return true if a reference is set
  */ 
  bool HasReference() const {return hasRefSet;}

  /**
      @return the first row of the range
  */
  unsigned short GetFirstRow() const {return m_FirstRow;}
  
  /**
      @return the last row of the range
  */
  unsigned short GetLastRow() const {return m_LastRow;}
  
  /**
      @return the first col of the range
  */
  unsigned char GetFirstCol() const {return m_FirstCol;}
  
  /**
      @return the last col of the range
  */
  unsigned char GetLastCol() const {return m_LastCol;}

  //-------------------------------------------------------------------------
  // Setter functions
  //-------------------------------------------------------------------------
  
  /**
    Set the sheet id

    @param idSheet the id of the sheet as given by Excel API
  */
  void SetSheetId(unsigned long idSheet)
  {
    m_SheetId = idSheet;
  }

  /**
    Set a range
    
    @param rwfirst indicates the first row
    @param colfirst indicates the first column
    @param nbrows indicates the number of rows conpounding the range (defalt is 1)
    @param nbrows indicates the number of columns conpounding the range (defalt is 1)
  */
  void SetRange(unsigned rwfirst,unsigned colfirst,
                unsigned nbrows = 1, unsigned nbcols = 1)
  {
    InitRange(rwfirst,colfirst,nbrows,nbcols);
  }
  
  /**
      the range reference become not valid     
  */
  void Clear(){ hasRefSet = false;}
private:
  
  // functions call by either the SetRange() function or the constructor
  void InitRange( unsigned rwfirst,unsigned colfirst,
                  unsigned nbrows = 1, unsigned nbcols = 1)
  {
    if ( !nbrows || !nbcols )
    {
      throw EXCEPTION_MSG(xlerrValue, "Attempt to create an empty range");
    }

    const unsigned rwlast = rwfirst + nbrows - 1,
                   collast = colfirst + nbcols - 1;

    if ( collast > 0xff || rwlast > 0xffff )
    {
      throw EXCEPTION_MSG(xlerrValue, "Cell coorindates out of range");
    }
    m_FirstRow = static_cast<WORD>(rwfirst);
    m_LastRow = static_cast<WORD>(rwlast);
    m_FirstCol = static_cast<BYTE>(colfirst);
    m_LastCol = static_cast<BYTE>(collast);

    hasRefSet = true;
  }

  unsigned long m_SheetId;
  unsigned short m_FirstRow;
  unsigned short m_LastRow;
  unsigned char m_FirstCol;
  unsigned char m_LastCol;

  // To distinguish if a range reference has been set or not 
  bool hasRefSet;
};

/**
    Convenient wrapper around XLOPER struct which is used to exchange the
    information between Excel and add-ins.

    Note that this class is a POD and has exactly the same layout as XLOPER
    struct which means that we can interpret pointers to XLOPER we get from
    Excel as pointers to the objects of this class and that we cannot add any
    members to it.

    Memory management: first of all, this only applies to operands of type
    xltypeStr, xltypeRef or xltypeMulti. Second, there are several different
    cases here:
      - Oper created in our code and destroyed in our code (in the meanwhile it
        was presumably passed as an argument to some XL function): in this case
        we allocate and free memory ourselves
      - Oper created in our code and returned to XL: you have to call Detach()
        to indicate that it shouldn't be freed or use static Return() method
      - Oper returned by XL: they are also freed by Excel in Clear()
 */

class Oper : private XLOPER
{
public:
  /**
      Typed synonyms for xltypeXXX constants.
   */
  enum Type
  {
    Type_Num      = xltypeNum,
    Type_Str      = xltypeStr,
    Type_Bool     = xltypeBool,
    Type_Ref      = xltypeRef,
    Type_Err      = xltypeErr,
    Type_Flow     = xltypeFlow,
    Type_Multi    = xltypeMulti,
    Type_Missing  = xltypeMissing,
    Type_Nil      = xltypeNil,
    Type_SRef     = xltypeSRef,
    Type_Int      = xltypeInt
  };

  /**
      Typed synonyms for xlerrXXX constants.
   */
  enum Err
  {
    Err_Null      = xlerrNull,
    Err_Div0      = xlerrDiv0,
    Err_Value     = xlerrValue,
    Err_Ref       = xlerrRef,
    Err_Name      = xlerrName,
    Err_Num       = xlerrNum,
    Err_NA        = xlerrNA
  };


  /**
      @name Static constructors.
   */
  //@{

  /**
      Returns an Oper of missing type, suitable for passing to XL::Call() for a
      parameter we don't specify.
   */
  static const Oper& Missing()
  {
    static const Oper s_operMissing(Type_Missing, static_cast<TypeTag *>(NULL));
    return s_operMissing;
  }

  /**
      Returns a value error object.

      This can be returned directly to Excel.
   */
  static XLOPER *ErrValue()
  {
    static Oper s_operErrValue(Err_Value, static_cast<ErrTag *>(NULL));
    return &s_operErrValue;
  }

  /**
      Creates a new Oper object and returns it to Excel.
   */
  template <typename T>
  static XLOPER *Return(T value)
  {
    Oper *oper = new Oper(value);
    return oper->Detach();
  }

  /**
      Creates a new array Oper object and returns it to Excel.

      @param w number of columns in the array
      @param h number of rows in the array
      @param data pointer to the array data, must be non-NULL and allocated
                  with new[]; we take ownership of it
   */
  static XLOPER *ReturnArray(unsigned w, unsigned h, Oper *data)
  {
    Oper *oper = new Oper(w, h, data);
    return oper->Detach();
  }

  //@}


  /**
      Default ctor initializes the operand to nil value.
   */
  Oper()
  {
    xltype = xltypeNil;
  }

  /**
      Copy constructor does deep copy.
   */
  Oper(const Oper& other)
  {
    Copy(other);
  }

  /**
      Constructor from a C string.

      Copies the string internally.

      May throw if string is too long (strings are limited to 255 characters in
      Excel C API).
   */
  Oper(const char *value)
  {
    InitWithString(value, strlen(value));
  }

  /**
      Constructor from a C++ string.

      Copies the string internally.

      May throw if string is too long (strings are limited to 255 characters in
      Excel C API).
   */
  Oper(const std::string& value)
  {
    InitWithString(value.c_str(), value.length());
  }

  /**
      Constructor from a int.

      Notice that Excel doesn't support ints greater than SHRT_MAX and this
      function throws unless the @a value is in range.
   */
  Oper(int value)
  {
    if ( value > SHRT_MAX || value < SHRT_MIN )
    {
      throw EXCEPTION_MSG(xlerrValue, "Integer number is too big");
    }

    xltype = xltypeInt;
    val.w = static_cast<short>(value);
  }

  /**
      Constructor from a double.
   */
  Oper(double value)
  {
    xltype = xltypeNum;
    val.num = value;
  }

  /**
      Constructor from a date.

      The date is converted into its Excel representation as double.
   */
  Oper(const Date& date)
  {
    xltype = xltypeNum;
    val.num = date.GetExcel();
  }

  /**
      Constructor for an array.

      Creates an object containing Excel array initialized with the given data.

      @param w number of columns in the array
      @param h number of rows in the array
      @param data pointer to the array data, must be non-NULL and allocated
                  with new[]; we take ownership of it
   */
  Oper(unsigned w, unsigned h, Oper *data)
  {
    const unsigned n = w*h;
    if ( !data || !n || w >= 0xff || h >= 0xffff )
    {
      throw EXCEPTION_MSG(xlerrValue, "Invalid array size");
    }

    InitAsArrayNoThrow(static_cast<WORD>(w), static_cast<WORD>(h));

    val.array.lparray = data;
  }

  /**
      Constructor for a simple reference.

      If any of the parameters is out of range, an exception is thrown.

      @param x starting column of the range
      @param y starting row of the range
      @param w number of columns in the range (1 by default, can't be 0)
      @param h number of rows in the range (1 by default, can't be 0)
   */
  Oper(unsigned x, unsigned y, unsigned w = 1, unsigned h = 1)
  {
    if ( !w || !h )
    {
      throw EXCEPTION_MSG(xlerrValue, "Attempt to create an empty range");
    }

    const unsigned x2 = x + w - 1,
                   y2 = y + h - 1;

    if ( x2 > UCHAR_MAX || y2 > USHRT_MAX )
    {
      throw EXCEPTION_MSG(xlerrValue, "Cell coorindates out of range");
    }

    xltype = xltypeSRef;

    val.sref.count = 1; // always 1

    XLREF& ref = val.sref.ref;
    ref.colFirst = static_cast<BYTE>(x);
    ref.rwFirst = static_cast<WORD>(y);
    ref.colLast = static_cast<BYTE>(x2);
    ref.rwLast = static_cast<WORD>(y2);
  }

  /**
      Constructor for a mref.

      @param ref the Ref object      
   */
  Oper(const Ref& ref)
  {
    InitAsRefNoThrow(ref);
  }

  /**
      Assignment operator does deep copy.
   */
  Oper& operator=(const Oper& other);

  /**
      Assignment operator from any of the primitive types supported by the
      constructors above.
   */
  template <typename T>
  Oper& operator=(T value)
  {
    return *this = Oper(value);
  }

  /**
      Destructor frees memory if necessary.
   */
  ~Oper();


  /**
      @name Accessors.
   */
  //@{

  /// Get the type of this object
  Type GetType() const
  {
    return static_cast<Type>(xltype & ~(xlbitXLFree | xlbitDLLFree));
  }

  /// Contains an error?
  bool IsError() const { return GetType() == Type_Err; }

  /// Is this object missing or NULL?
  bool IsMissing() const { return (GetType() & (Type_Missing | Type_Nil)) != 0; }

  /// Same as IsMissing() except it can also be used with NULL pointers
  static bool IsMissing(const Oper *o) { return !o || o->IsMissing(); }
  
  /// Return true if the object is missing, is an empty string or an empty range
  static bool IsEmpty(const Oper *o)
  {
    if ( IsMissing(o) ) return true;
    
    if ( o->GetType() & Type_Str )
    {
      std::string str;
      
      if ( o->As(&str) ) return str == ""?true:false;
      return false; //should never be there
    }
    else if ( o->GetType() & Type_Multi )
    {
      unsigned int w,h;
      o->GetAsArray(w,h);
      return h == 0?true:false;
    }
    return false;
  }

  /// Return the error code, only valid to call if IsError() returns true
  int GetErrorCode() const
  {
    ASSERT_MSG( IsError(), "invalid for non-error XLOPER" );

    return val.err;
  }


  /**
      Coerce the value to the given scalar type if possible.

      @param result if non-NULL, filled with the value converted to target type
                    if the conversion was successful
      @return true if coercion was successful, false otherwise
   */
  template <typename T>
  bool As(T *result) const
  {
    Oper xloRC;
    if ( CallNoThrow(xlCoerce, ByRef(xloRC), *this,
                     Oper(Traits<T>::Type)) != xlretSuccess || xloRC.IsError() )
      return false;

    if ( result )
      *result = Traits<T>::Extract(xloRC);

    return true;
  }

  /**
      Coerce the given object to an array.

      This is a low level method: if possible, use AsTable() instead.

      This should be used to represent ranges and real arrays as arrays
      independently of where do they come from.

      @param o the object to be filled with array data, GetAsArray() can be
               called on it later
      @return true if successful or false if coercion to array failed
   */
  bool AsArray(Oper *o) const;

  /**
      Access the value as array of XL::Oper pointers.

      This succeeds only for objects of type Type_Multi and returns @c NULL
      otherwise, i.e. unlike As() this method doesn't perform any coercions.
      It is a low level method and shouldn't be normally used directly at all,
      use AsTable() instead.

      @param w filled with the number of array columns if successful
      @param h filled with the number of array rows if successful
      @return the pointer to array data if successful or @c NULL
   */
  const Oper *GetAsArray(unsigned& w, unsigned& h) const;

#ifdef DOXYGEN
  /**
      Access the value as an N-column table.

      The value must be of type Type_Multi and have exactly N columns (but
      arbitrary number of rows). The values from the first column are stored in
      @a c1, from the second one -- in @a c2 and so on. If all conversions
      succeed, the function returns true, if any of them fail, it returns false.

      All parameters must provide std::vector-like interface.

      @param c1 receives the values of the first column in the table
      @param c2 receives the values of the second column in the table
      ...
      @param cN receives the values of the last column in the table
      @return true if all values were successfully stored, false otherwise
   */
  template <class C1, class C2, ..., class CN>
  bool AsTable(C1& c1, C2& c2, ..., CN& cN) const;
#else // !DOXYGEN
  #define MAKE_RESIZE(z, n, unused) c##n.resize(h);
  #define MAKE_AS_COLUMN(z, n, unused)                                        \
    if ( !data->As(&c##n[i]) )                                                \
      return false;                                                           \
    ++data;

  #define MAKE_AS_TABLE_n(z, n, unused)                                       \
  template < BOOST_PP_ENUM_PARAMS(n, class C) >                               \
  bool AsTable(BOOST_PP_ENUM_BINARY_PARAMS(n, C, &c)) const                   \
  {                                                                           \
    Oper array;                                                               \
    if ( !AsArray(&array) )                                                   \
      return false;                                                           \
    unsigned w, h;                                                            \
                                                                              \
    const Oper *data = array.GetAsArray(w, h);                                \
    if ( w != n )                                                             \
      return false;                                                           \
                                                                              \
    BOOST_PP_REPEAT(n, MAKE_RESIZE, ~)                                        \
    for ( unsigned i = 0; i < h; ++i )                                        \
    {                                                                         \
      BOOST_PP_REPEAT(n, MAKE_AS_COLUMN, ~)                                   \
    }                                                                         \
                                                                              \
    return true;                                                              \
  }

  BOOST_PP_REPEAT_FROM_TO(1, XL_MAX_TABLE_COLUMNS, MAKE_AS_TABLE_n, ~)

  #undef MAKE_AS_TABLE_n
#endif // DOXYGEN/!DOXYGEN


  /**
      Accesses the value as a reference.

      No coercion is done here, the object must be of type Type_SRef. If not,
      the function returns false (but doesn't throw).

      @param x receives the starting column of the reference
      @param y receives the starting row of the reference
      @param w reference the number of columns in the range if not NULL
      @param h reference the number of rows in the range if not NULL
      @return true if this is a simple reference, false otherwise
   */
  bool GetAsSRef(unsigned& x,
                 unsigned& y,
                 unsigned *w = NULL,
                 unsigned *h = NULL) const
  {
    if ( GetType() != Type_SRef )
      return false;

    const XLREF& ref = val.sref.ref;
    x = ref.colFirst;
    y = ref.rwFirst;
    if ( w )
      *w = ref.colLast - x + 1;
    if ( h )
      *h = ref.rwLast - y + 1;

    return true;
  }

  /**
      Accessees the sheet id of a Type_Ref object.
      The function doesn't throw but it return false if the operation is not valid 
      (the object is not a mref)

      @param idSheet receives the sheet id of the mref
      @return true if this is a mref, false otherwise
  */
  bool GetSheetId(DWORD& idSheet) const
  {
    if ( GetType() != Type_Ref )
      return false;
        
    idSheet = val.mref.idSheet;
    return true;
  }
  
  /**
      Accesses the value as a XLMREF.

      No coercion is done here, the object must be of type Type_Ref. If not,
      the function returns false.
      The function throw only if the object is a Type_Ref object with more than one range reference .

      @param ref receives the MREF Description ( range delimitation and the id of the sheet )
      @return true if this is a mreference, false otherwise
   */
  bool GetAsRef(Ref& ref) const
  {
    if ( GetType() != Type_Ref )
      return false;
    ref.Clear();
    ref.SetSheetId( val.mref.idSheet );

    if ( val.mref.lpmref )
    {
      const XLMREF& xlmref = *val.mref.lpmref;
      if ( xlmref.count != 1 )
        throw EXCEPTION_MSG(xlerrValue, "multiple references are not supported");

      ref.SetRange( xlmref.reftbl[0].rwFirst,
                    xlmref.reftbl[0].colFirst,
                    xlmref.reftbl[0].rwLast - xlmref.reftbl[0].rwFirst + 1,
                    xlmref.reftbl[0].colLast - xlmref.reftbl[0].colFirst + 1);
    }
    return true;
  }

  /**
      Generic comparison operator.

      This works for doubles, ints and booleans only currently.
   */
  template <typename T>
  bool operator==(T value) const
  {
    Oper xloRC;
    Call(xlCoerce, ByRef(xloRC), *this, Oper(Traits<T>::Type));
    if ( xloRC.IsError() )
      return false;

    return Traits<T>::IsSameAs(xloRC, value);
  }

  /// Comparison operator
  template <typename T>
  bool operator!=(T value) const { return !(*this == value); }

  //@}


  /**
      @name Passing to/from Excel.
   */
  //@{

  /**
      Explicit conversion to XLOPER.

      This only works for passing Oper by value to Call(), use ByRef for the
      return value and read/write parameters.
   */
  const XLOPER *ByVal() const { return this; }

  /**
      Implicit conversion to XLOPER.
   */
  operator const XLOPER *() const { return ByVal(); }

  /**
      Should be called when returning the object as a function result to Excel.

      This method sets the flag telling Excel to call us back to free this
      operand and also return the underlying XLOPER to be returned to Excel.
      Notice that the object itself should be allocated on the heap as it will
      be deleted later.
   */
  XLOPER *Detach() { xltype |= xlbitDLLFree; return this; }

  //@}

  /**
      Returns the symbolic name of the given error.
   */
  static std::string GetErrorName(WORD err);

  /**
      Debug helper: returns the string representation of this object.
   */
  std::string Dump() const;

private:
  // helper managing the global list of Oper objects allocated by Excel
  //
  // this helper class is not MT-safe because s_list is unprotected
  class XLAllocatedList
  {
  public:
    // add a new object to the list
    static void Add(Oper *oper)
    {
      if ( !Get().insert(oper).second )
      {
        FAIL( "XLAllocatedList already contained this object" );
      }
    }

    // check if the object is in the list and remove it from it if it is
    static bool RemoveIfPresent(Oper *oper)
    {
      return Get().erase(oper) != 0;
    }

  private:
    typedef std::set<Oper *> List;

    static List& Get()
    {
      static List s_list;
      return s_list;
    }
  };


  // constructor from type only used by Missing() static ctor
  Oper(Type type, struct TypeTag *)
  {
    xltype = static_cast<WORD>(type);
  }

  // constructor for error values used by ErrValue() and others
  Oper(Err err, struct ErrTag *)
  {
    xltype = Type_Err;
    val.err = static_cast<WORD>(err);
  }

  // initializes this object with the string of the given length
  //
  // does not free the existing data, doesn't throw
  void InitWithStringNoThrow(const char *s, size_t len)
  {
    xltype = xltypeStr;

    // no need for trailing NUL for Pascal string so don't use strdup()
    val.str = static_cast<char *>(malloc(len + 1));
    *val.str = static_cast<char>(len);
    memcpy(val.str + 1, s, len);
  }

  // same as InitWithString() but throws if the string is too long
  //
  // throws if the string is too long
  void InitWithString(const char *s, size_t len)
  {
    if ( len > 255 )
    {
      throw EXCEPTION_MSG(xlerrValue, "String is too long");
    }

    InitWithStringNoThrow(s, len);
  }

  // initializes this object to be an array of the specified size
  //
  // doesn't allocate memory for the array elements, use SetArrayData() for
  // this
  //
  // does not free the existing data, doesn't throw
  void InitAsArrayNoThrow(WORD w, WORD h)
  {
    xltype = xltypeMulti;
    val.array.columns = w;
    val.array.rows = h;
  }

  // fills in the array elements by copying from src (which must be non NULL)
  //
  // the array must have been initialized with InitAsArrayNoThrow() before
  void SetArrayData(const Oper *src)
  {
    ASSERT_MSG( xltype == xltypeMulti, "must be an array" );

    const WORD size = val.array.columns*val.array.rows;
    Oper *dst = new Oper[size];

    val.array.lparray = dst;
    for ( size_t i = 0; i < size; ++i )
      *dst++ = *src++;
  }

  void InitAsRefNoThrow(const Ref& ref)
  {
    xltype = xltypeRef;
    val.mref.idSheet = ref.GetSheetId();
            
    if ( ref.HasReference() )
    {
      val.mref.lpmref = new XLMREF;
    
      XLMREF& xlmref = *val.mref.lpmref; 
    
      xlmref.count = 1; // we don't support multiple ranges
      xlmref.reftbl[0].rwFirst = ref.GetFirstRow();
      xlmref.reftbl[0].rwLast = ref.GetLastRow();
      xlmref.reftbl[0].colFirst = ref.GetFirstCol();
      xlmref.reftbl[0].colLast = ref.GetLastCol();
    }
    else // set it to NULL
    {
      val.mref.lpmref = NULL;
    }
  }

  void InitAsRefNoThrow(const XLMREF* const lpmref,const DWORD idSheet )
  {
    xltype = xltypeRef;
    val.mref.idSheet = idSheet;

    if ( lpmref )
    {
      val.mref.lpmref = new XLMREF;
      XLMREF& xlmref = *val.mref.lpmref; 

      xlmref.count = lpmref->count; 
      xlmref.reftbl[0].rwFirst = lpmref->reftbl[0].rwFirst;
      xlmref.reftbl[0].rwLast = lpmref->reftbl[0].rwLast;
      xlmref.reftbl[0].colFirst = lpmref->reftbl[0].colFirst;
      xlmref.reftbl[0].colLast = lpmref->reftbl[0].colLast;
    }
    else
    {
      val.mref.lpmref = NULL;
    }
  }

  // return true if this object has any dynamically allocated memory
  bool NeedsToBeFreed() const
  {
    return (xltype & (xltypeStr | xltypeMulti | xltypeRef)) != 0;
  }

  // mark this object as being allocated by Excel
  void MarkFromExcel()
  {
    if ( !NeedsToBeFreed() )
      return;

    // we can't just set xlbitXLFree bit directly here as this would result in
    // a failure if m_oper is passed to another Excel function (as it would
    // free it and then fail the call) before it is destroyed, so instead put
    // it on XL-allocated list
    XLAllocatedList::Add(this);
  }

  // frees the existing data if necessary
  //
  // doesn't reset the type to xltypeNil nor anything else, the caller must do
  // it
  void Clear()
  {
    if ( NeedsToBeFreed() )
      DoClear();
  }

  // frees the existing data unconditionally, shouldn't be called for types
  // which don't have to be freed (such as numbers, bools, ints, ...)
  void DoClear();

  // overwrite the contents of this object with the other one doing deep copy
  // if necessary
  void Copy(const Oper& other);

  // does the deep copy, should be only called if it is really needed
  //
  // copies the data only, the type is set by Copy() in all cases
  void DeepCopy(const Oper& other);


  // it calls our MarkFromExcel()
  friend class ByRef;
};

// ----------------------------------------------------------------------------
// ByRef inline functions implementation
// ----------------------------------------------------------------------------

inline ByRef::~ByRef()
{
  m_oper.MarkFromExcel();
}

inline ByRef::operator XLOPER *() const
{
  return const_cast<XLOPER *>(m_oper.ByVal());
}

} // namespace XL

} // namespace ito33

// need to include this after Oper class is defined but before the inline
// functions below which use Call()
#include "ito33/XL/call.h"

// ----------------------------------------------------------------------------
// Oper inline functions implementation
// ---------------------------------------------------------------------------

inline bool ito33::XL::Oper::AsArray(Oper *o) const
{
  return CallNoThrow(xlCoerce, ByRef(*o), *this, Type_Multi) == xlretSuccess;
}

inline void ito33::XL::Oper::DoClear()
{
  if ( (xltype & xlbitXLFree) || XLAllocatedList::RemoveIfPresent(this) )
  {
    // this was allocated by XL and must be freed by it
    CallNoThrow(xlFree, NULL, *this);
  }
  else // we don't care if xlbitDLLFree is set, we must free memory anyhow
  {
    switch ( GetType() )
    {
      case Type_Str:
        free(val.str);
        break;

      case Type_Multi:
        delete [] static_cast<Oper *>(val.array.lparray);
        break;

      case Type_Ref:
        delete val.mref.lpmref;
        break;

      default:
        FAIL( "shouldn't be called for types not needing freeing" );
    }
  }
}

inline void ito33::XL::Oper::DeepCopy(const Oper& other)
{
  switch ( other.GetType() )
  {
    case Type_Str:
      InitWithStringNoThrow(other.val.str + 1, *other.val.str);
      break;

    case Type_Multi:
      InitAsArrayNoThrow(other.val.array.columns, other.val.array.rows);
      SetArrayData(static_cast<Oper *>(other.val.array.lparray));
      break;

    case Type_Ref:
      InitAsRefNoThrow(other.val.mref.lpmref,other.val.mref.idSheet);
      break;
      
    default:
      FAIL( "shouldn't be called for types not needing deep copy" );
  }
}

inline void ito33::XL::Oper::Copy(const Oper& other)
{
  xltype = other.xltype;

  if ( other.NeedsToBeFreed() )
  {
    DeepCopy(other);
  }
  else // simple value
  {
    memcpy(&val, &other.val, sizeof(val));
  }
}

inline ito33::XL::Oper& ito33::XL::Oper::operator=(const Oper& other)
{
  Clear();
  Copy(other);

  return *this;
}

inline ito33::XL::Oper::~Oper()
{
  Clear();
}

#endif // _ITO33_XL_OPER_H_

