//----------------------------------------------------------------------------
//
//   Group       : QR Credit
//
//   Filename    : OptimizerSpellucci.hpp
//
//   Description : Wrapper class for Spellucci's C code for the DONLP2 minimization algorithm
//                 This wrapper supports only unconstrained optimization
//
//   Date        : 31 August 2006
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_OPTIMIZERSPELLUCCI_CPP
#include "edginc/QSpellucci.hpp"
#include <cmath>
#include <valarray>
#include <set>
#include <limits>

#define TRY
#define CATCH(arg1)
#define ASSERT0(arg1, arg2) 
#define ASSERT1(arg1, arg2) 
#define ASSERT2(arg1, arg2) 

// typdefs
typedef std::map<int, double> DomainBound;
typedef std::map<int, double>::iterator DomainBoundIterator;

DRLIB_BEGIN_NAMESPACE

///////////////////////////////////////////////// MISCUTIL.H ///////////////////////////////////////

/*!
  \file miscutil.h

  \brief Miscellaneous low level utilities.
*/

//! Utility for packed matrices
inline unsigned int Triangle(unsigned int i, unsigned int j = 0)
{ 
TRY
  return (i * (i + 1)) / 2 + j;
CATCH(Triangle) 
}

//! Utility for packed matrices
inline unsigned int Triangle(unsigned int k, unsigned int l, unsigned int size)
{ 
TRY
  return size * k - (k * (k - 1)) / 2 + l;
CATCH(Triangle) 
}


/*! 
  \brief Return iterator to prior element. 

  Note that we pass the input interator by value so as to avoid
  modifying the original.

  Template argument \c Iter must be a ``forward iterator'', ie, \c operator++
  must be defined.
*/
template<class Iter>
Iter Next(Iter i)
{
  return ++i;
}
  
/*! 
  \brief Return iterator to next element. 

  Note that we pass the input interator by value so as to avoid
  modifying the original.

  Template argument \c Iter must be a ``backward iterator'', ie, \c operator--
  must be defined.
*/
template<class Iter>
Iter Prior(Iter i)
{
  return --i;
}

/*!
  \brief 
   Copy elements from container of type C1 into
   container of type C2.

   C2 must support member push_back().

template<class C1, class C2>
void ContainerToContainer(const C1& c1, C2& c2)
{
  for (C1::const_iterator iter = c1.begin(); iter != c1.end(); ++iter)
    c2.push_back(*iter);
}
*/

/*!
  \brief Work-around for the STL algorithm for_each

  STL doesn't support in-place replacement of elements of an STL 'set' 
  container as in-place replacement of elements could violate the requirement
  that elements are always sorted in ascending order.  Therefore the STL for_each
  algorithm cannot be used for STL sets. ForEachInSet provides this  functionality by making a local copy of the set, erasing the 
  elements in the original set, 
  and then applying Func to each element before inserting them in the 
  original set.
  Func must be a unary function returning void. Use an adaptor if Func is a
  member function!
*/    
template<class T, class Func>
void ForEachInSet(std::set<T> &st, Func func)
{
  std::set<T> _st(st);
  st.erase(st.begin(), st.end());

  typedef typename std::set<T>::iterator Iterator;

  for (Iterator iter = _st.begin(); iter != _st.end(); ++iter)
  {
    T t(*iter);
    func(t);
    st.insert(t);
  }
}

//! Utility for shifting values
template<class Arg>
void Shift(Arg &a,
           Arg &b,
           Arg &c,
           const Arg &d) 
{
  a = b;
  b = c;
  c = d;
}

//! Utility for shifting values
template<class Arg>
void Shift(Arg &a,
           Arg &b,
           const Arg &c)
{
  a = b;
  b = c;
}

//! Utility for swapping two values
template<class Arg>
void Shift(Arg &a,
           Arg &b)
{
  Arg tmp;
  Shift(tmp, a, b, tmp);
}

//}


///////////////////////////////////////////////// NORM.H ////////////////////////////////////////////

// namespace {

/*!

  \file norm.h

  \brief Vector norms

*/



/*! 
  \brief Calculates the OneNorm of subarray

  \arg V 'v's type

  \par Note: Template type V must suport size, value_type and operator[size_t i].
  \param v A vector
  \param begin Element to begin norm calculation.
  \param end Element after end norm calculation.

  \returns OneNorm

  \throw Argument Throw exception if start and end is outside of 'v's range.
*/

template<class V>
typename V::value_type OneNorm(const V &v, size_t begin, size_t end)
{
TRY
  ASSERT1(begin >= 0 && begin <= end && end <= v.size(),
    Argument("Illegal subarray range."));
  // check for numerical zero
  size_t i = 0;
  typename V::value_type absmax = Abs(v[begin]);
  for (i = begin + 1; i < end; ++i) 
    absmax = Max(absmax, Abs(v[i]));
  if (absmax == 0)
    return 0;
  typename V::value_type nv = 0.0;
  for (i = begin; i < end; ++i)
    nv += Abs(v[i]) / absmax;
  nv *= absmax;
  return nv;
CATCH(basnum::OneNorm)
}

/*! 
  \brief Calculates the OneNorm of 'v'

  \arg V 'v's type

  \par Note: Template type V must suport size, value_type and operator[size_t i].
  \param v A vector

  \returns OneNorm of 'v'
*/

template<class V>
inline typename V::value_type OneNorm(const V &v)
{
  return OneNorm(v, 0, v.size());
}

/*! 
  \brief Calculates the TwoNorm of subarray

  \arg V 'v's type

  \par Note: Template type V must suport size, value_type and operator[size_t i].
  \param v A vector
  \param begin Element to begin norm calculation.
  \param end Element after end norm calculation.

  \returns TwoNorm

  \throw Argument Throw exception if start and end is outside of 'v's range.
*/
  
template<class V>
typename V::value_type TwoNorm(const V &v, size_t begin, size_t end)
{
TRY
  ASSERT1(begin >= 0 && begin <= end && end <= v.size(),
    Argument("Illegal subarray range."));
  // check for numerical zero
  size_t i = 0;
  typename V::value_type absmax = Abs(v[begin]);
  for (i = begin + 1; i < end; ++i) 
    absmax = Max(absmax, Abs(v[i]));
  if (absmax == 0)
    return 0;
  typename V::value_type nv = 0.0;
  for (i = begin; i < end; ++i)
    nv += Square(Abs(v[i]) / absmax);
  nv = absmax * sqrt(nv);
  return nv;
CATCH(basnum::TwoNorm)
}
/*! 
  \brief Calculates the TwoNorm of 'v'

  \arg V 'v's type

  \par Note: Template type V must suport size, value_type and operator[size_t i].
  \param v A vector

  \returns TwoNorm of 'v'
*/  
template<class V>
inline typename V::value_type TwoNorm(const V &v)
{
  return TwoNorm(v, 0, v.size());
}
/*! 
  \brief Calculates the max norm of subarray

  \arg V 'v's type

  \par Note: Template type V must suport size, value_type and operator[size_t i].
  \param v A vector
  \param begin Element to begin norm calculation.
  \param end Element after end norm calculation.

  \returns max norm

  \throw Argument Throw exception if start and end is outside of 'v's range.
*/
template<class V>
typename V::value_type InfNorm(const V &v, size_t begin, size_t end)
{
TRY
  ASSERT1(begin >= 0 && begin <= end && end <= v.size(),
    Argument("Illegal subarray range."));
  typename V::value_type nv = 0.0;
  for (size_t i = begin; i < end; ++i)
    nv = Max(nv, Abs(v[i]));
  return nv;
CATCH(basnum::InfNorm)
}
/*! 
  \brief Calculates the max norm of 'v'

  \arg V 'v's type

  \par Note: Template type V must suport size, value_type and operator[size_t i].
  \param v A vector

  \returns max norm of 'v'
*/ 
template<class V>
inline typename V::value_type InfNorm(const V &v)
{
  return InfNorm(v, 0, v.size());
}


/*!
  \brief Abstract Norm class
  
  A Norm defines a norm on the space of valarray<double>'s.
*/
template<class V = valarray<double> >
class Norm
{
public:
// foundation
  //! Destructor
  virtual ~Norm() 
  {
  }
// inspection 
  /*! 
    \brief Return (vector-)norm of valarray.
    \param v Vector
    Note: for an empty vector the return is 0.0.
  */
  typename V::value_type operator()(const V &v) const
  {
    return operator()(v, 0, v.size());
  }
  /*! 
    \brief Return (vector-)norm of subarray 
    \param v Vector.
    \param begin Index of first subarray element.
    \param end One plus index of last subarray element.
    Note: if begin >= end the return is 0.0.
  */
  virtual typename V::value_type operator()(const V &v,
                                            size_t begin,
                                            size_t end) const = 0;
}; // Norm   

/*! 
  \brief A PowerNorm defines a norm on the space of valarray<double>'s.
*/
template<class V = valarray<double> >
class PowerNorm : public Norm<V>
{
public:
// foundation
  explicit PowerNorm(int power)
  : m_power(power),
    m_root(1.0 / static_cast<double>(m_power))
  {
  }
// inspection
  typename V::value_type operator()(const V &v,
                                    size_t begin,
                                    size_t end) const
  {
  TRY
    if (m_power == 1)
      return OneNorm(v, begin, end);
    if (m_power == 2)
      return TwoNorm(v, begin, end);
    ASSERT1(begin <= end && end <= v.size(),
        Argument("Illegal subarray range."));
    // check for numerical zero
    size_t i = 0;
    typename V::value_type absmax = Abs(v[begin]);
    for (i = begin + 1; i < end; ++i) 
      absmax = Max(absmax, Abs(v[i]));
    if (absmax == 0)
      return 0;
    typename V::value_type nv = 0.0;
    for (i = begin; i < end; ++i)
      nv += pow(Abs(v[i]) / absmax, m_power);
    nv = absmax * pow(nv, m_root);
    return nv;
  CATCH(basnum::PowerNorm::operator())
  }

private:
// data
  //! The positive integer power defining the norm
  int m_power;
  //! The inverse power
  double m_root;
}; // PowerNorm    

/*! 
  \brief MaxNorm designates the infinity/max norm.
*/
template<class V = valarray<double> >
class MaxNorm : public Norm<V>
{
public:
// inspection
  typename V::value_type operator()(const V &v,
                                    size_t begin,
                                    size_t end) const
  {
  TRY
    return InfNorm(v, begin, end);
  CATCH(basnum::MaxNorm::operator())
  }
}; // MaxNorm    

// }



//////////////////////////////////////////////// SAFEPTR.H //////////////////////////////////////////

// namespace {

/*!
  \file safeptr.h
  \brief Exception safe "smart pointer".
*/


/*! 
  \brief A ReferenceCount provides reference counting.

  This class can supply reference counting to a derived class which is 
  responsible for managing a resource to which multiple references
  may exist at the same time.
  By inheriting from ReferenceCount the derived class gets access 
  to the count and this may be shared with other derived types.

  The class has no public interface as it is vital that it be accessed
  only by the managing class.

  See the SafePtr class template for an instance of use.
*/
class ReferenceCount
{
protected:
// Foundation
  /*! 
    The destructor decrements the count and deletes
    counter iff count reaches zero.
    The destructor is pure virtual (only) in order to make 
    the class abstract,
  */
  virtual ~ReferenceCount() = 0;
  /*!
    \brief Default constructor
    \throw BadAlloc iff (the in-line call of) \c new() returned null pointer
  */
  ReferenceCount() 
  : m_refcount(new int(1))
  {
  }
  //! Copy constructor
  ReferenceCount(const ReferenceCount &spb) 
  : m_refcount(spb.m_refcount)  {
    ++(*m_refcount);
  }

  //! Copy assignment
  ReferenceCount& operator=(const ReferenceCount &spb) 
  {
    // self-assignment guard
    if (this != &spb)
    {
      // Warning: Do not call virtual destructor.
      // this->~ReferenceCount(); calls ~SafePtr first, if 'this' is a SafePtr.
      // Use private member function instead.
      mf_DecrementReferenceCount();
      m_refcount = spb.m_refcount;      ++(*m_refcount);
    }
    return *this;
  }

// Inspection
  int mf_Count() const 
  {
    return *m_refcount;
  }

private:
  /*! 
    \brief Decrements the count and deletes
     counter iff count reaches zero.
  */
  void mf_DecrementReferenceCount()
  {
  try
  {
    // decrement reference count and delete if zero  
    if (--(*m_refcount) == 0)
      delete m_refcount;
  }
  catch(...)
  {
  }
  }
  //! Pointer to reference count
  int* m_refcount;
}; // ReferenceCount


inline ReferenceCount::~ReferenceCount() 
{
  mf_DecrementReferenceCount();
}


/*! 
  \brief Safe pointers provide exception-safe automatic memory management
         while supporting all useful pointer functionality.

  This class guarantees that dynamic memory allocated to pointers
  is freed under all circumstances (exception throws, normal execution) 
  without the program explicitly calling operator delete.
  It also guarantees that no such memory is deleted while there are
  still 'live' references to it.

  A safe pointer can only be modified by the copy assignment operator.

  The idea is to never use 'naked' pointers in code; pointers should always 
  come 'dressed' as safe pointers. A safe pointer is constructed in one of  two ways: 
   - at the point of allocation using SafePtr<T>(new T(...))
   - by copying from another safe pointer

  Safe pointers support polymorphism, ie, a safe pointer to a derived class 
  can be used wherever a safe pointer to its base class is required.

  Safe pointers may be (copy) constructed from safe pointers to other classes
  if the corresponding cast of naked pointers is legal.

  Naturally, safe pointers support the standard operations for access through
  a pointer, ie, a safe pointer may be dereferenced by either '*' or '->' and 
  will return a reference or a pointer, respectively, to the underlying object.
  An attempt to dereference a null pointer throws an exception.

  A safe pointer may be queried for null-ness through the member function 
  IsNull().

  Helper functionality:

  Comparison is defined for safe pointers wrt equality as well as ordering.
  (See separate entry)

  Note: we intended to implement casting of safe pointers, but this idea was 
  abandoned due to lack of (MSVC 6.0) compiler support for friend templates
  (INTERNAL COMPILER ERROR C1001: msvc1.cpp, line 1794). Casting can be done
  only on the corresponding naked pointer by extracting the naked pointer 
  (be very careful here!) and performing the native cast on this. It is not
  possible to obtain a safe pointer of a different type pointing to the same
  object.*/
template<class T> 
class SafePtr : public ReferenceCount
{
public:
  // Foundation
  /*! 
    \brief The destructor decrements the reference count and deletes underlying 
    memory if and only if handle count becomes zero.
  */
  ~SafePtr();

  //! Default constructor
  SafePtr() 
  : m_pt(0)
  {
  }

  /*! 
    \brief Construction from new'ed pointer
    \throw BadAlloc iff (the in-line call of) \c new() returned null pointer

  */
  SafePtr(T *pt) 
  : m_pt(pt)
  {
  }
  /*!
    \brief Copy constructor taking SafePtr of any type.
     Construction will succeed only if DT is publicly derived from T.
  */
  template<class DT> 
  SafePtr(const SafePtr<DT> &spdt) 
  : ReferenceCount(spdt),
    m_pt(spdt.IsNull() ? 0 : &(*spdt))  {
  }
  /*!
    \brief Copy constructor taking SafePtr of same type

     The same-type copy constructor is a specialization 
     of the general copy constructor template above.
  */
  SafePtr(const SafePtr &spt) 
  : ReferenceCount(spt),
    m_pt(spt.m_pt)  {
  }
  /*!
    \brief Copy asignment from SafePtr of any type.
    
     The copy will succeed only if DT is publicly derived from T.
  */
  template<class DT> 
  SafePtr &operator=(const SafePtr<DT> &rhs)
  {
    // This function body would not compile (with MSVC 6.0) when moved    // outside class definition!

    // construct intermediate and copy to this
    *this = SafePtr(rhs);

    // return
    return *this;
  }
  /*!
    \brief Copy assigment to SafePtr of same type
    
     The same-type copy assignment operator is a specialization 
     of the general assignment operator template above.
  */
  SafePtr &operator=(const SafePtr &rhs)
  {
    // This function body would not compile (with MSVC 6.0) when moved    // outside class definition!

    // Self-assignment guard is more involved for SafePtr's as we want to
    // guard two situations:
    // 1) Assign the same physical SafePtr to itself, ie same physical 
    //    address on '*this' and on 'spt'
    // 2) Assign one SafePtr to another both pointing to the same object 
    //    instance, ie different physical addresses on '*this' and on 'spt'
    //    but referenced objects *(this->m_pt) and *(spt.m_spt) have same    //    physical address. This situation will not result in an error,
    //    avoiding it is only a matter of efficiency.
    if (this != &rhs && m_pt != rhs.m_pt)
    {
      // Delete underlying iff reference count is one
      // Note: we cannot call the destructor as this will invoke the 
      // ReferenceCount destructor which will also be called
      // by the call to 'ReferenceCount::operator=' in the next line.
      if (mf_Count() == 1)
        delete m_pt;
      // copy base class
      this->ReferenceCount::operator=(rhs);
      // copy naked pointer
      m_pt = rhs.m_pt;
    }
    // return
    return *this;
  }

  // Inspection
  /*!
    \brief Return reference to underlying object.

    \throw Pointer("Attempt to dereference null pointer")
           iff this is a null pointer.
  */
  T &operator*() const;

  /*!
    \brief Return member of underlying object

    \throw Pointer("Attempt to access through null pointer")
           iff this is a null pointer.
  */
  T *operator->() const;

  //! Return true iff this is a null pointer.
  bool IsNull() const 
  {
    return !m_pt;
  }

private:
  //! pointer to underlying object
  T* m_pt;
}; // SafePtr

template<class T> 
SafePtr<T>::~SafePtr() 
{
try
{
  // Check if reference count is one, ie, the safe pointer under destruction is
  // the last remaining reference.
  if (mf_Count() == 1)
    // This scope can only be reached during the life of the *m_pt object.
    delete m_pt;
}
catch(...){
}
}


template<class T> 
inline 
T& SafePtr<T>::operator*() const
{
  return *m_pt;
}

// return member of underlying object
template<class T> 
inline 
T* SafePtr<T>::operator->() const 
{
  return m_pt;
}

// Helper functions

/*!
  \page ordering_doc Definition of ordering operators for safe pointers

  We define the operators '<, >, <=, >=' for safe pointers such that
  the return value is the result of the corresponding comparisons of 
  the underlying objects (thus requiring that such comparison be suitably
  defined).

  The purpose of these definitions is to ease the use of safe pointers
  together with STL containers and algorithms. For example, in many STL
  algorithms the 'less' predicate object is the default. With the  definition given here for 'operator <' for safe pointers 'less' will compare
  underlying values and this is almost always what is desired.

  BEWARE that while operators == and != compare the IDENTITIES of underlying
  objects the ordering operators compare the VALUES of these. This means for  instance, that != is NOT EQUIVALENT to < || >.*/

/*! 
  \brief Return true iff safe pointers have same underlying object

  \relates SafePtr

  \ref ordering_doc "Goto ordering documentation"
*/
template<class T, class S>  
inline 
bool operator==(const SafePtr<T> &lhs, 
                  const SafePtr<S> &rhs)
{
  // check if underlying objects are identical
  return &(*lhs) == &(*rhs);
}

/*! 
  \brief Return true iff safe pointers have distinct underlying objects

  \relates SafePtr

  \ref ordering_doc "Goto ordering documentation"
*/
template<class T, class S>  
inline 
bool operator!=(const SafePtr<T> &lhs, 
                  const SafePtr<S> &rhs)
{
  return !(lhs == rhs);
}

/*! 
  \brief less than operator

  \relates SafePtr

  \ref ordering_doc "Goto ordering documentation"
*/
template<class T>
inline bool operator<(const SafePtr<T> &lhs, const SafePtr<T> &rhs)
{
  return *lhs < *rhs;
}
/*! 
  \brief greater than operator

  \relates SafePtr

  \ref ordering_doc "Goto ordering documentation"
*/
template<class T>
inline bool operator>(const SafePtr<T> &lhs, const SafePtr<T> &rhs)
{
  return *lhs > *rhs;
}
/*! 
  \brief less equal than operator

  \relates SafePtr

  \ref ordering_doc "Goto ordering documentation"
*/
template<class T>
inline bool operator<=(const SafePtr<T> &lhs, const SafePtr<T> &rhs)
{
  return !(*lhs > *rhs);
}
/*! 
  \brief greater equal than operator

  \relates SafePtr

  \ref ordering_doc "Goto ordering documentation"
*/
template<class T>
inline bool operator>=(const SafePtr<T> &lhs, const SafePtr<T> &rhs)
{
  return !(*lhs < *rhs);
}

//} // anmonymous namespace


//////////////////////////////////////////////// ADAPTOR.H //////////////////////////////////////////

// namespace{

/*!
  \file adaptor.h
  \brief Defines member function adaptors.

  \ref adaptor_doc "Goto documentation"
*/

/*!

  \page adaptor_doc Adaptor Documentation.

  We extend the set of STL adaptors with adaptors for member functions to be
  called on given objects.

  A number of algorithms take a function as input. Examples of this
  include root searching and root bracketing.  The function for which we want to perform a root search, may, however, come 
  in many different encapsulations e.g. non-member function, member function,
  object function, pointer to a function, etc.  The idea behind adaptors is to provide an interface between the 
  (encapsulated) function and the algorithm; that is, the adaptor unwraps the 
  function encapsulation and feeds the function to the algorithm in a format 
  proper for the algorithm.

  The adaptor is instantiated by a template function having as 
  return type a template class. It is this class that does the
  adaption of the non-conformant function to a simple function object.
  C++ compilers support deduction of template arguments for functions
  but not for classes. Adaptors are only used in setting up function calls
  and are there only used in-line and should therefore only
  exist as temporary objects.
  
  The generic algorithms, eg GoalSeek, requires
  a simple function object having a single input argument.
  The adaptor class wraps a member-function on which one wants to
  execute GoalSeek into a simple function object. To simplify
  the instantiation of these temporary adaptor classes adaptor
  functions are provided each returning an instance of the relevant 
  adaptor class.

  As an exemplification of the adaptor mechanism consider the following 
  example:
  \code
    // Define class containing a function
    class HasNonConstMember
    {
    public:
      HasNonConstMember(double c) : m_c(c) {}
      double func(double x) {return exp(x) - m_c;}

    private:
      double m_c;
    }; // HasNonConstMember

    // Initiate an object of the class
    HasNonConstMember sc(2.0);

    // Find root for the function using the adaptor function MemFunRef
    // as interface.
    GoalSeek(MemFunRef(sc,&HasNonConstMember::func), 
                       root, 0.0, 0.0, 2.0);

  \endcode
*/

/*! 
  \brief Internal adaptor for non-const member function of object 
  given by reference.

  The adaptor consists of a class template and a function template.
  The member function is assumed to take a single argument of type D
  and return void.
  \ref adaptor_doc "Goto adaptor documentation"
*/
template<class C, class D>
class CVoidMemFunRef
{
public:
  CVoidMemFunRef(C &c, void (C::* pmf) (D)) 
    : m_o(c), 
      m_pmf(pmf) 
  {
  }
  void operator() (D d) const 
  {
    (m_o.*m_pmf)(d);
  }

private:
  C &m_o;
  void (C::* m_pmf) (D);
}; // CVoidMemFunRef

/*! 
  \brief Adaptor for non-const member function of object given by reference.

  The adaptor consists of a class template and a function template
  The member function is assumed to take a single argument of type D
  and return void.
  \ref adaptor_doc "Goto adaptor documentation"
*/
template<class C, class D>
CVoidMemFunRef<C, D> 
MemFunRefVoid(C &c, void (C::* mf) (D))
{
  return CVoidMemFunRef<C, D>(c, mf);
}

/*! 
  \brief Internal adaptor for const member function of object given 
  by reference.

  The adaptor consists of a class template and a function template
  The member function is assumed to take a single argument of type D
  and return void.
  \ref adaptor_doc "Goto adaptor documentation"
*/
template<class C, class D>
class CConstVoidMemFunRef
{
public:
  CConstVoidMemFunRef(const C &c, void (C::* pcmf) (D) const) 
    : m_o(c), 
      m_pcmf(pcmf) 
  {
  }
  void operator() (D d) const 
  {
    (m_o.*m_pcmf)(d);
  }

private:
  const C &m_o;
  void (C::* m_pcmf) (D) const;
}; // CConstVoidMemFunRef

/*! 
  \brief Adaptor for const member function of object given by reference.

  The adaptor consists of a class template and a function template.
  The member function is assumed to take a single argument of type D
  and return void.
  \ref adaptor_doc "Goto adaptor documentation"
*/
template<class C, class D>
CConstVoidMemFunRef<C, D>
MemFunRefVoid(const C &c, void (C::* cmf) (D) const)
{
  return CConstVoidMemFunRef<C, D>(c, cmf);
}


/*! 
  \brief Internal adaptor for non-const member function 
  of object held in safe pointer

  The adaptor consists of a class template and a function template.
  The member function is assumed to take a single argument of type D
  and return type R.
  \ref adaptor_doc "Goto adaptor documentation"
*/
template<class C, class D, class R>
class CMemFun
{
public:
  // Constructor
  CMemFun(const SafePtr<C> &pc, R (C::* pmf) (D)) 
    : m_po(pc), 
      m_pmf(pmf) 
  {
  }

  R operator() (D d) const 
  {
    return (*m_po.*m_pmf)(d);
  }

private:
  const SafePtr<C> &m_po;
  R (C::* m_pmf) (D);
}; // CMemFun

/*! 
  \brief Adaptor for non-const member function of object held in safe pointer

  The adaptor consists of a class template and a function template.
  The member function is assumed to take a single argument of type D
  and return type R.
  \ref adaptor_doc "Goto adaptor documentation"
*/
template<class C, class D, class R>
CMemFun<C, D, R> 
MemFun(SafePtr<C> &pc, R (C::* mf) (D))
{
  return CMemFun<C, D, R>(pc, mf);
}

/*! 
  \brief Internal adaptor for const member function of 
  object held in safe pointer

  The adaptor consists of a class template and a function template.
  The member function is assumed to take a single argument of type D
  and return type R.
  \ref adaptor_doc "Goto adaptor documentation"
*/
template<class C, class D, class R>
class CConstMemFun
{
public:
  CConstMemFun(const SafePtr<C> &pc, R (C::* pmf) (D) const) 
    : m_po(pc), m_pmf(pmf) {}
  R operator() (D d) const 
  {
    return (*m_po.*m_pmf)(d);
  }

private:
  const SafePtr<C> &m_po;
  R (C::* m_pmf) (D) const;
}; // CConstMemFun

/*! 
  \brief Adaptor for const member function of object held in safe pointer

  The adaptor consists of a class template and a function template.
  The member function is assumed to take a single argument of type D
  and return type R.
  \ref adaptor_doc "Goto adaptor documentation"

*/
template<class C, class D, class R>
CConstMemFun<C, D, R> 
MemFun(const SafePtr<C> &pc, R (C::* mf) (D) const)
{
  return CConstMemFun<C, D, R>(pc, mf);
}

/*! 
  \brief Internal adaptor for non-const member 
  function of object given by reference.

  The adaptor consists of a class template and a function template.
  The member function is assumed to take a single argument of type D
  and return type R.
  \ref adaptor_doc "Goto adaptor documentation"

*/
template<class C, class D, class R>
class CMemFunRef
{
public:
  CMemFunRef(C &c, R (C::* pmf) (D)) 
    : m_o(c), 
      m_pmf(pmf) 
  {
  }
  R operator() (D d) const 
  {
    return (m_o.*m_pmf)(d);
  }

private:
  C &m_o;
  R (C::* m_pmf) (D);
}; // CMemFunRef

/*! 
  \brief Adaptor for non-const member function of object given by reference.

  The adaptor consists of a class template and a function template.
  The member function is assumed to take a single argument of type D
  and return type R.
  \ref adaptor_doc "Goto adaptor documentation"

*/
template<class C, class D, class R>
CMemFunRef<C, D, R> 
MemFunRef(C &c, R (C::* mf) (D))
{
  return CMemFunRef<C, D, R>(c, mf);
}

/*! 
  \brief Internal adaptor for const member function 
  of object given by reference.

  The adaptor consists of a class template and a function template.
  The member function is assumed to take a single argument of type D
  and return type R.
  \ref adaptor_doc "Goto adaptor documentation"
*/
template<class C, class D, class R>
class CConstMemFunRef
{
public:
  CConstMemFunRef(const C &c, R (C::* pcmf) (D) const) 
    : m_o(c), 
      m_pcmf(pcmf) 
  {
  }
  R operator() (D d) const 
  {
    return (m_o.*m_pcmf)(d);
  }

private:
  const C &m_o;
  R (C::* m_pcmf) (D) const;
}; // CConstMemFunRef

/*! 
  \brief Adaptor for const member function of object given by reference.

  The adaptor consists of a class template and a function template.
  The member function is assumed to take a single argument of type D
  and return type R.
  \ref adaptor_doc "Goto adaptor documentation"
*/
template<class C, class D, class R>
CConstMemFunRef<C, D, R>
MemFunRef(const C &c, R (C::* cmf) (D) const)
{
  return CConstMemFunRef<C, D, R>(c, cmf);
}

/*! 
  \brief Internal adaptor for const member function 
  of object given by reference.

  The adaptor consists of a class template and a function template.
  The member function is assumed to be a binary function with arguments 
  of type D1 and type D2 and return type R.  
  This adaptor is used to bind the first argument making the adapted 
  function a unary function.
  \ref adaptor_doc "Goto adaptor documentation"
*/
template<class C, class D1, class D2, class R>
class CConstMemFunRefBind1st
{
public:
  CConstMemFunRefBind1st(const C &c, D1 first, R (C::* pcmf) (D1, D2) const) 
    : m_first(first),
      m_o(c), 
      m_pcmf(pcmf) 
  {
  }
  R operator() (D2 d) const 
  {
    return (m_o.*m_pcmf)(m_first, d);
  }

private:
  D1 m_first;
  const C &m_o;
  R (C::* m_pcmf) (D1, D2) const;
}; // CConstMemFunRefBind1st

/*! 
  \brief Adaptor for const member function of object given by reference.

  The adaptor consists of a class template and a function template.
  The member function is assumed to be a binary function with arguments 
  of type D1 and type D2 and return type R.  
  This adaptor is used to bind the first argument making the adapted 
  function a unary function.
  \ref adaptor_doc "Goto adaptor documentation"
*/
template<class C, class D1, class D2, class R>
CConstMemFunRefBind1st<C, D1, D2, R>
MemFunRefBind1st(const C &c, D1 first, R (C::* cmf) (D1, D2) const)
{
  return CConstMemFunRefBind1st<C, D1, D2, R>(c, first, cmf);
}

/*! 
  \brief Internal adaptor for const member function 
  of object given by reference.

  The adaptor consists of a class template and a function template.
  The member function is assumed to be a binary function with arguments 
  of type D1 and type D2 and return type R.  
  This adaptor is used to bind the second argument making the adapted 
  function a unary function.
  \ref adaptor_doc "Goto adaptor documentation"
*/
template<class C, class D1, class D2, class R>
class CConstMemFunRefBind2nd
{
public:
  CConstMemFunRefBind2nd(const C &c, D2 second, R (C::* pcmf) (D1, D2) const) 
    : m_second(second),
      m_o(c), 
      m_pcmf(pcmf) 
  {
  }
  R operator() (D1 d) const 
  {
    return (m_o.*m_pcmf)(d, m_second);
  }

private:
  D2 m_second;
  const C &m_o;
  R (C::* m_pcmf) (D1, D2) const;
}; // CConstMemFunRefBind2nd

/*! 
  \brief Adaptor for const member function of object given by reference.

  The adaptor consists of a class template and a function template.
  The member function is assumed to be a binary function with arguments 
  of type D1 and type D2 and return type R.  
  This adaptor is used to bind the second argument making the adapted 
  function a unary function.
  \ref adaptor_doc "Goto adaptor documentation"
*/
template<class C, class D1, class D2, class R>
CConstMemFunRefBind2nd<C, D1, D2, R>
MemFunRefBind2nd(const C &c, D2 second, R (C::* cmf) (D1, D2) const)
{
  return CConstMemFunRefBind2nd<C, D1, D2, R>(c, second, cmf);
}

// } // anonymous namespace


//////////////////////////////////////////////// NUMUTIL.H //////////////////////////////////////////

// namespace{

/*!
  \file numutil.h

  \brief Definition of some low level numerical utilities.
*/

// Some people out there (eg Microsoft) like to define these as macros - let's get rid of that
#undef Max
#undef Min

//#pragma warning(disable : 4786)


/*!
  \brief Base of natural logarithm.
*/  
const double e = 2.71828182845905;


/*!
  \brief Ratio of circumference to diameter for flat circle.
*/
const double pi = 3.1415926535897932384626;

// one basis point
const double onebasispoint = 0.0001;

// days in a year
const double daysinoneyear = 365.0;

/*!
  \brief Normalization of standardized normal distribution (mean zero, unit variance)

  This is give as one divided by the square root of twice pi.
*/
const double oneoversqrttwopi = 0.398942280401433;

/*!
  \brief Maximum of two objects.

  Template function that the larger of two objects.
  
  Requires 'less than' (operator<) to be defined for the template class.

  Reimplements std::max since it is called _MAX in the Microsoft implemention.

  \return y if y > x and x otherwise
*/

template<class T> 
inline const T& Max(const T &x,
                    const T &y)
{
  return (x < y) ? y : x;
}

/*!
  \brief Minimum of two objects.
                 
  Template function to find the smaller of two objects.
  
  Requires 'less than' (operator<) to be defined for the template class.
  
  Reimplements std::min since it is called _MIN in the Microsoft implemention.

  \return x if x < y and y otherwise
*/

template<class T> 
inline const T& Min(const T &x,
                    const T &y)
{
  return (x < y) ? x  : y;
}



/*!

  \brief Absolute value of object.

  Template function to find the absolute value of an object.
  
  Requires 'less than' (operator<) and prefix unary minus (operator-) to  be
  defined for the template class.
  
  \return x if x > -x and -x otherwise

*/

template<class T> 
inline T Abs(const T &x)
{
  return Max(x, - x);
}

static inline
bool isZeroDouble(double x, double tolerance = 1e-6)
{
  return Abs(x) < tolerance;
}

//! Utility returning Absolute value of 'a' with sign of 'b'.
template<class Arg>
Arg Sign(Arg a, 
         Arg b) 
{
  return ((b) > 0.0) ? Abs(a) : -Abs(a);
}

/*!

  \brief Square of object.

  Template function to find the square value of an object.
  
  Requires multiplication (operator*) to be defined for the template class.
  
  \return x * x

*/

template<class T> 
inline T Square(T x)
{
  return x * x;
}

/*!
  \brief Return value of largest element in random access container.

  Complexity is linear.

  \par Template parameters
  \param C Type of random access container; must define \c value_type and support \c [] notation.

  \par Function arguments
  \param c Random access container.


template<class C>
typename C::value_type MaxElement(const C& c)
{
TRY
  C::value_type max;
  // Loop through container
  max = c[0];
  for (unsigned int i = 1; i < c.size(); ++i)
    if (c[i] > max)
      max = c[i];
  return max;
CATCH(Max)
}
*/

// } // anonymous namespace


//////////////////////////////////////////////// SLICEREF.H /////////////////////////////////////////

// namespace{ 

/*! 
  \file sliceref.h
  \brief Definition of matrix slice reference classes
*/

using std::valarray;
using std::slice;

// ********************* SliceRef classes *********************************

// forward declaration
template<typename T> class ConstSliceRef;

/*! 
  \brief A SliceRef provides random read/write access to a one-dimensional
  slice of a matrix.

  The SliceRef type is not intended to appear explicitly in client code---the
  intention is to leave all construction of SliceRef's to matrices and all 
  member function calls to helper and linear algebra functions.

  Note that despite appearances no operations actually act on a SliceRef---the
  action is always passed on to the matrix data referred to.

  Also note that a SliceRef is a reference type and that---as always with 
  references---it is up to the calling code to ensure that all references are valid.

  \par Template arguments
  \arg T The type of elements
*/
template<class T>
class SliceRef
{
public:
// friends
  //! Class ConstSliceRef is a friend in order to allow implicit conversion.
  friend class ConstSliceRef<T>;
// types 
  //! Typedef for the type of elements
  typedef T value_type;

// foundation
  /*! 
    \brief Construct from valarray of elements and a slice
    \param elem Vallaray of elements
    \param s The mapping (slice)
  */
  SliceRef(valarray<T> &elem, const slice &s);
  /*! 
    \brief Copy from any random access container 
    \par Template types
    \arg V Type of random access container to copy from; V must
    implement operator[](unsigned int) and size().
    \param v Vector
    \return A reference to SliceRef 
    \throw Argument Throw exception if sizes mismatch.
  */
  template<class V> SliceRef& operator=(const V &v)
  {
  TRY
    ASSERT1(size() == v.size(), Argument("Size mismatch."));
    for (unsigned int i=0; i < v.size(); ++i)
      operator[](i) = v[i];
    return *this;
  CATCH(SliceRef::operator=)
  }
  /*! 
    \brief Specialization of copy assignment to SliceRef 
    \param s Sliceref to be copied
    \return A reference to SliceRef 
    \throw Argument Throw exception if sizes mismatch.
  */
  //! Specialization of copy assignment to SliceRef
  SliceRef& operator=(const SliceRef &s)
  { 
  TRY
    if (this == &s)
      return *this;
    ASSERT1(size() == s.size(), Argument("SliceRef::operator=: Size mismatch."));
    for (unsigned int i=0; i < s.size(); ++i)
      operator[](i) = s[i];
    return *this;
    CATCH(SliceRef::operator=) 
  }

// inspection
  //! Return size of slice
  unsigned int size() const;
  //! Return value of i'th element (C notation)
  T operator[](unsigned int i) const;
  //! Return value of i'th element (Fortran notation)
  T operator()(unsigned int i) const;

// modification
  //! Return reference to i'th element (C notation)
  T& operator[](unsigned int i);
  //! Return reference to i'th element (Fortran notation)
  T& operator()(unsigned int i);
  //! Multiply slice elements by scalar
  SliceRef& operator*=(T a);

private:
// data
  valarray<T> &m_elem;
  slice m_s;
// functions
  //! Return reference to i'th element of slice
  T &mf_element(unsigned int i) const;
}; // SliceRef

/*!
  \brief A ConstSliceRef provides const (ie, read-only) access
  to a one-dimensional slice of a matrix

  The ConstSliceRef type is not intended to appear explicitly in client code---the
  intention is to leave all construction of ConstSliceRef's to matrices and all 
  member function calls to helper and linear algebra functions.

  Note that despite appearances no operations actually act on a ConstSliceRef---the
  action is always passed on to the matrix data referred to.

  Also note that a ConstSliceRef is a reference type and that---as always with 
  references---it is up to the calling code to ensure that all references are valid.

  \par Template arguments
  \arg T The type of elements
*/
template<class T>  
class ConstSliceRef
{
public:
// types 
  //! Typedef for the type of elements
  typedef T value_type;

// foundation
  /*! 
    \brief Construct from valarray of elements and a slice
    \param elem valarray of elements
    \param s The mapping (slice)
  */
  ConstSliceRef(const valarray<T> &elem, const slice &s);
  /*!
    \brief Construct from SliceRef.

    This constructor allows implicit conversion from SliceRef to ConstSliceRef.
    \param sr SliceRef to construct from.
  */
  ConstSliceRef(const SliceRef<T>& sr)
  : m_elem(sr.m_elem), m_s(sr.m_s)
  {
  }          

// inspection
  //! Return size of slice
  unsigned int size() const;
  //! Return value of i'th element (C notation)
  T operator[](unsigned int i) const;
  //! Return value of i'th element (Fortran notation)
  T operator()(unsigned int i) const;
  
private:
// data
  const valarray<T> &m_elem;
  slice m_s;
// functions
  //! Return value of i'th element of slice
  T mf_element(unsigned int i) const;
}; // ConstSliceRef

// *********************** SliceRef<> members *************************************
template<class T> inline
SliceRef<T>::SliceRef(valarray<T> &elem, const slice &s)
: m_elem(elem), 
  m_s(s)
{
}

template<class T> inline
unsigned int SliceRef<T>::size() const 
{
  return m_s.size();
}

template<class T> inline
T SliceRef<T>::operator[](unsigned int i) const 
{ 
TRY
  return mf_element(i);
CATCH(SliceRef::operator[]) 
}

template<class T> inline
T SliceRef<T>::operator()(unsigned int i) const 
{ 
TRY
  return mf_element(i);
CATCH(SliceRef::operator()) 
}

template<class T> inline 
T& SliceRef<T>::operator[](unsigned int i) 
{ 
TRY
  return mf_element(i);
CATCH(SliceRef::operator[]) }

template<class T> inline
T& SliceRef<T>::operator()(unsigned int i) 
{ 
TRY
  return mf_element(i);
CATCH(SliceRef::operator()) 
}

template<class T>
SliceRef<T>& SliceRef<T>::operator*=(T a)
{ 
TRY
  for (unsigned int i=0; i < size(); ++i)
    mf_element(i) *= a;
  return *this;
CATCH(SliceRef::operator*) 
}

template<class T> inline
T& SliceRef<T>::mf_element(unsigned int i) const
{ 
TRY
  ASSERT2(i < size(), Range("Index out of range."));
  return m_elem[m_s.start() + i * m_s.stride()];
CATCH(SliceRef::mf_element) 
}

// *********************** ConstSliceRef<> members *************************************
template<class T> inline  
ConstSliceRef<T>::ConstSliceRef(const valarray<T> &elem, const slice &s)
: m_elem(elem), 
  m_s(s)
{
}

template<class T> inline
unsigned int ConstSliceRef<T>::size() const 
{
  return m_s.size();
}

template<class T> inline
T ConstSliceRef<T>::operator[](unsigned int i) const 
{ 
TRY
  return mf_element(i);
CATCH(ConstSliceRef::operator[]) 
}

template<class T> inline 
T ConstSliceRef<T>::operator()(unsigned int i) const 
{ 
TRY
  return mf_element(i);
CATCH(ConstSliceRef::operator()) 
}
  
template<class T> inline 
T ConstSliceRef<T>::mf_element(unsigned int i) const
{ 
TRY
  ASSERT2(i < size(), Range("Index out of range."));
  return m_elem[m_s.start() + i * m_s.stride()];
CATCH(ConstSliceRef::mf_element) 
}

//} // anonymous namespace

//////////////////////////////////////////////// MATRIX.H ///////////////////////////////////////////

//namespace {

/*!
  \file matrix.h
  \brief Matrix type templates

  \par
  Design notes: 
  - layering of concrete types, not inheritance
  - matrix types defined by shapes; storage format 
    derived from shape.

  \par
  Conventions:
  - diagonals are numbered by integers witht the longest diagonal
    given the number 0 and diagonals below (above) given negative
    (positive) numbers
*/


// using declarations
using std::valarray;
using std::slice;

// ********************** Dense matrices *****************************

/*!
  \brief A Matrix is a general rectangular matrix.

  For this type of matrix there are no constraints on the elements
  and we store all (ncol x nrow) elements.
  \par Template types
  \arg T The element type
*/
template<class T>
class Matrix
{
public:
// types
  //! Type definition of element type
  typedef T value_type;

// foundation
  Matrix();
  /*! 
    \brief Construct from dimensions.
    \param nrow Number of rows - must be positive.
    \param ncol Number of columns - must be positive.
  */
  Matrix(unsigned int nrow, unsigned int ncol);
  /*!
    \brief Construct from value and dimensions.
    \param val Value assigned to all elements.
    \param nrow Number of rows - must be positive.
    \param ncol Number of columns - must be positive.
  */
  Matrix(T val, unsigned int nrow, unsigned int ncol);
  /*! 
    \brief Copy construct from given matrix of any type.
    \par Template arguments
    \arg M Matrix type; must implement members NRow(), NColumn() and
      operator()(unsigned int, unsigned int).
    \par Function arguments
    \param m Matrix to copy from
  */
  template<class M> explicit Matrix(const M &m)
  : m_nrow(m.NRow()),
    m_ncol(m.NColumn()),
    m_elem(m_nrow * m_ncol)
  { 
  TRY
    for (unsigned int i=0; i < m_nrow; ++i)
    {
      for (unsigned int j=0; j < m_nrow; ++j)
        m_elem[i * m_ncol + j] = m(i, j);
    }
  CATCH(Matrix) 
  }
  /*! 
    \brief Specialization of copy constructor to type Matrix
    \param m Matrix to copy from
  */
  Matrix(const Matrix &m);
  /*!
    \brief Copy assign from given matrix of any type.
    \par Template arguments
    \arg M Matrix type; must implement members NRow(), NColumn() and
      operator()(unsigned int, unsigned int).
    \par Function arguments
    \param rhs Matrix to copy from
  */
  template<class M> Matrix& operator=(const M &rhs)
  { 
  TRY
    m_nrow = rhs.NRow();
    m_ncol = rhs.NColumn();
    m_elem.resize(m_nrow * m_ncol);
    for (unsigned int i=0; i < m_nrow; ++i)
    {
      for (unsigned int j=0; j < m_nrow; ++j)
        m_elem[i * m_ncol + j] = rhs(i, j);
    }
    return *this;
  CATCH(Matrix::operator=) 
  }
  /*! 
    \brief Specialization of copy assignment to type Matrix
    \param rhs Matrix to copy from
  */
  Matrix& operator=(const Matrix &rhs);

// inspection
  //! Return # rows.
  unsigned int NRow() const;
  //! Return # columns.
  unsigned int NColumn() const;
  //! Return value of i'th element of j'th column.
  T operator()(unsigned int i, unsigned int j) const;
  //! Return i'th row for read-only access.
  ConstSliceRef<T> Row(unsigned int i) const;
  //! Return i'th column for read-only access.
  ConstSliceRef<T> Column(unsigned int i) const;

// modification
  //! Return reference to i'th element of j'th column.
  T& operator()(unsigned int i, unsigned int j);
  //! Return i'th row for read/write access
  SliceRef<T> Row(unsigned int i);
  //! Return i'th column for read/write access
  SliceRef<T> Column(unsigned int i);
  //! Multiply matrix by scalar
  Matrix& operator*=(T a);
  //! Resizing with initialization (defaulted)
  Matrix& Resize(unsigned int nrow, unsigned int ncol, T val=T())
  {
  TRY
    m_nrow = nrow;
    m_ncol = ncol;
    m_elem.resize(m_nrow * m_ncol, val);
    return *this;
  CATCH(Matrix::Resize) 
  }

  //! Clears Matrix contents
  void clear() { m_elem = valarray<T>(); }
private:
// data
  //! # rows, columns.
  unsigned int m_nrow, m_ncol;
  //! matrix elements.
  valarray<T> m_elem;
}; // class Matrix

/*!
  \brief A SquareMatrix is a generic square matrix

  For this type of matrix there are no constraints on the elements
  and we store all (ncol x nrow) elements.
  Note that diagonals are numbered such that the diagonal containing 
  element (i, 0) is numbered -i and the diagonal containing the element 
  (0, j) is numbered j.

  \par Template types
  \arg T The element type
*/
template<class T>
class SquareMatrix
{
public:
// types
  //! Type definition of element type
  typedef T value_type;

// foundation
  SquareMatrix();
  /*! 
    \brief Construct from size.
    \param size Number of rows/columns; must be positive.
  */
  explicit SquareMatrix(unsigned int size);
  /*! 
    \brief Construct from value and size.
    \param val Value assigned to all matrix elements.
    \param size Number of rows/columns; must be positive.
  */
  SquareMatrix(T val, unsigned int size);
  /*! 
    \brief Copy construct from given square matrix of any type.
    \par Template arguments
    \arg M Matrix type; must implement members NRow(), NColumn() and
      operator()(unsigned int, unsigned int).
    \par Function arguments
    \param m Matrix to copy from; must be square
    \throw Constructor Throw exception if m is not square.
  */
  template<class M> explicit SquareMatrix(const M &m)
  : m_size(m.NRow()),
    m_elem(m.NRow() * m.NColumn())
  { 
  TRY
    // make sure 'm' is square
    ASSERT1(m.NRow() == m.NColumn(), Constructor("Non-square input matrix."));
    for (unsigned int i=0; i < m_size; ++i)
    {
      for (unsigned int j=0; j < m_size; ++j)
        m_elem[i * m_size + j] = m(i, j);
    }
  CATCH(SquareMatrix) 
  }
  /*! 
    \brief Specialization of copy constructor to type SquareMatrix
    \param m SquareMatrix to copy from
  */
  SquareMatrix(const SquareMatrix<T> &m);
  /*!
    \brief Copy assign from given matrix of any type.
    \par Template arguments
    \arg M Matrix type; must implement members NRow(), NColumn() and
      operator()(unsigned int, unsigned int).
    \par Function arguments
    \param rhs Matrix to copy from; must be square
    \throw Argument Throw exception if rhs is not square.
  */
  template<class M> SquareMatrix& operator=(const M &rhs)
  { 
  TRY
    // make sure 'rhs' is square
    ASSERT1(rhs.NRow() == rhs.NColumn(), Argument("Non-square input matrix."));
    // copy data
    m_size = rhs.NRow();
    m_elem.resize(m_size * m_size);
    for (unsigned int i=0; i < m_size; ++i)
    {
      for (unsigned int j=0; j < m_size; ++j)
        m_elem[i * m_size + j] = rhs(i, j);
    }
    return *this;
  CATCH(SquareMatrix::operator=) 
  }
  /*! 
    \brief Specialization of copy asignment to type SquareMatrix.
    \param rhs SquareMatrix to copy from
  */
  SquareMatrix& operator=(const SquareMatrix &rhs);

// inspection
  //! Return size
  unsigned int Size() const;
  //! Return # rows.
  unsigned int NRow() const;
  //! Return # columns.
  unsigned int NColumn() const;
  //! Return value of i'th element of j'th column.
  T operator()(unsigned int i, unsigned int j) const;
  //! Return value of l'th element of k'th diagonal.
  T D(int k, unsigned int l) const;
  //! Return i'th row for read-only access.
  ConstSliceRef<T> Row(unsigned int i) const;
  //! Return i'th column for read-only access.
  ConstSliceRef<T> Column(unsigned int i) const;
  //! Return k'th diagonal for read-only access.
  ConstSliceRef<T> Diagonal(int k=0) const;

// modification
  //! Return reference to i'th element of j'th column.
  T& operator()(unsigned int i, unsigned int j);
  //! Return reference to l'th element of k'th diagonal.
  T& D(int k, unsigned int l);
  //! Return i'th row for read/write access.
  SliceRef<T> Row(unsigned int i);
  //! Return i'th column for read/write access.
  SliceRef<T> Column(unsigned int i);
  //! Return k'th diagonal for read/write access.
  SliceRef<T> Diagonal(int k=0);
  //! Multiply matrix by scalar
  SquareMatrix& operator*=(T a);
  //! Resizing
  SquareMatrix& Resize(unsigned int size, T val=T())
  {
  TRY
    m_size = size;
    m_elem.resize(Square(size), val);
    return *this;
  CATCH(SquareMatrix::Resize)
  }

private:
// data
  //! # size.
  unsigned int m_size;
  //! matrix elements.
  valarray<T> m_elem;
};  // SquareMatrix


/*!
  \brief A SymmetricMatrix is a symmetric matrix

  For this type of matrix we store only the lower triangle.
  All operations maintain symmetry.

  Note that diagonals are numbered such that the diagonal containing 
  element (i, 0) is numbered -i and the diagonal containing the element 
  (0, j) is numbered j.

  \par Template types
  \arg T The element type
*/
template<class T>
class SymmetricMatrix
{
public:
// typedef's
  //! Type definition of element type
  typedef T value_type;

// foundation
  SymmetricMatrix();
  /*! 
    \brief Construct from size.
    \param size Number of rows/columns; must be positive.
  */
  explicit SymmetricMatrix(unsigned int size);
  /*! 
    \brief Construct from value and size.
    \param val Value assigned to all matrix elements.
    \param size Number of rows/columns; must be positive.
  */
  SymmetricMatrix(T val, unsigned int size);
  /*! 
    \brief Copy constructor
    \param m SymmetricMatrix to copy from
  */
  SymmetricMatrix(const SymmetricMatrix<T> &m);
  /*!
    \brief Copy assignment
    \arg rhs Symmetric matrix assigned from
  */
  SymmetricMatrix& operator=(const SymmetricMatrix &rhs);

// inspection
  //! Return size.
  unsigned int Size() const;
  //! Return # rows
  unsigned int NRow() const;
  //! Return # columns
  unsigned int NColumn() const;
  //! Return value of j'th element of i'th row
  T operator()(unsigned int i, unsigned int j) const;
  //! Return value of l'th element of k'th diagonal
  T D(int k, unsigned int l) const;
  //! Return k'th diagonal for read-only access.
  ConstSliceRef<T> Diagonal(int k=0) const;

// modification
  //! Return reference to j'th element of i'th row
  T& operator()(unsigned int i, unsigned int j);
  //! Return reference to l'th element of k'th diagonal
  T& D(int k, unsigned int l);
  //! Return k'th diagonal for read/write access.
  SliceRef<T> Diagonal(int k=0);
  //! Multiply matrix by scalar
  SymmetricMatrix& operator*=(T a);
  //! Resizing
  SymmetricMatrix& Resize(unsigned int size, T val=T())
  {
  TRY
    m_size = size;
    m_elem.resize(Triangle(size, size, size) + 1, val);
    return *this;
  CATCH(SymmetricMatrix::Resize)
  }

private:
// data
  unsigned int m_size;
  valarray<T> m_elem;
}; // SymmetricMatrix

/*!
  \brief A LowerTriangleMatrix is a lower triangular matrix

  For this type of matrix we store only the lower triangle.
  All operations maintain shape.

  Note that diagonals are numbered such that the diagonal containing 
  element (i, 0) is numbered -i and the diagonal containing the element 
  (0, j) is numbered j.

  \par Template types
  \arg T The element type
*/
template<class T>
class LowerTriangleMatrix
{
public:
// typedef's
  //! Type definition of element type
  typedef T value_type;

// foundation
  LowerTriangleMatrix();
  /*! 
    \brief Construct from size.
    \param size Number of rows/columns; must be positive.
  */
  explicit LowerTriangleMatrix(unsigned int size);
  /*! 
    \brief Construct from value and size.
    \param val Value assigned to all matrix elements.
    \param size Number of rows/columns; must be positive.
  */
  LowerTriangleMatrix(T val, unsigned int size);
  /*! 
    \brief Copy constructor
    \param m LowerTriangleMatrix to copy from
  */
  LowerTriangleMatrix(const LowerTriangleMatrix<T> &m);
  /*!
    \brief Copy assignment
    \param rhs Lower triangular matrix assigned from
  */
  LowerTriangleMatrix<T>& operator=(const LowerTriangleMatrix &rhs);

// inspection
  //! Return size.
  unsigned int Size() const;
  //! Return # rows
  unsigned int NRow() const;
  //! Return # columns
  unsigned int NColumn() const;
  //! Return value of j'th element of i'th row
  T operator()(unsigned int i, unsigned int j) const;
  //! Return value of l'th element of k'th diagonal
  T D(int k, unsigned int l) const;
  //! Return k'th diagonal for read-only access.
  ConstSliceRef<T> Diagonal(int k=0) const;

// modification
  //! Return reference to j'th element of i'th row
  T& operator()(unsigned int i, unsigned int j);
  //! Return reference to l'th element of k'th diagonal
  T& D(int k, unsigned int l);
  //! Return k'th diagonal for read/write access.
  SliceRef<T> Diagonal(int k=0);
  //! Multiply matrix by scalar
  LowerTriangleMatrix& operator*=(T a);
  //! Resizing
  LowerTriangleMatrix& Resize(unsigned int size, T val=T())
  {
  TRY
    m_size = size;
    m_elem.resize(Triangle(size, size, size) + 1);
    return *this;
  CATCH(LowerTriangleMatrix::Resize)
  }

private:
// data
  unsigned int m_size;
  valarray<T> m_elem;
}; // LowerTriangleMatrix

/*!
  \brief A UpperTriangleMatrix is an upper triangular matrix

  For this type of matrix we store only the upper triangle.
  All operations maintain shape.

  Note that diagonals are numbered such that the diagonal containing 
  element (i, 0) is numbered -i and the diagonal containing the element 
  (0, j) is numbered j.

  \par Template types
  \arg T The element type
*/
template<class T>
class UpperTriangleMatrix
{
public:
// typedef's
  //! Type definition of element type
  typedef T value_type;

// foundation
  UpperTriangleMatrix();
  /*! 
    \brief Construct from size.
    \param size Number of rows/columns; must be positive.
  */
  explicit UpperTriangleMatrix(unsigned int size);
  /*! 
    \brief Construct from value and size.
    \param val Value assigned to all matrix elements.
    \param size Number of rows/columns; must be positive.
  */
  UpperTriangleMatrix(T val, unsigned int size);
  /*! 
    \brief Copy constructor
    \param m UpperTriangleMatrix to copy from
  */
  UpperTriangleMatrix(const UpperTriangleMatrix<T> &m);
  /*!
    \brief Copy assignment
    \param rhs Upper triangular matrix assigned from
  */
  UpperTriangleMatrix& operator=(const UpperTriangleMatrix &rhs);

// inspection
  //! Return size.
  unsigned int Size() const;
  //! Return # rows
  unsigned int NRow() const;
  //! Return # columns
  unsigned int NColumn() const;
  //! Return value of j'th element of i'th row
  T operator()(unsigned int i, unsigned int j) const;
  //! Return value of l'th element of k'th diagonal
  T D(int k, unsigned int l) const;
  //! Return k'th diagonal for read-only access.
  ConstSliceRef<T> Diagonal(int k=0) const;

// modification
  //! Return reference to j'th element of i'th row
  T& operator()(unsigned int i, unsigned int j);
  //! Return reference to l'th element of k'th diagonal
  T& D(int k, unsigned int l);
  //! Return k'th diagonal for read/write access.
  SliceRef<T> Diagonal(int k=0);
  //! Multiply matrix by scalar
  UpperTriangleMatrix& operator*=(T a);
  //! Resizing
  UpperTriangleMatrix& Resize(unsigned int size, T val=T())
  {
  TRY
    m_size = size;
    m_elem.resize(Triangle(size, size, size) + 1);
    return *this;
  CATCH(UpperTriangleMatrix::Resize)
  }

private:
// data
  unsigned int m_size;
  valarray<T> m_elem;
}; // UpperTriangleMatrix


// *************************** Band matrices **************************
/*!
  \brief A BandMatrix is a band matrix, ie, all elements are zero except in
  a number of diagonal bands around the diagonal.

  For this type of matrix we store only the bands.  All operations maintain shape.

  Note that bands (diagonals) are numbered such that the diagonal containing 
  element (i, 0) is numbered -i and the diagonal containing the element 
  (0, j) is numbered j.

  \par Template types
  \arg T The element type
  \arg L Number of 'lower' diagonals
  \arg U Number of 'upper' diagonals
*/
template<class T, unsigned int L, unsigned int U>
class BandMatrix
{
public:
// types
  //! Type definition of element type
  typedef T value_type;

// foundation
  BandMatrix();
  /*! 
    \brief Construct from size.
    \param size Number of rows/columns; must be positive.
  */
  explicit BandMatrix(unsigned int size);
  /*! 
    \brief Construct from value and size.
    \param val Value assigned to all matrix elements.
    \param size Number of rows/columns; must be positive.
  */
  BandMatrix(T val, unsigned int size);
  /*! 
    \brief Copy constructor
    \param m UpperTriangleMatrix to copy from
  */
  BandMatrix(const BandMatrix<T, L, U> &m);
  /*!
    \brief Copy assignment
    \param rhs Band matrix assigned from
  */
  BandMatrix& operator=(const BandMatrix &rhs);

// inspection
  //! Return size.
  unsigned int Size() const;
  //! Return # rows
  unsigned int NRow() const;
  //! Return # columns
  unsigned int NColumn() const;
  //! Return value of j'th element of i'th row
  T operator()(unsigned int i, unsigned int j) const;
  //! Return value of l'th element of k'th diagonal
  T D(int k, unsigned int l) const;
  //! Return k'th band (diagonal) for read-only access.
  ConstSliceRef<T> Diagonal(int k=0) const;

// modification
  //! Return reference to j'th element of i'th row
  T& operator()(unsigned int i, unsigned int j);
  //! Return reference to l'th element of k'th diagonal
  T& D(int k, unsigned int l);
  //! Return k'th band (diagonal) for read/write access.
  SliceRef<T> Diagonal(int k=0);
  //! Multiply matrix by scalar
  BandMatrix& operator*=(T a);
  //! Resizing
  BandMatrix& Resize(unsigned int size, T val=T())
  {
  TRY
    m_size = size;
    
    return *this;
  CATCH(BandMatrix::Resize)
  }

private:
// data
  unsigned int m_size;
  valarray<T> m_elem;
}; // BandMatrix



/*!
  \brief A TridiagonalMatrix is a band matrix with only three bands:
  the diagonal and the jo adjoining bands.

  For this type of matrix we store only the bands.  All operations maintain shape.

  Note that bands (diagonals) are numbered such that the diagonal containing 
  element (i, 0) is numbered -i and the diagonal containing the element 
  (0, j) is numbered j.

  \par Template types
  \arg T The element type
*/
template<class T>
class TridiagonalMatrix
{
public:
// types
  typedef T value_type;

// foundation
  TridiagonalMatrix();
  /*! 
    \brief Construct from size.
    \param size Number of rows/columns; must be positive.
  */
  explicit TridiagonalMatrix(unsigned int size);
  /*! 
    \brief Construct from value and size.
    \param val Value assigned to all matrix elements.
    \param size Number of rows/columns; must be positive.
  */
  TridiagonalMatrix(T val, unsigned int size);
  /*! 
    \brief Copy constructor
    \param m UpperTriangleMatrix to copy from
  */
  TridiagonalMatrix(const TridiagonalMatrix<T> &m);
  /*!
    \brief Copy assignment
    \param rhs Band matrix assigned from
  */
  TridiagonalMatrix& operator=(const TridiagonalMatrix &rhs);

// inspection
  //! Return size.
  unsigned int Size() const;
  //! Return # rows
  unsigned int NRow() const;
  //! Return # columns
  unsigned int NColumn() const;
  //! Return value of j'th element of i'th row
  T operator()(unsigned int i, unsigned int j) const;
  //! Return value of l'th element of k'th diagonal
  T D(int k, unsigned int l) const;
  //! Return k'th band (diagonal) for read-only access.
  ConstSliceRef<T> Diagonal(int k=0) const;

// modification
  //! Return reference to j'th element of i'th row
  T& operator()(unsigned int i, unsigned int j);
  //! Return reference to l'th element of k'th diagonal
  T& D(int k, unsigned int l);
  //! Return k'th band (diagonal) for read/write access.
  SliceRef<T> Diagonal(int k=0);
  //! Multiply matrix by scalar
  TridiagonalMatrix& operator*=(T a);

private:
// data
  BandMatrix<T, 1, 1> m_m;
}; // TridiagonalMatrix


/*!
  \brief A DiagonalMatrix is a band matrix with only one band: the diagonal.

  For this type of matrix we store only the diagonal.  All operations maintain shape.

  \par Template types
  \arg T The element type
*/
template<class T>
class DiagonalMatrix
{
public:
// types
  typedef T value_type;

// foundation
  DiagonalMatrix();
  /*! 
    \brief Construct from size.
    \param size Number of rows/columns; must be positive.
  */
  explicit DiagonalMatrix(unsigned int size);
  /*! 
    \brief Construct from value and size.
    \param val Value assigned to all matrix elements.
    \param size Number of rows/columns; must be positive.
  */
  DiagonalMatrix(T val, unsigned int size);
  /*! 
    \brief Copy constructor
    \param m UpperTriangleMatrix to copy from
  */
  DiagonalMatrix(const DiagonalMatrix<T> &m);
  /*!
    \brief Copy assignment
    \param rhs Band matrix assigned from
  */
  DiagonalMatrix& operator=(const DiagonalMatrix &rhs);

// inspection
  //! Return size.
  unsigned int Size() const;
  //! Return # rows
  unsigned int NRow() const;
  //! Return # columns
  unsigned int NColumn() const;
  //! Return value of j'th element of i'th row
  T operator()(unsigned int i, unsigned int j) const;
  //! Return value of l'th element of k'th diagonal
  T D(int k, unsigned int l) const;
  //! Return k'th band (diagonal) for read-only access.
  ConstSliceRef<T> Diagonal(int k=0) const;

// modification
  //! Return reference to j'th element of i'th row
  T& operator()(unsigned int i, unsigned int j);
  //! Return reference to l'th element of k'th diagonal
  T& D(int k, unsigned int l);
  //! Return k'th band (diagonal) for read/write access.
  SliceRef<T> Diagonal(int k=0);
  //! Multiply matrix by scalar
  DiagonalMatrix& operator*=(T a);

private:
// data
  BandMatrix<T, 0, 0> m_m;
}; // DiagonalMatrix


// ********************** Helper functions *****************************

/*!
  \brief Return transpose matrix.
  \par Template arguments
  \arg T element type
  \param m Object of class Matrix
  \return Object of class Matrix
*/
template<class T>
Matrix<T> Transpose(const Matrix<T> &m)
{ 
TRY
  Matrix<T> t(m.NColumn(), m.NRow());
  for (unsigned int i = 0; i < m.NRow(); ++i)
    t.Column(i) = m.Row(i);
  return t;
CATCH(Transpose) 
}

/*!
  \brief Return transpose matrix.
  \par Template arguments
  \arg T element type
  \param m Object of class SquareMatrix
  \return Object of class SquareMatrix
*/
template<class T>
SquareMatrix<T> Transpose(const SquareMatrix<T> &m)
{ 
TRY
  SquareMatrix<T> t(m.Size());
  for (unsigned int i = 0; i < m.Size(); ++i)
    t.Column(i) = m.Row(i);
  return t;
CATCH(Transpose) 
}

/*!
  \brief Return transpose matrix.
  \par Template arguments
  \arg T element type
  \param m Object of class LowerTriangleMatrix
  \return Object of class UpperTriangleMatrix
*/
template<class T>
UpperTriangleMatrix<T> Transpose(const LowerTriangleMatrix<T> &l)
{ 
TRY
  UpperTriangleMatrix<T> u(l.Size());
  for (unsigned int i = 0; i < l.Size(); ++i)
  {
    for (unsigned int j = 0; j <= i; ++j)
      u(j, i) = l(i, j);
  }
  return u;
CATCH(Transpose) 
}

/*!
  \brief Return transpose matrix.
  \par Template arguments
  \arg T element type
  \param m Object of class UpperTriangleMatrix
  \return Object of class LowerTriangleMatrix
*/
template<class T>
LowerTriangleMatrix<T> Transpose(const UpperTriangleMatrix<T> &u)
{ 
TRY
  LowerTriangleMatrix<T> l(u.Size());
  for (unsigned int i = 0; i < u.Size(); ++i)
  {
    for (unsigned int j = 0; j <= i; ++j)
      l(i, j) = u(j, i);
  }
  return l;
CATCH(Transpose) 
}


/*!
  \brief Set elements in the (zero'th) diagonal of matrix to elements 
  of vector.
  \par Template arguments
  \arg M Matrix type
  \arg V Vector type
  \param m Matrix
  \param v Vector
  \param k Diagonal to be set.   \par Note that diagonals are numbered such that the diagonal containing 
  element (i, 0) is numbered -i and the diagonal containing the element 
  (0, j) is numbered j.
*/
template<class M, class V>
void SetDiagonal(M &m, const V &v, unsigned int k = 0)
{ 
TRY
  for (unsigned int l=0; l < v.size(); ++l)
    m.D(k, l) = v[l];
CATCH(SetDiagonal) 
}

/*!
  \brief Copy vector u to vector v
  \par Template arguments
  \arg U Vector type
  \arg V Vector type
  \param u Vector
  \param v Vector
  \throw Argument Throw exception if size u and v mismatches.
  \return A reference to v
*/
template<class U, class V>
V& Copy(const U &u, V &v)
{ 
TRY
  ASSERT1(u.size() == v.size(), Argument("Size mismatch."));
  for (unsigned int i=0; i < u.size(); ++i)
    v[i] = u[i];
  return v;
CATCH(Copy)
}

/*!
  \brief Return square matrix of dimension 'size' with diagonal elements set to 'val'
  \par Template arguments
  \arg T Element type
  \param size Size of Matrix
  \param val Value
  \return SquareMatrix
*/
template<class T>
SquareMatrix<T> Diag(unsigned int size, T val)
{ 
TRY
  SquareMatrix<T> m(0.0, size);
  for (unsigned int i = 0; i < size; ++i)
    m(i, i) = val;
  return m;
CATCH(Diag) 
}

/********************* Matrix<> members *******************************/

template<class T> inline
Matrix<T>::Matrix()
{  
}

template<class T> inline
Matrix<T>::Matrix(unsigned int nrow, unsigned int ncol)
: m_nrow(nrow),
  m_ncol(ncol),
  m_elem(nrow * ncol)
{  
}

template<class T> inline
Matrix<T>::Matrix(T val, unsigned int nrow, unsigned int ncol)
: m_nrow(nrow),
  m_ncol(ncol),
  m_elem(val, nrow * ncol)
{
}

template<class T> inline
Matrix<T>::Matrix(const Matrix &m)
: m_nrow(m.NRow()),
  m_ncol(m.NColumn()),
  m_elem(m.m_elem)
{
}      

template<class T> inline
Matrix<T>& Matrix<T>::operator=(const Matrix &rhs)
{ 
TRY
  if (this != &rhs)
  {
    m_elem.resize(rhs.m_elem.size());
    m_elem = rhs.m_elem;
    m_nrow = rhs.m_nrow;
    m_ncol = rhs.m_ncol;
  }
  return *this;
CATCH(Matrix::operator=) 
}

template<class T> inline
unsigned int Matrix<T>::NRow() const
{
  return m_nrow;
}

template<class T> inline
unsigned int Matrix<T>::NColumn() const
{
  return m_ncol;
}

template<class T> inline
T Matrix<T>::operator()(unsigned int i, unsigned int j) const
{ 
TRY
  ASSERT2(i < NRow() && j < NColumn(), Range("Index out of range."));
  return m_elem[i * m_ncol + j];
CATCH(Matrix::operator()) 
}

template<class T> inline
ConstSliceRef<T> Matrix<T>::Row(unsigned int i) const 
{ 
TRY 
  ASSERT2(i < NRow(), Range("Row index out of range."));
  return ConstSliceRef<T>(m_elem, slice(i * m_ncol, m_ncol, 1));
CATCH(Matrix<T>::Row) 
}

template<class T> inline
ConstSliceRef<T> Matrix<T>::Column(unsigned int i) const 
{ 
TRY
  ASSERT2(i < NColumn(), Range("Column index out of range."));
  return ConstSliceRef<T>(m_elem, slice(i, m_nrow, m_ncol));
CATCH(Matrix::Column) 
}

template<class T> inline
T& Matrix<T>::operator()(unsigned int i, unsigned int j)
{ 
TRY
  ASSERT2(i < NRow() && j < NColumn(), Range("Index out of range."));
  return m_elem[i * m_ncol + j];
CATCH(Matrix::operator()) 
}

template<class T> inline
SliceRef<T> Matrix<T>::Row(unsigned int i)
{ 
TRY 
  ASSERT2(i < NRow(), Range("Row index out of range."));
  return SliceRef<T>(m_elem, slice(i * m_ncol, m_ncol, 1));
CATCH(Matrix::Row) 
}

template<class T> inline
SliceRef<T> Matrix<T>::Column(unsigned int i)
{ 
TRY
  ASSERT2(i < NColumn(), Range("Column index out of range."));
  return SliceRef<T>(m_elem, slice(i, m_nrow, m_ncol));
CATCH(Matrix::Column) 
}

template<class T> inline
Matrix<T>& Matrix<T>::operator*=(T a)
{ 
TRY
  m_elem *= a;
  return *this;
CATCH(Matrix::operator*=) 
}

/********************* SquareMatrix<> members *******************************/
template<class T> inline
SquareMatrix<T>::SquareMatrix()
{  
}

template<class T> inline
SquareMatrix<T>::SquareMatrix(unsigned int size)
: m_size(size),
  m_elem(size * size)
{  
}

template<class T> inline
SquareMatrix<T>::SquareMatrix(T val, unsigned int size)
: m_size(size),
  m_elem(val, size * size)
{
}

template<class T> inline
SquareMatrix<T>::SquareMatrix(const SquareMatrix<T> &m)
: m_size(m.Size()),
  m_elem(m.m_elem)
{
}

template<class T> inline
SquareMatrix<T>& SquareMatrix<T>::operator=(const SquareMatrix &rhs)
{ 
TRY
  if (this != &rhs)
  {
    m_elem.resize(rhs.m_elem.size());
    m_elem = rhs.m_elem;
    m_size = rhs.m_size;
  }
  return *this;
CATCH(SquareMatrix::operator=) 
}

template<class T> inline
unsigned int SquareMatrix<T>::Size() const
{
  return m_size;
}

template<class T> inline
unsigned int SquareMatrix<T>::NRow() const
{
  return m_size;
}

template<class T> inline
unsigned int SquareMatrix<T>::NColumn() const
{
  return m_size;
}

template<class T> inline
T SquareMatrix<T>::operator()(unsigned int i, unsigned int j) const
{ 
TRY
  ASSERT2(i < NRow() && j < NColumn(), Range("Index out of range."));
  return m_elem[i * m_size + j];
CATCH(SquareMatrix::operator()) 
}

template<class T> inline
T SquareMatrix<T>::D(int k, unsigned int l) const
{ 
TRY
  ASSERT2(Abs(k) < static_cast<int>(NRow()) && l < NRow() - Abs(k), Range("Index out of range."));
  return (k > 0) ? m_elem[l * m_size + k + l] : m_elem[(l - k) * m_size + l];
CATCH(SquareMatrix::D) 
}

template<class T> inline
ConstSliceRef<T> SquareMatrix<T>::Row(unsigned int i) const 
{ 
TRY 
  ASSERT2(i < NRow(), Range("Row index out of range."));
  return ConstSliceRef<T>(m_elem, slice(i * m_size, m_size, 1));
CATCH(SquareMatrix::Row) 
}

template<class T> inline
ConstSliceRef<T> SquareMatrix<T>::Column(unsigned int i) const 
{ 
TRY
  ASSERT2(i < NColumn(), Range("Column index out of range."));
  return ConstSliceRef<T>(m_elem, slice(i, m_size, m_size));
CATCH(SquareMatrix::Column) 
}

template<class T> inline
ConstSliceRef<T> SquareMatrix<T>::Diagonal(int k) const
{ 
TRY
  ASSERT2(Abs(k) < static_cast<int>(NRow()), Range("Diagonal index out of range."));
  return (k >= 0) ? ConstSliceRef<T>(m_elem, slice(k, m_size - k, m_size + 1))
                  : ConstSliceRef<T>(m_elem, slice(-k * m_size, m_size + k, m_size + 1));
CATCH(SquareMatrix::Diagonal) 
}

template<class T> inline
T& SquareMatrix<T>::operator()(unsigned int i, unsigned int j)
{ 
TRY
  ASSERT1(i < NRow() && j < NColumn(), Range("Index out of range."));
  return m_elem[i * m_size + j];
CATCH(SquareMatrix::operator()) 
}

template<class T> inline
T& SquareMatrix<T>::D(int k, unsigned int l)
{ 
TRY
  ASSERT1(Abs(k) < static_cast<int>(NRow()) && l < NRow() - Abs(k), Range("Index out of range."));
  return (k > 0) ? m_elem[l * m_size + k + l] : m_elem[(l - k) * m_size + l];
CATCH(SquareMatrix::D) 
}

template<class T> inline
SliceRef<T> SquareMatrix<T>::Row(unsigned int i)
{ 
TRY 
  ASSERT2(i < NRow(), Range("Row index out of range."));
  return SliceRef<T>(m_elem, slice(i * m_size, m_size, 1));
CATCH(SquareMatrix::Row) 
}

template<class T> inline
SliceRef<T> SquareMatrix<T>::Column(unsigned int i)
{ 
TRY
  ASSERT2(i < NColumn(), Range("Column index out of range."));
  return SliceRef<T>(m_elem, slice(i, m_size, m_size));
CATCH(SquareMatrix::Column) 
}

template<class T> inline
SliceRef<T> SquareMatrix<T>::Diagonal(int k)
{ 
TRY
  ASSERT1(Abs(k) < static_cast<int>(NRow()), Range("Diagonal index out of range."));
  return (k >= 0) ? SliceRef<T>(m_elem, slice(k, m_size - k, m_size + 1))
                  : SliceRef<T>(m_elem, slice(-k * m_size, m_size + k, m_size + 1));
CATCH(SquareMatrix::Diagonal) 
}

template<class T> inline
SquareMatrix<T>& SquareMatrix<T>::operator*=(T a)
{ 
TRY
  m_elem *= a;
  return *this;
CATCH(SquareMatrix::operator*) 
}

/********************* SymmetricMatrix<> members *******************************/
template<class T> inline
SymmetricMatrix<T>::SymmetricMatrix()
{
}

template<class T> inline
SymmetricMatrix<T>::SymmetricMatrix(unsigned int size)
: m_size(size),
  m_elem(Triangle(size, size, size) + 1)
  
{
}

template<class T> inline
SymmetricMatrix<T>::SymmetricMatrix(T val, unsigned int size)
  : m_size(size),
    m_elem(val, Triangle(size, size, size) + 1)  
{
}

template<class T> inline
SymmetricMatrix<T>::SymmetricMatrix(const SymmetricMatrix<T> &m)
: m_size(m.m_size),
  m_elem(m.m_elem)
{
}

template<class T> inline
SymmetricMatrix<T>& SymmetricMatrix<T>::operator=(const SymmetricMatrix &rhs)
{ 
TRY
  if (this != &rhs)
  {
    m_elem.resize(rhs.m_elem.size());
    m_elem = rhs.m_elem;
    m_size = rhs.m_size;
  }
  return *this;
CATCH(SymmetricMatrix::operator=) 
}

template<class T> inline
unsigned int SymmetricMatrix<T>::Size() const
{
  return m_size;
}

template<class T> inline
unsigned int SymmetricMatrix<T>::NRow() const
{
  return m_size;
}

template<class T> inline
unsigned int SymmetricMatrix<T>::NColumn() const
{
  return m_size;
}

template<class T> inline
T SymmetricMatrix<T>::operator()(unsigned int i, unsigned int j) const
{ 
TRY
  ASSERT2(i < NRow() && j < NColumn(), Range("Index out of range."));
  return (i < j) ? m_elem[Triangle(j - i, i, m_size)] : m_elem[Triangle(i - j, j, m_size)];
CATCH(SymmetricMatrix::operator()) 
}

template<class T> inline
T SymmetricMatrix<T>::D(int k, unsigned int l) const
{ 
TRY
  ASSERT1(Abs(k) < static_cast<int>(NRow()), 
                Range("Index out of range."));
  ASSERT1(l >= 0 && l < NRow() - Abs(k), 
                Range("Index out of range."));
  return m_elem[Triangle(Abs(k), l, m_size)];
CATCH(SymmetricMatrix::D) 
}

template<class T> inline
ConstSliceRef<T> SymmetricMatrix<T>::Diagonal(int k) const
{ 
TRY
  ASSERT1(Abs(k) < NRow(), Range("Index out of range."));
  return ConstSliceRef<T>(m_elem, slice(Triangle(Abs(k), 0, m_size), m_size - Abs(k), 1));
CATCH(SymmetricMatrix::Diagonal) 
}

template<class T> inline
T& SymmetricMatrix<T>::operator()(unsigned int i, unsigned int j)
{ 
TRY
  ASSERT1(i >= 0 && j >= 0 && i < NRow() && j < NColumn(), Range("Index out of range."));
  return (i < j) ? m_elem[Triangle(j - i, i, m_size)] : m_elem[Triangle(i - j, j , m_size)];
CATCH(SymmetricMatrix::operator()) 
}

template<class T> inline
T& SymmetricMatrix<T>::D(int k, unsigned int l)
{ 
TRY
  ASSERT1(Abs(k) < static_cast<int>(NRow()), 
                Range("Index out of range."));
  ASSERT1(l >= 0 && l < NRow() - Abs(k), 
                Range("Index out of range."));
  return m_elem[Triangle(Abs(k), l, m_size)];
CATCH(SymmetricMatrix::D) 
}

template<class T> inline
SliceRef<T> SymmetricMatrix<T>::Diagonal(int k)
{ 
TRY
  ASSERT1(Abs(k) < NRow(), Range("Index out of range."));
  return SliceRef<T>(m_elem, slice(Triangle(Abs(k), 0, m_size), m_size - Abs(k), 1));
CATCH(SymmetricMatrix::Diagonal) 
}

template<class T> inline
SymmetricMatrix<T>& SymmetricMatrix<T>::operator*=(T a)
{ 
TRY
  m_elem *= a;
  return *this;
CATCH(SymmetricMatrix::operator*) 
}

/********************* LowerTriangleMatrix<> members *******************************/
template<class T> inline
LowerTriangleMatrix<T>::LowerTriangleMatrix()
{
}

template<class T> inline
LowerTriangleMatrix<T>::LowerTriangleMatrix(unsigned int size)
: m_size(size),
  m_elem(Triangle(size, size, size) + 1)
  
{
}

template<class T> inline
LowerTriangleMatrix<T>::LowerTriangleMatrix(T val, unsigned int size)
: m_size(size),
  m_elem(val, Triangle(size, size, size) + 1)  
{
}

template<class T> inline
LowerTriangleMatrix<T>::LowerTriangleMatrix(const LowerTriangleMatrix<T> &m)
:   m_size(m.m_size),
    m_elem(m.m_elem)
{}

template<class T> inline
LowerTriangleMatrix<T>& LowerTriangleMatrix<T>::operator=(const LowerTriangleMatrix &rhs)
{ 
TRY
  if (this != &rhs)
  {
    m_elem.resize(rhs.m_elem.size());
    m_elem = rhs.m_elem;
    m_size = rhs.m_size;
  }
  return *this;
CATCH(LowerTriangleMatrix<T>::operator=) 
}

template<class T> inline
unsigned int LowerTriangleMatrix<T>::Size() const
{
  return m_size;
}

template<class T> inline
unsigned int LowerTriangleMatrix<T>::NRow() const
{
  return m_size;
}

template<class T> inline
unsigned int LowerTriangleMatrix<T>::NColumn() const
{
  return m_size;
}

template<class T> inline
T LowerTriangleMatrix<T>::operator()(unsigned int i, unsigned int j) const
{ 
TRY
  ASSERT1(i < NRow() && j < NColumn(), Range("Index out of range."));
  return (i < j) ? 0 : m_elem[Triangle(i - j, j, m_size)];
CATCH(LowerTriangleMatrix<T>::operator()) 
}

template<class T> inline
T LowerTriangleMatrix<T>::D(int k, unsigned int l) const
{ 
TRY
  ASSERT1(Abs(k) < static_cast<int>(NRow()), 
                Range("Index out of range."));
  ASSERT1(l >= 0 && l < NRow() - Abs(k), 
                Range("Index out of range."));
  return (k > 0) ? 0 : m_elem[Triangle(Abs(k), l, m_size)];
CATCH(LowerTriangleMatrix<T>::D) 
}

template<class T> inline
ConstSliceRef<T> LowerTriangleMatrix<T>::Diagonal(int k) const
{ 
TRY
  ASSERT1(k <= 0 && -k < static_cast<int>(NRow()), 
                Range("Index out of range."));
  return ConstSliceRef<T>(m_elem, slice(Triangle(Abs(k), 0, m_size), m_size - Abs(k), 1));
CATCH(LowerTriangleMatrix<T>::Diagonal) 
}

template<class T> inline
T& LowerTriangleMatrix<T>::operator()(unsigned int i, unsigned int j)
{ 
TRY
  ASSERT1(i < NRow() && j <= i, Range("Index out of range."));
  return m_elem[Triangle(i - j, j, m_size)];
CATCH(LowerTriangleMatrix<T>::operator()) 
}

template<class T> inline
T& LowerTriangleMatrix<T>::D(int k, unsigned int l)
{ 
TRY
  ASSERT1(k <= 0 && -k < static_cast<int>(NRow()), 
                Range("Index out of range."));
  ASSERT1(l >= 0 && l < NRow() - Abs(k), 
                Range("Index out of range."));
  return m_elem[Triangle(Abs(k), l, m_size)];
CATCH(LowerTriangleMatrix<T>::D) 
}

template<class T> inline
SliceRef<T> LowerTriangleMatrix<T>::Diagonal(int k)
{ 
TRY
  ASSERT1(k <= 0 && -k < static_cast<int>(NRow()), 
                Range("Index out of range."));
  return SliceRef<T>(m_elem, slice(Triangle(Abs(k), 0, m_size), m_size - Abs(k), 1));
CATCH(LowerTriangleMatrix<T>::Diagonal)
}

template<class T> inline
LowerTriangleMatrix<T>& LowerTriangleMatrix<T>::operator*=(T a)
{ 
TRY
  m_elem *= a;
  return *this;
CATCH(LowerTriangleMatrix<T>::operator*) 
}

/********************* UpperTriangleMatrix<> members *******************************/
template<class T> inline
UpperTriangleMatrix<T>::UpperTriangleMatrix()
{
}

template<class T> inline
UpperTriangleMatrix<T>::UpperTriangleMatrix(unsigned int size)
: m_elem(Triangle(size, size, size) + 1),
  m_size(size)
{
}

template<class T> inline
UpperTriangleMatrix<T>::UpperTriangleMatrix(T val, unsigned int size)
: m_elem(val, Triangle(size, size, size) + 1),
  m_size(size)
{
}

template<class T> inline
UpperTriangleMatrix<T>::UpperTriangleMatrix(const UpperTriangleMatrix<T> &m)
: m_elem(m.m_elem),
  m_size(m.m_size)
{
}

template<class T> inline
UpperTriangleMatrix<T>& UpperTriangleMatrix<T>::operator=(const UpperTriangleMatrix &rhs)
{ 
TRY
  if (this != &rhs)
  {
    m_elem.resize(rhs.m_elem.size());
    m_elem = rhs.m_elem;
    m_size = rhs.m_size;
  }
  return *this;
CATCH(UpperTriangleMatrix<T>::operator=)
}

template<class T> inline
unsigned int UpperTriangleMatrix<T>::Size() const
{
  return m_size;
}

template<class T> inline
unsigned int UpperTriangleMatrix<T>::NRow() const
{
  return m_size;
}

template<class T> inline
unsigned int UpperTriangleMatrix<T>::NColumn() const
{
  return m_size;
}

template<class T> inline
T UpperTriangleMatrix<T>::operator()(unsigned int i, unsigned int j) const
{ 
TRY
  ASSERT1(i >= 0 && j >= 0 && i < NRow() && j < NColumn(), Range("Index out of range."));
  return (i <= j) ? m_elem[Triangle(j - i, i, m_size)] : 0;
CATCH(UpperTriangleMatrix<T>::operator()) 
}

template<class T> inline
T UpperTriangleMatrix<T>::D(int k, unsigned int l) const
{ 
TRY
  ASSERT1(Abs(k) < static_cast<int>(NRow()), 
                Range("Index out of range."));
  ASSERT1(l >= 0 && l < NRow() - Abs(k), 
                Range("Index out of range."));
  return (k < 0) ? 0 : m_elem[Triangle(Abs(k), l, m_size)];
CATCH(UpperTriangleMatrix<T>::D) 
}

template<class T> inline
ConstSliceRef<T> UpperTriangleMatrix<T>::Diagonal(int k) const
{ 
TRY
  ASSERT1(k >= 0 && k < static_cast<int>(NRow()), 
                Range("Index out of range."));
  return ConstSliceRef<T>(m_elem, slice(Triangle(Abs(k), 0, m_size), m_size - Abs(k), 1));
CATCH(UpperTriangleMatrix<T>::Diagonal) 
}

template<class T> inline
T& UpperTriangleMatrix<T>::operator()(unsigned int i, unsigned int j)
{ 
TRY
  ASSERT1(i >= 0 && j >= 0 && i <= j && j < NColumn(), Range("Index out of range."));
  return m_elem[Triangle(j-i, i, m_size)];
CATCH(UpperTriangleMatrix<T>::operator()) 
}

template<class T> inline
T& UpperTriangleMatrix<T>::D(int k, unsigned int l)
{ 
TRY
  ASSERT1(k >= 0 && k < static_cast<int>(NRow()), 
                Range("Index out of range."));
  ASSERT1(l >= 0 && l < NRow() - Abs(k), 
                Range("Index out of range."));
  return m_elem[Triangle(Abs(k), l, m_size)];
CATCH(UpperTriangleMatrix<T>::D) 
}

template<class T> inline
SliceRef<T> UpperTriangleMatrix<T>::Diagonal(int k)
{ 
TRY
  ASSERT1(k >= 0 && k < static_cast<int>(NRow()), 
                Range("Index out of range."));
  return SliceRef<T>(m_elem, slice(Triangle(Abs(k), 0, m_size), m_size - Abs(k), 1));
CATCH(UpperTriangleMatrix<T>::Diagonal) 
}

template<class T> inline
UpperTriangleMatrix<T>& UpperTriangleMatrix<T>::operator*=(T a)
{ 
TRY
  m_elem *= a;
  return *this;
CATCH(UpperTriangleMatrix<T>::operator*) 
}

/********************* BandMatrix<> members *******************************/
template<class T, unsigned int L, unsigned int U> inline
BandMatrix<T, L, U>::BandMatrix()
{
}

template<class T, unsigned int L, unsigned int U> inline
BandMatrix<T, L, U>::BandMatrix(unsigned int size)
: m_elem((U + L + 1) * size),
  m_size(size)
{
}

template<class T, unsigned int L, unsigned int U> inline
BandMatrix<T, L, U>::BandMatrix(T val, unsigned int size)
: m_elem(val, (U + L + 1) * size),
  m_size(size)
{
}

template<class T, unsigned int L, unsigned int U> inline
BandMatrix<T, L, U>::BandMatrix(const BandMatrix<T, L, U> &m)
: m_elem(m.m_elem),
  m_size(m.m_size)
{
}

template<class T, unsigned int L, unsigned int U> inline
BandMatrix<T, L, U>& BandMatrix<T, L, U>::operator=(const BandMatrix &rhs)
{ 
TRY
  if (this != &rhs)
  {
    m_elem.resize(rhs.m_elem.size());
    m_elem = rhs.m_elem;
    m_size = rhs.m_size;
  }
  return *this;
CATCH(BandMatrix::operator=) 
}

template<class T, unsigned int L, unsigned int U> inline 
unsigned int BandMatrix<T, L, U>::Size() const
{
  return m_size;
}

template<class T, unsigned int L, unsigned int U> inline
unsigned int BandMatrix<T, L, U>::NRow() const
{
  return m_size;
}

template<class T, unsigned int L, unsigned int U> inline
unsigned int BandMatrix<T, L, U>::NColumn() const
{
  return m_size;
}

template<class T, unsigned int L, unsigned int U>
T BandMatrix<T, L, U>::operator()(unsigned int i, unsigned int j) const
{ 
TRY
  ASSERT1(i >= 0 && j >= 0 && i < NRow() && j < NColumn(), Range("Index out of range."));
  int k = static_cast<int>(j) - static_cast<int>(i);
  if (-k > static_cast<int>(L) || k > static_cast<int>(U))
    return 0;
  return (k > 0) ? m_elem[(k + L) * m_size + i] : m_elem[(k + L) * m_size + j];
CATCH(BandMatrix::operator()) 
}

template<class T, unsigned int L, unsigned int U>
T BandMatrix<T, L, U>::D(int k, unsigned int l) const
{ 
TRY
  ASSERT1(Abs(k) < NRow(), 
                Range("Index out of range."));
  ASSERT1(l >= 0 && l < NRow() - Abs(k), 
                Range("Index out of range."));
  return (-k > static_cast<int>(L) || k > static_cast<int>(U)) ? 
                                     0 : m_elem[(k + L) * m_size + l];
CATCH(BandMatrix::D) 
}

template<class T, unsigned int L, unsigned int U> inline
ConstSliceRef<T> BandMatrix<T, L, U>::Diagonal(int k) const
{ 
TRY
  ASSERT1(-k <= static_cast<int>(L) && k <= static_cast<int>(U), 
             Range("Index out of range."));     return ConstSliceRef<T>(m_elem, slice((k + L) * m_size, m_size - Abs(k), 1));
CATCH(BandMatrix::Diagonal) 
}

template<class T, unsigned int L, unsigned int U>
T& BandMatrix<T, L, U>::operator()(unsigned int i, unsigned int j)
{ 
TRY
  ASSERT1(i >= 0 && j >= 0 && i < NRow() && j < NColumn(), Range("Index out of range."));
  int k = static_cast<int>(j) - static_cast<int>(i);
  ASSERT1(-k <= static_cast<int>(L) && k <= static_cast<int>(U), 
             Range("Index out of range."));     return (k > 0) ? m_elem[(k + L) * m_size + i] : m_elem[(k + L) * m_size + j];
CATCH(BandMatrix::operator()) 
}

template<class T, unsigned int L, unsigned int U> inline
T& BandMatrix<T, L, U>::D(int k, unsigned int l)
{ 
TRY
  ASSERT1(-k <= static_cast<int>(L) && k <= static_cast<int>(U), 
             Range("Index out of range."));     ASSERT1(l >= 0 && l < NRow() - Abs(k), 
                Range("Index out of range."));
  return m_elem[(k + L) * m_size + l];
CATCH(BandMatrix::D) 
}

template<class T, unsigned int L, unsigned int U> inline
SliceRef<T> BandMatrix<T, L, U>::Diagonal(int k)
{ 
TRY
  ASSERT1(-k <= static_cast<int>(L) && k <= static_cast<int>(U), 
             Range("Index out of range."));     return SliceRef<T>(m_elem, slice((k + L) * m_size, m_size - Abs(k), 1));
CATCH(BandMatrix::Diagonal) 
}

template<class T, unsigned int L, unsigned int U> inline
BandMatrix<T, L, U>& BandMatrix<T, L, U>::operator*=(T a)
{ 
TRY
  m_elem *= a;
  return *this;
CATCH(BandMatrix::operator*) 
}

// ********************* TridiagonalMatrix<> members *******************************
template<class T> inline
TridiagonalMatrix<T>::TridiagonalMatrix()
{
}

template<class T> inline
TridiagonalMatrix<T>::TridiagonalMatrix(unsigned int size)
: m_m(size)
{
}

template<class T> inline
TridiagonalMatrix<T>::TridiagonalMatrix(T val, unsigned int size)
: m_m(val, size)
{
}

template<class T> inline
TridiagonalMatrix<T>::TridiagonalMatrix(const TridiagonalMatrix<T> &m)
: m_m(m.m_m)
{
}

template<class T> inline
TridiagonalMatrix<T>& TridiagonalMatrix<T>::operator=(const TridiagonalMatrix &rhs)
{ 
TRY
  if (this != &rhs)
    m_m = rhs.m_m;
  return *this;
CATCH(TridiagonalMatrix<T>::operator=) 
}

template<class T> inline
unsigned int TridiagonalMatrix<T>::Size() const
{
  return m_m.Size();
}

template<class T> inline
unsigned int TridiagonalMatrix<T>::NRow() const
{
  return m_m.NRow();
}

template<class T> inline
unsigned int TridiagonalMatrix<T>::NColumn() const
{
  return m_m.NColumn();
}

template<class T> inline
T TridiagonalMatrix<T>::operator()(unsigned int i, unsigned int j) const
{ 
TRY
  return m_m(i, j);
CATCH(TridiagonalMatrix::operator()) 
}

template<class T> inline
T TridiagonalMatrix<T>::D(int k, unsigned int l) const
{ 
TRY
  return m_m.D(k, l);
CATCH(TridiagonalMatrix::D) 
}

template<class T> inline
ConstSliceRef<T> TridiagonalMatrix<T>::Diagonal(int k) const
{ 
TRY
  return m_m.Diagonal(k);
CATCH(TridiagonalMatrix::Diagonal) 
}

template<class T> inline
T& TridiagonalMatrix<T>::operator()(unsigned int i, unsigned int j)
{ 
TRY
  return m_m(i, j);
CATCH(TridiagonalMatrix::operator()) 
}

template<class T> inline
T& TridiagonalMatrix<T>::D(int k, unsigned int l)
{ 
TRY
  return m_m.D(k, l);
CATCH(TridiagonalMatrix::D) 
}

template<class T> inline
SliceRef<T> TridiagonalMatrix<T>::Diagonal(int k)
{ 
TRY
  return m_m.Diagonal(k);
CATCH(TridiagonalMatrix::Diagonal) 
}

template<class T> inline
TridiagonalMatrix<T>& TridiagonalMatrix<T>::operator*=(T a)
{ 
TRY
  m_m *= a;
  return *this;
CATCH(TridiagonalMatrix::operator*) 
}

// ******************** DiagonalMatrix<> members *******************************
template<class T> inline
DiagonalMatrix<T>::DiagonalMatrix()
{
}

template<class T> inline
DiagonalMatrix<T>::DiagonalMatrix(unsigned int size)
: m_m(size)
{
}

template<class T> inline
DiagonalMatrix<T>::DiagonalMatrix(T val, unsigned int size)
: m_m(val, size)
{
}

template<class T> inline
DiagonalMatrix<T>::DiagonalMatrix(const DiagonalMatrix<T> &m)
: m_m(m.m_m)
{
}

template<class T> inline
DiagonalMatrix<T>& DiagonalMatrix<T>::operator=(const DiagonalMatrix &rhs)
{ 
TRY
  if (this != &rhs)
    m_m = rhs.m_m;
  return *this;
CATCH(DiagonalMatrix::operator=) 
}

template<class T> inline
unsigned int DiagonalMatrix<T>::Size() const
{
  return m_m.Size();
}

template<class T> inline
unsigned int DiagonalMatrix<T>::NRow() const
{
  return m_m.NRow();
}

template<class T> inline
unsigned int DiagonalMatrix<T>::NColumn() const
{
  return m_m.NColumn();
}

template<class T> inline
T DiagonalMatrix<T>::operator()(unsigned int i, unsigned int j) const
{ 
TRY
  return m_m(i, j);
CATCH(DiagonalMatrix::operator()) 
}

template<class T> inline
T DiagonalMatrix<T>::D(int k, unsigned int l) const
{ 
TRY
  return m_m.D(k, l);
CATCH(DiagonalMatrix::D) 
}

template<class T> inline
ConstSliceRef<T> DiagonalMatrix<T>::Diagonal(int k) const
{ 
TRY
  return m_m.Diagonal(k);
CATCH(DiagonalMatrix::Diagonal) 
}

template<class T> inline
T& DiagonalMatrix<T>::operator()(unsigned int i, unsigned int j)
{ 
TRY
  return m_m(i, j);
CATCH(DiagonalMatrix::operator()) 
}

template<class T> inline
T& DiagonalMatrix<T>::D(int k, unsigned int l)
{ 
TRY
  return m_m.D(k, l);
CATCH(DiagonalMatrix::D) 
}

template<class T> inline
SliceRef<T> DiagonalMatrix<T>::Diagonal(int k)
{ 
TRY
  return m_m.Diagonal(k);
CATCH(DiagonalMatrix::Diagonal) 
}

template<class T> inline
DiagonalMatrix<T>& DiagonalMatrix<T>::operator*=(T a)
{ 
TRY
  m_m *= a;
  return *this;
CATCH(DiagonalMatrix::operator*) 
}

//} // anonymous namespace

///////////////////////////////////////////////// LINALG.H + .CPP ///////////////////////////////////////

/*!
  \file linalg.h
  \brief Linear algebra functionality
*/

using std::vector;

/*! 
  \brief Return scalar product u.v

  \arg U 'u's type
  \arg V 'v's type

  \par Note: Template type U and V must suport size and operator[size_t i]. In
  additon, U operator[] multiplied V operator[] must be compatible with double.

  \param u A vector
  \param v A vector
  \param begin Element to begin scalar product.
  \param end Element after end of scalar product.


  \throw Argument Throw exception if start and end is outside of 'u's and 'v's range.

  \return scalar product of 'u' and 'v'.
*/
template<class U, class V>
inline double Dot(const U &u, const V &v, size_t begin, size_t end)
{
  ASSERT1(end <= Min(u.size(), v.size()), Argument("Dot:Illegal begin/end."));
  double ret = 0.0;
  for (size_t i = begin; i < end; ++i)
  {
    ret += u[i] * v[i];
  }
  return ret;
}

/*! 
  \brief Return scalar product u.v

  \arg U 'u's type
  \arg V 'v's type

  \par Note: Template type U and V must suport size and operator[size_t i]. In
  additon, U operator[] multiplied V operator[] must be compatible with double.

  \param u A vector
  \param v A vector

  \throw Argument Throw exception if u's size larger than 'v's size.

  \return scalar product of 'u' and 'v'.
*/
template<class U, class V>
inline double Dot(const U &u, const V &v)
{
  return Dot(u, v, 0, u.size());
}

/*! 
  \brief Adds Mu to v

  \arg M 'm's type
  \arg U 'u's type
  \arg V 'v's type

  \par Note: Template type M must suport NRow, NColumn and operator(size_t i, size_t j).  Template type U and V must suport size and operator[size_t i]. In
  additon, U operator[] multiplied M operator() must be compatible with V operator[] type.

  \param m A matrix
  \param u A vector
  \param v A vector 

  \throw Argument Throw exception if m.NColumn and u.size or m.NRow and v.size mismatches.
*/
template<class M, class U, class V>
inline void MVMultAdd(const M &m, const U &u, V &v)
{
TRY
  ASSERT1(m.NColumn() == u.size() && m.NRow() == v.size(), Argument("Size mismatch."));
  for (size_t i=0; i < v.size(); ++i)
  {
    for (size_t j=0; j < u.size(); ++j)
      v[i] += m(i, j) * u[j];
  }
CATCH(basnum::MVMultAdd);
}

/*! 
  \brief Returns Mu

  \arg M 'm's type
  \arg U 'u's type

  \par Note: Template type M must suport NRow, NColumn and operator(size_t i, size_t j).  Template type U must suport size, operator[size_t i] and value_type. In
  additon, M operator[] multiplied U operator[] must be compatible with M' value_type.

  \param m A matrix
  \param u A vector

  \return a valarray of M's valuetype.

  \throw Argument Throw exception if m.NColumn and u.size or m.NRow and v.size mismatches.
*/
template<class M, class U>
inline valarray<typename U::value_type> MVMult(const M &m, const U &u)
{
TRY
  // Here we use M::value_type to help the compiler catch any inconsistency
  valarray<typename M::value_type> v(0.0, m.NRow());
  MVMultAdd(m, u, v);
  return v;
CATCH(basnum::MVMult)
}

/*! 
  \brief Adds uM to v

  \arg M 'm's type
  \arg U 'u's type
  \arg V 'v's type

  \par Note: Template type M must suport NRow, NColumn and operator(size_t i, size_t j).  Template type U and V must suport size and operator[size_t i]. In
  additon, U operator[] multiplied M operator() must be compatible with V operator[] type.

  \param m A matrix
  \param u A vector
  \param v A vector 

  \throw Argument Throw exception if m.NColumn and v.size or m.NRow and u.size mismatches.
*/
template<class U, class M, class V>
inline void VMMultAdd(const U &u, const M &m, V &v)
{
TRY
  ASSERT1(m.NRow() == u.size() && m.NColumn() == v.size(), Argument("Size mismatch."));
  for (size_t i=0; i < v.size(); ++i)
  {
    for (size_t j=0; j < u.size(); ++j)
    {
      v[i] += u[j] * m(j, i);
    }
  }
CATCH(basnum::VMMultAdd);
}

/*! 
  \brief Returns uM

  \arg M 'm's type
  \arg U 'u's type

  \par Note: Template type M must suport NRow, NColumn, operator(size_t i, size_t j) 
  and value_type.  Template type U must suport size and operator[size_t i]. In
  additon, M operator[] multiplied U operator[] must be compatible with M' value_type.

  \param m A matrix
  \param u A vector

  \return a valarray of M's valuetype.

  \throw Argument Throw exception if m.NRow and u.size mismatches.
*/
template<class U, class M>
inline valarray<typename U::value_type> VMMult(const U &u, const M &m)
{
TRY
  // Here we use M::value_type to help the compiler catch any inconsistency
  valarray<typename M::value_type> v(0.0, m.NColumn());
  VMMultAdd(u, m, v);
  return v;
CATCH(basnum::VMMult)
}

/*! 
  \brief Evaluate quadratic form, ie, return uMv

  \arg M 'm's type
  \arg U 'u's type
  \arg V 'v's type

  \par Note: Template type M must suport NRow, NColumn and operator(size_t i, size_t j).  Template type U and V must suport size and operator[size_t i]. In
  additon, arithmetic operations (+,*) on V operator[] and M operator() must be 
  compatible with U' value_type.

  \param m A matrix
  \param u A vector
  \param v A vector 

  \returns uMv im 'U's value_type 

  \throw Argument Throw exception if m.NColumn and u.size or m.NRow and v.size mismatches.
*/
template<class U, class M, class V>
inline typename U::value_type UMVMult(const U &u, const M&m, const V&v)
{
TRY
  ASSERT1(m.NRow() == u.size() && m.NColumn() == v.size(), Argument("Size mismatch."));
  typename M::value_type res = 0;
  for (size_t i=0; i < u.size(); ++i)
  {
    typename M::value_type r = 0;
    for (size_t j=0; j < v.size(); ++j)
      r += m(i, j) * v[j];
    res += u[i] * r;
  }
  return res;
CATCH(UMVMult)
}

/*! 
  \brief Add ab to c

  \arg A 'a's type
  \arg B 'b's type
  \arg C 'c's type

  \par Note: Template type A, B and C must suport NRow, NColumn and 
  operator(size_t i, size_t j). In
  additon, arithmetic operations (+,*) on A operator() and B operator() must be 
  compatible with C operator().

  \param a A matrix
  \param b A matrix
  \param c A matrix 

  \returns uMv im 'U's value_type 

  \throw Argument Throw exception if a.NRow and c.NRow, b.NColumn and  c.NColumn or a.NColumn and b.NRow mismatches.
*/
template<class A, class B, class C>
inline void MMultAdd(const A &a, const B &b, C &c)
{
  // Consistency checks
  ASSERT1(a.NRow() == c.NRow() && b.NColumn() == c.NColumn(), Argument("Size mismatch."));
  ASSERT1(a.NColumn() == b.NRow(), Argument("Size mismatch."));

  // Compute
  for (size_t i=0; i < c.NRow(); ++i)
  {
    for (size_t j=0; j < c.NColumn(); ++j)
    {
      typename A::value_type t = 0;
      for (size_t k=0; k < a.NColumn(); ++k)
        t += a(i, k) * b(k, j);
      c(i, j) += t;
    }
  }
}

/*! 
  \brief Solves x in the equation \f$ tdm x = y \f$

  \arg T value_type of TridiagonalMatrix tdm
  \arg Y 'y's type
  \arg X 'x's type

  \par Note: Template type Y and X must suport size and 
  operator [size_t i]. In
  additon, Y's operator[], X's operator and T must suport all arithmetic 
  operations (-,/+,*) with return type value_type X.

  \param tdm TridiagonalMatrix with value_type x
  \param y rhs vector
  \param x lhs vector (the solution vector)

  \throw Numerical Throw exception if the equation has more than one solution or 
        first element in diagonal is 0.0.
*/
template<class T, class Y, class X>
inline void SolveMatrixEquation(const TridiagonalMatrix<T>& tdm,
                                const Y& y,
                                X& x,
                                valarray<typename X::value_type>& tmp)
{
TRY
  // consistency
  ASSERT1(tdm(0, 0) != 0.0, Numerical("Cannot solve."));

  // get the three "diagonals"
  ConstSliceRef<T> u = tdm.Diagonal(1);
  ConstSliceRef<T> d = tdm.Diagonal(0);
  ConstSliceRef<T> l = tdm.Diagonal(-1);

  // initialize
  T denom = d(0);
  x[0] = y[0] / denom;

  // decomposition and forward substitution
  size_t i;
  unsigned int sz = x.size();
  for (i = 1; i < sz; ++i) 
  { 
    tmp[i] = u(i-1) / denom;
    denom = d(i) - l(i-1) * tmp[i];
    ASSERT1(denom != 0.0, Numerical("Cannot solve."));
    x[i] = (y[i] - l(i-1) * x[i-1]) / denom;
  }

  // backsubstitution
  for (i = sz - 1; i > 0; --i)
    x[i-1] -= tmp[i] * x[i]; 
CATCH(basnum::SolveMatrixEquation)
}

/*! 
  \brief Solves x in the equation \f$ tdm x = y \f$

  \arg T value_type of TridiagonalMatrix tdm
  \arg Y 'y's type
  \arg X 'x's type

  \par Note: Template type Y and X must suport size and 
  operator [size_t i]. In
  additon, Y's operator[], X's operator and T must suport all arithmetic 
  operations (-,/+,*) with return type value_type X.

  \param tdm TridiagonalMatrix with value_type x
  \param y rhs vector
  \param x lhs vector (the solution vector)

  \throw Argument Throw exception if x.size outside 'y' or 'tdm' range
  \throw Numerical Throw exception if the equation has more than one solution or 
  first element in diagonal is 0.0.*/
template<class T, class Y, class X>
inline void SolveMatrixEquation(const TridiagonalMatrix<T> &tdm,
                                const Y &y,
                                X &x)
{
  // work space
  valarray<typename X::value_type> tmp(x.size());

  SolveMatrixEquation(tdm, y, x, tmp);
}

/*! 
  \brief Cholesky decomposition of positive definite, real, symmetric matrix.
  \param d The input matrix
  \param a The lower triagonal output 

  On output a is the left factor in the decomposition. Ie, the product of
  a with its (upper triangular) transpose is d.
  
  \throw Argument The function throws an argument exception if
  the inputs are of unequal sizes.
  \throw Numerical The function throws a numerical exception iff the input 
  matrix is not positive definite.
*/
void Cholesky(const SymmetricMatrix<double> &d, 
              LowerTriangleMatrix<double> &a);
void Cholesky(const SymmetricMatrix<double> &d, 
              LowerTriangleMatrix<float> &a);

/*! 
  \brief Orthogonal decomposition of symmetric matrix

  \param s      Symmetric matrix
  \param maxit  Maximal number of iterations
  \param o      Orthogonal matrix
  \param ev     Eigenvalues (output)
  \param tmp    Work space

  \par Note that we MUST have o.Size() == ev.size() == tmp.size() and o.Size() <= s.Size().
       If o.Size() < s.Size() only the submatrix in the upper left corner of s
       is decomposed.
*/
void OrthogonalDecomposition(const SymmetricMatrix<double>& s, 
                             size_t maxit, 
                             SquareMatrix<double>& o,
                             valarray<double>& ev,
                             valarray<double>& tmp);

/*! 
  \brief Orthogonal decomposition of symmetric matrix

  \param s      Symmetric matrix
  \param maxit  Maximal number of iterations
  \param o      Orthogonal matrix
  \param ev     Eigenvalues (output)

  \par Note that we MUST have o.Size() == ev.size() and o.Size() <= s.Size().
       If o.Size() < s.Size() only the submatrix in the upper left corner of s
       is decomposed.
*/
void OrthogonalDecomposition(const SymmetricMatrix<double>& s, 
                             size_t maxit, 
                             SquareMatrix<double>& o,
                             valarray<double>& ev);

void Matrix_Sym2Inv(SymmetricMatrix<double>& a, size_t maxit=30, size_t rank=0);

void OrthogonalDecomposition_Sort(const SymmetricMatrix<double>& s, 
                                  size_t maxit, 
                                  SquareMatrix<double>& o,
                                  valarray<double>& ev,
                                  valarray<double>& tmp);

void OrthogonalDecomposition_Sort(const SymmetricMatrix<double>& s, 
                                  size_t maxit, 
                                  SquareMatrix<double>& o,
                                  valarray<double>& ev);

void Matrix_TD2OD(valarray<double>& diag, valarray<double>& off, SquareMatrix<double>& a, size_t maxit);

/*!
  Householder reduction of a real, symmetric matrix 's'. On output, 'a' is replaced
  by the orthogonal matrix Q effecting the transformation. 'diag' returns the diagonal elements
  of the tridiagonal matrix, and 'off' the off-diagonal elements, with.
*/
void Symmetric2Tridagonal(const SymmetricMatrix<double>& s,
                          SquareMatrix<double>& a, 
                          valarray<double>& diag, 
                          valarray<double>& off);

void Symmetric2Tridagonal(SquareMatrix<double>& a, 
                          valarray<double>& diag, 
                          valarray<double>& off);

LowerTriangleMatrix<double>	LTInvert(const LowerTriangleMatrix<double>& lt);

void PCConstruct(SymmetricMatrix<double>& a, size_t n);


void Factors(const SymmetricMatrix<double> &c, 
                     Matrix<double> &a,
                     const double tol = 0.99,
                     const unsigned int maxIter = 100);

// Inversion of symmetric matrix by orthogonalisation and inversion of eigenvalues.
// For a singular matrix we use the given rank to determine how many eigenvalues to invert.
// If the given rank is '0' we invert all eigenvalues
void Matrix_Sym2Inv(SymmetricMatrix<double>& a, size_t maxit, size_t rank)
{
  size_t r = (rank == 0) ? a.NRow() : rank;

  // get work space
  valarray<double> ev(a.NRow()); // eigenvalues
  SquareMatrix<double> v(a);        // eigenvectors

  // find eigenvectors and -values
  OrthogonalDecomposition(a, maxit, v, ev);

  // construct inverse
  for (size_t k=0; k < r; k++)
    ev[k] = 1 / ev[k];
  for (size_t i=0; i < a.NRow(); i++)
    for (size_t j=0; j < a.NRow(); j++)         {
      a(i,j) = 0.0;
      for (size_t k=0; k < r; k++)
        a(i,j) += ev[k] * v(i,k) * v(j,k) ;
    }
}    

void OrthogonalDecomposition(const SymmetricMatrix<double>& s, 
                                     size_t maxit, 
                                     SquareMatrix<double>& o,
                                     valarray<double>& ev)
{
  // get work space
	valarray<double> tmp(ev.size());
	
  OrthogonalDecomposition(s, maxit, o, ev, tmp);
}

void OrthogonalDecomposition(const SymmetricMatrix<double>& s, 
                                     size_t maxit, 
                                     SquareMatrix<double>& o,
                                     valarray<double>& ev,
                                     valarray<double>& tmp)
{	
	// put matrix on tridiagonal form 
	Symmetric2Tridagonal(s, o, ev, tmp) ;
	
	// find eigenvalues and orthogonal matrix 
	Matrix_TD2OD(ev, tmp, o, maxit) ;
}

// orthogonal decomposition of symmetric matrix with sorting of eigenvalues and -vectors 
//	in order of decreasing eigenvalue
void OrthogonalDecomposition_Sort(const SymmetricMatrix<double>& s, 
                                          size_t maxit, 
                                          SquareMatrix<double>& o,
                                          valarray<double>& ev)
{
  // get work space
	valarray<double> tmp(ev.size());

  OrthogonalDecomposition_Sort(s, maxit, o, ev, tmp);
}


// orthogonal decomposition of symmetric matrix with sorting of eigenvalues and -vectors 
//	in order of decreasing eigenvalue
void OrthogonalDecomposition_Sort(const SymmetricMatrix<double>& s, 
                                          size_t maxit, 
                                          SquareMatrix<double>& o,
                                          valarray<double>& ev,
                                          valarray<double>& tmp)
{
	// do OD
  OrthogonalDecomposition(s, maxit, o, ev, tmp);
	
	  // rearrange eigenvectors  
	for (size_t l = 0; l < o.NRow(); l++)
  {
    size_t k=l;
    double v = ev[k] ;
    for (size_t j = l + 1; j < o.NRow(); j++)
      if (ev[j] >= v)
        v = ev[(k=j)] ;
    if (k != l)
    {
      ev[k] = ev[l] ;
      ev[l] = v ;
      for (size_t j = 0; j < o.NRow(); j++)
      {
        v = o(j,l) ;
        o(j,l) = o(j,k) ;
        o(j,k) = v ;
      }
    }
  }
}


// orthogonal decomposition of tri-diagonal matrix
void Matrix_TD2OD(valarray<double>& diag, valarray<double>& off, SquareMatrix<double>& a, size_t maxit)
{
TRY
  size_t dim = off.size();
	for (size_t j = 1; j < dim; j++)
		off[j - 1] = off[j] ;

	off[dim - 1] = 0.0 ;
	
	for (int k = 0; static_cast<size_t>(k) < dim; ++k)
	{
    size_t i;
    do
		{
			for (i = k; i < dim - 1; i++)
			{
        double y = Abs(diag[i]) + Abs(diag[i+1]) ;
				if (Abs(off[i]) + y == y) 
					break ;
			}
			// i == dim - 1  always !! 
			if (i != static_cast<size_t>(k))
			{
//				if (cnt++ == maxit)
  //        throw Numerical("Maximal number of iterations reached.");
				double ttmp = (diag[k + 1] - diag[k]) / (2.0 * off[k]) ;
				double temp = sqrt(Square(ttmp) + 1.0) ;				
        double zz = (ttmp < 0) ? -Abs(temp) : Abs(temp) ;
				ttmp = diag[i] - diag[k] + off[k] / (ttmp + zz) ;
				double tmp = 1.0 ;
        double z = 1.0;
				double ttemp = 0.0 ;
				for (int j = static_cast<size_t>(i - 1); j >= k; --j)
				{
					double x = tmp * off[j] ;
					double xx = z * off[j] ;
					if (Abs(x) >= Abs(ttmp))
					{
						z = ttmp / x ;
						temp = sqrt(Square(z) + 1.0) ;
						off[j + 1] = x * temp ;
						z *= (tmp = 1.0 / temp) ;
					}
					else
					{
						tmp = x / ttmp;
						temp = sqrt(Square(tmp) + 1.0) ;
						off[j + 1] = ttmp * temp ;
						tmp *= (z = 1.0 / temp) ;
					}
					ttmp = diag[j + 1] - ttemp ;
					temp = (diag[j] - ttmp) * tmp + 2.0 * z * xx ;
					ttemp = tmp * temp ;
					diag[j + 1] = ttmp + ttemp ;
					ttmp = z * temp - xx ;
					//
					for (size_t n = 0; n < dim; n++)
					{
						x = a(n,j + 1) ;
						a(n,j + 1) = tmp * a(n,j) + z * x ;
						a(n,j) = z * a(n,j) - tmp * x ;
					}
				}
				diag[k] = diag[k] - ttemp ;
				off[k] = ttmp ;
				off[i] = 0.0 ;
			}
		}	
		while (i != static_cast<size_t>(k)) ;
	}
CATCH(MatrixTD2OD)
}
					 				
		
		
void Symmetric2Tridagonal(const SymmetricMatrix<double>& s,
                                  SquareMatrix<double>& a, 
                                  valarray<double>& diag, 
                                  valarray<double>& off)
{
  a = s;
  Symmetric2Tridagonal(a, diag, off);
}

void Symmetric2Tridagonal(SquareMatrix<double>& a, 
                                  valarray<double>& diag, 
                                  valarray<double>& off)
{
  size_t dim = off.size();
	for (int j = static_cast<int>(dim - 1); j > 0; --j)
	{
		size_t n = static_cast<size_t>(j) - 1 ;
		double temp = 0.0;
    double fac = 0.0 ;
		if (n > 0)
		{
			for (size_t m = 0; m <= n; m++) 
				fac += Abs(a(j,m)) ;
			if (fac == 0.0)
				off[j] = a(j,n) ;
			else
			{
				for (size_t m = 0; m <= n; m++)
				{
					a(j,m) /= fac ;
					temp += a(j,m) * a(j,m);
				}
				double tmp = a(j,n);
				double ttmp = tmp > 0 ? -sqrt(temp) : sqrt(temp) ;
				off[j] = fac * ttmp ;
				temp -= tmp * ttmp ;
				a(j,n) = tmp - ttmp ;
				tmp = 0.0 ;
				for (size_t i = 0; i <= n; i++)
				{
					a(i,j) = a(j,i) / temp ;
					ttmp = 0.0 ;
					for (size_t m = 0; m <= i; m++)
						ttmp += a(i,m) * a(j,m) ;
					for (size_t m = i + 1; m <= n; m++)
						ttmp += a(m,i) * a(j,m) ;
					off[i] = ttmp / temp ;
					tmp += off[i] * a(j,i) ;	
				}
				double ttemp = tmp / (2.0 * temp) ;
				for (size_t i = 0; i <= n; i++)
				{
					tmp = a(j,i) ;
					off[i] = ttmp = off[i] - ttemp * tmp ;
					for (size_t m = 0; m <= i; m++)
						a(i,m) -= (tmp * off[m] + ttmp * a(j,m)) ;
				}
			}
		}	
		else
			off[j] = a(j,n) ;
		diag[j] = temp ;
	}

	diag[0] = 0.0 ;
	off[0] = 0.0 ;
	
	for (size_t j = 0; static_cast<size_t>(j) < dim; ++j)
	{
		int n = j - 1 ;
		if (diag[j])
		{
			for (int i = 0; i <= n; i++)
			{
				double ttmp = 0.0 ;
				for (int m = 0; m <= n; m++)
					ttmp += a(j,m) * a(m,i) ;
				for (int m = 0; m <= n; m++)
					a(m,i) -= ttmp * a(m,j) ;
			}
		}
		diag[j] = a(j,j) ;
		a(j,j) = 1.0 ;
		for (int i = 0; i <= n; i++)
			a(i,j) = a(j,i) = 0.0 ;
	}
}				


/**************************************************/

/*
  Householder reduction of a real, symmetric matrix a[1..n][1..n]. On output, a is replaced
  by the orthogonal matrix Q effecting the transformation. diag[1..n] returns the diagonal elements
  of the tridiagonal matrix, and off[1..n] the off-diagonal elements, with off[1]=0. Several
  statements, as noted in comments, can be omitted if only eigenvalues are to be found, in which
  case a contains no useful information on output. Otherwise they are to be included.

void Symmetric(SymmetricMatrix<double>& a, valarray<double> diag, valarray<double> off)
{
  int l,k,j,i;
  float scale,hh,h,g,f;

  for (i=n;i>=2;i--) 
  {
    l=i-1;
    h=scale=0.0;
    if (l > 1) 
    {
      for (k=1;k<=l;k++)
        scale += fabs(a[i][k]);
      if (scale == 0.0)
        off[i]=a[i][l];
      else 
      {
        for (k=1;k<=l;k++) 
        {
          a[i][k] /= scale;
          h += a[i][k]*a[i][k]; 
        }
        f=a[i][l];
        g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
        off[i]=scale*g;
        h -= f*g;
        a[i][l]=f-g;
        f=0.0;
        for (j=1;j<=l;j++) 
        {
          a[j][i]=a[i][j]/h;
          g=0.0;
          for (k=1;k<=j;k++)
            g += a[j][k]*a[i][k];
          for (k=j+1;k<=l;k++)
            g += a[k][j]*a[i][k];
          off[j]=g/h;
          f += off[j]*a[i][j];
        }
        hh=f/(h+h);
        for (j=1;j<=l;j++) 
        { 
          f=a[i][j];
          off[j]=g=off[j]-hh*f;
          for (k=1;k<=j;k++)
            a[j][k] -= (f*off[k]+g*a[i][k]);
        }
      }
    } 
    else
      off[i]=a[i][l];
    diag[i]=h;
  }
  diag[1]=0.0;
  off[1]=0.0;  for (i=1;i<=n;i++) 
  { 
    l=i-1;
    if (diag[i]) 
    { 
      for (j=1;j<=l;j++) 
      {
        g=0.0;
        for (k=1;k<=l;k++) 
          g += a[i][k]*a[k][j];
        for (k=1;k<=l;k++)
          a[k][j] -= g*a[k][i];
      }
    }
    diag[i]=a[i][i]; 
    a[i][i]=1.0;    for (j=1;j<=l;j++) 
      a[j][i]=a[i][j]=0.0;
  }
}


**************************************************/



void Cholesky(const SymmetricMatrix<double> &d, 
                      LowerTriangleMatrix<double> &a)
{
TRY
  // Consistency check
  ASSERT1(a.Size() == d.Size(), Argument("Unequal sizes of in- and output"));
    
  // Compute elements 
  for (size_t i=0; i < d.NRow(); ++i)
  {
    // do diagonal element
    a(i, i) = d(i, i);
    for (size_t j=0; j < i; ++j)   
    {
      // do non-diagonal elements
      a(i, j) = d(i, j);
      for (size_t k=0; k < j; ++k)
        a(i, j) -= a(i, k) * a(j, k);
//	    if (a(i,j) > 1.0e-12)
		    a(i, j) /= a(j, j);
      // do diagonal elements
      a(i, i) -= Square(a(i, j));
    }
//    ASSERT0(a(i, i) >= 0.0, Numerical("Input matrix not positive semi-definite."));
    a(i, i) = sqrt(a(i, i));
  }
CATCH(Cholesky)
}

void Cholesky(const SymmetricMatrix<double> &d, 
                      LowerTriangleMatrix<float> &a)
{
  LowerTriangleMatrix<double> ad(a.Size());
  Cholesky(d, ad);
  for (int k = 0; -k < static_cast<int>(a.Size()); --k)
    for (unsigned int l = 0; l < a.Size() + k; ++l)
      a.D(k,l) = static_cast<float>(ad.D(k,l));
}

// Invert lower triangular matrix
LowerTriangleMatrix<double>	LTInvert(const LowerTriangleMatrix<double>& lt)
{
	LowerTriangleMatrix<double> inv(lt.Size());
	
	for (size_t i=0; i < lt.Size(); ++i)
	{
//    if (lt(i,i) == 0.0)
  //    throw Logical("Cannot invert a LT matrix with zero on diagonal");

		// do diagonal elements
		inv(i,i) = 1.0 / lt(i,i);
		// do lower triangle
		for (size_t j=0; j < i; ++j)
		{
			inv(i,j) = 0.0;
			for (size_t k=j; k < i; ++k)
			{
				inv(i,j) += lt(i,k) * inv(k,j);
			}
      inv(i,j) /= -lt(i,i);
		}
	}
	return inv;
}

// Approximate a symmetric matrix by its first 'n' principal components
void PCConstruct(SymmetricMatrix<double>& s, size_t n)
{
  // get work space
   SquareMatrix<double> a(s);
	valarray<double> ev(a.NRow());
	
	// do orthogonal decomposition
	OrthogonalDecomposition_Sort(s, 30, a, ev);

	// build matrix from eigenvectors and -values
  for (size_t i=0; i < a.NRow(); ++i)
    for (size_t j=0; j < a.NRow(); ++j)         
    {
      s(i,j) = 0.0;
      for (size_t k=0; k < n; ++k)
        s(i,j) += ev[k] * a(i,k) * a(j,k) ;
    }
}


void Factors(const SymmetricMatrix<double> &c, 
                     Matrix<double> &a,
                     const double tol,
                     const unsigned int maxIter)
{
TRY
  // Consistency check
  ASSERT1(a.NRow() == c.Size(), Argument("Unequal sizes of in- and output"));

  // # of factors
  unsigned int numFactors = a.NColumn();

  // Work space
  SymmetricMatrix<double> cs(c);
  SquareMatrix<double> eigenVectors(c.Size());
  valarray<double> eigenValues(c.Size()), tmp(c.Size()), sqrtEVs(1.0, numFactors), evs(0.0, numFactors);

  // Initialize 
  a *= 0.0;

  // Recursion
  for (unsigned int i = 0; i < maxIter; ++i)
  {
    // Do orthogonal decomposition
    OrthogonalDecomposition_Sort(cs, 30, eigenVectors, eigenValues, tmp);

    // Compute change compared to previous 
    // and save relevant eigenvectors
    double d = 0.0;
    for (unsigned int j = 0; j < numFactors; ++j)
    {  
      d += Square(evs[j] - eigenValues[j]);
      evs[j] = eigenValues[j];
    }
    if (sqrt(d) < tol)
      break;

    // Replace diagonal by square of vector entries
    for (unsigned int k = 0; k < c.Size(); ++k)
    {
      cs(k,k) = 0.0;
      for (unsigned n = 0; n < numFactors; ++n)
        cs(k,k) += eigenValues[n] * Square(eigenVectors(k,n));
    }
  }
//  if (maxIter == i)
//    throw Numerical("Maximum number of iterations reached.");
//  else
    // Copy factors to output
    for (unsigned k = 0; k < a.NRow(); ++k)
      for (unsigned n = 0; n < numFactors; ++n)
        a(k,n) = sqrt(eigenValues[n]) * eigenVectors(k,n);

CATCH(Factors)
}

///////////////////////////////////////////////////////// GRADIENT.H ////////////////////////////////

/*!
  \file gradient.h

  \brief Calculation of numerical gradients
*/


/*! 
  \brief Numerical calculation of gradients.

  \param func     Function (object) for which a gradient is wanted.
  \param x        Calculation gradient at domain point x.
  \param difftype The numerical differention schema: Valid difftypes 
                  are 1, 2, and 3. Numerical differentation uses \f$n\f$,                  
                  \f$2n\f$, and \f$6n\f$ function evaluations for each 
                  gradient if difftype is 1, 2, 3 respectively.
                  - 1: Ordinary forward difference quotient with stepsize 
                      \f[ 0.1 epsfcn^{\frac{1}{2}}.\f]
                  - 2: Symmetric difference quotient with stepsize 
                      \f[ 0.1 epsfcn^{\frac{1}{3}}.\f]
                  - 3: Sixth order approximation based on Richardson 
                      extrapolation with 
                      \f[ 0.1 epsfcn^{\frac{1}{7}}.\f]
  \param epsfcn   The expected relative precision of the object function 
                  evaluation.
  \param taubnd   Amount by which bounds may be violated if numerical 
                  differentation is used.

  \return         The numerical gradient at point x.
*/
template<class Func>
valarray<double> Gradient(Func func,
                            valarray<double> &x,
                            int difftype = 1,
                            double epsfcn = 1.e-16,
                            double taubnd = 1.0){
TRY

  valarray<double> grad(x.size()); 
  // ordinary forward difference
  if (difftype == 1) 
  {
    double deldif = Min(0.1 * sqrt(epsfcn), 0.01);

    double gixloc = func(x); 
  
    for (size_t j = 0; j < x.size(); ++j)    {
      double xhelp = x[j];
      double xincr = Min(Min(0.01, deldif * Abs(xhelp) + deldif), taubnd);
      if (xhelp >= 0.0)        x[j] = xhelp + xincr;
      else 
        x[j] = xhelp - xincr;
      double fhelp = func(x); 
    
      grad[j] = (fhelp - gixloc) / (x[j] - xhelp);
      x[j] = xhelp;
    }
  } 
  // symmetric difference
  else if (difftype == 2) 
  {
    double deldif = Min(0.1 * pow(epsfcn, 1.0 / 3.0), 0.01);

    for (size_t j = 0; j < x.size(); ++j)    {
      double xhelp = x[j];
      double xincr = Min(Min(0.01, deldif * Abs(xhelp) + deldif), taubnd);
      x[j] = xhelp + xincr;
    
      double fhelp1 = func(x); 
    
      x[j] = xhelp - xincr;
    
      double fhelp2 = func(x); 
    
      grad[j] = (fhelp1 - fhelp2) / (xincr + xincr);
      x[j] = xhelp;
    }
  } 
  // sixth order approximation
  else if (difftype == 3) 
  {
    double deldif = Min(0.1 * pow(epsfcn, 1.0 / 7.0), 0.01);

    for (size_t j = 0; j < x.size(); ++j)    {
      double xhelp = x[j];
      double xincr = Min(Min(0.01, deldif * Abs(xhelp) + deldif), taubnd / 4.0);
      x[j] = xhelp - xincr;
    
      double fhelp1 = func(x);
    
      x[j] = xhelp + xincr;
    
      double fhelp2 = func(x);
    
      xincr  = xincr + xincr;
      double d1 = xincr;
      x[j] = xhelp - xincr;
    
      double fhelp3 = func(x);
    
      x[j] = xhelp + xincr;
    
      double fhelp4 = func(x);
    
      xincr = xincr + xincr;
      double d2 = xincr;
      x[j] = xhelp - xincr;
    
      double fhelp5 = func(x);
    
      x[j] = xhelp + xincr;
    
      double fhelp6 = func(x);
    
      x[j] = xhelp;
      double d3 = xincr + xincr;
      double sd1 = (fhelp2 - fhelp1) / d1;
      double sd2 = (fhelp4 - fhelp3) / d2;
      double sd3 = (fhelp6 - fhelp5) / d3;
      sd3 = sd2 - sd3;
      sd2 = sd1 - sd2;
      sd3 = sd2 - sd3;
      grad[j] = sd1 + 0.4 * sd2 + sd3 / 45.0;
    }
  }
  return grad;
CATCH(math::basnum::Gradient)
}

 



/////////////////////////////////////////////////////// QLIB STUFF ///////////////////////////////////

// namespace{

/*! 
  \brief An OptimizationProblem is a general optimization problem with
  equality and inequality constraints and domain bounds.

  Equality constraints are represented by functions which must evaluate
  to zero.  Similarly, inequality constraints are represented by functions which must
  evaluate non-negative.
  Domain bounds are implemented as two maps (one for upper and one for lower 
  bounds) from coordinate index (a int) to parameter bound (a double).  
  Iff a coordinate has no bound its index is absent from the map.
*/
class OptimizationProblem
{
public:
// foundation
  virtual ~OptimizationProblem() {}

// inspection
  //! Return (real manifold) dimension of solution space.
  virtual unsigned int SolutionDimension() const = 0;
  //! Return value of object function at domain point x.
  virtual double ObjectFunction(const valarray<double> &x) const = 0;
  //! Return number of equality constraints.
  virtual unsigned int NZeroConstraint() const;
  //! Return value of i'th zero constraint at domain point x.
  virtual double ZeroConstraint(int i, 
                                const valarray<double> &x) const;
  //! Return number of inequality constraints.
  virtual unsigned int NPositivityConstraint() const;
  //! Return value of i'th positivity constraint at domain point x.
  virtual double PositivityConstraint(int i, 
                                      const valarray<double> &x) const;
  //! Return upper bounds.
  virtual const DomainBound& UpperBound() const = 0;
  //! Return lower bounds.
  virtual const DomainBound& LowerBound() const = 0;
}; // OptimizationProblem



//////////////// INLINES ////////////////////////////

inline unsigned int OptimizationProblem::NZeroConstraint() const
{
  return 0;
}

inline double OptimizationProblem::ZeroConstraint(int i, 
                                const valarray<double> &x) const
{
  return 0.0;
}

inline unsigned int OptimizationProblem::NPositivityConstraint() const
{
  return 0;
}

inline double OptimizationProblem::PositivityConstraint(int i, 
                                      const valarray<double> &x) const
{
  return 0.0;
}

// Optimization problem class with QLib ctor (no support for constraints for now)
class QOptimizationProblem : public OptimizationProblem
{
public:
// foundation 
    //! Construct from MFunctionND and domain bounds
    QOptimizationProblem(const MFunctionND&  func);

// inspection
    //! Return (real manifold) dimension of solution space.
    unsigned int SolutionDimension() const { return static_cast<unsigned int>(m_func.getNbVars()); }
    //! Return value of object function at domain point x.
    double ObjectFunction(const valarray<double> &x) const;
    //! Return upper bounds.
    const DomainBound& UpperBound() const { return m_upperBound; }
    //! Return lower bounds.
    const DomainBound& LowerBound() const { return m_lowerBound; }

private:
// data
    const MFunctionND&  m_func;
    DomainBound m_upperBound, m_lowerBound;
// work space
    mutable CDoubleArray m_f, m_x;
}; // QOptimizationProblem


//} // unnamed namespace

/////////////////////////////////////////////// SPELLUCCI STUFF ////////////////////////////////////////

//DRLIB_BEGIN_NAMESPACE

/*! 
  \brief This class implements the Spellucci %optimization 
         algorithm called donlp2.

  The Spelluci %optimization algorithm is used to minimise an 
  object function, \f$ f(x) \f$, over an \f$ n \f$ dimensional 
  parameter space, given \f$ N_G \f$ inequality constraints 
  \f$ g_i(x) \f$, and \f$ N_H \f$ equality constraints 
  \f$ h_i(x) \f$. 

*/
class QSpellucci::Spellucci
{
public:
// foundation
/*! 
    Default construction

    \param del0     Parameter used to determine binding inequality constaints 
                    in the initial phase of the %optimization. A constraint is 
                    considered binding if
                    \f[ g_i(x) / \max \{ 1, \|\nabla g_i(x) \| \} 
                    \leq del0 \f]
    \param tau0     Bound determine how much the unscaled penalty-term may
                    derivate from zero. The algorithm assumes that within the
                    region described by
                    \f[ \sum_{i=1}^{N_H} |h_i(x)| - 
                    \sum_{i=1}^{N_G} \min\{ 0, g_i(x) \} \leq tau0 \f]
                    all functions may be evaluated safely.
    \param taubnd   Amount by which bounds may be violated if numerical 
                    differentation is used.
    \param epsdif   The expected relative precision of the gradient evaluation.
    \param epsfcn   The expected relative precision of the object function 
                    evaluation.
    \param difftype The numerical differention schema: Valid difftypes 
                    are 1, 2, and 3. Numerical differentation uses \f$n\f$, 
                    \f$2n\f$, and \f$6n\f$ function evaluations for each 
                    gradient if difftype is 1, 2, 3 respectively.
                    - 1: Ordinary forward difference quotient with step size 
                        \f[ 0.1 epsfcn^{\frac{1}{2}}.\f]
                    - 2: Symmetric difference quotient with step size 
                        \f[ 0.1 epsfcn^{\frac{1}{3}}.\f]
                    - 3: Sixth order approximation based on Richardson 
                        extrapolation with step size
                        \f[ 0.1 epsfcn^{\frac{1}{7}}.\f]
  \param iterma     The maximal number of main loop iterations.
  \param epsx       One of the termination criteria. 
                    \f[ \|\nabla L(x,\mu,\lambda) \| 
                    \leq epsx(1 + \|\nabla f(x) \|), \f]
                    with \f$ L \f$ being the Lagrange function corresponding to 
                    the constrained %optimization problem.
                    This is not the only termination criterion but the others
                    are purely internal.
  */
  Spellucci(double del0     = 0.2e0,  
            double tau0     = 1.e0,
            double taubnd   = 1.0,
            double epsdif   = 1.e-16,
            double epsfcn   = 1.e-16,
            int    difftype = 1,
            int    iterma   = 1000,
            double epsx     = 1.e-8) 
  : m_x(0.0, 0), 
    m_del0(del0),
    m_tau0(tau0),
    m_taubnd(taubnd),
    m_epsdif(epsdif),
    m_epsfcn(epsfcn),  
    m_difftype(difftype),
    m_iterma(iterma),
    m_epsx(epsx),
    gunit(1, 1),
    accinf(iterma + 1, 33),
    gres(1, 1),
    qr(1, 1),
    a(1, 1),
    fugrad(1, 1),
    xj(1, 1),
    r(1, 1),
    NSTEP(40)
  {
    mf_InitialiseOptimizer();
  }



// Inspection
  
  /*! 
  \brief Perform the minimization of a function \f$ n\f$ variables. 

  \param problem The %optimization problem

  \param y  The starting guess

  \throw Argument y.size() != problem.SolutionDimension()

  */
  void Optimize(const OptimizationProblem &problem, 
                const valarray<double> &y) const;
  //! Return optimal point.
  valarray<double> OptimalPoint() const;
  //! Return # evaluations of object function.
  long NObjectEvaluation() const;
  //! Return # evaluations of equality constraints.
  long NZeroConstraintEvaluation() const;
  //! Return # evaluations of inequality constraints.
  long NPositivityConstraintEvaluation() const;
private:
// types
  //! Wrapper for the constrained %optimization problem.
  class WrapperOptimizationProblem
  {
  public:
  // foundation
    WrapperOptimizationProblem(const OptimizationProblem &problem,
                               const Spellucci &spellucci)
    : m_op(problem),
      m_spellucci(spellucci)
    {
      mf_Initialise();
    }
  // inspection
    //! Straight evaluation on point.
    double ObjectFunction(const valarray<double> &x) const;
    /*! 
      \brief Evaluation of zero constraint on point.

      \param Zero constraint number.
      \param Point where constraint is evaluated.

      \return The constraint evaluated at x.
    */
    double ZeroConstraint(int i, const valarray<double> &x) const;
    /*! 
      \brief Evaluation of positivity constraint on point.

      \param i  Positivity constraint number.
      \param x  Point where constraint is evaluated

      \throw Argument i is less than zero or greater than the highest 
                      constraint number. 
      \return The constraint evaluated at x.

    */
    double PositivityConstraint(int i, const valarray<double> &x) const;
    //! Return (real manifold) dimension of solution space
    int SolutionDimension() const;
    //! Return number of equality constraints
    int NZeroConstraint() const;
    //! Return number of inequality constraints
    int NPositivityConstraint() const;
    //! Return upper bounds.
    const DomainBound& UpperBound() const;
    //! Return lower bounds.
    const DomainBound& LowerBound() const;
    //! Evaluation of lower bound constraint on point.
    double LowerBoundConstraint(int i, const valarray<double> &x) const;
    //! Evaluation of upper bound constraint on point.
    double UpperBoundConstraint(int i, const valarray<double> &x) const;
    
  private:
  // foundation
    void mf_Initialise() const;
  // data
    //! The constrained %optimization problem.
    const OptimizationProblem &m_op;
    //! The optimizer.
    const Spellucci &m_spellucci;
    //! Workspace for bound constraints.
    mutable std::vector<int> m_lower_index;
    mutable std::vector<double> m_lower_bound;
    mutable std::vector<int> m_upper_index;
    mutable std::vector<double> m_upper_bound;
  };

// friends
  friend class WrapperOptimizationProblem;

// Spellucci specific subroutines and variables 
  //! Initialise internal variables depending on input data.
  void mf_InitialiseOptimizer() const;
  //! Initialise the quasi Newton update with a multiple of the identity
  void mf_InitialiseQuasiNewtonUpdate(void) const;
  //! The driver routine. 
  void mf_Solve(const WrapperOptimizationProblem &problem) const;
  //! Initialise internal variables which does not depend on input data.
  void mf_InitialiseSolver(const WrapperOptimizationProblem &problem) const;
  //! Computation of new scaling factors for L1-penalty-function.
  void mf_ScalingFactors(void) const;
  /*! 
    \brief Computation of the Pantoja-Mayne BFGS-update of the Hessian matrix.

    The updating is done in the Cholesky decomposition of the Hessian
    and not in Hessian directly. 
    
    In cases with unconstrained %optimization problem a %Powell updating 
    schema is used instead.

    The Hessian is stored in the matrix a.
  */
  void mf_updateBFGS(void) const;
  /*! 
    Equality constrained recursive quadratic programming with multiple inactivation and
    superlinearly convergent projected BFGS-update.
  */
  void mf_RecursiveQuadraticSolver(const WrapperOptimizationProblem &problem) const;
  //! Compute the directional derivative of the L1-penalty-function
  void mf_DerivativePenaltyFunction(void) const;
  //! Cut d if appropriate and rescale.
  void mf_Cut_d(void) const;
  /*!
    Compute maximum stepsize \f$stmaxl\f$ such that projection on the box of lower and 
    upper bounds changes for sig in \f$\left[0, stmaxl\right]\f$, if such exists.
  */
  void mf_MaxStepSize(void) const;
  //! Restore the best point found so far to be the current new point
  void mf_RestoreCurrentMinimum(void) const;
  //! Save the best point found so far.
  void mf_SaveCurrentMinimum(void) const;
  //! Evaluate the object and constraint functions at some new point.
  double mf_Evaluate(double sigact,
                       const WrapperOptimizationProblem &problem) const;
  //! Determination of stepsize by an Armijo-like test for descent.
  void mf_StepSize(double sig1th, 
                   const WrapperOptimizationProblem &problem) const;
  //! Compute gradient of lagrangian.
  void mf_GradientLagrangian(valarray<double> &gphi) const;
  /*!
    QR-decomposition of matrix of gradients of binding constraints this set may be expanded using
    multiple calls to this method. No exit on singular r-factor here. Information on the 
    decompostion is stored in betaq and in and below the diagonal of qr. r-factor is stored in diag
    (diagonal) and above the diagonal of qr. cscal is the column scaling of the original matrix.
    Column pivoting is done here and is stored in colno.
  */
  void mf_QRDecompBindingConstraints(int nlow,
                                     int nrl) const;
  /*!  
    Application of Householder transformations stored in the lower or strict lower 
    (if incr = 0 or 1 resp.) triangle of a and in the vector beta on b giving c. Only columns is1
    to is2 are used in forward manner if id > 0, backwards otherwise. Rows is1 to m of c are 
    changed only
  */
  void mf_HouseholderTransformation(int id,
                                    int incr,
                                    int is1,
                                    int is2,
                                    int m,
                                    Matrix<double> &a,
                                    valarray<double> &beta,
                                    valarray<double> &b,
                                    valarray<double> &c) const;
  /*! 
    Solve triangular system \f$rx = b\f$, where \f$r\f$ is defined by the Householder-QR
    decomposition decomp (with column scaling).
  */
  void mf_SolveTriangularSystem(int nlow,
                                int nup,
                                valarray<double> &b,
                                valarray<double> &x) const;
  /*! 
    Solve triangular system \f$r^t = b\f$, where \f$r\f$ is defined by the Householder-QR
    decomposition decomp (with column scaling).
  */
  void mf_SolveTransposedTriangularSystem(int nlow,
                                          int nup,
                                          valarray<double> &b,
                                          valarray<double> &x) const;
// TODO: move outside class!
  //! Euclidian norm of vector \f$(a, b)\f$. 
  double norm(double a, 
                double b) const;
  /*! 
    Computes the upper triangular Cholesky factor \f$ r1 \f$ of 
    \f$ r^t r + z z^t -y y^t \f$ and restores it in \f$ r \f$. 
    The strict lower triangle of \f$ r \f$ remains unchanged.

    \throw Argument When input valarrays have unequal sizes.
    \throw Numerical When the decomposition does not exists.
  */
  void mf_UpperTriangularCholeskyFactor(SquareMatrix<double> &r,
                                        valarray<double> &z,
                                        valarray<double> &y) const;
  /*!
    Solves \f$ry = b\f$ where \f$r\f$ is the Cholesky factor of a stored in the upper half of a.
  
    \throw  Argument when dimensions of the matrix and vectors aren't compatible.
    \return The Euclidian norm of y.
  */
  double mf_SolveCholeskyFactorRight(SquareMatrix<double> &a,
                  valarray<double> &b,
                  valarray<double> &y,
                  int n) const;
  /*! 
    Solves \f$r^t y = b\f$ where \f$r\f$ is the Cholesky factor of a stored in the upper half of a.

    \throw  Argument when dimensions of the matrix and vectors aren't compatible.
    \return The Euclidian norm of y.
  */
  double mf_SolveCholeskyFactorLeft(SquareMatrix<double> &a,
                                      valarray<double> &b,
                                      valarray<double> &y,
                                      int n) const;
  //! Extended recursive quadratic solver.
  void mf_ExtendedRecursiveQuadraticSolver(void) const;
  //! Compute updated projected gradient (primal).
  void mf_UpdateProjGradient(valarray<double> &z) const;
  //! Compute correction of dual multipliers.
  void mf_CorrDualMultipliers(valarray<double> &rv) const;
  //! Delete constraint number l.
  void mf_DeleteConstraint(valarray<int> &ai, int l) const;
  //! Add constraint whose gradient is given by np.
  void mf_AddConstraint(void) const;
  /*! 
    Computes the inverse of the upper triangular matrix part of a and stores it in
    the upper triangle of the right lower part of x                                  

    \param n     a is a \f$n \times n\f$ matrix.
    \param ndual x is a \f$ndual \times ndual\f$ matrix.
  */
  void mf_InverseUpperTriangularMatrix(int n,
                                       SquareMatrix<double> &a,
                                       int ndual,
                                       SquareMatrix<double> &x) const;
// TODO: Remove after index offset refactorization.
  /*! 
    Object function. This member function performs a shift on x to overcome the problem with
    one index offset. 
  */
  double mf_ObjectFunction(valarray<double> &x,
               const WrapperOptimizationProblem &problem) const;
// TODO: Remove after index offset refactorization.
  /*! 
    Zero constraints. This member function performs a shift on x to overcome the problem with
    one index offset. 
  */
  double mf_ZeroConstraints(int i,
                              valarray<double> &x,
                              const WrapperOptimizationProblem &problem) const;
// TODO: Remove after index offset refactorization.
  /*! 
    Positivity constraints. This member function performs a shift on x to overcome the problem with
    one index offset. 
  */
  double mf_PositivityConstraints(int i,
                                    valarray<double> &x,
                                    const WrapperOptimizationProblem &problem) const;

  //! Initialise workspace, i.e. resizing of valarrays and matrices.
  void mf_InitialiseWorkSpace(const WrapperOptimizationProblem &problem) const;

// TODO: Remove after index offset refactorization.
  /*! 
    \brief Numerical calculation of gradients.

    \param func     Function (object) for which a gradient is wanted.
    \param x        Calculation gradient at domain point x.
    \param scale    Rescale x with scale before gradient calculation.
    \param difftype The numerical differention schema: Valid difftypes 
                    are 1, 2, and 3. Numerical differentation uses \f$n\f$, 
                    \f$2n\f$, and \f$6n\f$ function evaluations for each 
                    gradient if difftype is 1, 2, 3 respectively.
                    - 1: Ordinary forward difference quotient with stepsize 
                        \f[ 0.1 epsfcn^{\frac{1}{2}}.\f]
                    - 2: Symmetric difference quotient with stepsize 
                        \f[ 0.1 epsfcn^{\frac{1}{3}}.\f]
                    - 3: Sixth order approximation based on Richardson 
                        extrapolation with 
                        \f[ 0.1 epsfcn^{\frac{1}{7}}.\f]
    \param epsfcn   The expected relative precision of the object function 
                    evaluation.
    \param taubnd   Amount by which bounds may be violated if numerical 
                    differentation is used.

    \return         The numerical gradient at point x.

    \throw          Argument x.size() != scale.size()
  */
  template<class Func>
  valarray<double> NumericalGradient(Func func,
                                       valarray<double> &x,
                                       valarray<double> &scale,
                                       int difftype = 1,
                                       double epsfcn = 1.e-16,
                                       double taubnd = 1.0) const
  {
    // index shift as Gradient expects 0 index offset.
    valarray<double> xtr(x.size() - 1);
    std::copy(&x[1], &x[x.size()], &xtr[0]);

    valarray<double> scaletr(scale.size() - 1);
    std::copy(&scale[1], &scale[scale.size()], &scaletr[0]);

    // calculate gradient
    valarray<double> scaledpoint = xtr * scaletr;
    valarray<double> gradtr = Gradient(func, scaledpoint, difftype, epsfcn, taubnd);
    // rescale 
    gradtr *= scaletr;
    
    // index shift back to 1 index offset.
    valarray<double> grad(x.size());
    std::copy(&gradtr[0], &gradtr[gradtr.size()], &grad[1]);

    return grad;
  }

// data
  //! Current point
  mutable valarray<double> m_x;
  //! The number of object function evaluations.
  mutable long m_nobjeval;
  //! The number of zero constraint functions evaluations.
  mutable long m_nzeroconeval;
  //! The number of positivity constraint functions evaluations.
  mutable long m_nposconeval;

  double m_del0;
  double m_tau0;
  double m_taubnd;
  double m_epsdif;
  double m_epsfcn;
  int m_difftype;
  int m_iterma;
  double m_epsx; 

  // Internal variables
  mutable valarray<int> val;
  mutable valarray<int> gconst;
  mutable valarray<int> llow;
  mutable valarray<int> lup;
  mutable valarray<int> cres;
  mutable valarray<int> cgres;
  mutable valarray<int> aitr;
  mutable valarray<int> colno;
  mutable valarray<int> perm1;
  mutable valarray<int> perm;
  mutable valarray<int> violis;
  mutable valarray<int> alist;
  mutable valarray<int> bind;
  mutable valarray<int> bind0;
  mutable valarray<int> sort;

  mutable valarray<double> difx;
  mutable valarray<double> dd;
  mutable valarray<double> d0;
  mutable valarray<double> d;
  mutable valarray<double> xmin;
  mutable valarray<double> x1;
  mutable valarray<double> x0;
  mutable valarray<double> x;
  mutable valarray<double> gphi1;
  mutable valarray<double> gphi0;
  mutable valarray<double> qgf;
  mutable valarray<double> gradf;
  mutable valarray<double> diag0;
  mutable valarray<double> ug;
  mutable valarray<double> og;
  mutable valarray<double> xst;
  mutable valarray<double> xtr;
  mutable valarray<double> xsc;

  mutable valarray<double> resmin;
  mutable valarray<double> gresn;
  mutable valarray<double> diag;
  mutable valarray<double> cscal;
  mutable valarray<double> betaq;
  mutable valarray<double> colle;
  mutable valarray<double> u;
  mutable valarray<double> u0;
  mutable valarray<double> w;
  mutable valarray<double> w1;
  mutable valarray<double> res;
  mutable valarray<double> res0;
  mutable valarray<double> res1;
  mutable valarray<double> resst;
  mutable valarray<double> yu;
  mutable valarray<double> slack;
  mutable valarray<double> work;
  mutable valarray<double> delfac;
  mutable valarray<double> fu;

  mutable valarray<double> ddual;
  mutable valarray<double> np;

  mutable valarray<double> ud;
  mutable valarray<double> ud1;

  mutable Matrix<int> gunit;
  mutable Matrix<double> accinf;
  mutable Matrix<double> gres;
  mutable Matrix<double> qr;
  mutable SquareMatrix<double> a;
  mutable Matrix<double> fugrad;
  mutable SquareMatrix<double> xj;
  mutable Matrix<double> r;

  mutable int NX; 
  mutable int NRESM;
  mutable int NSTEP;
  mutable int NDUALM;
  mutable int MDUALM;

  mutable int rank;
  mutable int nreset;
  mutable int qpterm;
  mutable int ndual;
  mutable int iq;
  mutable int iptr;
  mutable int iqtr;
  mutable int n;
  mutable int nr;
  mutable int nres;
  mutable int nh;
  mutable int ng;
  mutable int inx;
  mutable int singul;
  mutable int ident;
  mutable int eqres;
  mutable int analyt;
  mutable int cold;
  mutable int icf;
  mutable int icgf;
  mutable int cfincr;
  mutable int itstep;
  mutable int phase; 
  mutable int iterma;
  mutable int ifill1; 
  mutable int lastdw;
  mutable int lastup;
  mutable int lastch; 
  mutable int bloc;
  mutable int valid;
  mutable int corr;
  mutable int difftype;

  mutable double gfn;
  mutable double scalm;
  mutable double scalm2;
  mutable double matsc;
  mutable double scf;
  mutable double scf0;
  mutable double infeas;
  mutable double rnorm;
  mutable double rlow;
  mutable double epsfcn;
  mutable double taubnd;
  mutable double epsmac;
  mutable double tolmac;
  mutable double deldif;
  mutable double epsdif;
  mutable double runtim;
  mutable double optite;
  mutable double upsi;
  mutable double upsi0;
  mutable double upsi1;
  mutable double upsist;
  mutable double psi;
  mutable double psi0;
  mutable double psi1;
  mutable double psist;
  mutable double psimin;
  mutable double phi;
  mutable double phi0;
  mutable double phi1;
  mutable double phimin;
  mutable double fx;
  mutable double fx0;
  mutable double fx1;
  mutable double fxst;
  mutable double fmin;
  mutable double b2n;
  mutable double b2n0;
  mutable double xnorm;
  mutable double x0norm;
  mutable double sig0;
  mutable double dscal;
  mutable double dnorm;
  mutable double d0norm;
  mutable double sig;
  mutable double sigmin;
  mutable double dirder;
  mutable double cosphi;
  mutable double upsim;
  mutable double del;
  mutable double del0;
  mutable double del01;
  mutable double delmin;
  mutable double tau0;
  mutable double tau;
  mutable double ny;
  mutable double smalld;
  mutable double smallw;
  mutable double rho;
  mutable double rho1;
  mutable double eta;
  mutable double epsx;
  mutable double c1d;
  mutable double scfmax;
  mutable double updmy0;
  mutable double tauqp;
  mutable double taufac;
  mutable double taumax;
  mutable double alpha;
  mutable double beta;
  mutable double theta;
  mutable double sigsm;
  mutable double sigla;
  mutable double delta;
  mutable double stptrm;
  mutable double delta1;
  mutable double stmaxl;
  mutable double level;
  mutable double clow;
  mutable double sstr;
  mutable double riitr;
};  // Spellucci

inline valarray<double> QSpellucci::Spellucci::OptimalPoint() const
{
  return m_x;
}

inline long QSpellucci::Spellucci::NObjectEvaluation() const
{ 
  return m_nobjeval;
}

inline long QSpellucci::Spellucci::NZeroConstraintEvaluation() const
{ 
  return m_nzeroconeval;
}

inline long QSpellucci::Spellucci::NPositivityConstraintEvaluation() const
{ 
  return m_nposconeval;
}

inline double 
QSpellucci::Spellucci::WrapperOptimizationProblem::ObjectFunction(const valarray<double> &x) const
{
  ++m_spellucci.m_nobjeval;
  return m_op.ObjectFunction(x);
}

inline double 
QSpellucci::Spellucci::WrapperOptimizationProblem::ZeroConstraint(int i, const valarray<double> &x) const
{
  ++m_spellucci.m_nzeroconeval;
  return m_op.ZeroConstraint(i, x);
}

inline double 
QSpellucci::Spellucci::WrapperOptimizationProblem::PositivityConstraint(int i, 
                                                            const valarray<double> &x) const
{
  // total number of positivity constraints, box constraints included
  int ntpc = NPositivityConstraint();
  // number of lower bound constraints
  int nlb = LowerBound().size();
  // number of upper bound constraints
  int nub = UpperBound().size();
  // number of ordinary positivity constraints.
  int nopc = ntpc - nlb - nub;

  ASSERT0(i > 0, Argument("Index of positivity less than 0!"));
  ASSERT0(i <= ntpc, 
    Argument("Index of positivity larger than total number of positivity constraint!"));

  ++m_spellucci.m_nposconeval;
  if (i <= nopc)
    return m_op.PositivityConstraint(i, x);
  if (i <= ntpc - nub)
    return LowerBoundConstraint(i - nopc, x);
  // if (i <= ntpc)
  return UpperBoundConstraint(i - nopc - nlb, x);
}

inline int QSpellucci::Spellucci::WrapperOptimizationProblem::SolutionDimension() const 
{ 
  return m_op.SolutionDimension();
} 

inline int QSpellucci::Spellucci::WrapperOptimizationProblem::NZeroConstraint() const
{ 
  return m_op.NZeroConstraint();
} 

inline int QSpellucci::Spellucci::WrapperOptimizationProblem::NPositivityConstraint() const
{ 
  return m_op.NPositivityConstraint() + LowerBound().size() + UpperBound().size();
} 

inline const DomainBound& QSpellucci::Spellucci::WrapperOptimizationProblem::UpperBound() const
{ 
  return m_op.UpperBound();
} 

inline const DomainBound& QSpellucci::Spellucci::WrapperOptimizationProblem::LowerBound() const
{ 
  return m_op.LowerBound();
} 

inline double 
QSpellucci::Spellucci::WrapperOptimizationProblem::LowerBoundConstraint(int i,
                                                            const valarray<double> &x) const
{
  return x[m_lower_index[i - 1]] - m_lower_bound[i - 1];
}

inline double 
QSpellucci::Spellucci::WrapperOptimizationProblem::UpperBoundConstraint(int i,
                                                            const valarray<double> &x) const
{
  return m_upper_bound[i - 1] - x[m_upper_index[i - 1]];
}

inline void QSpellucci::Spellucci::WrapperOptimizationProblem::mf_Initialise() const
{
  DomainBound d = LowerBound();
  DomainBoundIterator iter = d.begin();
  DomainBoundIterator end = d.end();
  for (; iter != end; ++iter)
  {
    m_lower_index.push_back(iter->first);
    m_lower_bound.push_back(iter->second);
  }
  d = UpperBound();
  iter = d.begin();
  end = d.end();
  for (; iter != end; ++iter)
  {
    m_upper_index.push_back(iter->first);
    m_upper_bound.push_back(iter->second);
  }
}

QOptimizationProblem::QOptimizationProblem(const MFunctionND&  func)
: m_func(func),
  m_f(1), m_x(func.getNbVars())
{ 
  if (func.getNbFuncs() != 1)
    throw ModelException("Spellucci", "Range dimension > 1");
  if (!func.getIntervals().empty())
    for (size_t i = 0; i < size_t(func.getIntervals().size()); ++i)
    {
      const Range& range = *func.getIntervals()[i];
      if (!range.getUpper().isInfinite())
        m_upperBound[i] = range.getUpper().getValue();
      if (!range.getLower().isInfinite())
        m_lowerBound[i] = range.getLower().getValue();
    }
} 

double QOptimizationProblem::ObjectFunction(const valarray<double> &x) const
{ 
  // Copy argument from valarray to CDoubleArray
  for (size_t i = 0; i < SolutionDimension(); ++i)
    m_x[i] = x[i];
  // Evaluate and return
  m_func(m_x, m_f);
  return m_f[0];
}

QSpellucci::~QSpellucci() 
{ 
  delete mp_spellucci; 
}

QSpellucci::QSpellucci(double del0,  
			 	               double tau0,
							         double taubnd,
                       double epsdif,
                       double epsfcn,
                       int    difftype,
                       int    iterma,
                       double epsx) 
: OptimizerND(TYPE),
  mp_spellucci(new Spellucci(del0, tau0, taubnd, epsdif, epsfcn, difftype, iterma, epsx))
{ 
  // empty
}            

void QSpellucci::minimize(const MFunctionND&  func,
                         const CDoubleArray& inxguess,      // initial guess
                         CDoubleArray&       outx) const    // result
{  
    CStringArray ids(0);
    minimize(func,
             inxguess,
             ids,
             outx);
}

void QSpellucci::minimize(const MFunctionND&  func,
                         const CDoubleArray& inxguess,      // initial guess
                         const CStringArray& ids,           // identifiers
                         CDoubleArray&       outx) const    // result
{  
    const string method = "QSpellucci::minimize";

/*  int m = func.getNbFuncs();
    int n = func.getNbVars();

        if (m <= 0){
            throw ModelException(method, 
                                 Format::toString(m) + " number of functions to minimize");
        }

        if (n <= 0){
            throw ModelException(method, 
                                 Format::toString(n) + " number of variables to minimize");
        }

        if (inxguess.size() != n){
            throw ModelException(method, 
                                 Format::toString("Length of xguess should be %ld; got %ld",
                                                  n, inxguess.size()));
        }

        if (outx.size() != n){
            throw ModelException(method, 
                                 Format::toString("Length of outx should be %ld; got %ld",
                                                  n, outx.size()));
        }
*/
        // Construct optimization problem 
        QOptimizationProblem prob(func);

        // Get work space
        valarray<double> y(prob.SolutionDimension());
        size_t i = 0;
        for (i = 0; i < prob.SolutionDimension(); ++i)
          y[i] = inxguess[i];

        // Optimize
        mp_spellucci->Optimize(prob, y);

        // Copy optimal point to output
        for (i = 0; i < prob.SolutionDimension(); ++i)
          outx[i] = mp_spellucci->OptimalPoint()[i];
}



/** Invoked when Class is 'loaded' */
void QSpellucci::load(CClassSP& clazz)
{
    REGISTER(QSpellucci, clazz);
    SUPERCLASS(OptimizerND);
    EMPTY_SHELL_METHOD(defaultCtor);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

IObject* QSpellucci::defaultCtor()
{
    return new QSpellucci();
}


CClassConstSP const QSpellucci::TYPE = CClass::registerClassLoadMethod(
    "QSpellucci", typeid(QSpellucci), load);



//////////////////// Original Spellucci class //////////////////////////////////////////////////////////////////

void QSpellucci::Spellucci::Optimize(const OptimizationProblem &problem,
                                     const valarray<double> &y) const
{
  // Consistency check
  ASSERT1(y.size() == problem.SolutionDimension(),
    Argument("Input start point has wrong size."));

  // Work space
  int dim = y.size();
/*
  try
  {
*/
    // Construct wrapper for optimization problem
    WrapperOptimizationProblem func(problem, *this);

    mf_InitialiseWorkSpace(func);

    // Set some dimensions
    m_x.resize(dim);

    // Save input point
    m_x = y;

    // Copy start guess to internal Spellucci representation
    std::copy(&m_x[0], &m_x[dim], &x[1]);
  
    // Call the Spellucci algorithm
    mf_Solve(func);

    // Copy reached point from internal Spellucci representation
    std::copy(&x[1], &x[dim + 1], &m_x[0]);
/*
  }
  catch(...)
  {
    // Copy reached point from internal Spellucci representation nomatter
    // what that might be!
    std::copy(&x[1], &x[dim + 1], &m_x[0]);
  }
*/
}

void QSpellucci::Spellucci::mf_InitialiseOptimizer() const
{
TRY
  m_nobjeval = 0;
  m_nzeroconeval = 0;
  m_nposconeval = 0;

  NX = 0;
  NRESM = 0;
  NDUALM = 0;
  MDUALM = 0;
  
  rank = 0;
  nreset = 0;
  qpterm = 0;
  ndual = 0;
  iq = 0;
  iptr = 0;
  iqtr = 0;
  n = 0;
  nr = 0;
  nres = 0;
  nh = 0;
  ng = 0;
  inx = 0;
  singul = 0;
  ident = 0;
  eqres = 0;
  analyt = 0;
  cold = 0;
  icf = 0;
  icgf = 0;
  cfincr = 0;
  itstep = 0;
  phase = 0; 
  iterma = 0;
  ifill1 = 0; 
  lastdw = 0;
  lastup = 0;
  lastch = 0; 
  bloc = 0;
  valid = 0;
  corr = 0;
  difftype = 0;

  gfn = 0.0;
  scalm = 0.0;
  scalm2 = 0.0;
  matsc = 0.0;
  scf = 0.0;
  scf0 = 0.0;
  infeas = 0.0;
  rnorm = 0.0;
  rlow = 0.0;
  epsfcn = 0.0;
  taubnd = 0.0;
  epsmac = 0.0;
  tolmac = 0.0;
  deldif = 0.0;
  epsdif = 0.0;
  runtim = 0.0;
  optite = 0.0;
  upsi = 0.0;
  upsi0 = 0.0;
  upsi1 = 0.0;
  upsist = 0.0;
  psi = 0.0;
  psi0 = 0.0;
  psi1 = 0.0;
  psist = 0.0;
  psimin = 0.0;
  phi = 0.0;
  phi0 = 0.0;
  phi1 = 0.0;
  phimin = 0.0;
  fx = 0.0;
  fx0 = 0.0;
  fx1 = 0.0;
  fxst = 0.0;
  fmin = 0.0;
  b2n = 0.0;
  b2n0 = 0.0;
  xnorm = 0.0;
  x0norm = 0.0;
  sig0 = 0.0;
  dscal = 0.0;
  dnorm = 0.0;
  d0norm = 0.0;
  sig = 0.0;
  sigmin = 0.0;
  dirder = 0.0;
  cosphi = 0.0;
  upsim = 0.0;
  del = 0.0;
  del0 = 0.0;
  del01 = 0.0;
  delmin = 0.0;
  tau0 = 0.0;
  tau = 0.0;
  ny = 0.0;
  smalld = 0.0;
  smallw = 0.0;
  rho = 0.0;
  rho1 = 0.0;
  eta = 0.0;
  epsx = 0.0;
  c1d = 0.0;
  scfmax = 0.0;
  updmy0 = 0.0;
  tauqp = 0.0;
  taufac = 0.0;
  taumax = 0.0;
  alpha = 0.0;
  beta = 0.0;
  theta = 0.0;
  sigsm = 0.0;
  sigla = 0.0;
  delta = 0.0;
  stptrm = 0.0;
  delta1 = 0.0;
  stmaxl = 0.0;
  level = 0.0;
  clow = 0.0;
  sstr = 0.0;
  riitr = 0.0;
CATCH(QSpellucci::Spellucci::mf_InitialiseOptimizer)
}


void QSpellucci::Spellucci::mf_InitialiseQuasiNewtonUpdate(void) const
{
TRY
a = Diag(NX + 1, matsc);
  diag0 = matsc;

  accinf(itstep, 27) = -1.0;
  accinf(itstep, 14) = 1.0;

  return;
CATCH(QSpellucci::Spellucci::mf_InitialiseQuasiNewtonUpdate)
}


void QSpellucci::Spellucci::mf_Solve(const WrapperOptimizationProblem &problem) const
{
  bloc = false;
  analyt = false;
  valid = false;
  epsfcn = m_epsfcn; 
  difftype = m_difftype;
  taubnd = m_taubnd; 

  xsc = 1.0;
  xtr = 0.0;

  epsdif = 1.e-16;

  gconst = false;
  val = false;
  gresn = 1.0;
  gresn[0] = 0.0;

  inx = false;
  cold = true;

  n = problem.SolutionDimension(); 
  nh = problem.NZeroConstraint(); 
  ng = problem.NPositivityConstraint(); 
  analyt = false;
  epsdif = m_epsdif; 
  
// TODO: make nreset part of public interface
  nreset = n;
  
  // del0 and tau0
  del0 = m_del0; 
  tau0 = m_tau0; 
  tau  = 0.1;

  // gunit-array. 
  int j = 0;
  for (j = 0; j <= nh + ng; ++j) 
  {
    gunit(1, j) = -1;
    gunit(2, j) = 0;
    gunit(3, j) = 0;
  }

  /*  
  DomainBound lower = problem.LowerBound();
  DomainBoundIterator iter = lower.begin();

  for (; iter != lower.end(); ++iter) 
  {
    gunit(1, iter->first) = 1;
    gunit(2, j) = 0;
    gunit(3, j) = 0;
  }
  */

  xst = x;
  x /= xsc;

  mf_InitialiseSolver(problem);

  int i = 1;
  for (i = 1; i <= n; i++)
  {
    // ug and og have been evaluted for the original variables here we use them internally
    // for the scaled ones
    if (llow[i]) 
      ug[i] /= xsc[i];
    if (lup[i]) 
      og[i] /= xsc[i];
  }

  // preevaluation of gradients of linear functions done only once.
  for (i = 0; i <= nres; i++)
  {
    if (gunit(1, i) != 1 && gconst[i])
    {
      // evaluate gradient once.
      if (i == 0)
      {
        val[0] = true;
        gradf = NumericalGradient(MemFunRef(problem, &WrapperOptimizationProblem::ObjectFunction), 
          x, xsc, difftype, epsfcn, taubnd);
      }
      else
      {
        valarray<double> yy(NX + 1);
        val[i] = true;
        if (i <= nh) 
          yy = NumericalGradient(MemFunRefBind1st(problem, i, 
                 &WrapperOptimizationProblem::ZeroConstraint), x, xsc, difftype, epsfcn, taubnd);
        else 
          yy =NumericalGradient(MemFunRefBind1st(problem, i - nh, 
                &WrapperOptimizationProblem::PositivityConstraint), 
                x, xsc, difftype, epsfcn, taubnd);
        int j = 1;
        for (j = 1; j <= n; j++)
          gres(j, i) = yy[j];
        gresn[i] = Max(1.0, TwoNorm(yy, 1, n + 1));
      }
    }
  }

//  clock_t start = clock();
  // Call the optimizer
  mf_RecursiveQuadraticSolver(problem);
//  clock_t finish = clock();

//  runtim = static_cast<double>((finish - start) / CLOCKS_PER_SEC);

  return;
}

void QSpellucci::Spellucci::mf_InitialiseSolver(
                                                  const WrapperOptimizationProblem &problem) const
{
TRY
  int     i, j;

  double term;
  epsmac = pow(2.0, -20);
  do
  {
    epsmac /= 2.0;
    term = 1.0 + epsmac;
  }
  while (term != 1.0);
  
  epsmac = epsmac + epsmac;
  tolmac = epsmac;
  
  double tol1;
  do
  {
    tol1 = tolmac;
    tolmac /= 16.0;
  }
  while (tolmac != 0.0);

  tolmac = tol1;

  // epsmac machine precision, tolmac smallest machine number larger than 0.0 (approximately
  // base 16 for exponent therefore division by 16 assumed)

  // warning
  // on some machines the computation of tolmac may result in an error because underflow is not
  // accepted as 0.0 as is assumed here

  if (tau0 == 0.0)
    tau0 = 1.0;
  if (del0 == 0.0)
    del0 = 0.5 * tau0;

  if (nreset > n)
    nreset = n;
  if (nreset <= 4)
    nreset = 4;

  // standard initialization
  lastch = 0;
  lastdw = 0;
  lastup = 0;
  level = 1.0;
  tau = 0.1;
  iterma = m_iterma;
  epsx = m_epsx;
  sigsm = sqrt(epsmac);
  smalld = 0.1;
  smallw = exp(2.0 * log(epsmac) / 3.0);
  rho = 1.e-6;
  rho1 = 1.e-10;
  del01 = del0 / 10.0;
  delmin = Min(del01, Max(1.e-6 * del0, smallw));
  if (!analyt) 
    delmin = Min(del01, Max(epsdif, delmin));
  c1d = 1.e-2;
  scfmax = 10000.0;
  taufac = 10.0;
  taumax = pow(scfmax, 4);
  updmy0 = 0.1;
  double infiny = epsmac / tolmac;
  fx = 0.0;
  b2n = 0.0;
  b2n0 = 0.0;
  nres = ng + nh;

  if (cold) 
  {
    a = Diag(NX + 1, 1.0);
    diag0 = 1.0;
  }

  diag = 0.0;
  // initialisation of qr and gres is done in mf_InitialiseWorkSpace.
  // qr = 0;
  // gres = 0;

  ug = -infiny;
  og = infiny;
  llow = false;
  lup = false;
  delfac = 1.0;   

  double gxi;
  valarray<double> xnull(NX + 1);
  for (i = nh + 1; i <= nres; i++) 
  {
    delfac[i] = 1.0;
    // scan for real lower or upper bounds
    if (gunit(1, i) == 1) 
    {
      gxi = mf_PositivityConstraints(i - nh, xnull, problem);
      if (gunit(3, i) > 0) 
      {
        llow[gunit(2, i)] = true;
        ug[gunit(2, i)] = -gxi / gunit(3, i);
      } 
      else 
      {
        lup[gunit(2, i)] = true;
        og[gunit(2, i)] = -gxi / gunit(3, i);
      }
    }
  }

  for (i = nh + 1; i <= nres; i++) 
  {
    // modify del0, such that lower and upper bound never become binding simultaneously
    if (gunit(1, i) == 1) 
    {
      j = gunit(2, i);
      if (og[j] < infiny && ug[j] > -infiny) 
        del0 = Min(del0, (og[j] - ug[j]) * 0.1 * abs(gunit(3, i)));
    }
  }

  for (i = nh + 1; i <= nres; i++) 
  {
    // delfac corresponds to an indirect primal scaling
    if (gunit(1, i) == 1) 
    {
      j = gunit(2, i);
      if (gunit(3, i) > 0) 
      {
        delfac[i] = Max(delfac[i], Abs(ug[j]) * 0.1);
        if (og[j] < infiny)
          delfac[i] = Min(delfac[i], (og[j] - ug[j]) / (10.0 * del0));
      } 
      else 
      {
        delfac[i] = Max(delfac[i], Abs(og[j]) * 0.1);
        if (ug[j] > -infiny)
          delfac[i] = Min(delfac[i], (og[j] - ug[j]) / (10.0 * del0));
      }
    }
  }
  double bd0 = infiny;
  for (i = 1; i <= n; i++)
  {
    if (ug[i] > 0.0) 
      bd0 = Min(bd0, og[i]);
    if (og[i] < 0.0) 
      bd0 = Min(bd0, -ug[i]);
  }
  
  corr = false;
  
  // evaluate gradients of supersimple functions only once
  // a function is said to be supersimple if it is of the form a*x[j]+b
  if (gunit(1, 0) == 1) 
  {
    gconst[0] = true;
    val[0] = true;
    gradf = 0.0;
    gradf[gunit(2, 0)] = gunit(3, 0) * xsc[gunit(2, 0)];
    gfn = abs(gunit(3, 0));
  } 
  else 
  {
    val[0] = false;
    gradf = 0.0;
  }

  double hxi;
  for (i = 1; i <= nh; i++) 
  {
    if (gunit(1, i) == 1) 
    {
      // a fixed variable. corrected if necessary
      val[i] = true;
      gconst[i] = true;
      gres(gunit(2, i), i) = gunit(3, i) * xsc[gunit(2, i)];
      gresn[i] = abs(gunit(3, i)) * xsc[gunit(2, i)];
      if (gresn[i] == 0.0) 
        exit(1);

      hxi = mf_ZeroConstraints(i, xnull, problem);
      
      term = -hxi / gunit(3, i);
      if (term != x[gunit(2, i)]) 
        corr = true;
      x[gunit(2, i)] = term;
    }
  }
  
  for (i = nh + 1; i <= nres; i++) 
  {
    if (gunit(1, i) == 1) 
    {
      if (gunit(3, i) == 0) 
        exit(1);

      gxi = mf_PositivityConstraints(i - nh, x, problem);
      
      gxi = 2.0 * delmin - gxi;
      if (gxi > 0.0) 
      {
        corr = true;
        x[gunit(2, i)] = x[gunit(2, i)] + gxi / gunit(3, i);
      }
      gres(gunit(2, i), i) = gunit(3, i) * xsc[gunit(2, i)];
      gresn[i] = abs(gunit(3, i)) * xsc[gunit(2, i)];
      val[i] = true;
      gconst[i] = true;
    }
  }

  bind = 0;
  bind0 = 0;
  u = 0.0;
  u0 = 0.0;
  cres = 0;
  cgres = 0;
  for (i = 1; i <= nres; i ++) 
  {
    // initial weights of the penalty-term
    if (cold) 
      w[i] = 1.0;
    sort[i] = i;
  }   
  clow = 1.0;
  ny = 2.0;

  // scf = weight factor for objective function
  // scf0 = damping factor for tangential direction

  scf = 1.0;
  scf0 = 1.0;
  sigla = 2048.0;
  beta = 4.0;
  alpha = 0.1;
  delta1 = 0.9;
  delta = 1.e-3;
  theta = 0.9;
  icf = 0;
  icgf = 0;

  return;
CATCH(QSpellucci::Spellucci::mf_InitialiseSolver)
}


void QSpellucci::Spellucci::mf_ScalingFactors(void) const
{
TRY
  double term;

  int wlow = false;
  int i;
  for (i = 1; i <= nres; ++i) 
  {
    // w1 tentative new weights
    term = ny * Abs(u[i]) + tau;
    if (term > w[i]) 
      w1[i] = term + tau;
    else 
    {
      w1[i] = w[i];
      if (term < w[i] * 0.5 && bind[i] == 1) 
        w1[i] = (term + w[i]) * 0.5;
    }
    if (w1[i] < w[i]) 
      wlow = true;
  }
  // wlow equals true if 1.0 tentative weight at least has been decreased 
  
  double s1 = 0.0;
  double s2 = 0.0;
  for (i = 1; i <= nres; ++i) 
  {
    if (i <= nh) 
    {
      s1 += w1[i] * Abs(resst[i]);
      s2 += w1[i] * Abs(res[i]);
    } 
    else 
    {
      s1 -= Min(0.0, resst[i]) *w1[i];
      s2 -= Min(0.0, res[i]) *w1[i];
    }
  }

  double diff0 = (fxst - fx) * scf + (s1 - s2);
  if (wlow && diff0 >= eta  *clow && itstep - lastdw > Max(5, Min(20, n / 10))) 
  {
    // accept new (diminished) weights 
    
    if (clow > itstep / 10) 
      eta = 1.3 * eta;
    lastch = itstep;
    lastdw = itstep;
    level = diff0 / iterma;
    psist = s1;
    psi = s2;
    for (i = 1; i <= nres; ++i) 
      w[i] = w1[i];
    clow += 1.0;
  } 
  else 
  {
    // increase individual weights if necessary. let weigths unchanged otherwise
    
    s1 = 0.0;
    s2 = 0.0;
    for (i = 1; i <= nres; ++i) 
    {
      if (w1[i] > w[i]) 
      {
        lastup = itstep;
        lastch = itstep;
      }
      w[i] = Max(w[i], w1[i]);
      if (i <= nh) 
      {
        s1 += w[i] * Abs(resst[i]);
        s2 += w[i] * Abs(res[i]);
      } 
      else 
      {
        s1 -= w[i] * Min(0.0, resst[i]);
        s2 -= w[i] * Min(0.0, res[i]);
      }
    }
    psist = s1;
    psi = s2;
  }   
  term = 0.0;
  if (nres >= 1)
    term = w[1];
  for (i = 2; i <= nres; i++) 
    term = Max(term, w[i]);

  accinf(itstep, 20) = term;
  accinf(itstep, 19) = clow;
  
CATCH(QSpellucci::Spellucci::QSpellucci::Spellucci::mf_ScalingFactors)
}

void QSpellucci::Spellucci::mf_updateBFGS(void) const
{
TRY
  const double p2 = 0.2;
  const double p8 = 0.8;
  const double tm3 = 0.001;
  const double p5 = 0.5;

  valarray<double> dg(NX + 1);
  valarray<double> adx(NX + 1);
  valarray<double> ltdx(NX + 1);
  valarray<double> gtdx(NRESM + 1);
  valarray<double> updz(NX + 1);
  valarray<double> updx(NX + 1);

  int i;
  // multiply dx = (s in the usual notation) by Cholesky factor stored in the upper half of a
  for (i = 1;  i <= n;  i++) 
    ltdx[i] = Dot(a.Row(i), difx, i, n + 1);
  dg = gphi1;
  dg -= gphi0;

  if (TwoNorm(dg, 1, n + 1) == 0.0) 
  {
    // suppress update
    accinf(itstep, 27) = 0.0;
    accinf(itstep, 28) = 0.0;
    accinf(itstep, 29) = 0.0;

    return;
  }
  // adx = a * (x-x0), x-x0 = difx 
  for (i = 1; i <= n; i++) 
    adx[i] = Dot(a.Column(i), ltdx, 1, i + 1); 
  
  // gtdx = grad(res)(transp)*(x-x0)
  for (i = 1; i <= alist[0]; i++) 
    gtdx[i] = Dot(gres.Column(alist[i]), difx, 1, n + 1) / gresn[alist[i]]; 
  
  double ndx = TwoNorm(difx, 1, n + 1);
  double tk = Min(p5, pow(dnorm, 2));
  double anorm = 0.0;
  double term1 = Abs(a(1, 1));
  anorm = 0.0;
  int j;
  for (i = 1; i <= n; i++) 
  {
    for (j = i; j <= n; j++) 
      anorm += pow(a(i, j), 2);
    term1 = Min(term1, Abs(a(i, i)));
  }
  double acond;
  if (term1 != 0.0) 
    acond = anorm / pow(term1, 2);
  else 
    acond = epsmac / tolmac;
  double den1 = TwoNorm(ltdx, 1, n + 1);
  den1 = pow(den1, 2);
  double den2 = Dot(dg, difx, 1, n + 1);
  if (den1 <= rho1 * anorm * pow(ndx, 2) || acond >= 1.0 / rho1) 
  {
    // take a restart step 
    mf_InitialiseQuasiNewtonUpdate();
    return;
  }
  double term;
  if (nres == 0) 
  {
    // in the unconstrained case we take the Powell update 
    double th = 1.0;
    if (den2 < p2 * den1) 
    {
      th = p8 * den1 / (den1 - den2);
      dg *= th;
      dg += (1.0 - th) * adx;
      den2 = Dot(dg, difx, 1, n + 1);
    }
    term = 1.0 / sqrt(den2);
    dg *= term;
    updz = dg;
    term = 1.0 / sqrt(den1);
    updx = term * adx;
    accinf(itstep, 28) = den2 / den1;
    accinf(itstep, 29) = th;
    accinf(itstep, 27) = 2.0;
    if (th != 1.0) 
      accinf(itstep, 27) = 3.0;
  } 
  else 
  {
    double ngtdx = TwoNorm(gtdx, 1, alist[0] + 1);
    term  = 1.0 / sqrt(den1);
    updx = term * adx;
    double xsik, den21;
    if (den2 >= rho1 * Dot(dg, dg, 1, n + 1) && TwoNorm(dg, 1, n + 1) >= sqrt(epsmac) * ndx) 
    {
      xsik = 0.0;
      updz = dg;
      den21 = den2;
    } 
    else 
    {
      // try Pantoja-Mayne modification
      double den3 = tk * pow(ndx, 2) + pow(ngtdx, 2);
      if (den2 >= rho1 * Dot(dg, dg, 1, n + 1))
        xsik = 1.0;
      else 
        xsik = 1.0 + (tk * pow(ndx, 2) + Abs(den2)) / den3;
      for (i = 1; i <= n; i++) 
      {
        term = 0.0;
        for (j = 1; j <= alist[0]; j++) 
        {
          term1 = gres(i, alist[j]) * gtdx[j];
          term1 /= gresn[alist[j]];
          term += term1;
        }
        updz[i] = dg[i] + xsik * (tk * difx[i] + term);
      }
      den21 = Dot(updz, difx, 1, n + 1);
    }
    term = 1.0 / sqrt(den21);
    updz *= term;
    double th = 1.0;
    if (den2 < p2 * den1) 
    {
      th = p8 * den1 / (den1 - den2);
      dg *= th;
      dg += (1.0 - th) * adx;
      den2 = Dot(dg, difx, 1, n + 1);
    }
    term = 1.0 / sqrt(den2);
    dg *= term;
    if (TwoNorm(dg, 1, n + 1) <= tm3 * TwoNorm(updz, 1, n + 1)) 
    {
      // the Powell update produces a smaller growth 
      updz = dg;
      accinf(itstep, 28) = den2 / den1;
      accinf(itstep, 29) = th;
      accinf(itstep, 27) = 2.0;
      if (th != 1.0) 
        accinf(itstep, 27) = 3.0;
    } 
    else 
    {
      // no update if strongly irregular 
      accinf(itstep, 27) = 1.0;
      accinf(itstep, 28) = tk;
      accinf(itstep, 29) = xsik;
    }
  }
  try
  {
    mf_UpperTriangularCholeskyFactor(a, updz, updx);
  
    // check illconditioning after updating
// TODO: take care of index off-set when the refactorizing is carried out
    valarray<double> v(NX);
    for (i = 1; i <= NX; ++i)
      v[i - 1] = Abs(a(i, i));
    term = v.max();
    term1 = v.min();

    if (pow(term1, 2) <= rho1 * pow(term, 2)) 
      mf_InitialiseQuasiNewtonUpdate(); 
  }
  catch(...)
  {
    mf_InitialiseQuasiNewtonUpdate(); 
  }

CATCH(QSpellucci::Spellucci::QSpellucci::Spellucci::mf_updateBFGS)
}


void QSpellucci::Spellucci::mf_RecursiveQuadraticSolver(
                                              const WrapperOptimizationProblem &problem) const
{
TRY
  const double seven = 7.e0;
  const double p7 = .7e0;
  const double tm4 = 1.e-4;
  const double tm3 = 1.e-3;
  const double tp1 = 1.e1;
  const double tp2 = 1.e2;
  const double tp3 = 1.e3;
  const double tp4 = 1.e4;
  const double tm2 = 1.e-2;
  const double tm1 = 1.e-1;
  const double p5 = .5e0;

  static double uminsc, slackn;
  
  double  delsig, delx, sum, term, umin, term1, scfh, unorm,
            del1, clwold, eps, delold, fac, tauqp0;

  int     iumin, rank0, nr0, csdifx, nrbas, nperm, qpnew,
            etaini, l, l0, i, j, k, csssig, csirup, csreg, 
            cschgx, csmdph;

  valarray<double> qtx(NX + 1);
  valarray<double> yy(NX + 1);
  valarray<double> yx(NX + 1);
  valarray<double> trvec(NX + 1);
  valarray<int> delist(NRESM + 1);
  valarray<int> bindba(NRESM + 1);

  // initialisation 
  d = 0.0;
  d0 = 0.0;
  itstep = 0;
  alist[0] = nh;
  delist[0] = 0;
  violis[0] = 0;
  upsi = 0.0;
  psi = 0.0;
  psi0 = 0.0;
  sig0 = 0.0;
  d0norm = 1.0;
  unorm = 1.0;

  // in order to have cosphi well defined for itstep = 1
  dnorm = 1.0;
  del = del0;

  // count successive regularization steps
  csreg  = 0;
  // count successive small changes in x
  cschgx = 0;
  // count small differences of fx
  csdifx = 0;
  // count irregular quasi-Newton-updates
  csirup = 0;
  // count successive small stepsizes
  csssig = 0;
  // count successive small differences of penalty-function
  csmdph = 0;
  matsc = 1.0;

  tauqp = 1.0;
  nperm = false;
  ident = false;
  etaini = false;

  for (i = 1; i <= n; ++i) 
  {
    perm[i] = i;
    perm1[i] = i;
  }
  bind0 = 1;
  bind = 1;
  for (i = 1; i <= nh; ++i) 
    alist[i] = i;

  if (analyt) 
    eps = Min(epsx, sqrt(epsmac));
  else 
  {
    eps = epsdif;
    if (epsx < pow(epsdif, 2))
      epsx = pow(epsdif, 2);
  }
  eps = Max(epsmac * tp3, Min(tm3, eps));

  // function and gradient values
  for (i = 1; i <= nres; i++) 
  {
    if (i <= nh)
    {
      res[i] = mf_ZeroConstraints(i, x, problem);
      term = Abs(res[i]);
      if (!gconst[i])
      {
        // we assume that the gradient can be evaluated whenever the function can
        yy = NumericalGradient(MemFunRefBind1st(problem, i, 
                &WrapperOptimizationProblem::ZeroConstraint), x, xsc, difftype, epsfcn, taubnd);
        val[i] = true;
        for (j = 1; j <= n; j++) 
          gres(j, i) = yy[j];
      }
    }
    else
    {
      bind[i] = 0;
      res[i] = mf_PositivityConstraints(i - nh, x, problem);
      term = -Min(0.0, res[i]);
      if (res[i] <= delmin) 
      {
        bind[i] = 1;
        alist[0] += 1;
        alist[alist[0]] = i;
        if (!gconst[i])
        {
          val[i] = true;
          // we assume that the gradient can be evaluated whenever the function can
          yy = NumericalGradient(MemFunRefBind1st(problem, i - nh, 
                 &WrapperOptimizationProblem::PositivityConstraint), 
                 x, xsc, difftype, epsfcn, taubnd);
          for (j = 1; j <= n; j++) 
            gres(j, i) = yy[j];
        }
      }
    }
    upsi += term;
    psi += term * w[i];
    if (val[i] && !gconst[i]) 
      gresn[i] = Max(1.0, TwoNorm(yy, 1, n + 1));
  }
  
  L100:

  // obtaining a point feasible within tau0 first
  if (upsi >= tau0)
  {
    scf = 0.0;
    phase = -1;
  } 
  else 
  {
    fx = mf_ObjectFunction(x, problem);
    if (!val[0]) 
    {
      // we assume that the gradient evaluation can be done whenever the function has been evaluated
      gradf = NumericalGradient(MemFunRef(problem, &WrapperOptimizationProblem::ObjectFunction), 
                x, xsc, difftype, epsfcn, taubnd);
      val[0] = true;
    }
    scf = 1.0;
    phase = 0;
    fxst = fx;
    psist = psi;
    upsist = upsi;
    resst = res;
    eta = 0.0;
  }

  L200:

  // main iteration loop: getting a better x
  if (!ident)
  {
    itstep += 1;
    if (itstep > iterma)
    {
      optite = -3.0;
      itstep = iterma;
      b2n = accinf(itstep, 8);
      b2n0 = accinf(itstep, 7);
      return;
    }
    qpnew = false;
    qpterm = 0;
    delold = del;
    del = 0.0;
    b2n0 = -1.0;
    b2n = -1.0;
    singul = false;
    nperm = false;
    for (i = 1; i <= n; i++)
    {
      nperm = nperm || (perm[i] != perm1[i]);
      perm[i] = perm1[i];
      diag[i] = 0.0;
    }
  }

  // current valid row permutation for QR-decomposition of
  // matrix of binding gradients in order to obtain continuity
  // of the QR-decomposition is given by perm

  nr = alist[0];
  nrbas = nr;
  bindba = bind;
  for (j = 1; j <= 32; j++) 
    accinf(itstep, j) = 0.0;
  gfn = TwoNorm(gradf, 1, n + 1);

  // compute new weight of objective function if useful
  if (nres > 0 && phase >= 0 && !ident && itstep > 1 &&
      ((accinf(itstep-1, 10) == -1. && scf0 == 1.0)
      || accinf(itstep-1, 10) == 1.))
  {
    // try rescaling the objective function
    term = 0.0;
    for (i = 1; i <= nres; i++) 
    {
      if (gunit(1, i) != 1)
        term = Max(term, gresn[i]);
    }
    scfh = term / Max(1.0 / scfmax, gfn);
    if (scfh < 1.0 / scfmax)
      scfh = 1.0 / scfmax;
    if (scfh > scfmax)
      scfh = scfmax;
    if ((fxst - fx) * scfh + scfh / scf * (psist - psi) >= 
        scfh / scf * eta * clow && lastch <= itstep-4 
        && scfh < tm1*scf || scfh > tp1*scf)
    {
      // rescale the objective function if this seems promising and the change is significant
      clow += 1;
      term = scfh / scf;
      psi *= term;
      psist *= term;
      u *= term;
      unorm *= term;
      scf = scfh;
      lastch = itstep;
      term = sqrt(term);
      for (i = 1; i <= n; i++) 
      {
        diag0[i] *= term;
        for (j = 1; j <= n; j++) 
          a(j, i) *= term;
      }
      matsc *= term;
    }
  }
  accinf(itstep, 1) = itstep;
  accinf(itstep, 2) = fx;
  accinf(itstep, 3) = scf;
  accinf(itstep, 4) = psi;
  accinf(itstep, 5) = upsi;

  // begin solver

  // QR-decomposition of matrix of binding gradients
  if (nr >= 1)
    mf_QRDecompBindingConstraints(1, nr);
  else
    rank = 0;

  term = mf_SolveCholeskyFactorLeft(a, gradf, yy, n);
    
  for (i = 1; i <= n; i++) 
    qgf[i] = yy[perm[i]];
  
  mf_HouseholderTransformation(1, 0, 1, rank, n, qr, betaq, qgf, trvec);
  
  qgf = trvec;
    
  // compute del as function of x (forcing infeasibility and the projected gradient to 0.0)
  b2n0 = TwoNorm(qgf, rank + 1, n + 1);
    
  sum = 0.0;
  for (i = 1; i <= nres; i++) 
  {
    if (i <= nh) 
      sum += Abs(res[i]) / gresn[i];
    else 
      sum -= Min(0.0, res[i]) / gresn[i];
  }
  if (itstep > 1 && accinf(itstep - 1, 8) >= 0.0 && !etaini && accinf(itstep-1, 18) >= 0)
  {
    etaini = true;
    eta = (accinf(itstep - 1, 8) / Max(1.0, gfn) + sum + Min(1.0, slackn) + 
            Min(1.0, Abs(uminsc))) / Min(30 * n, iterma);
    level = eta;
  }
  delx = delmin;
  term = scf * (fx0 - fx) + psi0 - psi;
  if (term > 0.0 && scf != 0.0) 
    delx = Max(delx, exp(p7 * p7 * log(term)));
  if (scf == 0.0) 
    delx = Min(del0 * tm4, Max(delx, upsi * tm2));
  delsig = delmin;
  
  // del should be large enough to include constraints hit in
  // step before violis comes from unimin
  
  for (i = 1; i <= violis[0]; i++) 
  {
    j = violis[i];
    delsig = Max(delsig, res[j] / gresn[j] / delfac[j] * (1.0 + tm1));
  }
  del = Min(del0, Max(Min(delsig, 5.0 * delx), delx));
  if (violis[0] == 0) 
    del = Min(del, del01);
  
  // if phase = 2 don't loose a binding constraint. phase = 2 implies delist[0] = 0
  if (phase == 2 && violis[0] == 0) 
  {
    for (i = nh + 1; i <= nres; i++) 
      if (bind0[i] == 1) 
        del = Min(del01, Max(del, Abs(res[i]) / gresn[i]));
  }
  // reduce del by a fixed factor (tm2) if useful, that is if delete-list in the 
  // previous step was not empty
  
  term = del;
  for (i = 1; i <= delist[0]; i++) 
  {
    j = delist[i];
    term1 = res[j] / gresn[j] * (1.0 - tm2) / delfac[j];
    if (term1 >= del * tm2) 
      term = Min(term, term1);
  }
  del = term;
  
  // if delta becomes too large, we may loose complementary slackness */
  if (itstep > 1 && ! ident && scf != 0.0) 
  {
    term = 0.0;
    for (i = nh + 1; i <= nres; i++) 
      term = term + Max(0.0, res[i] / gresn[i] - delmin) * Max(0.0, u[i] - smallw) / gresn[i];
    if (term > 0.0) 
    {
      for (i = nh + 1; i <= nres; i++) 
      {
        if (u[i] > smallw && res[i] / gresn[i] > delmin) 
          del = Max(delmin, Min(del, res[i] / gresn[i] * (1.0 - tm2) / delfac[i]));
      }
    }
  }

  // if the current step was singular and not successful, try a greater del the same, if qpterm 
  // in the last step did signal trouble, but stepsize selection was nevertheless successful
    
  if (itstep > 1 && accinf(itstep - 1, 30) < 0.0) del = Min(tp1 * delold, del0);
    
  // include nearly binding inequality constraints
  for (i = nh + 1; i <= nres; i++)
  {
    term = res[i] / gresn[i];
    if (bind[i] == 0 && term <= del*delfac[i])
    {
      // it may be useful to include constraints>0 if near its boundary
      // but avoid it, if the constraint was in the old delete-list
      bind[i] = 1;
      alist[0] += 1;
      alist[alist[0]] = i;
      if (!val[i]) 
      {
        val[i] = true;
        yy = NumericalGradient(MemFunRefBind1st(problem, i - nh, 
               &WrapperOptimizationProblem::PositivityConstraint), 
               x, xsc, difftype, epsfcn, taubnd);
        for (j = 1; j <= n; j++) 
          gres(j, i) = yy[j];
        gresn[i] = Max(1.0, TwoNorm(yy, 1, n + 1));
      }
    }
  }
  rank0 = rank;
  nr0 = nr;
  nr = alist[0];
  
  mf_QRDecompBindingConstraints(nr0 + 1, nr);
  
  mf_HouseholderTransformation(1, 0, rank0 + 1, rank, n, qr, betaq, qgf, trvec);
  
  qgf = trvec;
  yy = -qgf * scf;
  yu = 0.0;
  
  // first computation of lagrangian multipliers u
  
  mf_SolveTriangularSystem(1, rank, yy, yu);
  
  umin = 0.0;
  unorm = 0.0;
  u = 0.0;
  iumin = 0;
  uminsc = 0.0;
  for (i = 1; i <= rank; i++) 
  {
    unorm = Max(unorm, Abs(yu[i]));
    k = alist[colno[i]];
    u[k] = -yu[i];
    if (k > nh) 
    {
      if (-yu[i] / gresn[k] < uminsc) 
      {
        iumin = k;
        uminsc = -yu[i] / gresn[k];
      }
    }
  }
  if (scf != 0.0) 
  {
    for (i = 1; i <= n; i++) 
    {
      yx[i] = scf * gradf[i];
      for (j = 1; j <= nres; j++) 
        yx[i] = yx[i] - gres(i, j) * u[j];
    }
    b2n = TwoNorm(yx, 1, n + 1) / scf;
  } 
  else 
    b2n = -1.0;
  
  // compute new delta 
    
  del1 = del;
// TODO: is this really what we want? Needed to avoid exception in Griewank test example
  if (b2n > 0.0 /*b2n >= 0.0*/) del1 = Max(del, tm1 * Min(del0,
    exp(p7 * log(Abs(b2n) / (gfn + 1.0) + Max(0.0, -uminsc) + sum))));
                
  // exclude constraints which were candidates for inactivating in the previous step if useful
    
  for (i = 1; i <= delist[0]; i++) 
  {
    j = delist[i];
    term1 = res[j] / gresn[j] * (1.0 - tm2) / delfac[j];
    if (term1 >= del1 * tm2) 
      del1 = Max(delmin, Min(del1, term1));
  }
  slackn = 0.0;
  for (i = nh + 1; i <= nres; i++) 
    slackn += Max(0.0, res[i] / gresn[i] - delmin) * Max(0.0, u[i] - smallw) / gresn[i];
  if (upsi <= delmin && b2n <= epsx * (gfn + 1.0) && b2n != -1.0  && uminsc >= -smallw &&
      slackn <= delmin * smallw * nres) 
  {
      // sufficient accuracy in Kuhn-Tucker conditions 
      optite = 0.0;
      return;
  }
  // include additional constraints if necessary
  
  l0 = alist[0];
  for (i = 1; i <= nres; i++) 
  {
    term = res[i] / gresn[i] / delfac[i];
    if (term > del && term <= del1 && bind[i] == 0) 
    {
      bind[i] = 1;
      alist[0] += 1;
      alist[alist[0]] = i;
      if (! val[i]) 
      {
        val[i] = true;
        yy = NumericalGradient(MemFunRefBind1st(problem, i - nh, 
               &WrapperOptimizationProblem::PositivityConstraint), 
               x, xsc, difftype, epsfcn, taubnd);
        for (j = 1; j <= n; j++) 
          gres(j, i) = yy[j];
        gresn[i] = Max(1.0, TwoNorm(yy, 1, n + 1));
      }
    }
  }
  del = del1;
  accinf(itstep, 6) = del;
  accinf(itstep, 7) = b2n0;
  accinf(itstep, 9) = alist[0];
  accinf(itstep, 10) = -1.;
  nr = alist[0];
  if (l0 != nr) 
  {
    rank0 = rank;
    mf_QRDecompBindingConstraints(l0 + 1, nr);
    mf_HouseholderTransformation(1, 0, rank0 + 1, rank, n, qr, betaq, qgf, trvec);
    qgf = trvec;
  }
  if (rank != nr) 
    goto L400;

  // second solution for multipliers, rank may have changed!
  yy = -qgf * scf;
  yu = 0.0;

  mf_SolveTriangularSystem(1, rank, yy, yu);
  
  // remember the column interchanges in qr! yu[i] corresponds to u[alist[colno[i]]]
  umin  = 0.0;
  unorm = 0.0;
  u = 0.0;
  iumin = 0;
  uminsc = 0.0;
  for (i = 1; i <= rank; i++) 
  {
    unorm = Max(unorm, Abs(yu[i]));
    k = alist[colno[i]];
    u[k] = -yu[i];
    if (k > nh) 
    {
      umin = Min(umin, -yu[i]);
      if (-yu[i] / gresn[k] < uminsc) 
      {
        iumin = k;
        uminsc = -yu[i] / gresn[k];
      }
    }
  }
  if (scf != 0.0) 
  {
    for (i = 1; i <= n; i++) 
    {
      yx[i] = scf * gradf[i];
      for (j = 1; j <= nres; j++) 
        yx[i] = yx[i] - gres(i, j) * u[j];
    }
    b2n = TwoNorm(yx, 1, n + 1) / scf;
  }
  accinf(itstep, 8) = b2n;
  accinf(itstep, 11) = umin;
  
  delist[0] = 0;
  if (phase >= 0 && b2n != -1.0) 
  {
    if (Abs(uminsc) >= 
        Max(smallw, Abs(b2n) / (gfn + 1.0) * c1d)) 
    {
      for (i = nh + 1; i <= nr; i++) 
      {
        k = alist[colno[i]];
        if (-yu[i] / gresn[k] <= -smallw) 
        {
          delist[0] += 1;
          delist[delist[0]] = k;
        }
      }
    }
  }
  
  // the new delist doesn't influence the current d but only the  computation of the next del
  eqres = true;
  for (i = 1; i <= nres; i ++) 
    eqres = eqres && (bind[i] == bind0[i]);
  // compute condition number estimators of diag-r and diag of Cholesky-decomposition of b
  
  if (nr > 1) 
  {
    term = 0.0;
    term1 = 1.0;
    for (i = 1; i <= nr; i++) 
    {
      term  = Max(term, Abs(diag[i]));
      term1 = Min(term1, Abs(diag[i]));
    }
    accinf(itstep, 13) = term / term1;
  } 
  else if (nr == 1)
    accinf(itstep, 13) = 1.;
  else 
    accinf(itstep, 13) = -1.;
  term = Abs(a(1, 1));
  term1 = Abs(a(1, 1));
  i = 2;
  while (i <= n) 
  {
    term = Max(term, Abs(a(i, i)));
    term1 = Min(term1, Abs(a(i, i)));
    i += 1;
  }
  accinf(itstep, 14) = pow(term / term1, 2);
  
  // since a represents the Cholesky-factor, this square
  
  slackn = 0.0;
  for (i = nh + 1; i <= nres; i++) 
    slackn = slackn + Max(0.0, res[i] / gresn[i] - delmin) * Max(0.0, u[i] - smallw) / gresn[i];
  if (umin >= -smallw &&
        slackn <= delmin * smallw * nres &&
        upsi <= nres * delmin && upsi0 <= nres * delmin
        && Abs(fx - fx0) <= eps * (Abs(fx) + 1.0) &&
        b2n != -1.0 &&
        b2n <= tp2 * epsx * (gfn + 1.0)) 
  {
    csdifx += 1;
  } 
  else 
    csdifx = 0;
  
  if (phase >= 0 && (accinf(itstep, 14) > tp3 || !analyt) && csdifx > n) 
  {
    optite = 4.0;
    // to avoid possible slow convergence with singular projected hessian or 
    // inaccurate numerical gradients 
    return;
  }
  // compute damping factor for tangential component if upsi>tau0 / 2

  scf0 = 1.0;
  if (phase >= 0 && upsi > tau0*p5)
    scf0 = Max(1.0 / scfmax, (2.0 * (tau0 - upsi) / tau0) * upsi * tm1 / Max(1.0, gfn)) / scf;
  accinf(itstep, 15) = scf0;

  // compute tangential component
  
  for (i = nr + 1; i <= n; i++) 
    qtx[i] = yy[i] * scf0;
  // qtx[nr + 1],.., qtx[n] is s2
    
  // compute right hand side and vertical component use damping for inactivation direction if 
  // very large no indirect inactivation if infeasibility large and we are not almost stationary
  // on the current manifold 
  fac = 1.0;
  if (-umin * c1d > b2n + upsi && b2n != -1.0) 
    fac = c1d;
  if (upsi > tau0*p5) 
    fac = 0.0;
  for (i = 1; i <= nr; i++) 
  {
    k = alist[colno[i]];
    term = res[k];
    if (k > nh && -yu[i] < 0.0 && term > 0.0) 
      term = -term;
    if (k > nh && -yu[i] < 0.0) 
      term -= yu[i] * fac;
    yx[i] = -term;
  }
  mf_SolveTransposedTriangularSystem(1, nr, yx, qtx);
  
  // qtx is transformed direction of descent for phi
  
  mf_HouseholderTransformation(-1, 0, 1, nr, n, qr, betaq, qtx, yx);
  
  for (i = 1; i <= n; i++) 
    qtx[perm[i]] = yx[i];
  
  // solve l(transp)*d = qtx, l = a = Cholesky-factor of b
  
  term = mf_SolveCholeskyFactorRight(a, qtx, d, n);
    
  // end solver
    
  // compute new penalty weights : regular case
  clwold = clow;
  if (phase >= 0) 
    mf_ScalingFactors();
  if (clow > clwold) 
  {
    // tau_qp depends on the (new) weights
    term = w[1];
    for (i = 1; i <= nres; i++) 
      term = Max(term, w[i]);
    tauqp = Max(1.0, Min(tauqp, term));
  }
  // compute parameter phase and stopping criterion
  if (uminsc < -smallw) 
    phase = Min(1, phase);
  if (!eqres)
    phase = Min(0, phase);
  if (eqres && upsi < tau0) 
    phase = Max(1, phase);
  
  // rescale and project d if appropriate
  mf_Cut_d();
  
  // compute the directional derivative dirder
  mf_DerivativePenaltyFunction();
  
  // terminate if correction is small
  if (dnorm <= epsx * (xnorm + epsx) && upsi <= delmin
      && b2n != -1.0 && uminsc >= -smallw && b2n <= epsx * (gfn + 1.0)) 
  {
    optite = 1.0;
    return;
  }
  L350:

  // reenter from the singular case: dirder has been computed already
    
  accinf(itstep, 16) = xnorm;
  accinf(itstep, 17) = dnorm;
  accinf(itstep, 18) = phase;
    
  // compute stepsize
  cfincr = icf;

  // if no descent direction is obtained, check whether restarting the method might help
  if (dirder >= 0.0) 
  {
    // no direction of descent
    stptrm = -2.0;
    sig    = 0.0;
    
    goto L360;
  }
  // if directional derivative correct but very small, terminate
  // since no further progress might be possible
  
  if (-dirder <= epsmac * tp2 * (scf * Abs(fx) + psi + 1.0))
  {
    if (upsi > delmin * nres) 
    {
      optite = -1.0;
      stptrm = -1.0;
    } 
    else 
    {
      optite = 2.0;
      stptrm = 1.0;
    }
    sig = 0.0;
    
    goto L360;
  }
  // phase = 2 : we may hope to obtain superlinear convergence
  // switch to Maratos-correction is on then
  // return to phase = 1 if first order correction large
  
  if (phase >= 1 && dnorm <= smalld * (xnorm + smalld) && scf0 == 1.0 
    && uminsc >= -smallw && ! singul)
  {
    phase = 2;
  }   
  // return to phase 1 since correction large again
  
  if (phase == 2 && dnorm > (xnorm + smalld)) 
    phase = 1;

  mf_MaxStepSize();
  
  // stmaxl is the maximal stepsize such that point on projected ray changes with sigma, 
  // but sigla at most
  
  dd = 0.0;
  x1 = x;
  x1 += d;

  // compute second order correction of infeasibility if useful
  if (phase == 2 && dnorm > xnorm * sqrt(epsmac)) 
  {
    // only function values, from xtr = xsc*x1
    for (i = 1; i <= alist[0]; i++) 
    {
      yx[i] = 0.0;
      if (i <= nh && !gconst[alist[i]]) 
      {
        try
        {
          yx[i] = mf_ZeroConstraints(i, x1, problem);
        }
        catch(...)
        {
          goto L355;
        }
      } 
      else 
      {
        if (!gconst[alist[i]]) 
        {
          try
          {
            yx[i] = mf_PositivityConstraints(alist[i] - nh, x1, problem);
          }
          catch(...)
          {
            goto L355;
          }
        }
      }
      yx[i] = -yx[i];
    }
    for (i = 1; i <= alist[0]; i++) 
      yy[i] = yx[colno[i]];

    mf_SolveTransposedTriangularSystem(1, nr, yy, dd);
    
    for (i = nr + 1; i <= n; i++) 
      dd[i] = 0.0;

    mf_HouseholderTransformation(-1, 0, 1, nr, n, qr, betaq, dd, yx);
    
    for (i = 1; i <= n; i++) 
      dd[perm[i]] = yx[i];
      
    term = mf_SolveCholeskyFactorRight(a, dd, dd, n);
    
    if (sqrt(term) > p5*dnorm) 
    {
      // second order correction almost as large as first order 1.0: not useful
      dd = 0.0;
    }
  }
  L355:

  sig = Min(1.0, stmaxl);
  
  mf_StepSize(sig, problem);

  L360:

  cfincr = icf - cfincr;
  
  // count successive small steps
  term = scf * (fx0 - fx) + psi0 - psi;
  if (Abs(term) <= epsmac * tp3 * (scf * Abs(fx) + psi)) 
    csmdph += 1;
  else 
    csmdph = 0;
  // csmdph counts contiguous small differences of penalty function phi
  
  if (csmdph > Max(n, 10)) 
  {
    optite = seven;
    return;
  }
  if (sig <= 5.0*tm2) 
  {
    if (sig0 <= 5.0 * tm2) 
      csssig += 1;
  } 
  else 
    csssig = 0;
  
  // csssig counts the number of successive small sig's
  accinf(itstep, 21) = sig;
  accinf(itstep, 22) = cfincr;
  accinf(itstep, 23) = dirder;
  accinf(itstep, 24) = dscal;
  accinf(itstep, 25) = cosphi;
  accinf(itstep, 26) = violis[0];
    
  // no further significant progress possible
  if (sig == 0.0 && stptrm == 1.0 && optite == 2.0) 
    return;
  if (stptrm == 1.0 && sig <= tm4 && accinf(itstep, 13) > tp4 && ! singul && nres > 0) 
  {
    // try a regularized step, hopefully this will give a better d and larger sig
    if (accinf(itstep, 14) > tp4) 
      mf_InitialiseQuasiNewtonUpdate();
    ident = true;
    singul = true;
    
    goto L400;
  }
  if (stptrm < 0.0) 
  {
    // stepsize selection failed
    if (!ident) 
    {
      // try restart with a = identity scaled
      ident = true;
      delist[0] = 0;
      violis[0] = 0;
      csreg = 0;
      csssig = 0;
      csirup = 0;
      
      mf_InitialiseQuasiNewtonUpdate();
      
      alist[0] = nrbas;
      bind = bindba;
      if (upsi >= tau0) 
        goto L100;
      else 
        goto L200;
    }
    if (!singul && ident && accinf(itstep, 13) > tp4 && nres > 0) 
    {
      // try the full SQP-direction. this may be the third try for this point
      singul = true;
      ident = true;
      goto L400;
    }
    if (stptrm == -2.0) 
    {
      optite = -4.0;
      return;
    }
    if (sig == 0.0 && optite == -1.0) 
      return;
    
    // unidimensional search unsuccessfully terminated
    optite = -2.0;
    return;
  }
  if (singul && itstep > n && Abs(fx - fx0) <= 
      eps * (Abs(fx) + 1.0) && phase >= 0
      && upsi <= nres * delmin && upsi0 <= nres * delmin
      && slackn <= delmin * smallw * nres && infeas <= upsi && !ident) 
  {
    // since multipliers may be incorrect for infeas != 0.0 be careful
    optite = 4.0;
    // avoid slow progress in case of singular constraints
    return;
  }

  // relaxed termination criteria in the singular case
  if (singul && upsi <= delmin * nres && upsi0 <= delmin * nres
      && b2n != -1.0 && b2n <= (gfn + 1.0) * epsx * tp2 && phase >= 0
      && slackn <= delmin * smallw * nres && infeas <= upsi) 
  {
    // since multipliers may be incorrect for infeas != 0.0 be careful
    optite = 3.0;
    return;
  }
  
  k = 0;
  for (i = 1; i <= n; i++) 
    if (Abs(difx[i]) >= epsx * (Abs(x[i]) + tm2)) 
      k = 1;
  if (k == 0)
    cschgx += 1;
  else 
    cschgx = 0;
  if (cschgx > nreset && singul) 
  {
    // very slow progress in x in the singular case. terminate 
    optite = 5.0;
    return;
  }
  // new value of x has been accepted
    
  xnorm = TwoNorm(x, 1, n + 1);
  ident = false;
  
  mf_GradientLagrangian(gphi0);
  
  for (i = 0; i <= nres; i++) 
  {
    if (! gconst[i]) 
      val[i] = false;
  }
  
  // evaluate gradients only, since function values are already valid from unidimensional 
  // minimization argument is xtr = xsc*x
    
  if (phase >= 0 && !gconst[0]) 
  {
    val[0] = true;
    gradf = NumericalGradient(MemFunRef(problem, &WrapperOptimizationProblem::ObjectFunction), 
              x, xsc, difftype, epsfcn, taubnd);
  }
  for (i = 1; i <= alist[0]; i++) 
  {
    l = alist[i];
    if (!val[l]) 
    {
      val[l] = true;
      if (l <= nh) 
        yx = NumericalGradient(MemFunRefBind1st(problem, l, 
               &WrapperOptimizationProblem::ZeroConstraint), x, xsc, difftype, epsfcn, taubnd);
      else 
        yx = NumericalGradient(MemFunRefBind1st(problem, l - nh, 
               &WrapperOptimizationProblem::PositivityConstraint), 
               x, xsc, difftype, epsfcn, taubnd);
      for (j = 1; j <= n; j++) 
        gres(j, l) = yx[j];
      gresn[l] = Max(1.0, TwoNorm(yx, 1, n + 1));
    }
  }
  mf_GradientLagrangian(gphi1);
  
  yx = x;
  yx -= x0;
  yy = gphi1;
  yy -= gphi0;

  // since a represents the Cholesky-factor, this sqrt
  
  term = sqrt(TwoNorm(yy, 1, n + 1) / TwoNorm(yx, 1, n + 1));
  if (term != 0.0 && phase >= 0) 
    matsc = Max(1.0 / scfmax, Min(scfmax, term / 2));
  
  // current scaling of identity in case of restart
  alist[0] = 0;
  for (i = 1; i <= nres; i++) 
  {
    u0[i] = u[i];
    bind0[i] = bind[i];
    res0[i]  = res[i];
    if (i <= nh) 
      bind[i] = 1;
    else 
    {
      bind[i] = 0;
      if (res[i] / gresn[i] <= delmin) 
        bind[i] = 1;
      if (!val[i] && bind[i] == 1) 
      {
        val[i] = true;
        yx = NumericalGradient(MemFunRefBind1st(problem, i - nh, 
               &WrapperOptimizationProblem::PositivityConstraint), 
               x, xsc, difftype, epsfcn, taubnd);
        for (j = 1; j <= n; j++) 
          gres(j, i) = yx[j];
        gresn[i] = Max(1.0, TwoNorm(yx, 1, n + 1));
      }
    }
    if (bind[i] == 1) 
    {
      alist[0] += 1;
      alist[alist[0]] = i;
    }
  }

  if (scf != 0.0) 
  {
    if (csirup > nreset || csssig > nreset || csreg  > nreset) 
    {
      csreg  = 0;
      csssig = 0;
      csirup = 0;
   
      mf_InitialiseQuasiNewtonUpdate();
    } 
    else 
      mf_updateBFGS();
  }
  // proceed

  if (accinf(itstep, 27) == 1.0) 
  {
    if (itstep > 1 && accinf(itstep - 1, 29) != 0.0 && accinf(itstep, 29) != 0.0) 
        // count successive ir regular up dates 
      csirup += 1;
    else 
      csirup = 0;
  }
    
  if (phase == -1) 
    goto L100;
  else
    goto L200;

  L400:

  singul = true;
  phase  = Min(phase, 0);
  accinf(itstep, 10) = 1.0;
  
  // try to compute a descent direction using an extended quadratic program with individual
  // slack variable for any constraint compute damping factor for tangential component 
  // if upsi>tau0 / 2 by rescaling f if possible
    
  scf0 = 1.0;
  if (phase >= 0 && upsi > tau0 * p5) 
  {
    scfh = Max(1.0 / scfmax, 
               Min(scfmax, (2.0 * (tau0 - upsi) / tau0) * upsi * tau / Max(1.0, gfn)));
    if ((fxst - fx) * scfh + scfh / scf * (psist - psi) >= 
         scfh / scf * eta * clow && lastch <= itstep - 4
         && scfh < tm1 * scf || scfh > tp1 * scf)
    {
      // rescale the objective function if this seems promising and the change is significant
      clow += 1;
      term = scfh / scf;
      scf0 = term;
      psi *= term;
      psist *= term;
      u *= term;
      unorm *= term;
      scf = scfh;
      lastch = itstep;
      accinf(itstep, 15) = scf;
      term = sqrt(term);
      for (i = 1; i <= n; i++) 
      {
        diag0[i] *= term;
        for (j = 1; j <= n; j++) 
          a(j, i) *= term;
      }
      matsc *= term;
    }
  }
  // slack is upsi at most
  
  accinf(itstep, 32) = upsi;

  accinf(itstep, 13) = -1.0;
  term = Abs(a(1, 1));
  term1 = term;
  i = 2;
  while (i <= n) 
  {
    term = Max(term, Abs(a(i, i)));
    term1 = Min(term1, Abs(a(i, i)));
    i += 1;
  }
  accinf(itstep, 14) = pow(term / term1, 2);
  clwold = clow;
  
  // save for restart
  
  tauqp0 = tauqp;
  u = 0.0;
  dd = 0.0;
  mf_ExtendedRecursiveQuadraticSolver();
    
  if (dnorm == 0.0 && qpterm == 1 && optite == 3.0) 
    return;
  if (dnorm <= epsx * (Min(xnorm, 1.0) + epsx) && qpterm < 0) 
  {
    // may be it failed because of illconditioning 
    if (upsi >= nres * delmin && qpnew) 
    {
      // restarting the method has been done already: game is over
      optite = -1.0;
      return;
    }
    if (qpnew) 
    {
      optite = qpterm - 5.0;
      return;
    }
    // try a = id
    qpnew = true;
    w = 1.0;
    lastch = itstep;
    
    mf_InitialiseQuasiNewtonUpdate();
        
    ident = true;
    tauqp    = 1.0;
    alist[0] = nrbas;
    bind = bindba;
    if (scf == 0.0) 
      goto L100;
    else
      goto L200;
  }
  accinf(itstep, 11) = 0.0;
  delist[0] = 0;
  umin = 0.0;
  
  // b2n is defined also internally in mf_ExtendedRecursiveQuadraticSolver
  
  if ((qpterm >= 0 || qpterm == -3) && scf != 0.0) 
  {
    unorm = Abs(u[1]);
    for (i = 2; i <= nres; i++) 
      unorm = Max(unorm, Abs(u[i]));
    for (i = 1; i <= n; i++) 
    {
      yx[i] = scf * gradf[i];
      for (j = 1; j <= nres; j++) 
        yx[i] -= gres(i, j) * u[j];
    }
    b2n = TwoNorm(yx, 1, n + 1) / scf;
  } 
  else 
  {
    b2n = -1.0;
    // signals "undefined" here
  }
  
  // to avoid termination
  if (b2n == -1.0) 
    b2n = epsmac / tolmac;
  if (qpterm >= 0 && dnorm <= tm2 * epsx * (epsx + Min(1.0, xnorm))) 
  {
    if (upsi <= nres * delmin) 
      optite = 6.0;
    else 
      optite = -5.0;
    return;
  }
  // check whether QPsolver terminated unsuccessfully
  
  if (qpterm < 0) 
  {
    // we have a unfeasible solution for QP. try it but don't complain if it fails
  }
  if (clow > clwold) 
  {
    term = 1.0;
    for (i = 1; i <= nres; i++) 
      term = Max(term, w[i]);
    tauqp = Max(1.0, term);
  }
  if (tauqp > pow(taufac, 3) * tauqp0) 
    tauqp = pow(taufac, 3) * tauqp0;
  
  // no change of tauqp otherwise
  b2n0 = b2n;
  umin = 0.0;
  uminsc = 0.0;
  if (qpterm >= 1)
  {
    slackn = 0.0;
    for (i = nh + 1; i <= nres; i++) 
      slackn += Max(0.0, res[i] / gresn[i] - delmin) * Max(0.0, u[i] - smallw) / gresn[i];
  } 
  else 
  {
    // slack is undefined, since multipliers are undefined use this value to prevent 
    // premature termination
    slackn = 1.0;
  }
  goto L350;
CATCH(QSpellucci::Spellucci::mf_RecursiveQuadraticSolver)
}


void QSpellucci::Spellucci::mf_DerivativePenaltyFunction(void) const
{
TRY
  dirder = Dot(gradf, d, 1, n + 1) * scf;
  
  int i;
  for (i = 1; i <= nres; ++i) 
  {
    double term = Dot(gres.Column(i), d, 1, n + 1) * w[i];
    double term1 = res[i];
    if (i <= nh) 
    {
      if (term1 / gresn[i] <= -1000.0 * epsmac) 
        dirder -= term;
      else 
      {
        if (term1 / gresn[i] > 1000.0 * epsmac) 
          dirder += term;
        else 
          dirder += Abs(term);
      }
    } 
    else 
    {
      if (bind[i] == 1) 
      {
        if (Abs(term1) / gresn[i] <= 1000.0 * epsmac) 
          dirder -= Min(0.0, term);
        else 
        {
          if (term1 / gresn[i] < -1000.0 * epsmac) 
          {
            if (term > 0.0) 
              term = Min(term, -res[i] * w[i]);
            // only negative values of the constraints contribute to the penalty function
            dirder -= term;
          }
        }
      }
    }
  }
CATCH(QSpellucci::Spellucci::mf_DerivativePenaltyFunction)
}


void QSpellucci::Spellucci::mf_Cut_d(void) const
{
TRY
  xnorm = TwoNorm(x, 1, n + 1);
  double term = beta * (xnorm + 1.0);
  dnorm = TwoNorm(d, 1, n + 1);
  d0norm = TwoNorm(d0, 1, n + 1);
  dscal = 1.0;

  if (dnorm * d0norm != 0.0)
    cosphi = Dot(d, d0, 1, n + 1) / (d0norm * dnorm);
  else
    cosphi = 0.0;

  if (dnorm > term) 
  {
    // d too long: rescale
    double term1 = term / dnorm;
    dnorm = term;
    dscal = term1;
    int i;
    for (i = 1; i <= n; ++i)
    {
      d[i] *= term1;
      dd[i] *= pow(term1, 2);
    }
  }
  // since we project the ray with respect to the bounds, be sure to compute the directional
  // derivative correctly therefore correct d and dd appropriately

  int i;
  for (i = 1; i <= n; ++i)
  {
    if (llow[i] && x[i] + sigsm * d[i] <= ug[i])
    {
      d[i]  = 0.0;
      dd[i] = Max(0.0, dd[i]);
    }
    if (lup[i] && x[i] + sigsm * d[i] >= og[i])
    {
      d[i] = 0.0;
      dd[i] = Min(0.0, dd[i]);
    }
  }
  dnorm = TwoNorm(d, 1, n + 1);
CATCH(QSpellucci::Spellucci::mf_Cut_d)
}

void QSpellucci::Spellucci::mf_MaxStepSize(void) const
{
TRY
  int i;
  int exis = true;
  for (i = 1; i <= n; i++) 
    exis = exis && ((d[i] == 0.0) || (lup[i] && d[i] > 0.0) || (llow[i] && d[i] < 0.0));

  if (exis) 
  {
    stmaxl = sigsm;
    for (i = 1; i <= n; ++i)
    {
      if (llow[i] && d[i] < 0.0)
      {
        if (-d[i] * sigla >= x[i] - ug[i])
          stmaxl = Max(stmaxl, (x[i] - ug[i]) / (-d[i]));
        else 
          stmaxl = sigla;
      }
      if (lup[i] && d[i] > 0.0)
      {
        if (d[i] * sigla >= og[i] - x[i])
          stmaxl = Max(stmaxl, (og[i] - x[i]) / d[i]);
        else
          stmaxl = sigla;
      }
    }
  }
  else
    stmaxl = sigla;

  // but never use stepsize larger than sigla
  stmaxl = Min(sigla, stmaxl);
CATCH(QSpellucci::Spellucci::mf_MaxStepSize)
}

void QSpellucci::Spellucci::mf_RestoreCurrentMinimum(void) const
{
TRY
  phi1  = phimin;
  psi1  = psimin;
  upsi1 = upsim;
  sig   = sigmin;
  fx1   = fmin;

  x1 = xmin;
  res1 = resmin;

CATCH(QSpellucci::Spellucci::mf_RestoreCurrentMinimum)
}

void QSpellucci::Spellucci::mf_SaveCurrentMinimum(void) const
{
TRY
  phimin = phi1;
  upsim = upsi1;
  psimin = psi1;
  fmin = fx1;
  sigmin = sig;

  xmin = x1;
  resmin = res1;

CATCH(QSpellucci::Spellucci::mf_SaveCurrentMinimum)
}

double QSpellucci::Spellucci::mf_Evaluate(double sigact,
                                const WrapperOptimizationProblem &problem) const
{
TRY
  const double tp3 = 1.e3;
  
  sig = sigact;

  int i;
  for (i = 1; i <= n; ++i) 
  {
    x1[i] = x[i] + sig * (d[i] + sig * dd[i]);
    
    // project with respect to the box-constraints 
    if (llow[i]) 
      x1[i] = Max(x1[i], ug[i]);
    if (lup[i]) 
      x1[i] = Min(x1[i], og[i]);
  }

  double sigres = sig;
  upsi1 = 0.0;
  psi1 = 0.0;
  
  // only function values, from xtr = x1*xsc 
  int j;
  double term;
  for (j = 1; j <= nres; ++j) 
  {
    i = sort[j];
    if (i <= nh) 
    {
      res1[i] = mf_ZeroConstraints(i, x1 , problem);
      term = Abs(res1[i]);
    } 
    else 
    {
      res1[i] = mf_PositivityConstraints(i - nh, x1, problem);
      term = -Min(0.0, res1[i]);
      if (res1[i] < -delmin && bind[i] == 0) 
      {
        violis[0] += 1;
        violis[violis[0]] = i;
      }
    }
    // violis is the list of inequality-contraints currently       
    // not binding which have been hit during unidimensional search
    // sigres is the smallest 0.0 of secants through constraints   
    // which change sign along [x, x + d]                              

    upsi1 += term;
    if (upsi1 > tau0 && phase != -1) 
      throw 0; //Numerical("New point rejected");

    psi1 += term * w[i];
    if (res1[i] * res[i] < 0.0 && sig <= 1.0 && (bind[i] == 0 || (bind[i] == 1 && 
        (Abs(res[i]) / gresn[i] >= tp3 * epsmac 
        || Abs(res1[i]) / gresn[i] >= tp3 * epsmac))))
      sigres = Min(sigres, sig * res[i] / (res[i] - res1[i]));
  }
  if (phase != -1) 
    fx1 = mf_ObjectFunction(x1, problem);
  else 
    fx1 = 0.0;
  phi1 = scf * fx1 + psi1;
  
  return sigres;
CATCH(QSpellucci::Spellucci::mf_Evaluate)
}

void QSpellucci::Spellucci::mf_StepSize(double sig1th,
                            const WrapperOptimizationProblem &problem) const
{
TRY
  /*
    sig1th the first proposed stepsize for searching on the arc
    if sig = one did'nt work

    n = number of variables
    x = current point
    d = direction of descent. dd = second order correction
    x, d = input
    x0, d0 etc. information from previous step

    xnorm, dnorm = euclidean length of x and d
    stptrm = 1 on success , = -1 or = -2 otherwise
    sig    =  computed stepsize
    it is assumed that one is asymptotically optimal
    sigsm = smallest acceptable stepsize
    sigla = largest  acceptable stepsize
    alpha = smallest feasible reduction factor for stepsize
    delta = multiplier for derivative should be smaller than .25
    beta  = maximum feasible increase of x-norm for sig = 1
    theta = bound for cos(angle(current direction, previous direction))
            if overridden, stepsize larger than one is tried
  */

  int i, l, j, desc, descre, sminfe, lainfe;
  valarray<double> step(NSTEP + 1);

  step[1] = 0.5;
  step[2] = 0.25;
  for (i = 3; i <= NSTEP; i++) 
    step[i] = 0.1;

  // projection of d, rescaling and computing dirder has been done already
  l = 0;
  phi = scf * fx + psi;
  sig = sig1th;
  violis[0] = 0;
  
  L100:

  l += 1;
  if (l > NSTEP) 
  {
    stptrm = -1.0;
    sig = 0.0;
    return;
  }

  double sigres;
  // compute a new x and test for descent
  try
  {
    sigres = mf_Evaluate(sig, problem);
  }
  catch(...)
  {
    if (sig > 1.0) 
    {
      mf_RestoreCurrentMinimum();
      goto L200;
    } 
    else 
    {
      sig *= step[l];
      goto L100;
    }
  }
  
  // new function value
  if (sig > 1.0) 
  {
    if (phi1 >= phimin) 
    {
      // phi does'nt decrease further
      mf_RestoreCurrentMinimum();
      goto L200;
    } 
    else 
    {
      if (sig < stmaxl) 
      {
        mf_SaveCurrentMinimum();
        sig = Min(stmaxl, sig + sig);
        goto L100;
      } 
      else 
        goto L200;
    }
  }
  double diff;
  if (lastch >= itstep-3 || phase != 2 || singul) 
  {
    // require monotonic behaviour
    diff = phi - phi1;
  } 
  else 
  {
    double maxphi = phi;
    for (j = 1; j <= 3; j++) 
      maxphi = Max(scf * accinf(itstep - j, 2) + accinf(itstep - j, 4), maxphi);
    diff = maxphi - phi1;
  }
  desc = diff >= Min(-sig * delta * dirder, level);
  descre = upsi - upsi1 >= sig * pow(delta, 2) * upsi / tauqp;
  sminfe = upsi <= tau0 * 0.5 && upsi1 <= tau0;
  lainfe = upsi > tau0 * 0.5;
  
  if (desc && (sminfe || (lainfe && descre))) 
  {
    // Goldstein-Armijo descent test satisfied
    if (sig == 1.0 && ((cosphi >= theta && sig0 >= 1.0 
         && (phase + 1) * (phase - 2) != 0
         && !singul) || diff >= -sig * delta1 * dirder)
         && stmaxl > 1.0 && upsi < tau0 * 0.5) {
         
        // 1 >= delta1 >> delta > 0
        // Try stepsize larger than 1.0. Save the current point as the best 1.0
        mf_SaveCurrentMinimum();
        sig = Min(stmaxl, sig + sig);
        goto L100;
    }
    if (sig <= 1.0 && upsi > tau0*0.5 && upsi1 > upsi) 
      goto L300;
    
    goto L200;
  } 
  else 
    goto L300;

  L200:

  // accept new x, save old values 
  fx0 = fx;
  fx = fx1;
  upsi0  = upsi;
  upsi = upsi1;
  psi0 = psi;
  psi = psi1;
  stptrm = 1.0;
  sig0 = sig;
  
  x0 = x;
  d0 = d;
  x = x1;
  difx = x;
  difx -= x0;
  d0norm = dnorm;
  x0norm = xnorm;
  res = res1;

  return;
  
  // continue reducing sig
  L300:

  if (sigres < sig) 
    sig = Min(0.5 * sig, Max(step[l] * sig, sigres));
  else 
  {
    double term = (diff - dirder * sig) * 2.0;
    if (term > epsmac * (scf * Abs(fx) + psi)) 
      sig = Min(0.5 * sig, Max(step[l] * sig, -dirder * pow(sig, 2) / term));
    else 
      sig *= step[l];
  }
  
  if (sig * Max(1.0, dnorm) >= sigsm) 
    goto L100;
  stptrm = -1.0;
  sig = 0.0;
  
  return;
CATCH(QSpellucci::Spellucci::mf_StepSize)
}

void QSpellucci::Spellucci::mf_GradientLagrangian(valarray<double> &gphi) const
{
TRY
  int i;
  for (i = 1; i <= n; ++i) 
  {
    gphi[i] = gradf[i] * scf;
    int j;
    for (j = 1; j <= nh; ++j) 
      gphi[i] -= u[j] * gres(i, j);
    for (j = nh + 1; j <= alist[0]; ++j) 
    {
      int l = alist[j];
      if (u[l] > 0.0) 
        gphi[i] -= gres(i, l) * u[l];
     // include constraints, whose multipliers are of correct sign only
    }
  }
  return;
CATCH(QSpellucci::Spellucci::mf_GradientLagrangian)
}

void QSpellucci::Spellucci::mf_QRDecompBindingConstraints(int nlow,
                                                                       int nrl) const
{
TRY
  int     n1, n2, i, j, k, l, i1, i2, ipiv;

  valarray<double> qri(NX + 1);
  valarray<double> qri0(NX + 1);

  if (nlow > nrl) 
    return;
  if (nlow == 1) 
    rank = 0;
  double dl = 1.0 / (n + n + n);
  for (i = nlow; i <= nrl; i++) 
  {
    diag[i] = 0.0;
    betaq[i] = 0.0;
    colno[i] = i;
    for (j = 1; j <= n; j++) 
      qri[j] = gres(j, alist[i]);
    double sum = mf_SolveCholeskyFactorLeft(a, qri, qri0, n);
    
    if (sum == 0.0)
    {
      cscal[i] = 1.0;
      colle[i] = 0.0;
      for (j = 1; j <= n; j++)
        qr(j, i) = 0.0;
    } 
    else 
    {
      for (j = 1; j <= n; j++)
        qri[j] = qri0[perm[j]];
      double term = 1.0 / sqrt(Max(sum, pow(rho, 2)));
      cscal[i] = term;
      if (nlow > 1) 
      {
        mf_HouseholderTransformation(1, 0, 1, rank, n, qr, betaq, qri, qri0);
        qri = qri0;
      }
      for (j = 1; j <= n; j++)
        qr(j, i) = qri[j] * term;
      // colle : length of remaining column squared
      colle[i] = pow(TwoNorm(qri, rank + 1, n + 1) * term, 2);
    }
  }
  if (nlow > 1 && rank < nlow-1) 
  {
    // shift 0.0 block to the right
    i1 = nlow - 1 - rank;
    i2 = nrl - nlow + 1;
    for (i = 1; i <= Min(i1, i2); i++)
    {
      ipiv = rank + i;
      k = nrl - i + 1;
      double term = betaq[k];
      betaq[k] = betaq[ipiv];
      betaq[ipiv] = term;
      j = colno[k];
      colno[k] = colno[ipiv];
      colno[ipiv] = j;
      term = colle[k];
      colle[k] = colle[ipiv];
      colle[ipiv] = term;
      for (j = 1; j <= n; j++) 
      {
        term = qr(j, k);
        qr(j, k) = qr(j, ipiv);
        qr(j, ipiv) = term;
      }
    }
  }
  if (nlow > 1) 
  {
    n1 = rank + 1;
    n2 = n1 + nrl - nlow;
  } 
  else 
  {
    n1 = nlow;
    n2 = nrl;
  }
  for (i = n1; i <= n2; i++) 
  {
    // search for pivot column
    ipiv = i;
    double curle = colle[i];
    for (j = i + 1; j <= n2; j++) 
    {
      if (colle[j] > curle)
        curle = colle[j];
    }
    for (j = n2; j >= i; j--) 
    {
      if (colle[j] >= curle / 3.0) 
        ipiv = j;
    }
    // interchange columns explicitly, if necessary make interchanges continuous with respect to x
    if (ipiv != i) 
    {
      j = colno[i];
      colno[i] = colno[ipiv];
      colno[ipiv] = j;
      double term = colle[i];
      colle[i] = colle[ipiv];
      colle[ipiv] = term;
      for (k = 1; k <= n; k++) 
      {
        term = qr(k, i);
        qr(k, i) = qr(k, ipiv);
        qr(k, ipiv) = term;
      }
    }
    double sum = 0.0;
    for (j = i;  j <= n;  j++) 
    {
      double term = qr(j, i);
      qri[j] = term;
      sum += term*term;
    }
    if (sum <= pow(rho, 2))
    {
      // set tiny values to 0.0
      for (j = i; j <= n2; j++) 
      {
        colle[j] = 0.0;
        for (k = i; k <= n; k++) 
          qr(k, j) = 0.0;
      }
      rank = i - 1;
      
      return;
    }
    double qrii = qri[i];
    double dalpha = -sqrt(sum);
    if (Abs(qrii) <= -dalpha * dl) 
    {
      double term = 0.0;
      for (j = i + 1; j <= n; j++) 
      {
        if (Abs(qri[j]) > term) 
        {
          term = Abs(qri[j]);
          l = j;
        }
      }
      k = perm1[i];
      perm1[i] = perm1[l];
      perm1[l] = k;
    }
    if (qrii < 0.0) 
      dalpha = -dalpha;
    double dbeta = 1.0 / (sum - qrii * dalpha);
    diag[i] = dalpha;
    betaq[i] = dbeta;
    qri[i] = qrii-dalpha;
    qr(i, i) = qri[i];
    rank = i;
    i1 = i + 1;
    for (j = i1; j <= n2; j++) 
    {
      sum = dbeta * Dot(qr.Column(j), qri, i, n + 1); 
      for (k = i; k <= n; k++)
        qr(k, j) -= sum * qri[k];
      colle[j] -= pow(qr(i, j), 2);
    }
  }
  return;
CATCH(QSpellucci::Spellucci::mf_QRDecompBindingConstraints)
}

void QSpellucci::Spellucci::mf_HouseholderTransformation(int id,
                                                                      int incr,
                                                                      int is1,
                                                                      int is2,
                                                                      int m,
                                                                      Matrix<double> &a,
                                                                      valarray<double> &beta,
                                                                      valarray<double> &b,
                                                                      valarray<double> &c) const
{
TRY
  std::copy(&b[1], &b[m + 1], &c[1]);
  
  if (is1 > m) 
    return;
  if (is2 < is1) 
    return;
  
  int i, j, k, it;
  double  sum = 0.0;
  for (i = is1; i <= is2; ++i) 
  {
    it = i;
    if (id < 0) 
      it = is2 - it + is1;
    
      // it = index of transformation
    
    j = it + incr;
    sum = beta[it] * Dot(a.Column(it), c, j, m + 1); 
    for (k = j; k <= m; ++k) 
      c[k] -= sum * a(k, it);
  }
CATCH(QSpellucci::Spellucci::mf_HouseholderTransformation)
}

void QSpellucci::Spellucci::mf_SolveTriangularSystem(int nlow,
                                                                  int nup,
                                                                  valarray<double> &b,
                                                                  valarray<double> &x) const
{
TRY
  valarray<double> xl(0.0, NX + 1);

  int i;
  for (i = nup; i >= nlow; --i) 
  {
    int j;
    double sum = 0.0;
    for (j = i + 1; j <= nup; ++j) 
      sum += qr(i, j) * xl[j];
    xl[i] = (b[i] - sum) / diag[i];
  }
  for (i = nlow; i <= nup; ++i) 
    x[i] = xl[i] * cscal[colno[i]];
  // there must follow interchange of x as given by colno
  // e.g. xx(colno[i]) = x[i]                            
CATCH(QSpellucci::Spellucci::mf_SolveTriangularSystem)
}

void QSpellucci::Spellucci::mf_SolveTransposedTriangularSystem(int nlow,
                                                                      int nup,
                                                                      valarray<double> &b,
                                                                      valarray<double> &x) const
{
TRY
  int i;
  for (i = nlow; i <= nup; ++i) 
    // b has been permuted already!
    x[i] = b[i] * cscal[colno[i]];
  for (i = nlow; i <= nup; ++i) 
  {
    int j;
    double sum = 0.0;
    for (j = nlow; j <= i-1; ++j) 
      sum += qr(j, i) * x[j];
    x[i] = (x[i] - sum) / diag[i];
  }
CATCH(QSpellucci::Spellucci::mf_SolveTransposedTriangularSystem)
}

double QSpellucci::Spellucci::norm(double a, double b) const
{
  valarray<double> v(2);
  v[0] = a;
  v[1] = b;
  
  return TwoNorm(v, 0, 2);
}

void QSpellucci::Spellucci::mf_UpperTriangularCholeskyFactor(SquareMatrix<double> &r,
                                                                      valarray<double> &z,
                                                                      valarray<double> &y) const
{
TRY
  ASSERT0(z.size() == y.size(), 
    Argument("Dimension of update vectors does not match"));
  ASSERT0(a.Size() == y.size(), 
    Argument("Dimension of matrix and vectors does not match"));

  valarray<double> w(NX + 1);
  valarray<double> rn1(NX + 1);
  
  bool fail = false;
  
  int i;
  // save subdiagonal
  valarray<double> sdiag(NX + 1);
  for (i = 1; i <= n - 1; ++i)
  {
    sdiag[i] = r(i + 1, i);
    r(i + 1, i) = 0.0;
  }

  // step 1: include z*z^t
  double zl = 0.0;
  for (i = 1; i <= n; ++i)
    zl += pow(z[i], 2);
  
  if (zl != 0.0) 
  {
    // solve w^t * r = z^t
    double wl = mf_SolveCholeskyFactorLeft(r, z, w, n);
    
    wl = sqrt(wl + 1.0);
    
    // u[2]*u[3]*...*u[n]*w = (norm(w), 0,.., 0)(transpose)
    // u[i] rotations
    
    for (i = n; i >= 2; i--) 
    {
      if (w[i] != 0.0) 
      {
        int i1 = i - 1;
        double ai = w[i1];
        double bi = w[i];
        w[i1] = norm(ai, bi);
        ai /= w[i1];
        bi /= -w[i1];
        r(i, i1) = bi * r(i1, i1);
        r(i1, i1) *= ai;
        double h;
        int j;
        for (j = i; j <= n; j++) 
        {
          h = ai * r(i1, j) - bi * r(i, j);
          r(i, j) = bi * r(i1, j) + ai * r(i, j);
          r(i1, j) = h;
        }
      }
    }
    // r = d*r, d = diag(wl, 1,..., 1), r now Hessenberg
    
    for (i = 1; i <= n; i++) 
      r(1, i) *= wl;

    // r = u[n-1]*...*u[1]*r now upper triangular, u[i] givens-rotations                          
    for (i = 1; i <= n-1; i++) 
    {
      int i1 = i + 1;
      double ai = r(i, i);
      double bi = -r(i1, i);
      double h = norm(ai, bi);
      if (h != 0.0) 
      {
        ai /= h;
        bi /= h;
        r(i, i) = h;
        r(i1, i) = 0.0;
        int j;
        double h;
        for (j = i + 1; j <= n; j++) 
        {
          h = ai * r(i, j) - bi * r(i1, j);
          r(i1, j) = bi * r(i, j) + ai * r(i1, j);
          r(i, j)  = h;
        }
      }
    }
  }
  // step 2: include -y*y^t
  
  double yl = 0.0;
  for (i = 1; i <= n; i++) 
    yl += pow(y[i], 2);

  if (yl != 0.0) 
  {
    double wl = mf_SolveCholeskyFactorLeft(r, y, w, n);
    
    if (wl >= 1.0) 
      fail = true;
    else 
    {
      wl = sqrt(Abs(1.0 - wl));
      double wn1 = wl;
      
      //   (r(new), 0)                  (r  , w)
      //   (-----------)  =  u[1]*...u[n]*(-----------)
      //   (y(transp), 1)                  ((0,.., 0), wl)
      
      for (i = n; i >= 1; i--) 
      {
        double ai = wn1;
        double bi = w[i];
        wn1 = norm(ai, bi);
        if (wn1 != 0.0) 
        {
          ai /= wn1;
          bi /= wn1;
          rn1[i] = bi * r(i, i);
          r(i, i) *= ai;
          int j;
          for (j = i + 1; j <= n; j++) 
          {
            double h = ai * r(i, j) - bi * rn1[j];
            rn1[j] = bi * r(i, j) + ai * rn1[j];
            r(i, j) = h;
          }
        }
      }
    }
  }
  
  // restore subdiagonal 
  for (i = 1;  i <= n - 1; i++) 
    r(i + 1, i) = sdiag[i];

  ASSERT0(fail == false, Numerical("Decomposition does not exist"));

CATCH(QSpellucci::Spellucci::mf_UpperTriangularCholeskyFactor)
}

double QSpellucci::Spellucci::mf_SolveCholeskyFactorRight(SquareMatrix<double> &a,
                                                                         valarray<double> &b,
                                                                         valarray<double> &y,
                                                                         int n) const
{
TRY
  ASSERT0(b.size() == y.size(), 
    Argument("Dimension of vectors does not match"));
  ASSERT0(a.Size() == b.size(), 
    Argument("Dimension of matrix and vector does not match"));
  
  int i;
  double  h;
  double yl = 0.0;
  for (i = n; i >= 1; --i) 
  {
    h = b[i];
    int j;
    for (j = i + 1; j <= n; ++j) 
      h -= a(i, j) * y[j];
    h /= a(i, i);
    y[i] = h;
    yl += pow(h, 2);
  }
  return yl;
CATCH(QSpellucci::Spellucci::mf_SolveCholeskyFactorRight)
}

double QSpellucci::Spellucci::mf_SolveCholeskyFactorLeft(SquareMatrix<double> &a,
                                                                        valarray<double> &b,
                                                                        valarray<double> &y,
                                                                        int n) const
{
TRY
  ASSERT0(b.size() == y.size(), 
    Argument("Dimension of vectors does not match"));
  ASSERT0(a.NColumn() == b.size(), 
    Argument("Dimension of matrix and vector does not match"));

  int i;
  double  h;
  double yl = 0.0;

  for (i = 1; i <= n; ++i) 
  {
    h = b[i];
    int j;
    for (j = 1; j <= i-1; ++j) 
      h -= a(j, i) * y[j];
    h /= a(i, i);
    y[i] = h;
    yl += pow(h, 2);
  }
  return yl;
CATCH(QSpellucci::Spellucci::mf_SolveCholeskyFactorLeft)
}
    
//  scf * gradf(x) * d + (1 / 2) * d * a * d + summe(tauqp * d[i] + (my / 2) * pow(d[i], 2))
//  minimal subject to
//  d[i] >= 0, i = 1,..., nr
//  (gres[.][j] * d + res[j]) + vz * d[j] = 0, j = 1, ..., nh vz = -sign(res[j])
//  (gres[.][alist[j]]) * d + res[alist[j]]) + d[j] >= 0, j = nh + 1,...., nr
//  the weight tauqp is adapted during solution
//  the quasi-Newton-matrix a is taken from o8comm.h
//  a is regularized if not sufficiently well conditioned
//  the resulting d[1 + nr],..., d[n + nr] is a direction of descent for
//  the Zangwill function of the corresponding nonlinear
//  optimization problem
//  f(x) = Min, res[j] = 0, j = 1,..nh, res[j] >= 0 , j = nh + 1, nres
//  at the currrent point x if the weight tauqp is chosen appropriately
//  the quadratic programming problem is solved using the method
//  of Goldfarb and Idnani
//  variables are stored in xd (solution) and ud (multipliers)
//  in the following order xd = (d[j], j = 1, nr; d = direct. of desc.)
//  ud = (multipliers for d[i] >= 0 , i = 1,.., nr;
//  multipliers for the equality constraints ,
//  multipliers for the general inequality constraints)

void QSpellucci::Spellucci::mf_ExtendedRecursiveQuadraticSolver(void) const
{
TRY
    const double zero = 0.e0;
    const double one = 1.e0;    
    const double two = 2.e0;
    const double three = 3.e0;
    const double onep3 = 1.3e0;
    const double p8 = .8e0;
    const double tp1 = 1.e1;
    const double tp2 = 1.e2;
    const double tp3 = 1.e3;
    const double p5 = .5e0;

    int    mi, me, i, j, k, ip, l, incr, wlow;
    double infe1, s1, s2, tiny, my, zz, ss, su, t, t1, t2, f, psid, c1, c2, 
             cdiag, term, su1, su2, condr, infiny, term1, term2, diff0;

    valarray<double> y(NDUALM + 1);
    valarray<double> g0(NDUALM + 1);
    valarray<double> xd(NDUALM + 1);
    valarray<double> z(NDUALM + 1);
    valarray<double> vr(MDUALM + 1);
    valarray<double> cii(MDUALM + 1);

    valarray<double> ci0(MDUALM + 1);
    valarray<double> cei(MDUALM + 1);
    valarray<double> s(MDUALM + 1);

    valarray<double> xdold(NDUALM + 1);
    valarray<double> udold(MDUALM + 1);

    valarray<int> ai(MDUALM + 1);
    valarray<int> iai(MDUALM + 1);
    valarray<int> iaexcl(MDUALM + 1);
    valarray<int> aiold(MDUALM + 1);

    valarray<double> mult(NRESM + 1);

    valarray<int> qpdel(MDUALM + 1);
    
    ndual = n + nr;
    
    // number of equality constraints in QP-problem 
    
    me = nh;
    
    /* QP inequality constraints = active constraints - */
    /* equality-constraints + slack's                   */
    
    mi     = 2*nr-nh;
    infiny = epsmac / tolmac;
    if (analyt) 
    {
        tiny = 2*nr*epsmac*tp3;
    } 
    else 
    {
        tiny = 2*nr*Max(epsdif, epsmac*tp3);
    }
    qpterm = 0;
    for (i = 1; i <= nr; i++) {
    
        /* check gradients of active constraints against zero */
        
        for (j = 1; j <= n; j++) {
            y[j] = gres(j, alist[i]);
        }
        mult[i] = one;
        if (TwoNorm(y, 1, n + 1) == zero) {
            mult[i] = zero;
        }
    }
    /* restart point in case of increase of tauqp */
    
    L10:

    /* initialize matrices j and r */
    
    for (i = 1; i <= ndual; i++) {
        ddual[i] = zero;
        for (j = 1; j <= ndual; j++) {
            r(j, i)  = zero;
            xj(j, i) = zero;
        }
    }
    rnorm = one;
    rlow  = one;
    term1 = zero;
    for (i = 1; i <= nres; i++) {
        u[i] = zero;
        if (w[i] > term1) term1 = w[i];
    }
    accinf(itstep, 19) = clow;
    accinf(itstep, 20) = term1;
    accinf(itstep, 31) = tauqp;
    for (i = 1; i <= me + mi; i++) {
        ud[i] = zero;
    }
    c1 = Abs(a(1, 1));
    for (i = 1; i <= n; i++) {
        c1 = Max(c1, Abs(a(i, i)));
    }
    c1 = c1*tp1;
    
    /* we require much more regularity of a in the singular case */
    
    for (i = 1; i <= n; i++) {
        if (Abs(a(i, i)) < sqrt(rho1)*c1) a(i, i) = sqrt(rho1)*c1;
    }
    /* invert the Cholesky-factor and store in xj (Idnanis j-matrix) */
    
    mf_InverseUpperTriangularMatrix(n, a, ndual, xj);
    
    c1   = Abs(a(1, 1));
    incr = nr;
    c2   = Abs(xj(1 + incr, 1 + incr));
    for (i = 1; i <= n; i++) {
        c1 = Max(c1, Abs(a(i, i)));
        c2 = Max(c2, Abs(xj(i + incr, i + incr)));
    }
    my = zero;
    for (i = 1; i <= n; i++) {
        for (j = i; j <= n; j++) {
             my = my + pow(a(i, j), 2);
        }
    }
    my    = my / n;
    cdiag = one / sqrt(my);
    for (i = 1; i <= incr; i++) {
        xj(i, i) = cdiag;
    }
    for (i = 1; i <= ndual; i++) {
        if (i >  incr) g0[i] = gradf[i-incr]*scf;
        if (i <= incr) g0[i] = tauqp;
    }
    /* compute unconstrained solution */
    /* the Cholesky-factor of a is stored in the upper triangle */
    
    for (i = 1; i <= n; i++) {
        su = zero;
        for (j = 1; j <= i-1; j++) {
            su = su + a(j, i)*y[j + incr];
        }
        y[i + incr] = (g0[i + incr] - su) / a(i, i);
    }
    for (i = n; i >= 1; i--) {
        su = zero;
        for (j = i + 1; j <= n; j++) {
            su = su + a(i, j)*y[j + incr];
        }
        y[i + incr] = (y[i + incr] - su) / a(i, i);
    }
    for (i = 1; i <= incr; i++) {
    
        /* initially assume the slacks being zero */
        
        y[i] = zero;
    }
    for (i = 1; i <= ndual; i++) {
        xd[i] = -y[i];
    }
    /* unconstrained minimizer of the QP: slacks come first */
    
    f = p5*Dot(g0, xd, 1, ndual + 1);
    
    /* define the initial working set: all slacks are at their */
    /* lower bounds                                            */
    
    iq = nr;
    for (i = 1; i <= iq; i++) {
        ai[i]   = i;
        r(i, i) = one;
        ud[i]   = tauqp;
    }
    rnorm = one;
    rlow  = one;
    
    /* introduction of equality constraints */
    
    for (i = 1; i <= me; i++) {
        for (j = 1; j <= iq; j ++) {
            ud1[j] = ud[j];
        }
        L20:

        ud1[iq + 1] = zero;
        
        /* an equality constraint is indicated by the negative index */
    
        ai[iq + 1] = -nr-i;
        for (j = 1; j <= n; j++) {
        
            /* me = nh and alist[i] = i for i = 1,.., nh */
            
            cei[j + incr] = gres(j, i);
        }
        for (j = 1; j <= incr; j++) {
            cei[j] = zero;
        }
        cei[i] = one;
        if (res[i] > zero) cei[i] = -one;
        for (j = 1; j <= ndual; j++) {
            np[j] = cei[j];
        }
        mf_UpdateProjGradient(z);
        
        if (iq != 0) mf_CorrDualMultipliers(vr);
        
        /* |z| = 0? */
        
        /* old      zz = o8sc1(1, ndual, z, z)     */
            
        zz   = TwoNorm(z, 1, ndual + 1);
        term = Dot(z, np, 1, ndual + 1);
        
        /* old      if (zz != zero && term > zero)  */
        
        if (zz >= tiny*rnorm && term > zero) {
            t2 = (-Dot(np, xd, 1, ndual + 1) - res[i]) / term;
        } else if ((-Dot(np, xd, 1, ndual + 1)-res[i]) >= zero) {
            t2 = infiny;
        } else {
            t2 = -infiny;
        }
        /* in addition of an equality constraint, t2 may be positive or negative*/
        
        if (iq != 0) mf_CorrDualMultipliers(vr);
        l = 0;
        if (t2 > zero) {
            t1 = infiny;
            for (k = 1; k <= iq; k++) {
                if (vr[k] > zero && ai[k] > 0) {
                    if (ud1[k]/vr[k] < t1) {
                        t1 = ud1[k]/vr[k];
                    }
                }
            }
            t = Min(t1, t2);
        } else {
            t1 = infiny;
            for (k = 1; k <= iq; k++) {
                if (vr[k] < zero && ai[k] > 0) {
                    if (ud1[k]/Abs(vr[k]) < t1) {
                        t1 = ud1[k]/Abs(vr[k]);
                    }
                }
            }
            t1 = -t1;
            t  = Max(t1, t2);
            
            /* t now negative */
        }
        /* add constraint , otherwise we must first delete some */
        /* inequality constraint with negative multiplier       */
        /* first delete then add!                               */
        
        if (Abs(t)  >= infiny) goto L2000;
        
        if (Abs(t2) >= infiny) {
        
            /* purely dual step */
            
            for (k = 1; k <= iq; k++) {
                ud1[k] = ud1[k]+t*(-vr[k]);
                if (ud1[k] < zero && ai[k] > 0) ud1[k] = zero;
            }
            ud1[iq + 1] = ud1[iq + 1]+t;
            qpdel[0]  = 0;
            for (j = 1; j <= iq; j++) {
                if (ud1[j] <= tiny && ai[j] > 0) {
                    qpdel[0]        = qpdel[0]+1;
                    qpdel[qpdel[0]] = ai[j];
                }
            }
            for (k = 1; k <= qpdel[0]; k++) {
                l      = qpdel[k];
                iai[l] = l;
                
                mf_DeleteConstraint(ai, l);
            }
            goto L20;
        }
        for (k = 1; k <= ndual; k++) {
            xd[k] = xd[k]+t*z[k];
        }
        for (k = 1; k <= iq; k++) {
            ud1[k] = ud1[k]+t*(-vr[k]);
            if (ud1[k] < zero && ai[k] > 0) ud1[k] = zero;
        }
        ud1[iq + 1] = ud1[iq + 1] + t;
        
        f = f + t * Dot(z, np, 1, ndual + 1) * (p5 * t + ud1[iq + 1]);
        
        if (Abs(t2-t1) <= tiny) {
            qpdel[0] = 0;
            for (j = 1; j <= iq; j++) {
                if (ud1[j] <= tiny && ai[j] > 0) {
                    qpdel[0]        = qpdel[0]+1;
                    qpdel[qpdel[0]] = ai[j];
                }
            }
            for (k = 1; k <= qpdel[0]; k++) {
                l      = qpdel[k];
                iai[l] = l;
                
                mf_DeleteConstraint(ai, l);
            }
            ai[iq + 1] = -i-nr;
            
            mf_AddConstraint();
            
        } else if (t == t2) {
            ai[iq + 1] = -i-nr;
            
            mf_AddConstraint();
            
        } else {
            qpdel[0] = 0;
            for (j = 1; j <= iq; j++) {
                if (ud1[j] <= tiny && ai[j] > 0) {
                    qpdel[0]        = qpdel[0]+1;
                    qpdel[qpdel[0]] = ai[j];
                }
            }
            for (k = 1; k <= qpdel[0]; k++) {
                l      = qpdel[k];
                iai[l] = l;
                
                mf_DeleteConstraint(ai, l);
            }
            goto L20;
        }
        for (j = 1; j <= iq; j++) {
            ud[j] = ud1[j];
        }  
    }
    /* set iai = k\ai */
    
    for (i = 1; i <= mi; i++) {
        iai[i] = i;
    }
    /* step 1 */
    
    L50:

    /* ai = QP - working set , iai[i] = 0 if i in ai */
    
    for (i = 1; i <= iq; i++) {
        ip = ai[i];
        if (ip > 0) iai[ip] = 0;
    }
    /* s[xd] = ci(trans)*xd + ci0 >= 0 ? */
    
    psid = zero;
    
    /* psid : the measure of infeasibility */
    
    for (i = 1; i <= mi; i++) {
    
        /* iaexcl: if = 0, exclude from addition in this cycle */
        
        iaexcl[i] = 1;
        su        = zero;
        
        // numbers of inequality constraints:
        // i = 1,..., nr corresponds to the constraints v >= 0, u_a >= 0
        // i = nr + 1,...., mi to the regularized general inequalities
        
        if (i > nr) {
            k = alist[i + nh - incr];
            for (j = 1; j <= n; j++) {
                cii[j + incr] = gres(j, k);
            }
            for (j = 1; j <= incr; j++) {
                cii[j] = zero;
            }
            cii[nh + i - incr] = one;
            ci0[i]         = res[k];
        } else {
            for (j = 1; j <= ndual; j++) {
                cii[j] = zero;
            }
            ci0[i] = zero;
            cii[i] = one;
        }
        su   = Dot(cii, xd, 1, ndual + 1) + ci0[i];
        s[i] = su;
        psid = psid + Min(zero, su);
    }
    for (i = 1; i <= iq; i++) {
        udold[i] = ud[i];
        aiold[i] = ai[i];
    }
    for (i = 1; i <= ndual; i++) {
        xdold[i] = xd[i];
    }
    L60:

    ss = zero;
    ip = 0;
    
    /* introduce most violated inequality constraint */
    
    for (i = 1; i <= mi; i++) {
        if (s[i] < ss && iai[i] != 0 && iaexcl[i] != 0) {
            ss = s[i];
            ip = i;
        }
    }
    if (iq > 1) {
        condr = rnorm / rlow;
    } else {
        condr = one;
    }
    if (Abs(psid) <= tiny * (c1 * c2 + condr) || ip == 0) 
    {
    
        /* successful termination of QP-solver for current tauqp */
        
        qpterm = 1;
        accinf(itstep, 30) = one;
        accinf(itstep, 13) = condr;
        accinf(itstep, 14) = c1*c2;
        for (i = 1; i <= n; i++) {
            d[i] = xd[i + incr];
        }
        /* new : dnorm added */
        
        dnorm  = TwoNorm(d, 1, n + 1);
        infeas = zero;
        
        for (i = 1; i <= incr; i++) {
            infeas = infeas + Abs(xd[i]);
        }
        /* L1-norm of slack variables */
        
        accinf(itstep, 31) = tauqp;
        accinf(itstep, 32) = infeas;
        wlow = false;
        su1  = zero;
        su2  = zero;
        for (i = 1; i <= iq; i++) {
            if (ai[i] < 0) {
                u[-(ai[i]+nr)] = ud[i];
            } else {
                if (ai[i] > nr) u[alist[ai[i]+nh-nr]] = ud[i];
            }
        }
        term1 = zero;
        for (j = 1; j <= n; j++) {
            np[j] = gradf[j]*scf;
        }
        for (i = 1; i <= nres; i++) {
            for (j = 1; j <= n; j++) {
                np[j] = np[j]-gres(j, i)*u[i];
            }
        }
        b2n = TwoNorm(np, 1, n + 1);
        
        if (scf != zero) b2n = b2n / scf;
        
        /* correction in the original variables */
            
        infe1 = zero;
        for (i = 1; i <= nr; i++) {
            infe1 = infe1 + Abs(xd[i]) * mult[i];
        }
        if (upsi <= delmin*nres && upsi0 <= delmin*nres
            && b2n <= (gfn + one) * epsx * tp2 && phase >= 0
            && infeas <= delmin*nres) {
         
            /* since multipliers may be incorrect for infeas != zero be careful */
            /* we consider the problem as successfully solved with reduced      */
            /* requirements                                                     */

            for (i = 1; i <= n; i++) {
                d[i] = zero;
            }
            dnorm  = zero;
            optite = three;
            
            return;
        }
        /* there may be an additional increase of tauqp necessary again */

        if (infe1 >= (one - delta1 / tauqp) * upsi && 
            (TwoNorm(d, 1, n + 1) <= Min(infe1, pow(infe1, 2)) * tp1
            || upsi > tau0*p5)) {
            
            /* further increase tauqp ! */
            
            for (i = 1; i <= nres; i++) {
                u[i]     = zero;
                slack[i] = zero;
            }
            if (tauqp*taufac > taumax) {
                qpterm = -1;
                accinf(itstep, 30) = qpterm;
                accinf(itstep, 31) = tauqp;
                accinf(itstep, 32) = infeas;
                for (i = 1; i <= n; i++) {
                    d[i] = zero;
                }
                dnorm = zero;
                
                return;
                
            } else {
                tauqp = tauqp*taufac;
                
                goto L10;
            }
        }
        /* compute new weights for the penalty-function */
        
//        L500:

        for (i = 1; i <= nres; i++) {
            slack[i] = zero;
        }
        for (i = 1; i <= nr; i++) {
            slack[alist[i]] = xd[i];
        }
        wlow = false;
        for (i = 1; i <= nres; i++) {
            w1[i] = w[i];
            if (i <= nh) {
                if (Abs(slack[i]) > Abs(res[i])+tiny) {
                    w1[i] = Abs(u[i]);
                } else {
                    w1[i] = ny * Abs(u[i]) + tau;
                }
            } else {
                if (bind[i] == 0) {
                    w1[i] = Max(w[i] * p8, tau);
                } else {
                    if (res[i] >= zero && slack[i] <= tiny)
                      w1[i] = Max(ny * Abs(u[i]) + tau, 
                                  (Abs(u[i]) + w1[i]) * p5);
                    if (res[i] >= zero && slack[i] > tiny)
                      w1[i] = Abs(u[i]);
                    if (res[i] < zero && slack[i] <= -res[i] + tiny)
                      w1[i] = Max(ny * Abs(u[i]) + tau, 
                                  (Abs(u[i]) + w1[i]) * p5);
                    if (res[i] < zero && slack[i] > -res[i] + tiny)
                      w1[i] = Abs(u[i]);
                }
            } 
            if (w1[i] < w[i]) wlow = true;
        }
        if (wlow) {
            s1 = zero;
            s2 = zero;
            for (i = 1; i <= nres; i++) {
                if (i <= nh) {
                    s1 = s1 + w1[i] * Abs(resst[i]);
                    s2 = s2 + w1[i] * Abs(res[i]);
                } else {
                    s1 = s1 - Min(zero, resst[i]) * w1[i];
                    s2 = s2 - Min(zero, res[i]) * w1[i];
                }
            }
            diff0 = (fxst - fx) * scf + (s1 - s2);
            if (diff0 >= eta * clow && itstep - lastdw >= Max(5, Min(n / 10, 20))) {
            
                /* accept new (diminished) weights */
            
                lastdw = itstep;
                lastch = itstep;
                level  = diff0 / iterma;
                psist  = s1;
                psi    = s2;
                for (i = 1; i <= nres; i++) {
                    if (w1[i] != w[i]) lastch = itstep;
                    w[i] = w1[i];
                }
                clow = clow + one;
                if (clow > itstep / 10) {
                
                    /* additional increase of eta */
                    
                    eta = eta*onep3;
                }
                
                goto L1000;
            }
        }
        /* we cannot accept new weights */
        /* reset weights                */
        
        for (i = 1; i <= nres; i++) 
        {
          w1[i] = w[i];
          if (i <= nh) 
          {
            if (slack[i] > Abs(res[i])) 
              w1[i] = Abs(u[i]);
            if (slack[i] <= Abs(res[i])) 
            {
              if (w[i] <= Abs(u[i]) 
                  && Abs(u[i]) <= w[i] + tau)
              {
                w1[i] = w[i] + two * tau;
              } 
              else 
              {
                w1[i] = Max(w[i], ny*Abs(u[i]) + tau);
              }
            }
          } 
          else 
          {
            if (slack[i] > -Min(-tiny, res[i]) && bind[i] == 1) 
            {
              w1[i] = Abs(u[i]);
            } 
            else if (bind[i] == 1 && slack[i] <= -Min(-tiny, res[i]) 
                      && u[i] <= w[i] + tau && w[i] >= u[i])
            {
              w1[i] = w[i] + two * tau;
            } 
            else if (bind[i] == 1) 
            {
              w1[i] = Max(w[i], ny * Abs(u[i]) + tau);
            }
          }
        }
        term1 = zero;
        for (i = 1; i <= nres; i++) {
            if (w1[i] > w[i] || w1[i] < w[i]) lastch = itstep;
            if (w1[i] > w[i]) lastup = itstep;
            if (w1[i] < w[i]) lastdw = itstep;
            w[i]  = w1[i];
            term1 = Max(term1, w[i]);
        }
        s1 = zero;
        s2 = zero;
        for (i = 1; i <= nres; i++) {
            if (i <= nh) {
                s1 = s1 + w[i]*Abs(resst[i]);
                s2 = s2 + w[i]*Abs(res[i]);
            } else {
                s1 = s1 - w[i]*Min(zero, resst[i]);
                s2 = s2 - w[i]*Min(zero, res[i]);
            }
        }
        psist = s1;
        psi   = s2;
        accinf(itstep, 20) = term1;
        accinf(itstep, 19) = clow;
        
        goto L1000;
    }
    if (ip > nr) {
        k = alist[ip + nh - nr];
        for (j = 1; j <= n; j++) {
            cii[j + incr] = gres(j, k);
        }
        for (j = 1; j <= incr; j++) {
            cii[j] = zero;
        }
        cii[nh + ip - nr] = one;
        ci0[ip]       = res[k];
    } else {
        for (j = 1; j <= ndual; j++) {
            cii[j] = zero;
        }
        ci0[ip] = zero;
        cii[ip] = one;
    }
    for (i = 1; i <= ndual; i++) {
        np[i] = cii[i];
    }
    for (i = 1; i <= iq; i++) {
        ud1[i] = ud[i];
    }
    ud1[iq + 1] = zero;
    ai[iq + 1]  = ip;
    
    L100:

    /* step 2a */
    
    mf_UpdateProjGradient(z);
    
    if (iq != 0) mf_CorrDualMultipliers(vr);
    
    l  = 0;
    t1 = infiny;
    for (k = 1; k <= iq; k++) {
        if (ai[k] > 0 && vr[k] > zero) {
            if (ud1[k]/vr[k] < t1) {
                t1 = ud1[k]/vr[k];
            }
        }
    }
    /* |z| = 0? */
    
    /* old      zz = o8sc1(1, ndual, z, z) */
    
    zz   = TwoNorm(z, 1, ndual + 1);
    term = Dot(z, np, 1, ndual + 1);
    
    /* old      if (zz != zero && term > zero)  */
    
    if (zz >= tiny*rnorm && term > zero) {
        t2 = -s[ip]/term;
    } else {
        t2 = infiny;
    }
    t = Min(t1, t2);
    if (t >= infiny) goto L2000;

    if (t2 >= infiny) {
        for (k = 1; k <= iq; k++) {
            ud1[k] = ud1[k]+t*(-vr[k]);
            if (ud1[k] < zero && ai[k] > 0) ud1[k] = zero;
        }
        ud1[iq + 1] = ud1[iq + 1] + t;
        qpdel[0]  = 0;
        for (i = 1; i <= iq; i++) {
            if (ud1[i] <= tiny && ai[i] > 0) {
                qpdel[0]        = qpdel[0]+1;
                qpdel[qpdel[0]] = ai[i];
            }
        }
        for (k = 1; k <= qpdel[0]; k++) {
            l      = qpdel[k];
            iai[l] = l;
            
            mf_DeleteConstraint(ai, l);
        }
        goto L100;
    }
    for (k = 1; k <= ndual; k++) {
        xd[k] = xd[k]+t*z[k];
    }
    for (k = 1; k <= iq; k++) {
        ud1[k] = ud1[k]+t*(-vr[k]);
        if (ud1[k] < zero && ai[k] > 0) ud1[k] = zero;
    }
    ud1[iq + 1] = ud1[iq + 1] + t;
    
    f = f + t * Dot(z, np, 1, ndual + 1) * (p5 * t + ud1[iq + 1]);
    
    if (t2 <= t1-tiny) {
    
        /* ddual is computed by mf_UpdateProjGradient */
        
        if (TwoNorm(ddual, iq + 1, ndual + 1) < epsmac*rnorm) {
        
            /* degeneracy: adding this constraint gives a singular working set */
            /* theoretically impossible, but due to roundoff this may occur.   */
            /* mark this constraint and try to add another one                 */
            
            iptr = ip;
            iqtr = iq;
            for (i = 1; i <= iq; i++) {
                aitr[i] = ai[i];
            }
            sstr  = ss;
            riitr = TwoNorm(ddual, iq + 1, ndual + 1);
            
            iaexcl[ip] = 0;
            for (i = 1; i <= mi; i++) {
                iai[i] = i;
            }
            for (i = 1; i <= iq; i++) {
                ai[i] = aiold[i];
                if (ai[i] > 0) iai[ai[i]] = 0;
                ud1[i] = udold[i];
            }
            for (i = 1; i <= ndual; i++) {
                xd[i] = xdold[i];
            }
            goto L60;
        }
        /* add constraint, l-pair */
        
        mf_AddConstraint();
        
        iai[ip] = 0;
        for (i = 1; i <= iq; i++) {
            ud[i] = ud1[i];
        }
        goto L50;
    }
    su = zero;
    if (ip > nr) {
    
        /* a general linear inequality constraint */
        
        k = alist[ip + nh-nr];
        for (j = 1; j <= n; j++) {
            cii[j + incr] = gres(j, k);
        }
        for (j = 1; j <= incr; j++) {
            cii[j] = zero;
        }
        cii[nh + ip-nr] = one;
        ci0[ip]       = res[k];
        
        s[ip] = Dot(cii, xd, 1, ndual + 1) + ci0[ip];
    } else {
    
        /* a slack constraint */
        
        s[ip] = xd[ip];
    }
    /* now t = t1 */
    
    qpdel[0] = 0;
    for (i = 1; i <= iq; i++) {
        if (ud1[i] <= tiny && ai[i] > 0) {
            qpdel[0]        = qpdel[0]+1;
            qpdel[qpdel[0]] = ai[i];
        }
    }
    for (k = 1; k <= qpdel[0]; k++) {
        l      = qpdel[k];
        iai[l] = l;
        
        mf_DeleteConstraint(ai, l);
    }
    if (t2 <= t1 + tiny) {
        if (TwoNorm(ddual, iq + 1, ndual + 1) < epsmac*rnorm) {
        
            /* degeneracy */
            
            iptr = ip;
            iqtr = iq;
            for (i = 1; i <= iq; i++) {
                aitr[i] = ai[i];
            }
            sstr  = ss;
            riitr = TwoNorm(ddual, iq + 1, ndual + 1);
            
            iaexcl[ip] = 0;
            for (i = 1; i <= mi; i++) {
                iai[i] = i;
            }
            for (i = 1; i <= iq; i++) {
                ai[i] = aiold[i];
                if (ai[i] > 0) iai[ai[i]] = 0;
                ud1[i] = udold[i];
            }
            for (i = 1; i <= ndual; i++) {
                xd[i] = xdold[i];
            }
            goto L60;
        }
        /* add constraint, l-pair */
        
        mf_AddConstraint();
        
        iai[ip] = 0;
        for (i = 1; i <= iq; i++) {
            ud[i] = ud1[i];
        }
        goto L50;
        
    } else {
    
        goto L100;
    }
    /* this is the exit point of mf_ExtendedRecursiveQuadraticSolver        */
    /* we either may have successful or unsuccessful termination here       */
    /* the latter with qpterm = -2 or -3, in which case it may nevertheless */
    /* be possible to use the computed d. -2 or -3 exit is theoretically    */
    /* impossible but may occur due to roundoff effects.                    */
    /* we check the directional derivative of the penalty-function now      */
    
    L1000:

    /* cut and rescale d if appropriate */
    
    mf_Cut_d();
    
    /* compute the directional derivative dirder */
    
    mf_DerivativePenaltyFunction();
    
    if (dirder >= zero || (-dirder <= epsmac * tp2 * 
        (scf * Abs(fx) + psi + one) && infeas > Max(upsi, nres*delmin))) {
        if (tauqp <= taumax / taufac) {
            tauqp = tauqp*taufac;
            
            goto L10;
            
        } else {
            qpterm = -1;
            accinf(itstep, 30) = qpterm;
            accinf(itstep, 31) = tauqp;
            accinf(itstep, 32) = infeas;
            for (i = 1; i <= n; i++) {
                d[i] = zero;
            }
            dnorm = zero;
        }
    }
    return;
    
    L2000:

    /* QP infeasible (in this application impossible , theoretically) */
    
    qpterm = -2;
    accinf(itstep, 30) = -two;
    accinf(itstep, 13) = condr;
    accinf(itstep, 14) = c1*c2;
    for (i = 1; i <= n; i++) {
        d[i] = xd[i + incr];
    }
    dnorm = TwoNorm(d, 1, n + 1);
    su1   = zero;
    for (i = 1; i <= incr; i++) {
        su1 = su1 + Abs(xd[i]);
    }
    /* L1-norm of slack variables */
    
    accinf(itstep, 32) = su1;
    wlow = false;
    su1  = zero;
    su2  = zero;
    for (i = 1; i <= iq; i++) {
        if (ai[i] < 0) {
            u[-(ai[i]+nr)] = ud[i];
        } else {
            if (ai[i] > nr) u[alist[ai[i]-nr + nh]] = ud[i];
        }
    }
    term1 = zero;
    for (j = 1; j <= n; j++) {
        np[j] = gradf[j]*scf;
    }
    for (i = 1; i <= nres; i++) {
        for (j = 1; j <= n; j++) {
            np[j] = np[j]-gres(j, i)*u[i];
        }
        term  = Dot(gres.Column(i), d, 1, n + 1); 
        term2 = res[i]/Max(one, gresn[i]);
        w1[i] = w[i];
        if (i <= nh) {
            if (slack[i] > Abs(res[i])) w1[i] = Abs(u[i]);
            if (slack[i] <= Abs(res[i])) {
                if (w[i] <= Abs(u[i]) && Abs(u[i]) <= w[i]+tau) {
                    w1[i] = w[i]+two*tau;
                } else {
                    w1[i] = Max(w[i], ny * Abs(u[i]) + tau);
                }
            }
            su1 = su1 + Abs(res[i])*w1[i];
            su2 = su2 + Abs(resst[i])*w1[i];
        } else {
            if (slack[i] > -Min(-tiny, res[i]) && bind[i] == 1) {
                w1[i] = Abs(u[i]);
            } else if (bind[i] == 1 && slack[i] <= -Min(-tiny, res[i])
                        && u[i] <= w[i]+tau && w[i] >= u[i]) {
                w1[i] = w[i]+two*tau;
            } else if (bind[i] == 1) {
                w1[i] = Max(w[i], ny * Abs(u[i]) + tau);
            }
            su1 = su1-w1[i]*Min(zero, res[i]);
            su2 = su2-w1[i]*Min(zero, resst[i]);
        }
        if (w[i] != w1[i]) lastch = itstep;
        w[i]  = w1[i];
        term1 = Max(term1, w[i]);
    }
    psist = su2;
    psi   = su1;
    
    b2n   = sqrt(Dot(np, np, 1, n + 1));
    
    if (scf != zero) b2n = b2n / scf;
    if (wlow) {
        clow   = clow + one;
        lastch = itstep;
        lastdw = itstep;
    }
    accinf(itstep, 19) = clow;
    accinf(itstep, 20) = term1;
    accinf(itstep, 31) = tauqp;
    
    goto L1000;
CATCH(QSpellucci::Spellucci::mf_ExtendedRecursiveQuadraticSolver)
}

void QSpellucci::Spellucci::mf_UpdateProjGradient(valarray<double> &z) const
{
TRY
  int i;
  for (i = 1; i <= ndual; ++i) 
    ddual[i] = Dot(xj.Column(i), np, 1, ndual + 1);

  for (i = 1; i <= ndual; ++i) 
    z[i] = Dot(xj.Row(i), ddual, iq + 1, ndual + 1);

  return;
CATCH(QSpellucci::Spellucci::mf_UpdateProjGradient)
}

void QSpellucci::Spellucci::mf_CorrDualMultipliers(valarray<double> &rv) const
{
TRY
  int i;
  double s = 0.0;
  for (i = iq; i >= 1; i--) 
  {
    s = Dot(r.Row(i), rv, i + 1, iq + 1);
    rv[i] = (ddual[i] - s) / r(i, i);
  }
  return;
CATCH(QSpellucci::Spellucci::mf_CorrDualMultipliers)
}

void QSpellucci::Spellucci::mf_DeleteConstraint(valarray<int> &ai,
                                                             int l) const
{
TRY
  int i, qq;
  for (i = 1; i <= iq; ++i) 
  {
    if (ai[i] == l) 
    {
      qq = i;
      break;
    }
  }

  int j;
  for (i = qq; i <= iq - 1; ++i) 
  {
    ai[i] = ai[i + 1];
    ud1[i] = ud1[i + 1];
    for (j = 1; j <= ndual; ++j) 
      r(j, i) = r(j, i + 1);
  }

  ai[iq] = ai[iq + 1];
  ud1[iq] = ud1[iq + 1];
  ai[iq + 1] = 0;
  ud1[iq + 1] = 0.0;
  for (j = 1; j <= iq; ++j) 
    r(j, iq) = 0.0;
  iq -= 1;
  
  if (iq != 0) 
  {
    for (j = qq; j <= iq; ++j) 
    {
      double cc = r(j, j);
      double ss = r(j + 1, j);
      double h = norm(cc, ss);
    
      if (h == 0.0) continue; 
    
      double c1 = cc / h;
      double s1 = ss / h;
      r(j + 1, j) = 0.0;
      if (c1 < 0.0) 
      {
        r(j, j) = -h;
        c1 = -c1;
        s1 = -s1;
      } 
      else 
        r(j, j) = h;
      double t1, t2;
      int k;
      double xny = s1 / (1.0 + c1);
      for (k = j + 1; k <= iq; ++k) 
      {
        t1 = r(j, k);
        t2 = r(j + 1, k);
        r(j, k) = t1 * c1 + t2 * s1;
        r(j + 1, k) = xny * (t1 + r(j, k)) - t2;
      }
      for (k = 1; k <= ndual; ++k) 
      {
        t1 = xj(k, j);
        t2 = xj(k, j + 1);
        xj(k, j) = t1 * c1 + t2 * s1;
        xj(k, j + 1) = xny * (xj(k, j) + t1) - t2;
      }
    }
  }

  rnorm = 1.0;
  rlow = 1.0;
  
  if (iq >= 1) 
  {
    rnorm = Abs(r(1, 1));
    rlow = Abs(r(1, 1));
    i = 1;
    while (i < iq) 
    {
      i += 1;
      rnorm = Max(rnorm, Abs(r(i, i)));
      rlow = Min(rlow, Abs(r(i, i)));
    }
  }
CATCH(QSpellucci::Spellucci::mf_DeleteConstraint)
}

void QSpellucci::Spellucci::mf_AddConstraint(void) const
{
TRY
  int j;
  for (j = ndual; j >= iq + 2; --j) 
  {
    double cc = ddual[j - 1];
    double ss = ddual[j];
    double h = norm(cc, ss);
    
    if (h != 0.0)
    {
      ddual[j] = 0.0;
      double s1 = ss / h;
      double c1 = cc / h;
      if (c1 < 0.0) 
      {
        c1 = -c1;
        s1 = -s1;
        ddual[j - 1] = -h;
      } 
      else 
      {
        ddual[j - 1] = h;
      }
      double xny = s1 / (1.0 + c1);
      int k;
      for (k = 1; k <= ndual; k++) 
      {
        double t1 = xj(k, j - 1);
        double t2 = xj(k, j);
        xj(k, j - 1) = t1 * c1 + t2 * s1;
        xj(k, j) = xny * (t1 + xj(k, j - 1)) - t2;
      }
    }
  }
  iq = iq + 1;
  int i;
  for (i = 1; i <= iq; ++i) 
    r(i, iq) = ddual[i];

  rnorm = 1.0;
  rlow = 1.0;
  
  if (iq >= 1) 
  {
    rnorm = Abs(r(1, 1));
    rlow = Abs(r(1, 1));
    i = 1;
    while (i < iq) 
    {
      i += 1;
      rnorm = Max(rnorm, Abs(r(i, i)));
      rlow = Min(rlow, Abs(r(i, i)));
    }
  }
CATCH(QSpellucci::Spellucci::mf_AddConstraint)
}

void QSpellucci::Spellucci::mf_InverseUpperTriangularMatrix(int n,
                                                                    SquareMatrix<double> &a,
                                                                    int ndual,
                                                                    SquareMatrix<double> &x) const
{
TRY
  int incr = ndual - n;

  int l, j, k;
  for (j = n; j >= 1; j--) 
  {
    x(j + incr, j + incr) = 1.0 / a(j, j);
    for (k = j - 1; k >= 1; --k) 
    {
      double su = 0.0;
      for (l = k + 1; l <= j; ++l) 
        su += a(k, l) * x(l + incr, j + incr);
      x(k + incr, j + incr) = -su / a(k, k);
    }
  }
CATCH(QSpellucci::Spellucci::mf_InverseUpperTriangularMatrix)
}


double QSpellucci::Spellucci::mf_ObjectFunction(valarray<double> &x,
                                      const WrapperOptimizationProblem &problem) const
{
TRY

  xtr = x * xsc;
  valarray<double> point(problem.SolutionDimension());
  std::copy(&xtr[1], &xtr[problem.SolutionDimension() + 1], &point[0]);

  return problem.ObjectFunction(point);
CATCH(QSpellucci::Spellucci::mf_ObjectFunction)
}

double QSpellucci::Spellucci::mf_ZeroConstraints(int i,
                                       valarray<double> &x,
                                       const WrapperOptimizationProblem &problem) const
{
TRY

  xtr = x * xsc;
  valarray<double> point(problem.SolutionDimension());
  std::copy(&xtr[1], &xtr[problem.SolutionDimension() + 1], &point[0]);

  return problem.ZeroConstraint(i, point);
CATCH(QSpellucci::Spellucci::mf_ZeroConstraints)
}


double QSpellucci::Spellucci::mf_PositivityConstraints(int i,
                                             valarray<double> &x,
                                             const WrapperOptimizationProblem &problem) const
{
TRY
  
  xtr = x * xsc;
  valarray<double> point(problem.SolutionDimension());
  std::copy(&xtr[1], &xtr[problem.SolutionDimension() + 1], &point[0]);
  
  return problem.PositivityConstraint(i, point);
CATCH(QSpellucci::Spellucci::mf_PositivityConstraints)
}


void QSpellucci::Spellucci::mf_InitialiseWorkSpace(
                                                   const WrapperOptimizationProblem &problem) const
{
TRY
  NX = problem.SolutionDimension();
  NRESM = problem.NZeroConstraint() + problem.NPositivityConstraint();
  if (NRESM < NX) 
    NRESM = NX;

  NDUALM = NX + NRESM;
  MDUALM = NRESM * 2;

  val.resize(NRESM + 1);
  gconst.resize(NRESM + 1);
  llow.resize(NX + 1);
  lup.resize(NX + 1);
  cres.resize(NRESM + 1);
  cgres.resize(NRESM + 1);
  aitr.resize(MDUALM + 1);
  colno.resize(2*NRESM + 1);
  perm1.resize(NX + 1);
  perm.resize(NX + 1);
  violis.resize(NSTEP*NRESM + 1);
  alist.resize(NRESM + 1);
  bind.resize(NRESM + 1);
  bind0.resize(NRESM + 1);
  sort.resize(NRESM + 1);

  difx.resize(NX + 1);
  dd.resize(NX + 1);
  d0.resize(NX + 1);
  d.resize(NX + 1);
  xmin.resize(NX + 1);
  x1.resize(NX + 1);
  x0.resize(NX + 1);
  x.resize(NX + 1);
  gphi1.resize(NX + 1);
  gphi0.resize(NX + 1);
  qgf.resize(NX + 1);
  gradf.resize(NX + 1);
  diag0.resize(NX + 1);
  ug.resize(NX + 1);
  og.resize(NX + 1);
  xst.resize(NX + 1);
  xtr.resize(NX + 1);
  xsc.resize(NX + 1);

  resmin.resize(NRESM + 1);
  gresn.resize(NRESM + 1);
  diag.resize(NRESM + 1);
  cscal.resize(NRESM + 1);
  betaq.resize(NRESM + 1);
  colle.resize(NRESM + 1);
  u.resize(NRESM + 1);
  u0.resize(NRESM + 1);
  w.resize(NRESM + 1);
  w1.resize(NRESM + 1);
  res.resize(NRESM + 1);
  res0.resize(NRESM + 1);
  res1.resize(NRESM + 1);
  resst.resize(NRESM + 1);
  yu.resize(NRESM + 1);
  slack.resize(NRESM + 1);
  work.resize(NRESM + 1);
  delfac.resize(NRESM + 1);
  fu.resize(NRESM + 1);

  ddual.resize(NDUALM + 1);
  np.resize(NDUALM + 1);

  ud.resize(MDUALM + 1);
  ud1.resize(MDUALM + 1);

  gunit = Matrix<int>(4, NRESM + 1);
  gres = Matrix<double>(0.0, NX + 1, NRESM + 1);
  qr = Matrix<double>(0.0, NX + 1, NRESM + 1);
  a = SquareMatrix<double>(NX + 1, NX + 1);
  fugrad = Matrix<double>(NX + 1, NRESM + 1);
  xj = SquareMatrix<double>(NDUALM + 1, NDUALM + 1);
  r = Matrix<double>(NDUALM + 1, NDUALM + 1);
CATCH(QSpellucci::Spellucci::mf_InitialiseWorkSpace)
}

/* external symbol to allow class to be forced to be linked in */
bool QSpellucciLoad(){
	return (QSpellucci::TYPE != 0);
}

DRLIB_END_NAMESPACE
