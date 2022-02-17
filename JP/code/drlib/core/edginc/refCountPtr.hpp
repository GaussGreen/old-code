#ifndef EDG_REFCOUNTPTR_H
#define EDG_REFCOUNTPTR_H

#include "edginc/ModelException.hpp"

CORE_BEGIN_NAMESPACE

/**
   An alternative approach would be to use TR1::boost_shared_ptr, but so far only gcc-4.x is shippped with it.
   #include <tr1/memory>
   using std::tr1::boost_shared_ptr;

*/

#define refCountPtr shared_ptr
//#define refCountConstPtr shared_ptr

// NullDeleter can be passed as second argument of refcountPtr constructor:
// ptrSP = refCountPtr<T>( &some_obj , NullDeleter())
// this can be used to wrap a pointer into a shared one without taking ownership of it.

struct NullDeleter
{
	NullDeleter(){};
	~NullDeleter() {}
    void operator()(void const *) const {}
};


/// Some Qlib clients rely on lessThan template
template <class T>
bool lessThan( const refCountPtr<T>& first, const refCountPtr<T>& second )
{
    return (!first) || (!second) || (*first < *second);
}

#define DYNAMIC_POINTER_CAST boost::dynamic_pointer_cast
#define CONST_POINTER_CAST  boost::const_pointer_cast

CORE_END_NAMESPACE

#endif
