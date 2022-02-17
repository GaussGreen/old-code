/**
 * @file DECLARE.hpp
 */

#ifndef DRLIB_DECLARE_H
#define DRLIB_DECLARE_H

#include "edginc/Array.hpp"
#include <vector>

// switch off warning about template already instantiated - there appears to
// be a bug in MSVC 6 (at least) which means it warns about the template already
// being instantiated even though it hasn't see the the definition of the
// type(s) used by the template!
#if defined(DEBUG) && defined (_MSC_VER) && !(_MSC_VER >= 1300)
#pragma warning( disable : 4660 )
#endif

// NB You must use this macro inside a DRLIB_BEGIN_NAMESPACE etc
#define DECLARE(T)                                      \
    typedef smartConstPtr<T> T##ConstSP;                \
    typedef smartPtr<T> T##SP;                          \
    typedef array<T##SP, T> T##Array;                   \
    typedef smartPtr<T##Array> T##ArraySP;              \
    typedef smartConstPtr<T##Array> T##ArrayConstSP;

#define DECLARE_REF_COUNT(T)                    \
    typedef refCountPtr<const T> T##ConstSP;             \
    typedef refCountPtr<T> T##SP;                       \
    typedef std::vector<T##SP> T##Array;                     \
    typedef std::vector<T##ConstSP> T##ConstArray;            \
    typedef refCountPtr<T##Array> T##ArraySP;           \
    typedef refCountPtr<T##ConstArray> T##ConstArraySP; \
    typedef refCountPtr<const T##Array> T##ArrayConstSP;

#endif
