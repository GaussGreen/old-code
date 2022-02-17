//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Description : Used to instantiate any stl templates that we have 
//                 declared as 'extern' (useful for reducing size of .libs)
//
//   Author      : Mark A Robson
//
//   Date        : 30 Sep 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#ifndef QLIB_BUILD_DLL
// NB We don't use config.hpp as a precompiled header as templates there are
// defined as extern (and MSVC gets upset if you declare a template as extern
// and then instantiate it)
#if defined(DEBUG) && defined(_MSC_VER)
#include <string>
using namespace std;
#pragma warning( disable : 4660 )
template basic_string<char, char_traits<char>, allocator<char> >;
#endif

#else
// needed to ensure _List_global<bool> is instantiated (stlport bug?)
#include <list>

#endif
