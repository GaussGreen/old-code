/**
   This file contains various functions that are useful to be able to call
   within a debugger. Originally the functions were in the relevant .cpp file
   but Developer Studio can't cope with this when the application is built as
   a series of dlls. To fix this we include this file once per dll by
   including it each XXXXLib.cpp where XXXX is the directory name.
   Note that these functions are only for use within the debugger - which
   is why they're not inside a class etc.

   DO NOT INCLUDE THIS FILE FROM ANY OTHER FILE APART FROM THE XXXXLib.cpp FILES

   If you want to put some other functions in here which require access to
   other directories, then create eg MarketDebug.hpp and then include this file
   from there. Don't forget to replace the QLIB_TOOLKIT below with eg
   QLIB_MARKET. Then change eg MarketLib.cpp to include MarketDebug.hpp rather
   than this file.
 */
#ifndef QLIB_TREEDEBUG_HPP
#define QLIB_TREEDEBUG_HPP
#if defined (DEBUG) && (defined(QLIB_BUILD_DLL) || defined(QLIB_TREE))
// try to ensure that this only gets pulled in once per executable

#include "edginc/TreeSlice.hpp"

#endif
#endif
