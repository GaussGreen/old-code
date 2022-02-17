///////////////////////////////////////////////////////////////////////////////
// Filename: WarningsOn.h : header file
// Author: Adam Gladstone
// Date: 
// 
// Description
//
//    Based on John Robbins - Microsoft Systems Journal Bugslayer Column
//    Turns back on all the warnings that we had to turn off just to get
//    the MS headers to compile. Corresponding file is <Warningsoff.h>
//
///////////////////////////////////////////////////////////////////////////////


/* '' : signed/unsigned mismatch */
#pragma warning ( default : 4018 )

/* compiler limit : terminating line number emission */
#pragma warning ( default : 4049 )

/* Unreferenced formal parameter */
#pragma warning ( default : 4100 )

/* conditional expression is constant */
#pragma warning ( default : 4127 )

/* unary minus operator applied to unsigned type, result still unsigned */
#pragma warning ( default : 4146)

/* local variable is initialized but not referenced */
#pragma warning ( default : 4189)

/* Nonstandard extension : nameless struct/union */
#pragma warning ( default : 4201 )

/* 'return' : conversion from 'int' to 'unsigned short', possible loss of data */
#pragma warning ( default : 4244 )

/*  'initializing' : conversion from 'const int' to 'const unsigned int', signed/unsigned mismatch */
#pragma warning ( default : 4245 )

/* cast truncates constant value */
#pragma warning ( default : 4310 )

/* unreferenced local function has been removed */
#pragma warning ( default : 4505 )

/* default constructor could not be generated */
#pragma warning ( default : 4510 )

/* 'codecvt_base' : copy constructor could not be generated */
#pragma warning ( default : 4511 )

/* destructor could not be generated */
#pragma warning ( default : 4513 )

/* Unreferenced in inline function removed */
#pragma warning ( default : 4514 )

/*  C++ exception handler used, but unwind semantics are not enabled. Specify -GX */
#pragma warning ( default : 4530 )

/* struct '' can never be instantiated - user defined constructor */
/* required */
#pragma warning ( default : 4610 )

/* STL spews errors just including them! */

/* C++ language change: to explicitly specialize class template '' use the following syntax: */
#pragma warning ( default : 4663 )

/* local variable '...' may be used without having been initialized */
#pragma warning ( default : 4701 )

/* identifier was truncated to '255' characters in the browser information */
#pragma warning ( default : 4786 )