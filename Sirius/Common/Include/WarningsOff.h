///////////////////////////////////////////////////////////////////////////////
// Filename: WarningsOff.h : header file
// Author: Adam Gladstone
// Date: 
// 
// Description
//
//    Based on John Robbins - Microsoft Systems Journal Bugslayer Column
//    Turns off warnings so standard headers will compile.  Any warning
//    turned off, must be turned back in in <WarningsOn.h>
//
///////////////////////////////////////////////////////////////////////////////

// '' : signed/unsigned mismatch 
#pragma warning ( disable : 4018 )

/*// compiler limit : terminating line number emission 
#pragma warning ( disable : 4049 )

// Unreferenced formal parameter 
#pragma warning ( disable : 4100 )

// conditional expression is constant 
#pragma warning ( disable : 4127 )*/

// unary minus operator applied to unsigned type, result still unsigned 
#pragma warning ( disable : 4146)

/*// local variable is initialized but not referenced 
#pragma warning ( disable : 4189)

// Nonstandard extension : nameless struct/union 
#pragma warning ( disable : 4201 )*/

// 'return' : conversion from 'int' to 'unsigned short', possible loss of data 
#pragma warning ( disable : 4244 )

/*//  'initializing' : conversion from 'const int' to 'const unsigned int', signed/unsigned mismatch 
#pragma warning ( disable : 4245 )

// return type for '' is 'const double*' (ie; not a UDT or reference to a UDT. Will produce errors if applied using infix notation)
#pragma warning(disable : 4284)

// cast truncates constant value 
#pragma warning ( disable : 4310 )

// unreferenced local function has been removed 
#pragma warning ( disable : 4505 )

// default constructor could not be generated 
#pragma warning ( disable : 4510 )

// 'codecvt_base' : copy constructor could not be generated 
#pragma warning ( disable : 4511 )

// destructor could not be generated 
#pragma warning ( disable : 4513 )

// Unreferenced in inline function removed 
#pragma warning ( disable : 4514 )

//  C++ exception handler used, but unwind semantics are not enabled. Specify -GX 
#pragma warning ( disable : 4530 )

// struct '' can never be instantiated - user defined constructor 
// required 
#pragma warning ( disable : 4610 )

// STL spews errors just including them! 

// C++ language change: to explicitly specialize class template '' use the following syntax: 
#pragma warning ( disable : 4663 )

// local variable '...' may be used without having been initialized 
#pragma warning ( disable : 4701 )*/

// identifier was truncated to '255' characters in the browser information 
#pragma warning ( disable : 4786 )
