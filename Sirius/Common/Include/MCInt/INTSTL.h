
#ifndef __INTSTL_H__
#define __INTSTL_H__


// Include required library header files
#include <time.h>
#include <math.h>
#include <assert.h>
#include <stdarg.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <string>



using namespace std;


#define INFOPATH string("C:\\Temp\\")

// Constants
const static double PI = 3.14159265358979;

typedef std::vector< long			> STLLongVector;
typedef std::vector< double			> STLDoubleVector;
typedef std::vector< int				> STLIntegerVector;
typedef std::vector< bool				> STLBoolVector;
typedef std::vector< STLDoubleVector	> STLDoubleVectorVector;
typedef std::vector< STLIntegerVector > STLIntegerVectorVector;




#define EQSP_CC_STL_PAIR(T1,T2)           std::pair< T1, T2 >


#ifdef	ERROR
#undef	ERROR
#endif
#define ERROR(message) { string err; err = message; throw( err ); }



#endif 



