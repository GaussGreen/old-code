

#ifndef _TRACEARM_H
#define _TRACEARM_H

#include <stdio.h>
#include <string>
using namespace std;

class ARM_Object;

namespace ARM 
{

enum TraceLevelType
{
	NO_TRACE,
	LOW_LEVEL_TRACE,
	HIGH_LEVEL_TRACE
};

class TraceARM
{
	public:
		static void Print(int level, string& s);		
		
		static void Print(int level, const char* format, ...);
		
		static void Print(int level, ARM_Object& t);

		static void Init(string& filename, int level);
		static void Close();
		
		static string DefaultFilename();

	private:
		static int traceLevel;
		static FILE* traceFile;
};

}  //namespace

//
//
//
#ifndef	ARMLOG
	#define ARMLOG(logLevel,logMsg) \
	{ \
		std::stringstream _sstr; _sstr<<logMsg; \
		ARM::TraceARM::Print(logLevel,_sstr.str()); \
	} 
#endif // MERCLOG

class ARM_Date;

string DatetoString(ARM_Date& date);

#endif //_TRACEARM_H

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
