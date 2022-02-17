
#include <stdarg.h>
#include <time.h>
#include <sys/timeb.h>
#include <windows.h>
#include "TraceARM.h"
#include "armglob.h"
#include "dates.h"
//#include "Exception.h"

using namespace ARM;

int TraceARM::traceLevel = NO_TRACE;
FILE* TraceARM::traceFile = NULL;

void TraceARM::Close() 
{
	if ( traceFile ) 
	{
		fclose(traceFile);
	}
};

void TraceARM::Print(int level, string& s) 
{			
	if ( traceFile && (traceLevel >= level) ) 
	{
		fprintf(traceFile, s.c_str());
		fflush(traceFile);
		::OutputDebugString(s.c_str()); 
	}
}

void TraceARM::Print(int level, const char* format, ...) 
{			
	if ( traceFile && (traceLevel >= level) ) 
	{
		va_list ap;
		va_start(ap, format);

		char buffer[10240];
		vsprintf(buffer, format, ap);
		va_end(ap);

		fprintf(traceFile, buffer);
		fflush(traceFile);
		::OutputDebugString(buffer); 
	}
}

void TraceARM::Print(int level, ARM_Object& t) 
{ 
	if ( traceFile && (traceLevel >= level) ) 
	{
		t.View("", traceFile); 
		fflush(traceFile);
	}
}

void TraceARM::Init(string& filename, int level) 
{ 
	TraceARM::Close();
	traceLevel = level; 
	if ( (traceFile = fopen(filename.c_str(), "w")) == NULL ) 
	{
		char msg[100];
		sprintf( msg,"TraceARM::Init, unable to open Trace file %s", filename.c_str());
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,msg );
		//throw exception("Trace::Init, unable to open Trace file %s", filename.c_str());
	}
}

string TraceARM::DefaultFilename()
{
	char buffer[128];

	struct tm* datetime;
	time_t ltime;
	time( &ltime );
	datetime = localtime( &ltime );

	struct _timeb millisec;
	_ftime( &millisec );

	sprintf(buffer, "Traces %4.4d_%2.2d_%2.2d %2.2d.%2.2d.%2.2d.%3u.txt", 1900+datetime->tm_year, datetime->tm_mon+1, 
		datetime->tm_mday, datetime->tm_hour, datetime->tm_min, datetime->tm_sec, millisec.millitm);

	return string(buffer);
}

string DatetoString(ARM_Date& date)
{
	char dt[32];
	date.JulianToStrDate(dt);
	return string(dt);
}

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
