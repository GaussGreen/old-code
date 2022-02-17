/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: $
 *
 */

/*! \file Trace.cpp
 *
 *  \brief Class Trace enables tracing in a file using a trace level 
 *	\author  R. Anger
 *	\version 1.0
 *	\date December 2003
 */

#include <stdarg.h>
#include <time.h>
#include <sys/timeb.h>
#include <windows.h>
#include "Trace.h"
#include "armglob.h"
#include "dates.h"
#include "Exception.h"

using namespace mercure;

int Trace::traceLevel = NO_TRACE;
FILE* Trace::traceFile = NULL;

void Trace::Close() 
{
	if ( traceFile ) 
	{
		fclose(traceFile);
	}
};

void Trace::Print(int level, string& s) 
{			
	if ( traceFile && (traceLevel >= level) ) 
	{
		fprintf(traceFile, s.c_str());
		fflush(traceFile);
		::OutputDebugString(s.c_str()); 
	}
}

void Trace::Print(int level, const char* format, ...) 
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

void Trace::Print(int level, ARM_Object& t) 
{ 
	if ( traceFile && (traceLevel >= level) ) 
	{
		t.View("", traceFile); 
		fflush(traceFile);
	}
}

void Trace::Init(string& filename, int level) 
{ 
	Trace::Close();
	traceLevel = level; 
	if ( (traceFile = fopen(filename.c_str(), "w")) == NULL ) 
	{
		throw NestedException("Trace::Init, unable to open Trace file %s", filename.c_str());
	}
}

string Trace::DefaultFilename()
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

string toString(ARM_Date& date)
{
	char dt[32];
	date.JulianToStrDate(dt);
	return string(dt);
}

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
