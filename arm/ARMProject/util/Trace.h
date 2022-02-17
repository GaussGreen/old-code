/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: $
 *
 */

/*! \file Trace.h
 *
 *  \brief Class Trace enables tracing in a file using a trace level 
 *	\author  R. Anger
 *	\version 1.0
 *	\date December 2003
 */

#ifndef _TRACE_H
#define _TRACE_H

#include <stdio.h>
#include <string>
using namespace std;

class ARM_Object;

namespace mercure 
{

enum TraceLevelType
{
	NO_TRACE,
	LOW_LEVEL_TRACE,
	HIGH_LEVEL_TRACE
};

/*! \brief enables to print string, char* and ARM_Object in a trace file
	       depending on a trace level
*/
class Trace
{
	public:
		static void Print(int level, string& s);		
		//! \brief prints a formatted string as printf does
		static void Print(int level, const char* format, ...);
		//! \brief prints an ARM_Object using its method View
		static void Print(int level, ARM_Object& t);

		static void Init(string& filename, int level);
		static void Close();
		/*! \brief creates a default (unique) filename at date time millisec 
			       "Traces 2003_12_17 15.15.01.250.txt"
		*/
		static string DefaultFilename();

	private:
		static int traceLevel;
		static FILE* traceFile;
};

} // namespace

//
//
//
#ifndef	MERCLOG
	#define MERCLOG(logLevel,logMsg) \
	{ \
		std::stringstream _sstr; _sstr<<logMsg; \
		mercure::Trace::Print(logLevel,_sstr.str()); \
	} 
#endif // MERCLOG

class ARM_Date;

string toString(ARM_Date& date);

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
