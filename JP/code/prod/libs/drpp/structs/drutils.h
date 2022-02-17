#ifndef __drutils__h
#define __drutils__h

#include "drstring.h"
#include FSTREAM_H

extern "C" {
#include "macros.h"
#include "convert.h"
}

//f Converts string to double
double toNumber(const string&);

//f Returns a string with the current time (HH:MM:SS  M/D/Y)
DRString GetTimeStamp ();

//f Inline max function (didn't template cause of solaris)
inline double DRMax (double a, double b) {return (a>b) ? a : b;}
//f Inline min function
inline double DRMin (double a, double b) {return (a>b) ? b : a;}

//f Given f(low) = lowVal, f(high) = highVal, we
//f geometrically interpolate to get f(pos)
double GeometricInterp (double pos, double low, double high, double lowVal, double highVal);
//f Given f(low) = lowVal, f(high) = highVal, we
//f linearly interpolate to get f(pos)
double LinearInterp (double pos, double low, double high, double lowVal, double highVal);

//f (Spot-Strike)+
double Call (double spot, double strike); 
//f (Strike-Spot)+
double Put (double spot, double strike); 

//f (Spot - Strike)
double IntrinsicCall (double spot, double strike);
//f (Strike-Spot)
double IntrinsicPut (double spot, double strike); 

//f Reads from the stream until the first non-white character
void EatWhite(istream&);
//f Reads from the stream until the carriage return
void EatTillNewLine (istream&);
//f Reads from the stream till the first non-white char 
//f and returns the next char 
char EatWhiteAndPeek(istream&);



//f Converts day count name to GTO integer value
//f See web page for valid names
//int toDayCount(string);
int toDayCount(const string&);


//f Checks whether a char is line feed, carriage return,
//f space or tab
inline bool IsWhite(char c)
{
	return (c == '\n' || c == '\t' || c == ' ' || c == 13);
}

inline
void MyCountedOut (ostream& out, TDate* dates)
{
	int n = (int) dates[0];

	out << "[ ";
	for (int i = 1; i <= n; i++)
		out << GtoFormatDate(dates[i]) << " ";
	out << "]\n";
}

inline 
void MyCountedOut (ostream& out, double* values)
{
	int n = (int) values[0];

	out << "[ ";
	for (int i = 1; i <= n; i++)
		out << values[i] << " ";
	out << "]\n";
}

inline 
void MyCountedOut (ostream& out, int* values)
{
	int n = (int) values[0];

	out << "[ ";
	for (int i = 1; i <= n; i++)
		out << values[i] << " ";
	out << "]\n";
}

inline 
void MyCountedOut (ostream& out, char* values)
{
	int n = (int) values[0];

	out << "[ ";
	for (int i = 1; i <= n; i++)
		out << &values[GTO_STR_POS(i)] << " ";
	out << "]\n";
}



#endif
