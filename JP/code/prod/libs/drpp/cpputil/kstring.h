#ifndef __drstring__h
/**@#-*/
#define __drstring__h
/**@#+*/

#include "kplatdep.h"
#include IOSTREAM_H
#include STRING_H
#include STDLIB_H

#if 0 /* bmc */

#if defined (__SUNPRO_CC)
#if 0
#include "ksolaris_string.h"
#endif

#include <string>
using std;
typedef string DRString;
#endif

/**@#-*/
#if defined (_MSC_VER)
#include <string>
typedef string DRString;
#endif
/**@#+*/

#endif /* bmc*/



/* bmc: */

#include <string>
typedef string DRString;


/** String class where all characters are in upper case */
class upper_string : public string {
public:
	/** Construct an empty string */
	upper_string () {}
	/** Construct from a char* */
	upper_string (const char* x);
	/** Construct from a string */
	upper_string (const string&);
	
	friend istream& operator>>(istream& s, upper_string& a)
	{s >> (string&) a; a.make_upper(); return s;}

protected:
	void make_upper();
};


/** String class where all characters are in lower case */
class lower_string : public string {
public:
	/** Construct an empty string */
	lower_string (){}
	/** Construct from a char * */
	lower_string (const char* x);
	/** Construct from a string */
	lower_string (const string&);

	friend istream& operator>>(istream&, lower_string&);

protected:
	void make_lower ();
};

// Creates a left substring
// If len > 0, then substring has length len
// eg string x = "abc"; y = left_substring(x, 2); then y becomes "ab"
// If len < 0, then substring is shortened by len  
// eg string x = "abc"; y = left_substring(x, -2); then y becomes "a"

string left_substring (const string&, int len);	


// Same as left_substring, but creates a right substring

string right_substring (const string&, int len); 

string convertPathDelimiters(string s);		

inline 
string CheckPathDelimiters(string s) {return convertPathDelimiters(s);}		


//f Converts "TRUE" or "FALSE" to bool
//bool toBool(string);
bool toBool(const string&);

#ifndef __SUNPRO_CC
inline string toString (bool b)
{
	return (b) ? "TRUE" : "FALSE";
}
#endif

//f Converts double to string
inline string toString(double x)
{
	char temp[20];
	sprintf(temp, "%lf", x);
	return string(temp);
}

//f Converts integer to string
inline string toString(int x)
{
	char temp[20];
	sprintf(temp, "%d", x);
	return string(temp);
}

inline int toInteger (string s)
{
	return atoi (s.c_str());
}

inline double toDouble (string s)
{
	return atof (s.c_str());
}

#endif
