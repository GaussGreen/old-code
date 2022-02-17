#ifndef __kstring__h
/**@#-*/
#define __kstring__h
/**@#+*/

#include "kplatdep.h"
#include "kexception.h"
#include IOSTREAM_H
#include STRING_H
#include STDLIB_H

#if defined (__SUNPRO_CC)
#include "ksolaris_string.h"
#endif

//using namespace std;

/**@#-*/
#if defined (_MSC_VER)
#include <string>
typedef std::string KString;
#endif
/**@#+*/

/** String class where all characters are in upper case */
class upper_string : public KString {
public:
	/** Construct an empty KString */
	upper_string () {}
	/** Construct from a char* */
	upper_string (const char* x);
	/** Construct from a KString */
	upper_string (const KString&);
	
	friend std::istream& operator>>(std::istream& s, upper_string& a)
	{s >> (KString&) a; a.make_upper(); return s;}

protected:
	void make_upper();
};


/** String class where all characters are in lower case */
class lower_string : public KString {
public:
	/** Construct an empty KString */
	lower_string (){}
	/** Construct from a char * */
	lower_string (const char* x);
	/** Construct from a KString */
	lower_string (const KString&);

	friend std::istream& operator>>(std::istream&, lower_string&);

protected:
	void make_lower ();
};

// Creates a left substring
// If len > 0, then substring has length len
// eg KString x = "abc"; y = left_substring(x, 2); then y becomes "ab"
// If len < 0, then substring is shortened by len  
// eg KString x = "abc"; y = left_substring(x, -2); then y becomes "a"

KString left_substring (const KString&, int len);	


// Same as left_substring, but creates a right substring

KString right_substring (const KString&, int len); 

KString convertPathDelimiters(KString s);		

inline 
KString CheckPathDelimiters(KString s) {return convertPathDelimiters(s);}		


inline KString toString (bool b)
{
	return (b) ? "TRUE" : "FALSE";
}
inline KString toString(double x)
{
	char temp[20];
	sprintf(temp, "%lf", x);
	return KString(temp);
}

inline KString toString(int x)
{
	char temp[20];
	sprintf(temp, "%d", x);
	return KString(temp);
}



inline	bool toBool(const KString& s)
{
	upper_string temp (s);
	if (temp == "TRUE") return true;
	else if (temp == "FALSE") return false;
	else throw KException("Invalid boolean argument: ") << s;

	return false;
}


inline int toInteger(const KString &s)
{
	KString	tmp(s);
	return atoi (tmp.c_str());
}

inline double toDouble (const KString &s)
{
	KString	tmp(s);
	return atof (tmp.c_str());
}

#endif

