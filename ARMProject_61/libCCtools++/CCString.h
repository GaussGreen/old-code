/* -------------------------------------------------------------------------
 
   File: %M%
   Path: %P%
   Description: header CCString 
   Created: 98/06/02
   Author: Charles-Emmanuel MUSY
   Modified: %E% %U%
   Last maintained by: Charles-Emmanuel MUSY
   Revision: %I%
 
   -------------------------------------------------------------------------
 
   Note:
 
   ------------------------------------------------------------------------- */
 
#ifndef CCSTRING_H
#define CCSTRING_H

#include <CCcommon.h>
SCCS_ID (CCString_h_SccsId, "%W%, modified %E%");

#include <string.h>
#include <ctype.h>

#include <iosfwd>

#ifdef STL_WIN32
#include <vector>
#define VECTOR	std::vector
#else
#include <vector.h>
#define VECTOR	vector
#endif	// STL_WIN32



CCEXTERN_FUNCTION (char *CCtrim_right, (char* s));


class CCString
{
	friend std::ostream& operator<< (std::ostream& os, const CCString& s);

public:

	// constructeur void et chaine
	inline CCString (const char *n_str = NULL) { str = NULL; Set (n_str); }

    inline CCString (int n, const char *n_str = NULL) { str = NULL; Set (n, n_str); }


	// constructeur de copie
	inline CCString (const CCString& s) { str = NULL; Set (s); }

	// destructeur
	~CCString ();

	void Set (const char *n_str = NULL);
	void Set (int n, const char* n_str);
	void Set (const CCString& s) { Set ((const char *)s); }

	void Parser (char separator, VECTOR<CCString>& list); 
	void ParserWithEmpty (char separator, VECTOR<CCString>& list);

	void Replace (const char* n_str, char mask);
	void Replace (const char* strMask, const char* n_str);
	
	bool Contain (const CCString& chaine);

    double XL_value_convert ();

	inline void trim_right () { CCtrim_right (str); }

    // operateurs de recuperation de chaine allouee
	inline char *GetStr () const { return str ? strdup (str) : NULL; }
	
	/// interface C++
	char* c_str() const;
	operator char*() const { return c_str(); }

        // operateur de recuperation du buffer interne
	inline operator const char* () const { return str; }

	inline int GetLen () const { return len; }

	inline int operator! () const { return str == NULL; }

	CCString &operator= (const CCString& s);
	CCString &operator= (const char* s);
	CCString operator+ (const CCString &s) const;
	CCString &operator+= (const CCString& s);
	int operator== (const CCString& s) const;
	int operator== (const char*) const;
	int operator== (char*) const;
	int operator!= (const CCString& s) const;
	inline char operator[] (int index) const { return str[index]; }

	void toUpper ();
	
private:
	char *str;
	int len;
};

std::string CCSTringToSTLString( const CCString& s );


#ifdef CCString_c



#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include <CCmessage.h>

#define POINT	'.'
#define COMMA	','

#endif	// CCString_c

#endif	// STRING_H

// EOF %M%