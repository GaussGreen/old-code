/***************************************************************
 * Module:	
 * Submodule:	
 * File:	kdrstring.h
 * Function:	Standard Definition and Include files.
 * Author:	Andrew Chou
 * Revision:	$Header$
 ***************************************************************/
#ifndef _kdrstring_H
#define _kdrstring_H

#if defined (__SUNPRO_CC) && !defined(DQSTL)

#include "kstdinc.h"



//--------------------------------------------------------------
/**
 * This is a simple string class.  It's supposed to mimic STL string and is only used for solaris.
 */

class KDRString {
public:
	/** Constructor. */
	KDRString (const char* x = NULL); 

	/** Constructor. */
	~KDRString(); 

	KDRString (const KDRString&);

	/** Constructor. */
	KDRString& operator=(const KDRString&);

	/** Constructor. */
	KDRString& operator=(const char*);

	/** Gets length of string. */
	int size() const;

	/** Gets an individual char. */
	char& at(int loc) {return m_ptr[loc];}

	/** */
	const char* c_str () const {return m_ptr;}

	KDRString operator+(const KDRString&); // String concatenation

	KDRString operator+(const char*); // String concatenation

	friend KDRString operator+(const char*, const KDRString&); // String concatenation

	KDRString& operator+= (const char c); // Char concatenation

	KDRString& operator+= (const char* c); // String concatenation

	KDRString& operator+= (const KDRString&); // String concatenation

	KDRString substr(int pos, int n) const;

	friend bool operator==(const KDRString&, const KDRString&); //

	friend bool operator==(const KDRString&, const char*); //

	friend bool operator==(const char*, const KDRString&); //

	friend bool operator!=(const KDRString&, const KDRString&); //

	friend bool operator!=(const KDRString&, const char*); //

	friend bool operator!=(const char*, const KDRString&); //

	friend bool operator<(const KDRString&, const KDRString&); //

	friend bool operator<(const KDRString&, const char*); //

	friend bool operator<(const char*, const KDRString&); //

	friend ostream& operator<<(ostream&, const KDRString&);  //

	friend istream& operator>>(istream&, KDRString&);  //


private:
				/** String. */
	char* m_ptr;
				/** Length. */
	int m_length;
				/** Initializes KDRString to be length len. */
	void Alloc(int);
				/** Resizes a KDRString to a different size. */
	void resize(int);
				/** Private initializor. */
	KDRString (int len) {Alloc(len);}
};

//typedef KDRString string;


inline KDRString::KDRString (const KDRString &x)
{
	Alloc(x.m_length);
	for (int i=0; i<m_length; i++) m_ptr[i] = x.m_ptr[i];
}

inline KDRString::~KDRString() 
{
	if (m_ptr) delete [] m_ptr;
}

inline int KDRString::size() const
{
	return m_length;
}

inline bool operator==(const KDRString& a, const KDRString &b)
{
	return (strcmp(a.m_ptr, b.m_ptr) == 0);
}

inline bool operator==(const KDRString& a, const char* x)
{
	return (strcmp(a.m_ptr, x) == 0);
}

inline bool operator==(const char* x, const KDRString &a) 
{
	return (strcmp(x, a.m_ptr) == 0);
}

inline bool operator!=(const KDRString &a, const KDRString& b)
{
	return (strcmp(a.m_ptr, b.m_ptr) !=0);
}

inline bool operator!=(const KDRString& a, const char* x)
{
	return (strcmp(a.m_ptr, x) !=0);
}

inline bool operator!=(const char* x, const KDRString &a) 
{
	return (strcmp(x, a.m_ptr) != 0);
}

inline bool operator<(const KDRString& a, const KDRString& b)
{
	return (strcmp(a.m_ptr, b.m_ptr)<0);
}

inline bool operator<(const KDRString& a, const char* x)
{
	return (strcmp(a.m_ptr, x)<0);
}

inline bool operator<(const char* x, const KDRString &a) 
{
	return (strcmp(x, a.m_ptr)<0);
}

inline KDRString& KDRString::operator+=(const char c)
{	
	char* temp = m_ptr;
	Alloc(m_length+1);
	strcpy(m_ptr, temp);
	delete [] temp;
	m_ptr[m_length-1] = c;
	return *this;
}

inline KDRString& KDRString::operator+=(const char* c)
{	
	char* temp = m_ptr;
	int old_length = m_length;

	Alloc(m_length+strlen(c));
	strcpy(m_ptr, temp);
	delete [] temp;

	strcpy(m_ptr + old_length, c);
	return *this;
}

inline KDRString& KDRString::operator+=(const KDRString& s)
{	
	return (*this += s.c_str());
}

#endif

#endif // _kdrstring_H


/**@#+*/
