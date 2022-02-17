/**@#-*/
#ifndef solaris_string__h
#define solaris_string__h

#include "kplatform.h"
#if defined (__SUNPRO_CC)

#include "kplatdep.h"
#include STRING_H
#include SSTREAM_H

//k This is a simple string class.  
//k It's supposed to mimic STL string and is only used
//k for solaris.

class KString {
public:
	KString (const char* x = NULL); //
	KString (const std::string& other) 
	  {
	    operator+=( other.c_str() );
	  }
	~KString(); //

	operator std::string()
	  {
	    return std::string( m_ptr );
	  }

	KString (const KString&); //
	KString& operator=(const KString&); //
	KString& operator=(const char*); //

	int size() const; // Gets length of string
    char& at(int loc) {return m_ptr[loc];} // Gets an individual char

	const char* c_str () const {return m_ptr;} //

	KString operator+(const KString&); // String concatenation
	KString operator+(const char*); // String concatenation
	friend KString operator+(const char*, const KString&); // String concatenation

	KString& operator+= (const char c); // Char concatenation
	KString& operator+= (const char* c); // String concatenation
	KString& operator+= (const KString&); // String concatenation

	KString substr(int pos, int n) const;

	friend bool operator==(const KString&, const KString&); //
	friend bool operator==(const KString&, const char*); //
	friend bool operator==(const char*, const KString&); //

	friend bool operator!=(const KString&, const KString&); //
	friend bool operator!=(const KString&, const char*); //
	friend bool operator!=(const char*, const KString&); //

	friend bool operator<(const KString&, const KString&); //
	friend bool operator<(const KString&, const char*); //
	friend bool operator<(const char*, const KString&); //

	friend ostream& operator<<(ostream&, const KString&);  //
	friend istream& operator>>(istream&, KString&);  //


private:
	char* m_ptr;
	int m_length;
	void Alloc(int);		// initializes KString to be length len
	void resize(int);		// resizes a KString to a different size
	KString (int len) {Alloc(len);}	// private initializor
};

//typedef KString string;


inline KString::KString (const KString &x)
{ 
	Alloc(x.m_length);
	for (int i=0; i<m_length; i++) m_ptr[i] = x.m_ptr[i];
}

inline KString::~KString() 
{
	if (m_ptr) delete [] m_ptr;
}

inline int KString::size() const
{
	return m_length;
}

inline bool operator==(const KString& a, const KString &b)
{
	return (strcmp(a.m_ptr, b.m_ptr) == 0);
}

inline bool operator==(const KString& a, const char* x)
{
	return (strcmp(a.m_ptr, x) == 0);
}

inline bool operator==(const char* x, const KString &a) 
{
	return (strcmp(x, a.m_ptr) == 0);
}

inline bool operator!=(const KString &a, const KString& b)
{
	return (strcmp(a.m_ptr, b.m_ptr) !=0);
}

inline bool operator!=(const KString& a, const char* x)
{
	return (strcmp(a.m_ptr, x) !=0);
}

inline bool operator!=(const char* x, const KString &a) 
{
	return (strcmp(x, a.m_ptr) != 0);
}

inline bool operator<(const KString& a, const KString& b)
{
	return (strcmp(a.m_ptr, b.m_ptr)<0);
}

inline bool operator<(const KString& a, const char* x)
{
	return (strcmp(a.m_ptr, x)<0);
}

inline bool operator<(const char* x, const KString &a) 
{
	return (strcmp(x, a.m_ptr)<0);
}

inline KString& KString::operator+=(const char c)
{	
	char* temp = m_ptr;
	Alloc(m_length+1);
	strcpy(m_ptr, temp);
	delete [] temp;
	m_ptr[m_length-1] = c;
	return *this;
}

inline KString& KString::operator+=(const char* c)
{	
	char* temp = m_ptr;
	int old_length = m_length;

	Alloc(m_length+strlen(c));
	strcpy(m_ptr, temp);
	delete [] temp;

	strcpy(m_ptr + old_length, c);
	return *this;
}

inline KString& KString::operator+=(const KString& s)
{	
	return (*this += s.c_str());
}

#endif

#endif


/**@#+*/


