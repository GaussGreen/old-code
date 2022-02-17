/**@#-*/
#ifndef solaris_string__h
#define solaris_string__h

#if 0
#if defined (__SUNPRO_CC)

#include "kplatdep.h"
#include STRING_H
#include SSTREAM_H

//k This is a simple string class.  
//k It's supposed to mimic STL string and is only used
//k for solaris.

class DRString {
public:
	DRString (const char* x = NULL); //
	//~DRString(); //
	~DRString() { if(m_ptr) delete[] m_ptr;} //

	DRString (const DRString&); //
	DRString& operator=(const DRString&); //
	DRString& operator=(const char*); //

	int size() const; // Gets length of string
    char& at(int loc) {return m_ptr[loc];} // Gets an individual char

	const char* c_str () const {return m_ptr;} //

	DRString operator+(const DRString&); // String concatenation
	DRString operator+(const char*); // String concatenation
	friend DRString operator+(const char*, const DRString&); // String concatenation

	DRString& operator+= (const char c); // Char concatenation
	DRString& operator+= (const char* c); // String concatenation
	DRString& operator+= (const DRString&); // String concatenation

	DRString substr(int pos, int n) const;

	friend bool operator==(const DRString&, const DRString&); //
	friend bool operator==(const DRString&, const char*); //
	friend bool operator==(const char*, const DRString&); //

	friend bool operator!=(const DRString&, const DRString&); //
	friend bool operator!=(const DRString&, const char*); //
	friend bool operator!=(const char*, const DRString&); //

	friend bool operator<(const DRString&, const DRString&); //
	friend bool operator<(const DRString&, const char*); //
	friend bool operator<(const char*, const DRString&); //

	friend ostream& operator<<(ostream&, const DRString&);  //
	friend istream& operator>>(istream&, DRString&);  //


private:
	char* m_ptr;
	int m_length;
	void Alloc(int);		// initializes DRString to be length len
	void resize(int);		// resizes a DRString to a different size
	DRString (int len) {Alloc(len);}	// private initializor
};

typedef DRString string;


inline DRString::DRString (const DRString &x)
{
	Alloc(x.m_length);
	for (int i=0; i<m_length; i++) m_ptr[i] = x.m_ptr[i];
}

//inline DRString::~DRString() 
//{
//	if (m_ptr) delete [] m_ptr;
//}

inline int DRString::size() const
{
	return m_length;
}

inline bool operator==(const DRString& a, const DRString &b)
{
	return (strcmp(a.m_ptr, b.m_ptr) == 0);
}

inline bool operator==(const DRString& a, const char* x)
{
	return (strcmp(a.m_ptr, x) == 0);
}

inline bool operator==(const char* x, const DRString &a) 
{
	return (strcmp(x, a.m_ptr) == 0);
}

inline bool operator!=(const DRString &a, const DRString& b)
{
	return (strcmp(a.m_ptr, b.m_ptr) !=0);
}

inline bool operator!=(const DRString& a, const char* x)
{
	return (strcmp(a.m_ptr, x) !=0);
}

inline bool operator!=(const char* x, const DRString &a) 
{
	return (strcmp(x, a.m_ptr) != 0);
}

inline bool operator<(const DRString& a, const DRString& b)
{
	return (strcmp(a.m_ptr, b.m_ptr)<0);
}

inline bool operator<(const DRString& a, const char* x)
{
	return (strcmp(a.m_ptr, x)<0);
}

inline bool operator<(const char* x, const DRString &a) 
{
	return (strcmp(x, a.m_ptr)<0);
}

inline DRString& DRString::operator+=(const char c)
{	
	char* temp = m_ptr;
	Alloc(m_length+1);
	strcpy(m_ptr, temp);
	delete [] temp;
	m_ptr[m_length-1] = c;
	return *this;
}

inline DRString& DRString::operator+=(const char* c)
{	
	char* temp = m_ptr;
	int old_length = m_length;

	Alloc(m_length+strlen(c));
	strcpy(m_ptr, temp);
	delete [] temp;

	strcpy(m_ptr + old_length, c);
	return *this;
}

inline DRString& DRString::operator+=(const DRString& s)
{	
	return (*this += s.c_str());
}

#endif

#endif
#endif


/**@#+*/
