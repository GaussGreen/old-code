#if 0
#if defined (__SUNPRO_CC)

#include "ksolaris_string.h"
#include "kexception.h"


DRString::DRString (const char* x)
{
	if (x == NULL) {
		Alloc(0);
	}
	else {
		Alloc(strlen(x));
		for (int i=0; i<m_length; i++) m_ptr[i] = x[i];
	}
}

DRString& DRString::operator=(const DRString &x)
{
	if (this == &x) return *this;
	resize(x.m_length);
	for (int i=0; i<m_length; i++) m_ptr[i] = x.m_ptr[i];
	return *this;
}

DRString& DRString::operator=(const char* x)
{
	resize (strlen(x)); 
	for (int i=0; i<m_length; i++) m_ptr[i] = x[i];
	return *this;
}

	
DRString DRString::substr(int pos, int n) const
{
	DRString ans;
	while (pos < m_length && n-->0) {
		ans += m_ptr[pos++];
	}

	return ans;
}

DRString DRString::operator+(const DRString &a)
{
	int i;
	int len = a.m_length + m_length;

	DRString res(len);
	for (i=0; i<m_length; i++) res.m_ptr[i]=m_ptr[i];
	for (i=0; i<a.m_length; i++) res.m_ptr[i+m_length] = a.m_ptr[i];

	return res;
}

DRString DRString::operator+(const char* x)
{
	int i;
	int xlen = strlen(x);
	int len = m_length + xlen;

	DRString res (len);
	for (i=0; i<m_length; i++) res.m_ptr[i]=m_ptr[i];
	for (i=0; i<xlen; i++) res.m_ptr[i+m_length] = x[i];

	return res;
}

DRString operator+(const char* x, const DRString &a)
{
	int i;
	int xlen = strlen(x);
	int len = a.m_length + xlen;

	DRString res (len);
	for (i=0; i<xlen; i++) res.m_ptr[i]=x[i];
	for (i=0; i<a.m_length; i++) res.m_ptr[i+a.m_length] = a.m_ptr[i];

	return res;
}

ostream &operator<<(ostream &stream, const DRString &a) 
{
	for (int i=0; i<a.m_length; i++) stream << a.m_ptr[i];
	return stream;
}

static inline bool solaris_string_is_white (char c) 
{return (c==' ' || c=='\n' || c=='\t' || c==13);}
 
istream& operator>>(istream& s, DRString& a)
{
	char c;

	DRString ans;

	while (!s.eof() && solaris_string_is_white (s.peek())) s.get();

	while (!solaris_string_is_white(c = s.get())) ans += c;

	a = ans;

	while (!s.eof() && solaris_string_is_white (s.peek())) s.get();

	return s;
}

void DRString::resize (int len)
{
	if (m_length != len) {
		if (m_ptr) delete [] m_ptr;
		Alloc(len);
	}
}

void DRString::Alloc (int len)
{
	if (len < 0)
		throw DRException ("DRStrings cannot have negative lengths") << len;

	m_length = len;
	m_ptr = new char[m_length+1];
	if (!m_ptr)
		throw DRException ("Out of Memory in DRString");

	m_ptr[m_length] = '\0';
}


#endif
#endif
