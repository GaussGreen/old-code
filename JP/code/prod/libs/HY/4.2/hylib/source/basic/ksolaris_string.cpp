
#include "kplatform.h"
#if defined (__SUNPRO_CC)

#include "ksolaris_string.h"
#include "kexception.h"


KString::KString (const char* x)
{
	if (x == NULL) {
		Alloc(0);
	}
	else {
		Alloc(strlen(x));
		for (int i=0; i<m_length; i++) m_ptr[i] = x[i];
	}
}

KString& KString::operator=(const KString &x)
{
	if (this == &x) return *this;
	resize(x.m_length);
	for (int i=0; i<m_length; i++) m_ptr[i] = x.m_ptr[i];
	return *this;
}

KString& KString::operator=(const char* x)
{
	resize (strlen(x)); 
	for (int i=0; i<m_length; i++) m_ptr[i] = x[i];
	return *this;
}

	
KString KString::substr(int pos, int n) const
{
	KString ans;
	while (pos < m_length && n-->0) {
		ans += m_ptr[pos++];
	}

	return ans;
}

KString KString::operator+(const KString &a)
{
	int i;
	int len = a.m_length + m_length;

	KString res(len);
	for (i=0; i<m_length; i++) res.m_ptr[i]=m_ptr[i];
	for (i=0; i<a.m_length; i++) res.m_ptr[i+m_length] = a.m_ptr[i];

	return res;
}

KString KString::operator+(const char* x)
{
	int i;
	int xlen = strlen(x);
	int len = m_length + xlen;

	KString res (len);
	for (i=0; i<m_length; i++) res.m_ptr[i]=m_ptr[i];
	for (i=0; i<xlen; i++) res.m_ptr[i+m_length] = x[i];

	return res;
}

KString operator+(const char* x, const KString &a)
{
	int i;
	int xlen = strlen(x);
	int len = a.m_length + xlen;

	KString res (len);
	for (i=0; i<xlen; i++) res.m_ptr[i]=x[i];
	for (i=0; i<a.m_length; i++) res.m_ptr[i+a.m_length] = a.m_ptr[i];

	return res;
}

ostream &operator<<(ostream &stream, const KString &a) 
{
	for (int i=0; i<a.m_length; i++) stream << a.m_ptr[i];
	return stream;
}

static inline bool solaris_string_is_white (char c) 
{return (c==' ' || c=='\n' || c=='\t' || c==13);}
 
istream& operator>>(istream& s, KString& a)
{
	char c;

	KString ans;

	while (!s.eof() && solaris_string_is_white (s.peek())) s.get();

	while (!solaris_string_is_white(c = s.get())) ans += c;

	a = ans;

	while (!s.eof() && solaris_string_is_white (s.peek())) s.get();

	return s;
}

void KString::resize (int len)
{
	if (m_length != len) {
		if (m_ptr) delete [] m_ptr;
		Alloc(len);
	}
}

void KString::Alloc (int len)
{
	if (len < 0)
		throw KException ("KStrings cannot have negative lengths") << len;

	m_length = len;
	m_ptr = new char [m_length+1];
	if (!m_ptr)
		throw KException ("Out of Memory in KString");

	m_ptr[m_length] = '\0';
}


#endif
