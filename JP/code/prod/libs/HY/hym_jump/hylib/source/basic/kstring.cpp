#include "kstring.h"
#include "kexception.h"

upper_string::upper_string (const char* x) 
: KString (x)
{
	make_upper();
}


upper_string::upper_string (const KString& s) 
: KString (s)
{
	make_upper();
}

void upper_string::make_upper()
{
	for (int i = 0; i < size(); i++) 
		at(i) = toupper (at(i));
}


KString convertPathDelimiters(KString s)
{
	char c;
	KString ans;
	for (int i = 0; i < s.size(); i++) {
		c = s.at(i);
		if (c == '/' || c == 92) ans += PATH_DELIMITER;
		else ans += c;
	}

#ifdef __SUNPRO_CC
	return lower_string(ans);
#endif

#ifdef _MSC_VER
	return ans;
#endif
}

lower_string::lower_string (const char* x) 
: KString (x)
{
	make_lower();
}


lower_string::lower_string (const KString& s) 
: KString (s)
{
	make_lower();
}

void lower_string::make_lower()
{
	for (int i = 0; i < size(); i++) 
		at(i) = tolower (at(i));
}

std::istream& operator>>(std::istream& s, lower_string& a)
{
	s >> (KString&) a;
	a.make_lower();
	return s;
}

KString left_substring (const KString& s, int len)
{
	if (len < 0) len = s.size() + len;
	if (len > s.size()) 
		throw KException("Invalid Length ") << len << " for KString " << s;
	
	return s.substr(0, len);
}

KString right_substring (const KString& s, int len)
{
	int start;
	if (len < 0) start = -len;
	else start = s.size() - len;

	if (start >= s.size()) 
		throw KException("Invalid Length ") << len << " for KString " << s;
	
	return s.substr (start, s.size() - start); 
}






