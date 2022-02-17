#include "kstring.h"
#include "kexception.h"

upper_string::upper_string (const char* x) 
: string (x)
{
	make_upper();
}


upper_string::upper_string (const string& s) 
: string (s)
{
	make_upper();
}

void upper_string::make_upper()
{
	for (int i = 0; i < size(); i++) 
		at(i) = toupper (at(i));
}


string convertPathDelimiters(string s)
{
	char c;
	string ans;
	for (int i = 0; i < s.size(); i++) {
		c = s.at(i);
		if (c == '/' || c == 92) ans += PATH_DELIMITER;
		else ans += c;
	}


#ifdef _MSC_VER
	return ans;
#endif

	//#ifdef __SUNPRO_CC
	return lower_string(ans);
	//#endif

}

lower_string::lower_string (const char* x) 
: string (x)
{
	make_lower();
}


lower_string::lower_string (const string& s) 
: string (s)
{
	make_lower();
}

void lower_string::make_lower()
{
	for (int i = 0; i < size(); i++) 
		at(i) = tolower (at(i));
}

istream& operator>>(istream& s, lower_string& a)
{
	s >> (string&) a;
	a.make_lower();
	return s;
}

string left_substring (const string& s, int len)
{
	if (len < 0) len = s.size() + len;
	if (len > s.size()) 
		throw DRException("Invalid Length ") << len << " for string " << s;
	
	return s.substr(0, len);
}

string right_substring (const string& s, int len)
{
	int start;
	if (len < 0) start = -len;
	else start = s.size() - len;

	if (start >= s.size()) 
		throw DRException("Invalid Length ") << len << " for string " << s;
	
	return s.substr (start, s.size() - start); 
}

#if 0
bool toBool(string s)
{
	upper_string temp (s);
	if (temp == "TRUE") return true;
	else if (temp == "FALSE") return false;
	else throw DRException("Invalid boolean argument ") << s;

	return false;
}
#endif

bool toBool(const string& s)
{
	upper_string temp (s);
	if (temp == "TRUE") return true;
	else if (temp == "FALSE") return false;
	else throw DRException("Invalid boolean argument: ") << s;

	return false;
}





