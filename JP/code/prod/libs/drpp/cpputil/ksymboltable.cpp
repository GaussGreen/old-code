// DRSymbolTable.cpp: implementation of the DRSymbolTable class.
//
//////////////////////////////////////////////////////////////////////

#include "ksymboltable.h"
#include "kexception.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

inline static bool IsLetter (char c) {return (c >= 'A' && c <= 'Z') || (c>='a' && c <='z');}

string KSymbolTable::get (string s) 
{
	const char* f = s.c_str();
	string ans;

	int i = 0;
	char c;

	while ((c = f[i++]) != '\0') {
		if (c == '$') {
			DRString temp;

			while (IsLetter(c = f[i++]) || c == '_') temp += c;	

			iterator iter = find(temp);

			if (iter == end())
				throw DRException("Invalid path variable ") << temp;

			ans += get((*iter).second);
			ans += c;
		}
		else ans += c;
	}
	return ans;
}


KSymbolTable theSymbolTable;





