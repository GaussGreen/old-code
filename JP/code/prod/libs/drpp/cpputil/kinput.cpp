// mbsio.cpp: implementation of the mbsio class.
//
//////////////////////////////////////////////////////////////////////

#include "kinput.h"
#include "kexception.h"
#include "ksymboltable.h"

static inline bool IsWhite(char c) {return (c == '\n' || c == '\t' || c == ' ' || c == 13);}

static inline void EatWhite(istream& inFile)
{char first; while (!inFile.eof() && IsWhite(first = inFile.peek())) {inFile.get();}}

static inline void EatTillNewLine (istream& inFile) {char c; while (!inFile.eof() && (c=inFile.get()) != '\n') {}}

static inline char EatWhiteAndPeek(istream& inFile) {EatWhite(inFile); return inFile.peek();}


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


KInput::KInput() : m_type(EMPTY)
{}

KInput::KInput(const String& s) 
: m_type(STRING)
{
	m_stringPtr = StringPtr (new String (s));
}

KInput::KInput(const KInputVector& v) 
: m_type(VECTOR)
{
	m_vectorPtr = KInputVectorPtr (new KInputVector (v));
}

KInput::KInput(const KInputMap& m) 
: m_type(MAP)
{
	m_mapPtr = KInputMapPtr (new KInputMap (m));
}

KInputVector& KInput::get_vector()
{
	if (m_type != VECTOR) 
		throw KException("Illegal Vector Get on ") << *this;
	
	return *m_vectorPtr;
}

KInputMap& KInput::get_map()
{
	if (m_type != MAP) {
		throw KException("Illegal Map Get on ") << *this;
	}	
	return *m_mapPtr;
}

KInput::String& KInput::get_string()
{
	if (m_type != STRING) 
		throw KException("Illegal String Get on ") << *this;
	
	return *m_stringPtr;
}


void KInput::open (string fileName)
{
	ParseFile(fileName);
}

void KInput::ParseFile(string fileName)
{
	string temp = convertPathDelimiters(theSymbolTable.get(fileName));

	ifstream inFile(temp.c_str());
	if (!inFile) {
	        throw KException ("Failed to open ") << temp << " converted from " << fileName;
	}

	Parse(inFile);

	EatWhite(inFile);
	while (!inFile.eof() && inFile.peek() >= 0) {

		if (inFile.peek() != '#')  
			throw KException("Extra junk in the file ") << fileName;

		EatTillNewLine(inFile);		
		EatWhite(inFile);
	}
}

KInput KInput::Parse(istream& inFile)
{
	char first = PeekNextChar(inFile);

	if (first == '&') {
		inFile.get();
		string fileName = GetNextString(inFile);
		ParseFile(fileName);
	}
	else if (first == '{') {
		inFile.get();
		ParseMap(inFile);
	}
	else if (first == '[') {
		inFile.get();
		ParseVector(inFile);
	}
	else if (first == '\"') {
		inFile.get();
		ParseQuote (inFile);
	}
	else {
		ParseString(inFile);
	}

	return *this;
}

void KInput::ParseQuote (istream& inFile)
{
	string temp;
	char c;

	while (!inFile.eof() && ((c = inFile.get()) != '\"')) {
		temp += c;
	}

	if (c != '\"')
		throw KException ("Mismatched quotes");

	m_type = STRING;
	m_stringPtr = StringPtr(new String(temp));
}

void KInput::ParseVector (istream& inFile)
{
	m_vectorPtr = KInputVectorPtr(new KInputVector());
	m_type = VECTOR;

	KInputVector& vector = *m_vectorPtr;

	char first = PeekNextChar(inFile);

	while (first != ']') {
		if (inFile.eof())
			throw KException("Unmatched ]");

		KInput input;
		input.Parse(inFile);
		vector.push_back(input);
		first = PeekNextChar(inFile);
	}
	GetNextString(inFile);
}

void KInput::ParseMap (istream& inFile)
{
	m_mapPtr = KInputMapPtr(new KInputMap());
	m_type = MAP;
	
	KInputMap& map = *m_mapPtr;
	
	string input = GetNextString(inFile);
	
	if (map.find(input) != map.end())
		throw KException("Duplicate Name in Map: ") << input;
	
	while (input != "}") {
		if (inFile.eof())
			throw KException("Unmatched }");
		
		KInput drinput;
		drinput.Parse(inFile);
		map[input] = drinput;
		
		if (input == "PATH_VARIABLES")
		{
			KInputMap& pathMap = drinput.get_map();
			for (KInputMap::iterator iter = pathMap.begin(); iter != pathMap.end(); iter++) {
				theSymbolTable[(*iter).first] = (*iter).second.get_string();
			}
		}			
		input = GetNextString(inFile);
	}		
}

void KInput::ParseString (istream& inFile)
{
	upper_string tempString = GetNextString(inFile);

	m_stringPtr = StringPtr(new String(tempString));
	m_type = STRING;
}

char KInput::PeekNextChar (istream& inFile)
{
	char c = EatWhiteAndPeek(inFile);

	while (c == '#') {
		EatTillNewLine(inFile);
		c = EatWhiteAndPeek(inFile);
	}
	return c;
}

upper_string KInput::GetNextString (istream& inFile)
{
	static char temp[100];
	inFile >> temp;

	while (temp[0] == '#') {
		EatTillNewLine(inFile);
		inFile >> temp;
	}

	int len = strlen (temp);
	char c = temp[len-1];
	if (len != 1 && (c == '}' || c == ']'))  {
		temp[len-1] = '\0';
		inFile.putback(c);
	}

	return upper_string (temp);
}


ostream& operator<<(ostream& s, const KInput& a)
{
	a.print(s);
	s << endl;
	return s;
}


upper_string KInputMap::get (upper_string name, upper_string def)
{
	KInputMap::iterator loc = find(name);
	if (loc == end()) {
		return def;
	}
	else {
		return (*loc).second.get_string();
	}
}

upper_string KInputMap::get (upper_string name)
{
	KInputMap::iterator loc = find(name);

	if (loc == end()) {
		throw KException("Failed to find ") << name << " in " << *this;
	}

	return (*loc).second.get_string();
}

KInput& KInputMap::operator[](upper_string name)
{
	return KMap(KInput::String, KInput)::operator[](name);
}

KInput& KInputMap::get_KInput (upper_string name)
{
	KInputMap::iterator loc = find(name);

	if (loc == end()) {
		throw KException("Failed to find ") << name << " in " << *this;
	}

	return (*loc).second;
}

bool KInputMap::isKey (upper_string name)
{
	return (find(name) != end());
}

bool KInput::operator==(const KInput& a) const
{
	if (m_type != a.m_type) return false;

	if (m_type == EMPTY)
		return true;
	else if (m_type == STRING)
		return *m_stringPtr == *(a.m_stringPtr);
	else if (m_type == VECTOR)
		return *m_vectorPtr == *(a.m_vectorPtr);
	else if (m_type == MAP)
		return *m_mapPtr == *(a.m_mapPtr);

	return true;
}


KInputMap& KInputMap::store (upper_string name, upper_string value)
{
	KInput tempIO (value);
	operator[](name) = tempIO;

	return *this;
}

istream& operator>>(istream& s, KInput& a)
{
	a.Parse(s);
	return s;
}

ostream& operator<<(ostream& s, const KInputVector& a)
{
	a.print(s);
	return s;
}


ostream& operator<<(ostream& s, const KInputMap& a)
{
	a.print(s);
	return s;
}

static void PrintTabs(ostream& s, int level) 
{
	for (int i = 0; i < level; i++) s << "  ";
}
	
void KInput::print(ostream& s, int level) const
{
	if (m_type == EMPTY) return;
	else if (m_type == STRING) {
		s << *m_stringPtr << "  ";
	}
	else if (m_type == MAP) {
		KInputMap& map = *m_mapPtr;
		map.print(s, level);
	}
	else {
		KInputVector& vector = *m_vectorPtr;
		vector.print(s, level);
	}
}

void KInputMap::print(ostream& s, int level) const
{
	s << "{" << endl;
	
	const_iterator mapIter;
	for (mapIter = begin(); mapIter != end(); mapIter++) {
		PrintTabs(s, level+1);
		s << (*mapIter).first << "  ";
		(*mapIter).second.print(s, level+1);
		s << endl;
	}

	PrintTabs(s, level);
	s << "} ";
}

void KInputVector::print(ostream& s, int level) const
{
	s << "[" << endl;

	const_iterator vecIter;

	for (vecIter = begin(); vecIter != end(); vecIter++) {
		PrintTabs(s, level+1);
		vecIter->print(s, level+1);
		s << endl;
	}

	PrintTabs(s, level);
	s << "] ";
}
