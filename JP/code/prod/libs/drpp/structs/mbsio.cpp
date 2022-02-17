// mbsio.cpp: implementation of the mbsio class.
//
//////////////////////////////////////////////////////////////////////

#include "mbsio.h"
#include "drexception.h"
#include "drutils.h"

extern "C" {
#include "macros.h"
}

//void MakeCounted (DRInputVector& vec, TDate* &dates)
//{
//	dates = new TDate [vec.size() + 1];
//	dates[0] = vec.size();
//
//	int i = 1;
//	for (DRInputVector::iterator iter = vec.begin(); iter != vec.end();	iter++, i++) {
//		dates[i] = toDate((*iter).get_string());
//	}
// }

void MakeCounted (DRInputVector& vec, double* &values)
{
	values = new double [vec.size() + 1];
	values[0] = vec.size();

	int i = 1;
	for (DRInputVector::iterator iter = vec.begin(); iter != vec.end();	iter++, i++) {
		values[i] = toNumber((*iter).get_string());
	}
}

void MakeCounted (DRInputVector& vec, int* &values)
{
	values = new int [vec.size() + 1];
	values[0] = vec.size();

	int i = 1;
	for (DRInputVector::iterator iter = vec.begin(); iter != vec.end();	iter++, i++) {
		values[i] = toNumber((*iter).get_string());
	}
}

void MakeCounted (DRInputVector& vec, char* &values)
{
	values = new char [vec.size() * (GTO_MAX_STR_LEN + 1) + GTO_NUM_BYTES_FOR_STR_CNT];
	values[0] = vec.size();

	int i = 1;
	for (DRInputVector::iterator iter = vec.begin(); iter != vec.end();	iter++, i++) {
		strcpy (&values[GTO_STR_POS(i)],  ((*iter).get_string()).c_str());
	}
}

static double evalTable (DRInputMap& myMap)
{
	DRString key = myMap.get("Key");
	double add = toNumber(myMap.get("Add", "0"));
	double mult = toNumber (myMap.get("Mult", "1"));
	DRInputMap& table = myMap["TABLE"].get_map();
	double value = toNumber(table.get(key));	
	value = mult*value + add;
	return value;
}

double toNumber (DRInput& input)
{
	if (input.is_string()) 
		return toNumber(input.get_string());
	else if (input.is_map()) {
		return evalTable (input.get_map());
	}
	else
		throw DRException ("Can't convert to number");
	return 0;
}
