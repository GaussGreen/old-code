#include "drutils.h"
//#include "drmcprice_globals.h"
#include "drstack.h"
#include "drsymboltable.h"
#include "mbsio.h"

bool DR_LOGGING_ON = false;
bool DR_STDOUTPUT_ON = false;
ofstream DR_LOGGING_STREAM;

bool MID_OFFICE_STREAM_ON;
ofstream MID_OFFICE_STREAM;
DRString DEALID;

ofstream DEAL_OUTPUT_STREAM;

DRString GetTimeStamp ()
{
	tm* theTime;
	time_t long_time;
	
	time( &long_time );  
	theTime = localtime( &long_time );
	
	char tempString[80];
	sprintf (tempString, "%02d:%02d:%02d   %02d/%02d/%04d", 
		theTime->tm_hour, 
		theTime->tm_min,
		theTime->tm_sec,
		theTime->tm_mon + 1, 
		theTime->tm_mday,
		theTime->tm_year+1900);

	return DRString (tempString);
}


double GeometricInterp (double pos, double low, double high, double lowVal, double highVal)
{
	if (low == high) return lowVal;
	double w1 = (high - pos)/(high - low);
	double w2 = 1. - w1;

	double ans = w1 * log(lowVal) + w2 * log(highVal);
	ans = exp(ans);

	return ans;
}

double LinearInterp (double pos, double low, double high, double lowVal, double highVal)
{
	if (low == high) return lowVal;
	double w1 = (high - pos)/(high - low);
	double w2 = 1. - w1;

	double ans = w1 * lowVal + w2*highVal;
	return ans;
}


double Call (double spot, double strike) {return DRMax(spot-strike, 0);}

double Put (double spot, double strike) {return DRMax(strike-spot,0);}

double IntrinsicCall (double spot, double strike) {return spot-strike;}

double IntrinsicPut (double spot, double strike) {return strike-spot;}

void EatWhite(istream& inFile)
{
	char first;
	while (!inFile.eof() && IsWhite(first = inFile.peek())) {
		inFile.get();
	}
}

void EatTillNewLine (istream& inFile)
{
	char c;
	while (!inFile.eof() && (c=inFile.get()) != '\n') {}
}

char EatWhiteAndPeek(istream& inFile)
{
	EatWhite(inFile);
	return inFile.peek();
}


double toNumber(const string& s)
{
	drstack<double> myStack;
	
	ISTRINGSTREAM in(s.c_str());

	double num;
	char temp[100];
	while (!in.eof()) {
		in >> temp;

		if (strlen(temp) == 0)
			break;
		else if ((temp[0] >= '0' && temp[0] <= '9') || temp[0] == '.' || 
			(strlen(temp) > 1 && temp[0] == '-')) {
			myStack.push(atof(temp));
		}
		else if (temp[0] == '&') {
			DRString fileName = temp + 1;
			ifstream inFile (CheckPathDelimiters(theSymbolTable.get(fileName)).c_str());
			MbsIO input (inFile);
			num = toNumber(input.get_string());
			myStack.push(num);
		}
		else if (temp[0] == '$') {
			DRString variableName = temp + 1;
			upper_string temp (variableName);
			num = toNumber (theSymbolTable[temp]);
			myStack.push(num);
		}
		else if (strcmp(temp,"+")==0) {
			num = myStack.top(); myStack.pop();
			num += myStack.top(); myStack.pop();
			myStack.push(num);	
		}
		else if (strcmp(temp,"*")==0) {
			num = myStack.top(); myStack.pop();
			num *= myStack.top(); myStack.pop();
			myStack.push(num);	
		}
		else if (strcmp(temp,"-")==0) {
			num = -myStack.top(); myStack.pop();
			num += myStack.top(); myStack.pop();
			myStack.push(num);	
		}
		else if (strcmp(temp,"/")==0) {
			num = 1./myStack.top(); myStack.pop();
			num *= myStack.top(); myStack.pop();
			myStack.push(num);	
		}
		else 
			throw DRException ("Badly formed expression: ") << s;
	}
	double ans = myStack.top();
	myStack.pop();

	if (!myStack.empty())
		throw DRException ("Stack is not empty.  Suspect exp: ") << s;

	return ans;
}


int toDayCount(const string& t)
{
	upper_string s (t);

	if (s == "GTO_ACT_365") return 1;
	else if (s == "GTO_ACT_365F") return 2;
	else if (s == "GTO_ACT_360") return 3;
	else if (s == "GTO_B30_360") return 4;
	else if (s == "GTO_B30E_360") return 5;
	else throw DRException("Invalid day count: ") << s;
	
	return 0;
}

#if 0
int toDayCount(string t)
{
	upper_string s (t);

	if (s == "GTO_ACT_365") return 1;
	else if (s == "GTO_ACT_365F") return 2;
	else if (s == "GTO_ACT_360") return 3;
	else if (s == "GTO_B30_360") return 4;
	else if (s == "GTO_B30E_360") return 5;
	else throw DRException("Invalid day count: ") << s;
	
	return 0;
}
#endif
