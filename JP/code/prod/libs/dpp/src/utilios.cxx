/************************************************************************
 * Module:	PenGuin
 * File:
 * Function:    
 * Author:      C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include <ctype.h>
#include <errno.h>
#include "kutilios.h"
#include "kstlutil.h"


extern	"C" {
#include "cgeneral.h"
#include "cerror.h"
#include "bastypes.h"
#include "macros.h"             /* MAX */
#include "ldate.h"              /* GtoDayCountFraction */
#include "ratelist.h"		/* TRateList */
#include "tcurve.h"
#include "date_sup.h"

#include "drlio.h"		/* FScanStruct */
#include "drlstr.h"		/* StringLineVarScan */
#include "drltime.h"		/* TDatePrint */
#include "drloptio.h"		/* Black */
};

//---------------------------------------------------------------
//

bool checkSilent(ios& s) {

    if (s.bad()) {
		return false;
	}
	if (s.eof()) {
		return false;
	}
	if (s.fail()) {
		return false;
	}
	if (!s.good()) {
		return false;
	} else {
		return true;
	}
}

ios&
check(ios& s, const char* debugName)
{
	if (debugName) {
		if (s.bad()) {
			throw KFailure("ios::badbit on `%s'.\n"
				       "%s.\n",
					debugName,
					strerror(errno));
		}
		if (s.eof()) {
			throw KFailure("ios::eofbit on `%s'.\n"
				       "%s.\n",
					debugName,
					strerror(errno));
		}
		if (s.fail()) {
			throw KFailure("ios::failbit on `%s'.\n"
				       "%s.\n",
					debugName,
					strerror(errno));
		}
		if (!s.good()) {
			throw KFailure("ios::<unknown state> on `%s'.\n"
				       "%s.\n",
					debugName,
					strerror(errno));
		} else {
			//throw KFailure("ios::goodbit on %s.\n", debugName);
		}
	} else {
		if (s.bad()) {
			throw KFailure("ios::badbit.\n"
				       "%s.\n",
					strerror(errno));
		}
		if (s.eof()) {
			throw KFailure("ios::eofbit.\n"
				       "%s.\n",
					strerror(errno));
		}
		if (s.fail()) {
			throw KFailure("ios::failbit.\n"
				       "%s.\n",
					strerror(errno));
		}
		if (!s.good()) {
			throw KFailure("ios::<unknown state>.\n"
				       "%s.\n",
					strerror(errno));
		} else {
			//throw KFailure("ios::goodbit\n");
		}
	}
	return(s);
}

//---------------------------------------------------------------


static	inline
istream& eatwhite(istream &is)
{
	char	c;

	// skip over '#' commented lines 
	while (is.get(c)) {
	    if (!isspace(c)) {
		is.putback(c);
		break;
	    } else if (c == '#') {
		while ((is.get(c)) && (c != '\n'));
	    }
	}
	return is;
}


//---------------------------------------------------------------


static	inline
istream& getBuffer(
	istream &is,			// 
	char *s,			// 
	const char *debugName)		// 
{
	char    *q;
	char	c;

    try {

	q = s;

/*
	// skip over '#' commented lines 
	check(is, debugName); 
	while (is.get(c) && ((c == '#') || (isspace(c)))) {
            if (c == '#') {
		while ((is.get(c)) && (c != '\n'));
		check(is, debugName);
            }
        }


	// Scan string (possibly enclosed in "")
	check(is, debugName);
	if (c != '"') {
		*q++ = c;
	        while ((is.get(c)) && (!isspace(c)))
			*q++ = c;
		*q = '\0';
	} else {
	        while ((is.get(c)) && (c != '"'))
			*q++ = c;
		*q = '\0';
	}
	eatwhite(is);
*/

	// skip over '#' commented lines 
	check(is, debugName); 
	while (is.get(c) && ((c == '#') || (isspace(c)))) {
            if (c == '#') {
		while ((is.get(c)) && (c != '\n'));
		check(is, debugName);
            }
        }


	// Scan string (possibly enclosed in "")
	check(is, debugName);
	if (c != '"') {
		*q++ = c;
	        while ((is.get(c)) && (!isspace(c)))
			*q++ = c;
		*q = '\0';
	} else {
	        while ((is.get(c)) && (c != '"'))
			*q++ = c;
		*q = '\0';
	}

	eatwhite(is);


	return is;
    }
    catch (KFailure) {
	throw KFailure("getBuffer: failed on %s.\n",
			(debugName?debugName:"N/A")); 
    }
/*
    catch (...) {
	dppErr << "getBuffer: failed on " <<
		(debugName?debugName:"N/A") << " (uncaught)." << endl;
	throw KFailure();
    }
*/
}

//---------------------------------------------------------------



char*
getToken(istream &is, char *separators, const char *debugName)
{
	char	*q;
	char	c;
static	char	buf[256];

	q = buf;

	while ((is.get(c)) && (isspace(c)));

	if (c != '"') {
		*q++ = c;
	        while ((is.get(c)) && (!isspace(c)))
			*q++ = c;
		*q = '\0';
	} else {
	        while ((is.get(c)) && (c != '"'))
			*q++ = c;
		*q = '\0';
	}

	check(is, debugName);

	return buf;
}


//---------------------------------------------------------------

const char*
getString(istream& is, const char *debugName)
{
static	char	buf[1024];

	getBuffer(is, buf, debugName);
	return (buf);
}

//---------------------------------------------------------------

char
getChar(istream& is, const char *debugName)
{
	char	buf[128];

	getBuffer(is, buf, debugName);
	return (buf[0]);
}


//---------------------------------------------------------------

double
getDouble(istream& is, const char *debugName)
{
	char	buf[128];
	double	 value;

	getBuffer(is, buf, debugName);
	if (DrlCurScan(buf , &value) != SUCCESS) {
		is.clear(ios::badbit);
		throw KFailure("Failed scanning %s in `%s'.\n",
			debugName, buf);
	}
	return value;
}


//---------------------------------------------------------------
// Scans a int value in an istream and returns it.
// Throws a KFailure exception if fails.

int
getInt(istream& is, const char *debugName)
{
	char	buf[128];
	int	value;

	{
		getBuffer(is, buf, debugName);
		if (sscanf(buf, "%d", &value) != 1) {
			is.clear(ios::badbit);
			throw KFailure("Failed scanning %s in `%s'.\n",
				debugName, buf);
		}
	}


	return value;
}



//---------------------------------------------------------------

long
getLong(istream& is, const char *debugName)
{
	char	buf[128];
	long	value;

	getBuffer(is, buf, debugName);
	if (sscanf(buf, "%ld", &value) != 1) {
		is.clear(ios::badbit);
		throw KFailure("Failed scanning %s in `%s'.\n",
			debugName, buf);
	}

	return value;
}


//---------------------------------------------------------------

double
getDoubleInterval(istream& is, const char *debugName)
{
	char	buf[128];
	double	value;

	getBuffer(is, buf, debugName);
	if (DrlDoubleIntervalScan(buf, &value) != SUCCESS) {
		is.clear(ios::badbit);
		throw KFailure("Failed scanning %s in `%s'.\n",
			debugName, buf);
	}

	return value;
}



//---------------------------------------------------------------

TDate
getTDate(istream& is, const char *debugName)
{
	char	buf[128];
	TDate	value;

	getBuffer(is, buf, debugName);
	if (GtoStringToDate(buf, &value) != SUCCESS) {
		is.clear(ios::badbit);
		throw KFailure("Failed scanning %s in `%s'.\n",
			debugName, buf);
	}

	return value;
}



//---------------------------------------------------------------

TDateInterval
getTDateInterval(istream& is, const char *debugName)
{
	char	buf[128];
	TDateInterval	value;

	getBuffer(is, buf, debugName);
	if (DrlTDateIntervalScan(buf, &value) != SUCCESS) {
		is.clear(ios::badbit);
		throw KFailure("Failed scanning %s in `%s'.\n",
			debugName, buf);
	}

	return value;
}



//---------------------------------------------------------------

TDayCount
getTDayCount(istream& is, const char *debugName)
{
	char	buf[128];
	TDayCount	value;

	getBuffer(is, buf, debugName);
	if (GtoStringToDayCountConv(buf, &value) != SUCCESS) {
		is.clear(ios::badbit);
		throw KFailure("Failed scanning %s in `%s'.\n",
			debugName, buf);
	}

	return value;
}






//---------------------------------------------------------------


char*
getNextLine(istream &is, const char *debugName)
{
static	char	buf[1024];
	int	acceptFlag = 0;
	char	*p;

	do {
		// read next line 
		if (!is.getline(buf, sizeof(buf))) {
			throw KFailure("Failed scanning %s.\n",
				debugName);
		}
		check(is, debugName);

		// check for empty line (comment # ;) */
		p = buf;
		while ((isspace(*p)) && (*p != '\0')) p++;
		if ((*p != '#') && (*p != ';') && (*p != '\0'))
			acceptFlag = 1;
	} while (acceptFlag == 0);
	eatwhite(is);
	return(buf);
}











//---------------------------------------------------------------
//

istream&
advanceIs(istream& is, char *token)
{
	char	*p;
    try {
	for (;;) {
	    p = getToken(is, NULL);
	    if (!strcmp(p, token)) return (is);
	}
    }
    catch (...) {
	throw KFailure("advanceIs: can't find `%s'.\n",  token);
	check(is);
	throw KFailure();
    }
}





//---------------------------------------------------------------


char*
format(const char* fmt, ...)
{
static	char	buf[1024];
	va_list	ap;

	va_start(ap, fmt);
	vsprintf(buf, fmt, ap);
	va_end(ap);
	return(buf);
}


//---------------------------------------------------------------
//
string
formatDouble(const char* fmt, double c)
{
	string	outStr;
	
	string::size_type	pos = 0;

	outStr = string(format(fmt, c));
	
	//
	// Test for negative value and enclosed in parathesis
	//
	if((pos = outStr.find("-", pos)) != string::npos)
	{
		outStr.insert(pos, "(");
		outStr += ")";
	}

	return outStr;	

}


//---------------------------------------------------------------



KVector(char)
DppFileToVectChar(const char *fnam)
{
	KVector(char)	v;
	char		c;

	ifstream	is(fnam);
	check(is, fnam);

	while (is.get(c))
		v.insert(v.end(), c);

	return(v);

}



String
DppFileToString(const char *fnam)
{
	String	v;
	char	c;

	ifstream	is(fnam);
	check(is, fnam);

	while (is.get(c))
		v.insert(v.end(), c);

	return(v);

}



/*
 * Creates a string vector form a char-block style array.
 */
 
KVector(String)
DppStrVectorFromCharBlock(int n, char* array)
{
        int     i;
        KVector(String) strVector;
	String	str;

        for (i=0; i<n; i++)
	{
		str = &array[WRAP_STR_IDX(i+1)];
                strVector.insert(strVector.end(), str);
	}

        return(strVector);
}



/*
 * Creates a TDate vector form a CMLIB Date array.
 */
 
KVector(TDate)
DppDateArrayToTDateVector(const Array<Date> &a)
{
	int i;
	KVector(TDate) v;
 
	for (i=0; i<a.size(); i++)
		v.insert(v.end(), a[i]);
	return(v);
 
}

