/***************************************************************
 * Module:	PenGuin
 * File:	kutilios.h
 * Function:	Standard Definition and Include files.
 * Author:	Christian Daher
 * Revision:	$Header$
 ***************************************************************/
#ifndef	_kutilios_H
#define	_kutilios_H

#include "kstdinc.h"


//--------------------------------------------------------------
// IO utilities

	/**
	 * Stream checking utility.
	 *
	 * Checks that a stream is in a good state or throws
	 * and exception.
	 * @exception KFailure
	 */
	ios&		check(ios& s, const char* debugName = NULL);

    /**
	 * Stream checking utility.
	 *
	 * Checks that a stream is in a good state. If so,
     * return true, else return false
	 */
    bool checkSilent(ios& s);

	/**
	 * Read data from an input stream.
	 * 
	 * Get the next line skipping blank lines and comments.
	 * Returns static string.
	 * @exception KFailure
	 */
	char*		getNextLine(istream &is, const char *debugName = "");


	/**
	 * Read data from an input stream.
	 *
	 * Get the next token (a block of characters either enclosed
	 * in double brackets or not containing any space)
	 * and returns a static copy of it.
	 * (removes any double quotes).<br>
	 * WARNING: does not skip comments '\#'.
	 * @exception KFailure
	 */
	char*		getToken(istream &is, char *separators,
				const char *debugName = "");

	/**
	 * Read data from an input stream.
	 *
	 * Scans a string value in an istream
	 * and returns a static copy of it.
	 * Skips all commented lines (#) and ignores double quotes (")
	 * around strings.
	 * @exception KFailure
	 */
	const char*	getString(istream& is, const char *debugName = "");

	/**
	 * Read data from an input stream.
	 *
	 * Scans a single char in an istream and returns it.
	 * Skips all commented lines (#) and ignores double quotes (")
	 * around strings.
	 * @exception KFailure
	 */
	char		getChar(istream& is, const char *debugName = "");


	/**
	 * Read data from an input stream.
	 *
	 * Scans a double value in an istream as an interval
	 * (e.g. S/N, 1D, etc.)
	 * and returns it.
	 * Skips all commented lines (#) and ignores double quotes (")
	 * around strings.
	 * @exception KFailure
	 */
	double		getDoubleInterval(istream& is,
				const char *debugName = "");

	/**
	 * Read data from an input stream (supports currency format).
	 *
	 * Scans a double value in an istream and returns it.
	 * Skips all commented lines (#) and ignores double quotes (")
	 * around strings. Recognises currency format (e.g. 1,000.00)
	 * @exception KFailure
	 */
	double		getDouble(istream& is, const char *debugName = "");

	/**
	 * Read data from an input stream.
	 *
	 * Scans an int value in an istream and returns it.
	 * Skips all commented lines (#) and ignores double quotes (")
	 * around strings.
	 * @exception KFailure
	 */
	int		getInt(istream& is, const char *debugName = "");

	/**
	 * Read data from an input stream.
	 *
	 * Scans a long value in an istream and returns it.
	 * Skips all commented lines (#) and ignores double quotes (")
	 * around strings.
	 * @exception KFailure
	 */
	long		getLong(istream& is, const char *debugName = "");

	/**
	 * Read data from an input stream.
	 *
	 * Scans a TDate value in an istream and returns it.
	 * Skips all commented lines (#) and ignores double quotes (")
	 * around strings.
	 * @exception KFailure
	 */
	TDate		getTDate(istream& is, const char *debugName = "");

	/**
	 * Read data from an input stream.
	 *
	 * Scans a TDateInterval value in an istream and returns it.
	 * Skips all commented lines (#) and ignores double quotes (")
	 * around strings.
	 * @exception KFailure
	 */
	TDateInterval	getTDateInterval(istream& is,
				const char *debugName = "");

	/**
	 * Read data from an input stream.
	 *
	 * Scans a TDayCount value in an istream and returns it.
	 * Skips all commented lines (#) and ignores double quotes (")
	 * around strings.
	 * @exception KFailure
	 */
	TDayCount	getTDayCount(istream& is,
				const char *debugName = "");
	/**
	 * Advance to a token in an input stream.
	 *
	 * Advances in the stream util the token is encountered.
	 * (a block of characters either enclosed
	 * in double quotes or not containing any space).
	 * @exception KFailure
	 */
	istream&	advanceIs(istream& is, char *token);

	/**
	 * Utility to use printf style output to a stream.
	 *
	 * Format in printf style and returns pointer to static string.
	 * For example,
	 * <pre><dir>
	 * cout << format("value = %8.4f\n", 1.23);
	 * </pre></dir>
	 */
	char*		format(const char* fmt, ...);


	/**
	 * Print double and enclose in parenthesis () if negative
	 */
	string		formatDouble(const char* fmt, double c);

	/**
	 * Reads and puts the contents of a file "fnam"
	 * in a string (STL vector of char).
	 */
	String		DppFileToString(const char *fnam);


	/**
	 * Same, but vector of chars.
	 */
	KVector(char)	DppFileToVectChar(const char *fnam);



//---------------------------------------------------------------
//
//


/**
 * Reads an object from a file.
 */

template<class T> inline void
DppReadFromFile(const char* fnam, T& object)
{
	try {
		ifstream is(fnam);
		check(is, fnam);
		is >> object;
	}
	catch (...) {
		throw KFailure("ReadFromFile: failed on `%s'.\n", fnam);
	}
}


/**
 * Writes an object to a file (truncate file).
 */

template<class T> inline void
DppWriteToFile(const char* fnam, const T& object)
{
	ofstream os(fnam);
	os << object;
}


/**
 * Writes an object to a file (appends to file).
 */

template<class T> inline void
DppAppendToFile(const char* fnam, const T& object)
{
	ofstream os(fnam, ios::app);
	os << object;
}


/**
 * Reads an object from a char* string.
 */

template<class T> inline void
DppReadFromString(const char* str, T& object)
{
	try {
		istrstream is(str);
		check(is, str);
		is >> object;
	}
	catch (...) {
		throw KFailure("ReadFromString: failed on `%s'.\n", str);
	}
}



/**
 * Writes an object from a char* string.
 */

template<class T> inline String
DppWriteToString(T& object)
{
	ostrstream os;
	os << object;
	return String(os.str());
}






#endif


