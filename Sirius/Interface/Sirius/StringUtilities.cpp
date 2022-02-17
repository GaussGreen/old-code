/* -*-mode: C++; style: K&R; c-basic-offset: 4 ; -*- */

/*
 *  utility.cc - Simple portable string utility functions.
 *
 *  GNU MP3D - A portable(ish) MP3 server.
 *
 * Homepage:
 *   http://www.gnump3d.org
 *
 * Author:
 *  Steve Kemp <skx@tardis.ed.ac.uk>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
 *
 *  Steve Kemp
 *  ---
 *  http://www.steve.org.uk/
 *
 */


#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#ifndef ASSERT
# include <assert.h>
# define ASSERT assert
#endif /* ASSERT */

#include <algorithm>
#include <ctype.h>
#include <limits.h>
#include <stdio.h>

#include "debug.h"
#include "mutex.h"
//#include "utility.h"
#include "version.h"


/**
 *  Dummy object which will record the name and version of this file.
 */
InstallVersion dummyUtility( __FILE__,  __DATE__, "$Revision: 1.33 $" );


/**
 * Protected constructor - no instances of this class can be instantiated.
 */
CUtility::CUtility()
{
}

/**
 * Destructor.
 */
CUtility::~CUtility()
{
}


/**
 * Replace all occurences of one string with another in the given string.
 *
 * @param orig The string that should have the replacement applied to it.
 * @param from The substring which should be replaced.
 * @param to The replacement substring.
 * @param pos The starting position to make the replacement from.
 * @return The input string, with all occurrences of 'from' changed to 'to'.
 */
std::string CUtility::replaceStrInStr( std::string orig, std::string from, std::string to, std::string::size_type pos /* = 0 */ )
{
    CScopedMutex protection;

    size_t b = pos;
    for (;;)
    {
	b = orig.find(from, b);
	if (b == orig.npos) 
	    break;
	orig.replace(b, from.size(), to);
	b += to.size();
    }

    return orig;
}


/**
 * Convert part of a string to bold format, by wrapping it with html tags.
 * The segment of the string which is to be made bold is treated case 
 * insensitively.
 * @param orig The string which contains a substring to make bold
 * @param ebold The substring which should be made bold.
 * @return The potentially modified string.  
 */
std::string CUtility::boldString( std::string orig, std::string tobolden )
{
    CScopedMutex protection;

    std::string uo = CUtility::uc( orig );
    std::string ut = CUtility::uc( tobolden );

    size_t start = uo.find( ut );
    if ( start != std::string::npos )
    {
	std::string front = orig.substr( 0, start );
	std::string end   = orig.substr( start + tobolden.size() );

	/* TODO: Make these <b> literals parameters. */
	orig = front + "<b>" + tobolden + "</b>" + end ;
    }
    return( orig );
}

/**
 * Remove a section of text from a string.
 * Everything between the given start + end markers INCLUSIVE is removed.
 * Nothing will be removed unless BOTH tags are found.
 * @param orig The original text.
 * @param start The start tag to mark the place where removal should occur.
 * @param end The end tag to mark where removal should stop.
 * @return The modified text.
 */
/* static */ std::string CUtility::removeText( std::string orig, std::string start, std::string end )
{
    CScopedMutex protection;

    std::string::size_type sStart = orig.find( start );
    std::string::size_type sEnd   = orig.find( end );

    /*
     * return the original string if both markers weren't found.
     */
    if ( ( sStart == std::string::npos ) || ( sEnd == std::string::npos ) )
    {
	return( orig);
    }

    /*
     * Sanity check:
     * Make sure the `start` string occurs before the `end` string. 
     */
    ASSERT( sStart <= sEnd );

    std::string head = orig.substr( 0, sStart );
    std::string tail = orig.substr( sEnd + end.size() );

    return( head + tail );
}


/**
 * Normalize a path, by removing duplicate '/' characters.
 * @param path The path to normalize.
 * @return The normalized path.
 */
/* static */ std::string CUtility::normalizePath( std::string path )
{
    CScopedMutex protection;

    if ( path.find( "//" ) == std::string::npos )
    {
	return( path );
    }

    std::string::size_type length = path.size();
    char *p    = (char*)path.c_str();
    char c = '\0';
    std::string result;

    for( std::string::size_type i = 0; i < length ; i++ )
    {
	char x = p[i];
	if ( x == '/' )
	{
	    if ( x != c )
	    {
		result += x;
	    }
	}
	else
	{
	    result += x;
	}
	c = x;
    }
    return( result );
}


/**
 * Simple string splitter.
 * Return a list of all the substrings contained in the given
 * string, which are seperated by the given character.
 * @param str The string to process.
 * @param sep The character at which to split the string.
 * @return A vector containing the individual strings.
 */
std::vector< std::string > * CUtility::splitString( std::string str, char sep )
{
    CScopedMutex protection;

    std::vector< std::string > *results  = new std::vector< std::string >;

    std::string::size_type offset = 0;
    while ( ( offset = str.find( sep ) ) != std::string::npos )
    {
	std::string temp = str.substr( 0, offset  );
	if ( temp.size() > 0 )
	{
	    results->push_back( temp );
	}
	str = str.substr( offset+1 );
    }

    /* Make sure we get the last segment. */
    if ( str.size() > 0 )
    {
	results->push_back( str );
    }

    return( results );
}

/*
 * Convert a two byte hex string to an actual character.
 *
 * eg.  20 -> ' '
 */
unsigned int hex2int(const std::string& _s)
{
    CScopedMutex protection;

    unsigned int t = 0;
    std::string::const_iterator i;
    for(i=_s.begin();i!=_s.end();i++)
	if(*i>='0' && *i<='9')
	    t = (t<<4)+(*i-'0');
	else if(*i>='a' && *i<='f')
	    t = (t<<4)+(*i-'a'+10);
	else if(*i>='A' && *i<='F')
	    t = (t<<4)+(*i-'A'+10);
    else
	break;
    return t;
}


/**
 * Encode a string for use in an URL - as per RFC 1738
 * @param path The original string.
 * @return A string suitable for use as an URL.
 */
std::string CUtility::URLEncode( std::string path )
{
    CScopedMutex protection;

    unsigned char isAcceptable[96] =
	{/* 0x0 0x1 0x2 0x3 0x4 0x5 0x6 0x7 0x8 0x9 0xA 0xB 0xC 0xD 0xE 0xF */
	    0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0xF,0xE,0x0,0xF,0xF,0xC, /* 2x  !"#$%&'()*+,-./   */
	    0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0x8,0x0,0x0,0x0,0x0,0x0, /* 3x 0123456789:;<=>?   */
	    0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF, /* 4x @ABCDEFGHIJKLMNO   */
	    0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0x0,0x0,0x0,0x0,0xF, /* 5X PQRSTUVWXYZ[\]^_   */
	    0x0,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF, /* 6x `abcdefghijklmno   */
	    0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF,0x0,0x0,0x0,0x0,0x0  /* 7X pqrstuvwxyz{\}~DEL */
	};
    char *hexChars = "0123456789ABCDEF";

#define ACCEPTABLE(a)	( a>=32 && a<128 && ((isAcceptable[a-32]) & 0x8))

    std::string result;

    std::string::size_type length = path.size();

    for ( std::string::size_type i = 0; i < length; i++ )
    {
	unsigned char a = path[ i ];

	if ((!ACCEPTABLE(a))  || ( a == '+' ) )
	{
	    result += '%';	/* Means hex coming */
	    result += hexChars[a >> 4];
	    result += hexChars[a & 15];
	}
	else
	{
	    result += a;
	}
    }

    return( result );
}


/**
 * Decode a string from URL format - as per RFC 1738
 * @param path The original URL encoded string.
 * @return A normal string.
 */
std::string CUtility::URLDecode( std::string path )
{
    CScopedMutex protection;

    std::string decoded;
    std::string::size_type length = path.length();

    for( std::string::size_type i=0; i< length; i++)
    {
	if( path[i] == '+')
	{
	    decoded += ' ';
	}
	else if( ( path[i] == '%')  && ( i <= (length-2) ) )
        {
	    decoded += (char) hex2int(std::string(path,i+1,2));
	    i+=2;
	}
	else
	{
	    decoded += path[i];
	}
    }
    return decoded;
}

/**
 * Convert a string to upper case.
 * @param text The text to convert.
 * @return The given text converted to upper case.
 */
std::string CUtility::uc( std::string text )
{
    CScopedMutex protection;

    std::string s = text;
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    return( s );
}

/**
 * Convert a string to lower case.
 * @param text The text to convert.
 * @return The given text converted to lower case.
 */
std::string CUtility::lc( std::string text )
{
    std::string s = text;
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return( s );
}

/**
 * Compare two strings, case insensitively, like strcmp
 * @param a The first string.
 * @param b The second string.
 * @return A value describing the equality - like strcmp
 */
int CUtility::casecmp( std::string a, std::string b )
{
    CScopedMutex protection;

    std::string uA = CUtility::uc( a );
    std::string uB = CUtility::uc( b );

    return( strcmp( uA.c_str(), uB.c_str() ) );
}


/**
 * Remove all leading and trailing whitespace characters from a string
 * @param str The string to modify.
 * @return The modified string.
 */
std::string CUtility::stripWhiteSpace( std::string str )
{
    CScopedMutex protection;

    std::string trimchars = " \n\r\t\f\v";
    size_t b = str.find_first_not_of(trimchars);
    size_t e = str.find_last_not_of(trimchars);
	
    return b == str.npos ? "" : str.substr(b, e+1-b);
}


/**
 * Test to see if a string ends with a given string.
 * @param str The string to test.
 * @param suffix The suffix to look for.
 * @return 1 if the string has the given suffix, 0 otherwise.
 */
int CUtility::endsWith( std::string str, std::string suffix )
{
    CScopedMutex protection;

    /* Look for the string. */
    std::string::size_type i = str.rfind( suffix );
    if ( i == std::string::npos )
    {
	return 0;
    }

    /* We've found the suffix - but is it at the end? */
    if ( i == ( str.size() - suffix.size() ) )
    {
	/* Found the suffix at the end of the string. */
	return 1;
    }

    return 0;
}



/**
 * Test to see if the given string has a given suffix.
 * @param str The string to test.
 * @param prefix The prefix to test.
 * @return 1 if the string begins with the given suffix, 0 otherwise.
 */
int CUtility::startsWith( std::string str, std::string prefix )
{
    CScopedMutex protection;

    int sLen = str.size();
    int pLen = prefix.size();

    /* Early termination.  The string must be at least equal the the prefix
     * in size - otherwise it can't contain it. */
    if ( pLen > sLen )
    {
	return 0;
    }

    return ( str.substr(0,pLen) == prefix );
}


/**
 * Convert an integer to a string.
 * @param integer The integer to convert to a base 10 string.
 * @return The string value holding the character.
 */
std::string CUtility::int2str( int integer )
{
    CScopedMutex protection;

    std::string result = "";
    char *buf          = NULL;

    /* Maximum number of decimal digits the biggest integer could have. */
    int nDigits = INT_MAX % 10;

    /* Add one for '-' if negative, and one for the string terminator. */
    nDigits += 2;

    /* Create + format the string */
    buf = new char[ nDigits + 1];
    memset( buf, '\0', nDigits+1);

    if ( buf != NULL )
    {
	snprintf( buf, nDigits -1, "%d", integer );

	/* Copy results - ensure memory is free'd */
	result = buf;
	delete [] buf;
    }

    return( result );
}

