/****************************************************************************/
/*      Utility routines.                                                   */
/****************************************************************************/
/*      UTIL.c                                                              */
/****************************************************************************/


#include <stdio.h>
#include <stdlib.h>  
#include <string.h>

#include "fix123head.h"
#include "fix123plus.h"


int 
CompareManyStrings(
	const char			*str,				/* String to find */
	int					nbStrings,			/* Number of strings */
	const char			*stringList[])		/* Strings to find it in */
//////////////////////////////////////////////////////////////////////
// Finds a single string in a list of others.  Returns the position 
// in the list.
{
	int					i;					/* Loop index */

	/* Test each string */
	for (i=0; i<nbStrings; i++)
		if (strcmp(str, stringList[i]) == 0) return i;

	/* If we get this far, we didn't find it. */
	return -1;
}


void FindAndSkipComLine_CPP(
	FILE				*stream,			// (I) Stream
	const char			*Label,				// (I) Label
	const char			*Routine,			// (I) Calling routine
	const char			*FileName)			// (I) Filename
//////////////////////////////////////////////////////////////////////
// Finds a single string in a list of others.  Returns the position 
// in the list.
{
	int					rc;					// Result code

	if ((rc = FindAndSkipComLine(stream, (char*)Label, (char*)Routine, (char*)FileName)) != SUCCESS)
		throw "FindAndSkipComLine_CPP: Failed.";
}

	
char ReadCharFromFile(
	FILE				*stream,			// (I) Stream
	const char			*Label,				// (I) Label
	const char			*Routine,			// (I) Calling routine
	const char			*FileName)			// (I) Filename
//////////////////////////////////////////////////////////////////////
// Reads a single character from a file.
{
	char				ErrorData[256];		// Error message
	char				data;				// Data to read
	int					readerror;			// How many things read?

	// Read it.
    if ((readerror = fscanf (stream, "%c \n", &data)) != 1)
	{
		// Failed.
		sprintf(ErrorData, "%s: Failed to read %s from file '%s'.", Routine, Label, FileName);
		throw ErrorData;
	}

	// Return the new character
	return data;
}


int ReadIntFromFile(
	FILE				*stream,			// (I) Stream
	const char			*Label,				// (I) Label
	const char			*Routine,			// (I) Calling routine
	const char			*FileName)			// (I) Filename
//////////////////////////////////////////////////////////////////////
// Reads an integer from a file.
{
	char				ErrorData[256];		// Error message
	int					data;				// Data to read
	int					readerror;			// How many things read?

	// Read it.
    if ((readerror = fscanf (stream, "%d \n", &data)) != 1)
	{
		// Failed.
		sprintf(ErrorData, "%s: Failed to read %s from file '%s'.", Routine, Label, FileName);
		throw ErrorData;
	}

	// Return the new character
	return data;
}


long ReadLongFromFile(
	FILE				*stream,			// (I) Stream
	const char			*Label,				// (I) Label
	const char			*Routine,			// (I) Calling routine
	const char			*FileName)			// (I) Filename
//////////////////////////////////////////////////////////////////////
// Reads a long from a file.
{
	char				ErrorData[256];		// Error message
	long				data;				// Data to read
	int					readerror;			// How many things read?

	// Read it.
    if ((readerror = fscanf (stream, "%ld \n", &data)) != 1)
	{
		// Failed.
		sprintf(ErrorData, "%s: Failed to read %s from file '%s'.", Routine, Label, FileName);
		throw ErrorData;
	}

	// Return the new character
	return data;
}


double ReadDoubleFromFile(
	FILE				*stream,			// (I) Stream
	const char			*Label,				// (I) Label
	const char			*Routine,			// (I) Calling routine
	const char			*FileName)			// (I) Filename
//////////////////////////////////////////////////////////////////////
// Reads a double character from a file.
{
	char				ErrorData[256];		// Error message
	double				data;				// Data to read
	int					readerror;			// How many things read?

	// Read it.
    if ((readerror = fscanf (stream, "%lf \n", &data)) != 1)
	{
		// Failed.
		sprintf(ErrorData, "%s: Failed to read %s from file '%s'.", Routine, Label, FileName);
		throw ErrorData;
	}

	// Return the new character
	return data;
}












