/****************************************************************************/
/*      C++ Function templates for library.                                 */
/****************************************************************************/
/*      fix123head.h                                                        */
/****************************************************************************/


#ifndef _fix123plus_h
#define _fix123plus_h




/* Try to find a string str in a list of strings.  Return the */
/* number of the position where the string was found, or -1   */
/* if the string wasn't found.                                */
int 
CompareManyStrings(
	const char			*str,				/* String to find */
	int					nbStrings,			/* Number of strings */
	const char			*stringList[]);		/* Strings to find it in */


/* Read a comment line from file.  On error, throw an exception. */
void FindAndSkipComLine_CPP(
	FILE				*stream,			// (I) Stream
	const char	   		*Label,				// (I) Label
	const char	   		*Routine,			// (I) Calling routine
	const char	   		*FileName);			// (I) Filename


/* Read a char from file.  On error, throw an exception. */
char ReadCharFromFile(
	FILE				*stream,			// (I) Stream
	const char		   	*Label,				// (I) Label
	const char		   	*Routine,			// (I) Calling routine
	const char		   	*FileName);			// (I) Filename


/* Read a int from file.  On error, throw an exception. */
int ReadIntFromFile(
	FILE				*stream,			// (I) Stream
	const char			*Label,				// (I) Label
	const char			*Routine,			// (I) Calling routine
	const char			*FileName);			// (I) Filename


/* Read a long from file.  On error, throw an exception. */
long ReadLongFromFile(
	FILE				*stream,			// (I) Stream
	const char		   	*Label,				// (I) Label
 	const char			*Routine,			// (I) Calling routine
	const char		   	*FileName);			// (I) Filename


/* Read a double from file.  On error, throw an exception. */
double ReadDoubleFromFile(
	FILE				*stream,			// (I) Stream
	const char			*Label,				// (I) Label
	const char			*Routine,			// (I) Calling routine
	const char			*FileName);			// (I) Filename


/* Allocate a tree slice.  On error, throw an exception. */
double *Alloc_Slice_CPP(
	TREE_DATA			*tree_data);			// (I) Tree data


/* Allocate an array.  On error, throw an exception. */
void *DR_Array_CPP(
	int					type,				// (I) Type
    int					nl,					// (I) Lower bound
    int					nh);					// (I) Higher bound



#endif /* _fix123plus_h */
