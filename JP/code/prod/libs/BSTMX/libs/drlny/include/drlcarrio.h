/**************************************************************************
 * Module to support reading/writing of counted-array analytics input files.
 *************************************************************************/

#ifndef _drlcarrio_h
#define _drlcarrio_h

#include <stdio.h>

extern long drlReadArraySize(FILE* fp, const char* name);
extern int drlWriteArraySize(FILE* fp, const char* name, long size);
extern int drlReadDoubleArray(FILE* fp, const char* name, double** array);
extern int drlWriteDoubleArray(FILE* fp, const char* name, double array[]);
extern int drlReadLongArray(FILE* fp, const char* name, long** array);
extern int drlWriteLongArray(FILE* fp, const char* name, long array[]);
extern int drlReadStringArray(FILE* fp, const char* name, char*** array);
extern int drlWriteStringArray(FILE* fp, const char* name, char* array[]);

#endif /* _drlcarrio_h */
