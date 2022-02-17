/**************************************************************************
 * Module to support reading/writing of counted-array analytics input files.
 *************************************************************************/

#ifdef WIN_NT
#include <windows.h>
#else
#include <errno.h>
#endif
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#include <macros.h>
#include <cerror.h>
#include <drlcarrio.h>

static char linebuf[GTO_MAX_STR_LEN + sizeof("= \n")];


long drlReadArraySize(FILE* fp, const char* name)
{ 
  long size, c;
  char tag[256];

  sprintf(tag, "%s[%%ld]\n", name);

  while ((c = getc(fp)) != EOF)
    if (!isspace(c))
    {
      ungetc(c, fp);
      break;
    } 

  if (fgets(linebuf, sizeof(linebuf), fp) != NULL &&
      sscanf(linebuf, tag, &size) == 1)
    return size;
  else
  {
    GtoErrMsg("Could not read %s.\n", name); 
    return -1;
  }
}


int drlWriteArraySize(FILE* fp, const char* name, long arraySize)
{
  return fprintf(fp, "%s[%ld]\n", name, arraySize) > 0 ? 1 : 0; 
}


int drlReadDoubleArray(FILE* fp, const char* name, double** array)
{
  long size;
  long i = 0;
  double aValue;

  if ((size = drlReadArraySize(fp, name)) < 0)
    return 0;

  *array = (double*) malloc((size + 1) * sizeof(double));
  (*array)[0] = (double)size;

  while (i < size &&
	 fgets(linebuf, sizeof(linebuf), fp) != NULL &&
         sscanf(linebuf, "= %lf\n", &aValue) == 1)
    (*array)[++i] = aValue;

  if (i == size)
    return 1;
  else
  {
    GtoErrMsg("Not enough values for array: %s.\n", name);
    return 0;
  }
}


int drlWriteDoubleArray(FILE* fp, const char* name, double array[])
{
  long arraySize = (long)array[0];
  long i;

  if (drlWriteArraySize(fp, name, arraySize))
  {
    for (i = 1; i <= arraySize; ++i)
      /* Note that writing out 17 decimal digits guarantees no loss of precision
         when converting an IEEE double from binary to decimal and back again. */ 
      if (fprintf(fp, "= %.17lf\n", array[i]) < 0)
        return 0;
  
    return 1;
  }

  return 0;
}


int drlReadLongArray(FILE* fp, const char* name, long** array)
{
  long size;
  long i = 0;
  long aValue;

  if ((size = drlReadArraySize(fp, name)) < 0)
    return 0;

  *array = (long*) malloc((size + 1) * sizeof(long));
  (*array)[0] = size;

  while (i < size &&
	 fgets(linebuf, sizeof(linebuf), fp) != NULL &&
	 sscanf(linebuf, "= %d\n", &aValue) == 1)
    (*array)[++i] = aValue;

  if (i == size)
    return 1;
  else
  {
    GtoErrMsg("Not enough values for array: %s.\n", name);
    return 0;
  }
}


int drlWriteLongArray(FILE* fp, const char* name, long array[])
{
  long arraySize = (long)array[0];
  long i;

  if (drlWriteArraySize(fp, name, arraySize))
  {
    for (i = 1; i <= arraySize; ++i)
      if (fprintf(fp, "= %ld\n", array[i]) < 0)
        return 0;
  
    return 1;
  }

  return 0;
}


int drlReadStringArray(FILE* fp, const char* name, char*** array)
{
  long size, strLen;
  long i = 1;
  char* aValue = NULL;

  if ((size = drlReadArraySize(fp, name)) < 0)
  {
    *array = NULL;
    return 0;
  }

  *array = (char**)calloc(size + 1, sizeof(char*));
  if (*array == NULL)
    return 0;

  (*array)[0] = (char*)size;

  while (i <= size && fgets(linebuf, sizeof(linebuf), fp) != NULL)
    if (strncmp(linebuf, "= ", 2) == 0)
    {
      strLen = strcspn(linebuf + 2, "\n");
      if (strLen > GTO_MAX_STR_LEN)
      {
        GtoErrMsg("String element %d for %s too long (len = %d, max = %d).\n",
		  i, name, strLen, GTO_MAX_STR_LEN);
        return 0;
      }

      (*array)[i] = (char*)calloc(strLen + 1, sizeof(char));
      if ((*array)[i] == NULL)
        return 0;

      strncpy((*array)[i], linebuf + 2, strLen);
      ++i;
    }

  if (i == size + 1)
    return 1;
  else
  {
    GtoErrMsg("Not enough values for array: %s.\n", name);
    return 0;
  }
}


int drlWriteStringArray(FILE* fp, const char* name, char* stringArray[])
{
  long arraySize = (long)stringArray[0];
  long i;
  const char* s;

  if (drlWriteArraySize(fp, name, arraySize))
  {
    for (i = 1; i <= arraySize; ++i)
    {
      if (strlen(s = stringArray[i]) <= GTO_MAX_STR_LEN)
      {
        if (fprintf(fp, "= %s\n", s) < 0)
          return 0;
      }
      else
      {
        GtoErrMsg("String element %d for %s too long (len = %d, max = %d).\n",
		  i, name, strlen(s), GTO_MAX_STR_LEN);
	return 0;
      }
    }
  
    return 1;
  }

  return 0;
}
