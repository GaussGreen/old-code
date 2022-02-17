/*
***************************************************************************
** convert.h
**
** Various useful conversions routines which don't really have anywhere
** else to go.
***************************************************************************
*/
#ifndef IRX_CONVERT_H
#define IRX_CONVERT_H

#include "irxutils.h"

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * Formats a date. Can be called eight times from the same print
 * statement, but not more. Format is DD-MMM-YYYY.
 *
 * This function is not thread safe!
 */
extern char * irxFormatDate(IrxTDate date);

/**
 * Format a date interval.
 *
 * This function is not thread safe!
 */
extern char* irxFormatDateInterval(IrxTDateInterval interval);

/**
***************************************************************************
** Extracts an array from a structure array.
**
** For example, suppose you have a structure X which contains a date and
** a double, e.g. TDate aDate and double aDouble.
**
** Then given an array of structure X (not of pointers to structure X),
** then to extract the array of dates you would do something like this:
**
** size_t arraySize;
** X structArray[arraySize];
** TDate *dateArray;
** double *doubleArray;
**
** dateArray = (TDate*) irxArrayFromStructArray (sizeof(TDate),
**                                               sizeof(X),
**                                               &(structArray[0].aDate),
**                                               arraySize);
**
** doubleArray = (double*) irxArrayFromStructArray (sizeof(double),
**                                                  sizeof(X),
**                                                  &(structArray[0].aDouble),
**                                                  arraySize);
***************************************************************************
*/
void* irxArrayFromStructArray
(size_t elemSize,    /* (I) Size of elements of returned array: use sizeof */
 size_t structSize,  /* (I) Size of structure: use sizeof                  */
 void  *base,        /* (I) Address of field in first item in struct array */
 size_t numElems     /* (I) Size of the array                              */
);

#define EXTRACT_ARRAY(fieldType, structType, array, fieldName, arraySize)\
    (array == NULL) ? NULL : \
    (fieldType*)(irxArrayFromStructArray(sizeof(fieldType),\
                                         sizeof(structType),\
                                         &(array->fieldName),\
                                         (size_t)(arraySize)))


#ifdef __cplusplus
}
#endif

#endif    /* IRX_CONVERT_H */





