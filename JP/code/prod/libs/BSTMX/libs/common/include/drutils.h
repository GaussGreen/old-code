#ifndef _druril_h
#define _druril_h

#include <stdlib.h>
#include "drtypes.h"

#ifdef  __cplusplus
extern "C" {
#endif

/** Translate KDcc enum to DR character constant.

    @param dcc:  (I) EDcc enum

    @return          DR character constant
*/
char 
KDccToDRDcc(KDcc dcc);


/** Translate KFrequency enum to DR character constant.

    @param dcc:  (I) EFrequency enum

    @return          DR character constant
*/
char 
KFrequencyToDRFrequency(KFrequency frq);




/** Return a new copy of an array 

    @param src:  (I) pointer to the source data
    @param size: (I) size of the source data expressed
                     as number of elements times element size

    @return          pointer to the newly allocated array
*/
void* ArrayNewCopy(void const*  src,
                   size_t size);




/** Free heap allocated storage. Chacks pointer for null and 
    calls free.

    @param src:  (I) pointer to the data to be freed
*/
void SafeFree(void* ptr);

#ifdef  __cplusplus
}
#endif


#endif
