/*-----------------------------------------------------------------------
  SOURCE FILE:    errbuff.h
  
  PURPOSE:        Redirect error to buffer
  ---------------------------------------------------------------------- */
#ifndef DR_ERRBUFF_H
#define DR_ERRBUFF_H

#include "cgeneral.h"       /* GTO_EXPORT */

/* Internal buffer for error logging
 */
#define   DRL_MAX_BUFFER_LENGTH   8136L 


/* Initialize error logging into internal buffer
 */
GTO_EXPORT(void) drlErrorLogInit (void);


/* Retrieves error messages that are written to the 
 * global buffer by the error message call back function.
 * 
 * Returns FAILURE if outBufLen is less than the length
 * of the accumulated error message string
 */ 

GTO_EXPORT(int)   drlRetrieveErrorMessages 
(long    outBufLen,    /* (I) length of the input buffer. */ 
 char  * outBuffer     /* (I/O) output buffer. */ 
);

#endif



