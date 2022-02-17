/*----------------------------------------------------------------------------*
    bkshift.cpp
*----------------------------------------------------------------------------*/
 
 
 
 
#include "bkshift.h"




/*----------------------------------------------------------------------------*
    Constructor for bucket shift
*----------------------------------------------------------------------------*/
ARM_BucketShift::ARM_BucketShift(double bpShift, double bucketStartPeriod,
                       double bucketEndPeriod)
{
    Init();
 
    SetName(ARM_BUCKET_SHIFT);
 
    itsBPShift = bpShift;
 
    itsBucketStartPeriod = bucketStartPeriod;
 
    itsBucketEndPeriod = bucketEndPeriod;
}






/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
