/*
 $Log: calibration.cpp,v $
 Revision 1.4  2003/02/11 17:17:16  mab
 Added : Init()

 Revision 1.3  2002/11/21 16:37:45  mab
 Formatting

 Revision 1.2  2000/10/17 15:58:46  mab
 *** empty log message ***

 Revision 1.1  1999/11/10 08:54:46  nicolasm
 Initial revision

 */

/*----------------------------------------------------------------------------*/




#include "calibration.h"






ARM_Calibration::ARM_Calibration(ARM_Date& asOf, ARM_ZeroCurve* zc)
{
    Init();

    itsAsOf = asOf;
 
    itsZeroCurve = (ARM_ZeroCurve *) zc->Clone();
}



/*----------------------------------------------------------------------------*/
/*---- End Of File ----*/
