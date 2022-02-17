/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE         : smilecrv.cpp                                                */
/*                                                                            */
/* DESCRIPTION  : This file implements the ARM_SmileCurve class, a class for  */
/*                smile curves                                                */
/*                                                                            */
/* DATE         : Tue Apr 1 1997                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*---- System Include ----*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*---- Application Include ----*/

#include "dates.h"
#include "util.h"
#include "volcurv.h"
#include "smilecrv.h"
#include "expt.h"


 
/*----------------------------------------------------------------------------*/
/* Constructor(Copy)                                                          */
/*----------------------------------------------------------------------------*/

ARM_SmileCurve::ARM_SmileCurve(ARM_SmileCurve& volCurve) : ARM_VolCurve(volCurve)
{
    Init();

    BitwiseCopy(&volCurve);
}



/*----------------------------------------------------------------------------*/
/* Assignment operator                                                        */
/*----------------------------------------------------------------------------*/

ARM_SmileCurve& ARM_SmileCurve::operator = (ARM_SmileCurve& volCurve)
{
    (*this).ARM_VolCurve::operator =(volCurve);

    BitwiseCopy(&volCurve);

    return (*this);
}


 


// Rq: m the matruty is obsolete 

double ARM_SmileCurve::VolatilityFunction(double m, double K, double m2)
{
    double resVol;

    m2 = 0.0; // unused

    resVol = linInterpol(itsStrikes, K, itsVolats);

    return(resVol);
}



/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
