/*
 * $Log: volflat.cpp,v $
 * Revision 1.5  2004/03/05 16:49:12  mab
 * Use Set Method in constructor
 *
 * Revision 1.4  2004/02/19 13:30:00  mab
 * Added : SetVolatilities
 *
 * Revision 1.3  2003/08/18 09:52:46  mab
 *  improvement in the principal method :
 *  virtual double VolatilityFunction(double m1, double K, double m2)
 *  virtual double VolatilityFunction(double m1, double m2) :
 * No default parameter but 2 methods!
 *
 * Revision 1.2  2002/10/11 08:26:36  mab
 * Added : BumpVolatility
 *
 */

/*----------------------------------------------------------------------------*
 
    volflat.cpp
 
    This file implements the ARM_VolFlat class
*----------------------------------------------------------------------------*/




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "volflat.h"




/*----------------------------------------------------------------------------*/



void ARM_VolFlat::Set(ARM_Date& asOf, double vol)
{
    SetAsOfDate(asOf);

    itsVolatility = vol;

    ARM_Matrix* vols = new ARM_Matrix(2, 2, vol);

    SetVolatilities(vols);

    int expiriesSize = 2;

    ARM_Vector* tmpExp = new ARM_Vector(expiriesSize);

    tmpExp->Elt(0) = 0.0;
    tmpExp->Elt(1) = 1000.0;
   
    SetExpiryTerms(tmpExp);

 
    int tenorsSize = 2;

    ARM_Vector* tmpStrikes = new ARM_Vector(tenorsSize);

    tmpStrikes->Elt(0) = 0.0;
    tmpStrikes->Elt(1) = 1000.0;

    SetStrikes(tmpStrikes);
}



ARM_VolFlat::ARM_VolFlat(ARM_Date& asOf, double vol,
                         ARM_Currency* ccy)
{
    Init();
 
    Set(asOf, vol);

    SetCurrencyUnit(ccy);
}



void ARM_VolFlat::View(char* id, FILE* ficOut)
{
    FILE* fOut;

    char fOutName[200]; 
 
    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
 
       (void) unlink(fOutName);

       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

	fprintf(fOut, "\n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
    fprintf(fOut, "\n =================> ARM Flat Volatility \n\n");

    ARM_VolLInterpol::View(id, fOut);

    fprintf(fOut, "\n\n <================= ARM Flat Volatility \n\n");

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}


/*----------------------------------------------------------------------------*/
/*---- End of file ----*/