/*
 * $Log: ia3levrefval.cpp,v $
 * Revision 1.4  1998/12/02 13:52:53  nicolasm
 * Das le constructeur itsAmortRate0 etait affecté 3 fois au lieu de
 * itsAmortRate1 et itsAmortRate2
 *
 * Revision 1.3  1998/12/01 08:37:05  nicolasm
 * Getion des taux d'amortissement sous forme de pourcentage
 * plutot que decimal
 *
 * Revision 1.2  1998/11/30 15:35:26  nicolasm
 * Ajout du SetName dans le constructeur
 *
 * Revision 1.1  1998/11/27 14:58:06  nicolasm
 * Initial revision
 *
 */


#include "ia3levrefval.h"
#include "interpol.h"

ARM_IA3LevRefVal::ARM_IA3LevRefVal(double nomValue,
                                   double level0, double amort0, double level1, 
                                   double amort1, double level2, double amort2):
                            ARM_IARefVal(nomValue)

{
    SetName(ARM_IAREFVAL);

    itsLevel0 = level0;
    itsAmortRate0 = amort0;

    itsLevel1 = level1;
    itsAmortRate1 = amort1;

    itsLevel2 = level2;
    itsAmortRate2 = amort2;
}


double ARM_IA3LevRefVal::AmortNominal(double refIndex, double previousNominal)
{
    double amortRate;

    if (refIndex < itsLevel0)
    {
        amortRate = 0.0;
    }
    else if (itsLevel0 <= refIndex && refIndex < itsLevel1)
    {
        amortRate = linInterpol(refIndex, itsLevel0, itsAmortRate0, 
                                                     itsLevel1, itsAmortRate1);

    }
    else if (itsLevel1 <= refIndex && refIndex < itsLevel2)
    {
        amortRate = linInterpol(refIndex, itsLevel1, itsAmortRate1, 
                                                     itsLevel2, itsAmortRate2);
    }
    else
    {
        amortRate = itsAmortRate2;
    }

    double newNominal = (100.0 - amortRate)*previousNominal/100.0;

    return newNominal;

}
