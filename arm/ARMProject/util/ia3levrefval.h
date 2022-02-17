/*
 * $Log: ia3levrefval.h,v $
 * Revision 1.1  1998/11/27 14:58:49  nicolasm
 * Initial revision
 *
 */

/*----------------------------------------------------------------------------*

    ia3levrefval.h
 
    Header for the ARM_IIRefVal class, which implements
    calculation of reference values which depend on te value of an
    index
    
   

*----------------------------------------------------------------------------*/
#ifndef _IA3LREFVAL_H
#define _IA3LREFVAL_H





#include "iarefval.h"




/*----------------------------------------------------------------------------*/



class ARM_IA3LevRefVal : public ARM_IARefVal 
{
private:
    double itsLevel0;
    double itsAmortRate0;

    double itsLevel1;
    double itsAmortRate1;

    double itsLevel2;
    double itsAmortRate2;

public:
    ARM_IA3LevRefVal(void)
    {
        itsLevel0 = 0.0;
        itsAmortRate0 = 0.0;

        itsLevel1 = 0.0;
        itsAmortRate1 = 0.0;

        itsLevel2 = 0.0;
        itsAmortRate2 = 0.0;
    }

    ARM_IA3LevRefVal(double nominValue, double level0, double amort0, 
                     double level1, double amort1,
                     double level2, double amort2);

    ~ARM_IA3LevRefVal(void) {}

    double AmortNominal(double ref_index, double previous_nominal);

    double Nominal(void) {return CptReferenceValue((ARM_Date) "01.01.1990"); }



    // ***** Services de Copie *****

    void BitwiseCopy(const ARM_Object* srcRef)
    {
       ARM_IA3LevRefVal* ref = (ARM_IA3LevRefVal *) srcRef;

       itsLevel0 = ref->itsLevel0;
       itsAmortRate0 = ref->itsAmortRate0;

       itsLevel1 = ref->itsLevel1;
       itsAmortRate1 = ref->itsAmortRate1;

       itsLevel2 = ref->itsLevel2;
       itsAmortRate2 = ref->itsAmortRate2;
    }


    void Copy(const ARM_Object* srcRef)
    {
        ARM_IARefVal::Copy(srcRef);

        this->BitwiseCopy(srcRef);
    }


    ARM_Object* Clone(void)
    {
        ARM_IA3LevRefVal* theClone = new ARM_IA3LevRefVal();
 

        theClone->Copy(this);
  
        return(theClone);
    }

};

#endif
