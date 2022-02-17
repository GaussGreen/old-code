/*
 * $Log: iarefval.h,v $
 * Revision 1.1  1998/11/27 14:59:37  nicolasm
 * Initial revision
 *
 */

/*----------------------------------------------------------------------------*

    iarefval.h
 
    Header for the ARM_IIRefVal class, which implements
    calculation of reference values which depend on te value of an
    index
    
   

*----------------------------------------------------------------------------*/
#ifndef _IAREFVAL_H
#define _IAREFVAL_H


#include "refvalue.h"




class ARM_IARefVal : public ARM_ReferenceValue 
{
public:
    ARM_IARefVal(void) {}

    ARM_IARefVal(double value) : ARM_ReferenceValue(value)
    {
    }

    ~ARM_IARefVal(void) {}

    virtual double AmortNominal(double ref_index, double previous_nominal);

    virtual double Nominal(void)
    {
        return CptReferenceValue((ARM_Date) "01.01.1990"); 
    }


    // ***** Services de Copie *****

    void BitwiseCopy(const ARM_Object* srcRef)
    {
       ARM_IARefVal* ref = (ARM_IARefVal *) srcRef;
    }
 
    void Copy(const ARM_Object* srcRef)
    {
        ARM_ReferenceValue::Copy(srcRef);

        this->BitwiseCopy(srcRef);
    }

    ARM_Object* Clone(void)
    {
        ARM_IARefVal* theClone = new ARM_IARefVal();
 

        theClone->Copy(this);
  
        return(theClone);
    }

};

#endif
