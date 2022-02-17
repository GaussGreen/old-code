/*
 * $Log: zerointspreaded.h,v $
 * Revision 1.0  2006/08/09 16:33:00 mb
 * "Constification"
 *
 */

/*----------------------------------------------------------------------------*
 
    zerointspreaded.cpp
 
    This file implements the ARM_BasisCurve class, a class for 
         computing a ARM_ZeroCurve with spreads adjustments

*----------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "util.h"

#include "zerocurv.h"
#include "zeroint.h"
#include "zerointspreaded.h"
#include "currency.h"



/*----------------------------------------------------------------------------*/



ARM_Object* ARM_BasisCurve::Clone(void)
{
    ARM_BasisCurve* theClone = new ARM_BasisCurve();
            
    theClone->Copy(this);
            
    return(theClone);
}

 

ARM_BasisCurve& ARM_BasisCurve::operator = (const ARM_BasisCurve& zccurve)
{
	(*this).ARM_ZeroLInterpol::operator = (zccurve);
			
    BitwiseCopy(& zccurve);
			
    return(*this);
}



void ARM_BasisCurve::BitwiseCopy(const ARM_Object* srczintsp)
{
	Init();
    
    ARM_BasisCurve* zintsp = (ARM_BasisCurve*	) srczintsp;
	
    if (itsBasisCurve)	
    {	
        delete itsBasisCurve;	
    
        itsBasisCurve = NULL;	
    } 
	
    if (itsZeroCurve)	
    {	
        delete itsZeroCurve;	
        itsZeroCurve  = NULL;	
    } 

    if (zintsp->GetBasisCurve())
       itsBasisCurve = static_cast<ARM_ZeroCurve* >(zintsp->GetBasisCurve()->Clone());

	if (zintsp->GetZeroCurve())
       itsZeroCurve = static_cast<ARM_ZeroCurve* >(zintsp->GetZeroCurve()->Clone());

    itsMMFreq   = zintsp->itsMMFreq;

    itsSwapFreq = zintsp->itsSwapFreq;
}



void ARM_BasisCurve::Copy(const ARM_Object* srczintsp)
{
     ARM_ZeroLInterpol::Copy(srczintsp);
            
     BitwiseCopy(srczintsp);
}



void ARM_BasisCurve::Init() 
{	
	SetName(ARM_BASIS_CURVE);
		
    itsBasisCurve	= NULL;	
		
    itsZeroCurve	= NULL;
    
    itsMMFreq       = 2;

    itsSwapFreq     = 1;
}
   


void ARM_BasisCurve::GenerateBasisCurve(ARM_ZeroCurve* zcurve )	
{

	ARM_Currency*	ccy		 = itsZeroCurve->GetCurrencyUnit();
	ARM_Date		asofdate = itsZeroCurve->GetAsOfDate();
	

	ARM_BasisCurve tmp(asofdate,	
                       zcurve,	
                       itsZeroCurve, 
                       itsMMFreq, 	
                       itsSwapFreq,	
                       ccy);
	
    Copy(&tmp);	
}



void ARM_BasisCurve::GenerateZeroCurve(ARM_ZeroCurve* zcurve )	
{

	ARM_Currency*	ccy		 =	itsZeroCurve->GetCurrencyUnit();
	ARM_Date		asofdate =	itsZeroCurve->GetAsOfDate();

	ARM_BasisCurve tmp(asofdate, 
                       itsBasisCurve, 
                       zcurve, 
                       itsMMFreq, 	
                       itsSwapFreq,	
                       ccy);
	Copy(&tmp);	
}



ARM_BasisCurve*	ARM_BasisCurve::GenerateShiftCurve(ARM_CRV_TERMS& Term, ARM_Vector* epsilon, const char* curveToBump)
{
	ARM_BasisCurve*	vNewSpreadedCurve = NULL;
	ARM_ZeroCurve*	vBasisCurve       = NULL;
	ARM_ZeroCurve*	vZCCurve          = NULL;


	if ( strcmp(curveToBump, "BS") == 0 )
	{
	   vBasisCurve = itsBasisCurve->GenerateShiftCurve(Term, epsilon);
		
       vZCCurve = (ARM_ZeroCurve*) itsZeroCurve->Clone();
	}
	else // ZC
	{
	   vZCCurve = itsZeroCurve->GenerateShiftCurve(Term, epsilon);
	   
       vBasisCurve = (ARM_ZeroCurve*) itsBasisCurve->Clone();
	}

	vNewSpreadedCurve = new ARM_BasisCurve(itsZeroCurve->GetAsOfDate(), 
										   vBasisCurve, 
										   vZCCurve, 
										   itsMMFreq, 	
										   itsSwapFreq,	
										   itsZeroCurve->GetCurrencyUnit());

	delete vBasisCurve;
    delete vZCCurve;

    return(vNewSpreadedCurve);
}



ARM_BasisCurve*	ARM_BasisCurve::GenerateShiftCurveFwd(ARM_CRV_TERMS& Term, ARM_Vector* epsilon, const char* curveToBump)
{
	ARM_BasisCurve*	vNewSpreadedCurve = NULL;
	ARM_ZeroCurve*	vBasisCurve = NULL;
	ARM_ZeroCurve*	vZCCurve = NULL;

	if(strcmp(curveToBump, "BS") == 0)
	{
		vBasisCurve = itsBasisCurve->GenerateShiftCurveFwd(Term, epsilon);
		vZCCurve = (ARM_ZeroCurve*) itsZeroCurve->Clone();
	}
	else // ZC
	{
		vZCCurve = itsZeroCurve->GenerateShiftCurveFwd(Term, epsilon);
		vBasisCurve = (ARM_ZeroCurve*) itsBasisCurve->Clone();
	}

	vNewSpreadedCurve = new ARM_BasisCurve(itsZeroCurve->GetAsOfDate(), 
										   vBasisCurve, 
										   vZCCurve, 
										   itsMMFreq, 	
										   itsSwapFreq,	
										   itsZeroCurve->GetCurrencyUnit());

	return	vNewSpreadedCurve;
}



void ARM_BasisCurve::View(char* id, FILE* ficOut)
{	
    FILE* fOut;
    char  fOutName[200];
    
    
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

    fprintf(fOut, "\n ==============> ARM_BasisCurve <==============\n\n\n");

    fprintf(fOut, "\n Money Market Frequency: %d\n", itsMMFreq);
    fprintf(fOut, "\n Swap Frequency        : %d\n\n", itsSwapFreq);

    fprintf(fOut, "\n -------------------------------------------> The Initial ZC Curve:\n\n");

    itsZeroCurve->View(id, fOut);
			

    fprintf(fOut, "\n\n\n -------------------------------------------> The Spreads Curve:\n\n");

    itsBasisCurve->View(id, fOut);


    fprintf(fOut, "\n\n\n        >>>>>>>>>> The Generated basis adjusted ZC Curve <<<<<<<<<<\n\n");

    ARM_ZeroLInterpol::View(id, fOut);

    if ( ficOut == NULL )
	    fclose(fOut);
}


/*------------------------------------------------------------------------------------------*/
/*---- End Of File ----*/