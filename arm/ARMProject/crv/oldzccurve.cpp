
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : oldzccurve.cpp                                                   */
/*                                                                            */
/* DESCRIPTION : Header for the ARM_OldZcCurve class,							  */ 
/*               A curve keeping the discount from a day to another.		  */
/*                                                                            */
/* DATE        : Tue April  24 2007                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/


/*---- Application Include ----*/

#include "firsttoinc.h"
#include "oldzccurve.h"

ARM_OldZcCurve::ARM_OldZcCurve(ARM_Date& asOfDate, ARM_ZeroCurve* crv)
:
ARM_ZeroCurve(),
itsZeroCurve(static_cast<ARM_ZeroCurve*>(crv->Clone())),
itsAsOfDate(asOfDate)

{
	Init();
	SetCurrencyUnit(itsZeroCurve->GetCurrencyUnit());
	SetAsOfDate(asOfDate);
}
        
ARM_OldZcCurve::ARM_OldZcCurve(const ARM_OldZcCurve& rhs)
:
ARM_ZeroCurve(rhs),
itsZeroCurve(static_cast<ARM_ZeroCurve*>(rhs.itsZeroCurve->Clone())),
itsAsOfDate(rhs.itsAsOfDate)
{
}

ARM_OldZcCurve& ARM_OldZcCurve::operator=(const ARM_OldZcCurve& rhs)
{
	if (&rhs != this)
	{	this->~ARM_OldZcCurve ();
		new (this) ARM_OldZcCurve(rhs);
	}
	return *this;
}

ARM_Object* ARM_OldZcCurve::Clone(void)
{
	return new ARM_OldZcCurve(*this);
}

void ARM_OldZcCurve::View(char* id, FILE* ficOut)
{
	FILE* fOut;
    char fOutName[200];
 
    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

	fprintf(fOut, "\n xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx \n");
    fprintf(fOut, "\n\n\t          >>>>>>>>>>>>> Old Curve <<<<<<<<<<<<\n\n");
	fprintf(fOut, "\n\n New AsOfDate : %s \n", (const char *) itsAsOfDate.toString('.').c_str());
	fprintf(fOut, "\n\n");

	itsZeroCurve->View(id,fOut);

	if ( ficOut == NULL )
    {
       fclose(fOut);
    }

}

double ARM_OldZcCurve::DiscountFunction(double lag)
{
	ARM_Date prevAsOfDate = itsZeroCurve->GetAsOfDate();

	return itsZeroCurve->DiscountFunction(lag + (itsAsOfDate.GetJulian()-prevAsOfDate.GetJulian())/K_YEAR_LEN);
}


/*--------------------------------------------------------------------------*/
/*---- End of file ----*/
