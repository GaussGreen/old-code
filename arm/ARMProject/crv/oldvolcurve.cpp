
/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : oldvolcurve.cpp                                              */
/*                                                                            */
/* DESCRIPTION : Header for the ARM_OldCurve class,							  */ 
/*               A curve keeping the volatility from a day to another.		  */
/*                                                                            */
/* DATE        : Thu April  26 2007                                           */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/


/*---- Application Include ----*/

#include "firsttoinc.h"
#include "oldvolcurve.h"

ARM_OldVolCurve::ARM_OldVolCurve(ARM_Date& asOfDate, ARM_VolCurve* crv)
:
ARM_VolCurve(),
itsVolCurve(static_cast<ARM_VolCurve*>(crv->Clone())),
itsAsOfDate(asOfDate)

{
	Init();
	itsVolCurve->SetATMVolFX(NULL);
	SetCurrencyUnit(itsVolCurve->GetCurrency());
	SetAsOfDate(asOfDate);
	SetName(itsVolCurve->GetName());
}
        
ARM_OldVolCurve::ARM_OldVolCurve(const ARM_OldVolCurve& rhs)
:
ARM_VolCurve(rhs),
itsVolCurve(static_cast<ARM_VolCurve*>(rhs.itsVolCurve->Clone())),
itsAsOfDate(rhs.itsAsOfDate)
{
	SetName(itsVolCurve->GetName());
}

ARM_OldVolCurve& ARM_OldVolCurve::operator=(const ARM_OldVolCurve& rhs)
{
	if (&rhs != this)
	{	this->~ARM_OldVolCurve ();
		new (this) ARM_OldVolCurve(rhs);
	}
	return *this;
}

ARM_Object* ARM_OldVolCurve::Clone(void)
{
	return new ARM_OldVolCurve(*this);
}

void ARM_OldVolCurve::View(char* id, FILE* ficOut)
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
    fprintf(fOut, "\n\n\t          >>>>>>>>>>>>> Old Vol Curve <<<<<<<<<<<<\n\n");
	fprintf(fOut, "\n\n New AsOfDate : %s \n", (const char *) itsAsOfDate.toString('.').c_str());
	fprintf(fOut, "\n\n");

	itsVolCurve->View(id,fOut);

	if ( ficOut == NULL )
    {
       fclose(fOut);
    }

}

double ARM_OldVolCurve::computeVol(double moneyness, double maturity)
{
	ARM_Date prevAsOfDate = itsVolCurve->GetAsOfDate();

	return itsVolCurve->computeVol(moneyness, maturity + (itsAsOfDate.GetJulian()-prevAsOfDate.GetJulian())/K_YEAR_LEN);
}

double ARM_OldVolCurve::VolatilityFunction(double m1, double m2)
{
	ARM_Date prevAsOfDate = itsVolCurve->GetAsOfDate();

	return itsVolCurve->VolatilityFunction(m1 + (itsAsOfDate.GetJulian()-prevAsOfDate.GetJulian())/K_YEAR_LEN,m2);
}

ARM_VolLInterpol* ARM_OldVolCurve::ComputeATMFxVol(double lag) 
{
	ARM_Date prevAsOfDate = itsVolCurve->GetAsOfDate();
	return itsVolCurve->ComputeATMFxVol((itsAsOfDate.GetJulian()-prevAsOfDate.GetJulian())/K_YEAR_LEN);
}

int ARM_OldVolCurve::GetSmileFlag(void) const  
{ 
	return itsVolCurve->GetSmileFlag();   
}

ARM_Vector* ARM_OldVolCurve::GetExpiryTerms(void) const
{
	return itsVolCurve->GetExpiryTerms();
}


/*--------------------------------------------------------------------------*/
/*---- End of file ----*/
