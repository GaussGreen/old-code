/*
 * $Log: livretAcurv.cpp,v $
 * Revision 1.4  2004/06/09 16:48:42  mab
 * Correction when Fixings exist
 *
 * Revision 1.3  2004/06/09 13:57:19  mab
 * Correction
 *
 * Revision 1.2  2004/06/04 08:31:13  mab
 * Passing dates by ref. plus correction
 *
 */


/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE         : livretAcurv.cpp                                             */
/*                                                                            */
/* DESCRIPTION  : This file implements the ARM_LivretACurve class, a class    */
/*                for managing livret A rate curves                           */
/*                                                                            */
/* DATE         : Tue Mai 4 2004                                              */
/*                                                                            */
/*                                                                            */
/*----------------------------------------------------------------------------*/

#include "gpinflation/livretAcurv.h"
#include "gpinflation/resetmanager.h"

#define MONTH_FOR_AUGUST K_JUNE
#define MONTH_FOR_FEBRUARY K_NOVEMBER

CC_BEGIN_NAMESPACE( ARM )

/*!
 *    Default constructor
 */

ARM_LivretACurve::ARM_LivretACurve(void)
{
    Init();
}

/*!
 *     constructor
 */
ARM_LivretACurve::ARM_LivretACurve(ARM_Date& asOfDate, 
                                   ARM_InfCurv* infCurv, 
                                   ARM_ZeroCurve* euriborCurv,
                                   long flagArrondi, 
								   ARM_ResetManager* ResetManager,
                                   ARM_ReferenceValue* fixingLivretA,
                                   ARM_ReferenceValue* fixingRateEuribor,
								   long monthForAugustId,
								   long monthForFebruaryId) 
                                   : ARM_ZeroCurve(asOfDate)
{
    Init();

    SetFlagRounding(flagArrondi);
	
	ARM_InfCurv infCurvReset(*infCurv) ; 

	if(ResetManager)
		infCurvReset.SetResetManager(ResetManager);
    
	if (monthForAugustId == K_MONTH_DEFAULT)
		monthForAugustId = MONTH_FOR_AUGUST;

	if (monthForFebruaryId == K_MONTH_DEFAULT)
		monthForFebruaryId = MONTH_FOR_FEBRUARY;

	BuildCurve(&infCurvReset, euriborCurv, fixingLivretA, fixingRateEuribor, monthForAugustId, monthForFebruaryId);
}


/*!
 *     constructor with Euribor ResetManager
 */
ARM_LivretACurve::ARM_LivretACurve(ARM_Date& asOfDate, 
                                   ARM_InfCurv* infCurv, 
                                   ARM_ZeroCurve* euriborCurv,
                                   long flagArrondi, 
								   ARM_ResetManager* InflResetManager,
                                   ARM_ResetManager* LivretAResetManager,
                                   ARM_ResetManager* EuriborResetManager,
								   long monthForAugustId,
								   long monthForFebruaryId) 
                                   : ARM_ZeroCurve(asOfDate)
{
    Init();

    SetFlagRounding(flagArrondi);
	
	ARM_InfCurv infCurvReset(*infCurv) ; 

	if(InflResetManager)
		infCurvReset.SetResetManager(InflResetManager);
    
	if (monthForAugustId == K_MONTH_DEFAULT)
		monthForAugustId = MONTH_FOR_AUGUST;

	if (monthForFebruaryId == K_MONTH_DEFAULT)
		monthForFebruaryId = MONTH_FOR_FEBRUARY;

	BuildCurve(&infCurvReset, euriborCurv, LivretAResetManager, EuriborResetManager, monthForAugustId, monthForFebruaryId);
}


/*!
 *    Init function
 */
void ARM_LivretACurve::Init(void)
{

    SetName(ARM_LIVRETACURVE);

    SetFlagRounding(0); // no rounding
}


/*!
 *    Build Curve
 */
void ARM_LivretACurve::BuildCurve(ARM_InfCurv* infCurv, 
                                  ARM_ZeroCurve* euriborCurv, 
                                  ARM_ReferenceValue* fixingLivretA,
                                  ARM_ReferenceValue* fixingEuribor,
								  long monthForAugustId,
								  long monthForFebruaryId)
{    
    SetCurrencyUnit(euriborCurv->GetCurrencyUnit());

    // get min(end of inflation curve, end of euribor curve)
    double endInfCurve = infCurv->GetBucketEndPeriod();
    
    ARM_Vector* zcTerms = euriborCurv->GetDateTerms();
    
    int sizezc = zcTerms->GetSize();
    
    double endEurCurve = zcTerms->Elt(sizezc-1)+GetAsOfDateJul();

    double endEch = ( endInfCurve < endEurCurve ) 
                    ? endInfCurve : endEurCurve;

    double beginCalEch = GetBeginCalculSchedule();

    int nb = 0;  

    ARM_Date tmpDate = (ARM_Date) beginCalEch;

    // get the nb of the points to calculate for the curve livretA
    while ( tmpDate.GetJulian() < endEch )
    {
        nb++;

        tmpDate.AddMonths(6);
    }

    // get the nb of the correct points fixing 
    int sizelv = 0;

    if ( fixingLivretA != NULL )
    {
        for (int j = 0; j < fixingLivretA->GetSize(); j++)
        {
            if (fixingLivretA->GetDiscreteDates()->Elt(j) > beginCalEch)
                break;
            else 
               sizelv++;
        }

        nb += sizelv;
    }
    
    // build the curve
    ARM_Vector* rateValues = new ARM_Vector(nb);    
    ARM_Vector* dateTerms  = new ARM_Vector(nb);    
    ARM_Vector* yearTerms  = new ARM_Vector(nb);
	ARM_Vector* distFacts  = new ARM_Vector(nb);

	ARM_Date resetDate;//corriger 

    for (int i = 0; i < sizelv; i++)
    {
		resetDate = ARM_Date(fixingLivretA->GetDiscreteDates()->Elt(i));
		resetDate.AddMonths(-1);
		resetDate.AddDays(14);

		dateTerms->Elt(i)  = resetDate.GetJulian()-GetAsOfDate().GetJulian();

        rateValues->Elt(i) = fixingLivretA->GetDiscreteValues()->Elt(i);
        
        yearTerms->Elt(i)  = dateTerms->Elt(i)/K_YEAR_LEN;

		distFacts->Elt(i)  = exp(- rateValues->Elt(i)/100.00 * yearTerms->Elt(i));
    }

    for (i = sizelv; i < nb; i++)
    {
        tmpDate = ARM_Date(beginCalEch).AddMonths(6*(i-sizelv));
        
        ARM_Date dateBeginEuribor = tmpDate;

		if (dateBeginEuribor.GetMonth() == K_AUGUST)
			dateBeginEuribor.AddMonths(monthForAugustId - K_AUGUST);
		else
		{
			if (monthForFebruaryId > K_FEBRUARY)
				dateBeginEuribor.AddMonths(monthForFebruaryId - K_FEBRUARY - 12);
			else
				dateBeginEuribor.AddMonths(monthForFebruaryId - K_FEBRUARY);
		}
        
        rateValues->Elt(i)= RateLivretA(infCurv, euriborCurv, 
                                        dateBeginEuribor, fixingEuribor);

		resetDate = tmpDate;
		resetDate.AddMonths(-1);
		resetDate.AddDays(14);


        dateTerms->Elt(i) = resetDate.GetJulian()- GetAsOfDate().GetJulian();
        
        yearTerms->Elt(i) = dateTerms->Elt(i)/K_YEAR_LEN;

		distFacts->Elt(i)  = exp(- rateValues->Elt(i)/100.00 * yearTerms->Elt(i));
    }
    
    SetBucketEndPeriod(dateTerms->Elt(nb-1)+GetAsOfDate().GetJulian());

    SetBucketStartPeriod( dateTerms->Elt(0)+GetAsOfDate().GetJulian());
    
	SetZeroRates(rateValues);
    
    SetDateTerms(dateTerms);  
    
    SetYearTerms(yearTerms); 

	SetDiscountFactors(distFacts);
}


/*!
 *    Build Curve with Euribor ResetManager
 */
void ARM_LivretACurve::BuildCurve(ARM_InfCurv* infCurv, 
                                  ARM_ZeroCurve* euriborCurv, 
                                  ARM_ResetManager* LivretAResetManager,
                                  ARM_ResetManager* EuriborResetManager,
								  long monthForAugustId,
								  long monthForFebruaryId)
{    
    SetCurrencyUnit(euriborCurv->GetCurrencyUnit());

    // get min(end of inflation curve, end of euribor curve)
    double endInfCurve = infCurv->GetBucketEndPeriod();
    
    ARM_Vector* zcTerms = euriborCurv->GetDateTerms();
    
    int sizezc = zcTerms->GetSize();
    
    double endEurCurve = zcTerms->Elt(sizezc-1)+GetAsOfDateJul();

    double endEch = ( endInfCurve < endEurCurve ) 
                    ? endInfCurve : endEurCurve;

    double beginCalEch = GetBeginCalculSchedule();

    int nb = 0;  

    ARM_Date tmpDate = (ARM_Date) beginCalEch;

    // get the nb of the points to calculate for the curve livretA
    while ( tmpDate.GetJulian() < endEch )
    {
        nb++;

        tmpDate.AddMonths(6);
    }

    // get the nb of the correct points fixing 
    int sizelv = 0;
	OneResetHistory Reset = LivretAResetManager->GetResets();
	doubleDoubleMap listReset = Reset.GetResets();

    if ( LivretAResetManager != NULL )
    {
		doubleDoubleMap::const_iterator iter = listReset.begin();

        while (iter != listReset.end())
        {
            if (iter->first > beginCalEch)
                break;
            else 
               sizelv++;

			iter++;
        }

        nb += sizelv;
    }

    // build the curve
    ARM_Vector* rateValues = new ARM_Vector(nb);    
    ARM_Vector* dateTerms  = new ARM_Vector(nb);    
    ARM_Vector* yearTerms  = new ARM_Vector(nb);
	ARM_Vector* distFacts  = new ARM_Vector(nb);

	ARM_Date resetDate;//corriger 

	doubleDoubleMap::const_iterator iter = listReset.begin();
    for (int i = 0; i < sizelv; i++)
    {
		resetDate = ARM_Date(iter->first);
		resetDate.AddMonths(-1);
		resetDate.AddDays(14);

		dateTerms->Elt(i)  = resetDate.GetJulian()-GetAsOfDate().GetJulian();

        rateValues->Elt(i) = iter->second;
        
        yearTerms->Elt(i)  = dateTerms->Elt(i)/K_YEAR_LEN;

		distFacts->Elt(i)  = exp(- rateValues->Elt(i)/100.00 * yearTerms->Elt(i));

		iter++;
    }

    for (i = sizelv; i < nb; i++)
    {
        tmpDate = ARM_Date(beginCalEch).AddMonths(6*(i-sizelv));
        
        ARM_Date dateBeginEuribor = tmpDate;

		if (dateBeginEuribor.GetMonth() == K_AUGUST)
			dateBeginEuribor.AddMonths(monthForAugustId - K_AUGUST);
		else
		{
			if (monthForFebruaryId > K_FEBRUARY)
				dateBeginEuribor.AddMonths(monthForFebruaryId - K_FEBRUARY - 12);
			else
				dateBeginEuribor.AddMonths(monthForFebruaryId - K_FEBRUARY);
		}
        
        rateValues->Elt(i)= RateLivretA(infCurv, euriborCurv, 
                                        dateBeginEuribor, EuriborResetManager);

		resetDate = tmpDate;
		resetDate.AddMonths(-1);
		resetDate.AddDays(14);

        dateTerms->Elt(i) = resetDate.GetJulian()- GetAsOfDate().GetJulian();
        
        yearTerms->Elt(i) = dateTerms->Elt(i)/K_YEAR_LEN;

		distFacts->Elt(i)  = exp(- rateValues->Elt(i)/100.00 * yearTerms->Elt(i));
    }
    
    SetBucketEndPeriod(dateTerms->Elt(nb-1)+GetAsOfDate().GetJulian());

    SetBucketStartPeriod( dateTerms->Elt(0)+GetAsOfDate().GetJulian());
    
	SetZeroRates(rateValues);
    
    SetDateTerms(dateTerms);  
    
    SetYearTerms(yearTerms); 

	SetDiscountFactors(distFacts);
}



/*!
 *    Calcul Taux LivretA to build the curve
 */

// "dateIn" are the first days of the months;
// formule : eg.
// for the period from 1st August 2004 to 1st Febrary 2005
// Rate = {average of the Euribor3M from 1st mars to 31 mars 2004 
//           + (CPI on 1st mars 2004 / CPI on 1st mars 2003 -1)} / 2 + 0.25%    
double ARM_LivretACurve::RateLivretA(ARM_InfCurv* infCurv, 
                                     ARM_ZeroCurve* euriborCurv, 
                                     ARM_Date& dateIn, 
                                     ARM_ReferenceValue* fixingEuribor)
{    
    double rateEuribor = AverageZCRate(euriborCurv, dateIn, fixingEuribor);
    
    double ratioCPI = LivretACPIRatio(dateIn, infCurv);
    
    double rateLivretA = (rateEuribor+(ratioCPI-1)*100)/2.00 +0.25;
    
    if ( itsFlagRounding == 1 )
    {
       rateLivretA = Rounding(rateLivretA);
    }

    return(rateLivretA);
}



/*!
 *    Calcul Taux LivretA to build the curve
 */

// "dateIn" are the first days of the months;
// formule : eg.
// for the period from 1st August 2004 to 1st Febrary 2005
// Rate = {average of the Euribor3M from 1st mars to 31 mars 2004 
//           + (CPI on 1st mars 2004 / CPI on 1st mars 2003 -1)} / 2 + 0.25%    
double ARM_LivretACurve::RateLivretA(ARM_InfCurv* infCurv, 
                                     ARM_ZeroCurve* euriborCurv, 
                                     ARM_Date& dateIn, 
                                     ARM_ResetManager* EuriborResetManager)
{    
    double rateEuribor = AverageZCRate(euriborCurv, dateIn, EuriborResetManager);
    
    double ratioCPI = LivretACPIRatio(dateIn, infCurv);
    
    double rateLivretA = (rateEuribor+(ratioCPI-1)*100)/2.00 +0.25;
    
    if ( itsFlagRounding == 1 )
    {
       rateLivretA = Rounding(rateLivretA);
    }

    return(rateLivretA);
}


double ARM_LivretACurve::RateLivretA(ARM_InfCurv* infCurv, 
                                     ARM_ZeroCurve* euriborCurv, 
                                     double dateJulIn, 
                                     ARM_ReferenceValue* fixingEuribor)
{    
    ARM_Date dateIn = ARM_Date(dateJulIn);

    double rateEuribor = AverageZCRate(euriborCurv, dateIn, fixingEuribor);
    
    double ratioCPI = LivretACPIRatio(dateIn, infCurv);
    
    double rateLivretA = (rateEuribor+(ratioCPI-1)*100)/2.00+0.25;
    
    if ( itsFlagRounding == 1 )
    {
       rateLivretA = Rounding(rateLivretA);
    }

    return(rateLivretA);
}



/*!
 *    calculate the livretA CPI ratio
 */
double ARM_LivretACurve::LivretACPIRatio(ARM_Date& numDateIn,ARM_InfCurv* infCurv)
{
    long dailyInterpType = infCurv->GetDailyInterpType();

    double     CPInum, CPIdenom, CPIRatio;

    ARM_Date numDate,denomDate, tmpDate;

    tmpDate   = numDateIn;
    
    numDate   = tmpDate;
    
    denomDate = tmpDate.AddYears(-1);
    
    CPIdenom = ValueCPI(denomDate,infCurv);
    
    CPInum = ValueCPI(numDate, infCurv);
    
    CPIRatio = CPInum/CPIdenom;
    
    return(CPIRatio); 
}



double ARM_LivretACurve::LivretACPIRatio(double numJulDate, ARM_InfCurv* infCurv)
{
    long dailyInterpType = infCurv->GetDailyInterpType();

    double CPIdenom, CPInum, CPIRatio;

    ARM_Date numDate, denomDate, tmpDate;
    
    tmpDate      = ARM_Date(numJulDate);
    
    numDate   = ARM_Date(numJulDate);
    
    denomDate = tmpDate.AddYears(-1);
    
    CPIdenom  = ValueCPI(denomDate, infCurv);
    
    CPInum = ValueCPI(numDate, infCurv);
    
    CPIRatio = CPInum/CPIdenom;
    
    return(CPIRatio); 
}



/*!
 *    get the value of the CPI according to the inflation curve
 */

double ARM_LivretACurve::ValueCPI(const ARM_Date& resetDate, ARM_InfCurv* infCurv)
{
    long dailyInterpType = infCurv->GetDailyInterpType();
    
	ARM_Date TmpDate = resetDate;

	TmpDate.AddMonths(3);
    
    return(infCurv->CPIInterpolate(resetDate, TmpDate, dailyInterpType));
}



/*!
 *    Average Euribor Rate
 */
double ARM_LivretACurve::AverageZCRate(ARM_ZeroCurve* euriborCurv, 
                                       const ARM_Date& dateEuribor,
                                       ARM_ReferenceValue* fixingEuribor)
{
    ARM_Date asofdate = GetAsOfDate();        

    double sumValue = 0.0;
    double result   = 0.0;
    int i = 0;
    ARM_Date StartDateEurib = dateEuribor;

	ARM_Date EndDateEurib = StartDateEurib;

	char* calend = euriborCurv->GetCurrencyUnit()->GetCcyName();

    int compteur = 0;
    int nb = 0;

    ARM_Vector* vDates = NULL; 

    if ( fixingEuribor != NULL )
    {
       vDates = fixingEuribor->GetDiscreteDates();
    }

    if ( vDates == NULL ) // No Fixings exist
    {
        if (asofdate > StartDateEurib )
		    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
                            "check the Euribor fixing data");
		else
		{
			while ( EndDateEurib.GetMonth() == StartDateEurib.GetMonth() )
			{
				ARM_Date tmpDateConstatation(EndDateEurib);
				tmpDateConstatation.AdjustToBusDate(calend,K_PREVIOUS);

				ARM_Date tmpDateFwdStart(tmpDateConstatation);
				tmpDateFwdStart.NextBusinessDay(euriborCurv->GetCurrencyUnit()->GetSpotDays(),calend);

				nb++;

				sumValue += euriborCurv->ForwardYield(tmpDateFwdStart,
                              ARM_Date(tmpDateFwdStart.GetJulian()).AddMonths(3), -1, 1);

				EndDateEurib.AddDays(1);
			}
		}

		result = sumValue/nb ;

		return(result);
    }
    
    /// Fixings exist
	/// EndDateEurib = dateEuribor;
    EndDateEurib = StartDateEurib;

	nb=0;

    while ( EndDateEurib.GetMonth() == StartDateEurib.GetMonth() )
    {
		ARM_Date tmpDateConstatation(EndDateEurib);
		tmpDateConstatation.AdjustToBusDate(calend,K_PREVIOUS);

		ARM_Date tmpDateFwdStart(tmpDateConstatation);
		tmpDateFwdStart.NextBusinessDay(euriborCurv->GetCurrencyUnit()->GetSpotDays(),calend);

        nb++;

        if ( asofdate.GetJulian() > tmpDateConstatation.GetJulian() )
        {
            while (( compteur < vDates->GetSize() )
                    && (tmpDateConstatation > vDates->Elt(compteur)) )
                compteur++;

            if ( compteur >= vDates->GetSize() )
               throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
                        "check the Euribor fixing data");

            if ( tmpDateConstatation == vDates->Elt(compteur) )
               sumValue += fixingEuribor->GetDiscreteValues()->Elt(compteur);
            else
               throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
                        "check the Euribor fixing data");
        }
        else if ( asofdate.GetJulian() == tmpDateConstatation.GetJulian() )
        {
            if (( compteur == vDates->GetSize()-1 ) 
                || 
                ( tmpDateConstatation != vDates->Elt(compteur+1)))
					sumValue += euriborCurv->ForwardYield(tmpDateFwdStart,
                             ARM_Date(tmpDateFwdStart.GetJulian()).AddMonths(3),-1,1);

            else
                sumValue += fixingEuribor->GetDiscreteValues()->Elt(compteur+1);
        }
        else
		{
			sumValue += euriborCurv->ForwardYield(tmpDateFwdStart,
                                 ARM_Date(tmpDateFwdStart.GetJulian()).AddMonths(3),-1,1);
		}

        EndDateEurib.AddDays(1);
    }

    result = sumValue/nb;

	return(result);
}


/*!
 *    Average Euribor Rate
 */
double ARM_LivretACurve::AverageZCRate(ARM_ZeroCurve* euriborCurv, 
                                       const ARM_Date& dateEuribor,
                                       ARM_ResetManager* EuriborResetManager)
{
    ARM_Date asofdate = GetAsOfDate();

    double sumValue = 0.0;
    double result   = 0.0;
    int i = 0;
    ARM_Date StartDateEurib = dateEuribor;

	ARM_Date EndDateEurib = StartDateEurib;

    char* calend = euriborCurv->GetCurrencyUnit()->GetCcyName();

    int compteur = 0;
    int nb = 0;

    if ( EuriborResetManager == NULL ) // No Fixings exist
    {
        if (asofdate > StartDateEurib )
		    throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
                            "check the Euribor fixing data");
		else
		{
			while ( EndDateEurib.GetMonth() == StartDateEurib.GetMonth() )
			{
				ARM_Date tmpDateConstatation(EndDateEurib);
				tmpDateConstatation.AdjustToBusDate(calend,K_PREVIOUS);

				ARM_Date tmpDateFwdStart(tmpDateConstatation);
				tmpDateFwdStart.NextBusinessDay(euriborCurv->GetCurrencyUnit()->GetSpotDays(),calend);

				nb++;

				sumValue += euriborCurv->ForwardYield(tmpDateFwdStart,
                              ARM_Date(tmpDateFwdStart.GetJulian()).AddMonths(3), -1, 1);

				EndDateEurib.AddDays(1);
			}
		}

		result = sumValue/nb ;

		return(result);
    }
    
    /// Fixings exist
	/// EndDateEurib = dateEuribor;
    EndDateEurib = StartDateEurib;

	nb=0;

    while ( EndDateEurib.GetMonth() == StartDateEurib.GetMonth() )
    {
		ARM_Date tmpDateConstatation(EndDateEurib);
		tmpDateConstatation.AdjustToBusDate(calend,K_PREVIOUS);

		ARM_Date tmpDateFwdStart(tmpDateConstatation);
		tmpDateFwdStart.NextBusinessDay(euriborCurv->GetCurrencyUnit()->GetSpotDays(),calend);

        nb++;

        if ( asofdate.GetJulian() > tmpDateConstatation.GetJulian() )
        {
			sumValue += EuriborResetManager->GetReset(tmpDateConstatation.GetJulian());
        }
        else if ( asofdate.GetJulian() == tmpDateConstatation.GetJulian() )
        {
			try
			{
				sumValue += euriborCurv->ForwardYield(tmpDateFwdStart,
								 ARM_Date(tmpDateFwdStart.GetJulian()).AddMonths(3),-1,1);
			}
			catch(...)
			{
                sumValue += EuriborResetManager->GetReset(tmpDateConstatation.GetJulian());
			}
        }
        else
		{
			sumValue += euriborCurv->ForwardYield(tmpDateFwdStart,
                                 ARM_Date(tmpDateFwdStart.GetJulian()).AddMonths(3),-1,1);
		}
        EndDateEurib.AddDays(1);
    }

    result = sumValue/nb;

	return(result);
}


/*!
 *   Rounding Calculation 
 */
// rounded at the point closer to 1/4  eg. 2.208% ~ 2.25%
double ARM_LivretACurve::Rounding(double valueIn)
{
    double tmp = (valueIn-long(valueIn));
    
    long flag = long(tmp/0.125);
    
    switch(flag)
    {
        case 0: tmp = 0; break;
        case 1:
        case 2: tmp = 0.25; break;
        case 3: 
        case 4: tmp = 0.5; break;
        case 5: 
        case 6: tmp = 0.75; break;
        case 7:
        case 8: tmp = 1; 
    }

    double value = long(valueIn)+tmp;

    return value;
}



/*!
 *    Get the debut of echeancier
 */
double ARM_LivretACurve::GetBeginCalculSchedule(void)
{
    int month = GetAsOfDate().GetMonth();
    int day = 1;
    int year = GetAsOfDate().GetYear();


    if (( month >= 2 ) && ( month < 8 ))
    {
       month = 8;
    }
    else if (( month >= 8 ) && ( month <= 12 ))
    {
       month = 2;

       year = year+1;
    }
    else
       month = 2;
    
    ARM_Date beginEch(day, month, year);
    
    return(beginEch.GetJulian());
}



double ARM_LivretACurve::LivretAInterpolate(double yearTerm, int& nb)
{
	if (( yearTerm < GetYearTerms()->Elt(0) ) || 
			( yearTerm > GetYearTerms()->Elt(GetYearTerms()->GetSize()-1) ))
	{
		   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
						   "Date is out of the range" );
	}

    int i = 0;

    ARM_Vector* dateyearTerms = GetYearTerms();

    ARM_Vector* rateValues = GetZeroRates();

    while (!(( yearTerm >= dateyearTerms->Elt(i))
              && ( yearTerm < dateyearTerms->Elt(i+1) ) ))
    {
        i++;
    }
    
    nb = i;
    
    double value = rateValues->Elt(i);

    return(value);
}



double ARM_LivretACurve::LivretAInterpolate(const ARM_Date& rateDate, int& nb)
{
    double yearTerm = (rateDate.GetJulian()-GetAsOfDate().GetJulian())/K_YEAR_LEN;

    double value = LivretAInterpolate(yearTerm, nb);
      
    return(value);
}



// compMeth==0: no composing
double ARM_LivretACurve::DiscountFunctionLivretA(double yearTerm, int compMeth)
{
	if (( yearTerm < GetYearTerms()->Elt(0)) 
			||
			( yearTerm > GetYearTerms()->Elt(GetYearTerms()->GetSize()-1)))
	{
		   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						   "Date is out of the range" );
	}
    
    double T = yearTerm;

    int nbT0 = 0;
	int nbT = 0;

    double rateT0 = 0.0;
	double rateT = 0.0;

    rateT0 = LivretAInterpolate(0, nbT0);
    
    rateT  = LivretAInterpolate(T, nbT); 
    
    double z = 0.0;
 
    if (( compMeth == 0 ) || ( nbT0 == nbT )) 
    {
       z = exp(-rateT*T);
    }
    else
    {
       z = exp(-rateT0*(GetYearTerms()->Elt(nbT0+1)-1.00/K_YEAR_LEN))
            *exp(-rateT*(T-GetYearTerms()->Elt(nbT)+1.00/K_YEAR_LEN));

       int k = nbT;

       while ( k > nbT0+1 )
       {
           z *= exp(-GetZeroRates()->Elt(k-1)
                   *(GetYearTerms()->Elt(k)-GetYearTerms()->Elt(k-1)));

           k--;
       }
    }  
    
    return(z);
}



double ARM_LivretACurve::DiscountFunctionLivretA(const ARM_Date& dateIn, 
                                                 int compMeth)
{
    double T = (dateIn.GetJulian() - GetAsOfDate().GetJulian())/K_YEAR_LEN;

    double z = DiscountFunctionLivretA(T, compMeth);
    
    return(z);
}


double ARM_LivretACurve::DiscountFunction(double yearTerm)
{
	if (( yearTerm < GetYearTerms()->Elt(0))
			|| 
			( yearTerm > GetYearTerms()->Elt(GetYearTerms()->GetSize()-1)))
	{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
							"Date is out of the range" );
	}

    double T = yearTerm;

    int nbT;

    double rateT;    

    rateT  = LivretAInterpolate(T, nbT)/100.00; 
    
    double z = exp(-rateT*T);

    return(z);
}



double ARM_LivretACurve::DiscountYield(double yearTerm, int compMeth)
{
    double z = DiscountFunctionLivretA(yearTerm, compMeth);

    double value = -log(z)/yearTerm;
    
    return(value);
}



double ARM_LivretACurve::DiscountYield(ARM_Date& maturity, int compMeth)
{
    double z = DiscountFunctionLivretA(maturity, compMeth);

    double value = -log(z)/(maturity.GetJulian()-GetAsOfDate().GetJulian());

    return(value);
}



double ARM_LivretACurve::ForwardYield(double yearTerm1, double yearTerm2, 
                                      int compMeth)
{
    if ( yearTerm1 > yearTerm2 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                        "YearTerm2 must be after YearTerm2");
    }

    double z1 = DiscountFunctionLivretA(yearTerm1, compMeth);

    double z2 = DiscountFunctionLivretA(yearTerm2, compMeth);
    
    double value = -log(z2/z1)/(yearTerm2-yearTerm1 );

    return(value);
}



double ARM_LivretACurve::ForwardYield(ARM_Date& maturity1, 
                                      ARM_Date& maturity2, int compMeth)
{
    double yearTerm1 = (maturity1.GetJulian() - GetAsOfDateJul())/K_YEAR_LEN;
    
    double yearTerm2 = (maturity2.GetJulian() - GetAsOfDateJul())/K_YEAR_LEN;

    double value = ForwardYield(yearTerm1, yearTerm2, compMeth);
    
    return(value);
}



int ARM_LivretACurve::GetRateDateNb(ARM_Date& dateIn)
{
    long date = dateIn.GetJulian() - GetAsOfDateJul();

    if (( date < GetDateTerms()->Elt(0) )
        || 
        ( date > GetDateTerms()->Elt(GetYearTerms()->GetSize()-1 )))
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
                        "Date is out of the range" );
    }

    int i = 0;

    ARM_Vector* dateTerms = GetDateTerms();

    while (!(( date >= dateTerms->Elt(i))
              && ( date < dateTerms->Elt(i+1) ) ))
    {
        i++;
    }
    
    return(i);
}



void ARM_LivretACurve::SetFlagRounding(long flag)
{
    itsFlagRounding = flag;
}



long ARM_LivretACurve::GetFlagRounding() const
{
    return itsFlagRounding;
}


/*!
 * Destructor
 */
ARM_LivretACurve::~ARM_LivretACurve(void)
{
}



// Copy constructor
ARM_LivretACurve::ARM_LivretACurve(const ARM_LivretACurve& srcLivreACurv)
                                  :ARM_ZeroCurve(srcLivreACurv)
{
    Init();
    
    BitwiseCopy(&srcLivreACurv);
}


//    operator = 
ARM_LivretACurve& ARM_LivretACurve::operator = 
                                    (const ARM_LivretACurve& srcLivreACurv)
{
    if ( this == &srcLivreACurv ) 
    {
       return(*this);
    }

    (*this).ARM_ZeroCurve::operator = (srcLivreACurv);
    
    BitwiseCopy(&srcLivreACurv);

    return(*this);
}



//    function for ARM_Objet compatibility
void ARM_LivretACurve::BitwiseCopy(const ARM_Object* srcobj)
{
    ARM_LivretACurve* src = (ARM_LivretACurve *) srcobj;

    itsFlagRounding = src->itsFlagRounding;
}



void ARM_LivretACurve::Copy(const ARM_Object* src)
{
    ARM_ZeroCurve::Copy(src);

    BitwiseCopy(src);
}


ARM_Object* ARM_LivretACurve::Clone(void)
{
    ARM_LivretACurve* theClone = new ARM_LivretACurve();

    theClone->Copy(this);
    
    return(theClone);
}


void ARM_LivretACurve::View(char* id, FILE* ficOut)
{
    /// to be consistent with other interface
    /// need to use fprintf
    FILE* fOut;
 
    char fOutName[40];
    
    char strDate[20];
    
    /// do we have already a file opened?
    if ( ficOut == NULL )
    {
        ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
        
        fOut = fopen(fOutName, "w");
    }
    else
    {
        fOut = ficOut;
    }
    
    /// printing of the asofDate 
    fprintf(fOut, "\n\n\t =====> Taux Livret A Curve \n\n");

    ARM_Date tempDate = GetAsOfDate();
    
    tempDate.JulianToStrDateDay(strDate);

    /// AsOfDate
    fprintf (fOut, "\n\n\t AsOfDate  : %s \n\n", strDate);
    
    /// currency and calendar
    fprintf (fOut, "\t Currency  : %s \t\t Calendar  : %s\n", 
             GetCurrencyUnit()->GetCcyName(), GetCurrencyUnit()->GetCcyName() );
    
    /// points of the curve
    ARM_Vector* dateTerms = GetDateTerms();
    ARM_Vector* yearTerms = GetYearTerms();

    ARM_Vector* rate = GetZeroRates();

    int sz = dateTerms->size();
    
    fprintf (fOut, "\n\n ResetDate\t\t YearFrac \t\t Taux Livret A \n\n");
    
    int  i;

    for (i = 0; i < sz; ++i)
    {
        tempDate = ARM_Date( dateTerms->Elt(i)+GetAsOfDate().GetJulian() );

        tempDate.JulianToStrDateDay(strDate);
        
        fprintf (fOut, " %s \t %.4lf \t\t %.3lf \t\t \n", 
                strDate, yearTerms->Elt(i), rate->Elt(i) );
    }            

    if ( ficOut == NULL )
        fclose(fOut);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------*/
/*---- End Of File ----*/

