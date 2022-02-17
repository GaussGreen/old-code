/*
 * $Log: zeroint.cpp,v $
 * Revision 1.30  2004/01/21 14:31:26  jpriaudel
 * modif dans le bump
 *
 * Revision 1.29  2004/01/20 16:25:34  mab
 * Replace AddMonth by :
 * AddMonths(NbMonth, GOTO_END_OF_MONTH);
 *
 * Revision 1.28  2004/01/20 16:14:23  mab
 * dateCour.AddMonths(NbMonth); replaced by:
 * dateCour.AddMonths(NbMonth, GOTO_END_OF_MONTH);
 * in Spreaded curve
 *
 * Revision 1.27  2004/01/20 15:37:50  mab
 * Improvements in : Basis curves generation
 *
 * Revision 1.26  2003/09/25 09:16:48  mab
 * "Constification"
 *
 * Revision 1.25  2003/09/03 18:10:18  jpriaudel
 * ajout d'un constructeur plus le shift
 *
 * Revision 1.24  2003/09/01 14:00:05  mab
 * In double ARM_ZeroLInterpol::DiscountFunction(double yearTerm)
 * Return 0.0 when  yearTerm < 0.0
 * instead of exception
 *
 * Revision 1.23  2003/08/20 18:31:02  jpriaudel
 * memory leaked dans shiftcurve
 *
 * Revision 1.22  2003/06/16 13:21:22  jpriaudel
 * warning en moins
 *
 * Revision 1.21  2003/05/09 14:38:33  jpriaudel
 * Modif GenerateShiftCurve
 * inversement du shift dans le cas des futures
 *
 * Revision 1.20  2003/05/09 09:55:56  jpriaudel
 * Modification du GenerateShiftCurve
 * Inversion des 2 boucles For
 *
 * Revision 1.19  2003/04/09 17:11:49  jpriaudel
 * correction dans le GenerateShiftCurve
 * (Ajout des break dans les case)
 *
 * Revision 1.18  2003/03/19 16:25:15  sgasquet
 * Ajout de la ccy dans le premier constructeur
 *
 * Revision 1.17  2002/10/11 08:21:37  mab
 * Improvements
 *
 * Revision 1.16  2002/09/24 13:18:42  mab
 * Added : TOY management (See TOY section)
 *
 * Revision 1.15  2002/05/30 13:38:10  mab
 * char terms[ARM_NB_TERMS][6] : 6 replaced by 12
 *
 * Revision 1.14  2001/11/06 11:57:48  mab
 * Mise en forme
 *
 * Revision 1.13  2001/08/17 15:15:41  smysona
 * Methodes nv non vrituelles
 *
 * Revision 1.12  2001/07/02 08:40:17  abizid
 * Gestion des taux extrapoles amelioree
 *
 * Revision 1.11  2001/05/25 11:28:09  mab
 * Rajout de : delete sprdeeDates ; delete MMValues;  delete MMDates;
 *
 * Revision 1.10  2001/04/25 09:22:36  abizid
 * Debug creation (free memory)
 *
 * Revision 1.9  2001/04/03 12:00:48  nicolasm
 * unused variables et appel de linInterpol avec double*
 *
 * Revision 1.8  2001/02/02 10:31:14  mab
 * Suite correction precedente
 *
 * Revision 1.7  2001/02/01 19:27:51  mab
 * Modif MA,SM,AB, Si SpotDay=0 (GBP), utilisation de DiscountYield
 * a la place de ForwardYield pour courbe spreadee
 *
 * Revision 1.6  2001/02/01 16:31:29  mab
 * Remplacer : .AddMonths(NbMonth, GOTO_END_OF_MONTH);
 * par : .AddMonths(NbMonth)
 * Comme c'etait avant
 *
 * Revision 1.5  2000/12/12 12:15:15  smysona
 * Rajout d'un constructeur a partir d'une zerocurve qqconque
 *
 * Revision 1.4  2000/12/11 16:10:19  mab
 * Ds la generation de courbes : Pour la partie monetaire
 * >> Rajout d'un parametre a AddMonth
 *
 * Revision 1.3  1999/03/02 15:14:44  ypilchen
 * Rajout du comment. RCS
 *
 */


/*----------------------------------------------------------------------------*
 
    zeroint.cpp
 
    This file implements the ARM_ZeroLInterpol class, a class for 
         computing a ARM_ZeroCurve linear interpolation.

*----------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "util.h"
#include "zeroint.h"
#include "currency.h"
#include "expt.h"
#include "newton.h"
#include "paramview.h"




/*---------------------------------------------------------------------------*/


ARM_ZeroLInterpol::ARM_ZeroLInterpol(ARM_Date& asOf, 
                     ARM_Vector* yearTerms,
                     ARM_Vector* zeroYields, 
                     int compoundMeth, int lastBucketInt, 
                     int interpolMeth, 
                     ARM_Currency* ccy):ARM_ZeroCurve(asOf)
{
    // Check variables

    SetName(ARM_ZERO_LIN_INTERPOL);
    SetCurrencyUnit(ccy);
/*
    if ( yearTerms->GetSize() < K_MIN_NUM_YEAR_TERMS )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
        "The input year terms and yields must have a minimum of 3 elements");
    }
*/

    if ( yearTerms->GetSize() != zeroYields->GetSize() )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
        "The input year terms and yields must have the same size");
    }

    // Set variables

    if (lastBucketInt)
    {
       SetYearTerms((ARM_Vector*) yearTerms->Clone());

       SetZeroRates((ARM_Vector*) zeroYields->Clone());
    }
    else
    {
       int size = yearTerms->GetSize();

       ARM_Vector* terms = new ARM_Vector(yearTerms, size+1, 0, size-1);

       terms->Elt(size) = 1000.0;

       ARM_Vector* yields = new ARM_Vector(zeroYields, size+1, 0, size-1);

       yields->Elt(size) = zeroYields->Elt(size-1);

       SetYearTerms(terms);

       SetZeroRates(yields);
    }

    int Size = GetYearTerms()->GetSize();

    ARM_Vector* BPShifts = new ARM_Vector(Size, 0.0);

    for (int i = 0; i < Size; i++)
    {
        if ( GetYearTerms()->Elt(i) >= GetBucketStartPeriod()
                && GetYearTerms()->Elt(i) < GetBucketEndPeriod())
        {
           BPShifts->Elt(i) = GetBPShift()*0.01; // using yield in percent
        }
    }

    SetBPShifts(BPShifts);

    itsCompoundMeth = compoundMeth;

    itsLastBucketInt = lastBucketInt;

    itsInterpolMeth = interpolMeth;

    GenerateDateTerms();

    GenerateDiscountFactors(compoundMeth);
}



ARM_ZeroLInterpol::ARM_ZeroLInterpol(ARM_Date& asOf,
                                     ARM_Vector* zeroYields,
                                     ARM_CRV_TERMS& sTerms,
                                     int compMeth,
                                     int lastBucketInt, 
                                     int interpolMeth, 
                                     ARM_Currency* ccy):ARM_ZeroCurve(asOf)
{
    SetName(ARM_ZERO_LIN_INTERPOL);
    SetCurrencyUnit(ccy);

    // Set variables

    ARM_MarketData* Data = new ARM_MarketData(sTerms, zeroYields, compMeth, 
                                      lastBucketInt, K_PAR, interpolMeth, 3);

    SetMktData(Data);

    long nbSpotDays = ccy->GetSpotDays();

    ARM_Vector* yearterms = new ARM_Vector(zeroYields->GetSize());
    ARM_Date settleDate(asOf);

    settleDate.NextBusinessDay(nbSpotDays, ccy->GetCcyName());

    for (int i = 0; i < zeroYields->GetSize(); i++)
    {
        long isDouble = 0;

        int Nb;
        char cMatu;
        long freqId;

        sscanf(sTerms[i], "%d%c", &Nb, &cMatu);

        cMatu = toupper(cMatu);

        if ( cMatu == 'D' ) // Ex : "1D"
           freqId = K_DAILY;
        else if ( cMatu == 'W' )  
           freqId = K_WEEKLY;
        else if ( cMatu == 'M' ) 
           freqId = K_MONTHLY;
        else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
           freqId = K_ANNUAL;
        else
        {
           delete yearterms;
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                           "Bad Term for BS curve");
        }

        ARM_Date tmpDate;

        if ( freqId == K_DAILY )
        {
           tmpDate = asOf;
           tmpDate.NextBusinessDay(Nb, ccy->GetCcyName());
        }
        else
        {
           tmpDate = settleDate;
           tmpDate.AddPeriodMult(freqId, Nb, ccy->GetCcyName()).AdjustToBusDate(ccy->GetCcyName(), K_FOLLOWING);
        }

        yearterms->Elt(i) =  (tmpDate-asOf)/365.0;
    }

    if (lastBucketInt)
    {
       SetYearTerms((ARM_Vector*) yearterms->Clone());

       SetZeroRates((ARM_Vector*) zeroYields->Clone());
    }
    else
    {
       int size = yearterms->GetSize();

       ARM_Vector* terms = new ARM_Vector(yearterms, size+1, 0, size-1);

       terms->Elt(size) = 1000.0;

       ARM_Vector* yields = new ARM_Vector(zeroYields, size+1, 0, size-1);

       yields->Elt(size) = zeroYields->Elt(size-1);

       SetYearTerms(terms);

       SetZeroRates(yields);
    }

    delete yearterms;

    int Size = GetYearTerms()->GetSize();

    ARM_Vector* BPShifts = new ARM_Vector(Size, 0.0);

    for (int i = 0; i < Size; i++)
    {
        if ( GetYearTerms()->Elt(i) >= GetBucketStartPeriod()
                && GetYearTerms()->Elt(i) < GetBucketEndPeriod())
        {
           BPShifts->Elt(i) = GetBPShift()*0.01; // using yield in percent
        }
    }

    SetBPShifts(BPShifts);

    itsCompoundMeth = compMeth;

    itsLastBucketInt = lastBucketInt;

    itsInterpolMeth = interpolMeth;

    GenerateDateTerms();

    GenerateDiscountFactors(compMeth);
}



ARM_ZeroLInterpol::ARM_ZeroLInterpol(ARM_Date& asOf, 
                                     ARM_CRV_TERMS& terms, 
                                     ARM_Vector* mktData, 
                                     int MMVsFut, 
                                     int SwapVsFut, 
                                     int raw, 
                                     int Cont_Lin, 
                                     int lastBucketInt, 
                                     ARM_Currency* ccy,
									 int swapFrqId,
									 int fixDayCount) : ARM_ZeroCurve(asOf, 
                                                           terms, 
                                                           mktData, 
                                                           MMVsFut, 
                                                           SwapVsFut, 
                                                           raw, 
                                                           Cont_Lin, 
                                                           ccy,
														   swapFrqId,
														   fixDayCount)
{
    // Check variables

    SetName(ARM_ZERO_LIN_INTERPOL);

    char* smoothMeth = ARM_ParamView::GetMappingName(S_MISSING_MATURITY_FILLING, raw);

	string	strSmoothMeth(smoothMeth);

	SetSmoothingMethod(strSmoothMeth);

    // Set variables

    if (!lastBucketInt)
    {
       int size = GetYearTerms()->GetSize();

       ARM_Vector* terms = new ARM_Vector(GetYearTerms(), size+1, 0, size-1);

       terms->Elt(size) = 1000.0;

       ARM_Vector* yields = new ARM_Vector(GetZeroRates(), size+1, 0, size-1);

       yields->Elt(size) = GetZeroRates()->Elt(size-1);        
        
       ARM_Vector* df = new ARM_Vector(GetDiscountFactors(), size+1,
                                        0, size-1);

       df->Elt(size) = exp(-10.0*GetZeroRates()->Elt(size-1));        
        
       SetYearTerms(terms);

       SetZeroRates(yields);

       SetDiscountFactors(df);
    }

    int Size = GetYearTerms()->GetSize();
 
    ARM_Vector* BPShifts = new ARM_Vector(Size, 0.0);
 
    for (int i = 0; i < Size; i++)
    {
        if ( GetYearTerms()->Elt(i) >= GetBucketStartPeriod()
                && GetYearTerms()->Elt(i) < GetBucketEndPeriod())
        {
            BPShifts->Elt(i) = GetBPShift()*0.01; // using yield in percent
        }
    }

    SetBPShifts(BPShifts);

    itsCompoundMeth = K_COMP_CONT;

    itsInterpolMeth = Cont_Lin;

    itsLastBucketInt = lastBucketInt;
}



ARM_ZeroLInterpol::ARM_ZeroLInterpol(ARM_CRV_TERMS& terms,
                                     ARM_Date& asOf,
                                     ARM_Vector* mktData,
                                     int MMVsFut,
                                     int SwapVsFut,
                                     int raw,
                                     int Cont_Lin,
                                     int lastBucketInt,
                                     ARM_Currency* ccy) : ARM_ZeroCurve(terms,
                                                           asOf,
                                                           mktData,
                                                           MMVsFut,
                                                           SwapVsFut,
                                                           raw,
                                                           Cont_Lin,
                                                           ccy)
{
    // Check variables
 
    SetName(ARM_ZERO_LIN_INTERPOL);
 
    // Set variables
 
    if (!lastBucketInt)
    {
       int size = GetYearTerms()->GetSize();
 
       ARM_Vector* terms = new ARM_Vector(GetYearTerms(), size+1, 0, size-1);
 
       terms->Elt(size) = 1000.0;
 
       ARM_Vector* yields = new ARM_Vector(GetZeroRates(), size+1, 0, size-1);
 
       yields->Elt(size) = GetZeroRates()->Elt(size-1);
 
       ARM_Vector* df = new ARM_Vector(GetDiscountFactors(), size+1,
                                       0, size-1);
 
       df->Elt(size) = exp(-10.0*GetZeroRates()->Elt(size-1));
 
       SetYearTerms(terms);
 
       SetZeroRates(yields);
 
       SetDiscountFactors(df);
    }
 
    int Size = GetYearTerms()->GetSize();
 
    ARM_Vector* BPShifts = new ARM_Vector(Size, 0.0);
 
    for (int i = 0; i < Size; i++)
    {
        if ( GetYearTerms()->Elt(i) >= GetBucketStartPeriod()
                && GetYearTerms()->Elt(i) < GetBucketEndPeriod())
        {
            BPShifts->Elt(i) = GetBPShift()*0.01; // using yield in percent
        }
    }

    SetBPShifts(BPShifts);
 
    itsCompoundMeth = K_COMP_CONT;
 
    itsInterpolMeth = Cont_Lin;
 
    itsLastBucketInt = lastBucketInt;
}



ARM_ZeroLInterpol::ARM_ZeroLInterpol(ARM_Date& asOf, 
                                     ARM_CRV_TERMS& terms, 
                                     ARM_Vector* mktData, 
                                     ARM_Container* bonds, 
                                     ARM_Vector* yields,
                                     int MMVsFut, 
                                     int lastBucketInt, 
                                     ARM_Currency* ccy):ARM_ZeroCurve(asOf, 
                                                                      terms, 
                                                                      mktData, 
                                                                      bonds, 
                                                                      yields, 
                                                                      MMVsFut, 
                                                                      ccy)
{
    // Check variables

    SetName(ARM_ZERO_LIN_INTERPOL);

    // Set variables

    if (!lastBucketInt)
    {
       int size = GetYearTerms()->GetSize();

       ARM_Vector* terms = new ARM_Vector(GetYearTerms(), size+1, 0, size-1);

       terms->Elt(size) = 1000.0;

       ARM_Vector* yields = new ARM_Vector(GetZeroRates(), size+1, 0, size-1);

       yields->Elt(size) = GetZeroRates()->Elt(size-1);        
        
       ARM_Vector* df = new ARM_Vector(GetDiscountFactors(), size+1,
                                       0, size-1);

       df->Elt(size) = exp(-10.0*GetZeroRates()->Elt(size-1));        
                
       SetYearTerms(terms);

       SetZeroRates(yields);

       SetDiscountFactors(df);
    }

    int Size = GetYearTerms()->GetSize();
 
    ARM_Vector* BPShifts = new ARM_Vector(Size, 0.0);
 
    for (int i = 0; i < Size; i++)
    {
        if ( GetYearTerms()->Elt(i) >= GetBucketStartPeriod()
                && GetYearTerms()->Elt(i) < GetBucketEndPeriod())
        {
            BPShifts->Elt(i) = GetBPShift()*0.01; // using yield in percent
        }
    }

    SetBPShifts(BPShifts);

    itsCompoundMeth  = K_COMP_CONT;

    itsInterpolMeth = K_LINEAR;

    itsLastBucketInt = lastBucketInt;
}




// **********   Spreaded Curves ****************

ARM_ZeroLInterpol::ARM_ZeroLInterpol(ARM_Date asofdate,
                                     ARM_ZeroCurve* ZCSpread, 
                                     ARM_ZeroCurve* zcInit,
                                     int MMFreq, int SwapFreq,
                                     ARM_Currency* ccy)
{
    Init();
 
    SetName(ARM_ZERO_LIN_INTERPOL);
 
    SetAsOfDate(asofdate);
 
    SetCurrencyUnit(ccy);

    SetCompoundMeth(K_ZEROCOUPON);
 
    itsInterpolMeth = ((ARM_ZeroLInterpol *) zcInit)->itsInterpolMeth;

    int spotDay = ccy->GetSpotDays();
    double dayCount = FromDayCountToBasis(ccy->GetMMDayCount());
    
	double libDayCount = FromDayCountToBasis(ccy->GetLiborIndexDayCount());
    
	char ccyName[20];
	
	strcpy(ccyName, ccy->GetCcyName());
 
	int ccyMMdayCount = ccy->GetMMDayCount();

    if ( strcmp((const char *) ccyName, "JPY") == 0 )
	{
#ifdef WIN32

		if ( GetCALYPSOVersion == 0 )
		   strcpy(ccyName, "ZGJ");
		else
		   strcpy(ccyName, "GBPJPY");
#else
        strcpy(ccyName, "ZGJ");
#endif
	}

    int i;

    ARM_Date spotDate  = asofdate;
    spotDate.NextBusinessDay(spotDay, ccyName);
    
	ARM_Date dateCour = spotDate;
    
	ARM_Date BusDate, prevDate;
 
    double df_start = 0.;
    double df_cour = 0.;

 
    // ***********   Segment Monetaire   ****************
 
    double fwdRate, spread, matu1, matu2, fixFlow, cumFlow;
    double spotRate, term;

    int SizeM   = MMFreq;
	int SpotPoints = 0;
    
	int NbMonth = 12/MMFreq; // nb de mois entre chaque plot
 
    // PL : Like Summit : Spot, Spot+Spot, Spot+1W ?
    ARM_Vector MMDates_tmp(3);
    ARM_Vector MMValues_tmp(3);
 
    BusDate = spotDate;
    BusDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

    MMDates_tmp.Elt(0) = BusDate.GetJulian();

    matu1 = spotDate-asofdate;

	// SPOT
	SpotPoints++;
    
    if ( spotDay == 0 )
	{
		// spotDay = 0 => spotDate = next Business date !

		ARM_Date newSpotDate = spotDate;
		newSpotDate.NextBusinessDay(ccyName);
		
		fwdRate = zcInit->DiscountYield(newSpotDate, -1);

		term = CountYears(ccyMMdayCount, asofdate, spotDate);

		fwdRate *= dayCount/K_YEAR_LEN; // Only Fwd correction using basis!

		spread = ZCSpread->DiscountYield(BusDate, -1);

		double newMatu1 = newSpotDate-asofdate;

		fixFlow = (fwdRate/100.0-spread/10000.0)*(newMatu1/dayCount);
		
		MMValues_tmp.Elt(0) = (K_YEAR_LEN/newMatu1)*log(1.+fixFlow)*100.0; 
	}
	else
	{
		fwdRate = zcInit->ForwardYield(asofdate, spotDate, -1);

		term = CountYears(ccyMMdayCount, asofdate, spotDate);

		fwdRate *= dayCount/K_YEAR_LEN; // Only Fwd correction using basis!

		spread = ZCSpread->DiscountYield(BusDate, -1);

		fixFlow = (fwdRate/100.0-spread/10000.0)*term;

		MMValues_tmp.Elt(0) = (K_YEAR_LEN/matu1)*log(1.0+fixFlow)*100.0;
	}

	MMDates_tmp.Elt(0) = BusDate.GetJulian();

	// SPOT + SPOT
	double spotSpotDays = 2;

	SizeM++;
	SpotPoints++; // PL : add Spot+Spot
    ARM_Date spotSpotDate = spotDate;


	spotSpotDate.NextBusinessDay(spotSpotDays, ccyName);

	fwdRate = zcInit->ForwardYield(asofdate, spotSpotDate, -1);

	term = CountYears(ccyMMdayCount, asofdate, spotSpotDate);

	fwdRate *= dayCount/K_YEAR_LEN; // Only Fwd correction using basis!

	spread = ZCSpread->DiscountYield(spotSpotDate, -1);

	fixFlow = (fwdRate/100.0-spread/10000.0)*term;

	MMDates_tmp.Elt(1) = spotSpotDate.GetJulian();
	MMValues_tmp.Elt(1) = (K_YEAR_LEN/(spotSpotDate-asofdate))*log(1.0+fixFlow);
	MMValues_tmp.Elt(1) *= 100.0;

	if ( strcmp((const char *) ccy->GetCcyName(), "AUD") != 0 )
	{
		// SPOT + 1W (not in AUD)
		SizeM++;
		SpotPoints++; // PL : add Spot+1W
		ARM_Date weekSpotDate = spotDate;
		weekSpotDate.AddPeriod("1W", ccyName);

		fwdRate = zcInit->ForwardYield(asofdate, weekSpotDate, -1);

		term = CountYears(ccyMMdayCount, asofdate, weekSpotDate);

		fwdRate *= dayCount/K_YEAR_LEN; // Only Fwd correction using basis!

		spread = ZCSpread->DiscountYield(weekSpotDate, -1);

		fixFlow = (fwdRate/100.0-spread/10000.0)*term;

		MMDates_tmp.Elt(2) = weekSpotDate.GetJulian();
		MMValues_tmp.Elt(2) = (K_YEAR_LEN/(weekSpotDate-asofdate))*log(1.0+fixFlow);
		MMValues_tmp.Elt(2) *= 100.0;
	}

    ARM_Vector* MMDates  = new ARM_Vector(SizeM);
    ARM_Vector* MMValues = new ARM_Vector(SizeM);

	for (int k = 0; k < SpotPoints; k++)
	{
		MMDates->Elt(k) = MMDates_tmp.Elt(k);
		MMValues->Elt(k) = MMValues_tmp.Elt(k);
	}

    spotRate = MMValues->Elt(0);
    

	df_start = exp(-(spotRate/100.0)*(matu1/K_YEAR_LEN));

    cumFlow = 0.;

    for (i = SpotPoints; i < SizeM; i++)
    {
        // On calcule les plots de la nouvelle courbe
        // MMFreq est la frequence des plots (en 1/annee)
 
        prevDate = BusDate;

        dateCour.AddMonths(NbMonth, GOTO_END_OF_MONTH);

        BusDate = dateCour;
        BusDate.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);
 
        MMDates->Elt(i) = BusDate.GetJulian();
 
 
        // On calcule le taux fwd pour le plot a partir de zcInit
        if ( spotDay == 0 )
        {
           fwdRate = zcInit->DiscountYield(BusDate, -1);
        }
        else
        {
           fwdRate = zcInit->ForwardYield(spotDate, BusDate, -1);
        }

        fwdRate *= dayCount/K_YEAR_LEN;

        // On calcule le spread pour le plot a partir de spreads
 
        spread = ZCSpread->DiscountYield(BusDate, -1);
 
 
        // On calcule le flux fixe
 
		double term2 = CountYears(ccyMMdayCount, spotDate, BusDate);
        
		fixFlow = (fwdRate/100.0-spread/10000.0)*term2;

        matu1 = spotDate-asofdate;
        matu2 = BusDate-asofdate;
        

		df_cour = (df_start)/(1.+fixFlow);

        MMValues->Elt(i) = -(K_YEAR_LEN/matu2)*log(df_cour);

        MMValues->Elt(i) *= 100.;
    }
 
    // **************** Segment des Taux Longs *********************
 
    if ( ZCSpread->GetName() != ARM_ZERO_LIN_INTERPOL )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                       "ZCSpread must be Interpole Curve");
    }
 
	double sprdYearTermSize = zcInit->GetYearTerms()->GetSize();
    double lastDate = zcInit->GetYearTerms()->Elt(sprdYearTermSize-1);

    
    if ( lastDate > 900 )
       lastDate = zcInit->GetYearTerms()->Elt(sprdYearTermSize-2);
 

    // le vecteur SwDates ne commence qu'a 1Y d'ou le + 1Y sur asofdate

    ARM_Date oneYear = spotDate;
    oneYear.AddYears(1);

    ARM_Date resetDate, zeroTerm, unadjDate, swapMatu;

    int SizeS = floor((lastDate-1.)*SwapFreq)+1;
 
    ARM_Vector* sprdeeDates = new ARM_Vector(SizeS+SizeM);
    ARM_Vector* sprdeeValues = new ARM_Vector(SizeS+SizeM);
 
    // On recopie la partie monetaire
    for (i = 0; i < SizeM; i++)
    {
       sprdeeDates->Elt(i) = MMDates->Elt(i);
       sprdeeValues->Elt(i) = MMValues->Elt(i);
    }
 
    unadjDate = oneYear;
    dateCour = oneYear;
    dateCour.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);

    NbMonth = 12/SwapFreq; // nb de mois entre chaque plot

    for (i = SizeM ; i < SizeM+SizeS; i++)
    {
        sprdeeDates->Elt(i) = dateCour.GetJulian();
        sprdeeValues->Elt(i) = MMValues->Elt(SizeM-1);

        unadjDate.AddMonths(NbMonth, GOTO_END_OF_MONTH);

        dateCour = unadjDate;
        dateCour.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);
    }
 
 
    // Calcul du premier flux (a la date spotday)

    double flow0 = exp(-spotRate/100.*matu1/K_YEAR_LEN);
    double flowSum;
    double flowCour = 0;
    double discountRate, tmp;
    double term3;
 

    for (i = 0; i < SizeS; i++)
    {
        // Recherche du spread du swap

        swapMatu = ARM_Date(sprdeeDates->Elt(SizeM + i));

        spread = ZCSpread->DiscountYield((swapMatu-asofdate)/K_YEAR_LEN, -1);


        // Somme des flux intermediaires

        unadjDate = spotDate;
        resetDate = unadjDate;

        unadjDate.AddMonths(NbMonth, GOTO_END_OF_MONTH);
        zeroTerm  = unadjDate;
        zeroTerm.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);
 
        flowSum = 0;
 
        while (zeroTerm < swapMatu)
        {
            discountRate = linInterpol(sprdeeDates, zeroTerm.GetJulian(),
                                       sprdeeValues);

            fwdRate = zcInit->ForwardYield(resetDate, zeroTerm, -1);

            fwdRate *= (libDayCount/K_YEAR_LEN);

            term3 = CountYears(ccy->GetLiborIndexDayCount(), 
                               resetDate, zeroTerm);

            flowCour = (fwdRate-spread/100.0)/100.0
                       *term3*
                       exp(-discountRate/100.*(zeroTerm-asofdate)/K_YEAR_LEN);
 
            flowSum += flowCour;

            resetDate = zeroTerm;
            unadjDate.AddMonths(NbMonth, GOTO_END_OF_MONTH);
            zeroTerm = unadjDate;
            zeroTerm.GoodBusinessDay(K_MOD_FOLLOWING, ccyName);
        }
 
        tmp = flow0-flowSum;

        fwdRate = zcInit->ForwardYield(resetDate, zeroTerm, -1);
        fwdRate *= (libDayCount/K_YEAR_LEN);

        term3 = CountYears(ccy->GetLiborIndexDayCount(), 
                               resetDate, zeroTerm);

        tmp /= (1+((fwdRate-spread/100.0)/100.0
               *term3));
 
        sprdeeValues->Elt(SizeM+i) = -K_YEAR_LEN/(zeroTerm-asofdate)*log(tmp);
        sprdeeValues->Elt(SizeM+i) *= 100;
    }
 
    // Construction du vecteur DateTerm & YearTerm

    ARM_Vector* YearTerm = new ARM_Vector(SizeM+SizeS);
    ARM_Vector* DateTerm = new ARM_Vector(SizeM+SizeS);
 
    for (i = 0; i < SizeM+SizeS; i++)
    {
       DateTerm->Elt(i) = (sprdeeDates->Elt(i)-asofdate.GetJulian());
       YearTerm->Elt(i) = (sprdeeDates->Elt(i)-asofdate.GetJulian())/K_YEAR_LEN;
    }
 
    SetYearTerms(YearTerm);
    SetDateTerms(DateTerm);
    SetZeroRates(sprdeeValues);


    if (sprdeeDates)
    {
       delete sprdeeDates;
    } 

    if (MMValues)
    {
       delete MMValues;
    }
 
    if (MMDates)
    {
       delete MMDates; 
    }

    // to be consistent with the naive constructor
    if ( itsLastBucketInt == 0 )
    {
       ARM_Vector* yearTerms  = GetYearTerms();

       ARM_Vector* zeroYields = GetZeroRates(); 

       int size = yearTerms->GetSize();

       ARM_Vector* terms = new ARM_Vector(yearTerms, size+1, 0, size-1);

       terms->Elt(size) = 1000.0;

       ARM_Vector* yields = new ARM_Vector(zeroYields, size+1, 0, size-1);

       yields->Elt(size) = zeroYields->Elt(size-1);

       SetYearTerms(terms);

       SetZeroRates(yields);

     
       ARM_Vector* dateTerms = GetDateTerms(); 

       ARM_Vector* dates = new ARM_Vector(dateTerms, size+1, 0, size-1);

       dates->Elt(size) = dateTerms->Elt(size-1);

       SetDateTerms(dates);
    }

    GenerateDiscountFactors(K_ZEROCOUPON);
}



ARM_ZeroLInterpol::ARM_ZeroLInterpol(const ARM_ZeroLInterpol& zeroLInterpol)
                   :ARM_ZeroCurve(zeroLInterpol)
{
    Init();

    SetName(ARM_ZERO_LIN_INTERPOL);

    BitwiseCopy(&zeroLInterpol);
}



ARM_ZeroLInterpol& ARM_ZeroLInterpol::operator= (const ARM_ZeroLInterpol&
                                                   zeroLInterpol)
{
    (*this).ARM_ZeroCurve::operator = (zeroLInterpol);

    BitwiseCopy(&zeroLInterpol);

    return(*this);
}



/*----------------------------------------------------------------------------*
   Returns the discount price with maturity yearTerm years from settlement, 
    computed from splines.
   We assume that the interpolation is constant, out of the range of YearTerms
*----------------------------------------------------------------------------*/

double ARM_ZeroLInterpol::DiscountFunction(double yearTerm)
{
    double z, zeroShift=0.0, intYield;
    
    int    zrIsCloned = 0;

    ARM_Vector* zeroRates = NULL;

    if ( yearTerm < 0.0 )
    {
       // throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
       //     "YearTerm must be non negative");
        
       return 0.0; // Past flows have a PV of zero as they are in the
                   // accrued bucket.
    }

    if (GetBucketEndPeriod() > K_DOUBLE_TOL
        && yearTerm >= GetBucketStartPeriod() 
        && yearTerm < GetBucketEndPeriod())
    {
       zeroShift = 0.0001 * GetBPShift();
    }
    

    if ( GetBPShift() > 0 )
    {
       zeroRates = (ARM_Vector *) GetZeroRates()->Clone();

       zrIsCloned = 1;
    }
    else
    {
       zeroRates = (ARM_Vector *) GetZeroRates();
    }

/*  
    On ajoute la perturbation avant l'interpolation

    *zeroRates += *GetBPShifts();
*/


/*  interpolate */

    double* zeroRatesElt = zeroRates->GetElt();
    double* yearTermsElt = GetYearTerms()->GetElt();
    int lastIndx = GetYearTerms()->GetSize()-1;

    if ( itsLastBucketInt == 0 )
       lastIndx--;

    if ( itsInterpolMeth == K_LIN_EXTRAPOL )
    {
       intYield = linInterpol(yearTermsElt, lastIndx+1, yearTerm, zeroRatesElt)
                              /100.0;
    }
    else if ( yearTerm < yearTermsElt[0] - K_DOUBLE_TOL ) 
    {
       intYield = zeroRatesElt[0]/100.0;
    }
    else if ( yearTerm >= yearTermsElt[lastIndx] - K_DOUBLE_TOL) 
    {
       intYield = zeroRatesElt[lastIndx]/100.0;
    }
    else if (itsInterpolMeth != K_CONTINUOUS)
    {
       intYield = linInterpol(yearTermsElt, lastIndx+1, yearTerm, zeroRatesElt)
                   / 100.0;
    }
    else
    {
       int index;

       index = indexBeforeValue(GetYearTerms(), yearTerm);

       double t1, t2, DF1, DF2;

       t1 = yearTerm - yearTermsElt[index];

       t2 = yearTermsElt[index+1] - yearTermsElt[index];

       DF1 = GetDiscountFactors()->Elt(index);

       DF2 = GetDiscountFactors()->Elt(index+1);

       z = pow(DF1, 1.0-t1/t2)*pow(DF2, t1/t2);

       intYield = -log(z)/yearTerm;
    }

    if ( itsCompoundMeth == 0 )
    {
       z = exp(-yearTerm*(intYield + zeroShift));
    }
    else if ( itsCompoundMeth == -1 )
    {
       z = 1.0/(1.0+yearTerm*(intYield + zeroShift));
    }
    else if ( itsCompoundMeth > 0 )
    {
       z = pow(1.0+(zeroShift + intYield)/itsCompoundMeth,
                 -yearTerm*itsCompoundMeth);
    }

    if (zeroRates && zrIsCloned)
       delete zeroRates;

    return(z);
}



/*----------------------------------------------------------------------------*
   Returns the discount price with maturity yearTerm years from settlement, 
    computed from splines.
   We assume that the interpolation is constant, out of the range of YearTerms
*----------------------------------------------------------------------------*/

double ARM_ZeroLInterpol::ComputeSlope(double yearTerm)
{
    double z, zeroShift=0.0, intYield;

    int    zrIsCloned = 0;

    ARM_Vector* zeroRates = NULL;

    if ( yearTerm < 0.0 )
    {
       // throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
       //     "YearTerm must be non negative");
        
       return 0.0; // Past flows have a PV of zero as they are in the
                   // accrued bucket.
    }

    if (GetBucketEndPeriod() > K_DOUBLE_TOL
        && yearTerm >= GetBucketStartPeriod() 
        && yearTerm < GetBucketEndPeriod())
    {
       zeroShift = 0.0001 * GetBPShift();
    }
    

    if ( GetBPShift() > 0 )
    {
       zeroRates = (ARM_Vector *) GetZeroRates()->Clone();

       zrIsCloned = 1;
    }
    else
    {
       zeroRates = (ARM_Vector *) GetZeroRates();
    }

/*  
    On ajoute la perturbation avant l'interpolation

    *zeroRates += *GetBPShifts();
*/


/*  interpolate */

    double* zeroRatesElt = zeroRates->GetElt();
    double* yearTermsElt = GetYearTerms()->GetElt();
    int lastIndx = GetYearTerms()->GetSize()-1;

    if ( itsLastBucketInt == 0 )
       lastIndx--;

    if ( yearTerm < yearTermsElt[0] - K_DOUBLE_TOL ) 
    {
       intYield = SlopelinInterpol(yearTermsElt, lastIndx+1, yearTerm+0.001, zeroRatesElt);
    }
    else if ( yearTerm >= yearTermsElt[lastIndx] - K_DOUBLE_TOL) 
    {
       intYield = 0.;
    }
    else if (itsInterpolMeth != K_CONTINUOUS)
    {
       intYield = SlopelinInterpol(yearTermsElt, lastIndx+1, yearTerm, zeroRatesElt);
    }
    else
    {
       int index;

       index = indexBeforeValue(GetYearTerms(), yearTerm);

       double t1, t2, DF1, DF2;

       t1 = yearTerm - yearTermsElt[index];
       t2 = yearTermsElt[index+1] - yearTermsElt[index];

       DF1 = GetDiscountFactors()->Elt(index);
       DF2 = GetDiscountFactors()->Elt(index+1);

       z = pow(DF1, 1.0-t1/t2)*pow(DF2, t1/t2);

       intYield = -log(z)/yearTerm;
    }

    if (zeroRates && zrIsCloned)
       delete zeroRates;

    return(intYield);
}



/*----------------------------------------------------------------------------*
    Returns the d(discount price) / d(yearTerm) with maturity yearTerm 
    years from settlement, computed from linear interpolation.
*----------------------------------------------------------------------------*/

double ARM_ZeroLInterpol::D1DiscountFunction(double yearTerm)
{
    double zp;

    if ( yearTerm < 0.0 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                       "YearTerm must be non negative");
    }

    zp = (DiscountFunction(yearTerm+0.000001)
        - DiscountFunction(yearTerm-0.000001))*500000.0; 
    
    return (zp);
}


    
/*----------------------------------------------------------------------------*
    Returns the d2(discount price) / d(yearTerm)2 with maturity yearTerm 
    years from settlement, computed from splines.
*----------------------------------------------------------------------------*/
double ARM_ZeroLInterpol::D2DiscountFunction(double yearTerm)
{
    double zp;

    if ( yearTerm < 0.0 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "YearTerm must be non negative");
    }

    zp = (D1DiscountFunction(yearTerm+0.000001)
        - D1DiscountFunction(yearTerm-0.000001))*500000.0; 

    return (zp);
}



/*------------------------ TAM/T4M Curve ------------------------------------*/


ARM_ZeroLInterpol::ARM_ZeroLInterpol(ARM_Date& asOf,
                                     char* terms[ARM_NB_TERMS],
                                     ARM_Vector* mktData,
                                     double mean_rates,
                                     int raw,
                                     int Cont_Lin,
                                     int lastBucketInt,
                                     ARM_Currency* ccy) : ARM_ZeroCurve(asOf,
                                                           terms,
                                                           mktData,
                                                           mean_rates,
                                                           raw,
                                                           Cont_Lin,
                                                           ccy)
 
{
    // Check variables

    SetName(ARM_ZERO_LIN_INTERPOL);

    // Set variables

    if (!lastBucketInt)
    {
       int size = GetYearTerms()->GetSize();

       ARM_Vector* terms = new ARM_Vector(GetYearTerms(), size+1, 0, size-1);

       terms->Elt(size) = 1000.0;

       ARM_Vector* yields = new ARM_Vector(GetZeroRates(), size+1, 0, size-1);

       yields->Elt(size) = GetZeroRates()->Elt(size-1);

       ARM_Vector* df = new ARM_Vector(GetDiscountFactors(), size+1, 0,
                                        size-1);

       df->Elt(size) = exp(-10.0*GetZeroRates()->Elt(size-1));

       SetYearTerms(terms);

       SetZeroRates(yields);

       SetDiscountFactors(df);
    }

    int Size = GetYearTerms()->GetSize();

    ARM_Vector* BPShifts = new ARM_Vector(Size, 0.0);

    for (int i = 0; i < Size; i++)
    {
        if ( GetYearTerms()->Elt(i) >= GetBucketStartPeriod()
                && GetYearTerms()->Elt(i) < GetBucketEndPeriod())
        {
            BPShifts->Elt(i) = GetBPShift()*0.01; // using yield in percent
        }
    }

    SetBPShifts(BPShifts);

    itsCompoundMeth  = K_COMP_CONT;

    itsInterpolMeth = Cont_Lin;

    itsLastBucketInt = lastBucketInt;
}


ARM_ZeroLInterpol::ARM_ZeroLInterpol(ARM_ZeroCurve& inZeroCurve)
{
    int i;

    Init();

    ARM_ZeroCurve::Copy((const ARM_Object*)&inZeroCurve);


    switch(inZeroCurve.GetName())
    {
        case ARM_ZERO_FLAT:
        {
            double dYearTerms[]={1,100,500};
            ARM_Vector *yearTerms = new ARM_Vector(3,dYearTerms);
            ARM_Vector* zeroYields = new ARM_Vector(3,0.0);

            for (i=2; i>=0;i--)
            {
                (*zeroYields)[i] = inZeroCurve.DiscountYield((double)
                                                             (*yearTerms)[i],
                                                             itsCompoundMeth);
            }

            if (itsLastBucketInt)
            {
               SetYearTerms(yearTerms);

               SetZeroRates(zeroYields);
            }
            else
            {
               int size = yearTerms->GetSize();

               ARM_Vector* terms = new ARM_Vector(yearTerms, size+1,
                                                   0, size-1);

               terms->Elt(size) = 1000.0;

               ARM_Vector* yields = new ARM_Vector(zeroYields, size+1,
                                                    0, size-1);

               yields->Elt(size) = zeroYields->Elt(size-1);

               SetYearTerms(terms);

               SetZeroRates(yields);
            }

            int Size = GetYearTerms()->GetSize();

            ARM_Vector* BPShifts = new ARM_Vector(Size, 0.0);

            for ( i = Size-1; i >=0; i--)
            {
                if ( GetYearTerms()->Elt(i) >= GetBucketStartPeriod()
                     && GetYearTerms()->Elt(i) < GetBucketEndPeriod())
                {
                    BPShifts->Elt(i)=GetBPShift()*0.01; //using yield in percent
                }
            }

            SetBPShifts(BPShifts);

            GenerateDateTerms();

            GenerateDiscountFactors(itsCompoundMeth);
        }
        break;


        case ARM_ZERO_VASICEK:
            break;


        default :
        {}
    }

}


double ARM_ZeroLInterpol::nvDiscountPrice(ARM_Date& maturity)
{
    double yearTerm;

    double res = 0.0;

    
    yearTerm = (maturity.GetJulian()-GetAsOfDate().GetJulian())/K_YEAR_LEN;

    res = ARM_ZeroLInterpol::DiscountFunction(yearTerm);

    return(res);
}


double ARM_ZeroLInterpol::nvDiscountPrice(double yearTerm)
{
    return(ARM_ZeroLInterpol::DiscountFunction(yearTerm));
}




/*----------------- TOY METHODS -----------------*/


ARM_ZeroLInterpol::ARM_ZeroLInterpol(ARM_Date& asOf,
                                     ARM_CRV_TERMS& terms,
                                     ARM_Vector* mktData,
                                     int MMVsFut,
                                     int SwapVsFut,
                                     int raw,
                                     int Cont_Lin,
                                     int lastBucketInt,
                                     int Frequency,
                                     ARM_Currency* ccy) : ARM_ZeroCurve(asOf,
                                                           terms,
                                                           mktData,
                                                           MMVsFut,
                                                           SwapVsFut,
                                                           raw,
                                                           Cont_Lin,
                                                           Frequency,
                                                           ccy)
{
    Init();

    // Set variables

    if (!lastBucketInt)
    {
       int size = GetYearTerms()->GetSize();

       ARM_Vector* terms = new ARM_Vector(GetYearTerms(), size+1, 0, size-1);

       terms->Elt(size) = 1000.0;

       ARM_Vector* yields = new ARM_Vector(GetZeroRates(), size+1, 0, size-1);

       yields->Elt(size) = GetZeroRates()->Elt(size-1);

       ARM_Vector* df = new ARM_Vector(GetDiscountFactors(), size+1,
                                        0, size-1);

       df->Elt(size) = exp(-10.0*GetZeroRates()->Elt(size-1));

       SetYearTerms(terms);

       SetZeroRates(yields);

       SetDiscountFactors(df);
    }

    int Size = GetYearTerms()->GetSize();

    ARM_Vector* BPShifts = new ARM_Vector(Size, 0.0);

    for (int i = 0; i < Size; i++)
    {
        if ( GetYearTerms()->Elt(i) >= GetBucketStartPeriod()
                && GetYearTerms()->Elt(i) < GetBucketEndPeriod())
        {
           BPShifts->Elt(i) = GetBPShift()*0.01; // using yield in percent
        }
    }

    SetBPShifts(BPShifts);

    itsCompoundMeth = K_COMP_CONT;

    itsInterpolMeth = Cont_Lin;

    itsLastBucketInt = lastBucketInt;
}

/*------------- END OF TOY METHODS --------------*/


ARM_ZeroCurve* ARM_ZeroLInterpol::GenerateShiftCurve(ARM_CRV_TERMS& Term, ARM_Vector* epsilon)
{
   if (GetMktData() == NULL)
   {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "The input curve cant be shifted");
   }

   ARM_ZeroLInterpol* zc=NULL;
   ARM_Date tmpdat;
   ARM_MarketData* MktData = (ARM_MarketData*) GetMktData()->Clone();

   // verification de l'existence des datamarkets
   if (!MktData)
      return NULL;

   int szinp =0;
   szinp = epsilon->GetSize();

   int k=0;
   int l=0;

   for (k=0; k<ARM_NB_TERMS; k++)
   {
      for (l=0; l<szinp; l++)
      {
         if (!strcmp(Term[l],MktData->itsMktTerms[k]))
         {
            // Cas Futur
            if ( MktData->itsMktValue->Elt(k) > 50.0 )
               MktData->itsMktValue->Elt(k) = 
                           MktData->itsMktValue->Elt(k) - epsilon->Elt(l);
            else
               MktData->itsMktValue->Elt(k) = 
                           MktData->itsMktValue->Elt(k) + epsilon->Elt(l);

            l=szinp;
         }
      }
   }

   switch (MktData->itsConstructionMeth)
   {
       case 0 :
               zc = new ARM_ZeroLInterpol(GetAsOfDate(),
                        MktData->itsMktTerms,
                        MktData->itsMktValue,
                        MktData->itsMMVsFut,
                        MktData->itsSwapVsFut,
                        MktData->itsraw,
                        MktData->itsCont_Lin,
                        this->itsLastBucketInt,
                        GetCurrencyUnit());
       break;

       case 1 :
               zc = new ARM_ZeroLInterpol(MktData->itsMktTerms,
                        GetAsOfDate(),
                        MktData->itsMktValue,
                        MktData->itsMMVsFut,
                        MktData->itsSwapVsFut,
                        MktData->itsraw,
                        MktData->itsCont_Lin,
                        this->itsLastBucketInt,
                        GetCurrencyUnit());
       break;

       case 2 :
               zc = new ARM_ZeroLInterpol(GetAsOfDate(), 
                        MktData->itsMktTerms, 
                        MktData->itsMktValue, 
                        MktData->itsMMVsFut, 
                        MktData->itsSwapVsFut, 
                        MktData->itsraw, 
                        MktData->itsCont_Lin,
                        this->itsLastBucketInt,
                        0,
                        GetCurrencyUnit());
       break;

       case 3 :
               zc = new ARM_ZeroLInterpol(GetAsOfDate(),
                                          MktData->itsMktValue,     
                                          MktData->itsMktTerms, 
                                          MktData->itsMMVsFut, 
                                          MktData->itsSwapVsFut, 
                                          MktData->itsCont_Lin,
                                          GetCurrencyUnit());
       break;

       default:
       {
           if (zc)
              delete zc;

           zc = NULL;

           delete MktData;

           return(NULL);
       }
   }

   delete MktData;

   return(zc);
}






/*---------------------------------------------------------------------------*/
/*---- End of file ----*/