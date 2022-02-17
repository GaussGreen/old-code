/*
 * $Log: fromto.cpp,v $
 * Revision 1.23  2004/03/19 19:37:10  jpriaudel
 * include swap.h in fromto.cpp
 * (compilation pb)
 *
 * Revision 1.22  2004/03/18 15:32:37  mab
 * added:
 * FromAbsoluteVegaToRelativeVega, MakeGenSwap
 * ApproximatedForwardRate
 *
 * Revision 1.21  2004/02/02 10:29:18  jpriaudel
 * GetFrqAndNbFromString added
 *
 * Revision 1.20  2003/12/08 15:58:30  jpriaudel
 * modif StringToGenMatu
 *
 * Revision 1.19  2003/11/28 18:18:40  jpriaudel
 * modif dans FromStringToIndexType
 *
 * Revision 1.18  2003/11/14 11:12:48  mab
 * Added : CalcFwdFXSpot(..
 *
 * Revision 1.17  2003/10/23 17:00:15  ykhlif
 *  ajout de ConvertCouponAsSummit()
 *
 * Revision 1.16  2003/09/30 10:23:08  emezzine
 * change the precision for 1.0e-2 to 1.0e-1 in line 757
 *
 * Revision 1.15  2003/09/17 18:04:12  ebenhamou
 * added reverse operation of StringMatuToYearTerm
 *
 * Revision 1.14  2003/09/16 11:10:53  jpriaudel
 * Ajout de FromSummitGapToARMGap
 *
 * Revision 1.13  2003/09/09 12:51:06  jpriaudel
 * ajout de FromSummitFreqToARMFreq et FromSummitDaycountToARMDaycount
 *
 * Revision 1.12  2003/09/02 15:43:02  jpriaudel
 * ajout de GetSummitCurveIndexFromCurrency
 *
 * Revision 1.11  2003/08/25 09:02:37  jpriaudel
 * Added FromStringToIndexType
 *
 * Revision 1.10  2003/07/22 18:08:15  jpriaudel
 * correction dans FromRateToRate
 *
 * Revision 1.9  2002/03/13 10:15:07  mab
 * *** empty log message ***
 *
 * Revision 1.8  2001/06/25 15:24:41  vberger
 * ajout de FromLiborTypeToFrequency
 *
 * Revision 1.7  2001/02/15 18:36:50  nicolasm
 * Ajout FromFrequencyToXiborType
 *
 * Revision 1.6  1999/11/19 16:49:54  nicolasm
 * Ajout de fonctions de converstion de date generique et de Product type
 *
 * Revision 1.5  1999/09/10 14:24:34  nicolasm
 * Ajout FromTOYCodeToDates
 *
 * Revision 1.4  1999/02/18 18:45:52  nicolasm
 * Ajout fonction FromFutureCodeToDates
 *
 * Revision 1.3  1999/02/08 16:14:47  nicolasm
 * Ajout fonction FromReuterCodeToMonth
 *
 */

/*----------------------------------------------------------------------------*
 
    fromto.cpp
 
    Sources of functions computing the equivalent Rate or zero price

    using various compounding methods
 
    Copyright (c) 1997 

    term = rate or zero maturity in year term already computed
            using the right basis

    compMeth =   Rate Compounding methods 

 
    Continuous        K_COMP_CONT         0
    Proportional      K_COMP_PROP         -1
    Annual            K_COMP_ANNUAL       1
    Semiannual        K_COMP_SEMIANNUAL   2
    Quarterly         K_COMP_QUARTERLY    4
    Monthly           K_COMP_MONTHLY      12



*----------------------------------------------------------------------------*/

#include <math.h>

#include "fromto.h"
#include "armglob.h"
//#include "swap.h"

#include "gpbase\gpvector.h"


/*----------------------------------------------------------------------------*/

double FromStrMatuToDouble(const char* matu, ARM_Date* asOfDate)
{
    double resMatu = 0.0;
    int inMatu;

    char buf[20];
  


    buf[0] = ' ';

    int    month;
    int    year;
    ARM_Date    vFutureMaturity;

    if( isalpha(matu[0]) )    // Maturité d'un Future
    {
        if( ! asOfDate )
        {
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "asOfDate should not be null");
        }

        GetMonthYearFromExpiryDate((char*)matu, &month, &year);
        vFutureMaturity.ChangeDate(1, month, year);
        vFutureMaturity.PiborDelivery();

        resMatu = ( vFutureMaturity.GetJulian() - asOfDate->GetJulian() ) / 365.0;
    }
    else
    {
        sscanf(matu, "%d%s", &inMatu, buf);

        if (( toupper(buf[0]) == 'J' ) || ( toupper(buf[0]) == 'D' ))
        {
           resMatu = inMatu/365.0;
        }
        else if (( toupper(buf[0]) == 'S' ) || ( toupper(buf[0]) == 'W' ))
        {
           resMatu = inMatu*7.0/365.0;
        }
        else if ( toupper(buf[0]) == 'M' )
        {
           resMatu = inMatu/12.0;
        }
        else if (( toupper(buf[0]) == 'A' ) || ( toupper(buf[0]) == 'Y' ))
        {
           resMatu = inMatu;
        }
        else
        {
           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                 "===> Expiry or Maturity must look like : 3D, 2M, 5Y, ...");
        }
    }

    return(resMatu);
}



ARM_Swap* MakeGenSwap(ARM_Date& curDate,
                      double expiry,
                      double matu,
                      ARM::ARM_Currency* ccy,
					  long daycountId)
{
 //   ARM_Date swapStartDate;
 //   ARM_Date swapEndDate;
 //   ARM_Date tmpDate;


 //   ARM_Date stDate = curDate;

 //   char* ISOccy = ccy->GetCcyName();

 //   int spotDays = ccy->GetSpotDays();

 //   char payCal[50];

 //   ccy->CalcFloatPayCal(payCal);

 //   swapStartDate = curDate;

 //   swapStartDate.NextBusinessDay(spotDays, payCal);
	//
	//swapStartDate.AddDays(int(expiry*365.0));
	//
	//tmpDate = swapStartDate;
	//
	//swapEndDate = tmpDate.AddDays(int(matu*365.0));

 //   double fixedRate = -1.0; // 10.0;
 //   double spread    = 0.0;

 //   int resetFreq = -1;
 //   int payFreq   = -1;

 //   double swapPrice = 0.0;
 //   double swapRate  = 0.0;


 //   ARM_INDEX_TYPE CCY_INDEX = GetDefaultIndexFromCurrency(ISOccy);

 //   ARM_Swap* curSwap = new ARM_Swap(swapStartDate, swapEndDate,
 //                                    (ARM_INDEX_TYPE) CCY_INDEX,
 //                                    spread, fixedRate, K_RCV,
 //                                    resetFreq, payFreq, ccy, daycountId);
 //   return(curSwap);
	return NULL;
}



double ApproximatedForwardRate(double Mat,
                               double Tenor,
                               ARM_ZeroCurve* zc,
                               int isSwapRate,
							   long daycountId)
{

	//double curMatu   = Tenor;
 //   double curExpiry = Mat;
	//
	//int spotDays = zc->GetCurrencyUnit()->GetSpotDays();
	//
	//ARM_Date d_Matu = zc->GetAsOfDate();
	//ARM_Date d_Expiry = zc->GetAsOfDate();
	//ARM_Date tmpDate;
	//
	//char payCal[50];
 //   zc->GetCurrencyUnit()->CalcFloatPayCal(payCal);


 //   if (( curMatu == 0.0 ) || ( curMatu < 0.0 ))
 //   {
 //      throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
 //        "===>ApproximatedForwardRate: Tenor must be not null and > 0");
 //   }

 //   if ( curExpiry < 0.0 )
 //   {
 //      throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
 //            "ApproximatedForwardRate: ===> Maturity must be null or > 0");
 //   }

 //   char* ISOccy      = zc->GetCurrencyUnit()->GetCcyName();

 //   double fwdYield;

 //   if (!(isSwapRate)
 //       ||
 //       ( curMatu < 1 )
 //      )
 //   {
 //      /*--- Compute the Fwd Rate ---*/
	//			
	//   if ( daycountId < 0 )
	//   {
 // 		   if ((strcmp(ISOccy,"ZAR") != 0)
	//			 && (strcmp(ISOccy,"PLN") != 0)
	//			 && (strcmp(ISOccy,"BEF") != 0)
	//			 && (strcmp(ISOccy,"GBP") != 0)
	//			 && (strcmp(ISOccy,"PTE") != 0)
	//			 && (strcmp(ISOccy,"AUD") != 0)
	//			)
	//		{
	//			fwdYield = zc->ForwardYield(curExpiry,
	//									  curExpiry+curMatu, -1)*360.0/365.0;
	//		}
	//		else
	//		{
	//			fwdYield = zc->ForwardYield(curExpiry,
	//									  curExpiry+curMatu, -1);
	//		}
	//   }
	//   else
	//   {

	//	d_Expiry.NextBusinessDay(spotDays, payCal);
	//
	//	d_Expiry.AddDays(int(curExpiry*365.0));
	//
	//	tmpDate = d_Expiry;
	//
	//	d_Matu = tmpDate.AddDays(int(curMatu*365.0));

	//	fwdYield = zc->ForwardYieldWithDayCount(d_Expiry, d_Matu, -1, (int) daycountId);
	//   }
 //   }
 //   else
 //   {
 //      /*--- Compute the Swap Fwd Rate ---*/

 //      ARM_Swap* curSwap = MakeGenSwap(zc->GetAsOfDate(),
 //                                      curExpiry,
 //                                      curMatu,
 //                                      zc->GetCurrencyUnit(),
	//								   daycountId);

 //      /*double swapPrice = 0.0;
 //      double swapRate  = 0.0;

 //      ARM_YCModel theModel(zc);

 //      curSwap->SetModelVariable(NULL);
 //      curSwap->SetModel(&theModel);

 //      swapRate = curSwap->PriceToRate(zc->GetAsOfDate(), swapPrice);

 //      fwdYield = swapRate;*/

 //      delete curSwap;
 //   }

 //   return(fwdYield);
	return 0;
}



double FromAbsoluteVegaToRelativeVega(char* expiry,
                                      char* matu,
                                      ARM_VolCurve* vol,
                                      ARM_ZeroCurve* zc,
                                      double absoluteVega)
{
    double relativeVega = 0.0;

    double convExpiry;
    double convMatu;

    convExpiry = FromStrMatuToDouble(expiry);
    convMatu   = FromStrMatuToDouble(matu);

    double fwd = ApproximatedForwardRate(convExpiry,
                                         convMatu,
                                         zc,
                                         1)/100.0;

    double volBS = 0;//vol->ComputeVolatility(convExpiry, convMatu)/100.0;

    relativeVega = absoluteVega*fwd*exp(-volBS*volBS*convExpiry/8.0); 

    return(relativeVega);
}



double FromRateToRate(double rate, double term,
                      int fromCompMeth, int toCompMeth)
{
    double newRate, Rate;

    Rate = rate/100.0;

    switch(fromCompMeth) 
    {
        case K_COMP_ANNUAL :
        case K_COMP_SEMIANNUAL :
        case K_COMP_QUARTERLY :
        case K_COMP_BIMONTHLY :
        case K_COMP_MONTHLY :
        {
            switch (toCompMeth) 
            {
                case K_COMP_ANNUAL :
                case K_COMP_SEMIANNUAL :
                case K_COMP_QUARTERLY :
                case K_COMP_BIMONTHLY :
                case K_COMP_MONTHLY :
                {
                    newRate = 
                        pow(1.0+(Rate / double(fromCompMeth)), double(fromCompMeth)/double(toCompMeth));

                    newRate -= 1.0;

                    newRate *= double(toCompMeth);
                    break;
                };

                case K_COMP_CONT :
                {
                    newRate = log(1.0+Rate/double(fromCompMeth));
                    newRate *= double(toCompMeth);
                    break;
                };

                case K_COMP_PROP :
                {
                    newRate = 
                        pow(1.0+Rate/double(fromCompMeth), 
                            double(fromCompMeth)*term);

                    newRate -= 1.0;

                    newRate /= term;
                    break;
                };

                default :
                    newRate = Rate;
            }
            break;
        }

        case K_COMP_CONT :
        {
            switch (toCompMeth) 
            {
                case K_COMP_ANNUAL :
                case K_COMP_SEMIANNUAL :
                case K_COMP_QUARTERLY :
                case K_COMP_BIMONTHLY :
                case K_COMP_MONTHLY :
                {
                    newRate = exp(Rate/double(toCompMeth)) - 1.0;

                    newRate *= double(toCompMeth);
                    break;
                };

                case K_COMP_PROP :
                {
                    newRate = exp(Rate*term) - 1.0;

                    newRate /= term;
                    break;
                };
                default :
                    newRate = Rate;
            };
            break;
        };

        case K_COMP_PROP :
        {
            switch (toCompMeth) 
            {
                case K_COMP_ANNUAL :
                case K_COMP_SEMIANNUAL :
                case K_COMP_QUARTERLY :
                case K_COMP_BIMONTHLY :
                case K_COMP_MONTHLY :
                {
                    newRate = 
                        pow(1.0+Rate*term, 1.0/double(term*toCompMeth));

                    newRate -= 1.0;

                    newRate *= double(toCompMeth);
                    break;
                };

                case K_COMP_CONT :
                {
                    newRate = log(1.0+Rate*term);

                    newRate /= term;
                    break;
                };
                default :
                    newRate = Rate;
            };
            break;
        };
        default :
            newRate = Rate;
    }
    return(100.0*newRate);
}



/*-----------------------------------------------------------------------------*

  Compute the equivalent rate in toRateType term from fromRateType term

  possible values of toRateType and fromRateType are  

    Rate Type               dayCount
    
    K_MM_RATE   =  -1       K30_360         interest = d*R
    K_ACT_RATE   =   1      KACTUAL_ACTUAL  interest = (1+R)^d - 1
    K_CONT_RATE =   0       KACTUAL_365     interest = exp(d*R)- 1

  where interest period are computed using the dayCount column Basis 

*-----------------------------------------------------------------------------*/

double FromRateToRate(double rate, 
                      ARM_Date& startDate, ARM_Date& endDate,
                      int fromRateType, int toRateType)
{
    double newRate, Rate, fromTerm, toTerm;

    Rate = rate/100.0;
    
    switch (fromRateType) 
    {
        case K_MM_RATE :
        {
            fromTerm = CountYears(KACTUAL_360, startDate, endDate);
            
            switch (toRateType) 
            {
                case K_ACT_RATE :
                {
                    toTerm = CountYears(KACTUAL_ACTUAL, startDate, endDate);
                    
                    newRate = pow(1.0+fromTerm*Rate, 1.0/double(toTerm))-1.0;

                    newRate *= 100.0;
                    break;
                };

                case K_CONT_RATE :
                {
                    toTerm = CountYears(KACTUAL_365, startDate, endDate);
                    
                    newRate = log(1.0+fromTerm*Rate)/double(toTerm);
                    
                    newRate *= 100.0;
                    break;
                };

                default :
                    newRate = rate;
            }
            break;
        }

        case K_ACT_RATE :
        {
            fromTerm = CountYears(KACTUAL_ACTUAL, startDate, endDate);
            
            switch (toRateType) 
            {
                case K_MM_RATE :
                {
                    toTerm = CountYears(KACTUAL_360, startDate, endDate);
                    
                    newRate = 
                        pow(1.0+Rate, double(fromTerm))-1.0;

                    newRate *= 100.0/double(toTerm);
                    break;
                };

                case K_CONT_RATE :
                {
                    toTerm = CountYears(KACTUAL_365, startDate, endDate);
                    
                    newRate = fromTerm*log(1.0+Rate);
                    
                    newRate *= 100.0/double(toTerm);
                    break;
                };
                default :
                    newRate = rate;
            }
            break;
        }

        case K_CONT_RATE :
        {
            fromTerm = CountYears(KACTUAL_365, startDate, endDate);
            
            switch (toRateType) 
            {
                case K_MM_RATE :
                {
                    toTerm = CountYears(KACTUAL_360, startDate, endDate);
                    
                    newRate = exp(fromTerm*Rate)-1.0;
                    
                    newRate *= 100.0/double(toTerm);
                    break;
                };

                case K_ACT_RATE :
                {
                    toTerm = CountYears(KACTUAL_ACTUAL, startDate, endDate);
                    
                    newRate = exp(fromTerm*Rate/toTerm)-1.0;
                    
                    newRate *= 100.0;
                    break;
                };
                default :
                    newRate = rate;
            }
            break;
        }

        default :
            newRate = rate;
    }

    return(newRate);
}



double FromZeroToRate(double zero, double term, int RateCompMeth)
{
    double Rate;

    switch (RateCompMeth) 
    {
        case K_COMP_ANNUAL :
        case K_COMP_SEMIANNUAL :
        case K_COMP_QUARTERLY :
        case K_COMP_BIMONTHLY :
        case K_COMP_MONTHLY :
        {
            Rate = pow(1.0/zero, 1.0/double(term*RateCompMeth));
            Rate -= 1.0;
            Rate *= RateCompMeth;
            break;
        }

        case K_COMP_PROP :
        {
            Rate = 1.0/zero - 1.0;
            Rate /= term;
            break;
        }
 
        case K_COMP_CONT :
        {
            Rate = -log(zero)/term;
            break;
        }

        default :
            Rate = -log(zero)/term;
    }
    return(100.0*Rate); 
}



double FromRateToZero(double rate, double term, int RateCompMeth)
{
    double Zero, Rate;

    Rate = rate/100.0;

    switch (RateCompMeth) 
    {
        case K_COMP_ANNUAL :
        case K_COMP_SEMIANNUAL :
        case K_COMP_QUARTERLY :
        case K_COMP_BIMONTHLY :
        case K_COMP_MONTHLY :
        {
            Zero = 
                pow(1.0+ Rate/double(RateCompMeth), -double(RateCompMeth*term));
            break;
        }

        case K_COMP_CONT :
        {
            Zero = exp(-Rate*term);
            break;
        }
 
        case K_COMP_PROP :
        {
            Zero = 1.0/(1.0+Rate*term);
            break;
        }

        default :
            Zero = exp(-Rate*term);
    }

    return(Zero); 
}



int FromFrequencyToLiborType(int freq)
{
    int libType;

    switch (freq) 
    {
        case K_MONTHLY :
        {
            libType = K_LIBOR1M;
            break;
        }

        case K_BIMONTHLY :
        {
            libType = K_LIBOR2M;
            break;
        }

        case K_QUARTERLY :
        {
            libType = K_LIBOR3M;
            break;
        }
 
        case K_SEMIANNUAL :
        {
            libType = K_LIBOR6M;
            break;
        }

        case K_ANNUAL :
        {
            libType = K_LIBOR1Y;
            break;
        }
        default :
            libType = K_LIBOR6M;
    }
    return(libType); 
}



int FromFrequencyToPiborType(int freq)
{
    int libType;

    switch (freq) 
    {
        case K_MONTHLY :
        {
            libType = K_PIBOR1M;
            break;
        }

        case K_BIMONTHLY :
        {
            libType = K_PIBOR2M;
            break;
        }

        case K_QUARTERLY :
        {
            libType = K_PIBOR3M;
            break;
        }
 
        case K_SEMIANNUAL :
        {
            libType = K_PIBOR6M;
            break;
        }

        case K_ANNUAL :
        {
            libType = K_PIBOR1Y;
            break;
        }
        default :
            libType = K_PIBOR6M;
    }

    return(libType); 
}



int FromFrequencyToEuriborType(int freq)
{
    int libType;
 
    switch (freq)
    {
        case K_MONTHLY :
        {
            libType = K_EURIBOR1M;
            break;
        }
 
        case K_BIMONTHLY :
        {
            libType = K_EURIBOR2M;
            break;
        }
 
        case K_QUARTERLY :
        {
            libType = K_EURIBOR3M;
            break;
        }
 
        case K_SEMIANNUAL :
        {
            libType = K_EURIBOR6M;
            break;
        }
 
        case K_ANNUAL :
        {
            libType = K_EURIBOR1Y;
            break;
        }
        default :
            libType = K_EURIBOR6M;
    }

    return(libType);
}


 
int FromFrequencyToXiborType(int freq, char *ccy)
{
    if (!strcmp("EUR",ccy))
       return(FromFrequencyToEuriborType(freq));
    else
       return(FromFrequencyToLiborType(freq));
}
 


int FromReuterCodeToMonth(char code)
{
     switch(code)
     {
         case 'F' : return 1; break;
         case 'G' : return 2; break;
         case 'H' : return 3; break;
         case 'J' : return 4; break;
         case 'K' : return 5; break;
         case 'M' : return 6; break;
         case 'N' : return 7; break;
         case 'Q' : return 8; break;
         case 'U' : return 9; break;
         case 'V' : return 10; break;
         case 'X' : return 11; break;
         case 'Z' : return 12; break;

         default :
         {
             throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                 "Invalid Reuter code");

             return 0;
         }
     }
}



void FromFutureCodeToDates(char* terms, ARM_Date& dateContrat,
                           ARM_Date& dateUnder)
{
    char mthChr, yearChr, matuChr;
 
 
    // Decodage du terms
 
    sscanf(terms, "%c%c%c", &mthChr, &yearChr, &matuChr);
 
    int mthInt  = FromReuterCodeToMonth(mthChr);
 
    int yearInt = atoi(&yearChr);
 
    if (yearInt == 9)
       yearInt = 1999;
    else
       yearInt += 2000;
 
 
    ARM_Date dateEch;
 
    dateEch.ChangeDate(1, mthInt, yearInt);
 
    // Si on retire GetDayOfWeek au 1er du mois on se place
    // au dimanche qui precede le premier du mois.
    //

    if ( dateEch.GetDayOfWeek() <= 3 )
    {
       dateEch.AddDays(17 - dateEch.GetDayOfWeek());
    }
    else
    {
       dateEch.AddDays(24 - dateEch.GetDayOfWeek());
    }
 
    dateContrat = dateEch;
 
 
    ARM_Date dateMatu = dateEch;
 
    if ( matuChr == 'M' )
       dateMatu.AddMonths(1);
    else
       dateMatu.AddMonths(3);
 
    dateMatu.ChangeDate(1, dateMatu.GetMonth(), dateMatu.GetYear());
 
 
    if ( dateMatu.GetDayOfWeek() <= 3 )
    {
       dateMatu.AddDays(17 - dateMatu.GetDayOfWeek());
    }
    else
    {
       dateMatu.AddDays(24 - dateMatu.GetDayOfWeek());
    }
 
 
    dateUnder = dateMatu;
}



void FromTOYCodeToDates(char* terms, ARM_Date& dateContrat,
                        ARM_Date& dateUnder)
{
    char car1, car2;
    int yearInt;
 
 
    // Decodage du terms (TY99 signifie passge 99 / 2000)
 
    sscanf(terms, "%c%c%d", &car1, &car2, &yearInt);
 
 
    if ( yearInt == 99 )
       yearInt = 1999;
    else
       yearInt += 2000;
 
    ARM_Date dateEch(31, 12, yearInt);
 
    dateEch.GoodBusinessDay(-1);

    dateContrat = dateEch;
 
    dateUnder = dateEch;
 
    dateUnder.NextBusinessDay();
}



ARM_PRODUCT_TYPE StringToProductType(char *str)
{
    ARM_PRODUCT_TYPE pt;
 
    if (!strcmp(str, "IRG"))
        pt = PT_IRG;
    else if (!strcmp(str,"SWOPT" ))
        pt = PT_SWOPT;
    else
        pt = PT_IRG;
 
    return pt;
}



ARM_GEN_MATU StringToGenMatu(char *str)
{
    ARM_GEN_MATU gm;

 
    if (!strcmp(str, "3M"))
       gm = GM_3M;
    else if (!strcmp(str, "6M"))
       gm = GM_6M;
    else if (!strcmp(str, "9M"))
       gm = GM_9M;
    else if (!strcmp(str, "12M"))
       gm = GM_1Y;
    else if (!strcmp(str, "1Y"))
       gm = GM_1Y;
    else if (!strcmp(str, "2Y"))
       gm = GM_2Y;
    else if (!strcmp(str, "5Y"))
       gm = GM_5Y;
    else if (!strcmp(str, "10Y"))
       gm = GM_10Y;
    else if (!strcmp(str, "15Y"))
       gm = GM_15Y;
    else if (!strcmp(str, "20Y"))
       gm = GM_20Y;
    else if (!strcmp(str, "25Y"))
       gm = GM_25Y;
    else if (!strcmp(str, "30Y"))
       gm = GM_30Y;
    else if (!strcmp(str, "40Y"))
       gm = GM_40Y;
    else
       gm = GM_40Y;
 
 
    return gm;
}


double StringMatuToYearTerm(char *str)
{
    double yt = 0.0;
    int nb;
    char unit;
 
    sscanf(str, "%d%c", &nb, &unit);
 
    unit = toupper(unit);
 
    switch (unit)
    {
        case 'D' : yt = nb / 365.0; break;
        case 'W' : yt = nb * 7.0 / 365.0; break;
        case 'M' : yt = nb / 12.0; break;
        case 'Y' : yt = nb * 1.0; break;
        default  : throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "invalid maturity string should be like nb[d/w/m/y]");
    }
 
    return yt;
}



/// returns true if manage to convert
/// else return false
bool YearTermToStringMatuHelper( double d, char* res, double period, const char* periodString )
{

    const double precision = 1e-1;

    /// priority to year then to month then week then to day to day
    if( fabs( period * d - int( period * d ) ) < precision )
    {
        sprintf( res, "%d%s", int( period  *d ), periodString );
        return true;
    }
    else
        return false;
}



string YearTermToStringMatu( double d)
{
	char cres[10];

    /// because of potential ambiguity 
    /// priority to year then to month then week then to day to day
	
	// except for 0D :
	if ( d == 0 )
        return string("0D");

    /// choice is that for more than 2 year
    /// we round to the corresponding year
    if( d > 2.0 )
    {
        sprintf( cres, "%dY", int( d ) );
        return (string) cres;
    }

    if( YearTermToStringMatuHelper( d, cres, 1.0, "Y" ) )
	{
        return (string) cres;
	}
    
    if( YearTermToStringMatuHelper( d, cres, 12.0, "M" ) )
	{
        return (string) cres;
	}
    
    if( YearTermToStringMatuHelper( d, cres, 365.0/7.0, "W" ) )
	{
        return (string) cres;
	}
    
    if( YearTermToStringMatuHelper( d, cres, 365.0, "D" ) )
	{
        return (string) cres;
	}
    
    /// if did not manage to convert throw an exception
    throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "Invalid double period");
}


string ConvertYearTermToStringMatu(double d)
{
	char cres[10];

	double	dNbMonths = d * 12.0;
	int		iNbMonths = d * 12.0;

	if ( fabs(dNbMonths - (double)iNbMonths) <= 15./30. )
	{
		if ( (iNbMonths % 12 == 0) )
			sprintf(cres,"%dY",iNbMonths / 12);
		else
			sprintf(cres,"%dM",iNbMonths);
	}
	else if (fabs(dNbMonths - (double)(iNbMonths+1)) < 15./30.)
	{
		if ( ((iNbMonths+1) % 12 == 0) )
			sprintf(cres,"%dY",(iNbMonths+1) / 12);
		else
			sprintf(cres,"%dM",iNbMonths+1);
	}

	return (string) cres;
}




/////////////////////////////////////////////////////////////////
/// split a given string into a vector of string according to a given delimiter
/////////////////////////////////////////////////////////////////
vector< string> splitStringIntoPieces( const string& s, const string& delimiter )
{
	vector<string> subTags;
	string::size_type pos = 0, prev_pos = 0;
	while( (pos = s.find_first_of(delimiter, pos ) ) != string::npos )
	{
		subTags.push_back( s.substr( prev_pos,pos-prev_pos ) );
		prev_pos = ++pos;
	}

	/// last one not to be forgotten
	subTags.push_back( s.substr( prev_pos,pos-prev_pos ) );
	return subTags;
}


double GenMatuToYearTerm(ARM_GEN_MATU gm)
{
    double yt = 0.0;
 
    switch (gm)
    {
        case GM_3M : yt = 0.25; break;
        case GM_6M : yt = 0.5;  break;
        case GM_9M : yt = 0.75; break;
        case GM_1Y : yt = 1.0;  break;
        case GM_2Y : yt = 2.0;  break;
        case GM_3Y : yt = 3.0;  break;
        case GM_4Y : yt = 4.0;  break;
        case GM_5Y : yt = 5.0;  break;
        case GM_6Y : yt = 6.0;  break;
        case GM_7Y : yt = 7.0;  break;
        case GM_8Y : yt = 8.0;  break;
        case GM_9Y : yt = 9.0;  break;
        case GM_10Y : yt = 10.0;  break;
        case GM_12Y : yt = 12.0;  break;
        case GM_15Y : yt = 15.0;  break;
        case GM_20Y : yt = 20.0;  break;
        case GM_25Y : yt = 25.0;  break;
        case GM_30Y : yt = 30.0;  break;
        case GM_40Y : yt = 40.0;  break;
        default : yt = 40.0;
    }
 
    return yt;
}



int FromLiborTypeToFrequency(int index)
{
    switch (index) 
    {
        case K_LIBOR1M :
        case K_PIBOR1M :
        case K_EURIBOR1M :
                return 12 ;
            break;

        case K_LIBOR3M :
        case K_PIBOR3M :
        case K_EURIBOR3M :
        case K_EUR3M :
                return 4 ;
            break;

        case K_LIBOR6M :
        case K_PIBOR6M :
        case K_EURIBOR6M :
                return 2 ;
            break;

        case K_LIBOR1Y :
        case K_PIBOR1Y :
        case K_EURIBOR1Y :
        case K_EUR12 :
                return 1 ;
            break;
        default : return 2;
     }
       
}

// Methode qui renvoie le type de maturite (de type ARM_INDEX_TYPE)
// a partir du term et de l'index (sous forme de char*)

int FromStringToIndexType(const char* term, const char* index)
{
    if (!strcmp(index, "FIXED"))
        return K_FIXED;

    if (!strcmp(index, "PIBOR"))
    {
        if (!strcmp(term, "1M"))
            return K_PIBOR1M;
        if (!strcmp(term, "3M"))
            return K_PIBOR3M;
        if (!strcmp(term, "6M"))
            return K_PIBOR6M;
        if (!strcmp(term, "1Y"))
            return K_PIBOR1Y;

        throw Exception(__LINE__,__FILE__,ERR_INVALID_ARGUMENT,
            "Not a good term for Pibor");
    }

    if (!strncmp(index, "EUR", 3))
    {
        if (!strcmp(term, "1M"))
            return K_EURIBOR1M;
        if (!strcmp(term, "3M"))
            return K_EURIBOR3M;
        if (!strcmp(term, "6M"))
            return K_EURIBOR6M;
        if (!strcmp(term, "1Y"))
            return K_EURIBOR1Y;
        if (!strcmp(term, "12M"))
            return K_EURIBOR1Y;

        throw Exception(__LINE__,__FILE__,ERR_INVALID_ARGUMENT,
            "Not a good term for Euribor");
    }

    // Pour tous les autres index (LIBOR, AUBB...) -> LIBOR dans ARM
    
    if (!strcmp(term, "1M"))
        return K_LIBOR1M;
    if (!strcmp(term, "3M"))
        return K_LIBOR3M;
    if (!strcmp(term, "6M"))
        return K_LIBOR6M;
    if (!strcmp(term, "1Y"))
        return K_LIBOR1Y;
    if (!strcmp(term, "12M"))
        return K_LIBOR1Y;

        throw Exception(__LINE__,__FILE__,ERR_INVALID_ARGUMENT,
            "Not a good term for Libor");
}



void GetSummitCurveIndexFromCurrency(char* ccy, char* index)
{
    static int firstCall = 1;

    char* country = NULL;


    if (firstCall)
    {
       country = getenv("ARM_COUNTRY");

       firstCall = 0;
    }

    if ( country != NULL )
    {
       if ( strcmp(country, "JPY") == 0 )
       {
          if ( strcmp(ccy, "EUR") == 0 )
          {
             strcpy(index, "LIBOR");

             return;
          }
       }
    }


    if ( strcmp(ccy, "FRF") == 0 )
    {
       strcpy(index, "PIBOR");
    }
    else if ( strcmp(ccy, "EUR") == 0 )
    {
       strcpy(index, "EURIB");
    }
    else if ( strcmp(ccy, "SEK") == 0 )
    {
       strcpy(index, "STIBO");
    }
    else if ( strcmp(ccy, "DKK") == 0 )
    {
       strcpy(index, "CIBOR");
    }
    else if ( strcmp(ccy, "AUD") == 0 )
    {
       strcpy(index, "AUBB");
    }
    else if ( strcmp(ccy, "NZD") == 0 )
    {
       strcpy(index, "BBR");
    }
    else if ( strcmp(ccy, "ZAR") == 0 )
    {
       strcpy(index, "JIBAR");
    }
    else if ( strcmp(ccy, "FIM") == 0 )
    {
       strcpy(index, "HELIB");
    }
    else if ( strcmp(ccy, "PTE") == 0 )
    {
       strcpy(index, "LISBO");
    }
    else if ( strcmp(ccy, "ATS") == 0 )
    {
       strcpy(index, "VIBOR");
    }
    else if ( strcmp(ccy, "NOK") == 0 )
    {
       strcpy(index, "OIBOR");
    }
    else if ( strcmp(ccy, "IEP") == 0 )
    {
       strcpy(index, "DIBOR");
    }
    else if ( strcmp(ccy, "NLG") == 0 )
    {
       strcpy(index, "AIBOR");
    }
    else if ( strcmp(ccy, "BEF") == 0 )
    {
       strcpy(index, "BIBOR");
    }
/**    else if ( strcmp(ccy, "THB") == 0 )
    {
       strcpy(index, "SOR");
    }
	*/ 
    else if ( strcmp(ccy, "SGD") == 0 )
    {
       strcpy(index, "SOR");
    }
    else
    {
       strcpy(index, "LIBOR");
    }
}

//	-----------------------------------------------------------------------------------------
std::string GetSummitCurveIndexFromCurrency(const std::string& ccy)
{
	char index[16] ; // just like mercure
	GetSummitCurveIndexFromCurrency((char*)ccy.c_str(),index); 
	return index; 
}

int FromSummitDaycountToARMDaycount(const char* daycount)
{
    if ( strcmp(daycount,"30/360") == 0)
        return K30_360;
    else if ( (strcmp(daycount,"ACT") == 0) ||
		(strcmp(daycount,"ACT/ACT")==0) )
        return KACTUAL_ACTUAL;
    else if ( (strcmp(daycount,"A365F") == 0) ||
			(strcmp(daycount,"ACT/365") == 0) )
        return KACTUAL_365;
    else if ( (strcmp(daycount,"A360") == 0)||
			(strcmp(daycount,"ACT/360") == 0) )
        return KACTUAL_360;
    else if (strcmp(daycount,"A30") == 0)
        return K30_360;
    else if ( (strcmp(daycount,"ACT29") == 0)||
			(strcmp(daycount,"ACT/ACT29") == 0)) 
        return KACTUAL_FEB29;
    else if (strcmp(daycount,"COUP") == 0)
        return KCOUPON;
	else if ( strcmp(daycount,"30E") == 0 )
        return K30_360E;

    throw Exception(__LINE__,__FILE__,ERR_INVALID_ARGUMENT,
            "Not a good Summit Daycount");

}

int FromSummitFreqToARMFreq(const char* freq)
{
    if ( (strcmp(freq,"A") == 0) || (strcmp(freq,"1") == 0) 
		|| (strcmp(freq,"PA") == 0) )
        return K_ANNUAL;
    else if ( (strcmp(freq,"S") == 0) || (strcmp(freq,"2") == 0) 
		|| (strcmp(freq,"SA") == 0) )
        return K_SEMIANNUAL;
    else if ( (strcmp(freq,"Q") == 0) || (strcmp(freq,"3") == 0) 
			|| (strcmp(freq,"QTR") == 0) )
        return K_QUARTERLY;
    else if ( (strcmp(freq,"M") == 0) || (strcmp(freq,"4") == 0) 
		|| (strcmp(freq,"MTH") == 0) )
        return K_MONTHLY;
    else if ( (strcmp(freq,"B") == 0)
		|| (strcmp(freq,"BIM") == 0) )
        return K_BIMONTHLY;
    else if ( (strcmp(freq,"W") == 0) || (strcmp(freq,"5") == 0) 
		|| (strcmp(freq,"WK") == 0) )
        return K_WEEKLY;
    else if ( (strcmp(freq,"D") == 0) || (strcmp(freq,"6") == 0) 
		|| (strcmp(freq,"1D") == 0) )
        return K_DAILY;
    else if ( (strcmp(freq,"Z") == 0) || (strcmp(freq,"7") == 0) 
		|| (strcmp(freq,"ZC") == 0) )
        return K_ZEROCOUPON;

    throw Exception(__LINE__,__FILE__,ERR_INVALID_ARGUMENT,
            "Not a good Summit Frequency");

}


int FromSummitGapToARMGap(const char* dayGap)
{	
	std::string str(dayGap);
	
	if (str[str.size() - 1] != 'D')
		if (str[str.size() - 1] == 'Y' || str[str.size() - 1] == 'C')
			throw Exception(__LINE__ , __FILE__ , ERR_INVALID_ARGUMENT , "Not a summit DayGap taken in account by ARM");

    return atoi(str.substr(0, str.size()-1).c_str());
}


double ConvertCouponAsSummit(double coupon, int Freq, int DayCount, int newFreq, int newDayCount)
{
    double period = 1.0/Freq;
    switch (DayCount)
    {
        case K30_360 : 
        case K30_360E :
        {
            period *= 1.0;
        };
        break;

        case KACTUAL_360 : 
        {
            period *= 365.25/360.0;
        };
        break;

        case KACTUAL_365 : 
        {
            period *= 365.25/365.0;
        };
        break;

        case KACTUAL_ACTUAL : 
        {
            period *= 1.0;
        };
        break;

        default :
        {
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "FixedDayCount not treated In CalibrationPF ");
        };
        break;
    }        
    
    double newperiod = (double) (1.0/newFreq);
    switch (newDayCount)
    {
        case K30_360 : 
        case K30_360E :
        {
            newperiod *= 1.0;
        };
        break;

        case KACTUAL_360 : 
        {
            newperiod *= 365.25/360.0;
        };
        break;

        case KACTUAL_365 : 
        {
            newperiod *= 365.25/365.0;
        };
        break;

        case KACTUAL_ACTUAL : 
        {
            newperiod *= 1.0;
        };
        break;

        default :
        {
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "FixedDayCount not treated In CalibrationPF ");
        };
        break;
    }

    double puissance = ((double) (Freq)/ (double) (newFreq));
    double newCoupon = 100.0*(pow((1.0+period*coupon/100.0),puissance)-1.0)/newperiod;
    return newCoupon;
}



double CalcFwdFXSpot(ARM_Date& AsofDate,
                     double Spot,
                     ARM_Date& aFwdDate,
                     ARM_ZeroCurve* NumDiscountCurve,
                     ARM_ZeroCurve* UndDiscountCurve)
{
    double SpotMat, FwdMat, Denom;
    double Fwd = Spot;
    double UndDf, NumDf;
   
    ARM_Date FwdDate = aFwdDate;
    ARM_Date SpotDate = AsofDate;
   
    /*int SpotDays = MAX(NumDiscountCurve->GetCurrencyUnit()->GetSpotDays(),
                       UndDiscountCurve->GetCurrencyUnit()->GetSpotDays());

    char PayCal[20];
    strcpy(PayCal, NumDiscountCurve->GetCurrencyUnit()->GetCcyName());
    strcat(PayCal, UndDiscountCurve->GetCurrencyUnit()->GetCcyName());*/

//    SpotDate.NextBusinessDay(SpotDays, PayCal);
    //FwdDate.NextBusinessDay(SpotDays, PayCal);

    //SpotMat = CountYears((int) KACTUAL_365, AsofDate, SpotDate);
    //FwdMat  = CountYears((int) KACTUAL_365, AsofDate, FwdDate);

    //UndDf = UndDiscountCurve->DiscountPrice(SpotMat);
    //Denom = UndDiscountCurve->DiscountPrice(FwdMat);
   
    //// filtering the zero-divide:
    //if ( fabs(Denom) < 1.0e-16 )
    //{
    //   throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
    //                   "CalcFwdFXSpot: zero-divide");

    //   return(0.0);
    //}

    //UndDf /= Denom;

    //NumDf = NumDiscountCurve->DiscountPrice(SpotMat);
    //Denom = NumDiscountCurve->DiscountPrice(FwdMat);

    //// filtering the zero-divide:
    //if ( fabs(Denom) < 1.0e-16 )
    //{
    //   throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
    //      "CalcFwdFXSpot: zero-divide");

    //   return(0.0);
    //}

    //NumDf /= Denom;

    //Fwd *= NumDf;
    //Fwd /= UndDf;

    //return(Fwd);
	return 0;
}


void GetFrqAndNbFromString(const char* period, int& Nb, long& freqId)
{
    char matu;

    sscanf(period, "%d%c", &Nb, &matu);

    switch (matu)
    {
    case 'D':
        freqId = K_DAILY;
        break;

    case 'W':
        freqId = K_WEEKLY;
        break;

    case 'M':
        freqId = K_MONTHLY;
        break;

    case 'Y':
        freqId = K_ANNUAL;
        break;

    default:
        throw Exception(__LINE__,__FILE__,ERR_INVALID_ARGUMENT,
                   "Not a good gap for adding period");
    }
}


long ARM_ConvPayResetRule (const char* aRule)
{
    
    char temp[20];
    strcpy(temp,aRule);
    char* tmp = ARM_UPPERCASE(temp);

    if(strcmp(tmp, "ADV")==0)
    {
        return K_ADVANCE;
    }
    else if(strcmp(tmp, "ADVANCE")==0)
    {
        return K_ADVANCE;
    }
    else if(strcmp(tmp, "ARR")==0)
    {
        return K_ARREARS;
    }
    else if(strcmp(tmp, "ARREARS")==0)
    {
        return K_ARREARS;
    }
    else if(strcmp(tmp, "1")==0)
    {
        return K_ADVANCE;
    }
    else if(strcmp(tmp, "-1")==0)
    {
        return K_ARREARS;
    }
    throw Exception(__LINE__, __FILE__, ERR_INVALID_INPUT,
                        "ARM_ERR: Invalid Payment or Reset Rule - Valid are  ADV, ADVANCE, 1 or ARR, ARREARS, -1");
    return ARM_DEFAULT_ERR;
}

long ARM_ConvCvMethod (const char* aCvMeth)
{
	char temp[20];
	strcpy(temp,aCvMeth);
	char* tmp = ARM_UPPERCASE(temp);

	if ( tmp[0] == 'P' )
	{
	   return(K_PAR);
	}

    if ( tmp[0] == 'R' )
	{
       return(K_RAW);
	}

	if ( tmp[0] == 'F' )
	{
       return(K_FORWARD);
	}

    if ( tmp[0] == '0' )
	{
       return(K_PAR);
	}

    if ( tmp[0] == '1' )
	{
       return(K_RAW);
	}

	if ( tmp[0] == '2' )
	{
       return(K_FORWARD);
	}

	// pour le CreateZcFromSummit
    if ( (strcmp(aCvMeth,"DEFAULT") == 0) || (strcmp(aCvMeth,"") == 0) )
	{
        return K_DEFAULT_CURVEMOD;
	}

	throw Exception(__LINE__, __FILE__, ERR_INVALID_INPUT,
                        "ARM_ERR: Invalid Cv Method - Valid are: (PAR,0) - (RAW,1) - (FORWARD,2)");
	return ARM_DEFAULT_ERR;
}

long ARM_ConvIrType (const char*  aIrType)
{
    long irType = ARM_ConvIrTypeWithoutException(aIrType);
    
    if (irType != ARM_DEFAULT_ERR)
        return irType;
    else
        throw Exception(__LINE__, __FILE__, ERR_INVALID_INPUT,
                            "ARM_ERR: Invalid Type index - Valid are Libor3m Pibor6m CMS10 etc ...");
}


long ARM_ConvIrTypeWithoutException (const char*  aIrType)
{
    char temp[20];
    strcpy(temp,aIrType);
    char* tmp = ARM_UPPERCASE(temp);

    
    if ( strcmp(tmp, "FIXED") == 0 )
    {
        return K_FIXED;
    }

    if ( strcmp(tmp, "LIBOR1M") == 0 )
    {
        return K_LIBOR1M;
    }

    if ( strcmp(tmp, "LIBOR2M") == 0 )
    {
        return K_LIBOR2M;
    }

    if ( strcmp(tmp, "LIBOR3M") == 0 )
    {
        return K_LIBOR3M;
    }

    if ( strcmp(tmp, "LIBOR6M") == 0 )
    {
        return K_LIBOR6M;
    }

    if ( strcmp(tmp, "LIBOR1Y") == 0 )
    {
        return K_LIBOR1Y;
    }

    if ( strcmp(tmp, "LIBOR12M") == 0 )
    {
        return K_LIBOR1Y;
    }

    if ( strcmp(tmp, "PIBOR1M") == 0 )
    {
        return K_PIBOR1M;
    }
    
    if ( strcmp(tmp, "PIBOR2M") == 0 )
    {
        return K_PIBOR2M;
    }

    if ( strcmp(tmp, "PIBOR3M") == 0 )
    {
        return K_PIBOR3M;
    }

    if ( strcmp(tmp, "PIBOR6M") == 0 )
    {
        return K_PIBOR6M;
    }

    if ( strcmp(tmp, "PIBOR12M") == 0 )
    {
        return K_PIBOR1Y;
    }

    if ( strcmp(tmp, "PIBOR1Y") == 0 )
    {
        return K_PIBOR1Y;
    }

    if ( strcmp(tmp, "EURIBOR1M") == 0 )
    {
        return K_EURIBOR1M;
    }

    if ( strcmp(tmp, "EURIBOR2M") == 0 )
    {
        return K_EURIBOR2M;
    }

    if ( strcmp(tmp, "EURIBOR3M") == 0 )
    {
        return K_EURIBOR3M;
    }

    if ( strcmp(tmp, "EURIBOR6M") == 0 )
    {
        return K_EURIBOR6M;
    }

    if ( strcmp(tmp, "EURIBOR1Y") == 0 )
    {
        return K_EURIBOR1Y;
    }

    if ( strcmp(tmp, "EURIBOR12M") == 0 )
    {
        return K_EURIBOR1Y;
    }

    if ( strcmp(tmp, "CMT1") == 0 )
    {
        return K_CMT1;
    }

    if ( strcmp(tmp, "CMT2") == 0 )
    {
        return K_CMT2;
    }

    if ( strcmp(tmp, "CMT5") == 0 )
    {
        return K_CMT5;
    }

    if ( strcmp(tmp, "CMT10") == 0 )
    {
        return K_CMT10;
    }

    if ( strcmp(tmp, "CMT15") == 0 )
    {
        return K_CMT15;
    }

    if ( strcmp(tmp, "CMT20") == 0 )
    {
        return K_CMT20;
    }

    if ( strcmp(tmp, "CMT30") == 0 )
    {
        return K_CMT30;
    }

    if ( strcmp(tmp, "CMS1") == 0 )
    {
        return K_CMS1;
    }

    if ( strcmp(tmp, "CMS2") == 0 )
    {
        return K_CMS2;
    }

    if ( strcmp(tmp, "CMS3") == 0 )
    {
        return K_CMS3;
    }

    if ( strcmp(tmp, "CMS4") == 0 )
    {
        return K_CMS4;
    }

    if ( strcmp(tmp, "CMS5") == 0 )
    {
        return K_CMS5;
    }

    if ( strcmp(tmp, "CMS6") == 0 )
    {
        return K_CMS6;
    }   

    if ( strcmp(tmp, "CMS7") == 0 )
    {
        return K_CMS7;
    }   

    if ( strcmp(tmp, "CMS8") == 0 )
    {
        return K_CMS8;
    }   

    if ( strcmp(tmp, "CMS9") == 0 )
    {
        return K_CMS9;
    }   

 
    if ( strcmp(tmp, "CMS10") == 0 )
    {
        return K_CMS10;
    }

    if ( strcmp(tmp, "CMS11") == 0 )
    {
        return K_CMS11;
    }

    if ( strcmp(tmp, "CMS12") == 0 )
    {
        return K_CMS12;
    }

    if ( strcmp(tmp, "CMS13") == 0 )
    {
        return K_CMS13;
    }

    if ( strcmp(tmp, "CMS14") == 0 )
    {
        return K_CMS14;
    }

    if ( strcmp(tmp, "CMS15") == 0 )
    {
        return K_CMS15;
    }

    if ( strcmp(tmp, "CMS16") == 0 )
    {
        return K_CMS16;
    }

    if ( strcmp(tmp, "CMS17") == 0 )
    {
        return K_CMS17;
    }

    if ( strcmp(tmp, "CMS18") == 0 )
    {
        return K_CMS18;
    }

    if ( strcmp(tmp, "CMS19") == 0 )
    {
        return K_CMS19;
    }

    if ( strcmp(tmp, "CMS20") == 0 )
    {
        return K_CMS20;
    }

    if ( strcmp(tmp, "CMS21") == 0 )
    {
        return K_CMS21;
    }

    if ( strcmp(tmp, "CMS22") == 0 )
    {
        return K_CMS22;
    }

    if ( strcmp(tmp, "CMS23") == 0 )
    {
        return K_CMS23;
    }

    if ( strcmp(tmp, "CMS24") == 0 )
    {
        return K_CMS24;
    }

    if ( strcmp(tmp, "CMS25") == 0 )
    {
        return K_CMS25;
    }

    if ( strcmp(tmp, "CMS26") == 0 )
    {
        return K_CMS26;
    }

    if ( strcmp(tmp, "CMS27") == 0 )
    {
        return K_CMS27;
    }

    if ( strcmp(tmp, "CMS28") == 0 )
    {
        return K_CMS28;
    }

    if ( strcmp(tmp, "CMS29") == 0 )
    {
        return K_CMS29;
    }

    if ( strcmp(tmp, "CMS30") == 0 )
    {
        return K_CMS30;
    }

    if ( strcmp(tmp, "TEC5") == 0 )
    {
        return K_TEC5;
    }

    if ( strcmp(tmp, "TEC10") == 0 )
    {
        return K_TEC10;
    }

    if ( strcmp(tmp, "T4M") == 0 )
    {
        return K_T4M;
    }

    if ( strcmp(tmp, "T4M_FIXED") == 0 )
    {
        return K_T4M_FIXED;
    }

    if ( strcmp(tmp, "TAM") == 0 )
    {
        return K_TAM;
    }

    if ( strcmp(tmp, "TAG") == 0 )
    {
        return K_TAG;
    }

    if ( strcmp(tmp, "TMP") == 0 )
    {
        return K_TMP;
    }

    if ( strcmp(tmp, "EONIA") == 0 )
    {
        return K_EONIA;   // In fact K_EONIA == K_TMP
    }

    if ( strcmp(tmp, "STD") == 0 )
    {
        return K_STD;
    }

    if ( strcmp(tmp, "EUR3M") == 0 )
    {
        return K_EUR3M;
    }

	if ( strcmp(tmp, "EUR1M") == 0 )
    {
        return K_EURIBOR1M;
    }

    if(strcmp(tmp, "EUR12") ==0)
    {
        return K_EUR12;
    }

    return(ARM_DEFAULT_ERR);
}


long ARM_ConvDayCount (const char* aDayCount)
{
    char temp[20];
    strcpy(temp,aDayCount);
    char* tmp = ARM_UPPERCASE(temp);

    if(strcmp(tmp, "ACTUAL") ==0 || strcmp(tmp, "1") ==0)
    {
        return KACTUAL_ACTUAL;
    }
    if(strcmp(tmp, "A365") ==0 || strcmp(tmp, "2") ==0)
    {
        return KACTUAL_365;
    }
    if(strcmp(tmp, "A360") ==0 || strcmp(tmp, "3") ==0)
    {
        return KACTUAL_360;
    }
    if(strcmp(tmp, "30/360") ==0 || strcmp(tmp, "4") ==0)
    {
        return K30_360;
    }
    if(strcmp(tmp, "ACTREAL") ==0 || strcmp(tmp, "5") ==0)
    {
        return KACTUAL_REAL;
    }
    if(strcmp(tmp, "ACT29") ==0 || strcmp(tmp, "6") ==0)
    {
        return KACTUAL_FEB29;
    }

    if(strcmp(tmp, "ISMA") ==0 || strcmp(tmp, "7") ==0)
    {
        return KACTUAL_ISMA;
    }

    if(strcmp(tmp, "NOBASE") ==0 || strcmp(tmp, "8") ==0)
    {
        return KNOBASE;
    }
	if(strcmp(tmp, "30/360E") ==0 || strcmp(tmp, "9") ==0)
    {
        return K30_360E;
    }

    if(strcmp(tmp, "-1") ==0)
    {
        return -1;
    }

    throw Exception(__LINE__, __FILE__, ERR_INVALID_INPUT,
                        "ARM_ERR: Invalid DayCount - Valid are ACTUAL: A365, A360, 30/360, ACTREAL, ISMA, ACT29, NOBASE, 30/360E");
    return ARM_DEFAULT_ERR;

}


long ARM_ConvDecompFrequency (const char* aFrequency)
{
    char temp[20];
    strcpy(temp,aFrequency);
    char* tmp = ARM_UPPERCASE(temp);

    if(strcmp(tmp, "C") ==0)
    {
        return K_COMP_CONT;
    }
    if(strcmp(tmp, "P") ==0)
    {
        return K_COMP_PROP;
    }
    if(strcmp(tmp, "A") ==0)
    {
        return K_COMP_ANNUAL;
    }
    if(strcmp(tmp, "S") ==0)
    {
        return K_COMP_SEMIANNUAL;
    }
    if(strcmp(tmp, "Q") ==0)
    {
        return K_COMP_QUARTERLY;
    }
    if(strcmp(tmp, "M") ==0)
    {
        return K_COMP_MONTHLY;
    }
    if(strcmp(tmp, "B") ==0)
    {
        return K_COMP_BIMONTHLY;
    }
    if(strcmp(tmp, "D360") ==0)
    {
        return K_COMP_DAILY_360;
    }
    if(strcmp(tmp, "D365") ==0)
    {
        return K_COMP_DAILY_365;
    }
    if(strcmp(tmp, "0") ==0)
    {
        return K_COMP_CONT;
    }
    if(strcmp(tmp, "-1") ==0)
    {
        return K_COMP_PROP;
    }
    if(strcmp(tmp, "1") ==0)
    {
        return K_COMP_ANNUAL;
    }
    if(strcmp(tmp, "2") ==0)
    {
        return K_COMP_SEMIANNUAL;
    }
    if(strcmp(tmp, "4") ==0)
    {
        return K_COMP_QUARTERLY;
    }
    if(strcmp(tmp, "12") ==0)
    {
        return K_COMP_MONTHLY;
    }
    if(strcmp(tmp, "6") ==0)
    {
        return K_COMP_BIMONTHLY;
    }
    if(strcmp(tmp, "360") ==0)
    {
        return K_COMP_DAILY_360;
    }
    if(strcmp(tmp, "365") ==0)
    {
        return K_COMP_DAILY_365;
    }

    throw Exception(__LINE__, __FILE__, ERR_INVALID_INPUT,
                        "ARM_ERR: Invalid Type Name - Valid are A, S, Q, B(imonth), M, D, Z");
    return ARM_DEFAULT_ERR;
    
}


long ARM_ConvIntRule (const char* aRule)
{
    char temp[20];
    strcpy(temp,aRule);
    char* tmp = ARM_UPPERCASE(temp);

    
    if(strcmp(tmp, "ADJ") ==0)
    {
        return K_ADJUSTED;
    }
    if(strcmp(tmp, "A") ==0)
    {
        return K_ADJUSTED;
    }
    if(strcmp(tmp, "NOADJ") ==0)
    {
        return K_UNADJUSTED;
    }
    if(strcmp(tmp, "UNADJ") ==0)
    {
        return K_UNADJUSTED;
    }
    if(strcmp(tmp, "MATUNADJUSTED") ==0)
    {
        return K_MATUNADJUSTED;
    }
    if(strcmp(tmp, "MATUNADJ") ==0)
    {
        return K_MATUNADJUSTED;
    }
    if(strcmp(tmp, "M") ==0)
    {
        return K_MATUNADJUSTED;
    }
    if(strcmp(tmp, "N") ==0)
    {
        return K_UNADJUSTED;
    }
    if(strcmp(tmp, "U") ==0)
    {
        return K_UNADJUSTED;
    }
    if(strcmp(tmp, "1") ==0)
    {
        return K_ADJUSTED;
    }
    if(strcmp(tmp, "0") ==0)
    {
        return K_UNADJUSTED;
    }
    if(strcmp(tmp, "Y") ==0)
    {
        return K_ADJUSTED;
    }
    if(strcmp(tmp, "N") ==0)
    {
        return K_UNADJUSTED;
    }
    
    throw Exception(__LINE__, __FILE__, ERR_INVALID_INPUT,
                        "ARM_ERR: Invalid Type Name - Valid are A(ADJ,1,Y),N(UNADJ,0),M(MATUNADJ)");
    return ARM_DEFAULT_ERR;

}

long ARM_ConvStartAdjRule(const string& aRule)
{
	char temp[20];
    strcpy(temp,aRule.c_str());
    char* tmp = ARM_UPPERCASE(temp);

    
    if(strcmp(tmp, "ADJ") ==0)
    {
        return K_ADJUSTED;
    }
    if(strcmp(tmp, "A") ==0)
    {
        return K_ADJUSTED;
    }
    if(strcmp(tmp, "NOADJ") ==0)
    {
        return K_UNADJUSTED;
    }
    if(strcmp(tmp, "UNADJ") ==0)
    {
        return K_UNADJUSTED;
    }
    if(strcmp(tmp, "N") ==0)
    {
        return K_UNADJUSTED;
    }
    if(strcmp(tmp, "U") ==0)
    {
        return K_UNADJUSTED;
    }
    if(strcmp(tmp, "1") ==0)
    {
        return K_ADJUSTED;
    }
    if(strcmp(tmp, "0") ==0)
    {
        return K_UNADJUSTED;
    }
    if(strcmp(tmp, "Y") ==0)
    {
        return K_ADJUSTED;
    }
    if(strcmp(tmp, "N") ==0)
    {
        return K_UNADJUSTED;
    }
    
    throw Exception(__LINE__, __FILE__, ERR_INVALID_INPUT,
                        "ARM_ERR: Invalid Type Name - Valid are A(ADJ,1,Y),N(UNADJ,0)");
    return ARM_DEFAULT_ERR;

}

long ARM_ConvStubRule (const char* aRule)
{
    char temp[20];
    strcpy(temp,aRule);
    char* tmp = ARM_UPPERCASE(temp);


    if( (strcmp(tmp, "SS") ==0) || (strcmp(tmp, "") ==0) )
    {
        return K_SHORTSTART;
    }
    if(strcmp(tmp, "LS") ==0)
    {
        return K_LONGSTART;
    }
    if(strcmp(tmp, "SE") ==0)
    {
        return K_SHORTEND;
    }
    if(strcmp(tmp, "LE") ==0)
    {
        return K_LONGEND;
    }
    if(strcmp(tmp, "1") ==0)
    {
        return K_SHORTSTART;
    }
    if(strcmp(tmp, "2") ==0)
    {
        return K_LONGSTART;
    }
    if(strcmp(tmp, "3") ==0)
    {
        return K_SHORTEND;
    }
    if(strcmp(tmp, "4") ==0)
    {
        return K_LONGEND;
    }

    throw Exception(__LINE__, __FILE__, ERR_INVALID_INPUT,
                        "ARM_ERR: Invalid Stub Rule - Valid are  SS , LS, SE, LE or 1, 2, 3, 4");
    return ARM_DEFAULT_ERR;
}



long ARM_NotionalExchange(const char* aRule)
{
    char temp[20];
    strcpy(temp,aRule);
    char* tmp = ARM_UPPERCASE(temp);

    if ( strcmp(tmp, "NXNONE") ==0)
    {
       return(K_NX_NONE);
    }
    if ( strcmp(tmp, "NXSTART") ==0)
    {
        return(K_NX_START);
    }
    if ( strcmp(tmp, "NXEND") ==0)
    {
        return(K_NX_END);
    }
    if ( strcmp(tmp, "NXBOTH") ==0)
    {
        return(K_NX_BOTH);
    }
    if ( strcmp(tmp, "START") ==0)
    {
        return(K_NX_START);
    }
    if ( strcmp(tmp, "END") ==0)
    {
        return(K_NX_END);
    }
    if ( strcmp(tmp, "BOTH") ==0)
    {
        return(K_NX_BOTH);
    }
    if ( strcmp(tmp, "NEITHER") ==0)
    {
        return(K_NX_NONE);
    }

    throw Exception(__LINE__, __FILE__, ERR_INVALID_INPUT,
                        "ARM_ERR: Invalid Notional Exchange Rule - Valid are  NXNONE, NXSTART, NXEND or NXBOTH");
    return ARM_DEFAULT_ERR;
}



long ARM_ConvPricCorridorLD(const  char* param)
{
    char temp[20];

    strcpy(temp,param);
    char* tmp = ARM_UPPERCASE(temp);

    if ((strcmp(tmp, "DIG") ==0) || (strcmp(tmp, "DIGIT") ==0))
        return K_DIGITALE;

    if ((strcmp(tmp, "SPR") ==0) || (strcmp(tmp, "SPREAD") ==0))
        return K_SPREAD;

    throw Exception(__LINE__, __FILE__, ERR_INVALID_INPUT,
                        "ARM_ERR: Invalid Parameter - Valid are DIG or SPR");
    return ARM_DEFAULT_ERR;
}



ARM_Vector* CreateARMVectorFromVECTOR(const std::vector<double>& param, int newSize)
{
    long size = param.size();

    if (size == 0)
        return NULL;

    int mySize;

    if (( newSize < 0 ) || ( newSize > size ))
       mySize = size;
    else
       mySize = newSize;

    ARM_Vector* res = new ARM_Vector(mySize);

    for (int i = 0; i < mySize; i++)
        res->Elt(i) = param[i];

    return res;
}


int IsFixedIndex(ARM_INDEX_TYPE idx)
{
    switch(idx)
    {
        case IDXFIXED :
        {
            return(1);
        };
        break;
 
        default :
        {
            return(0);
        };
        break;
    }
}



int IsLiborIndex(ARM_INDEX_TYPE idx)
{
    switch(idx)
    {
        case PIBOR1M :
        case PIBOR3M :
        case PIBOR6M :
        case PIBOR1Y :
        case LIBOR1M :
        case LIBOR3M :
        case LIBOR6M :
        case LIBOR1Y :
        case EURIBOR1M :
        case EURIBOR3M :
        case EURIBOR6M :
        case EURIBOR1Y :
		case EUR3M:
		case EUR12:
        {
            return(1);
        };
        break;

        default :
        {
            return(0);
        };
        break;
    }
}



int IsCMSIndex(ARM_INDEX_TYPE idx)
{
    switch(idx)
    {
        case CMS1 :
        case CMS2 :
        case CMS3 :
        case CMS4 :
        case CMS5 :
        case CMS6 :
        case CMS7 :
        case CMS8 :
        case CMS9 :
        case CMS10:
        case CMS11:
        case CMS12:
        case CMS13:
        case CMS14:
        case CMS15:
        case CMS16:
        case CMS17:
        case CMS18:
        case CMS19:
        case CMS20:
        case CMS21:
        case CMS22:
        case CMS23:
        case CMS24:
        case CMS25:
        case CMS26:
        case CMS27:
        case CMS28:
        case CMS29:
        case CMS30:
        {
            return(1);
        };
        break;

        default :
        {
            return(0);
        };
        break;
    }
}



int IsCMTIndex(ARM_INDEX_TYPE idx)
{
    switch(idx)
    {
        case CMT1 :
        case CMT2 :
        case CMT5 :
        case CMT10:
        case CMT15:
        case CMT20:
        case CMT30:
        {
            return(1);
        };
        break;

        default :
        {
            return(0);
        };
        break;
    }
}



long ARM_ConvWhatIsInterp (const char* aIntMeth)
{
    char tmp[40];

    strcpy(tmp, aIntMeth);

    ARM_UPPERCASE(tmp);

    if( strcmp(tmp, "DEFAULT") == 0 )
    {
        return FXINTERP_STRIKE;
    }

    if( tmp[0] == 'S' )
    {
        return FXINTERP_STRIKE;
    }
    if( tmp[0] == 'D' )
    {
        return FXINTERP_DELTA;
    }
    if( tmp[0] == 'E' )
    {
        return FXINTERP_DELTA_SMOOTH;
    }
    if( tmp[0] == '0' )
    {
        return FXINTERP_STRIKE;
    }
    if( tmp[0] == '1' )
    {
        return FXINTERP_DELTA;
    }
    if( tmp[0] == '2' )
    {
        return FXINTERP_DELTA_SMOOTH;
    }
    
    throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA, "Invalid what is interpolated - Valid are: (S,0) - (D,1) - (E,2)");

    return ARM_DEFAULT_ERR;
}


long ARM_ConvMonth (const char* aMonth)
{
    char temp[20];
    strcpy(temp,aMonth);
    char* tmp = ARM_UPPERCASE(temp);

    if(strcmp(tmp, "JAN") ==0)
    {
        return K_JANUARY;
    }
    if(strcmp(tmp, "FEB") ==0)
    {
        return K_FEBRUARY;
    }
    if(strcmp(tmp, "MAR") ==0)
    {
        return K_MARCH;
    }
    if(strcmp(tmp, "APR") ==0)
    {
        return K_APRIL;
    }
    if(strcmp(tmp, "MAY") ==0)
    {
        return K_MAY;
    }
    if(strcmp(tmp, "JUN") ==0)
    {
        return K_JUNE;
    }
    if(strcmp(tmp, "JUL") ==0)
    {
        return K_JULY;
    }
    if(strcmp(tmp, "AUG") ==0)
    {
        return K_AUGUST;
    }
    if(strcmp(tmp, "SEP") ==0)
    {
        return K_SEPTEMBER;
    }
    if(strcmp(tmp, "OCT") ==0)
    {
        return K_OCTOBER;
    }
    if(strcmp(tmp, "NOV") ==0)
    {
        return K_NOVEMBER;
    }
    if(strcmp(tmp, "DEC") ==0)
    {
        return K_DECEMBER;
    }
    if(strcmp(tmp, "DEFAULT") ==0)
    {
        return K_MONTH_DEFAULT;
    }

    throw Exception(__LINE__, __FILE__, ERR_INVALID_INPUT,
                        "ARM_ERR: Invalid MONTH Name - Valid are MAY, JUN, NOV, DEC, DEFAULT");
    return ARM_DEFAULT_ERR;
    
} 


void FromRRSTRtoVol(ARM_Matrix* RR, ARM_Matrix* STR, ARM_Vector* pivotVol)
{
    int nbCols = RR->GetNumCols();

    ARM_Matrix* newRR  = new ARM_Matrix(RR->GetNumLines(), nbCols);
    ARM_Matrix* newSTR = new ARM_Matrix(RR->GetNumLines(), nbCols);

    for (int i = 0; i < RR->GetNumLines(); i++)
    {
        for (int j = 0; j < nbCols; j++)
        {
            newRR->Elt(i, j)  = pivotVol->Elt(i)+STR->Elt(i, nbCols-j-1)-RR->Elt(i, j)/2.0;
        
            newSTR->Elt(i, j) = pivotVol->Elt(i)+STR->Elt(i, j)+RR->Elt(i, nbCols-j-1)/2.0;
        }
    }

    *RR  = *newRR;
    *STR = *newSTR;

    delete newRR;
    delete newSTR;
}

long ARM_ConvFwdRule (const char* aFwdRule)
{
    char tmp[20];

    strcpy(tmp, aFwdRule);

    (void) ARM_UPPERCASE(tmp);

    if( strcmp(tmp, "F") == 0 )
    {
        return K_FOLLOWING;
    }

    if( strcmp(tmp,"MF") == 0 )
    {
        return K_MOD_FOLLOWING;
    }

    if( strcmp(tmp, "P") == 0 )
    {
        return K_PREVIOUS;
    }

    if( strcmp(tmp, "MP") == 0 )
    {
        return K_MOD_PREVIOUS;
    }

    if( strcmp(tmp, "1") == 0 )
    {
        return K_FOLLOWING;
    }

    if( strcmp(tmp, "2") == 0 )
    {
        return K_MOD_FOLLOWING;
    }

    if(strcmp(tmp , "-1")==0 )
    {
        return K_PREVIOUS;
    }

    if(strcmp(tmp , "-2") == 0)
    {
        return K_MOD_PREVIOUS;
    }

    return ARM_DEFAULT_ERR;
}



char* ARM_STRDUP(char* inStr)
{
    if ( inStr == NULL )
    {
       return(NULL);
    }

    char* resStr = new char [strlen(inStr)+1];

    strcpy(resStr, inStr);

    return(resStr);
}

ARM_ReferenceValue* ARM_FromStartRefToResetRef(ARM_ReferenceValue* valstartRef, 
											   ARM_SwapLeg* legRef)
{
	//ARM_Vector* refdates = valstartRef->GetDiscreteDates(); 
	//ARM_Vector* refvalues = valstartRef->GetDiscreteValues();

	//ARM_Vector *schedStartDates   = legRef->GetFlowStartDates();
	//ARM_Vector *schedResetDates   = legRef->GetResetDates();

	//if ( !schedStartDates || !schedResetDates )
	//{
	//	throw Exception (__LINE__, __FILE__, ERR_OBJECT_NULL, "ARM_FromStartRefToResetRef: Error generating schedule");
	//}		

	//ARM_ReferenceValue* valresetRef = NULL;
	//ARM_Vector* resetDates = new ARM_Vector(refdates->GetSize());

	//if (refdates)
	//{
	//	// on construit le vecteur avec les dates de reset
	//	int idx;

	//	for (int i = 0; i < refdates->GetSize(); i++)
	//	{
	//		if (refdates->Elt(i) > schedStartDates->Elt(schedStartDates->GetSize()-1))
	//		{
	//			resetDates->Elt(i) = refdates->Elt(i);
	//		}
	//		else
	//		{
	//			idx = schedStartDates->find(refdates->Elt(i));
	//			// on ne la trouve pas : problème
	//			if (idx != -1)
	//			{
	//				resetDates->Elt(i) = schedResetDates->Elt(idx);
	//			}
	//			else
	//			{
	//				// eventuellement, on ajuste la date et on reteste
	//				int rule = legRef->GetIRIndex()->GetFwdRule();
	//				char* calendar = legRef->GetPayCalName();					
	//				double adjustedRef = ARM_Date(refdates->Elt(i)).GoodBusinessDay(rule, calendar).GetJulian();

	//				idx = schedStartDates->find(adjustedRef);
	//				if (idx != -1)
	//				{
	//					resetDates->Elt(i) = schedResetDates->Elt(idx);
	//				}
	//				else
	//				{
	//					char msg[200];
	//					char strDate[11];
	//					ARM_Date(refdates->Elt(i)).JulianToStrDate(strDate);
	//					sprintf(msg, "ARM_FromStartRefToResetRef: Expiry[%d] (%s) does not match any start date", i, strDate);
	//					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg);
	//				}
	//			}
	//		}
	//	}

	//	valresetRef = new ARM_ReferenceValue((ARM_Vector *)resetDates->Clone(),(ARM_Vector *)refvalues->Clone());

	//	valresetRef->SetCalcMethod(K_STEPUP_LEFT);

	//	if (resetDates)
	//		delete resetDates;
	//	resetDates = NULL;
	//}

	return  NULL;//valresetRef;
}



ARM_ReferenceValue* ARM_FromStartDSToResetDS(ARM_ReferenceValue* valstartRef, 
											   const ARM_DateStrip& DS)
{
	ARM_Vector* refdates = valstartRef->GetDiscreteDates(); 
	ARM_Vector* refvalues = valstartRef->GetDiscreteValues();

	std::vector<double> *schedStartDates   = DS.GetFlowStartDates();
	std::vector<double> *schedResetDates   = DS.GetResetDates();

	if ( !schedStartDates || !schedResetDates )
	{
		throw Exception (__LINE__, __FILE__, ERR_OBJECT_NULL, "ARM_FromStartRefToResetRef: Error generating schedule");
	}		

	ARM_ReferenceValue* valresetRef = NULL;

	if (refdates)
	{
		ARM_Vector* resetDates = new ARM_Vector(refdates->GetSize());

		// on construit le vecteur avec les dates de reset
		for (int i = 0; i < refdates->GetSize(); i++)
		{
			if ((*refdates)[i] > (*schedStartDates)[schedStartDates->size()-1])
			{
				(*resetDates)[i]= (*refdates)[i];
			}
			else
			{
				int idx = -1;
				size_t k;
				for (k=0; k<schedStartDates->size(); ++k)
				{

					if (fabs((*refdates)[i]- (*schedStartDates)[k]) < K_NEW_DOUBLE_TOL)
					{
						idx = k;
						break;
					}
				}

				// on ne la trouve pas : problème
				if (idx != -1)
				{
					(*resetDates)[i] = (*schedResetDates)[idx];
				}
				else
				{
					// eventuellement, on ajuste la date et on reteste
/*					int rule = legRef->GetIRIndex()->GetFwdRule();
					char* calendar = legRef->GetPayCalName();					
					double adjustedRef = ARM_Date(refdates->Elt(i)).GoodBusinessDay(rule, calendar).GetJulian();

					idx = schedStartDates->find(adjustedRef);
					if (idx != -1)
					{
						resetDates->Elt(i) = schedResetDates->Elt(idx);
					}
					else
					{
						char msg[200];
						char strDate[11];
						ARM_Date(refdates->Elt(i)).JulianToStrDate(strDate);
						sprintf(msg, "ARM_FromStartRefToResetRef: Expiry[%d] (%s) does not match any start date", i, strDate);
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg);
					}
*/				}
			}
		}

		valresetRef = new ARM_ReferenceValue((ARM_Vector *)resetDates->Clone(),(ARM_Vector *)refvalues->Clone());

		valresetRef->SetCalcMethod(K_STEPUP_LEFT);

		if (resetDates)
			delete resetDates;
		resetDates = NULL;
	}
	else
	{
		valresetRef = (ARM_ReferenceValue*) valstartRef->Clone();
	}

	return valresetRef;
}

char* ARM_FromIndexTypeToStringIndex(const ARM_INDEX_TYPE aIndexType)
{
	char* indexType = new char [20];

    switch(aIndexType)
    {
        case CMS1:
        {
            strcpy(indexType, "CMS1Y");
        };
        break;
 
        case CMS2:
        {
            strcpy(indexType, "CMS2Y");
        };
        break;
 
        case CMS3:
        {
            strcpy(indexType, "CMS3Y");
        };
        break;
 
        case CMS4:
        {
            strcpy(indexType, "CMS4Y");
        };
        break;
 
        case CMS5:
        {
            strcpy(indexType, "CMS5Y");
        };
        break;
 
        case CMS6:
        {
            strcpy(indexType, "CMS6Y");
        };
        break;
 
        case CMS7:
        {
            strcpy(indexType, "CMS7Y");
        };
        break;
 
        case CMS8:
        {
            strcpy(indexType, "CMS8Y");
        };
        break;
 
        case CMS9:
        {
            strcpy(indexType, "CMS9Y");
        };
        break;
 
        case CMS10:
        {
            strcpy(indexType, "CMS10Y");
        };
        break;
 
        case CMS11:
        {
            strcpy(indexType, "CMS11Y");
        };
        break;
 
        case CMS12:
        {
            strcpy(indexType, "CMS12Y");
        };
        break;
 
        case CMS13:
        {
            strcpy(indexType, "CMS13Y");
        };
        break;
 
        case CMS14:
        {
            strcpy(indexType, "CMS14Y");
        };
        break;
 
        case CMS15:
        {
            strcpy(indexType, "CMS15Y");
        };
        break;
 
 
        case CMS16:
        {
            strcpy(indexType, "CMS16Y");
        };
        break;
 
        case CMS17:
        {
            strcpy(indexType, "CMS17Y");
        };
        break;
 
        case CMS18:
        {
            strcpy(indexType, "CMS18Y");
        };
        break;
 
        case CMS19:
        {
            strcpy(indexType, "CMS19Y");
        };
        break;
 
        case CMS20:
        {
            strcpy(indexType, "CMS20Y");
        };
        break;
 
        case CMS21:
        {
            strcpy(indexType, "CMS21Y");
        };
        break;
 
        case CMS22:
        {
            strcpy(indexType, "CMS22Y");
        };
        break;
 
        case CMS23:
        {
            strcpy(indexType, "CMS23Y");
        };
        break;
 
        case CMS24:
        {
            strcpy(indexType, "CMS24Y");
        };
        break;
 
        case CMS25:
        {
            strcpy(indexType, "CMS25Y");
        };
        break;
 
        case CMS26:
        {
            strcpy(indexType, "CMS26Y");
        };
        break;
 
        case CMS27:
        {
            strcpy(indexType, "CMS27Y");
        };
        break;
 
        case CMS28:
        {
            strcpy(indexType, "CMS28Y");
        };
        break;
 
        case CMS29:
        {
            strcpy(indexType, "CMS29Y");
        };
        break;
 
        case CMS30:
        {
            strcpy(indexType, "CMS30Y");
        };
        break;
 
        case CMT2:
        {
            strcpy(indexType, "CMT2Y");
        };
        break;
 
        case CMT5:
        {
            strcpy(indexType, "CMT5Y");
        };
        break;
 
        case CMT10:
        {
            strcpy(indexType, "CMT10Y");
        };
        break;
 
        case CMT15:
        {
            strcpy(indexType, "CMT15Y");
        };
        break;
 
        case CMT20:
        {
            strcpy(indexType, "CMT20Y");
        };
        break;
 
        case CMT30:
        {
            strcpy(indexType, "CMT30Y");
        };
        break;
 
        case TEC5:
        {
            strcpy(indexType, "TEC5");
        };
        break;
 
        case TEC10:
        {
            strcpy(indexType, "TEC10");
        };
        break;
 
        case LIBOR1M:
        {
            strcpy(indexType, "LIBOR1M");
        };
        break;
 
        case PIBOR1M:
        {
            strcpy(indexType, "PIBOR1M");
        };
        break;
 
        case EURIBOR1M:
        {
            strcpy(indexType, "EURIBOR1M");
        };
        break;
 
        case LIBOR3M:
        {
            strcpy(indexType, "LIBOR3M");
        };
        break;
 
        case PIBOR3M:
        {
            strcpy(indexType, "PIBOR3M");
        };
        break;
 
        case EURIBOR3M:
        {
            strcpy(indexType, "EURIBOR3M");
        };
        break;
 
        case LIBOR6M:
        {
            strcpy(indexType, "LIBOR6M");
        };
        break;
 
        case PIBOR6M:
        {
            strcpy(indexType, "PIBOR6M");
        };
        break;
 
        case EURIBOR6M:
        {
            strcpy(indexType, "EURIBOR6M");
        };
        break;
 
        case LIBOR1Y:
        {
            strcpy(indexType, "LIBOR1Y");
        };
        break;
 
        case PIBOR1Y:
        {
            strcpy(indexType, "PIBOR1Y");
        };
        break;
 
        case EURIBOR1Y:
        {
            strcpy(indexType, "EURIBOR1Y");
        };
        break;
 
        case T4M:
        {
            strcpy(indexType, "T4M");
        };
        break;
 
        case T4M_FIXED:
        {
            strcpy(indexType, "T4M_FIXED");
        };
        break;
 
        case TAM:
        {
            strcpy(indexType, "TAM");
        };
        break;
 
        case TAG:
        {
            strcpy(indexType, "TAG");
        };
        break;
 
        case EONIA: // idem TMP:
        {
            strcpy(indexType, "TMP");
        };
        break;

        case IDXFIXED:
        {
            strcpy(indexType, "FIXED");
        };
        break;
 
        default:
        {
            strcpy(indexType, "  ");
        };
        break;
    }
	return indexType;
}

// Calypso function

int calypso2ARMDaycount(const char* daycount)
{
    if ( strcmp(daycount,"30/360") == 0)
        return K30_360;
    else if (strcmp(daycount,"ACT/ACT")==0)
        return KACTUAL_ACTUAL;
    else if (strcmp(daycount,"ACT/365") == 0)
        return KACTUAL_365;
    else if (strcmp(daycount,"ACT/360") == 0)
        return KACTUAL_360;
    else if (strcmp(daycount,"30/360") == 0)
        return K30_360;
    else if ( strcmp(daycount,"ACT/ACT29") == 0)
        return KACTUAL_FEB29;
    else if (strcmp(daycount,"COUP") == 0)
        return KCOUPON;
	else if ( strcmp(daycount,"30E/360") == 0 )
        return K30_360E;

	char errorMsg[128];
	sprintf(errorMsg, "Invalid Calypso Daycount : %s", daycount);
    throw Exception(__LINE__,__FILE__,ERR_INVALID_ARGUMENT,errorMsg);

}

int calypso2ARMFreq(const char* freq)
{
	if ( strcmp(freq,"WK") == 0 )
        return K_WEEKLY;
   	else if ( strcmp(freq,"MTH") == 0 )
        return K_MONTHLY;
	else if ( strcmp(freq,"BIM") == 0 )
        return K_BIMONTHLY;
	else if ( strcmp(freq,"QTR") == 0 )
        return K_QUARTERLY;
	else if ( strcmp(freq,"SA") == 0 )
        return K_SEMIANNUAL;
	else if ( strcmp(freq,"PA") == 0 )
        return K_ANNUAL;
	else if ( strcmp(freq,"ZC") == 0 )
        return K_ZEROCOUPON;
	else if ( strcmp(freq,"DLY") == 0 )
        return K_DAILY;

	char errorMsg[128];
	sprintf(errorMsg, "Invalid Calypso frequency : %s", freq);
    throw Exception(__LINE__,__FILE__,ERR_INVALID_ARGUMENT,errorMsg);

}

long calypso2ARMStubRule (const char* aRule)
{
    char temp[20];
    strcpy(temp,aRule);
    char* tmp = ARM_UPPERCASE(temp);

    if(strcmp(tmp, "LONG FIRST") ==0)
        return K_LONGSTART;
    
	if(strcmp(tmp, "SHORT LAST") ==0)
        return K_SHORTEND;
    
	if(strcmp(tmp, "LONG LAST") ==0)
        return K_LONGEND;
    
	return K_SHORTSTART;	// "SHORT FIRST","FULL COUPON", "NONE"
}

int calypso2ARMPayOrRec(int payRec)
{
	if (payRec == 1)
		return K_PAY;
	else
		return K_RCV;
}

int calypso2ARMTiming(const char* calypsoTiming) 
{
	if (strcmp(calypsoTiming, "BEG_PER") == 0)
		return K_ADVANCE;
	if (strcmp(calypsoTiming, "END_PER") == 0)
		return K_ARREARS;
	
	return K_ADVANCE;
	char errorMsg[128];
	sprintf(errorMsg, "Invalid Calypso timing '%s'", calypsoTiming);
	throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, errorMsg);
}


int calypso2ARMInterestRule(string calypsoRule)
{
	if (calypsoRule =="ADJUSTED")
		return K_ADJUSTED;
	else if (calypsoRule == "MAT_UNADJUSTED")
		return K_MATUNADJUSTED;
	else
		return K_UNADJUSTED;
}



string calypso2ARMNotionalExchange(bool initExchange, bool finalExchange) {
			string res = "NXNONE";
			if(initExchange && finalExchange) {
				res = "NXBOTH";
			} else if(initExchange) {
				res = "NXSTART";
			} else if(finalExchange) {
				res = "NXEND";
			}
			return res;
}


int calypso2ARMIndexType(string indexName,string indexTenor) {
		string indexType = "";
		//string indexName = NULL;
		//if(index!=NULL) {
		//	indexName = index.getName();
			if(!((indexName == "CMS")||(indexName=="LIBOR")||(indexName=="EURIB")||(indexName =="EURIBOR"))) {
				indexName = "LIBOR";
			}
			indexType = indexName+indexTenor;			
		//}
		return calypso2ARMIndexType(indexType);
}	



int calypso2ARMIndexType(string calypsoIndex)
{
	if (calypsoIndex == "LIBOR1M")
		return K_LIBOR1M;
	if (calypsoIndex == "LIBOR3M")
		return K_LIBOR3M;
	if (calypsoIndex == "LIBOR6M")
		return K_LIBOR6M;
	if (calypsoIndex == "LIBOR1Y")
		return K_LIBOR1Y;
	if (calypsoIndex == "EURIB1M" || calypsoIndex == "EURIBOR1M")
		return K_EURIBOR1M;
	if (calypsoIndex == "EURIB3M" || calypsoIndex == "EURIBOR3M")
		return K_EURIBOR3M;
	if (calypsoIndex == "EURIB6M" || calypsoIndex == "EURIBOR6M")
		return K_EURIBOR6M;
	if (calypsoIndex == "EURIB1Y" || calypsoIndex == "EURIBOR1Y")
		return K_EURIBOR1Y;
	if (calypsoIndex == "CMS1Y")
		return K_CMS1;
	if (calypsoIndex == "CMS2Y")
		return K_CMS2;
	if (calypsoIndex == "CMS3Y")
		return K_CMS3;
	if (calypsoIndex == "CMS4Y")
		return K_CMS4;
	if (calypsoIndex == "CMS5Y")
		return K_CMS5;
	if (calypsoIndex == "CMS10Y")
		return K_CMS10;
	if (calypsoIndex == "CMS15Y")
		return K_CMS15;
	if (calypsoIndex == "CMS20Y")
		return K_CMS20;
	if (calypsoIndex == "CMS25Y")
		return K_CMS25;
	if (calypsoIndex == "CMS30Y")
		return K_CMS30;

	char errorMsg[128];
	sprintf(errorMsg, "Invalid rate index name : %s", calypsoIndex);
	throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, errorMsg);
}


int calypso2ARMDateRoll (string dateRollStr)
{
	if (dateRollStr=="PRECEDING")
	{
        return K_PREVIOUS;
	}
	if (dateRollStr=="MOD_PRECED")
	{
        return K_MOD_PREVIOUS;
	}
	if (dateRollStr=="FOLLOWING")
	{
        return K_FOLLOWING;
	}
	if (dateRollStr =="MOD_FOLLOW")
	{
        return K_MOD_FOLLOWING;
	}
	return K_MOD_FOLLOWING;
}



double FromIndexTypeToTerm(ARM_INDEX_TYPE indexType)
{
	/*if(IsLiborIndex(indexType))
	{
       ARM_IRIndex index(indexType);
	   return index.GetYearTerm();
	}
	else
	{
       ARM_IRIndex index(indexType,(ARM_INDEX_TYPE)K_LIBOR3M);
	   return index.GetCMYearTerm();
	}*/
	return 0;
}


/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
