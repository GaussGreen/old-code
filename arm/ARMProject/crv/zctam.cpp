/*
 * $Log: zctam.cpp,v $
 * Revision 1.3  2001/04/23 09:22:33  smysona
 * unused variables
 *
 * Revision 1.2  1999/10/28 10:14:44  mab
 * Rajout de LOG pour RCS
 *
 */


/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE         : zctam.cpp                                                   */
/*                                                                            */
/* DESCRIPTION  : This file implements the ARM_ZeroCurve class, a class for   */
/*                managing zero curves                                        */
/* DATE         : Thu Aug  1 1996                                             */
/*                                                                            */
/*                                                                            */
/*----------------------------------------------------------------------------*/


/*---- System Include ----*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/*
#include <iostream.h>
#include <fstream.h>
*/


/*---- Application Include ----*/

#include "firsttoinc.h"
#include "dates.h"
#include "interpol.h"
#include "zerocurv.h"
#include "expt.h"
//#include "ycmodel.h"
#include "containr.h"
#include "security.h"
//#include "bond.h"

#include "armglob.h"





/*----------------------------------------------------------------------------*/
/* Construction de la courbe TAM T4M                                          */
/*----------------------------------------------------------------------------*/

// Constructeur de la Courbe TAM T4M

ARM_ZeroCurve::ARM_ZeroCurve(ARM_Date& asOf, 
                             char* terms[ARM_NB_TERMS], 
                             ARM_Vector* mktData,
                             double mean_rates,
                             int raw, 
                             int Cont_Lin, 
                             ARM_Currency* ccy)
{
    Init();

    SetName(ARM_ZERO_CURVE);

    itsAsOfDate = asOf;

    // Set variables

    SetCurrencyUnit(ccy);

    ZCFromMarketRates(terms, mktData, mean_rates,  raw, Cont_Lin, ccy);
}



// Fonction de calcul d'un taux de swap equivalent

double ARM_ZeroCurve::EquivalentRate(double MktSwapRates, int nber_of_days)
{
    double new_rate;
    double base360 = 360.0 ;

    new_rate = (pow((1+MktSwapRates/100.0*nber_of_days/base360), 
                   1/(double)nber_of_days )-1 )*base360;

    return(new_rate*100.0);
}



double ARM_ZeroCurve::InvertEquivalentRate(double MktSwapRates, 
                                           int nber_of_days)
{
    double new_rate;
    double base360 = 360.0 ;

    new_rate = (pow((1+MktSwapRates/100.0/base360 ), 
                 nber_of_days)-1)*(base360/(double) nber_of_days );

    return(new_rate*100.0);
}



double ARM_ZeroCurve::DecapitalizedTAMRate(double MktSwapRates, 
                           int nber_of_days, double baseACTUAL)
{
    double new_rate;


    new_rate = (pow((1+MktSwapRates/100.0), 
                (double)nber_of_days/baseACTUAL)-1)
                *baseACTUAL/(double)nber_of_days ;

    return(new_rate*100.0);
}




// Fonction de calcul d'un pilier de reference : premier jour ouvre du mois

double ARM_ZeroCurve::RefDate(double julian_date, ARM_Currency* ccy)
{
    SetCurrencyUnit(ccy);
    char* ccyName = ccy->GetCcyName();

    ARM_Date datej(julian_date);

    datej.GoodBusinessDay(K_FOLLOWING, ccyName);

    return(datej.GetJulian());
}



double ARM_ZeroCurve::KeepFourFigures(double rate)
{
    double tmp = rate*10000.0;

    int tmp2= tmp;

    return((double) tmp2 / 10000.0 );
}



// Fonction de calcul d'une date en j+1 ouvr‹

double ARM_ZeroCurve::AddOneDayBis(double julian_date, ARM_Currency* ccy)
{
    SetCurrencyUnit(ccy);
    char* ccyName = ccy->GetCcyName();

    ARM_Date datej(julian_date);
    datej.GoodBusinessDay(K_FOLLOWING, ccyName);
    datej.AddDays(1);
    datej.GoodBusinessDay(K_FOLLOWING, ccyName);
    return datej.GetJulian();
}




// Fonction de calcul d'une jambe variable pour un swap T4M de maturit‹ 1M

double ARM_ZeroCurve::FirstFloatingT4MNPV(double last_ref_df,
                                          double start_date,
                                          double first_df,
                                          double first_date,
                                          double mean,
                                          int fwd_bwd,
                                          int Cont_Lin, 
                                          ARM_Currency* ccy,
                                          double *first_ref,
                                          double *first_fixing,
                                          double *last_df)
{
    double fwd, first_ref_df=first_df, first_ref_date, 
           start_date_df=0.0, end_date_df, act_df;

    double VAN = 0.0 ;


    SetCurrencyUnit(ccy);
    char* ccyName = ccy->GetCcyName();

    ARM_Date end_date(start_date);

    end_date.AddMonths(1);

    int known_days, remaining_fwd_days, total_days, base360=360 ;


    if ( fwd_bwd == K_FWD )        // depart le 01 du mois suivant
    {    
       double nb_days_swap = start_date - itsAsOfDate.GetJulian() ;
       double nb_days_MM = first_date - itsAsOfDate.GetJulian() ;

       start_date_df = pow(first_df, nb_days_swap/nb_days_MM); 

       double nb_days_ref = RefDate(start_date, ccy)-itsAsOfDate.GetJulian();

       first_ref_df = pow(first_df, nb_days_ref/nb_days_MM);
       first_ref_date = RefDate(start_date, ccy);

        end_date_df = Interpol(first_ref_date , 
                               RefDate(end_date.GetJulian(), ccy), 
                               end_date.GetJulian(), first_ref_df, 
                               last_ref_df, Cont_Lin, KACTUAL_360 );

        act_df = Interpol(first_ref_date, RefDate(end_date.GetJulian(),ccy), 
                          AddOneDayBis(end_date.GetJulian(),ccy), 
                          first_ref_df, last_ref_df, Cont_Lin, KACTUAL_360);
        
        fwd = ( start_date_df / end_date_df - 1.0 )  ;

    }
    else // depart le 01 du mois en cours
    {
        remaining_fwd_days =  end_date.GetJulian() - first_date ;
        known_days = itsAsOfDate.GetJulian() - start_date ;
        total_days = end_date.GetJulian() - start_date;
        
        // calcul du fwd

        end_date_df = Interpol(first_date, 
                               RefDate(end_date.GetJulian(), ccy), 
                               end_date.GetJulian(), 
                               first_df, last_ref_df, Cont_Lin, KACTUAL_360);

        fwd = ( 1.0 / end_date_df  - 1.0 ) ;

        fwd += mean/100.0 * (known_days / (double) base360 );

        //os2<< "fwd" << end_date.GetJulian()-start_date <<"/t"<< fwd<<endl;

        act_df = Interpol(first_date, RefDate(end_date.GetJulian(),ccy), 
                          AddOneDayBis(end_date.GetJulian(),ccy),
                          first_df, last_ref_df, Cont_Lin, KACTUAL_360);
    }

    VAN = fwd * act_df ;

    *first_fixing = start_date_df; // on retourne le premier fixing
    *last_df = end_date_df;        // on retourne le dernier fixing
    *first_ref = first_ref_df ;    // on retourne le premier pilier

    return(VAN);
}



// calcul de la branche variable d'un swap TAM 1Y

double ARM_ZeroCurve::FirstFloatingTAMNPV(double last_ref_df,
                                          double start_date,
                                          double first_swap_fixing,
                                          double lastT4M_ref_df,
                                          double lastT4M_ref_date,
                                          double mean,
                                          int fwd_bwd,
                                          int Cont_Lin,
                                          ARM_Currency* ccy,
                                          double *last_fixing_df)
{
    double VAN, fwd, end_date_df, act_df;

    SetCurrencyUnit(ccy);
    char* ccyName = ccy->GetCcyName();

    ARM_Date end_date(start_date);

    end_date.AddYears(1);


    if ( fwd_bwd == K_FWD ) // depart le 01 du mois suivant
    {    
       end_date_df = Interpol(lastT4M_ref_date, 
                              RefDate(end_date.GetJulian(), ccy), 
                              end_date.GetJulian(), 
                              lastT4M_ref_df, last_ref_df, Cont_Lin, 
                              KACTUAL_ACTUAL);

       act_df = Interpol(lastT4M_ref_date, RefDate(end_date.GetJulian(),ccy), 
                         AddOneDayBis(end_date.GetJulian(),ccy), 
                         lastT4M_ref_df, last_ref_df, 
                         Cont_Lin, KACTUAL_ACTUAL );

       fwd = ( first_swap_fixing / end_date_df - 1.0 );
    }
    else  // depart le 01 du mois en cours
    {
       // a verifier

       end_date_df = Interpol(lastT4M_ref_date , 
                              RefDate(end_date.GetJulian(), ccy), 
                              end_date.GetJulian(), 
                              lastT4M_ref_df, last_ref_df, 
                              Cont_Lin, KACTUAL_ACTUAL);

       double capi_fact = ( 1 + mean * (itsAsOfDate.GetJulian()
                            -start_date)/36000.0 ) ;
                                    
       fwd = (capi_fact*(1.0/end_date_df))-1.0;

       act_df = Interpol(lastT4M_ref_date, 
                         RefDate(end_date.GetJulian(),ccy), 
                         AddOneDayBis(end_date.GetJulian(),ccy), 
                         lastT4M_ref_df, last_ref_df, 
                         Cont_Lin, KACTUAL_ACTUAL);
    }
    
    VAN = fwd * act_df ;

    *last_fixing_df = end_date_df ;  // on retourne le dernier fixing

    return(VAN);
}



double ARM_ZeroCurve::FloatingTAMNPVWithRompu(
                       double last_ref_df, // pilier cherche
                       double start_date,
                       double first_swap_fixing,
                       double first_ref_df, // premier pilier de la courbe TAM
                       double first_ref_date, 
                       ARM_Matrix SwapT4M,
                       int SizeT4M,
                       double mean,
                       int nb_months,
                       int fwd_bwd,
                       int Cont_Lin,
                       ARM_Currency* ccy)
{
    double VAN, base360 = 360.0, baseACTUAL, fwd,
           end_date_reset_df, end_date_reset_df2, act_df;
    int k;

    SetCurrencyUnit(ccy);
    char* ccyName = ccy->GetCcyName();

    ARM_Date end_date(start_date);
    ARM_Date end_date_for_base(start_date);
    end_date_for_base.AddYears(1);

    end_date.AddMonths(nb_months); // on considere le premier pilier 
                                   // du swap (rompu)

    baseACTUAL = DaysBetweenDates(KACTUAL_ACTUAL, (ARM_Date) start_date, 
                                  (ARM_Date) end_date_for_base.GetJulian());

    // on cherche les deux piliers T4M encadrant cette maturite

    if ( SwapT4M.Elt(SizeT4M, 1) > end_date.GetJulian() )
    {
       k = 1;

       while ( SwapT4M.Elt(k, 1) < end_date.GetJulian() ) 
             k++;

       end_date_reset_df = Interpol(SwapT4M.Elt(k-1, 1), 
                                    SwapT4M.Elt(k, 1), 
                                    end_date.GetJulian(), 
                                    SwapT4M.Elt(k-1, 0), 
                                    SwapT4M.Elt(k, 0), 
                                    Cont_Lin, KACTUAL_ACTUAL);

       act_df = Interpol(SwapT4M.Elt(k-1, 1), SwapT4M.Elt(k, 1), 
                          AddOneDayBis (end_date.GetJulian(), ccy), 
                          SwapT4M.Elt(k-1, 0), SwapT4M.Elt(k, 0), 
                          Cont_Lin, KACTUAL_ACTUAL);
    }
    else 
    {
       end_date_reset_df = Interpol(SwapT4M.Elt(SizeT4M, 1), first_ref_date,
                                    end_date.GetJulian(), 
                                    SwapT4M.Elt(SizeT4M, 0),
                                    first_ref_df,  Cont_Lin, 
                                    KACTUAL_ACTUAL);

       act_df = Interpol(SwapT4M.Elt(SizeT4M, 1),
                         first_ref_date,
                         AddOneDayBis (end_date.GetJulian(), ccy), 
                         SwapT4M.Elt(SizeT4M, 0), first_ref_df,  
                         Cont_Lin, KACTUAL_ACTUAL);
    }


    if ( fwd_bwd == K_FWD )  // depart le 01 du mois suivant
       fwd = ( first_swap_fixing / end_date_reset_df - 1.0 )  ;
    else // depart le 01 du mois en cours
    {
       double capi_fact = ( 1 + mean * (itsAsOfDate.GetJulian()
                          -start_date)/base360/100.0) ;

       fwd = (capi_fact / end_date_reset_df) - 1;
    }
    
    VAN = fwd * act_df ;

    // on considere le second pilier du swap (on ne considere que les swaps 
    // avec rompu de mat comprises entre 1Y et 2Y)

    end_date.AddYears(1);    

    end_date_reset_df2 = Interpol(first_ref_date, 
                                  RefDate(end_date.GetJulian(), ccy), 
                                  end_date.GetJulian(), 
                                  first_ref_df, last_ref_df, 
                                  Cont_Lin, KACTUAL_ACTUAL);

    act_df = Interpol(first_ref_date, RefDate(end_date.GetJulian(), ccy), 
                      AddOneDayBis (end_date.GetJulian(), ccy), 
                      first_ref_df, last_ref_df, Cont_Lin, KACTUAL_ACTUAL);

    VAN += ( end_date_reset_df / end_date_reset_df2 - 1.0 ) * act_df ;

    return(VAN);
}



// Fonction de calcul d'une jambe variable pour un swap T4M ou TAM 

double ARM_ZeroCurve::FloatingNPV(double last_ref_df,
                                  double last_ref_date,
                                  double first_ref_df,
                                  double first_ref_date,
                                  double first_fixing_df,
                                  double first_fixing_date,
                                  double last_NPV,
                                  int nb_missing_flows,
                                  int T4M_TAM,
                                  int Cont_Lin, 
                                  ARM_Currency* ccy)
{
    double fwdk, dfk_fixing, dfk_act;
    double previous_df = first_fixing_df;
    double VAN = last_NPV ;

    ARM_Date datek(first_fixing_date);
    SetCurrencyUnit(ccy);
    char* ccyName = ccy->GetCcyName();


    for (int k = 1; k <= nb_missing_flows + 1; k++)
    {
        if ( T4M_TAM == K_T4M_INDEX )
        {
           datek.AddMonths(1);

           dfk_fixing = Interpol(RefDate(first_fixing_date, ccy), 
                                 last_ref_date, datek.GetJulian(), 
                                 first_ref_df, last_ref_df, Cont_Lin, 
                                 KACTUAL_360);

           dfk_act = Interpol(RefDate(first_fixing_date, ccy), 
                              last_ref_date, 
                              AddOneDayBis(datek.GetJulian(),ccy), 
                              first_ref_df, last_ref_df, Cont_Lin, 
                              KACTUAL_360);
        }
        else 
        {
           datek.AddYears(1);        

           dfk_fixing = Interpol(first_ref_date, last_ref_date, 
                                 datek.GetJulian(), 
                                 first_ref_df, last_ref_df, 
                                 Cont_Lin, KACTUAL_ACTUAL);

           dfk_act = Interpol(first_ref_date, last_ref_date, 
                              AddOneDayBis(datek.GetJulian(),ccy), 
                              first_ref_df, last_ref_df, 
                              Cont_Lin, KACTUAL_ACTUAL);

        }

        fwdk = previous_df / dfk_fixing - 1 ;

        VAN += fwdk * dfk_act;

        previous_df = dfk_fixing ;
    }

    return(VAN);
}



// Fonction de calcul d'une jambe fixe pour un swap T4M 

double ARM_ZeroCurve::FirstFixedNPV(double last_ref_df,
                double last_ref_date,
                double first_df,// soit df monetaire soit dernier df T4M
                double first_date,// soit date monetaire soit derniere date T4M
                double swap_start_date,
                int fwd_bwd, 
                int T4M_TAM,
                int Cont_Lin,
                ARM_Currency* ccy)
{
    double first_ref_df, first_ref_date;
    double df_act, base360=360.0, previous_date = swap_start_date;
    double VAN = 0.0 ;
    ARM_Date end_date(swap_start_date);
    int nb_days;


    if ( T4M_TAM == K_T4M_INDEX )
    {
       end_date.AddMonths(1);
       nb_days = DaysBetweenDates(KACTUAL_360, (ARM_Date) previous_date, 
                                  (ARM_Date) end_date.GetJulian());
    }
    else
    {
       end_date.AddYears(1);

       nb_days = base360;
    }
    
    if ( fwd_bwd == K_FWD )  // depart le 01 du mois suivant
    {        
       double nb_days_MM = first_date - itsAsOfDate.GetJulian() ;

       double nb_days_ref = RefDate(swap_start_date, ccy) 
                            - itsAsOfDate.GetJulian();
        
       if ( T4M_TAM == K_T4M_INDEX)
       {
          first_ref_df = pow(first_df, nb_days_ref/nb_days_MM);

          first_ref_date = RefDate(swap_start_date, ccy);
       }
       else
       {
          first_ref_df = first_df;

          first_ref_date = first_date;
       }
        
       df_act = Interpol(first_ref_date, RefDate(end_date.GetJulian(),ccy), 
                     AddOneDayBis(end_date.GetJulian(),ccy), first_ref_df, 
                     last_ref_df, Cont_Lin, (T4M_TAM == K_T4M_INDEX ? 
                            KACTUAL_360 :  KACTUAL_ACTUAL) );
    }
    else // depart le 01 du mois en cours
    {
       if ( T4M_TAM == K_T4M_INDEX)
          df_act = Interpol(first_date, RefDate(end_date.GetJulian(),ccy), 
                            AddOneDayBis(end_date.GetJulian(),ccy),
                            first_df, last_ref_df, Cont_Lin, KACTUAL_360);    
        else
           df_act = Interpol(first_date, RefDate(end_date.GetJulian(),ccy), 
                             AddOneDayBis(end_date.GetJulian(),ccy),
                             first_df, last_ref_df, Cont_Lin, 
                             KACTUAL_ACTUAL);    
    }
        
    VAN += df_act * nb_days/base360;

    return(VAN);
}



double ARM_ZeroCurve::FixedTAMNPVWithRompu(
                     double last_ref_df,  // pilier cherche
                     double start_date,
                     double first_ref_df, // premier pilier de la courbe TAM
                     double first_ref_date, 
                     double swap_rate,
                     ARM_Matrix SwapT4M,
                     int SizeT4M,
                     int nb_months,
                     int Cont_Lin,
                     ARM_Currency* ccy)
{
    double VAN, act_df, nb_days, baseACTUAL;
    int k;

    SetCurrencyUnit(ccy);
    char* ccyName = ccy->GetCcyName();

    ARM_Date end_date(start_date);
    ARM_Date end_date_for_base(start_date);
    end_date_for_base.AddYears(1);

    end_date.AddMonths(nb_months); // on considere le premier pilier 
                                   // du swap (rompu)

    // on cherche les deux piliers T4M encadrant cette maturite

    if ( SwapT4M.Elt(SizeT4M, 1) > end_date.GetJulian() )
    {
       k = 1;

       while ( SwapT4M.Elt(k, 1) < end_date.GetJulian() ) 
             k++;

       act_df = Interpol(SwapT4M.Elt(k-1, 1) , SwapT4M.Elt(k, 1), 
                          AddOneDayBis (end_date.GetJulian(), ccy), 
                          SwapT4M.Elt(k-1, 0), SwapT4M.Elt(k, 0), 
                          Cont_Lin, KACTUAL_ACTUAL);
    }
    else 
       act_df = Interpol(SwapT4M.Elt(SizeT4M, 1), first_ref_date,
                          AddOneDayBis (end_date.GetJulian(), ccy), 
                          SwapT4M.Elt(SizeT4M, 0), first_ref_df,
                          Cont_Lin, KACTUAL_ACTUAL);

    nb_days = DaysBetweenDates(KACTUAL_ACTUAL, (ARM_Date) start_date, 
                               (ARM_Date) end_date.GetJulian());

    baseACTUAL = DaysBetweenDates(KACTUAL_ACTUAL, (ARM_Date) start_date, 
                                  (ARM_Date) end_date_for_base.GetJulian());

    VAN = act_df*DecapitalizedTAMRate(swap_rate, nb_days, 
                             baseACTUAL)/100.0 * nb_days / baseACTUAL;

    // on considere le second pilier du swap (on ne considere 
    // que les swaps avec rompu de mat comprises entre 1Y et 2Y)

    ARM_Date end_date2(end_date);
    end_date2.AddYears(1);    

    act_df = Interpol(first_ref_date, RefDate(end_date2.GetJulian(), ccy), 
                      AddOneDayBis (end_date2.GetJulian(), ccy), 
                      first_ref_df, last_ref_df,  Cont_Lin, KACTUAL_ACTUAL);
        
    VAN += act_df * swap_rate /100.0;

    return(VAN);
}



// Fonction de calcul d'une jambe fixe pour un swap T4M ou TAM

double ARM_ZeroCurve::FixedNPV(double last_ref_df,
                               double last_ref_date,
                               double first_ref_df,
                               double first_ref_date, 
                               double last_swap_reset_date,
                               double last_NPV,
                               int nb_missing_flows,
                               int T4M_TAM,
                               int Cont_Lin, 
                               ARM_Currency* ccy)
{
    double dfk_plus_un, base, base360 = 360.0, base365 = 365.0, 
           previous_date = last_swap_reset_date;

    double VAN = last_NPV ;

    ARM_Date datek(last_swap_reset_date);
    ARM_Date datek_moins_un;
    int nb_days;

    for (int k = 1; k <= nb_missing_flows + 1; k++)
    {
        if ( T4M_TAM == K_T4M_INDEX )
        {
           datek.AddMonths(1);

           dfk_plus_un = Interpol(RefDate(last_swap_reset_date, ccy), 
                                  last_ref_date, 
                                  AddOneDayBis(datek.GetJulian(), ccy),  
                                  first_ref_df, last_ref_df, 
                                  Cont_Lin, KACTUAL_360);    

           nb_days = DaysBetweenDates(KACTUAL_360, (ARM_Date) previous_date, 
                     (ARM_Date) datek.GetJulian());

           base = base360;
        }
        else
        {
           datek.AddYears(1);

           nb_days = base360; 

           base = base360;

           dfk_plus_un = Interpol(first_ref_date, last_ref_date, 
                                  AddOneDayBis(datek.GetJulian(), ccy),  
                                  first_ref_df, last_ref_df, Cont_Lin, 
                                  KACTUAL_ACTUAL);    
        }

        previous_date = datek.GetJulian() ;

        VAN += dfk_plus_un * nb_days / base;
    }

    return(VAN);
}
                                  



// Construction de la Courbe TAM_T4M : calcul des discount facteurs

void ARM_ZeroCurve::ZCFromMarketRates(char* Terms[ARM_NB_TERMS], 
                                      ARM_Vector* data, 
                                      double mean_rates, 
                                      int raw, int Cont_Lin, 
                                      ARM_Currency* ccy)
{
    int SizeM = 0, SizeT4M = 1, SizeTAM = 1, i = 0, Nb;
    int Tx_Spot = 0;
    ARM_Date matDate;

    int CurrentMonth;
    int CurrentDay;

    ARM_Matrix MMRates(30, 2, 0.0);
    ARM_Matrix T4MSwapRates(50, 2, 0.0);
    ARM_Matrix TAMSwapRates(80, 3, 0.0);

    SetCurrencyUnit(ccy);

    char* ccyName = ccy->GetCcyName();
    int fwdRule = ccy->GetFwdRule();
    int spotDays = ccy->GetSpotDays();

    char matu;
    int month = 0, k=0;
    double first_fixing_df;
    double base360 = 360.0;
    double baseACTUAL;
    double pseudo_T4M_swap_value, pseudo_T4M_swap_mat;
    
    
    CurrentMonth = itsAsOfDate.GetMonth();
    CurrentDay = itsAsOfDate.GetDay();
    ARM_Date startSwapDate(1, CurrentMonth, itsAsOfDate.GetYear());

    if ( CurrentDay >= MID_MONTH_DAY )
       startSwapDate.AddMonths(1);

        
    matu = 'X';

    while ( Terms[i][0] != 'X' )
    {          
        sscanf(Terms[i], "%d%c", &Nb, &matu);

        matu = toupper(matu);
        
       if (data->Elt(i) < 1.0E-6)
       {
          throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Non zero rates are required for building zero curve");
       }

        // Conversion des Termes char -> Termes double

        if ( matu == 'D' )            // on renseigne le taux "1D"
        {    
           matDate = itsAsOfDate;
           matDate.NextBusinessDay(Nb, ccyName);    
           MMRates.Elt(SizeM, 0) = data->Elt(i);            
           MMRates.Elt(SizeM, 1) = matDate.GetJulian();            
           Tx_Spot = 1;

           SizeM++;
        }
        //  ->implicitement ce sont des taux de swap T4M
        else if (( matu == 'M' ) && ( Nb < TAM_PERIOD_IN_MONTH ) )        
        {   //  Ex : "9M"    

           matDate = startSwapDate;          

           matDate.AddMonths(Nb);

           T4MSwapRates.Elt(Nb, 0) = data->Elt(i) ;                                              
           T4MSwapRates.Elt(Nb, 1) = matDate.GetJulian();

           SizeT4M++;
        }
        //  ->implicitement ce sont des taux de swap TAM
        else if (( matu == 'Y')    || (( matu == 'M' ) 
                && ( Nb > TAM_PERIOD_IN_MONTH ) && 
                   ( Nb < 2*TAM_PERIOD_IN_MONTH )))    
        {   
            //  Ex : "15Y"
            matDate = startSwapDate;

            if ( matu == 'Y')
            {
               matDate.AddYears(Nb);

               TAMSwapRates.Elt(Nb + k, 1) = matDate.GetJulian();
               TAMSwapRates.Elt(Nb + k, 0) = data->Elt(i);
            }
            else 
            {
               int nb_years = (int)Nb/12.0;

               matDate.AddYears(nb_years);

               matDate.AddMonths(Nb - nb_years*12);

               // on decale d'autant l'indicage des swaps et 
               // leur maturite en nombre d'annees

               k++; 

               // on indique qu'il s'agit d'une p‹riode bris‹e

               TAMSwapRates.Elt(nb_years + k, 2) = (int)(Nb - 12 * nb_years) ; 
  
               TAMSwapRates.Elt(nb_years + k, 1) = matDate.GetJulian();
    
               TAMSwapRates.Elt(nb_years + k, 0) = data->Elt(i);
            }

            SizeTAM++;
        }
        else  
           break;

        i++;
    }

    SizeM = 1;            // on ne prend en compte que le taux JJ

    if ( Tx_Spot == 0 ) // pas de taux spot
    {
       char BUF[100];


       sprintf(BUF,"Spot rate not found ");

       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,BUF);
    }

    // on utilise le taux JJ pour le calcul des T4M
            
    // il faut le taux swap T4M 1M 
    if (( T4MSwapRates.Elt(1, 0) < 0.0001 ))
    {
       char BUF[100];

       sprintf(BUF, "1 month Swap rate not found");

       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, BUF);
    }

    // il faut le taux swap TAM 1Y 
    if (( TAMSwapRates.Elt(1, 0) < 0.0001 ))
    {
       char BUF[100];

       sprintf(BUF, "1 year Swap rate not found");

       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, BUF);
    }

    // calcul du pseudo T4M 12M 
    pseudo_T4M_swap_mat = TAMSwapRates.Elt(1, 1);

    baseACTUAL = pseudo_T4M_swap_mat
                  -startSwapDate.GetJulian(); // base ACTUAL vaut 365 ou 366

    pseudo_T4M_swap_value =  base360/baseACTUAL * TAMSwapRates.Elt(1, 0);


    // Les dates affectees aux instruments sont celles correspondantes 
    // aux dates de reset des taux var 
    // (premier jour du mois)
    // et non aux dates des piliers
    // c'est a la fin du calcul qu'elles seront converties en dates 
    // piliers (j ouvr‹)

    for (i = 1; i < T4MSwapRates.GetNumLines(); i++) 
    {
        if (T4MSwapRates.Elt(i, 1) > 0.0) 
           SizeT4M = i; // SizeT4M = maturite max de la courbe
                        // indexation entre 1 et n

    }

    for (i = 1; i <= SizeT4M; i++)                       
    {
        if ( T4MSwapRates.Elt(i, 1) < 0.0001 ) 
        {
           matDate = startSwapDate;

           matDate.AddMonths(i);// Flux mensuels necessairement

           T4MSwapRates.Elt(i, 1) = matDate.GetJulian();// On reserve la 
                                                        // place d'un swap 
                                                        // non sur un pilier
         }
    }

    for (i = 1; i < TAMSwapRates.GetNumLines(); i++) 
    {
        if ( TAMSwapRates.Elt(i, 1) > 0.0 )
           SizeTAM = i; // SizeTAM = maturite max de la courbe + nbre de 
                        // swaps brises
    }                   // indexation entre 1 et n

    k = 0;

    for (i = 1; i <= SizeTAM; i++)
    {
        if ( TAMSwapRates.Elt(i, 2) > 0.0001 ) 
           k++; // indicage des swaps a periode brisee

        if ( TAMSwapRates.Elt(i, 1) < 0.0001 ) 
        {
           matDate = startSwapDate;
           matDate.AddYears(i-k);

           // on reserve la place d'un swap non sur un pilier
           TAMSwapRates.Elt(i, 1) = matDate.GetJulian(); 
        }
    }

    // Calcul des taux ZC implicites dans chaque segment du marche
    
    CptMMZeroRates(&MMRates, SizeM);

    // prise en compte du df JJ pour le calcul du permier 
    // pilier des swaps T4M

    T4MSwapRates.Elt(0, 0) = MMRates.Elt( 0, 0);
    T4MSwapRates.Elt(0, 1) = MMRates.Elt( 0, 1);

    CptSwapT4MZeroRates(mean_rates, startSwapDate, &first_fixing_df,
                        &T4MSwapRates, 
                        SizeT4M, pseudo_T4M_swap_value, pseudo_T4M_swap_mat, 
                        raw, Cont_Lin, ccy);

    CptSwapTAMZeroRates(mean_rates, startSwapDate, first_fixing_df, 
                        T4MSwapRates, 
                        SizeT4M, &TAMSwapRates, SizeTAM, 
                        raw, Cont_Lin, ccy);

    // Preparation du tableau des taux ZC correspondant a GenerateCurve

    int zcSize = 0;

    for (i = 0; i < SizeM; i++)
        zcSize++;

    if ( MMRates.Elt(SizeM-1, 1) == T4MSwapRates.Elt(1, 1))
    {
       for (i = 2; i <= SizeT4M; i++)
           zcSize++;
    }
    else
    {
       for (i = 1; i <= SizeT4M; i++)
           zcSize++;
    }

    for (i = 1; i <= SizeTAM; i++)
        zcSize++;

    itsDateTerms = new ARM_Vector (zcSize, 0.0);
    itsYearTerms = new ARM_Vector (zcSize, 0.0);
    itsZeroRates = new ARM_Vector (zcSize, 0.0);
    itsDiscountFactors = new ARM_Vector(zcSize, 0.0);

    zcSize = 0;

    for (i = 0; i < SizeM; i++)
    {
        itsDateTerms->Elt(zcSize) = MMRates.Elt(i, 1) - itsAsOfDate.GetJulian();     
        itsDiscountFactors->Elt(zcSize) = MMRates.Elt(i, 0);     

        zcSize++;
    }

    // on ne prend pas en compte le premier pilier des swaps T4M

    if ( MMRates.Elt(SizeM-1, 1) == T4MSwapRates.Elt(1, 1))
    {
       for (i = 2; i <= SizeT4M; i++)
       {
           itsDateTerms->Elt(zcSize) = T4MSwapRates.Elt(i, 1) 
                                       - itsAsOfDate.GetJulian();     

           itsDiscountFactors->Elt(zcSize) = T4MSwapRates.Elt(i, 0);     

           zcSize++;
        }
    }
    else
    {
       for (i = 1; i <= SizeT4M; i++)
       {
           itsDateTerms->Elt(zcSize) = T4MSwapRates.Elt(i, 1) 
                                       - itsAsOfDate.GetJulian();     

           itsDiscountFactors->Elt(zcSize) = T4MSwapRates.Elt(i, 0);     

           zcSize++;
       }
    }


    for (i = 1; i <= SizeTAM; i++)
    {
        itsDateTerms->Elt(zcSize) = TAMSwapRates.Elt(i, 1) 
                                    - itsAsOfDate.GetJulian();

        itsDiscountFactors->Elt(zcSize) = TAMSwapRates.Elt(i, 0);

        zcSize++;
    }

    for (i = 0; i < zcSize; i++)
    {    
        itsYearTerms->Elt(i) = itsDateTerms->Elt(i)/365.0; 

        itsZeroRates->Elt(i) = -100.0*log(itsDiscountFactors->Elt(i))
                               /itsYearTerms->Elt(i);
    }
}



/*----------------------------------------------------------------------------*
 Construction des piliers sur la courbe T4M
 cette fonction renvoie la valeur du dernier DF T4M et de la date du pilier 
 correspondante cette fonction renvoie aussi la valeur du premier fixing 
 qui sera necessaire pour le calcul des swaps TAM
*----------------------------------------------------------------------------*/

void ARM_ZeroCurve::CptSwapT4MZeroRates(double mean_rates, 
                                        ARM_Date startSwapDate, 
                                        double* first_fixing_df, 
                                        ARM_Matrix *SwapT4M, int &SizeT4M, 
                                        double pseudo_T4M_swap_value, 
                                        double pseudo_T4M_swap_mat,
                                        int raw, int Cont_Lin, 
                                        ARM_Currency* ccy)
{
    int i, j;
    int base360 = 360;
    int fwd_bwd;
    double newton_result, last_df, floatNPV, fixedNPV, 
           der_floatNPV, der_fixedNPV;

    double TMP1, TMP2, TMP3, tmp1, tmp2, tmp3, first_fixing, first_ref_df, 
           last_fixing, lastFloatNPV=0.0, lastFixedNPV=0.0;

    double Df_Ctr;

    char* ccyName = ccy->GetCcyName();

    int missing_flows, last_missing_flows; 

    int m1=startSwapDate.GetMonth();
    int m2=itsAsOfDate.GetMonth();
    
    if ( m1 == m2 )
        fwd_bwd = K_BWD;
    else
       fwd_bwd = K_FWD;


    // si methode PAR, on calcule les taux de swap manquant par 
    // interpolation lineaire
    // des taux de swap renseignes 
    // au delÌ des taux de march‹ on utilise le pseudo taux T4M 12 M 

    if ( raw == K_PAR )
    {
       for (i = 1; i <= SizeT4M; i++)
       {
           if ( SwapT4M->Elt(i, 0) < ISPRESENT )
           {
              j = i+1;

              while ( SwapT4M->Elt(j, 0) < ISPRESENT ) 
                    j++;

              double first_swap = SwapT4M->Elt(i-1,0);
              double scd_swap = SwapT4M->Elt(j,0); 

              double tmp_swap = first_swap+(scd_swap-first_swap)  
                            *DaysBetweenDates(K30_360, 
                 (ARM_Date) SwapT4M->Elt(i-1, 1),(ARM_Date) SwapT4M->Elt(i, 1))
                 /DaysBetweenDates(K30_360, 
                (ARM_Date) SwapT4M->Elt(i-1, 1),(ARM_Date) SwapT4M->Elt(j, 1));

              SwapT4M->Elt(i, 0) = tmp_swap;
           }
       }

       j = SizeT4M +1;

       while ( j < TAM_PERIOD_IN_MONTH )
       {
           double first_swap = SwapT4M->Elt(j-1,0) ; 
           double scd_swap = pseudo_T4M_swap_value ; 
           ARM_Date date(SwapT4M->Elt(j-1,1));

           date.AddMonths(1);
           SwapT4M->Elt(j, 1) = date.GetJulian();

           double tmp_swap = first_swap + (scd_swap-first_swap)
                    *DaysBetweenDates(K30_360, 
                  (ARM_Date) SwapT4M->Elt(j-1, 1),(ARM_Date) SwapT4M->Elt(j, 1))
               /DaysBetweenDates(K30_360, (ARM_Date) SwapT4M->Elt(j-1, 1),
               (ARM_Date) pseudo_T4M_swap_mat);

            SwapT4M->Elt(j, 0) = tmp_swap ;

            j++;
            SizeT4M++;
       }
    }
    // Conversion des taux en taux ‹quivalents
    for (i=1 ; i <= SizeT4M ; i++)
        SwapT4M->Elt(i, 0) = EquivalentRate(SwapT4M->Elt(i, 0), 
                  (int) SwapT4M->Elt(i, 1) - startSwapDate.GetJulian());

    // Calcul des discount facteurs
        
    for (i = 1 ; i <= SizeT4M ; i++)
    {
        if ( (SwapT4M->Elt(i, 0) < ISPRESENT) && (raw == K_RAW) )
        {    // Cas des swaps absents dans la methode raw
             // on part du principe que les taux sont connus entre i-1 et j-1

            j = i+1;

            while ( SwapT4M->Elt(j, 0) < ISPRESENT ) 
                  j++;    

            missing_flows = j-i;
        }
        else
        {
           j = i;

           missing_flows = 0.0;
        }

        // on cherche a egaliser les branches fixes et variables des swaps

        Df_Ctr = FIRST_DF_FOR_NEWTON;

        do
        {
            // cas du premier flux

            if (i == 1)
            {
               ARM_Date end_date=startSwapDate;

               end_date.AddMonths(1);

               floatNPV = FirstFloatingT4MNPV(Df_Ctr, 
                        startSwapDate.GetJulian(),  SwapT4M->Elt(0,0), 
                        SwapT4M->Elt(0,1), 
                        mean_rates, fwd_bwd, Cont_Lin, ccy, 
                        &first_ref_df, &first_fixing, &last_fixing);

               TMP2 =last_fixing ;
               TMP1 =first_ref_df ;
               TMP3 = first_fixing ;

               der_floatNPV = (FirstFloatingT4MNPV(Df_Ctr+H_FOR_DERIVATIVE, 
                             startSwapDate.GetJulian(), SwapT4M->Elt(0,0), 
                             SwapT4M->Elt(0,1), mean_rates, fwd_bwd, Cont_Lin, 
                      ccy, &tmp1, &tmp2, &tmp3) - floatNPV )/H_FOR_DERIVATIVE;
                        

               fixedNPV = FirstFixedNPV(Df_Ctr, RefDate(SwapT4M->Elt(j, 1),
                                                     ccy), SwapT4M->Elt(0, 0),
                           SwapT4M->Elt(0, 1), startSwapDate.GetJulian(), 
                           fwd_bwd,
                           K_T4M_INDEX, Cont_Lin, ccy);
                    
               der_fixedNPV = (FirstFixedNPV(Df_Ctr+H_FOR_DERIVATIVE, 
                        RefDate(SwapT4M->Elt(j, 1),ccy), SwapT4M->Elt(0, 0),
                        SwapT4M->Elt(0, 1), startSwapDate.GetJulian(), 
                        fwd_bwd, K_T4M_INDEX, 
                        Cont_Lin,  ccy)-fixedNPV)/H_FOR_DERIVATIVE;
            }
            else // cas des autres flux
            {
               if (i==2)
                  last_df = TMP2;  // valeur du dernier fixing du premier flux
               else 
                  last_df = Interpol(RefDate(SwapT4M->Elt(i-2
                                     -last_missing_flows, 1), ccy),  
                     RefDate(SwapT4M->Elt(i-1, 1), ccy), SwapT4M->Elt(i-1, 1), 
                     SwapT4M->Elt(i-2-last_missing_flows, 0), 
                         SwapT4M->Elt(i-1, 0), Cont_Lin, KACTUAL_360);


               floatNPV = FloatingNPV(Df_Ctr, 
                         RefDate(SwapT4M->Elt(j, 1),ccy), SwapT4M->Elt(i-1, 0),
                         RefDate(SwapT4M->Elt(i-1, 1),ccy),
                         last_df, SwapT4M->Elt(i-1, 1), lastFloatNPV, 
                         missing_flows, 
                         K_T4M_INDEX, Cont_Lin, ccy);

                der_floatNPV = ( FloatingNPV (Df_Ctr+H_FOR_DERIVATIVE, 
                                  RefDate(SwapT4M->Elt(j, 1),ccy), 
                               SwapT4M->Elt(i-1, 0), 
                              RefDate(SwapT4M->Elt(i-1, 1),ccy),last_df, 
                              SwapT4M->Elt(i-1, 1), lastFloatNPV, 
                              missing_flows, 
                              K_T4M_INDEX, Cont_Lin, ccy)
                              - floatNPV )/H_FOR_DERIVATIVE;
 

                fixedNPV = FixedNPV( Df_Ctr, RefDate(SwapT4M->Elt(j, 1),ccy), 
                           SwapT4M->Elt(i-1, 0),
                           RefDate(SwapT4M->Elt(i-1, 1),ccy),
                           SwapT4M->Elt(i-1, 1), lastFixedNPV, 
                           missing_flows,  K_T4M_INDEX, Cont_Lin,  ccy);
                    
                der_fixedNPV =(FixedNPV(Df_Ctr+H_FOR_DERIVATIVE,
                               RefDate(SwapT4M->Elt(j, 1),ccy), 
                               SwapT4M->Elt(i-1, 0), 
                               RefDate(SwapT4M->Elt(i-1, 1),ccy),
                               SwapT4M->Elt(i-1, 1), lastFixedNPV, 
                               missing_flows, K_T4M_INDEX, 
                               Cont_Lin,  ccy) - fixedNPV)/H_FOR_DERIVATIVE;
            }
                    
            newton_result = ( SwapT4M->Elt(j, 0)/100.0*fixedNPV 
                             - floatNPV)/(SwapT4M->Elt(j, 0)
                             /100.0*der_fixedNPV - der_floatNPV ) ;

            Df_Ctr = Df_Ctr - newton_result;
         }
         while ( fabs(newton_result) > PRECISION );

         lastFloatNPV = floatNPV;
         lastFixedNPV = fixedNPV;
         i=i+missing_flows;
         SwapT4M->Elt(i, 0) = Df_Ctr;
         SwapT4M->Elt(0, 0) = TMP1;
        last_missing_flows=missing_flows;
    }
    
    SwapT4M->Elt(0,0) = TMP1; // valeur du premier pilier correspondant 
                              // au depart fwd des swaps

    SwapT4M->Elt(0,1) = RefDate(startSwapDate.GetJulian(),ccy); 

    *first_fixing_df = TMP3; // valeur du premier fixing necessaire 
                             // aux calculs des swaps TAM

    // si meth = RAW, virer les taux manquant initialement 
    if ( raw == K_RAW )
    {
       for (i=1; i<=SizeT4M; i++)
       {
           if (SwapT4M->Elt(i, 0)<0.0001) // ils n'ont pas ete calcule
           {
              for (j=i; j<SizeT4M; j++)
              {
                  SwapT4M->Elt(j, 0) = SwapT4M->Elt(j+1, 0);
                  SwapT4M->Elt(j, 1) = SwapT4M->Elt(j+1, 1);
              }

              SizeT4M--;
              i--;
            }
        }
    }


    // on rajoute un pilier correspondant Ì la date de depart des swaps
    if (fwd_bwd == K_FWD)
    {
        for (j=SizeT4M; j>=0; j--){
                SwapT4M->Elt(j+1, 0) = SwapT4M->Elt(j, 0);
                SwapT4M->Elt(j+1, 1) = SwapT4M->Elt(j, 1);
        }
        SizeT4M+=1;
    }


    // on calcule la date de maturit‹ du pilier 
    // (premier jour ouvre du mois)

    for (i=1; i<=SizeT4M; i++)
        SwapT4M->Elt(i, 1) = RefDate(SwapT4M->Elt(i, 1), ccy);

}
    


/*----------------------------------------------------------------------------*/
/* Construction des piliers sur la courbe TAM                                 */
/*----------------------------------------------------------------------------*/

void ARM_ZeroCurve::CptSwapTAMZeroRates(double mean_rates, 
                                        ARM_Date startSwapDate, 
                                        double first_fixing_df,
                                        ARM_Matrix SwapT4M, int SizeT4M,
                                        ARM_Matrix *SwapTAM, int &SizeTAM, 
                                        int raw, int Cont_Lin, 
                                        ARM_Currency* ccy)
{
    int i, j, k;
    int base360 = 360;
    int fwd_bwd;
    double newton_result, last_df, last_fixing_df, floatNPV, 
           fixedNPV, der_floatNPV, der_fixedNPV;
    double TMP, tmp, lastFloatNPV=0.0, lastFixedNPV=0.0;
    double last_T4M_df = SwapT4M.Elt(SizeT4M, 0);
    double last_T4M_date = SwapT4M.Elt(SizeT4M, 1);

    double Df_Ctr;
    char* ccyName = ccy->GetCcyName();

    int missing_flows, last_missing_flows; 

    int m1=startSwapDate.GetMonth();
    int m2=itsAsOfDate.GetMonth();
    
    if (m1 == m2) 
       fwd_bwd = K_BWD;
    else
       fwd_bwd = K_FWD;


    // Si methode PAR, on calcule les taux de swap manquant 
    // par interpolation lineaire
    // des taux de swap renseignes 

    if ( raw == K_PAR )
    {
       for (i=1; i<=SizeTAM; i++)
       {
           if ( SwapTAM->Elt(i, 0) < ISPRESENT )
           {
              j = i+1;
              while ( SwapTAM->Elt(j, 0) < ISPRESENT ) 
                    j++;

              double first_swap =SwapTAM->Elt(i-1,0);
              double scd_swap = SwapTAM->Elt(j,0);

              double tmp = (first_swap + (scd_swap-first_swap)
                *DaysBetweenDates(KACTUAL_ACTUAL, 
             (ARM_Date) SwapTAM->Elt(i-1, 1),(ARM_Date) SwapTAM->Elt(i, 1))
             /DaysBetweenDates(KACTUAL_ACTUAL, 
             (ARM_Date) SwapTAM->Elt(i-1, 1),(ARM_Date) SwapTAM->Elt(j, 1)));

              SwapTAM->Elt(i, 0) = (tmp);
            }
       }
    }
    
    // Calcul des discount facteurs

    for (i = 1; i <= SizeTAM ; i++)
    {
        if ( (SwapTAM->Elt(i, 0) < ISPRESENT) && (raw == K_RAW) )
        { 
           // Cas des swaps absents dans la methode raw
           // on part du principe que les taux sont connus entre i-1 et j-1

           j = i+1;

           while ( SwapTAM->Elt(j, 0) < ISPRESENT ) 
               j++;

           missing_flows = j-i;
        }
        else
        {
           j = i;

           missing_flows = 0.0;
        }

        //  on cherche a egaliser les branches fixes et variables des swaps

        Df_Ctr = FIRST_DF_FOR_NEWTON;
            
        do
        {    // premier flux

             if ( i == 1 )
             {
                floatNPV = FirstFloatingTAMNPV ( Df_Ctr, 
                                    startSwapDate.GetJulian(), 
                                    first_fixing_df, last_T4M_df, 
                                    last_T4M_date, mean_rates, fwd_bwd, 
                                    Cont_Lin, ccy, &TMP);

                last_fixing_df = TMP; // on recupere la valeur du dernier fixing

                der_floatNPV = (FirstFloatingTAMNPV(Df_Ctr+H_FOR_DERIVATIVE, 
                                startSwapDate.GetJulian(),  
                                first_fixing_df, last_T4M_df, 
                                last_T4M_date, mean_rates, fwd_bwd, 
                                Cont_Lin, ccy, &tmp) - 
                                floatNPV )/H_FOR_DERIVATIVE;
                        

                fixedNPV = FirstFixedNPV(Df_Ctr, 
                             RefDate(SwapTAM->Elt(j, 1),ccy), last_T4M_df ,
                             last_T4M_date, startSwapDate.GetJulian(), 
                             fwd_bwd,
                             K_TAM_INDEX, Cont_Lin,  ccy);
                    
                der_fixedNPV = (FirstFixedNPV( Df_Ctr + H_FOR_DERIVATIVE, 
                                RefDate(SwapTAM->Elt(j, 1),ccy), last_T4M_df ,
                                last_T4M_date, startSwapDate.GetJulian(),
                                fwd_bwd,
                                K_TAM_INDEX, Cont_Lin,  ccy)
                                - fixedNPV)/H_FOR_DERIVATIVE;

             }
             else // autres flux
             {
                if ( (SwapTAM->Elt(j, 2) < ISPRESENT )
                    && ( SwapTAM->Elt(i-1, 2) 
                         < ISPRESENT) ) // il n'y a pas de rompu
                {
                    if ( i == 2 )
                       last_df = last_fixing_df; 
                    else 
                       last_df = Interpol( RefDate(SwapTAM->Elt(i-2
                                      -last_missing_flows, 1), ccy),  
                                      RefDate(SwapTAM->Elt(i-1, 1), ccy), 
                                      SwapTAM->Elt(i-1, 1), 
                                      SwapTAM->Elt(i-2-last_missing_flows, 0), 
                                      SwapTAM->Elt(i-1, 0), Cont_Lin,
                                      KACTUAL_ACTUAL);

                    k = i-1;

                    floatNPV = FloatingNPV(Df_Ctr, 
                                           RefDate(SwapTAM->Elt(j, 1),ccy), 
                                           SwapTAM->Elt(i-1, 0),
                                           RefDate(SwapTAM->Elt(i-1, 1),ccy), 
                                           last_df, SwapTAM->Elt(k, 1), 
                                           lastFloatNPV, missing_flows, 
                                           K_TAM_INDEX, 
                                           Cont_Lin, ccy);

                    der_floatNPV = (FloatingNPV(Df_Ctr+H_FOR_DERIVATIVE, 
                                    RefDate(SwapTAM->Elt(j, 1),ccy), 
                                    SwapTAM->Elt(i-1, 0), 
                                    RefDate(SwapTAM->Elt(i-1, 1),ccy), last_df, 
                                    SwapTAM->Elt(k, 1), lastFloatNPV,
                                    missing_flows, K_TAM_INDEX, Cont_Lin, ccy)
                                    -floatNPV)/H_FOR_DERIVATIVE;

                    fixedNPV = FixedNPV(Df_Ctr, 
                                        RefDate(SwapTAM->Elt(j, 1),ccy), 
                                        SwapTAM->Elt(i-1, 0),
                                        RefDate(SwapTAM->Elt(i-1, 1), ccy), 
                                        SwapTAM->Elt(k, 1), 
                                        lastFixedNPV, missing_flows, 
                                        K_TAM_INDEX, Cont_Lin,  ccy);
                    
                    der_fixedNPV = (FixedNPV( Df_Ctr+H_FOR_DERIVATIVE, 
                                    RefDate(SwapTAM->Elt(j, 1),ccy), 
                                    SwapTAM->Elt(i-1, 0), 
                                    RefDate( SwapTAM->Elt(i-1, 1), ccy), 
                                    SwapTAM->Elt(k, 1), lastFixedNPV, 
                                    missing_flows,  
                                    K_TAM_INDEX, Cont_Lin, ccy)- fixedNPV)
                                    /H_FOR_DERIVATIVE;

                }
                else
                {
                   // le premier taux de swap est brise
                   if ( SwapTAM->Elt(i-1, 2) > ISPRESENT )  
                   {    
                      k = i-2;

                      while (SwapTAM->Elt(k, 2) >ISPRESENT) 
                          k--; // k correspond au dernier pilier 
                               // non brise renseigne
                   }
                                

                   // le dernier taux de swap n'est pas brise
                   if  ( SwapTAM->Elt(j, 2) < ISPRESENT )
                   {
                       if (( i == 2 )||( k-1-last_missing_flows == 0 )) 
                          last_df = last_fixing_df; 
                       else
                          last_df = Interpol(RefDate(SwapTAM->Elt(k-1
                                        -last_missing_flows, 1), ccy),  
                                        RefDate(SwapTAM->Elt(k, 1), ccy), 
                                        SwapTAM->Elt(k, 1), 
                                        SwapTAM->Elt(k-1-last_missing_flows, 
                                         0), SwapTAM->Elt(k, 0), Cont_Lin,
                                        KACTUAL_ACTUAL);


                       floatNPV = FloatingNPV(Df_Ctr, 
                                              RefDate(SwapTAM->Elt(j, 1),ccy), 
                                              SwapTAM->Elt(i-1, 0),
                                              RefDate(SwapTAM->Elt(i-1, 1),
                                              ccy), last_df, SwapTAM->Elt(k, 1),
                                              lastFloatNPV, missing_flows, 
                                              K_TAM_INDEX, 
                                              Cont_Lin, ccy);

                       der_floatNPV = (FloatingNPV(Df_Ctr+H_FOR_DERIVATIVE, 
                                       RefDate(SwapTAM->Elt(j, 1),ccy), 
                                       SwapTAM->Elt(i-1, 0),
                                       RefDate(SwapTAM->Elt(i-1, 1),ccy), 
                                       last_df, 
                                       SwapTAM->Elt(k, 1), lastFloatNPV,
                                       missing_flows, K_TAM_INDEX, 
                                       Cont_Lin, ccy)-floatNPV)
                                       /H_FOR_DERIVATIVE;

                       fixedNPV = FixedNPV(Df_Ctr, RefDate(SwapTAM->Elt(j, 1),
                                       ccy), SwapTAM->Elt(i-1, 0),
                                   RefDate( SwapTAM->Elt(i-1, 1) , ccy), 
                                   SwapTAM->Elt(k, 1), 
                                   lastFixedNPV, missing_flows, K_TAM_INDEX, 
                                   Cont_Lin, ccy);
                    
                       der_fixedNPV = (FixedNPV(Df_Ctr+H_FOR_DERIVATIVE, 
                                        RefDate(SwapTAM->Elt(j, 1),ccy), 
                                        SwapTAM->Elt(i-1, 0), 
                                        RefDate( SwapTAM->Elt(i-1, 1) , ccy), 
                                        SwapTAM->Elt(k, 1), lastFixedNPV, 
                                        missing_flows,  
                                        K_TAM_INDEX, Cont_Lin, ccy)
                                        -fixedNPV)/ H_FOR_DERIVATIVE;
                   }
                   else // (SwapTAM->Elt(i-1, 2) >ISPRESENT)  
                   {    // le premier taux de swap est bris‹

                      floatNPV = FloatingTAMNPVWithRompu(Df_Ctr, 
                                       startSwapDate.GetJulian(), 
                                       first_fixing_df,
                                       SwapTAM->Elt(1, 0), 
                                       SwapTAM->Elt(1, 1), SwapT4M, 
                                       SizeT4M, mean_rates, 
                                       SwapTAM->Elt(j, 2), fwd_bwd, 
                                       Cont_Lin, ccy);

                      der_floatNPV = (FloatingTAMNPVWithRompu(Df_Ctr
                               +H_FOR_DERIVATIVE, startSwapDate.GetJulian(), 
                                     first_fixing_df,
                                     SwapTAM->Elt(1, 0), SwapTAM->Elt(1, 1), 
                                     SwapT4M, SizeT4M, mean_rates, 
                                     SwapTAM->Elt(j, 2), 
                                     fwd_bwd, Cont_Lin, ccy)
                                     - floatNPV)/H_FOR_DERIVATIVE;
 

                      fixedNPV = FixedTAMNPVWithRompu(Df_Ctr, 
                                   startSwapDate.GetJulian(), 
                                   SwapTAM->Elt(1, 0),
                                   SwapTAM->Elt(1, 1), SwapTAM->Elt(j, 0), 
                                   SwapT4M, SizeT4M, SwapTAM->Elt(j, 2),
                                   Cont_Lin, ccy);

                    
                      der_fixedNPV = (FixedTAMNPVWithRompu(Df_Ctr
                              +H_FOR_DERIVATIVE, startSwapDate.GetJulian(), 
                              SwapTAM->Elt(1, 0), SwapTAM->Elt(1, 1), 
                              SwapTAM->Elt(j, 0), SwapT4M, 
                              SizeT4M, SwapTAM->Elt(j, 2) , 
                              Cont_Lin,  ccy)- fixedNPV)    
                              /H_FOR_DERIVATIVE;
                            }
                }
             }

             if ( SwapTAM->Elt(j, 2) > ISPRESENT )
             {
                newton_result= ( fixedNPV - floatNPV)/
                                    (der_fixedNPV - der_floatNPV);
             }
             else
             {
                newton_result = ( SwapTAM->Elt(j, 0)/100.0*fixedNPV 
                                  - floatNPV)/
                                  ( SwapTAM->Elt(j, 0)/100.0*der_fixedNPV 
                                  - der_floatNPV);
             }

             Df_Ctr = Df_Ctr - newton_result;
        }
        while ( fabs(newton_result) > PRECISION );
                
        if ( SwapTAM->Elt(j, 2) < ISPRESENT )
        {
           lastFloatNPV = floatNPV;
           lastFixedNPV = fixedNPV;
        }

        i = i + missing_flows;
        SwapTAM->Elt(i, 0) = Df_Ctr;
        last_missing_flows=missing_flows;
    }
                    
    // si meth = RAW, virer les taux manquant initialement 
    if ( raw == K_RAW )
    {  
       for (i=1; i<=SizeTAM; i++)
       {
           if (SwapTAM->Elt(i, 0)<ISPRESENT) // ils n'ont pas ete calcules
           {
              for (j=i; j<SizeTAM; j++)
              {
                  SwapTAM->Elt(j, 0) = SwapTAM->Elt(j+1, 0);
                  SwapTAM->Elt(j, 1) = SwapTAM->Elt(j+1, 1);
              }

              SizeTAM--;

              i--;
           }
        }
    }

    // on recalcule les dates de mat des piliers (premier jour ouvre du mois)

    for (i=1; i<=SizeTAM; i++)
        SwapTAM->Elt(i, 1) = RefDate(SwapTAM->Elt(i, 1), ccy);


}




/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
