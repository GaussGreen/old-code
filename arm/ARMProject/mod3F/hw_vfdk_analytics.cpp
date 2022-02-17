/*===========================================================
  Name    : hw_vfdk_analytics.cpp
  Owner   : DK
  Created : DK
  Comment : Don't touch if your name is not DK !
  Time    : Last update 13:36 06/28/2004
=============================================================*/
#include <iostream>
#include <strstream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifndef WIN32
#include <sys/param.h>
#endif 
#define ACC 1.0e-04
#define JMAX 5000
#define MAX(x,y)    (x >= y ? x : y)
#define g_dPi 3.141592653589790 

unsigned int g_SP=1;

#include "Swaption.h"
#include "dk_utils.h"
#include "DKMaille.h"
#include "DKMaille2D.h"
#include "hw_vfdk_analytics.h"
#include "hw_vfdk_LDHD_lattice.h"


double Analytical_VFDK_HW1F(double dSpotDate,
                            DKMaille<double> dDiscountDates,
                            DKMaille<double> dDiscountRates,
                            DKMaille<double> dSwaptionDates,
                            double dNoticePeriod,
                            double dMeanReversion,
                            double dBlackScholesVol,
                            double dOutput)
{
    double dBlackScholesPrice =0.;

    DKMaille<double> dDiscountFactors(dSwaptionDates.entries());
    DKMaille<double> dPeriods(dSwaptionDates.entries());
    DKMaille<double> dSigma(dSwaptionDates.entries());
    DKMaille<double> dZCRates(dDiscountRates.entries());

    dPeriods[0]=0.;
    double dAnnuity = 0.;
    double dSwapRate = 0.;

    unsigned int ui = 0;
    for(ui=0;ui<dDiscountDates.entries();ui++)
        dDiscountDates[ui]=(dDiscountDates[ui]-dSpotDate)/365.;

    for(ui=0;ui<dSwaptionDates.entries();ui++)
        dSwaptionDates[ui]=(dSwaptionDates[ui]-dSpotDate)/365.;

    for(ui=0;ui<dDiscountDates.entries();ui++)
    {
        if(dDiscountDates.at(ui)>0.) dZCRates.at(ui)=-log(dDiscountRates.at(ui))/dDiscountDates.at(ui);
        else dZCRates.at(ui)=-log(dDiscountRates.at(ui+1))/dDiscountDates.at(ui+1);
    }

    for(ui=0;ui<dSwaptionDates.entries();ui++)
    {
        dDiscountFactors[ui]=exp(-rateinterpolation_dk_maille(2,dSwaptionDates[ui],dDiscountDates,
                                 dZCRates,dDiscountDates.entries()-1)*dSwaptionDates[ui]);
        if(ui>0)
        {
            dPeriods[ui]=dSwaptionDates[ui]-dSwaptionDates[ui-1];
            dAnnuity+=dPeriods[ui]*dDiscountFactors[ui];
        }
    }
    dSwapRate = (dDiscountFactors[0]-dDiscountFactors[dSwaptionDates.entries()-1])/dAnnuity;

    if(dOutput==0.) return dAnnuity;
    else if(dOutput==1.)	return dSwapRate;
    else if(dOutput==2.)
    {
        double dObservationDate = 0.;
        double dTerminalDate    = dSwaptionDates[0]-dNoticePeriod;
        double dPi = 3.141592653589790 ;
        DKMaille<double> dVector(dSwaptionDates.entries());

        for(ui=0;ui<dSwaptionDates.entries();ui++)
        {
            dVector[ui]=0.;
            dSigma[ui]=1.;
        }

        dVector=dKacVector(dSwaptionDates.entries(),dSwapRate,dPeriods,dDiscountFactors);
        double dForwardStdDev=0.;
        double dIntermediateCalc=0.;

        for(unsigned int iK=0;iK<dSwaptionDates.entries();iK++) {
            for(unsigned int iM=0;iM<dSwaptionDates.entries();iM++) {
                dIntermediateCalc=dVector[iK]*dVector[iM]*
                                  dKacDeterminant(dMeanReversion,dObservationDate,dTerminalDate,dSigma[iK],
                                                  dSigma[iM],dSwaptionDates[iK],dSwaptionDates[iM]);
                dForwardStdDev += dIntermediateCalc;
            }
        }
        dBlackScholesPrice = dAnnuity*dSwapRate*ATMpricevol_DK(dTerminalDate,dBlackScholesVol)*10000.;
        double dxx=0.;
        int ii=0;
        double xprice = 0.;
        double xprice_derivative = 0.;
        double dVol=dBlackScholesVol*dSwapRate;
        while(ii<=100)
        {
            xprice = 10000. * (dVol/fabs(dMeanReversion)) * sqrt(dForwardStdDev) / sqrt ( 2. * dPi ) - dBlackScholesPrice;
            xprice_derivative = 10000. * (1./fabs(dMeanReversion)) * sqrt(dForwardStdDev) / sqrt ( 2. * dPi );
            dxx=xprice/xprice_derivative;
            dVol-=dxx;
            if(fabs(dxx)<1.e-05) ii=100;
            ii+=1;
        }
        return dVol;
    }
    else
    {
        return 0.;
    }
}// Analytical_VFDK_HW1F


double ATMnormal(double dSigma,double dTimeToExpiry)
{

    return dSigma*sqrt(dTimeToExpiry/(2.*3.141592653589790));
}

double ATMlognormal(double dSigma,double dRate,double dTimeToExpiry)
{
    double dArgument=(dSigma)*sqrt(dTimeToExpiry)/2.;
    double N1=0.; double N2=0.;
    normalDK(&N1,dArgument);
    normalDK(&N2,-dArgument);
    return (dRate)*(N1-N2);
}

double ATMq(double dSigma,double dRate,double dTimeToExpiry,double dQ)
{
    double dArgument=dQ*(dSigma/dRate)*sqrt(dTimeToExpiry)/2.;
    double N1=0.; double N2=0.;
    normalDK(&N1,dArgument);
    normalDK(&N2,-dArgument);
    if(dQ<0.) throw("Q negative is not a option");
    if(dQ>0.) return (dRate/dQ)*(N1-N2);
    else return dSigma*sqrt(dTimeToExpiry/(2.*3.141592653589790));
}


double Analytical_QModel(double dSpotDate,
                         DKMaille<double> dDiscountDates,
                         DKMaille<double> dDiscountRates,
                         DKMaille<double> dSwaptionDates,
                         double dNoticePeriod,
                         double dStrike,
                         double dLognormalNormalSmileModelType,
                         double dQVol,
                         double dQParameter,
                         double dCallPut)
{
    DKMaille<double> dDiscountFactors(dSwaptionDates.entries());
    DKMaille<double> dPeriods(dSwaptionDates.entries());
    DKMaille<double> dZCRates(dDiscountRates.entries());
    dPeriods[0]=0.;
    double dAnnuity = 0.;
    double dSwapRate = 0.;
    unsigned int ui = 0;
    for(ui=0;ui<dDiscountDates.entries();ui++)
        dDiscountDates[ui]=(dDiscountDates[ui]-dSpotDate)/365.;
    for(ui=0;ui<dSwaptionDates.entries();ui++)
        dSwaptionDates[ui]=(dSwaptionDates[ui]-dSpotDate)/365.;
    for(ui=0;ui<dDiscountDates.entries();ui++)
    {
        if(dDiscountDates.at(ui)>0.) dZCRates.at(ui)=-log(dDiscountRates.at(ui))/dDiscountDates.at(ui);
        else dZCRates.at(ui)=-log(dDiscountRates.at(ui+1))/dDiscountDates.at(ui+1);
    }
    for(ui=0;ui<dSwaptionDates.entries();ui++)
    {
        dDiscountFactors[ui]=exp(-rateinterpolation_dk_maille(2,dSwaptionDates[ui],dDiscountDates,
                                 dZCRates,dDiscountDates.entries()-1)*dSwaptionDates[ui]);
        if(ui>0)
        {
            dPeriods[ui]=dSwaptionDates[ui]-dSwaptionDates[ui-1];
            dAnnuity+=dPeriods[ui]*dDiscountFactors[ui];
        }
    }
    dSwapRate = (dDiscountFactors[0]-dDiscountFactors[dSwaptionDates.entries()-1])/dAnnuity;
    double dTerminalDate=dSwaptionDates[0]-dNoticePeriod;
    double dPi=3.141592653589790 ;
    double x0=0.;
    double argument1=0.;
    double argument2=0.;
    double variance=0.;

    variance=(dQVol)*sqrt(dTerminalDate);


    if(dLognormalNormalSmileModelType==2.)
    {
        // Calibration of Q-Model to normal model for ATM option
        double dGuess=dQVol;
        double dPrint=dQVol;
        double dRatio=0.;
        int iCounter=0;
        double dPrice=ATMnormal(dQVol,dTerminalDate);
        while(iCounter<100)
        {
            iCounter++;
            double dFunction=ATMq(dGuess,dSwapRate,dTerminalDate,dQParameter)-dPrice;
            double dFunctionUp=ATMq(dGuess+0.00001,dSwapRate,dTerminalDate,dQParameter)-dPrice;
            dRatio=0.00001*dFunction/(dFunctionUp-dFunction);
            dGuess-=dRatio;
            if(fabs(dRatio)<.5e-5)
            {
                iCounter=100;
            }
        } // while iCounter < 100
        dQVol=dGuess;
        if(dGuess<0.||dGuess>1.) dQVol=dPrint;
        variance=(dQVol/dSwapRate)*sqrt(dTerminalDate);
    }

    if(dLognormalNormalSmileModelType==0.)
    {
        // Calibration of Q-Model to normal model for ATM option
        double dGuess=dQVol/dSwapRate;
        double dRatio=0.;
        int iCounter=0;
        double dPrice=ATMnormal(dQVol,dTerminalDate);
        while(iCounter<100)
        {
            iCounter++;
            double dFunction=ATMlognormal(dGuess,dSwapRate,dTerminalDate)-dPrice;
            double dFunctionUp=ATMlognormal(dGuess+0.00001,dSwapRate,dTerminalDate)-dPrice;
            dRatio=0.00001*dFunction/(dFunctionUp-dFunction);
            dGuess-=dRatio;
            if(fabs(dRatio)<.5e-5)
            {
                iCounter=100;
            }
        } // while iCounter < 100
        dQVol=dGuess;
        variance=(dQVol)*sqrt(dTerminalDate);
    }

    double argument3=0.;
    double argument4=0.;
    double dPrice=0.;
    // dSwapRate=0.02;
    if(dLognormalNormalSmileModelType==2.)
    {
        if(dQParameter>=1.e-9)
        {
            x0=(1./dQParameter)*log(1.+dQParameter*(-1.+(dStrike/dSwapRate)));
            argument1=(0.5*dQParameter*variance*variance-x0)/variance;
            argument2=(-0.5*dQParameter*variance*variance-x0)/variance;
            normalDK(&argument3,argument1);
            normalDK(&argument4,argument2);
            dPrice=dAnnuity*((dSwapRate/dQParameter)*argument3+(dSwapRate-dStrike-(dSwapRate/dQParameter))*argument4);
        }
        else
        {
            x0=(dStrike/dSwapRate-1.);
            argument1=-x0/variance;
            normalDK(&argument3,argument1);
            argument4=exp(-0.5*(argument1)*(argument1))*(1./sqrt(2.*dPi));
            dPrice=dAnnuity*((dSwapRate-dStrike)*argument3+dSwapRate*variance*argument4);
        }
    }
    if(dLognormalNormalSmileModelType==0.)
    {
        argument1=(0.5*variance*variance+log(dSwapRate/dStrike))/variance;
        argument2=(-0.5*variance*variance+log(dSwapRate/dStrike))/variance;
        normalDK(&argument3,argument1);
        normalDK(&argument4,argument2);
        dPrice=dAnnuity*((dSwapRate)*argument3-dStrike*argument4);
    }
    if(dLognormalNormalSmileModelType==1.)
    {
        dPrice=GetOptionPrice(0.,dSwapRate,dAnnuity,dStrike,dQVol,dTerminalDate)/10000.;
    }
    if(dCallPut==1.) dPrice=dPrice-dAnnuity*(dSwapRate-dStrike);
    return dPrice*10000.;
}// Analytical_QModel

double Delta_QVol(double dSpotDate,
                  DKMaille<double> dDiscountDates,
                  DKMaille<double> dDiscountRates,
                  DKMaille<double> dSwaptionDates,
                  double dNoticePeriod,
                  double dLognormalNormalSmileModelType,
                  double dQVol,
                  double dQParameter,
                  double dShift)
{
    DKMaille<double> dPeriods(dSwaptionDates.entries());
    dPeriods.at(0)=0.;
    for(unsigned int ui=0;ui<dSwaptionDates.entries();ui++)
    {
        if(ui>0)
        {
            dPeriods.at(ui)=(dSwaptionDates[ui]-dSwaptionDates[ui-1])/365.;
        }
    }
    DKMaille<double> dShiftedDiscountRates(dDiscountRates.entries());
    for(ui=0;ui<dDiscountRates.entries();ui++)
    {
        dShiftedDiscountRates.at(ui)=dDiscountRates.at(ui)
                                     -dShift*(dDiscountRates.at(ui))*(dDiscountDates.at(ui)-dSpotDate)/365.;
    }

    double dATMSwapRate=GetSwapRate(dSpotDate,dSwaptionDates,dPeriods,dDiscountDates,dDiscountRates,dDiscountRates);
    double dBase=Analytical_QModel(dSpotDate,dDiscountDates,dShiftedDiscountRates,dSwaptionDates,dNoticePeriod,dATMSwapRate,1.,dQVol,dQParameter,0.)+
                 Analytical_QModel(dSpotDate,dDiscountDates,dShiftedDiscountRates,dSwaptionDates,dNoticePeriod,dATMSwapRate,1.,dQVol,dQParameter,1.)-
                 Analytical_QModel(dSpotDate,dDiscountDates,dDiscountRates,dSwaptionDates,dNoticePeriod,dATMSwapRate,1.,dQVol,dQParameter,0.)-
                 Analytical_QModel(dSpotDate,dDiscountDates,dDiscountRates,dSwaptionDates,dNoticePeriod,dATMSwapRate,1.,dQVol,dQParameter,1.);

    // Calibration of Q-Model to normal model for ATM option
    double dGuess=0.;
    double dRatio=0.;
    int iCounter=0;
    while(iCounter<100)
    {
        iCounter++;
        double dFunction=Analytical_QModel(dSpotDate,dDiscountDates,dShiftedDiscountRates,dSwaptionDates,dNoticePeriod,dATMSwapRate,dLognormalNormalSmileModelType-3.,dQVol+dGuess,dQParameter,0.)+
                         Analytical_QModel(dSpotDate,dDiscountDates,dShiftedDiscountRates,dSwaptionDates,dNoticePeriod,dATMSwapRate,dLognormalNormalSmileModelType-3.,dQVol+dGuess,dQParameter,1.)-
                         Analytical_QModel(dSpotDate,dDiscountDates,dDiscountRates,dSwaptionDates,dNoticePeriod,dATMSwapRate,dLognormalNormalSmileModelType-3.,dQVol,dQParameter,0.)-
                         Analytical_QModel(dSpotDate,dDiscountDates,dDiscountRates,dSwaptionDates,dNoticePeriod,dATMSwapRate,dLognormalNormalSmileModelType-3.,dQVol,dQParameter,1.)
                         -dBase;
        double dFunctionUp=Analytical_QModel(dSpotDate,dDiscountDates,dShiftedDiscountRates,dSwaptionDates,dNoticePeriod,dATMSwapRate,dLognormalNormalSmileModelType-3.,dQVol+dGuess+0.00001,dQParameter,0.)+
                           Analytical_QModel(dSpotDate,dDiscountDates,dShiftedDiscountRates,dSwaptionDates,dNoticePeriod,dATMSwapRate,dLognormalNormalSmileModelType-3.,dQVol+dGuess+0.00001,dQParameter,1.)-
                           Analytical_QModel(dSpotDate,dDiscountDates,dDiscountRates,dSwaptionDates,dNoticePeriod,dATMSwapRate,dLognormalNormalSmileModelType-3.,dQVol,dQParameter,0.)-
                           Analytical_QModel(dSpotDate,dDiscountDates,dDiscountRates,dSwaptionDates,dNoticePeriod,dATMSwapRate,dLognormalNormalSmileModelType-3.,dQVol,dQParameter,1.)
                           -dBase;
        dRatio=0.00001*dFunction/(dFunctionUp-dFunction);
        dGuess-=dRatio;
        if(fabs(dRatio)<.5e-5)
        {
            iCounter=100;
        }
    } // while iCounter < 100
    return dGuess;
}// Delta_QVol

double GetAdjAnnuity(double dSpotDate,
                     const DKMaille<double> &dSwaptionDates,
                     const DKMaille<double> &dAccrualPeriods,
                     const DKMaille<double> &dDiscountDates,
                     const DKMaille<double> &dDiscountRates,
                     const DKMaille<double> &dAdjustedDiscountRates,
                     double dJulianObservationDate)
{
    DKMaille<double> dDiscountFactors(dSwaptionDates.entries());
    DKMaille<double> dPeriods(dAccrualPeriods.entries());

    dPeriods[0]=dAccrualPeriods[0];

    DKMaille<double> dZCRates(dDiscountDates.entries());
    DKMaille<double> dZCAdjRates(dDiscountDates.entries());
    FromDiscountToZero(dSpotDate,dZCRates,dZCAdjRates,dDiscountDates,dDiscountRates,dAdjustedDiscountRates);
    double dDiscountObservation=0.;
    if(dDiscountObservation!=dSpotDate) dDiscountObservation=exp(-rateinterpolation_dk_maille(2,dJulianObservationDate,dDiscountDates,
                dZCAdjRates,dDiscountDates.entries()-1)*(dJulianObservationDate-dSpotDate)/365.);
    else dDiscountObservation=1.;

    double dAnnuity = 0.;

    unsigned ui=0;
    for(;ui<dSwaptionDates.entries();ui++)
    {
        dDiscountFactors[ui]=exp(-rateinterpolation_dk_maille(2,dSwaptionDates[ui],dDiscountDates,
                                 dZCAdjRates,dDiscountDates.entries()-1)*(dSwaptionDates[ui]-dSpotDate)/365.)/dDiscountObservation;
        if(ui>0)
        {
            dPeriods[ui]=dAccrualPeriods[ui];
            dAnnuity+=dPeriods[ui]*dDiscountFactors[ui];
        }
    }

    return dAnnuity;
} // GetAdjAnnuity(...)


double BasisConversion(double dSpotDate,
                       DKMaille<double> &dFundingDiscountCurveDates,
                       DKMaille<double> &dFundingDiscountCurve,
                       DKMaille<double> &dFundingDiscountCurveAdjusted,
                       DKMaille<double> &dDomesticDiscountCurveDates,
                       DKMaille<double> &dDomesticDiscountCurve,
                       DKMaille<double> &dDomesticDiscountCurveAdjusted,
                       DKMaille<double> &dFundingStartDates,
                       DKMaille<double> &dFundingEndDates,
                       DKMaille<double> &dDomesticStartDates,
                       DKMaille<double> &dDomesticEndDates,
                       double dActualisationFunding,
                       double dActualisationDomestic,
                       double dFundingNotional,
                       double dDomesticNotional,
                       double dFundingSpread,
                       double dFXSpotRate)
{

    DKMaille<double> dFundingZCRates(dFundingDiscountCurveDates.entries());
    DKMaille<double> dFundingZCAdjRates(dFundingDiscountCurveDates.entries());
    DKMaille<double> dDomesticZCRates(dDomesticDiscountCurveDates.entries());
    DKMaille<double> dDomesticZCAdjRates(dDomesticDiscountCurveDates.entries());
    FromDiscountToZero(dSpotDate,dFundingZCRates,dFundingZCAdjRates,dFundingDiscountCurveDates,dFundingDiscountCurve,dFundingDiscountCurveAdjusted);
    FromDiscountToZero(dSpotDate,dDomesticZCRates,dDomesticZCAdjRates,dDomesticDiscountCurveDates,dDomesticDiscountCurve,dDomesticDiscountCurveAdjusted);

    // Calc Funding Currency PV
    double dFundingPV=0.;
    for(unsigned int uiM=0;uiM<dFundingStartDates.entries();uiM++)
    {
        // Start Date ZC
        double dZCDiscountStart=exp(-rateinterpolation_dk_maille(2,
                                    dFundingStartDates.at(uiM),
                                    dFundingDiscountCurveDates,
                                    dFundingZCAdjRates,
                                    dFundingDiscountCurveDates.entries()-1)
                                    *(dFundingStartDates.at(uiM)-dSpotDate)/365.);

        double dZCDiscountEnd=exp(-rateinterpolation_dk_maille(2,
                                  dFundingEndDates.at(uiM),
                                  dFundingDiscountCurveDates,
                                  dFundingZCAdjRates,
                                  dFundingDiscountCurveDates.entries()-1)
                                  *(dFundingEndDates.at(uiM)-dSpotDate)/365.);

        double dZCDiscountStartCash=exp(-rateinterpolation_dk_maille(2,
                                        dFundingStartDates.at(uiM),
                                        dFundingDiscountCurveDates,
                                        dFundingZCRates,
                                        dFundingDiscountCurveDates.entries()-1)
                                        *(dFundingStartDates.at(uiM)-dSpotDate)/365.);

        double dZCDiscountEndCash=exp(-rateinterpolation_dk_maille(2,
                                      dFundingEndDates.at(uiM),
                                      dFundingDiscountCurveDates,
                                      dFundingZCRates,
                                      dFundingDiscountCurveDates.entries()-1)
                                      *(dFundingEndDates.at(uiM)-dSpotDate)/365.);


        double dPeriod=(dFundingEndDates.at(uiM)-dFundingStartDates.at(uiM))/dActualisationFunding;

        double dFundingLIBOR=(dZCDiscountStart/dZCDiscountEnd-1.)/dPeriod;

        double dBasisAdjustment=((dZCDiscountStartCash/dZCDiscountEndCash)/(dZCDiscountStart/dZCDiscountEnd)-1.)/dPeriod;


        dFundingPV+=dFundingNotional*(/*dFundingLIBOR*/+dFundingSpread+dBasisAdjustment)*dZCDiscountEnd*dPeriod;

        /*
        if(uiM==dFundingStartDates.entries()-1)
        	dFundingPV+=dFundingNotional*dZCDiscountEnd;

        if(uiM==0) dFundingPV-=dFundingNotional*dZCDiscountStart;
        */
    }

    // Calc Domestic Currency PV
    double dDomesticPV=0.;
    double dDomesticDV01=0.;
    for(uiM=0;uiM<dDomesticStartDates.entries();uiM++)
    {
        // Start Date ZC
        double dZCDiscountStart=exp(-rateinterpolation_dk_maille(2,
                                    dDomesticStartDates.at(uiM),
                                    dDomesticDiscountCurveDates,
                                    dDomesticZCAdjRates,
                                    dDomesticDiscountCurveDates.entries()-1)
                                    *(dDomesticStartDates.at(uiM)-dSpotDate)/365.);

        double dZCDiscountEnd=exp(-rateinterpolation_dk_maille(2,
                                  dDomesticEndDates.at(uiM),
                                  dDomesticDiscountCurveDates,
                                  dDomesticZCAdjRates,
                                  dDomesticDiscountCurveDates.entries()-1)
                                  *(dDomesticEndDates.at(uiM)-dSpotDate)/365.);

        double dZCDiscountStartCash=exp(-rateinterpolation_dk_maille(2,
                                        dDomesticStartDates.at(uiM),
                                        dDomesticDiscountCurveDates,
                                        dDomesticZCRates,
                                        dDomesticDiscountCurveDates.entries()-1)
                                        *(dDomesticStartDates.at(uiM)-dSpotDate)/365.);

        double dZCDiscountEndCash=exp(-rateinterpolation_dk_maille(2,
                                      dDomesticEndDates.at(uiM),
                                      dDomesticDiscountCurveDates,
                                      dDomesticZCRates,
                                      dDomesticDiscountCurveDates.entries()-1)
                                      *(dDomesticEndDates.at(uiM)-dSpotDate)/365.);


        double dPeriod=(dDomesticEndDates.at(uiM)-dDomesticStartDates.at(uiM))/dActualisationDomestic;

        double dDomesticLIBOR=(dZCDiscountStart/dZCDiscountEnd-1.)/dPeriod;

        double dBasisAdjustment=((dZCDiscountStartCash/dZCDiscountEndCash)/(dZCDiscountStart/dZCDiscountEnd)-1.)/dPeriod;

        dDomesticDV01+=dDomesticNotional*dZCDiscountEnd*dPeriod/dFXSpotRate;

        dDomesticPV+=dDomesticNotional*(/*dDomesticLIBOR*/+dBasisAdjustment)*dZCDiscountEnd*dPeriod/dFXSpotRate;

        /*
        if(uiM==dDomesticStartDates.entries()-1)
        	dDomesticPV+=dDomesticNotional*dZCDiscountEnd/dFXSpotRate;

        if(uiM==0)
        	dDomesticPV-=dZCDiscountStart*dDomesticNotional/dFXSpotRate;
        */
    }





    return (dFundingPV-dDomesticPV)/dDomesticDV01;

} // BasisConversion




double GetSwapRate(double dSpotDate,
                   const DKMaille<double> &dSwaptionDates,
                   const DKMaille<double> &dAccrualPeriods,
                   const DKMaille<double> &dDiscountDates,
                   const DKMaille<double> &dDiscountRates,
                   const DKMaille<double> &dAdjustedDiscountRates)
{
    DKMaille<double> dDiscountFactors(dSwaptionDates.entries());
    DKMaille<double> dPeriods(dAccrualPeriods.entries());
    DKMaille<double> dZCRates(dDiscountDates.entries());
    DKMaille<double> dZCAdjRates(dDiscountDates.entries());

    FromDiscountToZero(dSpotDate,dZCRates,dZCAdjRates,dDiscountDates,dDiscountRates,dAdjustedDiscountRates);
    dPeriods[0]=dAccrualPeriods[0];

    double dAnnuity = 0.;
    double dSwapRate = 0.;

    unsigned int ui=0;
    for(;ui<dSwaptionDates.entries();ui++)
    {
        dDiscountFactors[ui]=exp(-rateinterpolation_dk_maille(2,dSwaptionDates[ui],dDiscountDates,
                                 dZCRates,dDiscountDates.entries()-1)*(dSwaptionDates[ui]-dSpotDate)/365.);
        if(ui>0)
        {
            dPeriods[ui]=dAccrualPeriods[ui];
            dAnnuity+=dPeriods[ui]*dDiscountFactors[ui];
        }
    }
    dSwapRate = (dDiscountFactors[0]-dDiscountFactors[dSwaptionDates.entries()-1])/dAnnuity;

    return dSwapRate;
} // dGetSwapRate

DKMaille<double> GetItoVector(unsigned int uiNumberOfSwaptionDates,
                              unsigned int uiNumberOfBasisSwaptionDates,
                              const DKMaille<double> &dCash,
                              const DKMaille<double> &dDiscounts)
{
    DKMaille<double> dVector(uiNumberOfSwaptionDates+uiNumberOfBasisSwaptionDates);
    unsigned int uiIndex;
    unsigned int uiSize=uiNumberOfSwaptionDates+uiNumberOfBasisSwaptionDates;
    // fill the first element
    dVector[0]=dDiscounts[0];
    // fill the last fixed pay element
    dVector[uiNumberOfSwaptionDates-1]=-(1.+dCash[uiNumberOfSwaptionDates-1])*dDiscounts[uiNumberOfSwaptionDates-1];
    // fill the last float p[ay element
    dVector[uiSize-1]=dCash[uiSize-1]*dDiscounts[uiSize-1];
    if(uiNumberOfSwaptionDates>2)
    {
        for(uiIndex=1; uiIndex<uiNumberOfSwaptionDates-1; uiIndex++)
        {
            dVector[uiIndex]=-dCash[uiIndex]*dDiscounts[uiIndex];
        }
    }
    if(uiNumberOfBasisSwaptionDates>1)
    {
        for(uiIndex=uiNumberOfSwaptionDates; uiIndex<uiSize-1; uiIndex++)
        {
            dVector[uiIndex]=dCash[uiIndex]*dDiscounts[uiIndex];
        }
    }

    return dVector;
} // GetItoVector

DKMaille<double> GetItoVectorDK(unsigned int uiNumberOfSwaptionDates,
                                const DKMaille<double> &dCash,
                                const DKMaille<double> &dDiscounts)
{
    DKMaille<double> dVector(uiNumberOfSwaptionDates);
    unsigned int uiIndex;
    unsigned int uiSize=uiNumberOfSwaptionDates;
    // fill the first element
    dVector.at(0)=dDiscounts[0];
    // fill the last element
    dVector.at(uiNumberOfSwaptionDates-1)=-(1.+dCash[uiNumberOfSwaptionDates-1])*dDiscounts[uiNumberOfSwaptionDates-1];
    if(uiNumberOfSwaptionDates>2)
    {
        for(uiIndex=1; uiIndex<uiNumberOfSwaptionDates-1; uiIndex++)
        {
            dVector.at(uiIndex)=-dCash[uiIndex]*dDiscounts[uiIndex];
        }
    }
    return dVector;
} // GetItoVectorDK




double GetItoDeterminant(double dMeanFieldDecay, double dObservationDate, double dTerminalDate,
                         double dSigmaDate1, double dSigmaDate2, double dDate1, double dDate2)
{
    double dTerm1 = (dTerminalDate-dObservationDate);

    double dTerm2 = -(1. /dMeanFieldDecay)*(exp(-dMeanFieldDecay*(dDate1-dTerminalDate))
                                            -exp(-dMeanFieldDecay*(dDate1-dObservationDate)));

    double dTerm3 = -(1. /dMeanFieldDecay)*(exp(-dMeanFieldDecay*(dDate2-dTerminalDate))
                                            -exp(-dMeanFieldDecay*(dDate2-dObservationDate)));

    double dTerm4 = +(0.5/dMeanFieldDecay)*(exp(-dMeanFieldDecay*(dDate1+dDate2-2.*dTerminalDate))
                                            -exp(-dMeanFieldDecay*(dDate1+dDate2-2.*dObservationDate)));

    return (dTerm1+dTerm2+dTerm3+dTerm4)*dSigmaDate1*dSigmaDate2;
} // GetItoDeterminant


double GetItoDeterminant_HW1F_TV(double dMeanFieldDecay,
                                 double dObservationDate,
                                 double dTerminalDate,
                                 const DKMaille<double> &dDatesStrip,
                                 const DKMaille<double> &dSigmaStrip,
                                 double dDate1,
                                 double dDate2)

{
    int i=0;
    double dSum=0.;
    while(dDatesStrip[i]<dTerminalDate-0.00001)
    {
        dSum+=GetItoDeterminant(dMeanFieldDecay,dDatesStrip[i],dDatesStrip[i+1],
                                dSigmaStrip[i],dSigmaStrip[i],dDate1,dDate2);
        i++;
    }
    return dSum;
} // GetItoDeterminant_HW1F_TV(...)


void datesForSlices(DKMaille<double> inputDates,
                    int howManySteps, // initial number of steps
                    DKMaille<double> &timeArray)
{
    unsigned int howManyDates=inputDates.entries();
    double startdate;
    double stepsize;
    double newstepsize;
    int stepgap;
    int counter;
    register int j;
    register int i;
    timeArray.insert(inputDates.at(0));
    counter=0;
    startdate=inputDates.at(0);
    stepsize=(inputDates.at(howManyDates-1)-startdate)/(howManySteps);

    for(i=1;i<howManyDates;i++)
    {
        // Check dates are not equal
        if(inputDates.at(i)!=inputDates.at(i-1))
        {
            if((inputDates.at(i)-inputDates.at(i-1))/stepsize<1.) {
                newstepsize=(inputDates.at(i)-inputDates.at(i-1));
                stepgap=1;
            }
            else
            {
                stepgap=(int)(1.*(inputDates.at(i)-inputDates.at(i-1))/stepsize);
                if((double)stepgap<(inputDates.at(i)-inputDates.at(i-1))/stepsize-0.5) stepgap+=1;
                newstepsize=(inputDates.at(i)-inputDates.at(i-1))/(double)stepgap;
            }
            for (j=1;j<=stepgap;j++)
            {
                counter+=1;
                timeArray.insert(timeArray.at(counter-1)+newstepsize/365.);
            }
        }
    }
} // datesForSlices(...)





void Slicing(double dIntegrationStart,
             double dIntegrationEnd,
             double dBondDate1,
             double dBondDate2,
             double dt,
             DKMaille<double> &dSlices)
{

    if(dBondDate1<dIntegrationEnd) throw("Unauthorised input in GetItoDeterminant_DK1F::Slicing");
    if(dBondDate2<dIntegrationEnd) throw("Unauthorised input in GetItoDeterminant_DK1F::Slicing");
    if(dIntegrationEnd<dIntegrationStart+1.e-09) throw("Unauthorised input in GetItoDeterminant_DK1F::Slicing");


    double dFirstBondDate=0.; double dSecondBondDate=0.;
    if(dBondDate1>dBondDate2)
    {
        dFirstBondDate=dBondDate2;
        dSecondBondDate=dBondDate1;
    }
    else if(dBondDate2>dBondDate1)
    {
        dFirstBondDate=dBondDate1;
        dSecondBondDate=dBondDate2;
    }
    else
    {
        dFirstBondDate=dBondDate1;
        dSecondBondDate=dBondDate1;
    }

    double dSliceLine=dIntegrationStart;
    dSlices.insert(dSliceLine);
    double dSliceEnd=dSecondBondDate;
    unsigned int uiOut=0; unsigned int uiFilled1=0; unsigned int uiFilled2=0; unsigned int uiFilled3=0;
    while(dSliceLine<dSliceEnd-1.e-9)
    {
        uiOut=0;
        if(dSliceLine+dt>=dIntegrationEnd+1.e-9&&uiFilled1==0)
        {
            dSliceLine=dIntegrationEnd;
            uiOut=1;
        }
        if(dSliceLine+dt>=dFirstBondDate+1.e-9&&dFirstBondDate!=dIntegrationEnd&&dFirstBondDate!=dSecondBondDate&&uiOut==0&&uiFilled2==0)
        {
            dSliceLine=dFirstBondDate;
            uiOut=1;
        }
        if(dSliceLine+dt>=dSecondBondDate+1.e-9&&uiOut==0&&uiFilled3==0)
        {
            dSliceLine=dSecondBondDate;
            uiOut=1;
        }
        if(uiOut==0) dSliceLine=dSliceLine+dt;
        dSlices.insert(dSliceLine);
        // Checking
        if(dSliceLine==dIntegrationEnd) uiFilled1=1;
        if(dSliceLine==dFirstBondDate) uiFilled2=1;
        if(dSliceLine==dSecondBondDate) uiFilled3=1;
    }

    /*
    	int howManySteps=int(4.*(dSecondBondDate-dIntegrationStart));
    	if(howManySteps==0) howManySteps=1;
    	DKMaille<double> inputDates;
    	inputDates.insert(dIntegrationStart);
    	inputDates.insert(dIntegrationEnd);
    	if(dFirstBondDate>dIntegrationEnd+1.e-9) inputDates.insert(dFirstBondDate);
    	if(dSecondBondDate>dIntegrationEnd+1.e-9&&fabs(dSecondBondDate-dFirstBondDate)>1.e-9) inputDates.insert(dSecondBondDate);
      datesForSlices(inputDates, 
                     howManySteps, 
                     dSlices);
    */
}


void dIntermediateCalcSwap(double dSpotDate,
                           const DKMaille<double> &dDiscountDates,
                           const DKMaille<double> &dDiscountRates,
                           const DKMaille<double> &dAdjustedDiscountRates,
                           const DKMaille<double> &dSwaptionDates,
                           const DKMaille<double> &dBasisSwaptionDates,
                           const DKMaille<double> &dAccrualPeriods,
                           const DKMaille<double> &dBasis,
                           const DKMaille<double> &dBasisAccrualPeriods,
                           DKMaille<double> &dDiscountFactors,
                           DKMaille<double> &dDiscountFactorsOnBasisDates,
                           double *dSwapRate,
                           double *dAnnuity,
                           double *dBasisAnnuity,
                           double dJulianObservationDate)
{
    double dPrice=0.;
    DKMaille<double> dZCRates(dDiscountDates.entries());
    DKMaille<double> dZCAdjRates(dDiscountDates.entries());
    FromDiscountToZero(dSpotDate,dZCRates,dZCAdjRates,dDiscountDates,dDiscountRates,dAdjustedDiscountRates);
    // DKMaille<double> dDiscountFactors(dSwaptionDates.entries());
    // DKMaille<double> dDiscountFactorsOnBasisDates(dBasisSwaptionDates.entries());
    // Build SwapRate analytics

    *dAnnuity = 0.;
    *dSwapRate = 0.;
    for(unsigned int ui=0;ui<dSwaptionDates.entries();ui++)
    {
        dDiscountFactors.at(ui)=exp(-rateinterpolation_dk_maille(2,dSwaptionDates[ui],dDiscountDates,
                                    dZCAdjRates,dDiscountDates.entries()-1)*(dSwaptionDates[ui]-dSpotDate)/365.);
        if(ui>0)
            *dAnnuity+=dAccrualPeriods[ui]*dDiscountFactors.at(ui);
    }
    *dSwapRate = (dDiscountFactors.at(0)-dDiscountFactors.at(dSwaptionDates.entries()-1))/(*dAnnuity);
    // Calculate New Swap Rate
    *dBasisAnnuity=0.;
    for(ui=0;ui<dBasisSwaptionDates.entries();ui++)
    {
        dDiscountFactorsOnBasisDates.at(ui)=exp(-rateinterpolation_dk_maille(2,dBasisSwaptionDates[ui],dDiscountDates,
                                                dZCAdjRates,dDiscountDates.entries()-1)*(dBasisSwaptionDates[ui]-dSpotDate)/365.);
        *dBasisAnnuity+=dBasis[ui]*dBasisAccrualPeriods[ui]*dDiscountFactorsOnBasisDates.at(ui);
    }
    *dSwapRate+=(*dBasisAnnuity/(*dAnnuity));
}

double GetSwaptionAbsVol_VFDK_HW1F(double dSpotDate,
                                   DKMaille<double> dSwaptionDates,
                                   const DKMaille<double> &dAccrualPeriods,
                                   DKMaille<double> dBasisSwaptionDates,
                                   const DKMaille<double> &dBasis,
                                   const DKMaille<double> &dBasisAccrualPeriods,
                                   double dNoticePeriod,
                                   double dMeanReversion,
                                   DKMaille<double> dVolStripDates,
                                   DKMaille<double> &dVolStrip,
                                   const DKMaille<double> &dDiscountFactors,
                                   const DKMaille<double> &dDiscountFactorsOnBasisDates,
                                   double dSwapRate,
                                   double dBasisAnnuity,
                                   double dAnnuity,
                                   double dJulianObservationDate)
{

    double dPrice=0.;
    double dObservationDate=(dJulianObservationDate-dSpotDate)/365.;

    double dTerminalDate    = (dSwaptionDates[0]-dSpotDate)/365.-dNoticePeriod;
    for(unsigned int ui=0;ui<dVolStripDates.entries(); ui++)
        dVolStripDates[ui]=(dVolStripDates[ui]-dSpotDate)/365.;

    for(ui=0;ui<dSwaptionDates.entries(); ui++)
        dSwaptionDates[ui]=(dSwaptionDates[ui]-dSpotDate)/365.;

    for(ui=0;ui<dBasisSwaptionDates.entries(); ui++)
        dBasisSwaptionDates[ui]=(dBasisSwaptionDates[ui]-dSpotDate)/365.;

    DKMaille<double> dVector(dSwaptionDates.entries()+dBasisSwaptionDates.entries());
    DKMaille<double> dDiscounts(dDiscountFactors.entries()+dDiscountFactorsOnBasisDates.entries());
    DKMaille<double> dCash(dDiscountFactors.entries()+dDiscountFactorsOnBasisDates.entries());
    DKMaille<double> dDates(dSwaptionDates.entries()+dBasisSwaptionDates.entries());

    if(
        (dSwaptionDates.entries()!=dDiscountFactors.entries())
        ||(dSwaptionDates.entries()!=dAccrualPeriods.entries())
        ||(dBasisSwaptionDates.entries()!=dBasisAccrualPeriods.entries())
        ||(dBasisSwaptionDates.entries()!=dDiscountFactorsOnBasisDates.entries())
        ||(dBasisSwaptionDates.entries()!=dBasis.entries())
    )
        throw("Input deal dates are inconsistent");

    for(ui=0;ui<dVector.entries();ui++)
    {
        if(ui<dAccrualPeriods.entries())
        {
            dCash[ui]=dAccrualPeriods[ui]*dSwapRate;
            dDiscounts[ui]=dDiscountFactors[ui];
            dDates[ui]=dSwaptionDates[ui];
        }
        else
        {
            dCash[ui]=dBasisAccrualPeriods[ui-dAccrualPeriods.entries()]*dBasis[ui-dAccrualPeriods.entries()];
            dDiscounts[ui]=dDiscountFactorsOnBasisDates[ui-dAccrualPeriods.entries()];
            dDates[ui]=dBasisSwaptionDates[ui-dAccrualPeriods.entries()];
        }
    }
    dVector=GetItoVector(dSwaptionDates.entries(),dBasisSwaptionDates.entries(),dCash,dDiscounts);

    // New functionality
    // Make sure integration slices are flexible to allow for fast pricing of callable reverse floaters and callable CMS spread options
    // Checks on dates
    DKMaille<double> dIntegrationStrip=dCheckDates(dVolStripDates,dObservationDate,dTerminalDate);
    DKMaille<double> dIntegrationStripVols=dNewVols(dVolStripDates,dVolStrip,dIntegrationStrip);


    // Date tests
    bool bCheckDis=true;
    unsigned int uiM=0;
    for(;uiM<dIntegrationStrip.entries();uiM++)
    {
        if(fabs(dIntegrationStrip[uiM]-dTerminalDate)<0.00001)
            bCheckDis=false;
    }
    if(bCheckDis)
        throw("Dates of Std Dev Strip Are Inconsistent with Option Expiry");

    double dForwardStdDev=0.;
    double dIntermediateCalc=0.;

    // This is the core calculation of the VFDK approximation
    // VFDK approximation for time dependent piece-wise constant standard deviation
    // Two dimensional Kac Determinant applied through the entire space of dDates
    unsigned int uiK=0;
    for(;uiK<dDates.entries();uiK++)
    {
        unsigned int uiM=0;
        for(;uiM<dDates.entries();uiM++)
        {
            dIntermediateCalc=dVector[uiK]
                              * dVector[uiM]
                              * GetItoDeterminant_HW1F_TV(dMeanReversion,
                                                          dObservationDate,
                                                          dTerminalDate,
                                                          dIntegrationStrip,
                                                          dIntegrationStripVols,
                                                          dDates[uiK],
                                                          dDates[uiM]);

            dForwardStdDev += dIntermediateCalc;
        }
    }
    // This is the end of the calculation
    if(dForwardStdDev<0.) throw("Implied Negative Standard Deviation");
    dForwardStdDev = (1./dAnnuity) * fabs(1./dMeanReversion) * sqrt(dForwardStdDev/(dTerminalDate-dObservationDate));
    return dForwardStdDev;
} // GetSwaptionAbsVol_VFDK_HW1F






double GetItoDeterminant2F(double dMeanFieldDecay, double dMeanFieldDecay2,
                           double dObservationDate, double dTerminalDate,
                           double dSigma1, double dSigma2, double dDate1,
                           double dDate2)
{
    double dTerm1 = (dTerminalDate-dObservationDate);
    double dTerm2 = -(1. /dMeanFieldDecay)*(exp(-dMeanFieldDecay*(dDate1-dTerminalDate))
                                            -exp(-dMeanFieldDecay*(dDate1-dObservationDate)));
    double dTerm3 = -(1. /dMeanFieldDecay2)*(exp(-dMeanFieldDecay2*(dDate2-dTerminalDate))
                    -exp(-dMeanFieldDecay2*(dDate2-dObservationDate)));
    double dTerm4 = +(1./(dMeanFieldDecay+dMeanFieldDecay2))
                    *(
                        exp(-dMeanFieldDecay*(dDate1-dTerminalDate)-dMeanFieldDecay2*(dDate2-dTerminalDate))
                        -exp(-dMeanFieldDecay*(dDate1-dObservationDate)-dMeanFieldDecay2*(dDate2-dObservationDate))
                    );
    return (dTerm1+dTerm2+dTerm3+dTerm4)*dSigma1*dSigma2;
} // double GetItoDeterminant2F


double GetItoDeterminant_HW2F_TV(double dMeanFieldDecay,
                                 double dMeanFieldDecay2,
                                 double dObservationDate,
                                 double dTerminalDate,
                                 const DKMaille<double> &dDatesStrip,
                                 const DKMaille<double> &dSigmaStrip,
                                 const DKMaille<double> &dSigmaStrip2,
                                 double dDate1,
                                 double dDate2)
{
    int i=0;
    double dSum=0.;
    while(dDatesStrip[i]<dTerminalDate-0.00001)
    {
        dSum+=GetItoDeterminant2F(dMeanFieldDecay,dMeanFieldDecay2,dDatesStrip[i],dDatesStrip[i+1],
                                  dSigmaStrip[i],dSigmaStrip2[i],dDate1,dDate2);
        i++;
    }
    return dSum;
} // GetItoDeterminant_HW2F_TV

double GetSwaptionAbsVol_VFDK_HW2F(double dSpotDate,
                                   DKMaille<double> dSwaptionDates,
                                   const DKMaille<double> &dAccrualPeriods,
                                   DKMaille<double> dBasisSwaptionDates,
                                   const DKMaille<double> &dBasis,
                                   const DKMaille<double> &dBasisAccrualPeriods,
                                   double dNoticePeriod,
                                   double dMeanReversion1,
                                   double dMeanReversion2,
                                   double dRelativeFactor,
                                   double dCorrelation,
                                   DKMaille<double> dVolStripDates,
                                   DKMaille<double> &dVolStrip,
                                   const DKMaille<double> &dDiscountFactors,
                                   const DKMaille<double> &dDiscountFactorsOnBasisDates,
                                   double dSwapRate,
                                   double dBasisAnnuity,
                                   double dAnnuity,
                                   double dJulianObservationDate)
{

    double dPrice=0.;
    double dObservationDate=(dJulianObservationDate-dSpotDate)/365.;
    double dTerminalDate    = (dSwaptionDates[0]-dSpotDate)/365.-dNoticePeriod;

    for(unsigned int ui=0;ui<dVolStripDates.entries(); ui++)
        dVolStripDates[ui]=(dVolStripDates[ui]-dSpotDate)/365.;

    for(ui=0;ui<dSwaptionDates.entries(); ui++)
        dSwaptionDates[ui]=(dSwaptionDates[ui]-dSpotDate)/365.;
    for(ui=0;ui<dBasisSwaptionDates.entries(); ui++)
        dBasisSwaptionDates[ui]=(dBasisSwaptionDates[ui]-dSpotDate)/365.;

    DKMaille<double> dVector(dSwaptionDates.entries()+dBasisSwaptionDates.entries());
    DKMaille<double> dDiscounts(dDiscountFactors.entries()+dDiscountFactorsOnBasisDates.entries());
    DKMaille<double> dCash(dDiscountFactors.entries()+dDiscountFactorsOnBasisDates.entries());
    DKMaille<double> dDates(dSwaptionDates.entries()+dBasisSwaptionDates.entries());

    if(
        (dSwaptionDates.entries()!=dDiscountFactors.entries())
        ||(dSwaptionDates.entries()!=dAccrualPeriods.entries())
        ||(dBasisSwaptionDates.entries()!=dBasisAccrualPeriods.entries())
        ||(dBasisSwaptionDates.entries()!=dDiscountFactorsOnBasisDates.entries())
        ||(dBasisSwaptionDates.entries()!=dBasis.entries())
    )
        throw("Input deal dates are inconsistent");

    for(ui=0;ui<dVector.entries();ui++)
    {
        if(ui<dAccrualPeriods.entries())
        {
            dCash[ui]=dAccrualPeriods[ui]*dSwapRate;
            dDiscounts[ui]=dDiscountFactors[ui];
            dDates[ui]=dSwaptionDates[ui];
        }
        else
        {
            dCash[ui]=dBasisAccrualPeriods[ui-dAccrualPeriods.entries()]*dBasis[ui-dAccrualPeriods.entries()];
            dDiscounts[ui]=dDiscountFactorsOnBasisDates[ui-dAccrualPeriods.entries()];
            dDates[ui]=dBasisSwaptionDates[ui-dAccrualPeriods.entries()];
        }
    }
    dVector=GetItoVector(dSwaptionDates.entries(),dBasisSwaptionDates.entries(),dCash,dDiscounts);

    // New functionality
    // Make sure integration slices are flexible to allow for fast pricing of callable reverse floaters and callable CMS spread options
    // Checks on dates
    DKMaille<double> dIntegrationStrip=dCheckDates(dVolStripDates,dObservationDate,dTerminalDate);
    DKMaille<double> dIntegrationStripVols=dNewVols(dVolStripDates,dVolStrip,dIntegrationStrip);


    // Date tests
    bool bCheckDis=true;
    unsigned int uiM=0;
    for(;uiM<dIntegrationStrip.entries();uiM++)
    {
        if(fabs(dIntegrationStrip[uiM]-dTerminalDate)<0.00001)
            bCheckDis=false;
    }
    if(bCheckDis)
        throw("Dates of Std Dev Strip Are Inconsistent with Option Expiry");

    double dForwardStdDev=0.;
    double dIntermediateCalc=0.;
    // This is the core calculation of the VFDK approximation
    // Create second index
    DKMaille<double> dVolStrip2(dIntegrationStrip.entries());
    unsigned int uiK=0;
    for(;uiK<dVolStrip2.entries();uiK++) dVolStrip2[uiK]=dRelativeFactor*dIntegrationStripVols[uiK];
    // iterate through Kac space
    for(uiK=0;uiK<dDates.entries();uiK++)
    {
        unsigned int uiM=0;
        for(;uiM<dDates.entries();uiM++)
        {
            dIntermediateCalc=dVector[uiK]*dVector[uiM]*
                              (
                                  (1./(dMeanReversion1*dMeanReversion1))*
                                  GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion1,
                                                            dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dIntegrationStripVols,dDates[uiK],dDates[uiM])
                                  +(1./(dMeanReversion2*dMeanReversion2))*
                                  GetItoDeterminant_HW2F_TV(dMeanReversion2,dMeanReversion2,
                                                            dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip2,dVolStrip2,dDates[uiK],dDates[uiM])
                                  +2.*dCorrelation*
                                  (1./(dMeanReversion1*dMeanReversion2))*
                                  GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion2,
                                                            dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dVolStrip2,dDates[uiK],dDates[uiM])
                              );
            dForwardStdDev += dIntermediateCalc;
        }
    }
    // This is the end of the calculation
    if(dForwardStdDev<0.) throw("Implied Negative Standard Deviation");
    dForwardStdDev = (1./dAnnuity) * sqrt(dForwardStdDev/(dTerminalDate-dObservationDate));
    return dForwardStdDev;
} // GetSwaptionAbsVol_VFDK_HW2F





double GetSwaptionAbsVol_VFDK_HW3F(double dSpotDate,
                                   DKMaille<double> dSwaptionDates,
                                   const DKMaille<double> &dAccrualPeriods,
                                   DKMaille<double> dBasisSwaptionDates,
                                   const DKMaille<double> &dBasis,
                                   const DKMaille<double> &dBasisAccrualPeriods,
                                   double dNoticePeriod,
                                   double dMeanReversion1,
                                   double dMeanReversion2,
                                   double dMeanReversion3,
                                   double dRelativeFactor12,
                                   double dRelativeFactor13,
                                   double dCorrelation12,
                                   double dCorrelation13,
                                   double dCorrelation23,
                                   DKMaille<double> dVolStripDates,
                                   DKMaille<double> &dVolStrip,
                                   const DKMaille<double> &dDiscountFactors,
                                   const DKMaille<double> &dDiscountFactorsOnBasisDates,
                                   double dSwapRate,
                                   double dBasisAnnuity,
                                   double dAnnuity,
                                   double dJulianObservationDate)
{

    double dPrice=0.;
    double dObservationDate=(dJulianObservationDate-dSpotDate)/365.;

    for(unsigned int ui=0;ui<dVolStripDates.entries(); ui++)
        dVolStripDates[ui]=(dVolStripDates[ui]-dSpotDate)/365.;

    double dTerminalDate    = (dSwaptionDates[0]-dSpotDate)/365.-dNoticePeriod;


    for(ui=0;ui<dSwaptionDates.entries(); ui++)
        dSwaptionDates[ui]=(dSwaptionDates[ui]-dSpotDate)/365.;
    for(ui=0;ui<dBasisSwaptionDates.entries(); ui++)
        dBasisSwaptionDates[ui]=(dBasisSwaptionDates[ui]-dSpotDate)/365.;

    DKMaille<double> dVector(dSwaptionDates.entries()+dBasisSwaptionDates.entries());
    DKMaille<double> dDiscounts(dDiscountFactors.entries()+dDiscountFactorsOnBasisDates.entries());
    DKMaille<double> dCash(dDiscountFactors.entries()+dDiscountFactorsOnBasisDates.entries());
    DKMaille<double> dDates(dSwaptionDates.entries()+dBasisSwaptionDates.entries());

    if(
        (dSwaptionDates.entries()!=dDiscountFactors.entries())
        ||(dSwaptionDates.entries()!=dAccrualPeriods.entries())
        ||(dBasisSwaptionDates.entries()!=dBasisAccrualPeriods.entries())
        ||(dBasisSwaptionDates.entries()!=dDiscountFactorsOnBasisDates.entries())
        ||(dBasisSwaptionDates.entries()!=dBasis.entries())
    )
        throw("Input deal dates are inconsistent");

    for(ui=0;ui<dVector.entries();ui++)
    {
        if(ui<dAccrualPeriods.entries())
        {
            dCash[ui]=dAccrualPeriods[ui]*dSwapRate;
            dDiscounts[ui]=dDiscountFactors[ui];
            dDates[ui]=dSwaptionDates[ui];
        }
        else
        {
            dCash[ui]=dBasisAccrualPeriods[ui-dAccrualPeriods.entries()]*dBasis[ui-dAccrualPeriods.entries()];
            dDiscounts[ui]=dDiscountFactorsOnBasisDates[ui-dAccrualPeriods.entries()];
            dDates[ui]=dBasisSwaptionDates[ui-dAccrualPeriods.entries()];
        }
    }
    dVector=GetItoVector(dSwaptionDates.entries(),dBasisSwaptionDates.entries(),dCash,dDiscounts);

    // New functionality
    // Make sure integration slices are flexible to allow for fast pricing of callable reverse floaters and callable CMS spread options
    // Checks on dates
    DKMaille<double> dIntegrationStrip=dCheckDates(dVolStripDates,dObservationDate,dTerminalDate);
    DKMaille<double> dIntegrationStripVols=dNewVols(dVolStripDates,dVolStrip,dIntegrationStrip);


    // Date tests
    bool bCheckDis=true;
    unsigned int uiM=0;
    for(;uiM<dIntegrationStrip.entries();uiM++)
    {
        if(fabs(dIntegrationStrip[uiM]-dTerminalDate)<0.00001)
            bCheckDis=false;
    }
    if(bCheckDis)
        throw("Dates of Std Dev Strip Are Inconsistent with Option Expiry");

    double dForwardStdDev=0.;
    double dIntermediateCalc=0.;
    // This is the core calculation of the VFDK approximation
    // Create second index
    // Create third index
    DKMaille<double> dVolStrip2(dIntegrationStrip.entries());
    DKMaille<double> dVolStrip3(dIntegrationStrip.entries());

    for(unsigned int uiK=0;uiK<dIntegrationStrip.entries();uiK++) dVolStrip2[uiK]=dRelativeFactor12*dIntegrationStripVols[uiK];
    for(uiK=0;uiK<dIntegrationStrip.entries();uiK++) dVolStrip3[uiK]=dRelativeFactor13*dIntegrationStripVols[uiK];

    // iterate through Kac space
    for(uiK=0;uiK<dDates.entries();uiK++)
    {
        unsigned int uiM=0;
        for(;uiM<dDates.entries();uiM++)
        {
            dIntermediateCalc=dVector[uiK]*dVector[uiM]*
                              (
                                  (1./(dMeanReversion1*dMeanReversion1))*
                                  GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion1,
                                                            dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dVolStrip,dDates[uiK],dDates[uiM])
                                  +(1./(dMeanReversion2*dMeanReversion2))*
                                  GetItoDeterminant_HW2F_TV(dMeanReversion2,dMeanReversion2,
                                                            dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip2,dVolStrip2,dDates[uiK],dDates[uiM])
                                  +(1./(dMeanReversion3*dMeanReversion3))*
                                  GetItoDeterminant_HW2F_TV(dMeanReversion3,dMeanReversion3,
                                                            dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip3,dVolStrip3,dDates[uiK],dDates[uiM])
                                  +2.*dCorrelation12*
                                  (1./(dMeanReversion1*dMeanReversion2))*
                                  GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion2,
                                                            dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dVolStrip2,dDates[uiK],dDates[uiM])
                                  +2.*dCorrelation13*
                                  (1./(dMeanReversion1*dMeanReversion3))*
                                  GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion3,
                                                            dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dVolStrip3,dDates[uiK],dDates[uiM])
                                  +2.*dCorrelation23*
                                  (1./(dMeanReversion2*dMeanReversion3))*
                                  GetItoDeterminant_HW2F_TV(dMeanReversion2,dMeanReversion3,
                                                            dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip2,dVolStrip3,dDates[uiK],dDates[uiM])
                              );
            dForwardStdDev += dIntermediateCalc;
        }
    }
    // This is the end of the calculation

    if(dForwardStdDev<0.) throw("Implied Negative Standard Deviation");

    dForwardStdDev = (1./dAnnuity) * sqrt(dForwardStdDev/(dTerminalDate-dObservationDate));
    return dForwardStdDev;
} // GetSwaptionAbsVol_VFDK_HW3F

double GetOptionPriceBS(double dCallPut,
                        double dSwapRate,
                        double dAnnuity,
                        double dStrike,
                        double dForwardVol,
                        double dTimeToExpiry)
{
    double d1=(log(dSwapRate/dStrike)+0.5*dForwardVol*dForwardVol*dTimeToExpiry)/(dForwardVol*sqrt(dTimeToExpiry));
    double d2=(log(dSwapRate/dStrike)-0.5*dForwardVol*dForwardVol*dTimeToExpiry)/(dForwardVol*sqrt(dTimeToExpiry));
    double dND1=0.;
    double dND2=0.;
    normalDK(&dND1,d1);
    normalDK(&dND2,d2);
    double dCall=dSwapRate*dND1-dStrike*dND2;
    double dResult;
    if(dCallPut==0.) dResult=dCall;
    if(dCallPut==1.) dResult=dCall-dSwapRate+dStrike;
    return dResult;
} // GetOptionPrice


double GetOptionPrice(double dCallPut,
                      double dSwapRate,
                      double dAnnuity,
                      double dStrike,
                      double dForwardStdDev,
                      double dTimeToExpiry)
{
    double denominator = (dForwardStdDev*sqrt(dTimeToExpiry));
    double dArgument = (dSwapRate-dStrike)/denominator;
    if(dCallPut==1.) dArgument*=-1.;
    double dNormsDist;
    normalDK(&dNormsDist,dArgument);
    double dPi = 3.141592653589790 ;
    double dDerivativeTerm = dForwardStdDev * sqrt((dTimeToExpiry)/(2.*dPi)) *
                             exp(-0.5*(dArgument*dArgument));
    double dPrice;
    if(dCallPut==0.) dPrice = 10000. * dAnnuity * ((dSwapRate-dStrike)*dNormsDist+dDerivativeTerm);
    if(dCallPut==1.) dPrice = 10000. * dAnnuity * (-(dSwapRate-dStrike)*dNormsDist+dDerivativeTerm);
    return dPrice;
} // GetOptionPrice




double GetAdjSwapRate(double dSpotDate,
                      const DKMaille<double> &dSwaptionDates,
                      const DKMaille<double> &dAccrualPeriods,
                      const DKMaille<double> &dBasisSwaptionDates,
                      const DKMaille<double> &dBasis,
                      const DKMaille<double> &dBasisAccrualPeriods,
                      const DKMaille<double> &dDiscountDates,
                      const DKMaille<double> &dDiscountRates,
                      const DKMaille<double> &dAdjustedDiscountRates)
{
    DKMaille<double> dDiscountFactors(dSwaptionDates.entries());
    DKMaille<double> dDiscountFactorsOnBasisDates(dBasisSwaptionDates.entries());
    DKMaille<double> dPeriods(dAccrualPeriods.entries());
    DKMaille<double> dZCRates(dDiscountDates.entries());
    DKMaille<double> dZCAdjRates(dDiscountDates.entries());
    FromDiscountToZero(dSpotDate,dZCRates,dZCAdjRates,dDiscountDates,dDiscountRates,dAdjustedDiscountRates);
    dPeriods[0]=dAccrualPeriods[0];
    double dAnnuity = 0.;
    double dSwapRate = 0.;
    for(unsigned ui=0;ui<dSwaptionDates.entries();ui++)
    {
        dDiscountFactors[ui]=exp(-rateinterpolation_dk_maille(2,dSwaptionDates[ui],dDiscountDates,
                                 dZCAdjRates,dDiscountDates.entries()-1)*(dSwaptionDates[ui]-dSpotDate)/365.);
        if(ui>0)
        {
            dPeriods[ui]=dAccrualPeriods[ui];
            dAnnuity+=dPeriods[ui]*dDiscountFactors[ui];
        }
    }
    dSwapRate = (dDiscountFactors[0]-dDiscountFactors[dSwaptionDates.entries()-1])/dAnnuity;

    // Calculate New Swap Rate
    double dBasisAnnuity=0.;
    for(ui=0;ui<dBasisSwaptionDates.entries();ui++)
    {
        dDiscountFactorsOnBasisDates[ui]=exp(-rateinterpolation_dk_maille(2,dBasisSwaptionDates[ui],dDiscountDates,
                                             dZCAdjRates,dDiscountDates.entries()-1)*(dBasisSwaptionDates[ui]-dSpotDate)/365.); ;
        dBasisAnnuity+=dBasis[ui]*dBasisAccrualPeriods[ui]*dDiscountFactorsOnBasisDates[ui];
    }

    dSwapRate+=(dBasisAnnuity/dAnnuity);

    return dSwapRate;
} // double GetAdjSwapRate

double GetPriceFromAbsVol(double dSpotDate,
                          double dJulianObservationDate,
                          const DKMaille<double> &dSwaptionDates,
                          const DKMaille<double> &dAccrualPeriods,
                          double dNoticePeriod,
                          double dAbsoluteVol,
                          const DKMaille<double> &dDiscountDates,
                          const DKMaille<double> &dDiscountRates,
                          const DKMaille<double> &dAdjustedDiscountRates)
{
    double dBlackScholesPrice =0.;
    DKMaille<double> dDiscountFactors(dSwaptionDates.entries());
    DKMaille<double> dZCRates(dDiscountDates.entries());
    DKMaille<double> dZCAdjRates(dDiscountDates.entries());
    FromDiscountToZero(dSpotDate,dZCRates,dZCAdjRates,dDiscountDates,dDiscountRates,dAdjustedDiscountRates);
    double dAnnuity = 0.;
    double dSwapRate = 0.;
    unsigned int ui=0;
    for(;ui<dSwaptionDates.entries();ui++)
    {
        dDiscountFactors[ui]=exp(-rateinterpolation_dk_maille(2,dSwaptionDates[ui],dDiscountDates,
                                 dZCRates,dDiscountDates.entries()-1)*(dSwaptionDates[ui]-dSpotDate)/365.);
        if(ui>0)
        {
            dAnnuity+=dAccrualPeriods[ui]*dDiscountFactors[ui];
        }
    }
    dSwapRate = (dDiscountFactors[0]-dDiscountFactors[dSwaptionDates.entries()-1])/dAnnuity;
    double dObservationDate   =    0.;
    double dTerminalDate      =    (dSwaptionDates[0]-dJulianObservationDate)/365.-dNoticePeriod;
    double dPi = 3.141592653589790 ;
    double dFactor = sqrt(dTerminalDate)/sqrt(2.*dPi);
    return dAnnuity*dFactor*dAbsoluteVol*10000.;
} // GetPriceFromAbsVol()




double  Bootstrapping1F(double dSwapStart,
                        double dSpotDate,
                        unsigned int k,                                       // notice index
                        const DKMaille2D<double> &dVolInfo,										// all info regarding vols
                        double dAbsVol,                                       // AbsVol
                        const DKMaille<double> &swaptionDates,                // notice dates from the fixed leg
                        const DKMaille<double> &accrualPeriods,               // accrual periods from fixed leg
                        const DKMaille<double> &basisSwaptionDates,           // notice dates from the float leg
                        const DKMaille<double> &basisAccrualPeriods,          // accrual periods from the float leg
                        const DKMaille<double> &basis,                        // basis from float leg
                        const DKMaille<double> &dModelParameters,
                        const DKMaille<double> &dDiscountDates,
                        const DKMaille<double> &dDiscountRates,
                        const DKMaille<double> &dAdjustedDiscountRates,
                        double dJulianObservationDate
                       )
{

    DKMaille<double> dVolStripDates(k+2);
    DKMaille<double> dVolStrip(k+2);

    dVolStripDates.at(0)=dSpotDate;

    // fill past vols and dates
    for(unsigned int ii=0;ii<k+1;ii++)
    {
        dVolStripDates.at(ii+1)=dVolInfo.at(ii,0);
        if(k>0)
        {
            dVolStrip.at(ii)=dVolInfo.at(ii,1);
        }
    } // for ii

    // Get Model parameters
    double dMeanReversion1 = dModelParameters.at(1);

    // Get AdjSwap Rate
    double dAdjSwapRate=GetAdjSwapRate(dSpotDate,
                                       swaptionDates,
                                       accrualPeriods,
                                       basisSwaptionDates,
                                       basis,
                                       basisAccrualPeriods,
                                       dDiscountDates,
                                       dDiscountRates,
                                       dAdjustedDiscountRates);

    // Get AdjSwap Annuity
    double dAdjAnnuity=GetAdjAnnuity(dSpotDate,
                                     swaptionDates,
                                     accrualPeriods,
                                     dDiscountDates,
                                     dDiscountRates,
                                     dAdjustedDiscountRates,
                                     dJulianObservationDate);

    // Get Swap Rate
    double dSwapRate=GetSwapRate(dSpotDate,
                                 swaptionDates,
                                 accrualPeriods,
                                 dDiscountDates,
                                 dDiscountRates,
                                 dAdjustedDiscountRates);

    double dNoticeDate = dVolInfo.at(k,0);

    double dNoticePeriod = (swaptionDates.at(0)-dNoticeDate)/365.;

    // Catch on dates
    if(dNoticePeriod<0.)
    {
        throw("Notice Date After Swap Start");
    }

    // Get straddle price from AbsoluteVol
    double dMarketPrice=2.*GetPriceFromAbsVol(dSpotDate,
                        dJulianObservationDate,
                        swaptionDates,
                        accrualPeriods,
                        dNoticePeriod,
                        dAbsVol,
                        dDiscountDates,
                        dDiscountRates,
                        dAdjustedDiscountRates);

    double dGuess;

    if(0==k)
    {
        dGuess = dAbsVol;
    }
    else
    {
        dGuess = dVolStrip.at(k-1);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Look to see if there is NO solution i.e. we don't need Newton-Raphson
    /////////////////////////////////////////////////////////////////////////////
    dVolStrip.at(k)=0.0001;

    DKMaille<double> dDiscountFactors(swaptionDates.entries());
    DKMaille<double> dDiscountFactorsOnBasisDates(basisSwaptionDates.entries());
    double dSwapRate_=0.; double dAnnuity_=0.; double dBasisAnnuity_=0.;
    dIntermediateCalcSwap(dSpotDate,
                          dDiscountDates,
                          dDiscountRates,
                          dAdjustedDiscountRates,
                          swaptionDates,
                          basisSwaptionDates,
                          accrualPeriods,
                          basis,
                          basisAccrualPeriods,
                          dDiscountFactors,
                          dDiscountFactorsOnBasisDates,
                          &dSwapRate_,
                          &dAnnuity_,
                          &dBasisAnnuity_,
                          dJulianObservationDate);

    double dFwdStdDev = GetSwaptionAbsVol_VFDK_HW1F(dSpotDate,
                        swaptionDates,
                        accrualPeriods,
                        basisSwaptionDates,
                        basis,
                        basisAccrualPeriods,
                        dNoticePeriod,
                        dMeanReversion1,
                        dVolStripDates,
                        dVolStrip,
                        dDiscountFactors,
                        dDiscountFactorsOnBasisDates,
                        dSwapRate_,
                        dBasisAnnuity_,
                        dAnnuity_,
                        dJulianObservationDate);

    double dFunctionValue = GetOptionPrice(0.,dAdjSwapRate,dAdjAnnuity,dSwapRate,dFwdStdDev,(dNoticeDate-dJulianObservationDate)/365.)
                            +GetOptionPrice(1.,dAdjSwapRate,dAdjAnnuity,dSwapRate,dFwdStdDev,(dNoticeDate-dJulianObservationDate)/365.)
                            -dMarketPrice;

    if(dFunctionValue >= 1.) //if positive, tolerance = 1 BP
    {
        return 0.; // failure state - no solution ever (???)
    }

    /////////////////////////////////////////////////////////////////////////////
    // End of check pre-NR check ...
    /////////////////////////////////////////////////////////////////////////////

    // Newton-Raphson loop to get fwd std dev
    int iCounter=0;

    while(iCounter<100)
    {
        iCounter++;
        dVolStrip.at(k)=dGuess;

        // Function
        double dFwdStdDev = GetSwaptionAbsVol_VFDK_HW1F(dSpotDate,
                            swaptionDates,
                            accrualPeriods,
                            basisSwaptionDates,
                            basis,
                            basisAccrualPeriods,
                            dNoticePeriod,
                            dMeanReversion1,
                            dVolStripDates,
                            dVolStrip,
                            dDiscountFactors,
                            dDiscountFactorsOnBasisDates,
                            dSwapRate_,
                            dBasisAnnuity_,
                            dAnnuity_,
                            dJulianObservationDate);

        if(dAdjSwapRate!=dSwapRate)
        {
            dFunctionValue=GetOptionPrice(0.,dAdjSwapRate,dAdjAnnuity,dSwapRate,dFwdStdDev,(dNoticeDate-dJulianObservationDate)/365.)
                           +GetOptionPrice(1.,dAdjSwapRate,dAdjAnnuity,dSwapRate,dFwdStdDev,(dNoticeDate-dJulianObservationDate)/365.)
                           -dMarketPrice;
        }
        else
        {
            dFunctionValue=2.*10000.*dAdjAnnuity*dFwdStdDev*sqrt((dNoticeDate-dJulianObservationDate)/365./(2.*g_dPi))-dMarketPrice;
        }

        // FunctionUp
        dVolStrip.at(k)=dGuess+0.0001;

        double dFwdStdDevUp=GetSwaptionAbsVol_VFDK_HW1F(dSpotDate,
                            swaptionDates,
                            accrualPeriods,
                            basisSwaptionDates,
                            basis,
                            basisAccrualPeriods,
                            dNoticePeriod,
                            dMeanReversion1,
                            dVolStripDates,
                            dVolStrip,
                            dDiscountFactors,
                            dDiscountFactorsOnBasisDates,
                            dSwapRate_,
                            dBasisAnnuity_,
                            dAnnuity_,
                            dJulianObservationDate);
        double dFunctionValueUp=0.;

        if(dAdjSwapRate!=dSwapRate)
        {
            dFunctionValueUp=GetOptionPrice(0.,dAdjSwapRate,dAdjAnnuity,dSwapRate,dFwdStdDevUp,(dNoticeDate-dJulianObservationDate)/365.)
                             +GetOptionPrice(1.,dAdjSwapRate,dAdjAnnuity,dSwapRate,dFwdStdDevUp,(dNoticeDate-dJulianObservationDate)/365.)
                             -dMarketPrice;
        }
        else
        {
            dFunctionValueUp=2.*10000.*dAdjAnnuity*dFwdStdDevUp*sqrt((dNoticeDate-dJulianObservationDate)/365./(2.*g_dPi))-dMarketPrice;
        }

        double dRatio=0.0001*dFunctionValue/(dFunctionValueUp-dFunctionValue);
        dGuess-=dRatio;

        // Force -ve guess to v.small +ve
        if (dGuess < 0.0)
        {
            dGuess = 0.0000001;
        }

        if(fabs(dRatio)<2.e-05)
        {
            dVolStrip.at(k)=dGuess;
            break;
        }
    } // while iCounter < 100

    if(dVolStrip.at(k)<0.0001)
        throw("Zero or Negative Vol in bootstrapping");
    return dVolStrip.at(k);
} // Bootstrapping1F(...)



double  Bootstrapping2F(double dSwapStart,
                        double dSpotDate,
                        unsigned int k,                                       // notice index
                        const DKMaille2D<double> &dVolInfo,										// all info regarding vols
                        double dAbsVol,                                       // AbsVol
                        const DKMaille<double> &swaptionDates,                // notice dates from the fixed leg
                        const DKMaille<double> &accrualPeriods,               // accrual periods from fixed leg
                        const DKMaille<double> &basisSwaptionDates,           // notice dates from the float leg
                        const DKMaille<double> &basisAccrualPeriods,          // accrual periods from the float leg
                        const DKMaille<double> &basis,                        // basis from float leg
                        const DKMaille<double> &dModelParameters,
                        const DKMaille<double> &dDiscountDates,
                        const DKMaille<double> &dDiscountRates,
                        const DKMaille<double> &dAdjustedDiscountRates,
                        double dJulianObservationDate
                       )
{

    DKMaille<double> dVolStripDates(k+2);
    DKMaille<double> dVolStrip(k+2);

    dVolStripDates.at(0)=dSpotDate;

    // fill past vols and dates
    for(unsigned int ii=0;ii<k+1;ii++)
    {
        dVolStripDates.at(ii+1)=dVolInfo.at(ii,0);
        if(k>0)
        {
            dVolStrip.at(ii)=dVolInfo.at(ii,1);
        }
    } // for ii

    // Get Model parameters
    double dMeanReversion1 = dModelParameters.at(1);
    double dMeanReversion2 = dModelParameters.at(2);
    double dRelativeFactor = dModelParameters.at(4);
    double dCorrelation    = dModelParameters.at(6);


    // Get AdjSwap Rate
    double dAdjSwapRate=GetAdjSwapRate(dSpotDate,
                                       swaptionDates,
                                       accrualPeriods,
                                       basisSwaptionDates,
                                       basis,
                                       basisAccrualPeriods,
                                       dDiscountDates,
                                       dDiscountRates,
                                       dAdjustedDiscountRates);

    // Get AdjSwap Annuity
    double dAdjAnnuity=GetAdjAnnuity(dSpotDate,
                                     swaptionDates,
                                     accrualPeriods,
                                     dDiscountDates,
                                     dDiscountRates,
                                     dAdjustedDiscountRates,
                                     dJulianObservationDate);

    // Get Swap Rate
    double dSwapRate=GetSwapRate(dSpotDate,
                                 swaptionDates,
                                 accrualPeriods,
                                 dDiscountDates,
                                 dDiscountRates,
                                 dAdjustedDiscountRates);

    double dNoticeDate = dVolInfo.at(k,0);

    double dNoticePeriod = (swaptionDates.at(0)-dNoticeDate)/365.;

    // Catch on dates
    if(dNoticePeriod<0.)
    {
        throw("Notice Date After Swap Start");
    }

    // Get straddle price from AbsoluteVol
    double dMarketPrice=2.*GetPriceFromAbsVol(dSpotDate,
                        dJulianObservationDate,
                        swaptionDates,
                        accrualPeriods,
                        dNoticePeriod,
                        dAbsVol,
                        dDiscountDates,
                        dDiscountRates,
                        dAdjustedDiscountRates);

    double dGuess;

    if(0==k)
    {
        dGuess = dAbsVol;
    }
    else
    {
        dGuess = dVolStrip.at(k-1);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Look to see if there is NO solution i.e. we don't need Newton-Raphson
    /////////////////////////////////////////////////////////////////////////////
    dVolStrip.at(k)=0.0001;

    DKMaille<double> dDiscountFactors(swaptionDates.entries());
    DKMaille<double> dDiscountFactorsOnBasisDates(basisSwaptionDates.entries());
    double dSwapRate_=0.; double dAnnuity_=0.; double dBasisAnnuity_=0.;
    dIntermediateCalcSwap(dSpotDate,
                          dDiscountDates,
                          dDiscountRates,
                          dAdjustedDiscountRates,
                          swaptionDates,
                          basisSwaptionDates,
                          accrualPeriods,
                          basis,
                          basisAccrualPeriods,
                          dDiscountFactors,
                          dDiscountFactorsOnBasisDates,
                          &dSwapRate_,
                          &dAnnuity_,
                          &dBasisAnnuity_,
                          dJulianObservationDate);

    double dFwdStdDev = GetSwaptionAbsVol_VFDK_HW2F(dSpotDate,
                        swaptionDates,
                        accrualPeriods,
                        basisSwaptionDates,
                        basis,
                        basisAccrualPeriods,
                        dNoticePeriod,
                        dMeanReversion1,
                        dMeanReversion2,
                        dRelativeFactor,
                        dCorrelation,
                        dVolStripDates,
                        dVolStrip,
                        dDiscountFactors,
                        dDiscountFactorsOnBasisDates,
                        dSwapRate_,
                        dBasisAnnuity_,
                        dAnnuity_,
                        dJulianObservationDate);

    double dFunctionValue = GetOptionPrice(0.,dAdjSwapRate,dAdjAnnuity,dSwapRate,dFwdStdDev,(dNoticeDate-dJulianObservationDate)/365.)
                            +GetOptionPrice(1.,dAdjSwapRate,dAdjAnnuity,dSwapRate,dFwdStdDev,(dNoticeDate-dJulianObservationDate)/365.)
                            -dMarketPrice;

    if(dFunctionValue >= 0.)
    {
        return 0.; // failure state - no solution ever (???)
    }

    /////////////////////////////////////////////////////////////////////////////
    // End of check pre-NR check ...
    /////////////////////////////////////////////////////////////////////////////

    // Newton-Raphson loop to get fwd std dev
    int iCounter=0;

    while(iCounter<100)
    {
        iCounter++;
        dVolStrip.at(k)=dGuess;

        // Function
        dFwdStdDev=GetSwaptionAbsVol_VFDK_HW2F(dSpotDate,
                                               swaptionDates,
                                               accrualPeriods,
                                               basisSwaptionDates,
                                               basis,
                                               basisAccrualPeriods,
                                               dNoticePeriod,
                                               dMeanReversion1,
                                               dMeanReversion2,
                                               dRelativeFactor,
                                               dCorrelation,
                                               dVolStripDates,
                                               dVolStrip,
                                               dDiscountFactors,
                                               dDiscountFactorsOnBasisDates,
                                               dSwapRate_,
                                               dBasisAnnuity_,
                                               dAnnuity_,
                                               dJulianObservationDate);


        dFunctionValue=GetOptionPrice(0.,dAdjSwapRate,dAdjAnnuity,dSwapRate,dFwdStdDev,(dNoticeDate-dJulianObservationDate)/365.)
                       +GetOptionPrice(1.,dAdjSwapRate,dAdjAnnuity,dSwapRate,dFwdStdDev,(dNoticeDate-dJulianObservationDate)/365.)
                       -dMarketPrice;

        // FunctionUp
        dVolStrip.at(k)=dGuess+0.0001;

        double dFwdStdDevUp=GetSwaptionAbsVol_VFDK_HW2F(dSpotDate,
                            swaptionDates,
                            accrualPeriods,
                            basisSwaptionDates,
                            basis,
                            basisAccrualPeriods,
                            dNoticePeriod,
                            dMeanReversion1,
                            dMeanReversion2,
                            dRelativeFactor,
                            dCorrelation,
                            dVolStripDates,
                            dVolStrip,
                            dDiscountFactors,
                            dDiscountFactorsOnBasisDates,
                            dSwapRate_,
                            dBasisAnnuity_,
                            dAnnuity_,
                            dJulianObservationDate);

        double dFunctionValueUp=GetOptionPrice(0.,dAdjSwapRate,dAdjAnnuity,dSwapRate,dFwdStdDevUp,(dNoticeDate-dJulianObservationDate)/365.)
                                +GetOptionPrice(1.,dAdjSwapRate,dAdjAnnuity,dSwapRate,dFwdStdDevUp,(dNoticeDate-dJulianObservationDate)/365.)
                                -dMarketPrice;

        double dRatio=0.0001*dFunctionValue/(dFunctionValueUp-dFunctionValue);
        dGuess-=dRatio;

        // Force -ve guess to v.small +ve
        if (dGuess < 0.0)
        {
            dGuess = 0.0000001;
        }

        if(fabs(dRatio)<2.e-05)
        {
            dVolStrip.at(k)=dGuess;
            break;
        }
    } // while iCounter < 100
    if(dVolStrip.at(k)<0.0001) throw("Zero or Negative Vol in bootstrapping");
    return dVolStrip.at(k);
} // Bootstrapping2F(...)


double  Bootstrapping3F(double dSwapStart,
                        double dSpotDate,
                        unsigned int k,                                       // notice index
                        const DKMaille2D<double> &dVolInfo,										// all info regarding vols
                        double dAbsVol,                                       // AbsVol
                        const DKMaille<double> &swaptionDates,                // notice dates from the fixed leg
                        const DKMaille<double> &accrualPeriods,               // accrual periods from fixed leg
                        const DKMaille<double> &basisSwaptionDates,           // notice dates from the float leg
                        const DKMaille<double> &basisAccrualPeriods,          // accrual periods from the float leg
                        const DKMaille<double> &basis,                        // basis from float leg
                        const DKMaille<double> &dModelParameters,
                        const DKMaille<double> &dDiscountDates,
                        const DKMaille<double> &dDiscountRates,
                        const DKMaille<double> &dAdjustedDiscountRates,
                        double dJulianObservationDate
                       )
{

    DKMaille<double> dVolStripDates(k+2);
    DKMaille<double> dVolStrip(k+2);

    dVolStripDates.at(0)=dSpotDate;

    // fill past vols and dates
    for(unsigned int ii=0;ii<k+1;ii++)
    {
        dVolStripDates.at(ii+1)=dVolInfo.at(ii,0);
        if(k>0)
        {
            dVolStrip.at(ii)=dVolInfo.at(ii,1);
        }
    } // for ii

    // Get Model parameters
    double dMeanReversion1 = dModelParameters.at(1);
    double dMeanReversion2 = dModelParameters.at(2);
    double dMeanReversion3 = dModelParameters.at(3);
    double dRelativeFactor12 = dModelParameters.at(4);
    double dRelativeFactor13 = dModelParameters.at(5);
    double dCorrelation12    = dModelParameters.at(6);
    double dCorrelation13		 = dModelParameters.at(7);
    double dCorrelation23		 = dModelParameters.at(8);


    // Get AdjSwap Rate
    double dAdjSwapRate=GetAdjSwapRate(dSpotDate,
                                       swaptionDates,
                                       accrualPeriods,
                                       basisSwaptionDates,
                                       basis,
                                       basisAccrualPeriods,
                                       dDiscountDates,
                                       dDiscountRates,
                                       dAdjustedDiscountRates);

    // Get AdjSwap Annuity
    double dAdjAnnuity=GetAdjAnnuity(dSpotDate,
                                     swaptionDates,
                                     accrualPeriods,
                                     dDiscountDates,
                                     dDiscountRates,
                                     dAdjustedDiscountRates,
                                     dJulianObservationDate);

    // Get Swap Rate
    double dSwapRate=GetSwapRate(dSpotDate,
                                 swaptionDates,
                                 accrualPeriods,
                                 dDiscountDates,
                                 dDiscountRates,
                                 dAdjustedDiscountRates);

    double dNoticeDate = dVolInfo.at(k,0);

    double dNoticePeriod = (swaptionDates.at(0)-dNoticeDate)/365.;

    // Catch on dates
    if(dNoticePeriod<0.)
    {
        throw("Notice Date After Swap Start");
    }

    // Get straddle price from AbsoluteVol
    double dMarketPrice=2.*GetPriceFromAbsVol(dSpotDate,
                        dJulianObservationDate,
                        swaptionDates,
                        accrualPeriods,
                        dNoticePeriod,
                        dAbsVol,
                        dDiscountDates,
                        dDiscountRates,
                        dAdjustedDiscountRates);

    double dGuess;

    if(0==k)
    {
        dGuess = dAbsVol;
    }
    else
    {
        dGuess = dVolStrip.at(k-1);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Look to see if there is NO solution i.e. we don't need Newton-Raphson
    /////////////////////////////////////////////////////////////////////////////
    dVolStrip.at(k)=0.0001;

    DKMaille<double> dDiscountFactors(swaptionDates.entries());
    DKMaille<double> dDiscountFactorsOnBasisDates(basisSwaptionDates.entries());
    double dSwapRate_=0.; double dAnnuity_=0.; double dBasisAnnuity_=0.;
    dIntermediateCalcSwap(dSpotDate,
                          dDiscountDates,
                          dDiscountRates,
                          dAdjustedDiscountRates,
                          swaptionDates,
                          basisSwaptionDates,
                          accrualPeriods,
                          basis,
                          basisAccrualPeriods,
                          dDiscountFactors,
                          dDiscountFactorsOnBasisDates,
                          &dSwapRate_,
                          &dAnnuity_,
                          &dBasisAnnuity_,
                          dJulianObservationDate);

    double dFwdStdDev = GetSwaptionAbsVol_VFDK_HW3F(dSpotDate,
                        swaptionDates,
                        accrualPeriods,
                        basisSwaptionDates,
                        basis,
                        basisAccrualPeriods,
                        dNoticePeriod,
                        dMeanReversion1,
                        dMeanReversion2,
                        dMeanReversion3,
                        dRelativeFactor12,
                        dRelativeFactor13,
                        dCorrelation12,
                        dCorrelation13,
                        dCorrelation23,
                        dVolStripDates,
                        dVolStrip,
                        dDiscountFactors,
                        dDiscountFactorsOnBasisDates,
                        dSwapRate_,
                        dBasisAnnuity_,
                        dAnnuity_,
                        dJulianObservationDate);


    double dFunctionValue = GetOptionPrice(0.,dAdjSwapRate,dAdjAnnuity,dSwapRate,dFwdStdDev,(dNoticeDate-dJulianObservationDate)/365.)
                            +GetOptionPrice(1.,dAdjSwapRate,dAdjAnnuity,dSwapRate,dFwdStdDev,(dNoticeDate-dJulianObservationDate)/365.)
                            -dMarketPrice;

    if(dFunctionValue >= 0.)
    {
        return 0.; // failure state - no solution ever (???)
    }

    /////////////////////////////////////////////////////////////////////////////
    // End of check pre-NR check ...
    /////////////////////////////////////////////////////////////////////////////

    // Newton-Raphson loop to get fwd std dev
    int iCounter=0;

    while(iCounter<100)
    {
        iCounter++;
        dVolStrip.at(k)=dGuess;

        // Function
        dFwdStdDev=GetSwaptionAbsVol_VFDK_HW3F(dSpotDate,
                                               swaptionDates,
                                               accrualPeriods,
                                               basisSwaptionDates,
                                               basis,
                                               basisAccrualPeriods,
                                               dNoticePeriod,
                                               dMeanReversion1,
                                               dMeanReversion2,
                                               dMeanReversion3,
                                               dRelativeFactor12,
                                               dRelativeFactor13,
                                               dCorrelation12,
                                               dCorrelation13,
                                               dCorrelation23,
                                               dVolStripDates,
                                               dVolStrip,
                                               dDiscountFactors,
                                               dDiscountFactorsOnBasisDates,
                                               dSwapRate_,
                                               dBasisAnnuity_,
                                               dAnnuity_,
                                               dJulianObservationDate);

        dFunctionValue=GetOptionPrice(0.,dAdjSwapRate,dAdjAnnuity,dSwapRate,dFwdStdDev,(dNoticeDate-dJulianObservationDate)/365.)
                       +GetOptionPrice(1.,dAdjSwapRate,dAdjAnnuity,dSwapRate,dFwdStdDev,(dNoticeDate-dJulianObservationDate)/365.)
                       -dMarketPrice;

        // FunctionUp
        dVolStrip.at(k)=dGuess+0.0001;

        double dFwdStdDevUp=GetSwaptionAbsVol_VFDK_HW3F(dSpotDate,
                            swaptionDates,
                            accrualPeriods,
                            basisSwaptionDates,
                            basis,
                            basisAccrualPeriods,
                            dNoticePeriod,
                            dMeanReversion1,
                            dMeanReversion2,
                            dMeanReversion3,
                            dRelativeFactor12,
                            dRelativeFactor13,
                            dCorrelation12,
                            dCorrelation13,
                            dCorrelation23,
                            dVolStripDates,
                            dVolStrip,
                            dDiscountFactors,
                            dDiscountFactorsOnBasisDates,
                            dSwapRate_,
                            dBasisAnnuity_,
                            dAnnuity_,
                            dJulianObservationDate);

        double dFunctionValueUp=GetOptionPrice(0.,dAdjSwapRate,dAdjAnnuity,dSwapRate,dFwdStdDevUp,(dNoticeDate-dJulianObservationDate)/365.)
                                +GetOptionPrice(1.,dAdjSwapRate,dAdjAnnuity,dSwapRate,dFwdStdDevUp,(dNoticeDate-dJulianObservationDate)/365.)
                                -dMarketPrice;

        double dRatio=0.0001*dFunctionValue/(dFunctionValueUp-dFunctionValue);
        dGuess-=dRatio;

        // Force -ve guess to v.small +ve
        if (dGuess < 0.0)
        {
            dGuess = 0.0000001;
        }

        if(fabs(dRatio)<2.e-05)
        {
            dVolStrip.at(k)=dGuess;
            break;
        }
    } // while iCounter < 100

    if(dVolStrip.at(k)<0.0001) throw("Zero or Negative Vol in bootstrapping");
    return dVolStrip.at(k);
} // Bootstrapping3F(...)

DKMaille<double> Bootstrapping_VFDK_HW1To3F(DKMaille2D<double> dAbsVol,
        DKMaille<double> dSwaptionExpiry,
        DKMaille<double> dSwaptionTenor,
        DKMaille<double> dDiscountDates,
        DKMaille<double> dDiscountRates,
        DKMaille<double> dAdjustedDiscountRates,
        DKMaille<double> dNoticeDates,
        DKMaille<double> dSwapStartDates,
        DKMaille<double> dSwapEndDates,
        DKMaille<double> dModelParameters,
        double dJulianObservationDate)
{

    unsigned int uiFixedIncomeFactors = (unsigned int)dModelParameters.at(0);

    if(dSwapStartDates.entries()!=dNoticeDates.entries()-2)
        throw("Wrong Input Dates");

    if(dSwapStartDates.entries()!=dSwapEndDates.entries())
        throw("Wrong Input Dates");

    static DKMaille<double> vol(dNoticeDates.entries()-2);
    vol.resize(dNoticeDates.entries()-2);
    // for all notice dates
    DKMaille<double> swaptionDates;      // swap start date, then all n payment dates of the fixed leg.
    DKMaille<double> basisSwaptionDates; // all m payment dates of the floating leg.
    DKMaille<double> accrualPeriods;     // bond basis adjusted period length of fixed leg in basis years (model years) array of size n+1 first element is zero.
    DKMaille<double> basisAccrualPeriods;// bond basis adjusted period length of floating leg in basis years (model years)
    DKMaille<double> basis;              // basis adjustments

    // Transpose dAbsVol
    // Update the function when more time
    DKMaille2D<double> dInputSurfaceZArray(dAbsVol.columns(),dAbsVol.rows());
    for(int ui=0;ui<dAbsVol.rows();ui++)
        for(int uk=0;uk<dAbsVol.columns();uk++)
            dInputSurfaceZArray.at(uk,ui) = dAbsVol.at(ui,uk);

    double dSpotDate=dNoticeDates.at(0);
    if(dJulianObservationDate!=dSpotDate)
        throw("Observation Date Different From Spot Date in Boostrapping. Feature will not be supported");

    DKMaille2D<double> dVolInfo(dNoticeDates.entries(),2);

    double dVol=0.;

    for(unsigned int uiN=0;uiN<dNoticeDates.entries()-2;uiN++)
    {
        // Create the underlying vanilla swap for relevant notice date
        CreateVanillaSwap(dSpotDate,																						// SpotDate input
                          (dSwapEndDates.at(uiN)-dSwapStartDates.at(uiN))/365.,	// dSwapTenorInYears input
                          (dNoticeDates.at(uiN+1)-dSpotDate)/365.,			// dOptionExpiryInYears input
                          swaptionDates,
                          accrualPeriods,
                          basisSwaptionDates,
                          basisAccrualPeriods,
                          basis,
                          (dSwapStartDates.at(uiN)-dNoticeDates.at(uiN+1))/365., // NoticePeriodInDays
                          dDiscountDates,
                          dDiscountRates,
                          dAdjustedDiscountRates);

        // Interpolate Absolute Volatility
        double dAbsoluteVolFromSurface=InterpolateMatrix(dInputSurfaceZArray,
                                       (dSwapEndDates.at(uiN)-dSwapStartDates.at(uiN))/365.,
                                       (dNoticeDates.at(uiN+1)-dNoticeDates.at(0))/365.,
                                       dSwaptionTenor,
                                       dSwaptionExpiry);

        dVolInfo.at(uiN,0)=dNoticeDates.at(uiN+1);
        dVolInfo.at(uiN,1)=0.01;


        switch (uiFixedIncomeFactors)
        {
        case 3:
            {
                dVol  = Bootstrapping3F(dSwapStartDates.at(uiN),
                                        dSpotDate,
                                        uiN,																// notice index
                                        dVolInfo,														// vol info
                                        dAbsoluteVolFromSurface,            // ABSVOL
                                        swaptionDates,											// notice dates from the fixed leg
                                        accrualPeriods,											// accrual periods from fixed leg
                                        basisSwaptionDates,									// notice dates from the float leg
                                        basisAccrualPeriods,								// accrual periods from the float leg
                                        basis,															// basis from float leg
                                        dModelParameters,
                                        dDiscountDates,
                                        dDiscountRates,
                                        dAdjustedDiscountRates,
                                        dJulianObservationDate
                                       );
                vol.at(uiN) = dVol;
            }
            break;
        case 2:
            {
                dVol  = Bootstrapping2F(dSwapStartDates.at(uiN),
                                        dSpotDate,
                                        uiN,																// notice index
                                        dVolInfo,														// vol info
                                        dAbsoluteVolFromSurface,            // ABSVOL
                                        swaptionDates,											// notice dates from the fixed leg
                                        accrualPeriods,											// accrual periods from fixed leg
                                        basisSwaptionDates,									// notice dates from the float leg
                                        basisAccrualPeriods,								// accrual periods from the float leg
                                        basis,															// basis from float leg
                                        dModelParameters,
                                        dDiscountDates,
                                        dDiscountRates,
                                        dAdjustedDiscountRates,
                                        dJulianObservationDate
                                       );
                vol.at(uiN) = dVol;
            }
            break;
        case 1:
            {
                dVol  = Bootstrapping1F(dSwapStartDates.at(uiN),
                                        dSpotDate,
                                        uiN,																// notice index
                                        dVolInfo,														// vol information
                                        dAbsoluteVolFromSurface,            // ABSVOL
                                        swaptionDates,											// notice dates from the fixed leg
                                        accrualPeriods,											// accrual periods from fixed leg
                                        basisSwaptionDates,									// notice dates from the float leg
                                        basisAccrualPeriods,								// accrual periods from the float leg
                                        basis,															// basis from float leg
                                        dModelParameters,
                                        dDiscountDates,
                                        dDiscountRates,
                                        dAdjustedDiscountRates,
                                        dJulianObservationDate
                                       );
                // the
                vol.at(uiN) = dVol;
            } // case 1-factor
            break;
        default:
            {
                throw("You should never be here");
            }
        } // switch

        if(dVol<0.0001)
            throw("Zero or Negative Vol in bootstrapping");
        dVolInfo.at(uiN,1)=dVol;

    }
    return vol;
} // Bootstrapping_VFDK_HW1To3F()


void CalcShortRates(double dSpotDate,
                    DKMaille<double> &dDates,
                    DKMaille<double> &dRates,
                    DKMaille<double> &dAdjRates,
                    DKMaille<double> &dShortRateDates,
                    DKMaille<double> &dShortRates)
{
    DKMaille<double> dZCRates(dRates);
    DKMaille<double> dZCAdjRates(dAdjRates);
    DKMaille<double> t(dDates);
    FromDiscountToZero(dSpotDate,dZCRates,dZCAdjRates,dDates,dRates,dAdjRates);
    for(unsigned int uiK=0;uiK<dDates.entries();uiK++)
    {
        t.at(uiK)=(dDates.at(uiK)-dSpotDate)/365.;
    }
    double t1=0.;
    double tN=61.;
    double dt=0.25;
    while(t1<tN)
    {
        double t2=t1+dt;
        double d1=exp(-rateinterpolation_dk_maille(2,t1,t,dZCAdjRates,t.entries()-1)*t1);
        double d2=exp(-rateinterpolation_dk_maille(2,t2,t,dZCAdjRates,t.entries()-1)*t2);
        dShortRates.insert((d1/d2-1.)/dt);
        dShortRateDates.insert(t1);
        t1+=dt;
    }
}

void CalcXCCYShortRates(double dSpotDate,
                        DKMaille<double> &dDates,
                        DKMaille<double> &dRates,
                        DKMaille<double> &dAdjRates,
                        DKMaille<double> &dForeignDates,
                        DKMaille<double> &dForeignRates,
                        DKMaille<double> &dForeignAdjRates,
                        DKMaille<double> &dShortRateDates,
                        DKMaille<double> &dShortRates,
                        DKMaille<double> &dShortRatesForeign)
{
    DKMaille<double> dZCRates(dRates);
    DKMaille<double> dZCAdjRates(dAdjRates);
    DKMaille<double> dZCForeignRates(dForeignRates);
    DKMaille<double> dZCForeignAdjRates(dForeignAdjRates);
    DKMaille<double> t(dDates);
    DKMaille<double> tf(dForeignDates);
    FromDiscountToZero(dSpotDate,dZCRates,dZCAdjRates,dDates,dRates,dAdjRates);
    FromDiscountToZero(dSpotDate,dZCForeignRates,dZCForeignAdjRates,dForeignDates,dForeignRates,dForeignAdjRates);
    for(unsigned int uiK=0;uiK<dDates.entries();uiK++)
    {
        t.at(uiK)=(dDates.at(uiK)-dSpotDate)/365.;
    }
    for(uiK=0;uiK<dForeignDates.entries();uiK++)
    {
        tf.at(uiK)=(dForeignDates.at(uiK)-dSpotDate)/365.;
    }
    double t1=0.;
    double tN=61.;
    double dt=0.25;
    while(t1<tN)
    {
        double t2=t1+dt;
        double d1=exp(-rateinterpolation_dk_maille(2,t1,t,dZCAdjRates,t.entries()-1)*t1);
        double d1_=exp(-rateinterpolation_dk_maille(2,t1,tf,dZCForeignAdjRates,tf.entries()-1)*t1);
        double d2=exp(-rateinterpolation_dk_maille(2,t2,t,dZCAdjRates,t.entries()-1)*t2);
        double d2_=exp(-rateinterpolation_dk_maille(2,t2,tf,dZCForeignAdjRates,tf.entries()-1)*t2);
        dShortRates.insert((d1/d2-1.)/dt);
        dShortRatesForeign.insert((d1_/d2_-1.)/dt);
        // dShortRates.insert(1.);
        // dShortRatesForeign.insert(1.);
        dShortRateDates.insert(t1);
        t1+=dt;
    }
}





double SwaptionPrice_VFDK_HW1To3F(double dSwaptionExpiryInYears,
                                  double dSwaptionTenorInYears,
                                  double dNoticePeriodInDays,
                                  double dStrike,
                                  double dCallPut,
                                  DKMaille<double> dDiscountCurveDates,
                                  DKMaille<double> dDiscountCurveRates,
                                  DKMaille<double> dAdjustedDiscountCurveRates,
                                  DKMaille<double> dNoticeDates,
                                  DKMaille<double> dTimeDependentStandardDeviation,
                                  DKMaille<double> dModelParameters,
                                  double dJulianObservationDate)
{

    unsigned int uiFixedIncomeFactors = (unsigned int)dModelParameters.at(0);

    if(dTimeDependentStandardDeviation.at(dTimeDependentStandardDeviation.entries()-1)!=
            dTimeDependentStandardDeviation.at(dTimeDependentStandardDeviation.entries()-3))
        throw("Wrong Input Model Vols");

    if(dTimeDependentStandardDeviation.at(dTimeDependentStandardDeviation.entries()-1)!=
            dTimeDependentStandardDeviation.at(dTimeDependentStandardDeviation.entries()-2))
        throw("Wrong Input Model Vols");

    if(dTimeDependentStandardDeviation.entries()!=dNoticeDates.entries())
        throw("Wrong Input Dates in Model Vols");

    // for all notice dates
    DKMaille<double> swaptionDates;      // swap start date, then all n payment dates of the fixed leg.
    DKMaille<double> basisSwaptionDates; // all m payment dates of the floating leg.
    DKMaille<double> accrualPeriods;     // bond basis adjusted period length of fixed leg in basis years (model years) array of size n+1 first element is zero.
    DKMaille<double> basisAccrualPeriods;// bond basis adjusted period length of floating leg in basis years (model years)
    DKMaille<double> basis;              // basis adjustments

    double dSpotDate=dNoticeDates.at(0);

    DKMaille2D<double> dVolInfo(dNoticeDates.entries(),2);

    // Create the underlying vanilla swap for relevant notice date
    CreateVanillaSwap(dSpotDate,									// SpotDate input
                      dSwaptionTenorInYears,			// dSwapTenorInYears input
                      dSwaptionExpiryInYears+(dJulianObservationDate-dSpotDate)/365.,			// swaption ExpiryInYears corresponds always to the real option expiry
                      swaptionDates,
                      accrualPeriods,
                      basisSwaptionDates,
                      basisAccrualPeriods,
                      basis,
                      dNoticePeriodInDays/365.,		// NoticePeriodInYears
                      dDiscountCurveDates,
                      dDiscountCurveRates,
                      dAdjustedDiscountCurveRates);

    // Get AdjSwap Rate
    double dAdjSwapRate=GetAdjSwapRate(dSpotDate,
                                       swaptionDates,
                                       accrualPeriods,
                                       basisSwaptionDates,
                                       basis,
                                       basisAccrualPeriods,
                                       dDiscountCurveDates,
                                       dDiscountCurveRates,
                                       dAdjustedDiscountCurveRates);

    // Get AdjSwap Annuity
    double dAdjAnnuity=GetAdjAnnuity(dSpotDate,
                                     swaptionDates,
                                     accrualPeriods,
                                     dDiscountCurveDates,
                                     dDiscountCurveRates,
                                     dAdjustedDiscountCurveRates,
                                     dJulianObservationDate);

    // Get Swap Rate
    double dSwapRate=GetSwapRate(dSpotDate,
                                 swaptionDates,
                                 accrualPeriods,
                                 dDiscountCurveDates,
                                 dDiscountCurveRates,
                                 dAdjustedDiscountCurveRates);


    double dFwdStdDev=0.;

    // Get Absolute Volatility
    switch (uiFixedIncomeFactors)
    {
    case 2:
        {

            // Get Model parameters
            double dMeanReversion1 = dModelParameters.at(1);
            double dMeanReversion2 = dModelParameters.at(2);
            double dRelativeFactor = dModelParameters.at(4);
            double dCorrelation    = dModelParameters.at(6);
            DKMaille<double> dDiscountFactors(swaptionDates.entries());
            DKMaille<double> dDiscountFactorsOnBasisDates(basisSwaptionDates.entries());
            double dSwapRate_=0.; double dAnnuity_=0.; double dBasisAnnuity_=0.;
            dIntermediateCalcSwap(dSpotDate,
                                  dDiscountCurveDates,
                                  dDiscountCurveRates,
                                  dAdjustedDiscountCurveRates,
                                  swaptionDates,
                                  basisSwaptionDates,
                                  accrualPeriods,
                                  basis,
                                  basisAccrualPeriods,
                                  dDiscountFactors,
                                  dDiscountFactorsOnBasisDates,
                                  &dSwapRate_,
                                  &dAnnuity_,
                                  &dBasisAnnuity_,
                                  dJulianObservationDate);

            dFwdStdDev=GetSwaptionAbsVol_VFDK_HW2F(dSpotDate,
                                                   swaptionDates,
                                                   accrualPeriods,
                                                   basisSwaptionDates,
                                                   basis,
                                                   basisAccrualPeriods,
                                                   dNoticePeriodInDays/365.,
                                                   dMeanReversion1,
                                                   dMeanReversion2,
                                                   dRelativeFactor,
                                                   dCorrelation,
                                                   dNoticeDates,
                                                   dTimeDependentStandardDeviation,
                                                   dDiscountFactors,
                                                   dDiscountFactorsOnBasisDates,
                                                   dSwapRate_,
                                                   dBasisAnnuity_,
                                                   dAnnuity_,
                                                   dJulianObservationDate);

        } // case 2-factor
        break;
    case 1:
        {
            // Get Model parameters
            double dMeanReversion1 = dModelParameters.at(1);
            DKMaille<double> dDiscountFactors(swaptionDates.entries());
            DKMaille<double> dDiscountFactorsOnBasisDates(basisSwaptionDates.entries());
            double dSwapRate_=0.; double dAnnuity_=0.; double dBasisAnnuity_=0.;
            dIntermediateCalcSwap(dSpotDate,
                                  dDiscountCurveDates,
                                  dDiscountCurveRates,
                                  dAdjustedDiscountCurveRates,
                                  swaptionDates,
                                  basisSwaptionDates,
                                  accrualPeriods,
                                  basis,
                                  basisAccrualPeriods,
                                  dDiscountFactors,
                                  dDiscountFactorsOnBasisDates,
                                  &dSwapRate_,
                                  &dAnnuity_,
                                  &dBasisAnnuity_,
                                  dJulianObservationDate);

            dFwdStdDev=GetSwaptionAbsVol_VFDK_HW1F(dSpotDate,
                                                   swaptionDates,
                                                   accrualPeriods,
                                                   basisSwaptionDates,
                                                   basis,
                                                   basisAccrualPeriods,
                                                   dNoticePeriodInDays/365.,
                                                   dMeanReversion1,
                                                   dNoticeDates,
                                                   dTimeDependentStandardDeviation,
                                                   dDiscountFactors,
                                                   dDiscountFactorsOnBasisDates,
                                                   dSwapRate_,
                                                   dBasisAnnuity_,
                                                   dAnnuity_,
                                                   dJulianObservationDate);

        } // case 1-factor
        break;
    case 3:
        {
            // Get Model parameters
            double dMeanReversion1 = dModelParameters.at(1);
            double dMeanReversion2 = dModelParameters.at(2);
            double dMeanReversion3 = dModelParameters.at(3);

            double dRelativeFactor12 = dModelParameters.at(4);
            double dRelativeFactor13 = dModelParameters.at(5);

            double dCorrelation12    = dModelParameters.at(6);
            double dCorrelation13    = dModelParameters.at(7);
            double dCorrelation23    = dModelParameters.at(8);

            DKMaille<double> dDiscountFactors(swaptionDates.entries());
            DKMaille<double> dDiscountFactorsOnBasisDates(basisSwaptionDates.entries());
            double dSwapRate_=0.; double dAnnuity_=0.; double dBasisAnnuity_=0.;
            dIntermediateCalcSwap(dSpotDate,
                                  dDiscountCurveDates,
                                  dDiscountCurveRates,
                                  dAdjustedDiscountCurveRates,
                                  swaptionDates,
                                  basisSwaptionDates,
                                  accrualPeriods,
                                  basis,
                                  basisAccrualPeriods,
                                  dDiscountFactors,
                                  dDiscountFactorsOnBasisDates,
                                  &dSwapRate_,
                                  &dAnnuity_,
                                  &dBasisAnnuity_,
                                  dJulianObservationDate);

            dFwdStdDev=GetSwaptionAbsVol_VFDK_HW3F(dSpotDate,
                                                   swaptionDates,
                                                   accrualPeriods,
                                                   basisSwaptionDates,
                                                   basis,
                                                   basisAccrualPeriods,
                                                   dNoticePeriodInDays/365.,
                                                   dMeanReversion1,
                                                   dMeanReversion2,
                                                   dMeanReversion3,
                                                   dRelativeFactor12,
                                                   dRelativeFactor13,
                                                   dCorrelation12,
                                                   dCorrelation13,
                                                   dCorrelation23,
                                                   dNoticeDates,
                                                   dTimeDependentStandardDeviation,
                                                   dDiscountFactors,
                                                   dDiscountFactorsOnBasisDates,
                                                   dSwapRate_,
                                                   dBasisAnnuity_,
                                                   dAnnuity_,
                                                   dJulianObservationDate);
        }
        break;
    } // switch

    double dPrice=GetOptionPrice(dCallPut,dAdjSwapRate,dAdjAnnuity,dStrike,dFwdStdDev,dSwaptionExpiryInYears);

    return dPrice;

} // SwaptionPrice_VFDK_HW1To3F()




double SwaptionAbsVol_VFDK_HW1To3F(double dSwaptionExpiryInYears,
                                   double dSwaptionTenorInYears,
                                   double dNoticePeriodInDays,
                                   DKMaille<double> dDiscountCurveDates,
                                   DKMaille<double> dDiscountCurveRates,
                                   DKMaille<double> dAdjustedDiscountCurveRates,
                                   DKMaille<double> dNoticeDates,
                                   DKMaille<double> dTimeDependentStandardDeviation,
                                   DKMaille<double> dModelParameters,
                                   double dJulianObservationDate)
{

    unsigned int uiFixedIncomeFactors = (unsigned int)dModelParameters.at(0);

    if(dTimeDependentStandardDeviation.at(dTimeDependentStandardDeviation.entries()-1)!=
            dTimeDependentStandardDeviation.at(dTimeDependentStandardDeviation.entries()-3))
        throw("Wrong Input Model Vols");

    if(dTimeDependentStandardDeviation.at(dTimeDependentStandardDeviation.entries()-1)!=
            dTimeDependentStandardDeviation.at(dTimeDependentStandardDeviation.entries()-2))
        throw("Wrong Input Model Vols");

    if(dTimeDependentStandardDeviation.entries()!=dNoticeDates.entries())
        throw("Wrong Input Dates in Model Vols");

    // for all notice dates
    DKMaille<double> swaptionDates;      // swap start date, then all n payment dates of the fixed leg.
    DKMaille<double> basisSwaptionDates; // all m payment dates of the floating leg.
    DKMaille<double> accrualPeriods;     // bond basis adjusted period length of fixed leg in basis years (model years) array of size n+1 first element is zero.
    DKMaille<double> basisAccrualPeriods;// bond basis adjusted period length of floating leg in basis years (model years)
    DKMaille<double> basis;              // basis adjustments

    double dSpotDate=dNoticeDates.at(0);

    DKMaille2D<double> dVolInfo(dNoticeDates.entries(),2);

    // Create the underlying vanilla swap for relevant notice date
    CreateVanillaSwap(dSpotDate,									// SpotDate input
                      dSwaptionTenorInYears,			// dSwapTenorInYears input
                      dSwaptionExpiryInYears+(dJulianObservationDate-dSpotDate)/365.,			// swaption ExpiryInYears corresponds always to the real option expiry
                      swaptionDates,
                      accrualPeriods,
                      basisSwaptionDates,
                      basisAccrualPeriods,
                      basis,
                      dNoticePeriodInDays/365.,		// NoticePeriodInYears
                      dDiscountCurveDates,
                      dDiscountCurveRates,
                      dAdjustedDiscountCurveRates);

    // Get AdjSwap Rate
    double dAdjSwapRate=GetAdjSwapRate(dSpotDate,
                                       swaptionDates,
                                       accrualPeriods,
                                       basisSwaptionDates,
                                       basis,
                                       basisAccrualPeriods,
                                       dDiscountCurveDates,
                                       dDiscountCurveRates,
                                       dAdjustedDiscountCurveRates);

    // Get AdjSwap Annuity
    double dAdjAnnuity=GetAdjAnnuity(dSpotDate,
                                     swaptionDates,
                                     accrualPeriods,
                                     dDiscountCurveDates,
                                     dDiscountCurveRates,
                                     dAdjustedDiscountCurveRates,
                                     dJulianObservationDate);

    // Get Swap Rate
    double dSwapRate=GetSwapRate(dSpotDate,
                                 swaptionDates,
                                 accrualPeriods,
                                 dDiscountCurveDates,
                                 dDiscountCurveRates,
                                 dAdjustedDiscountCurveRates);


    double dFwdStdDev=0.;

    // Get Absolute Volatility
    switch (uiFixedIncomeFactors)
    {
    case 2:
        {

            // Get Model parameters
            double dMeanReversion1 = dModelParameters.at(1);
            double dMeanReversion2 = dModelParameters.at(2);
            double dRelativeFactor = dModelParameters.at(4);
            double dCorrelation    = dModelParameters.at(6);
            DKMaille<double> dDiscountFactors(swaptionDates.entries());
            DKMaille<double> dDiscountFactorsOnBasisDates(basisSwaptionDates.entries());
            double dSwapRate_=0.; double dAnnuity_=0.; double dBasisAnnuity_=0.;
            dIntermediateCalcSwap(dSpotDate,
                                  dDiscountCurveDates,
                                  dDiscountCurveRates,
                                  dAdjustedDiscountCurveRates,
                                  swaptionDates,
                                  basisSwaptionDates,
                                  accrualPeriods,
                                  basis,
                                  basisAccrualPeriods,
                                  dDiscountFactors,
                                  dDiscountFactorsOnBasisDates,
                                  &dSwapRate_,
                                  &dAnnuity_,
                                  &dBasisAnnuity_,
                                  dJulianObservationDate);

            dFwdStdDev=GetSwaptionAbsVol_VFDK_HW2F(dSpotDate,
                                                   swaptionDates,
                                                   accrualPeriods,
                                                   basisSwaptionDates,
                                                   basis,
                                                   basisAccrualPeriods,
                                                   dNoticePeriodInDays/365.,
                                                   dMeanReversion1,
                                                   dMeanReversion2,
                                                   dRelativeFactor,
                                                   dCorrelation,
                                                   dNoticeDates,
                                                   dTimeDependentStandardDeviation,
                                                   dDiscountFactors,
                                                   dDiscountFactorsOnBasisDates,
                                                   dSwapRate_,
                                                   dBasisAnnuity_,
                                                   dAnnuity_,
                                                   dJulianObservationDate);

        } // case 2-factor
        break;
    case 1:
        {
            // Get Model parameters
            double dMeanReversion1 = dModelParameters.at(1);
            DKMaille<double> dDiscountFactors(swaptionDates.entries());
            DKMaille<double> dDiscountFactorsOnBasisDates(basisSwaptionDates.entries());
            double dSwapRate_=0.; double dAnnuity_=0.; double dBasisAnnuity_=0.;
            dIntermediateCalcSwap(dSpotDate,
                                  dDiscountCurveDates,
                                  dDiscountCurveRates,
                                  dAdjustedDiscountCurveRates,
                                  swaptionDates,
                                  basisSwaptionDates,
                                  accrualPeriods,
                                  basis,
                                  basisAccrualPeriods,
                                  dDiscountFactors,
                                  dDiscountFactorsOnBasisDates,
                                  &dSwapRate_,
                                  &dAnnuity_,
                                  &dBasisAnnuity_,
                                  dJulianObservationDate);

            dFwdStdDev=GetSwaptionAbsVol_VFDK_HW1F(dSpotDate,
                                                   swaptionDates,
                                                   accrualPeriods,
                                                   basisSwaptionDates,
                                                   basis,
                                                   basisAccrualPeriods,
                                                   dNoticePeriodInDays/365.,
                                                   dMeanReversion1,
                                                   dNoticeDates,
                                                   dTimeDependentStandardDeviation,
                                                   dDiscountFactors,
                                                   dDiscountFactorsOnBasisDates,
                                                   dSwapRate_,
                                                   dBasisAnnuity_,
                                                   dAnnuity_,
                                                   dJulianObservationDate);

        } // case 1-factor
        break;
    case 3:
        {
            // Get Model parameters
            double dMeanReversion1 = dModelParameters.at(1);
            double dMeanReversion2 = dModelParameters.at(2);
            double dMeanReversion3 = dModelParameters.at(3);

            double dRelativeFactor12 = dModelParameters.at(4);
            double dRelativeFactor13 = dModelParameters.at(5);

            double dCorrelation12    = dModelParameters.at(6);
            double dCorrelation13    = dModelParameters.at(7);
            double dCorrelation23    = dModelParameters.at(8);

            DKMaille<double> dDiscountFactors(swaptionDates.entries());
            DKMaille<double> dDiscountFactorsOnBasisDates(basisSwaptionDates.entries());
            double dSwapRate_=0.; double dAnnuity_=0.; double dBasisAnnuity_=0.;
            dIntermediateCalcSwap(dSpotDate,
                                  dDiscountCurveDates,
                                  dDiscountCurveRates,
                                  dAdjustedDiscountCurveRates,
                                  swaptionDates,
                                  basisSwaptionDates,
                                  accrualPeriods,
                                  basis,
                                  basisAccrualPeriods,
                                  dDiscountFactors,
                                  dDiscountFactorsOnBasisDates,
                                  &dSwapRate_,
                                  &dAnnuity_,
                                  &dBasisAnnuity_,
                                  dJulianObservationDate);

            dFwdStdDev=GetSwaptionAbsVol_VFDK_HW3F(dSpotDate,
                                                   swaptionDates,
                                                   accrualPeriods,
                                                   basisSwaptionDates,
                                                   basis,
                                                   basisAccrualPeriods,
                                                   dNoticePeriodInDays/365.,
                                                   dMeanReversion1,
                                                   dMeanReversion2,
                                                   dMeanReversion3,
                                                   dRelativeFactor12,
                                                   dRelativeFactor13,
                                                   dCorrelation12,
                                                   dCorrelation13,
                                                   dCorrelation23,
                                                   dNoticeDates,
                                                   dTimeDependentStandardDeviation,
                                                   dDiscountFactors,
                                                   dDiscountFactorsOnBasisDates,
                                                   dSwapRate_,
                                                   dBasisAnnuity_,
                                                   dAnnuity_,
                                                   dJulianObservationDate);
        }
        break;
    } // switch

    return dFwdStdDev;

} // SwaptionAbsVol_VFDK_VFDK_HW1To3F()



double dIntegratedCorrelation_VFDK_HW2F_TV(double dSpotDate,
        DKMaille<double> dDiscountDates,
        const DKMaille<double> &dDiscountRates,
        const DKMaille<double> &dAdjustedDiscountRates,
        DKMaille<double> dSwaptionDates,
        const DKMaille<double> &dAccrualPeriods,
        DKMaille<double> dBasisSwaptionDates,
        DKMaille<double> dBasis,
        const DKMaille<double> &dBasisAccrualPeriods,
        DKMaille<double> dSwaptionDates2,
        const DKMaille<double> &dAccrualPeriods2,
        DKMaille<double> dBasisSwaptionDates2,
        DKMaille<double> dBasis2,
        const DKMaille<double> &dBasisAccrualPeriods2,
        double dNoticePeriod,
        double dMeanReversion1,
        double dMeanReversion2,
        double dRelativeFactor,
        double dCorrelation,
        DKMaille<double> dVolStripDates,
        DKMaille<double> &dVolStrip,
        double dJulianObservationDate)
{
    DKMaille<double> dZCRates(dDiscountDates.entries());
    DKMaille<double> dZCAdjRates(dDiscountDates.entries());
    FromDiscountToZero(dSpotDate,dZCRates,dZCAdjRates,dDiscountDates,dDiscountRates,dAdjustedDiscountRates);

    double dPrice=0.;

    unsigned int ui=0;
    for(ui=0;ui<dDiscountDates.entries(); ui++)
        dDiscountDates[ui]=(dDiscountDates[ui]-dSpotDate)/365.;

    for(ui=0;ui<dSwaptionDates.entries(); ui++)
        dSwaptionDates[ui]=(dSwaptionDates[ui]-dSpotDate)/365.;

    for(ui=0;ui<dSwaptionDates2.entries(); ui++)
        dSwaptionDates2[ui]=(dSwaptionDates2[ui]-dSpotDate)/365.;

    for(ui=0;ui<dBasisSwaptionDates.entries(); ui++)
        dBasisSwaptionDates[ui]=(dBasisSwaptionDates[ui]-dSpotDate)/365.;

    for(ui=0;ui<dBasisSwaptionDates2.entries(); ui++)
        dBasisSwaptionDates2[ui]=(dBasisSwaptionDates2[ui]-dSpotDate)/365.;

    for(ui=0;ui<dVolStripDates.entries(); ui++)
        dVolStripDates[ui]=(dVolStripDates[ui]-dSpotDate)/365.;

    DKMaille<double> dDiscountFactors(dSwaptionDates.entries());
    DKMaille<double> dDiscountFactorsOnBasisDates(dBasisSwaptionDates.entries());
    DKMaille<double> dDiscountFactors2(dSwaptionDates2.entries());
    DKMaille<double> dDiscountFactorsOnBasisDates2(dBasisSwaptionDates2.entries());

    // Build SwapRate1 analytics
    double dAnnuity = 0.;
    double dSwapRate = 0.;
    for(ui=0;ui<dSwaptionDates.entries();ui++)
    {
        dDiscountFactors[ui]=
            exp(-rateinterpolation_dk_maille(2,dSwaptionDates[ui],dDiscountDates,dZCRates,dDiscountDates.entries()-1)*dSwaptionDates[ui]);
        if(ui>0) dAnnuity+=dAccrualPeriods[ui]*dDiscountFactors[ui];
    }
    dSwapRate = (dDiscountFactors[0]-dDiscountFactors[dSwaptionDates.entries()-1])/dAnnuity;
    double dBasisAnnuity=0.;
    for(ui=0;ui<dBasisSwaptionDates.entries();ui++)
    {
        dDiscountFactorsOnBasisDates[ui]=
            exp(-rateinterpolation_dk_maille(2,dBasisSwaptionDates[ui],dDiscountDates,
                                             dZCRates,dDiscountDates.entries()-1)*dBasisSwaptionDates[ui]);
        dBasisAnnuity+=dBasis[ui]*dBasisAccrualPeriods[ui]*dDiscountFactorsOnBasisDates[ui];
    }
    dSwapRate+=(dBasisAnnuity/dAnnuity);

    // Build SwapRate2 analytics
    double dAnnuity2 = 0.;
    double dSwapRate2 = 0.;
    for(ui=0;ui<dSwaptionDates2.entries();ui++)
    {
        dDiscountFactors2[ui]=
            exp(-rateinterpolation_dk_maille(2,dSwaptionDates2[ui],dDiscountDates,dZCRates,dDiscountDates.entries()-1)*dSwaptionDates2[ui]);
        if(ui>0) dAnnuity2+=dAccrualPeriods2[ui]*dDiscountFactors2[ui];
    }
    dSwapRate2 = (dDiscountFactors2[0]-dDiscountFactors2[dSwaptionDates2.entries()-1])/dAnnuity2;
    double dBasisAnnuity2=0.;
    for(ui=0;ui<dBasisSwaptionDates2.entries();ui++)
    {
        dDiscountFactorsOnBasisDates2[ui]=
            exp(-rateinterpolation_dk_maille(2,dBasisSwaptionDates2[ui],dDiscountDates,
                                             dZCRates,dDiscountDates.entries()-1)*dBasisSwaptionDates2[ui]);
        dBasisAnnuity2+=dBasis2[ui]*dBasisAccrualPeriods2[ui]*dDiscountFactorsOnBasisDates2[ui];
    }
    dSwapRate2+=(dBasisAnnuity2/dAnnuity2);


    double dObservationDate=(dJulianObservationDate-dSpotDate)/365.;
    double dTerminalDate    = dSwaptionDates[0]-dNoticePeriod;
    if(dSwaptionDates[0]!=dSwaptionDates2[0]) throw("Start dates of underlying swaps are inconsistent. System does not yet support calculation of correlation between rates having different reset date.");

    double dPi = 3.141592653589790 ;

    DKMaille<double> dVector(dSwaptionDates.entries()+dBasisSwaptionDates.entries());
    DKMaille<double> dDiscounts(dDiscountFactors.entries()+dDiscountFactorsOnBasisDates.entries());
    DKMaille<double> dCash(dDiscountFactors.entries()+dDiscountFactorsOnBasisDates.entries());
    DKMaille<double> dDates(dSwaptionDates.entries()+dBasisSwaptionDates.entries());
    if(
        (dSwaptionDates.entries()!=dDiscountFactors.entries())
        ||(dSwaptionDates.entries()!=dAccrualPeriods.entries())
        ||(dBasisSwaptionDates.entries()!=dBasisAccrualPeriods.entries())
        ||(dBasisSwaptionDates.entries()!=dDiscountFactorsOnBasisDates.entries())
        ||(dBasisSwaptionDates.entries()!=dBasis.entries())
    )
        throw("Input deal dates are inconsistent");

    DKMaille<double> dVector2(dSwaptionDates2.entries()+dBasisSwaptionDates2.entries());
    DKMaille<double> dDiscounts2(dDiscountFactors2.entries()+dDiscountFactorsOnBasisDates2.entries());
    DKMaille<double> dCash2(dDiscountFactors2.entries()+dDiscountFactorsOnBasisDates2.entries());
    DKMaille<double> dDates2(dSwaptionDates2.entries()+dBasisSwaptionDates2.entries());
    if(
        (dSwaptionDates2.entries()!=dDiscountFactors2.entries())
        ||(dSwaptionDates2.entries()!=dAccrualPeriods2.entries())
        ||(dBasisSwaptionDates2.entries()!=dBasisAccrualPeriods2.entries())
        ||(dBasisSwaptionDates2.entries()!=dDiscountFactorsOnBasisDates2.entries())
        ||(dBasisSwaptionDates2.entries()!=dBasis2.entries())
    )
        throw("Input deal dates are inconsistent");

    // Build ItoVector for rate 1
    for(ui=0;ui<dVector.entries();ui++)
    {
        if(ui<dAccrualPeriods.entries())
        {
            dCash[ui]=dAccrualPeriods[ui]*dSwapRate;
            dDiscounts[ui]=dDiscountFactors[ui];
            dDates[ui]=dSwaptionDates[ui];
        }
        else
        {
            dCash[ui]=dBasisAccrualPeriods[ui-dAccrualPeriods.entries()]*dBasis[ui-dAccrualPeriods.entries()];
            dDiscounts[ui]=dDiscountFactorsOnBasisDates[ui-dAccrualPeriods.entries()];
            dDates[ui]=dBasisSwaptionDates[ui-dAccrualPeriods.entries()];
        }
    }
    dVector=GetItoVector(dSwaptionDates.entries(),dBasisSwaptionDates.entries(),dCash,dDiscounts);

    // Build ItoVector for rate 2
    for(ui=0;ui<dVector2.entries();ui++)
    {
        if(ui<dAccrualPeriods2.entries())
        {
            dCash2[ui]=dAccrualPeriods2[ui]*dSwapRate2;
            dDiscounts2[ui]=dDiscountFactors2[ui];
            dDates2[ui]=dSwaptionDates2[ui];
        }
        else
        {
            dCash2[ui]=dBasisAccrualPeriods2[ui-dAccrualPeriods2.entries()]*dBasis2[ui-dAccrualPeriods2.entries()];
            dDiscounts2[ui]=dDiscountFactorsOnBasisDates2[ui-dAccrualPeriods2.entries()];
            dDates2[ui]=dBasisSwaptionDates2[ui-dAccrualPeriods2.entries()];
        }
    }
    dVector2=GetItoVector(dSwaptionDates2.entries(),dBasisSwaptionDates2.entries(),dCash2,dDiscounts2);





    // New functionality
    // Make sure integration slices are flexible to allow for fast pricing of callable reverse floaters and callable CMS spread options
    // Checks on dates
    DKMaille<double> dIntegrationStrip=dCheckDates(dVolStripDates,dObservationDate,dTerminalDate);
    DKMaille<double> dIntegrationStripVols=dNewVols(dVolStripDates,dVolStrip,dIntegrationStrip);


    // Date tests
    bool bCheckDis=true;
    unsigned int uiM=0;
    for(;uiM<dIntegrationStrip.entries();uiM++)
    {
        if(fabs(dIntegrationStrip[uiM]-dTerminalDate)<0.00001)
            bCheckDis=false;
    }
    if(bCheckDis)
        throw("Dates of Std Dev Strip Are Inconsistent with Option Expiry");



    double dForwardStdDev=0.;
    double dIntermediateCalc=0.;
    double dForwardStdDev2=0.;
    double dIntermediateCalc2=0.;
    double dForwardStdDev3=0.;
    double dIntermediateCalc3=0.;



    // This is the core calculation of the VFDK approximation
    // Create second index
    DKMaille<double> dVolStrip2(dIntegrationStrip.entries());
    unsigned int uiK = 0;
    for(uiK=0;uiK<dVolStrip2.entries();uiK++) dVolStrip2[uiK]=dRelativeFactor*dIntegrationStripVols[uiK];


    // iterate through Kac space
    // calculation of absolute standard deviation rate 1
    for(uiK=0;uiK<dDates.entries();uiK++)
    {
        for(unsigned int uiM=0;uiM<dDates.entries();uiM++)
        {
            dIntermediateCalc=dVector[uiK]*dVector[uiM]*
                              (
                                  (1./(dMeanReversion1*dMeanReversion1))*
                                  GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion1,
                                                            dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dIntegrationStripVols,dDates[uiK],dDates[uiM])
                                  +(1./(dMeanReversion2*dMeanReversion2))*
                                  GetItoDeterminant_HW2F_TV(dMeanReversion2,dMeanReversion2,
                                                            dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip2,dVolStrip2,dDates[uiK],dDates[uiM])
                                  +2.*dCorrelation*
                                  (1./(dMeanReversion1*dMeanReversion2))*
                                  GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion2,
                                                            dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dVolStrip2,dDates[uiK],dDates[uiM])
                              );
            dForwardStdDev += dIntermediateCalc;
        }
    }
    dForwardStdDev = (1./dAnnuity) * sqrt(dForwardStdDev/(dTerminalDate-dObservationDate));

    // calculation of absolute standard deviation rate 2
    for(uiK=0;uiK<dDates2.entries();uiK++)
    {
        for(unsigned int uiM=0;uiM<dDates2.entries();uiM++)
        {
            dIntermediateCalc2=dVector2[uiK]*dVector2[uiM]*
                               (
                                   (1./(dMeanReversion1*dMeanReversion1))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion1,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dIntegrationStripVols,dDates2[uiK],dDates2[uiM])
                                   +(1./(dMeanReversion2*dMeanReversion2))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion2,dMeanReversion2,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip2,dVolStrip2,dDates2[uiK],dDates2[uiM])
                                   +2.*dCorrelation*
                                   (1./(dMeanReversion1*dMeanReversion2))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion2,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dVolStrip2,dDates2[uiK],dDates2[uiM])
                               );
            dForwardStdDev2 += dIntermediateCalc2;
        }
    }
    dForwardStdDev2 = (1./dAnnuity2) * sqrt(dForwardStdDev2/(dTerminalDate-dObservationDate));

    // calculation of mixed term in the variance equation
    for(uiK=0;uiK<dDates.entries();uiK++)
    {
        for(unsigned int uiM=0;uiM<dDates2.entries();uiM++)
        {
            dIntermediateCalc3=dVector[uiK]*dVector2[uiM]*
                               (
                                   (1./(dMeanReversion1*dMeanReversion1))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion1,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dIntegrationStripVols,dDates[uiK],dDates2[uiM])
                                   +(1./(dMeanReversion2*dMeanReversion2))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion2,dMeanReversion2,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip2,dVolStrip2,dDates[uiK],dDates2[uiM])
                                   +dCorrelation*
                                   (1./(dMeanReversion1*dMeanReversion2))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion2,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dVolStrip2,dDates[uiK],dDates2[uiM])
                                   +dCorrelation*
                                   (1./(dMeanReversion2*dMeanReversion1))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion2,dMeanReversion1,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip2,dIntegrationStripVols,dDates[uiK],dDates2[uiM])
                               );
            dForwardStdDev3 += dIntermediateCalc3;
        }
    }
    dForwardStdDev3 = (1./(dAnnuity*dAnnuity2) * dForwardStdDev3/(dTerminalDate-dObservationDate));

    return dForwardStdDev3/(dForwardStdDev*dForwardStdDev2);

} // dIntegratedCorrelation_VFDK_HW2F_TV

double dIntegratedCorrelation_VFDK_HW3F_TV(double dSpotDate,
        DKMaille<double> dDiscountDates,
        const DKMaille<double> &dDiscountRates,
        const DKMaille<double> &dAdjustedDiscountRates,
        DKMaille<double> dSwaptionDates,
        const DKMaille<double> &dAccrualPeriods,
        DKMaille<double> dBasisSwaptionDates,
        DKMaille<double> dBasis,
        const DKMaille<double> &dBasisAccrualPeriods,
        DKMaille<double> dSwaptionDates2,
        const DKMaille<double> &dAccrualPeriods2,
        DKMaille<double> dBasisSwaptionDates2,
        DKMaille<double> dBasis2,
        const DKMaille<double> &dBasisAccrualPeriods2,
        double dNoticePeriod,
        double dMeanReversion1,
        double dMeanReversion2,
        double dMeanReversion3,
        double dRelativeFactor12,
        double dRelativeFactor13,
        double dCorrelation12,
        double dCorrelation13,
        double dCorrelation23,
        DKMaille<double> dVolStripDates,
        DKMaille<double> &dVolStrip,
        double dJulianObservationDate)
{
    DKMaille<double> dZCRates(dDiscountDates.entries());
    DKMaille<double> dZCAdjRates(dDiscountDates.entries());
    FromDiscountToZero(dSpotDate,dZCRates,dZCAdjRates,dDiscountDates,dDiscountRates,dAdjustedDiscountRates);

    double dPrice=0.;

    unsigned int ui=0;
    for(ui=0;ui<dDiscountDates.entries(); ui++)
        dDiscountDates[ui]=(dDiscountDates[ui]-dSpotDate)/365.;

    for(ui=0;ui<dSwaptionDates.entries(); ui++)
        dSwaptionDates[ui]=(dSwaptionDates[ui]-dSpotDate)/365.;

    for(ui=0;ui<dSwaptionDates2.entries(); ui++)
        dSwaptionDates2[ui]=(dSwaptionDates2[ui]-dSpotDate)/365.;

    for(ui=0;ui<dBasisSwaptionDates.entries(); ui++)
        dBasisSwaptionDates[ui]=(dBasisSwaptionDates[ui]-dSpotDate)/365.;

    for(ui=0;ui<dBasisSwaptionDates2.entries(); ui++)
        dBasisSwaptionDates2[ui]=(dBasisSwaptionDates2[ui]-dSpotDate)/365.;

    for(ui=0;ui<dVolStripDates.entries(); ui++)
        dVolStripDates[ui]=(dVolStripDates[ui]-dSpotDate)/365.;

    DKMaille<double> dDiscountFactors(dSwaptionDates.entries());
    DKMaille<double> dDiscountFactorsOnBasisDates(dBasisSwaptionDates.entries());
    DKMaille<double> dDiscountFactors2(dSwaptionDates2.entries());
    DKMaille<double> dDiscountFactorsOnBasisDates2(dBasisSwaptionDates2.entries());

    // Build SwapRate1 analytics
    double dAnnuity = 0.;
    double dSwapRate = 0.;
    for(ui=0;ui<dSwaptionDates.entries();ui++)
    {
        dDiscountFactors[ui]=
            exp(-rateinterpolation_dk_maille(2,dSwaptionDates[ui],dDiscountDates,dZCRates,dDiscountDates.entries()-1)*dSwaptionDates[ui]);
        if(ui>0) dAnnuity+=dAccrualPeriods[ui]*dDiscountFactors[ui];
    }
    dSwapRate = (dDiscountFactors[0]-dDiscountFactors[dSwaptionDates.entries()-1])/dAnnuity;
    double dBasisAnnuity=0.;
    for(ui=0;ui<dBasisSwaptionDates.entries();ui++)
    {
        dDiscountFactorsOnBasisDates[ui]=
            exp(-rateinterpolation_dk_maille(2,dBasisSwaptionDates[ui],dDiscountDates,
                                             dZCRates,dDiscountDates.entries()-1)*dBasisSwaptionDates[ui]);
        dBasisAnnuity+=dBasis[ui]*dBasisAccrualPeriods[ui]*dDiscountFactorsOnBasisDates[ui];
    }
    dSwapRate+=(dBasisAnnuity/dAnnuity);

    // Build SwapRate2 analytics
    double dAnnuity2 = 0.;
    double dSwapRate2 = 0.;
    for(ui=0;ui<dSwaptionDates2.entries();ui++)
    {
        dDiscountFactors2[ui]=
            exp(-rateinterpolation_dk_maille(2,dSwaptionDates2[ui],dDiscountDates,dZCRates,dDiscountDates.entries()-1)*dSwaptionDates2[ui]);
        if(ui>0) dAnnuity2+=dAccrualPeriods2[ui]*dDiscountFactors2[ui];
    }
    dSwapRate2 = (dDiscountFactors2[0]-dDiscountFactors2[dSwaptionDates2.entries()-1])/dAnnuity2;
    double dBasisAnnuity2=0.;
    for(ui=0;ui<dBasisSwaptionDates2.entries();ui++)
    {
        dDiscountFactorsOnBasisDates2[ui]=
            exp(-rateinterpolation_dk_maille(2,dBasisSwaptionDates2[ui],dDiscountDates,
                                             dZCRates,dDiscountDates.entries()-1)*dBasisSwaptionDates2[ui]);
        dBasisAnnuity2+=dBasis2[ui]*dBasisAccrualPeriods2[ui]*dDiscountFactorsOnBasisDates2[ui];
    }
    dSwapRate2+=(dBasisAnnuity2/dAnnuity2);


    double dObservationDate=(dJulianObservationDate-dSpotDate)/365.;
    double dTerminalDate    = dSwaptionDates[0]-dNoticePeriod;
    if(dSwaptionDates[0]!=dSwaptionDates2[0]) throw("Start dates of underlying swaps are inconsistent. System does not yet support calculation of correlation between rates having different reset date.");

    double dPi = 3.141592653589790 ;

    DKMaille<double> dVector(dSwaptionDates.entries()+dBasisSwaptionDates.entries());
    DKMaille<double> dDiscounts(dDiscountFactors.entries()+dDiscountFactorsOnBasisDates.entries());
    DKMaille<double> dCash(dDiscountFactors.entries()+dDiscountFactorsOnBasisDates.entries());
    DKMaille<double> dDates(dSwaptionDates.entries()+dBasisSwaptionDates.entries());
    if(
        (dSwaptionDates.entries()!=dDiscountFactors.entries())
        ||(dSwaptionDates.entries()!=dAccrualPeriods.entries())
        ||(dBasisSwaptionDates.entries()!=dBasisAccrualPeriods.entries())
        ||(dBasisSwaptionDates.entries()!=dDiscountFactorsOnBasisDates.entries())
        ||(dBasisSwaptionDates.entries()!=dBasis.entries())
    )
        throw("Input deal dates are inconsistent");

    DKMaille<double> dVector2(dSwaptionDates2.entries()+dBasisSwaptionDates2.entries());
    DKMaille<double> dDiscounts2(dDiscountFactors2.entries()+dDiscountFactorsOnBasisDates2.entries());
    DKMaille<double> dCash2(dDiscountFactors2.entries()+dDiscountFactorsOnBasisDates2.entries());
    DKMaille<double> dDates2(dSwaptionDates2.entries()+dBasisSwaptionDates2.entries());
    if(
        (dSwaptionDates2.entries()!=dDiscountFactors2.entries())
        ||(dSwaptionDates2.entries()!=dAccrualPeriods2.entries())
        ||(dBasisSwaptionDates2.entries()!=dBasisAccrualPeriods2.entries())
        ||(dBasisSwaptionDates2.entries()!=dDiscountFactorsOnBasisDates2.entries())
        ||(dBasisSwaptionDates2.entries()!=dBasis2.entries())
    )
        throw("Input deal dates are inconsistent");

    // Build ItoVector for rate 1
    for(ui=0;ui<dVector.entries();ui++)
    {
        if(ui<dAccrualPeriods.entries())
        {
            dCash[ui]=dAccrualPeriods[ui]*dSwapRate;
            dDiscounts[ui]=dDiscountFactors[ui];
            dDates[ui]=dSwaptionDates[ui];
        }
        else
        {
            dCash[ui]=dBasisAccrualPeriods[ui-dAccrualPeriods.entries()]*dBasis[ui-dAccrualPeriods.entries()];
            dDiscounts[ui]=dDiscountFactorsOnBasisDates[ui-dAccrualPeriods.entries()];
            dDates[ui]=dBasisSwaptionDates[ui-dAccrualPeriods.entries()];
        }
    }
    dVector=GetItoVector(dSwaptionDates.entries(),dBasisSwaptionDates.entries(),dCash,dDiscounts);

    // Build ItoVector for rate 2
    for(ui=0;ui<dVector2.entries();ui++)
    {
        if(ui<dAccrualPeriods2.entries())
        {
            dCash2[ui]=dAccrualPeriods2[ui]*dSwapRate2;
            dDiscounts2[ui]=dDiscountFactors2[ui];
            dDates2[ui]=dSwaptionDates2[ui];
        }
        else
        {
            dCash2[ui]=dBasisAccrualPeriods2[ui-dAccrualPeriods2.entries()]*dBasis2[ui-dAccrualPeriods2.entries()];
            dDiscounts2[ui]=dDiscountFactorsOnBasisDates2[ui-dAccrualPeriods2.entries()];
            dDates2[ui]=dBasisSwaptionDates2[ui-dAccrualPeriods2.entries()];
        }
    }
    dVector2=GetItoVector(dSwaptionDates2.entries(),dBasisSwaptionDates2.entries(),dCash2,dDiscounts2);





    // New functionality
    // Make sure integration slices are flexible to allow for fast pricing of callable reverse floaters and callable CMS spread options
    // Checks on dates
    DKMaille<double> dIntegrationStrip=dCheckDates(dVolStripDates,dObservationDate,dTerminalDate);
    DKMaille<double> dIntegrationStripVols=dNewVols(dVolStripDates,dVolStrip,dIntegrationStrip);


    // Date tests
    bool bCheckDis=true;
    unsigned int uiM=0;
    for(;uiM<dIntegrationStrip.entries();uiM++)
    {
        if(fabs(dIntegrationStrip[uiM]-dTerminalDate)<0.00001)
            bCheckDis=false;
    }
    if(bCheckDis)
        throw("Dates of Std Dev Strip Are Inconsistent with Option Expiry");


    double dForwardStdDev=0.;
    double dIntermediateCalc=0.;
    double dForwardStdDev2=0.;
    double dIntermediateCalc2=0.;
    double dForwardStdDev3=0.;
    double dIntermediateCalc3=0.;



    // This is the core calculation of the VFDK approximation
    // Create second and third index
    DKMaille<double> dVolStrip2(dIntegrationStrip.entries());
    DKMaille<double> dVolStrip3(dIntegrationStrip.entries());

    for(unsigned int uiK=0;uiK<dVolStrip2.entries();uiK++) dVolStrip2[uiK]=dRelativeFactor12*dIntegrationStripVols[uiK];

    for(uiK=0;uiK<dVolStrip3.entries();uiK++) dVolStrip3[uiK]=dRelativeFactor13*dIntegrationStripVols[uiK];

    // iterate through Kac space
    // calculation of absolute standard deviation rate 1
    for(uiK=0;uiK<dDates.entries();uiK++)
    {
        for(unsigned int uiM=0;uiM<dDates.entries();uiM++)
        {
            dIntermediateCalc=dVector[uiK]*dVector[uiM]*
                              (
                                  (1./(dMeanReversion1*dMeanReversion1))*
                                  GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion1,
                                                            dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dIntegrationStripVols,dDates[uiK],dDates[uiM])
                                  +(1./(dMeanReversion2*dMeanReversion2))*
                                  GetItoDeterminant_HW2F_TV(dMeanReversion2,dMeanReversion2,
                                                            dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip2,dVolStrip2,dDates[uiK],dDates[uiM])
                                  +(1./(dMeanReversion3*dMeanReversion3))*
                                  GetItoDeterminant_HW2F_TV(dMeanReversion3,dMeanReversion3,
                                                            dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip3,dVolStrip3,dDates[uiK],dDates[uiM])
                                  +2.*dCorrelation12*
                                  (1./(dMeanReversion1*dMeanReversion2))*
                                  GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion2,
                                                            dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dVolStrip2,dDates[uiK],dDates[uiM])
                                  +2.*dCorrelation13*
                                  (1./(dMeanReversion1*dMeanReversion3))*
                                  GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion3,
                                                            dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dVolStrip3,dDates[uiK],dDates[uiM])
                                  +2.*dCorrelation23*
                                  (1./(dMeanReversion2*dMeanReversion3))*
                                  GetItoDeterminant_HW2F_TV(dMeanReversion2,dMeanReversion3,
                                                            dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip2,dVolStrip3,dDates[uiK],dDates[uiM])
                              );
            dForwardStdDev += dIntermediateCalc;
        }
    }
    dForwardStdDev = (1./dAnnuity) * sqrt(dForwardStdDev/(dTerminalDate-dObservationDate));

    // calculation of absolute standard deviation rate 2
    for(uiK=0;uiK<dDates2.entries();uiK++)
    {
        for(unsigned int uiM=0;uiM<dDates2.entries();uiM++)
        {
            dIntermediateCalc2=dVector2[uiK]*dVector2[uiM]*
                               (
                                   (1./(dMeanReversion1*dMeanReversion1))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion1,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dIntegrationStripVols,dDates2[uiK],dDates2[uiM])
                                   +(1./(dMeanReversion2*dMeanReversion2))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion2,dMeanReversion2,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip2,dVolStrip2,dDates2[uiK],dDates2[uiM])
                                   +(1./(dMeanReversion3*dMeanReversion3))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion3,dMeanReversion3,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip3,dVolStrip3,dDates2[uiK],dDates2[uiM])
                                   +2.*dCorrelation12*
                                   (1./(dMeanReversion1*dMeanReversion2))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion2,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dVolStrip2,dDates2[uiK],dDates2[uiM])
                                   +2.*dCorrelation13*
                                   (1./(dMeanReversion1*dMeanReversion3))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion3,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dVolStrip3,dDates2[uiK],dDates2[uiM])
                                   +2.*dCorrelation23*
                                   (1./(dMeanReversion2*dMeanReversion3))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion2,dMeanReversion3,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip2,dVolStrip3,dDates2[uiK],dDates2[uiM])
                               );
            dForwardStdDev2 += dIntermediateCalc2;
        }
    }
    dForwardStdDev2 = (1./dAnnuity2) * sqrt(dForwardStdDev2/(dTerminalDate-dObservationDate));

    // calculation of mixed term in the variance equation
    for(uiK=0;uiK<dDates.entries();uiK++)
    {
        for(unsigned int uiM=0;uiM<dDates2.entries();uiM++)
        {
            dIntermediateCalc3=dVector[uiK]*dVector2[uiM]*
                               (
                                   (1./(dMeanReversion1*dMeanReversion1))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion1,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dIntegrationStripVols,dDates[uiK],dDates2[uiM])
                                   +(1./(dMeanReversion2*dMeanReversion2))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion2,dMeanReversion2,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip2,dVolStrip2,dDates[uiK],dDates2[uiM])
                                   +(1./(dMeanReversion3*dMeanReversion3))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion3,dMeanReversion3,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip3,dVolStrip3,dDates[uiK],dDates2[uiM])
                                   +dCorrelation12*
                                   (1./(dMeanReversion1*dMeanReversion2))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion2,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dVolStrip2,dDates[uiK],dDates2[uiM])
                                   +dCorrelation12*
                                   (1./(dMeanReversion2*dMeanReversion1))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion2,dMeanReversion1,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip2,dIntegrationStripVols,dDates[uiK],dDates2[uiM])
                                   +dCorrelation13*
                                   (1./(dMeanReversion1*dMeanReversion3))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion1,dMeanReversion3,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dIntegrationStripVols,dVolStrip3,dDates[uiK],dDates2[uiM])
                                   +dCorrelation13*
                                   (1./(dMeanReversion3*dMeanReversion1))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion3,dMeanReversion1,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip3,dIntegrationStripVols,dDates[uiK],dDates2[uiM])
                                   +dCorrelation23*
                                   (1./(dMeanReversion2*dMeanReversion3))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion2,dMeanReversion3,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip2,dVolStrip3,dDates[uiK],dDates2[uiM])
                                   +dCorrelation23*
                                   (1./(dMeanReversion3*dMeanReversion2))*
                                   GetItoDeterminant_HW2F_TV(dMeanReversion3,dMeanReversion2,
                                                             dObservationDate,dTerminalDate,dIntegrationStrip,dVolStrip3,dVolStrip2,dDates[uiK],dDates2[uiM])
                               );
            dForwardStdDev3 += dIntermediateCalc3;
        }
    }
    dForwardStdDev3 = (1./(dAnnuity*dAnnuity2) * dForwardStdDev3/(dTerminalDate-dObservationDate));

    return dForwardStdDev3/(dForwardStdDev*dForwardStdDev2);

} // dIntegratedCorrelation_VFDK_HW3F_TV

double ImpliedFwdCorrelation_VFDK_HW1To3F(double dSwaptionExpiryInYears,
        double dSwaptionTenorInYears,
        double dSwaptionTenor2InYears,
        double dNoticePeriodInDays,
        DKMaille<double> dDiscountCurveDates,
        DKMaille<double> dDiscountCurveRates,
        DKMaille<double> dAdjustedDiscountCurveRates,
        DKMaille<double> dNoticeDates,
        DKMaille<double> dTimeDependentStandardDeviation,
        DKMaille<double> dModelParameters,
        double dJulianObservationDate)
{

    unsigned int uiFixedIncomeFactors = (unsigned int)dModelParameters.at(0);

    if(dTimeDependentStandardDeviation.at(dTimeDependentStandardDeviation.entries()-1)!=
            dTimeDependentStandardDeviation.at(dTimeDependentStandardDeviation.entries()-3))
        throw("Wrong Input Model Vols");

    if(dTimeDependentStandardDeviation.at(dTimeDependentStandardDeviation.entries()-1)!=
            dTimeDependentStandardDeviation.at(dTimeDependentStandardDeviation.entries()-2))
        throw("Wrong Input Model Vols");

    if(dTimeDependentStandardDeviation.entries()!=dNoticeDates.entries())
        throw("Wrong Input Dates in Model Vols");


    DKMaille<double> swaptionDates;      // swap start date, then all n payment dates of the fixed leg.
    DKMaille<double> basisSwaptionDates; // all m payment dates of the floating leg.
    DKMaille<double> accrualPeriods;     // bond basis adjusted period length of fixed leg in basis years (model years) array of size n+1 first element is zero.
    DKMaille<double> basisAccrualPeriods;// bond basis adjusted period length of floating leg in basis years (model years)
    DKMaille<double> basis;              // basis adjustments

    DKMaille<double> swaptionDates2;      // swap start date, then all n payment dates of the fixed leg.
    DKMaille<double> basisSwaptionDates2; // all m payment dates of the floating leg.
    DKMaille<double> accrualPeriods2;     // bond basis adjusted period length of fixed leg in basis years (model years) array of size n+1 first element is zero.
    DKMaille<double> basisAccrualPeriods2;// bond basis adjusted period length of floating leg in basis years (model years)
    DKMaille<double> basis2;              // basis adjustments

    double dSpotDate=dNoticeDates.at(0);

    DKMaille2D<double> dVolInfo(dNoticeDates.entries(),2);

    // Create the underlying vanilla swap for relevant notice date
    CreateVanillaSwap(dSpotDate,									// SpotDate input
                      dSwaptionTenorInYears,			// dSwapTenorInYears input
                      dSwaptionExpiryInYears+(dJulianObservationDate-dSpotDate)/365.,			// dOptionExpiryInYears input
                      swaptionDates,
                      accrualPeriods,
                      basisSwaptionDates,
                      basisAccrualPeriods,
                      basis,
                      dNoticePeriodInDays/365.,		// NoticePeriodInYears
                      dDiscountCurveDates,
                      dDiscountCurveRates,
                      dAdjustedDiscountCurveRates);

    // Create the underlying vanilla swap for relevant notice date
    CreateVanillaSwap(dSpotDate,									// SpotDate input
                      dSwaptionTenor2InYears,			// dSwapTenorInYears input
                      dSwaptionExpiryInYears+(dJulianObservationDate-dSpotDate)/365.,			// dOptionExpiryInYears input
                      swaptionDates2,
                      accrualPeriods2,
                      basisSwaptionDates2,
                      basisAccrualPeriods2,
                      basis2,
                      dNoticePeriodInDays/365.,		// NoticePeriodInYears
                      dDiscountCurveDates,
                      dDiscountCurveRates,
                      dAdjustedDiscountCurveRates);


    double dFwdCorrelation=0.;

    // Get Absolute Volatility
    switch (uiFixedIncomeFactors)
    {
    case 2:
        {

            // Get Model parameters
            double dMeanReversion1 = dModelParameters.at(1);
            double dMeanReversion2 = dModelParameters.at(2);
            double dRelativeFactor = dModelParameters.at(4);
            double dCorrelation    = dModelParameters.at(6);

            dFwdCorrelation=dIntegratedCorrelation_VFDK_HW2F_TV(dSpotDate,
                            dDiscountCurveDates,
                            dDiscountCurveRates,
                            dAdjustedDiscountCurveRates,
                            swaptionDates,
                            accrualPeriods,
                            basisSwaptionDates,
                            basis,
                            basisAccrualPeriods,
                            swaptionDates2,
                            accrualPeriods2,
                            basisSwaptionDates2,
                            basis2,
                            basisAccrualPeriods2,
                            ((swaptionDates.at(0)-dSpotDate)/365.-(dSwaptionExpiryInYears+(dJulianObservationDate-dSpotDate)/365.)),
                            dMeanReversion1,
                            dMeanReversion2,
                            dRelativeFactor,
                            dCorrelation,
                            dNoticeDates,
                            dTimeDependentStandardDeviation,
                            dJulianObservationDate);

        } // case 2-factor
        break;
    case 1:
        {
            // Get Model parameters
            double dMeanReversion1 = dModelParameters.at(1);

            dFwdCorrelation = 1.0;

        } // case 1-factor
        break;
    case 3:
        {
            // Get Model parameters
            double dMeanReversion1 = dModelParameters.at(1);
            double dMeanReversion2 = dModelParameters.at(2);
            double dMeanReversion3 = dModelParameters.at(3);
            double dRelativeFactor12 = dModelParameters.at(4);
            double dRelativeFactor13 = dModelParameters.at(5);
            double dCorrelation12    = dModelParameters.at(6);
            double dCorrelation13    = dModelParameters.at(7);
            double dCorrelation23    = dModelParameters.at(8);

            dFwdCorrelation=dIntegratedCorrelation_VFDK_HW3F_TV(dSpotDate,
                            dDiscountCurveDates,
                            dDiscountCurveRates,
                            dAdjustedDiscountCurveRates,
                            swaptionDates,
                            accrualPeriods,
                            basisSwaptionDates,
                            basis,
                            basisAccrualPeriods,
                            swaptionDates2,
                            accrualPeriods2,
                            basisSwaptionDates2,
                            basis2,
                            basisAccrualPeriods2,
                            ((swaptionDates.at(0)-dSpotDate)/365.-(dSwaptionExpiryInYears+(dJulianObservationDate-dSpotDate)/365.)),
                            dMeanReversion1,
                            dMeanReversion2,
                            dMeanReversion3,
                            dRelativeFactor12,
                            dRelativeFactor13,
                            dCorrelation12,
                            dCorrelation13,
                            dCorrelation23,
                            dNoticeDates,
                            dTimeDependentStandardDeviation,
                            dJulianObservationDate);
        }
        break;
    } // switch


    return dFwdCorrelation;

} // ImpliedFwdCorrelation_VFDK_VFDK_HW1To3F()



double ForwardFXVolatility3FConstant(double dTimeToReset,
                                     double dDomesticStdDev,
                                     double dForeignStdDev,
                                     double dSpotFXVol,
                                     double dDomesticMeanReversion,
                                     double dForeignMeanReversion,
                                     double dCorrelationSpotFXDomestic,
                                     double dCorrelationSpotFXForeign,
                                     double dCorrelationForeignDomestic)
{
    double dTerm1 = dSpotFXVol * dSpotFXVol;

    double dTerm2 = pow(dForeignStdDev / dForeignMeanReversion,2.)*(1.
                    -2.*(1.-exp(-dForeignMeanReversion*dTimeToReset))/(dForeignMeanReversion*dTimeToReset)
                    +(1./(2.*dForeignMeanReversion*dTimeToReset))*(1.-exp(-2.*dForeignMeanReversion*dTimeToReset)));

    double dTerm3 = pow(dDomesticStdDev / dDomesticMeanReversion,2.)*(1.
                    -2.*(1.-exp(-dDomesticMeanReversion*dTimeToReset))/(dDomesticMeanReversion*dTimeToReset)
                    +(1./(2.*dDomesticMeanReversion*dTimeToReset))*(1.-exp(-2.*dDomesticMeanReversion*dTimeToReset)));

    double dTerm4 = - 2.*dCorrelationSpotFXForeign*dSpotFXVol*(dForeignStdDev/dForeignMeanReversion)*
                    (1.-(1.-exp(-dForeignMeanReversion*dTimeToReset))/(dForeignMeanReversion*dTimeToReset));

    double dTerm5 = 2.*dCorrelationSpotFXDomestic*dSpotFXVol*(dDomesticStdDev/dDomesticMeanReversion)*
                    (1.-(1.-exp(-dDomesticMeanReversion*dTimeToReset))/(dDomesticMeanReversion*dTimeToReset));

    double dTerm6 = -2.* dCorrelationForeignDomestic * (dDomesticStdDev/dDomesticMeanReversion) * (dForeignStdDev/dForeignMeanReversion)
                    *(1.
                      -(1.-exp(-dDomesticMeanReversion*dTimeToReset))/(dDomesticMeanReversion*dTimeToReset)
                      -(1.-exp(-dForeignMeanReversion*dTimeToReset))/(dForeignMeanReversion*dTimeToReset)
                      +(1./((dForeignMeanReversion+dDomesticMeanReversion)*dTimeToReset))
                      *(1.-exp(-(dForeignMeanReversion+dDomesticMeanReversion)*dTimeToReset)));

    double dForwardFXVol = sqrt(dTerm1+dTerm2+dTerm3+dTerm4+dTerm5+dTerm6);

    return dForwardFXVol;
} // LDHD_3F_FX_ANALYTICS(...)






double ForwardFXVolatility3FTD(double dSpotDate,
                               double dStartDate,
                               double dEndDate,
                               double dForwardMaturityDate,
                               DKMaille<double> dDomesticStrip,
                               DKMaille<double> dForeignStrip,
                               DKMaille<double> dSpotFXStrip,
                               DKMaille<double> dDomesticStdDev,
                               DKMaille<double> dForeignStdDev,
                               DKMaille<double> dSpotFXVol,
                               double dDomesticMeanReversion,
                               double dForeignMeanReversion,
                               double dCorrelationSpotFXDomestic,
                               double dCorrelationSpotFXForeign,
                               double dCorrelationForeignDomestic)
{

    DKMaille<double> dStrip=CreateDateStrip(dSpotDate,dDomesticStrip,dForeignStrip,dSpotFXStrip);

    DKMaille<double> dStripDomesticStdDev(dStrip.entries());
    DKMaille<double> dStripForeignStdDev(dStrip.entries());
    DKMaille<double> dStripSpotFXVol(dStrip.entries());

    AggregateVolStrip(dStrip,dDomesticStrip,dDomesticStdDev,dStripDomesticStdDev);
    AggregateVolStrip(dStrip,dForeignStrip,dForeignStdDev,dStripForeignStdDev);
    AggregateVolStrip(dStrip,dSpotFXStrip,dSpotFXVol,dStripSpotFXVol);

    double dResult=0.;

    dResult=GetForwardFXVol(dSpotDate,
                            dStartDate,
                            dEndDate,
                            dForwardMaturityDate,
                            dStrip,
                            dStripDomesticStdDev,
                            dStripForeignStdDev,
                            dStripSpotFXVol,
                            dDomesticMeanReversion,
                            dForeignMeanReversion,
                            dCorrelationSpotFXDomestic,
                            dCorrelationSpotFXForeign,
                            dCorrelationForeignDomestic);

    return dResult;

}


double ForwardFXVolatility3FTD_DK(double dSpotDate,
                                  double dStartDate,
                                  double dEndDate,
                                  double dForwardMaturityDate,
                                  DKMaille<double> dDomesticStrip,
                                  DKMaille<double> dForeignStrip,
                                  DKMaille<double> dSpotFXStrip,
                                  DKMaille<double> dDomesticStdDev,
                                  DKMaille<double> dForeignStdDev,
                                  DKMaille<double> dSpotFXVol,
                                  double dDomesticMeanReversion,
                                  double dForeignMeanReversion,
                                  double dCorrelationSpotFXDomestic,
                                  double dCorrelationSpotFXForeign,
                                  double dCorrelationForeignDomestic,
                                  DKMaille<double> shortRateDates,
                                  DKMaille<double> shortRates,
                                  DKMaille<double> shortRatesForeign)
{

    DKMaille<double> dStrip=CreateDateStrip(dSpotDate,dDomesticStrip,dForeignStrip,dSpotFXStrip);

    DKMaille<double> dStripDomesticStdDev(dStrip.entries());
    DKMaille<double> dStripForeignStdDev(dStrip.entries());
    DKMaille<double> dStripSpotFXVol(dStrip.entries());

    AggregateVolStrip(dStrip,dDomesticStrip,dDomesticStdDev,dStripDomesticStdDev);
    AggregateVolStrip(dStrip,dForeignStrip,dForeignStdDev,dStripForeignStdDev);
    AggregateVolStrip(dStrip,dSpotFXStrip,dSpotFXVol,dStripSpotFXVol);

    double dResult=0.;

    dResult=GetForwardFXVol_DK(dSpotDate,
                               dStartDate,
                               dEndDate,
                               dForwardMaturityDate,
                               dStrip,
                               dStripDomesticStdDev,
                               dStripForeignStdDev,
                               dStripSpotFXVol,
                               dDomesticMeanReversion,
                               dForeignMeanReversion,
                               dCorrelationSpotFXDomestic,
                               dCorrelationSpotFXForeign,
                               dCorrelationForeignDomestic,
                               shortRateDates,
                               shortRates,
                               shortRatesForeign);

    return dResult;

}


void TransformVolatilitiesSingle( DKMaille<double> dStripDates,
                                  DKMaille<double> &dNewStripDates,
                                  DKMaille<double> &dVolStrip)
{
    DKMaille<double> dVolStripPrime(dNewStripDates.entries());

    for(unsigned int uf=0;uf<dNewStripDates.entries();uf++)
    {
        dVolStripPrime.at(uf)=numerical_function(dNewStripDates.at(uf),dStripDates,dVolStrip);
    }

    dVolStrip.resize(dNewStripDates.entries());

    for(uf=0;uf<dNewStripDates.entries();uf++)
    {
        dVolStrip.at(uf)=dVolStripPrime.at(uf);
    }
}





void TransformVolatilities(DKMaille<double> dDomesticVolStripDates,
                           DKMaille<double> dForeignVolStripDates,
                           DKMaille<double> dSpotFXStrip,
                           DKMaille<double> &dNewStripDates,
                           DKMaille<double> &dDomesticVolStrip,
                           DKMaille<double> &dForeignVolStrip,
                           DKMaille<double> &dSpotFXVol)
{
    DKMaille<double> dDomesticVolStripPrime(dNewStripDates.entries());
    DKMaille<double> dForeignVolStripPrime(dNewStripDates.entries());
    DKMaille<double> dSpotFXVolPrime(dNewStripDates.entries());

    for(unsigned int uf=0;uf<dNewStripDates.entries();uf++)
    {
        dDomesticVolStripPrime.at(uf)=numerical_function(dNewStripDates.at(uf),dDomesticVolStripDates,dDomesticVolStrip);
        dForeignVolStripPrime.at(uf)=numerical_function(dNewStripDates.at(uf),dForeignVolStripDates,dForeignVolStrip);
        dSpotFXVolPrime.at(uf)=numerical_function(dNewStripDates.at(uf),dSpotFXStrip,dSpotFXVol);
    }
    dDomesticVolStrip.resize(dNewStripDates.entries());
    dForeignVolStrip.resize(dNewStripDates.entries());
    dSpotFXVol.resize(dNewStripDates.entries());

    for(uf=0;uf<dNewStripDates.entries();uf++)
    {
        dDomesticVolStrip.at(uf)=dDomesticVolStripPrime.at(uf);
        dForeignVolStrip.at(uf)=dForeignVolStripPrime.at(uf);
        dSpotFXVol.at(uf)=dSpotFXVolPrime.at(uf);
    }
}


void TransformSwaptionVolatilities(DKMaille<double> dDomesticVolStripDates,DKMaille<double> dForeignVolStripDates,
                                   DKMaille<double> dSpotFXStrip,DKMaille<double> &dNewStripDates,DKMaille<double>
                                   &dDomesticVolStrip,DKMaille<double> &dForeignVolStrip)
{
    unsigned int uiCount=0;

    for(unsigned int ui=0;ui<dSpotFXStrip.entries();ui++)
    {
        if(dSpotFXStrip.at(ui)<1.000) uiCount++;
    }
    unsigned int uiAfter=49*4+1;

    // Modify interpolation to match the lattice
    dNewStripDates.resize(uiAfter+uiCount*4);

    int uf1=0;
    for(unsigned int uf=0;uf<uiCount;uf++)
    {
        dNewStripDates.at(uf1)=dSpotFXStrip.at(uf);
        uf1++;
        for(int i=1;i<4;i++)
        {
            dNewStripDates.at(uf1)=dSpotFXStrip.at(uf)+0.25*i*(dSpotFXStrip.at(uf+1)-dSpotFXStrip.at(uf));
            uf1++;
        }
    }
    for(uf=0;uf<uiAfter;uf++)
    {
        dNewStripDates.at(uf1+uf)=1.+(((double)(uf))*0.25);
    }

    uiCount*=4;

    DKMaille<double> dDomesticVolStripPrime(uiCount+uiAfter); DKMaille<double> dForeignVolStripPrime(uiCount+uiAfter);
    for(uf=0;uf<uiAfter+uiCount;uf++)
    {
        dDomesticVolStripPrime.at(uf)=numerical_function(dNewStripDates.at(uf),dDomesticVolStripDates,dDomesticVolStrip);
        dForeignVolStripPrime.at(uf)=numerical_function(dNewStripDates.at(uf),dForeignVolStripDates,dForeignVolStrip);
    }
    dDomesticVolStrip.resize(uiAfter+uiCount); dForeignVolStrip.resize(uiAfter+uiCount);
    for(uf=0;uf<uiAfter+uiCount;uf++)
    {
        dDomesticVolStrip.at(uf)=dDomesticVolStripPrime.at(uf);
        dForeignVolStrip.at(uf)=dForeignVolStripPrime.at(uf);
    }
}

void TransformFXVolatilities(DKMaille<double> dSpotFXStrip,
                             DKMaille<double> &dNewStripDates,
                             DKMaille<double> &dSpotFXVol)
{
    unsigned int uiCount=0;

    for(unsigned int ui=0;ui<dSpotFXStrip.entries();ui++)
    {
        if(dSpotFXStrip.at(ui)<1.000) uiCount++;
    }
    unsigned int uiAfter=49*4+1;

    // Modify interpolation to match the lattice
    dNewStripDates.resize(uiAfter+uiCount*4);

    int uf1=0;
    for(unsigned int uf=0;uf<uiCount;uf++)
    {
        dNewStripDates.at(uf1)=dSpotFXStrip.at(uf);
        uf1++;
        for(int i=1;i<4;i++)
        {
            dNewStripDates.at(uf1)=dSpotFXStrip.at(uf)+0.25*i*(dSpotFXStrip.at(uf+1)-dSpotFXStrip.at(uf));
            uf1++;
        }
    }
    for(uf=0;uf<uiAfter;uf++)
    {
        dNewStripDates.at(uf1+uf)=1.+(((double)(uf))*0.25);			}

    uiCount*=4;

    DKMaille<double> dSpotFXVolPrime(uiCount+uiAfter);
    for(uf=0;uf<uiAfter+uiCount;uf++)
    {
        dSpotFXVolPrime.at(uf)=numerical_function(dNewStripDates.at(uf),dSpotFXStrip,dSpotFXVol);
    }
    dSpotFXVol.resize(uiAfter+uiCount);
    for(uf=0;uf<uiAfter+uiCount;uf++)
    {
        dSpotFXVol.at(uf)=dSpotFXVolPrime.at(uf);
    }
}


DKMaille<double>  &FromSpotToForwardVol3F(	double dSpotDate,
        DKMaille<double> dStartDate,
        DKMaille<double> dEndDate,
        DKMaille<double> dForwardMaturityDate,
        DKMaille<double> dBaseDiscountCurveDates,
        DKMaille<double> dBaseDiscountCurveRates,
        DKMaille<double> dBaseAdjustedDiscountCurveRates,
        DKMaille<double> dForeignDiscountCurveDates,
        DKMaille<double> dForeignDiscountCurveRates,
        DKMaille<double> dForeignAdjustedDiscountCurveRates,
        DKMaille<double> dDomesticVolX,
        DKMaille<double> dDomesticVolY,
        DKMaille2D<double> dDomesticVolZ,
        DKMaille<double> dForeignVolX,
        DKMaille<double> dForeignVolY,
        DKMaille2D<double> dForeignVolZ,
        DKMaille<double> dSpotFXStrip,
        DKMaille<double> dSpotFXVol,
        double dDomesticMeanReversion,
        double dForeignMeanReversion,
        double dCorrelationSpotFXDomestic,
        double dCorrelationSpotFXForeign,
        double dCorrelationForeignDomestic,
        double dCurrencyPair1, // 0=jpy 1=usd 2=aud 3=eur
        double dCurrencyPair2, // 0=jpy 1=usd 2=aud 3=eur
        DKMaille<double> dNoticeDates,
        DKMaille<double> dFXCouponResetDates,
        DKMaille<double> dFXCouponPaymentDates,
        double dIsSwaptionCalibrationWithBasis) 
{
    DKMaille<double> dStripDatesVolFX;
    DKMaille<double> dStripDatesVolBase;
    DKMaille<double> dStripDatesVolForeign;
    DKMaille<double> dVolFX;
    DKMaille<double> dVolBase;
    DKMaille<double> dVolForeign;

    int maxLength = MAX(MAX(dNoticeDates.entries(), dFXCouponResetDates.entries()),
                        dFXCouponPaymentDates.entries());
    DKMaille2D<double> dBoosterData(maxLength, 3);
    for (int i = 0; i < maxLength; i++)
    {
        if ( i < dNoticeDates.entries() )
            dBoosterData.at(i,0) = dNoticeDates.at(i);
        else
            dBoosterData.at(i,0)=0.0;

        if ( i < dFXCouponResetDates.entries() )
            dBoosterData.at(i,1) = dFXCouponResetDates.at(i);
        else
            dBoosterData.at(i,1) = 0.;

        if ( i < dFXCouponPaymentDates.entries() )
            dBoosterData.at(i,2) = dFXCouponPaymentDates.at(i);
        else
            dBoosterData.at(i,2) = 0.;
    }


    Processing_3F_Data( dStripDatesVolFX,
                        dStripDatesVolBase,
                        dStripDatesVolForeign,
                        dVolFX,
                        dVolBase,
                        dVolForeign,
                        4,  /////////////// a determiner
                        3,  /////////////// a determiner
                        1,  /////////////// a determiner
                        dSpotDate,
                        dBaseDiscountCurveDates,
                        dBaseAdjustedDiscountCurveRates,
                        dForeignDiscountCurveDates,
                        dForeignAdjustedDiscountCurveRates,
                        dDomesticVolX,
                        dDomesticVolY,
                        dDomesticVolZ,
                        dDomesticMeanReversion,
                        dForeignVolX,
                        dForeignVolY,
                        dForeignVolZ,
                        dForeignMeanReversion,
                        dSpotFXStrip,
                        dSpotFXVol,
                        dCorrelationForeignDomestic,
                        dCorrelationSpotFXDomestic,
                        dCorrelationSpotFXForeign,
                        10, /////////////// a determiner
                        dBoosterData,
                        dBaseDiscountCurveRates,
                        dForeignDiscountCurveRates,
                        dIsSwaptionCalibrationWithBasis);


    // les stripdates sont renvoyes en date julienne, on les remplace par leurs valeurs en annees
    for(i=0;i<dStripDatesVolFX.entries();i++)
        dStripDatesVolFX.at(i)=(dStripDatesVolFX.at(i)-dSpotDate)/365.;
    for(i=0;i<dStripDatesVolBase.entries();i++)
        dStripDatesVolBase.at(i)=(dStripDatesVolBase.at(i)-dSpotDate)/365.;
    for(i=0;i<dStripDatesVolForeign.entries();i++)
        dStripDatesVolForeign.at(i)=(dStripDatesVolForeign.at(i)-dSpotDate)/365.;


    // Calculate fx forward vol
    DKMaille<double> dForwardFXVol;
    dForwardFXVol.resize(dForwardMaturityDate.entries());

    for(unsigned int uk=0;uk<dForwardMaturityDate.entries();uk++)
    {
        if(dForwardMaturityDate.at(uk)<=0.)
            dForwardFXVol.at(uk)=dSpotFXVol.at(0);
        else
        {
            dForwardFXVol.at(uk)=ForwardFXVolatility3FTD(0.,
                                 dStartDate.at(uk),
                                 dEndDate.at(uk),
                                 dForwardMaturityDate.at(uk),
                                 dStripDatesVolBase,
                                 dStripDatesVolForeign,
                                 dStripDatesVolFX,
                                 dVolBase,
                                 dVolForeign,
                                 dVolFX,
                                 dDomesticMeanReversion,
                                 dForeignMeanReversion,
                                 dCorrelationSpotFXDomestic,
                                 dCorrelationSpotFXForeign,
                                 dCorrelationForeignDomestic);
        }
    }
    static DKMaille<double> dResult;
    dResult.resize(dForwardMaturityDate.entries());
    for(uk=0;uk<dForwardMaturityDate.entries();uk++) dResult.at(uk)=dForwardFXVol.at(uk);
    return dResult;
}
// FromSpotToForwardVol3F


DKMaille<double>  &FromSpotToForwardVolDK3F(double dSpotDate,
        DKMaille<double> dStartDate,
        DKMaille<double> dEndDate,
        DKMaille<double> dForwardMaturityDate,
        DKMaille<double> dBaseDiscountCurveDates,
        DKMaille<double> dBaseDiscountCurveRates,
        DKMaille<double> dBaseAdjustedDiscountCurveRates,
        DKMaille<double> dForeignDiscountCurveDates,
        DKMaille<double> dForeignDiscountCurveRates,
        DKMaille<double> dForeignAdjustedDiscountCurveRates,
        DKMaille<double> dDomesticVolX,
        DKMaille<double> dDomesticVolY,
        DKMaille2D<double> dDomesticVolZ,
        DKMaille<double> dForeignVolX,
        DKMaille<double> dForeignVolY,
        DKMaille2D<double> dForeignVolZ,
        DKMaille<double> dSpotFXStrip,
        DKMaille<double> dSpotFXVol,
        double dDomesticMeanReversion,
        double dForeignMeanReversion,
        double dCorrelationSpotFXDomestic,
        double dCorrelationSpotFXForeign,
        double dCorrelationForeignDomestic,
        double dCurrencyPair1, // 0=jpy 1=usd 2=aud 3=eur
        double dCurrencyPair2, // 0=jpy 1=usd 2=aud 3=eur
        DKMaille<double> dNoticeDates,
        DKMaille<double> dFXCouponResetDates,
        DKMaille<double> dFXCouponPaymentDates,
        double dIsSwaptionCalibrationWithBasis) 
{

    DKMaille<double> dShortRates;
    DKMaille<double> dShortRateDates;
    DKMaille<double> dShortRatesForeign;
    CalcXCCYShortRates(dSpotDate,
                       dBaseDiscountCurveDates,
                       dBaseAdjustedDiscountCurveRates,
                       dBaseAdjustedDiscountCurveRates,
                       dForeignDiscountCurveDates,
                       dForeignAdjustedDiscountCurveRates,
                       dForeignAdjustedDiscountCurveRates,
                       dShortRateDates,
                       dShortRates,
                       dShortRatesForeign);



    DKMaille<double> dStripDatesVolFX;
    DKMaille<double> dStripDatesVolBase;
    DKMaille<double> dStripDatesVolForeign;
    DKMaille<double> dVolFX;
    DKMaille<double> dVolBase;
    DKMaille<double> dVolForeign;

    int maxLength = MAX(MAX(dNoticeDates.entries(), dFXCouponResetDates.entries()),
                        dFXCouponPaymentDates.entries());
    DKMaille2D<double> dBoosterData(maxLength, 3);
    for (int i = 0; i < maxLength; i++)
    {
        if ( i < dNoticeDates.entries() )
            dBoosterData.at(i,0) = dNoticeDates.at(i);
        else
            dBoosterData.at(i,0)=0.0;

        if ( i < dFXCouponResetDates.entries() )
            dBoosterData.at(i,1) = dFXCouponResetDates.at(i);
        else
            dBoosterData.at(i,1) = 0.;

        if ( i < dFXCouponPaymentDates.entries() )
            dBoosterData.at(i,2) = dFXCouponPaymentDates.at(i);
        else
            dBoosterData.at(i,2) = 0.;
    }


    Processing_3F_Data_DK(dStripDatesVolFX,
                          dStripDatesVolBase,
                          dStripDatesVolForeign,
                          dVolFX,
                          dVolBase,
                          dVolForeign,
                          4,  /////////////// a determiner
                          3,  /////////////// a determiner
                          1,  /////////////// a determiner
                          dSpotDate,
                          dBaseDiscountCurveDates,
                          dBaseAdjustedDiscountCurveRates,
                          dForeignDiscountCurveDates,
                          dForeignAdjustedDiscountCurveRates,
                          dDomesticVolX,
                          dDomesticVolY,
                          dDomesticVolZ,
                          dDomesticMeanReversion,
                          dForeignVolX,
                          dForeignVolY,
                          dForeignVolZ,
                          dForeignMeanReversion,
                          dSpotFXStrip,
                          dSpotFXVol,
                          dCorrelationForeignDomestic,
                          dCorrelationSpotFXDomestic,
                          dCorrelationSpotFXForeign,
                          10, /////////////// a determiner
                          dBoosterData,
                          dBaseDiscountCurveRates,
                          dForeignDiscountCurveRates,
                          dIsSwaptionCalibrationWithBasis);


    // les stripdates sont renvoyes en date julienne, on les remplace par leurs valeurs en annees
    for(i=0;i<dStripDatesVolFX.entries();i++)
        dStripDatesVolFX.at(i)=(dStripDatesVolFX.at(i)-dSpotDate)/365.;
    for(i=0;i<dStripDatesVolBase.entries();i++)
        dStripDatesVolBase.at(i)=(dStripDatesVolBase.at(i)-dSpotDate)/365.;
    for(i=0;i<dStripDatesVolForeign.entries();i++)
        dStripDatesVolForeign.at(i)=(dStripDatesVolForeign.at(i)-dSpotDate)/365.;



    // Calculate fx forward vol
    DKMaille<double> dForwardFXVol;
    dForwardFXVol.resize(dForwardMaturityDate.entries());


    for(unsigned int uk=0;uk<dForwardMaturityDate.entries();uk++)
    {
        if(dForwardMaturityDate.at(uk)<=0.)
            dForwardFXVol.at(uk)=dSpotFXVol.at(0);
        else
        {
            dForwardFXVol.at(uk)=ForwardFXVolatility3FTD_DK(0.,
                                 dStartDate.at(uk),
                                 dEndDate.at(uk),
                                 dForwardMaturityDate.at(uk),
                                 dStripDatesVolBase,
                                 dStripDatesVolForeign,
                                 dStripDatesVolFX,
                                 dVolBase,
                                 dVolForeign,
                                 dVolFX,
                                 dDomesticMeanReversion,
                                 dForeignMeanReversion,
                                 dCorrelationSpotFXDomestic,
                                 dCorrelationSpotFXForeign,
                                 dCorrelationForeignDomestic,
                                 dShortRateDates,
                                 dShortRates,
                                 dShortRatesForeign);
        }
    }
    static DKMaille<double> dResult;
    dResult.resize(dForwardMaturityDate.entries());
    for(uk=0;uk<dForwardMaturityDate.entries();uk++) dResult.at(uk)=dForwardFXVol.at(uk);
    return dResult;
}
// FromSpotToForwardVolDK3F




// middle bounds function
double numerical_derivative(double dDate,
                            DKMaille<double> dDates,
                            DKMaille<double> dFunction)
{
    int iFlag=0.;
    unsigned int uiK=0;
    double dAnchorDate=dDates.at(0);

    if(dDates.at(0)!=0.)
        throw("Error Input Vol Function");

    if(dDates.entries()<3)
        throw("Error Input Vol Function");

    // if(dFunction.at(dFunction.entries()-1)!=dFunction.at(dFunction.entries()-2))
    //	throw("Error Input Vol Function");

    // if((dDate<=dAnchorDate+(dDates.at(1)-dDates.at(0))/2.)||(dDate>=dDates.at(dDates.entries()-2)+(dDates.at(dDates.entries()-1)-dDates.at(dDates.entries()-2))/2.))
    //		return 0.;

    if((dDate>=dDates.at(dDates.entries()-2)+(dDates.at(dDates.entries()-1)-dDates.at(dDates.entries()-2))/2.))
        return 0.;

    // European option
    if(dDates.entries()==3) return 0.;

    // Initial derivative
    double dX10=(dDates.at(1)-dDates.at(0))/2.;
    double dX20=(dDates.at(2)-dDates.at(1))/2.;
    if(dDate<dDates.at(0)+dX10)
    {
        double dXX1=dDates.at(0)+dX10;
        double dXX2=dDates.at(0+1)+dX20;
        double dYY1=dFunction.at(0);
        double dYY2=dFunction.at(0+1);
        // return 0.;
        return (dYY1-dYY2)/(dXX1-dXX2);
    }


    while(iFlag==0&&uiK<dDates.entries()-2)
    {
        double dX1=(dDates.at(uiK+1)-dDates.at(uiK))/2.;
        double dX2=(dDates.at(uiK+2)-dDates.at(uiK+1))/2.;
        if(dDate>=dDates.at(uiK)+dX1&&dDate<dDates.at(uiK+1)+dX2)
        {
            iFlag=1;
            double dXX1=dDates.at(uiK)+dX1;
            double dXX2=dDates.at(uiK+1)+dX2;
            double dYY1=dFunction.at(uiK);
            double dYY2=dFunction.at(uiK+1);
            return (dYY1-dYY2)/(dXX1-dXX2);
        }
        uiK++;
    }

    throw("ERROR: Out of bounds in derivative interpolation");

}

double numerical_function(double dDate,
                          DKMaille<double> dDates,
                          DKMaille<double> dFunction)
{

    // return rateinterpolation_dk_maille(0,dDate,dDates,dFunction,dFunction.entries()-1);

    if(dDate>dDates.at(dDates.entries()-1)) return dFunction.at(dFunction.entries()-1);

    if(dDates.at(0)!=0.)
        throw("Error Input Vol Function");

    if(dDates.entries()<3)
        throw("Error Input Vol Function");


    if(dFunction.at(dFunction.entries()-1)!=dFunction.at(dFunction.entries()-2))
        throw("Error Input Vol Function");

    // European option
    if(dDates.entries()==3) return dFunction.at(0);

    double a=0.;
    double b=0.;
    // Left edge
    if((dDate>=dDates.at(0))&&(dDate<(dDates.at(1)+0.5*(dDates.at(2)-dDates.at(1)))))
    {
        //    if(dDate<dDates.at(1)/2.)
        //		{
        //      a=0.;
        //      b=dFunction.at(0);
        //		}
        //	else
        {
            double dT01=0.5*(dDates.at(1));
            a=numerical_derivative(dDate,dDates,dFunction);
            b=dFunction.at(0)-a*dT01;
        }
    }
    // middle
    unsigned int uiK=1;
    int iFlag=0;
    while(iFlag==0&&uiK<dDates.entries()-2)
    {
        double dX1=0.5*(dDates.at(uiK+1)-dDates.at(uiK));
        double dX2=0.5*(dDates.at(uiK+2)-dDates.at(uiK+1));
        if(dDate>=dDates.at(uiK)+dX1&&dDate<dDates.at(uiK+1)+dX2)
        {
            iFlag=1;
            double dT01=dDates.at(uiK)+dX1;
            a=numerical_derivative(dDate,dDates,dFunction);
            b=dFunction.at(uiK)-a*dT01;
        }
        uiK++;
    }
    // Right edge
    if(dDate>=dDates.at(dDates.entries()-2)+0.5*(dDates.at(dDates.entries()-1)-dDates.at(dDates.entries()-2)))
    {
        a=0.0;
        b=dFunction.at(dFunction.entries()-1);
    }

    return a*dDate+b;
}

DKMaille<double> BootstrappingSpotFXVolatility3FTD(double dSpotDate,
        DKMaille<double> dStartDate,
        DKMaille<double> dEndDate,
        DKMaille<double> dForwardMaturityDate,
        DKMaille<double> dForwardVol,
        DKMaille<double> dDomesticStrip,
        DKMaille<double> dForeignStrip,
        DKMaille<double> dDomesticStdDev,
        DKMaille<double> dForeignStdDev,
        double dDomesticMeanReversion,
        double dForeignMeanReversion,
        double dCorrelationSpotFXDomestic,
        double dCorrelationSpotFXForeign,
        double dCorrelationForeignDomestic)

{
    static DKMaille<double> dResult;
    dResult.resize(dForwardMaturityDate.entries());
    DKMaille<double> dSpotFXDateStrip(dForwardMaturityDate.entries()+1);
    dSpotFXDateStrip.at(0)=dSpotDate;
    for(unsigned int ui=0;ui<dForwardMaturityDate.entries();ui++)
        dSpotFXDateStrip.at(ui+1)=dForwardMaturityDate.at(ui);
    for(unsigned int uiK=0;uiK<dForwardMaturityDate.entries();uiK++)
    {
        int iCounter=0;
        double dGuess=dForwardVol.at(uiK);
        dResult.at(uiK)=dGuess;
        while(iCounter<100)
        {
            iCounter++;
            dResult.at(uiK)=dGuess;
            // Function
            double dFunction=ForwardFXVolatility3FTD(dSpotDate,dStartDate.at(uiK),dEndDate.at(uiK),dForwardMaturityDate.at(uiK),
                             dDomesticStrip,dForeignStrip,dSpotFXDateStrip,dDomesticStdDev,dForeignStdDev,
                             dResult,dDomesticMeanReversion,dForeignMeanReversion,dCorrelationSpotFXDomestic,
                             dCorrelationSpotFXForeign,dCorrelationForeignDomestic)-dForwardVol.at(uiK);
            // FunctionUp
            dResult.at(uiK)=dGuess+0.0001;
            double dFunctionUp=ForwardFXVolatility3FTD(dSpotDate,dStartDate.at(uiK),dEndDate.at(uiK),dForwardMaturityDate.at(uiK),
                               dDomesticStrip,dForeignStrip,dSpotFXDateStrip,dDomesticStdDev,dForeignStdDev,
                               dResult,dDomesticMeanReversion,dForeignMeanReversion,dCorrelationSpotFXDomestic,
                               dCorrelationSpotFXForeign,dCorrelationForeignDomestic)-dForwardVol.at(uiK);
            double dRatio=0.0001*dFunction/(dFunctionUp-dFunction);
            dGuess-=dRatio;
            dResult.at(uiK)=dGuess;
            if(fabs(dRatio)<.1e-4)
            {
                iCounter=100;
            }
        } // while iCounter < 100
    }
    return dResult;
}

DKMaille<double> BootstrappingSpotFXVolatility3FTD_DK(double dSpotDate,
        DKMaille<double> dStartDate,
        DKMaille<double> dEndDate,
        DKMaille<double> dForwardMaturityDate,
        DKMaille<double> dForwardVol,
        DKMaille<double> dDomesticStrip,
        DKMaille<double> dForeignStrip,
        DKMaille<double> dDomesticStdDev,
        DKMaille<double> dForeignStdDev,
        double dDomesticMeanReversion,
        double dForeignMeanReversion,
        double dCorrelationSpotFXDomestic,
        double dCorrelationSpotFXForeign,
        double dCorrelationForeignDomestic,
        DKMaille<double> dCurveDates,
        DKMaille<double> dCurve,
        DKMaille<double> dCurveForeignDates,
        DKMaille<double> dCurveForeign)

{
    static DKMaille<double> dResult;
    dResult.resize(dForwardMaturityDate.entries());
    DKMaille<double> dSpotFXDateStrip(dForwardMaturityDate.entries()+1);
    dSpotFXDateStrip.at(0)=dSpotDate;
    for(unsigned int ui=0;ui<dForwardMaturityDate.entries();ui++)
        dSpotFXDateStrip.at(ui+1)=dForwardMaturityDate.at(ui);
    DKMaille<double> dShortRates;
    DKMaille<double> dShortRateDates;
    DKMaille<double> dShortRatesForeign;
    CalcXCCYShortRates(dCurveDates.at(0),dCurveDates,dCurve,dCurve,dCurveForeignDates,dCurveForeign,dCurveForeign,dShortRateDates,dShortRates,dShortRatesForeign);

    for(unsigned int uiK=0;uiK<dForwardMaturityDate.entries();uiK++)
    {
        int iCounter=0;
        double dGuess=dForwardVol.at(uiK);
        dResult.at(uiK)=dGuess;
        while(iCounter<100)
        {
            iCounter++;
            dResult.at(uiK)=dGuess;
            // Function
            double dFunction=ForwardFXVolatility3FTD_DK(dSpotDate,dStartDate.at(uiK),dEndDate.at(uiK),dForwardMaturityDate.at(uiK),
                             dDomesticStrip,dForeignStrip,dSpotFXDateStrip,dDomesticStdDev,dForeignStdDev,
                             dResult,dDomesticMeanReversion,dForeignMeanReversion,dCorrelationSpotFXDomestic,
                             dCorrelationSpotFXForeign,dCorrelationForeignDomestic,dShortRateDates,dShortRates,dShortRatesForeign)
                             -dForwardVol.at(uiK);
            // FunctionUp
            dResult.at(uiK)=dGuess+0.0001;
            double dFunctionUp=ForwardFXVolatility3FTD_DK(dSpotDate,dStartDate.at(uiK),dEndDate.at(uiK),dForwardMaturityDate.at(uiK),
                               dDomesticStrip,dForeignStrip,dSpotFXDateStrip,dDomesticStdDev,dForeignStdDev,
                               dResult,dDomesticMeanReversion,dForeignMeanReversion,dCorrelationSpotFXDomestic,
                               dCorrelationSpotFXForeign,dCorrelationForeignDomestic,dShortRateDates,dShortRates,dShortRatesForeign)
                               -dForwardVol.at(uiK);
            double dRatio=0.0001*dFunction/(dFunctionUp-dFunction);
            dGuess-=dRatio;
            dResult.at(uiK)=dGuess;
            if(fabs(dRatio)<.1e-4)
            {
                iCounter=100;
            }
        } // while iCounter < 100
    }
    return dResult;
}




double BondFromVFDKAnalytics_HW1F(double dNorm,
                                  double dDecay,
                                  double dBrownian)
{
    return dNorm*exp(-dDecay*dBrownian);
} // BondFromVFDKAnalytics_HW1F(...)

double BondFromVFDKAnalytics_HW2F(double dNorm,
                                  double dDecay1,
                                  double dBrownian1,
                                  double dDecay2,
                                  double dBrownian2)
{
    return dNorm*exp(-dDecay1*dBrownian1-dDecay2*dBrownian2);
} // BondFromVFDKAnalytics_HW2F(...)

double BondFromVFDKAnalytics_HW3F(double dNorm,
                                  double dDecay1,
                                  double dBrownian1,
                                  double dDecay2,
                                  double dBrownian2,
                                  double dDecay3,
                                  double dBrownian3)
{
    return dNorm*exp(-dDecay1*dBrownian1-dDecay2*dBrownian2-dDecay3*dBrownian3);
} // BondFromVFDKAnalytics_HW3F(...)


double ShuffleVol(double dDate,DKMaille<double> &dInputVolDates,DKMaille<double> &dInputVol)
{

    unsigned int uk=0;
    int iFlag=0;
    double dReturn;
    while(iFlag==0&&uk<dInputVol.entries()-1)
    {
        if((dDate<dInputVolDates.at(uk+1))&&(dDate>=dInputVolDates.at(uk)))
        {
            dReturn=dInputVol.at(uk);
            iFlag=1;
        }
        uk++;
    }
    if(dDate>=dInputVolDates.at(dInputVolDates.entries()-1)) dReturn=dInputVol.at(dInputVolDates.entries()-1);

    return dReturn;

}

double  Analytics_HWVFDK_3F_(double dSpotDate,
                             DKMaille<double> dBaseDates,
                             DKMaille<double> dBaseRates,
                             DKMaille<double> dForeignDates,
                             DKMaille<double> dForeignRates,
                             DKMaille<double> dStdDevBaseX,
                             DKMaille<double> dStdDevBaseZ,
                             double dMeanReversionBase,
                             DKMaille<double> dStdDevForeignX,
                             DKMaille<double> dStdDevForeignZ,
                             double dMeanReversionForeign,
                             double dSpotFX,
                             DKMaille<double> dSpotFXVolDatesTD,
                             DKMaille<double> dSpotFXVolTD,
                             double dBaseForeignCorrelation,
                             double dBaseSpotFXCorrelation,
                             double dForeignSpotFXCorrelation,
                             DKMaille2D<double> dBoosterData,
                             DKMaille<double> dRedemptionData)
{

    unsigned int uiPeriods=dBoosterData.rows();

    // Re-shuffle vols
    // FX Reset Dates

    DKMaille<double> _dStdDevBaseZ(dSpotFXVolDatesTD.entries());
    DKMaille<double> _dStdDevForeignZ(dSpotFXVolDatesTD.entries());
    DKMaille<double> _dSpotFXVolTD(dSpotFXVolDatesTD.entries());

    DKMaille<double> dZCBaseRates(dBaseDates.entries());
    DKMaille<double> dZCForeignRates(dForeignDates.entries());
    FromDiscountToZero(dSpotDate,dZCBaseRates,dZCBaseRates,dBaseDates,dBaseRates,dBaseRates);
    FromDiscountToZero(dSpotDate,dZCForeignRates,dZCForeignRates,dForeignDates,dForeignRates,dForeignRates);


    for(unsigned int ui=0;ui<dSpotFXVolDatesTD.entries();ui++)
    {
        _dStdDevBaseZ.at(ui)=ShuffleVol(dSpotFXVolDatesTD.at(ui),dStdDevBaseX,dStdDevBaseZ);
        _dStdDevForeignZ.at(ui)=ShuffleVol(dSpotFXVolDatesTD.at(ui),dStdDevForeignX,dStdDevForeignZ);
        _dSpotFXVolTD.at(ui)=ShuffleVol(dSpotFXVolDatesTD.at(ui),dSpotFXVolDatesTD,dSpotFXVolTD);
    }

    double dPV=0.;
    for(ui=0;ui<uiPeriods;ui++)
    {
        // FX Coupon
        // Calculate FX Vol
        double dForwardFXVol=GetForwardFXVol(dSpotDate,
                                             dSpotDate,
                                             dBoosterData.at(ui,1),
                                             dBoosterData.at(ui,1),
                                             dSpotFXVolDatesTD,
                                             _dStdDevBaseZ,
                                             _dStdDevForeignZ,
                                             _dSpotFXVolTD,
                                             dMeanReversionBase,
                                             dMeanReversionForeign,
                                             dBaseSpotFXCorrelation,
                                             dForeignSpotFXCorrelation,
                                             dBaseForeignCorrelation);

        // FXCoupon
        double dFXResetDate=dBoosterData.at(ui,1);
        double dFXPaymentDate=dBoosterData.at(ui,2);
        double dStrike=dBoosterData.at(ui,8);
        double dCap=dBoosterData.at(ui,9);
        if(dFXResetDate!=dFXPaymentDate)
            throw("Analytics for convexity not implemented yet");
        double dDomesticZC=rateinterpolation_dk_maille(2,dFXResetDate,dBaseDates,dZCBaseRates,dBaseDates.entries()-1);
        double dForeignZC=rateinterpolation_dk_maille(2,dFXResetDate,dForeignDates,dZCForeignRates,dForeignDates.entries()-1);
        double dFXForward=dSpotFX*exp((dDomesticZC-dForeignZC)*(dFXResetDate-dSpotDate)/365.);
        double dDomesticZCPayment=rateinterpolation_dk_maille(2,dFXPaymentDate,dBaseDates,dZCBaseRates,dBaseDates.entries()-1);
        double dDiscount=exp(-dDomesticZCPayment*(dFXPaymentDate-dSpotDate)/365.);
        double dOptionFX=dBlackScholes((dFXResetDate-dSpotDate)/365.,dForwardFXVol,dStrike,dFXForward,0.)*dDiscount;
        double dOptionCap=dBlackScholes((dFXResetDate-dSpotDate)/365.,dForwardFXVol,dCap,dFXForward,0.)*dDiscount;
        double dFXAccrualPeriod=dBoosterData.at(ui,10);
        double dNotionalMultiplier=dBoosterData.at(ui,7);
        dOptionFX=dOptionFX-dOptionCap;
        double dFXCoupon=10000.*dOptionFX*dFXAccrualPeriod*dNotionalMultiplier;

        // Funding
        double dLIBORStartDate=dBoosterData.at(ui,3);
        double dLIBOREndDate=dBoosterData.at(ui,4);
        double dLIBORPaymentDate=dBoosterData.at(ui,5);
        double dZCStart=rateinterpolation_dk_maille(2,dLIBORStartDate,dBaseDates,dZCBaseRates,dBaseDates.entries()-1);
        double dZCEnd=rateinterpolation_dk_maille(2,dLIBOREndDate,dBaseDates,dZCBaseRates,dBaseDates.entries()-1);
        double dZCPayment=rateinterpolation_dk_maille(2,dLIBORPaymentDate,dBaseDates,dZCBaseRates,dBaseDates.entries()-1);
        double dDiscountStart=exp(-dZCStart*(dLIBORStartDate-dSpotDate)/365.);
        double dDiscountEnd=exp(-dZCEnd*(dLIBOREndDate-dSpotDate)/365.);
        double dDiscountPayment=exp(-dZCPayment*(dLIBORPaymentDate-dSpotDate)/365.);
        double dLIBORAccrualPeriod=dBoosterData.at(ui,11);
        double dLIBOR=(dDiscountStart/dDiscountEnd-1)/dLIBORAccrualPeriod;
        double dSpread=dBoosterData.at(ui,6);
        double dFunding=(dLIBOR+dSpread)*dLIBORAccrualPeriod*10000.*dDiscountPayment;
        dPV+=dFunding-dFXCoupon;
    }

    if(dRedemptionData.at(0)>0)
    {
        if(dRedemptionData.at(1)==0)
        {
            double dBaseDiscount=rateinterpolation_dk_maille(	2, dRedemptionData.at(3), dBaseDates, dBaseRates, dBaseRates.entries()-1);
            double dForeignDiscount=rateinterpolation_dk_maille(	2, dRedemptionData.at(3), dForeignDates, dForeignRates, dForeignDates.entries()-1);
            double dForward=dSpotFX*dForeignDiscount/dBaseDiscount;
            dPV-=10000*dBaseDiscount*(dForward-dRedemptionData.at(2))/dRedemptionData.at(2);
        }
        else
        {
            double dBaseDiscount=rateinterpolation_dk_maille(	2, dRedemptionData.at(4), dBaseDates, dBaseRates, dBaseRates.entries()-1);
            double dBaseDiscount1=rateinterpolation_dk_maille(	2, dRedemptionData.at(3), dBaseDates, dBaseRates, dBaseRates.entries()-1);
            double dForeignDiscount=rateinterpolation_dk_maille(	2, dRedemptionData.at(4), dForeignDates, dForeignRates, dForeignDates.entries()-1);
            double dForward=dSpotFX*dForeignDiscount/dBaseDiscount;
            double dTimeToExpiry=(dRedemptionData.at(4)-dSpotDate)/365;
            double dForwardFXVol=GetForwardFXVol(	dSpotDate,
                                                  dSpotDate,
                                                  dRedemptionData.at(4),
                                                  dRedemptionData.at(4),
                                                  dSpotFXVolDatesTD,
                                                  _dStdDevBaseZ,
                                                  _dStdDevForeignZ,
                                                  _dSpotFXVolTD,
                                                  dMeanReversionBase,
                                                  dMeanReversionForeign,
                                                  dBaseSpotFXCorrelation,
                                                  dForeignSpotFXCorrelation,
                                                  dBaseForeignCorrelation);
            double* dNormDist1=new double;
            double* dNormDist2=new double;

            double aa=-(log(dForward/dRedemptionData.at(2))+0.5*dForwardFXVol*dForwardFXVol*dTimeToExpiry)/dForwardFXVol/sqrt(dTimeToExpiry);
            double bb=-(log(dForward/dRedemptionData.at(2))-0.5*dForwardFXVol*dForwardFXVol*dTimeToExpiry)/dForwardFXVol/sqrt(dTimeToExpiry);

            int i1=normalDK(dNormDist1,aa);
            int i2=normalDK(dNormDist2,bb);

            dPV+=10000*(dBaseDiscount1/dRedemptionData.at(2))*(-dForward**dNormDist1+dRedemptionData.at(2)**dNormDist2);
        }
    }


    return dPV;

}

double  Analytics_DK_3F_(double dSpotDate,
                         DKMaille<double> dBaseDates,
                         DKMaille<double> dBaseRates,
                         DKMaille<double> dForeignDates,
                         DKMaille<double> dForeignRates,
                         DKMaille<double> dStdDevBaseX,
                         DKMaille<double> dStdDevBaseZ,
                         double dMeanReversionBase,
                         DKMaille<double> dStdDevForeignX,
                         DKMaille<double> dStdDevForeignZ,
                         double dMeanReversionForeign,
                         double dSpotFX,
                         DKMaille<double> dSpotFXVolDatesTD,
                         DKMaille<double> dSpotFXVolTD,
                         double dBaseForeignCorrelation,
                         double dBaseSpotFXCorrelation,
                         double dForeignSpotFXCorrelation,
                         DKMaille2D<double> dBoosterData,
                         DKMaille<double> dRedemptionData)
{

    unsigned int uiPeriods=dBoosterData.rows();

    // Re-shuffle vols
    // FX Reset Dates

    DKMaille<double> _dStdDevBaseZ(dSpotFXVolDatesTD.entries());
    DKMaille<double> _dStdDevForeignZ(dSpotFXVolDatesTD.entries());
    DKMaille<double> _dSpotFXVolTD(dSpotFXVolDatesTD.entries());

    DKMaille<double> dZCBaseRates(dBaseDates.entries());
    DKMaille<double> dZCForeignRates(dForeignDates.entries());
    FromDiscountToZero(dSpotDate,dZCBaseRates,dZCBaseRates,dBaseDates,dBaseRates,dBaseRates);
    FromDiscountToZero(dSpotDate,dZCForeignRates,dZCForeignRates,dForeignDates,dForeignRates,dForeignRates);


    for(unsigned int ui=0;ui<dSpotFXVolDatesTD.entries();ui++)
    {
        _dStdDevBaseZ.at(ui)=ShuffleVol(dSpotFXVolDatesTD.at(ui),dStdDevBaseX,dStdDevBaseZ);
        _dStdDevForeignZ.at(ui)=ShuffleVol(dSpotFXVolDatesTD.at(ui),dStdDevForeignX,dStdDevForeignZ);
        _dSpotFXVolTD.at(ui)=ShuffleVol(dSpotFXVolDatesTD.at(ui),dSpotFXVolDatesTD,dSpotFXVolTD);
    }

    DKMaille<double> dShortRates;
    DKMaille<double> dShortRateDates;
    DKMaille<double> dShortRatesForeign;
    CalcXCCYShortRates(dSpotDate,dBaseDates,dBaseRates,dBaseRates,dForeignDates,dForeignRates,dForeignRates,dShortRateDates,dShortRates,dShortRatesForeign);

    double dPV=0.;
    for(ui=0;ui<uiPeriods;ui++)
    {
        // FX Coupon
        // Calculate FX Vol
        double dForwardFXVol=GetForwardFXVol_DK(dSpotDate,
                                                dSpotDate,
                                                dBoosterData.at(ui,1),
                                                dBoosterData.at(ui,1),
                                                dSpotFXVolDatesTD,
                                                _dStdDevBaseZ,
                                                _dStdDevForeignZ,
                                                _dSpotFXVolTD,
                                                dMeanReversionBase,
                                                dMeanReversionForeign,
                                                dBaseSpotFXCorrelation,
                                                dForeignSpotFXCorrelation,
                                                dBaseForeignCorrelation,
                                                dShortRateDates,
                                                dShortRates,
                                                dShortRatesForeign);

        // FXCoupon
        double dFXResetDate=dBoosterData.at(ui,1);
        double dFXPaymentDate=dBoosterData.at(ui,2);
        double dStrike=dBoosterData.at(ui,8);
        double dCap=dBoosterData.at(ui,9);
        if(dFXResetDate!=dFXPaymentDate)
            throw("Analytics for convexity not implemented yet");
        double dDomesticZC=rateinterpolation_dk_maille(2,dFXResetDate,dBaseDates,dZCBaseRates,dBaseDates.entries()-1);
        double dForeignZC=rateinterpolation_dk_maille(2,dFXResetDate,dForeignDates,dZCForeignRates,dForeignDates.entries()-1);
        double dFXForward=dSpotFX*exp((dDomesticZC-dForeignZC)*(dFXResetDate-dSpotDate)/365.);
        double dDomesticZCPayment=rateinterpolation_dk_maille(2,dFXPaymentDate,dBaseDates,dZCBaseRates,dBaseDates.entries()-1);
        double dDiscount=exp(-dDomesticZCPayment*(dFXPaymentDate-dSpotDate)/365.);
        double dOptionFX=dBlackScholes((dFXResetDate-dSpotDate)/365.,dForwardFXVol,dStrike,dFXForward,0.)*dDiscount;
        double dOptionCap=dBlackScholes((dFXResetDate-dSpotDate)/365.,dForwardFXVol,dCap,dFXForward,0.)*dDiscount;
        double dFXAccrualPeriod=dBoosterData.at(ui,10);
        double dNotionalMultiplier=dBoosterData.at(ui,7);
        dOptionFX=dOptionFX-dOptionCap;
        double dFXCoupon=10000.*dOptionFX*dFXAccrualPeriod*dNotionalMultiplier;

        // Funding
        double dLIBORStartDate=dBoosterData.at(ui,3);
        double dLIBOREndDate=dBoosterData.at(ui,4);
        double dLIBORPaymentDate=dBoosterData.at(ui,5);
        double dZCStart=rateinterpolation_dk_maille(2,dLIBORStartDate,dBaseDates,dZCBaseRates,dBaseDates.entries()-1);
        double dZCEnd=rateinterpolation_dk_maille(2,dLIBOREndDate,dBaseDates,dZCBaseRates,dBaseDates.entries()-1);
        double dZCPayment=rateinterpolation_dk_maille(2,dLIBORPaymentDate,dBaseDates,dZCBaseRates,dBaseDates.entries()-1);
        double dDiscountStart=exp(-dZCStart*(dLIBORStartDate-dSpotDate)/365.);
        double dDiscountEnd=exp(-dZCEnd*(dLIBOREndDate-dSpotDate)/365.);
        double dDiscountPayment=exp(-dZCPayment*(dLIBORPaymentDate-dSpotDate)/365.);
        double dLIBORAccrualPeriod=dBoosterData.at(ui,11);
        double dLIBOR=(dDiscountStart/dDiscountEnd-1)/dLIBORAccrualPeriod;
        double dSpread=dBoosterData.at(ui,6);
        double dFunding=(dLIBOR+dSpread)*dLIBORAccrualPeriod*10000.*dDiscountPayment;
        dPV+=dFunding-dFXCoupon;
    }

    if(dRedemptionData.at(0)>0)
    {
        if(dRedemptionData.at(1)==0)
        {
            double dBaseDiscount=rateinterpolation_dk_maille(	2, dRedemptionData.at(3), dBaseDates, dBaseRates, dBaseRates.entries()-1);
            double dForeignDiscount=rateinterpolation_dk_maille(	2, dRedemptionData.at(3), dForeignDates, dForeignRates, dForeignDates.entries()-1);
            double dForward=dSpotFX*dForeignDiscount/dBaseDiscount;
            dPV-=10000*dBaseDiscount*(dForward-dRedemptionData.at(2))/dRedemptionData.at(2);
        }
        else
        {
            double dBaseDiscount=rateinterpolation_dk_maille(	2, dRedemptionData.at(4), dBaseDates, dBaseRates, dBaseRates.entries()-1);
            double dBaseDiscount1=rateinterpolation_dk_maille(	2, dRedemptionData.at(3), dBaseDates, dBaseRates, dBaseRates.entries()-1);
            double dForeignDiscount=rateinterpolation_dk_maille(	2, dRedemptionData.at(4), dForeignDates, dForeignRates, dForeignDates.entries()-1);
            double dForward=dSpotFX*dForeignDiscount/dBaseDiscount;
            double dTimeToExpiry=(dRedemptionData.at(4)-dSpotDate)/365;
            double dForwardFXVol=GetForwardFXVol_DK(dSpotDate,
                                                    dSpotDate,
                                                    dRedemptionData.at(4),
                                                    dRedemptionData.at(4),
                                                    dSpotFXVolDatesTD,
                                                    _dStdDevBaseZ,
                                                    _dStdDevForeignZ,
                                                    _dSpotFXVolTD,
                                                    dMeanReversionBase,
                                                    dMeanReversionForeign,
                                                    dBaseSpotFXCorrelation,
                                                    dForeignSpotFXCorrelation,
                                                    dBaseForeignCorrelation,
                                                    dShortRateDates,
                                                    dShortRates,
                                                    dShortRatesForeign);
            double* dNormDist1=new double;
            double* dNormDist2=new double;

            double aa=-(log(dForward/dRedemptionData.at(2))+0.5*dForwardFXVol*dForwardFXVol*dTimeToExpiry)/dForwardFXVol/sqrt(dTimeToExpiry);
            double bb=-(log(dForward/dRedemptionData.at(2))-0.5*dForwardFXVol*dForwardFXVol*dTimeToExpiry)/dForwardFXVol/sqrt(dTimeToExpiry);

            int i1=normalDK(dNormDist1,aa);
            int i2=normalDK(dNormDist2,bb);

            dPV+=10000*(dBaseDiscount1/dRedemptionData.at(2))*(-dForward**dNormDist1+dRedemptionData.at(2)**dNormDist2);
        }
    }


    return dPV;

}


DKMaille2D<double> dItoTensor(DKMaille<double> &dDiscounts,
                              DKMaille<double> &dPeriods,
                              DKMaille<double> &dCashBasis,
                              unsigned int uiFixedCouponsExcludingLast,
                              unsigned int uiFloatCoupons,
                              double dSwapRate,
                              double dAnnuity)
{
    DKMaille2D<double> dTensor(dDiscounts.entries(),dDiscounts.entries());
    unsigned int uiSize=dDiscounts.entries();
    //////////////////////////////
    //  diagonal elements       //
    //////////////////////////////
    // forward date
    dTensor.at(0,0)=0.;
    // payment date
    dTensor.at(1,1)=0.;
    // fixed coupon dates excluding last
    unsigned int ui=0;
    for(ui=2;ui<2+uiFixedCouponsExcludingLast;ui++)
    {
        dTensor.at(ui,ui)=2.*dSwapRate*dDiscounts.at(1)*pow(dPeriods.at(ui)/dAnnuity,2.);
    }
    // float coupon dates
    for(ui=2+uiFixedCouponsExcludingLast;ui<2+uiFixedCouponsExcludingLast+uiFloatCoupons;ui++)
    {
        dTensor.at(ui,ui)=0.;
    }
    // last fixed coupon payment date
    dTensor.at(uiSize-1,uiSize-1)=2.*dDiscounts.at(1)*dPeriods.at(uiSize-1)*(1.+dPeriods.at(uiSize-1)*dSwapRate)/pow(dAnnuity,2.);
    //////////////////////////////
    //  non-diagonal elements   //
    //////////////////////////////
    // forward date against payment date
    dTensor.at(0,1)=1./dAnnuity;
    dTensor.at(1,0)=dTensor.at(0,1);
    // forward date against fixed coupon dates excluding last
    for(ui=2;ui<2+uiFixedCouponsExcludingLast;ui++)
    {
        dTensor.at(0,ui)=-dPeriods.at(ui)*dDiscounts.at(1)/pow(dAnnuity,2.);
        dTensor.at(ui,0)=dTensor.at(0,ui);
    }
    // forward date against float coupon dates
    for(ui=2+uiFixedCouponsExcludingLast;ui<2+uiFixedCouponsExcludingLast+uiFloatCoupons;ui++)
    {
        dTensor.at(0,ui)=0.;
        dTensor.at(ui,0)=0.;
    }
    // forward date against final fixed coupon date
    dTensor.at(0,uiSize-1)=-dPeriods.at(ui)*dDiscounts.at(1)/pow(dAnnuity,2.);
    dTensor.at(uiSize-1,0)=dTensor.at(0,uiSize-1);
    // payment date against fixed coupon dates excluding last
    for(ui=2;ui<2+uiFixedCouponsExcludingLast;ui++)
    {
        dTensor.at(1,ui)=-dPeriods.at(ui)*dSwapRate/dAnnuity;
        dTensor.at(ui,1)=dTensor.at(1,ui);
    }
    // payment date against float coupon dates
    for(ui=2+uiFixedCouponsExcludingLast;ui<2+uiFixedCouponsExcludingLast+uiFloatCoupons;ui++)
    {
        dTensor.at(1,ui)=+dCashBasis.at(ui)/dAnnuity;
        dTensor.at(ui,1)=dTensor.at(1,ui);
    }
    // payment date against final fixed coupon date
    dTensor.at(1,uiSize-1)=-1./dAnnuity-dPeriods.at(uiSize-1)*dSwapRate/dAnnuity;;
    dTensor.at(uiSize-1,1)=dTensor.at(1,uiSize-1);

    // fixed coupon dates excluding last against fixed coupon dates excluding last except when different
    for(ui=2;ui<2+uiFixedCouponsExcludingLast;ui++)
    {
        for(unsigned int uj=2;uj<2+uiFixedCouponsExcludingLast;uj++)
        {
            if(uj!=ui)
            {
                dTensor.at(ui,uj)=2.*dSwapRate*dDiscounts.at(1)*dPeriods.at(ui)*dPeriods.at(uj)*pow(1./dAnnuity,2.);
            }
        }
    }

    // fixed coupon dates excluding last against floating coupon dates
    for(ui=2;ui<2+uiFixedCouponsExcludingLast;ui++)
    {
        for(unsigned int uj=2+uiFixedCouponsExcludingLast;uj<2+uiFixedCouponsExcludingLast+uiFloatCoupons;uj++)
        {
            dTensor.at(ui,uj)=-dPeriods.at(ui)*dCashBasis.at(uj)*dDiscounts.at(1)/pow(dAnnuity,2.);
            dTensor.at(uj,ui)=dTensor.at(ui,uj);
        }
    }
    // fixed coupon dates excluding last against last fixed coupon date
    for(ui=2;ui<2+uiFixedCouponsExcludingLast;ui++)
    {
        dTensor.at(ui,uiSize-1)=(dPeriods.at(ui)*dDiscounts.at(1)/pow(dAnnuity,2.))*(1.+2.*dPeriods.at(uiSize-1)*dSwapRate);
        dTensor.at(uiSize-1,ui)=dTensor.at(ui,uiSize-1);
    }
    // float coupon dates against float coupon date when different
    for(ui=2+uiFixedCouponsExcludingLast;ui<2+uiFixedCouponsExcludingLast+uiFloatCoupons;ui++)
    {
        for(unsigned int uj=2+uiFixedCouponsExcludingLast;uj<2+uiFixedCouponsExcludingLast+uiFloatCoupons;uj++)
        {
            if(ui!=uj)
            {
                dTensor.at(ui,uj)=0.;
            }
        }
    }
    // float coupon dates against last fixed coupon date
    unsigned int uj=0;
    for(uj=2+uiFixedCouponsExcludingLast;uj<2+uiFixedCouponsExcludingLast+uiFloatCoupons;uj++)
    {
        dTensor.at(uj,uiSize-1)=-dPeriods.at(uiSize-1)*dCashBasis(uj)*dDiscounts.at(1)/pow(dAnnuity,2.);
        dTensor.at(uiSize-1,uj)=dTensor.at(uj,uiSize-1);
    }

    // Multiply by discount factors

    for(ui=0;ui<uiSize;ui++)
    {
        for(uj=0;uj<uiSize;uj++)
        {
            dTensor.at(ui,uj)*=dDiscounts.at(ui)*dDiscounts.at(uj);
        }
    }

    // Done

    return dTensor;
}


double ki_cap_analytics_price(double dSpotDate,
                              double dDateOne,
                              double dDateTwo,
                              double dStdDev1,
                              double dStdDev2,
                              double dRate1,
                              double dRate2,
                              double dStrike1,
                              double dStrike2,
                              double dCorrelation,
                              double dOutput)
{
    double dPi = 3.141592653589790 ;

    dDateOne=(dDateOne-dSpotDate)/365.;
    dDateTwo=(dDateTwo-dSpotDate)/365.;

    double dArgument1=(dRate1-dStrike1)/(dStdDev1*sqrt(dDateOne));
    double dArgument2=(dRate2-dStrike2)/(dStdDev2*sqrt(dDateTwo));

    double dIntegral = CumulativeNormalBivariateIntegral(dArgument1,dArgument2,dCorrelation);

    double dIntegralDerivative1 = (CumulativeNormalBivariateIntegral(dArgument1+0.0001,dArgument2,dCorrelation)-
                                   CumulativeNormalBivariateIntegral(dArgument1,dArgument2,dCorrelation))/0.0001;
    double dIntegralDerivative2 = (CumulativeNormalBivariateIntegral(dArgument1,dArgument2+0.0001,dCorrelation)-
                                   CumulativeNormalBivariateIntegral(dArgument1,dArgument2,dCorrelation))/0.0001;

    double dN1=0.; double dN2=0.;

    // if 1 is a lognormal variable then dN1=N(d2)
    normalDK(&dN1,dArgument1);
    normalDK(&dN2,dArgument2);

    double dKnockInCap=(dRate2-dStrike2)*dIntegral + dStdDev2*sqrt(dDateTwo)*(dIntegralDerivative2+dCorrelation*dIntegralDerivative1);
    double dKnockInForward=dRate2*dN1 + dCorrelation*dStdDev2*sqrt(dDateTwo)*(1./sqrt(2.*dPi))*exp(-0.5*(dArgument1*dArgument1));

    double dXX=0.;

    if(dOutput==1.) dXX = dKnockInForward;
    if(dOutput==2.) dXX = dKnockInCap;

    return dXX;

}// ki_cap_analytics_price


void Processing_3F_Data( DKMaille<double> &dStripDatesVolFX,
                         DKMaille<double> &dStripDatesVolBase,
                         DKMaille<double> &dStripDatesVolForeign,
                         DKMaille<double> &dVolFX,
                         DKMaille<double> &dVolBase,
                         DKMaille<double> &dVolForeign,
                         double dNumTimeLinesBeforeFirstNotice,
                         double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                         double dNumTimeLinesPerYear,
                         double dSpotDate,
                         DKMaille<double> &dBaseDates,
                         DKMaille<double> &dBaseRates,
                         DKMaille<double> &dForeignDates,
                         DKMaille<double> &dForeignRates,
                         DKMaille<double> &dStdDevBaseX,
                         DKMaille<double> &dStdDevBaseY,
                         DKMaille2D<double> &dStdDevBaseZ,
                         double dMeanReversionBase,
                         DKMaille<double> &dStdDevForeignX,
                         DKMaille<double> &dStdDevForeignY,
                         DKMaille2D<double> &dStdDevForeignZ,
                         double dMeanReversionForeign,
                         DKMaille<double> dSpotFXVolDatesTD,
                         DKMaille<double> dSpotFXVolTD,
                         double dBaseForeignCorrelation,
                         double dBaseSpotFXCorrelation,
                         double dForeignSpotFXCorrelation,
                         double dOptimal,
                         DKMaille2D<double> &dBoosterData,
                         DKMaille<double> &dBaseRatesNoBasis,
                         DKMaille<double> &dForeignRatesNoBasis,
                         double dIsSwaptionCalibrationWithBasis)
{

    DKMaille<double> NoticeDates(dBoosterData.rows());

    for(unsigned int i=0;i<dBoosterData.rows();i++)
    {
        NoticeDates.at(i)=(dBoosterData.at(i,0)-dSpotDate)/365.;
    }


    DKMaille<double> NoticeDatesDomestic;
    DKMaille<double> VolDomestic;
    DKMaille<double> NoticeDatesForeign;
    DKMaille<double> VolForeign;

    DKMaille<double> dNoticeDatesBis;
    DKMaille<double> dNoticeDates_(dBoosterData.rows());
    for(unsigned int uiBD=0;uiBD<dBoosterData.rows();uiBD++)
        dNoticeDates_.at(uiBD)=(dBoosterData.at(uiBD,0)-dSpotDate)/365.;
    DKMaille<double> dLocalNoticeDates;
    unsigned int uiNoticeSize=0;
    for(uiBD=0;uiBD<dNoticeDates_.entries();uiBD++)
    {
        if(dNoticeDates_.at(uiBD)==0.)
            throw("Notice time is the same as evaluation time. Not authorised");
        // only value the future portion
        if(dNoticeDates_.at(uiBD)>0.)
        {
            if(dNoticeDates_.at(uiBD)!=(0.-dSpotDate)/365.)
            {
                dLocalNoticeDates.insert(dNoticeDates_.at(uiBD));
                uiNoticeSize++;
            }
        }
    }
    dNoticeDatesBis.resize(uiNoticeSize);
    for(uiBD=0;uiBD<uiNoticeSize;uiBD++)
    {
        dNoticeDatesBis.at(uiBD)=dLocalNoticeDates.at(uiBD);
    }



    if(dIsSwaptionCalibrationWithBasis==1.0) 
    {
        BootstrappedVolsSingleCurrencyFixedIncome1F( NoticeDatesDomestic,
                VolDomestic,
                dNoticeDatesBis,
                dBaseDates,
                dBaseRates,
                dBaseRatesNoBasis,
                dStdDevBaseX,
                dStdDevBaseY,
                dStdDevBaseZ,
                dSpotDate,
                dBoosterData,
                dMeanReversionBase);
        BootstrappedVolsSingleCurrencyFixedIncome1F( NoticeDatesForeign,
                VolForeign,
                dNoticeDatesBis,
                dForeignDates,
                dForeignRates,
                dForeignRatesNoBasis,
                dStdDevForeignX,
                dStdDevForeignY,
                dStdDevForeignZ,
                dSpotDate,
                dBoosterData,
                dMeanReversionForeign);
    }
    else if(dIsSwaptionCalibrationWithBasis==0.0) // SwaptionCalibrationWithoutBasis // CDC ??? // Does not make sense but looks like this is the way they do it
    {
        BootstrappedVolsSingleCurrencyFixedIncome1F( NoticeDatesDomestic,
                VolDomestic,
                dNoticeDatesBis,
                dBaseDates,
                dBaseRatesNoBasis,
                dBaseRatesNoBasis,
                dStdDevBaseX,
                dStdDevBaseY,
                dStdDevBaseZ,
                dSpotDate,
                dBoosterData,
                dMeanReversionBase);
        BootstrappedVolsSingleCurrencyFixedIncome1F( NoticeDatesForeign,
                VolForeign,
                dNoticeDatesBis,
                dForeignDates,
                dForeignRatesNoBasis,
                dForeignRatesNoBasis,
                dStdDevForeignX,
                dStdDevForeignY,
                dStdDevForeignZ,
                dSpotDate,
                dBoosterData,
                dMeanReversionForeign);
    }
    else
    {
        throw("Unauthorised Input in Swaption Calibration Method");
    }


    DKMaille<double> dStripAsLattice;
    double LastDate=(dBoosterData.at(dBoosterData.rows()-1,1)-dSpotDate)/365.;

    dStripAsLattice=Create_Strip_Analytics( dNoticeDatesBis,
                                            dNumTimeLinesBeforeFirstNotice,
                                            dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                                            dNumTimeLinesPerYear,
                                            dOptimal,
                                            LastDate);

    TransformVolatilitiesSingle(NoticeDatesDomestic, dStripAsLattice, VolDomestic);
    TransformVolatilitiesSingle(NoticeDatesForeign, dStripAsLattice, VolForeign);


    SetSpotFXVolAnalytics( dSpotFXVolDatesTD,
                           dSpotFXVolTD,
                           dSpotDate,
                           dStripAsLattice,
                           dStripAsLattice,
                           VolDomestic,
                           VolForeign,
                           dMeanReversionBase,
                           dMeanReversionForeign,
                           dBaseForeignCorrelation,
                           dBaseSpotFXCorrelation,
                           dForeignSpotFXCorrelation,
                           dBoosterData);

    dStripDatesVolBase.resize(dStripAsLattice.entries());
    dStripDatesVolForeign.resize(dStripAsLattice.entries());
    dVolBase.resize(dStripAsLattice.entries());
    dVolForeign.resize(dStripAsLattice.entries());

    for(unsigned int ui=0;ui<dStripAsLattice.entries();ui++)
    {
        dStripDatesVolBase.at(ui)=365.*dStripAsLattice.at(ui)+dSpotDate;
        dStripDatesVolForeign.at(ui)=365.*dStripAsLattice.at(ui)+dSpotDate;
        dVolBase.at(ui)=VolDomestic.at(ui);
        dVolForeign.at(ui)=VolForeign.at(ui);
    }

    dStripDatesVolFX.resize(dSpotFXVolDatesTD.entries());
    dVolFX.resize(dSpotFXVolTD.entries());

    for(ui=0;ui<dSpotFXVolTD.entries();ui++)
    {
        dStripDatesVolFX.at(ui)=365.*dSpotFXVolDatesTD.at(ui)+dSpotDate;
        dVolFX.at(ui)=dSpotFXVolTD.at(ui);
    }

}

void Processing_3F_Data_DK(DKMaille<double> &dStripDatesVolFX,
                           DKMaille<double> &dStripDatesVolBase,
                           DKMaille<double> &dStripDatesVolForeign,
                           DKMaille<double> &dVolFX,
                           DKMaille<double> &dVolBase,
                           DKMaille<double> &dVolForeign,
                           double dNumTimeLinesBeforeFirstNotice,
                           double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                           double dNumTimeLinesPerYear,
                           double dSpotDate,
                           DKMaille<double> &dBaseDates,
                           DKMaille<double> &dBaseRates,
                           DKMaille<double> &dForeignDates,
                           DKMaille<double> &dForeignRates,
                           DKMaille<double> &dStdDevBaseX,
                           DKMaille<double> &dStdDevBaseY,
                           DKMaille2D<double> &dStdDevBaseZ,
                           double dMeanReversionBase,
                           DKMaille<double> &dStdDevForeignX,
                           DKMaille<double> &dStdDevForeignY,
                           DKMaille2D<double> &dStdDevForeignZ,
                           double dMeanReversionForeign,
                           DKMaille<double> dSpotFXVolDatesTD,
                           DKMaille<double> dSpotFXVolTD,
                           double dBaseForeignCorrelation,
                           double dBaseSpotFXCorrelation,
                           double dForeignSpotFXCorrelation,
                           double dOptimal,
                           DKMaille2D<double> &dBoosterData,
                           DKMaille<double> &dBaseRatesNoBasis,
                           DKMaille<double> &dForeignRatesNoBasis,
                           double dIsSwaptionCalibrationWithBasis)
{

    DKMaille<double> NoticeDates(dBoosterData.rows());

    for(unsigned int i=0;i<dBoosterData.rows();i++)
    {
        NoticeDates.at(i)=(dBoosterData.at(i,0)-dSpotDate)/365.;
    }


    DKMaille<double> NoticeDatesDomestic;
    DKMaille<double> VolDomestic;
    DKMaille<double> NoticeDatesForeign;
    DKMaille<double> VolForeign;

    DKMaille<double> dNoticeDatesBis;
    DKMaille<double> dNoticeDates_(dBoosterData.rows());
    for(unsigned int uiBD=0;uiBD<dBoosterData.rows();uiBD++)
        dNoticeDates_.at(uiBD)=(dBoosterData.at(uiBD,0)-dSpotDate)/365.;
    DKMaille<double> dLocalNoticeDates;
    unsigned int uiNoticeSize=0;
    for(uiBD=0;uiBD<dNoticeDates_.entries();uiBD++)
    {
        if(dNoticeDates_.at(uiBD)==0.)
            throw("Notice time is the same as evaluation time. Not authorised");
        // only value the future portion
        if(dNoticeDates_.at(uiBD)>0.)
        {
            if(dNoticeDates_.at(uiBD)!=(0.-dSpotDate)/365.)
            {
                dLocalNoticeDates.insert(dNoticeDates_.at(uiBD));
                uiNoticeSize++;
            }
        }
    }
    dNoticeDatesBis.resize(uiNoticeSize);
    for(uiBD=0;uiBD<uiNoticeSize;uiBD++)
    {
        dNoticeDatesBis.at(uiBD)=dLocalNoticeDates.at(uiBD);
    }



    if(dIsSwaptionCalibrationWithBasis==1.0) 
    {
        BootstrappedVolsSingleCurrencyFixedIncome1F_DK(NoticeDatesDomestic,
                VolDomestic,
                dNoticeDatesBis,
                dBaseDates,
                dBaseRates,
                dBaseRatesNoBasis,
                dStdDevBaseX,
                dStdDevBaseY,
                dStdDevBaseZ,
                dSpotDate,
                dBoosterData,
                dMeanReversionBase);

        BootstrappedVolsSingleCurrencyFixedIncome1F_DK(NoticeDatesForeign,
                VolForeign,
                dNoticeDatesBis,
                dForeignDates,
                dForeignRates,
                dForeignRatesNoBasis,
                dStdDevForeignX,
                dStdDevForeignY,
                dStdDevForeignZ,
                dSpotDate,
                dBoosterData,
                dMeanReversionForeign);
    }
    else if(dIsSwaptionCalibrationWithBasis==0.0) // SwaptionCalibrationWithoutBasis // CDC ??? // Does not make sense but looks like this is the way they do it
    {
        BootstrappedVolsSingleCurrencyFixedIncome1F_DK(NoticeDatesDomestic,
                VolDomestic,
                dNoticeDatesBis,
                dBaseDates,
                dBaseRatesNoBasis,
                dBaseRatesNoBasis,
                dStdDevBaseX,
                dStdDevBaseY,
                dStdDevBaseZ,
                dSpotDate,
                dBoosterData,
                dMeanReversionBase);

        BootstrappedVolsSingleCurrencyFixedIncome1F_DK(NoticeDatesForeign,
                VolForeign,
                dNoticeDatesBis,
                dForeignDates,
                dForeignRatesNoBasis,
                dForeignRatesNoBasis,
                dStdDevForeignX,
                dStdDevForeignY,
                dStdDevForeignZ,
                dSpotDate,
                dBoosterData,
                dMeanReversionForeign);
    }
    else
    {
        throw("Unauthorised Input in Swaption Calibration Method");
    }


    DKMaille<double> dStripAsLattice;
    double LastDate=(dBoosterData.at(dBoosterData.rows()-1,1)-dSpotDate)/365.;

    dStripAsLattice=Create_Strip_Analytics( dNoticeDatesBis,
                                            dNumTimeLinesBeforeFirstNotice,
                                            dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                                            dNumTimeLinesPerYear,
                                            dOptimal,
                                            LastDate);

    TransformVolatilitiesSingle(NoticeDatesDomestic, dStripAsLattice, VolDomestic);
    TransformVolatilitiesSingle(NoticeDatesForeign, dStripAsLattice, VolForeign);

    SetSpotFXVolAnalytics_DK(dSpotFXVolDatesTD,
                             dSpotFXVolTD,
                             dSpotDate,
                             dStripAsLattice,
                             dStripAsLattice,
                             VolDomestic,
                             VolForeign,
                             dMeanReversionBase,
                             dMeanReversionForeign,
                             dBaseForeignCorrelation,
                             dBaseSpotFXCorrelation,
                             dForeignSpotFXCorrelation,
                             dBoosterData,
                             dBaseDates,
                             dBaseRates,
                             dForeignDates,
                             dForeignRates);

    dStripDatesVolBase.resize(dStripAsLattice.entries());
    dStripDatesVolForeign.resize(dStripAsLattice.entries());
    dVolBase.resize(dStripAsLattice.entries());
    dVolForeign.resize(dStripAsLattice.entries());

    for(unsigned int ui=0;ui<dStripAsLattice.entries();ui++)
    {
        dStripDatesVolBase.at(ui)=365.*dStripAsLattice.at(ui)+dSpotDate;
        dStripDatesVolForeign.at(ui)=365.*dStripAsLattice.at(ui)+dSpotDate;
        dVolBase.at(ui)=VolDomestic.at(ui);
        dVolForeign.at(ui)=VolForeign.at(ui);
    }

    dStripDatesVolFX.resize(dSpotFXVolDatesTD.entries());
    dVolFX.resize(dSpotFXVolTD.entries());

    for(ui=0;ui<dSpotFXVolTD.entries();ui++)
    {
        dStripDatesVolFX.at(ui)=365.*dSpotFXVolDatesTD.at(ui)+dSpotDate;
        dVolFX.at(ui)=dSpotFXVolTD.at(ui);
    }

}







void SetSpotFXVolAnalytics( DKMaille<double> &dFXVolDates,
                            DKMaille<double> &dFXVol,
                            double dSpotDate,
                            DKMaille<double> &dStdDevBaseX,
                            DKMaille<double> &dStdDevForeignX,
                            DKMaille<double> &dStdDevBaseZ,
                            DKMaille<double> &dStdDevForeignZ,
                            double dMeanReversionBase,
                            double dMeanReversionForeign,
                            double dBaseForeignCorrelation,
                            double dBaseSpotFXCorrelation,
                            double dForeignSpotFXCorrelation,
                            DKMaille2D<double> &dBoosterData)
{

    if(dFXVolDates.at(0)==0.)
        throw("Unauthorised first date in spot FX vol");

    if(dFXVol.entries()!=dFXVolDates.entries()+2)
        throw("Unauthorised input for FX volatility");

    double dFinalMaturity=(dBoosterData.at(dBoosterData.rows()-1,2)-dSpotDate)/365.;

    double dCutOff=dFXVol.at(dFXVol.entries()-2);

    DKMaille<double> dNewVolStripDates;
    DKMaille<double> dNewVols;

    if(dCutOff==0.)
    {

        // no bootstrapping
        // Add anchor to the date strip

        DKMaille<double> dRawSpotFXVol;
        DKMaille<double> dStrip;

        dStrip.insert(0.);

        dRawSpotFXVol.insert(dFXVol.at(0));

        // Long dated historical Spot FX Vol

        double dLongDatedSpotFXVol=dFXVol.at(dFXVol.entries()-1);

        for(unsigned int ui=0;ui<dFXVolDates.entries();ui++)
        {
            dStrip.insert(dFXVolDates.at(ui));
            if(ui<dFXVolDates.entries()-1) dRawSpotFXVol.insert(dFXVol.at(ui+1));
            else dRawSpotFXVol.insert(dLongDatedSpotFXVol);
        }

        dStrip.insert(40.);
        dRawSpotFXVol.insert(dLongDatedSpotFXVol);
        dStrip.insert(50.);

        dRawSpotFXVol.insert(dLongDatedSpotFXVol);
        dNewVolStripDates=dStrip;

        dNewVols=dRawSpotFXVol;
    }
    else
    {
        double dCutOffBis=0.;

        //verification du positionnement du CutOff
        for(unsigned int ui=0;ui<dFXVolDates.entries();ui++)
        {
            if(fabs(dCutOff-dFXVolDates.at(ui))<=0.1)
            {
                dCutOff=dFXVolDates.at(ui);
                if(ui<dFXVolDates.entries()-1)
                    dCutOffBis=dFXVolDates.at(ui+1);
                else
                    dCutOffBis=40.;
                ui=dFXVolDates.entries();
            }
        }

        // si le CutOff n'est pas dans la liste des dates on sort
        if(dCutOffBis==0.)
            throw("Mauvais CutOff");

        DKMaille<double> dFXVolDatesPrime;
        DKMaille<double> dFXVolPrime;

        for(ui=0;ui<dFXVolDates.entries();ui++)
        {
            if(dFXVolDates.at(ui)<=dCutOff)
            {
                dFXVolDatesPrime.insert(dFXVolDates.at(ui));
                dFXVolPrime.insert(dFXVol.at(ui));
            }
        }

        DKMaille<double> dStartDate(dFXVolDatesPrime.entries());
        for(ui=0;ui<dFXVolDatesPrime.entries();ui++)
        {
            dStartDate.at(ui)=0.;
        }

        dNewVols=BootstrappingSpotFXVolatility3FTD( 0.,
                 dStartDate,
                 dFXVolDatesPrime,
                 dFXVolDatesPrime,
                 dFXVolPrime,
                 dStdDevBaseX,
                 dStdDevForeignX,
                 dStdDevBaseZ,
                 dStdDevForeignZ,
                 dMeanReversionBase,
                 dMeanReversionForeign,
                 dBaseSpotFXCorrelation,
                 dForeignSpotFXCorrelation,
                 dBaseForeignCorrelation);

        // Modification par D Kalafatis & S Pannetier.
        // si la vol spot FX sort negative a cause des donnees de marche
        // renvoyer a l'utilisateur un message d'erreur
        for(unsigned int uiJ=0;uiJ<dNewVols.entries();uiJ++)
        {
			if(dNewVols.at(uiJ)<=0.)
			{
				char msg[100];
				sprintf(msg,"Implied Spot FX Vol Is Negative. Check Input FX Vols and Input Swaption Vols %d %lf",uiJ,dNewVols.at(uiJ));
				throw(msg);
			}
		}

        dNewVolStripDates.insert(0.);
        for(ui=0;ui<dFXVolDatesPrime.entries();ui++)
            dNewVolStripDates.insert(dFXVolDatesPrime.at(ui));
        dNewVolStripDates.insert(dCutOffBis);
        dNewVolStripDates.insert(50);

        if(fabs(dCutOff-30.)<=0.1||fabs(dFXVol.at(dFXVol.entries()-1))<0.00001)
        {
            double dLastVol=dNewVols.at(dNewVols.entries()-1);
            dNewVols.insert(dLastVol);
            dNewVols.insert(dLastVol);
            dNewVols.insert(dLastVol);
        }
        else if(fabs(dFXVol.at(dFXVol.entries()-1))<0.011)
        {
            double dLastVol=SpotFXVolIntegral(dNewVolStripDates,dNewVols,0.,dCutOff);
            dLastVol=sqrt(dLastVol/dCutOff);
            dNewVols.insert(dLastVol);
            dNewVols.insert(dLastVol);
            dNewVols.insert(dLastVol);
        }
        else
        {
            dNewVols.insert(dFXVol.at(dFXVol.entries()-1));
            dNewVols.insert(dFXVol.at(dFXVol.entries()-1));
            dNewVols.insert(dFXVol.at(dFXVol.entries()-1));
        }
    }

    dFXVolDates.resize(dNewVolStripDates.entries());
    dFXVol.resize(dNewVolStripDates.entries());

    for(unsigned int ui=0;ui<dNewVolStripDates.entries();ui++)
    {
        dFXVolDates.at(ui)=dNewVolStripDates.at(ui);
        dFXVol.at(ui)=dNewVols.at(ui);
    }

}


void SetSpotFXVolAnalytics_DK(DKMaille<double> &dFXVolDates,
                              DKMaille<double> &dFXVol,
                              double dSpotDate,
                              DKMaille<double> &dStdDevBaseX,
                              DKMaille<double> &dStdDevForeignX,
                              DKMaille<double> &dStdDevBaseZ,
                              DKMaille<double> &dStdDevForeignZ,
                              double dMeanReversionBase,
                              double dMeanReversionForeign,
                              double dBaseForeignCorrelation,
                              double dBaseSpotFXCorrelation,
                              double dForeignSpotFXCorrelation,
                              DKMaille2D<double> &dBoosterData,
                              DKMaille<double> &dCurveDates,
                              DKMaille<double> &dCurve,
                              DKMaille<double> &dCurveForeignDates,
                              DKMaille<double> &dCurveForeign)
{

    if(dFXVolDates.at(0)==0.)
        throw("Unauthorised first date in spot FX vol");

    if(dFXVol.entries()!=dFXVolDates.entries()+2)
        throw("Unauthorised input for FX volatility");

    double dFinalMaturity=(dBoosterData.at(dBoosterData.rows()-1,2)-dSpotDate)/365.;

    double dCutOff=dFXVol.at(dFXVol.entries()-2);

    DKMaille<double> dNewVolStripDates;
    DKMaille<double> dNewVols;

    if(dCutOff==0.)
    {

        // no bootstrapping
        // Add anchor to the date strip

        DKMaille<double> dRawSpotFXVol;
        DKMaille<double> dStrip;

        dStrip.insert(0.);

        dRawSpotFXVol.insert(dFXVol.at(0));

        // Long dated historical Spot FX Vol

        double dLongDatedSpotFXVol=dFXVol.at(dFXVol.entries()-1);

        for(unsigned int ui=0;ui<dFXVolDates.entries();ui++)
        {
            dStrip.insert(dFXVolDates.at(ui));
            if(ui<dFXVolDates.entries()-1) dRawSpotFXVol.insert(dFXVol.at(ui+1));
            else dRawSpotFXVol.insert(dLongDatedSpotFXVol);
        }

        dStrip.insert(40.);
        dRawSpotFXVol.insert(dLongDatedSpotFXVol);
        dStrip.insert(50.);

        dRawSpotFXVol.insert(dLongDatedSpotFXVol);
        dNewVolStripDates=dStrip;

        dNewVols=dRawSpotFXVol;
    }
    else
    {
        double dCutOffBis=0.;

        //verification du positionnement du CutOff
        for(unsigned int ui=0;ui<dFXVolDates.entries();ui++)
        {
            if(fabs(dCutOff-dFXVolDates.at(ui))<=0.1)
            {
                dCutOff=dFXVolDates.at(ui);
                if(ui<dFXVolDates.entries()-1)
                    dCutOffBis=dFXVolDates.at(ui+1);
                else
                    dCutOffBis=40.;
                ui=dFXVolDates.entries();
            }
        }

        // si le CutOff n'est pas dans la liste des dates on sort
        if(dCutOffBis==0.)
            throw("Mauvais CutOff");

        DKMaille<double> dFXVolDatesPrime;
        DKMaille<double> dFXVolPrime;

        for(ui=0;ui<dFXVolDates.entries();ui++)
        {
            if(dFXVolDates.at(ui)<=dCutOff)
            {
                dFXVolDatesPrime.insert(dFXVolDates.at(ui));
                dFXVolPrime.insert(dFXVol.at(ui));
            }
        }

        DKMaille<double> dStartDate(dFXVolDatesPrime.entries());
        for(ui=0;ui<dFXVolDatesPrime.entries();ui++)
        {
            dStartDate.at(ui)=0.;
        }

        dNewVols=BootstrappingSpotFXVolatility3FTD_DK(0.,
                 dStartDate,
                 dFXVolDatesPrime,
                 dFXVolDatesPrime,
                 dFXVolPrime,
                 dStdDevBaseX,
                 dStdDevForeignX,
                 dStdDevBaseZ,
                 dStdDevForeignZ,
                 dMeanReversionBase,
                 dMeanReversionForeign,
                 dBaseSpotFXCorrelation,
                 dForeignSpotFXCorrelation,
                 dBaseForeignCorrelation,
                 dCurveDates,
                 dCurve,
                 dCurveForeignDates,
                 dCurveForeign
                                                     );

        // Modification par D Kalafatis & S Pannetier.
        // si la vol spot FX sort negative a cause des donnees de marche
        // renvoyer a l'utilisateur un message d'erreur
        for(unsigned int uiJ=0;uiJ<dNewVols.entries();uiJ++)
            if(dNewVols.at(uiJ)<=0.) throw("Implied Spot FX Vol Is Negative. Check Input FX Vols and Input Swaption Vols");


        dNewVolStripDates.insert(0.);
        for(ui=0;ui<dFXVolDatesPrime.entries();ui++)
            dNewVolStripDates.insert(dFXVolDatesPrime.at(ui));
        dNewVolStripDates.insert(dCutOffBis);
        dNewVolStripDates.insert(50);


        if(fabs(dCutOff-30.)<=0.1||fabs(dFXVol.at(dFXVol.entries()-1))<0.00001)
        {
            double dLastVol=dNewVols.at(dNewVols.entries()-1);
            dNewVols.insert(dLastVol);
            dNewVols.insert(dLastVol);
            dNewVols.insert(dLastVol);
        }
        else if(fabs(dFXVol.at(dFXVol.entries()-1))<0.011)
        {
            double dLastVol=SpotFXVolIntegral(dNewVolStripDates,dNewVols,0.,dCutOff);
            dLastVol=sqrt(dLastVol/dCutOff);
            dNewVols.insert(dLastVol);
            dNewVols.insert(dLastVol);
            dNewVols.insert(dLastVol);
        }
        else
        {
            dNewVols.insert(dFXVol.at(dFXVol.entries()-1));
            dNewVols.insert(dFXVol.at(dFXVol.entries()-1));
            dNewVols.insert(dFXVol.at(dFXVol.entries()-1));
        }

    }

    dFXVolDates.resize(dNewVolStripDates.entries());
    dFXVol.resize(dNewVolStripDates.entries());

    for(unsigned int ui=0;ui<dNewVolStripDates.entries();ui++)
    {
        dFXVolDates.at(ui)=dNewVolStripDates.at(ui);
        dFXVol.at(ui)=dNewVols.at(ui);
    }

}












DKMaille<double> Create_Strip_Analytics( DKMaille<double> NoticeDates,
        int dNumTimeLinesBeforeFirstNotice,
        int dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
        int dNumTimeLinesPerYear,
        double dOptimalDate,
        double LastDate)
{
    // determination du nombre de pas avant la premiere notice
    int uiSlices1;
    if(NoticeDates.at(0)>1.1)
        uiSlices1=(int)(dNumTimeLinesBeforeFirstNotice*NoticeDates.at(0));
    else
        uiSlices1=dNumTimeLinesBeforeFirstNotice;
    if(uiSlices1<4) uiSlices1=4;


    // determination de la NoticeDate la plus proche de l'OptimalDate
    int Optimal=-1;
    for(int ui=0;ui<NoticeDates.entries();ui++)
    {
        if(NoticeDates.at(ui)>dOptimalDate)
        {
            if(ui==0)
            {
                Optimal=0;
            }
            else
            {
                if((NoticeDates.at(ui)-dOptimalDate)>(dOptimalDate-NoticeDates.at(ui-1)))
                    Optimal=ui-1;
                else
                    Optimal=ui;
            }
            ui=NoticeDates.entries();
        }

    }

    // si Optimal est encore a -1  c que dOptimalDate est superieure a la derniere notice
    if(Optimal==-1)
    {
        Optimal=NoticeDates.entries()-1;
    }


    bool annuel=true;
    if(NoticeDates.entries()>1 && NoticeDates.at(1)-NoticeDates.at(0)<0.75)
        annuel=false;

    // uiSlices2 est le nombre de slices entre chaque notice avant la date optimale
    int uiSlices2=0;
    if(Optimal!=0)
    {
        if(annuel)
        {
            uiSlices2=dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal;
        }
        else
        {
            if(dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal%2 == 0)
            {
                uiSlices2=dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal/2;
            }
            else
            {
                uiSlices2=(dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal+1)/2;
            }
        }
    }

    // uiSlices3 est le nombre de slices entre chaque notice apres la date optimale
    int uiSlices3=0;
    if(Optimal!=NoticeDates.entries()-1)
    {
        if(annuel)
        {
            uiSlices3=dNumTimeLinesPerYear;
        }
        else
        {
            if(dNumTimeLinesPerYear%2 == 0)
            {
                uiSlices3=dNumTimeLinesPerYear/2;
            }
            else
            {
                uiSlices3=(dNumTimeLinesPerYear+1)/2;
            }
        }
    }

    int Slices2=uiSlices2*Optimal;
    int Slices3=uiSlices3*(NoticeDates.entries()-Optimal-1);
    int NbPasTotal=uiSlices1+Slices2+Slices3;

    // on ajoute a NbPasTotal le nombre de pas apres la derniere notice
    // Pour stabiliser le reseau on en met dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal par annee si Optimal est a la derniere notice
    // Sinon on en met dNumTimeLinesPerYear par annee
    int NbPasAfterLastNotice=0;
    if(Optimal!=NoticeDates.entries()-1)
    {
        NbPasAfterLastNotice=(int) (dNumTimeLinesPerYear*(LastDate-NoticeDates.at(NoticeDates.entries()-1))+0.5);
    }
    else
    {
        NbPasAfterLastNotice=(int) (dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal*(LastDate-NoticeDates.at(NoticeDates.entries()-1))+0.5);
    }
    NbPasAfterLastNotice=MAX(NbPasAfterLastNotice,1);
    NbPasTotal+=NbPasAfterLastNotice;



    // on remplit le dStrip
    DKMaille<double> dStrip;
    dStrip.resize(1+NbPasTotal);

    dStrip.at(0)=0.;
    for(ui=1;ui<=uiSlices1;ui++)
    {
        dStrip.at(ui)=(double)(ui)/(double)(uiSlices1)*NoticeDates.at(0);
    }

    for(ui=0;ui<Optimal;ui++)
    {
        for(int j=1;j<=uiSlices2;j++)
        {
            dStrip.at(uiSlices1+ui*uiSlices2+j)
            =NoticeDates.at(ui)+(double)(j)/(double)(uiSlices2)*(NoticeDates.at(ui+1)-NoticeDates.at(ui));
        }
    }

    for(ui=Optimal;ui<NoticeDates.entries()-1;ui++)
    {
        for(int j=1;j<=uiSlices3;j++)
        {
            dStrip.at(uiSlices1+Slices2+(ui-Optimal)*uiSlices3+j)
            =NoticeDates.at(ui)+(double)(j)/(double)(uiSlices3)*(NoticeDates.at(ui+1)-NoticeDates.at(ui));
        }
    }

    for(ui=1;ui<=NbPasAfterLastNotice;ui++)
    {
        dStrip.at(uiSlices1+Slices2+Slices3+ui)=NoticeDates.at(NoticeDates.entries()-1)
                                                +(double)(ui)/(double)(NbPasAfterLastNotice)*(LastDate-NoticeDates.at(NoticeDates.entries()-1));
    }

    return dStrip;
}




void BootstrappedVolsSingleCurrencyFixedIncome1F(DKMaille<double> &dNewVolStripDates,
        DKMaille<double> &dOutputStrip,
        DKMaille<double> &dNoticeDatesBis,
        DKMaille<double> &dDates,
        DKMaille<double> &dRates,
        DKMaille<double> &dRatesNoBasis,
        DKMaille<double> &dStdDevX,
        DKMaille<double> &dStdDevY,
        DKMaille2D<double> &dStdDevZ,
        double dSpotDate,
        DKMaille2D<double> &dBoosterData,
        double dMeanReversion)

{
    // Sort Notice Dates for domestic currency bootstrapping
    DKMaille<double> dVolStripDates;
    SelectBootstrappingDates(dVolStripDates,dNoticeDatesBis,dStdDevX);
    dNewVolStripDates.resize(dVolStripDates.entries()+2);
    DKMaille<double> dSwapStartDates(dVolStripDates.entries());
    DKMaille<double> dSwapEndDates(dVolStripDates.entries());
    dNewVolStripDates.at(0)=dSpotDate;
    double dFinalMaturity=dBoosterData.at(dBoosterData.rows()-1,2);
    for(unsigned int uiD=0;uiD<dVolStripDates.entries();uiD++)
    {
        dNewVolStripDates.at(uiD+1)=dSpotDate+dVolStripDates.at(uiD)*365.;
        // dNoticePeriod is set to 0 days for cash settlement in the bootstrapping procedure. It will be questionable to change it
        dSwapStartDates.at(uiD)=dNewVolStripDates.at(uiD+1)+0.;
        // Only final maturity diagonal is allowed for calibration purposes for the moment
        // Change the below line to make it accept any date
        dSwapEndDates.at(uiD)=dFinalMaturity;
    }
    dNewVolStripDates.at(dNewVolStripDates.entries()-1)=dFinalMaturity;
    DKMaille<double> dModelParameters(2);
    dModelParameters.at(0)=1.;
    dModelParameters.at(1)=dMeanReversion;


    DKMaille<double> dVolStrip;
    double dIsConstant=isConstant(dStdDevZ);
    if(dIsConstant==0.)
    {
        dVolStrip=Bootstrapping_VFDK_HW1To3F(dStdDevZ,
                                             dStdDevX,
                                             dStdDevY,
                                             dDates,
                                             dRatesNoBasis,
                                             dRates,
                                             dNewVolStripDates,
                                             dSwapStartDates,
                                             dSwapEndDates,
                                             dModelParameters,
                                             dNewVolStripDates.at(0));
    }
    else
    {
        for(unsigned int uiVol=0;uiVol<dNewVolStripDates.entries()-2;uiVol++)
        {
            dVolStrip.insert(dIsConstant);
        }
    }
    if(dVolStrip.entries()!=dNewVolStripDates.entries()-2)
        throw("Error in input bootstrapping dates");
    dOutputStrip.resize(dVolStrip.entries()+2);
    for(uiD=0;uiD<dVolStrip.entries();uiD++)
        dOutputStrip.at(uiD)=dVolStrip.at(uiD);
    // fill the last two holders
    dOutputStrip.at(dVolStrip.entries())=dVolStrip.at(dVolStrip.entries()-1);
    dOutputStrip.at(dVolStrip.entries()+1)=dVolStrip.at(dVolStrip.entries()-1);

    // Write on top of the vol holders
    for(uiD=0;uiD<dNewVolStripDates.entries();uiD++)
        dNewVolStripDates.at(uiD)=(dNewVolStripDates.at(uiD)-dSpotDate)/365.;
}

void BootstrappedVolsSingleCurrencyFixedIncome1F_DK(DKMaille<double> &dNewVolStripDates,
        DKMaille<double> &dOutputStrip,
        DKMaille<double> &dNoticeDatesBis,
        DKMaille<double> &dDates,
        DKMaille<double> &dRates,
        DKMaille<double> &dRatesNoBasis,
        DKMaille<double> &dStdDevX,
        DKMaille<double> &dStdDevY,
        DKMaille2D<double> &dStdDevZ,
        double dSpotDate,
        DKMaille2D<double> &dBoosterData,
        double dMeanReversion)

{
    // Sort Notice Dates for domestic currency bootstrapping
    DKMaille<double> dVolStripDates;
    SelectBootstrappingDates(dVolStripDates,dNoticeDatesBis,dStdDevX);
    dNewVolStripDates.resize(dVolStripDates.entries()+2);
    DKMaille<double> dSwapStartDates(dVolStripDates.entries());
    DKMaille<double> dSwapEndDates(dVolStripDates.entries());
    dNewVolStripDates.at(0)=dSpotDate;
    double dFinalMaturity=dBoosterData.at(dBoosterData.rows()-1,2);
    for(unsigned int uiD=0;uiD<dVolStripDates.entries();uiD++)
    {
        dNewVolStripDates.at(uiD+1)=dSpotDate+dVolStripDates.at(uiD)*365.;
        // dNoticePeriod is set to 0 days for cash settlement in the bootstrapping procedure. It will be questionable to change it
        dSwapStartDates.at(uiD)=dNewVolStripDates.at(uiD+1)+0.;
        // Only final maturity diagonal is allowed for calibration purposes for the moment
        // Change the below line to make it accept any date
        dSwapEndDates.at(uiD)=dFinalMaturity;
    }
    dNewVolStripDates.at(dNewVolStripDates.entries()-1)=dFinalMaturity;
    DKMaille<double> dModelParameters(2);
    dModelParameters.at(0)=1.;
    dModelParameters.at(1)=dMeanReversion;


    DKMaille<double> dVolStrip;
    double dIsConstant=isConstant(dStdDevZ);
    if(dIsConstant==0.)
    {
        if(g_SP==0)
        {
                throw("Unauthorised input !");

        }
        else
        {
            dVolStrip=Bootstrapping_DK3F_Numerical_SP(dStdDevZ,
                      dStdDevX,
                      dStdDevY,
                      dDates,
                      dRatesNoBasis,
                      dRates,
                      dNewVolStripDates,
                      dSwapStartDates,
                      dSwapEndDates,
                      dModelParameters,
                      dNewVolStripDates.at(0),
                      4,
                      3,
                      1,
                      10.,
                      0.01);
        }
    }
    else
    {
        for(unsigned int uiVol=0;uiVol<dNewVolStripDates.entries()-2;uiVol++)
        {
            dVolStrip.insert(dIsConstant);
        }
    }
    if(dVolStrip.entries()!=dNewVolStripDates.entries()-2)
        throw("Error in input bootstrapping dates");
    dOutputStrip.resize(dVolStrip.entries()+2);
    for(uiD=0;uiD<dVolStrip.entries();uiD++)
        dOutputStrip.at(uiD)=dVolStrip.at(uiD);
    // fill the last two holders
    dOutputStrip.at(dVolStrip.entries())=dVolStrip.at(dVolStrip.entries()-1);
    dOutputStrip.at(dVolStrip.entries()+1)=dVolStrip.at(dVolStrip.entries()-1);

    // Write on top of the vol holders
    for(uiD=0;uiD<dNewVolStripDates.entries();uiD++)
        dNewVolStripDates.at(uiD)=(dNewVolStripDates.at(uiD)-dSpotDate)/365.;
}




void FromForwardToSpotVol3F(double dSpotDate,
                            DKMaille<double> dBaseDiscountCurveDates,
                            DKMaille<double> dBaseDiscountCurveRates,
                            DKMaille<double> dBaseAdjustedDiscountCurveRates,
                            DKMaille<double> dForeignDiscountCurveDates,
                            DKMaille<double> dForeignDiscountCurveRates,
                            DKMaille<double> dForeignAdjustedDiscountCurveRates,
                            DKMaille<double> dDomesticVolX,
                            DKMaille<double> dDomesticVolY,
                            DKMaille2D<double> dDomesticVolZ,
                            DKMaille<double> dForeignVolX,
                            DKMaille<double> dForeignVolY,
                            DKMaille2D<double> dForeignVolZ,
                            DKMaille<double> dForwardFXStrip,
                            DKMaille<double> dForwardFXVol,
                            double dDomesticMeanReversion,
                            double dForeignMeanReversion,
                            double dCorrelationSpotFXDomestic,
                            double dCorrelationSpotFXForeign,
                            double dCorrelationForeignDomestic,
                            double dCurrencyPair1, // 0=jpy 1=usd 2=aud 3=eur
                            double dCurrencyPair2, // 0=jpy 1=usd 2=aud 3=eur
                            DKMaille<double> dNoticeDates,
                            DKMaille<double> dFXCouponResetDates,
                            DKMaille<double> dFXCouponPaymentDates,
                            double dIsSwaptionCalibrationWithBasis, 
                            DKMaille<double> &dSpotFXStrip,
                            DKMaille<double> &dSpotFXVol)
{
    DKMaille<double> dStripDatesVolBase;
    DKMaille<double> dStripDatesVolForeign;
    DKMaille<double> dVolBase;
    DKMaille<double> dVolForeign;

    int maxLength = MAX(MAX(dNoticeDates.entries(), dFXCouponResetDates.entries()),
                        dFXCouponPaymentDates.entries());
    DKMaille2D<double> dBoosterData(maxLength, 3);
    for (int i = 0; i < maxLength; i++)
    {
        if ( i < dNoticeDates.entries() )
            dBoosterData.at(i,0) = dNoticeDates.at(i);
        else
            dBoosterData.at(i,0)=0.0;

        if ( i < dFXCouponResetDates.entries() )
            dBoosterData.at(i,1) = dFXCouponResetDates.at(i);
        else
            dBoosterData.at(i,1) = 0.;

        if ( i < dFXCouponPaymentDates.entries() )
            dBoosterData.at(i,2) = dFXCouponPaymentDates.at(i);
        else
            dBoosterData.at(i,2) = 0.;
    }


    Processing_3F_Data( dSpotFXStrip,
                        dStripDatesVolBase,
                        dStripDatesVolForeign,
                        dSpotFXVol,
                        dVolBase,
                        dVolForeign,
                        4,  /////////////// a determiner
                        3,  /////////////// a determiner
                        1,  /////////////// a determiner
                        dSpotDate,
                        dBaseDiscountCurveDates,
                        dBaseAdjustedDiscountCurveRates,
                        dForeignDiscountCurveDates,
                        dForeignAdjustedDiscountCurveRates,
                        dDomesticVolX,
                        dDomesticVolY,
                        dDomesticVolZ,
                        dDomesticMeanReversion,
                        dForeignVolX,
                        dForeignVolY,
                        dForeignVolZ,
                        dForeignMeanReversion,
                        dForwardFXStrip,
                        dForwardFXVol,
                        dCorrelationForeignDomestic,
                        dCorrelationSpotFXDomestic,
                        dCorrelationSpotFXForeign,
                        10, /////////////// a determiner
                        dBoosterData,
                        dBaseDiscountCurveRates,
                        dForeignDiscountCurveRates,
                        dIsSwaptionCalibrationWithBasis);

}
// FromForwardToSpotVol3F




