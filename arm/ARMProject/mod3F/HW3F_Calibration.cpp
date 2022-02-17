/*
 * $Log: HW3F_Calibration.cpp,v $
 * Revision 1.9  2004/05/28 12:35:23  mab
 * version without printfun Print out function
 *
 * Revision 1.8  2004/05/27 13:40:17  mab
 * Correction
 *
 * Revision 1.5  2004/05/25 13:30:41  mab
 * Added : Calibration control
 *
 * Revision 1.4  2004/05/18 13:59:54  mab
 * corrections in Calibration
 *
 * Revision 1.3  2004/05/04 13:28:59  mcampet
 * MC from SP modif options e04jbc
 *
 * Revision 1.2  2004/04/29 07:29:22  mcampet
 *  MC Updates from SP
 *
 * Revision 1.1  2004/04/22 12:08:20  mcampet
 * Initial revision
 *
 */

/*----------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef WIN32

#include <float.h>

#include <sys/timeb.h>
#include "nag.h"
#include "nage04.h"

#endif


#include "hw_vfdk_analytics.h"

#include "HW3F_Calibration.h"

#include "expt.h"

#ifndef TRUE
#define TRUE 1
#endif


/*----------------------------------------------------------------------------*/

HW3F_Calibration_Params dHW3F_Calibration_Params;


#ifdef WIN32

typedef void (NAG_CALL *OBJFUN)(long n, double* x, double* f, 
                                double* g, Nag_Comm* comm);


typedef void (NAG_CALL *CONFUN)(long a, long b, long* c, double* d,
                                double* e, double* f, Nag_Comm* comm);



void NAG_CALL confun(long n, long ncnlin, long* needc, double* x,
                     double* conf, double* conjac, Nag_Comm* comm)
{
    throw("confun: You should never come here in e04ucc");
}



double ComputeNorm(long size, double* X1, double* X2)
{
    double norm = 0.0;


	for (long i; i < size; i++)
	{
		norm += fabs(X1[i]-X2[i])*fabs(X1[i]-X2[i]);
	}

	norm = sqrt(norm);

	return(norm);
}



void NAG_CALL printfun(const Nag_Search_State* st, Nag_Comm* comm)
{

     dHW3F_Calibration_Params.iterCounter++; // "Global" iterations counter

     if ( dHW3F_Calibration_Params.iterCounter >= CONSTANT_X_ITER )
     {
        double curNorm = ComputeNorm(st->n,
                                     st->x,
                                     dHW3F_Calibration_Params.X_PREC);

        if ( curNorm <= THE_3F_PRECISION )
        {
           // We already "converged" so:
           // generate a user exit from NAg

           comm->flag = -1000; // Thus NE_USER_STOP will be raised!
        }
        else
        {
           memcpy(dHW3F_Calibration_Params.X_PREC, st->x,
                  sizeof(double)*st->n);
        }
     }
     else
     {
        memcpy(dHW3F_Calibration_Params.X_PREC, st->x,
               sizeof(double)*st->n);
     }
}



void NAG_CALL objfun(long n, double* x, double* objf, double* g, 
					 Nag_Comm* comm)
{
    *objf = HW3F_Calibration_Function_TO_Minimize(x);
}

#endif



double GetSwaptionData(double dXX, 
                       double dYY,
                       DKMaille<double> &dSwaptionDates,
                       DKMaille<double> &dAccrualPeriods,
                       DKMaille<double> &dBasisSwaptionDates,
                       DKMaille<double> &dBasisAccrualPeriods,
                       DKMaille<double> &dBasis,
                       double dFixedPaymentFrequency,
                       double dFloatPaymentFrequency)
{
    // find the size of the arrays
    unsigned int uiStart=0; unsigned int uiEnd=-1; unsigned int ui=0; unsigned um=0;
    unsigned int uiSize= (unsigned int)(dFixedPaymentFrequency *(dXX)+1);
    dSwaptionDates.resize(uiSize);
    dAccrualPeriods.resize(uiSize);

    // Artificially set basis data

    dBasisSwaptionDates.resize(uiSize-1);
    dBasisAccrualPeriods.resize(uiSize-1);
    dBasis.resize(uiSize-1);

    for(ui=0;ui<uiSize;ui++)
    {
        dSwaptionDates.at(ui)=dHW3F_Calibration_Params.SpotDate+2.+dYY*365.
                               +(double)(ui)*365./2.;

        if(ui==0) 
          dAccrualPeriods.at(ui)=0.0;

        else dAccrualPeriods.at(ui)=(dSwaptionDates.at(ui)-dSwaptionDates.at(ui-1))/365.;

        if(ui!=0)
        {
            dBasisSwaptionDates.at(ui-1)=dSwaptionDates.at(ui);
            dBasisAccrualPeriods.at(ui-1)=0.0;
            dBasis.at(ui-1)=0.0;
        }
    }

    return 1.0;
}



double HW3F_Calibration_Function_TO_Minimize(double* dLocalModelParameters)
{

    DKMaille<double> dSwaptionDates;
    DKMaille<double> dAccrualPeriods;
    DKMaille<double> dBasisSwaptionDates;
    DKMaille<double> dBasisAccrualPeriods;
    DKMaille<double> dBasis;
    DKMaille<double> dSwaptionDates2;
    DKMaille<double> dAccrualPeriods2;
    DKMaille<double> dBasisSwaptionDates2;
    DKMaille<double> dBasisAccrualPeriods2;
    DKMaille<double> dBasis2;
    DKMaille<double> dModelParameters(dHW3F_Calibration_Params.ModelParameters.entries());
    DKMaille<double> dVolStrip(dHW3F_Calibration_Params.NoticeDatesBis.entries());
    DKMaille<double> dVolStripDates(dHW3F_Calibration_Params.NoticeDatesBis.entries());

    double dFunctionValue=0.;
    double dAbsoluteVolFromSurface;
    double dStrike;
    double dStraddle;
    double dMarketStraddle;
    double dImpliedCorrelation;
    double dReturn;
    double dSurfaceWeight;


    unsigned int ui=0;
    for(ui=0;ui<dHW3F_Calibration_Params.NoticeDatesBis.entries();ui++)
        dVolStripDates.at(ui)=dHW3F_Calibration_Params.NoticeDatesBis.at(ui);

    dHW3F_Calibration_Params.uiIterations+=1;

    unsigned uk=0;

    for(ui=0;ui<dModelParameters.entries();ui++)
    {
        if(dHW3F_Calibration_Params.WhichModelParametersToCalibrate.at(ui)==1.)
        {
            dModelParameters.at(ui)=dLocalModelParameters[uk];
            uk++;
        }
        else
        {
            dModelParameters.at(ui)=dHW3F_Calibration_Params.ModelParameters.at(ui);
        }
    }

    DKMaille<double> dOutputVolStrip;

    try
    {
        dOutputVolStrip=Bootstrapping_VFDK_HW1To3F(dHW3F_Calibration_Params.InputSurface_AbsVols,
                                                    dHW3F_Calibration_Params.InputArray_X_AbsVols,
                                                    dHW3F_Calibration_Params.InputArray_Y_AbsVols,
                                                    dHW3F_Calibration_Params.ZCDates,
                                                    dHW3F_Calibration_Params.ZCRates,
                                                    dHW3F_Calibration_Params.ZCRates,
                                                    dHW3F_Calibration_Params.NoticeDatesBis,
                                                    dHW3F_Calibration_Params.SwapStartDatesBis,
                                                    dHW3F_Calibration_Params.SwapEndDatesBis,
                                                    dModelParameters,
                                                    dHW3F_Calibration_Params.SpotDate);
    }

    catch(char* pszMsg)
    {
        dHW3F_Calibration_Params.Plantage++;

        if (strcmp(pszMsg,"Zero or Negative Vol in bootstrapping")==0.)
           return 1.e10;
        else if(strcmp(pszMsg,"Implied Negative Standard Deviation")==0.)
            return 1.e10;
        else
            throw(pszMsg);
    }
   
	catch(...)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
                           "Bootstrapping_VFDK_HW1To3F: Unrecognized failure");
    }

    for(ui=0;ui<dOutputVolStrip.entries();ui++)
        dVolStrip.at(ui)=dOutputVolStrip.at(ui);
    dVolStrip.at(dOutputVolStrip.entries())=dVolStrip.at(dOutputVolStrip.entries()-1);
    dVolStrip.at(dOutputVolStrip.entries()+1)=dVolStrip.at(dOutputVolStrip.entries()-1);

    unsigned int uiK=0;
    unsigned int uiM=0;
    double dXX=0.;
    double dYY=0.;

    double dFixedPaymentFrequency=2.;
    double dFloatPaymentFrequency=2.;

    
    double dSum=0.;

    try 
    {
        //  Create swaption matrix loop
        for(uiK=0;uiK<dHW3F_Calibration_Params.InputArray_X_AbsVols.entries();uiK++)
        {
            for(uiM=0;uiM<dHW3F_Calibration_Params.InputArray_Y_AbsVols.entries();uiM++)
            {
                // Create swaption dates
                // dXX is option expiry in years
                dXX = dHW3F_Calibration_Params.InputArray_X_AbsVols[uiK];
                // dYY is swap tenor in years
                dYY = dHW3F_Calibration_Params.InputArray_Y_AbsVols[uiM];

                dReturn = GetSwaptionData(dYY,dXX,dSwaptionDates,dAccrualPeriods,
                                          dBasisSwaptionDates,dBasisAccrualPeriods,
                                          dBasis,dFixedPaymentFrequency,
                                          dFloatPaymentFrequency);

                dStrike = GetSwapRate(dHW3F_Calibration_Params.SpotDate,dSwaptionDates,
                                      dAccrualPeriods,dHW3F_Calibration_Params.ZCDates,
                                      dHW3F_Calibration_Params.ZCRates,
                                      dHW3F_Calibration_Params.ZCRates);

                dAbsoluteVolFromSurface = dHW3F_Calibration_Params.InputSurface_AbsVols.at(uiK,uiM) ;

                dMarketStraddle=GetPriceFromAbsVol(dHW3F_Calibration_Params.SpotDate,
                                                   dHW3F_Calibration_Params.SpotDate,
                                                   dSwaptionDates,
                                                   dAccrualPeriods,
                                                   2./365.,
                                                   dAbsoluteVolFromSurface,
                                                   dHW3F_Calibration_Params.ZCDates,
                                                   dHW3F_Calibration_Params.ZCRates,
                                                   dHW3F_Calibration_Params.ZCRates);
                

                dSurfaceWeight = dHW3F_Calibration_Params.InputSurface_AbsVolsWeights.at(uiK,uiM);
                dStraddle=0.;

                if(dSurfaceWeight>0.0)
                {
                    dStraddle = 0.5*(SwaptionPrice_VFDK_HW1To3F( dXX,
                                                                 dYY,
                                                                 2.,
                                                                 dStrike,
                                                                 0.,
                                                                 dHW3F_Calibration_Params.ZCDates,
                                                                 dHW3F_Calibration_Params.ZCRates,
                                                                 dHW3F_Calibration_Params.ZCRates,
                                                                 dVolStripDates,
                                                                 dVolStrip,
                                                                 dModelParameters,
                                                                 dHW3F_Calibration_Params.SpotDate)
                                  +SwaptionPrice_VFDK_HW1To3F(dXX,
                                                              dYY,
                                                              2.,
                                                              dStrike,
                                                              1.,
                                                              dHW3F_Calibration_Params.ZCDates,
                                                              dHW3F_Calibration_Params.ZCRates,
                                                              dHW3F_Calibration_Params.ZCRates,
                                                              dVolStripDates,
                                                              dVolStrip,
                                                              dModelParameters,
                                                              dHW3F_Calibration_Params.SpotDate));
                }

                dFunctionValue+= dSurfaceWeight*(dStraddle-dMarketStraddle)
                                 *(dStraddle-dMarketStraddle);

                dSum+=dSurfaceWeight;
            }
        }

        dFunctionValue/=dSum;

        dSum=0.;

        // Add fwd correlation constraint
        double dFwdCorrelationFunctionValue=0.;
        //  Create fwd correlation array loop
        for(uiK=0;uiK<dHW3F_Calibration_Params.InputArray_X_FwdCorrelation_2_20.entries();uiK++)
        {
            // Create swaption dates
            // dXX is swap tenor of rate 1 in years
            dXX = 2.;
            // dYY is swap tenor of rate 2 in years
            dYY = 20.;

            double dSurfaceWeight = dHW3F_Calibration_Params.InputArray_FwdCorrelation_2_20Weights.at(uiK);

            // option expiry
            double dZZ=dHW3F_Calibration_Params.InputArray_X_FwdCorrelation_2_20.at(uiK);


            if(dSurfaceWeight>0.0)
            {
                double dMarketCorrelation = 
                     dHW3F_Calibration_Params.InputArray_FwdCorrelation_2_20.at(uiK);

                dImpliedCorrelation=ImpliedFwdCorrelation_VFDK_HW1To3F(dZZ,
                                           dXX,
                                           dYY,
                                           2.,
                                           dHW3F_Calibration_Params.ZCDates,
                                           dHW3F_Calibration_Params.ZCRates,
                                           dHW3F_Calibration_Params.ZCRates,
                                           dVolStripDates,
                                           dVolStrip,
                                           dModelParameters,
                                           dHW3F_Calibration_Params.SpotDate);               

                dFwdCorrelationFunctionValue+= 10000*dSurfaceWeight
                                               *(dImpliedCorrelation-dMarketCorrelation)
                                               *(dImpliedCorrelation-dMarketCorrelation);

                dSum+=dSurfaceWeight;
            } // if on weights non-zero
        }

        if (dSum!=0.) 
           dFwdCorrelationFunctionValue/=dSum;

        dFunctionValue=dFunctionValue+dFwdCorrelationFunctionValue;

    }

    catch(...)
    {
        dHW3F_Calibration_Params.Plantage++; 
        dFunctionValue=1.e10;
    }

    if(dHW3F_Calibration_Params.uiIterations>1 
           && finite(dFunctionValue) && dFunctionValue<dHW3F_Calibration_Params.BestResult)
    {
        dHW3F_Calibration_Params.BestResult=dFunctionValue;
        dHW3F_Calibration_Params.BestParams=dModelParameters;
    }

    return dFunctionValue;
}



void NAG_Function_e04UCC(int NbParams, 
                         double* x, 
                         double* bl, 
                         double* bu, 
                         double* delta, 
                         double* Minimum)
{
#ifdef WIN32

    int IFAIL;

    double objf;
    Nag_BoundType bound;
    Nag_E04_Opt options;
    static NagError fail, fail2;
    double g[200];
   
	// Useful for managing slow convergence cases

	Nag_Comm comm;

	// Init.

    e04xxc(&options);

    options.optim_tol = THE_3F_PRECISION;

    options.local_search = FALSE;
    options.obj_deriv = FALSE;
    options.con_deriv = FALSE;
    
	options.f_prec = 0.0001;

   //  options.verify_grad = Nag_NoCheck;

    char fOutName[200];

    ARM_GetTmpAbsFile("nag.txt", fOutName);

    strcpy(options.outfile, fOutName);

    options.max_iter = 200;

    bound = Nag_Bounds;

    // Print options

    options.print_level = Nag_Soln_Iter_Const;

    options.minor_print_level = Nag_Iter;

    fail.print = FALSE;

	// For managing special slow convergence cases!
    // Init. the needed variables

    // comm.it_maj_prt = TRUE;

    // options.print_fun = printfun;

    // Previous iteration Point

    memset(dHW3F_Calibration_Params.X_PREC, 0.0, sizeof(double)*MAX_3F_PARAM);

    dHW3F_Calibration_Params.iterCounter = 0; // "Global" iterations counter

    double* a = NULL;

    fail.code = -10000;

	try
    {
       e04ucc(NbParams, 0, 0, a, 0, bl, bu,
              objfun, confun, x, &objf, g, &options, 
		      &comm,
              &fail);
	}

	catch(char* a3FException)
    {
         char buf[255];

         strcpy(buf, "HW Calibration(NAg): ");

         strcat(buf, a3FException);

         throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL, buf);
    }

    catch(Exception& theExpt)
    {
        throw theExpt;
    }

    catch(...)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
                           "HW Calibration(NAg): Unrecognized failure");
    }


    if ( fail.code == NE_USER_STOP )
    {
       // a monitored user exit

       IFAIL = 0; // RET_OK
    }
    else if ( fail.code == NW_COND_MIN )
    {
       // A Solution has been found, but may be it's not the best one!!"

       IFAIL = 0; // RET_OK
    }
    else if ( fail.code == NE_NOERROR )
    {
       IFAIL = 0;  // RET_OK
    }
    else if ( fail.code == NW_KT_CONDITIONS )
	{
	   IFAIL = 0; // RET_OK
	}
	else
    {
       IFAIL = -1; // RET_KO
    }

    // Clean Up
 
    e04xzc(&options,"all", &fail2);

    *Minimum = objf;

#endif
}



void Fill(double dSpotDate,
          DKMaille<double> &dZCDates,
          DKMaille<double> &dZCRates,
          DKMaille<double> &dModelParameters,
          DKMaille<double> &dWhichModelParametersToCalibrate,
          DKMaille<double> &dInputArray_X_AbsVols,
          DKMaille<double> &dInputArray_Y_AbsVols,
          DKMaille2D<double> &dInputSurface_AbsVols,
          DKMaille<double> &dInputArray_X_FwdCorrelation_2_20,
          DKMaille<double> &dInputArray_FwdCorrelation_2_20,
          DKMaille2D<double> &dInputSurface_AbsVolsWeights,
          DKMaille<double> &dInputArray_FwdCorrelation_2_20Weights,
          DKMaille<double> &dNoticesDates,
          DKMaille<double> &dSwapStartDates,
          DKMaille<double> &dSwapEndDates)
{
    dHW3F_Calibration_Params.SpotDate=dSpotDate;
    dHW3F_Calibration_Params.ZCDates=dZCDates;
    dHW3F_Calibration_Params.ZCRates=dZCRates;
    dHW3F_Calibration_Params.ModelParameters=dModelParameters;
    dHW3F_Calibration_Params.WhichModelParametersToCalibrate=dWhichModelParametersToCalibrate;
    dHW3F_Calibration_Params.InputArray_X_AbsVols=dInputArray_X_AbsVols;
    dHW3F_Calibration_Params.InputArray_Y_AbsVols=dInputArray_Y_AbsVols;
    dHW3F_Calibration_Params.InputSurface_AbsVols=dInputSurface_AbsVols;
    dHW3F_Calibration_Params.InputArray_X_FwdCorrelation_2_20=dInputArray_X_FwdCorrelation_2_20;
    dHW3F_Calibration_Params.InputArray_FwdCorrelation_2_20=dInputArray_FwdCorrelation_2_20;
    dHW3F_Calibration_Params.InputSurface_AbsVolsWeights=dInputSurface_AbsVolsWeights;
    dHW3F_Calibration_Params.InputArray_FwdCorrelation_2_20Weights
                                  =dInputArray_FwdCorrelation_2_20Weights;
    dHW3F_Calibration_Params.NoticeDatesBis=dNoticesDates;
    dHW3F_Calibration_Params.SwapStartDatesBis=dSwapStartDates;
    dHW3F_Calibration_Params.SwapEndDatesBis=dSwapEndDates;
    dHW3F_Calibration_Params.PreviousModelParameters.resize(dModelParameters.entries()+1);
    dHW3F_Calibration_Params.BestResult=1.e15;
    dHW3F_Calibration_Params.BestParams.resize(dModelParameters.entries());
    dHW3F_Calibration_Params.uiIterations=0;
    dHW3F_Calibration_Params.Plantage=0;    
}




DKMaille<double> HW3F_CALIBRATION(  double dSpotDate,
                                    DKMaille<double> &dZCDates,
                                    DKMaille<double> &dZCRates,
                                    DKMaille<double> &dModelParameters,
                                    DKMaille<double> &dInputArray_X_AbsVols,
                                    DKMaille<double> &dInputArray_Y_AbsVols,
                                    DKMaille2D<double> &dInputSurface_AbsVols,
                                    DKMaille<double> &dInputArray_X_FwdCorrelation_2_20,
                                    DKMaille<double> &dInputArray_FwdCorrelation_2_20,
                                    DKMaille2D<double> &dInputSurface_AbsVolsWeights,
                                    DKMaille<double> &dInputArray_FwdCorrelation_2_20Weights,
                                    DKMaille<double> &dNoticesDates,
                                    DKMaille<double> &dSwapStartDates,
                                    DKMaille<double> &dSwapEndDates
                                    //double MaxTime
                                        )
{

    DKMaille<double> dWhichModelParametersToCalibrate;
    dWhichModelParametersToCalibrate.resize(dModelParameters.entries());
    
    bool recherche_continue=0;

    if(dModelParameters.at(0)>3)
    {
        recherche_continue=1;
        dModelParameters.at(0)=3;
    }

    double BestResult=1.e15;
    DKMaille<double> BestParams(dModelParameters.entries());
    
    double LowerBound[MAX_3F_PARAM];
    double UpperBound[MAX_3F_PARAM];
    double Delta[MAX_3F_PARAM];

    for(int i=0;i<dModelParameters.entries();i++)
    {
        dWhichModelParametersToCalibrate.at(i)=0.;
    }

    if(dModelParameters.at(0)==1)
    {
        dWhichModelParametersToCalibrate.at(1)=1.;
        LowerBound[0]=0.0001;
        UpperBound[0]=2.;
        Delta[0]=0.0001;
    }
    else if(dModelParameters.at(0)==2)
    {
        dWhichModelParametersToCalibrate.at(1)=1.;
        dWhichModelParametersToCalibrate.at(2)=1.;
        dWhichModelParametersToCalibrate.at(4)=1.;
        dWhichModelParametersToCalibrate.at(6)=1.;
        LowerBound[0]=0.0001;
        LowerBound[1]=0.0001;
        LowerBound[2]=0.0001;
        LowerBound[3]=-0.99;
        UpperBound[0]=1.e10;
        UpperBound[1]=1.e10;
        UpperBound[2]=1.e10;
        UpperBound[3]=.99;
        Delta[0]=0.0001;
        Delta[1]=0.0001;
        Delta[2]=0.0003;
        Delta[3]=0.0005;
    }
    else if(dModelParameters.at(0)==3)
    {
        dWhichModelParametersToCalibrate.at(1)=1.;
        dWhichModelParametersToCalibrate.at(2)=1.;
        dWhichModelParametersToCalibrate.at(3)=1.;
        dWhichModelParametersToCalibrate.at(4)=1.;
        dWhichModelParametersToCalibrate.at(5)=1.;
        dWhichModelParametersToCalibrate.at(6)=1.;
        dWhichModelParametersToCalibrate.at(7)=1.;
        dWhichModelParametersToCalibrate.at(8)=1.;
        LowerBound[0]=0.0001;
        LowerBound[1]=0.0001;
        LowerBound[2]=0.0001;
        LowerBound[3]=0.0001;
        LowerBound[4]=0.0001;
        LowerBound[5]=-0.99;
        LowerBound[6]=-0.99;
        LowerBound[7]=-0.99;
        UpperBound[0]=1.e10;
        UpperBound[1]=1.e10;
        UpperBound[2]=1.e10;
        UpperBound[3]=1.e10;
        UpperBound[4]=1.e10;
        UpperBound[5]=.99;
        UpperBound[6]=.99;
        UpperBound[7]=.99;
        Delta[0]=0.0001;
        Delta[1]=0.0001;
        Delta[2]=0.0001;
        Delta[3]=0.0003;
        Delta[4]=0.0003;
        Delta[5]=0.0005;
        Delta[6]=0.0005;
        Delta[7]=0.0005;
    }
    else
        throw("Wrong Model Parameters");

    Fill( dSpotDate,
          dZCDates,
          dZCRates,
          dModelParameters,
          dWhichModelParametersToCalibrate,
          dInputArray_X_AbsVols,
          dInputArray_Y_AbsVols,
          dInputSurface_AbsVols,
          dInputArray_X_FwdCorrelation_2_20,
          dInputArray_FwdCorrelation_2_20,
          dInputSurface_AbsVolsWeights,
          dInputArray_FwdCorrelation_2_20Weights,
          dNoticesDates,
          dSwapStartDates,
          dSwapEndDates);

    unsigned int uiHowManyToOptimize=0;

    for(unsigned int uM=0;uM<dWhichModelParametersToCalibrate.entries();uM++)
    {
        if(dWhichModelParametersToCalibrate.at(uM)==1.)
            uiHowManyToOptimize++;
    }

    // parameter0 = Nb Factors
    // parameter1 = mean rev 1
    // parameter2 = mean rev 2
    // parameter3 = mean rev 3
    // parameter4 = relative factor 1-2
    // parameter5 = relative factor 1-3
    // parameter6 = correlation 1-2
    // parameter7 = correlation 1-3
    // parameter8 = correlation 2-3



    double X0[MAX_3F_PARAM]; // In fact [uiHowManyToOptimize] need;

    // Fill the depart array
    int uk=0;
    for(int ui=0;ui<dModelParameters.entries();ui++)
    {
        if(dWhichModelParametersToCalibrate.at(ui)==1.)
        {

            X0[uk]=dModelParameters.at(ui);
            uk++;
        }
    }

    double FinalResult=0.;
    
    try
    {
        NAG_Function_e04UCC(uiHowManyToOptimize, X0, LowerBound, UpperBound,
                            Delta, &FinalResult);
    }

	catch(char* a3FException)
    {
         char buf[255];

         strcpy(buf, "HW3F_CALIBRATION(NAG_Function_e04UCC): ");

         strcat(buf, a3FException);

         throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL, buf);
    }

    catch(Exception& theExpt)
    {
        throw theExpt;
    }

    catch(...)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
                           "HW3F_CALIBRATION(NAG_Function_e04UCC): Unrecognized failure");
    }
    
    if(recherche_continue)
    {
        //on sauve le resultat precedent
        if(finite(FinalResult) && FinalResult<BestResult)
        {
            BestResult=FinalResult;
            uk=0;
            for(ui=0;ui<dModelParameters.entries();ui++)
            {
                if(dWhichModelParametersToCalibrate.at(ui)==1.)
                {
                    BestParams.at(ui)=X0[uk];
                    uk++;
                }
                else
                {
                    BestParams.at(ui)=dModelParameters.at(ui);
                }
            }
        }
        
        double StartValue;
        FILE* fic;

        char fOutName[200];

        ARM_GetTmpAbsFile("nagResult.txt", fOutName);

        fic = fopen(fOutName, "w");
        
        //  Initiate the Random generator
        srand( (unsigned)time( NULL ) );

        while(TRUE)
        {
            // On tire les points de  depart au hasard
            dModelParameters.at(1)=1.*rand()/RAND_MAX;
            dModelParameters.at(2)=1.*rand()/RAND_MAX;
            dModelParameters.at(3)=1.*rand()/RAND_MAX;
            dModelParameters.at(4)=3.*rand()/RAND_MAX;
            dModelParameters.at(5)=3.*rand()/RAND_MAX;
            dModelParameters.at(6)=-0.99+1.98*rand()/RAND_MAX;
            dModelParameters.at(7)=-0.99+1.98*rand()/RAND_MAX;
            dModelParameters.at(8)=-0.99+1.98*rand()/RAND_MAX;


            Fill(dSpotDate,
                 dZCDates,
                 dZCRates,
                 dModelParameters,
                 dWhichModelParametersToCalibrate,
                 dInputArray_X_AbsVols,
                 dInputArray_Y_AbsVols,
                 dInputSurface_AbsVols,
                 dInputArray_X_FwdCorrelation_2_20,
                 dInputArray_FwdCorrelation_2_20,
                 dInputSurface_AbsVolsWeights,
                 dInputArray_FwdCorrelation_2_20Weights,
                 dNoticesDates,
                 dSwapStartDates,
                 dSwapEndDates);

            // Fill the depart array
            uk=0;
            for(ui=0;ui<dModelParameters.entries();ui++)
            {
                if(dWhichModelParametersToCalibrate.at(ui)==1.)
                {

                    X0[uk]=dModelParameters.at(ui);
                    uk++;
                }
            }

            StartValue=HW3F_Calibration_Function_TO_Minimize(X0);

            if ( finite(StartValue) && StartValue < 1.e7 )
            {
                try
                {
                    NAG_Function_e04UCC(uiHowManyToOptimize, X0, LowerBound, 
                                        UpperBound, Delta, &FinalResult);
                }

				catch(char* a3FException)
				{
                    char buf[255];

                    strcpy(buf, "HW3F_CALIBRATION(NAG_Function_e04UCC): ");

                    strcat(buf, a3FException);

                    throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL, buf);
				}

                catch(Exception& theExpt)
				{
                    throw theExpt;
				}

                catch(...)
				{
                     throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
                           "HW3F_CALIBRATION(NAG_Function_e04UCC): Unrecognized failure");
				}

                if(finite(FinalResult) && FinalResult<BestResult)
                {
                    BestResult=FinalResult;
                    uk=0;
                    for(ui=0;ui<dModelParameters.entries();ui++)
                    {
                        if(dWhichModelParametersToCalibrate.at(ui)==1.)
                        {
                            BestParams.at(ui)=X0[uk];
                            uk++;
                        }
                        else
                        {
                            BestParams.at(ui)=dModelParameters.at(ui);
                        }
                    }
                }
                
                fprintf(fic, "\n");
                fprintf(fic, "////////////////////////////////////////////////\n");
                fprintf(fic, "Meilleur Resultat\n");
                fprintf(fic, "Ecart = ");
                fprintf(fic, "%lf",BestResult);
                fprintf(fic, "\n");
                fprintf(fic, "Param[0] = ");
                fprintf(fic, "%lf",BestParams.at(0));
                fprintf(fic, "\n");
                fprintf(fic, "Param[1] = ");
                fprintf(fic, "%lf",BestParams.at(1));
                fprintf(fic, "\n");
                fprintf(fic, "Param[2] = ");
                fprintf(fic, "%lf",BestParams.at(2));
                fprintf(fic, "\n");
                fprintf(fic, "Param[3] = ");
                fprintf(fic, "%lf",BestParams.at(3));
                fprintf(fic, "\n");
                fprintf(fic, "Param[4] = ");
                fprintf(fic, "%lf",BestParams.at(4));
                fprintf(fic, "\n");
                fprintf(fic, "Param[5] = ");
                fprintf(fic, "%lf",BestParams.at(5));
                fprintf(fic, "\n");
                fprintf(fic, "Param[6] = ");
                fprintf(fic, "%lf",BestParams.at(6));
                fprintf(fic, "\n");
                fprintf(fic, "Param[7] = ");
                fprintf(fic, "%lf",BestParams.at(7));
                fprintf(fic, "\n");
                fprintf(fic, "Param[8] = ");
                fprintf(fic, "%lf",BestParams.at(8));
                fprintf(fic, "\n");
                fprintf(fic, "////////////////////////////////////////////////\n");
            }
        }

        // On remet le meilleur resultat dans X0 et FinalResult
        FinalResult=BestResult;
        uk=0;
        for(ui=0;ui<dModelParameters.entries();ui++)
        {
            if(dWhichModelParametersToCalibrate.at(ui)==1.)
            {
                X0[uk]=BestParams.at(ui);
                uk++;
            }
        }

        fclose(fic);
    }
    

    DKMaille<double> dModelVolsResult(dModelParameters.entries()+2);

    // Fill the result array
    uk=0;
    for(ui=0;ui<dModelParameters.entries();ui++)
    {
        if (dWhichModelParametersToCalibrate.at(ui)==1.)
        {
           dModelVolsResult.at(ui)=X0[uk];
           uk++;
        }
        else
        {
           dModelVolsResult.at(ui)=dModelParameters.at(ui);
        }
    }

    dModelVolsResult.at(dModelParameters.entries())=(double) dHW3F_Calibration_Params.Plantage;
    dModelVolsResult.at(dModelParameters.entries()+1) = FinalResult;
    
    return dModelVolsResult;
}



/*----------------------------------------------------------------------------*/
/*---- End Of File ----*/
