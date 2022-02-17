/*
 $Log: calibrationhwsv.cpp,v $
 Revision 1.28  2003/09/08 11:58:36  mab
 the input portfolio is deeply cloned now in constructors!

 Revision 1.27  2003/09/02 16:47:35  sgasquet
  Ajout constructeur prenant le modele de marche en input

 Revision 1.26  2003/07/04 07:24:11  ebenhamou
 move the destructor to the cpp file since we only know the declaration of the ARM_Portfolio and not its definition

 Revision 1.25  2003/06/30 16:42:41  ebenhamou
 remove unused var

 Revision 1.24  2003/02/28 18:13:50  sgasquet
 Ajout controle du vega pour les flexaccretswaption

 Revision 1.23  2003/01/20 19:08:31  sgasquet
 Ajout modif pour ARM_FLEXACCRETSWAPTION

 Revision 1.22  2002/12/10 10:38:29  sgasquet
 Modif vol par defaut dans constructeurs; Ajout fonctions pour ARM_OPTIONALACCRUALZC

 Revision 1.21  2002/11/21 12:39:19  mab
 Init added to some constructors

 Revision 1.20  2002/10/09 09:36:38  mab
 Formatting

 Revision 1.19  2002/02/21 11:09:31  vberger
 correction de la non prise en compte de itsvolmin et itsvolmax
 ds la fct void ARM_CalibrationHWSV::Calibrate(char * msg, ARM_Calibration_PORTTYPE porttype, ARM_CalibrationHWSV_OUTPUT& Output)

 Revision 1.18  2002/01/25 15:30:47  sgasquet
 Augmentation precision calage

 Revision 1.17  2002/01/22 15:28:42  sgasquet
 Externationalisation Calcul Vol Implicite

 Revision 1.16  2002/01/18 16:24:32  sgasquet
 Initialisation indicateur d'erreur

 Revision 1.15  2002/01/16 14:47:53  sgasquet
 Ajout infos pour calage stickys

 Revision 1.12  2001/11/22 14:06:31  nicolasm
 Initialisation volmin et volmax

 Revision 1.11  2001/11/22 14:01:20  nicolasm
 Initialisation des mb itsvolmin et itsvolmax

 Revision 1.10  2001/11/20 10:56:59  vberger
 allˆgement de la fct calaibrate

 Revision 1.9  2001/10/24 09:06:40  vberger
 modifs dues aux modifs de l'enum portype
 modifs dues aux modifs de la structure d'ouput

 Revision 1.8  2001/06/28 15:21:48  vberger
 ajout gTrace pour ‰viter printf ou cout

 Revision 1.7  2001/06/25 11:53:24  vberger
 modif de calib pour la calibration des MCPRCS

 Revision 1.6  2001/05/07 15:57:17  vberger
 ajout de void ARM_CalibrationHWSV::Calibrate(char * msg,ARM_Calibration_PORTTYPE porttype,ARM_CalibrationHWSV_OUTPUT& Output)

 Revision 1.5  2001/01/23 09:12:11  nicolasm
 Ds Calibrate suppression des instruments

 Revision 1.4  2000/11/14 08:03:30  nicolasm
 Division de l'indicateur prix  par PFSize

 Revision 1.3  2000/10/03 17:17:55  nicolasm
 Ajout du try dans fonction Calibrate
 Prise en compte des pds nuls
 Passage de 5 a 10 jour sur les plots cv sigma

 Revision 1.2  2000/08/03 12:23:15  nicolasm
  Ajout indicateurs

 Revision 1.1  1999/11/10 08:55:12  nicolasm
 Initial revision

 */


/*-----------------------------------------------------------------*/


#include "calibrationhwsv.h"
#include "security.h"
#include "bootstrapcalibration.h"

#include "bssmiled.h"

#include "swaption.h"

#include "reverse.h"

#include "capfloor.h"


/*-----------------------------------------------------------------*/



ARM_ConfigHWSVCalibration* ARM_CalibrationHWSV::itsConfigCalib = NULL;



/// because the destructor calls the one of ARM_Portfolio
/// it has to know the exact definition of the class
/// otherwise it will issue a warning

ARM_CalibrationHWSV::~ARM_CalibrationHWSV(void)
{
    if (itsPf)
    {
       itsPf->FreePortfolioAndAssets();

       delete itsPf; 
       itsPf = NULL;
    }

    if (itsSigmaDates)
       delete itsSigmaDates;

    if (itsMktModelToFit)
    {
       delete itsMktModelToFit;

       itsMktModelToFit = NULL;
    }
}


ARM_CalibrationHWSV::ARM_CalibrationHWSV(ARM_Date& asOf, ARM_ZeroCurve* zc, 
                                         ARM_VolCurve* vol, ARM_Security* sec,
                                         double amin ,double amax , 
                                         ARM_Vector* dates, ARM_Portfolio* pf)
                                        : ARM_Calibration(asOf, zc)
{
    Init();

    itsVolCurve = vol;

    itsSecurity = sec;

    itsamin = amin;
    itsamax = amax;
    itsvolmin = 0.001;
    itsvolmax = 5.0;

    itsIndicPrix = 0.0;
    itsIndicVol = 0.0; 

    if (dates)
       itsSigmaDates = (ARM_Vector *) dates->Clone();

    if (pf)
       itsPf = (ARM_Portfolio *) pf->DeepClone();
}



ARM_CalibrationHWSV::ARM_CalibrationHWSV(ARM_Date& asOf, ARM_ZeroCurve* zc, 
                                         ARM_VolCurve* vol, ARM_Security* sec,
                                         double amin , double amax,
                                         double volmin , double volmax,
                                         ARM_Vector* dates , ARM_Portfolio* pf,
                                         ARM_VolCurve* RhoSwopt,
                                         ARM_VolCurve* NuSwopt,
                                         ARM_VolCurve* BetaSwopt,
										 int isAlphaOrSigma) 
                                         : ARM_Calibration(asOf, zc)
{
    Init();

    itsVolCurve = vol;

    itsSecurity = sec;

    itsamin = amin;

    itsamax = amax;

    itsvolmin = volmin;

    itsvolmax = volmax;

    itsIndicPrix = 0.0;

    itsIndicVol = 0.0; 


    if (dates)
       itsSigmaDates = (ARM_Vector *) dates->Clone();

    if (pf)
       itsPf = (ARM_Portfolio *) pf->DeepClone();

    
    // SABR Model

    ARM_BSModel* swoptModel = NULL;

	if (RhoSwopt && NuSwopt)
    { 
       if (BetaSwopt)
       {
          int SABR_Flag;

          if (BetaSwopt->IsEqualToOne())
          {
             SABR_Flag = K_SABR_ARITH;
          }
          else
          {
             SABR_Flag = K_SABR_IMPLNVOL;
          }

          swoptModel = new ARM_BSSmiledModel(asOf, 
                                             0,
								             zc,
								             zc,
								             vol,
								             K_YIELD,
								             RhoSwopt,
								             NuSwopt,
       						                 SABR_Flag,
                                             BetaSwopt,
											 0.5,				// Default value: not used in this case
											 isAlphaOrSigma);
       }
       else
       {
	      swoptModel = new ARM_BSSmiledModel(asOf, 
                                             0,
									         zc,
									         zc,
									         vol,
									         K_YIELD,
									         RhoSwopt,
									         NuSwopt,
       								         K_SABR_ARITH,
											 NULL,				// Default value: not used in this case
											 0.5,				// Default value: not used in this case
											 isAlphaOrSigma);
       }
    }
    else
    {
		swoptModel = new ARM_BSModel(asOf,
									 0, 
									 zc,
									 zc, 
									 vol,
									 K_PRICE);
	}

    itsMktModelToFit = swoptModel;
}




ARM_CalibrationHWSV::ARM_CalibrationHWSV(ARM_Date& asOf, 
                                         ARM_Model* model, ARM_Security* sec,
                                         double amin, double amax, 
                                         double volmin, double volmax, 
                                         ARM_Vector* dates, ARM_Portfolio* pf)
                                         : ARM_Calibration(asOf, model->GetZeroCurve())
{
    Init();

    itsMktModelToFit = (ARM_Model *) model->Clone();

    itsSecurity = sec;

    itsamin = amin;

    itsamax = amax;

    itsvolmin = volmin;

    itsvolmax = volmax;

    itsIndicPrix = 0.0;

    itsIndicVol = 0.0; 

    if (dates)
       itsSigmaDates = (ARM_Vector *) dates->Clone();

    if (pf)
       itsPf = (ARM_Portfolio *) pf->DeepClone();
}



void ARM_CalibrationHWSV::Calibrate(char* msg,
                                    ARM_Calibration_PORTTYPE porttype,
                                    ARM_CalibrationHWSV_OUTPUT& Output)

{
    ARM_GYCSigVarModel* mod = NULL;


    // fabrication du portfeuille par l'instrument

    // -------------------------------------------

    int diagsize;

    try 
    {
        if (!itsPf)
        {
           itsPf = ((ARM_Reverse *) itsSecurity)->CalibrationPF(GetAsOfDate(),
                                                                porttype, diagsize);
        }

        if (!itsPf)
        {
           ARM_Date asof = GetAsOfDate();

           Output.meanrev = 0.03;

           Output.HWVolSize = 10;

           Output.PFSize = 10;

           Output.HWVectorStability = NOJUMP ;

           Output.calagequality = 0.0;

           Output.Asof = GetAsOfDate(); 

           for (int i = 0; i < 10; i++)
           {
               Output.datesched[i] = asof.AddYears(1+i);

               Output.volsched[i] = 0.6 ;

               Output.SecurityExeDate[i] = asof.AddYears(1+i);

               Output.SecurityStartDate[i] = asof.AddYears(1+i);

               Output.SecurityEndDate[i] = asof.AddYears(5+i);

               Output.SecurityStrike[i] = -1000000;

               Output.InputPrice[i]= 0.1 ;

               Output.InputVol[i]= 0.1;

               Output.OutputPrice[i]= 0.1;

               Output.ErrorPrice[i]= 0.1; 

               Output.OutputVol[i]= 0.1;

               Output.ErrorVol[i]= 0.1 ;
           }

           return;
        }

        // fabrication du model

        ARM_Date testdate = GetAsOfDate() ;

        ARM_BSModel* bsmod = NULL;

        if ( itsMktModelToFit == NULL )
        {
           bsmod = new ARM_BSModel(GetAsOfDate(), 0.0, 
                                   GetZeroCurve(),
                                   GetZeroCurve(), 
                                   itsVolCurve,
                                   K_YIELD);
   
           itsMktModelToFit = bsmod;
        }
        else
        {
           bsmod = (ARM_BSModel *) itsMktModelToFit;
        }

        itsPf->SetModel(bsmod);

        double price;

        int PFSize = itsPf->GetSize();

        ARM_Vector* mktprice = itsPf->GetMktPrices();

        ARM_Vector* weights  = itsPf->GetWeights();

        int i;

        for (i = 0; i < PFSize; i++)
        {
            price = itsPf->GetAsset(i)->ComputePrice();

            mktprice->Elt(i) = price;

            if ( price < 1e-08 )
               weights->Elt(i) = 0.0;
            else
               weights->Elt(i) = 1.0/(price*price);
        }

        // recuperation de la liste de dates

        // ---------------------------------

        double X[200];

        double Y[200];

        double vals_curve[200];

        MEMSET(Y, 0.0, sizeof(double)*200);

        MEMSET(vals_curve, 0.0, sizeof(double)*200);


        Output.Asof = GetAsOfDate(); 

        if (gTrace)
        {
           cout<<"PF asset exercising date   "<<endl;
        }

        for (i = 0; i < PFSize; i++)
        {
            Output.SecurityExeDate[i] = itsPf->GetAsset(i)->GetExpiryDate() ;

            Output.SecurityStartDate[i] = 
                 ((ARM_Swap*)(itsPf->GetAsset(i)))->GetStartDate() ;

            Output.SecurityEndDate[i] = 
             ((ARM_Swap*)(itsPf->GetAsset(i)))->GetEndDate();

            Output.SecurityStrike[i] = 
                      ((ARM_Swaption*)(itsPf->GetAsset(i)))->GetStrike();
        }

        for (i = 0; i < diagsize; i++)
        {
            Output.datesched[i] = itsPf->GetAsset(i)->GetExpiryDate()  ;

            Output.datesched[i].AddDays(1);

            X[i] = (itsPf->GetAsset(i)->GetExpiryDate()
                    -GetAsOfDate())/365.0 + 1.0/365.0;

        }

        // Appel de la fonction de minimisation

        // ------------------------------------

        ARM_GYCSigVarModel* sigVarModel = new ARM_GYCSigVarModel(itsamin, 
                                                                 diagsize, 
                                                                 X,
                                                                 Y,
                                                                 GetZeroCurve());
        for (i = 0; i < PFSize; i++)
        {
            itsPf->GetAsset(i)->SetModel(sigVarModel);
        }

        double precision_meanRev = 0.01,

        precision_vol = 0.001;

        long nbMaxIter = 100; //1000

        double optimal_meanrev = ComputeMeanRevAndParasCurve(sigVarModel,
                                                             itsPf,
                                                             diagsize,
                                                             X,
                                                             precision_meanRev,
                                                             precision_vol,
                                                             itsamin,
                                                             itsvolmin,
                                                             itsamax,
                                                             itsvolmax,
                                                             nbMaxIter,
                                                             vals_curve);


        Output.meanrev = optimal_meanrev;

        // HW volatility array

        for (i = 0; i < diagsize; i++)
        {
            Output.volsched[i] = vals_curve[i];   
        }

        delete sigVarModel;

        mod = new ARM_GYCSigVarModel(optimal_meanrev,
                                     diagsize,
                                     X, 
                                     vals_curve,
                                     GetZeroCurve());

        // calcul des indicateurs

        double ecart, prixHW, prixBS;

        ecart = 0.0;

        Output.HWVolSize = diagsize;

        Output.PFSize = PFSize;

        Output.HWVectorStability = NOJUMP ;

        Output.calagequality = 0.0;

        for (i = 0; i < PFSize; i++)
        {
            prixBS = mktprice->Elt(i);

            itsPf->GetAsset(i)->SetModel(mod);

            prixHW = itsPf->GetAsset(i)->ComputePrice();

            ecart = prixHW / prixBS - 1;

            ecart *= ecart;

            itsIndicPrix += ecart ;

            Output.InputPrice[i] = prixBS ;

            Output.InputVol[i] = itsPf->GetAsset(i)->ComputeImpliedVol(prixBS);

            Output.OutputPrice[i] = prixHW;

            Output.ErrorPrice[i] = prixHW-prixBS; 

            Output.OutputVol[i] = itsPf->GetAsset(i)->ComputeImpliedVol(prixHW);

            Output.ErrorVol[i] = Output.OutputVol[i]-Output.InputVol[i] ;

            Output.calagequality += Output.ErrorVol[i]*Output.ErrorVol[i] ;
        }

        Output.calagequality = sqrt(Output.calagequality/PFSize)*100;

        double minvol =  1000000;

        double maxvol = -1000000;

        double voli ;

        for (i = 0; i < diagsize-1; i++)
        {
            voli = Output.InputPrice[i];

            if ( voli > maxvol )
            {
               maxvol = voli;
            }

            if ( voli<minvol )
            {
               minvol = voli;
            }

            if ( fabs (Output.InputPrice[i+1]-voli)
                 > MIN(voli,Output.InputPrice[i+1]) 
               )
               Output.HWVectorStability = PREVIOUSTONEXTJUMP ;
        }

        voli = Output.InputPrice[diagsize-1];

        if ( voli > maxvol )
        {
           maxvol = voli;
        }

        if ( voli < minvol )
        {
           minvol = voli;
        }

        if ( maxvol-minvol > minvol )
        {
           if ( Output.HWVectorStability == PREVIOUSTONEXTJUMP )
              Output.HWVectorStability = BOTHJUMP ;
           else
              Output.HWVectorStability = MAXMINJUMP ;

        }

        itsIndicPrix = sqrt(itsIndicPrix);

        itsIndicPrix /= PFSize;

        // if (bsmod)
        // delete bsmod ;

        // Portfolio DEEP! deletion 
        for (i = 0; i < PFSize; i++)
        {
            if (itsPf->GetAsset(i))
               delete itsPf->GetAsset(i);
        }

        delete itsPf;

        itsPf = NULL;
    } // end try

    catch(Exception& m)
    {
        // printf("Exception caught in Calibrate in calibrationhwsv.cpp \n");

		std::string e=m.GetErrorString(); 
		if (msg) strncpy(msg,e.c_str(),ARM_EXCEPTION_MSG_MAX_SIZE); 
        
    }
}



ARM_GYCSigVarModel* ARM_CalibrationHWSV::Calibrate(char* msg)
{
    ARM_GYCSigVarModel* mod   = NULL;
    ARM_Model*          bsmod = NULL;



    // fabrication du portfeuille par l'instrument
    // -------------------------------------------

    try 
    {
        if (!itsPf)
           itsPf = itsSecurity->CalibrationPF(GetAsOfDate());

        // fabrication du model si pas donne en input
        if (itsMktModelToFit)
        {
           bsmod = (ARM_Model*) itsMktModelToFit->Clone();
        }
        else
        {
           bsmod = new ARM_BSModel(GetAsOfDate(), 0.0, 
                                   GetZeroCurve(),
                                   GetZeroCurve(), 
                                   itsVolCurve,
                                   K_YIELD);
        }
        
        itsPf->SetModel(bsmod);

        double price;

        int PFSize = itsPf->GetSize();

        ARM_Vector* mktprice = itsPf->GetMktPrices();
        ARM_Vector* weights  = itsPf->GetWeights();

        int i;

        for (i = 0; i < PFSize; i++)
        {
            price = itsPf->GetAsset(i)->ComputePrice();
         
            mktprice->Elt(i) = price;

            if ( price < 1e-08 )
               weights->Elt(i) = 0.0;
            else
               weights->Elt(i) = 1/(price*price);
        }
    
        // recuperation de la liste de dates
        // ---------------------------------

        // calcul de la maturity de l'instrument

        ARM_Date matDate = itsSecurity->GetExpiryDate();

        matDate.AddDays(5);

        double maturity = (matDate - GetAsOfDate())/K_YEAR_LEN;

        ARM_GEN_MATU secgm;

        if ( maturity < 5.0 )
           secgm = GM_5Y;
        else if (maturity < 10.0)
           secgm = GM_10Y;
        else if (maturity < 15.0)
           secgm = GM_15Y;
        else if (maturity < 20.0)
           secgm = GM_20Y;
        else if (maturity < 30.0)
           secgm = GM_30Y;
        else
           throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                           "security maturity too great !");

       // Calcul du Product_type

       ARM_PRODUCT_TYPE pt = itsSecurity->CalibrationPT();

       // Calcul de la maturite du sous-jacent

       ARM_GEN_MATU undergm = itsSecurity->CalibrationMatuIndex();

       ARM_SigmaDatesElt *sigDatesElt = 
               GetConfigCalib()->GetSigmaDates(pt, secgm, undergm);

       ARM_Vector* dates = sigDatesElt->GetSigmaDates();

       int sz = dates->GetSize();
    
       double X[200];
       double Y[200];
       double vals_curve[200];

       MEMSET(Y, 0.0, sizeof(double)*200);
       MEMSET(vals_curve, 0.0, sizeof(double)*200);

       for (i = 0; i < sz; i++)
       {
           X[i] = dates->Elt(i) + 10.0/365.0;
       }

       // Appel de la fonction de minimisation
       // ------------------------------------

       ARM_GYCSigVarModel* sigVarModel = new ARM_GYCSigVarModel(itsamin, 
                                                                sz, 
                                                                X,
                                                                Y,
                                                                GetZeroCurve());

       for (i = 0; i < PFSize; i++)
       {
           itsPf->GetAsset(i)->SetModel(sigVarModel);
       }

       double precision_meanRev = 0.01,
              precision_vol = 0.001;
           
       long nbMaxIter = 100; //1000

       double res = ComputeMeanRevAndParasCurve(sigVarModel,
                                                itsPf,
                                                sz,
                                                X,
                                                precision_meanRev,
                                                precision_vol,
                                                itsamin,
                                                itsvolmin,
                                                itsamax,
                                                itsvolmax,
                                                nbMaxIter,
                                                vals_curve);

       double optimal_meanrev = res;

       delete sigVarModel;

       mod = new ARM_GYCSigVarModel(optimal_meanrev,
                                    sz,
                                    X, 
                                    vals_curve,
                                    GetZeroCurve());

       // calcul des indicateurs
       double ecart, prixHW, prixBS;
       ecart = 0.0;

       for (i = 0; i < PFSize; i++)
       {
           if ( weights->Elt(i) != 0.0 )
           {
              prixBS = mktprice->Elt(i);
              itsPf->GetAsset(i)->SetModel(mod);
              prixHW = itsPf->GetAsset(i)->ComputePrice();
              ecart = prixHW / prixBS - 1;
              ecart *= ecart;
              itsIndicPrix += ecart;
           }
       }

       itsIndicPrix = sqrt(itsIndicPrix);

       itsIndicPrix /= PFSize;

       if (bsmod)
          delete bsmod;

       for (i = 0; i < PFSize; i++)
       {
           if (itsPf->GetAsset(i))
              delete itsPf->GetAsset(i);
       }

       delete itsPf;

       itsPf = NULL;

    } // end try

    catch(Exception& m)
    {
        // printf("Exception caught!! \n");
		std::string e=m.GetErrorString(); 
		if (msg) strncpy(msg,e.c_str(),ARM_EXCEPTION_MSG_MAX_SIZE); 
        return(NULL);
    }

    return(mod);
}



ARM_GYCSigVarModel* ARM_CalibrationHWSV::Calibrate(char* msg,
                                                   ARM_Calibration_TYPE porttype,
                                                   ARM_CalibrationHWSV_OUTPUT& Output)

{
    ARM_GYCSigVarModel* mod = NULL;
    ARM_Model*        bsmod = NULL;
    int IsPfComputed = 0;


    try 
    {
        // fabrication du portfeuille par l'instrument

        if (!itsPf)
        {
           itsPf = itsSecurity->CalibrationPF(GetAsOfDate(), porttype);
           IsPfComputed = 1;
        }

        // fabrication du model
        if (itsMktModelToFit)
        {
           bsmod = (ARM_Model *) itsMktModelToFit->Clone();
        }
        else
        {
           bsmod = new ARM_BSModel(GetAsOfDate(), 0.0, 
                                   GetZeroCurve(),
                                   GetZeroCurve(), 
                                   itsVolCurve,
                                   K_YIELD);
        }

        itsPf->SetModel(bsmod);

        double price;
        
        int PFSize = itsPf->GetSize();
        int nbPeriods;
        
        Output.PFSize = PFSize;

        Output.Asof = GetAsOfDate(); 

        ARM_Vector* mktprice = itsPf->GetMktPrices();

        ARM_Vector* weights  = itsPf->GetWeights();

        int i;

        for (i = 0; i < PFSize; i++)
        {
            price = itsPf->GetAsset(i)->ComputePrice();

            mktprice->Elt(i) = price;

            if ( price < 1e-08 )
               weights->Elt(i) = 0.0;
            else
               weights->Elt(i) = 1.0/(price*price);

        }

        double X[200];
        double Y[200];
        double vals_curve[200];

        MEMSET(Y, 0.0, sizeof(double)*200);

        MEMSET(vals_curve, 0.0, sizeof(double)*200);


        if (( itsSecurity->GetName() == ARM_OPTIONALACCRUALZC ) 
            || 
            ( itsSecurity->GetName() == ARM_FLEXACCRETSWAPTION ))
        {
            int j =0;
            int n =0;
            nbPeriods = (long) ((double)PFSize) * 0.5;

            // we check if the first security is not to ITM / OTM
            double vega = ((ARM_Swaption*)(itsPf->GetAsset(0)))->GetBSVega();

            if ( fabs(vega) < 0.0005 )
            {
               weights->Elt(0) = 0.0;
               weights->Elt(1) = 0.0;
            }

            for (i = 0; i < nbPeriods; i++)
            {
                if ( weights->Elt(j) > 0.0 )
                {
                   X[n] = (itsPf->GetAsset(j)->GetExpiryDate().AddDays(2)
                          -GetAsOfDate())/365.0;

                   n = n + 1;                
                }

                j = j+2;
            }

            nbPeriods = n;

            if ( nbPeriods == 0 ) // on garde au moins la maturite la plus eloignee
            {
               X[0] = (itsPf->GetAsset(PFSize-1)->GetExpiryDate().AddDays(2)
                      -GetAsOfDate())/365.0;

               nbPeriods = 1;
            }
        }
        else
        {
            nbPeriods = PFSize;
            for (i = 0; i < PFSize; i++)
            {
                X[i] = 
                 (itsPf->GetAsset(i)->GetExpiryDate().AddDays(2)
                                        -GetAsOfDate())/365.0;
            }
        }


        // Appel de la fonction de minimisation
        ARM_GYCSigVarModel* sigVarModel = new ARM_GYCSigVarModel(itsamin, 
                                                                 nbPeriods, 
                                                                 X,
                                                                 Y,
                                                                 GetZeroCurve());

        for (i = 0; i < PFSize; i++)
        {
            itsPf->GetAsset(i)->SetModel(sigVarModel);
        }

        double precision_meanRev = 0.0001, precision_vol = 0.00001; //, 
               //min_vol = 0.001,  max_vol = 5.0;    

        long nbMaxIter = 1000;

        double optimal_meanrev = ComputeMeanRevAndParasCurve(sigVarModel,
                                                             itsPf,
                                                             nbPeriods,
                                                             X,
                                                             precision_meanRev,
                                                             precision_vol,
                                                             itsamin,
                                                             itsvolmin,
                                                             itsamax,
                                                             itsvolmax,
                                                             nbMaxIter,
                                                             vals_curve);
        delete sigVarModel;

        mod = new ARM_GYCSigVarModel(optimal_meanrev, nbPeriods, X, 
                                     vals_curve, GetZeroCurve());

        double prixHW, prixBS;

        Output.calagequality  = 0.0;


        for (i = 0; i < PFSize; i++)
        {
            if ( itsPf->GetAsset(i)->GetName() == ARM_CAPFLOOR )
            {
                Output.SecurityStartDate[i] = 
                 ((ARM_CapFloor*)(itsPf->GetAsset(i)))->GetSwapLeg()->GetStartDate() ;

                Output.SecurityEndDate[i] = 
                    ((ARM_CapFloor*)(itsPf->GetAsset(i)))->GetSwapLeg()->GetEndDate();

                Output.SecurityStrike[i] = 
                      ((ARM_CapFloor*)(itsPf->GetAsset(i)))->GetStrike();
            }
            else
            {
                if ( itsPf->GetAsset(i)->GetName() == ARM_SWAPTION )            
                {
                    Output.SecurityStartDate[i] = 
                     ((ARM_Swaption*)(itsPf->GetAsset(i)))->GetStartDate() ;

                    Output.SecurityEndDate[i] = 
                        ((ARM_Swaption*)(itsPf->GetAsset(i)))->GetEndDate();

                    Output.SecurityStrike[i] = 
                          ((ARM_Swaption*)(itsPf->GetAsset(i)))->GetStrike();
                }
                else
                {
                    throw Exception(__LINE__, __FILE__, ERR_OBJECT_NULL,
                    "Unknown security in calibration portfolio");
                }
            }

            Output.meanrev = optimal_meanrev;
            Output.datesched[i] = itsPf->GetAsset(i)->GetExpiryDate().AddDays(2);
            Output.volsched[i] = vals_curve[i];  
            prixBS = mktprice->Elt(i);
            itsPf->GetAsset(i)->SetModel(mod);
            prixHW = itsPf->GetAsset(i)->ComputePrice();
            Output.InputPrice[i]= prixBS;
            Output.InputVol[i]= -100.0;
            Output.OutputPrice[i]= prixHW;
            Output.ErrorPrice[i]= prixHW - prixBS; 
            Output.OutputVol[i]= -100.0;
            Output.ErrorVol[i] = 0.0;
        }

        for (i = 0; i < PFSize; i++)
        {
            Output.InputVol[i] = 
                itsPf->GetAsset(i)->ComputeImpliedVol(Output.InputPrice[i]);

            Output.OutputVol[i] = 
                itsPf->GetAsset(i)->ComputeImpliedVol(Output.OutputPrice[i]);

            Output.ErrorVol[i]= Output.OutputVol[i] - Output.InputVol[i] ;
            Output.calagequality += Output.ErrorVol[i] * Output.ErrorVol[i] ;
        }

        if (bsmod)
           delete bsmod ;

        for (i = 0; i < PFSize; i++)
        {
            if (itsPf->GetAsset(i))
               delete itsPf->GetAsset(i);
        }

        delete itsPf;

        itsPf = NULL;
    } // end try

    catch(Exception& m)
    {
        if (IsPfComputed)
        {
           delete itsPf;
           itsPf = NULL;
        }

        // printf("Exception caught!!!\n");
		std::string e=m.GetErrorString(); 
		if (msg) strncpy(msg,e.c_str(),ARM_EXCEPTION_MSG_MAX_SIZE); 


        return(NULL);
    }

    if (IsPfComputed)
    {
       delete itsPf;
       itsPf = NULL;
    }

    return(mod);
}



/*------------------------------------------------------------------*/
/*----- End Of File ----*/
