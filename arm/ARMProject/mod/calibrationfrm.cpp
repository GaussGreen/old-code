/*
 $Log: calibrationfrm.cpp,v $
 Revision 1.7  2003/09/08 12:37:03  mab
 Added : composite objects Cloning()+Formatting

 Revision 1.6  2003/09/02 16:07:07  sgasquet
 Ajout constructeur prenant le modele de marche en input

 Revision 1.5  2002/01/25 15:31:30  sgasquet
 Amelioration Vol implicites

 Revision 1.4  2002/01/18 16:25:22  sgasquet
 Ajout parametre de filtration du portefeuille

 Revision 1.3  2002/01/17 14:23:30  sgasquet
  Ajout parametres de calage

 Revision 1.2  2002/01/16 14:49:07  sgasquet
  Ajout infos pour calage stickys

*/


/*------------------------------------------------------------------*/
#include "calibrationfrm.h"
#include "security.h"
#include "bootstrapcalibration.h"



ARM_CalibrationFRM::ARM_CalibrationFRM(ARM_Date& asOf, ARM_ZeroCurve *zc,
                                       ARM_VolCurve *VolCurve,
                                       ARM_Security *sec,
                                       int mfine , int shape_type  , 
                                       double decay, double slope, 
                                       double asymptote, int NbFactor , 
                                       double precisionVol, double minVols,
                                       double maxVols, long nbMaxIter,
                                       ARM_Vector* CorrelatedIndexes,  
                                       ARM_Vector* Indexes , 
                                       ARM_Matrix* Correlations,
                                       ARM_Portfolio* Pf)
                    : ARM_Calibration(asOf, zc)
{
    Init();

    itsSecurity = sec;
    itsVolCurve = VolCurve;

    itsMfine = mfine;
    itsShape_type = shape_type;
    itsDecay = decay;
    itsSlope = slope;
    itsAsymptote = asymptote;
    itsNbFactor = NbFactor;

    if (CorrelatedIndexes)
       itsCorrelatedIndexes = (ARM_Vector *) CorrelatedIndexes->Clone();

    if (Indexes)
       itsIndexes = (ARM_Vector *) Indexes->Clone();

    if (Correlations)
       itsCorrelations = (ARM_Matrix *) Correlations->Clone();

    itsIndicPrix = 0.0;
    itsIndicVol = 0.0; 

    itsPf = Pf;

    itsprecisionVol = precisionVol;
    itsminVols = minVols;
    itsmaxVols = maxVols;
    itsnbMaxIter = nbMaxIter;
}



ARM_CalibrationFRM::ARM_CalibrationFRM(ARM_Date& asOf, ARM_Model *model,
                                       ARM_Security* sec,
                                       int mfine, int shape_type  , 
                                       double decay, double slope, 
                                       double asymptote, int NbFactor , 
                                       double precisionVol, double minVols,
                                       double maxVols, long nbMaxIter,
                                       ARM_Vector* CorrelatedIndexes,  
                                       ARM_Vector* Indexes , 
                                       ARM_Matrix* Correlations,
                                       ARM_Portfolio* Pf)
                    : ARM_Calibration(asOf, model->GetZeroCurve())
{
    Init();

    itsSecurity = sec;
    itsMktModelToFit = model;

    itsMfine = mfine;
    itsShape_type = shape_type;
    itsDecay = decay;
    itsSlope = slope;
    itsAsymptote = asymptote;
    itsNbFactor = NbFactor;

    if (CorrelatedIndexes)
       itsCorrelatedIndexes = (ARM_Vector *) CorrelatedIndexes->Clone();

    if (Indexes)
       itsIndexes = (ARM_Vector *) Indexes->Clone();

    if (Correlations)
       itsCorrelations = (ARM_Matrix *) Correlations->Clone();

    itsIndicPrix = 0.0;
    itsIndicVol = 0.0; 

    itsPf = Pf;

    itsprecisionVol = precisionVol;
    itsminVols = minVols;
    itsmaxVols = maxVols;
    itsnbMaxIter = nbMaxIter;
}



ARM_CalibrationFRM::~ARM_CalibrationFRM(void) 
{
    if (itsPf)
    {
       delete itsPf; 
       itsPf = NULL;
    }

    if (itsCorrelatedIndexes)
    {
       delete itsCorrelatedIndexes;
       itsCorrelatedIndexes =NULL;
    }

    if (itsIndexes)
    {
       delete itsIndexes;
       itsIndexes =NULL;
    }

    if (itsCorrelations)
    {
       delete itsCorrelations;
       itsCorrelations = NULL;
    }
}



ARM_FrmAna* ARM_CalibrationFRM::Calibrate(char* msg,
                                          ARM_Calibration_TYPE caltype,
                                          ARM_CalibrationHWSV_OUTPUT& Output)
{
    ARM_FrmAna* AnaModel = NULL;
    ARM_Model*  bsmod    = NULL;

    // fabrication du portfeuille par l'instrument
    // -------------------------------------------
      try 
    {     
          if (!itsPf)
          {
             //itsSecurity->PrepareToCalibrate(bsmod);
             itsPf = itsSecurity->CalibrationPF(GetAsOfDate(), caltype);
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
          
          // fabrication du model
          itsPf->SetModelVariable(NULL);
          itsPf->SetModel(bsmod);

          double price;
          int PFSize = itsPf->GetSize();
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
                    
    
          // recuperation de la liste de dates
           // ---------------------------------
          double X[200];
          double Y[200];
          double vals_curve[200];

          MEMSET(Y, 0.0, sizeof(double)*200);
          MEMSET(vals_curve, 0.0, sizeof(double)*200);

          if ( itsShape_type == K_ROW )
          {
               for (i = 0; i < PFSize; i++)
               {
                   X[i] = (itsPf->GetAsset(i)->GetExpiryDate()-GetAsOfDate())
                          /365.0+1.0/365.0;

                   Output.datesched[i] = itsPf->GetAsset(i)->GetExpiryDate().AddDays(1);
               }
          }
          else
          {
               for (i = 0; i < PFSize; i++)
               {
                   ARM_Security* curAsset = itsPf->GetAsset(i);
                   ARM_CapFloor* cap = NULL;

                   SmartCast(curAsset, cap);

                   if (cap)
                   {
                      int nbFlows = cap->GetSwapLeg()->GetNumFlows();

                      X[i] = (cap->GetSwapLeg()->GetFlowEndDates()->Elt(nbFlows-1)
                             -GetAsOfDate().GetJulian())/365.0;                   

                      Output.datesched[i] = 
                         cap->GetSwapLeg()->GetFlowEndDates()->Elt(nbFlows-1);               
                    }
                    else
                    {
                       X[i] = (curAsset->GetExpiryDate() - GetAsOfDate())/365.0 + 1.0/365.0;               
                       Output.datesched[i] = curAsset->GetExpiryDate().AddDays(1);
                    }
               }
          }

          // Appel de la fonction de minimisation
          // ------------------------------------    
          AnaModel = new ARM_FrmAna(GetZeroCurve(), itsPf, itsMfine, itsShape_type,
                                    itsDecay, itsSlope, itsAsymptote, 
                                    itsNbFactor , itsCorrelatedIndexes, itsIndexes,
                                    itsCorrelations);


          for (i = 0; i < PFSize; i++)
          {
              itsPf->GetAsset(i)->SetModelVariable(NULL);
              itsPf->GetAsset(i)->SetModel(AnaModel);
          }
   
          int filtreType = 0;

          if ( itsShape_type == K_DIAG )
             filtreType = K_DIAG_TYPE;               

  
          ComputeParasCurve(AnaModel, itsPf, PFSize,
                            X, itsprecisionVol,
                            itsminVols, itsmaxVols, 
                            itsnbMaxIter, vals_curve, filtreType);

          // calcul des indicateurs
          double ecart, prixFRM, prixBS;
          Output.calagequality = 0.0;
          ecart = 0.0;

          for (i = 0; i < PFSize; i++)
          {
               prixBS = mktprice->Elt(i);

               itsPf->GetAsset(i)->SetModelVariable(NULL);
               itsPf->GetAsset(i)->SetModel(AnaModel);

               prixFRM = itsPf->GetAsset(i)->ComputePrice();
               ecart = prixFRM/prixBS-1;
               ecart *= ecart;
               itsIndicPrix += ecart ;
     
               Output.SecurityStartDate[i] = ((ARM_CapFloor *)
                      (itsPf->GetAsset(i)))->GetSwapLeg()->GetStartDate() ;

               Output.SecurityEndDate[i] = ((ARM_CapFloor *)
                      (itsPf->GetAsset(i)))->GetSwapLeg()->GetEndDate();

               Output.SecurityStrike[i] = ((ARM_CapFloor *)
                                           (itsPf->GetAsset(i)))->GetStrike();

               Output.meanrev = itsDecay;
               Output.volsched[i] = vals_curve[i];  

               prixBS = mktprice->Elt(i);

               itsPf->GetAsset(i)->SetModelVariable(NULL);
               itsPf->GetAsset(i)->SetModel(AnaModel);

               prixFRM = itsPf->GetAsset(i)->ComputePrice();
          
               Output.InputPrice[i]  = prixBS;
               Output.InputVol[i]    = -100.0;
               Output.OutputPrice[i] = prixFRM;
               Output.ErrorPrice[i]  = prixFRM-prixBS; 
               Output.OutputVol[i]   = -100.0;
               Output.ErrorVol[i]    = 0.0;
          }

          for (i = 0; i < PFSize; i++)
          {
               Output.InputVol[i]= 
                   itsPf->GetAsset(i)->ComputeImpliedVol(Output.InputPrice[i]);

               Output.OutputVol[i]= 
                   itsPf->GetAsset(i)->ComputeImpliedVol(Output.OutputPrice[i]);

               Output.ErrorVol[i]= Output.OutputVol[i] - Output.InputVol[i] ;
               Output.calagequality += Output.ErrorVol[i] * Output.ErrorVol[i] ;
          }

          itsIndicPrix  = sqrt(itsIndicPrix);
          itsIndicPrix /= PFSize;

          if (bsmod)       
             delete bsmod ;

          for (i = 0; i < PFSize; i++)
          {
              if (itsPf->GetAsset(i))
                 delete itsPf->GetAsset(i);
          }

          delete itsPf;
          itsPf = NULL;
          
     }// end try
	catch(Exception&e)
	{
		std::string m=e.GetErrorString(); 
		if (msg) strncpy(msg,m.c_str(),ARM_EXCEPTION_MSG_MAX_SIZE); 
	}

     return(AnaModel);
}


/*-----------------------------------------------------------------------*/
/*----- End Of File ----*/
