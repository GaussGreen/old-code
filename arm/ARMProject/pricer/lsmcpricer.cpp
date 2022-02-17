/*
 * $Log: lsmcpricer.cpp,v $
 * Revision 1.19  2003/07/01 21:09:57  arm
 * long i added instead of for (long i;
 *
 * Revision 1.18  2003/06/30 08:57:27  ebenhamou
 * remove unused var
 *
 * Revision 1.17  2003/02/03 10:20:28  sgasquet
 * Modifs pour prise en compte OptionOnPortfolio et modele LogDecale
 *
 * Revision 1.16  2002/03/01 16:51:50  mab
 * MODIF : YLI Corrections
 *
 * Revision 1.15  2001/09/06 17:15:11  smysona
 * correction discount basis pour variable de controle
 *
 * Revision 1.14  2001/08/09 18:15:43  smysona
 * sort the calibPF by price
 *
 * Revision 1.13  2001/07/30 09:11:36  smysona
 * modifs pour basis discount.
 *
 * Revision 1.12  2001/05/23 17:44:19  smysona
 * Ajour de end / begin pricing
 *
 * Revision 1.11  2001/03/26 15:09:19  nicolasm
 * Initialisation du modele au debut
 * Corrections des leaks
 * Optimisation des perfs
 * makeup general
 *
 * Revision 1.10  2001/03/22 15:36:42  nicolasm
 * Modif de la boucle d'aggregation des prix
 *
 * Revision 1.9  2001/03/21 16:46:46  mab
 * Remontee chgt CN pour les CRF
 *
 * Revision 1.8  2001/03/20 16:32:42  sgasquet
 * Modif Prise en compte var de controle
 *
 * Revision 1.7  2001/03/15 19:01:01  smysona
 * Correction adressage si pas de variable de controle
 *
 * Revision 1.6  2001/03/12 19:37:31  smysona
 * Acceleration
 *
 * Revision 1.5  2001/03/05 10:14:30  smysona
 * Control variables are freed only if they have been allocated
 *
 * Revision 1.4  2001/02/22 19:42:20  smysona
 * SpeedUp + ajout ControleVariate
 *
 * Revision 1.3  2000/11/13 13:31:22  sgasquet
 * Retrait headers
 *
 * Revision 1.2  2000/10/24 12:49:13  sgasquet
 * First revision
 *
 * Revision 1.1  2000/10/24 12:20:58  sgasquet
 * Initial revision
 *
 */


#include "lsmcpricer.h"
#include "portfolio.h"
#include "frmana.h"
#include "smc_frm.h"
// #include "stats.h"
// #include "linalgpp.h"
// #include "mcentrop.h"

void ARM_LSMonteCarloPricer::Init(void)
{
    itsMode = 0;
}




void ARM_LSMonteCarloPricer::EndPricing(void)
{
    
    ARM_Frm* frm = NULL;

    SmartCast(itsModel, frm);

    if(frm)
    {
        if (frm->GetPriceOutpoutFilename())
        {
            ofstream fout(frm->GetPriceOutpoutFilename(), ios::app);
            time_t end;

            time(&end);
            double elapsed = difftime(end, GetTime());

            fout << endl << endl;
            fout << "Pricing time : " << elapsed << " (s)" << endl;

            fout.close();
        }
    }
    
}



// if returns 1: error
// Y: FittedVector   X = Value    B = OutParams   Y = X B + eps
// B = (X'X)^(-1) X'Y 

long ARM_LSMonteCarloPricer::ComputeRegression(	ARM_Vector* OutParams, ARM_Vector* FittedVector, 
			ARM_Matrix* Value, long nbRelevantInputs, long nbRelevantFactors)
{
	if (OutParams->GetSize() != nbRelevantFactors)
		return 1;
	if (FittedVector->GetSize() < nbRelevantInputs)
		return 1;
	if ( (Value->GetNumCols() < nbRelevantFactors) || (Value->GetNumLines() < nbRelevantInputs) )
		return 1;

	double det;

	// We compute (X'X)
	ARM_Matrix* Xprime = new ARM_Matrix(nbRelevantFactors, nbRelevantInputs);
	ARM_Matrix* X = new ARM_Matrix(nbRelevantInputs, nbRelevantFactors);
		
        long i;
	for (i = 0 ; i < nbRelevantFactors; i++)		// l colonnes
	{
		for (long j = 0 ; j < nbRelevantInputs; j++)	// p lignes
		{
			Xprime->Elt(i, j) = Value->Elt(j, i);
			X->Elt(j, i) = Value->Elt(j, i);
		}
	}
	ARM_Matrix XPrimeX = (*Xprime) * (*X) ;

	// We compute X'Y
	ARM_Matrix* Y = new ARM_Matrix(nbRelevantInputs, 1);

	for ( i = 0 ; i < nbRelevantInputs; i++)			// l lignes
	{
		Y->Elt(i, 0) = FittedVector->Elt(i);
	}

	ARM_Matrix XprimeY = (*Xprime) * (*Y) ;
		
	for ( i = 0 ; i < nbRelevantFactors; i++)			// l lignes
	{
		OutParams->Elt(i) = XprimeY.Elt(i, 0);
	}

	// We solve X'Y = X'X B
	XPrimeX.LinSolve(OutParams, det);

    cleanit(Xprime);
    cleanit(X);
    cleanit(Y);

	if (det!=0)
		return 0;
	else
		return 1;
}





/***************************************************************/
/*                                                             */
/*            THE price method                                 */
/*                                                             */
/***************************************************************/


double ARM_LSMonteCarloPricer::Price(void)
{


try
{
    StartPricing();


    //Initial Setup
    ARM_SMCFrm* frmmc;
    SmartCast(itsModel, frmmc);
	long IsFRMModel = 0;

    if (frmmc)
    {
        frmmc->InitTrajectory(0);
		itsSecurity->SetModel(itsModel);
		IsFRMModel = 1;
    }

   long i, j, k, l, p, kMax;
   ARM_Date calcDate = itsModel->GetStartDate();
   

   itsSecurity->PrepareToPrice(calcDate);                    
   ARM_ExerciseStyle* securExerciseStyle = itsSecurity->GetExerciseStyle();

   ARM_MCModel *mcmod = dynamic_cast<ARM_MCModel *> (itsModel);    

   mcmod->BeFittedTo(itsSecurity);                         

   // YAN 02/2002
   itsSecurity->ComputeUnderlyingPrice(calcDate, mcmod);

   long nbTimeSteps = mcmod->NbSteps();
   long nbPaths = mcmod->NbIters();


//////////////////////////////////////////////////////////////
// We define and allocate the matrices if possible
// We compute the cash flows at each date k
// We compute the regressors value only at the exercise date
// If there is no option exercise, the cash flows should
// be backwarded until the first period in the function
// ComputeBackwardCashFlowsGrid
// In this case, nbExercisesDates=0 and nbExercises=1
//////////////////////////////////////////////////////////////
	long RankForFirstDate = mcmod->GetRankForFirstTimeStep(); // 1 for FRM et 0 for LD

    ARM_Matrix* CashFlows = new ARM_Matrix( nbPaths, nbTimeSteps+1-RankForFirstDate);
    ARM_Matrix* DiscountPeriod = new ARM_Matrix( nbPaths, nbTimeSteps+1);

    // YAN 02/2002
	ARM_Matrix* CvDiscountPeriod = NULL;
	if (mcmod->IsBasisDiscountModel())
		CvDiscountPeriod = new ARM_Matrix( nbPaths, nbTimeSteps+1);
	else
		CvDiscountPeriod = DiscountPeriod;
    // Fin YAN 02/2002

    ARM_Vector* ValueOption = new ARM_Vector( nbPaths );
    ARM_Vector* IntrinValue = new ARM_Vector( nbPaths );
    int         nbExercisesDates = 0;
    int         nbExercises = 0;
    int         nbRegressionFactors;
    int         nbComputedExercises;
    double***   RegressorValues; 
    int*        ExerciseRank;
    int***      OptimalExDates; 
    ARM_Date*   currDate = NULL;


    if (securExerciseStyle)
       nbExercisesDates = securExerciseStyle->GetExerciseDatesSize();
    else
       nbExercisesDates = 0;


    k = 0; 
    kMax = 0;


    while ( ( k < nbExercisesDates )
            &&
            ( calcDate.GetJulian() 
              > securExerciseStyle->GetExerciseEndDates()->Elt(k) ) 
          )
    {    
        nbExercisesDates --;
        k++;
    }

    nbExercises = itsSecurity->NbExercises(); // nb cartouches = 1 by default


    /////////////////////////////////////////////////////////////////////
    // We allocate RegressorValues if there is at least one exerciseDates
    /////////////////////////////////////////////////////////////////////

    RegressorValues = NULL;

    if (nbExercisesDates > 0)
    {
        RegressorValues = new double**[nbExercisesDates];
        l = 0;

        for ( k = 0; k < nbTimeSteps+1-RankForFirstDate; k++)
        {          
			if (k == 0 ) 
				currDate = new ARM_Date(calcDate);
			else
				currDate = mcmod->FromTimeStepToDate(k-1+RankForFirstDate);    

            if ( (securExerciseStyle)
                 &&
                 (securExerciseStyle->IsExerciseDate(*currDate)))
            {
                long nb = itsSecurity->nbRegressionFactors( *currDate );

                if (nb > 0)
                {
                    RegressorValues[l] = InitializeDoubleMatrix(nbPaths, nb);
                    l++;
                }
                else
                {
                    RegressorValues[l] = NULL;
                }
            }
            else
            {        
                RegressorValues[l] = NULL;
            }
			if ( (currDate) && (k==0) ) delete currDate;
        }
    }


    /////////////////////////////////////////////////////////////////////
    // Controle variate 
    /////////////////////////////////////////////////////////////////////


    //Temporary code : assume we use an FRM model
	ARM_Portfolio *controleVariates = NULL;
	ARM_FrmAna* frmana = NULL;

    if (IsFRMModel)
    {
		frmana = new ARM_FrmAna(); 
		frmana->DummyCopy(mcmod);
	

		controleVariates = 
        itsSecurity->CreateControleVariate(frmana, CALIB_PORTFOLIO_MODE::Sort);
	}
    // DEBUG
    // POUR SUPPRIMER LE CONTROLE
    // cleanit(controleVariates);


    ARM_Security *controleVariate = NULL;
    ARM_Matrix *CvCashFlows = NULL;
    int CvExercise = /*-1*/0;

    if (controleVariates)
    {
        controleVariate = controleVariates->GetAsset(0);
    }


    if (controleVariate)
    {
        controleVariate->PrepareToPrice(calcDate);
        CvCashFlows = new ARM_Matrix(nbPaths, nbTimeSteps);
    }



    /////////////////////////////////////////////////////////////////////
    // We allocate and compute rank exercises
    /////////////////////////////////////////////////////////////////////

    if (securExerciseStyle)
    {
        ExerciseRank = new int[nbExercisesDates];    
        l = 0;
        for ( k = nbTimeSteps-RankForFirstDate; k >=0 ; k--)    
        {
			if (k == 0 ) 
				currDate = new ARM_Date(calcDate);
			else
				currDate = mcmod->FromTimeStepToDate(k-1+RankForFirstDate);    

            if (securExerciseStyle->IsExerciseDate( currDate->GetJulian() ))                 
            {
                ExerciseRank[l] = k;
                l++;                
            }        
            if (controleVariate)
            {
                if (controleVariate->GetExerciseStyle()->IsExerciseDate( currDate->GetJulian() ))                 
                {
                    CvExercise = k;
                }

            }
			if ( (currDate) && (k==0) ) delete currDate;
        }
    }
    else
    {
        ExerciseRank = new int[1];    
        ExerciseRank[0] = 0; // we take only the first cash flows 
    }
    

    /////////////////////////////////////////////////////////////////////
    // We allocate and initialize optimal exercises
    /////////////////////////////////////////////////////////////////////

    OptimalExDates = new int**[nbExercises];

    for ( l = 0; l < nbExercises; l++)
    {
        OptimalExDates[l] = InitializeIntMatrix(l+1, nbPaths);

        for ( i = 0; i <= l; i++)
        {
            for ( p = 0; p < nbPaths; p++)
            {
                OptimalExDates[l][i][p] = ExerciseRank[i];
            }
        }
    }



    /////////////////////////////////////////////////////////////////////
    // We build the matrix of cash flows
	// We compute DF and CashFlow at startdate 
	// As a result DF(0) = Df between asOf and first date for schedule
	// For FRM as the first date of the schedule is indeed the asOf, no pbs
	// For LD, we need to create this date and shift all the date by -1 with 
	// the parameter RankForFirstDate
     /////////////////////////////////////////////////////////////////////
    nbComputedExercises = 0;

    ARM_Date * pDate = mcmod->FromTimeStepToDate(RankForFirstDate);
    double firstDf = mcmod->GetDiscountCurve()->DiscountPrice(*pDate);

    // YAN 02/2002
	mcmod->SetBasisDiscountFlag(0);
	double CvfirstDf = mcmod->GetDiscountCurve()->DiscountPrice(*pDate);
	mcmod->SetBasisDiscountFlag(1);
   // FIN YAN 02/2002

    //Speed up
    double* DiscountPeriodElt = NULL;   
    // YAN 02/2002
	double* CvDiscountPeriodElt = NULL;
    // Fin YAN 2002

    double* CashFlowsElt = NULL;
    double* CvCashFlowsElt = NULL;

    for ( p = 0; p < nbPaths; p++)
    {
        DiscountPeriodElt = DiscountPeriod->GetDbleLine(p);
        // YAN 02/2002
		CvDiscountPeriodElt = CvDiscountPeriod->GetDbleLine(p);
       // Fin YAN 02/2002
        CashFlowsElt      = CashFlows->GetDbleLine(p);
		mcmod->InitTrajectory(p);

        if (controleVariate)
           CvCashFlowsElt = CvCashFlows->GetDbleLine(p);    

		if (IsFRMModel)
			DiscountPeriodElt[0] =  firstDf;
		else
			DiscountPeriodElt[0] =  mcmod->FirstPeriodDiscountFactor(p);  //firstDf;

        // YAN 02/2002
		CvDiscountPeriodElt[0] = DiscountPeriodElt[0];
       // Fin YAN

        for ( k = 0; k < nbTimeSteps+1-RankForFirstDate; k++)
        {
			if (k == 0 ) 
				currDate = new ARM_Date(calcDate);
			else
				currDate = mcmod->FromTimeStepToDate(k-1+RankForFirstDate);    

			if (!itsSecurity->IsComputationRegressorsAutomatic(*currDate))				
			{
				CashFlowsElt[k] = itsSecurity->ComputePayOffForGrid(*currDate, mcmod);
			}
			else
			{
				CashFlowsElt[k] = itsSecurity->ComputeUnderlyingPayOffForGrid(*currDate, mcmod);
			}

            if (controleVariate)
            {
               // YAN 02/2002
			   mcmod->SetBasisDiscountFlag(0);
               CvCashFlowsElt[k] = controleVariate->ComputePayOffForGrid(*currDate, mcmod);
               // YAN 02/2002
			   mcmod->SetBasisDiscountFlag(1);
            }

            nbRegressionFactors = itsSecurity->nbRegressionFactors(*currDate);
            
            if ( (securExerciseStyle) 
                 && (securExerciseStyle->IsExerciseDate( *currDate )) 
                 && ( nbRegressionFactors > 0) ) //k < ExerciseRank[0]))
            {
				if (!itsSecurity->IsComputationRegressorsAutomatic(*currDate))
				{
					itsSecurity->ComputeRegressorPayOffForGrid(*currDate, mcmod,
									 RegressorValues[nbComputedExercises][p]);
				}
				else
				{
					RegressorValues[nbComputedExercises][p][0] = 1.0;
					for (long n = 1 ; n < nbRegressionFactors; n++)
						RegressorValues[nbComputedExercises][p][n] = pow(CashFlowsElt[k] , n);					
				}
                nbComputedExercises++;
            }

            // CashFlows should be computed for all k and ranked from nbTimeSteps-1 to 0    
            // RegressorValues should be computed only for exercise dates from nbExercisesDates-1 to 0
			if (itsSecurity->IsComputationRegressorsAutomatic(*currDate))
			{
				CashFlowsElt[k] = MAX(CashFlowsElt[k], 0);
			}

            if (k < nbTimeSteps-1)
            {
                mcmod->PropagateCurve(p, k+1); //1 premiere diffusion

                DiscountPeriodElt[k+1] =  mcmod->FirstPeriodDiscountFactor(p);        
 
                // YAN 02/2002
				mcmod->SetBasisDiscountFlag(0);
				CvDiscountPeriodElt[k+1] =  mcmod->FirstPeriodDiscountFactor(p);
				mcmod->SetBasisDiscountFlag(1);
                // Fin YAN 
            }
			if ( (currDate) && (k==0) ) delete currDate;
        }

        itsSecurity->ComputeBackwardCashFlowsGrid(p, k, mcmod, CashFlows, DiscountPeriod);
        if (controleVariate)
        {
            // YAN 02/2002
            controleVariate->ComputeBackwardCashFlowsGrid(p, k, mcmod, CvCashFlows, CvDiscountPeriod);


            // OLDcontroleVariate->ComputeBackwardCashFlowsGrid(p, k, mcmod, CvCashFlows, DiscountPeriod);
        }


        nbComputedExercises = 0;
    }




////////////////////////////////////////////////////////////////////////
// We compute the exercise strategy only if there is an exercise choice
////////////////////////////////////////////////////////////////////////

if ( nbExercises < nbExercisesDates)
{
    long nbComputedExercises = 0;

    // We start the regression for each exercise date and we go until date1
    // The last actualisation is done at date 1 for date 0
    for  ( k = nbTimeSteps  - RankForFirstDate ; k > 0 ; k--)
    {    

        currDate = mcmod->FromTimeStepToDate(k-1+RankForFirstDate);

        if (  securExerciseStyle->IsExerciseDate(*currDate) ) 
        // it's an Exercise date and we need to regress    
        // because all the cartouche have been potentially
        // exercised
        {
            // We regress independently each right to postpone the lth exercise
            for ( l = MIN(nbExercises, nbComputedExercises)-1; l >= 0; l--)        
            {
                /////////////////////////////////////////////
                // We compute the regression of the pay offs
                /////////////////////////////////////////////

                long NbRelevantCashFlows = 0;

                nbRegressionFactors = itsSecurity->nbRegressionFactors(*currDate);

                ARM_Matrix* FactorVector = new ARM_Matrix(nbPaths, nbRegressionFactors );
                ARM_Vector* Value = new ARM_Vector(nbPaths); 

                // we will only use the NbRelevantCashFlows 1st elts
                long* PathRank = new long[nbPaths];

                for (p = 0; p < nbPaths; p++) 
                {
                    if ( itsSecurity->IsRelevantForRegression(CashFlows->Elt(p, k ) ) ) 
                         // Positive by default
                    {
                        PathRank[NbRelevantCashFlows] = p;
                                
                        Value->Elt(NbRelevantCashFlows) = 0.0;
                        // Value is the value of the option if all the exercises
                        for ( i = 0; i <= l; i++)
                        {
                            Value->Elt(NbRelevantCashFlows) +=
                                CashFlows->Elt( p, OptimalExDates[l][i][p]); 
                        }
                
                        // Regression variables                        
                        for ( j = 0; j < nbRegressionFactors; j++)
                        {
                            FactorVector->Elt(NbRelevantCashFlows, j) =
                            RegressorValues[nbExercisesDates-nbComputedExercises-1][p][j]; 
                            
                            // we use the same regressor for the different pay-offs
                            // nbComputedExercises: rank of the Regressor value
                        }    
            
                        NbRelevantCashFlows++;
                    }
                }

                if (NbRelevantCashFlows > 0) 
                // Else: we don't exercise; nothing changes
                {
                    ARM_Vector* Params = new ARM_Vector(nbRegressionFactors);

                    int test = ComputeRegression(Params, Value, FactorVector, NbRelevantCashFlows,
                           nbRegressionFactors);
                
                    if (test)
                    {
                        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                                "Pb in Regression process")    ;
                    }
                    

            ///////////////////////////////////////////////////////////////////////
            // We compute the value of the option and the intrinseque value
            // We check if we exercise or not only on the path relevant for
            // the regression
            ///////////////////////////////////////////////////////////////////////

                    for (j = 0; j < NbRelevantCashFlows; j++)
                    {            
                        long p = PathRank[j];
                        ValueOption->Elt(p) = .0 ; 
                        IntrinValue->Elt(p) = CashFlows->Elt( p, k); 

                        for ( i = 0; i < nbRegressionFactors; i++)
                        {                    
                            ValueOption->Elt(p) +=    Params->Elt(i)
                            *RegressorValues[nbExercisesDates-nbComputedExercises-1][p][i];
                        }

                        for ( i = 0; i <= l-1; i++)
                        {
                            long tmp = OptimalExDates[l-1][i][p];
                            IntrinValue->Elt(p) += CashFlows->Elt( p, OptimalExDates[l-1][i][p]); 
                        }

                        if ( IntrinValue->Elt(p) > ValueOption->Elt(p))  
                        // We exercise
                        {
                            OptimalExDates[l][0][p] = k;
                            for ( i = 0; i <= l-1; i++)
                                OptimalExDates[l][i+1][p] = OptimalExDates[l-1][i][p];

                        }    
                    }

					cleanit(Params);
                }

                cleanit(FactorVector);
				cleanit(Value);
                cleanit(PathRank);
                
            }  // end loop on l

            nbComputedExercises++;
        } // end condition is exercise date or not


    ///////////////////////////////////////////////////////////////////////
    // We actualize the payoff of the option until the next period
    // from k-1 to k
    ///////////////////////////////////////////////////////////////////////

        for ( p = 0; p < nbPaths; p++)
        {
            CashFlowsElt      = CashFlows->GetDbleLine(p);
        if (controleVariate)
          CvCashFlowsElt    = CvCashFlows->GetDbleLine(p);
            DiscountPeriodElt = DiscountPeriod->GetDbleLine(p);
            // YAN 02/2002
			CvDiscountPeriodElt = CvDiscountPeriod->GetDbleLine(p);

            for (long t = nbTimeSteps - RankForFirstDate ; t > k-1; t--)
            {
                //CashFlows->Elt( p, t) = CashFlows->Elt( p, t) * DiscountPeriod->Elt(p, k-1);
                CashFlowsElt[t] = CashFlowsElt[t] * DiscountPeriodElt[k-1];
            }

            if (controleVariate)
            {
                for (long t = nbTimeSteps - RankForFirstDate ; t > k-1; t--)
                {
                    //CvCashFlows->Elt( p, t) = CvCashFlows->Elt( p, t) * DiscountPeriod->Elt(p, k-1);
                    // YAN 02/2002 remplacement DiscountPeriodElt par CvDiscountPeriodElt
                    CvCashFlowsElt[t] = CvCashFlowsElt[t] * CvDiscountPeriodElt[k-1];
                }                
            }

        }
    } // end loop on time steps
}

 ///////////////////////////////////////////////////////////////////////////////////////
 //  There is as many cartouche as exercise dates or no exercise dates
 //  we only actualize the pay off from time n-2 (period n-2 to n-1) to 0 (period 0 to 1)
 //  n-1 is the rank of the last date
 ///////////////////////////////////////////////////////////////////////////////////////
else
{
    for  ( k = nbTimeSteps - 1 - RankForFirstDate ; k >=0 ; k--)
    {
        for ( p = 0; p < nbPaths; p++)
        {
            CashFlowsElt      = CashFlows->GetDbleLine(p);
            DiscountPeriodElt = DiscountPeriod->GetDbleLine(p);
            // YAN 02/2002
			CvDiscountPeriodElt = CvDiscountPeriod->GetDbleLine(p);

            for (long t = nbTimeSteps - RankForFirstDate ; t > k; t--)
            {
                //CashFlows->Elt( p, t) = CashFlows->Elt( p, t) * DiscountPeriod->Elt(p, k);
                CashFlowsElt[t] = CashFlowsElt[t] * DiscountPeriodElt[k];
            }

            if (controleVariate)
            {
            CvCashFlowsElt    = CvCashFlows->GetDbleLine(p);

                for (long t = nbTimeSteps - RankForFirstDate ; t > k; t--)
                {
                    // YAN 02/2002 remplacement de DiscountPeriodElt par CvDiscountPeriodElt
                    CvCashFlowsElt[t] = CvCashFlowsElt[t] * CvDiscountPeriodElt[k];
                }

            }

        }
    }

}
            
            
///////////////////////////////////////////////////////////////////////
// We can now compute the option value and average it
///////////////////////////////////////////////////////////////////////
    
    ARM_Vector *Price = new ARM_Vector(nbPaths, 0.0);
    double* PriceElt = Price->GetElt();
    double Average = 0.0 , StdDev = 0.0;


    for (p = 0 ; p < nbPaths ; p++)
    {
        for  (i = nbExercises - 1 ; i >=0 ; i--)
        {
            PriceElt[p] += CashFlows->Elt(p, OptimalExDates[nbExercises-1][i][p] );
        }
	}

        if (controleVariate)
        {
        for (p = 0 ; p < nbPaths ; p++)
        {
            PriceElt[p] -= CvCashFlows->Elt(p,CvExercise);
        }
    }

    for (p = 0 ; p < nbPaths ; p++)
    {   
        Average += PriceElt[p];
        StdDev  += PriceElt[p] * PriceElt[p];

    }



    Average /= nbPaths;
    StdDev  /= nbPaths;
    StdDev  = sqrt(StdDev - Average * Average ); 
    SetEstimatedSdev(StdDev);

    double cvPrice = 0;
    

    if (controleVariate)
    {
        //DIRTY, should be put in a method of the security
        //Shift the price for basis discount 
        //Remember : cashflows are computed using the real curve, discountCurve is only used
        //to discount the CFs
        cvPrice = controleVariates->GetMktPrices()->Elt(0);
        ARM_Date expiryDate(controleVariates->GetAsset(0)->GetExpiryDate());    
 // YAN 02/2002 Mise en commentaire    
/*
        cvPrice *=  mcmod->GetDiscountCurve()->DiscountPrice(expiryDate)
                  / mcmod->GetZeroCurve()->DiscountPrice(expiryDate);
*/
// Fin YAN
        Average += cvPrice;
    }

	((ARM_CapFloor*)itsSecurity)->GetSwapLeg()->SetCashFlowValue(0,StdDev);


///////////////////////////////////////////////////////////////////////
// We free the matrices
///////////////////////////////////////////////////////////////////////

    delete CashFlows;
    delete DiscountPeriod;
    // YAN 02/2002
	if (mcmod->IsBasisDiscountModel())
		delete CvDiscountPeriod;
     // Fin YAN 02/2002
    delete ValueOption;
    delete IntrinValue;
    delete ExerciseRank;

    l = 0;


    if (RegressorValues)
    {
        for ( k = 0; k < nbTimeSteps+1-RankForFirstDate; k++)
        {
			if (k == 0 ) 
				currDate = new ARM_Date(calcDate);
			else
				currDate = mcmod->FromTimeStepToDate(k-1+RankForFirstDate);  

			if ( (securExerciseStyle) && securExerciseStyle->IsExerciseDate(*currDate) )
            {
                long nb = itsSecurity->nbRegressionFactors(*currDate);            
                if (RegressorValues[l] != NULL)
                {
                    FreeDoubleMatrix(RegressorValues[l] , nbPaths);
                    l++;
                }
            }
			if ( (currDate) && (k==0) ) delete currDate;
        }
        delete RegressorValues;
    }

    
    for ( l = 0; l < nbExercises; l++)
    {
        FreeIntMatrix(OptimalExDates[l], l+1);
    }


    delete OptimalExDates;
    delete Price;

    cleanit(CvCashFlows);
    if (controleVariates)
      controleVariates->FreePortfolioAndAssets();
    cleanit(controleVariates);


    if (frmana)
    {
		frmana->CleanDummyCopy();
		cleanit(frmana);
	}

    EndPricing();


    return Average;
}

catch (Exception& m)
{
    //One day we will include #ifdef to check wether we are
    //compiled for summit or ARM
    cerr << "Exception caught in ARM pricing !" << endl;

    // char errMsg[50];
    // m.GetErrorMessage(errMsg);
    cerr << m.GetErrorString() << endl;

    return(0);
}
catch (...)
{
    cerr << "Something caught in ARM pricing !" << endl;
    return(0);
}

    return(0);
}


