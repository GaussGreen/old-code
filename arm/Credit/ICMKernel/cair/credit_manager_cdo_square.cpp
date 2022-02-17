#include "ARMKernel\glob\firsttoinc.h"
/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		CREDIT_MANAGER_CDO_SQUARE.CPP
	PROJECT:	CAIR
	
	DESCRIPTION:	implementation based on Monte-Carlo for CDO^2 Pricing

   -----------------------------------------------------------------
   
	CAIR Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */


#include "ICMKernel\cair\credit_manager.h"
#include "ICMKernel\glob\icm_maths.h"

#include <numeric>  // accumulate


// -------------------------------------------------------
// PREPARE CDO SQUARE DATA from Excel Interface
// -------------------------------------------------------
void	CreditManager::PrepareCDOSquareData()
{
	int	i, j;
	int	ii, jj;

	double	CDO_Notional;

	// only 0 or 1!
	ii = 0;
	its_NbRelevant_CDO_Underlyings	=	accumulate(its_CDO_Square_CDO_Includes.begin(), its_CDO_Square_CDO_Includes.end(), 0);

	its_CDO_Square_Relevant_CDO_Ids.clear();
	its_CDO_Square_Relevant_CDO_MaturitiesInYearFraction.clear();
	its_CDO_Square_Relevant_CDO_NbCredits.clear();
	its_CDO_Square_Relevant_CDO_Notionals.clear();
	its_CDO_Square_Relevant_CDO_Rescaling_Factor.clear();

	// to be optimized
	its_CDO_Square_Relevant_Credit_Ids 		= new ICM_QMatrix<int>(Get_NbCredits(), its_NbRelevant_CDO_Underlyings);
	its_CDO_Square_Relevant_Credit_Includes = new ICM_QMatrix<bool>(Get_NbCredits(), its_NbRelevant_CDO_Underlyings);
	its_CDO_Square_Relevant_Credit_Losses 	= new ICM_QMatrix<double>(Get_NbCredits(), its_NbRelevant_CDO_Underlyings);

	int		*Current_Ids;
	double	*Current_Losses;
	int		Current_Nb_Credits;
	int		TheId;

	itsCDOSquareNotional	=	0.0;

	for (i=0; i<its_CDO_Square_Nb_Underlyings; i++)
	{
		if (its_CDO_Square_CDO_Includes[i]) // ??? CHECK EQUAL ???
		{
			// included
			its_CDO_Square_Relevant_CDO_Ids.push_back(i);
			its_CDO_Square_Relevant_CDO_MaturitiesInYearFraction.push_back(MATHTIME(its_CDO_Square_CDO_Maturities[i]));

			its_CDO_Square_Relevant_CDO_Attachments_Pct.push_back(its_CDO_Square_CDO_Attachments[i]);
			its_CDO_Square_Relevant_CDO_Detachments_Pct.push_back(its_CDO_Square_CDO_Detachments[i]);

			Current_Nb_Credits	=	its_CDO_Square_CDO_NbCredits[i];
			its_CDO_Square_Relevant_CDO_NbCredits.push_back(Current_Nb_Credits);

			// Now deal with Credits
			Current_Ids		=	new	int [Current_Nb_Credits];
			Current_Losses	=	new	double [Current_Nb_Credits];

			// loop
			CDO_Notional	=	0.0;

			for (j=0, jj=0; j<Get_NbCredits(); j++)
			{
				TheId	=	(*its_CDO_Square_Credit_Ids)(j, i);
				
				// discarded credits have -1 Flag
				if (TheId > 0)
				{
					Current_Ids[jj]		=	(*its_CDO_Square_Credit_Ids)(j, i) - 1;
					Current_Losses[jj]	=	(*its_CDO_Square_Credit_Losses)(j, i);
					jj++;
					(*its_CDO_Square_Relevant_Credit_Includes)(j, i)	=	true;
				}
				else
					(*its_CDO_Square_Relevant_Credit_Includes)(j, i)	=	false;

				// Notional
				CDO_Notional	+= (*its_CDO_Square_Credit_Notionals)(j, i);
			}

			// CDO Notional, and Attachment, Detachment in Notionals
			its_CDO_Square_Relevant_CDO_Notionals.push_back(CDO_Notional);
			its_CDO_Square_Relevant_CDO_Rescaling_Factor.push_back(its_CDO_Square_CDO_Notionals[i] / CDO_Notional);

			its_CDO_Square_Relevant_CDO_Attachments_Notio.push_back(its_CDO_Square_CDO_Attachments[i] * CDO_Notional);
			its_CDO_Square_Relevant_CDO_Detachments_Notio.push_back(its_CDO_Square_CDO_Detachments[i] * CDO_Notional);

			// Add this Column in the QMatrix Structure
			its_CDO_Square_Relevant_Credit_Ids->SetCol(ii, Current_Ids);
			its_CDO_Square_Relevant_Credit_Losses->SetCol(ii, Current_Losses);

			// One more Underlying CDO 
			ii++; // --> accumulate works fine

			itsCDOSquareNotional	+=	its_CDO_Square_CDO_Notionals[i];	//CDO_Notional;

			// free memory
			delete[] Current_Ids;
			delete[] Current_Losses;
		}
	}
}


// -------------------------------------------------------
// DISPLAY CDO SQUARE DATA
// -------------------------------------------------------
void	CreditManager::DisplayCDOSquareData()
{
	
	int	i, j;
//    FILE*	fOut;
	int	NbCredits;

    char	fOutName[100]	=	"C:\\test\\CDO_Square_Data.txt";
 
    if ((its_fOut = fopen(fOutName, "w")) == NULL)
		ICMTHROW(ERR_INVALID_DATA,"Unable to write CDO Square in file! Check C:\test directory.");

	fprintf(its_fOut, "\n\n------------------------------------------------------------------\n");
	fprintf(its_fOut, "---------------- Credit Manager - CDO Square Data -----------------------------\n");
	fprintf(its_fOut, "------------------------------------------------------------------\n\n\n");

	fprintf(its_fOut, "---------------- For Each Underlying CDO: CDO Included Id -----------------------------\n");
	
	for (i=0; i<its_NbRelevant_CDO_Underlyings; i++)
		fprintf(its_fOut, "%lu\t", its_CDO_Square_Relevant_CDO_Ids[i]);
	fprintf(its_fOut, "\n");

	fprintf(its_fOut, "---------------- For Each Underlying CDO: Maturities of CDO -----------------------------\n");
	
	for (i=0; i<its_NbRelevant_CDO_Underlyings; i++)
		fprintf(its_fOut, "%f\t", its_CDO_Square_Relevant_CDO_MaturitiesInYearFraction[i]);
	fprintf(its_fOut, "\n");

	fprintf(its_fOut, "---------------- For Each Underlying CDO: Attachments of CDO -----------------------------\n");
	
	for (i=0; i<its_NbRelevant_CDO_Underlyings; i++)
		fprintf(its_fOut, "%f\t", its_CDO_Square_Relevant_CDO_Attachments_Pct[i]);
	fprintf(its_fOut, "\n");

	fprintf(its_fOut, "---------------- For Each Underlying CDO: Detachments of CDO -----------------------------\n");
	
	for (i=0; i<its_NbRelevant_CDO_Underlyings; i++)
		fprintf(its_fOut, "%f\t", its_CDO_Square_Relevant_CDO_Detachments_Pct[i]);
	fprintf(its_fOut, "\n");

	fprintf(its_fOut, "---------------- For Each Underlying CDO: Number of Credits -----------------------------\n");
	
	for (i=0; i<its_NbRelevant_CDO_Underlyings; i++)
		fprintf(its_fOut, "%u\t", its_CDO_Square_Relevant_CDO_NbCredits[i]);
	fprintf(its_fOut, "\n");

	fprintf(its_fOut, "---------------- For Each Underlying CDO: Rescaling Factor of CDO -----------------------------\n");
	
	for (i=0; i<its_NbRelevant_CDO_Underlyings; i++)
		fprintf(its_fOut, "%f\t", its_CDO_Square_Relevant_CDO_Rescaling_Factor[i]);
	fprintf(its_fOut, "\n");

	fprintf(its_fOut, "---------------- For Each Underlying CDO: Credits Included Ids -----------------------------\n");
	for (i=0; i<its_NbRelevant_CDO_Underlyings; i++)
	{
		NbCredits	=	its_CDO_Square_Relevant_CDO_NbCredits[i];
		fprintf(its_fOut, "CDO %u\t\t", i);
		for (j=0; j<NbCredits; j++)
			fprintf(its_fOut, "%u\t", (*its_CDO_Square_Relevant_Credit_Ids)(j,i));
		fprintf(its_fOut, "\n");
	}
	
	int Flag;
	fprintf(its_fOut, "---------------- For Each Underlying CDO: Credits Included Include Flag -----------------------------\n");
	for (i=0; i<its_NbRelevant_CDO_Underlyings; i++)
	{
		NbCredits	=	its_CDO_Square_Relevant_CDO_NbCredits[i];
		fprintf(its_fOut, "CDO %u\t\t", i);
		for (j=0; j<NbCredits; j++)
		{
			if ((*its_CDO_Square_Relevant_Credit_Includes)(j,i))
				Flag	=	1;
			else
				Flag	=	0;
			fprintf(its_fOut, "%u\t", Flag);
		}
		fprintf(its_fOut, "\n");
	}

	fprintf(its_fOut, "---------------- For Each Underlying CDO: Credits Included Losses -----------------------------\n");
	for (i=0; i<its_NbRelevant_CDO_Underlyings; i++)
	{
		NbCredits	=	its_CDO_Square_Relevant_CDO_NbCredits[i];
		fprintf(its_fOut, "CDO %u\t\t", i);
		for (j=0; j<NbCredits; j++)
			fprintf(its_fOut, "%.2lf\t", (*its_CDO_Square_Relevant_Credit_Losses)(j,i));
		fprintf(its_fOut, "\n");
	}

	fclose(its_fOut);

}



// ----------------------------------------------------------------------------
// ComputeCumulativeLossesAndNbDefaults for CDO SQUARE
// ----------------------------------------------------------------------------

void CreditManager::ComputeCumulativeLossesAndNbDefaults_CDO_SQUARE()
{
	int		j;
	int		SumNbDef	=	0;
	double	SumLosses	=	0.0;
	double	TheLossAmount	=	0.0;

	int	The_id;

	// -------------------------------------------------------------
	DoubleVector	TheCumulativeLoss_CDOs;
	DoubleVector	TheCumulativeNbDef_CDOs;

	TheCumulativeLoss_CDOs.resize(its_NbRelevant_CDO_Underlyings);
	TheCumulativeNbDef_CDOs.resize(its_NbRelevant_CDO_Underlyings);

	// reset 0.0
	fill(TheCumulativeLoss_CDOs.begin(), TheCumulativeLoss_CDOs.end(), 0.0);
	fill(TheCumulativeNbDef_CDOs.begin(), TheCumulativeNbDef_CDOs.end(), 0.0);
	// -------------------------------------------------------------
	
	// I want to find all items between Tmin and Tmax
	set<Tau_Item>::iterator iter;

	iter = itsSortedDefaultTimes.begin();

	// ----------------------------------------------------------------------------
	// THE LOSS FOR THE CDO SQUARE
	// ----------------------------------------------------------------------------
	double	Current_CDO_Cumul_Loss;
	double	Current_CDO_Loss;
	double	CDO_Square_Loss;
	double	Prev_CDO_Square_Loss;
	
	double	Current_CDO_Attachment_Notio;
	double	Current_CDO_Detachment_Notio;

	double	TheLossAmountUsed;
	int		TheNbDefUsed;

	double	TheRescalingFactor;

	// ----------------------------------------------------------------------------
	Prev_CDO_Square_Loss	=	0.0;

	while (iter != itsSortedDefaultTimes.end())
	{
		iter->cumul_loss	=	SumLosses;		// excluded.

		The_id			=	iter->id;
		TheLossAmount	=	its_LossesAmount[The_id];
		
		// standard CDO
		SumLosses		+=	TheLossAmount;
		
		iter->cumul_nbdef	=	SumNbDef;		// excluded.
		SumNbDef++;

		iter->loss_at_that_time	=	TheLossAmount;

		// CDO Square		
		iter->cumul_loss_CDOs.resize(its_NbRelevant_CDO_Underlyings);
		iter->cumul_nbdef_CDOs.resize(its_NbRelevant_CDO_Underlyings);
		iter->tranche_losses_CDOs.resize(its_NbRelevant_CDO_Underlyings);

		CDO_Square_Loss	=	0.0;

		for  (j=0; j<its_NbRelevant_CDO_Underlyings; j++)
		{
			// improve the handling of 
			if ((*its_CDO_Square_Relevant_Credit_Includes)(The_id, j))
			{
				TheLossAmountUsed	=	TheLossAmount;
				TheNbDefUsed		=	1;
			}
			else
			{
				TheLossAmountUsed	=	0.0;
				TheNbDefUsed		=	0;
			}

			// Portfolio Losses --> EXCLUDED???
			Current_CDO_Cumul_Loss	=	TheCumulativeLoss_CDOs[j];
			
			Current_CDO_Cumul_Loss		+=	TheLossAmountUsed;
			TheCumulativeNbDef_CDOs[j]	+=	TheNbDefUsed;

			TheCumulativeLoss_CDOs[j]	=	Current_CDO_Cumul_Loss;

			// Portfolios' Losses --> INCLUDED
			(iter->cumul_loss_CDOs)[j]	=	Current_CDO_Cumul_Loss;
			(iter->cumul_nbdef_CDOs)[j]	=	TheCumulativeNbDef_CDOs[j];
			
			// Tranches' Losses --> INCLUDED
			// required for CDO^2 Losses computation at Default Dates
			Current_CDO_Attachment_Notio	=	its_CDO_Square_Relevant_CDO_Attachments_Notio[j];
			Current_CDO_Detachment_Notio	=	its_CDO_Square_Relevant_CDO_Detachments_Notio[j];
			Current_CDO_Loss	=	FMAX(Current_CDO_Cumul_Loss - Current_CDO_Attachment_Notio, 0.0) - FMAX(Current_CDO_Cumul_Loss - Current_CDO_Detachment_Notio, 0.0);

			(iter->tranche_losses_CDOs)[j]	=	Current_CDO_Loss;

			// Takes into account Rescaling Factors
			TheRescalingFactor	=	its_CDO_Square_Relevant_CDO_Rescaling_Factor[j];
			CDO_Square_Loss		+=	Current_CDO_Loss * TheRescalingFactor;
		}

		// Keep CDO Square Loss at Current Default Time
		iter->cumul_loss_square					=	CDO_Square_Loss;
		iter->prev_cumul_loss_square			=	Prev_CDO_Square_Loss;
		iter->equiv_current_loss_square			=	CDO_Square_Loss - Prev_CDO_Square_Loss;		
		
		Prev_CDO_Square_Loss	=	CDO_Square_Loss;

		// NEXT ONE
		iter++;

	}

	// ----------------------------------------------------------------------------
	// DISPLAY
	// ----------------------------------------------------------------------------
/*	
	// ----------------------------------------------------------------------------
	if ((its_fOut = fopen("c:\\test\\mc_ComputeCumulativeLossesAndNbDefaults_CDO_SQUARE.txt", "a+")) == NULL) return;
	fprintf(its_fOut, " ----------------- itsSortedDefaultTimes\t\tSimul\t%u ----------------- \n", its_SimulId);
	// ----------------------------------------------------------------------------
	for (iter = itsSortedDefaultTimes.begin();
				iter != itsSortedDefaultTimes.end();
				++iter)
	{
		fprintf(its_fOut, "Id:\t%u\t\ttau:\t%.2lf\t\tCumul loss:\t%.2lf\t\tCumul NbDef:\t%u\t\tCDO Square Cumul Loss:\t%.2lf\n", iter->id, iter->tau, iter->cumul_loss, iter->cumul_nbdef, iter->cumul_loss_square);
		for (j=0; j<its_NbRelevant_CDO_Underlyings; j++)
		{
			fprintf(its_fOut, "CDO %u:\t%.2lf\t\t%u\t\t", j, (iter->cumul_loss_CDOs)[j], (iter->cumul_nbdef_CDOs)[j]);
		}
		fprintf(its_fOut, "Tranche Losses\n");
		for (j=0; j<its_NbRelevant_CDO_Underlyings; j++)
		{
			fprintf(its_fOut, "CDO %u:\t%.2lf\t\t", j, (iter->tranche_losses_CDOs)[j]);
		}
		fprintf(its_fOut, "\n");
	}
	fclose(its_fOut);
	// ----------------------------------------------------------------------------
*/	
}


// ----------------------------------------------------------------------------
// ComputeLossesAmountBeforeCreditWindowsDates for CDO SQUARE
// ----------------------------------------------------------------------------

void CreditManager::ComputeLossesAmountBeforeCreditWindowsDates_CDO_SQUARE()
{
	int	j;

	// to be optimized
	// for each Credit Date, find the cumulative loss
	set<Tau_Item>::iterator iter;
	set<Tau_Item>::iterator iter_DefTime;
	set<Tau_Item>::iterator iter_prev_end;
	set<Tau_Item>::iterator iter_DefTimeBis;

	double	CreditDate, tau_date, PrevCreditDate;
	double	TheCumulativeLoss	=	0.0;
	int		TheCumulativeNbDef	=	0;

	int		TheId;
	double	tmpDate;

	DoubleVector	PrevNbDef(its_NbRelevant_CDO_Underlyings);
	DoubleVector	PrevLoss(its_NbRelevant_CDO_Underlyings);

	bool exit_loop;

	// -------------------------------------------------------------

	iter	= its_SortedCreditObservationDatesAndLosses.begin();
	while (iter != its_SortedCreditObservationDatesAndLosses.end())
	{
		// current Credit Date
		CreditDate	=	iter->tau;

		// RESET!!!
		iter->loss_at_that_time	=	0.0;

		// RESET for CDO SQUARE
		iter->cumul_loss_CDOs.resize(its_NbRelevant_CDO_Underlyings);
		iter->cumul_nbdef_CDOs.resize(its_NbRelevant_CDO_Underlyings);
		iter->tranche_losses_CDOs.resize(its_NbRelevant_CDO_Underlyings);

		// Credit Id
		iter_DefTime	=	itsSortedDefaultTimes.upper_bound(*iter);
		
		// test
		if (iter_DefTime !=	itsSortedDefaultTimes.end())
		{
			TheId		=	iter_DefTime->id;
			tau_date	=	iter_DefTime->tau;

			// CHECK_EQUAL ?
			if ((iter_DefTimeBis	=	itsSortedDefaultTimes.find(*iter)) != itsSortedDefaultTimes.end()) // tau_date == CreditDate
			{
				tmpDate	=	iter_DefTimeBis->tau;
				iter->loss_at_that_time	=	iter_DefTimeBis->loss_at_that_time;
			}
			
			// Cumulative Loss + Current Loss
			TheCumulativeLoss	=	iter_DefTime->cumul_loss;
			TheCumulativeNbDef	=	iter_DefTime->cumul_nbdef;

			// Set the Cumulative Loss
			iter->cumul_loss	=	TheCumulativeLoss;
			iter->cumul_nbdef	=	TheCumulativeNbDef;	// just posterior

			// -------------------------------------------------------------------------
			// CDO SQUARE
			if (iter_DefTime == itsSortedDefaultTimes.begin())
			{
				// Except if Instantaneous or imposed defaults
				TheCumulativeLoss	=	0.0;
				TheCumulativeNbDef	=	0;

				fill(PrevNbDef.begin(), PrevNbDef.end(), 0.0);
				fill(PrevLoss.begin(), PrevLoss.end(), 0.0);

			}
			else
			{
				iter_prev_end	=	iter_DefTime;
				--iter_prev_end;
				PrevCreditDate		=	iter_prev_end->tau;
				
				TheCumulativeLoss	=	TheCumulativeLoss;		//iter_prev_end->cumul_loss;
				TheCumulativeNbDef	=	TheCumulativeNbDef;		//iter_prev_end->cumul_nbdef;

				PrevNbDef	=	iter_prev_end->cumul_nbdef_CDOs;
				PrevLoss	=	iter_prev_end->cumul_loss_CDOs;
			}

			// MAY DO BETTER
			exit_loop	=	false;

			while ((iter != its_SortedCreditObservationDatesAndLosses.end()) && (! exit_loop))
			{
				// current Credit Date
				CreditDate	=	iter->tau;

				if (CreditDate >= tau_date)
					exit_loop	=	true;
				else
				{
					iter->cumul_loss	=	TheCumulativeLoss;
					iter->cumul_nbdef	=	TheCumulativeNbDef;		// to be checked!

					// -------------------------------------------------------------------------
					// CDO SQUARE
	/*				for (j=0; j<its_NbRelevant_CDO_Underlyings; j++)
					{
						// do I have to test? otherwise, fill!!!
	//					if ((*its_CDO_Square_Relevant_Credit_Includes)(TheId, j))
	//					{
							(iter->cumul_loss_CDOs)[j]	=	TheCumulativeLoss;
							(iter->cumul_nbdef_CDOs)[j]	=	TheCumulativeNbDef;
	//					}				
					}
	*/
					iter->cumul_nbdef_CDOs	=	PrevNbDef;
					iter->cumul_loss_CDOs	=	PrevLoss;

					// -------------------------------------------------------------------------

					// NEXT ONE
					++iter;
				}
			}
/*
			for (j=0; j<its_NbRelevant_CDO_Underlyings; j++)
			{
				// do I have to test?
				if ((*its_CDO_Square_Relevant_Credit_Includes)(TheId, j))
				{
					(iter->cumul_loss_CDOs)[j]	=	(iter_DefTime->cumul_loss_CDOs)[j];
					(iter->cumul_nbdef_CDOs)[j]	=	(iter_DefTime->cumul_loss_CDOs)[j];
				}				
			}
			// -------------------------------------------------------------------------
*/			
			// NEXT ONE
//			++iter;

		}
		else
		{
			// I am already out of the range, so that for all posterior Credit Dates (the current included)
			// the Loss Amount, is given by the last computed CumulativeLoss
			if (itsSortedDefaultTimes.size() > 0)
			{
				iter_prev_end	=	itsSortedDefaultTimes.end();
				--iter_prev_end;
				PrevCreditDate		=	iter_prev_end->tau;
				
				TheCumulativeLoss	=	iter_prev_end->cumul_loss + its_LossesAmount[iter_prev_end->id];
				TheCumulativeNbDef	=	iter_prev_end->cumul_nbdef + 1;

				PrevNbDef	=	iter_prev_end->cumul_nbdef_CDOs;
				PrevLoss	=	iter_prev_end->cumul_loss_CDOs;
			}
			else
			{
				TheCumulativeLoss	=	0.0;
				TheCumulativeNbDef	=	0;

				fill(PrevNbDef.begin(), PrevNbDef.end(), 0.0);
				fill(PrevLoss.begin(), PrevLoss.end(), 0.0);
			}

			while (iter != its_SortedCreditObservationDatesAndLosses.end())
			{
				// current Credit Date
				CreditDate	=	iter->tau;

				iter->cumul_loss	=	TheCumulativeLoss;
				iter->cumul_nbdef	=	TheCumulativeNbDef;		// to be checked!

				// -------------------------------------------------------------------------
				// CDO SQUARE
/*				for (j=0; j<its_NbRelevant_CDO_Underlyings; j++)
				{
					// do I have to test? otherwise, fill!!!
//					if ((*its_CDO_Square_Relevant_Credit_Includes)(TheId, j))
//					{
						(iter->cumul_loss_CDOs)[j]	=	TheCumulativeLoss;
						(iter->cumul_nbdef_CDOs)[j]	=	TheCumulativeNbDef;
//					}				
				}
*/
				iter->cumul_nbdef_CDOs	=	PrevNbDef;
				iter->cumul_loss_CDOs	=	PrevLoss;

				// -------------------------------------------------------------------------

				// NEXT ONE
				++iter;
			}
		}
	}

	// --------------------------------------------------------------------------------------
	// Not optimum but more readable
	// --------------------------------------------------------------------------------------
	double	CDO_Square_Loss;
	double	Current_CDO_Loss;
	double	Current_CDO_Cumul_Loss;
	double	Current_CDO_Attachment_Notio;
	double	Current_CDO_Detachment_Notio;
	double	TheRescalingFactor;

	iter	= its_SortedCreditObservationDatesAndLosses.begin();

	while (iter != its_SortedCreditObservationDatesAndLosses.end())
	{
		// CDO Square
		CDO_Square_Loss	=	0.0;

		for  (j=0; j<its_NbRelevant_CDO_Underlyings; j++)
		{
			Current_CDO_Cumul_Loss	=	(iter->cumul_loss_CDOs)[j];

			Current_CDO_Attachment_Notio	=	its_CDO_Square_Relevant_CDO_Attachments_Notio[j];
			Current_CDO_Detachment_Notio	=	its_CDO_Square_Relevant_CDO_Detachments_Notio[j];

			Current_CDO_Loss	=	FMAX(Current_CDO_Cumul_Loss - Current_CDO_Attachment_Notio, 0.0) - FMAX(Current_CDO_Cumul_Loss - Current_CDO_Detachment_Notio, 0.0);

			(iter->tranche_losses_CDOs)[j]	=	Current_CDO_Loss;

			// Takes into account Rescaling Factors
			TheRescalingFactor	=	its_CDO_Square_Relevant_CDO_Rescaling_Factor[j];
			CDO_Square_Loss		+=	Current_CDO_Loss * TheRescalingFactor;
		}

		// Keep CDO Square Loss at Current Default Time
		iter->cumul_loss_square	=	CDO_Square_Loss;
		
		// NEXT ONE
		iter++;
	}
	// --------------------------------------------------------------------------------------

/*
	//	set<Tau_Item>::iterator iter;
	if (its_temp_index <= 0)		// either Price or First hedge Credit
	{
		// ----------------------------------------------------------------------------
		if ((its_fOut = fopen("c:\\test\\mc_ComputeLossesAmountBeforeCreditWindowsDates_CDO_SQUARE.txt", "a+")) == NULL) return;
		fprintf(its_fOut, " ----------------- ComputeLossesAmountBeforeCreditWindowsDates_CDO_SQUARE:\t\t%u ----------------- \n", its_SimulId);
		// ----------------------------------------------------------------------------
		for (iter = its_SortedCreditObservationDatesAndLosses.begin();
					iter != its_SortedCreditObservationDatesAndLosses.end();
					++iter)
		{
			fprintf(its_fOut, "Id:\t%u\t\tCredit Date:\t%.2lf\t\tCumul loss:\t%.2lf\t\tCumul NbDef:\t%u\t\tCumul Loss Square:\t%.2lf\n", iter->id, iter->tau, iter->cumul_loss, iter->cumul_nbdef, iter->cumul_loss_square);
			
			for (j=0; j<its_NbRelevant_CDO_Underlyings; j++)
			{
				fprintf(its_fOut, "CDO %u:\t%.2lf\t\t%u\t\t", j, (iter->cumul_loss_CDOs)[j], (iter->cumul_nbdef_CDOs)[j]);
			}
			fprintf(its_fOut, "\n");
			fprintf(its_fOut, "Tranche Losses\n");
			for (j=0; j<its_NbRelevant_CDO_Underlyings; j++)
			{
				fprintf(its_fOut, "CDO %u:\t%.2lf\t\t", j, (iter->tranche_losses_CDOs)[j]);
			}
			fprintf(its_fOut, "\n");
		}
		fclose(its_fOut);
	}
*/
}


// ----------------------------------------------------------------------------
// GetHowMuchLossBeforeACreditObservation for CDO SQUARE
// ----------------------------------------------------------------------------

void CreditManager::GetHowMuchLossBeforeACreditObservation_CDO_SQUARE(int LowOrUpFlag, double& TheLoss)
{
	// hypothesis: TheDate must be in 'its_SortedCreditObservationDatesAndLosses'
	// uses the current id: its_CurrentIndex
	// otherwise:	implementation of a search in its_Default

	set<Tau_Item>::iterator	iter;

//	double	TheLossAtThatTime;

	// LowOrUpFlag = 0 --> Low
	if (LowOrUpFlag == 0)
	{
		iter	=	its_CreditObervationDatesSortedLowsIds[its_CurrentIndex];
//		TheLossAtThatTime	=	iter->loss_at_that_time;
	}
	else
	{
		iter	=	its_CreditObervationDatesSortedUpsIds[its_CurrentIndex];
//		TheLoss	=	iter->cumul_loss;
//		TheLossAtThatTime	=	0.0;
	}

	// -------------------------------------------------------------------------
	// CDO SQUARE
	TheLoss	=	iter->cumul_loss_square; // - TheLossAtThatTime;
	// -------------------------------------------------------------------------
}


// USELESS ?

void CreditManager::GetHowMuchLossBeforeACreditObservation_CDO_SQUARE(int LowOrUpFlag, DoubleVector& TheLosses)
{
	int j;
	// hypothesis: TheDate must be in 'its_SortedCreditObservationDatesAndLosses'
	// uses the current id: its_CurrentIndex
	// otherwise:	implementation of a search in its_Default

	// TheLosses must have the good size.

	set<Tau_Item>::iterator	iter;

//	double	TheLossAtThatTime;

	// LowOrUpFlag = 0 --> Low
	if (LowOrUpFlag == 0)
	{
		iter	=	its_CreditObervationDatesSortedLowsIds[its_CurrentIndex];
//		TheLossAtThatTime	=	iter->loss_at_that_time;
	}
	else
	{
		iter	=	its_CreditObervationDatesSortedUpsIds[its_CurrentIndex];
//		TheLoss	=	iter->cumul_loss;
//		TheLossAtThatTime	=	0.0;
	}

	// -------------------------------------------------------------------------
	// CDO SQUARE
	for (j=0; j<its_NbRelevant_CDO_Underlyings; j++)
		TheLosses[j]	=	(iter->cumul_loss_CDOs)[j]; // - TheLossAtThatTime;
	// -------------------------------------------------------------------------
}


void CreditManager::GetHowMuchLossBeforeADate_CDO_SQUARE(RelativeDate TheDate, double& TheLoss)
{
	if (itsSortedDefaultTimes.size() == 0)
	{
		TheLoss	=	0.0;
		return;
	}

	set<Tau_Item>::iterator	iter			=	itsSortedDefaultTimes.begin();
	set<Tau_Item>::iterator	iter_prev_end	=	itsSortedDefaultTimes.end();	
	--iter_prev_end;

	// the Item to be searched
	Tau_Item	SearchItem(0, (double) TheDate, 0.0);

	// iter points at the first item which is above or equal TheDate
	iter	=	itsSortedDefaultTimes.lower_bound(SearchItem);
	
	if (iter == itsSortedDefaultTimes.end())
	{		
		// -------------------------------------------------------------------------
		// CDO SQUARE
		TheLoss	=	iter_prev_end->cumul_loss_square;
//		TheLoss	+=	its_LossesAmount[iter_prev_end->id];
		// -------------------------------------------------------------------------
	}
	else
	{
		// -------------------------------------------------------------------------
		// CDO SQUARE
		TheLoss	=	iter->cumul_loss_square;
		// -------------------------------------------------------------------------
	}

}

// USELESS ?

void CreditManager::GetHowMuchLossBeforeADate_CDO_SQUARE(RelativeDate TheDate, DoubleVector& TheLosses)
{
	// TheLosses must have the good size.
	int j;
	double	TheLoss;

	if (itsSortedDefaultTimes.size() == 0)
	{
		fill(TheLosses.begin(), TheLosses.end(), 0.0);
		return;
	}

	set<Tau_Item>::iterator	iter			=	itsSortedDefaultTimes.begin();
	set<Tau_Item>::iterator	iter_prev_end	=	itsSortedDefaultTimes.end();	
	--iter_prev_end;

	// the Item to be searched
	Tau_Item	SearchItem(0, (double) TheDate, 0.0);

	// iter points at the first item which is above or equal TheDate
	iter	=	itsSortedDefaultTimes.lower_bound(SearchItem);
	
	if (iter == itsSortedDefaultTimes.end())
	{		
		// -------------------------------------------------------------------------
		// CDO SQUARE
		for (j=0; j<its_NbRelevant_CDO_Underlyings; j++)
		{
			TheLoss	=	(iter_prev_end->cumul_loss_CDOs)[j];
//			TheLoss	+=	its_LossesAmount[iter_prev_end->id];

			TheLosses[j]	=	TheLoss;
		}
		// -------------------------------------------------------------------------
	}
	else
	{
		// -------------------------------------------------------------------------
		// CDO SQUARE
		for (j=0; j<its_NbRelevant_CDO_Underlyings; j++)
		{
			TheLoss	=	(iter->cumul_loss_CDOs)[j];
			TheLosses[j]	=	TheLoss;
		}
		// -------------------------------------------------------------------------
	}

}


void CreditManager::GetHowManyDefaultsBeforeACreditObservation_CDO_SQUARE(int LowOrUpFlag, int& TheNbDef, set<Tau_Item>::iterator&	iter)
{
	// hypothesis: TheDate must be in 'its_SortedCreditObservationDatesAndLosses'
	// uses the current id: its_CurrentIndex
	// otherwise:	implementation of a search in its_Default

	// LowOrUpFlag = 0 --> Low
	if (LowOrUpFlag == 0)
		iter	=	its_CreditObervationDatesSortedLowsIds[its_CurrentIndex];
	else
		iter	=	its_CreditObervationDatesSortedUpsIds[its_CurrentIndex];


	// -------------------------------------------------------------------------
	// CDO SQUARE
	TheNbDef	=	iter->cumul_nbdef_square;
	// -------------------------------------------------------------------------

}

// USELESS ?

void CreditManager::GetHowManyDefaultsBeforeACreditObservation_CDO_SQUARE(int LowOrUpFlag, IntVector& TheNbDefs, set<Tau_Item>::iterator&	iter)
{
	int j;
	// TheNbDefs must have the good size.

	// hypothesis: TheDate must be in 'its_SortedCreditObservationDatesAndLosses'
	// uses the current id: its_CurrentIndex
	// otherwise:	implementation of a search in its_Default

	// LowOrUpFlag = 0 --> Low
	if (LowOrUpFlag == 0)
	{
		iter	=	its_CreditObervationDatesSortedLowsIds[its_CurrentIndex];
	}
	else
	{
		iter	=	its_CreditObervationDatesSortedUpsIds[its_CurrentIndex];
	}
	// -------------------------------------------------------------------------
	// CDO SQUARE
	for (j=0; j<its_NbRelevant_CDO_Underlyings; j++)
		TheNbDefs[j]	=	(iter->cumul_nbdef_CDOs)[j];
	// -------------------------------------------------------------------------

}


void CreditManager::GetHowManyDefaultsBeforeADate_CDO_SQUARE(RelativeDate TheDate, int& TheNbDef, set<Tau_Item>::iterator&	iter)
{
	// TheNbDefs must have the good size.

	iter	=	itsSortedDefaultTimes.begin();

	if (itsSortedDefaultTimes.size() == 0)
	{
		TheNbDef	=	0;
		return;
	}

	set<Tau_Item>::iterator	iter_prev_end	=	itsSortedDefaultTimes.end();
	--iter_prev_end;

	// the Item to be searched
	Tau_Item	SearchItem(0, (double) TheDate, 0.0);

	// iter points at the first item which is above or equal TheDate
	iter	=	itsSortedDefaultTimes.lower_bound(SearchItem);
	
	if (iter == itsSortedDefaultTimes.end())
	{				
		// -------------------------------------------------------------------------
		// CDO SQUARE
		TheNbDef	=	iter_prev_end->cumul_nbdef_square;
//		TheNbDef++;
		// -------------------------------------------------------------------------
	}
	else
		// -------------------------------------------------------------------------
		// CDO SQUARE
		TheNbDef	=	iter->cumul_nbdef_square;
		// -------------------------------------------------------------------------

}


// USELESS ?

void CreditManager::GetHowManyDefaultsBeforeADate_CDO_SQUARE(RelativeDate TheDate, IntVector& TheNbDefs, set<Tau_Item>::iterator&	iter)
{
	int		j, TheNbDef;

	// TheNbDefs must have the good size.

	iter	=	itsSortedDefaultTimes.begin();

	if (itsSortedDefaultTimes.size() == 0)
	{
		TheNbDef	=	0;
		return;
	}

	set<Tau_Item>::iterator	iter_prev_end	=	itsSortedDefaultTimes.end();
	--iter_prev_end;

	// the Item to be searched
	Tau_Item	SearchItem(0, (double) TheDate, 0.0);

	// iter points at the first item which is above or equal TheDate
	iter	=	itsSortedDefaultTimes.lower_bound(SearchItem);
	
	if (iter == itsSortedDefaultTimes.end())
	{				
		// -------------------------------------------------------------------------
		// CDO SQUARE
		for (j=0; j<its_NbRelevant_CDO_Underlyings; j++)
		{
			TheNbDef	=	(iter_prev_end->cumul_nbdef_CDOs)[j];
//			TheNbDef++;

			TheNbDefs[j]	=	TheNbDef;
		}
		// -------------------------------------------------------------------------
	}
	else
		// -------------------------------------------------------------------------
		// CDO SQUARE
		for (j=0; j<its_NbRelevant_CDO_Underlyings; j++)
		{
			TheNbDef	=	(iter->cumul_nbdef_CDOs)[j];
			TheNbDefs[j]	=	TheNbDef;
		}
		// -------------------------------------------------------------------------

}


// ************************************************************************************************************
// ************************************************************************************************************

// ************************************************************************************************************
// ************************************************************************************************************

// ----------------------------------------------------------------------------
// Default Leg Pricing for a Simulation NbDef
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Premium Leg Pricing for a Simulation LOSSES _CDO_SQUARE
// ----------------------------------------------------------------------------

/*
void CreditManager::PricePremiumLegForThisSimulation_Losses_CDO_SQUARE(DoubleVector& Outputs)
{
}
*/

// ----------------------------------------------------------------------------
// Premium Leg Pricing for a Simulation NBDEF for CDO SQUARE
// ----------------------------------------------------------------------------
/*
void CreditManager::PricePremiumLegForThisSimulation_NbDef_CDO_SQUARE(DoubleVector& Outputs)
{
}
*/
