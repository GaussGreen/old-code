
#include "ARMKernel\glob\firsttoinc.h"

#include "ICMKernel\cair\credit_manager.h"
#include "ICMKernel\glob\icm_maths.h"

#include "ICMKernel\util\icm_gengauss1F.h"
#include "ICMKernel\util\icm_gen_student1f.h"
#include "ICMKernel\util\icm_gengauss_1F_factorloading_2points.h"

#include "ICMKernel\crv\ICM_Constant_Piecewise.h"

#include "ICMKernel\util\icm_integrator.h"

#include <nags.h>		//	s15abc
#include <nagg01.h>		//	g01fac

#include <algorithm>
#include <functional>
#include <numeric>

#include "ICMKernel\util\icm_cholesky.h"
#include "ICMKernel\util\icm_pgcd.h"

#include "ICMKernel\cair\ICM_Probability_Density_Standard_Gaussian.h"
#include "ICMKernel\cair\ICM_Probability_Density_Student.h"

#include "ICMKernel\glob\icm_smile_correlation.h"
#include "ICMKernel\glob\icm_corrmatrix.h"

#define		NB_ITER_NEWTON	20

//------------------------------------------------
// Initialization of datas members
//------------------------------------------------

using namespace std ; 

void CreditManager::Init(void)
{
	ICM_Pricer_Security::Init();

	SetName(ICM_CREDIT_MANAGER);

	itsValDateAsDouble	=	0.0;
	itsValDate	=	"01/01/2001";
	itsMaxCDSDate = "01/01/2001";
	itsZCValuesForCalibration = NULL;
	itsElapsed_Time	=	0.0;

	// ZC curve
	its_ZeroCurve = NULL;
	its_ImposedIRCurveFlag	=	false;

	its_IR_Lags.clear();
	its_IR_ZC_Yields.clear();

	// Flags
	its_RollDateFlag		=	true;
	its_KeepCalibrationFlag	=	false;
	its_Bump_Choice		=	CHB_NO;
	its_BumpSpread_Type	=	BT_ADD;

	// Values
	its_BumpSpread		=	0.0;
	its_BumpRecovery	=	0.0;
	its_BumpCorrelation	=	0.0;

	its_CreditDataSpreads	=	NULL;

	its_Categories.clear();
	its_Currencies.clear();
	its_Accrueds.clear();
	its_Recoveries.clear();
	its_Notionals.clear();
	its_Input_Losses.clear();
	its_DefaultDates.clear();
	its_AmortizationDates.clear();

	// Maturities
	its_Maturities.clear();

	memset(its_Maturities_AsChar,'\0', sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS); 	
//	its_Maturities_AsChar	=	NULL;

	its_ModelMultiCurves	=	NULL;

	// Useful
	its_LastPricingDate	=	0.0;

	// Pricer Data
	its_NbSimul			=	0;
	its_CopulaType		=	CCT_GAUSSIAN;
	its_CorrelationType	=	CT_FLAT;
	its_CreditModelType	=	CMT_MONTECARLO;
	its_FreedomDegree	=	4;
	its_CorrelationValue	=	0.5;

	its_NIntegration_1F	=	51;
	its_N_FFT			=	512;

	// Default Leg Data
	its_DL_CreditWindowLow		=	0.0;
	its_DL_CreditWindowUp		=	0.0;
	its_DL_LossMin				=	0.0;
	its_DL_LossMax				=	0.0;
	its_DL_NbDefMin				=	0;
	its_DL_NbDefMax				=	0;
	its_DL_PaymentType			=	CEP_ATDEFAULTDATE;
	its_DL_PaymentDate			=	0.0;
	its_DL_PaymentLag			=	0;

	// Premium Leg Data
	its_PL_CreditWindowLows.clear();
	its_PL_CreditWindowUps.clear();
	its_PL_StartDates.clear();
	its_PL_EndDates.clear();
	its_PL_PaymentDates.clear();

	its_PL_LossMins.clear();
	its_PL_LossMaxs.clear();
	its_PL_NbDefMins.clear();
	its_PL_NbDefMaxs.clear();
	its_PL_Ratios.clear();
	its_PL_Notios.clear();
	its_PL_Spreads.clear();

	its_PL_CreditFlags.clear();
	its_PL_CreditSpreadCaps.clear();
	its_PL_Redemptions.clear();

	its_PL_CreditWindowLow	=	0.0;
	its_PL_CreditWindowUp	=	0.0;
	its_PL_StartDate		=	0.0;
	its_PL_EndDate			=	0.0;
	its_PL_PaymentDate		=	0.0;

	its_PL_LossMin			=	0.0;
	its_PL_LossMax			=	0.0;
	its_PL_Ratio			=	0.0;
	its_PL_Notio			=	0.0;
	its_PL_Spread			=	0.0;

	its_PL_CreditFlag		=	CPLT_NOTIONAL;
	its_PL_CreditSpreadCap	=	0.0;
	
	its_PL_PaymentType			=	CEP_ATDEFAULTDATE;

	// internal data
	its_PL_NbFlows		=	0;
	its_CurrentIndex	=	0;
	its_SimulId			=	0;
	its_iFlow			=	0;

	its_PL_CreditWindowLow	=	0.0;
	its_PL_CreditWindowUp	=	0.0;
	its_PL_StartDate	=	0.0;
	its_PL_EndDate	=	0.0;
	its_PL_PaymentDate	=	0.0;

	its_PL_LossMin	=	0.0;
	its_PL_LossMax	=	0.0;
	its_PL_NbDefMin	=	0.0;
	its_PL_NbDefMax	=	0.0;
	its_PL_Ratio	=	0.0;
	its_PL_Notio	=	0.0;
	its_PL_Spread	=	0.0;

	its_PL_CreditFlag		=	CPLT_OUTSTANDING;
	its_PL_CreditSpreadCap	=	0.0;
	its_PL_Redemption		=	0.0;

	its_Cumulative_Coupon	=	0.0;
	its_Last_Coupon			=	0.0;
	its_Last_FlowForATMMargin	=	0.0;;
	its_Toggle_Redemption	=	false;
	its_Final_Redemption	=	false;

	its_Useful_MaxNbDef		=	0;
	its_Useful_MaxLoss		=	0.0;

	its_DF.clear();

	itsSortedDefaultTimes.clear();
	itsBarriers_Standard.clear();
	itsBarriers_Shifted.clear();

	itsDef_Prob.clear();

	itsGenerator = NULL;

	its_CreditsLabels.clear();

	its_PricingLegsType			=	CBN_PREMIUMMINUSDEFAULT;
	its_CreditObservationType	=	CO_LOSSES;
	its_CreditPremiumLegAccrued	=	CAP_NON_PRORATA;
	its_ATMDataFlag				=	CPLADT_PURESPREAD;
	its_NPV_Type				=	CNPV_STANDARD;

	DefNPVFlag	=	0.0;
	PremNPVFlag	=	0.0;

	its_DF_Flag	=	true;
	its_MaxDate	=	0.0;

	its_MaxCreditDate	=	0.0;

	its_PL_CreditWindowUps_DF.clear();
	its_PL_PaymentDates_DF.clear();

	its_PL_Sorted_PaymentDates.clear();
	its_PL_Sorted_CreditWindowUps.clear();
	its_PL_Sorted_CreditWindowUps_DF.clear();
	its_PL_Sorted_PaymentDates_DF.clear();

	its_SortedCreditObservationDates.clear();
	its_CreditObervationDatesSortedLowsIds.clear();
	its_CreditObervationDatesSortedUpsIds.clear();

	its_LossesAmount.clear();
	its_SortedCreditObservationDatesAndLosses.clear();

	its_CreditCalibrator_Lags.clear();
	its_CreditCalibrator_Prob_Default.clear();
	its_CreditCalibrator_Std_Error.clear();

	its_Data_Matrix_NbRows		=	0;
	its_Data_Matrix_NbColumns	=	0;


	its_HedgesRunning	=	false;
	its_CurrentHedgesIndex	=	-1;
	its_AllNPVS.clear();

	its_AllMatrixNPVS	=	NULL;
	its_Shift_Results	=	NULL;
	
	its_AllDefLegs.clear();
	its_AllPremLegs.clear();

	its_DefaultTimes.clear();
	its_DefaultTimesShifted.clear();
	itsSortedDefaultTimes_Keep.clear();


	its_Hedges_Delta_NPV.clear();
	its_Hedges_Hedge_Ratio.clear();
	its_Hedges_Delta_CDS.clear();
	its_Hedges_CDS_ATM_Margin.clear();
	its_Hedges_Carry.clear();

	its_Hedges_MatrixOutput.clear();

	// Correlation
	its_Beta.clear();
	its_Base_Correlation_Strikes.clear();
	its_Base_Correlation_Values.clear();
	its_FL_Alpha.clear();
	its_FL_Beta1.clear();
	its_FL_Beta2.clear();
	its_FL_Beta1_Complement.clear();
	its_FL_Beta2_Complement.clear();

	itsCorrelation			=	NULL;
	its_CorrelationMatrix	=	NULL;
	its_CholeskyMatrix		=	NULL;
	its_ArrayDefaultCrv		=	NULL;
	its_MatrixShiftedDefaultCrv	=	NULL;

	its_CC_Lags_String	=	"";
	its_CC_Prob_String	=	"";
	its_CC_Std_Error_String	=	"";

	its_MC_Variance_Reduction = CMCVR_NONE;

	its_BasketNotional	=	0.0;
	its_TotalLoss		=	0.0;
	its_LossesAmount.clear();
	TheDFAtAFixedDate	=	0.0;

	// Outputs
	its_UpFront_RunningSpread	=	0.0;
	its_UpFront_Premium	=	0.0;

	NPV	=	0.0;
	DefaultLegPV	=	0.0;
	PremiumLegPV	=	0.0;
	ATMPremiumLeg	=	0.0;
	ATMPremiumLegWithoutNotio	=	0.0;
	ATMSpread	=	0.0;
	ATMUpFront	=	0.0;
	ATMRunningSpread	=	0.0;
	NPVSquare	=	0.0;
	StdError	=	0.0;
	Central_NPV	=	0.0;

	its_Cash	=	0.0;
	its_Accrued	=	0.0;

	// CALIBRATION
	its_CreditCalibratorChoice	=	CCS_SINGLENAME;
	its_CreditCalibrator_NameId	=	0;
	its_CreditCalibrator_NTD	=	1;
	its_CreditCalibrator_Loss	=	0.0;

	its_CreditCalibratorNbMaturities	=	0;

	its_CreditCalibrator_Density_NTD	=	1;
	its_CreditCalibrator_Density_Loss_Down	=	0.0;
	its_CreditCalibrator_Density_Loss_Up	=	0.0;
	its_CreditCalibrator_Density_Loss_Step	=	0.0;

	size_DENSITY_LOSS	=	0;

	its_Data_Matrix_NbRows		=	0;
	its_Data_Matrix_NbColumns	=	0;

	// VARIANCE REDUCTION
	its_MC_Variance_Reduction	=	CMCVR_NONE;
	its_IS_Theta	=	0.0;
	its_IS_Mu		=	0.0;
	its_IS_Loss_Level	=	0.0;
	its_IS_Loss_Maturity	=	0.0;
	IS_Common_exp_twists	=	0.0;
	IS_Common_FactorShift	=	0.0;
	IS_Individual_FactorShift	=	0.0;

	IS_Common_exp_twistss.clear();
	IS_individual_exp_twists.clear();
	IS_individual_exp_twists_Kept.clear();
	TheCond_DefProb_Vect.clear();

	// HEDGES
	its_temp_index	=	0;
	its_HedgesDefault	=	false;

	its_HedgesCDSMaturity	=	"";
	HedgesCDSToMaturityFlag	=	false;

	its_fOut	=	NULL;

	// 1 FACTOR
	its_Time_Step_Prorata_1F	=	30.0;

	// RECURSIVE
	its_Recursive_1F_Loss_Unit_Choice	=	COFRLUC_PGCD;
	its_LossRate.clear();

	its_LossUnit	=	0.0;
	its_Recursive_1F_Loss_Unit_Min	=	0.0;
	its_Recursive_1F_NbLossStep	=	0;

	// 1 FACTOR
	its_ProbCond			=	NULL;
	its_ProbCond_Perturb	=	NULL;
	its_lossdistrib_perturb	=	NULL;
	its_barrier_perturb		=	NULL;
	its_taildistrib_perturb	=	NULL;

	its_pdef_perturb.clear();

	// VARIANCE REDUCTION
	its_Target_Default_Prob	=	0.0;

	// CDO SQUARE
	its_CDO_Type	=	CPT_STD_CDO;

	its_CDO_Square_Credit_Ids		=	NULL;
	its_CDO_Square_Credit_Losses	=	NULL;
	its_CDO_Square_Credit_Notionals	=	NULL;

	its_CDO_Square_CDO_Includes.clear();
	its_CDO_Square_CDO_Maturities.clear();
	its_CDO_Square_CDO_NbCredits.clear();
	its_CDO_Square_CDO_Attachments.clear();
	its_CDO_Square_CDO_Detachments.clear();
	its_CDO_Square_CDO_Notionals.clear();

	its_CDO_Square_Relevant_CDO_Ids.clear();
	its_CDO_Square_Relevant_CDO_MaturitiesInYearFraction.clear();
	its_CDO_Square_Relevant_CDO_NbCredits.clear();
	its_CDO_Square_Relevant_CDO_Attachments_Pct.clear();
	its_CDO_Square_Relevant_CDO_Detachments_Pct.clear();
	its_CDO_Square_Relevant_CDO_Attachments_Notio.clear();
	its_CDO_Square_Relevant_CDO_Detachments_Notio.clear();

	its_CDO_Square_Relevant_CDO_Notionals.clear();
	its_CDO_Square_Relevant_CDO_Rescaling_Factor.clear();

	its_CDO_Square_Relevant_Credit_Ids	=	NULL;
	its_CDO_Square_Relevant_Credit_Losses	=	NULL;
	its_CDO_Square_Relevant_Credit_Includes	=	NULL;

	// LARGE PORTFOLIO: HERMITE EXPANSIONS
	its_LP_Degree	=	0;		// STANDARD GAUSSIAN APPROXIMATION

	its_DefCurve	=	NULL;

	its_CurrentHermiteCoeffs.clear();
	
	its_CommonFactors.clear();
	its_BarrierDerivatives_FastHedge.clear();

	// NIG
	its_NIG_Alpha	=	0.0;
	its_NIG_Beta	=	0.0;
	its_NIG_Mu		=	0.0;
	its_NIG_Delta	=	0.0;
	its_NIG_Rho		=	0.0;

	// STOCHASTIC CORRELATION
	its_SC_Rho1		=	0.0;
	its_SC_Rho2		=	0.0;
	its_SC_Prob		=	0.0;

	its_Used_Beta_SC_Vector.clear();
	its_Used_SQRT_OneMinusBetaSquare_SC_Vector.clear();
	its_SC_Coefficients.clear();

	its_Correlation_Calibrator	=	NULL;

	// ------------------------------------------------------------
	its_Credit_Product	=	NULL;
	its_Credit_Market_Data	=	NULL;

	its_Hedges_To_Compute	=	true;
	// ------------------------------------------------------------

	its_Algo_Type = 1;

	// ------------------------------------------------------------
	// DISPLAY
	its_Display_DefaultTimes_Flag	=	false;
	its_Display_LossesAmount_Before_CreditDates_Flag	=	false;
	its_Display_Sorted_DefaultTimes_Flag	=	false;
	its_Display_CumulativeLossesAndDefaults_Flag	=	false;

	its_Loss_Distrib_Map.clear();

//	Test_Integration_Dim();
}

/*
void CreditManager::BitwiseCopy(const ARM_Object* src)
{
	CreditManager* dc = (CreditManager*) src;

	itsValDate = dc->itsValDate;
	itsMaxCDSDate = dc->itsMaxCDSDate;

	if (dc->itsZCValuesForCalibration)
		itsZCValuesForCalibration = (ARM_Vector*) dc->itsZCValuesForCalibration->Clone();

	itsElapsed_Time = dc->itsElapsed_Time;

	if (dc->its_ZeroCurve)
		its_ZeroCurve = (ARM_ZeroCurve*) dc->its_ZeroCurve->Clone();

	its_ImposedIRCurveFlag = dc->its_ImposedIRCurveFlag;

	its_NbCredits = dc->its_NbCredits;

	if (dc->its_CreditDataSpreads)
		its_CreditDataSpreads = (ICM_QMatrix*) dc->its_CreditDataSpreads->Clone();

}


void CreditManager::Copy(const ARM_Object* src)
{
   ARM_Object::Copy(src);

   BitwiseCopy(src);
}
*/

void CreditManager::Reset()
{
	ICM_Pricer_Security::Reset();

	if (itsZCValuesForCalibration)
		delete itsZCValuesForCalibration;
	itsZCValuesForCalibration	=	NULL;

	if (its_ZeroCurve)
		delete its_ZeroCurve;
	its_ZeroCurve	=	NULL;

	if (its_CreditDataSpreads)
		delete its_CreditDataSpreads;
	its_CreditDataSpreads	=	NULL;

}

void CreditManager::SetZeroCurve(ARM_ZeroCurve* crv) 
{ 
	if (its_ZeroCurve) delete its_ZeroCurve;	
	its_ZeroCurve = (ARM_ZeroCurve*) crv->Clone();
}


// some tests to be added?
double	CreditManager::DiscountPrice(int i)
{
	return (*itsZCValuesForCalibration)[i];
}

/*
// CHANGE DEFAULT CURVE SET...
void	CreditManager::Set_CreditDataMaturitiesAsChar(char** data)
{
	int	i;

	its_Maturities_AsChar = new char*[its_NbMaturities];

	for (i=0;i<its_NbMaturities;i++)
	{
		its_Maturities_AsChar[i] = new char[_SIZE_CREDIT_MATURITIES_];
		strcpy(its_Maturities_AsChar[i], data[i]);
	}
}
*/

void	CreditManager::Set_CreditDataMaturitiesAsChar(char data[ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS])
{
	memcpy(its_Maturities_AsChar, data, sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);
}


/*void	CreditManager::Set_CreditsLabelsAsChar(char** data)
{
	int	i;

	its_CreditsLabelsAsChar = new char*[its_NbCredits];

	for (i=0;i<its_NbCredits;i++)
	{
		its_CreditsLabelsAsChar[i] = new char[_size_creditlabel_];		// 60
		strcpy(its_CreditsLabelsAsChar[i], data[i]);
	}
}

void	CreditManager::Set_CreditsLabelsAsChar(const std::vector<std::string>&data)
{
	int	i;

	its_CreditsLabelsAsChar = new char*[its_NbCredits];

	for (i=0;i<its_NbCredits;i++)
	{
		its_CreditsLabelsAsChar[i] = new char[_size_creditlabel_];		// 60
		strcpy(its_CreditsLabelsAsChar[i], data[i].c_str());
	}
}
*/
void	CreditManager::Set_CreditCalibratorMaturitiesAsChar(int NbMat, char data[ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS])
{
	its_CreditCalibratorNbMaturities	=	NbMat;
	memcpy(its_CreditCalibratorMaturities_AsChar, data, sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);
}


void	CreditManager::Set_CreditDataHedgesCDSMaturity(string data)
{
	its_HedgesCDSMaturity = data;

	string	s1	= "TO_MATURITY";

	if (its_HedgesCDSMaturity == s1)
		SetHedgesCDSToMaturity();
	else
		SetHedgesCDSToSpreadCurve();
}

void CreditManager::CholeskyComputation()
{
	if (its_CholeskyMatrix)
		delete its_CholeskyMatrix;
	its_CholeskyMatrix	=	NULL;

	if (its_CorrelationMatrix == NULL) return;

	its_CholeskyMatrix	=	new ICM_QMatrix<double>(its_CreditsLabels.size(), its_CreditsLabels.size());

	int	i,j;

	const ICM_QMatrix<double>&	TheCorrelMatrix = its_CorrelationMatrix->GetMatrix();

    
    {
		for (i=0;i<its_CreditsLabels.size();i++)
			for (j=0;j<its_CreditsLabels.size();j++)
				its_CholeskyMatrix->SetValue(i,j,TheCorrelMatrix.Getvalue(i,j));

		double tmp=0.;

		cholesky Obj=cholesky();

		Obj.init(*its_CholeskyMatrix,its_CreditsLabels.size());
		Obj.choldc(its_CreditsLabels.size());

		for (i=0;i<its_CreditsLabels.size();i++)
		{
			its_CholeskyMatrix->SetValue(i,i,Obj.m_eigenvalues[i]);

			for (j=i+1;j<its_CreditsLabels.size();j++)
				its_CholeskyMatrix->SetValue(i,j,0.);
		}

		for (i=0;i<its_CreditsLabels.size();i++)
		{
			for (j=0;j<i;j++)
			{
				tmp=Obj.m_matrix(i,j);
				its_CholeskyMatrix->SetValue(i,j,Obj.m_matrix(i,j));
				tmp=its_CholeskyMatrix->Getvalue(i,j);
			}
		}

    }
   
}

//--------------------------------------------
//	Display Method: TEMPORARY
//--------------------------------------------
void CreditManager::Display()
{
    char fOutName[100]	=	"C:\\Credit\\Excel\\display_Credit_Manager.txt";
 

    if ((its_fOut = fopen(fOutName, "w")) == NULL)
		// pb.
		return ;

	fprintf(its_fOut, "\n\n------------------------------------------------------------------\n");
	fprintf(its_fOut, "---------------- Credit Manager -----------------------------\n");
	fprintf(its_fOut, "------------------------------------------------------------------\n\n\n");

	fprintf(its_fOut, "Val Date\t:\t%ld.%ld.%ld\n\n", itsValDate.GetDay(), itsValDate.GetMonth(), itsValDate.GetYear());
	fprintf(its_fOut, "Max CDS Date\t:\t%ld.%ld.%ld\n\n", itsMaxCDSDate.GetDay(), itsMaxCDSDate.GetMonth(), itsMaxCDSDate.GetYear());
	fprintf(its_fOut, "Elapsed Time (in s)\t:\t%.lf\n\n", itsElapsed_Time);

	fprintf(its_fOut, "\n\nIndex\t\tZC Value\n");

	int NbPoints;
	NbPoints	=	itsZCValuesForCalibration->GetSize();

	for (int i = 0; i < NbPoints ; i++)
	{
		fprintf(its_fOut, "\n");	
		fprintf(its_fOut, "\t%lu\t\t%f", i, (*itsZCValuesForCalibration)[i]);
	}		 

	fclose(its_fOut);
 
}


/**********************************************************************************************************
****************	SOLVER TEST	*******************
**********************************************************************************************************/
/*
static Index TotalHitsCount;

extern "C" ReturnCode TheFunction(void* Param, DoubleVector& CalibrationVector, double& Result)
{
	ReturnCode Err;
	Result = 0;

	// From The X vector --> go to the 'ThePoint'
	Err = ((CreditManager*)Param)->ExtractCurrentCalibration(CalibrationVector);
	if (Err != RetOk) return Err;

	Err = ((CreditManager*)Param)->CalibrationFunction(Result);
	if (Err != RetOk) return Err;

	TotalHitsCount++;
	return RetOk;
}

ReturnCode CreditManager :: ExtractCurrentCalibration(DoubleVector& CalibrationVector)
{
	ThePoint	=	CalibrationVector;
	return RetOk;
}


ReturnCode CreditManager :: CalibrationFunction(double& Result)
{
	Index	TheDimCol;
	Index	TheDimRow;
	
	The_Coeffs.GetCount(TheDimRow,TheDimCol);

	double	ModelPrice, TheValue;

	Result	=	0.0;
	int i,j;

	for (i=0; i<TheDimRow; i++)
	{
		ModelPrice	=	0.0;

		for (j=0; j<TheDimCol; j++)
		{
			TheValue	=	0.0;

			switch (j)
			{
			case 0:
				// coeff x2
				TheValue	=	The_Coeffs[i][j]	*	ThePoint[0] * ThePoint[0];
				break;
			case 1:
				// coeff log y
				TheValue	=	The_Coeffs[i][j]	*	log(ThePoint[1]);
				break;
			case 2:
				// coeff exp z
//				TheValue	=	The_Coeffs[i][j]	*	exp(ThePoint[2]);
				break;
			case 3:
				// coeff exp-x
				TheValue	=	The_Coeffs[i][j]	*	exp(-ThePoint[0]);
				break;
			case 4:
				// coeff y
				TheValue	=	The_Coeffs[i][j]	*	ThePoint[1];
				break;

			case 5:
				// coeff 1/z
//				TheValue	=	The_Coeffs[i][j]	*	1.0 / ThePoint[2];
				break;

			}
			ModelPrice	+=	TheValue;

		}
		Result	+=	(ModelPrice - The_Targets[i]) * (ModelPrice - The_Targets[i]);
	}

	return RetOk;
}


/**********************************************************************************************************
****************	SOLVER TEST	*******************
**********************************************************************************************************/
/*
void CreditManager::TestSolver()
{
	QuasiNewton		TheSolver;

	CairString	Info;
	double		Norm;
	Index		NbIter;

	TheSolver.SetNbMaxIter(200);
	TheSolver.SetTarget(1e-6);

	DoubleVector Steps;

	Index	TheDimCol	=	6;
	Index	TheDimRow	=	2;

	// Values --> X vector
	The_Coeffs.Resize(TheDimRow,TheDimCol);
	The_Coeffs[0][0]	=	1.0;
	The_Coeffs[0][1]	=	1.0;
	The_Coeffs[0][2]	=	0.0;
	The_Coeffs[0][3]	=	0.0;
	The_Coeffs[0][4]	=	0.0;
	The_Coeffs[0][5]	=	0.0;
	The_Coeffs[1][0]	=	0.0;
	The_Coeffs[1][1]	=	0.0;
	The_Coeffs[1][2]	=	1.0;
	The_Coeffs[1][3]	=	1.0;
	The_Coeffs[1][4]	=	1.0;
	The_Coeffs[1][4]	=	1.0;

	The_Targets.Resize(TheDimRow);
	The_Targets[0]	=	16.0;
	The_Targets[1]	=	10.0;

	ThePoint.Resize(TheDimRow);

	// Init
	ThePoint[0]		=	1.0;
	ThePoint[1]		=	1.0;

	Steps.Resize(TheDimRow);

	for (int i=0; i<TheDimRow;i++)
		Steps[i]	=	ThePoint[i]	*	0.01;

	ReturnCode	Err;
	double	Result;

	Err = TheSolver.Solve(this, ThePoint, Steps, TheFunction, Info, Result, Norm, NbIter);
	Result = sqrt(Result);

}

*/
//--------------------------------------------
//	MARKET DATA PARAMETERS
//--------------------------------------------

void	CreditManager::SetMarketDataParameters(ICM_Matrix<ARM_Vector>* parameters)
{
	if (parameters == NULL) return;

	ARM_Vector* THE_VECT=NULL;

	// ------------------------------------------------------
	// DOUBLE TYPE: ValDate from Excel
	THE_VECT	= parameters->GetColVect("EXCEL_VALDATE");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  no Excel Date vector found");

	itsValDateAsDouble = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// BOOLEAN TYPE: Imposed Curve or not
	THE_VECT	= parameters->GetColVect("IMPOSED_ZC_CURVE");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters:  no Imposed_ZC_Curve vector found");

	its_ImposedIRCurveFlag = (THE_VECT->Elt(0) ? false : true);
	// ------------------------------------------------------

}


//--------------------------------------------
//	IMPOSED IR CURVES
//--------------------------------------------

void	CreditManager::SetImposedIRCurves(ICM_Matrix<ARM_Vector>* parameters)
{
	if (parameters == NULL) return;
	
	ARM_Vector* THE_VECT=NULL;
	int	i, size;

	// ------------------------------------------------------
	// IR_LAGS:	VECTOR OF INT

	THE_VECT	= parameters->GetColVect("IR_LAGS");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no IR_LAGS vector found");

	size =	THE_VECT->GetSize();	
	its_IR_Lags.resize(size);

	for (i=0;i<size;i++)
		its_IR_Lags[i]	=	(int) (*THE_VECT)[i];
	// ------------------------------------------------------

	// ------------------------------------------------------
	// IR_ZC_YIELDS:	VECTOR OF DOUBLE

	THE_VECT	= parameters->GetColVect("IR_ZC_YIELDS");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no IR_ZC_YIEDLS vector found");

	if (size !=	THE_VECT->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad IR_ZC_YIEDLS vector found");
	
	its_IR_ZC_Yields.resize(size);

	for (i=0;i<size;i++)
		its_IR_ZC_Yields[i]	=	(double) (*THE_VECT)[i];		
	// ------------------------------------------------------

}


//--------------------------------------------
//	CREDIT DATA PARAMETERS
//--------------------------------------------

void	CreditManager::SetCreditDataParameters(ICM_Matrix<ARM_Vector>* parameters)
{
	if (parameters == NULL) return;
	
	ARM_Vector* THE_VECT=NULL;
	int TheValue;

	// ------------------------------------------------------
	// ROLL_DATE:	BOOLEAN

	THE_VECT	= parameters->GetColVect("ROLL_DATE");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad ROLL_DATE");

	its_RollDateFlag = (THE_VECT->Elt(0) ? false : true);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// HEDGES_TYPE:	INT --> ENUM

	THE_VECT	= parameters->GetColVect("HEDGES_TYPE");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad HEDGES_TYPE");

	TheValue = (int) THE_VECT->Elt(0);		// enum type, must check the values
	if ((TheValue < CHB_NO) || (TheValue > CHB_CORRELATION)) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad HEDGES_TYPE");
	
	its_Bump_Choice = (CreditHedgeBump) TheValue;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// BUMP_SPREAD:	DOUBLE

	THE_VECT	= parameters->GetColVect("BUMP_SPREAD");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad BUMP_SPREAD");

	its_BumpSpread = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// BUMP_SPREAD_TYPE:	int

	THE_VECT	= parameters->GetColVect("BUMP_SPREAD_TYPE");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad BUMP_SPREAD_TYPE");

	TheValue = (int) THE_VECT->Elt(0);		// enum type, must check the values
	if ((TheValue < BT_ADD) || (TheValue > BT_MULT)) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad BUMP_SPREAD_TYPE");
	
	its_BumpSpread_Type = (BumpType) TheValue;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// BUMP_RECOVERY:	DOUBLE

	THE_VECT	= parameters->GetColVect("BUMP_RECOVERY");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad BUMP_RECOVERY");

	its_BumpRecovery = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// BUMP_CORRELATION:	DOUBLE

	THE_VECT	= parameters->GetColVect("BUMP_CORRELATION");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad BUMP_CORRELATION");

	its_BumpCorrelation = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CDO_SQUARE_NB_UNDERLYINGS:	INT

	THE_VECT	= parameters->GetColVect("CDO_SQUARE_NB_UNDERLYINGS");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CDO_SQUARE_NB_UNDERLYINGS");

	its_CDO_Square_Nb_Underlyings = (int) THE_VECT->Elt(0);
	// ------------------------------------------------------

}


//--------------------------------------------
//	DATA DESCRIPTION
//--------------------------------------------

void	CreditManager::SetCreditDataDescription(ICM_Matrix<ARM_Vector>* parameters)
{
	if (parameters == NULL) return;
	
	ARM_Vector* THE_VECT=NULL;
	int	i, size;

	// ------------------------------------------------------
	// CATEGORY:	VECTOR OF INT --> TO MATCH WITH ENUM TYPE

	THE_VECT	= parameters->GetColVect("CATEGORY");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no CATEGORY vector found");

	size =	THE_VECT->GetSize();	
	its_Categories.resize(size);

	for (i=0;i<size;i++)
		its_Categories[i]	=	(CreditCategory) (int) (*THE_VECT)[i];
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CURRENCY:	VECTOR OF INT --> TO MATCH WITH ENUM TYPE

	THE_VECT	= parameters->GetColVect("CURRENCY");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no CURRENCY vector found");

	if (size !=	THE_VECT->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CURRENCY vector found");
	
	its_Currencies.resize(size);

	for (i=0;i<size;i++)
		its_Currencies[i]	=	(CurrencyName) (int) (*THE_VECT)[i];		
	// ------------------------------------------------------

	// ------------------------------------------------------
	// ACCRUED:	VECTOR OF INT

	THE_VECT	= parameters->GetColVect("ACCRUED");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no ACCRUED vector found");

	if (size !=	THE_VECT->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad ACCRUED vector found");
	
	its_Accrueds.resize(size);

	for (i=0;i<size;i++)
		its_Accrueds[i]	=	(int) (*THE_VECT)[i];		
	// ------------------------------------------------------

	// ------------------------------------------------------
	// RECOVERY:	VECTOR OF DOUBLE

	THE_VECT	= parameters->GetColVect("RECOVERY");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no RECOVERY vector found");

	if (size !=	THE_VECT->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad RECOVERY vector found");
	
	its_Recoveries.resize(size);

	for (i=0;i<size;i++)
		its_Recoveries[i]	=	(double) (*THE_VECT)[i];		
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NOTIONAL:	VECTOR OF DOUBLE

	THE_VECT	= parameters->GetColVect("NOTIONAL");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no NOTIONAL vector found");

	if (size !=	THE_VECT->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad NOTIONAL vector found");
	
	its_Notionals.resize(size);

	for (i=0;i<size;i++)
		its_Notionals[i]	=	(double) (*THE_VECT)[i];		
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS:	VECTOR OF DOUBLE

	THE_VECT	= parameters->GetColVect("LOSS");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no LOSS vector found");

	if (size !=	THE_VECT->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad LOSS vector found");
	
	its_Input_Losses.resize(size);

	for (i=0;i<size;i++)
		its_Input_Losses[i]	=	(double) (*THE_VECT)[i];		
	// ------------------------------------------------------

	// ------------------------------------------------------
	// DEFAULT_DATE:	VECTOR OF DOUBLE

	THE_VECT	= parameters->GetColVect("DEFAULT_DATE");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no DEFAULT_DATE vector found");

	if (size !=	THE_VECT->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad DEFAULT_DATE vector found");
	
	its_DefaultDates.resize(size);

	for (i=0;i<size;i++)
		its_DefaultDates[i]	=	(RelativeDate) ((*THE_VECT)[i] - itsValDateAsDouble);	//	relative to AsOf
	// ------------------------------------------------------

	// ------------------------------------------------------
	// AMORTIZATION_DATE:	VECTOR OF DOUBLE

	THE_VECT	= parameters->GetColVect("AMORTIZATION_DATE");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no AMORTIZATION_DATE vector found");

	if (size !=	THE_VECT->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad AMORTIZATION_DATE vector found");
	
	its_AmortizationDates.resize(size);

	for (i=0;i<size;i++)
		its_AmortizationDates[i]	=	(RelativeDate) ((*THE_VECT)[i] - itsValDateAsDouble);	//	relative to AsOf
	// ------------------------------------------------------
}


// DO I HAVE TO THROW ERRORS or RATHER JUST SET DEFAULT VALUES?

void	CreditManager::SetCreditModelParameters(ICM_Matrix<ARM_Vector>* parameters)
{
	if (parameters == NULL) return;
	
	ARM_Vector* THE_VECT=NULL;
	int TheValue;

	// ------------------------------------------------------
	// MODEL_TYPE:	INT --> ENUM

	THE_VECT	= parameters->GetColVect("MODEL_TYPE");
	if (!THE_VECT)
		its_CreditModelType = CMT_ANALYTIC_LHP;
	else
	{
		TheValue = (int) THE_VECT->Elt(0);	
		if ((TheValue < CMT_MONTECARLO) || (TheValue > CMT_ANALYTIC_FFT_1F)) 
			its_CreditModelType = CMT_ANALYTIC_LHP;
		else		
			its_CreditModelType = (CreditModelType) TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// KEEP_CALIBRATION:	BOOLEAN

	THE_VECT	= parameters->GetColVect("KEEP_CALIBRATION");
	if (!THE_VECT)
		its_KeepCalibrationFlag = true;
	else
		its_KeepCalibrationFlag = (THE_VECT->Elt(0) ? true : false);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CORRELATION_TYPE:	INT --> ENUM

	THE_VECT	= parameters->GetColVect("CORRELATION_TYPE");
	if (!THE_VECT)
		its_CorrelationType = CT_FLAT;
	else
	{
		TheValue = (int) THE_VECT->Elt(0);
		if ((TheValue < CT_FLAT) || (TheValue > CT_BASE_STRIKE)) 			
			its_CorrelationType = CT_FLAT;
		else		
			its_CorrelationType = (CorrelationType) TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// COPULA_TYPE:	INT --> ENUM

	THE_VECT	= parameters->GetColVect("COPULA_TYPE");
	if (!THE_VECT)
		its_CopulaType = CCT_GAUSSIAN;
	else
	{
		TheValue = (int) THE_VECT->Elt(0);		// enum type, must check the values
		if ((TheValue < CCT_GAUSSIAN) || (TheValue > CCT_NIG)) 
			its_CopulaType = CCT_GAUSSIAN;
		else
			its_CopulaType = (CreditCopulaType) TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NB_SIMUL:	DOUBLE

	THE_VECT	= parameters->GetColVect("NB_SIMUL");
	if (!THE_VECT)
		its_NbSimul = 1000;
	else
		its_NbSimul = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------


	// ------------------------------------------------------
	// its_FreedomDegree:	INT

	THE_VECT	= parameters->GetColVect("FREEDOM_DEGREE");
	if (!THE_VECT)
		its_FreedomDegree = 4;
	else
	{
		TheValue = (int) THE_VECT->Elt(0);		// must be > 0
		if (TheValue <= 0)
			its_FreedomDegree = 4;
		else
			its_FreedomDegree = TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// its_NIntegration_1F:	INT

	THE_VECT	= parameters->GetColVect("N_INTEGRATION_1F");
	if (!THE_VECT)
		its_NIntegration_1F = 51;
	else
	{
		TheValue = (int) THE_VECT->Elt(0);		// must be > 0
		if (TheValue <= 0)
			its_NIntegration_1F = 51;
		else
			its_NIntegration_1F = TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// its_N_FFT:	INT

	THE_VECT	= parameters->GetColVect("N_FFT_1F");
	if (!THE_VECT)
		its_N_FFT = 512;
	else
	{
		TheValue = (int) THE_VECT->Elt(0);		// must be > 0 and even
		if (TheValue <= 0)
			its_N_FFT = 512;
		else
			its_N_FFT = TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// VARIANCE_REDUCTION:	INT --> ENUM

	THE_VECT	= parameters->GetColVect("VARIANCE_REDUCTION");
	if (!THE_VECT)
		its_MC_Variance_Reduction	= CMCVR_NONE;
	else
	{
		TheValue = (int) THE_VECT->Elt(0);		// enum type, must check the values
		if ((TheValue < CMCVR_NONE) || (TheValue > CMCVR_IS_FACTORS)) 
			its_MC_Variance_Reduction	= CMCVR_NONE;
		else		
			its_MC_Variance_Reduction	= (CreditMonteCarloVarianceReduction) TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// IS_THETA:	DOUBLE

	THE_VECT	= parameters->GetColVect("IS_THETA");
	if (!THE_VECT)
		its_IS_Theta = 0.0;
	else
		its_IS_Theta = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// IS_MU:	DOUBLE

	THE_VECT	= parameters->GetColVect("IS_MU");
	if (!THE_VECT)
		its_IS_Mu = 0.0;
	else
		its_IS_Mu = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// IS_LOSS_LEVEL:	DOUBLE

	THE_VECT	= parameters->GetColVect("IS_LOSS_LEVEL");
	if (!THE_VECT)
		its_IS_Loss_Level = 0.0;
	else
		its_IS_Loss_Level = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// IS_LOSS_MATURITY:	DOUBLE

	THE_VECT	= parameters->GetColVect("IS_LOSS_MATURITY");
	if (!THE_VECT)
		its_IS_Loss_Maturity = 5.0;
	else
		its_IS_Loss_Maturity = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS_UNIT_CHOICE:	INT --> ENUM

	THE_VECT	= parameters->GetColVect("LOSS_UNIT_CHOICE");
	if (!THE_VECT)
		its_Recursive_1F_Loss_Unit_Choice	= COFRLUC_PGCD;
	else
	{
		TheValue = (int) THE_VECT->Elt(0);		// enum type, must check the values
		if ((TheValue < COFRLUC_PGCD) || (TheValue > COFRLUC_NB_LOSS_STEP)) 
			its_Recursive_1F_Loss_Unit_Choice	= COFRLUC_PGCD;
		else		
			its_Recursive_1F_Loss_Unit_Choice	= (CreditOneFactorRecursiveLossUnitChoice) TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS_UNIT_MIN:	DOUBLE

	THE_VECT	= parameters->GetColVect("LOSS_UNIT_MIN");
	if (!THE_VECT)
		its_Recursive_1F_Loss_Unit_Min = 0.0;
	else
		its_Recursive_1F_Loss_Unit_Min = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS_UNIT_NB_STEP:	INT

	THE_VECT	= parameters->GetColVect("LOSS_UNIT_NB_STEP");
	if (!THE_VECT)
		its_Recursive_1F_NbLossStep = 0;
	else
	{
		TheValue = (int) THE_VECT->Elt(0);		// must be > 0 and even
		if (TheValue <= 0)
			its_Recursive_1F_NbLossStep = 0;
		else
			its_Recursive_1F_NbLossStep = TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// TARGET_DEFAULT_PROB:	DOUBLE

	THE_VECT	= parameters->GetColVect("TARGET_DEFAULT_PROB");
	if (!THE_VECT)
		its_Target_Default_Prob = 0.0;
	else
		its_Target_Default_Prob = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// VARIANCE_REDUCTION_THETA_CHOICE:	INT --> ENUM

	THE_VECT	= parameters->GetColVect("VARIANCE_REDUCTION_THETA_CHOICE");
	if (!THE_VECT)
		its_Theta_Choice	= CVRI_IMPOSED;
	else
	{
		TheValue = (int) THE_VECT->Elt(0);		// enum type, must check the values
		if ((TheValue < CVRI_IMPOSED) || (TheValue > CVRI_OPTIM)) 
			its_Theta_Choice	= CVRI_IMPOSED;
		else		
			its_Theta_Choice	= (CreditVarianceReductionIdiosyncratic) TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// VARIANCE_REDUCTION_MU_CHOICE:	INT --> ENUM

	THE_VECT	= parameters->GetColVect("VARIANCE_REDUCTION_MU_CHOICE");
	if (!THE_VECT)
		its_Mu_Choice	= CVRF_TARGET_DEFPROB;
	else
	{
		TheValue = (int) THE_VECT->Elt(0);		// enum type, must check the values
		if ((TheValue < CVRF_IMPOSED) || (TheValue > CVRF_TARGET_DEFPROB)) 
			its_Mu_Choice	= CVRF_TARGET_DEFPROB;
		else
			its_Mu_Choice	= (CreditVarianceReductionFactors) TheValue;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LHP_MATURITY:	DOUBLE

	THE_VECT	= parameters->GetColVect("LHP_MATURITY");
	if (!THE_VECT)
		its_LHP_Maturity = 5.0 * 365.0;
	else
		its_LHP_Maturity = (RelativeDate) (THE_VECT->Elt(0) - itsValDateAsDouble);	//	relative to AsOf
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LHP_SPREAD:	DOUBLE

	THE_VECT	= parameters->GetColVect("LHP_SPREAD");
	if (!THE_VECT)
		its_LHP_Spread = 0.0;
	else
		its_LHP_Spread = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LHP_RECOVERY:	DOUBLE

	THE_VECT	= parameters->GetColVect("LHP_RECOVERY");
	if (!THE_VECT)
		its_LHP_Recovery = 0.0;
	else
		its_LHP_Recovery = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LP_DEGREE:	INT

	THE_VECT	= parameters->GetColVect("LP_DEGREE");
	if (!THE_VECT)
		its_LP_Degree = 0;
	else
		its_LP_Degree = (int) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NIG_ALPHA:	DOUBLE

	THE_VECT	= parameters->GetColVect("NIG_ALPHA");
	if (!THE_VECT)
		its_NIG_Alpha	= 0.0;
	else
		its_NIG_Alpha	= (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NIG_BETA:	DOUBLE

	THE_VECT	= parameters->GetColVect("NIG_BETA");
	if (!THE_VECT)
		its_NIG_Beta	= 0.0;
	else
	{
		its_NIG_Beta	= (double) THE_VECT->Elt(0);

		if (its_NIG_Beta > its_NIG_Alpha )
			its_NIG_Beta	= its_NIG_Alpha;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NIG_MU:	DOUBLE

	THE_VECT	= parameters->GetColVect("NIG_MU");
	if (!THE_VECT)
		its_NIG_Mu = 0.0;
	else
		its_NIG_Mu = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NIG_DELTA:	DOUBLE

	THE_VECT	= parameters->GetColVect("NIG_DELTA");
	if (!THE_VECT)
		its_NIG_Delta = 0.0;
	else
	{
		its_NIG_Delta = (double) THE_VECT->Elt(0);

		if (its_NIG_Delta < 0. )
			its_NIG_Delta	=	0.0;
	}
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NIG_RHO:	DOUBLE

	THE_VECT	= parameters->GetColVect("NIG_RHO");
	if (!THE_VECT)
		its_NIG_Rho = 0.0;
	else
		its_NIG_Rho = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// SC_RHO_1:	DOUBLE

	THE_VECT	= parameters->GetColVect("SC_RHO_1");
	if (!THE_VECT)
		its_SC_Rho1	= 0.0;
	else
		its_SC_Rho1	= (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// SC_RHO_2:	DOUBLE

	THE_VECT	= parameters->GetColVect("SC_RHO_2");
	if (!THE_VECT)
		its_SC_Rho2	= 0.0;
	else
		its_SC_Rho2	= (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// SC_RHO_PROB:	DOUBLE

	THE_VECT	= parameters->GetColVect("SC_RHO_PROB");
	if (!THE_VECT)
		its_SC_Prob	= 0.0;
	else
		its_SC_Prob	= (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// DISPLAY_DEFAULT_TIMES OUTPUT:	BOOLEAN

	THE_VECT	= parameters->GetColVect("DISPLAY_DEFAULT_TIMES");
	if (!THE_VECT)
		its_Display_DefaultTimes_Flag = true;
	else
		its_Display_DefaultTimes_Flag = (THE_VECT->Elt(0) ? true : false);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// DISPLAY_LOSSES_AMOUNT OUTPUT:	BOOLEAN

	THE_VECT	= parameters->GetColVect("DISPLAY_LOSSES_AMOUNT");
	if (!THE_VECT)
		its_Display_LossesAmount_Before_CreditDates_Flag = true;
	else
		its_Display_LossesAmount_Before_CreditDates_Flag = (THE_VECT->Elt(0) ? true : false);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// DISPLAY_SORTED_DEFAULT_TIMES OUTPUT:	BOOLEAN

	THE_VECT	= parameters->GetColVect("DISPLAY_SORTED_DEFAULT_TIMES");
	if (!THE_VECT)
		its_Display_Sorted_DefaultTimes_Flag = true;
	else
		its_Display_Sorted_DefaultTimes_Flag = (THE_VECT->Elt(0) ? true : false);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// DISPLAY_CUMULATIVE_LOSSES OUTPUT:	BOOLEAN

	THE_VECT	= parameters->GetColVect("DISPLAY_CUMULATIVE_LOSSES");
	if (!THE_VECT)
		its_Display_CumulativeLossesAndDefaults_Flag = true;
	else
		its_Display_CumulativeLossesAndDefaults_Flag = (THE_VECT->Elt(0) ? true : false);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CORRELATION_RFL:	INT --> ENUM

	THE_VECT	= parameters->GetColVect("CORRELATION_RFL");
	if (!THE_VECT)
		its_Correlation_RFL_Type = CRFL_2F_INDIVIDUAL;
	else
	{
		TheValue = (int) THE_VECT->Elt(0);
		if ((TheValue < CRFL_2F_CONSTANT) || (TheValue > CRFL_TANH)) 			
			its_Correlation_RFL_Type = CRFL_2F_INDIVIDUAL;
		else		
			its_Correlation_RFL_Type = (Correlation_RFL) TheValue;
	}
	// ------------------------------------------------------


}

void	CreditManager::SetCreditProductDefaultMatrix(ICM_Matrix<ARM_Vector>* parameters)
{
	int	TheValue;
	if (parameters == NULL) return;

	ARM_Vector* DATA	=	NULL;
	
	// ------------------------------------------------------
	// CREDIT_WINDOW_LOW:	RELATIVE DATE

	DATA	= parameters->GetColVect("CREDIT_WINDOW_LOW");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no CREDIT_WINDOW_LOW vector found");

	its_DL_CreditWindowLow = (RelativeDate) (DATA->Elt(0) - itsValDateAsDouble);	//	relative to AsOf
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CREDIT_WINDOW_UP:	RELATIVE DATE

	DATA	= parameters->GetColVect("CREDIT_WINDOW_UP");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no CREDIT_WINDOW_UP vector found");

	its_DL_CreditWindowUp = (RelativeDate) (DATA->Elt(0) - itsValDateAsDouble);	//	relative to AsOf
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS_MIN: DOUBLE

	DATA	= parameters->GetColVect("LOSS_MIN");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no LOSS_MIN vector found");

	its_DL_LossMin = (double) DATA->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS_MAX: DOUBLE

	DATA	= parameters->GetColVect("LOSS_MAX");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no LOSS_MAX vector found");

	its_DL_LossMax = (double) DATA->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PAYMENT_TYPE:	RELATIVE DATE

	DATA	= parameters->GetColVect("PAYMENT_TYPE");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no PAYMENT_TYPE vector found");

	TheValue	=	(int)	DATA->Elt(0);
	if ((TheValue < CEP_ATDEFAULTDATE) || (TheValue > CEP_ATFIXEDDATE))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad PAYMENT_TYPE item found");

	its_DL_PaymentType = (CreditEventPayment) TheValue;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PAYMENT_DATE:	RELATIVE DATE

	DATA	= parameters->GetColVect("PAYMENT_DATE");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no PAYMENT_DATE vector found");

	its_DL_PaymentDate = (RelativeDate) (DATA->Elt(0) - itsValDateAsDouble);	//	relative to AsOf
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PAYMENT LAG: INT

	DATA	= parameters->GetColVect("PAYMENT_LAG");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no PAYMENT_LAG vector found");

	its_DL_PaymentLag = (int) DATA->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NB_MIN : INT

	DATA	= parameters->GetColVect("NB_MIN");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no NB_MIN vector found");

	its_DL_NbDefMin = (int) DATA->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// NB_MAX : INT

	DATA	= parameters->GetColVect("NB_MAX");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no NB_MAX vector found");

	its_DL_NbDefMax = (int) DATA->Elt(0);
	// ------------------------------------------------------
}


void	CreditManager::SetCreditProductPremiumMatrix(ICM_Matrix<ARM_Vector>* parameters)
{
	if (parameters == NULL) return;

	ARM_Vector* DATA	=	NULL;
	int	i, size;
	int	TheValue;

	// ------------------------------------------------------
	// CREDIT_WINDOW_LOW:	VECTOR OF RELATIVE DATE

	DATA	= parameters->GetColVect("CREDIT_WINDOW_LOW");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no CREDIT_WINDOW_LOW vector found");

	size	=	DATA->GetSize();
	its_PL_NbFlows	=	size;
	its_PL_CreditWindowLows.resize(size);

	for (i=0;i<size;i++)
		its_PL_CreditWindowLows[i]	=	(RelativeDate) ((*DATA)[i] - itsValDateAsDouble);	//	relative to AsOf

	// ------------------------------------------------------

	// ------------------------------------------------------
	// CREDIT_WINDOW_UP:	VECTOR OF RELATIVE DATE

	DATA	= parameters->GetColVect("CREDIT_WINDOW_UP");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no CREDIT_WINDOW_UP vector found");

	if (size !=	DATA->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CREDIT_WINDOW_UP vector found");
	
	its_PL_CreditWindowUps.resize(size);

	for (i=0;i<size;i++)
		its_PL_CreditWindowUps[i]	=	(RelativeDate) ((*DATA)[i] - itsValDateAsDouble);	//	relative to AsOf

	// ------------------------------------------------------

	// ------------------------------------------------------
	// START DATES:	VECTOR OF RELATIVE DATE

	DATA	= parameters->GetColVect("START_DATE");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no START_DATE vector found");

	if (size !=	DATA->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad START_DATE vector found");
	
	its_PL_StartDates.resize(size);

	for (i=0;i<size;i++)
		its_PL_StartDates[i]	=	(RelativeDate) ((*DATA)[i] - itsValDateAsDouble);	//	relative to AsOf

	// ------------------------------------------------------

	// ------------------------------------------------------
	// END DATES:	VECTOR OF RELATIVE DATE

	DATA	= parameters->GetColVect("END_DATE");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no END_DATE vector found");

	if (size !=	DATA->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad END_DATE vector found");
	
	its_PL_EndDates.resize(size);

	for (i=0;i<size;i++)
		its_PL_EndDates[i]	=	(RelativeDate) ((*DATA)[i] - itsValDateAsDouble);	//	relative to AsOf

	// ------------------------------------------------------

	// ------------------------------------------------------
	// PAYMENT DATES:	VECTOR OF RELATIVE DATE

	DATA	= parameters->GetColVect("PAYMENT_DATE");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no PAYMENT_DATE vector found");

	if (size !=	DATA->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad PAYMENT_DATE vector found");
	
	its_PL_PaymentDates.resize(size);

	for (i=0;i<size;i++)
		its_PL_PaymentDates[i]	=	(RelativeDate) ((*DATA)[i] - itsValDateAsDouble);	//	relative to AsOf
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LOSS MINS:	VECTOR OF DOUBLE

	DATA	= parameters->GetColVect("LOSS_MIN");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no LOSS MINS vector found");

	if (size !=	DATA->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad LOSS MINS vector found");
	
	its_PL_LossMins.resize(size);

	for (i=0;i<size;i++)
		its_PL_LossMins[i]	=	(*DATA)[i];		

	// ------------------------------------------------------
	// LOSS MAXS:	VECTOR OF DOUBLE

	DATA	= parameters->GetColVect("LOSS_MAX");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no LOSS MAXS vector found");

	if (size !=	DATA->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad LOSS MAXS vector found");
	
	its_PL_LossMaxs.resize(size);

	for (i=0;i<size;i++)
		its_PL_LossMaxs[i]	=	(*DATA)[i];		

	// ------------------------------------------------------
	// NB MINS:	VECTOR OF INT

	DATA	= parameters->GetColVect("NB_MINS");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no NB_MINS vector found");

	if (size !=	DATA->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad NB MINS vector found");
	
	its_PL_NbDefMins.resize(size);

	for (i=0;i<size;i++)
		its_PL_NbDefMins[i]	=	(int) (*DATA)[i];		

	// ------------------------------------------------------
	// NB MAXS:	VECTOR OF INT

	DATA	= parameters->GetColVect("NB_MAXS");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no NB_MAXS vector found");

	if (size !=	DATA->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad NB MAXS vector found");
	
	its_PL_NbDefMaxs.resize(size);

	for (i=0;i<size;i++)
		its_PL_NbDefMaxs[i]	=	(int) (*DATA)[i];		

	// ------------------------------------------------------
	// RATIOS:	VECTOR OF DOUBLE

	DATA	= parameters->GetColVect("RATIOS");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no RATIOS vector found");

	if (size !=	DATA->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad RATIOS vector found");
	
	its_PL_Ratios.resize(size);

	for (i=0;i<size;i++)
		its_PL_Ratios[i]	=	(*DATA)[i];		

	// ------------------------------------------------------

	// ------------------------------------------------------
	// NOTIOS:	VECTOR OF DOUBLE

	DATA	= parameters->GetColVect("NOTIOS");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no NOTIOS vector found");

	if (size !=	DATA->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad NOTIOS vector found");
	
	its_PL_Notios.resize(size);

	for (i=0;i<size;i++)
		its_PL_Notios[i]	=	(*DATA)[i];		

	// ------------------------------------------------------

	// ------------------------------------------------------
	// SPREADS:	VECTOR OF DOUBLE

	DATA	= parameters->GetColVect("SPREADS");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no SPREADS vector found");

	if (size !=	DATA->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad SPREADS vector found");
	
	its_PL_Spreads.resize(size);

	for (i=0;i<size;i++)
		its_PL_Spreads[i]	=	(*DATA)[i];		

	// ------------------------------------------------------

	// ------------------------------------------------------
	// SPREAD CAPS:	VECTOR OF DOUBLE

	DATA	= parameters->GetColVect("SPREAD_CAP");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no SPREAD CAPS vector found");

	if (size !=	DATA->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad SPREAD CAPS vector found");
	
	its_PL_CreditSpreadCaps.resize(size);

	for (i=0;i<size;i++)
		its_PL_CreditSpreadCaps[i]	=	(*DATA)[i];		

	// ------------------------------------------------------

	// ------------------------------------------------------
	// REDEMPTION:	VECTOR OF DOUBLE

	DATA	= parameters->GetColVect("REDEMPTION");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no REDEMPTION vector found");

	if (size !=	DATA->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad REDEMPTION vector found");
	
	its_PL_Redemptions.resize(size);

	for (i=0;i<size;i++)
		its_PL_Redemptions[i]	=	(*DATA)[i];		

	// ------------------------------------------------------

	// ------------------------------------------------------
	// CREDIT FLAGS:	VECTOR OF DOUBLE

	DATA	= parameters->GetColVect("CREDIT_FLAG");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no CREDIT FLAGS vector found");

	if (size !=	DATA->GetSize())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CREDIT FLAGS vector found");
	
	its_PL_CreditFlags.resize(size);

	for (i=0;i<size;i++)
	{
		TheValue	=	(int) (*DATA)[i];
		if ((TheValue < CPLT_GUARANTEED) || (TheValue > CPLT_CUMULREDEMPTION)) 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
				"Parameters :  bad Credit Flag for Premium Leg");

		its_PL_CreditFlags[i]	=	(CreditPremiumLegType) TheValue;		
	}

	// ------------------------------------------------------
}


void	CreditManager::SetCreditProductPricingParametersMatrix(ICM_Matrix<ARM_Vector>* parameters)
{
	if (parameters == NULL) return;
	int	TheValue;

	ARM_Vector* DATA	=	NULL;
	
	// ------------------------------------------------------
	// PRICING_LEGS:	INT --> ENUM

	DATA	= parameters->GetColVect("PRICING_LEGS");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no PRICING_LEGS item found");

	TheValue = (int) DATA->Elt(0);
	if ((TheValue < CBN_DEFAULTMINUSPREMIUM) || (TheValue > CBN_DEFAULTLEGONLY))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad PRICING_LEGS item found");
	
	its_PricingLegsType	=	(CreditBasketNPV) TheValue;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CREDIT_OBSERVATION:	INT --> ENUM

	DATA	= parameters->GetColVect("CREDIT_OBSERVATION");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no PRICING_LEGS item found");

	TheValue = (int) DATA->Elt(0);
	if ((TheValue < CO_NBDEFAULTS) || (TheValue > CO_LOSSES))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CREDIT_OBSERVATION item found");
	
	its_CreditObservationType	=	(CreditObservation) TheValue;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PRICE_FORMAT:	INT --> ENUM

	DATA	= parameters->GetColVect("PRICE_FORMAT");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no PRICING_LEGS item found");

	TheValue = (int) DATA->Elt(0);
	if ((TheValue < CPLADT_PURESPREAD) || (TheValue > CPLADT_UF_SPREAD))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad PRICE_FORMAT item found");
	
	its_ATMDataFlag	=	(CreditPremiumLegATMData) TheValue;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// UP_FRONT_PREMIUM:	DOUBLE

	DATA	= parameters->GetColVect("UP_FRONT_PREMIUM");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no PRICING_LEGS UP_FRONT_PREMIUM item found");	
	
	its_UpFront_Premium	=	(double) DATA->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// UP_FRONT_RUNNING_SPREAD:	DOUBLE

	DATA	= parameters->GetColVect("UP_FRONT_RUNNING_SPREAD");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no PRICING_LEGS UP_FRONT_RUNNING_SPREAD item found");	
	
	its_UpFront_RunningSpread	=	(double) DATA->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PREMIUM_ACCRUED:	INT --> ENUM

	DATA	= parameters->GetColVect("PREMIUM_ACCRUED");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no PREMIUM_ACCRUED item found");

	TheValue = (int) DATA->Elt(0);
	if ((TheValue < CAP_PRORATA) || (TheValue > CAP_NON_PRORATA))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad PREMIUM_ACCRUED item found");
	
	its_CreditPremiumLegAccrued	=	(CreditAccruedPayment) TheValue;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PREMIUM_ACCRUED_PAYMENT:	INT --> ENUM

	DATA	= parameters->GetColVect("PREMIUM_ACCRUED_PAYMENT");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no PREMIUM_ACCRUED_PAYMENT item found");

	TheValue = (int) DATA->Elt(0);
	if ((TheValue < CEP_ATDEFAULTDATE) || (TheValue > CEP_ATFIXEDDATE))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad PREMIUM_ACCRUED_PAYMENT item found");
	
	its_PL_PaymentType	=	(CreditEventPayment) TheValue;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// TIME_STEP_1F_PRORATA_NB_DAYS:	DOUBLE

	DATA	= parameters->GetColVect("TIME_STEP_1F_PRORATA_NB_DAYS");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no TIME_STEP_1F_PRORATA_NB_DAYS item found");

	its_Time_Step_Prorata_1F	=	(double) DATA->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// COMPUTE_HEDGES:	INT: 0 or 1

	DATA	= parameters->GetColVect("COMPUTE_HEDGES");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no COMPUTE_HEDGES item found");

	TheValue = (int) DATA->Elt(0);
	if (TheValue)
		its_HedgesRunning	=	true;
	else
		its_HedgesRunning	=	false;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PRODUCT_TYPE:	INT --> ENUM

	DATA	= parameters->GetColVect("PRODUCT_TYPE");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no PRODUCT_TYPE item found");

	TheValue = (int) DATA->Elt(0);
	if ((TheValue < CPT_STD_CDO) || (TheValue > CPT_CDO_SQUARE))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad PRODUCT_TYPE item found");
	
	its_CDO_Type	=	(CreditProductType) TheValue;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PRICE_TYPE:	INT --> ENUM

	DATA	= parameters->GetColVect("PRICE_TYPE");
	if (!DATA)
//		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
 //          "Parameters :  no PRICING_LEGS item found");
		TheValue =	CNPV_STANDARD;
	else
	{
		TheValue = (int) DATA->Elt(0);

		if ((TheValue < CNPV_STANDARD) || (TheValue > CNPV_WITH_RUNNING_SPREAD_AND_UP_FRONT_PREMIUM))
			throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
				"Parameters :  bad PRICE_TYPE item found");
		
		its_NPV_Type	=	(CreditNPV) TheValue;
	}
	// ------------------------------------------------------

}


void	CreditManager::SetFactorLoadingsParameters(ICM_Matrix<ARM_Vector>* parameters)
{
	if (parameters == NULL) return;

	ARM_Vector* DATA	=	NULL;
	int	i, size;

	// ------------------------------------------------------
	// ALPHA coeff

	DATA	= parameters->GetColVect("ALPHA");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no ALPHA vector found");

	size	=	DATA->GetSize();
//	if (size != its_NbCredits)
//		ICMTHROW(ERR_INVALID_DATA,"Bad Size for Factor Loadings Coefficients! " << size << " instead of " << its_NbCredits);

	its_FL_Alpha.resize(size);

	for (i=0;i<size;i++)
		its_FL_Alpha[i]	=	(double) (*DATA)[i];

	// ------------------------------------------------------

	// ------------------------------------------------------
	// BETA_1 coeff

	DATA	= parameters->GetColVect("BETA_1");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no BETA_1 vector found");

	if (size != DATA->GetSize())
		ICMTHROW(ERR_INVALID_DATA,"Bad Size for Factor Loadings Coefficients! " << size << " instead of " << its_CreditsLabels.size());

	its_FL_Beta1.resize(size);

	for (i=0;i<size;i++)
		its_FL_Beta1[i]	=	(double) (*DATA)[i];

	// ------------------------------------------------------

	// ------------------------------------------------------
	// BETA_2 coeff

	DATA	= parameters->GetColVect("BETA_2");
	if (!DATA)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  no BETA_2 vector found");

	if (size != DATA->GetSize())
		ICMTHROW(ERR_INVALID_DATA,"Bad Size for Factor Loadings Coefficients! " << size << " instead of " << its_CreditsLabels.size());

	its_FL_Beta2.resize(size);

	for (i=0;i<size;i++)
		its_FL_Beta2[i]	=	(double) (*DATA)[i];

	// ------------------------------------------------------

}


//--------------------------------------------
//	CREDIT DATA DATA: for each Underlying Credit
//--------------------------------------------

void	CreditManager::SetCreditDataCDOSquare_Data(ICM_Matrix<ARM_Vector>* parameters)
{
	if (parameters == NULL) return;
	
	ARM_Vector* THE_VECT=NULL;
	int	i, size;

	if ((its_CDO_Square_Nb_Underlyings <= 0) || (its_CreditsLabels.size() <= 0))
		return;

	// ------------------------------------------------------------------------------------------------------
	// ALLOCATIONS
	its_CDO_Square_Credit_Ids		= new ICM_QMatrix<int>(its_CreditsLabels.size(), its_CDO_Square_Nb_Underlyings);
	its_CDO_Square_Credit_Losses	= new ICM_QMatrix<double>(its_CreditsLabels.size(), its_CDO_Square_Nb_Underlyings);
	its_CDO_Square_Credit_Notionals	= new ICM_QMatrix<double>(its_CreditsLabels.size(), its_CDO_Square_Nb_Underlyings);
	
	its_CDO_Square_CDO_NbCredits.resize(its_CDO_Square_Nb_Underlyings);
	// ------------------------------------------------------------------------------------------------------

	// ------------------------------------------------------------------------------------------------------
	string	TheLabel;
	char	buffer[20];
	int j;
	int	NbCreditsForCurrentCDO;
	// ------------------------------------------------------------------------------------------------------

	for (j=0; j<its_CDO_Square_Nb_Underlyings; j++)
	{
		// ------------------------------------------------------
		// INCLUDE:	VECTOR OF INT
		TheLabel = "INCLUDE_";
		itoa(j+1, buffer, 10);
		TheLabel.append(buffer);

		THE_VECT	= parameters->GetColVect(TheLabel);
		if (!THE_VECT)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
				"Parameters :  no INCLUDE vector found");

		size =	THE_VECT->GetSize();	
		if (size != its_CreditsLabels.size())
			ICMTHROW(ERR_INVALID_DATA,"Bad Size for CDO Square INCLUDE no " << j << ".");

		NbCreditsForCurrentCDO	=	0;

		for (i=0;i<size;i++)
		{

			if ((*THE_VECT)[i])  // ??? CHECK_EQUAL ???
			{
				(*its_CDO_Square_Credit_Ids)(i, j)	=	i+1;		// i is the id+1
				NbCreditsForCurrentCDO++;
			}
			else
				(*its_CDO_Square_Credit_Ids)(i,j)	=	0;		// flag, excluded
			
		}
		its_CDO_Square_CDO_NbCredits[j]	=	NbCreditsForCurrentCDO;
		// ------------------------------------------------------

		// ------------------------------------------------------
		// NOTIONAL:	VECTOR OF DOUBLE
		TheLabel = "NOTIONAL_";
		itoa(j+1, buffer, 10);
		TheLabel.append(buffer);

		THE_VECT	= parameters->GetColVect(TheLabel);
		if (!THE_VECT)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
				"Parameters :  no NOTIONAL vector found");

		size =	THE_VECT->GetSize();	
		if (size != its_CreditsLabels.size())
			ICMTHROW(ERR_INVALID_DATA,"Bad Size for CDO Square NOTIONAL no " << j << ".");

		for (i=0;i<size;i++)
		{
			(*its_CDO_Square_Credit_Notionals)(i, j)	=	(double) (*THE_VECT)[i];
			(*its_CDO_Square_Credit_Losses)(i, j)	=	(double) (*THE_VECT)[i];
		}
		// ------------------------------------------------------

		// ------------------------------------------------------
		// LOSSES:	VECTOR OF DOUBLE
		TheLabel = "LOSS_";
		itoa(j+1, buffer, 10);
		TheLabel.append(buffer);

		THE_VECT	= parameters->GetColVect(TheLabel);
		if (!THE_VECT)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
				"Parameters :  no LOSS vector found");

		size =	THE_VECT->GetSize();	
		if (size != its_CreditsLabels.size())
			ICMTHROW(ERR_INVALID_DATA,"Bad Size for CDO Square LOSS no " << j << ".");

		for (i=0;i<size;i++)
			(*its_CDO_Square_Credit_Losses)(i, j)	*=	(double) (*THE_VECT)[i];
		// ------------------------------------------------------
	}

}


//--------------------------------------------
//	CREDIT DATA PARAMETERS: Underlying CDO descriptions
//--------------------------------------------

void	CreditManager::SetCreditDataCDOSquare_Parameters(ICM_Matrix<ARM_Vector>* parameters)
{
	if (parameters == NULL) return;
	
	ARM_Vector* THE_VECT=NULL;
	int	i, size;

	if ((its_CDO_Square_Nb_Underlyings <= 0) || (its_CreditsLabels.size() <= 0))
		return;

	// ALLOCATION	
	its_CDO_Square_CDO_Includes.resize(its_CDO_Square_Nb_Underlyings);
	its_CDO_Square_CDO_Maturities.resize(its_CDO_Square_Nb_Underlyings);
	its_CDO_Square_CDO_Attachments.resize(its_CDO_Square_Nb_Underlyings);
	its_CDO_Square_CDO_Detachments.resize(its_CDO_Square_Nb_Underlyings);
	its_CDO_Square_CDO_Notionals.resize(its_CDO_Square_Nb_Underlyings);

	string	TheLabel;

	// ------------------------------------------------------
	// CDO_SQR_INCLUDE:	VECTOR OF INT
	TheLabel = "CDO_SQR_INCLUDE";

	THE_VECT	= parameters->GetColVect(TheLabel);
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			"Parameters :  no CDO_SQR_INCLUDE vector found");

	size =	THE_VECT->GetSize();	
	if (size != its_CDO_Square_Nb_Underlyings)
		ICMTHROW(ERR_INVALID_DATA,"Bad Size for CDO Square CDO_SQR_INCLUDE!");

	for (i=0;i<size;i++)
		its_CDO_Square_CDO_Includes[i]	=	(int) (*THE_VECT)[i];
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CDO_SQR_MATURITIES:	VECTOR OF DOUBLE
	TheLabel = "CDO_SQR_MATURITIES";

	THE_VECT	= parameters->GetColVect(TheLabel);
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			"Parameters :  no CDO_SQR_MATURITIES vector found");

	size =	THE_VECT->GetSize();	
	if (size != its_CDO_Square_Nb_Underlyings)
		ICMTHROW(ERR_INVALID_DATA,"Bad Size for CDO Square CDO_SQR_MATURITIES!");

	for (i=0;i<size;i++)
		its_CDO_Square_CDO_Maturities[i]	=	((double) (*THE_VECT)[i] - itsValDateAsDouble);	//	relative to AsOf
	// ------------------------------------------------------


	// ------------------------------------------------------
	// CDO_SQR_NOTIONALS:	VECTOR OF DOUBLE
	TheLabel = "CDO_SQR_NOTIONALS";

	THE_VECT	= parameters->GetColVect(TheLabel);
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			"Parameters :  no CDO_SQR_NOTIONALS vector found");

	size =	THE_VECT->GetSize();	
	if (size != its_CDO_Square_Nb_Underlyings)
		ICMTHROW(ERR_INVALID_DATA,"Bad Size for CDO Square CDO_SQR_NOTIONALS!");

	for (i=0;i<size;i++)
		its_CDO_Square_CDO_Notionals[i]	=	(double) (*THE_VECT)[i];
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CDO_SQR_ATTACHMENTS:	VECTOR OF DOUBLE
	TheLabel = "CDO_SQR_ATTACHMENTS";

	THE_VECT	= parameters->GetColVect(TheLabel);
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			"Parameters :  no CDO_SQR_ATTACHMENTS vector found");

	size =	THE_VECT->GetSize();	
	if (size != its_CDO_Square_Nb_Underlyings)
		ICMTHROW(ERR_INVALID_DATA,"Bad Size for CDO Square CDO_SQR_ATTACHMENTS!");

	for (i=0;i<size;i++)
		its_CDO_Square_CDO_Attachments[i]	=	(double) (*THE_VECT)[i];
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CDO_SQR_DETACHMENTS:	VECTOR OF DOUBLE
	TheLabel = "CDO_SQR_DETACHMENTS";

	THE_VECT	= parameters->GetColVect(TheLabel);
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			"Parameters :  no CDO_SQR_DETACHMENTS vector found");

	size =	THE_VECT->GetSize();	
	if (size != its_CDO_Square_Nb_Underlyings)
		ICMTHROW(ERR_INVALID_DATA,"Bad Size for CDO Square CDO_SQR_DETACHMENTS!");

	for (i=0;i<size;i++)
		its_CDO_Square_CDO_Detachments[i]	=	(double) (*THE_VECT)[i];
	// ------------------------------------------------------
}

// ***********************************************************
// ***********************************************************

// ***********************************************************
// ***********************************************************

// -----------------------------------------------------------
// Useful DATA
// -----------------------------------------------------------

void	CreditManager::ComputeMaxDFDate()
{
	int	i;
	RelativeDate	CurrentDate, MaxDate;
	
	MaxDate	=	0.0;
	
	// I have to loop on Premium Leg (not sorted) and take into account the last Payment Date
	for (i=0;i<its_PL_NbFlows;i++)
	{
		CurrentDate	=	its_PL_PaymentDates[i];
		if (CurrentDate > MaxDate)
			MaxDate	=	CurrentDate;
	}

	// Default Leg, according to Payment Type
	switch (its_DL_PaymentType)
	{
	case CEP_ATNEXTCREDITDATE:
		break;
	case CEP_ATNEXTPAYMENTDATE:
		break;
	case CEP_ATFIXEDDATE:
		CurrentDate	=	its_DL_PaymentDate;
		break;
	case CEP_ATDEFAULTDATE:
		CurrentDate	=	its_DL_CreditWindowUp;
		break;
	}
	// Add the Lag
	CurrentDate	+=	its_DL_PaymentLag;

	its_MaxDate	=	FMAX(MaxDate, CurrentDate);

}

// ----------------------------------------------------
// Computes all DF
// ----------------------------------------------------

void	CreditManager::ComputeAllDF()
{
	int i;

	ComputeMaxDFDate();

	int	size	=	(int) its_MaxDate;
	
	if (size <= 0)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
            "Parameters :  bad Max DF Date");

	its_DF.resize(size+1);	// in order to take into account the date T=0

	if (!its_ImposedIRCurveFlag)
		for (i=0;i<=size;i++)
			its_DF[i] = its_ZeroCurve->DiscountPrice(MATHTIME(i));
	else
	{
		// do the interpolation in the vector of imposed values
		for (i=0;i<size;i++)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Interpolation in a vector<double> must be implemented!");
//			its_DF[i] = its_IR_ZC_Yieldsits_ZeroCurve->DiscountPrice((double) i) ;
		}
	}

	// PREMIUM LAG --> no specific lag
	// Pre Compute DF for the following given dates
	// its_PL_PaymentDates
	// its_PL_CreditWindowUps
	its_PL_Sorted_CreditWindowUps.clear();
	its_PL_CreditWindowUps_DF.clear();
	its_PL_Sorted_CreditWindowUps_DF.clear();

	its_PL_Sorted_PaymentDates.clear();
	its_PL_PaymentDates_DF.clear();
	its_PL_Sorted_PaymentDates_DF.clear();
	
	RelativeDate	TheDate;
	size	=	its_PL_CreditWindowUps.size();

	for (i=0;i<size;i++)
	{
		// SORTED AND NON SORTED ONES
		TheDate	=	its_PL_CreditWindowUps[i];

//		its_PL_Sorted_PaymentDates.insert(TheDate);
		its_PL_CreditWindowUps_DF.push_back(its_DF[(int)TheDate]);
//		its_PL_Sorted_CreditWindowUps_DF.insert(its_DF[(int)TheDate]);

		TheDate	=	its_PL_PaymentDates[i];

//		its_PL_Sorted_PaymentDates.insert(TheDate);
		its_PL_PaymentDates_DF.push_back(its_DF[(int)TheDate]);
//		its_PL_Sorted_PaymentDates_DF.insert(its_DF[(int)TheDate]);
	}

	its_PL_Sorted_CreditWindowUps	= its_PL_CreditWindowUps;
	sort(its_PL_Sorted_CreditWindowUps.begin(), its_PL_Sorted_CreditWindowUps.end());
	its_PL_Sorted_PaymentDates	= its_PL_PaymentDates;
	sort(its_PL_Sorted_PaymentDates.begin(), its_PL_Sorted_PaymentDates.end());

	// DEFAULT LEG
	if (its_DL_PaymentType == CEP_ATFIXEDDATE)
		TheDFAtAFixedDate	=	its_DF[(int)its_DL_PaymentDate];	// no lag
}


// ----------------------------------------------------
// Get Max relevant Credit Date
// ----------------------------------------------------

void	CreditManager::ComputeMaxCreditDate()
{
	int	i;
	RelativeDate	CurrentDate, MaxDate;
	
	MaxDate	=	0.0;
	
	// I have to loop on Premium Leg (not sorted) and take into account the last Payment Date
	for (i=0;i<its_PL_NbFlows;i++)
	{
		CurrentDate	=	its_PL_CreditWindowUps[i]; // should check WindowLow is <= WindowHigh
		if (CurrentDate > MaxDate)
			MaxDate	=	CurrentDate;
	}

	MaxDate	=	FMAX(MaxDate, CurrentDate);

	its_MaxCreditDate	=	MaxDate;
}


// ----------------------------------------------------
// Get Max relevant Credit Date
// ----------------------------------------------------

void	CreditManager::SortAllCreditObservationDates()
{
	its_SortedCreditObservationDates.clear();
	
	merge(its_PL_Sorted_CreditWindowUps.begin(), its_PL_Sorted_CreditWindowUps.end(),
		its_PL_CreditWindowLows.begin(), its_PL_CreditWindowLows.end(),
		back_insert_iterator<DateVector>(its_SortedCreditObservationDates));

	// move consecutive duplicates past the end; store new end
	DateVector::iterator	the_end	=	its_SortedCreditObservationDates.end();
	DateVector::iterator	new_end	=	unique(its_SortedCreditObservationDates.begin(), the_end);
	
	// delete all elements past new_end
	its_SortedCreditObservationDates.erase(new_end, the_end);

	int	nb_sorted_dates	=	its_SortedCreditObservationDates.size();
	Tau_Item	DummyItem(0, 0.0, 0.0, 0);

	its_SortedCreditObservationDatesAndLosses.clear();
	
	// I want to find all items between Tmin and Tmax
	set<Tau_Item>::iterator iter;
	int	the_id;	
	int	i;

	for (i=0;i<nb_sorted_dates;i++)
	{		
		DummyItem.id			=	i;
		DummyItem.tau			=	MATHTIME(its_SortedCreditObservationDates[i]);		// in year fractions
		DummyItem.cumul_loss	=	0.0;
		DummyItem.cumul_nbdef	=	0;

		if (IsCDOSquare())
		{
			DummyItem.cumul_loss_CDOs.resize(its_NbRelevant_CDO_Underlyings);
			DummyItem.cumul_nbdef_CDOs.resize(its_NbRelevant_CDO_Underlyings);
			DummyItem.tranche_losses_CDOs.resize(its_NbRelevant_CDO_Underlyings);

			DummyItem.cumul_loss_square		=	0.0;
			DummyItem.cumul_nbdef_square	=	0;
		}

		its_SortedCreditObservationDatesAndLosses.insert(DummyItem);
	}

	// iterator in order to skip the redondant items
	set<Tau_Item>::iterator	iter_found;

	its_CreditObervationDatesSortedLowsIds.clear();
	its_CreditObervationDatesSortedUpsIds.clear();

	for (i=0;i<its_PL_CreditWindowLows.size();i++)
	{
		// find
		DummyItem.tau	=	MATHTIME(its_PL_CreditWindowLows[i]);

		iter_found	=	find(its_SortedCreditObservationDatesAndLosses.begin(), its_SortedCreditObservationDatesAndLosses.begin(), DummyItem);

		// sure to be in the list (otherwise, test iter_found != 
		the_id	=	iter_found->id;
		its_CreditObervationDatesSortedLowsIds.push_back(iter_found);

	}

	if (its_CreditObervationDatesSortedLowsIds.size() != its_PL_CreditWindowLows.size())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
            "Parameters:  Pb in SortAllCreditObservationDates with its_CreditObervationDatesSortedLowsIds");

	for (i=0;i<its_PL_CreditWindowUps.size();i++)
	{
		// find
		DummyItem.tau	=	MATHTIME(its_PL_CreditWindowUps[i]);

		iter_found	=	its_SortedCreditObservationDatesAndLosses.find(DummyItem);

		// sure to be in the list (otherwise, test iter_found != 
		the_id	=	iter_found->id;
		its_CreditObervationDatesSortedUpsIds.push_back(iter_found);
	}

	if (its_CreditObervationDatesSortedUpsIds.size() != its_PL_CreditWindowLows.size())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
            "Parameters:  Pb in SortAllCreditObservationDates with its_CreditObervationDatesSortedUpsIds");

}


// ----------------------------------------------------
// Check the Consistency
// ----------------------------------------------------
void	CreditManager::CheckPremiumLegSchedule()
{
	// TO DO
}

// ----------------------------------------------------
// Get the DF at the Next Credit Date (Window Up)
// ----------------------------------------------------
void	CreditManager::Get_DFAtNextCreditDate(RelativeDate	TheDate, double& NextDate, double& TheDF)
{
	int	i;

	// Sorted
	NextDate	=	-1.0;
	TheDF		=	-1.0;

	for (i=0;i<its_PL_NbFlows;i++)
	{
		NextDate = its_PL_Sorted_CreditWindowUps[i];
		if (NextDate >= TheDate)
		{
			// Ok
			TheDF	=	its_DF[(int)NextDate];
			break;
		}
	}
	// ------------------------------------------------------------
	// this test should be avoided, according to Schedule Chekings
	if (TheDF == -1.0)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
            "Parameters:  Can't find next Credit Date");
	// ------------------------------------------------------------
}


void	CreditManager::Get_DFAtNextPaymentDate(RelativeDate TheDate, double& NextDate, double& TheDF)
{
	int	i;

	// Sorted
	NextDate	=	-1.0;
	TheDF		=	-1.0;

	for (i=0;i<its_PL_NbFlows;i++)
	{
		NextDate = its_PL_Sorted_PaymentDates[i];
		if (NextDate >= TheDate)
		{
			// Ok
			TheDF	=	its_DF[(int)NextDate];	
			break;
		}
	}

	// ------------------------------------------------------------
	// this test should be avoided, according to Schedule Chekings
	if (TheDF == -1.0)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
            "Parameters:  Can't find next Payment Date");
	// ------------------------------------------------------------

}


// ----------------------------------------------------
// Generate All Default Curves and Multi Curves Model
// ----------------------------------------------------

void	CreditManager::GenerateAllDefaultCurves()
{
	if ((KeepCalibration()) && (its_ModelMultiCurves != NULL)) return;

	int i;
	ARM_Vector*			Spread		=	NULL;
	ICM_DefaultCurve*	DefCurve	=	NULL;
	ICM_DefaultCurve*	DefCurve_Clone	=	NULL;

	// ------------------------------------------------------------
	// Default Curves
	its_ArrayDefaultCrv = new ICM_DefaultCurve*[its_CreditsLabels.size()];

	ComputeBasketNotional();
	// ------------------------------------------------------------
	// Here, according to the Spread -999., I should adjust the maturities

	// Roll Dates
	qCDS_ADJ	CDS_Roll_Type;
	CDS_Roll_Type	=	(its_RollDateFlag) ? qCredit_Adjust20 : qCredit_Default;

	long PayFreq = K_QUARTERLY;

	// bool	BrentFlag	=	true;
	qDEFCURVE_CALIB_ALGO calibAlgo= qDEFCURVE_BRENT ; 


	std::vector<const ICM_DefaultCurve*> items(its_CreditsLabels.size()); 
	// I have to get.. The Curve
	for (i=0;i<its_CreditsLabels.size();i++)
	{
		// in order to get an 'ARM_Vector'
		Spread = its_CreditDataSpreads->TruncRowAsVector(i,its_Maturities.size()); 
		ARM_Security* psec = GetSecurity();
		string ccy=("");
		if ( psec == NULL) ccy = ARM_DEFAULT_COUNTRY;
		else ccy = string(GetSecurity()->GetCurrencyUnit()->GetCcyName());
		// Test Roll Dates ?
		DefCurve = new ICM_Constant_Piecewise(
												itsValDate,
												its_Maturities,
												Spread,
												its_Recoveries[i],
												its_ZeroCurve,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
												CDS_Roll_Type,					// / Roll Dates
												ccy/*ARM_DEFAULT_COUNTRY*/,			//	ccy
												its_CreditsLabels[i],		//	label
												false,							//	issummitcurve
												NULL,
												//2 NULL,
												PayFreq,
												calibAlgo,"STD",
												ARM_Currency(ccy.c_str()).GetCreditStartDateLag(),
												ICM_Parameters());					


		// ARM_OK = 0 ???
		items[i]=its_ArrayDefaultCrv[i] = DefCurve;

		// Delete
		if (Spread)
			delete Spread;
		Spread = NULL;

		if (IsHomogeneousBasket())
		{
			i++;
			while (i<its_CreditsLabels.size())
			{
				DefCurve_Clone	=	(ICM_DefaultCurve*) (DefCurve->Clone());
				DefCurve_Clone->SetLabel(its_CreditsLabels[i]);

				items[i]=its_ArrayDefaultCrv[i] = DefCurve_Clone;

				i++;
			}
		}
	}

	ARM_Vector theLosses;

	// theLosses	=	(double *)malloc(its_NbCredits*sizeof(double));
	for (i=0;i<its_CreditsLabels.size();i++)
		theLosses[i]	=	its_Input_Losses[i] * its_Notionals[i];

	// Create MultiCurveModel
	if (its_ModelMultiCurves != NULL)
		delete its_ModelMultiCurves;

	its_ModelMultiCurves = new ICM_ModelMultiCurves(
											// its_CreditsLabels.size(),
											items ,
											its_ZeroCurve, 
											theLosses,
											itsCorrelation);

	// if (theLosses != NULL)
		// free(theLosses);
}


// ----------------------------------------------------
// BASKET USEFUL DATA
// ----------------------------------------------------

void	CreditManager::ComputeBasketNotional()
{
	int	i;

	its_LossesAmount.clear();
	its_LossesAmount.resize(its_CreditsLabels.size());

	its_BasketNotional	=	0.0;
	its_TotalLoss		=	0.0;
	
	double	CurrentLossAmount;
	double	CurrentLoss, CurrentNotional;

	double	FirstLossAmount;
	FirstLossAmount	=	its_Notionals[0];
	FirstLossAmount	*=	its_Input_Losses[0];
	
	its_IsHomogeneousBasket	=	true;
	its_IsHomogeneousBasketLossUnit	=	true;

	//Estimation de la vraie taille du basket (en fonction de la taille
	//de la tranche)
	double rapport = 0., first_total=0.;
	for (i=0;i<its_CreditsLabels.size();i++)
		first_total += MAX(its_Notionals[i],0.);

	//Renormalisation
	if (!CHECK_NULL(first_total * (its_PL_LossMaxs[0] - its_PL_LossMins[0])))
		rapport = its_PL_Notios[0] / (first_total * (its_PL_LossMaxs[0] - its_PL_LossMins[0]));
	
	for (i=0;i<its_CreditsLabels.size();i++)
		its_Notionals[i] = its_Notionals[i] * rapport;

	for (i=0;i<its_CreditsLabels.size();i++)
	{
		CurrentNotional		=	its_Notionals[i];
		CurrentLoss			=	its_Input_Losses[i];
		
		its_BasketNotional	+=	MAX(CurrentNotional,0.);
		CurrentLossAmount	=	CurrentLoss * CurrentNotional;

		// LOSS UNIT COMPUTATION
		if ((its_CreditModelType == CMT_ANALYTIC_RECURSIVE_1F) && (its_Recursive_1F_Loss_Unit_Choice == COFRLUC_LOSS_UNIT_MIN))
		{
			// just hope there will not have any side effect
			if (!CHECK_EQUAL(round(CurrentLoss/LOSS_UNIT_MIN),round(CurrentLoss)/LOSS_UNIT_MIN))
				CurrentLossAmount = round(CurrentLoss/LOSS_UNIT_MIN)*LOSS_UNIT_MIN;		
		}

		its_LossesAmount[i]	=	CurrentLossAmount;
		its_TotalLoss		+=	CurrentLossAmount;

		// Work on Loss Amount rather than on Notionals
		if ((its_IsHomogeneousBasket) && (CurrentLossAmount != FirstLossAmount))
		{
			its_IsHomogeneousBasket	=	false;
			its_IsHomogeneousBasketLossUnit	=	false;
		}
	}

	TheBasketNotional	=	its_BasketNotional;

	// LOSS UNIT COMPUTATION
	if ((its_CreditModelType == CMT_ANALYTIC_RECURSIVE_1F) || (its_CreditModelType == CMT_ANALYTIC_RECURSIVE_INTERPOLATION_1F))
	{
		switch (its_Recursive_1F_Loss_Unit_Choice)
		{
		case COFRLUC_PGCD:
		case COFRLUC_LOSS_UNIT_MIN:

			if (its_CreditsLabels.size() > 1)
				its_LossUnit	=	fabs(pgcd::solve(its_LossesAmount));
			else
				its_LossUnit	=	fabs(its_LossesAmount[0]);

			its_LossRate.clear();
			its_LossRate.resize(its_CreditsLabels.size());

			for (i=0;i<its_CreditsLabels.size();i++)
				its_LossRate[i]	=	floor(its_LossesAmount[i] / its_LossUnit + 0.5);

			its_Recursive_NbLossSteps	=	its_TotalLoss	/	its_LossUnit;
			
			break;

		case COFRLUC_NB_LOSS_STEP:
			its_Recursive_NbLossSteps	=	its_Recursive_1F_NbLossStep;
			break;
		}
	}

	// Check with Spreads
	if (its_IsHomogeneousBasket)
		its_IsHomogeneousBasket	=	AreSpreadsHomogeneous();
}


bool	CreditManager::AreSpreadsHomogeneous()
{
	int	i, j;
	const ARM_Vector*	Spread_Init	=	NULL;
	const ARM_Vector*	Spread		=	NULL;

	const ICM_DefaultCurve*	TheDefaultCurve=NULL;

	// in order to get an 'ARM_Vector'
	if (its_CreditDataSpreads)
	{
		Spread_Init = its_CreditDataSpreads->TruncRowAsVector(0, its_Maturities.size());

		for (i=1;i<its_CreditsLabels.size();i++)
		{
			// in order to get an 'ARM_Vector'
			Spread = its_CreditDataSpreads->TruncRowAsVector(i, its_Maturities.size());

			for (j=0; j<its_Maturities.size(); j++)
			{
				if ((*Spread)[j] != (*Spread_Init)[j])
					return false;
			}
		}
	}
	else
	{
		// try with Model Multi Curves
		if (its_ModelMultiCurves)
		{
			TheDefaultCurve	=	its_ModelMultiCurves->GetDefaultCurve(0);
			Spread_Init		=	TheDefaultCurve->GetRates();

			for (i=1;i<its_CreditsLabels.size();i++)
			{
				// in order to get an 'ARM_Vector'
				TheDefaultCurve	=	its_ModelMultiCurves->GetDefaultCurve(i);
				Spread			= TheDefaultCurve->GetRates();

				// recall that at time t=0 a copy has been made
				for (j=1; j<its_Maturities.size(); j++)
				{
					if ((*Spread)[j] != (*Spread_Init)[j])
						return false;
				}
			}
		}
		else
			return false;	// enable to say
	}
	
	return true;
}

// ----------------------------------------------------
// BASKET USEFUL DATA
// ----------------------------------------------------
/*
void	CreditManager::ComputeBasketNotional_CDO_Square()
{
	int	i;

	its_LossesAmount.clear();
	its_LossesAmount.resize(its_NbCredits);

	its_BasketNotional	=	0.0;
	its_TotalLoss		=	0.0;
	
	double	CurrentLossAmount;
	double	CurrentLoss, CurrentNotional;

	double	FirstLossAmount;
	FirstLossAmount	=	its_Notionals[0];
	FirstLossAmount	*=	its_Input_Losses[0];
	
	its_IsHomogeneousBasket	=	true;

	for (i=0;i<its_NbCredits;i++)
	{
		CurrentNotional		=	its_Notionals[i];
		CurrentLoss			=	its_Input_Losses[i];
		
		its_BasketNotional	+=	CurrentNotional;
		CurrentLossAmount	=	CurrentLoss * CurrentNotional;

		// LOSS UNIT COMPUTATION
		if ((its_CreditModelType == CMT_ANALYTIC_RECURSIVE_1F) && (its_Recursive_1F_Loss_Unit_Choice == COFRLUC_LOSS_UNIT_MIN))
		{
			// just hope there will not have any side effect
			if (!CHECK_EQUAL(round(CurrentLoss/LOSS_UNIT_MIN),round(CurrentLoss)/LOSS_UNIT_MIN))
				CurrentLossAmount = round(CurrentLoss/LOSS_UNIT_MIN)*LOSS_UNIT_MIN;		
		}

		its_LossesAmount[i]	=	CurrentLossAmount;
		its_TotalLoss		+=	CurrentLossAmount;

		// Work on Loss Amount rather than on Notionals
		if ((its_IsHomogeneousBasket) && (CurrentLossAmount != FirstLossAmount))
			its_IsHomogeneousBasket	=	false;
	}

	// LOSS UNIT COMPUTATION
	if ((its_CreditModelType == CMT_ANALYTIC_RECURSIVE_1F) || (its_CreditModelType == CMT_ANALYTIC_RECURSIVE_INTERPOLATION_1F))
	{
		switch (its_Recursive_1F_Loss_Unit_Choice)
		{
		case COFRLUC_PGCD:
		case COFRLUC_LOSS_UNIT_MIN:

			if (its_NbCredits > 1)
				its_LossUnit	=	fabs(pgcd::solve(its_LossesAmount));
			else
				its_LossUnit	=	fabs(its_LossesAmount[0]);

			its_LossRate.clear();
			its_LossRate.resize(its_NbCredits);

			for (i=0;i<its_NbCredits;i++)
				its_LossRate[i]	=	its_LossesAmount[i] / its_LossUnit;

			its_Recursive_NbLossSteps	=	its_TotalLoss	/	its_LossUnit;
			
			break;

		case COFRLUC_NB_LOSS_STEP:
			its_Recursive_NbLossSteps	=	its_Recursive_1F_NbLossStep;
			break;
		}
	}

}
*/

// ----------------------------------------------------
// Price
// ----------------------------------------------------

void	CreditManager::PriceOrHedge()
{
	its_temp_index	=	0;

	// no shift
	DeActivateShift();

	Price();

	// ---------------
	if (its_HedgesRunning)
	{
		// Keep Central NPV
		Central_NPV	=	NPV;
		ComputeHedges();
		return;
	}

}

void	CreditManager::Price()
{
	DoubleVector	Outputs;

	// -------------------------------------------------------
	// Prepare Outputs
	// -------------------------------------------------------
	PrepareOutputs(Outputs);

	// -------------------------------------------------------
	// Schedule Checkings
	// -------------------------------------------------------
	ScheduleCheckings();

	// -------------------------------------------------------
	// Compute Useful DF Data
	// -------------------------------------------------------
	ComputeAllDF();

	// -------------------------------------------------------
	// Compute Useful Data for the Basket
	// -------------------------------------------------------
	ComputeBasketNotional();

	// -------------------------------------------------------
	// CDO Square Specific Data
	// -------------------------------------------------------
	if (IsCDOSquare())
	{
		if (its_CreditModelType != CMT_MONTECARLO)
			ICMTHROW(ERR_INVALID_DATA,"CDO Square only available with MONTE-CARLO Numerical Method!");

		PrepareCDOSquareData();
//		DisplayCDOSquareData();

		TheBasketNotional	=	itsCDOSquareNotional;
	}
	else
		TheBasketNotional	=	its_BasketNotional;

	// -------------------------------------------------------
	// Compute Max Nb Def and Loss Levels
	// -------------------------------------------------------
	ComputeMaxLevels();

	// -------------------------------------------------------
	// Prepare for Pricing
	// -------------------------------------------------------
	SortAllCreditObservationDates();

	// -------------------------------------------------------
	// All Credit Curves have already been generated
	// -------------------------------------------------------
//	GenerateAllDefaultCurves();

	// -------------------------------------------------------
	AllocateStructures_1F();

	if (its_CreditModelType == CMT_MONTECARLO)
		Price_MonteCarlo(Outputs);
	else
		Price_Analytic_1F(Outputs);
	// -------------------------------------------------------

	// -------------------------------------------------------
	// Get Cash and Accrued 
	Compute_Cash_And_Accrued();

	SetOutputData();
	// -------------------------------------------------------


}



void	CreditManager::Price_MonteCarlo(DoubleVector& Outputs)
{
		// ----------------------------------------------------------------------------
//		if ((its_fOut = fopen("c:\\test\\mc_Price_Monte_Carlo.txt", "a+")) == NULL) return;
//		fprintf(its_fOut, " ----------------- NPVs ----------------- \n");
		// ----------------------------------------------------------------------------
	// -------------------------------------------------------
	// Compute Barriers for Default Times
	// -------------------------------------------------------
	ComputeMaxCreditDate();
	ComputeBarriers(its_MaxCreditDate);

	// -------------------------------------------------------
	// Variance Reduction useful data computation
	// -------------------------------------------------------
	ComputeVarianceReductionData();

	// -------------------------------------------------------
	// Initialization of the Random Generator
	// -------------------------------------------------------
	SetRandomGenerator();

	// -------------------------------------------------------
	// According to Correlation
	// -------------------------------------------------------
	switch (its_CorrelationType)
	{
		case CT_FLAT:
		case CT_BETA:	
			break;

		case CT_MATRIX:
			// Cholesky computation
			CholeskyComputation();

			break;

		case CT_FACTOR_LOADING_2:
			break;
	}

//	tmpPrem.clear();
//	tmpPrem.resize(its_PL_NbFlows);

	// -------------------------------------------------------
	// The Pricing Loop  
	// -------------------------------------------------------
	for (its_SimulId=0;its_SimulId<its_NbSimul;its_SimulId++)
	{
		// -------------------------------------------------------
		// According to Correlation
		// -------------------------------------------------------
		switch (its_CorrelationType)
		{
			case CT_FLAT:
			case CT_BETA:
			case CT_FACTOR_LOADING_2:
				
				GenerateDefaultTimes_Betas();
				break;
			case CT_MATRIX:
				GenerateDefaultTimes_Matrix();		// Cholevsky
				break;
			
		}

		// -------------------------------------------------------
		// DISPLAY
		// -------------------------------------------------------
//		Display_DefaultTimes();		// and shifted if there any

		// -------------------------------------------------------
		// ComputeCumulativeLossesAndNbDefaults
		// -------------------------------------------------------
		if (IsCDOStandard())
			ComputeCumulativeLossesAndNbDefaults();
		else // if (IsCDOSquare())
			ComputeCumulativeLossesAndNbDefaults_CDO_SQUARE();

		// -------------------------------------------------------
		// Price Basket For This Simul
		// -------------------------------------------------------
		PriceBasketForThisSimulation(Outputs);

		NPV				+=	Outputs[0];
		DefaultLegPV	+=	Outputs[1];
		PremiumLegPV	+=	Outputs[2];
		ATMPremiumLeg	+=	Outputs[3];
		ATMPremiumLegWithoutNotio		+=	Outputs[4];
		NPVSquare		+=	Outputs[0] * Outputs[0];

//		fprintf(its_fOut, "%u\t\t%.4lf\n", its_SimulId, Outputs[0]);

	}

	// ------------------------------------
	// VARIANCE REDUCTION
	// ------------------------------------
	double	TheValue;
	double	TheSquareValue;

	if (its_CreditObservationType == CO_LOSSES)
	{
		TheValue	=	1.0;					
		TheSquareValue	=	1.0;
			
		switch (its_MC_Variance_Reduction)
		{					
		case CMCVR_IS_PURE_FACTORS:
		case CMCVR_IS_FACTORS:
			TheValue		=	IS_Common_FactorShift;
			TheSquareValue	=	TheValue * TheValue;

			NPV				*=	TheValue;
			DefaultLegPV	*=	TheValue;
			PremiumLegPV	*=	TheValue;
			ATMPremiumLeg	*=	TheValue;
			ATMPremiumLegWithoutNotio		*=	TheValue;
			NPVSquare		*=	TheSquareValue;
			
		case CMCVR_IS_IDIOSYNCRATIC:

			break;

		case CMCVR_NONE:
		default:
		
			break;
		}
	}

//	fclose(its_fOut);

//	for (int i=0; i<its_PL_NbFlows; i++)
//		tmpPrem[i] /= its_NbSimul;

//	string	tmpLabel("All Premium Flows");
	
//	Display_Vector(tmpLabel, tmpPrem);
}


void	CreditManager::Price_Analytic_1F(DoubleVector& Outputs)
{
	// According to the Model Type
	// I may test
	its_Loss_Distrib_Map.clear();

	// Pricing
	PriceBasket_1F(Outputs);

	// Outputs
	NPV				=	Outputs[0];
	DefaultLegPV	=	Outputs[1];
	PremiumLegPV	=	Outputs[2];
	ATMPremiumLeg	=	Outputs[3];
	ATMPremiumLegWithoutNotio		=	Outputs[4];
	NPVSquare		=	0.0;
}


void	CreditManager::SetRandomGenerator()
{
	// TO BE REVIEWED, IN ORDER TO TAKE INTO ACCOUNT GAUSSIAN OR STUDENT
	// Must depend on the Correlation Type
	// Remark, TheSeed, must depend on the Uniform Random Generator.
	
	int i;

	int	TheSeed	=	0;
	vector<double> betas;
	betas.clear();
	vector<double> betas_2;
	betas_2.clear();
	vector<double> alphas;
	alphas.clear();

	switch (its_CorrelationType)
	{
	case CT_FLAT:

		if (fabs(its_CorrelationValue - 1.0) <= DB_TOL)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
				"ERROR: Bad Correlation Value must lie in ]-1.0 ; 1.0[");

		betas.resize(its_CreditsLabels.size(), sqrt(its_CorrelationValue));

		// A Gaussian Variable
		// seed - nb Factors - Beta vector
		switch (its_CopulaType)
		{
			case CCT_GAUSSIAN:
				SetGenerator((ICM_Root_Generator *) new ICM_Gengauss1F(TheSeed, its_CreditsLabels.size(), betas));
			break;

			case CCT_STUDENT:
				SetGenerator((ICM_Root_Generator *) new ICM_Gen_Student1F(TheSeed, its_CreditsLabels.size(), its_FreedomDegree, betas));
				break;
		}

		break;

	case CT_BETA:

		for (i=0;i<its_CreditsLabels.size();i++)
			betas.push_back(its_Beta[i]);

			// A Gaussian Variable
			// seed - nb Factors - Beta vector
		switch (its_CopulaType)
		{
			case CCT_GAUSSIAN:
				SetGenerator((ICM_Root_Generator *) new ICM_Gengauss1F(TheSeed, its_CreditsLabels.size(), betas));
			break;

			case CCT_STUDENT:
				SetGenerator((ICM_Root_Generator *) new ICM_Gen_Student1F(TheSeed, its_CreditsLabels.size(), its_FreedomDegree, betas));
				break;
		}
		break;

	case CT_MATRIX:
		for (i=0;i<its_CreditsLabels.size();i++)
			betas.push_back(0.0);

			// A Gaussian Variable
			// seed - nb Factors - Beta vector
		switch (its_CopulaType)
		{
			case CCT_GAUSSIAN:
				SetGenerator((ICM_Root_Generator *) new ICM_Gengauss1F(TheSeed, its_CreditsLabels.size(), betas));
			break;

			case CCT_STUDENT:
				SetGenerator((ICM_Root_Generator *) new ICM_Gen_Student1F(TheSeed, its_CreditsLabels.size(), its_FreedomDegree, betas));
				break;
		}
		break;


	case CT_FACTOR_LOADING_2:

			// A Gaussian Variable
		for (i=0;i<its_CreditsLabels.size();i++)
		{
			betas.push_back(its_FL_Beta1[i]);
			betas_2.push_back(its_FL_Beta2[i]);
			alphas.push_back(its_FL_Alpha[i]);
		}

			// seed - nb Factors - Beta vector
		switch (its_CopulaType)
		{
			case CCT_GAUSSIAN:
				SetGenerator((ICM_Root_Generator *) new ICM_Gengauss1F_FactorLoading_2Points(TheSeed, its_CreditsLabels.size(), betas, betas_2, alphas));
			break;

			case CCT_STUDENT:
				ICMTHROW(ERR_INVALID_DATA,"Factor Loading with Student Law not yet implemented!");
				break;
		}
		break;


	// take into account CorrelationObject
	default:
		ICMTHROW(ERR_INVALID_DATA,"Bad Correlation Type for Random Generation!");
	}

	// --------------------------------------
	// VARIANCE REDUCTION
	// --------------------------------------
	double	TheDrift;
	ICM_Root_Generator* Gen = GetGenerator();
	
	switch (its_MC_Variance_Reduction)
	{
		case CMCVR_IS_FACTORS:
		case CMCVR_IS_PURE_FACTORS:
			TheDrift	=	its_IS_Mu;
			IS_Common_FactorShift	=	exp(0.5 * its_IS_Mu * its_IS_Mu);
			break;

		case CMCVR_IS_IDIOSYNCRATIC:
		case CMCVR_NONE:
		default:
			TheDrift	=	0.0;
			IS_Common_FactorShift	=	1.0;
			break;
	}

	if (Gen)
		Gen->setCommonFactor_Drift(TheDrift);
	// --------------------------------------

	// --------------------------------------
	betas.clear();
	betas_2.clear();
	alphas.clear();
	// --------------------------------------

}



// ----------------------------------------------------------------------------
// Compute Barriers in complete case
// ----------------------------------------------------------------------------
void CreditManager::ComputeBarriers(double	TimeT)
{
	int i;
	double DefProb;
	double	TheValue;
//	double	TheCoef;

	itsDef_Prob.resize(its_CreditsLabels.size());
	itsBarriers_Standard.resize(its_CreditsLabels.size());
	
	// Default Proba expects only yearterms
/*	for (i=0;i<its_NbCredits;i++)
		DefProbs[i] = (its_ArrayDefaultCrv[i])->DefaultProba(MATHTIME(its_MaxCreditDate));

	switch (its_MC_Variance_Reduction)
	{
		case CMCVR_IS_IDIOSYNCRATIC:
		case CMCVR_IS_FACTORS:

			for (i=0;i<its_NbCredits;i++)
			{
				DefProb	=	DefProbs[i];
				TheCoef	=	IS_individual_exp_twists[i];

				TheValue	=	DefProb * TheCoef / (1.0 + DefProb * (TheCoef - 1.0));

				DefProbs[i]	=	TheValue;
			}

		case CMCVR_NONE:
		case CMCVR_IS_PURE_FACTORS:
		default:

			break;
	}
*/

	switch (its_CopulaType)
	{
		case CCT_GAUSSIAN:

			if (IsHomogeneousBasket())
			{
				itsDef_Prob[0] = (its_ArrayDefaultCrv[0])->DefaultProba(MATHTIME(TimeT));	// in year fraction MATHTIME(its_MaxCreditDate));
				DefProb = itsDef_Prob[0];

				// should test if DefProb = 0.0 or = 1.0
				if (DefProb == 0.0)
					TheValue	=	_MINUS_INFINITY_;	//	minus infinity
				else if (DefProb == 1.0)
					TheValue	=	_PLUS_INFINITY_;	//	plus infinity
				else
				{
					TheValue	=	its_CopulaChoice->Inverse_Cumulative_Density_Function(DefProb);
//					TheValue	=	g01fac(Nag_LowerTail, DefProb, NAGERR_DEFAULT);
				}

				for (i=0;i<its_CreditsLabels.size();i++)
					itsBarriers_Standard[i] = TheValue;

			}
			else
				for (i=0;i<its_CreditsLabels.size();i++)
				{
					itsDef_Prob[i] = (its_ArrayDefaultCrv[i])->DefaultProba(MATHTIME(TimeT));	// in year fraction MATHTIME(its_MaxCreditDate));
					DefProb = itsDef_Prob[i];

					// should test if DefProb = 0.0 or = 1.0
					if (DefProb == 0.0)
						TheValue	=	_MINUS_INFINITY_;	//	minus infinity
					else if (DefProb == 1.0)
						TheValue	=	_PLUS_INFINITY_;	//	plus infinity
					else
					{
						TheValue	=	its_CopulaChoice->Inverse_Cumulative_Density_Function(DefProb);
	//					TheValue	=	g01fac(Nag_LowerTail, DefProb, NAGERR_DEFAULT);
					}

					itsBarriers_Standard[i] = TheValue;
				}
		break;

		case CCT_STUDENT:

			ICMTHROW(ERR_INVALID_DATA,"Student Copula not yet implemented for 1F Model!");
		break;

		default:
			ICMTHROW(ERR_INVALID_DATA,"Unknown Copula Type!");
		break;
	}

}


// ----------------------------------------------------------------------------
// Compute Barrier for a Shifted Curve in complete case
// ----------------------------------------------------------------------------
void CreditManager::ComputeBarriersShifted(double TimeT, int col_id)
{
	// col_id represents the hedges' scenario id

	int i;
	double DefProb;
	double	TheValue;

	itsBarriers_Shifted.resize(its_CreditsLabels.size());

	if (IsHomogeneousBasket())
	{
		DefProb = ((*its_MatrixShiftedDefaultCrv)(0, col_id))->DefaultProba(MATHTIME(TimeT)); // year fraction MATHTIME (its_MaxCreditDate));	// Default Proba expects only yearterms

		// should test if DefProb = 0.0 or = 1.0
		if (DefProb == 0.0)
			TheValue	=	_MINUS_INFINITY_;	//	minus infinity
		else if (DefProb == 1.0)
			TheValue	=	_PLUS_INFINITY_;	//	plus infinity
		else
		{
			TheValue	=	its_CopulaChoice->Inverse_Cumulative_Density_Function(DefProb);
			TheValue=NAG_deviates_normal(DefProb); 
		}
		
		for (i=0;i<its_CreditsLabels.size();i++)
			itsBarriers_Shifted[i] = TheValue;

	}
	else
		for (i=0;i<its_CreditsLabels.size();i++)
		{
			DefProb = ((*its_MatrixShiftedDefaultCrv)(i, col_id))->DefaultProba(MATHTIME(TimeT)); // year fraction MATHTIME (its_MaxCreditDate));	// Default Proba expects only yearterms

			// should test if DefProb = 0.0 or = 1.0
			if (DefProb == 0.0)
				TheValue	=	_MINUS_INFINITY_;	//	minus infinity
			else if (DefProb == 1.0)
				TheValue	=	_PLUS_INFINITY_;	//	plus infinity
			else
			{
				TheValue	=	its_CopulaChoice->Inverse_Cumulative_Density_Function(DefProb);
				TheValue=NAG_deviates_normal(DefProb); 
			}

			itsBarriers_Shifted[i] = TheValue;
		}

}

// ----------------------------------------------------------------------------
// Generation of default times with Betas
// ----------------------------------------------------------------------------
void CreditManager::GenerateDefaultTimes_Betas()
{
	itsSortedDefaultTimes.clear();
	its_DefaultTimes.clear();
	its_DefaultTimesShifted.clear();

	ICM_Root_Generator* Gen = GetGenerator();
	Gen->CommonFactors();
	
	double	bari;
	double aux = 0.;
	double ri, t;
	unsigned int i;

	DoubleVector	All_ri;
	All_ri.resize(its_CreditsLabels.size());

	int j;
	double	TheExpectedLoss_z;
	double	TheCommonFactor;
	
	TheCond_DefProb_Vect.clear();
	TheCond_DefProb_Vect.resize(its_CreditsLabels.size());

	double	Theta, Psi_Prime, Psi_Second;
	double	NextTheta;
		
	double	rj, bar_j, den;
	double	exp_twist_j;
	double	cond_def_prob_j;
	double	t_old, new_aux;

	DoubleVector	Cond_Def_Prob;

	// -------------------------------------------------------
	// VARIANCE REDUCTION
	// -------------------------------------------------------
//	if (its_CreditModelType == CMT_MONTECARLO)		// always MC in this function
//	{
		// ----------------------------------------------------------------------------
//		if ((its_fOut = fopen("c:\\test\\mc_Gaussian_Random.txt", "a+")) == NULL) return;
//		fprintf(its_fOut, " ----------------- Theta ----------------- \n");
		// ----------------------------------------------------------------------------

		// required in order to Compute Barriers
		IS_Individual_FactorShift	=	1.0;
		TheCommonFactor	=	Gen->GetCommonFactor();

		// Common Factor
//		fprintf(its_fOut, "\n%u\t%.6lf", its_SimulId, TheCommonFactor);

		switch (its_MC_Variance_Reduction)
		{
			case CMCVR_IS_FACTORS:
				IS_Individual_FactorShift	=	exp(- its_IS_Mu * TheCommonFactor);

			case CMCVR_IS_IDIOSYNCRATIC:

				its_IS_Theta	=	its_Kept_Theta;

				// FIRST STEP - computes Sum pk(Z) * ck
				TheExpectedLoss_z	=	0.0;
				for (j=0; j<its_CreditsLabels.size(); j++)
				{
					bar_j	=	itsBarriers_Standard[j];
					cond_def_prob_j				=	Gen->getConditionalDefaultProbability(bar_j, j);
					// keep it for later loops
					TheCond_DefProb_Vect[j]	=	cond_def_prob_j;

					TheExpectedLoss_z	+=	cond_def_prob_j * its_LossesAmount[j]; // / TheBasketNotional;
				}
				
				i=0;
				if (TheExpectedLoss_z <= its_IS_Loss_Level) // / TheBasketNotional)
				{					
					switch (its_Theta_Choice)
					{
					case CVRI_OPTIM:

						Theta		=	1.0;
						NextTheta	=	0.0;

						// NEWTON method in order to get a guess of Theta

						while ((i<NB_ITER_NEWTON) && (FABS(NextTheta / Theta - 1.0) > 1e-6))
						{
							Theta		=	NextTheta;
							Psi_Prime	=	Psi_First_Deriv(Theta) - its_IS_Loss_Level / TheBasketNotional;	// CDO!
							Psi_Second	=	Psi_Second_Deriv(Theta);
							NextTheta	=	Theta - Psi_Prime / Psi_Second;
							i++;
							if (NextTheta < 0.0)
								NextTheta	=	i*10.0;	// ???
						}

						if ((i == NB_ITER_NEWTON) && (FABS(NextTheta / Theta - 1.0) > 1e-6))
							its_IS_Theta	=	its_Kept_Theta;
						else
							its_IS_Theta	=	Theta;

						for (j=0; j<its_CreditsLabels.size(); j++)
						{
							exp_twist_j		=	exp(its_IS_Theta * its_LossesAmount[j] / TheBasketNotional);
							IS_individual_exp_twists[j]			=	exp_twist_j;
							IS_individual_exp_twists_Kept[j]	=	exp_twist_j;
						}

					case CVRI_IMPOSED:

						IS_Common_exp_twists	=	1.0;
						for (j=0; j<its_CreditsLabels.size(); j++)
						{
//							exp_twist_j		=	exp(its_IS_Theta * its_LossesAmount[j] / TheBasketNotional);
//							IS_individual_exp_twists[j]	=	exp_twist_j;
							// IMPOSED CASE --> already computed, Kept
							exp_twist_j	=	IS_individual_exp_twists_Kept[j];

							cond_def_prob_j	=	TheCond_DefProb_Vect[j];

							den		=	1.0 + cond_def_prob_j * (exp_twist_j - 1.0);
							IS_Common_exp_twists	*=	den;
						}

						break;
					}
				}
				else
				{
					its_IS_Theta	=	0.0;
					IS_Common_exp_twists	=	1.0;

					fill(IS_individual_exp_twists.begin(), IS_individual_exp_twists.end(), 1.0);
				}						
					
//				ICMLOG("For Simulation " << its_SimulId << " Theta is " << its_IS_Theta << " - in " << i << " iterations.");
						
				// Now compute New Default Times and IS_Common_exp_twists

				// Computes the Factor CGF
				// in order to have some indication about Theta

//				IS_Common_exp_twists	=	1.0;
				for (j=0; j<its_CreditsLabels.size(); j++)
				{
					exp_twist_j	=	IS_individual_exp_twists[j];
					// barrier
					bar_j	=	itsBarriers_Standard[j];
					rj		=	Gen->generateRandom(j);
//							cond_def_prob_j	=	TheCond_DefProb_Vect[j];

//					fprintf(its_fOut, "\t%.6lf", Gen->GetIdiosyncraticValue());
					
					// Shifted Default Times
					if (bar_j == _PLUS_INFINITY_)
					{
						t		=	0.0;
						itsSortedDefaultTimes.insert(Tau_Item(j,t));
					}
					else if (bar_j == _MINUS_INFINITY_)
					{
						t		=	_PLUS_INFINITY_;
					}
					else // if (rj < bar_j)
					{
						if (its_CorrelationType == CT_FACTOR_LOADING_2)
							aux	=	Gen->CumulativeDistribution(j, rj);
						else
							switch (its_CopulaType)
							{
							case CCT_GAUSSIAN:
								aux	=	NAG_cumul_normal(rj);
								break;
							case CCT_STUDENT:
								aux =	NAG_prob_students_t( rj, its_FreedomDegree);
								break;
							}

						// use aux in order to compute the adjusted exponential twist
						// transform aux = N(X(j)), with X(j) = beta(j) * Z + sqrt(1-beta(j)^2) * Z(j)
						// to be commented
						t_old = (its_ArrayDefaultCrv[j])->DefProbInverse(aux);

						new_aux	=	aux;
						new_aux	/=	exp_twist_j * (1.0 - aux) + aux;						

						t = (its_ArrayDefaultCrv[j])->DefProbInverse(new_aux);
						if (t <= MATHTIME(its_MaxCreditDate))
							itsSortedDefaultTimes.insert(Tau_Item(j,t));
					}						

					// The default time is always stored
					if (CHECK_EQUAL(t, 0.0))
						t	=	1e-6;
					its_DefaultTimes.push_back(t);
				}

//					fprintf(its_fOut, "Simul Id:\t%u\t\tFactor Value:\t%.6lf\t\tTheta:\t%.4lf\n", its_SimulId, TheCommonFactor, Theta);	
//					fclose(its_fOut);

//					IS_Common_exp_twistss.push_back(IS_Common_exp_twists);

				break;

			case CMCVR_IS_PURE_FACTORS:
				IS_Individual_FactorShift	=	exp(- its_IS_Mu * TheCommonFactor);

			case CMCVR_NONE:

				for (i=0;i<its_CreditsLabels.size();i++)
				{	
					ri = Gen->generateRandom(i);
					// stored for Variance Reduction
//					All_ri[i]	=	ri;

					bari	=	itsBarriers_Standard[i];
					if (bari == _PLUS_INFINITY_)
					{
						t	=	0.0;
						itsSortedDefaultTimes.insert(Tau_Item(i,t));
					}
					else if (bari == _MINUS_INFINITY_)
						t	=	_PLUS_INFINITY_;
					else if (ri < bari)
					{
						if (its_CorrelationType == CT_FACTOR_LOADING_2)
							aux	=	Gen->CumulativeDistribution(i, ri);
						else
							switch (its_CopulaType)
							{
							case CCT_GAUSSIAN:
								aux = NAG_cumul_normal(ri);
								break;
							case CCT_STUDENT:
								aux =	NAG_prob_students_t( ri, its_FreedomDegree);
								break;
							}

						t = (its_ArrayDefaultCrv[i])->DefProbInverse(aux);
						if (CHECK_EQUAL(t, 0.0))
							t	=	1e-6;		// some T_EPS ???

						itsSortedDefaultTimes.insert(Tau_Item(i,t));
					}
					else
						t	=	-1.0;		// if I want to keep track

					its_DefaultTimes.push_back(t);
				}

			default:

				break;
		}

//		fclose(its_fOut);

//	}

	// ----------------------------------------------------------------------------
	if (its_Display_Sorted_DefaultTimes_Flag)
	{
		set<Tau_Item>::iterator iter;

		// ----------------------------------------------------------------------------
		if ((its_fOut = fopen("c:\\test\\mc_GenerateDefaultTimes_Betas.txt", "a+")) == NULL) return;
		fprintf(its_fOut, " ----------------- GenerateDefaultTimes_Betas ----------------- \n");
		// ----------------------------------------------------------------------------
		for (iter = itsSortedDefaultTimes.begin();
					iter != itsSortedDefaultTimes.end();
					++iter)
		{
			fprintf(its_fOut, "Id:\t%u\t\ttau:\t%.2lf\t\tCumul loss:\t%.2lf\t\tCumul NbDef:\t%u\n", iter->id, iter->tau, iter->cumul_loss, iter->cumul_nbdef);
		}
		fclose(its_fOut);
	}

}	

double	CreditManager::Psi_First_Deriv(double theta)
{
	int j;
	double	TheValue;
	double	Sum	=	0.0;
	double	The_individual_exp_twist;


	for (j=0; j<its_CreditsLabels.size(); j++)
	{
		The_individual_exp_twist	=	exp(theta * its_LossesAmount[j] / TheBasketNotional);
		TheValue	=	TheCond_DefProb_Vect[j] * its_LossesAmount[j] / TheBasketNotional * The_individual_exp_twist;
		TheValue	/=	(1 + TheCond_DefProb_Vect[j] * (The_individual_exp_twist - 1.0));
		Sum +=	TheValue;
	}

	return Sum;
}


double	CreditManager::Psi_Second_Deriv(double theta)
{
	int j;
	double	TheValue;
	double	Sum	=	0.0;
	double	The_individual_exp_twist;

	for (j=0; j<its_CreditsLabels.size(); j++)
	{
		The_individual_exp_twist	=	exp(theta * its_LossesAmount[j] / TheBasketNotional);
		TheValue	=	TheCond_DefProb_Vect[j] * its_LossesAmount[j] / TheBasketNotional * its_LossesAmount[j] / TheBasketNotional * The_individual_exp_twist;
		TheValue	*=	(1 - TheCond_DefProb_Vect[j]);
		TheValue	/=	(1 + TheCond_DefProb_Vect[j] * (The_individual_exp_twist - 1.0));
		TheValue	/=	(1 + TheCond_DefProb_Vect[j] * (The_individual_exp_twist - 1.0));
		Sum +=	TheValue;
	}

	return Sum;
}


// ----------------------------------------------------------------------------
// Generation of default times with Betas
// ----------------------------------------------------------------------------
void CreditManager::GenerateDefaultTimes_Betas_ForHedge(int col_id, bool FastHedgeFlag)
{
	if (its_Bump_Choice == CHB_RECOVERY_LOSS)
	{
		GenerateDefaultTimes_Betas();
		return;
	}

	itsSortedDefaultTimes.clear();
	its_DefaultTimes.clear();
	its_DefaultTimesShifted.clear();

	ICM_Root_Generator* Gen = GetGenerator();
	Gen->CommonFactors();

	double aux = 0.;
	double ri, t, t_shift, bari;
	unsigned int i;


	for (i=0; i<its_CreditsLabels.size(); i++)
	{
		aux	=	-1.0;

		ri		=	Gen->generateRandom(i);
		bari	=	itsBarriers_Standard[i];

		if (bari == _PLUS_INFINITY_)
		{
			t	=	0.0;
			itsSortedDefaultTimes.insert(Tau_Item(i,t));
		}
		else if (bari == _MINUS_INFINITY_)
			t	=	_PLUS_INFINITY_;
		else if (ri < bari)
		{
			if (its_CorrelationType == CT_FACTOR_LOADING_2)
				aux	=	Gen->CumulativeDistribution(i, ri);
			else
				switch (its_CopulaType)
				{
				case CCT_GAUSSIAN:
					aux = NAG_cumul_normal(ri);
					break;
				case CCT_STUDENT:
					aux = NAG_prob_students_t( ri, its_FreedomDegree);
					break;
				}

			t = (its_ArrayDefaultCrv[i])->DefProbInverse(aux);
			itsSortedDefaultTimes.insert(Tau_Item(i,t));
		}
		else
			t	=	-1.0;		// if I want to keep track

		its_DefaultTimes.push_back(t);

		if (IsHedgesDefault())
		{
			// I already know, that the New Tau is 0.0 (Today)
			t_shift	=	0.0;
		}
		else
			if (ri < itsBarriers_Shifted[i])
			{
				if (aux == -1.0)
					if (its_CorrelationType == CT_FACTOR_LOADING_2)
						aux	=	Gen->CumulativeDistribution(i, ri);
					else
						switch (its_CopulaType)
						{
						case CCT_GAUSSIAN:
							aux = NAG_cumul_normal(ri);
							break;
						case CCT_STUDENT:
							aux = NAG_prob_students_t( ri, its_FreedomDegree );
							break;
						}

				t_shift = ((*its_MatrixShiftedDefaultCrv)(i, col_id))->DefProbInverse(aux);
			}
			else
				t_shift	=	-1.0;		// if I want to keep track

		its_DefaultTimesShifted.push_back(t_shift);
/*
		if (IsHomogeneousBasket())
		{
			for (i=1; i<its_NbCredits; i++)
			{
				its_DefaultTimes.push_back(t);
				if (itsSortedDefaultTimes.size() > 0)
					itsSortedDefaultTimes.insert(Tau_Item(i,t));
				its_DefaultTimesShifted.push_back(t_shift);
			}
		}
*/	}

	if (FastHedgeFlag)
	{
		double			TheValue;

		DoubleVector	TheBetas;
		DoubleVector	TheCours;
		double			TheCommonFactor;	// 1F !!!
		double			TheCour;
		double			DefProb;

		bari		=	itsBarriers_Standard[i];
		TheBetas	=	GetBetas();
		TheCours	=	GetComplementaryBetas();
		TheCommonFactor	=	Gen->GetCommonFactor();

		for (i=0; i<its_CreditsLabels.size(); i++)
		{
			TheCour	=	TheCours[i];

			// I should consider the case Correlation = 1.0
			if (CHECK_EQUAL(TheCour, 0.0))
				its_BarrierDerivatives_FastHedge[i]	=	0.0;
			else
			{
				DefProb	=	itsDef_Prob[i];
				
				aux	=	(bari - TheBetas[i] * TheCommonFactor) / TheCour;
				TheValue	=	StandardGaussianDensity(aux);

				TheValue	/=	TheCour;
				TheValue	*=	(1.0 - DefProb);

				aux	=	StandardGaussianDensity(bari);
//				if (CHECK_EQUAL(aux, 0.0))		// should not happen

				TheValue	/=	aux;

				// hypothesis Epsilon = 1.0
				its_BarrierDerivatives_FastHedge[i]	=	TheValue;
			}
		}
	}

/*
	set<Tau_Item>::iterator iter;

	// ----------------------------------------------------------------------------
	if ((its_fOut = fopen("c:\\test\\mc_GenerateDefaultTimes_Betas_ForHedge.txt", "a+")) == NULL) return;
	fprintf(its_fOut, " ----------------- GenerateDefaultTimes_Betas_ForHedge ----------------- \n");
	// ----------------------------------------------------------------------------
	for (iter = itsSortedDefaultTimes.begin();
				iter != itsSortedDefaultTimes.end();
				++iter)
	{
		fprintf(its_fOut, "Id:\t%u\t\ttau:\t%.2lf\t\tCumul loss:\t%.2lf\t\tCumul NbDef:\t%u\n", iter->id, iter->tau, iter->cumul_loss, iter->cumul_nbdef);
	}
	fclose(its_fOut);
*/	
}	


// ----------------------------------------------------------------------------
// Generation of default times with Matrix via Cholesky
// ----------------------------------------------------------------------------
void CreditManager::GenerateDefaultTimes_Matrix()
{
	itsSortedDefaultTimes.clear();

	ICM_Root_Generator* Gen = GetGenerator();
	Gen->CommonFactors();

	std::vector<double>	its_indep_normal_vector; 
	std::vector<double> its_corr_normal_vector; 
	
	its_indep_normal_vector.resize(its_CreditsLabels.size());
	its_corr_normal_vector.resize(its_CreditsLabels.size());

	double aux = 0.;
	double ri, t;
	unsigned int i,j;

	// Independant variables generation
	for (i=0;i<its_CreditsLabels.size();i++)
	{	
		its_corr_normal_vector[i]	=	0.;
		its_indep_normal_vector[i]	=	NAG_random_normal(0.,1.);
	}

	for (i=0;i<its_CreditsLabels.size();i++)
		for (j=0;j<its_CreditsLabels.size();j++)
			its_corr_normal_vector[i]	+=	its_indep_normal_vector[j]	*	its_CholeskyMatrix->Getvalue(i,j);

	if (its_CopulaType == CCT_STUDENT)
	{
		double	Chi2Factor	=	((ICM_Gen_Student1F*) Gen)->GetChi2Factor();
		for (i=0;i<its_CreditsLabels.size();i++)
		{
			its_corr_normal_vector[i]	=	NAG_prob_students_t( its_corr_normal_vector[i], its_FreedomDegree );		
		}
	}


	for (i=0;i<its_CreditsLabels.size();i++)
	{	
		ri = its_corr_normal_vector[i];		
		if (ri < itsBarriers_Standard[i])
		{
			if (its_CorrelationType == CT_FACTOR_LOADING_2)
				aux	=	Gen->CumulativeDistribution(i, ri);
			else
				switch (its_CopulaType)
				{
				case CCT_GAUSSIAN:
					aux = NAG_cumul_normal(ri);
					break;
				case CCT_STUDENT:
					aux = NAG_prob_students_t( ri, its_FreedomDegree);
					break;
				}

			t = (its_ArrayDefaultCrv[i])->DefProbInverse(aux);
			itsSortedDefaultTimes.insert(Tau_Item(i,t));
		}
	}


/*	set<Tau_Item>::iterator iter;

	// ----------------------------------------------------------------------------
	if ((its_fOut = fopen("c:\\test\\mc_GenerateDefaultTimes_Betas.txt", "a+")) == NULL) return;
	fprintf(its_fOut, " ----------------- GenerateDefaultTimes_Betas ----------------- \n");
	// ----------------------------------------------------------------------------
	for (iter = itsSortedDefaultTimes.begin();
				iter != itsSortedDefaultTimes.end();
				++iter)
	{
		fprintf(its_fOut, "Id:\t%u\t\ttau:\t%.2lf\t\tCumul loss:\t%.2lf\t\tCumul NbDef:\t%u\n", iter->id, iter->tau, iter->cumul_loss, iter->cumul_nbdef);
	}
	fclose(its_fOut);
*/	
}	




void CreditManager::ComputeCumulativeLossesAndNbDefaults()
{
	int		SumNbDef	=	0;
	double	SumLosses	=	0.0;
	double	TheLossAmount	=	0.0;

	int	The_id;

	// I want to find all items between Tmin and Tmax
	set<Tau_Item>::iterator iter;


	iter = itsSortedDefaultTimes.begin();

	while (iter != itsSortedDefaultTimes.end())
	{
		iter->cumul_loss	=	SumLosses;		// excluded.
		The_id	=	iter->id;
		TheLossAmount	=	its_LossesAmount[The_id];
		SumLosses	+=	TheLossAmount;
		
		iter->cumul_nbdef	=	SumNbDef;		// excluded.
		SumNbDef++;

		iter->loss_at_that_time	=	TheLossAmount;

		iter++;

		// TO OPTIMIZE: early EXIT (but more tests)
/*		if (its_CreditObservationType == CO_LOSSES)
		{
			if (SumLosses >= its_Useful_MaxLoss)
			{
				while (iter != itsSortedDefaultTimes.end())
				{
					iter->cumul_loss	=	SumLosses;
					iter++;
				}
			}
		}
		else	// NB DEF
		{
			if (SumNbDef >= its_Useful_MaxNbDef)
			{
				while (iter != itsSortedDefaultTimes.end())
				{
					iter->cumul_nbdef	=	SumNbDef;
					iter++;
				}
			}
		}
*/
	}

	// ----------------------------------------------------------------------------
	// ----------------------------------------------------------------------------
	
	if (its_Display_CumulativeLossesAndDefaults_Flag)
	{
		// ----------------------------------------------------------------------------
		if ((its_fOut = fopen("c:\\test\\mc_ComputeCumulativeLossesAndNbDefaults.txt", "a+")) == NULL) return;
		fprintf(its_fOut, " ----------------- itsSortedDefaultTimes\t\t%u ----------------- \n", its_SimulId);
		// ----------------------------------------------------------------------------
		for (iter = itsSortedDefaultTimes.begin();
					iter != itsSortedDefaultTimes.end();
					++iter)
		{
			fprintf(its_fOut, "Id:\t%u\t\ttau:\t%.4lf\t\tCumul loss:\t%.2lf\t\tCumul NbDef:\t%u\n", iter->id, iter->tau, iter->cumul_loss, iter->cumul_nbdef);
		}
		fclose(its_fOut);
	}
	
}


void CreditManager::UpdateCumulativeLossesAndNbDefaults(int	Num)
{
	int		SumNbDef	=	0;
	double	SumLosses	=	0.0;
	double	tau, tau_shifted;
	double	TheLossAmount	=	0.0;

	int	The_id;

	// I want to find all items between Tmin and Tmax
	set<Tau_Item>::iterator iter;

	if (its_Bump_Choice != CHB_RECOVERY_LOSS)
	{

		// Num is the id of the Credit
		if (its_DefaultTimes[Num] == -1.0)
		{
			// the Credit Num default time was not inside the range
			if (its_DefaultTimesShifted[Num] == -1.0)
			{
				// nor the shifted, so nothing else to do
				// REMARK --> no PRICE change, I should avoid a repricing
				return;		
			}
			else
			{
				// The shifted one is inside.
				tau_shifted	=	its_DefaultTimesShifted[Num];
				itsSortedDefaultTimes.insert(Tau_Item(Num,tau_shifted));
			}
		}
		else
		{
			// the Credit Num default time is inside the range
			// remove it
			tau	=	its_DefaultTimes[Num];
			// to be optimized
			iter	=	itsSortedDefaultTimes.find(Tau_Item(Num,tau));
			itsSortedDefaultTimes.erase(iter);

			if (its_DefaultTimesShifted[Num] == -1.0)
			{
				// but the shifted is not
				;
			}
			else
			{
				// insert the shifted one
				tau_shifted	=	its_DefaultTimesShifted[Num];
				itsSortedDefaultTimes.insert(Tau_Item(Num,tau_shifted));
			}
		}
	}

	for (iter = itsSortedDefaultTimes.begin();
				iter != itsSortedDefaultTimes.end();
				++iter)
	{
		The_id	=	iter->id;

		iter->cumul_loss	=	SumLosses;		// excluded.
		SumLosses	+=	its_LossesAmount[The_id];
		
		iter->cumul_nbdef	=	SumNbDef;		// excluded.
		SumNbDef++;

		TheLossAmount	=	its_LossesAmount[The_id];
		iter->loss_at_that_time	=	TheLossAmount;
	}
/*
if (Num == 0)
{
	// ----------------------------------------------------------------------------
	if ((its_fOut = fopen("c:\\test\\mc_UpdateCumulativeLossesAndNbDefaults.txt", "a+")) == NULL) return;
	fprintf(its_fOut, " ----------------- itsSortedDefaultTimes\t\t%u\t\t%u ----------------- \n", its_SimulId,Num);
	// ----------------------------------------------------------------------------
	for (iter = itsSortedDefaultTimes.begin();
				iter != itsSortedDefaultTimes.end();
				++iter)
	{
		fprintf(its_fOut, "Id:\t%u\t\ttau:\t%.2lf\t\tCumul loss:\t%.2lf\t\tCumul NbDef:\t%u\n", iter->id, iter->tau, iter->cumul_loss, iter->cumul_nbdef);
	}
	fclose(its_fOut);
}
*/
}

// ----------------------------------------------------------------------------
// Prepare Outputs
// ----------------------------------------------------------------------------

void CreditManager::ResetOutputs()
{
	NPV	=	0.0;
	DefaultLegPV	=	0.0;
	PremiumLegPV	=	0.0;
	ATMPremiumLeg	=	0.0;
	ATMPremiumLegWithoutNotio	=	0.0;
	ATMSpread		=	0.0;
	ATMUpFront		=	0.0;
	ATMRunningSpread	=	0.0;
	NPVSquare	=	0.0;

	StdError		=	0.0;

	switch (its_PricingLegsType)
	{
	case CBN_PREMIUMLEGONLY:
		PremNPVFlag	=	1.0;
		DefNPVFlag	=	0.0;
		break;

	case CBN_DEFAULTLEGONLY:
		PremNPVFlag	=	0.0;
		DefNPVFlag	=	1.0;
		break;

	case CBN_DEFAULTMINUSPREMIUM:
		PremNPVFlag	=	-1.0;
		DefNPVFlag	=	1.0;
		break;

	case CBN_PREMIUMMINUSDEFAULT:
		PremNPVFlag	=	1.0;
		DefNPVFlag	=	-1.0;
		break;

	}
}

void CreditManager::PrepareOutputs(DoubleVector& Outputs)
{
	Outputs.clear();
	ResetOutputs();
}


void CreditManager::SetOutputData()
{
	// AVERAGES
	if (its_CreditModelType == CMT_MONTECARLO)
	{
		NPV	/=	its_NbSimul;
		DefaultLegPV	/=	its_NbSimul;
		PremiumLegPV	/=	its_NbSimul;
		ATMPremiumLeg	/=	its_NbSimul;
		ATMPremiumLegWithoutNotio	/=	its_NbSimul;
		
		// Standard Error
		NPVSquare	/=	its_NbSimul;
		if (NPVSquare - NPV * NPV >= 0.0)
			StdError	=	sqrt(NPVSquare - NPV * NPV);
		else
			StdError	=	_MINUS_INFINITY_;	// Pb!
	}

	double	UpFront_PV	=	0.0;

	switch (its_NPV_Type)
	{
	case CNPV_STANDARD:
		break;
		
	case CNPV_WITH_RUNNING_SPREAD:
		UpFront_PV	=	its_UpFront_RunningSpread * ATMPremiumLeg;
		break;
		
	case CNPV_WITH_UP_FRONT_PREMIUM:
		UpFront_PV	=	its_UpFront_Premium * its_PL_Notios[0];

		break;
		
	case CNPV_WITH_RUNNING_SPREAD_AND_UP_FRONT_PREMIUM:
		UpFront_PV	=	its_UpFront_Premium * its_PL_Notios[0];
		UpFront_PV	+=	its_UpFront_RunningSpread * ATMPremiumLeg;
		
		break;
	}


	if (UpFront_PV)
	{
		PremiumLegPV	+=	UpFront_PV;

		NPV	+=	PremNPVFlag * UpFront_PV;
	}
	if (ATMPremiumLeg)
		ATMSpread	=	DefaultLegPV / ATMPremiumLeg * 10000.0;	// in bps.
	else
		ATMPremiumLeg	=	_MINUS_INFINITY_;	


	switch (its_ATMDataFlag)
	{
	case CPLADT_PURESPREAD:
		break;

	case CPLADT_UF_UPFRONT:
		// find the ATM UpFront
		// I have some Running Spread
		// ATMSpread and its_UpFront_RunningSpread (=500) in bps.
		// Output is ATMUpFront in % (30)
		ATMUpFront	=	(ATMSpread - its_UpFront_RunningSpread) / 100.0 * ATMPremiumLegWithoutNotio;

		break;

	case CPLADT_UF_SPREAD:
		// find the ATM Running Spread 
		// I have an UpFrontValue
		ATMRunningSpread	= ATMSpread - its_UpFront_Premium * 100.0 / ATMPremiumLegWithoutNotio;

		break;
	}


}


void CreditManager::ScheduleCheckings()
{
	// ------------------------------------------------------------------------------
	// if Default Leg is paid at Fixed Lag
	// ensure that the Payment date + Lag is >= TMax.
	if (its_DL_PaymentType == CEP_ATFIXEDDATE)
	{
		if (its_DL_PaymentDate + its_DL_PaymentLag < its_DL_CreditWindowLow)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
		        "Parameters:  Mismatch for Default Leg Fixed Payment specification");
	}
	// ------------------------------------------------------------------------------
}


void CreditManager::ComputeMaxLevels()
{
	// ------------------------------------------------------------------------------
	// Get MaxLoss % or MaxNbDef
	//	this information will be used later on, in order to short MC computation
	// ------------------------------------------------------------------------------
	its_Useful_MaxLoss	=	0.0;
	its_Useful_MaxNbDef	=	0.0;

	double	TheValue;
	DoubleVector::iterator iter;

	if (its_CreditObservationType == CO_LOSSES)
	{
		// first init with Default Leg
		if (its_DL_LossMax < its_DL_LossMin)
			ICMTHROW(ERR_INVALID_DATA,"Default Leg: Max Loss must be > Min Loss!");

		its_Useful_MaxLoss	=	its_DL_LossMax;

		// then scan Premium Leg
		for (iter = its_PL_LossMaxs.begin();
					iter != its_PL_LossMaxs.end();
					++iter)
		{
			TheValue	=	(*iter);
			if (TheValue > its_Useful_MaxLoss)
				its_Useful_MaxLoss	=	TheValue;
		}

		its_Useful_MaxLoss	*=	TheBasketNotional;
	}
	else
	{
		// NbDef
		// first init with Default Leg
		if (its_DL_NbDefMax < its_DL_NbDefMin)
			ICMTHROW(ERR_INVALID_DATA,"Default Leg: Nb Def Max must be > Nb Def Loss!");

		its_Useful_MaxNbDef	=	its_DL_NbDefMax;

		// then scan Premium Leg
		for (iter = its_PL_NbDefMaxs.begin();
					iter != its_PL_NbDefMaxs.end();
					++iter)
		{
			TheValue	=	(*iter);
			if (TheValue > its_Useful_MaxNbDef)
				its_Useful_MaxNbDef	=	TheValue;
		}

	}

}


void	CreditManager::GetDataFromLabel(string TheLabel, double& TheValue)
{
	// transform to Upper
	transform(TheLabel.begin(), TheLabel.end(), TheLabel.begin(), toupper);

	if (! TheLabel.compare("NPV"))	// return 0 if they are the same. must I imposed Upper???
		TheValue	=	NPV;
	else if (! TheLabel.compare("PREMIUM_PV"))
		TheValue	=	PremiumLegPV;
	else if (! TheLabel.compare("DEFAULT_PV"))
		TheValue	=	DefaultLegPV;
	else if (! TheLabel.compare("RPV01"))
		TheValue	=	ATMPremiumLegWithoutNotio;
	else if (! TheLabel.compare("ATM_SPREAD"))
		TheValue	=	ATMSpread;
	else if (! TheLabel.compare("ATM_UPFRONT"))
		TheValue	=	ATMUpFront;
	else if (! TheLabel.compare("ATM_RUNNINGSPREAD"))
		TheValue	=	ATMRunningSpread;
	else if (! TheLabel.compare("STD_ERROR"))
		TheValue	=	StdError;
	else if (! TheLabel.compare("DATA_MATRIX_NB_ROWS"))
		TheValue	=	its_Data_Matrix_NbRows;
	else if (! TheLabel.compare("DATA_MATRIX_NB_COLUMNS"))
		TheValue	=	its_Data_Matrix_NbColumns;
	else if (! TheLabel.compare("CASH"))
		TheValue	=	its_Cash;
	else if (! TheLabel.compare("ACCRUED"))
		TheValue	=	its_Accrued;
	else
		TheValue	=	0.0;

}


void	CreditManager::GetDataMatrixFromLabel(string TheLabel, vector<double*>& OutputMatrix, vector<string>& OutputLabels, int& OutputNbRows, int& OutputNbCols)
{
	// transform to Upper
	transform(TheLabel.begin(), TheLabel.end(), TheLabel.begin(), toupper);

	int i, j, the_size;

	// ---------------------------------------------------------------
	// compare: return 0 if they are the same. must I imposed Upper???
	// ---------------------------------------------------------------
	OutputLabels.clear();
	OutputMatrix.clear();

	ICM_QMatrix<double>*	CC_Outputs;
	CC_Outputs	=	NULL;

	if (! TheLabel.compare("CREDIT_CALIBRATOR_LAGS"))
	{
		the_size	=	its_CreditCalibrator_Lags.size();

		OutputLabels.push_back(its_CC_Lags_String);
//		OutputLabels.push_back("DAY LAG");

		OutputNbCols	=	1;	
		OutputMatrix.resize(OutputNbCols);

		for (j=0;j<OutputNbCols;j++)
			OutputMatrix[j] = new double[the_size];

		for (i=0;i<the_size;i++)
			OutputMatrix[0][i]	=	its_CreditCalibrator_Lags[i];

	}
	else if (! TheLabel.compare("CREDIT_CALIBRATOR_PROB_DEFAULT"))
	{
		OutputLabels.push_back(its_CC_Prob_String);

//		OutputLabels.push_back("DEFAULT PROB");
//		its_CC_Prob_String	=	"LOSS DENSITY";
		
		the_size	=	its_CreditCalibrator_Prob_Default.size();

		OutputNbCols	=	1;	
		OutputMatrix.resize(OutputNbCols);

		for (j=0;j<OutputNbCols;j++)
			OutputMatrix[j] = new double[the_size];

		for (i=0;i<the_size;i++)
			OutputMatrix[0][i]	=	its_CreditCalibrator_Prob_Default[i];
	}
	else if (! TheLabel.compare("CREDIT_CALIBRATOR_STD_ERROR"))
	{
		OutputLabels.push_back(its_CC_Std_Error_String);

//		OutputLabels.push_back("DEFAULT PROB");
//		its_CC_Prob_String	=	"LOSS DENSITY";
		
		the_size	=	its_CreditCalibrator_Std_Error.size();

		OutputNbCols	=	1;	
		OutputMatrix.resize(OutputNbCols);

		for (j=0;j<OutputNbCols;j++)
			OutputMatrix[j] = new double[the_size];

		for (i=0;i<the_size;i++)
			OutputMatrix[0][i]	=	its_CreditCalibrator_Std_Error[i];
	}
	else if ((! TheLabel.compare("CREDIT_HEDGES_PARALLEL")) || (! TheLabel.compare("CREDIT_HEDGES_DEFAULT"))
					|| (! TheLabel.compare("CREDIT_HEDGES_PARALLEL_FAST")))
	{
		// outputs will be 
		//	-	Labels
		//	-	its_Hedges_Delta_NPV;
		//	-	its_Hedges_Hedge_Ratio;
		//	-	its_Hedges_Delta_CDS;
		//	-	its_Hedges_CDS_ATM_Margin;
		//	-	its_Hedges_Carry;
//		OutputLabels.push_back("LABELS");
		OutputLabels.push_back("DELTA NPV");
		OutputLabels.push_back("HEDGE RATIO");
		OutputLabels.push_back("DELTA CDS");
		OutputLabels.push_back("CDS ATM MARGIN");
		OutputLabels.push_back("CARRY");

		the_size	=	its_CreditsLabels.size();

//		OutputNbCols	=	6;	
		OutputNbCols	=	5;
		OutputMatrix.resize(OutputNbCols);

		for (j=0;j<OutputNbCols;j++)
			OutputMatrix[j] = new double[the_size];

		for (i=0;i<the_size;i++)
		{
			j=0;
//			OutputMatrix[j][i]	=	its_CreditsLabelsAsChar[i];
//			j++;
			OutputMatrix[j][i]	=	its_Hedges_Delta_NPV[i];
			j++;
			OutputMatrix[j][i]	=	its_Hedges_Hedge_Ratio[i];
			j++;
			OutputMatrix[j][i]	=	its_Hedges_Delta_CDS[i];
			j++;
			OutputMatrix[j][i]	=	its_Hedges_CDS_ATM_Margin[i];
			j++;
			OutputMatrix[j][i]	=	its_Hedges_Carry[i];
		}
	}
	else if ((! TheLabel.compare("CREDIT_HEDGES_RECOVERY_SENS")) || (! TheLabel.compare("CREDIT_HEDGES_RECOVERY_LOSS")))
	{
		// outputs will be 
		//	-	Labels
		//	-	its_Hedges_Delta_NPV;
		//	-	its_Hedges_Carry;
//		OutputLabels.push_back("LABELS");
		OutputLabels.push_back("DELTA NPV");
		OutputLabels.push_back("SENSITIVITY");

		the_size	=	its_CreditsLabels.size();

//		OutputNbCols	=	3;
		OutputNbCols	=	2;
		OutputMatrix.resize(OutputNbCols);

		for (j=0;j<OutputNbCols;j++)
			OutputMatrix[j] = new double[the_size];

		for (i=0;i<the_size;i++)
		{
			j=0;
//			OutputMatrix[j][i]	=	its_CreditsLabelsAsChar[i];
//			j++;
			OutputMatrix[j][i]	=	its_Hedges_Delta_NPV[i];
			j++;
			OutputMatrix[j][i]	=	its_Hedges_Carry[i];
		}
	}
	else if (! TheLabel.compare("CREDIT_HEDGES_CORRELATION"))
	{
		// outputs will be 
		//	-	Labels
		//	-	its_Hedges_Delta_NPV;
		//	-	its_Hedges_Carry;
//		OutputLabels.push_back("LABELS");
		OutputLabels.push_back("DELTA NPV");
		OutputLabels.push_back("SENSITIVITY");

		switch (its_CorrelationType)
		{
			case CT_FLAT:
				the_size	=	1;
				break;
			case CT_BETA:		// TO BE DONE, IF BETA FLAT, A SINGLE OUTPUT?
				the_size	=	its_CreditsLabels.size();
				break;
			case CT_FACTOR_LOADING_2:	// TO BE DONE, IF BETA FLAT, A SINGLE OUTPUT?
				the_size	=	its_CreditsLabels.size();
				break;

		}		

//		OutputNbCols	=	3;
		OutputNbCols	=	2;
		OutputMatrix.resize(OutputNbCols);

		for (j=0;j<OutputNbCols;j++)
			OutputMatrix[j] = new double[the_size];

		for (i=0;i<the_size;i++)
		{
			j=0;
//			OutputMatrix[j][i]	=	its_CreditsLabelsAsChar[i];
//			j++;
			OutputMatrix[j][i]	=	its_Hedges_Delta_NPV[i];
			j++;
			OutputMatrix[j][i]	=	its_Hedges_Carry[i];
		}
	}
	else if (! TheLabel.compare("CREDIT_HEDGES_RECOVERY_GLOBAL"))
	{
		// outputs will be 
		//	-	Labels
		//	-	its_Hedges_Delta_NPV;
		//	-	its_Hedges_Carry;
//		OutputLabels.push_back("LABELS");
		OutputLabels.push_back("DELTA NPV");
		OutputLabels.push_back("SENSITIVITY");

		the_size	=	1;

//		OutputNbCols	=	3;
		OutputNbCols	=	2;
		OutputMatrix.resize(OutputNbCols);

		for (j=0;j<OutputNbCols;j++)
			OutputMatrix[j] = new double[the_size];

		for (i=0;i<the_size;i++)
		{
			j=0;
//			OutputMatrix[j][i]	=	its_CreditsLabelsAsChar[i];
//			j++;
			OutputMatrix[j][i]	=	its_Hedges_Delta_NPV[i];
			j++;
			OutputMatrix[j][i]	=	its_Hedges_Carry[i];
		}
	}
	else if (! TheLabel.compare("CORRELATION_CALIBRATOR"))
	{
		// outputs will be
		// Market Implied Correlation
		// Market Up-Front Premium (0.0)
		// Market Implied Spreads (no up-Front right now)
		
		// Model Dependant:
		// - Base Correlation
		
		OutputLabels.push_back("MARKET IMPLIED CORRELATION");
		OutputLabels.push_back("MARKET IMPLIED UP-FRONT PREMIUM");
		OutputLabels.push_back("MARKET IMPLIED SPREAD");

		// TO IMPROVE
		OutputLabels.push_back("BASE CORRELATION");

		its_Correlation_Calibrator->Get_Correlation_Calibrator_Outputs(CC_Outputs);

		the_size	=	CC_Outputs->Getnbrows();
		// test nbrows?

		OutputNbCols	=	CC_Outputs->Getnbcols();;
		OutputMatrix.resize(OutputNbCols);

		for (j=0; j<OutputNbCols; j++)
		{
			OutputMatrix[j] = new double[the_size];

			for (i=0; i<the_size; i++)
				OutputMatrix[j][i]	=	(*CC_Outputs)(i,j);
		}

	}
	
	OutputNbRows	=	the_size;

	its_Data_Matrix_NbRows		=	OutputNbRows+1;	// label included!
	its_Data_Matrix_NbColumns	=	OutputNbCols;
}

//--------------------------------------------
//	CREDIT CALIBRATOR PARAMETERS
//--------------------------------------------

void	CreditManager::SetCreditCalibrator_Description(ICM_Matrix<ARM_Vector>* parameters)
{
	if (parameters == NULL) return;
	
	int TheValue;
	ARM_Vector* THE_VECT=NULL;

	// ------------------------------------------------------
	// CALIBRATOR_SELECTION:	INT --> ENUM

	THE_VECT	= parameters->GetColVect("CALIBRATOR_SELECTION");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CALIBRATOR_SELECTION");

	TheValue = (int) THE_VECT->Elt(0);	
	if ((TheValue < CCS_SINGLENAME) || (TheValue > CCS_EXPECTED_LOSS)) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CALIBRATOR_SELECTION");
	
	its_CreditCalibratorChoice = (CreditCalibratorSelection) TheValue;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CC_LABEL:	INT

	THE_VECT	= parameters->GetColVect("CC_LABEL");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CC_LABEL");
	
	its_CreditCalibrator_NameId = (int) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CC_NTD:	INT

	THE_VECT	= parameters->GetColVect("CC_NTD");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CC_NTD");

	its_CreditCalibrator_NTD = (int) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CC_LOSS:	DOUBLE

	THE_VECT	= parameters->GetColVect("CC_LOSS");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CC_LOSS");

	its_CreditCalibrator_Loss = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------


	// ------------------------------------------------------
	// CC_DENSITY_NTD:	INT

	THE_VECT	= parameters->GetColVect("CC_DENSITY_NTD");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CC_DENSITY_NTD");

	its_CreditCalibrator_Density_NTD = (int) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CC_DENSITY_LOSS_DOWN:	DOUBLE

	THE_VECT	= parameters->GetColVect("CC_DENSITY_LOSS_DOWN");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CC_DENSITY_LOSS_DOWN");

	its_CreditCalibrator_Density_Loss_Down = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CC_DENSITY_LOSS_UP:	DOUBLE

	THE_VECT	= parameters->GetColVect("CC_DENSITY_LOSS_UP");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CC_DENSITY_LOSS_UP");

	its_CreditCalibrator_Density_Loss_Up = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CC_DENSITY_LOSS_STEP:	DOUBLE

	THE_VECT	= parameters->GetColVect("CC_DENSITY_LOSS_STEP");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CC_DENSITY_LOSS_STEP");

	its_CreditCalibrator_Density_Loss_Step = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CC_EXPECTED_LOSS_DOWN:	DOUBLE

	THE_VECT	= parameters->GetColVect("CC_EXPECTED_LOSS_DOWN");
	if (!THE_VECT)
//		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
  //          "Parameters :  bad CC_EXPECTED_LOSS_DOWN");
		its_CreditCalibrator_Density_Loss_Up = (double) 0.0;
	else
		its_CreditCalibrator_Expected_Loss_Down = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CC_EXPECTED_LOSS_UP:	DOUBLE

	THE_VECT	= parameters->GetColVect("CC_EXPECTED_LOSS_UP");
	if (!THE_VECT)
//		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
//            "Parameters :  bad CC_EXPECTED_LOSS_UP");
		its_CreditCalibrator_Density_Loss_Step = (double) 0.0;
	else
		its_CreditCalibrator_Expected_Loss_Up = (double) THE_VECT->Elt(0);
	// ------------------------------------------------------
}

void	CreditManager::Credit_Calibrate(const StringVector& Maturities)
{
	// ------------------------------------------------------
	// Choice has to be made between
	// Calibrate Basket - Single Name / NTD / NLOSS / All
	// ------------------------------------------------------

	CreditCalibratorSelection	TheCreditCalibratorChoice;
	const ICM_DefaultCurve*	TheDefaultCurve=NULL;
	const ARM_Vector*			TheNbDays=NULL;
	const ARM_Vector*	TheSurvivalProba=NULL;
	
	int i;
	int size;

	GetCreditCalibratorChoice(TheCreditCalibratorChoice);

	// CALIBRATION HAS NOT BEEN CARRIED OUT.
	if (its_ModelMultiCurves == NULL)
		GenerateAllDefaultCurves();
//		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
//			"Parameters:  Model Multi Curve is NULL!");

	AllocateStructures_1F();

	switch (TheCreditCalibratorChoice)
	{
	case CCS_SINGLENAME:
		
		//	Retrieve Name from Id
		TheDefaultCurve	=	its_ModelMultiCurves->GetDefaultCurve(its_CreditCalibrator_NameId);
		if (TheDefaultCurve	== NULL)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Parameters:  invalid Single Name");

		// Calibrate Lambda
		/** JLA: assumed calibrated, or calibrated on demand.
		TheDefaultCurve->Calibrate();
		TheDefaultCurve->CptTermsSurvivalProba(); **/ 

		// Set the Curve
		TheNbDays	=	TheDefaultCurve->GetNbDays();
		if (TheNbDays == NULL)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Parameters:  invalid Single Name NbDays");
		
		TheSurvivalProba	=	TheDefaultCurve->GetSurvivalProba();
		if (TheSurvivalProba == NULL)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Parameters:  invalid Single Name TheSurvivalProba");

		its_CreditCalibrator_Lags.clear();
		its_CreditCalibrator_Prob_Default.clear();

		size = TheNbDays->GetSize();
		for (i=0; i<size; i++)
		{
			its_CreditCalibrator_Lags.push_back(TheNbDays->Elt(i));
			its_CreditCalibrator_Prob_Default.push_back(TheSurvivalProba->Elt(i));
		}

		its_CC_Lags_String	=	"LAGS (days)";
		its_CC_Prob_String	=	"DEFAULT PROBA";

		break;

	case CCS_NTD:
	case CCS_NLOSS:
	case CCS_DENSITY_NTD:
	case CCS_DENSITY_LOSS:
	case CCS_EXPECTED_LOSS:
		// Single Calibration for every name
		// second test (its_ModelMultiCurves == NULL) is already checked
		GenerateAllDefaultCurves();	// only if required

		BasketCalibration(TheCreditCalibratorChoice);

		break;
	}

}


void	CreditManager::BasketCalibration(CreditCalibratorSelection TheCreditCalibratorChoice)
{
	double	TotalRef;

	DoubleVector	TheLocalMaturities;

	// -------------------------------------------------------
	// Compute Useful Data for the Basket
	// -------------------------------------------------------
	ComputeBasketNotional();

	// ----------------------------------------------------------------------------------
	if (TheCreditCalibratorChoice == CCS_NTD)
	{
		if (its_CreditCalibrator_NTD > its_CreditsLabels.size())
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Parameters:  NTD Integer should be less than Nb Credits");

		TotalRef = (double) its_CreditCalibrator_NTD;

		// -------------------------------------------------------
		// Compute Max Nb Def and Loss Levels
		// -------------------------------------------------------
		its_CreditObservationType	=	CO_NBDEFAULTS;
		its_Useful_MaxNbDef	=	TotalRef;

	}
	else if (TheCreditCalibratorChoice == CCS_NLOSS)
	{
		if ((its_CreditCalibrator_Loss <= 0.0) || (its_CreditCalibrator_Loss > 1.0))
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Parameters:  Loss Pct should be in ]0;1]!");

		TotalRef = (double) its_CreditCalibrator_Loss;
		TotalRef *= TheBasketNotional;

		// -------------------------------------------------------
		// Compute Max Nb Def and Loss Levels
		// -------------------------------------------------------
		its_CreditObservationType	=	CO_LOSSES;
		its_Useful_MaxLoss	=	TotalRef;

	}
	else if (TheCreditCalibratorChoice == CCS_DENSITY_NTD)
	{
		if (its_CreditCalibrator_Density_NTD > its_CreditsLabels.size())
			ICMTHROW(ERR_INVALID_DATA,"Parameters:  NTD Integer should be less than Nb Credits for NTD Density Calibration!");

		// this represents the maximum of events to be taken into account
		TotalRef = (double) its_CreditCalibrator_Density_NTD;

		// -------------------------------------------------------
		// Compute Max Nb Def and Loss Levels
		// -------------------------------------------------------
		its_CreditObservationType	=	CO_NBDEFAULTS;
		its_Useful_MaxNbDef	=	TotalRef;

	}
	else if (TheCreditCalibratorChoice == CCS_DENSITY_LOSS)
	{
/*
		// DO THE TESTS HERE?
		
		if ((its_CreditCalibrator_Loss <= 0.0) || (its_CreditCalibrator_Loss > 1.0))
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Parameters:  Loss Pct should be in ]0;1]!");
*/
		TotalRef = (double) its_CreditCalibrator_Density_Loss_Up;
		TotalRef *= TheBasketNotional;

		// -------------------------------------------------------
		// Compute Max Nb Def and Loss Levels
		// -------------------------------------------------------
		its_CreditObservationType	=	CO_LOSSES;
		its_Useful_MaxLoss	=	TotalRef;

	}
	else if (TheCreditCalibratorChoice == CCS_EXPECTED_LOSS)
	{
		;
	}
	else
		return;
	// ----------------------------------------------------------------------------------

	// ----------------------------------------------------------------------------------
	// What about Maturities
	// ----------------------------------------------------------------------------------
	// Roll Dates
	int	i;
	qCDS_ADJ	CDS_Roll_Type;
	CDS_Roll_Type	=	(its_RollDateFlag) ? qCredit_Adjust20 : qCredit_Default;
	
	ARM_Date	MaturityDate;
	double		Maturity_YTi;
	bool	AdjorNot	=	true;	// Business Days

	TheLocalMaturities.clear();

	for (i=0; i<its_CreditCalibratorNbMaturities; i++)
	{
		// from Maturities to Dates
		// TO BE CHANGED, ARM_DEFAULT_CURRENCY
		ARM_Security* psec = GetSecurity();
		string ccy("");
		if (psec == NULL) ccy = ARM_DEFAULT_COUNTRY;
		else ccy = string(GetSecurity()->GetCurrencyUnit()->GetCcyName());
		MaturityDate	= AddPeriod(itsValDate, its_CreditCalibratorMaturities_AsChar[i] , 
									ccy, /*ARM_DEFAULT_COUNTRY,*/
									AdjorNot, CDS_Roll_Type);
		Maturity_YTi	= (MaturityDate.GetJulian() - itsValDate.GetJulian())/K_YEAR_LEN;
		TheLocalMaturities.push_back(Maturity_YTi);
	}

	switch (TheCreditCalibratorChoice)
	{
	case CCS_NTD:
	case CCS_NLOSS:

		BasketCalibrationMaturities(TheCreditCalibratorChoice, TheLocalMaturities, TotalRef);

		break;

	case CCS_DENSITY_NTD:
	case CCS_DENSITY_LOSS:

		BasketCalibrationDensity(TheCreditCalibratorChoice, Maturity_YTi);

		break;

	case CCS_EXPECTED_LOSS:

//		BasketCalibrationExpectedLoss(TheCreditCalibratorChoice, Maturity_YTi);

		break;
	}

	return;
}


void	CreditManager::BasketCalibrationMaturities(CreditCalibratorSelection TheCreditCalibratorChoice, DoubleVector& TheLocalMaturities, double TotalRef)
{
	double	ProbTi;
	double	EventBeforeT;
	int		i;
	double	Maturity_YTi;
	
	// -------------------------------------------------------
	// Compute Barriers for Default Times
	// -------------------------------------------------------
	// ComputeMaxCreditDateForCalibration();
	its_MaxCreditDate	=	DAYSTIME(TheLocalMaturities[its_CreditCalibratorNbMaturities-1]);

	ComputeBarriers(its_MaxCreditDate);

	// -------------------------------------------------------
	// According to Correlation
	// -------------------------------------------------------
	switch (its_CorrelationType)
	{
		case CT_FLAT:
		case CT_BETA:			
		case CT_FACTOR_LOADING_2:
			break;

		case CT_MATRIX:
			// Cholesky computation
			CholeskyComputation();
			break;

	}

	// -------------------------------------------------------
	// The Pricing Loop  
	// -------------------------------------------------------
	DoubleVector	Probas;
	DoubleVector	TheProbas;

	DoubleVector	StdErrors;

	double	TheValue, TheSquareValue;

	Probas.resize(its_CreditCalibratorNbMaturities);
	TheProbas.resize(1);
	StdErrors.resize(its_CreditCalibratorNbMaturities);

	switch (its_CreditModelType)
	{
	case CMT_MONTECARLO:

		// -------------------------------------------------------
		// Variance Reduction useful data computation
		// -------------------------------------------------------
		ComputeVarianceReductionData();

		// -------------------------------------------------------
		// Initialization of the Random Generator
		// -------------------------------------------------------
		SetRandomGenerator();

		// -------------------------------------------------------
		// VARIANCE REDUCTION
		// -------------------------------------------------------

		for (its_SimulId=0;its_SimulId<its_NbSimul;its_SimulId++)
		{
			// -------------------------------------------------------
			// According to Correlation
			// -------------------------------------------------------
			switch (its_CorrelationType)
			{
				case CT_FLAT:
				case CT_BETA:				
				case CT_FACTOR_LOADING_2:
					GenerateDefaultTimes_Betas();
					break;

				case CT_MATRIX:
					GenerateDefaultTimes_Matrix();		// Cholevsky
					break;

			}

			// -------------------------------------------------------
			// ComputeCumulativeLossesAndNbDefaults
			// -------------------------------------------------------
			if (IsCDOStandard())
				ComputeCumulativeLossesAndNbDefaults();
			else // if (IsCDOSquare())
				ComputeCumulativeLossesAndNbDefaults_CDO_SQUARE();

			// -------------------------------------------------------
			// For all Maturities, get the Proba
			// -------------------------------------------------------
			for (i=0; i<its_CreditCalibratorNbMaturities; i++)
			{
				Maturity_YTi = TheLocalMaturities[i];

				// may use a pointer.
				if (TheCreditCalibratorChoice == CCS_NLOSS)
					GetHowMuchLossBeforeADate(Maturity_YTi, EventBeforeT);
				else
					GetHowManyDefaultsBeforeADate(Maturity_YTi, EventBeforeT);

				// Not optimal, when dealing with NTD or first Loss
				if (EventBeforeT >= TotalRef)
				{
					TheValue	=	1.0;					
					switch (its_MC_Variance_Reduction)
					{					
					case CMCVR_IS_PURE_FACTORS:
						TheValue	=	IS_Individual_FactorShift;
						
						Probas[i] += TheValue;
						StdErrors[i] += TheValue * TheValue;
						
						break;

					case CMCVR_IS_FACTORS:
						TheValue	=	IS_Individual_FactorShift;

					case CMCVR_IS_IDIOSYNCRATIC:
						TheValue		*=	exp(- its_IS_Theta * EventBeforeT / TheBasketNotional);		// scale because it is in %
						TheValue		*=	IS_Common_exp_twists;
						
						Probas[i] += TheValue;
						StdErrors[i] += TheValue * TheValue;
						
						break;

					case CMCVR_NONE:
					default:
						Probas[i] += 1.0;

						break;

					}
				}
			}
		}

		//// Apply a common
		//DoubleVector::iterator	iter_StdErr;
		//DoubleVector::iterator	iter;

		// ------------------------------------
		// VARIANCE REDUCTION
		// ------------------------------------
		TheValue	=	1.0;
		TheSquareValue	=	1.0;

		switch (its_MC_Variance_Reduction)
		{					
		case CMCVR_IS_PURE_FACTORS:
		case CMCVR_IS_FACTORS:
			TheValue		=	IS_Common_FactorShift;
			TheSquareValue	=	TheValue * TheValue;

// FIXMEFRED: mig.vc8 (28/05/2007 14:58:10):portée
		case CMCVR_IS_IDIOSYNCRATIC:
			{
				DoubleVector::iterator	iter_StdErr;
				DoubleVector::iterator	iter;
				for (iter = Probas.begin(), iter_StdErr = StdErrors.begin(); iter != Probas.end(); ++iter, iter_StdErr++)
				{
					(*iter) /= (double) its_NbSimul;

					// StdErrors for each Maturity
					(*iter_StdErr)	/= (double) its_NbSimul;
					(*iter_StdErr)	-= (*iter) * (*iter);
				}
			}
			break;

		case CMCVR_NONE:
		default:
			{
				DoubleVector::iterator	iter_StdErr;
				DoubleVector::iterator	iter;
				for (iter = Probas.begin(), iter_StdErr = StdErrors.begin(); iter != Probas.end(); ++iter, iter_StdErr++)
				{
					(*iter) /= (double) its_NbSimul;

					// StdErrors for each Maturity
					(*iter_StdErr)	=	(*iter) * (1 - (*iter));
				}

				break;
			}
		}
		
		// ------------------------------------
		
		its_CC_Std_Error_String	=	"STANDARD ERROR";

			break;

	case CMT_ANALYTIC_RECURSIVE_1F:

	case CMT_ANALYTIC_RECURSIVE_INTERPOLATION_1F:

		CreditLossComputation	TheCreditLossComputation;

		TheCreditLossComputation	=	(its_CreditModelType == CMT_ANALYTIC_RECURSIVE_1F) ? CLC_PROBABILITY : CLC_PROBABILITY_INTERPOLATION;

		// -------------------------------------------------------
		// For all Maturities, get the Proba, and sum it up
		// -------------------------------------------------------
		for (i=0; i<its_CreditCalibratorNbMaturities; i++)
		{
			Maturity_YTi = TheLocalMaturities[i];

			// may use a pointer.
			if (TheCreditCalibratorChoice == CCS_NTD)
				ProbabilityOfAtLeastNDefaultBeforeT_1F(Maturity_YTi, TotalRef, ProbTi);
			else
			{
				// I do want to compute Loss Probabilities > x%
				// may be improved, taking into account a Tranche [x% ; y%]
				ComputeProbabilityOrExpectationLossesBeforeT_RECURSIVE_1F(Maturity_YTi, TheCreditLossComputation, TotalRef, TheBasketNotional, ProbTi);
			}

			Probas[i] = ProbTi;
		}
		
		its_CC_Std_Error_String	=	"";

		break;

	case CMT_ANALYTIC_RECURSIVE_STOCHASTIC_CORRELATION_BERNOULLI_1F:
		ICMTHROW(ERR_INVALID_ARGUMENT,"BasketCalibrationMaturities not implemented for CMT_ANALYTIC_RECURSIVE_STOCHASTIC_CORRELATION_BERNOULLI_1F yet!");
		break;

	case CMT_ANALYTIC_LARGE_PORTFOLIO_1F:
		ICMTHROW(ERR_INVALID_ARGUMENT,"BasketCalibrationMaturities not implemented for CMT_ANALYTIC_LARGE_PORTFOLIO_1F yet!");
		break;
	
	case CMT_ANALYTIC_FFT_1F:
		ICMTHROW(ERR_INVALID_ARGUMENT,"BasketCalibrationMaturities not implemented for CMT_ANALYTIC_FFT_1F yet!");

		// -------------------------------------------------------
		// For all Maturities, get the Proba, and sum it up
		// -------------------------------------------------------
		for (i=0; i<its_CreditCalibratorNbMaturities; i++)
		{
			Maturity_YTi = TheLocalMaturities[i];

			// may use a pointer.
			if (TheCreditCalibratorChoice == CCS_NTD)
				// to be changed?
				ProbabilityOfAtLeastNDefaultBeforeT_1F(Maturity_YTi, TotalRef, ProbTi);
			else
			{
				// I do want to compute Loss Probabilities > x%
				// may be improved, taking into account a Tranche [x% ; y%]
				ComputeProbabilityOrExpectationLossesBeforeT_FFT_1F(Maturity_YTi, CLC_PROBABILITY, TotalRef, its_TotalLoss, TheProbas);
				ProbTi	=	TheProbas[0];
			}

			Probas[i] = ProbTi;
		}
		
		its_CC_Std_Error_String	=	"";


		break;

	case CMT_ANALYTIC_LHP:
	case CMT_ANALYTIC_LHP_JPM:
	case CMT_ANALYTIC_LHP_PLUS:
		ICMTHROW(ERR_INVALID_ARGUMENT,"BasketCalibrationMaturities not implemented for CMT_ANALYTIC_LHP yet!");
		break;
	
	}

	// Here, one has [1 - exp(-Integrale(Lambda))] over all Ti maturities
	// reste plus qu'à construire Lambda

	DoubleVector	Values;
	double	LastDate, NextDate; 
	double	LastIntegral, /*NextIntegral, TheValue, */Proba, StdErr;

	LastDate = 0;
	LastIntegral = 0.0;

	its_CreditCalibrator_Lags.clear();
	its_CreditCalibrator_Prob_Default.clear();
	its_CreditCalibrator_Std_Error.clear();
	
	for (i=0; i<its_CreditCalibratorNbMaturities; i++)
	{
		NextDate = TheLocalMaturities[i];
		Proba = Probas[i];
		StdErr	=	StdErrors[i];

		its_CreditCalibrator_Lags.push_back(NextDate*365);		// from year fraction
		its_CreditCalibrator_Prob_Default.push_back(Proba);
		its_CreditCalibrator_Std_Error.push_back(StdErr);
	}

	its_CC_Lags_String	=	"LAGS (days)";
	its_CC_Prob_String	=	"DEFAULT PROBA";
	its_CC_Std_Error_String	=	"STD ERROR";
}


void	CreditManager::BasketCalibrationDensity(CreditCalibratorSelection TheCreditCalibratorChoice, double& Maturity_YTi)
{
	double	EventBeforeT;
	
	// -------------------------------------------------------
	// Initialization of the Random Generator
	// -------------------------------------------------------
	SetRandomGenerator();

	// -------------------------------------------------------
	// Compute Barriers for Default Times
	// -------------------------------------------------------
	// ComputeMaxCreditDateForCalibration();
	its_MaxCreditDate	=	DAYSTIME(Maturity_YTi);

	ComputeBarriers(its_MaxCreditDate);

	// -------------------------------------------------------
	// According to Correlation
	// -------------------------------------------------------
	switch (its_CorrelationType)
	{
		case CT_FLAT:
		case CT_BETA:			
		case CT_FACTOR_LOADING_2:
			break;

		case CT_MATRIX:
			// Cholesky computation
			CholeskyComputation();
			break;

	}
	
	DoubleVector	Probas;
	DoubleVector	StdErrors;

	if (TheCreditCalibratorChoice == CCS_DENSITY_NTD)
	{
		Probas.resize(its_CreditCalibrator_Density_NTD+2);
		StdErrors.resize(its_CreditCalibrator_Density_NTD+2);
	}
	else
	{
		// LOSS
		size_DENSITY_LOSS	=	(int) (1.0 / its_CreditCalibrator_Density_Loss_Step);

		if ((1.0 / its_CreditCalibrator_Density_Loss_Step) != size_DENSITY_LOSS)
			size_DENSITY_LOSS++;

		Probas.resize(size_DENSITY_LOSS+2);
		StdErrors.resize(size_DENSITY_LOSS+2);
	}

	fill(Probas.begin(), Probas.end(), 0.0);

	int	k;
	double	PortfolioLossInPct;
	double	TrancheLossInPct;
	double	TrancheSizeInPct;

	// -------------------------------------------------------
	// SOME CHECKINGS

	if (its_CreditCalibrator_Density_Loss_Up == its_CreditCalibrator_Density_Loss_Down)
		ICMTHROW(ERR_INVALID_DATA,"Bad Tranche Definition: NULL size!");

	if ((its_CreditCalibrator_Density_Loss_Step <= 0.0) || (its_CreditCalibrator_Density_Loss_Step > 1.0))
		ICMTHROW(ERR_INVALID_DATA,"Size must lie between ]0.0%;100.0%[ !");
	// -------------------------------------------------------

	TrancheSizeInPct	=	(its_CreditCalibrator_Density_Loss_Up - its_CreditCalibrator_Density_Loss_Down);

	// -------------------------------------------------------
	double	Result, LastResult;
	// -------------------------------------------------------

	// -------------------------------------------------------
	// The Pricing Loop  
	// -------------------------------------------------------

	switch (its_CreditModelType)
	{
	case CMT_MONTECARLO:
		
		for (its_SimulId=0;its_SimulId<its_NbSimul;its_SimulId++)
		{
			// -------------------------------------------------------
			// According to Correlation
			// -------------------------------------------------------
			switch (its_CorrelationType)
			{
				case CT_FLAT:
				case CT_BETA:				
				case CT_FACTOR_LOADING_2:
					GenerateDefaultTimes_Betas();
					break;

				case CT_MATRIX:
					GenerateDefaultTimes_Matrix();		// Cholevsky
					break;

			}

			// -------------------------------------------------------
			// ComputeCumulativeLossesAndNbDefaults
			// -------------------------------------------------------
			if (IsCDOStandard())
				ComputeCumulativeLossesAndNbDefaults();
			else // if (IsCDOSquare())
				ComputeCumulativeLossesAndNbDefaults_CDO_SQUARE();

			// -------------------------------------------------------
			// For the Maturity, get the Proba
			// -------------------------------------------------------

			if (TheCreditCalibratorChoice == CCS_DENSITY_NTD)
			{
				GetHowManyDefaultsBeforeADate(Maturity_YTi, EventBeforeT);
				Probas[FMIN(EventBeforeT, its_CreditCalibrator_Density_NTD+1)] += 1;
			}
			else
			{
				// CCS_DENSITY_LOSS
				GetHowMuchLossBeforeADate(Maturity_YTi, EventBeforeT);

				PortfolioLossInPct	=	EventBeforeT / TheBasketNotional;
				TrancheLossInPct	=	(FMAX(its_CreditCalibrator_Density_Loss_Up - PortfolioLossInPct, 0.0) - FMAX(its_CreditCalibrator_Density_Loss_Down - PortfolioLossInPct, 0.0));
				TrancheLossInPct	/=	TrancheSizeInPct;
				TrancheLossInPct	=	1.0 - TrancheLossInPct;

				// Integer Division
				if (CHECK_EQUAL(TrancheLossInPct, 0.0))
					k	=	0;
				else if (CHECK_EQUAL(TrancheLossInPct, 1.0))
					k	=	size_DENSITY_LOSS + 1;
				else
					k	=	(int) (TrancheLossInPct / its_CreditCalibrator_Density_Loss_Step) + 1;
				if ((k < 0) || (k >= its_CreditCalibrator_Density_NTD+1))
					ICMTHROW(ERR_INVALID_DATA,"Bad Index for Loss Density Computation with Monte-Carlo!");

				Probas[k]	+=	1;
			}
		}
		{
		// Apply a common
		DoubleVector::iterator	iter;
		DoubleVector::iterator	iter_StdErr;
		
		for (iter = Probas.begin(), iter_StdErr = StdErrors.begin(); iter != Probas.end(); ++iter, iter_StdErr++)
		{
			(*iter) /= (double) its_NbSimul;
			// StdErrors for each Maturity
			(*iter_StdErr)	=	(*iter) * (1 - (*iter));
		}
		}
		its_CC_Std_Error_String	=	"STANDARD ERROR";
			
		break;

	case CMT_ANALYTIC_RECURSIVE_1F:
	case CMT_ANALYTIC_RECURSIVE_INTERPOLATION_1F:

			CreditLossComputation	TheCreditLossComputation;

			TheCreditLossComputation	=	(its_CreditModelType == CMT_ANALYTIC_RECURSIVE_1F) ? CLC_PROBABILITY : CLC_PROBABILITY_INTERPOLATION;

			if (TheCreditCalibratorChoice == CCS_DENSITY_NTD)
				ProbabilityOfNDefaultsBeforeT_1F(Maturity_YTi, its_CreditCalibrator_Density_NTD, Probas);
			else
			{
				// I do want to compute Loss Probabilities > x%
				// may be improved, taking into account a Tranche [x% ; y%]
				
				// begin with Prob(Loss > its_CreditCalibrator_Density_Loss_Down (%))

				int k = 0;

				ComputeProbabilityOrExpectationLossesBeforeT_RECURSIVE_1F(Maturity_YTi, TheCreditLossComputation, its_CreditCalibrator_Density_Loss_Down * TheBasketNotional, 
							TheBasketNotional, Result);

				Probas[k] = 1.0 - Result;
				LastResult	=	Result;

				double	StepSizeInPct	=	(its_CreditCalibrator_Density_Loss_Up - its_CreditCalibrator_Density_Loss_Down) * its_CreditCalibrator_Density_Loss_Step;

				while (its_CreditCalibrator_Density_Loss_Step * k <= 1.0)
				{
					k++;

					ComputeProbabilityOrExpectationLossesBeforeT_RECURSIVE_1F(Maturity_YTi, TheCreditLossComputation,
							(its_CreditCalibrator_Density_Loss_Down + StepSizeInPct * k) * TheBasketNotional,
								TheBasketNotional, Result);

					Probas[k]	=	LastResult - Result;
					LastResult	=	Result;

				}
				
				// and the last one, Proba to be > at 100% of the interval
				Probas[k]	=	LastResult;

				its_CC_Std_Error_String	=	"";
			}
		
			break;

	case CMT_ANALYTIC_RECURSIVE_STOCHASTIC_CORRELATION_BERNOULLI_1F:
		ICMTHROW(ERR_INVALID_ARGUMENT,"BasketCalibrationDensity not implemented for CMT_ANALYTIC_RECURSIVE_STOCHASTIC_CORRELATION_BERNOULLI_1F yet!");
		break;
	
	case CMT_ANALYTIC_LARGE_PORTFOLIO_1F:
		ICMTHROW(ERR_INVALID_ARGUMENT,"BasketCalibrationDensity not implemented for CMT_ANALYTIC_LARGE_PORTFOLIO_1F yet!");
		break;
	
	case CMT_ANALYTIC_FFT_1F:
		ICMTHROW(ERR_INVALID_ARGUMENT,"BasketCalibrationDensity not implemented for CMT_ANALYTIC_FFT_1F yet!");

//			if (TheCreditCalibratorChoice == CCS_DENSITY_NTD)
				// to be changed?
//				ProbabilityOfNDefaultsBeforeT_1F(Maturity_YTi, its_CreditCalibrator_Density_NTD, Probas);
//			else
				// I do want to compute Loss Probabilities > x%
				// may be improved, taking into account a Tranche [x% ; y%]
//				ProbabilityOrExpectationLossesBeforeT_FFT_1F(Maturity_YTi, CLC_PROBABILITY, its_CreditCalibrator_Density_Loss_Down * TheBasketNotional, 
//							its_CreditCalibrator_Density_Loss_Up * TheBasketNotional, Probas);
		
			its_CC_Std_Error_String	=	"";
		
		break;

	case CMT_ANALYTIC_LHP:
	case CMT_ANALYTIC_LHP_JPM:
	case CMT_ANALYTIC_LHP_PLUS:
		ICMTHROW(ERR_INVALID_ARGUMENT,"BasketCalibrationDensity not implemented for CMT_ANALYTIC_LHP yet!");
		break;
	

	}

	double	Proba, StdErr;

	its_CreditCalibrator_Lags.clear();
	its_CreditCalibrator_Prob_Default.clear();
	its_CreditCalibrator_Std_Error.clear();
	
	int	i;

	its_CC_Std_Error_String	=	"STD ERROR";

	if (TheCreditCalibratorChoice == CCS_DENSITY_NTD)
	{
		its_CC_Lags_String	=	"NB DEFAULTS";
		its_CC_Prob_String	=	"LOSS DENSITY";

		for (i=0; i<its_CreditCalibrator_Density_NTD+2; i++)
		{			
			Proba = Probas[i];
			StdErr	=	StdErrors[i];

			its_CreditCalibrator_Lags.push_back((double) i);
			its_CreditCalibrator_Prob_Default.push_back(Proba);
			its_CreditCalibrator_Std_Error.push_back(StdErr);
		}
	}
	else // CCS_DENSITY_LOSS
	{
		its_CC_Lags_String	=	"LOSS %";
		its_CC_Prob_String	=	"LOSS DENSITY";

		for (i=0; i<=size_DENSITY_LOSS + 1; i++)
		{			
			Proba = Probas[i];
			StdErr	=	StdErrors[i];

			its_CreditCalibrator_Lags.push_back((double) i*its_CreditCalibrator_Density_Loss_Step * 100.0);
			its_CreditCalibrator_Prob_Default.push_back(Proba);
			its_CreditCalibrator_Std_Error.push_back(StdErr);
		}
	}
}



// ----------------------------------------------------
// HEDGES
// ----------------------------------------------------

// ----------------------------------------------------
// Hedges for standard parallel shift of spreads
// ----------------------------------------------------

void	CreditManager::Hedges(double& dummyNPV)
{
	switch (its_CreditModelType)
	{
	case CMT_MONTECARLO:
		Hedges_Monte_Carlo(dummyNPV);
		break;

	case CMT_ANALYTIC_LHP:
		Hedges_LHP(dummyNPV);
		break;

	case CMT_ANALYTIC_RECURSIVE_1F:
	case CMT_ANALYTIC_RECURSIVE_INTERPOLATION_1F:
		Hedges_Recursive(dummyNPV);
		break;

	case CMT_ANALYTIC_RECURSIVE_STOCHASTIC_CORRELATION_BERNOULLI_1F:
	case CMT_ANALYTIC_LARGE_PORTFOLIO_1F:
	case CMT_ANALYTIC_LHP_JPM:
	case CMT_ANALYTIC_FFT_1F:
		
	case CMT_ANALYTIC_LHP_PLUS:

		ICMTHROW(ERR_INVALID_DATA,"Hedges not implemented for this model yet!");

	}

	return;
}


void	CreditManager::Hedges_Monte_Carlo(double& dummyNPV)
{
	if (!its_HedgesRunning) return;

	int	i;
	double	Loss_i;
	DoubleVector	Outputs;
	
	ICM_DefaultCurve*	DefCurve	=	NULL;
	ICM_DefaultCurve*	ShiftedDefCurve	=	NULL;

	DoubleVector	TheLossesAmountShifted;

//	double	PremLeg;
//	double	DefLeg;

	// First passage
	if (its_CurrentHedgesIndex	==	-1)
	{
		fill(its_AllNPVS.begin(), its_AllNPVS.end(), 0.0);
		fill(its_AllDefLegs.begin(), its_AllDefLegs.end(), 0.0);
		fill(its_AllPremLegs.begin(), its_AllPremLegs.end(), 0.0);

		// -------------------------------------------------------
		// Reset Outputs
		// -------------------------------------------------------
		ResetOutputs();	// NPV to 0.0

		// -------------------------------------------------------
		// Schedule Checkings
		// -------------------------------------------------------
//		ScheduleCheckings();	// already done in Central NPV

		// -------------------------------------------------------
		// Compute Useful DF Data
		// -------------------------------------------------------
//		ComputeAllDF();			// already done in Central NPV

		// -------------------------------------------------------
		// Compute Useful Data for the Basket
		// -------------------------------------------------------
//		ComputeBasketNotional(); // already done in Central NPV

		// -------------------------------------------------------
		// Prepare for Pricing
		// -------------------------------------------------------
//		SortAllCreditObservationDates(); // already done in Central NPV

		// -------------------------------------------------------
		// All Credit Curves have already been generated
		// -------------------------------------------------------
	//	GenerateAllDefaultCurves();

		// -------------------------------------------------------
		// Initialization of the Random Generator
		// -------------------------------------------------------
		SetRandomGenerator();

		// -------------------------------------------------------
		// Compute Barriers for Default Times
		// -------------------------------------------------------
//		ComputeMaxCreditDate(); // already done in Central NPV
		
		ComputeBarriers(its_MaxCreditDate);
		if ((its_Bump_Choice != CHB_DEFAULT) && (its_Bump_Choice != CHB_RECOVERY_LOSS))
			ComputeBarriersShifted(its_MaxCreditDate);

		// -------------------------------------------------------
		// According to Correlation
		// -------------------------------------------------------
		switch (its_CorrelationType)
		{
			case CT_FLAT:
			case CT_BETA:				
			case CT_FACTOR_LOADING_2:
				break;

			case CT_MATRIX:
				// Cholesky computation
				CholeskyComputation();
				break;

		}

		// Computes Shifted Losses Amounts if required
		if ((its_Bump_Choice == CHB_RECOVERY_SENS) || (its_Bump_Choice == CHB_RECOVERY_LOSS) || (its_Bump_Choice == CHB_RECOVERY_GLOBAL))
		{
			TheLossesAmountShifted.resize(its_CreditsLabels.size());
			for (i=0; i<its_CreditsLabels.size(); i++)
				TheLossesAmountShifted[i]	=	(its_Input_Losses[i] - its_BumpRecovery) * its_Notionals[i];

		}

		// -------------------------------------------------------
		// The Pricing Loop  
		// -------------------------------------------------------
		for (its_SimulId=0;its_SimulId<its_NbSimul;its_SimulId++)
		{
			// -------------------------------------------------------
			// According to Correlation
			// -------------------------------------------------------
			switch (its_CorrelationType)
			{
				case CT_FLAT:
					// take care of correl = 1
				case CT_BETA:					
				case CT_FACTOR_LOADING_2:
					GenerateDefaultTimes_Betas_ForHedge();
					break;

				case CT_MATRIX:
					// Cholevsky
					break;
			}

			// -------------------------------------------------------
			// ComputeCumulativeLossesAndNbDefaults
			// with initial Curves (non shifted)
			// -------------------------------------------------------
			if (IsCDOStandard())
				ComputeCumulativeLossesAndNbDefaults();
			else // if (IsCDOSquare())
				ComputeCumulativeLossesAndNbDefaults_CDO_SQUARE();

			// Keep SortedDefaultTimes
			if (its_Bump_Choice != CHB_RECOVERY_LOSS)
				Keep_SortedDefaultTimes();

			// -------------------------------------------------------
			// Here is the difference for the Hedge
			// -------------------------------------------------------

//			Display_DefaultTimes();		// and shifted

			for (i=0;i<its_CreditsLabels.size();i++)
			{
				its_temp_index	=	i;

				// Specific Hedges Default is handled previously
				if ((its_Bump_Choice == CHB_RECOVERY_SENS) || (its_Bump_Choice == CHB_RECOVERY_LOSS) || (its_Bump_Choice == CHB_RECOVERY_LOSS))
				{
					Loss_i	=	its_LossesAmount[i];
					its_LossesAmount[i]	=	TheLossesAmountShifted[i];
				}

				// -------------------------------------------------------
				// UpdateCumulativeLossesAndNbDefaults
				// with initial Curves and the shifted one
				// -------------------------------------------------------
				UpdateCumulativeLossesAndNbDefaults(i);
				
				// -------------------------------------------------------
				// Price Basket For This Simul
				// -------------------------------------------------------
				PriceBasketForThisSimulation(Outputs);

				its_AllNPVS[i]		+=	Outputs[0];
				its_AllDefLegs[i]	+=	Outputs[1];
				its_AllPremLegs[i]	+=	Outputs[2];

				// Restore losses Amount
				if ((its_Bump_Choice == CHB_RECOVERY_SENS) || (its_Bump_Choice == CHB_RECOVERY_LOSS) || (its_Bump_Choice == CHB_RECOVERY_GLOBAL))
					its_LossesAmount[i]	=	Loss_i;

				if (its_Bump_Choice != CHB_RECOVERY_LOSS)
					// restore Sorted Default Times
					Restore_SortedDefaultTimes();

				if (IsHomogeneousBasket())
				{
					for (i=1; i<its_CreditsLabels.size(); i++)
					{
						its_AllNPVS[i]		+=	Outputs[0];
						its_AllDefLegs[i]	+=	Outputs[1];
						its_AllPremLegs[i]	+=	Outputs[2];
					}
				}

			}

		}

		dummyNPV	=	0.0;

	}
	else
	{
		// simply get the right value
		dummyNPV	=	its_AllNPVS[its_CurrentHedgesIndex];
		dummyNPV	/=	its_NbSimul;

//		PremLeg	=	its_AllPremLegs[its_CurrentHedgesIndex];
//		PremLeg	/=	its_NbSimul;

//		DefLeg	=	its_AllDefLegs[its_CurrentHedgesIndex];
//		DefLeg	/=	its_NbSimul;
	}

}



void	CreditManager::ComputeHedges()
{
/*
	if (IsHedgesCDSToMaturity())
	{
		// its_MaxDate
	}
*/
	switch (its_Bump_Choice)
	{
	case CHB_SPREAD_ALL_MATURITIES:
		ComputeHedgeSpreadAllMaturities();
		break;
	case CHB_SPREAD_PARALLEL:
	case CHB_FAST_SPREAD_PARALLEL:
		ComputeHedgeSpreadParallel();
		break;
	case CHB_DEFAULT:
		ComputeHedgeDefault();
		break;
	case CHB_RECOVERY_SENS:
		ComputeHedgeRecoverySens();
		break;
	case CHB_RECOVERY_LOSS:
		ComputeHedgeRecoveryLoss();
		break;
	case CHB_RECOVERY_GLOBAL:
		ComputeHedgeRecoveryGlobal();
		break;
	case CHB_CORRELATION:
		ComputeHedgeCorrelation();
		break;
	}
	
	if (its_MatrixShiftedDefaultCrv)
		delete its_MatrixShiftedDefaultCrv;
	its_MatrixShiftedDefaultCrv	=	NULL;

	return;
}


void	CreditManager::ComputeHedgeSpreadParallel()
{
	int i;

	double	dummy_NPV, NPV_Shifted;
	double	DiffCDS, HedgeRatio, ATM_Margin;

	// Reset All Outputs
	ResetAndResizeAllHedgesOutputs();

	//	CDS Prices
	//	Generates Shifted Default Curves
	//	CDS Shifted Prices
	PrepareHedges();		// its_nb_scenarii has been computed

	// All Prices will be computed
	its_HedgesRunning = true;
	its_CurrentHedgesIndex = -1;
	
	its_AllNPVS.resize(its_CreditsLabels.size());
	its_AllDefLegs.resize(its_CreditsLabels.size());
	its_AllPremLegs.resize(its_CreditsLabels.size());

	// Price for Hedges
	UnSetHedgesDefault();		// just to be sure
	
	// ------------------------------
	// To be enhanced

		switch (its_Bump_Choice)
		{
		case CHB_SPREAD_PARALLEL:
			Hedges(dummy_NPV);
			break;

		case CHB_FAST_SPREAD_PARALLEL:
			// maybe temporary code
			HedgesFastSpreadParallel(dummy_NPV);
			break;
		}
	// ------------------------------	

	// Add one more loop with 'its_nb_scenarii' for ALL MATURITIES
	its_CurrentHedgesScenario	=	0;	// PARALLEL

	// Get All Prices
	for (i=0; i<its_CreditsLabels.size(); i++)
	{
		// in order to get 
		its_CurrentHedgesIndex =	i;

		// ------------------------------	
		switch (its_Bump_Choice)
		{
		case CHB_SPREAD_PARALLEL:
			// Set NPV
			Hedges(NPV_Shifted);

			// Delta NPV
			its_Hedges_Delta_NPV[i]	=	NPV_Shifted - Central_NPV;
			
			break;

		case CHB_FAST_SPREAD_PARALLEL:
			// Set NPV
			HedgesFastSpreadParallel(NPV_Shifted);

			// Delta NPV
			its_Hedges_Delta_NPV[i]	=	NPV_Shifted;

			break;
		}
		// ------------------------------	


		DiffCDS	=	its_Hedges_Delta_CDS[i];

		// Hedge Ratio
		if (fabs(DiffCDS) < 1e-12)
			HedgeRatio	=	0.0;
		else
			HedgeRatio	=	(NPV_Shifted - Central_NPV) / DiffCDS;
		its_Hedges_Hedge_Ratio[i]	=	HedgeRatio;

		// Delta CDS NPV
		its_Hedges_Delta_CDS[i]	=	DiffCDS;

		// ATM_Margin

		// Carry
		ATM_Margin	=	its_Hedges_CDS_ATM_Margin[i];
		its_Hedges_Carry[i]	=	HedgeRatio * ATM_Margin;	// is it in bps?
	}

	return;
}


void	CreditManager::HedgesFastSpreadParallel(double& NPV)
{
	switch (its_CreditModelType)
	{
	case CMT_MONTECARLO:
		HedgesFastSpreadParallel_Monte_Carlo(NPV);
		break;

	case CMT_ANALYTIC_RECURSIVE_STOCHASTIC_CORRELATION_BERNOULLI_1F:
	case CMT_ANALYTIC_RECURSIVE_INTERPOLATION_1F:
		HedgesFastSpreadParallel_Recursive(NPV);
		break;
	}
}


void	CreditManager::ComputeHedgeSpreadAllMaturities()
{
	return;
}

void	CreditManager::ComputeHedgeDefault()
{
	int i;

	double	dummy_NPV, NPV_Shifted;
	double	DiffCDS, HedgeRatio, ATM_Margin;

	// Reset All Outputs
	ResetAndResizeAllHedgesOutputs();

	//	CDS Prices
	//	Generates Shifted Default Curves
	//	CDS Shifted Prices
	PrepareHedgesDefault();

	// All Prices will be computed
	its_HedgesRunning = true;
	its_CurrentHedgesIndex = -1;
	
	its_AllNPVS.resize(its_CreditsLabels.size());
	its_AllDefLegs.resize(its_CreditsLabels.size());
	its_AllPremLegs.resize(its_CreditsLabels.size());

	// Price for Hedges
	SetHedgesDefault();
	Hedges(dummy_NPV);

	// Get All Prices
	for (i=0; i<its_CreditsLabels.size(); i++)
	{
		// in order to get 
		its_CurrentHedgesIndex =	i;

		// Set NPV
		Hedges(NPV_Shifted);

		// Delta NPV
		its_Hedges_Delta_NPV[i]	=	NPV_Shifted - Central_NPV;

		DiffCDS	=	its_Hedges_Delta_CDS[i];

		// Hedge Ratio
		if (fabs(DiffCDS) < 1e-12)
			HedgeRatio	=	0.0;
		else
			HedgeRatio	=	(NPV_Shifted - Central_NPV) / DiffCDS;
		its_Hedges_Hedge_Ratio[i]	=	HedgeRatio;

		// Delta CDS NPV
		its_Hedges_Delta_CDS[i]	=	DiffCDS;

		// ATM_Margin

		// Carry
		ATM_Margin	=	its_Hedges_CDS_ATM_Margin[i];
		its_Hedges_Carry[i]	=	HedgeRatio * ATM_Margin;	// is it in bps?
	}

	return;
}

void	CreditManager::ComputeHedgeRecoverySens()
{
	int i;

	double	dummy_NPV, NPV_Shifted;

	// Reset All Outputs
	ResetAndResizeAllHedgesOutputs();

	//	CDS Prices
	//	Generates Shifted Default Curves
	//	CDS Shifted Prices
	PrepareHedges();

	// All Prices will be computed
	its_HedgesRunning = true;
	its_CurrentHedgesIndex = -1;
	
	its_AllNPVS.resize(its_CreditsLabels.size());
	its_AllDefLegs.resize(its_CreditsLabels.size());
	its_AllPremLegs.resize(its_CreditsLabels.size());

	// Price for Hedges
	UnSetHedgesDefault();		// just to be sure
	Hedges(dummy_NPV);

	// Get All Prices
	for (i=0; i<its_CreditsLabels.size(); i++)
	{
		// in order to get 
		its_CurrentHedgesIndex =	i;

		// Set NPV
		Hedges(NPV_Shifted);

		// Delta NPV
		its_Hedges_Delta_NPV[i]	=	NPV_Shifted - Central_NPV;

		// Carry (ugly!)
		its_Hedges_Carry[i]	=	(NPV_Shifted - Central_NPV) / its_BumpRecovery;
	}

	return;
}

void	CreditManager::ComputeHedgeRecoveryLoss()
{
	int i;

	double	dummy_NPV, NPV_Shifted;

	// Reset All Outputs
	ResetAndResizeAllHedgesOutputs();

	Copy_MatrixShiftedDefaultCrv_From_DefaultCrv();

	//	CDS Prices
	//	Generates Shifted Default Curves
	//	CDS Shifted Prices

	// All Prices will be computed
	its_HedgesRunning = true;
	its_CurrentHedgesIndex = -1;
	
	its_AllNPVS.resize(its_CreditsLabels.size());
	its_AllDefLegs.resize(its_CreditsLabels.size());
	its_AllPremLegs.resize(its_CreditsLabels.size());

	// Price for Hedges
	UnSetHedgesDefault();		// just to be sure
	Hedges(dummy_NPV);

	// Get All Prices
	for (i=0; i<its_CreditsLabels.size(); i++)
	{
		// in order to get 
		its_CurrentHedgesIndex =	i;

		// Set NPV
		Hedges(NPV_Shifted);

		// Delta NPV
		its_Hedges_Delta_NPV[i]	=	NPV_Shifted - Central_NPV;

		// Carry (ugly!)
		its_Hedges_Carry[i]	=	(NPV_Shifted - Central_NPV) / its_BumpRecovery;
	}

	return;
}

void	CreditManager::ComputeHedgeRecoveryGlobal()
{
	double	NPV_Shifted;

	// Reset All Outputs
	ResetAndResizeAllHedgesOutputs();

	//	CDS Prices
	//	Generates Shifted Default Curves
	//	CDS Shifted Prices
	PrepareHedges();

	ICM_DefaultCurve**	Copy_ArrayDefaultCrv;
	Copy_ArrayDefaultCrv	= new ICM_DefaultCurve* [its_CreditsLabels.size()];

	int	i, iscenario;
	iscenario	=	0;

	for (i=0; i<its_CreditsLabels.size(); i++)
	{
		Copy_ArrayDefaultCrv[i]	=	its_ArrayDefaultCrv[i];	//(ICM_DefaultCurve*) its_ArrayDefaultCrv[i]->Clone();
		its_ArrayDefaultCrv[i]	=	(*its_MatrixShiftedDefaultCrv)(i, iscenario);
	}

	// All Prices will be computed
	its_HedgesRunning = true;
	its_CurrentHedgesIndex = -1;
	
	// ------------------------------------

	// A single output
	its_AllNPVS.resize(1);
	its_AllDefLegs.resize(1);
	its_AllPremLegs.resize(1);


	// ------------------------------------
	// LOSSES AMOUNT HAVE TO BE UPDATED

	double	TheLoss;
	DoubleVector	Kept_Losses;

	Kept_Losses.resize(its_CreditsLabels.size());

	for (i=0; i<its_CreditsLabels.size(); i++)
	{
		TheLoss			=	its_Input_Losses[i];
		Kept_Losses[i]	=	TheLoss;
		its_Input_Losses[i]	-=	its_BumpRecovery;
	}
	// ------------------------------------

	// Price will computes its_Input_Losses_Amount via ComputeBasketNotional
	Price();

	NPV_Shifted	=	NPV;

	// Delta NPV
	its_Hedges_Delta_NPV[0]	=	NPV_Shifted - Central_NPV;

	// Carry (ugly!)	// sensitivity
	its_Hedges_Carry[0]	=	(NPV_Shifted - Central_NPV) / its_BumpRecovery;

	// ------------------------------------
	// LOSSES AMOUNT HAVE TO BE RESTORED

	for (i=0; i<its_CreditsLabels.size(); i++)
		its_Input_Losses[i]	=	Kept_Losses[i];
	// ------------------------------------

	for (i=0; i<its_CreditsLabels.size(); i++)
		its_ArrayDefaultCrv[i]	=	Copy_ArrayDefaultCrv[i];

	delete Copy_ArrayDefaultCrv;
	Copy_ArrayDefaultCrv	=	NULL;
}


void	CreditManager::ComputeHedgeCorrelation()
{
	int i;

	double	dummy_NPV, NPV_Shifted;

	// Reset All Outputs
	ResetAndResizeAllHedgesOutputs();

	Copy_MatrixShiftedDefaultCrv_From_DefaultCrv();
	
	//	CDS Prices
	//	Generates Shifted Default Curves
	//	CDS Shifted Prices

	// All Prices will be computed
	its_HedgesRunning = true;
	its_CurrentHedgesIndex = -1;
	
	// ------------------------------------
	switch (its_CorrelationType)
	{
		case CT_FLAT:

			// A single output
			its_AllNPVS.resize(1);
			its_AllDefLegs.resize(1);
			its_AllPremLegs.resize(1);

			// Bump The Correlation
			if (fabs(its_CorrelationValue + its_BumpCorrelation) > 1.0)
				ICMTHROW(ERR_INVALID_DATA,"Shifted Correlation is out of range!");

			its_CorrelationValue	+=	its_BumpCorrelation;

			Price();

			NPV_Shifted	=	NPV;

			// Delta NPV
			its_Hedges_Delta_NPV[0]	=	NPV_Shifted - Central_NPV;

			// Carry (ugly!)	// sensitivity
			its_Hedges_Carry[0]	=	(NPV_Shifted - Central_NPV) / its_BumpCorrelation;

			// restore Correlation
			its_CorrelationValue	-=	its_BumpCorrelation;

			return;

		case CT_BETA:			
			break;

		case CT_MATRIX:
				ICMTHROW(ERR_INVALID_ARGUMENT,"Unable to compute Correlation Sensitivities with a Full Matrix!");
			break;
		case CT_FACTOR_LOADING_2:
				ICMTHROW(ERR_INVALID_ARGUMENT,"Not yet implemented! Correlation Sensitivities with a 2 Factor Loading!");
			break;
	}
	// ------------------------------------

	// ------------------------------------------
	// BETA cases
	// ------------------------------------------

	its_AllNPVS.resize(its_CreditsLabels.size());
	its_AllDefLegs.resize(its_CreditsLabels.size());
	its_AllPremLegs.resize(its_CreditsLabels.size());

	// Price for Hedges
	UnSetHedgesDefault();		// just to be sure
	Hedges(dummy_NPV);

	// Get All Prices
	double	TheSQRT_BumpCorrelation;

	TheSQRT_BumpCorrelation	=	sqrt(its_BumpCorrelation);

	for (i=0; i<its_CreditsLabels.size(); i++)
	{
		// Bump The Correlation (not the beta...)

		if (fabs(its_Beta[i] + TheSQRT_BumpCorrelation) > 1.0)
			ICMTHROW(ERR_INVALID_DATA,"Shifted Correlation is out of range!");

		its_Beta[i]	+=	TheSQRT_BumpCorrelation;

		Price();

		NPV_Shifted	=	NPV;

		// Delta NPV
		its_Hedges_Delta_NPV[i]	=	NPV_Shifted - Central_NPV;

		// Carry (ugly!)	// sensitivity
		its_Hedges_Carry[i]	=	(NPV_Shifted - Central_NPV) / its_BumpCorrelation;

		// restore Correlation
		its_Beta[i]	-=	TheSQRT_BumpCorrelation;

	}

	return;
}


void	CreditManager::PrepareHedges()
{
	double	ATM_Margin;
	double	CDS_NPV, CDS_NPV_Shifted;
	double	TheCDSNotional;

	DoubleVector	CDSDiffNPVS;
	
	CDSDiffNPVS.resize(its_CreditsLabels.size());

	const ICM_DefaultCurve*	TheDefaultCurve=NULL;
	const ICM_DefaultCurve*	TheShiftedDefaultCurve=NULL;
	
	int i;
	int	nocurve = -1;
	
	// -------------------------------------------
	ARM_Date	Maturity;				// 20/06/2010
	string		Tenor	=	"NULL";		// '5Y'
	char* plot = new char[10];

	// -------------------------------------------
	if (IsHedgesCDSToMaturity())
	{
		// its_MaxDate
		Maturity	=	itsValDate;
		Maturity.AddDays(its_MaxDate);
	}
	else
		strcpy(plot, its_HedgesCDSMaturity.c_str());
	// -------------------------------------------


	its_nb_scenarii	=	1;
	qSENSITIVITY_TYPE	shift_mode;
	double	epsilon;

	const ARM_Vector*	TmpVect;

	switch (its_Bump_Choice)
	{
	case CHB_SPREAD_ALL_MATURITIES:
		// get the number of scenarii from the first Default Curve (hypothesis, they all have the same plots)
		TheDefaultCurve	=	its_ModelMultiCurves->GetDefaultCurve(0);
		
		if (TheDefaultCurve == NULL)
		{
			ICMTHROW(ERR_INVALID_DATA,"NULL Default Curve in PrepareHedges!");
		}
		else
		{
			TmpVect	=	TheDefaultCurve->GetRates();
			if (TmpVect == NULL)
				ICMTHROW(ERR_INVALID_DATA,"NULL Default Curve Data in PrepareHedges!");

			its_nb_scenarii	=	TmpVect->GetSize();

			if (its_nb_scenarii <= 0)
				ICMTHROW(ERR_INVALID_DATA,"Bad Default Curve Data in PrepareHedges, number of points <= 0!");
		}

		break;

	case CHB_SPREAD_PARALLEL:
		// PARALLEL SHIFT
		its_nb_scenarii	=	1;

		// discriminate additive or multiplicative shift 
		if (its_BumpSpread_Type == BT_ADD)
		{
			shift_mode	=	ICMSPREAD_TYPE;			
			epsilon		=	its_BumpSpread / 100.0;	// in pct... s = s + epsilon / 100 ; epsilon = 10 for a 10 additive bump
		}
		else // BT_MULT
		{
			shift_mode	=	ICMSPRELSHIFT_TYPE;			
			epsilon		=	its_BumpSpread;	// in pct... s = s * (1+epsilon) ; epsilon = 10 for a 10% multiplicative bump
		}
		strcpy(plot,"NONE");			// parallel shift

		break;

	case CHB_DEFAULT:

		PrepareHedgesDefault();
		if (plot)
			delete[] plot;
		plot = NULL;
		return;

		break;
		
	case CHB_RECOVERY_SENS:
		// RECOVERY ON BOTH CALIBRATION AND PRICING

		shift_mode	=	ICMRECOVERY_TYPE;
		epsilon		=	its_BumpRecovery;
		strcpy(plot,"NONE");			// parallel shift

		break;
		
	case CHB_RECOVERY_LOSS:
		break;
		
	case CHB_RECOVERY_GLOBAL:
		// RECOVERY ON BOTH CALIBRATION AND PRICING

		shift_mode	=	ICMRECOVERY_TYPE;
		epsilon		=	its_BumpRecovery;
		strcpy(plot,"NONE");			// parallel shift

		break;

	case CHB_CORRELATION:
		break;

	case CHB_NO:
		if (plot)
			delete[] plot;
		plot = NULL;
		return;
	}

	// -------------------------------------------

	if (its_ModelMultiCurves == NULL)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
			"Parameters:  Model Multi Curve is NULL!");

	// -------------------------------------------------------------------------
	// Bump a package of Default Curves according to a given bump profile
	// -------------------------------------------------------------------------
	ICM_ModelMultiCurves* ModelDef_Shifted = its_ModelMultiCurves->GenerateShiftModel(
																	shift_mode,
																	plot, 
																	"NONE",		// all curves
																	nocurve,	// useless
																	epsilon);

	// ------------------------------------------------------------
	// Default Curves
	its_MatrixShiftedDefaultCrv = new	ICM_QMatrix<ICM_DefaultCurve*>(its_CreditsLabels.size(), its_nb_scenarii);
	// ------------------------------------------------------------

	if (its_Bump_Choice != CHB_CORRELATION)
	{
		double	DefLeg	=	0.0;
		double	RPV01	=	0.0;

		int	iscenario	=	0;

		for (i=0; i<its_CreditsLabels.size(); i++)
		{
			//	Retrieve Default Curve from Id
			TheDefaultCurve	=	its_ModelMultiCurves->GetDefaultCurve(i);
			if (TheDefaultCurve == NULL)
				ICMTHROW(ERR_INVALID_DATA,"Parameters:  The Default Curve is NULL!");

			// CDS_NPV
			// waiting for Damien's delivery
			CDS_NPV	=	0.0;

			// ATM_Margin in bps. i.e. 50 for 0.5%
			ATM_Margin	=	TheDefaultCurve->NPV_Implicit_Cds(Maturity, Tenor, qCMPSPREAD);
	//		ATM_Margin	=	0.0;

			// Now deal with shifted Curve
			TheShiftedDefaultCurve	=	ModelDef_Shifted->GetDefaultCurve(i);
			if (TheShiftedDefaultCurve == NULL)
				ICMTHROW(ERR_INVALID_DATA,"Parameters:  The Shifted Default Curve is NULL!");

			// store it (to be clones?)
			(*its_MatrixShiftedDefaultCrv)(i, iscenario)	=	(ICM_DefaultCurve*) TheShiftedDefaultCurve->Clone();

			// CDS_NPV_Shifted
			// waiting for Damien's delivery
			DefLeg	=	TheDefaultCurve->NPV_Implicit_Cds(Maturity, Tenor, qCMPDEFLEGPV);
			RPV01	=	TheDefaultCurve->NPV_Implicit_Cds(Maturity, Tenor, qCMPDURATION);
			TheCDSNotional	=	TheDefaultCurve->GetNotional();

			CDS_NPV_Shifted	=	DefLeg - ATM_Margin * RPV01 / 10000.0;
			CDS_NPV_Shifted	/=	TheCDSNotional;


			// Set in vectors
			its_Hedges_Delta_CDS[i]			=	CDS_NPV_Shifted - CDS_NPV;
			its_Hedges_CDS_ATM_Margin[i]	=	ATM_Margin;
		}
	}

	if (plot)
		delete[] plot;
	plot = NULL;
}


void	CreditManager::Copy_MatrixShiftedDefaultCrv_From_DefaultCrv()
{
	// ------------------------------------------------------------
	// Default Curves
	its_MatrixShiftedDefaultCrv = new	ICM_QMatrix<ICM_DefaultCurve*>(its_CreditsLabels.size(), 1);
	// ------------------------------------------------------------

	for (int i=0; i<its_CreditsLabels.size(); i++)
	{
		// store it (to be clones?)
		(*its_MatrixShiftedDefaultCrv)(i, 0)	=	(ICM_DefaultCurve*) ((its_ArrayDefaultCrv[i])->Clone());
	}
}


// ------------------------------------------
// much much simplier
// ------------------------------------------
void	CreditManager::PrepareHedgesDefault()
{
	if (its_Bump_Choice != CHB_DEFAULT)
	{
		PrepareHedges();
		return;
	}

	double	TheNotional	=	1.0;
	double	ATM_Margin;
	double	CDS_NPV, CDS_NPV_Shifted;

	DoubleVector	CDSDiffNPVS;
	
	CDSDiffNPVS.resize(its_CreditsLabels.size());

	const ICM_DefaultCurve*	TheDefaultCurve=NULL;
	const ICM_DefaultCurve*	TheShiftedDefaultCurve=NULL;
	
	int i;
	int	nocurve = -1;
	
	// -------------------------------------------
	ARM_Date	Maturity;				// 20/06/2010
	string		Tenor	=	"NULL";		// '5Y'
	char* plot = new char[10];

	// -------------------------------------------
	if (IsHedgesCDSToMaturity())
	{
		// its_MaxDate
		Maturity	=	itsValDate;
		Maturity.AddDays(its_MaxDate);
	}
	else
		strcpy(plot, its_HedgesCDSMaturity.c_str());
	// -------------------------------------------

	if (its_ModelMultiCurves == NULL)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
			"Parameters:  Model Multi Curve is NULL!");

	for (i=0; i<its_CreditsLabels.size(); i++)
	{
		//	Retrieve Default Curve from Id
		TheDefaultCurve	=	its_ModelMultiCurves->GetDefaultCurve(i);
		if (TheDefaultCurve == NULL)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
				"Parameters:  The Default Curve is NULL!");

		// CDS_NPV
		// waiting for Damien's delivery
		CDS_NPV	=	0.0;

		// ATM_Margin
		ATM_Margin	=	TheDefaultCurve->NPV_Implicit_Cds(Maturity, Tenor, qCMPSPREAD);
//		ATM_Margin	=	0.0;

		// CDS_NPV_Shifted
		CDS_NPV_Shifted	=	(1.0 - TheDefaultCurve->GetRecovery()) * TheNotional;

		// Set in vectors
		its_Hedges_Delta_CDS[i]			=	CDS_NPV_Shifted - CDS_NPV;
		its_Hedges_CDS_ATM_Margin[i]	=	ATM_Margin;
	}

}


void	CreditManager::	ResetAndResizeAllHedgesOutputs()
{
	its_Hedges_Delta_NPV.clear();
	its_Hedges_Hedge_Ratio.clear();
	its_Hedges_Delta_CDS.clear();
	its_Hedges_CDS_ATM_Margin.clear();
	its_Hedges_Carry.clear();

	its_Hedges_MatrixOutput.clear();

	its_Hedges_Delta_NPV.resize(its_CreditsLabels.size());
	its_Hedges_Hedge_Ratio.resize(its_CreditsLabels.size());
	its_Hedges_Delta_CDS.resize(its_CreditsLabels.size());
	its_Hedges_CDS_ATM_Margin.resize(its_CreditsLabels.size());
	its_Hedges_Carry.resize(its_CreditsLabels.size());

//	its_Hedges_MatrixOutput.resize(its_CreditsLabels.size(),);

}


void	CreditManager::Display_DefaultTimes()
{
	int	i;
	DoubleVector::iterator iter;
	DoubleVector::iterator iter_shifted;

	if (its_Display_DefaultTimes_Flag)
	{
		// ----------------------------------------------------------------------------
		if ((its_fOut = fopen("c:\\test\\mc_Display_DefaultTimes.txt", "a+")) == NULL) return;
		fprintf(its_fOut, " ----------------- Display_DefaultTimes ----------------- \n");
		// ----------------------------------------------------------------------------
		i	=	0;

		if (its_DefaultTimesShifted.empty())
			for (iter = its_DefaultTimes.begin();
						iter != its_DefaultTimes.end();
						++iter,i++)
			
				fprintf(its_fOut, "Id:\t%u\t\ttau:\t%.2lf\n", i, *iter);
			
		else
			for (iter = its_DefaultTimes.begin(), iter_shifted = its_DefaultTimesShifted.begin();
						iter != its_DefaultTimes.end();
						++iter,i++,++iter_shifted)
			
				fprintf(its_fOut, "Id:\t%u\t\ttau:\t%.2lf\t\ttau shifted:\t%.2lf\n", i, *iter, *iter_shifted);
			

		fclose(its_fOut);
	}

}


void	CreditManager::Display_Vector(string& label, DoubleVector&	data)
{
	int	i;
	DoubleVector::iterator iter;

	// ----------------------------------------------------------------------------
	if ((its_fOut = fopen("c:\\test\\mc_Display_Vector.txt", "a+")) == NULL) return;
	fprintf(its_fOut, " ----------------- Display_Vector %s----------------- \n", label.c_str());
	// ----------------------------------------------------------------------------
	
	i =	0;
	for (iter = data.begin(); iter != data.end(); ++iter,i++)	
		fprintf(its_fOut, "Id:\t%u\t\t:\t%.6lf\n", i, *iter);		

	fclose(its_fOut);

}


// -------------------------------------------------------
// VARIANCE REDUCTION
// -------------------------------------------------------
void	CreditManager::ComputeVarianceReductionData()
{
	int		j;
	double	TheCorrelation;
	double	TheValue, Target_Inv_Prob;

	if (its_CreditModelType == CMT_MONTECARLO)
	{
		switch (its_MC_Variance_Reduction)
		{
			case CMCVR_IS_FACTORS:
			case CMCVR_IS_PURE_FACTORS:

				// mu computation
				its_IS_Mu	=	0.0;

				// N-1(its_Target_Default_Prob)				
				if (CHECK_EQUAL(its_Target_Default_Prob, 0.0))
					Target_Inv_Prob	=	_MINUS_INFINITY_;	//	minus infinity
				else if (CHECK_EQUAL(its_Target_Default_Prob, 1.0))
					Target_Inv_Prob	=	_PLUS_INFINITY_;	//	plus infinity
				else
				{
					Target_Inv_Prob =	its_CopulaChoice->Inverse_Cumulative_Density_Function(its_Target_Default_Prob);
					Target_Inv_Prob	=	NAG_deviates_normal(its_Target_Default_Prob );
				}
			
				// N-1(DefaultProb at Maturity)
				for (j=0; j<its_CreditsLabels.size(); j++)
				{
					switch (its_CorrelationType)
					{
						case CT_FLAT:
							TheCorrelation	=	its_CorrelationValue;
							break;
						
						case CT_BETA:
							TheCorrelation	=	its_Beta[j] * its_Beta[j];
							break;

						default:
							ICMTHROW(ERR_INVALID_DATA,"Variance Reduction is only available with Flat or Beta Correlation Structures!");

							break;
					}

					if (CHECK_EQUAL(TheCorrelation, 0.0))
						TheValue	=	0.0;
					else
					{
						TheValue	= itsBarriers_Standard[j];
						TheValue	-= sqrt(1.0 - TheCorrelation) * Target_Inv_Prob;
						TheValue	/=	sqrt(TheCorrelation);
					}

					its_IS_Mu	+=	TheValue;
				}
				its_IS_Mu	/=	its_CreditsLabels.size();
				
			// does it work?
			if (its_MC_Variance_Reduction == CMCVR_IS_PURE_FACTORS)
				break;

			case CMCVR_IS_IDIOSYNCRATIC:

				IS_individual_exp_twists.clear();
				IS_individual_exp_twists_Kept.clear();

				if (its_Theta_Choice == CVRI_OPTIM)
				{
					// completely ad-hoc choice
					its_IS_Theta	=	1.0;
					its_Kept_Theta	=	1.0;
				}
				else
					its_Kept_Theta	=	its_IS_Theta;

				// Normalization in order to keep a tractable Theta (always the same order)
				for (j=0; j<its_CreditsLabels.size(); j++)
				{
					IS_individual_exp_twists.push_back(exp(its_IS_Theta * its_LossesAmount[j] / TheBasketNotional));
					IS_individual_exp_twists_Kept.push_back(exp(its_IS_Theta * its_LossesAmount[j] / TheBasketNotional));
				}
				break;

			case CMCVR_NONE:
				break;
			
			default:
				break;
		}
	}
}


// -------------------------------------------------------
// STRUCTURE ALLOCATION FOR 1F or Monte-Carlo
// -------------------------------------------------------
void	CreditManager::AllocateStructures_1F()
{
	ICM_Integrator	TheIntegrator;

	ARM_Vector LHP_Spreads;
	ARM_Vector LHP_Maturities_AsYF;

	// DEFAULT CURVE, parameters
	long	PayFreq		=	K_QUARTERLY;
	qDEFCURVE_CALIB_ALGO algo =	qDEFCURVE_BRENT;
	// Roll Dates
//	qCDS_ADJ	CDS_Roll_Type;
//	CDS_Roll_Type	=	(its_RollDateFlag) ? qCredit_Adjust20 : qCredit_Default;

	switch (its_CreditModelType)
	{
	case CMT_MONTECARLO:

		break;

	case CMT_ANALYTIC_RECURSIVE_1F:
	case CMT_ANALYTIC_RECURSIVE_INTERPOLATION_1F:
	case CMT_ANALYTIC_RECURSIVE_STOCHASTIC_CORRELATION_BERNOULLI_1F:

		// ProbCond: a QCube (cf. required only for hedges purposes)
		// LossDistrib: a vector

		if (its_NIntegration_1F <=0)
			ICMTHROW(ERR_INVALID_DATA,"Numerical integration number must be strictly positive! Current value: " << its_NIntegration_1F);

		if (its_ProbCond) delete its_ProbCond;
		its_ProbCond	=	new ICM_QCubix<double>(its_CreditsLabels.size()+1, its_NIntegration_1F, 1, 0.);

		its_lossdistrib.clear();
		its_lossdistrib.resize(1);

	case CMT_ANALYTIC_LARGE_PORTFOLIO_1F:

		// Integration
		its_Xi.resize(its_NIntegration_1F+1);
		its_Wi.resize(its_NIntegration_1F+1);

		// Gauss-Legendre Integration
		// TO BE DIVERSIFIED!!!
		TheIntegrator.ComputeAbscissasAndWeightsGaussLegendre(0.0, 1.0, its_NIntegration_1F);
		
		TheIntegrator.GetAbscissaVector(its_Xi);
		TheIntegrator.GetWeightVector(its_Wi);

		break;

	case CMT_ANALYTIC_LHP_JPM:

		// Integration
		its_Xi.resize(its_NIntegration_1F+1);
		its_Wi.resize(its_NIntegration_1F+1);

		// Gauss-Legendre Integration
		// TO BE DIVERSIFIED!!!
		TheIntegrator.ComputeAbscissasAndWeightsGaussLegendre(0.0, 1.0, its_NIntegration_1F);
		
		TheIntegrator.GetAbscissaVector(its_Xi);
		TheIntegrator.GetWeightVector(its_Wi);

	case CMT_ANALYTIC_LHP:
	case CMT_ANALYTIC_LHP_PLUS:	// needs that part + the Heterogeneous part
		{
		// FLAT CORRELATION
		if (its_CorrelationType != CT_FLAT)
			ICMTHROW(ERR_INVALID_DATA,"FLAT CORRELATION EXPECTED for LHP Model!");

		// Roll Dates
		qCDS_ADJ	CDS_Roll_Type;
		CDS_Roll_Type	=	(its_RollDateFlag) ? qCredit_Adjust20 : qCredit_Default;

		LHP_Spreads.push_back(its_LHP_Spread / 10000.0);		// in bps.
		LHP_Maturities_AsYF.push_back(MATHTIME(its_LHP_Maturity));
		ARM_Security* psec = GetSecurity();
		string currency("");
		if (psec == NULL) currency = ARM_DEFAULT_COUNTRY;
		else currency = string(GetSecurity()->GetCurrencyUnit()->GetCcyName());
		// DEFAULT CURVE
		its_DefCurve = new ICM_Constant_Piecewise(
												itsValDate,
												LHP_Maturities_AsYF,
												LHP_Spreads,
												its_LHP_Recovery,
												its_ZeroCurve,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
												qCredit_Default,//qCredit_Special_None_YF,  useless
												currency, //ARM_DEFAULT_COUNTRY,	//	ccy
												"DUMMY",				//	label
												false,					//	issummitcurve
												NULL,					//	Vol Curve
												false,
												//2 NULL,
												PayFreq,
												algo,
												"STD",
												ARM_Currency(currency.c_str()).GetCreditStartDateLag());					

		}
		break;

	case CMT_ANALYTIC_FFT_1F:
		break;

	}

	AllocateCopula_1F();

	// CORRELATION
	int	j;
	double	tmpValue;

	if (its_CreditModelType == CMT_ANALYTIC_RECURSIVE_STOCHASTIC_CORRELATION_BERNOULLI_1F)
		Allocation_Stochastic_Correlation_Data();
	else
	{
		switch (its_CorrelationType)
		{
			case CT_FLAT:
				its_Used_Beta.resize(its_CreditsLabels.size());
				its_Used_SQRT_OneMinusBetaSquare.resize(its_CreditsLabels.size());

				tmpValue	=	sqrt(1.0 - its_CorrelationValue);
				fill(its_Used_Beta.begin(), its_Used_Beta.end(), sqrt(its_CorrelationValue));
				fill(its_Used_SQRT_OneMinusBetaSquare.begin(), its_Used_SQRT_OneMinusBetaSquare.end(), tmpValue);

				break;

			case CT_BETA:			
				its_Used_Beta.resize(its_CreditsLabels.size());
				its_Used_SQRT_OneMinusBetaSquare.resize(its_CreditsLabels.size());

				for (j=0; j<its_CreditsLabels.size(); j++)
				{
					tmpValue			=	its_Beta[j];
					its_Used_Beta[j]	=	tmpValue;
					its_Used_SQRT_OneMinusBetaSquare[j]	=	sqrt(1.0 - tmpValue * tmpValue);
				}

				break;

			case CT_MATRIX:
				break;

			case CT_FACTOR_LOADING_2:
				break;
		}
	}
}


// -------------------------------------------------------
// COPULA ALLOCATION FOR 1F or Monte-Carlo
// -------------------------------------------------------
void	CreditManager::AllocateCopula_1F()
{
	switch (its_CopulaType)
	{
	case CCT_GAUSSIAN:
		its_CopulaChoice	=	new ICM_Probability_Density_Standard_Gaussian();
		break;
		
	case CCT_STUDENT:
		its_CopulaChoice	=	new ICM_Probability_Density_Student(its_FreedomDegree);
		break;
		
	case CCT_NIG:
		// to do
		break;
	}
}


// ----------------------------------------------------
// HEDGE FOR ANALYTIC LHP PLUS
// ----------------------------------------------------

void	CreditManager::Hedges_LHP(double& dummyNPV)
{
	if (!its_HedgesRunning) return;

	int	i;
	double	PremLeg, DefLeg;
	DoubleVector	Outputs;
	
	// First passage
	if (its_CurrentHedgesIndex	==	-1)
	{
		fill(its_AllNPVS.begin(), its_AllNPVS.end(), 0.0);
		fill(its_AllDefLegs.begin(), its_AllDefLegs.end(), 0.0);
		fill(its_AllPremLegs.begin(), its_AllPremLegs.end(), 0.0);

		// -------------------------------------------------------
		// The model will have to be set to LHP PLUS...
		// -------------------------------------------------------
		if (its_CreditModelType	!= CMT_ANALYTIC_LHP)
			ICMTHROW(ERR_INVALID_DATA,"LHP Hedges, only with LHP Model! Be consistent.");

		its_CreditModelType	=	CMT_ANALYTIC_LHP_PLUS;

		try
		{
			for (i=0;i<its_CreditsLabels.size();i++)
			{
				// -------------------------------------------------------
				// Reset Outputs
				// -------------------------------------------------------
				ResetOutputs();	// NPV to 0.0
				
				// Pricing
				PriceBasket_1F(Outputs);

				// Outputs
				its_AllNPVS[i]		=	Outputs[0];
				its_AllDefLegs[i]	=	Outputs[1];
				its_AllPremLegs[i]	=	Outputs[2];
			}
		}

		catch(...)
		{
			// restore
			its_CreditModelType	= CMT_ANALYTIC_LHP;
		}

		dummyNPV	=	0.0;

	}
	else
	{
		// simply get the right value
		dummyNPV	=	its_AllNPVS[its_CurrentHedgesIndex];

		PremLeg	=	its_AllPremLegs[its_CurrentHedgesIndex];

		DefLeg	=	its_AllDefLegs[its_CurrentHedgesIndex];
	}

}


void	CreditManager::Set_CorrelationCalibrator(ICM_Matrix<ARM_Vector>* Matrix_Description, ICM_Matrix<ARM_Vector>* Matrix_Parameters, vector<ICM_Mez*>& Vector_Mezz)
{
	if (Matrix_Description == NULL) return;
	
	ARM_Vector* THE_VECT=NULL;
	int	i, size;

	DoubleVector	The_Attachment;
	DoubleVector	The_Detachment;
	DoubleVector	The_Compound_Correlation;
	DoubleVector	The_Market_UpFrontPremiums;
	DoubleVector	The_Market_Spreads;

	// ------------------------------------------------------
	// ATTACHMENT, DETACHMENT, COMPOUND CORRELATION:	VECTOR OF DOUBLE

	THE_VECT	= Matrix_Description->GetColVect("CB_ATTACHMENT");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Matrix_Description :  no CB_ATTACHMENT vector found");

	size =	THE_VECT->GetSize();
	
	The_Attachment.resize(size);

	for (i=0;i<size;i++)
		The_Attachment[i]	=	(double) (*THE_VECT)[i];		

	THE_VECT	= Matrix_Description->GetColVect("CB_DETACHMENT");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Matrix_Description :  no CB_DETACHMENT vector found");

	The_Detachment.resize(size);

	for (i=0;i<size;i++)
		The_Detachment[i]	=	(double) (*THE_VECT)[i];		

	THE_VECT	= Matrix_Description->GetColVect("CB_COMPOUND_CORRELATION");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Matrix_Description :  no CB_COMPOUND_CORRELATION vector found");

	The_Compound_Correlation.resize(size);

	for (i=0;i<size;i++)
		The_Compound_Correlation[i]	=	(double) (*THE_VECT)[i];		

	THE_VECT	= Matrix_Description->GetColVect("CB_MARKET_UPFRONT_PREMIUM");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Matrix_Description :  no CB_MARKET_UPFRONT_PREMIUM vector found");

	The_Market_UpFrontPremiums.resize(size);

	for (i=0;i<size;i++)
		The_Market_UpFrontPremiums[i]	=	(double) (*THE_VECT)[i];		

	THE_VECT	= Matrix_Description->GetColVect("CB_MARKET_SPREADS");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Matrix_Description :  CB_MARKET_SPREADS vector found");

	The_Market_Spreads.resize(size);

	for (i=0;i<size;i++)
		The_Market_Spreads[i]	=	(double) (*THE_VECT)[i];		
	// ------------------------------------------------------
	
	if (Matrix_Parameters == NULL) return;
	
	int TheValue;

	// ------------------------------------------------------
	// MODEL_TYPE:	INT --> ENUM
	CreditManager_Correlation_Calibration	The_Correlation_Calibration_Type;

	THE_VECT	= Matrix_Parameters->GetColVect("CORREL_CALIBRATION_TYPE");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad MODEL_TYPE");

	TheValue = (int) THE_VECT->Elt(0);	
	if ((TheValue < CMCC_FROM_LHP_COMPOUND_TO_LHP_BASE) || (TheValue > CMCC_FROM_LHP_JPM_COMPOUND_TO_LHP_JPM_NIG)) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CORREL_CALIBRATION_TYPE");
	
	The_Correlation_Calibration_Type = (CreditManager_Correlation_Calibration) TheValue;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// MODEL_TYPE:	INT --> ENUM
	CreditManager_Correlation_Calibration_Input	The_Correlation_Input_Type;

	THE_VECT	= Matrix_Parameters->GetColVect("CB_INPUT_TYPE");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CB_INPUT_TYPE");

	TheValue = (int) THE_VECT->Elt(0);	
	if ((TheValue < CMCCI_COMPOUND_CORRELATION) || (TheValue > CMCCI_SPREADS)) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CB_INPUT_TYPE");
	
	The_Correlation_Input_Type = (CreditManager_Correlation_Calibration_Input) TheValue;
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CB_LHP_SPREAD:	DOUBLE

	THE_VECT	= Matrix_Parameters->GetColVect("CB_LHP_SPREAD");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CB_LHP_SPREAD");

	double	CB_LHP_Spread	= (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// LHP_MATURITY:	DOUBLE

	THE_VECT	= Matrix_Parameters->GetColVect("CB_LHP_MATURITY");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CB_LHP_MATURITY");

	RelativeDate	CB_LHP_Maturity = (RelativeDate) (THE_VECT->Elt(0) - itsValDateAsDouble);	//	relative to AsOf
	// ------------------------------------------------------

	// ------------------------------------------------------
	// CB_LHP_RECOVERY:	DOUBLE

	THE_VECT	= Matrix_Parameters->GetColVect("CB_LHP_RECOVERY");
	if (!THE_VECT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "Parameters :  bad CB_LHP_RECOVERY");

	double	CB_LHP_Recovery	= (double) THE_VECT->Elt(0);
	// ------------------------------------------------------

	// ------------------------------------------------------
	// PARAMETERS
	// ------------------------------------------------------
	if (its_Correlation_Calibrator)
		delete its_Correlation_Calibrator;

	its_Correlation_Calibrator	=	new ICM_Credit_Manager_Calibrator();


	// ATTACHEMENT POINTS
	its_Correlation_Calibrator->Set_Attach_Points(The_Attachment);
	its_Correlation_Calibrator->Set_Detach_Points(The_Detachment);

	its_Correlation_Calibrator->Set_Market_Correlations(The_Compound_Correlation);
	its_Correlation_Calibrator->Set_Market_UpFrontPremiums(The_Market_UpFrontPremiums);
	its_Correlation_Calibrator->Set_Market_Spreads(The_Market_Spreads);

	its_Correlation_Calibrator->Set_Calibration_Type(The_Correlation_Calibration_Type);
	its_Correlation_Calibrator->Set_Input_Type(The_Correlation_Input_Type);

	its_Correlation_Calibrator->Set_LHP_Recovery(CB_LHP_Recovery);
	its_Correlation_Calibrator->Set_LHP_Spread(CB_LHP_Spread);
	its_Correlation_Calibrator->Set_LHP_Maturity(CB_LHP_Maturity);

	// Mezz Objects
	its_Correlation_Calibrator->Set_Credit_Tranches(Vector_Mezz);
}


void	CreditManager::Correlation_Calibrate()
{
	its_Correlation_Calibrator->Set_CreditManager(this);
	its_Correlation_Calibrator->Calibrate();
}


void	CreditManager::Compute_Cash_And_Accrued()
{
	int	i;

	double	TheCash;
	double	TheAccrued;

	TheCash		=	0.0;
	TheAccrued	=	0.0;

	if (its_PricingLegsType != CBN_DEFAULTLEGONLY)
	{
		for (i=0;i<its_PL_NbFlows;i++)
		{
			// Flows are paid if at the Credit Fixing Date, the event has occurred
			// using 'static' variables
			its_PL_StartDate	=	its_PL_StartDates[i];
			its_PL_EndDate		=	its_PL_EndDates[i];

			// Ratio
			its_PL_Ratio			=	its_PL_Ratios[i];
				
			// LOSS OBSERVATIONS
			its_PL_Notio			=	its_PL_Notios[i];
			its_PL_Spread			=	its_PL_Spreads[i];
			
			// ---------------------------------------------------
			// is there Cash or Accrued
			// ---------------------------------------------------
			if (its_PL_StartDate < 0)
			{
				if (its_PL_EndDate <= 0.0)
					TheCash		+=	its_PL_Spread * its_PL_Ratio * its_PL_Notio;
				else
				{
					TheAccrued	=	its_PL_Spread * its_PL_Ratio * its_PL_Notio * its_PL_StartDate / (its_PL_StartDate - its_PL_EndDate);
					// go out
					break;
				}
			}
		}
	}

	// Cash
	its_Cash	=	TheCash;
	// Accrued
	its_Accrued	=	TheAccrued;

}



// *******************************************************************************************
// Set method
//*******************************************************************************************

void CreditManager::Set(ARM_Security *option, ARM_Model *mod, const ICM_Parameters& parameters,const ARM_Date&asof)
{
	SetName(ICM_CREDIT_MANAGER);

	// ------------------------------------------------------------
	// PRODUCT
	if (option != NULL)
	{
		its_Credit_Product		=	(ICM_Customized_CDO*) option->Clone();

		// Retrieve from Stack
		// Set Cash Flows Matrix
		its_Credit_Product->FillCreditProductPricingParametersMatrix();	// do it first, because of ValDate
		its_Credit_Product->FillCreditProductDefaultMatrix();
		its_Credit_Product->FillCreditProductPremiumMatrix();

		SetSecurity(its_Credit_Product);
	}
	// ------------------------------------------------------------


	// ------------------------------------------------------------
	// MODEL
	ARM_CLASS_NAME	name;

	if (mod != NULL)
	{
		name = mod->GetName();
		if (name == ICM_MODELMULTICURVES)
		{
			Fill_DataFrom_ModelMultiCurves(option, mod);
		}
		else	//  ICM_CUSTOMIZED_CREDIT_MULTICURVES
		{
			its_Credit_Market_Data	=	(ICM_Customized_Credit_MultiCurves*) mod->Clone();

			// Retrieve from Stack
			its_Credit_Market_Data->FillCreditDataParameters();
			its_Credit_Market_Data->FillCreditDataDescription();
		}
	}
	// ------------------------------------------------------------

	// ------------------------------------------------------------
	// Parameters:
	// for Hedges purposes: SetCreditDataParameters
	ICM_Matrix<ARM_Vector>* tmp_Parameters = NULL;

	// if (parameters != NULL)
	if (!parameters .empty())
	{
		tmp_Parameters	=	unconst(parameters).GetDblParams();
		
		SetCreditModelParameters(tmp_Parameters);
	}
	// ------------------------------------------------------------

	// ------------------------------------------------------------
	// PRODUCT
	// ------------------------------------------------------------
	GetDataFromCreditProduct();

	// ------------------------------------------------------------
	// ZERO CURVE
	// ------------------------------------------------------------
	its_ZeroCurve	=	(ARM_ZeroCurve*) (((ICM_ModelMultiCurves*) mod)->GetZeroCurve()->Clone());

	// ------------------------------------------------------------
	// MODEL
	// ------------------------------------------------------------
	GetDataFromCreditMultiCurves();

	// ------------------------------------------------------------
	// CORRELATION
	// ------------------------------------------------------------

	bool	PriceToDo	=	true;

	double	New_KMin;
	double	New_KMax;
	double	Beta_KMin;
	double	Beta_KMax;
	double	Kept_DL_LossMin;
	double	Kept_DL_LossMax;

	double	Kept_NPV;
	double	Kept_Central_NPV;
	double	Kept_PremiumLegPV;
	double	Kept_DefaultLegPV;
	double	Kept_ATMPremiumLegWithoutNotio;
	double	Kept_NPVSquare;
	double	Kept_ATMPremiumLeg;

	double	UpFront_PV;

	ICM_Correlation*	The_Correlation	=	NULL;
	ARM_Vector*		Betas_ARM_Vector	=	NULL;

	ICM_QMatrix<double>*	The_CorrelationMatrix	=	NULL;

	switch (its_CorrelationType)
	{
	case CT_FLAT:
	case CT_BETA:

		The_Correlation	=	((ICM_Beta_Correlation*) ((ICM_Customized_Credit_MultiCurves*) mod)->GetCorrelation());

		// just two options
		double	TheValue;

		TheValue	=	((ICM_Beta_Correlation*) The_Correlation)->GetFixedBeta();
		
		if (TheValue != CREDIT_DEFAULT_VALUE)
		{
			its_CorrelationType		=	CT_FLAT;
			its_CorrelationValue	=	TheValue * TheValue;
		}
		else
		{
			its_CorrelationType		=	CT_BETA;
			
			its_Beta.clear();
			its_Beta.resize(its_CreditsLabels.size());

			Betas_ARM_Vector	=	((ICM_Beta_Correlation*) The_Correlation)->GetBetas();

			if (Betas_ARM_Vector->GetSize() != its_CreditsLabels.size())
				ICMTHROW(ERR_INVALID_DATA,"Bad Size for Betas from ICM Correlation Object!");

			for (int i=0; i<its_CreditsLabels.size(); i++)
				its_Beta[i]	=	(*Betas_ARM_Vector)[i];

		}

		break;

	case CT_MATRIX:
	
		// MATRIX case
		its_CorrelationType	=	CT_MATRIX;

		The_Correlation	=	((ICM_CorrMatrix*) ((ICM_Customized_Credit_MultiCurves*) mod)->GetCorrelation());

		if (its_CorrelationMatrix)
			delete	its_CorrelationMatrix;

		its_CorrelationMatrix	=	(ICM_CorrMatrix*) (The_Correlation->Clone());

		break;

	case CT_FACTOR_LOADING_2:

		// TO DO
		ICMTHROW(ERR_INVALID_DATA,"Factor Loading Dim 2 not done yet!");

		break;
	case CT_BASE_STRIKE:

		// CT_BASE_STRIKE case
		The_Correlation	=	((ICM_Smile_Correlation*) ((ICM_Customized_Credit_MultiCurves*) mod)->GetCorrelation());

		PriceToDo	=	false;

		if (its_PricingLegsType != CBN_PREMIUMLEGONLY)
		{
			New_KMin	=	its_DL_LossMin;
			New_KMax	=	its_DL_LossMax;
		}
		else
		{
			// take some Average...
			New_KMin	=	accumulate(its_PL_LossMins.begin(), its_PL_LossMins.end(), 0.0);
			New_KMin	/=	its_PL_NbFlows;
			New_KMax	=	accumulate(its_PL_LossMaxs.begin(), its_PL_LossMaxs.end(), 0.0);
			New_KMax	/=	its_PL_NbFlows;
		}

		// Get Correlation for these two new Strikes
		// in order to get its_MaxCreditDate
		ComputeMaxCreditDate();

		Beta_KMin	=	((ICM_Smile_Correlation*) The_Correlation)->GetBeta("", MATHTIME(its_MaxCreditDate), New_KMin);
		Beta_KMax	=	((ICM_Smile_Correlation*) The_Correlation)->GetBeta("", MATHTIME(its_MaxCreditDate), New_KMax);

		DoubleVector	Kept_PL_LossMins;
		DoubleVector	Kept_PL_LossMaxs;
		DoubleVector	Kept_PL_Notios;

		Kept_PL_LossMins	=	its_PL_LossMins;
		Kept_PL_LossMaxs	=	its_PL_LossMaxs;
		Kept_PL_Notios		=	its_PL_Notios;

		Kept_DL_LossMin	=	its_DL_LossMin;
		Kept_DL_LossMax	=	its_DL_LossMax;

		its_CorrelationType	=	CT_FLAT;
				
		int	i;

		if (CHECK_EQUAL(New_KMin, 0.0))
		{
			// EQUITY TRANCHE
			its_DL_LossMin	=	0.0;
			its_DL_LossMax	=	New_KMax;

			fill(its_PL_LossMins.begin(), its_PL_LossMins.end(), 0.0);
			fill(its_PL_LossMaxs.begin(), its_PL_LossMaxs.end(), New_KMax);

			for (i=0; i < its_PL_NbFlows; i++)
				its_PL_Notios[i]	*=	New_KMax / (Kept_PL_LossMaxs[i] - Kept_PL_LossMins[i]);

			its_CorrelationValue	=	Beta_KMax * Beta_KMax;

			// ------------------------------------------------------------
			// READY FOR PRICE
			// ------------------------------------------------------------
			ResetRootPricer();

			Price();
		}
		else
		{
			// MEZZANINE TRANCHE
			its_DL_LossMin	=	0.0;
			fill(its_PL_LossMins.begin(), its_PL_LossMins.end(), 0.0);

			// FIRST PRICE
			its_DL_LossMax	=	New_KMin;
			fill(its_PL_LossMaxs.begin(), its_PL_LossMaxs.end(), New_KMin);
			its_CorrelationValue	=	Beta_KMin * Beta_KMin;

			for (i=0; i < its_PL_NbFlows; i++)
				its_PL_Notios[i]	*=	New_KMin / (Kept_PL_LossMaxs[i] - Kept_PL_LossMins[i]);

			// ------------------------------------------------------------
			// READY FOR PRICE
			// ------------------------------------------------------------
			ResetRootPricer();

			Price();

			// Just keep relevant data
			Kept_NPV			=	NPV;
			Kept_Central_NPV	=	Central_NPV;
			Kept_PremiumLegPV	=	PremiumLegPV;
			Kept_DefaultLegPV	=	DefaultLegPV;
			Kept_ATMPremiumLegWithoutNotio	=	ATMPremiumLegWithoutNotio;
			Kept_NPVSquare		=	NPVSquare;
			Kept_ATMPremiumLeg	=	ATMPremiumLeg;

			double	Notio1	=	its_PL_Notios[0];

			// SECOND PRICE
			its_DL_LossMax	=	New_KMax;
			fill(its_PL_LossMaxs.begin(), its_PL_LossMaxs.end(), New_KMax);
			its_CorrelationValue	=	Beta_KMax * Beta_KMax;

			for (i=0; i < its_PL_NbFlows; i++)
				its_PL_Notios[i]	=	Kept_PL_Notios[i] * New_KMax / (Kept_PL_LossMaxs[i] - Kept_PL_LossMins[i]);

			// ------------------------------------------------------------
			// READY FOR PRICE
			// ------------------------------------------------------------
			ResetRootPricer();

			Price();

			// update NPVs and all
			NPV				-=	Kept_NPV;
			Central_NPV		-=	Kept_Central_NPV;
			PremiumLegPV	-=	Kept_PremiumLegPV;
			DefaultLegPV	-=	Kept_DefaultLegPV;
//			ATMPremiumLegWithoutNotio	-=	Kept_ATMPremiumLegWithoutNotio;
			ATMPremiumLeg	-=	Kept_ATMPremiumLeg;

			double	Notio2	=	its_PL_Notios[0];

			ATMPremiumLegWithoutNotio	*=	Notio2;
			ATMPremiumLegWithoutNotio	-=	Kept_ATMPremiumLegWithoutNotio * Notio1;
			ATMPremiumLegWithoutNotio	/=	(Notio2 - Notio1);
					
			// Standard Error
/*
			NPVSquare	-=	Kept_NPVSquare;
			if (NPVSquare - NPV * NPV >= 0.0)
				StdError	=	sqrt(NPVSquare - NPV * NPV);
			else
				StdError	=	_MINUS_INFINITY_;	// Pb!
*/
			StdError	=	0.0;
			UpFront_PV	=	0.0;

			switch (its_NPV_Type)
			{
			case CNPV_STANDARD:
				break;
				
			case CNPV_WITH_RUNNING_SPREAD:
				UpFront_PV	=	its_UpFront_RunningSpread * ATMPremiumLeg;
				break;
				
			case CNPV_WITH_UP_FRONT_PREMIUM:
				UpFront_PV	=	its_UpFront_Premium * its_PL_Notios[0];

				break;
				
			case CNPV_WITH_RUNNING_SPREAD_AND_UP_FRONT_PREMIUM:
				UpFront_PV	=	its_UpFront_Premium * its_PL_Notios[0];
				UpFront_PV	+=	its_UpFront_RunningSpread * ATMPremiumLeg;
				
				break;
			}


			if (UpFront_PV)
			{
				PremiumLegPV	+=	UpFront_PV;

				NPV	+=	PremNPVFlag * UpFront_PV;
			}

			double	TheAvgSpread;

			TheAvgSpread	=	accumulate(its_PL_Spreads.begin(), its_PL_Spreads.end(), 0.0);
			TheAvgSpread	/=	its_PL_NbFlows;

			if (ATMPremiumLeg)
				ATMSpread	=	DefaultLegPV / ATMPremiumLeg * 10000.0;	// in bps.
			else
				ATMPremiumLeg	=	_MINUS_INFINITY_;	
/*
			if (PremiumLegPV)
				ATMSpread	=	DefaultLegPV / PremiumLegPV * 10000.0 * TheAvgSpread;	// in bps.
			else
				ATMSpread	=	_MINUS_INFINITY_;	
*/

			switch (its_ATMDataFlag)
			{
			case CPLADT_PURESPREAD:
				break;

			case CPLADT_UF_UPFRONT:
				// find the ATM UpFront
				// I have some Running Spread
				// ATMSpread and its_UpFront_RunningSpread (=500) in bps.
				// Output is ATMUpFront in % (30)
				ATMUpFront	=	(ATMSpread - its_UpFront_RunningSpread) / 100.0 * ATMPremiumLegWithoutNotio;

				break;

			case CPLADT_UF_SPREAD:
				// find the ATM Running Spread 
				// I have an UpFrontValue
				ATMRunningSpread	= ATMSpread - its_UpFront_Premium * 100.0 / ATMPremiumLegWithoutNotio;

				break;
			}

		}

		// RESTORE
		its_DL_LossMin	=	Kept_DL_LossMin;
		its_DL_LossMax	=	Kept_DL_LossMax;

		its_PL_LossMins	=	Kept_PL_LossMins;
		its_PL_LossMaxs	=	Kept_PL_LossMaxs;
		
		its_PL_Notios	=	Kept_PL_Notios;

		break;
	}



	if (PriceToDo)
	{
		// ------------------------------------------------------------
		// READY FOR PRICE
		// ------------------------------------------------------------
		ResetRootPricer();

		Price();

	}

	// ------------------------------------------------------------
	// SET OUTPUTS
	// ------------------------------------------------------------
	SetICMOutputsFromCreditManager();
}



void CreditManager::SetICMOutputsFromCreditManager()
{
	//Net Present Value FeeLegPV - DefLegPV (we are short protection)
	SetPrice(NPV);

	// not sure
	//Accrued PV
	// SetAccruedPV(its_Cash);

	// not sure
	//Accrued
	SetAccrued(its_Accrued);

	//Inital price (for sensitivity computing)
	SetInitialPrice(Central_NPV);	// NPV?

	//Premium Leg PV
	SetFeeLegPrice(PremiumLegPV);

	// not sure
	SetFeeLegPrice_Unity(ATMPremiumLegWithoutNotio);

	 //Defaultable Leg PV
	SetDefLegPrice(DefaultLegPV);

	//Breakeven spread 
	SetSpread(ATMSpread);

	//Risky Duration
	SetDuration(ATMPremiumLegWithoutNotio);

	// FLAGS????
}



// ------------------------------------------------------------
// GET DATA FROM PRODUCT TO CREDIT MANAGER
// ------------------------------------------------------------

void CreditManager::GetDataFromCreditProduct()
{
	if (its_Credit_Product == NULL)
		return;

	its_Credit_Product->Get_ValDate(itsValDate);
	its_Credit_Product->Get_ValDateAsDouble(itsValDateAsDouble);
	
	its_Credit_Product->Get_PL_NbFlows(its_PL_NbFlows);

	its_Credit_Product->Get_DL_CreditWindowLow(its_DL_CreditWindowLow);
	its_Credit_Product->Get_DL_CreditWindowUp(its_DL_CreditWindowUp);
	its_Credit_Product->Get_DL_PaymentDate(its_DL_PaymentDate);
	its_Credit_Product->Get_DL_LossMin(its_DL_LossMin);
	its_Credit_Product->Get_DL_LossMax(its_DL_LossMax);
//		int				its_DL_NbDefMin;					// Credit Nb Def Min
//		int				its_DL_NbDefMax;					// Credit Nb Def Max
	its_Credit_Product->Get_DL_PaymentLag(its_DL_PaymentLag);
	its_Credit_Product->Get_DL_PaymentType(its_DL_PaymentType);


	its_Credit_Product->Get_PL_CreditWindowLows(its_PL_CreditWindowLows);
	its_Credit_Product->Get_PL_CreditWindowUps(its_PL_CreditWindowUps);
	its_Credit_Product->Get_PL_StartDates(its_PL_StartDates);
	its_Credit_Product->Get_PL_EndDates(its_PL_EndDates);
	its_Credit_Product->Get_PL_PaymentDates(its_PL_PaymentDates);
	its_Credit_Product->Get_PL_LossMins(its_PL_LossMins);
	its_Credit_Product->Get_PL_LossMaxs(its_PL_LossMaxs);
	its_Credit_Product->Get_PL_Ratios(its_PL_Ratios);
	its_Credit_Product->Get_PL_Notios(its_PL_Notios);
	its_Credit_Product->Get_PL_Spreads(its_PL_Spreads);
	its_Credit_Product->Get_PL_CreditFlags(its_PL_CreditFlags);
	its_Credit_Product->Get_PL_CreditSpreadCaps(its_PL_CreditSpreadCaps);
	its_Credit_Product->Get_PL_Redemptions(its_PL_Redemptions);
		
//		DoubleVector	its_PL_NbDefMins;
//		DoubleVector	its_PL_NbDefMaxs;
		
	its_Credit_Product->Get_PL_PaymentType(its_PL_PaymentType);


	its_Credit_Product->Get_UpFront_RunningSpread(its_UpFront_RunningSpread);
	its_Credit_Product->Get_UpFront_Premium(its_UpFront_Premium);

	its_Credit_Product->Get_PricingLegsType(its_PricingLegsType);

	its_Credit_Product->Get_DefNPVFlag(DefNPVFlag);
	its_Credit_Product->Get_PremNPVFlag(PremNPVFlag);

	its_Credit_Product->Get_CreditObservationType(its_CreditObservationType);
	its_Credit_Product->Get_CreditAccruedPayment(its_CreditPremiumLegAccrued);
	its_Credit_Product->Get_ATMDataFlag(its_ATMDataFlag);
	its_Credit_Product->Get_NPV_Type(its_NPV_Type);
	its_Credit_Product->Get_CDO_Type(its_CDO_Type);
	its_Credit_Product->Get_TimeStep_ProrataF(its_Time_Step_Prorata_1F);

}


// ------------------------------------------------------------
// GET DATA FROM PRODUCT TO CREDIT MANAGER
// ------------------------------------------------------------

void CreditManager::GetDataFromCreditMultiCurves()
{
	if (its_Credit_Market_Data == NULL)
		return;

	int nbCredit = its_Credit_Market_Data->Get_NbCredits();
	int NbMaturities  = its_Credit_Market_Data->Get_NbMaturities();
	
	its_Credit_Market_Data->Get_RollDateFlag(its_RollDateFlag);
	its_Credit_Market_Data->Get_CDO_Square_Nb_Underlyings(its_CDO_Square_Nb_Underlyings);

	its_Credit_Market_Data->Get_BumpSpread(its_BumpSpread);
	its_Credit_Market_Data->Get_BumpType(its_BumpSpread_Type);
	its_Credit_Market_Data->Get_BumpRecovery(its_BumpRecovery);
	its_Credit_Market_Data->Get_BumpCorrelation(its_BumpCorrelation);


	its_Credit_Market_Data->Get_Categories(its_Categories);
	its_Credit_Market_Data->Get_Currencies(its_Currencies);
	its_Credit_Market_Data->Get_Accrueds(its_Accrueds);
	its_Credit_Market_Data->Get_Recoveries(its_Recoveries);
	its_Credit_Market_Data->Get_Notionals(its_Notionals);
	its_Credit_Market_Data->Get_Losses(its_Input_Losses);
	its_Credit_Market_Data->Get_DefaultDates(its_DefaultDates);
	its_Credit_Market_Data->Get_AmortizationDates(its_AmortizationDates);

	Set_CreditsLabels(its_Credit_Market_Data->Get_CreditsLabels());

	Set_CreditDataMaturities(its_Credit_Market_Data->Get_CreditDataMaturities());

	// Spreads
	its_CreditDataSpreads	=	its_Credit_Market_Data->Get_CreditDataSpreads();

	GenerateAllDefaultCurves();
}


// ------------------------------------------------------------
// COMPUTE SENSITIVITY METHOD FOR BASKETS
// ------------------------------------------------------------

double	CreditManager::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
										  const std::string&  plot, 
										  const std::string&label, 
											double  epsvalue,  double epsilonGamma //useless 
											)
{
	int	credit_Id;

	if (AreHedgesToCompute())
	{
		its_HedgesRunning	=	true;

		switch (typesensi)
		{
			case ICM_CM_SPREAD_PARALLEL:

				its_Bump_Choice	=	CHB_SPREAD_PARALLEL;

				if (epsvalue != CREDIT_DEFAULT_VALUE)
					its_BumpSpread	=	epsvalue;

				break;

			case ICM_CM_DEFAULT:

				its_Bump_Choice	=	CHB_DEFAULT;

				break;

			case ICM_CM_RECOVERY_SENS:

				its_Bump_Choice	=	CHB_RECOVERY_SENS;

				if (epsvalue != CREDIT_DEFAULT_VALUE)
					its_BumpRecovery	=	epsvalue;

				break;

			case ICM_CM_RECOVERY_LOSS:

				its_Bump_Choice	=	CHB_RECOVERY_LOSS;

				if (epsvalue != CREDIT_DEFAULT_VALUE)
					its_BumpRecovery	=	epsvalue;

				break;

			case ICM_CM_RECOVERY_GLOBAL:

				its_Bump_Choice	=	CHB_RECOVERY_GLOBAL;

				if (epsvalue != CREDIT_DEFAULT_VALUE)
					its_BumpRecovery	=	epsvalue;

				break;

			case ICM_CM_CORRELATION:

				its_Bump_Choice	=	CHB_CORRELATION;

				if (epsvalue != CREDIT_DEFAULT_VALUE)
					its_BumpCorrelation	=	epsvalue;

				break;

			default:

				its_Bump_Choice	=	CHB_NO;
		}

		// -----------------------
		PriceOrHedge();
		// -----------------------
	
		SetHedges_Computed();
		its_HedgesRunning	=	false;
	}
	
	// -----------------------
	// OUTPUTS
	// -----------------------
	switch (typesensi)
	{
		case ICM_CM_SPREAD_PARALLEL:

			credit_Id	=	Get_LabelId(label);

			break;

		case ICM_CM_DEFAULT:

			credit_Id	=	Get_LabelId(label);

			break;

		case ICM_CM_RECOVERY_SENS:

			credit_Id	=	Get_LabelId(label);

			break;

		case ICM_CM_RECOVERY_LOSS:

			credit_Id	=	Get_LabelId(label);

			break;

		case ICM_CM_RECOVERY_GLOBAL:

			credit_Id	=	0;

			break;

		case ICM_CM_CORRELATION:

			credit_Id	=	0;

			break;

		default:

			credit_Id	=	0;

			break;
	}
	
	if ((credit_Id < 0) || (credit_Id >= its_CreditsLabels.size()))
		ICMTHROW(ERR_INVALID_DATA,"Bad Id For Hedges!");

	
	return	(its_Hedges_Delta_NPV[credit_Id]);

}


int	CreditManager::Get_LabelId(const std::string&  data)
{
	int	i;

	for (i=0;i<its_CreditsLabels.size();i++)
	{
		if (!strcmp(its_CreditsLabels[i].c_str(), data.c_str()))
			return	i;
	}

	return -1;
}

void	
CreditManager::Set_CorrelationMatrix(ICM_CorrMatrix* correlmatrix)
{
	if (its_CorrelationMatrix)
		delete its_CorrelationMatrix;
	its_CorrelationMatrix = (ICM_CorrMatrix*) correlmatrix->Clone();
}



CreditManager::~CreditManager(void)
	{
		if (its_ZeroCurve)
			delete	its_ZeroCurve;
		its_ZeroCurve	=	NULL;

		if (itsZCValuesForCalibration)
			delete itsZCValuesForCalibration;
		itsZCValuesForCalibration = NULL;

		its_IR_Lags.clear();
		its_IR_ZC_Yields.clear();
		
		if (its_CreditDataSpreads)
			delete its_CreditDataSpreads;
		its_CreditDataSpreads	=	NULL;

		its_Categories.clear();
		its_Currencies.clear();
		its_Accrueds.clear();
		its_Recoveries.clear();
		its_Notionals.clear();
		its_Input_Losses.clear();
		its_DefaultDates.clear();
		its_AmortizationDates.clear();

		its_CreditsLabels.clear();

		if (its_ModelMultiCurves)
			delete its_ModelMultiCurves;
		its_ModelMultiCurves	=	NULL;

		its_Maturities.clear();

		its_Beta.clear();
		its_Base_Correlation_Strikes.clear();
		its_Base_Correlation_Values.clear();

		its_FL_Alpha.clear();
		its_FL_Beta1.clear();
		its_FL_Beta2.clear();
		its_FL_Beta1_Complement.clear();
		its_FL_Beta2_Complement.clear();

		its_Used_Beta.clear();
		its_Used_SQRT_OneMinusBetaSquare.clear();

		// Premium Leg Data
		its_PL_CreditWindowLows.clear();
		its_PL_CreditWindowUps.clear();
		its_PL_StartDates.clear();
		its_PL_EndDates.clear();
		its_PL_PaymentDates.clear();

		its_PL_CreditWindowUps_DF.clear();
		its_PL_PaymentDates_DF.clear();

		its_PL_Sorted_PaymentDates.clear();
		its_PL_Sorted_CreditWindowUps.clear();
		its_PL_Sorted_CreditWindowUps_DF.clear();
		its_PL_Sorted_PaymentDates_DF.clear();

		its_PL_LossMins.clear();
		its_PL_LossMaxs.clear();
		its_PL_NbDefMins.clear();
		its_PL_NbDefMaxs.clear();
		its_PL_Ratios.clear();
		its_PL_Notios.clear();
		its_PL_Spreads.clear();

		its_PL_CreditFlags.clear();
		its_PL_CreditSpreadCaps.clear();
		its_PL_Redemptions.clear();

		if (itsCorrelation)
			delete itsCorrelation;
		itsCorrelation	=	NULL;

		if (its_CorrelationMatrix)
			delete its_CorrelationMatrix;
		its_CorrelationMatrix	=	NULL;

		if (its_CholeskyMatrix)
			delete its_CholeskyMatrix;
		its_CholeskyMatrix	=	NULL;

		if (its_ArrayDefaultCrv)
			for(int i=0; i<its_CreditsLabels.size(); i++){
				if(its_ArrayDefaultCrv[i]) 
					delete its_ArrayDefaultCrv[i];
			}
			delete [] its_ArrayDefaultCrv;
		its_ArrayDefaultCrv	=	NULL;

		if (its_MatrixShiftedDefaultCrv)
			delete its_MatrixShiftedDefaultCrv;
		its_MatrixShiftedDefaultCrv	=	NULL;

		its_DF.clear();

		itsSortedDefaultTimes.clear();
		itsBarriers_Standard.clear();
		itsBarriers_Shifted.clear();

		itsDef_Prob.clear();
		
		its_DefaultTimes.clear();
		its_DefaultTimesShifted.clear();
		itsSortedDefaultTimes_Keep.clear();

		if (itsGenerator)
			delete itsGenerator;
		itsGenerator = NULL;

//			memset(its_CreditCalibratorMaturities_AsChar,'\0', sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS); 

//			memset(its_Maturities_AsChar,'\0', sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS); 
//			FreePointerTabChar(its_Maturities_AsChar,	its_NbMaturities);	

		//FreePointerTabChar(its_CreditsLabelsAsChar,	its_NbCredits);	

		its_SortedCreditObservationDates.clear();
		its_CreditObervationDatesSortedLowsIds.clear();
		its_CreditObervationDatesSortedUpsIds.clear();

		its_SortedCreditObservationDatesAndLosses.clear();

		its_CreditCalibrator_Lags.clear();
		its_CreditCalibrator_Prob_Default.clear();
		its_CreditCalibrator_Std_Error.clear();

		its_AllNPVS.clear();

		if (its_AllMatrixNPVS)
			delete its_AllMatrixNPVS;
		its_AllMatrixNPVS = NULL;

		if (its_Shift_Results)
			delete its_Shift_Results;
		its_Shift_Results = NULL;

		its_AllDefLegs.clear();
		its_AllPremLegs.clear();

		its_Hedges_Delta_NPV.clear();
		its_Hedges_Hedge_Ratio.clear();
		its_Hedges_Delta_CDS.clear();
		its_Hedges_CDS_ATM_Margin.clear();
		its_Hedges_Carry.clear();

		its_Hedges_MatrixOutput.clear();

		its_CC_Lags_String	=	"";
		its_CC_Prob_String	=	"";
		its_CC_Std_Error_String	=	"";
		its_HedgesCDSMaturity	=	"";

		IS_Common_exp_twistss.clear();
		IS_individual_exp_twists.clear();
		IS_individual_exp_twists_Kept.clear();
		TheCond_DefProb_Vect.clear();

		its_LossRate.clear();

		its_lossdistrib.clear();
		its_Current_DefProb.clear();

		if (its_ProbCond)
			delete its_ProbCond;
		its_ProbCond = NULL;

		if (its_ProbCond_Perturb)
			delete its_ProbCond_Perturb;
		its_ProbCond_Perturb = NULL;

		if (its_taildistrib_perturb)
			delete its_taildistrib_perturb;
		its_taildistrib_perturb = NULL;

		if (its_barrier_perturb)
			delete its_barrier_perturb;
		its_barrier_perturb = NULL;

		if (its_CDO_Square_Credit_Ids)
			delete its_CDO_Square_Credit_Ids;
		its_CDO_Square_Credit_Ids	=	NULL;
		if (its_CDO_Square_Credit_Losses)
			delete its_CDO_Square_Credit_Losses;
		its_CDO_Square_Credit_Losses	=	NULL;
		if (its_CDO_Square_Credit_Notionals)
			delete its_CDO_Square_Credit_Notionals;
		its_CDO_Square_Credit_Notionals	=	NULL;

		its_CDO_Square_CDO_Includes.clear();
		its_CDO_Square_CDO_Maturities.clear();
		its_CDO_Square_CDO_NbCredits.clear();
		its_CDO_Square_CDO_Attachments.clear();
		its_CDO_Square_CDO_Notionals.clear();

		its_CDO_Square_Relevant_CDO_Ids.clear();
		its_CDO_Square_Relevant_CDO_MaturitiesInYearFraction.clear();
		its_CDO_Square_Relevant_CDO_NbCredits.clear();
		its_CDO_Square_Relevant_CDO_Attachments_Pct.clear();
		its_CDO_Square_Relevant_CDO_Detachments_Pct.clear();
		its_CDO_Square_Relevant_CDO_Attachments_Notio.clear();
		its_CDO_Square_Relevant_CDO_Detachments_Notio.clear();

		its_CDO_Square_Relevant_CDO_Notionals.clear();
		its_CDO_Square_Relevant_CDO_Rescaling_Factor.clear();

		if (its_CDO_Square_Relevant_Credit_Ids)
			delete its_CDO_Square_Relevant_Credit_Ids;
		its_CDO_Square_Relevant_Credit_Ids	=	NULL;

		if (its_CDO_Square_Relevant_Credit_Includes)
			delete its_CDO_Square_Relevant_Credit_Includes;
		its_CDO_Square_Relevant_Credit_Includes	=	NULL;

		if (its_CDO_Square_Relevant_Credit_Losses)
			delete its_CDO_Square_Relevant_Credit_Losses;
		its_CDO_Square_Relevant_Credit_Losses	=	NULL;

		if (its_DefCurve)
			delete its_DefCurve;
		its_DefCurve	=	NULL;

		its_CurrentHermiteCoeffs.clear();

		its_CommonFactors.clear();
		its_BarrierDerivatives_FastHedge.clear();

		its_CopulaChoice	=	NULL;

//			if (its_fOut)
//				delete its_fOut;
		its_fOut	=	NULL;

		its_IsActivateShift	=	false;
		LossTShifts.clear();

		its_Used_Beta_SC_Vector.clear();
		its_Used_SQRT_OneMinusBetaSquare_SC_Vector.clear();
		its_SC_Coefficients.clear();

//			if (its_Correlation_Calibrator)
//				delete its_Correlation_Calibrator;
//			its_Correlation_Calibrator	=	NULL;

		// ------------------------------------------------------------
//			if (its_Credit_Product)
//				delete its_Credit_Product;
//			its_Credit_Product	=	NULL;

//			if (its_Credit_Market_Data)
//				delete its_Credit_Market_Data;
//			its_Credit_Market_Data	=	NULL;
		// ------------------------------------------------------------

		its_Loss_Distrib_Map.clear();
	}