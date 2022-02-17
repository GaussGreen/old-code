#pragma warning(disable :4786)

#include "firstToBeIncluded.h"

#include "ICM_local_Str.h"

#include <ARM\libarm_local\ARM_local_mod.h>
#include <ARM\libarm_local\ARM_local_persistent.h>
#include <ARM\libarm_local\ARM_local_glob.h>

/// remove warnings on va_start and va_end
/// headers to remove definition of va_start and va_end as this is redefined later on!
/// handle this with care after sorting out why this is so!!
/// and do not change the order as this one has to be precisely here
/// if you do not know what you are doing, please ask someone who does!
#include <ARM\libarm_local\undef_va_vars.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>


#include "ICMKernel\crv\icm_defaultcurve.h"
#include "ICMKernel\inst\icm_gen.h"

#include "ICMKernel\cair\credit_manager.h"
#include "ICMKernel\mod\modelmulticurves.h"
#include "ICMKernel/glob/icm_corrmatrix.h"

long ICMLOCAL_SetCreditManager_MarketData(
							long			CreditManagerId,
							const double&	AsOfDate,
							CCString		CurrencyLabel,
							int				IR_SummitCurveId,
							long			MarketData_ParametersId,
							long			Imposed_IR_ValuesId,
							ARM_result&		result)
{
	double dResult=0.;

	CreditManager*	ACreditManager	=	NULL;

	ARM_ZeroCurve* zc = NULL;
	ICM_GenCF* TheCashFlows = NULL;

	ICM_Matrix<ARM_Vector>* Matrix_Parameters = NULL;
	ICM_Matrix<ARM_Vector>* Matrix_Imposed_IR = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	ARM_Date NewDate;
	double ModDate = 0.;

	// Get AsOfDate
	char* sDate = new char[11];
	Local_XLDATE2ARMDATE(AsOfDate,sDate);

	CCString msg ("");

	try
	{
		// ---------------------------------------------------------------------------------
		ACreditManager = (CreditManager*) LOCAL_PERSISTENT_OBJECTS->GetObject(CreditManagerId);
		// ICM_CREDIT_MANAGER
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ACreditManager, ICM_PRICER) == 0)
		{
			result.setMsg ("ARM_ERR: CreditManager is not of a good type");
			
			return ARM_KO;
		}

		// ---------------------------------------------------------------------------------
		// Summit ZC Curve
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(IR_SummitCurveId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Default Zc Curve is not of a good type");
			return ARM_KO;
		}
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// Market Data, series of Cash Flows
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(MarketData_ParametersId); 

		if (TheCashFlows)	Matrix_Parameters = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// Market Data, series of Cash Flows
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(Imposed_IR_ValuesId); 

		if (TheCashFlows)	Matrix_Imposed_IR = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		ACreditManager->SetValDate((ARM_Date)sDate);
		ACreditManager->SetZeroCurve(zc);
		ACreditManager->SetMarketDataParameters(Matrix_Parameters);
		ACreditManager->SetImposedIRCurves(Matrix_Imposed_IR);

		// ---------------------------------------------------------------------------------
		result.setDouble(dResult);

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		return ARM_OK;
	}

	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	if (sDate)
		delete [] sDate;
	sDate = NULL;

	return ARM_KO;
}


// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


long ICMLOCAL_SetCreditManager_CreditData(
							long			CreditManagerId,
							VECTOR<CCString>	labels,
							long				CD_DescriptionId,
							vector<double>		MatIssuersSpread,
							VECTOR<CCString>	maturities,
							long				CD_ParametersId,
//							const VECTOR<long>& DefCurvesID,
							CCString			HedgesCDSMaturity,
							long				CD_CDOSquare_ParametersId,
							long				CD_CDOSquare_DataId,
							ARM_result&			result)
{
	double dResult=0.;

	CreditManager*	ACreditManager	=	NULL;

	ICM_DefaultCurve* pCdsZC = NULL;
	ICM_GenCF* TheCashFlows = NULL;

	ICM_Matrix<ARM_Vector>* Matrix_Description = NULL;
	ICM_Matrix<ARM_Vector>* Matrix_Parameters = NULL;
	ICM_Matrix<ARM_Vector>* Matrix_CDOSquare_Parameters = NULL;
	ICM_Matrix<ARM_Vector>* Matrix_CDOSquare_Data = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	int i,j,il;
	int nbissuers = 0;
	int nbmat= 0;
	char** pmaturities = NULL;

//	int nbdefcurves = DefCurvesID.size();

	vector<string> AllCreditLabels;
	vector<string> AllMaturities;

	try
	{
		// ---------------------------------------------------------------------------------
		ACreditManager = (CreditManager*) LOCAL_PERSISTENT_OBJECTS->GetObject(CreditManagerId);
		// ICM_CREDIT_MANAGER
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ACreditManager, ICM_PRICER) == 0)
		{
			result.setMsg ("ARM_ERR: CreditManager is not of a good type");
			
			return ARM_KO;
		}

		// ---------------------------------------------------------------------------------
		// Labels
		nbissuers = labels.size();
		AllCreditLabels.resize(nbissuers);

		if (nbissuers)
		{
			for (j=0; j<nbissuers; j++)
				AllCreditLabels[j] = (string)(labels[j]);
		}

		/*char** psLabels = new char*[nbissuers];

		for (j=0;j<nbissuers;j++)
		{
			psLabels[j] = new char[60];
			sprintf(psLabels[j],(const char*)labels[j]);
		}*/

		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// Credit Data Description
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(CD_DescriptionId); 

		if (TheCashFlows)	Matrix_Description = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// Credit Data Parameters
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(CD_ParametersId); 

		if (TheCashFlows)	Matrix_Parameters = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// CDO^2 Credit Data Parameters
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(CD_CDOSquare_ParametersId); 

		if (TheCashFlows)	Matrix_CDOSquare_Parameters = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// CDO^2 Credit Data Data
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(CD_CDOSquare_DataId); 

		if (TheCashFlows)	Matrix_CDOSquare_Data = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// Maturities
		nbmat = maturities.size();
		AllMaturities.resize(nbmat);

		if (nbmat)
		{
			for (il=0; il<nbmat; il++)
			{
				AllMaturities[il]	=	(string)(maturities[il]);
			}
		}

		/*char psMatu [ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS];

		for (j=0;j<ARM_NB_TERMS;j++)	//	400!	x ARM_NB_MAX_CHAR_TERMS (30)
			sprintf(psMatu[j],"X");
		
		for (j=0;j<nbmat;j++)
			sprintf(psMatu[j], (const char*)maturities[j]); 
		*/
		// in order to skip some pillars (the same for every credit!)
/*		i2=0;
		for (j=0;j<nbmat;j++)
		{
			if (Rates[j] != -999.)
			{
				sprintf(psMatu[i2], (const char*)maturities[j]);
				i2++;
			}
		}
*/
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// Matrice de Spread	
		// should check the size
		ICM_QMatrix<double>* CreditSpreads = new ICM_QMatrix<double>(nbissuers,nbmat,CREDIT_DEFAULT_VALUE);
			
		for (i=0;i<nbissuers;i++)
			for (j=0;j<nbmat;j++)
				CreditSpreads->SetValue(i,j,MatIssuersSpread[i*nbmat+j] / 10000.0);		// in bps.

		// in order to skip -999.
/*
		pdRate = new double[real_size2];
		
		i2=0;
		for(i = 0; i < real_size; i++)
		{ if (Rates[i] != -999.) {pdRate[i2] = Rates[i] / 10000.;i2++;} }

		VRates  = new ARM_Vector(real_size2,pdRate);
*/			
		// ---------------------------------------------------------------------------------
/*
		if (nbdefcurves)
		{
			for (j=0; j < nbdefcurves; j++) {if (DefCurvesID[j] == -1) sizereal-=1;}

			ICM_DefaultCurve** pZCCDS = new ICM_DefaultCurve*[sizereal];
			
			modif = 0;
			for (j = 0; j < nbdefcurves; j++)
			{
				if (DefCurvesID[j] == -1) continue;

				// on recupere les courbes ****************************************************************
				pCdsZC = (ICM_DefaultCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DefCurvesID[j]);

				if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pCdsZC, ICM_DEFAULTCURVE) == 0)
				{
					result.setMsg ("ARM_ERR: one CDS Zc Curve is not of a good type");
				
					return ARM_KO;
				}

				pZCCDS[modif] = pCdsZC;
				modif++;
			}
		}
*/
		// ---------------------------------------------------------------------------------
		//ACreditManager->Set_NbCredits(nbissuers);
		ACreditManager->Set_CreditsLabels(AllCreditLabels);
//		ACreditManager->Set_CreditsLabelsAsChar(psLabels);
		ACreditManager->SetCreditDataDescription(Matrix_Description);
		ACreditManager->SetCreditDataParameters(Matrix_Parameters);
		ACreditManager->Set_CreditDataSpreads(CreditSpreads);
		ACreditManager->Set_CreditDataMaturities(AllMaturities);
//		ACreditManager->Set_CreditDataMaturitiesAsChar(psMatu);
		ACreditManager->Set_CreditDataHedgesCDSMaturity(string (HedgesCDSMaturity));
		// ---------------------------------------------------------------------------------
		
		ACreditManager->SetCreditDataCDOSquare_Parameters(Matrix_CDOSquare_Parameters);
		ACreditManager->SetCreditDataCDOSquare_Data(Matrix_CDOSquare_Data);

		// ---------------------------------------------------------------------------------
		// All Default Curves to be created
		ACreditManager->GenerateAllDefaultCurves();

//		ACreditManager->GenerateAllDefaultCurves();

		// ---------------------------------------------------------------------------------
		result.setDouble(dResult);

		return ARM_OK;
	}

	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	return ARM_KO;
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


long ICMLOCAL_SetCreditManager_CreditModel(
							long			CreditManagerId,
							long				CM_ParameterId,
							double				correlation_value,
							VECTOR<double>		beta_vector,
							VECTOR<double>		base_correlation_strikes,
							VECTOR<double>		base_correlation_values,
							long				CorrelationMatrixId,
							long				CorrelationId,
							long				FactorLoading_ParameterId,
							ARM_result&			result)
{
	double dResult=0.;

	CreditManager*	ACreditManager	=	NULL;

	ICM_DefaultCurve* pCdsZC = NULL;
	ICM_GenCF* TheCashFlows = NULL;

	ICM_Matrix<ARM_Vector>* Matrix_Parameters = NULL;
	ICM_Matrix<ARM_Vector>* FactorLoadings_Parameters = NULL;
	ICM_Correlation* pCorrelation = NULL;
	ICM_CorrMatrix* pCorrelationMatrix = NULL;

	double** CorrelationMtx = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		// ---------------------------------------------------------------------------------
		ACreditManager = (CreditManager*) LOCAL_PERSISTENT_OBJECTS->GetObject(CreditManagerId);
		// ICM_CREDIT_MANAGER
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ACreditManager, ICM_PRICER) == 0)
		{
			result.setMsg ("ARM_ERR: CreditManager is not of a good type");
			
			return ARM_KO;
		}

		// ---------------------------------------------------------------------------------
		// Credit Model Description
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(CM_ParameterId); 

		if (TheCashFlows)	Matrix_Parameters = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// Correlation Value
		// ---------------------------------------------------------------------------------
		
		// ---------------------------------------------------------------------------------
		// BETA vector
		// ---------------------------------------------------------------------------------		

		// ---------------------------------------------------------------------------------
		// Base Correlation Strikes
		// ---------------------------------------------------------------------------------		

		// ---------------------------------------------------------------------------------
		// Base Correlation Values
		// ---------------------------------------------------------------------------------		

		// ---------------------------------------------------------------------------------
		// Correlation Matrix Object
		// ---------------------------------------------------------------------------------		
		pCorrelationMatrix = (ICM_CorrMatrix *) LOCAL_PERSISTENT_OBJECTS->GetObject(CorrelationMatrixId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pCorrelationMatrix, ICM_CORRELATION) == 0) 
		{
			result.setMsg ("ARM_ERR: Correlation Matrix Object is not of a good type");
			return ARM_KO;
		}

		// ---------------------------------------------------------------------------------
		// Correlation Object
		// ---------------------------------------------------------------------------------		
		pCorrelation = (ICM_Correlation *) LOCAL_PERSISTENT_OBJECTS->GetObject(CorrelationId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pCorrelation, ICM_CORRELATION) == 0) 
		{
			result.setMsg ("ARM_ERR: Correlation Object is not of a good type");
//			return ARM_KO;
		}
		else
			ACreditManager->SetCorrelation(pCorrelation);

		// ---------------------------------------------------------------------------------
		// Factor Loadings Data
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(FactorLoading_ParameterId); 

		if (TheCashFlows)	FactorLoadings_Parameters = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		ACreditManager->SetCreditModelParameters(Matrix_Parameters);
		ACreditManager->Set_CorrelationValue(correlation_value);
		ACreditManager->Set_Beta(beta_vector);
		ACreditManager->Set_Base_Correlation_Strikes(base_correlation_strikes);
		ACreditManager->Set_Base_Correlation_Values(base_correlation_values);
		ACreditManager->Set_CorrelationMatrix(pCorrelationMatrix);
//		ACreditManager->SetCorrelation(pCorrelation);
		ACreditManager->SetFactorLoadingsParameters(FactorLoadings_Parameters);
		// ---------------------------------------------------------------------------------
		result.setDouble(dResult);

		return ARM_OK;
	}

	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	return ARM_KO;
}



extern long ICMLOCAL_SetCreditManager_CreditProduct(
							long				CreditManagerId,
							long				CP_DefaultLegId,
							long				CP_PremiumLegId,
							long				CP_PricingParametersId,
							ARM_result&			result)
{
	double dResult=0.;
	CCString msg ("");

	CreditManager*	ACreditManager	=	NULL;

	ICM_GenCF* TheCashFlows = NULL;

	ICM_Matrix<ARM_Vector>* Matrix_Default = NULL;
	ICM_Matrix<ARM_Vector>* Matrix_Premium = NULL;
	ICM_Matrix<ARM_Vector>* Matrix_PricingParameters = NULL;

	try
	{
		// ---------------------------------------------------------------------------------
		ACreditManager = (CreditManager*) LOCAL_PERSISTENT_OBJECTS->GetObject(CreditManagerId);
		// ICM_CREDIT_MANAGER
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ACreditManager, ICM_PRICER) == 0)
		{
			result.setMsg ("ARM_ERR: CreditManager is not of a good type");
			
			return ARM_KO;
		}

		// ---------------------------------------------------------------------------------
		// Market Data, series of Cash Flows
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(CP_DefaultLegId); 

		if (TheCashFlows)	Matrix_Default = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// Market Data, series of Cash Flows
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(CP_PremiumLegId); 

		if (TheCashFlows)	Matrix_Premium = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// Pricing Parameters, series of Cash Flows
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(CP_PricingParametersId); 

		if (TheCashFlows)	Matrix_PricingParameters = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		ACreditManager->SetCreditProductDefaultMatrix(Matrix_Default);
		ACreditManager->SetCreditProductPremiumMatrix(Matrix_Premium);
		ACreditManager->SetCreditProductPricingParametersMatrix(Matrix_PricingParameters);

		// temporary
		ACreditManager->PriceOrHedge();

		result.setDouble(dResult);

		return ARM_OK;
	}

	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	return ARM_KO;
}


// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


extern long ICMLOCAL_GetCreditManager_DataFromLabel(
							long				CreditManagerId,
							CCString			DataLabel,
							ARM_result&			result)
{
	double dResult=0.;

	CreditManager*	ACreditManager	=	NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		// ---------------------------------------------------------------------------------
		ACreditManager = (CreditManager*) LOCAL_PERSISTENT_OBJECTS->GetObject(CreditManagerId);
		// ICM_CREDIT_MANAGER
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ACreditManager, ICM_PRICER) == 0)
		{
			result.setMsg ("ARM_ERR: CreditManager is not of a good type");
			
			return ARM_KO;
		}

		// ---------------------------------------------------------------------------------
		ACreditManager->GetDataFromLabel(string (DataLabel), dResult);
		// ---------------------------------------------------------------------------------
		
		result.setDouble(dResult);

		return ARM_OK;
	}

	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	return ARM_KO;

}


// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


extern long ICMLOCAL_GetCreditManager_DataMatrixFromLabel(
							long				CreditManagerId,
							CCString			DataLabel,
							VECTOR<double*>&	OutputMatrix,
							VECTOR<CCString>&	OutputLabels,
							int&				OutputNbRows,
							int&				OutputNbCols,
							ARM_result&			result)

{
	CreditManager*	ACreditManager	=	NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	int	j;
	vector<string> AllCreditLabels;

	try
	{
		// ---------------------------------------------------------------------------------
		ACreditManager = (CreditManager*) LOCAL_PERSISTENT_OBJECTS->GetObject(CreditManagerId);
		// ICM_CREDIT_MANAGER
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ACreditManager, ICM_PRICER) == 0)
		{
			result.setMsg ("ARM_ERR: CreditManager is not of a good type");
			
			return ARM_KO;
		}


		// ---------------------------------------------------------------------------------
		ACreditManager->GetDataMatrixFromLabel(string (DataLabel), OutputMatrix, AllCreditLabels, OutputNbRows, OutputNbCols);
		// ---------------------------------------------------------------------------------
		
		// ---------------------------------------------------------------------------------
		// Labels
		OutputLabels.clear();
		for (j=0; j<OutputNbCols; j++)
			OutputLabels.push_back((AllCreditLabels[j]).c_str());

		// ---------------------------------------------------------------------------------
		return ARM_OK;
	}

	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	return ARM_KO;

}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


extern long ICMLOCAL_SetCreditManager_CreditCalibrator(
							long			CreditManagerId,
							long				CC_DescriptionId,
							VECTOR<CCString>	maturities,
							ARM_result&			result)
{
	double dResult=0.;

	CreditManager*	ACreditManager	=	NULL;
	ICM_Matrix<ARM_Vector>* Matrix_Description = NULL;

	ICM_GenCF* TheCashFlows = NULL;
	char** pmaturities = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");
	
	int j,il;
	int nbmat= 0;
	vector<string> AllMaturities;

	try
	{
		// ---------------------------------------------------------------------------------
		ACreditManager = (CreditManager*) LOCAL_PERSISTENT_OBJECTS->GetObject(CreditManagerId);
		// ICM_CREDIT_MANAGER
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ACreditManager, ICM_PRICER) == 0)
		{
			result.setMsg ("ARM_ERR: CreditManager is not of a good type");
			
			return ARM_KO;
		}

		// ---------------------------------------------------------------------------------
		// Credit Data Description
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(CC_DescriptionId); 

		if (TheCashFlows)	Matrix_Description = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// Maturities
		nbmat = maturities.size();
		AllMaturities.resize(nbmat);

		if (nbmat)
		{
			for (il=0; il<nbmat; il++)
			{
				AllMaturities[il]	=	(string)(maturities[il]);
			}
		}

		char psMatu [ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS];

		for (j=0;j<ARM_NB_TERMS;j++)	//	400!	x ARM_NB_MAX_CHAR_TERMS (30)
			sprintf(psMatu[j],"X");
		
		for (j=0;j<nbmat;j++)
			sprintf(psMatu[j], (const char*)maturities[j]);

		// ---------------------------------------------------------------------------------
		// ---------------------------------------------------------------------------------
		ACreditManager->SetCreditCalibrator_Description(Matrix_Description);
		ACreditManager->Set_CreditCalibratorMaturitiesAsChar(nbmat,psMatu);
		ACreditManager->Credit_Calibrate(AllMaturities);
		// ---------------------------------------------------------------------------------
		
		result.setDouble(dResult);

		return ARM_OK;
	}

	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	return ARM_KO;

}


// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


extern long ICMLOCAL_SetCreditManager_CorrelationCalibrator(
							long				CreditManagerId,
							long				CC_DescriptionId,
							long				CC_ParametersId,
							VECTOR<long>		TranchesId,
							ARM_result&			result)
{
	int	j;

	double dResult=0.;

	CreditManager*	ACreditManager	=	NULL;
	ICM_Matrix<ARM_Vector>* Matrix_Description = NULL;
	ICM_Matrix<ARM_Vector>* Matrix_Parameters = NULL;

	ICM_GenCF* TheCashFlows = NULL;
	char** pmaturities = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");
	
	try
	{
		// ---------------------------------------------------------------------------------
		ACreditManager = (CreditManager*) LOCAL_PERSISTENT_OBJECTS->GetObject(CreditManagerId);
		// ICM_CREDIT_MANAGER
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ACreditManager, ICM_PRICER) == 0)
		{
			result.setMsg ("ARM_ERR: CreditManager is not of a good type");
			
			return ARM_KO;
		}

		// ---------------------------------------------------------------------------------
		// Calibrator Data Description (attachment, detachment, compound correlation, etc.)
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(CC_DescriptionId); 

		if (TheCashFlows)	Matrix_Description = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// Calibrator Data Description (attachment, detachment, compound correlation, etc.)
		// ---------------------------------------------------------------------------------
		TheCashFlows = (ICM_GenCF*) LOCAL_PERSISTENT_OBJECTS->GetObject(CC_ParametersId); 

		if (TheCashFlows)	Matrix_Parameters = TheCashFlows->GetMatrix();
		// ---------------------------------------------------------------------------------

		// ---------------------------------------------------------------------------------
		// The Tranches
		// ---------------------------------------------------------------------------------

		int	size	=	TranchesId.size();
		vector<ICM_Mez*>	arraypMez;
		ICM_Mez*	pMez;
		
		for (j=0; j<size; j++)
		{
			pMez = (ICM_Mez *) LOCAL_PERSISTENT_OBJECTS->GetObject(TranchesId[j]);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pMez, ARM_SECURITY) == 0)	// must implement GetRootName for 
			{
				result.setMsg ("ARM_ERR: one Mezz is not of a good type");
			
				return ARM_KO;
			}

			// do I have to clone?
			arraypMez.push_back(pMez);
		}
		
		// ---------------------------------------------------------------------------------
		ACreditManager->Set_CorrelationCalibrator(Matrix_Description, Matrix_Parameters, arraypMez);
		ACreditManager->Correlation_Calibrate();
		// ---------------------------------------------------------------------------------
		
		result.setDouble(dResult);

		return ARM_OK;
	}

	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	return ARM_KO;

}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//	OBJECT CREATION
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


extern long ICMLOCAL_CreditManager(ARM_result&	result,
								   long objId)
{
	ICM_Pricer*	PrevPricer	= NULL;
	ICM_Pricer*	APricer	= NULL;

	long PricerId	=	0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{

		// ---------------------------------------------------------------------------------
		APricer	=	(ICM_Pricer*) (new CreditManager());
		// ---------------------------------------------------------------------------------

		if (APricer == NULL)
		{
			result.setMsg ("ARM_ERR: Credit Manager is null");
			return ARM_KO;
		}

		
		// ---------------------------------------------------------------------------------
		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			PricerId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)APricer);

			if (PricerId == RET_KO)
			{
				if (APricer)
					delete APricer;
				APricer = NULL;
	
				result.setMsg ("ARM_ERR: Pb with inserting object - APricer");				
				return ARM_KO;
			}

			result.setLong(PricerId);

			return ARM_OK;
		}
		else
		{
			PrevPricer = (ICM_Pricer*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);
			// ICM_CREDIT_MANAGER
//			if ((LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(PrevPricer, ICM_CREDIT_MANAGER) == 1) ||
//					(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(PrevPricer, ICM_PRICER) == 1))
			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(PrevPricer, ICM_PRICER) == 1)
			{
				if (PrevPricer)
				{
					delete PrevPricer;
					PrevPricer = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)APricer, objId);

				return ARM_OK;
			}
			else
			{
				if (APricer)
					delete APricer;
				APricer = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type - APricer");
				return ARM_KO;
			}
		}
	}

	catch (Exception& x)
	{
		x.DebugPrint();

		if (APricer)
			delete APricer;
		APricer = NULL;

		ARM_RESULT();
	}

	catch (...)
	{
		result.setMsg ("ARM_ERR: Credit Manager: unrecognized failure");
		return ARM_KO;
	}

}