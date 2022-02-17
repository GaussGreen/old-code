#pragma warning(disable : 4541)
#pragma warning(disable : 4250)


#include <GP_Base\gpbase\curve.h>
#include "firstToBeIncluded.h"
#include <libCCdate\CCdate.h>
#include <libCCtools++\CCstring.h>

#include <GP_Base\gpbase\checkarg.h>

#include <ARM\local_xlarm\ARM_local_interglob.h>

#include <ARM\libarm\ARM_result.h>
#include "ARM_local_class.h"
#include "ARM_local_wrapper.h"

#include "ARM_local_persistent.h"
#include "ARM_local_glob.h"
#include <inst\portfolio.h>
#include <inst\swaption.h>

extern "C"
{
	#include <util\ournag.h>
	#include <nag_stdlib.h>
	#include <nage04.h>
	#include <inst\optqn.h>
	#include <inst\minimize.h>
}


#include <crv\zerovsk.h>
#include <crv\zerospli.h>
#include <mod\x1fhw.h>
#include <mod\x1fhwls.h>
#include <mod\x2fhw.h>
#include <mod\x1fhwsigvar.h>
#include <mod\irtbk.h>
#include <mod\bootstrapcalibration.h>
#include <mod\calibrationlogdec.h>

#include <util\fromto.h>

/*---  VSK Initial Params ---*/
 
#define R_INI           0.1
#define S_INI           0.0
#define G_INI           -0.001
 
#define A_INI           0.25
 
/*---  SPL Initial Params ---*/
 
#define SP1_INI         -0.033067178
#define SP2_INI         -0.001159509
#define SP3_INI         -0.000374861
#define SP4_INI         0.0007462389
#define SP5_INI         -0.000498036
#define SP6_INI         0.0001187675


/*---  HW Initial Model Params ---*/

#define HW_A_INI        0.03
#define HW_S_INI        1.2

#define HW2F_A_INI        0.53
#define HW2F_S_INI        1.35
#define HW2F_DA_INI        0.504
#define HW2F_DS_INI        2.29
#define HW2F_DC_INI        -0.717
 
#define BK_A_INI        0.22
#define BK_S_INI        0.25

 
#define GYLS_A_INI        0.53
#define GYLS_S_INI        1.35
#define GYLS_X_INI        2.0
#define GYLS_XV_INI        20.0
#define GYLS_XC_INI        0.7
#define GYLS_LR_INI        0.0


double* SV_VECT_ABS(int nbPar, double* res)
{
    int k;


    for (k = 1; k < nbPar; k++)
    {
        res[k] = fabs(res[k]);
    }

    return(res);
} 


long ARMLOCAL_PF  (const VECTOR<long>& insts,
				   const VECTOR<double>& coeffs,
				   const VECTOR<double>& marketPrices,
				   const VECTOR<double>& precisions,
				   ARM_result& result,
				   long objId)
{
	long portId;
	int i;
	ARM_StdPortfolio* createdPF=NULL;
	ARM_StdPortfolio* pf=NULL;
	ARM_Security*  curAsset=NULL;

	double* dCoefs;
	double* dMarketPrices;
	double* dPrecisions = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		int n = insts.size ();

        if(!((n == coeffs.size ()) && (n == marketPrices.size ())))
        {
            result.setMsg ("ARM_ERR: instruments, coefficients and prices array must have same size");
            return ARM_KO;
        }

		if ( (precisions.size() > 0) && (!(n == precisions.size())) )
        {
            result.setMsg ("ARM_ERR: precisions, instruments, coefficients and prices array must have same size");
            return ARM_KO;
        }

		dCoefs = new double[n];
		dMarketPrices = new double[n];

		if (precisions.size() > 0)
			dPrecisions  = new double[n];

		for (i=0;i<n;i++)
		{
			dCoefs[i] = coeffs[i];
			dMarketPrices[i] = marketPrices[i];
			if (dPrecisions)
				dPrecisions[i] = precisions[i];
		}

		createdPF = new ARM_StdPortfolio(n, dCoefs, dMarketPrices,dPrecisions);

		for (i = 0; i < n; i++)
		{
			curAsset = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(insts[i]);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(curAsset, ARM_SECURITY) == 0)
			{
				if (dCoefs)
					delete [] dCoefs;
				dCoefs = NULL;

				if (dMarketPrices)
					delete [] dMarketPrices;
				dMarketPrices = NULL;

				if (dPrecisions)
					delete [] dPrecisions;
				dPrecisions = NULL;

				if (createdPF)
					delete createdPF;
				createdPF = NULL;

				result.setMsg("ARM_ERR: instruments have to be a security to constitute a portfolio");
				return ARM_KO;
			}

			createdPF->SetAsset((ARM_Security*)curAsset->Clone(), i);
		}

		if (dCoefs)
			delete [] dCoefs;
		dCoefs = NULL;

		if (dMarketPrices)
			delete [] dMarketPrices;
		dMarketPrices = NULL;

		if (dPrecisions)
			delete [] dPrecisions;
		dPrecisions = NULL;

		if (createdPF == NULL)
		{
			result.setMsg ("ARM_ERR: Portfolio is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
 
			portId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(createdPF);

			if (portId == RET_KO)
			{
				if (createdPF)
					delete createdPF;
				createdPF = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(portId);

			return ARM_OK;
		}
		else
		{
			pf = (ARM_StdPortfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf, ARM_PORTFOLIO) == 1)
			{
				if (pf)
				{
					delete pf;
					pf = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdPF, objId);

				return ARM_OK;
			}

			else
			{
				if (createdPF)
					delete createdPF;
				createdPF = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}
	}
    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (pf)
			delete pf;
		pf = NULL;

		if (dCoefs)
			delete [] dCoefs;
		dCoefs = NULL;

		if (dMarketPrices)
			delete [] dMarketPrices;
		dMarketPrices = NULL;

		if (dPrecisions)
			delete [] dPrecisions;
		dPrecisions = NULL;

		ARM_RESULT();
		return ARM_KO;
    }
}


////////////////////////////////////////////
//// Function to build a portfolio
////////////////////////////////////////////
extern long ARMLOCAL_PFFILL_Create(const VECTOR<long >& assetsIdVec,
		 VECTOR<double >& weightsVec,
         VECTOR<double >& mktpricesVec,
         VECTOR<double >& vegasVec,
		long portfolioId,
         ARM_result&	result, 
         long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_StdPortfolio* createdPF=NULL;
	try
	{   
        CC_NS(ARM_Check,CheckSameArgSize)(assetsIdVec, weightsVec, "AssetsArray", "WeightsArray");
        CC_NS(ARM_Check,CheckSameArgSize)(assetsIdVec, mktpricesVec, "AssetsArray", "WeightsArray");
        CC_NS(ARM_Check,CheckSameArgSize)(assetsIdVec, vegasVec, "AssetsArray", "PrecisionsArray");

        int size = assetsIdVec.size();
        vector<ARM_Security*> assets(size, NULL);
        int i;
        for (i = 0; i < size; i++)
		{
            if( !GetObjectFromId( &assets[i], assetsIdVec[i], ARM_SECURITY) )
	        {
		        result.setMsg ("ARM_ERR: Products Vector is not of a good type, please check all products");
		        return ARM_KO;
	        }
        }

        
        if (portfolioId != ARM_NULL_OBJECT)
        {
            ARM_StdPortfolio* PFToFill=NULL;  
            if( !GetObjectFromId( &PFToFill, portfolioId, ARM_PORTFOLIO) )
	        {
		        result.setMsg ("ARM_ERR: Potfolio to fill is not of a good type");
		        return ARM_KO;
	        }
            if(!PFToFill->GetPrecision())
            {
		        result.setMsg ("ARM_ERR: PotfolioId has a null precision Vector");
		        return ARM_KO;
	        }

            for (i = 0; i < PFToFill->size(); i++)
		    {
                assets.push_back((ARM_Security*)PFToFill->GetAsset(i)); /// will be cloned by the constructor !
                weightsVec.push_back((*(PFToFill->GetWeights()))[i]);
                mktpricesVec.push_back((*PFToFill->GetMktPrices())[i]);
                vegasVec.push_back((*PFToFill->GetPrecision())[i]);            
            }
        }

        createdPF = new ARM_StdPortfolio(assets,weightsVec,mktpricesVec,vegasVec);
        createdPF->sort();


        /// assign object
		if( !assignObject( createdPF, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete createdPF;
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete createdPF;
		result.setMsg ("ARM_ERR: unrecognized failure in Creating A portfolio");
		return ARM_KO;
	}
}


////////////////////////////////////////////
//// Function to build a portfolio
////////////////////////////////////////////
extern long ARMLOCAL_PF_Merge(
		 VECTOR<long> pfIds,
         ARM_result&	result,
         long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_StdPortfolio* createdPF=NULL;
	try
	{
		vector<ARM_Security*> assets;
		vector<double> weights;
		vector<double> mktPrices;
		vector<double> precisions;
		ARM_StdPortfolio* pf=NULL;
		size_t pfIdx,i;
		for(pfIdx=0;pfIdx<pfIds.size();++pfIdx)
		{
            if( pfIds[pfIdx] != ARM_NULL_OBJECT && !GetObjectFromId( &pf, pfIds[pfIdx], ARM_PORTFOLIO) )
	        {
		        result.setMsg ("ARM_ERR: Potfolio is not of a good type");
		        return ARM_KO;
	        }
			for(i=0;i<pf->size();++i)
			{
				assets.push_back((ARM_Security*)pf->GetAsset(i)); /// will be cloned by the constructor !
				weights.push_back((*(pf->GetWeights()))[i]);
				mktPrices.push_back((*pf->GetMktPrices())[i]);
				precisions.push_back((*pf->GetPrecision())[i]);            
			}
		}

        createdPF = new ARM_StdPortfolio(assets,weights,mktPrices,precisions);
        createdPF->sort();


        /// assign object
		if( !assignObject( createdPF, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete createdPF;
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete createdPF;
		result.setMsg ("ARM_ERR: unrecognized failure in merging portfolios");
		return ARM_KO;
	}
}








long ARMLOCAL_PFGYCSIGVARFIT (long pfId,
							  long curveId,
							  const VECTOR<double>& matCurve,
							  double in_min_meanRev,
							  double in_max_meanRev,
							  double in_min_vol,
							  double in_max_vol,
							  double in_precision_meanRev,
							  double in_precision_vol,
							  long nbMaxIter,
							  ARM_result& result,
							  long objId)
{
	long modId;

	VECTOR<CCString> matCurve_str;

	ARM_Portfolio* pf = NULL;
	ARM_ZeroCurve* zc = NULL;
	ARM_Security* curAsset = NULL;
	ARM_GYCSigVarModel* mod = NULL;
	ARM_GYCSigVarModel* prevMod = NULL;

	ARM_GYCSigVarModel* SigVarModel = NULL;

	int sz = matCurve.size ();
	double res, optimal_meanrev;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* sDate = new char[11];

	CCString msg ("");

	try
	{
		for (int i = 0; i < matCurve.size (); i++)
		{
			Local_XLDATE2ARMDATE(matCurve[i],sDate);
			matCurve_str.push_back(sDate);
		}

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		pf = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(pfId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf, ARM_PORTFOLIO) == 0)
		{
			result.setMsg ("ARM_ERR: Portfolio is not of a good type");
			return ARM_KO;
		}
		
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(curveId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		double X[200];
		double Y[200];
		double vals_curve[200];

		MEMSET(X, 0.0, sizeof(double)*200);   
		MEMSET(Y, 0.0, sizeof(double)*200);   
		MEMSET(vals_curve, 0.0, sizeof(double)*200);   

		for (i = 0; i < sz; i++)
		{
			ARM_Date tmpDate(matCurve_str[i]);

			double dDate = (tmpDate.GetJulian() - zc->GetAsOfDate().GetJulian())/K_YEAR_LEN;

			X[i] = dDate;
		}

		SigVarModel = new ARM_GYCSigVarModel(in_min_meanRev, sz, X, Y, zc); 

		long PFSize = pf->GetSize();

		for (i = 0; i < PFSize; i++)
		{
		   curAsset = pf->GetAsset(i);

		   curAsset->SetModel(SigVarModel);
		}

		res = ComputeMeanRevAndParasCurve(SigVarModel,
										   pf,
										   sz,
										   X,
										   in_precision_meanRev,
										   in_precision_vol,
										   in_min_meanRev,
										   in_min_vol,
										   in_max_meanRev,
										   in_max_vol,
										   nbMaxIter,
										   vals_curve);
		optimal_meanrev = res;

		if (SigVarModel)
			delete SigVarModel;
		SigVarModel = NULL;

		mod = new ARM_GYCSigVarModel(optimal_meanrev, 
                             sz, X, vals_curve, zc);

		if (mod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
 
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(mod);

			if (modId == RET_KO)
			{
				if (mod)
					delete mod;
				mod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			prevMod = (ARM_GYCSigVarModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevMod, ARM_GYCSIGVARMODEL) == 1)
			{
				if (prevMod)
				{
					delete prevMod;
					prevMod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)mod, objId);

				return ARM_OK;
			}

			else
			{
				if (mod)
					delete mod;
				mod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (sDate)
			delete [] sDate;
		sDate = NULL;

		if (SigVarModel)
			delete SigVarModel;
		SigVarModel = NULL;

		if (mod)
			delete mod;
		mod = NULL;

		ARM_RESULT();
		return ARM_KO;
	}
}



long ARMLOCAL_PFLOGDECVOLFIT(long pfId,
							 long curveId,
							 const VECTOR<double>& proba,
							 double theAccuracy,
							 long shapeType,
							 double decay,
							 double slope,
							 double asymptote,
							 const VECTOR<double>& matCurve,
							 const VECTOR<double>& volinit_vect,
							 const VECTOR<double>& coeff_vect,
							 ARM_result& result,
							 long objId)
{
	long modId;

    /*--- parameters checking ---*/

	if (shapeType)
	{
		if (volinit_vect.size () != proba.size ())
		{
			result.setMsg ("ARM_ERR: Sizes of Initial Vols and Proba Vectors must be equal");
			return(ARM_KO);
		}
	}
	else
	{
		if (volinit_vect.size () != matCurve.size ())
		{
			result.setMsg ("ARM_ERR: Sizes of Initial Vols and Maturities Vectors must be equal");
			return(ARM_KO);
		}
	}

	if ( volinit_vect.size () != coeff_vect.size () )
	{
		result.setMsg ("ARM_ERR: Sizes of Initial Vols and Coeffs. Vectors must be equal");
		return(ARM_KO);
	}

    ARM_Portfolio* pf = NULL;
    ARM_ZeroCurve* zc = NULL;
    double dx[ARM_NB_TERMS];

    ARM_CalibrationLOGDEC* CalibMod = NULL;

	ARM_Vector* terms=NULL;
	ARM_Vector* probas=NULL;
	ARM_Vector* mVol=NULL;
	ARM_Vector* mCoefVol=NULL;

	ARM_LogDecalANA* createdLDCANA=NULL;
	ARM_LogDecalANA* newMod = NULL;
    
	int i, volsize;

	int sz = matCurve.size ();

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* sDate = new char[11];

	CCString msg ("");

	try
	{
        long realMatCurveSize = 0;
		VECTOR<CCString> matCurve_str;

        if ( matCurve.size () != 0 )
        {
			realMatCurveSize = matCurve.size ();

            for (int i = 0; i < matCurve.size (); i++)
            {
				Local_XLDATE2ARMDATE(matCurve[i],sDate);
				matCurve_str.push_back(sDate);
            }
        }
        else
        {
            matCurve_str.push_back("NULLDATE");  
        }

		if (sDate)
			delete [] sDate;
		sDate = NULL;

		pf = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(pfId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf, ARM_PORTFOLIO) == 0)
		{
			result.setMsg ("ARM_ERR: Portfolio is not of a good type");
			return ARM_KO;
		}
		
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(curveId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		ARM_Date dateDeb = zc->GetAsOfDate();

		if (!realMatCurveSize)
		{
			terms = NULL;
			mVol = NULL;
			mCoefVol = NULL;
		}
		else
		{
			for (i = 0; i < realMatCurveSize; i++)
			{
				sDate = (char*)matCurve_str[i];
				ARM_Date tmpDate(sDate);
				delete sDate;

				double dDate = (tmpDate.GetJulian()-dateDeb.GetJulian())
							   / K_YEAR_LEN;

				dx[i] = dDate;
			}
			terms = new ARM_Vector(realMatCurveSize, dx);
		}

		probas = CreateARMVectorFromVECTOR(proba);
		
		// test vol initiale positive (config : proba geree)
		if (volinit_vect.size() > 0)
		{
			if (volinit_vect[0] > 0.0)
			{
				if (shapeType)
					volsize = proba.size();
				else
					volsize = realMatCurveSize;

				mVol = new ARM_Vector(volsize);
				mCoefVol = new ARM_Vector(volsize);

				for (i = 0; i < volsize; i++)
				{
					mVol->Elt(i) = volinit_vect[i];
					mCoefVol->Elt(i) = coeff_vect[i];
				}
			}
			else
			{
				mVol = NULL;
				mCoefVol = NULL;
			}
		}
		else
		{
			mVol = NULL;
			mCoefVol = NULL;
		}

		CalibMod = new ARM_CalibrationLOGDEC(zc, terms, probas, pf, theAccuracy, 
											 shapeType, decay, slope, asymptote,
											 mVol, mCoefVol);

		// ce modele issu du calibrage ne sert que pour recuperer les parametres
		// ne pas utiliser tel quel car la ZeroCurve est detruite dans CalibMod
		createdLDCANA = CalibMod->Calibrate();

		if (terms)
		{
			delete terms;
			terms = NULL;
		}
		terms = (ARM_Vector *) createdLDCANA->GetResetYearTerms()->Clone();

		ARM_Vector *shift = NULL;
		shift = (ARM_Vector *) createdLDCANA->Shifts()->Clone();

		ARM_Matrix *Vols = NULL;
		Vols = (ARM_Matrix *) createdLDCANA->FwdVol()->Clone();

		ARM_ShapeVector *newShape = NULL;
		newShape = (ARM_ShapeVector *) createdLDCANA->ShapeVector()->Clone();

		newMod = (ARM_LogDecalANA *) new ARM_LogDecalANA(zc, 1, terms, shift,
														 Vols);
		newMod->SetShapeVector(newShape);

/*		if (terms)
			delete terms;
		terms = NULL;
*/
		if (probas)
			delete probas;
		probas = NULL;

		if (mVol)
			delete mVol;
		mVol = NULL;

		if (mCoefVol)
			delete mCoefVol;
		mCoefVol = NULL;

		if (createdLDCANA)
			delete createdLDCANA;
		createdLDCANA = NULL;

		if (CalibMod)
			delete CalibMod;
		CalibMod = NULL;

		if (newMod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
 
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(newMod);

			if (modId == RET_KO)
			{
				if (newMod)
					delete newMod;
				newMod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			ARM_LogDecalANA* prevMod = dynamic_cast<ARM_LogDecalANA *> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevMod, ARM_LOGDECALANA) == 1)
			{
				if (prevMod)
				{
					delete prevMod;
					prevMod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newMod, objId);

				return ARM_OK;
			}

			else
			{
				if (newMod)
					delete newMod;
				newMod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (sDate)
			delete sDate;
		sDate = NULL;

		if (terms)
			delete terms;
		terms = NULL;

		if (probas)
			delete probas;
		probas = NULL;

		if (mVol)
			delete mVol;
		mVol = NULL;

		if (mCoefVol)
			delete mCoefVol;
		mCoefVol = NULL;

		if (createdLDCANA)
			delete createdLDCANA;
		createdLDCANA = NULL;

		if (CalibMod)
			delete CalibMod;
		CalibMod = NULL;

		if (newMod)
			delete newMod;
		newMod = NULL;

		ARM_RESULT();

		return ARM_KO;

	}
}



long ARMLOCAL_PFGYCSIGVARPENALFIT(long pfId,
								  long curveId,
								  const VECTOR<double>& matCurve,
								  double start_meanRev,
								  const VECTOR<double>& start_vol_vect,
								  const VECTOR<double>& penal_vect,
								  const VECTOR<double>& coeff_vect,
								  double theAccuracy,
								  ARM_result& result,
								  long objId)
{
    ARM_Portfolio* pf=NULL;
    ARM_GYCSigVarModel* mod = NULL;
    ARM_ZeroCurve* zc = NULL;
    ARM_Security* curAsset=NULL;
 
    double initVect[200];
 
    long modId;
 
	if ( matCurve.size () != start_vol_vect.size () )
	{
		result.setMsg ("ARM_ERR: Sizes of Maturities and Volatilities Vectors must be equal");
		return(ARM_KO);
	}

	if ( (penal_vect.size ()-1) != start_vol_vect.size () )
	{
		result.setMsg ("ARM_ERR: Sizes of Penalty and Volatilities Vectors must be equal");
		return(ARM_KO);
	}

	if ( penal_vect.size () != coeff_vect.size () )
	{
		result.setMsg ("ARM_ERR: Sizes of Penalty and Coeffs. Vectors must be equal");
		return(ARM_KO);
	}

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	VECTOR<CCString> matCurve_str;
	char* sDate = new char[11];

	int szVol = start_vol_vect.size();
	int sz = matCurve.size();

	CCString msg ("");
	
	try
	{
		pf = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(pfId);
 
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf, ARM_PORTFOLIO) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			result.setMsg ("ARM_ERR: Portfolio is not of a good type");
			return ARM_KO;
		}
 
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(curveId);
 
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			if (sDate)
				delete [] sDate;
			sDate = NULL;

			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		for (int i = 0; i < matCurve.size (); i++)
		{
			Local_XLDATE2ARMDATE(matCurve[i],sDate);
			matCurve_str.push_back(sDate);
		}

		double X[200];
		double Y[200];

		MEMSET(X, 0.0, sizeof(double)*200);

		if ( szVol == 0 )
			MEMSET(Y, HW_S_INI, sizeof(double)*200);
		else if ( szVol == 1 )
			MEMSET(Y, start_vol_vect[0], sizeof(double)*200);
		else
		{
			if ( sz != szVol )
			{
				throw Exception(__LINE__, __FILE__, ERR_NUMERICAL_CALL,
				"Problem in PFHWSigVarPenalFit : Nb Start. Vol. Values Different from NbPlots!! ");
			}

			for (i = 0; i < szVol; i++)
				Y[i] = start_vol_vect[i];
		}

		ARM_Date beginDate = zc->GetAsOfDate();

		for (i = 0; i < sz; i++)
		{
			sDate = matCurve_str[i];
			ARM_Date tmpDate(sDate);
			delete sDate;

			double dDate = (tmpDate.GetJulian()-beginDate.GetJulian())/K_YEAR_LEN;

			X[i] = dDate;
		}
		sDate = NULL;

		mod = new ARM_GYCSigVarModel(start_meanRev, sz, X, Y, zc);

		long PFSize = pf->GetSize();

		for (i = 0; i < PFSize; i++)
		{
			curAsset = pf->GetAsset(i);
			curAsset->SetModel(mod);
		}

		char settl[20];
		beginDate.JulianToStrDate(settl);

		initVect[0] = start_meanRev;

		for ( i = 1; i < sz+1; i++)
			initVect[i] = Y[i-1];

		mod->SetParams(initVect);

		double* dPenalVect = new double[szVol];
		double* dCoefVect = new double[szVol];

		for (i=0;i<szVol;i++)
		{
			dPenalVect[i] = penal_vect[i];
			dCoefVect[i] = coeff_vect[i];
		}

		long retCode = quasi_N((void *) pf, mod->GetNbParams(), settl,
							   initVect, dPenalVect, dCoefVect, theAccuracy);

		if (dPenalVect)
			delete [] dPenalVect;
		dPenalVect = NULL;

		if (dCoefVect)
			delete [] dCoefVect;
		dCoefVect = NULL;

		if (retCode)
		{
		   if (mod)
			  delete mod;

		   mod = NULL;

		   throw Exception(__LINE__, __FILE__, ERR_NUMERICAL_CALL,
					"Problem in PFHWSigVarPenalFit");
		}

		mod->SetParams(SV_VECT_ABS(mod->GetNbParams(), initVect));

		if (mod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
 
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(mod);

			if (modId == RET_KO)
			{
				if (mod)
					delete mod;
				mod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			ARM_GYCSigVarModel* prevMod = (ARM_GYCSigVarModel*) (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevMod, ARM_GYCSIGVARMODEL) == 1)
			{
				if (prevMod)
				{
					delete prevMod;
					prevMod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)mod, objId);

				return ARM_OK;
			}

			else
			{
				if (mod)
					delete mod;
				mod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

    catch(Exception& x)
    {
        x.DebugPrint();
 
		if (sDate)
			delete [] sDate;
		sDate = NULL;

		if (mod)
			delete mod;
		mod = NULL;

		ARM_RESULT();
		return ARM_KO;
	}
}



long ARMLOCAL_PFINSTLOGDECVOLFIT(long pfId,
								 long secId,
								 long curveId,
								 const VECTOR<double>& proba,
								 long UsePFResetDates,
								 double theAccuracy,
								 long shapeType,
								 double decay,
								 double slope,
								 double asymptote,
								 long VolBSId,
								 const VECTOR<double>& matCurve,
								 const VECTOR<double>& volinit_vect,
								 const VECTOR<double>& coeff_vect,
								 ARM_result& result,
								 long objId)
{
	long modId;

	ARM_VolCurve* volBSInit = NULL;
	ARM_Portfolio* pf = NULL;
	ARM_Security* Sec = NULL;
	ARM_ZeroCurve* zc = NULL;
	double dx[ARM_NB_TERMS];

	ARM_CalibrationLOGDEC* CalibMod = NULL;
	ARM_LogDecalANA* newMod = NULL;

	if (shapeType)
	{
		if (volinit_vect.size () != proba.size ())
		{
			result.setMsg ("ARM_ERR: Sizes of Initial Vols and Proba Vectors must be equal");
        
			return(ARM_KO);
		}
	}
	else
	{
		if (volinit_vect.size () != matCurve.size ())
		{
			result.setMsg ("ARM_ERR: Sizes of Initial Vols and Maturities Vectors must be equal");
        
			return(ARM_KO);
		}
	}

	if ( volinit_vect.size () != coeff_vect.size () )
	{
		result.setMsg ("ARM_ERR: Sizes of Initial Vols and Coeffs. Vectors must be equal");

		return(ARM_KO);
	}

	VECTOR<CCString> matCurve_str;

	long realMatCurveSize = 0;

	char* strDate = new char[11];

	if ( matCurve.size () != 0 )
	{
		realMatCurveSize = matCurve.size ();

		for (int i = 0; i < matCurve.size (); i++)
		{
			Local_XLDATE2ARMDATE(matCurve[i],strDate);
			matCurve_str.push_back(strDate);
		}
	}
	else
	{
		matCurve_str.push_back("NULLDATE");  
	}

	if (strDate)
		delete [] strDate;
	strDate = NULL;

	int i(0);

	CCString msg(" ");

	try
	{
		pf = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(pfId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf, ARM_PORTFOLIO) == 0)
		{
			result.setMsg ("ARM_ERR: Portfolio is not of a good type");
			return ARM_KO;
		}

		Sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(Sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: Security is not of a good type");
			return ARM_KO;
		}

		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(curveId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: ZC is not of a good type");
			return ARM_KO;
		}

		if (VolBSId != -1)
		{
			volBSInit = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(VolBSId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(volBSInit, ARM_VOL_CURVE) == 0)
			{
				result.setMsg ("ARM_ERR: Vol BS is not of a good type");
				return ARM_KO;
			}
		}

		ARM_Date dateDeb = zc->GetAsOfDate();

		ARM_Vector* terms;
		ARM_Vector* probas;
		ARM_Vector* mVol;
		ARM_Vector* mCoefVol;

		if (!realMatCurveSize)
		{
			terms = NULL;
			mVol = NULL;
			mCoefVol = NULL;
		}
		else
		{
			for (i = 0; i < realMatCurveSize; i++)
			{
				ARM_Date tmpDate(matCurve_str[i]);

				double dDate = (tmpDate.GetJulian()-dateDeb.GetJulian())
							   / K_YEAR_LEN;

				dx[i] = dDate;
			}
			terms = new ARM_Vector(realMatCurveSize, dx);
		}

		long probasize (proba.size());
		long volsize(0);

		probas = new ARM_Vector(probasize, 0.0);
		for (i=0;i<proba.size();i++)
			probas->Elt(i) = proba[i];

		// test vol initiale positive (config : proba geree)
		if (volinit_vect.size() > 0)
		{
			if (shapeType)
				volsize = probasize;
			else
				volsize = realMatCurveSize;

			mVol = new ARM_Vector(volsize, 0.);
			mCoefVol = new ARM_Vector(volsize, 0.);

			for (i=0;i<volsize;i++)
			{
				mVol->Elt(i) = volinit_vect[i];
				mCoefVol->Elt(i) = coeff_vect[i];
			}
		}
		else
		{
			mVol = NULL;
			mCoefVol = NULL;
		}

        CalibMod = new ARM_CalibrationLOGDEC(zc, terms, probas, pf, Sec,
                                             UsePFResetDates, theAccuracy,
                                             shapeType, decay, slope, asymptote,
                                             mVol, mCoefVol, volBSInit);

        ARM_LogDecalANA* createdLDCANA = CalibMod->Calibrate();

        if (terms)
        {
            delete terms;
            terms = NULL;
        }
        terms = (ARM_Vector *) createdLDCANA->GetResetYearTerms()->Clone();

        ARM_Vector *shift = NULL;
        shift = (ARM_Vector *) createdLDCANA->Shifts()->Clone();

        ARM_Matrix *Vols = NULL;
        Vols = (ARM_Matrix *) createdLDCANA->FwdVol()->Clone();

        ARM_ShapeVector *newShape = NULL;
        newShape = (ARM_ShapeVector *) createdLDCANA->ShapeVector()->Clone();

        newMod = (ARM_LogDecalANA *) new ARM_LogDecalANA(zc, 1, terms, shift,
                                                         Vols);
        newMod->SetShapeVector(newShape);

		if (probas)
			delete probas;
		probas = NULL;

		if (mVol)
			delete mVol;
		mVol = NULL;

		if (mCoefVol)
			delete mCoefVol;
		mCoefVol = NULL;

		if (createdLDCANA)
			delete createdLDCANA;
		createdLDCANA = NULL;

		if (CalibMod)
			delete CalibMod;
		CalibMod = NULL;

		if (newMod == NULL)
		{
			result.setMsg ("ARM_ERR: Model is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
 
			modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(newMod);

			if (modId == RET_KO)
			{
				if (newMod)
					delete newMod;
				newMod = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(modId);

			return ARM_OK;
		}
		else
		{
			ARM_LogDecalANA* prevMod = dynamic_cast<ARM_LogDecalANA *> (LOCAL_PERSISTENT_OBJECTS->GetObject(objId));

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevMod, ARM_LOGDECALANA) == 1)
			{
				if (prevMod)
				{
					delete prevMod;
					prevMod = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newMod, objId);

				return ARM_OK;
			}

			else
			{
				if (newMod)
					delete newMod;
				newMod = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}

    catch(Exception& x)
    {
        x.DebugPrint();

		if (newMod)
			delete newMod;
		newMod = NULL;
 
		ARM_RESULT();
		return ARM_KO;
	}
}



// ************************************************************************************
// ARM_PFMODFIT functions                                                          
// ************************************************************************************

long CC_LOCAL_PFModEstimateAlgo1(char* mod_name,
								 int pfId,
								 char* settl,
								 int zcId,
								 int size,
								 double* initVect,
								 int* flagEstim,
								 int step,
								 char* horizon,
								 long& retour )
{
    static char* __ARM_FUNC__ = "CC_PFModEstimateAlgo1";

    ARM_Portfolio* pf = NULL;
    ARM_Model* mod = NULL;
    ARM_ZeroCurve* zc = NULL;
    ARM_Security* curAsset = NULL;
    double iniParam[20];
    ARM_Vector* vParam = NULL;
    int CST[20];
    int idMod;
    int i;
    long retCode=0;
    char msg[100];
    ARM_Container* Cont = NULL;


 
    //ARM_BEGIN(pfTrace);

	ARM_result result; 

    if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
     
	MEMSET(CST, 0, sizeof(CST)); 

    pf = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(pfId);
 
    
	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf, ARM_PORTFOLIO) == 0)
	{
		result.setMsg ("ARM_ERR: PORTFOLIO is not of a good type");
		return ARM_KO;
	}

    if ( zcId != ARM_NULL_OBJECT )
    {
       zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

      
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: ZERO_CURVE is not of a good type");
			return ARM_KO;
		}
    }

    try
    {
        if (( strcmp(mod_name, "GYCM") == 0 )
            || 
            ( strcmp(mod_name, "HW2F") == 0 )
            ||
            ( strcmp(mod_name, "BKTR") == 0 )
            ||
            ( strcmp(mod_name, "GYLS") == 0 )
           )
        {
           if ( zc == NULL )
           {
              throw Exception(__LINE__, __FILE__, ERR_NUMERICAL_CALL,
                           "Need a Yield Curve as input");
           }
        }
 
        if ( strcmp(mod_name, "GYCM") == 0 )
        {
           iniParam[0] = HW_A_INI;
           iniParam[1] = HW_S_INI;

           for (i = 0; i < 2; i++)
           {
               if ( i < size )
               {
                  iniParam[i] = initVect[i];

                  CST[i] = flagEstim[i];
               }
               else
                  CST[i] = 0;
           }
 
           mod = (ARM_Model *) new ARM_GYCModel(iniParam[0], iniParam[1], zc);
        }
        else if ( strcmp(mod_name, "BKTR") == 0 )
        {
           iniParam[0] = BK_A_INI;
           iniParam[1] = BK_S_INI;
 
           for (i = 0; i < 2; i++)
           {
               if ( i < size )
               {
                  iniParam[i] = initVect[i];

                  CST[i] = flagEstim[i];
               }
               else
                  CST[i] = 0;
            }
 
            mod = (ARM_Model *) new ARM_BKIRTree((ARM_Date) settl ,
                            (ARM_Date) horizon,
                            step,
                            iniParam[0], iniParam[1], zc);
 
        }
        else if ( strcmp(mod_name, "HW2F") == 0 )
        {
           iniParam[0] = HW2F_A_INI;
           iniParam[1] = HW2F_S_INI;
           iniParam[2] = HW2F_DA_INI;
           iniParam[3] = HW2F_DS_INI;
           iniParam[4] = HW2F_DC_INI;

           for (i = 0; i < 5; i++)
           {
               if ( i < size ) 
               {
                  iniParam[i] = initVect[i];

                  CST[i] = flagEstim[i];
               }
               else
                  CST[i] = 0;
           }

           mod = new ARM_HWTwoFactorModel(iniParam[0], iniParam[1], 
                        iniParam[2], iniParam[3], iniParam[4], zc);
        }
        else if ( strcmp(mod_name, "GYLS") == 0 )
        {
           iniParam[0] = GYLS_A_INI;
           iniParam[1] = GYLS_S_INI;
           iniParam[2] = GYLS_X_INI;
           iniParam[3] = GYLS_XV_INI;
           iniParam[4] = GYLS_XC_INI;;
           iniParam[5] = GYLS_LR_INI;;

           for (i = 0; i < 6; i++)
           {
               if ( i < size ) 
               {
                  iniParam[i] = initVect[i];

                  CST[i] = flagEstim[i];
               }
               else
                  CST[i] = 0;
           } 

           mod = new ARM_GYCLSModel(iniParam[0], iniParam[1], zc,
                        iniParam[2], iniParam[3], iniParam[4], iniParam[5]);
        }
        else if ( strcmp(mod_name, "ZCVS") == 0 )
        {
           iniParam[0] = A_INI;
           iniParam[1] = R_INI;
           iniParam[2] = S_INI;
           iniParam[3] = G_INI;

           for (i = 0; i < 4; i++)
           {
               if ( i < size )
               {
                  iniParam[i] = initVect[i];

                  CST[i] = flagEstim[i];
               }
               else
                  CST[i] = 0;
            }

            ARM_ZeroVasicek* lzc = new ARM_ZeroVasicek((ARM_Date) settl, 
                                                       iniParam);

            mod = (ARM_Model *) new ARM_YCModel(lzc);

            //mod->SetNoSharedObjects();
			//idMod=LOCAL_PERSISTENT_OBJECTS->SetPersistent(mod);
        }
        else if ( strcmp(mod_name, "ZCSP") == 0 )
        {
           iniParam[0] = SP1_INI;
           iniParam[1] = SP2_INI;
           iniParam[2] = SP3_INI;
           iniParam[3] = SP4_INI;
           iniParam[4] = SP5_INI;
           iniParam[5] = SP6_INI;
 
           for (i = 0; i < 6; i++)
           {
               if ( i < size )
               {
                  iniParam[i] = initVect[i];
                  CST[i] = flagEstim[i];
               }
               else
                  CST[i] = 0;
           }

           ARM_ZeroSplines* lzc = new ARM_ZeroSplines((ARM_Date) settl, 
                                        6, iniParam);

           mod = (ARM_Model *) new ARM_YCModel(lzc);

           mod->SetNoSharedObjects();
		   //idMod=LOCAL_PERSISTENT_OBJECTS->SetPersistent(mod);
        }
        else
        {
           if (mod)
           {
              delete mod;

              mod = NULL;
           }

           throw Exception(__LINE__, __FILE__, ERR_OBJECT_NULL,
                           "Problem in Model Name");
        }

        int sz = pf->GetSize();
 
 
        ARM_Date date = (ARM_Date) settl;
 
        Cont = pf->PfToCont();

        double OP[19];

        //MEMSET(OP, 0, sizeof(OP)*sizeof(double));
		MEMSET(OP, 0, sizeof(OP));

        ARM_Matrix* Data = pf->PfToMarketData();

        retCode = mod->FitToMarketData(Cont, Data, CST, 0, 0, OP);

        if (Cont)
        {
           delete Cont;
           Cont = NULL;
        }

        if (Data) 
        {
           if (Data)
              delete Data;

           Data = NULL;
        }
 
        if ( retCode >= 3 )
        {
           if (mod)
           {
              delete mod;

              mod = NULL;
           }

           throw Exception(__LINE__, __FILE__, ERR_OBJECT_NULL,
                           "Problem in Fitting");
        }
        else
        {
           if ( retCode == 0 )
              strcpy(msg, "Optimum Reached !");
 
           if ( retCode == 1 )
              strcpy(msg, "Maximum Iterations Reached !");
 
           if (retCode == 2)
              strcpy(msg,
                     "Optimum conditions not satisfied but no better solution found locally !");
        }

        /*---- The result ----*/
 
		/*
        FILE* MOD = ARM_GetOutputHomeFile("MOD", "a+");
 

        fprintf(MOD, "Model de Type %s \n", mod_name);

        ARM_Vector* Parametre;

        Parametre = mod->GetParameters();

        for (i=0; i < mod->GetNbParams(); i++)
        {
            fprintf(MOD, "Param[%d] = %15.10f\n", i, Parametre->Elt(i));
        }

        fclose(MOD);

		*/

        CREATE_GLOBAL_OBJECT();
 
        idMod = LOCAL_PERSISTENT_OBJECTS->SetPersistent(mod);
		retour = idMod;

    }
 
    catch(Exception& x)
    {
        if (mod)
        {
           delete mod;
   
           mod = NULL;
        }
 
        x.DebugPrint();
 
 		//ARM_RESULT();
		return ARM_KO;
    }
 
 
    //ARM_END(pfTrace);
    //rc = RC_GetIntValue(idMod);
 
    //strcpy(rc.msg, msg);
 
    return(ARM_OK);
} 



long CC_LOCAL_PFModEstimate(char* mod_name,
							int pfId,
							char* settl,
							int zcId,
                            int size,
							double* initVect,
							int* flagEstim,
                            int step,
							char* horizon,
							int NAG_ALGO,
							int MOD_ID,
							long & retour)
{
    static char* __ARM_FUNC__ = "CC_LOCAL_PFModEstimate";


    //ARM_BEGIN(pfTrace);

    if ( NAG_ALGO == 0 )
    {
       if (CC_LOCAL_PFModEstimateAlgo1(mod_name, pfId, settl, zcId, 
                                  size, initVect, flagEstim,
                                  step, horizon,retour) == ARM_OK)
       return ARM_OK;
	   else
	   return ARM_KO;
    }

    ARM_Portfolio* pf = NULL;
    ARM_Model* mod = NULL;
    ARM_ZeroCurve* zc = NULL;
    ARM_Security* curAsset = NULL;
    double res[10];
    double iniParam[20];
    ARM_Vector* vParam = NULL;
    int CST[20];
    int idMod;
    int i;


 	ARM_result result; 

 	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

    pf = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(pfId);
 
  
 	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf, ARM_PORTFOLIO) == 0)
	{
		result.setMsg ("ARM_ERR: PORTFOLIO is not of a good type");
		return ARM_KO;
	}

    if ( zcId != ARM_NULL_OBJECT )
    {
       zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);
 
      
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
		   result.setMsg ("ARM_ERR: ZERO_CURVE is not of a good type");
		   return ARM_KO;
		}
    }

    try
    {
        if (( strcmp(mod_name, "GYCM") == 0 )
            ||
            ( strcmp(mod_name, "HW2F") == 0 )
            ||
            ( strcmp(mod_name, "BKTR") == 0 )
            ||
            ( strcmp(mod_name, "GYLS") == 0 )
           )
        {
           if ( zc == NULL )
           {
              throw Exception(__LINE__, __FILE__, ERR_NUMERICAL_CALL,
                           "Fiiting Needs a Yield Curve as input");
           }
        }

        if ( strcmp(mod_name, "GYCM") == 0 )
        {
           iniParam[0] = HW_A_INI;
           iniParam[1] = HW_S_INI;
 
           for (i = 0; i < 2; i++)
           {
               if ( i < size )
               {
                  iniParam[i] = initVect[i];

                  CST[i] = flagEstim[i];
               }
               else
                  CST[i] = 0;
           }
 
           mod = (ARM_Model *) new ARM_GYCModel(iniParam[0], iniParam[1], zc);
        }
        else if ( strcmp(mod_name, "BKTR") == 0 )
        {
           iniParam[0] = BK_A_INI;
           iniParam[1] = BK_S_INI;
 
           for (i = 0; i < 2; i++)
           {
               if ( i < size )
               {
                  iniParam[i] = initVect[i];
 
                  CST[i] = flagEstim[i];
               }
               else
                  CST[i] = 0;
            }
 
            mod = (ARM_Model *) new ARM_BKIRTree((ARM_Date) settl ,
                                                 (ARM_Date) horizon,
                                                 step,
                                                 iniParam[0], iniParam[1], zc);
 
        }
        else if ( strcmp(mod_name, "HW2F") == 0 )
        {
           iniParam[0] = HW2F_A_INI;
           iniParam[1] = HW2F_S_INI;
           iniParam[2] = HW2F_DA_INI;
           iniParam[3] = HW2F_DS_INI;
           iniParam[4] = HW2F_DC_INI;
 
           for (i = 0; i < 5; i++)
           {
               if ( i < size )
               {
                  iniParam[i] = initVect[i];
 
                  CST[i] = flagEstim[i];
               }
               else
                  CST[i] = 0;
           }
 
           mod = new ARM_HWTwoFactorModel(iniParam[0], iniParam[1],
                        iniParam[2], iniParam[3], iniParam[4], zc);
        }
        else if ( strcmp(mod_name, "GYLS") == 0 )
        {
           iniParam[0] = GYLS_A_INI;
           iniParam[1] = GYLS_S_INI;
           iniParam[2] = GYLS_X_INI;
           iniParam[3] = GYLS_XV_INI;
           iniParam[4] = GYLS_XC_INI;;
           iniParam[5] = GYLS_LR_INI;;
 
           for (i = 0; i < 6; i++)
           {
               if ( i < size )
               {
                  iniParam[i] = initVect[i];
 
                  CST[i] = flagEstim[i];
               }
               else
                  CST[i] = 0;
           }
 
           mod = new ARM_GYCLSModel(iniParam[0], iniParam[1], zc,
                        iniParam[2], iniParam[3], iniParam[4], iniParam[5]);
        }
        else if ( strcmp(mod_name, "ZCVS") == 0 )
        {
           iniParam[0] = A_INI;
           iniParam[1] = R_INI;
           iniParam[2] = S_INI;
           iniParam[3] = G_INI;
 
           for (i = 0; i < 4; i++)
           {
               if ( i < size )
               {
                  iniParam[i] = initVect[i];
 
                  CST[i] = flagEstim[i];
               }
               else
                  CST[i] = 0;
            }
 
            ARM_ZeroVasicek* lzc = new ARM_ZeroVasicek((ARM_Date) settl,
                                                       iniParam);
 
            mod = (ARM_Model *) new ARM_YCModel(lzc);
 
            mod->SetNoSharedObjects();
        }
        else if ( strcmp(mod_name, "ZCSP") == 0 )
        {
           iniParam[0] = SP1_INI;
           iniParam[1] = SP2_INI;
           iniParam[2] = SP3_INI;
           iniParam[3] = SP4_INI;
           iniParam[4] = SP5_INI;
           iniParam[5] = SP6_INI;
 
           for (i = 0; i < 6; i++)
           {
               if ( i < size )
               {
                  iniParam[i] = initVect[i];
                  CST[i] = flagEstim[i];
               }
               else
                  CST[i] = 0;
           }
 
           ARM_ZeroSplines* lzc = new ARM_ZeroSplines((ARM_Date) settl,
                                        6, iniParam);
 
           mod = (ARM_Model *) new ARM_YCModel(lzc);
 
           mod->SetNoSharedObjects();
		   //idMod=LOCAL_PERSISTENT_OBJECTS->SetPersistent(mod);
        }
        else
        {
           if (mod)
           {
              delete mod;
 
              mod = NULL;
           }
 
           throw Exception(__LINE__, __FILE__, ERR_OBJECT_NULL,
                           "Problem in Model Name");
        }

        int sz = pf->GetSize();
 
        for (i = 0; i < sz; i++)
        {
           curAsset = pf->GetAsset(i);
 
           curAsset->SetModel(mod);
        }
 
        ARM_Date date = (ARM_Date) settl;
 
        long retCode = EstimateModelWithFlags(mod->GetNbParams(), (void *) pf,
                                     iniParam, settl, res,flagEstim);
 
        if (retCode)
        {
           if (mod)
              delete mod;
 
           mod = NULL;
 
           throw Exception(__LINE__, __FILE__, ERR_NUMERICAL_CALL,
                           "Problem in Model Fitting");
        }

        /*---- The result ----*/

		
 
        mod->SetParams(res);

		/*
 
        FILE* MOD = ARM_GetOutputHomeFile("MOD", "a+");
 
 
        fprintf(MOD, "Model de Type %s \n", mod_name);
 
        for (i = 0; i < mod->GetNbParams(); i++)
        {
            fprintf(MOD, "Param[%d] = %15.10f\n", i, res[i]);
        }
 
        fclose(MOD);
		
		*/

        CREATE_GLOBAL_OBJECT();
 
        if ( MOD_ID < 0 )
        {
           idMod = LOCAL_PERSISTENT_OBJECTS->SetPersistent(mod);
		   retour = idMod;
        }
        else
        {
           idMod = MOD_ID;
 
           (void) LOCAL_PERSISTENT_OBJECTS->SetPersistent(mod, MOD_ID);
		   retour = idMod;
        }

		
    }

    catch(Exception& x)
    {
        if (mod)
        {
           delete mod;
 
           mod = NULL;
        }
 
        x.DebugPrint();
 
		//ARM_RESULT();
		return ARM_KO;
    }
 
 
    //ARM_END(pfTrace);
 
    //return(RC_GetIntValue(idMod));
	return(ARM_OK);
}
 

long CC_LOCAL_SetPFModEstimatedAlgo1(int modId,
									 char* mod_name,
									 int pfId,
									 char* settl,
									 int zcId,
									 int size,
									 double* initVect,
									 int* flagEstim,
									 int step,
									 char* horizon,
									 long & retour)
{
    static char* __ARM_FUNC__ = "CC_LOCAL_SetPFModEstimatedAlgo1";

    ARM_Portfolio* pf=NULL;
    ARM_Model* mod = NULL;
    ARM_Model* tmpmod = NULL;
    ARM_ZeroCurve* zc = NULL;
    ARM_Security* curAsset=NULL;
    double iniParam[20];
    ARM_Vector* vParam = NULL;
    int CST[20];
    int i;
    ARM_Container* Cont = NULL;
    long retCode=0;
    char msg[100];
    ARM_RET_CODE rc;
    
 
     //ARM_BEGIN(pfTrace);
 	ARM_result result; 
    //CHECK_GLOBAL_OBJECT();

 	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

    pf = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(pfId);
 
	
    MEMSET(CST, 0, sizeof(CST));

    mod = (ARM_Model *) LOCAL_PERSISTENT_OBJECTS->GetObject(modId);

   
 	if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(mod, ARM_MODEL) == 0)
	{
		result.setMsg ("ARM_MODEL: MODEL is not of a good type");
		return ARM_KO;
	}

    pf = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(pfId);
 
    
 	if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pf, ARM_PORTFOLIO) == 0)
	{
		result.setMsg ("ARM_ERR: PORTFOLIO is not of a good type");
		return ARM_KO;
	}

    if ( zcId != ARM_NULL_OBJECT )
    {
       zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

       
	   if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
       {
		  result.setMsg ("ARM_ERR: ZERO_CURVE is not of a good type");
		  return ARM_KO;
       }
    }

    try
    {

        if (strcmp(mod_name, "GYCM")==0 || strcmp(mod_name, "HW2F")==0 ||
            strcmp(mod_name, "BKTR")==0)
        {
            if ( zc == NULL )
            {
                throw Exception(__LINE__, __FILE__, ERR_NUMERICAL_CALL,
                           "Need a curve as input");
            }
        }

        if (strcmp(mod_name, "GYCM")==0)
        {
               iniParam[0] = HW_A_INI;
               iniParam[1] = HW_S_INI;

            for (i = 0; i < 2; i++)
            {
                if ( i < size )
                {
                   iniParam[i] = initVect[i];
                   CST[i] = flagEstim[i];
                }
                else
                {
                   CST[i] = 0;
                }
            }

            if ( mod->GetName() != ARM_GYCMODEL )
            {
                tmpmod = new ARM_GYCModel(iniParam[0], iniParam[1], zc);

                if (mod)
                {
                   delete mod;

                   mod = NULL;
                }

                mod = tmpmod;

                LOCAL_PERSISTENT_OBJECTS->SetPersistent(mod, modId);

            }
            else
            {
                if ( initVect != NULL )
                {
                    vParam = new ARM_Vector(2, iniParam);
                    mod->SetParameters(vParam);

                    mod->SetZeroCurve(zc);
                }
                else
                {
                    mod->SetZeroCurve(zc);
                }
            }
        }
        else if (strcmp(mod_name, "BKTR")==0)
        {
               iniParam[0] = BK_A_INI;
               iniParam[1] = BK_S_INI;

            for (i = 0; i < 2; i++)
            {
                if ( i < size ) 
                {
                   iniParam[i] = initVect[i];

                   CST[i] = flagEstim[i];
                }
                else
                   CST[i] = 0;
            }

            if (mod->GetName() != ARM_BKIRTREE
                || ((ARM_PDESolver *) mod)->GetNumSteps() != step
                || ((ARM_PDESolver *) mod)->GetHorizon()  != (ARM_Date) horizon)
            {
                tmpmod = new ARM_BKIRTree((ARM_Date) settl ,
                                          (ARM_Date) horizon, step,
                                          iniParam[0], iniParam[1], zc);

                if (mod)
                {
                   delete mod;

                   mod = NULL;
                }

                mod = tmpmod;

                LOCAL_PERSISTENT_OBJECTS->SetPersistent(mod, modId);

            }
            else
            {
                if (initVect != NULL)
                {
                    vParam = new ARM_Vector(2, iniParam);
                    mod->SetParameters(vParam);
                    // delete vParam; Not Cloned !

                    mod->SetZeroCurve(zc);
                }
                else
                    mod->SetZeroCurve(zc);
            }

        }
        else if (strcmp(mod_name, "HW2F")==0)
        {

            iniParam[0] = HW2F_A_INI;
            iniParam[1] = HW2F_S_INI;
            iniParam[2] = HW2F_DA_INI;
            iniParam[3] = HW2F_DS_INI;
            iniParam[4] = HW2F_DC_INI;

            for (i = 0; i < 5; i++)
            {
                if (i<size) 
                {
                   iniParam[i] = initVect[i];
                   CST[i] = flagEstim[i];
                }
                else
                   CST[i] = 0;
            }

 
            if (mod->GetName() != ARM_HWTWOFACTORMODEL)
            {
                tmpmod = (ARM_Model *) new ARM_HWTwoFactorModel(iniParam[0], 
                 iniParam[1], iniParam[2], iniParam[3], iniParam[4], zc);
 
                if (mod)
                {
                   delete mod;
                   mod = NULL;
                }
 
                mod = tmpmod;
 
                LOCAL_PERSISTENT_OBJECTS->SetPersistent(mod, modId);
            }
            else
            {
                if (initVect != NULL)
                {
                    vParam = new ARM_Vector(5, iniParam);
                    mod->SetParameters(vParam);
                    mod->SetZeroCurve(zc);
                    // delete vParam; Not Cloned !
                }
                else
                    mod->SetZeroCurve(zc);
            }
        }

        else if (strcmp(mod_name, "GYLS")==0)
        {

            iniParam[0] = GYLS_A_INI;
            iniParam[1] = GYLS_S_INI;
            iniParam[2] = GYLS_X_INI;
            iniParam[3] = GYLS_XV_INI;
            iniParam[4] = GYLS_XC_INI;
            iniParam[5] = GYLS_LR_INI;

            for (i=0;i<6;i++)
            {
                if (i<size) 
                {
                    iniParam[i] = initVect[i];
                    CST[i] = flagEstim[i];
                }
                else
                    CST[i] = 0;
            }

 
            if (mod->GetName() != ARM_GYCLSMODEL)
            {
               tmpmod = (ARM_Model *) new ARM_GYCLSModel(iniParam[0], 
                         iniParam[1],
                         zc, 
                         iniParam[2], 
                         iniParam[3], iniParam[4], iniParam[5]);
 
               if (mod)
               {
                  delete mod;
                  mod = NULL;
               }
 
               mod = tmpmod;
 
               LOCAL_PERSISTENT_OBJECTS->SetPersistent(mod, modId);
            }
            else
            {
               if (initVect != NULL)
               {
                  vParam = new ARM_Vector(6, iniParam);
                  mod->SetParameters(vParam);
                  mod->SetZeroCurve(zc);
                  // delete vParam; Not Cloned !
                }
                else
                   mod->SetZeroCurve(zc);
            }
        }
        else if (strcmp(mod_name, "ZCVS")==0)
        {
            iniParam[0] = A_INI;
            iniParam[1] = R_INI;
            iniParam[2] = S_INI;
            iniParam[3] = G_INI;

            for (i=0;i<4;i++)
            {
                if (i<size) 
                {
                   iniParam[i] = initVect[i];
                   CST[i] = flagEstim[i];
                }
                else
                   CST[i] = 0;
            }

 
            if (mod->GetName() != ARM_YCMODEL || 
                mod->GetZeroCurve()->GetName() != ARM_ZERO_VASICEK)
            {
               ARM_ZeroVasicek* zc = new ARM_ZeroVasicek((ARM_Date) settl, 
                                                          iniParam);
 
               tmpmod = (ARM_Model *) new ARM_YCModel(zc);

               if (mod)
               {
                  delete mod;
                  mod = NULL;
               }
 
               mod = tmpmod;
 
               LOCAL_PERSISTENT_OBJECTS->SetPersistent(mod, modId);
            }
            else
            {
               if (initVect != NULL)
               {
                  vParam = new ARM_Vector(4, iniParam);
                  mod->SetParameters(vParam);
                  // delete vParam; Not Cloned !
               }
            }
        }
        else if (strcmp(mod_name, "ZCSP")==0)
        {
            iniParam[0] = SP1_INI;
            iniParam[1] = SP2_INI;
            iniParam[2] = SP3_INI;
            iniParam[3] = SP4_INI;
            iniParam[4] = SP5_INI;
            iniParam[5] = SP6_INI;
 
            for (i = 0; i < 6; i++)
            {
                if (i<size) 
                {
                   iniParam[i] = initVect[i];
                   CST[i] = flagEstim[i];
                }
                else
                   CST[i] = 0;
            }

            if (mod->GetName() != ARM_YCMODEL || 
                mod->GetZeroCurve()->GetName() != ARM_ZERO_SPLINES)
            {
 
                ARM_ZeroSplines* zc = new ARM_ZeroSplines((ARM_Date) settl, 6, 
                                                iniParam);

                tmpmod = (ARM_Model *) new ARM_YCModel(zc);

                if (mod)
                {
                   delete mod;
                   mod = NULL;
                }
 
                mod = tmpmod;
 
                LOCAL_PERSISTENT_OBJECTS->SetPersistent(mod, modId);
            }
            else
            {
                if (initVect != NULL)
                {
                    vParam = new ARM_Vector(6, iniParam);
                    mod->SetParameters(vParam);
                    // delete vParam; Not Cloned !
                }
            }
        }
        else
        {
            throw Exception(__LINE__, __FILE__, ERR_NUMERICAL_CALL,
                           "Problem in Model Name");
        }

        int sz = pf->GetSize();
 
        for (i = 0; i < sz; i++)
        {
           curAsset = pf->GetAsset(i);
 
           curAsset->SetModel(mod);
        }
 
        ARM_Date date = (ARM_Date) settl;
 
        Cont = pf->PfToCont();


        double OP[19];
 
    
		MEMSET(OP, 0, sizeof(OP));

        ARM_Matrix* Data = pf->PfToMarketData();
 
        retCode = mod->FitToMarketData(Cont, Data, CST, 0, 0, OP);
 
        if (Cont)
        {
           delete Cont;
           Cont = NULL;
        }

        if (Data)
        {
           if (Data)
              delete Data;

           Data = NULL;
        }
 
        if ( retCode >= 3 )
        {
           throw Exception(__LINE__, __FILE__, ERR_NUMERICAL_CALL,
                           "Problem in Fitting !");
        }
        else
        {
            if (retCode == 0) 
               strcpy(msg, "Optimum Reached !");

            if (retCode == 1) 
               strcpy(msg, "Maximum Iterations Reached !");

            if (retCode == 2) 
               strcpy(msg, 
         "Optimum conditions not satisfied but no better solution locally !");
        }

        /*---- The result ----*/
	
	    /*
        FILE* MOD = ARM_GetOutputHomeFile("MOD", "a+");
 
 
        ARM_Vector* Parametre;
 
        Parametre = mod->GetParameters();
 

        fprintf(MOD, "Model de Type %s \n", mod_name);
 
        for (i=0; i<mod->GetNbParams(); i++)
        {
            fprintf(MOD, "Param[%d] = %10.10f\n", i, Parametre->Elt(i));
        }
 
        fclose(MOD);
		*/

        CREATE_GLOBAL_OBJECT();
 
        //modId = LOCAL_PERSISTENT_OBJECTS->SetPersistent(mod);
		retour= modId;
    }
 
    catch(Exception& x)
    {
        x.DebugPrint();
 
		return ARM_KO;
		//ARM_RESULT();
    }
 
    //ARM_END(pfTrace);
 
    retour = modId;

    strcpy(rc.msg, msg);

    return(ARM_OK);
}



long CC_LOCAL_SetPFModEstimated(int modId,
								char* mod_name,
								int pfId,
								char* settl,
								int zcId,
								int size,
								double* initVect,
								int* flagEstim,
								int step,
								char* horizon,
								int NAG_ALGO,
								long & retour)
{
    static char* __ARM_FUNC__ = "CC_LOCAL_SetPFModEstimated";
 
 
    //ARM_BEGIN(pfTrace);
 
    if ( NAG_ALGO == 0 )
    {
       if (CC_LOCAL_SetPFModEstimatedAlgo1(modId, mod_name, pfId, settl, zcId,
                                      size, initVect, flagEstim,
                                      step, horizon,retour)== ARM_OK)
        return ARM_OK;
	   else
		return ARM_KO;
    }


    try
    {
       if (CC_LOCAL_PFModEstimate(mod_name, pfId, settl, zcId,
                              size, initVect, flagEstim,
                              step, horizon, NAG_ALGO, modId,retour) == ARM_OK)
		  return ARM_OK;
	   else
		  return ARM_KO;
 
    }

    catch(Exception& x)
    {
        x.DebugPrint();
 		return ARM_KO;
    }
 
	return ARM_OK;
} 







// ************************************************************************************
// ARM_PFMODFIT                                                             ***********
// ************************************************************************************

long ARMLOCAL_PFMODFIT (const CCString& modName,
						long pfId,
						double settlement,
						long zcId,
						VECTOR<double>& vList,
						VECTOR<double>& fList,
						long step,
						double dhorizon,
						long nag_algo,
						ARM_result& result,
						long objId)
{
        long rc;
 
        int idMod, idZc, idPF;
        int i, size;
        char mname[50];
        char settl[50];
        char horizon[50];
        int NAG_ALGO;
	    double* VVal = NULL;
	    int*   VFlag = NULL;

	// ************************************************************************
	// Conversion des parametres 
	// ************************************************************************

        long tmp_size = vList.size();

        if ( tmp_size == 0 )
        {
           if ( modName == "GYCM" )
           {
              
              tmp_size = 2;

              for (int j = 0; j < tmp_size; j++)
              {
                  vList.push_back(0.0);
                  fList.push_back(0.0);
              }
           }
           else
           {
              vList.push_back(0.0);
              fList.push_back(0.0);
           }
        }

        VECTOR<long> fListTmp;
        for (i = 0; i < fList.size (); i++)
        {
            fListTmp.push_back ((long)fList[i]);
        }



    try
    {
        /*--- parameters checking ---*/
        if ( vList.size () != fList.size ()) 
        {
            result.setMsg ("ARM_ERR: values and flags array must have same size");
            return ARM_KO;
        }

        if ( objId != -1 )
        {
			// SET

		        strcpy(mname, modName);
		        idPF  = pfId;
		        
                Local_XLDATE2ARMDATE(settlement, settl);
		    
                idZc = zcId;
		        size = tmp_size;

				Local_XLDATE2ARMDATE(dhorizon, horizon);
		        
                NAG_ALGO = nag_algo;
				idMod    = objId;

		        if ( size > 0 )
				{
					VVal = (double *) malloc(sizeof(double)*size);
					VFlag = (int *) malloc(sizeof(int)*size);
					
					
		            for ( i = 0; i < size ; i++)
				    {
						VVal[i] = vList[i];
		                VFlag[i] = fList[i];
				    }
				}
 
				if ( CC_LOCAL_SetPFModEstimated(idMod, mname, idPF, settl, idZc, size,
                                                VVal, VFlag, step, horizon, NAG_ALGO,rc) == ARM_KO)
				{
				   if (VVal)
					  free(VVal);
   
				   if (VFlag)
					  free(VFlag);

					return ARM_KO;
				}	
 
				if (VVal)
				   free(VVal);
   
				if (VFlag)
				   free(VFlag);

				result.setLong(rc);

				return(ARM_OK);            
		}
        else
        {	
		     strcpy(mname, modName);
		     idPF  = pfId;
		        
             Local_XLDATE2ARMDATE(settlement,settl);
		        
             idZc = zcId;
		        
             size = tmp_size;
				
             Local_XLDATE2ARMDATE(dhorizon,horizon);
		        
             NAG_ALGO   = nag_algo;

		     if ( size > 0 )
			 {
				VVal = (double *) malloc(sizeof(double)*size);
				VFlag = (int *) malloc(sizeof(int)*size);

				for ( i = 0; i < size ; i++)
				{
	                VVal[i] = vList[i];
	                VFlag[i] = fList[i];
			    }
			}
 
			if ( CC_LOCAL_PFModEstimate(mname, idPF, settl, idZc, size, VVal, VFlag,
                                        step, horizon, NAG_ALGO, -1,rc) != ARM_OK)

			{	

				if (VVal)
				   free(VVal);
   
				if (VFlag)
				   free(VFlag);

				return(ARM_KO);
			}


				if (VVal)
				   free(VVal);
   
				if (VFlag)
				   free(VFlag);

				result.setLong(rc);
                
				return(ARM_OK);
        }    
     }

     catch(Exception& x)
     {
        x.DebugPrint();

	/*	if (newMod)
			delete newMod;
		newMod = NULL;
	*/
	//	ARM_RESULT();
		return ARM_KO;
     }

     return ARM_KO;
}



/////////////MRS Calibration PF//////////////////////////////////////////////////
///////////////Portfolio for the calibration of the mean reversion //////////////
long ARMLOCAL_MRSCalibrationPF(long zcId,
							long volId,
							long secId,
							long portfolioId,
							int freq,
							bool ATMStrikes,
							ARM_result& result,
							long objId = -1)
{
	long pfId;

	ARM_Portfolio* newPortfolio = NULL;
	ARM_Portfolio* oldPortfolio = NULL;
	ARM_Portfolio* firstPortfolio = NULL;

	ARM_ZeroCurve* zc = NULL;
	ARM_BSModel * bsModel = NULL;
	ARM_Swaption* sec = NULL;
	
	CCString msg (" ");

	try
	{
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);
		firstPortfolio = (ARM_Portfolio*) LOCAL_PERSISTENT_OBJECTS->GetObject(portfolioId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(firstPortfolio, ARM_PORTFOLIO) == 0)
		{
			result.setMsg ("ARM_ERR:Portfolio is not of a good type");
			return ARM_KO;
		}

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR:Zc Curve is not of a good type");
			return ARM_KO;
		}

		bsModel = (ARM_BSModel *) LOCAL_PERSISTENT_OBJECTS->GetObject(volId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(bsModel, ARM_BSMODEL) == 0)
		{
			result.setMsg ("ARM_ERR: VolCurve is not of a good type");
			return ARM_KO;
		}

		sec = (ARM_Swaption *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: security is not of a good type");
			return ARM_KO;
		}

		newPortfolio = sec->MRSCalibrationPF(zc,bsModel,firstPortfolio,freq,ATMStrikes);

		if (newPortfolio == NULL)
		{
			result.setMsg ("ARM_ERR: Portfolio is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			pfId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newPortfolio);

			if (pfId == RET_KO)
			{				
				if (newPortfolio)
				{
					newPortfolio->FreePortfolioAndAssets();
					delete newPortfolio;
					newPortfolio = NULL;
				}

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(pfId);

			return ARM_OK;
		}
		else
		{
			oldPortfolio = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldPortfolio, ARM_PORTFOLIO) == 1)
			{
				if (oldPortfolio)
				{
					oldPortfolio->FreePortfolioAndAssets();
					delete oldPortfolio;
					oldPortfolio = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newPortfolio, objId);

				return ARM_OK;
			}
			else
			{
				if (newPortfolio)
				{
					newPortfolio->FreePortfolioAndAssets();
					delete newPortfolio;
					newPortfolio = NULL;
				}

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newPortfolio)
		{
			newPortfolio->FreePortfolioAndAssets();
			delete newPortfolio;
			newPortfolio = NULL;
		}
		
		ARM_RESULT();
	}
}



// YK CalibrationPF//////////////////////////////////////////////////////////////

long ARMLOCAL_CalibrationPF(long zcId,
							long volId,
							long secId,
							long modeType,
							long portfolioType,
							double shift,
							ARM_result& result,
							long objId)

{
	long pfId;

	ARM_Portfolio* newPortfolio = NULL;
	ARM_Portfolio* oldPortfolio = NULL;

	ARM_ZeroCurve* zc = NULL;
	ARM_VolCurve* vol = NULL;
	ARM_Security* sec = NULL;
	
	CCString msg (" ");

	try
	{
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(zcId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR:Zc Curve is not of a good type");
			return ARM_KO;
		}

		vol = (ARM_VolCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(volId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(vol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: VolCurve is not of a good type");
			return ARM_KO;
		}

		sec = (ARM_Security *) LOCAL_PERSISTENT_OBJECTS->GetObject(secId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(sec, ARM_SECURITY) == 0)
		{
			result.setMsg ("ARM_ERR: security is not of a good type");
			return ARM_KO;
		}
		if (2 == portfolioType);
			//newPortfolio = sec->MRSCalibrationPF(zc,vol,modeType);
		else if (1 == portfolioType)
			newPortfolio = sec->CalibrationPF(zc,vol,4);
		else
			newPortfolio = sec->CalibrationPF(zc,vol,modeType);

		if (newPortfolio == NULL)
		{
			result.setMsg ("ARM_ERR: Portfolio is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
		
			pfId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newPortfolio);

			if (pfId == RET_KO)
			{				
				if (newPortfolio)
				{
					newPortfolio->FreePortfolioAndAssets();
					delete newPortfolio;
					newPortfolio = NULL;
				}

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(pfId);

			return ARM_OK;
		}
		else
		{
			oldPortfolio = (ARM_Portfolio *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(oldPortfolio, ARM_PORTFOLIO) == 1)
			{
				if (oldPortfolio)
				{
					oldPortfolio->FreePortfolioAndAssets();
					delete oldPortfolio;
					oldPortfolio = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newPortfolio, objId);

				return ARM_OK;
			}
			else
			{
				if (newPortfolio)
				{
					newPortfolio->FreePortfolioAndAssets();
					delete newPortfolio;
					newPortfolio = NULL;
				}

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (newPortfolio)
		{
			newPortfolio->FreePortfolioAndAssets();
			delete newPortfolio;
			newPortfolio = NULL;
		}
		
		ARM_RESULT();
	}
}

extern string SecurityTypeToClass(long portfolioId, long assetIdIndex)
{
    /// Find the model name for interface
    ARM_Portfolio* portfolio = dynamic_cast<ARM_Portfolio *>(LOCAL_PERSISTENT_OBJECTS->GetObject(portfolioId));
	if (!portfolio)
		return LOCAL_ANY_CLASS;

    if(assetIdIndex > portfolio->GetSize()-1)
		return LOCAL_ANY_CLASS;

    switch(portfolio->GetAsset(assetIdIndex)->GetName())
    {
        case ARM_CAPFLOOR :
		case ARM_BARRIER :
            return LOCAL_CAPFLOOR_CLASS;

        case ARM_SWAPTION :
            return LOCAL_SWAPTION_CLASS;

        case ARM_SWAP :
            return LOCAL_SWAP_CLASS;

        case ARM_SWAPLEG :
        case ARM_FIXEDLEG :
            return LOCAL_SWAPLEG_CLASS;

        case ARM_CORRIDORLEG :
            return LOCAL_CORRIDORLEG_CLASS;

        case ARM_DIGITAL :
            return LOCAL_DIGITAL_CLASS;

        case ARM_CMSLEG :
            return LOCAL_CMSLEG_CLASS;
			
        case ARM_SPREADOPTION :
        case ARM_CORRIDORDBLCONDITION :
            return LOCAL_SPREAD_OPTION_CLASS;

        case ARM_OPTION :
		case ARM_STRIPOPTION :
		case ARM_STRIPDIGITALOPTION :
            return LOCAL_OPTION_CLASS;

        case ARM_PORTFOLIO :
            return LOCAL_PF_CLASS;

		case ARM_POWERREVERSE :
			return LOCAL_POWER_REVERSE_CLASS;

		case ARM_GLOBALCAP :
			return LOCAL_GLOBALCAP_CLASS;

		case ARM_FXSPREADSTRIPOPTION :
			return LOCAL_FXSTRIP_CLASS;

		default :
        return LOCAL_ANY_CLASS;
    }

    return LOCAL_ANY_CLASS;
 }

ARM_CLASS_NAME  GetClassNameFromSecurity(ARM_Security* security)
{
    switch(security->GetName())
        {
            case ARM_SWAPLEG :
            case ARM_FIXEDLEG :
                return ARM_SWAPLEG;

            case ARM_CAPFLOOR :
            case ARM_SWAPTION :
            case ARM_SWAP :
            case ARM_CORRIDORLEG :
            case ARM_DIGITAL :
            case ARM_CMSLEG :
            case ARM_PORTFOLIO :
			case ARM_POWERREVERSE :
			case ARM_GLOBALCAP :
			case ARM_FXSPREADSTRIPOPTION :
            case ARM_CORRIDORDBLCONDITION :
                return security->GetName();

            case ARM_OPTION :
            case ARM_STRIPOPTION :
            case ARM_STRIPDIGITALOPTION :
                return ARM_OPTION;

            default :
				return ARM_SECURITY;
        }
}

extern long ARMLOCAL_GetAssetFromPF(long portfolioId, 
                                    long assetIdIndex,
                                    ARM_result& result,
                                    long        objId )
{
    /// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;
    /// used in the MACRO ARM_RESULT
	CCString msg ("");

    ARM_Portfolio* portfolio = NULL;
    ARM_Object* object=NULL;
	ARM_CLASS_NAME className;

	try
	{
		if( !GetObjectFromId( &portfolio, portfolioId, ARM_PORTFOLIO ) )
		{
			result.setMsg ("ARM_ERR: Portfolio is not of a good type");
			return ARM_KO;
		};
        /// test
        if(assetIdIndex > portfolio->GetSize()-1)
        {
			result.setMsg ("ARM_ERR: Index is out of range, please advise!");
			return ARM_KO;
		};

        ARM_Security* security = (ARM_Security*)(portfolio->GetAsset(assetIdIndex))->Clone();
        object = (ARM_Object*)security;
        className = GetClassNameFromSecurity(security);        

		/// Assign the object in ARM cache
		if( !assignObject( object, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

}

extern long ARMLOCAL_GetSizeOfPF(long portfolioId, 
                                    long& sizeofPf,
                                    ARM_result& result)
{
    /// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;
    /// used in the MACRO ARM_RESULT
	CCString msg ("");

    ARM_Portfolio* portfolio = NULL;

	try
	{
		if( !GetObjectFromId( &portfolio, portfolioId, ARM_PORTFOLIO ) )
		{
			result.setMsg ("ARM_ERR: Portfolio is not of a good type");
			return ARM_KO;
		};

        sizeofPf = portfolio->GetSize(); 
        
        return ARM_OK; 
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

}


extern long ARMLOCAL_ComputeAll(long pfId, long modelId,ARM_result& result)
{
	ARM_Date* settlement=NULL;
	ARM_Portfolio* pf = NULL;
	ARM_Model* model= NULL;

	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	try
	{
		if( !GetObjectFromId( &pf, pfId, ARM_PORTFOLIO) )
		{
			result.setMsg ("ARM_ERR: portfolio is not of a good type");
			return ARM_KO;
		};

		if( !GetObjectFromId( &model, modelId, ARM_MODEL) )
		{
			result.setMsg ("ARM_ERR: model is not of a good type");
			return ARM_KO;
		};

		ARM_Vector* cptRes = pf->ComputeAll(model, settlement);

		int nb =  cptRes->GetSize();

		for ( int i=0; i<nb; i++)
		{
			result.setArray(cptRes->Elt(i),i);	
		}
		
		result.setDouble(nb);

		if(cptRes)
			delete cptRes;
		cptRes = NULL;

		return ARM_OK;

	}
	
	catch(Exception& x)
	{	
		x.DebugPrint();
		
		ARM_RESULT();
	}
}


 
