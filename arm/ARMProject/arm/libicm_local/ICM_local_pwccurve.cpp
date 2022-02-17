

#include "firstToBeIncluded.h"
#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_persistent.h>

#include <ARM\libarm\ARM_result.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <ARM\local_xlarm\XL_local_xlarm_common.h>
#include <ARM\libarm_local\arm_local_init.h>

#include <ARM\libarm_frometk\arm_local_etoolkit_for_ICM.h>
#include <ARM\libarm_frometk\arm_local_parsexml_for_ICM.h>
#include <ARM\libarm_frometk\arm_local_parsexml_nt_curves.h>

#include <ccy\currency.h>
// #include <ICMKernel\crv\icm_defcurve_index_pwc.h>
#include <ICMKernel\crv\icm_constant_piecewise.h>
#include <ICMKernel\crv\icm_linear_piecewise.h>
#include <ICMKernel\crv\icm_interpoldefcrv.h>
#include <ICMKernel\crv\icm_volinterpol.h>
#include <ICMKernel\crv\icm_implied_loss_tree.h>
#include <ICMKernel\crv\icm_defaultcurveSimple.h>
#include <ICMKernel\crv\icm_fixing_curve.h>
#include <ICMKernel\inst\icm_credit_index.h>

#include <ARM\libarm_frometk\arm_local_paesexml_calypso.h> 

#include "icm_local_pwccurve.h" 

#include <ios>

#include <streambuf>

 
long ICMLOCAL_SurvivalProba (long idCurve, 
							 double matu,
							 ARM_result& result)
{
	double dResult;
	ARM_ZeroCurve* zc = NULL;
	ICM_DefaultCurve* DC = NULL;

	ARM_Date startdate;

	char* sDate = new char[11];
	Local_XLDATE2ARMDATE(matu,sDate);
	ARM_Date maturity(sDate);

	double fracyear = 0;
	ARM_CLASS_NAME objname = ARM_ZERO_CURVE;


	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		zc = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ARM_ZERO_CURVE) == 0)
		{
			objname = ICM_DEFAULTCURVE;

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ICM_DEFAULTCURVE) == 0)
			{
			result.setMsg ("ARM_ERR: Default Curve is not of a good type");
			return ARM_KO;
			}

			DC = (ICM_DefaultCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);
		}

		startdate = zc->GetAsOfDate();
		fracyear = (maturity.GetJulian() - startdate.GetJulian())/360.;

		if (objname == ARM_ZERO_CURVE)
			dResult = zc->DiscountPrice(fracyear);
		else 
			dResult = DC->SurvivalProba(maturity);

		if (sDate)
			delete[] sDate;
		sDate = NULL;

		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		if (sDate)
			delete[] sDate;
		sDate = NULL;

		x.DebugPrint();

		ARM_RESULT();
	}
}


long ICMLOCAL_DefaultProba (long idCurve, 
							 double matu,
							 ARM_result& result)
{
	double dResult;
	ICM_DefaultCurve* DC = NULL;

	ARM_Date startdate;

	char* sDate = new char[11];
	Local_XLDATE2ARMDATE(matu,sDate);
	ARM_Date maturity(sDate);

	double fracyear = 0;
	ARM_CLASS_NAME objname = ARM_ZERO_CURVE;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		DC = (ICM_DefaultCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(DC, ICM_DEFAULTCURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Default Curve is not of a good type");
			return ARM_KO;
		}

		dResult = DC->DefaultProba(maturity);

		if (sDate)
			delete[] sDate;
		sDate = NULL;

		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{

		if (sDate)
			delete[] sDate;
		sDate = NULL;

		x.DebugPrint();

		ARM_RESULT();
	}
}



long ICMLOCAL_DefaultIntensity (long idCurve, 
								double matu, 
								long meth,
							    ARM_result& result)
{
	double dResult;
	ICM_DefaultCurve* zc;

	ARM_Date startdate;

	char* sDate = new char[11];
	Local_XLDATE2ARMDATE(matu,sDate);
	ARM_Date maturity(sDate);

	double fracyear = 0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		zc = (ICM_DefaultCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ICM_DEFAULTCURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Zc Curve is not of a good type");
			return ARM_KO;
		}

		dResult = zc->GetPWCIntensity(maturity);

		if (sDate)
			delete[] sDate;
		sDate = NULL;

		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{

		if (sDate)
			delete[] sDate;
		sDate = NULL;

		x.DebugPrint();

		ARM_RESULT();
	}
}


double ICMLOCAL_CptImplicitSpreadInterpol (long defcurveId, 
										   const std::string& plot) 
{

	ICM_DefaultCurve* defCurve; 
	LocalPersistent::get().convert(defcurveId,defCurve); 
	if (!defCurve) ICMTHROW(ERR_INVALID_ARGUMENT,"Not a default curve"); 
	return defCurve->ImpliedSpreadInterpol(plot); 
}	

double ICMLOCAL_CptImplicitSpreadInterpol (long defcurveId, 
										   const ARM_Date& date) 
{

	ICM_DefaultCurve* defCurve; 
	LocalPersistent::get().convert(defcurveId,defCurve); 
	if (!defCurve) ICMTHROW(ERR_INVALID_ARGUMENT,"Not a default curve"); 
	return defCurve->ImpliedSpreadInterpol(date); 
}
long  ICMLOCAL_createFlatCurve (long defcurveId, const std::string& plot,long prevId) 
{
	ICM_DefaultCurve* defCurve; 
	LocalPersistent::get().convert(defcurveId,defCurve); 
	if (!defCurve) ICMTHROW(ERR_INVALID_ARGUMENT,"Not a default curve"); 
	ICM_DefaultCurve * item = defCurve->createFlatCurve(plot); 
	return LocalPersistent::get().adopt(item,prevId); 
}	
long ICMLOCAL_createFlatCurve (long defcurveId, const ARM_Date& date,long prevId) 
{
	ICM_DefaultCurve* defCurve; 
	LocalPersistent::get().convert(defcurveId,defCurve); 
	if (!defCurve) ICMTHROW(ERR_INVALID_ARGUMENT,"Not a default curve"); 
	ICM_DefaultCurve * item = defCurve->createFlatCurve(date); 
	return LocalPersistent::get().adopt(item,prevId); 
}

long ICMLOCAL_createDefCurveFromBase(long defcurveIdCDS,long defcurveIdIndex , const ARM_Vector& vBase,long prevId)
{
	ICM_DefaultCurve* defCurveCDS;
	ICM_DefaultCurve* defCurveIndex; 
	LocalPersistent::get().convert(defcurveIdIndex,defCurveIndex); 
	LocalPersistent::get().convert(defcurveIdCDS,defCurveCDS); 

	if (!defcurveIdCDS) ICMTHROW(ERR_INVALID_ARGUMENT,"Not a default curve"); 
	if (!defCurveIndex) ICMTHROW(ERR_INVALID_ARGUMENT,"Not a default curve Index");
	ICM_DefaultCurve * item = defCurveCDS->createDefCurveFromBase(defCurveIndex, vBase); 
	return LocalPersistent::get().adopt(item,prevId); 

}
long ICMLOCAL_GetDPFromSummit (const ARM_Date&AsOf,
							   const CCString& issuer,
							   const CCString& CurveName,
							   const long& PWCcurveId,	
							   const CCString& label,
							   ARM_result& result,
							   long objId)
{

	if (GetDataRetrieverVersion() == FFRETRIEVER)
	{
		result.setMsg ("ARM_ERR: Function not implemented without ETK");
		return ARM_KO;
	}

	long curveId = 0;

	ICM_DefaultCurve* pwcdefault = NULL;
	ICM_DefaultCurve* prevpwcdefault = NULL;

	ARM_ZeroCurve* pwcshort = NULL;
	ARM_Currency* ccy = NULL;

	CCString msg (" ");

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}


	try
	{
		CCString xmlResponse_Credit;
		CCString messageList;

		long out = etoolkit_GetDefProbCurve(issuer,AsOf, CurveName, xmlResponse_Credit,messageList);

		if (out == ARM_KO)  
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"etoolkit_GetDefProbCurve fail");

		pwcshort = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(PWCcurveId);


		pwcdefault = ARMLOCAL_ParseDefProbCurve(xmlResponse_Credit,AsOf,pwcshort,(const char*)label);

		if (pwcdefault == NULL)
		{
			result.setMsg("Object is Null");
			return ARM_KO;
		}
	}
	catch (Exception& x)
	{
		x.DebugPrint();

		if (pwcdefault)
			delete pwcdefault;
		pwcdefault = NULL;

		ARM_RESULT();
	}

	if(objId == -1)
	{
		CREATE_GLOBAL_OBJECT();

		curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)pwcdefault);

		if (curveId == RET_KO)
		{
			if (pwcdefault)
				delete pwcdefault;
			pwcdefault = NULL;

			result.setMsg ("ARM_ERR: Pb with inserting object");				
			return ARM_KO;
		}

		result.setLong(curveId);

		return ARM_OK;
	}
	else
	{
		prevpwcdefault = (ICM_DefaultCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(prevpwcdefault, ICM_DEFAULTCURVE) == 1)
		{
			if (prevpwcdefault)
			{
				delete prevpwcdefault;
				prevpwcdefault = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)pwcdefault, objId);

			return ARM_OK;
		}
		else
		{
			result.setMsg ("ARM_ERR: previous object is not of a good type");
			return ARM_KO;
		}
	}

}

//---------------------------------------------------------------------------
// Recuperation DefProb From calypso
//--------------------------------------------------------------------------

long
ICMLOCAL_GetDPFromCalypso (const ARM_Date& AsOf,
							const std::string& issuer,
							const std::string& Seniority,
							const std::string& Currency,
							const std::string& PricingEnv,
							const std::string& ForceCurveName,
							const long& PWCcurveId,	
							const std::string& label,
							const std::string& xmlFileName,
							long prevId)
{		
	ARM_ZeroCurve* pwcshort = NULL;
	
	std::string xmlContent ;
	ARM_CalypsoToolkit::GetProbabilityCurve(issuer,Currency,Seniority,ForceCurveName,PricingEnv,AsOf,xmlFileName,xmlContent); 
	
	LocalPersistent::get().convert(PWCcurveId,pwcshort); 

	// creates the AMR object & acquire ownership  
	std::auto_ptr<ICM_DefaultCurve> pwcdefault (  ARMLOCAL_ParseDefProbCurveCalypso(xmlContent ,AsOf,pwcshort,issuer,Seniority,Currency,ForceCurveName) ) ;

	// overwrites the label ... 
	if ( (label!="") && (label!="NONE")) pwcdefault->SetLabel(label); 

	long retId = LocalPersistent::get().adopt(pwcdefault.get(),prevId); 
	pwcdefault.release(); 
	return retId; 
}

long ICMLOCAL_DPMktDataFromSummit (double AsOf,
										  const CCString& issuer,
										  const CCString& CurveName,
										  VECTOR<CCString>&	Matu,
										  VECTOR<double>& Spread,
										  VECTOR<double>& Recovery,
										  CCString& Currency,
										  CCString& IndexName,
										  long& AIMMADJ, 
										  ARM_result& result)
{

	if (GetDataRetrieverVersion() == FFRETRIEVER)
	{
		result.setMsg ("ARM_ERR: Function not implemented without ETK");
		return ARM_KO;
	}

	long retour = 0;

	CCString msg (" ");

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* sDate = new char[11];
	Local_XLDATE2ARMDATE(AsOf,sDate);
	ARM_Date CurveDate(sDate);

	try
	{
		CCString xmlResponse;
		CCString messageList;

		char sDate[11];
		Local_XLDATE2ARMDATE(AsOf,sDate);
		ARM_Date myDate(sDate);

		long out = etoolkit_GetDefProbCurve(issuer,CurveDate, CurveName, xmlResponse,messageList);

		if (out == ARM_KO)  
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"etoolkit_GetDefProbCurve fail");

		retour = ARMLOCAL_ParseDefProbCurveMktData(xmlResponse,Matu,Spread,Recovery,Currency,IndexName,AIMMADJ);

		return ARM_OK;

	}
	catch (Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

	return ARM_KO;
}

void  ICMLOCAL_DPMktDataFromCalypso(const ARM_Date& AsOf,
									const std::string & pricingEnv,
									const std::string & issuer,
									const std::string & seniority,
									const std::string & ccy,
									const std::string & forceCurveName,
									const std::string & xmlFileName,
									std::vector<std::string>& Matu,
									std::vector<double>& Spread,
									double & Recovery,
									std::string & Currency,
									std::string & IndexName,
									qCDS_ADJ & AIMMADJ )  
{
	

	std::string xmlOutput ;
	ARM_CalypsoToolkit::GetProbabilityCurve(issuer,ccy,seniority,forceCurveName,pricingEnv,AsOf,xmlFileName,xmlOutput); 
	ARMLOCAL_ParseDefProbCurveMktDataCalypso(xmlOutput,Matu,Spread,Recovery,Currency,IndexName,AIMMADJ);
}

long ICMLOCAL_CorrFromSummit(double C_AsOfDate,
							const CCString& Name,
							const CCString& Name2,
							const CCString& CurveId,
							ARM_result& result)
{

	if (GetDataRetrieverVersion() == FFRETRIEVER)
	{
		result.setMsg ("ARM_ERR: Function not implemented without ETK");
		return ARM_KO;
	}

	long retour = 0;

	ARM_Currency* ccy = NULL;

	CCString msg (" ");

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	char* sDate = new char[11];
	Local_XLDATE2ARMDATE(C_AsOfDate,sDate);
	ARM_Date CurveDate(sDate);

	try
	{
		CCString xmlResponse;
		CCString messageList;

		long out = etoolkit_GetCorrFromSummit(CurveDate,Name,Name2,CurveId, xmlResponse,messageList);

		result.setDouble(0.); 

		if (sDate)
			delete[] sDate;
		sDate = NULL;

		return ARM_OK;

	}
	catch (Exception& x)
	{
		x.DebugPrint();

		if (sDate)
			delete[] sDate;
		sDate = NULL;

		ARM_RESULT();
	}

	return ARM_KO;
}


long ICMLOCAL_ConstantDefaultCurve (const  ARM_Date& AsOf, 
								const std::vector<std::string> & Tenors, 
							   VECTOR<double>& Rates,
							   double recovery,
							   long ircurveid,
							   int intRule,
							   int adjStartRule,
							   const std::string & ccy, 
							   const std::string & Label, 
							   qCDS_ADJ AdjCalType,
							   bool IsSummitCurve,
								qDEFCURVE_CALIB_ALGO calibrationAlgo,
							   long VolCurveId,
							   const std::string& calibrationData,
							   int Lag,
							   long paramId,
							   long prevId)
{
    ARM_Vector VRates; //= NULL;

	// non optional ZeroCurve 
	ARM_ZeroCurve* ZC =NULL;
	LocalPersistent::get().convert(ircurveid,ZC); 
	if (!ZC) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_ConstantDefaultCurve: ZeroCurve not provided."); 
	
	// optional VolCurve
	ARM_VolCurve* VolCurve = NULL;
	LocalPersistent::get().convert(VolCurveId,VolCurve); 

	// optional Parameter
	ICM_Parameters* params =NULL; 
	LocalPersistent::get().convert(paramId,params); 

	int i=0;
	int real_size = Tenors.size ();
	int real_size2 = real_size;

	
	vector<double> vRate; 

	int nbdefvalues = 0;
	int i2 = 0;

	

	for(i = 0; i < real_size; i++) { if (Rates[i] == -999.) nbdefvalues+=1; }
	real_size2 -=nbdefvalues;

	vRate.resize(real_size2);

		
	i2=0;
	for(i = 0; i < real_size; i++)
	{ if (Rates[i] != -999.) {vRate[i2] = Rates[i] / 10000.;i2++;} }

	VRates  = ARM_Vector(vRate);

	vector<string> terms;

	i2=0;
	for(int j = 0; j < real_size; j++) 
	{
		if (Rates[j] != -999.) {terms.push_back(Tenors[j]);}
	}


	long PayFreq = K_QUARTERLY;


	ICM_DefaultCurve* createdDC = new ICM_Constant_Piecewise(  AsOf,
											terms, 
											&VRates, 
											recovery,
											ZC,
											intRule,	// intRule
											adjStartRule,	// adjStartDate
											 AdjCalType, 
											ccy, 
											Label,
											IsSummitCurve,
											//2 NULL,
											VolCurve,
											PayFreq,
											calibrationAlgo,
											calibrationData,
											Lag,
											params?*params:ICM_Parameters());

	if (!createdDC) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Can't create ICM_Constant_Piecewise"); 

	long objId = LocalPersistent::get().adopt(createdDC,prevId); 

	return objId; 
}


long ICMLOCAL_LinearDefaultCurve (const ARM_Date& AsOf, 
								  const std::vector<std::string>&  Tenors, 
							   VECTOR<double>& Rates,
							   double recovery,
							   long ircurveid,
							   const CCString& ccy, 
							   const CCString& Label, 
							   qCDS_ADJ AdjCalType,
							   bool IsSummitCurve,
							   ARM_result& result, 
							   long objId)
{
	long curveId;

	ICM_DefaultCurve* createdDC = NULL;
	ICM_DefaultCurve* prevDC = NULL;
	ARM_ZeroCurve* ZC =NULL;
	// ARM_Currency* aCcy = NULL;
    ARM_Vector* VRates = NULL;

	int i=0;
	int real_size = Tenors.size ();

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	double * pdRate = NULL;

	CCString msg ("");

    try
	{
		ZC = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(ircurveid);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ZC, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: previous object IRcurve is not of a good type");
			return ARM_KO;
		}

		pdRate = new double[real_size];
		
		for(i = 0; i < real_size; i++)
			pdRate[i] = Rates[i] / 10000.;

		VRates  = new ARM_Vector(real_size,pdRate);

		vector<string> terms;
		for ( int i= 0; i< real_size ; i++)	
			terms.push_back(Tenors[i]);

		createdDC = new ICM_Linear_Piecewise(  AsOf,
												terms, 
												VRates, 
												recovery,
												ZC,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
												AdjCalType, 
												CCSTringToSTLString(ccy), 
												CCSTringToSTLString(Label),
												false /*IsSummitCurve*/);

		if (pdRate)
			delete[] pdRate;
		pdRate = NULL;

		if (VRates)
			delete VRates;
		VRates = NULL;


		if (createdDC == NULL)
		{
			result.setMsg ("ARM_ERR: Curve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdDC);

			if (curveId == RET_KO)
			{
				if (createdDC)
					delete createdDC;
				createdDC = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			prevDC = (ICM_DefaultCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevDC, ICM_LINEAR_PIECEWISE) == 1)
			{
				if (prevDC)
				{
					delete prevDC;
					prevDC = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdDC, objId);

				return ARM_OK;
			}
			else
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (pdRate)
			delete[] pdRate;
		pdRate = NULL;

		if (VRates)
			delete VRates;
		VRates = NULL;


		ARM_RESULT();
    }
}


long ICMLOCAL_FwdSpread (long DefCurveId, 
						 double matu1,
						 double matu2,
						 double fwdstart,
						 double fwdend,
						 long VolId, 
						 ARM_result& result)
{
	double dResult;
	ICM_DefaultCurve* DC = NULL;
	ARM_VolCurve* vol = NULL;

	char* sDate1 = new char[11];
	Local_XLDATE2ARMDATE(matu1,sDate1);
	ARM_Date date1(sDate1);

	char* sDate2 = new char[11];
	Local_XLDATE2ARMDATE(matu2,sDate2);
	ARM_Date date2(sDate2);

	char* sfwdstart = new char[11];
	ARM_Date dfwdstart;

	char* sfwdend = new char[11];
	ARM_Date dfwdend;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		if (fwdstart>0)
		{
		Local_XLDATE2ARMDATE(fwdstart,sfwdstart);
		dfwdstart = (ARM_Date) sfwdstart;
		}

		if (fwdend>0)
		{
		Local_XLDATE2ARMDATE(fwdend,sfwdend);
		dfwdend = (ARM_Date) sfwdend;
		}

		DC = (ICM_DefaultCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(DefCurveId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(DC, ICM_DEFAULTCURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Default Curve is not of a good type");
			return ARM_KO;
		}


		if (VolId != -1)
		{
		vol = (ARM_VolCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(VolId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(vol, ARM_VOL_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Volatility Curve is not of a good type");
			return ARM_KO;
		}
		}

		dResult = DC->FwdSpread(date1,date2);
		if (VolId != -1)
		{
			double yt1 = (dfwdstart-DC->GetAsOfDate())/365.;
			double yt2 = (dfwdend-DC->GetAsOfDate())/365.;
			//double volatility = vol->ComputeVolatility(yt1,0.)/100.;
			dResult+= DC->AjustConvexity(yt1,yt2,dResult,vol); 
		}

		if(sfwdstart) delete sfwdstart;
			sfwdstart = NULL;
		if(sfwdend) delete sfwdend;
			sfwdend = NULL;

		if (sDate1)
			delete[] sDate1;
		sDate1 = NULL;

		if (sDate2)
			delete[] sDate2;
		sDate2 = NULL;

		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{

		if(sfwdstart) delete sfwdstart;
			sfwdstart = NULL;
		if(sfwdend) delete sfwdend;
			sfwdend = NULL;
		if (sDate1)
			delete[] sDate1;
		sDate1 = NULL;

		if (sDate2)
			delete[] sDate2;
		sDate2 = NULL;

		x.DebugPrint();

		ARM_RESULT();
	}
}


long ICMLOCAL_GetZC_DP_FromSummit (const CCString& Issuer,
								   const CCString& currency,
								   const CCString& cvName,
								   double aSdate,
								   long ircurveid,
								   const CCString& Label,
								   ARM_result& result,
								   long objId)
{
	long curveId;

	ICM_DefaultCurve* createdZcLin = NULL;
	ICM_DefaultCurve* prevZcLin = NULL;
	ARM_ZeroCurve* zpy =NULL;

	char sDate[11];
	Local_XLDATE2ARMDATE(aSdate,sDate);

	CCString msg (" ");

	if (GetDataRetrieverVersion() == FFRETRIEVER)
	{
		result.setMsg ("ARM_ERR: Function not implemented without ETK");
		return ARM_KO;
	}

	try
	{
		CCString xmlResponse;
		CCString xmlMsgList;

		ARM_Date myDate(sDate);

		zpy = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(ircurveid);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zpy, ARM_ZERO_CURVE) != 1)
		{
			result.setMsg ("ARM_ERR: previous object IrCurve is not of a good type");
			return ARM_KO;
		}

		int retour = etoolkit_getXML_ZC_DP_FromSummit(Issuer,currency,cvName,myDate,xmlResponse,xmlMsgList);
		if (retour == ARM_KO)
		{
			result.setMsg("Unable to extract Zero Coupon for DefProb");
			return ARM_KO;
		}

		createdZcLin = ARMLOCAL_XML_DEFPROB_with_summit_stripper_Etk(xmlResponse,zpy,Label);

		if (createdZcLin == NULL)
		{
			result.setMsg("Object is Null");
			return ARM_KO;
		}
	}
	catch (Exception& x)
	{
		x.DebugPrint();

		if (createdZcLin)
			delete createdZcLin;
		createdZcLin = NULL;

		ARM_RESULT();
	}

	if(objId == -1)
	{
		CREATE_GLOBAL_OBJECT();

		curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZcLin);

		if (curveId == RET_KO)
		{
			if (createdZcLin)
				delete createdZcLin;
			createdZcLin = NULL;

			result.setMsg ("ARM_ERR: Pb with inserting object");				
			return ARM_KO;
		}

		result.setLong(curveId);

		return ARM_OK;
	}
	else
	{
		prevZcLin = (ICM_DefaultCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(prevZcLin, ICM_DEFAULTCURVE) == 1)
		{
			if (prevZcLin)
			{
				delete prevZcLin;
				prevZcLin = NULL;
			}

			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdZcLin, objId);

			return ARM_OK;
		}
		else
		{
			result.setMsg ("ARM_ERR: previous object is not of a good type");
			return ARM_KO;
		}
	}
}


long ICMLOCAL_GetZC_DP_FromCalypso(const ARM_Date& AsOfDate,
								   const std::string& pricingEnv,
								   const std::string& issuer,
								   const std::string& seniority,
								   const std::string& ccy,
								   const std::string& forceCurveName,
								   const std::string& xmlFileName,
								   long ircurveid,
								   const std::string& forceLabel, long objId) 
{
	// 
	ARM_ZeroCurve* irCurve =0 ; 
	LocalPersistent::get().convert(ircurveid,irCurve); 

	std::string xmlOutput ;
	ARM_CalypsoToolkit::GetProbabilityCurve(issuer,ccy,seniority,forceCurveName,pricingEnv,AsOfDate,xmlFileName,xmlOutput); 
	std::auto_ptr<ICM_DefaultCurve> item ( ARMLOCAL_XML_DEFPROB_with_calypso_stripper(xmlOutput,*irCurve,forceLabel) ); 
	
	long id= LocalPersistent::get().adopt(item.get(),objId) ;
	item.release(); 
	return id; 

}

/** 

long ICMLOCAL_CstDefCurve (const ARM_Date& AsOfDate, 
						   const std::vector<double>& Yearfractions,
						   const std::vector<double>& Inputs,
								  double Type,
								  double recovery,
								  long ircurveid,
								  int intRule,
								  int adjStartRule,
								  const std::string& ccy, 
								  const std::string& Label,
								  long VolCurveId,
								  const std::string & calibrationMethod,
								  int lag,
								  long prevId)
{

	
	
    ARM_Vector* VInputs = NULL;
	
	double* pYearfractions = NULL;


	// 
	ARM_ZeroCurve* ZC =NULL;
	LocalPersistent::get().convert(ircurveid,ZC); 
	if (!ZC) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_CstDefCurve: IR Curve missing"); 
	//
	ARM_VolCurve* VolCurve = NULL;
	LocalPersistent::get().convert(VolCurveId,VolCurve); 

	int i=0;
	int real_size = Yearfractions.size ();
	int real_size2 = real_size;

	pYearfractions = new double[real_size];

	for (i=0; i<real_size; i++)
		pYearfractions[i] = Yearfractions[i];

	double * pdInputs = NULL;
	// char* myCurveDate=NULL;

	int nbdefvalues = 0;
	int i2 = 0;

	CCString msg ("");


		for(i = 0; i < real_size; i++) { if (Inputs[i] == -999.) nbdefvalues+=1; }
		real_size2 -=nbdefvalues;


	 
		pdInputs = new double[real_size2];
		
		i2=0;
		for(i = 0; i < real_size; i++)
		{ if (Inputs[i] != -999.) {pdInputs[i2] = Inputs[i] ;i2++;} }

		VInputs  = new ARM_Vector(real_size2,pdInputs);



		ICM_DefaultCurve* createdDC  =  new ICM_Constant_Piecewise(AsOfDate,
												pYearfractions, 
												VInputs, 
												Type,
												recovery,
												ZC,
												intRule,	// intRule
												adjStartRule,	// adjStartDate
												qCredit_Default,
												ccy, 
												Label,
												VolCurve,
												qDEFCURVE_DICHO,
												calibrationMethod, 
												lag);

		if (!createdDC) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_CstDefCurve:Can't create"); 

		if (pdInputs)
			delete[] pdInputs;
		pdInputs = NULL;

		if (VInputs)
			delete VInputs;
		VInputs = NULL;


		long objId = LocalPersistent::get().adopt(createdDC,prevId); 


	return objId; 
}
**/ 
long ICMLOCAL_CstDefCurve_Dates (const ARM_Date& AsOfDate, 
										VECTOR<double>& Dates,
										VECTOR<double>& Inputs,
										double Type,
										double recovery,
										long ircurveid,
										const CCString& ccy, 
										const CCString& Label,
										long VolCurveId,
										const std::string& calibrationMethod,
										int lag,
										long paramId,
										long prevId)
{
// 	long curveId;

    ARM_Vector  VInputs  ;
	ARM_Vector  VDates  ;
	
	bool flagdate = true;

	ARM_ZeroCurve* ZC =NULL;
	LocalPersistent::get().convert(ircurveid,ZC); 
	if (!ZC) ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_CstDefCurve_Dates: mising IR Curve"); 
	ICM_Parameters* params=NULL; 
	LocalPersistent::get().convert(paramId,params); 
	
	ARM_VolCurve* VolCurve = NULL;
	LocalPersistent::get().convert(VolCurveId,VolCurve); 

	int i=0;
	int real_size = Dates.size ();
	int real_size2 = real_size;
	int real_size3 = real_size;


	char* tmp = NULL;
	double * pdInputs = NULL;
	double * pdDates = NULL;
	char* AuxDate=NULL;

	int nbdefvalues = 0;
	int i2 = 0, i3= 0;

	CCString msg ("");


		for(i = 0; i < real_size; i++) { if (Inputs[i] == -999.) nbdefvalues+=1; }
		real_size2 -=nbdefvalues;


		pdInputs = new double[real_size2];
		pdDates = new double[real_size3];
		
		for(i = 0; i < real_size; i++)
		{ if (Inputs[i] != -999.) {pdInputs[i2] = Inputs[i] ;i2++;} }

		for(i = 0; i < real_size; i++)
		{
			if ((Dates[i] != -999.)&&(Dates[i] >100.))
			{
				AuxDate = new char[11] ;
				Local_XLDATE2ARMDATE(Dates[i],AuxDate);

				pdDates[i3] = ((ARM_Date)AuxDate).GetJulian() ;

				if (AuxDate)
					delete [] AuxDate;
				AuxDate = NULL;

				i3++;
			}
			else
			{ pdDates[i3] = Dates[i]; flagdate=false;}
		}

		VInputs =ARM_Vector(real_size2,pdInputs);
		VDates = ARM_Vector(real_size3,pdDates);

		// char psMatu [ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS];
		// for (int j = 0; j < ARM_NB_TERMS; j++)
		// 	sprintf(psMatu[j],"X");


	ICM_DefaultCurve* createdDC = new ICM_Constant_Piecewise(AsOfDate,
												VDates, 
												VInputs, 
												(int)Type,
												recovery,
												ZC,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
												qCredit_Default,
												CCSTringToSTLString(ccy), 
												CCSTringToSTLString(Label),
												VolCurve,qDEFCURVE_DICHO,
												calibrationMethod, lag,
												params?*params:ICM_Parameters());
	
	long objId = LocalPersistent::get().adopt(createdDC,prevId); 

		if (pdInputs)
			delete[] pdInputs;
		pdInputs = NULL;


		if (pdDates)
			delete[] pdDates;
		pdDates = NULL;



		return objId; 
}




long ICMLOCAL_InputDefCurve_Dates (const ARM_Date& AsOfDate, 
								   const std::vector<ARM_Date>& Dates,
										const ARM_Vector  & Inputs, 
										double recovery,
										long ircurveid,
										const std::string & ccy, 
										const std::string & Label,
										// qINTERPOL_TYPE InterpolType, 
										// ARM_result& result, 
										long prevId)
{

	ARM_ZeroCurve* zeroCurve =NULL ;
	LocalPersistent::get().convert(ircurveid,zeroCurve); 
	if (!zeroCurve) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_InputDefCurve_Dates: ZeroCurve not provided"); 

	ARM_Vector dates(Dates.size()); 
	for(int i=0;i<Dates.size();i++) dates[i]=Dates[i].GetJulian();
	ICM_InterpolDefCrv* curve = new ICM_InterpolDefCrv(  AsOfDate,
												dates, 
												Inputs, 
												recovery,
												zeroCurve,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
												qCredit_Adjust20, 
												ccy, 
												Label
												// InterpolType
												); 
	long objId = LocalPersistent::get().adopt(curve,prevId); 
	return objId; 
/** 
	ICM_DefaultCurve* createdDC = NULL;
	ICM_DefaultCurve* prevDC = NULL;
	ARM_ZeroCurve* ZC =NULL;
	// ARM_Currency* aCcy = NULL;
    ARM_Vector* VInputs = NULL;
	ARM_Vector* VDates = NULL;
	bool flagdate = true;

	int i=0;
	int real_size = Dates.size ();
	int real_size2 = real_size;
	int real_size3 = real_size;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	//char* tmp = NULL;
	double * pdInputs = NULL;
	double * pdDates = NULL;
	char* myCurveDate=NULL;
	char* AuxDate=NULL;

	int nbdefvalues = 0;
	int i2 = 0, i3= 0;

	CCString msg ("");

    try
	{

		for(i = 0; i < real_size; i++) { if (Inputs[i] == -999.) nbdefvalues+=1; }
		real_size2 -=nbdefvalues;

		myCurveDate = new char[11];

		Local_XLDATE2ARMDATE(AsOfDate,myCurveDate);
		ZC = (ARM_ZeroCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(ircurveid);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(ZC, ARM_ZERO_CURVE) == 0)
		{
			result.setMsg ("ARM_ERR: previous object IRcurve is not of a good type");
			return ARM_KO;
		}

	 

		pdInputs = new double[real_size2];
		pdDates = new double[real_size3];
		
		for(i = 0; i < real_size; i++)
		{ if (Inputs[i] != -999.) {pdInputs[i2] = Inputs[i] ;i2++;} }

		for(i = 0; i < real_size; i++)
		{
			if ((Dates[i] != -999.)&&(Dates[i] >100.))
			{
				AuxDate = new char[11] ;
				Local_XLDATE2ARMDATE(Dates[i],AuxDate);

				pdDates[i3] = ((ARM_Date)AuxDate).GetJulian() ;

				if (AuxDate)
					delete [] AuxDate;
				AuxDate = NULL;

				i3++;
			}
			else
			{ pdDates[i3] = Dates[i]; flagdate=false;}
		}

		VInputs  = new ARM_Vector(real_size2,pdInputs);
		VDates  = new ARM_Vector(real_size3,pdDates);

		// char psMatu [ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS];
		// for (int j = 0; j < ARM_NB_TERMS; j++)
		// 	sprintf(psMatu[j],"X");



	
	
	createdDC = new ICM_InterpolDefCrv((ARM_Date) myCurveDate,
												VDates, 
												VInputs, 
												recovery,
												ZC,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
												qCredit_Adjust20, 
												CCSTringToSTLString(ccy), 
												CCSTringToSTLString(Label)
												// InterpolType
												); 

		if (pdInputs)
			delete[] pdInputs;
		pdInputs = NULL;

		if (VInputs)
			delete VInputs;
		VInputs = NULL;

		if (pdDates)
			delete[] pdDates;
		pdDates = NULL;

		if (VDates)
			delete VDates;
		VDates = NULL;

		if (myCurveDate)
			delete [] myCurveDate;
		myCurveDate = NULL;


		if (createdDC == NULL)
		{
			result.setMsg ("ARM_ERR: Curve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdDC);

			if (curveId == RET_KO)
			{
				if (createdDC)
					delete createdDC;
				createdDC = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			prevDC = (ICM_DefaultCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevDC, ICM_INTERPOL_DEFCRV) == 1)
			{
				if (prevDC)
				{
					delete prevDC;
					prevDC = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdDC, objId);

				return ARM_OK;
			}
			else
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (pdInputs)
			delete[] pdInputs;
		pdInputs = NULL;

		if (VInputs)
			delete VInputs;
		VInputs = NULL;

		if (myCurveDate)
			delete [] myCurveDate;
		myCurveDate = NULL;

		ARM_RESULT();
    }
**/ 	
}

/**
long ICMLOCAL_DefaultCurveIndexPWC (const ARM_Date& AsOf, 
									const std::vector<std::string>& Tenors, 
							   VECTOR<double>& Rates,
							   VECTOR<double>& RefRates,
							   double recovery,
							   long ircurveid,
							   const CCString& ccy, 
							   const CCString& Label, 
							   qCDS_ADJ AdjCalType,
							   bool IsSummitCurve,
							   qDEFCURVE_CALIB_ALGO calibrationAlgo,
							   long VolCurveId,
							   const std::string& calibrationData,
							   int lag,
							   long prevId)
{
	

	
	
	ARM_ZeroCurve* ZC =NULL;
	LocalPersistent::get().convert(ircurveid,ZC); 
	if (!ZC) 
		ICMTHROW(ERR_INVALID_ARGUMENT,""); 

    ARM_Vector* VRates = NULL;
	ARM_Vector* VRefRates = NULL;
	ARM_VolCurve* VolCurve = NULL;
	LocalPersistent::get().convert(VolCurveId,VolCurve); 

	int i=0;
	int real_size = Tenors.size ();
	int real_size2 = real_size;


	double * pdRate = NULL;
	double * pdRefRate = NULL;

	int nbdefvalues = 0;
	int i2 = 0;

	CCString msg ("");

		for(i = 0; i < real_size; i++) { if (Rates[i] == -999.) nbdefvalues+=1; }
		real_size2 -=nbdefvalues;

	

	
		pdRate = new double[real_size2];
		pdRefRate = new double[real_size2];
		
		i2=0;
		for(i = 0; i < real_size; i++)
		{ if (Rates[i] != -999.) 
			{pdRate[i2] = Rates[i] / 10000.;
		     pdRefRate[i2] = RefRates[i] / 10000.;
			 i2++;}
		}

		VRates  = new ARM_Vector(real_size2,pdRate);
		VRefRates  = new ARM_Vector(real_size2,pdRefRate);

		vector<string> terms;

		for(int j = 0; j < real_size; j++) {
			if (Rates[j] != -999.) {terms.push_back( Tenors[j]  ); } 
		}

		long PayFreq = K_QUARTERLY;

		ICM_DefaultCurve* createdDC = new ICM_DefcurveIndex(  AsOf,
												terms, 
												VRates, 
												VRefRates,
												recovery,
												ZC,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
												 AdjCalType, 
												CCSTringToSTLString(ccy), 
												CCSTringToSTLString(Label),
												IsSummitCurve,
												VolCurve,
												PayFreq,
												calibrationAlgo,
												calibrationData, lag);

		long objId = LocalPersistent::get().adopt(createdDC,prevId); 
		if (pdRate)
			delete[] pdRate;
		pdRate = NULL;

		if (VRates)
			delete VRates;
		VRates = NULL;

		if (pdRefRate)
			delete[] pdRefRate;
		pdRefRate = NULL;

		if (VRefRates)
			delete VRefRates;
		VRefRates = NULL;

		return objId; 
}

**/ 

long ICMLOCAL_ImpliedLossTree(double AsOfDate, 
							  long DefCurveId,
							  long VolAsCorrelId,
							  long correltype,
							  bool adjusted,
							  int step,
							  ARM_result& result, 
							  long objId)
{
	long curveId;

	ICM_DefaultCurve* pDefCurve = NULL;
	ICM_Smile_Correlation* pVolCurve = NULL;
	ICM_ImpLossTree* pImpLossTree = NULL;
	ICM_ImpLossTree* pPrevImpLossTree = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

    try
	{

		pDefCurve = (ICM_DefaultCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(DefCurveId);
		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pDefCurve, ICM_DEFAULTCURVE) == 0)
		{
			result.setMsg ("ARM_ERR: previous object DefCurve is not of a good type");
			return ARM_KO;
		}

		if (VolAsCorrelId != -1)
		{
			pVolCurve = (ICM_Smile_Correlation *) LOCAL_PERSISTENT_OBJECTS->GetObject(VolAsCorrelId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pVolCurve, ICM_CORRELATION) == 0)
			{
				result.setMsg ("ARM_ERR: previous object VolCurve is not of a good type");
				return ARM_KO;
			}
		}

		pImpLossTree = new ICM_ImpLossTree(pDefCurve,pVolCurve,(qCAL_INDEX_CORR_TYPE)correltype,adjusted,step);


		if (pImpLossTree == NULL)
		{
			result.setMsg ("ARM_ERR: Tree is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)pImpLossTree);

			if (curveId == RET_KO)
			{
				if (pImpLossTree)
					delete pImpLossTree;
				pImpLossTree = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);

			return ARM_OK;
		}
		else
		{
			pPrevImpLossTree = (ICM_ImpLossTree*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pPrevImpLossTree,  ICM_IMPLIED_LOSS_TREE) == 1)
			{
				if (pPrevImpLossTree)
				{
					delete pPrevImpLossTree;
					pPrevImpLossTree = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)pImpLossTree, objId);

				return ARM_OK;
			}
			else
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
    }
}


long ICMLOCAL_ModifyLossTree (long TreeId,
							  long indexno,
							  const VECTOR<double>& yearterms,
							  const VECTOR<double>& losses,
							  const VECTOR<double>& transprobas,
							  const VECTOR<double>& statelosses)
{
	ICM_ImpLossTree* tree = NULL;

	CCString msg ("");

	try
	{
		tree = (ICM_ImpLossTree*) LOCAL_PERSISTENT_OBJECTS->GetObject(TreeId);
		
		if (tree->IsCalibrated()==false)
		{
			ICM_QMatrix<double> mat(1,yearterms.size());

			for (int l=0;l<yearterms.size();l++)
				{mat(0,l) = yearterms[l];}

			tree->SetTimeStep(mat);
			tree->CreateTrees();
		}

		for (int i=0;i<yearterms.size();i++)
			for (int j=0;j<MIN(i+1,yearterms.size());j++)
			{
			int idx = -1;
			for (int k=0;k<tree->GetTreeVector(indexno).depth();k++)
				{if (fabs(yearterms[i]-tree->GetTreeVector(indexno).Time(k))<1.e-3) 
				{idx = k;break;}}

				if ((statelosses.size()>0)&&(idx>=0))
				{
				tree->GetTreeVector(indexno).data(idx,(int)losses[j]).state=statelosses[i*losses.size()+j];
				tree->CptTransitionProba();
				}

				if ((transprobas.size()>0)&&(idx>=0))
				{
				tree->GetTreeVector(indexno).data(idx,(int)losses[j]).p_def=transprobas[i*losses.size()+j];
				tree->GetTreeVector(indexno).data(idx,(int)losses[j]).p_nodef=1.-transprobas[i*losses.size()+j];
				}
			}

		tree->SetCalibrated(true);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
	}

	return ARM_KO;
}
long ICMLOCAL_DefProbInverse (long idCurve, 
							 double DefaultProba,
							 ARM_result& result)
{
	double dResult;
	ICM_DefaultCurve* DC = NULL;

	ARM_CLASS_NAME objname = ARM_ZERO_CURVE;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		DC = (ICM_DefaultCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(idCurve);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(DC, ICM_DEFAULTCURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Default Curve is not of a good type");
			return ARM_KO;
		}

		dResult = DC->DefProbInverse(DefaultProba);

		result.setDouble(dResult);
		return ARM_OK;
	}

	catch(Exception& x)
	{

		x.DebugPrint();
		ARM_RESULT();
	}
}


long ICMLOCAL_QuickDefaultCurve (double spread,
								 double recovery,
								 CCString label,
							     ARM_result& result, 
							     long objId)
{
	long curveId;

	ICM_DefaultCurve* createdDC = NULL;
	ICM_DefaultCurve* prevDC = NULL;

	int i=0;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}
	
	CCString msg ("");

    try
	{

		createdDC = new ICM_DefCurveSimple(spread,recovery,CCSTringToSTLString(label));


		if (createdDC == NULL)
		{
			result.setMsg ("ARM_ERR: Curve is null");
			return ARM_KO;
		}

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			curveId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdDC);

			if (curveId == RET_KO)
			{
				if (createdDC)
					delete createdDC;
				createdDC = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(curveId);
			return ARM_OK;
		}
		else
		{
			prevDC = (ICM_DefaultCurve*) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(prevDC,ICM_DEFCURVE_QUICK) == 1)
			{
				if (prevDC)
				{
					delete prevDC;
					prevDC = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdDC, objId);
				return ARM_OK;
			}
			else
			{
				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}

	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
    } 
	catch(std::exception& se)
	{
		// msg erreur
		CCString msg (" std exception : ");
		msg+= se.what();
		Exception x(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (char*) msg.c_str());
		x.DebugPrint();
		ARM_RESULT();
	}
	catch(...)
	{
		// msg erreur
		CCString msg (" unknown exception ");
		Exception x(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (char*) msg.c_str());
		x.DebugPrint();
		ARM_RESULT();
	}

}

long ICMLOCAL_Fixing_Curve (const std::vector<ARM_Date>& vDates,
								   const std::vector<double>& vValues,
								   const ARM_Date& asOf,
								   const string& indexName,
								   long indexNameID,
							       ARM_result& result, 
								   long objId)
{
	try {
		ARM_IRIndex* pIndex = NULL;
		std::map<ARM_Date, double> lMap;
		if (indexName.empty())
		{
			if (indexNameID == -1) {
				ICMTHROW(ERR_INVALID_ARGUMENT, "ICMLOCAL_Fixing_Curve IndexId is no good");
			} else {
				pIndex = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(indexNameID);
				if ((LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pIndex, ARM_IRINDEX) == 0 ) && 
					(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pIndex, ICM_CREDIT_INDEX) == 0))
				{
					ICMTHROW(ERR_INVALID_ARGUMENT, "ICMLOCAL_Fixing_Curve : Pay Index is not of a type ARM_IRINDEX nor ICM_CREDIT_INDEX");
				}
			}
		}
		if (vDates.size() != vValues.size()) {
				ICMTHROW(ERR_INVALID_ARGUMENT, "ICMLOCAL_Fixing_Curve : Dates vector and values vector don't have the same size");
		}	
		for (int i=0; i<vDates.size(); i++){
			lMap[vDates[i]] = vValues[i];
		}
		ICM_Fixing_Curve*  pNewFixingCurve = new ICM_Fixing_Curve(indexName, asOf, lMap, pIndex);
		long QId = LocalPersistent::get().adopt(pNewFixingCurve,objId );
		return QId; 
	}
	catch(Exception& x)
    {
		CCString msg ("");
		x.DebugPrint();
		ARM_RESULT();
	}
	catch(...){
		result.setMsg ("Unknown exception");
		return ARM_KO;
	}
}
	


long ICMLOCAL_CptImplicitSpreadInterpolOLD (long defcurveId, 
									CCString plot,
									double date,
									double slope,
									double Exactdate,
									VECTOR<double>& output,
									ARM_result& result)
{
	double dResult;
	double Recovery;

	ICM_DefaultCurve* zc = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	ARM_Date NewDate;
	double ModDate = 0.;

	if (date)
	{
	char* sDate = new char[11];
	Local_XLDATE2ARMDATE(date,sDate);
	NewDate = (ARM_Date)sDate;
	ModDate = NewDate.GetJulian();
	}

	ARM_Date ExactNewDate;
	double ExactModDate = 0.;

	if (Exactdate)
	{
	char* sDate = new char[11];
	Local_XLDATE2ARMDATE(Exactdate,sDate);
	ExactNewDate = (ARM_Date)sDate;
	ExactModDate = ExactNewDate.GetJulian();
	}
	
	CCString msg ("");

	try
	{
		zc = (ICM_DefaultCurve *) LOCAL_PERSISTENT_OBJECTS->GetObject(defcurveId);

		if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(zc, ICM_DEFAULTCURVE) == 0)
		{
			result.setMsg ("ARM_ERR: Default Zc Curve is not of a good type");
			return ARM_KO;
		}

		dResult = zc->ImpliedSpreadInterpolOLD((char*)plot,slope,ExactModDate,ModDate);
		Recovery = 	zc->GetRecovery();

		output.push_back(dResult);
		output.push_back(Recovery);

		return ARM_OK;
	}

	catch(Exception& x)
	{

		x.DebugPrint();

		ARM_RESULT();
	}

	return ARM_KO;
}	


long ICMLOCAL_ABS_PWCDefaultCurve (const  ARM_Date& AsOf, 
							   const std::vector<std::string> & Tenors, 
							   VECTOR<double>& Rates,
							   double recovery,
							   long ircurveid,
							   const std::string& ccy, 
							   const std::string& Label, 
							   qCDS_ADJ AdjCalType,
							   bool IsSummitCurve,
							   qDEFCURVE_CALIB_ALGO calibrationAlgo,
							   long VolCurveId,
							   const std::string& calibrationData,
							   int Lag,
							   VECTOR<double>& Upfront, 
							   VECTOR<long>& RefValId, 
							   long prevId)
{
    ARM_Vector VRates; //= NULL;
    ARM_Vector* VUpfront= NULL;
	vector<ARM_ReferenceValue*> VRefId;VRefId.resize(RefValId.size());

	// non optional ZeroCurve 
	ARM_ZeroCurve* ZC =NULL;
	LocalPersistent::get().convert(ircurveid,ZC); 
	if (!ZC) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICMLOCAL_ConstantDefaultCurve: ZeroCurve not provided."); 
	
	// optional VolCurve
	ARM_VolCurve* VolCurve = NULL;
	LocalPersistent::get().convert(VolCurveId,VolCurve); 

	int i=0;
	int real_size = Tenors.size ();
	int real_size2 = real_size;

	vector<double> vRate; 

	int nbdefvalues = 0;
	int i2 = 0;

	CCString msg ("");

	if (Upfront.size()>0) 
	{VUpfront = new ARM_Vector(Upfront.size(),Upfront.begin());}

	for(i = 0; i < RefValId.size(); i++) 
	{ LocalPersistent::get().convert(RefValId[i],VRefId[i]); }
	
	for(i = 0; i < real_size; i++) { if (Rates[i] == -999.) nbdefvalues+=1; }
	real_size2 -=nbdefvalues;

	vRate.resize(real_size2);

	i2=0;
	for(i = 0; i < real_size; i++)
	{ if (Rates[i] != -999.) {vRate[i2] = Rates[i] / 10000.;i2++;} }

	VRates  = ARM_Vector(vRate);

	vector<string> terms;

	i2=0;
	for(int j = 0; j < real_size; j++) {
		if (Rates[j] != -999.) {terms.push_back(Tenors[j]);}
	}

	long PayFreq = K_QUARTERLY;


	ICM_DefaultCurve* createdDC = new ICM_Constant_Piecewise(  AsOf,
												terms, 
												&VRates, 
												recovery,
												ZC,
												K_ADJUSTED,	// intRule
												K_ADJUSTED,	// adjStartDate
												 AdjCalType, 
												ccy, 
												Label,
												IsSummitCurve,
												//2 NULL,
												VolCurve,
												PayFreq,
												calibrationAlgo,
												calibrationData,
												Lag,
												ICM_Parameters(),
												VUpfront,
												VRefId);
	if (!createdDC) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Can't create ICM_Constant_Piecewise"); 

	long objId = LocalPersistent::get().adopt(createdDC,prevId); 

	return objId; 
}