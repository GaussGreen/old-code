#include "firsttoinc.h"

#include "model.h"
#include "security.h"
#include "ICMKernel\pricer\icm_pricer.h"
#include "ICMKernel\crv\icm_distriblosscalculator.h"
#include "ICMKernel/glob/icm_mktdatamng.h"


double FeesPricing(ARM_ReferenceValue* ref,ARM_ZeroCurve* zc);
//	only called by ctors. 
void ICM_Pricer::Init()
{
	SetName(ICM_PRICER);

	itsSecurity = NULL;
	itsModel = NULL;
	itsInitialPrice = 0.;
	itsInitialPriceFlg = false;
	itsFaster = false;
	itsFeeLegPrice_UnityFlg = false;

	// 17783 itsSensiManager = NULL;

	//itsParameters = NULL;
	
	// 
	for(int i=0;i<qCMPLAST;i++) itsFlags[i]=false; 

	//	---------------------------------------
	//	This is the ARM_Vector* cache 
	for(i=0;i<qVECTLAST;i++) itsObjectValues[i]=0; 
}

ICM_Pricer::~ICM_Pricer(void) 
{
	// 17783 if (itsSensiManager)
	// 17783 	delete itsSensiManager;
	// 17783 itsSensiManager = NULL;

	// if (itsParameters)
	//	delete itsParameters;
	// itsParameters = NULL;
	for(int i=0;i<qVECTLAST;i++) 
	{
		if (itsObjectValues[i]) delete itsObjectValues[i]; 
	}

	std::vector<ARM_Vector*>::iterator it = itsVectorMeasureValues.begin(); 
	for(;it!=itsVectorMeasureValues.end();++it) delete *it; 
}
void ICM_Pricer::SetPrice(const double& price)  
{
	setValue(qCMPPRICE,price); 

	if (!itsInitialPriceFlg)
	{	itsInitialPrice = price;
		itsInitialPriceFlg = true;}
}
// 17783 void ICM_Pricer::SetSensiManager(ICM_SensiManager*item) 
// 17783 {
// 17783 	if (itsSensiManager==item) return; 
// 17783 	if (itsSensiManager) delete itsSensiManager; 
// 17783 	itsSensiManager=item; 
// 17783 }
void ICM_Pricer::SetMktDataMng(ICM_MktDataMng* model) 
{
	if (itsModel == model)
		return;
	itsModel = model;
}
void ICM_Pricer::ResetPricer(void)  
{
	Reset();
	ResetRootPricer();
}
//	------------------------------------------------------------------------------------------------------
void 
ICM_Pricer::unsetFlgs()
{
	for(int i=0;i<qCMPLAST;i++) itsFlags[i]=false; 
	for(i=0;i<qVECTLAST;i++) { 
	if (itsObjectValues[i])
		delete itsObjectValues[i] ; 
	itsObjectValues[i]= NULL ; }
	itsMeasureNames.clear(); 
	itsMeasureValues.clear(); 
}
//	------------------------------------------------------------------------------------------------------
ARM_Object* 
ICM_Pricer::Clone() 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Cloning is not available for pricers"); 
} 

//	JLA Added
void ICM_Pricer::Set(ARM_Security*sec,ARM_Object*mod,const ICM_Parameters&params,const ARM_Date*asof)
{
	itsSecurity = sec;
	itsModel = mod;
	SetParameters(params); 

	if (asof) 
		SetAsOfDate(*asof); 
	else 
	{
		ARM_Model*mod_=dynamic_cast<ARM_Model*>(mod); 
		if (!mod_) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"Can't find AsOf date... "); 
		SetAsOfDate(mod_->GetStartDate());
	}
// 17783 	SetSensiManager(new ICM_SensiManager()); 
}

//	JLA - Accessors for the pricer parameters
//		returns false on 
//	-----------------------------------------------------------------------------------
bool 
ICM_Pricer::getParam(const std::string&paramName,double&ret,bool throwOnError) const
{
	ARM_Vector local; 
	if (!getParam(paramName,local,throwOnError)) return false ; 
	if (local.GetNumLines()==0) 
	{
		if (throwOnError) ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer::getParam "<<paramName<<": empty value! "); 
		return false ; 
	}
	ret=local.Elt(0); 
	return true; 
}
//	-----------------------------------------------------------------------------------
bool 
ICM_Pricer::getParam(const std::string&paramName,long&ret,bool throwOnError) const
{
	ARM_Vector local; 
	if (!getParam(paramName,local,throwOnError)) return false ; 
	if (local.GetNumLines()==0) 
	{
		if (throwOnError) ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer::getParam "<<paramName<<": empty value! "); 
		return false ; 
	}
	ret=local.Elt(0); 
	return true; 
}
//	-----------------------------------------------------------------------------------
bool 
ICM_Pricer::getParam(const std::string&paramName,ARM_Vector&ret,bool throwOnError) const
{
	// We might not have any parameter defined... 
	// if (!itsParameters)
	// {
	// 	ret.Resize(0); 
	// 	if (throwOnError) ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer::getParam "<<paramName<<": no parameters !"); 
	// 	return false ;
	// }
	// We might not found the given parameter name
	/// UGLY FIX TO MAKE COMPIL WORK
	ARM_Vector* item = itsParameters.GetColVect( const_cast<char*>(paramName.c_str() ) ); 
	if (item==0) 
	{
		ret.Resize(0); 
		if (throwOnError) ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer::getParam "<<paramName<<": not found !"); 
		return false ;
	}
	ret=*item; 
	return true; 
}
//	-------------------------------------------------------------------------------------------------
double ICM_Pricer::Price(qCMPMETH measure)
{
	return ComputePrice(measure); 
}
//	-------------------------------------------------------------------------------------------------
double ICM_Pricer::Price(const std::string& measureName)
{
	// -- cache not managerd here.
	std::string valueList ;
	qCMPMETH measure ;
	bool isOk; 
	ICM_EnumsCnv::cnv(measureName,measure,isOk,valueList ); 
	if (isOk) return Price(measure); 
	double tmp ;
	if (!getFlg(measureName)) 
		tmp = DoPrice(measureName); 
	return getValue(measureName); 
}
//	-------------------------------------------------------------------------------------------------
double 
ICM_Pricer::ComputePrice(qCMPMETH mode)  
{	double Price = CheckForPrice(mode);
	return (Price);  
}
//	-------------------------------------------------------------------------------------------------
const ARM_Vector& 
ICM_Pricer::PriceVector(qVECTMETH measure)
{
	if (!getFlgVector(measure)) DoPriceVector(measure); 
	return getVectorValue(measure); 
}
//	-------------------------------------------------------------------------------------------------
const ARM_Vector& 
ICM_Pricer::PriceVector(const std::string& measureName)
{
	std::string valueList ;
	qVECTMETH measure ;
	bool isOk; 
	ICM_EnumsCnv::cnv(measureName,measure,isOk,valueList ); 
	if (isOk) return PriceVector(measure); 
	if (!getFlgVector(measureName)) DoPriceVector(measureName); 
	return getVectorValue(measureName); 
}
//	-------------------------------------------------------------------------------------------------
double 
ICM_Pricer::DoPrice(qCMPMETH measure) 
{
	ICMTHROW (ERR_UNIMP_METHOD_CALL,"Unimplemented <DoPrice> method "<<measure); return (-99999999999.); 
}
//	-------------------------------------------------------------------------------------------------
double 
ICM_Pricer::DoPrice(const std::string& measure) 
{
	ICMTHROW (ERR_UNIMP_METHOD_CALL,"Unimplemented <DoPrice> method " << measure); return (-99999999999.); 
}
//	-------------------------------------------------------------------------------------------------
void
ICM_Pricer::DoPriceVector(const std::string& measure) 
{
	ICMTHROW (ERR_UNIMP_METHOD_CALL,"Unimplemented <DoPriceVector> method " << measure);  
}
//	-------------------------------------------------------------------------------------------------
void 
ICM_Pricer::DoPriceVector(qVECTMETH measure) 
{
	ICMTHROW (ERR_UNIMP_METHOD_CALL,"Unimplemented <DoPriceVector> method " << measure); 
}
//	-------------------------------------------------------------------------------------------------
double ICM_Pricer::CheckForPrice(qCMPMETH pricingmode)
{
	double result = 0.;

	switch (pricingmode)
	{
	case qCMPPRICE: // case default avant !!
		if (GetPriceFlg()) result = GetPrice() ;
		else 
		{ 
		  result = FeeLegPV()-DefLegPV();	
		  SetPrice(result);
		}
		break;
	case qCMPSPREAD:
		if (GetSpreadFlg())	result = GetSpread();
		else {result = ComputeSpread();}
		break;
	case qCMPFEELEGPV:
		if (GetFeeLegPriceFlg()) result = GetFeeLegPrice();
		else {result = FeeLegPV();}
		break;
	case qCMPDEFLEGPV:
		if (GetDefLegPriceFlg()) result = GetDefLegPrice();
		else {result = DefLegPV();}
		break;
	case qCMPACCRUED:
		/** if (GetAccruedFlg()) result = GetAccrued();
		else {result = Accrued();}
		**/ 
		if (getFlg(qCMPACCRUED)) result = itsValues[qCMPACCRUED]; 
		else result = Accrued();
		break;
	case qCMPDURATION:
		if (getFlg(qCMPDURATION)) result=getValue(qCMPDURATION);
		else result = ComputeDuration() ;
		break;
	case qCMPFEES :
		if (GetSecurity()->GetFee()) {result=FeesPricing(GetSecurity()->GetFee(),
														 GetModel()->GetDiscountCurve());}
		break;
	case qCMPEL:
		if (getFlg(qCMPEL)) result=getValue(qCMPEL);; 
		result = CptPtf_ELoss(this) ;
		break;
	/*case qCMPFLATCORR:
		computeFlatCorrel();
		break;*/
	default:
		result = DoPrice(pricingmode);
	}

	return (result);
}

void ICM_Pricer::SetPriceFollowMode(qCMPMETH pricingmode,const double& value)
{
	double result = 0.;

	switch (pricingmode)
	{
	case qCMPSPREAD:
		SetSpread(value);
		break;
	case qCMPFEELEGPV:
		SetFeeLegPrice(value);
		break;
	case qCMPDEFLEGPV:
		SetDefLegPrice(value);
		break;
	case qCMPACCRUED:
		/** SetAccrued(value); **/ 
		setValue(qCMPACCRUED,value); 
		break;
	case qCMPDURATION:
		SetDuration(value);
		break;
	case qCMPEL:
		/** SetAccrued(value); **/ 
		setValue(qCMPEL,value); 
		break;
	default:
		SetPrice(value);
	}
}

double FeesPricing(ARM_ReferenceValue* ref,ARM_ZeroCurve* zc){
	ARM_Vector* dates = ref->GetDiscreteDates();
	ARM_Vector* values = ref->GetDiscreteValues();
	double MAX_DAT = 1000.;
	double FeePv = 0.;

	for (int i=0; i<dates->GetSize(); i++)
	{
		if (dates->Elt(i)<MAX_DAT) {
			ICMTHROW(ERR_INVALID_ARGUMENT,"Dates must'nt be yearfractions !"); }
		FeePv += zc->DiscountPrice((ARM_Date)dates->Elt(i))*values->Elt(i);	}

	return (FeePv);
}
//	-----------------------------------------------------------------------------------
void 
ICM_Pricer::ResetRootPricer()  
{
	itsFeeLegPrice_UnityFlg = false;
	unsetFlgs(); 
}
//	-----------------------------------------------------------------------------------
void ICM_Pricer::View(char* id, FILE* ficOut)
{
	FILE* fOut;
	char  fOutName[200];

	if ( ficOut == NULL )
	{
	ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

	(void) unlink(fOutName);

	fOut = fopen(fOutName, "w"); 
	}
	else
	{
	fOut = ficOut;
	} 

	int size =0;

	fprintf(fOut, "\t\t\t ----------------- Basis Pricer ----------------- \n\n");
	fprintf(fOut, "\tType\t\t\t\t\tStatus\n\n");

	
	if (getFlg(qCMPPRICE)) fprintf(fOut, "\tPrice: %.2lf\t\t\t%s\n",getValue(qCMPPRICE),"Ok");
	if (itsInitialPriceFlg) fprintf(fOut, "\tInitial Price: %.2lf\t\t%s\n",itsInitialPrice,"Ok");
	if (getFlg(qCMPFEELEGPV)) fprintf(fOut, "\tFee Leg Price: %.2lf\t\t%s\n",getValue(qCMPFEELEGPV),"Ok");
	if (itsFeeLegPrice_UnityFlg) fprintf(fOut, "\tFee Leg Unity Price: %.2lf\t\t%s\n",itsFeeLegPrice_Unity,"Ok");
	if (getFlg(qCMPDEFLEGPV)) fprintf(fOut, "\tDefault Leg Price: %.2lf\t%s\n",getValue(qCMPDEFLEGPV),"Ok");
	if (getFlg(qCMPSPREAD)) fprintf(fOut, "\tBreak even spread: %.2lf\t\t\t%s\n",getValue(qCMPSPREAD),"Ok");
	if (getFlg(qCMPDURATION)) fprintf(fOut, "\tDuration: %.2lf\t\t\t%s\n",getValue(qCMPDURATION),"Ok");
	
	fprintf(fOut, "\n");


	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}

//virtual 
double ICM_Pricer::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon , double  epsilonGamma)
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer::ComputeSensitivity: Not implemented for "
		<<typeid(*this).name()); 
}

// virtual 
void ICM_Pricer::SetModel(ARM_Model* model)
{
	if (itsModel == model) return;
	itsModel = model;
	if (itsModel) PropagateModel((ARM_Model*)itsModel);
}
