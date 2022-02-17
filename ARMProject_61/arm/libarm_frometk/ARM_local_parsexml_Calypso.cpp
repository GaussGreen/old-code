#include <ARM\libarm_local\firstToBeIncluded.h>

#include <crv\zeroint.h>

#include <ARM\libarm_frometk\arm_local_parsexml_common.h> 
#include <ARM\libarm_frometk\arm_local_parsexml_util.h> 
#include <ARM\libarm_frometk\arm_local_paesexml_calypso.h> 
#include "PaserManagerUtilities.h"

#include <fstream>
#include <libCCatl\CCatl.h>
#include <ARM\libarm_frometk\XMLTools.h>
#include <ARM\libarm_frometk\VariantTools.h>


#import <JavaExcelCom.dll>	// 
using namespace MSXML2;

// this is part of ICM Kernel, but should be located on ARM_Date 
extern std::ostream& operator<<(std::ostream&,const ARM_Date&); 

//	------------------------------------------------------------------------------------
//
//		This is the singleton class holding a JavaExcelCom instance
//
class JavaExcelCom
{
public:	
	static JAVAEXCELCOMLib::IAddinPtr& get()
	{
		if (itsInstance.get()==0) 
		{
			itsInstance=std::auto_ptr<JavaExcelCom>(new JavaExcelCom) ;
		}
		return itsInstance->addinptr(); 
	}
public:
	JAVAEXCELCOMLib::IAddinPtr& addinptr() { return *itsAddinPtr; }
private:
	JavaExcelCom()
	{
		ICMLOG("JavaExcelCom::JavaExcelCom: Creating the COM instance"); 
		itsAddinPtr = new JAVAEXCELCOMLib::IAddinPtr ; 
		itsAddinPtr->CreateInstance(__uuidof(JAVAEXCELCOMLib::Addin)); 
	}
	~JavaExcelCom()
	{
		ICMLOG("JavaExcelCom::~JavaExcelCom"); 
		try 
		{ 
			delete itsAddinPtr ; 
		} 
		catch(...)
		{
		} 
	}
	JavaExcelCom(const JavaExcelCom&); //NA
	JavaExcelCom& operator=(const JavaExcelCom&); // NA 
private:
	static std::auto_ptr<JavaExcelCom> itsInstance; 
	JAVAEXCELCOMLib::IAddinPtr * itsAddinPtr; 
	friend class std::auto_ptr<JavaExcelCom>; 
}; 
std::auto_ptr<JavaExcelCom> JavaExcelCom::itsInstance(0);
//	------------------------------------------------------------------------------------
//	static 
void 
ARM_CalypsoToolkit::readTextFile(const std::string& xmlFileName,std::string& xmlContent)
{
	// ICMLOG("ARM_CalypsoToolkit::readTextFile:"<<xmlFileName<<" ...") ; 
	PreciseChrono c;	
	c.Start(); 
	std::ifstream its(xmlFileName.c_str());
	std::stringstream sstr ; sstr<<its.rdbuf();
	xmlContent=sstr.str(); 
	c.Stop() ;
	// ICMLOG("ARM_CalypsoToolkit::readTextFile:"<<xmlFileName<<":"<<c.GetDurationInMilliseconds()<<" ms") ;
}
//	------------------------------------------------------------------------------------
void 
ARM_CalypsoToolkit::init(const std::string& env)
{
	ICMLOG("ARM_CalypsoToolkit::init... "); 
	PreciseChrono c; 
	c.Start() ;
	try {
		JavaExcelCom::get().CreateInstance(__uuidof(JAVAEXCELCOMLib::Addin)); 
		JavaExcelCom::get()->Calypso_Init(env.c_str()); 
		JavaExcelCom::get()->Calypso_SetLogging(10); 
	}
	catch(_com_error &e)
	{
		std::string desc (e.Description()); 
		std::string src (e.Source()); 
		ICMTHROW(ERR_INVALID_ARGUMENT,src<<":"<<desc); 
	}
	c.Stop() ;
	ICMLOG("ARM_CalypsoToolkit::init:"<<c.GetDurationInMilliseconds()<<" ms"); 
}
//	------------------------------------------------------------------------------------
void 
ARM_CalypsoToolkit::GetProbabilityCurve(const std::string& issuer,
										const std::string& ccy,
										const std::string& seniority,
										const std::string& forceCurveName,
										const std::string& pricingEnv,
										const ARM_Date& date,
										const std::string& xmlFileName, 
										std::string& xmlOutput) 
{
	// ICMLOG("ARM_CalypsoToolkit::GetProbabilityCurve:"<<issuer<<","<<ccy<<","<<seniority<<","<<forceCurveName
	// 	<<","<<pricingEnv<<","<<date<<","<<xmlFileName<<" ..."); 
	ICMELAPSER("ARM_CalypsoToolkit::GetProbabilityCurve:"<<issuer<<","<<ccy<<","<<seniority<<","<<forceCurveName
		<<","<<pricingEnv<<","<<date<<","<<xmlFileName<<" ..."); 
	PreciseChrono c;	
	c.Start(); 

	if (!xmlFileName.empty()) 
	{
		ARM_CalypsoToolkit::readTextFile(xmlFileName,xmlOutput); 
		return ; 
	}

	std::string curveName(forceCurveName) ; 
	if (curveName.empty()) 
		curveName = pricingEnv+"_"+issuer+"_"+seniority+"_ANY_ANY" ;
	
	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetMonth()<<"/"
		<<setw(2)<<setfill('0')<<date.GetDay()<<"/"<<date.GetYear(); 
	_variant_t probaCurve; 
	_variant_t err; 
	try {
		JavaExcelCom::get()->GetProbabilityCurve(curveName.c_str(),pricingEnv.c_str(),sstr.str().c_str(),&probaCurve,&err); 
	}
	catch(_com_error &e)
	{
		std::string desc (e.Description()); 
		std::string src (e.Source()); 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_CalypsoToolkit::GetProbabilityCurve:"<<curveName<<":"<<src<<":"<<desc); 
	}
	xmlOutput=_bstr_t(probaCurve); 
	if (xmlOutput.empty()) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_CalypsoToolkit::GetProbabilityCurve:"<<curveName<<": xml is empty"); 
	c.Stop(); 
	// ICMLOG("ARM_CalypsoToolkit::GetProbabilityCurve:"<<curveName<<":"<<c.GetDurationInMilliseconds()<<" ms");  
}
//	------------------------------------------------------------------------------------
void 
ARM_CalypsoToolkit::GetCurveZero(const std::string& index,
								const std::string& ccy,
								const std::string& term,
								const std::string& forceCurveName,
								const std::string& pricingEnv,
								const ARM_Date& date,
								const std::string& xmlFileName,
								std::string& xmlOutput )
{
	ICMELAPSER("ARM_CalypsoToolkit::GetCurveZero:"<<index<<","<<ccy<<","<<term<<","<<forceCurveName
		<<","<<pricingEnv<<","<<date<<","<<xmlFileName<<" ..."); 

	PreciseChrono c;
	c.Start(); 

	if (!xmlFileName.empty()) 
	{
		ARM_CalypsoToolkit::readTextFile(xmlFileName,xmlOutput); 
		return ; 
	}

	std::string curveName(forceCurveName) ; 
	if (curveName.empty()) 
		curveName = pricingEnv+"_"+ccy+"_"+index+"_"+term ;

	std::stringstream sstr; 
	//sstr.width(2);
	//sstr.fill('0');
	sstr<<setw(2)<<setfill('0')<<date.GetMonth()<<"/"
		<<setw(2)<<setfill('0')<<date.GetDay()<<"/"<<date.GetYear(); 
	_variant_t probaCurve; 
	_variant_t err ; 
	try {
		JavaExcelCom::get()->GetCurveZero(curveName.c_str(),pricingEnv.c_str(),sstr.str().c_str(),&probaCurve,&err); 
	}
	catch(_com_error &e)
	{
		std::string desc (e.Description()); 
		std::string src (e.Source()); 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_CalypsoToolkit::GetCurveZero:"<<curveName<<":"<<src<<":"<<desc); 
	}
	xmlOutput=_bstr_t(probaCurve); 
	if (xmlOutput.empty())
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_CalypsoToolkit::GetCurveZero:"<<curveName<<":xml is empty"); 
	c.Stop();
	// ICMLOG("ARM_CalypsoToolkit::GetCurveZero:"<<curveName<<":"<<c.GetDurationInMilliseconds()<<" ms"); 
}
//	------------------------------------------------------------------------------------
// 	static 
void ARM_CalypsoToolkit::GetBasketCorrel(const std::string& forceCurveName,
										 const std::string& pricingEnv,
										 const ARM_Date& date,
										 const std::string& xmlFileName,
										 std::string& xmlOutput)
{
	ICMELAPSER("ARM_CalypsoToolkit::GetBasketCorrel:"<<forceCurveName
		<<","<<pricingEnv<<","<<date<<","<<xmlFileName<<" ..."); 

	PreciseChrono c;
	c.Start(); 

	if (!xmlFileName.empty()) 
	{
		ARM_CalypsoToolkit::readTextFile(xmlFileName,xmlOutput); 
		return ; 
	}


	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetMonth()<<"/"
		<<setw(2)<<setfill('0')<<date.GetDay()<<"/"<<date.GetYear(); 
	_variant_t basketCorrel; 
	_variant_t err ; 
	try {
		JavaExcelCom::get()->GetBasketCorrel(forceCurveName.c_str(),pricingEnv.c_str(),sstr.str().c_str(),&basketCorrel,&err); 
	}
	catch(_com_error &e)
	{
		std::string desc (e.Description()); 
		std::string src (e.Source()); 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_CalypsoToolkit::GetCurveZero:"<<forceCurveName<<":"<<src<<":"<<desc); 
	}
	xmlOutput=_bstr_t(basketCorrel); 
	if (xmlOutput.empty())
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_CalypsoToolkit::GetBasketCorrel:"<<forceCurveName<<":xml is empty"); 
	c.Stop();
	// ICMLOG("ARM_CalypsoToolkit::GetBasketCorrel:"<<forceCurveName<<":"<<c.GetDurationInMilliseconds()<<" ms"); 
}

//	------------------------------------------------------------------------------------
// 	static 
void ARM_CalypsoToolkit::GetTrade(const std::string& calypsoId,
								  const ARM_Date& date,
								  std::string& xmlOutput)
{
	ICMELAPSER("ARM_CalypsoToolkit::GetTrade:"<<calypsoId); 

	PreciseChrono c;
	c.Start(); 

	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetMonth()<<"/"
		<<setw(2)<<setfill('0')<<date.GetDay()<<"/"<<date.GetYear(); 

	_variant_t calypsoTrade; 
	_variant_t err ; 
	try {
		JavaExcelCom::get()->GetTradeCashFlow(calypsoId.c_str(),"MO",sstr.str().c_str(),0,&calypsoTrade,&err); 
	}
	catch(_com_error &e)
	{
		std::string desc (e.Description()); 
		std::string src (e.Source()); 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_CalypsoToolkit::GetTrade:"<<src<<":"<<desc); 
	}
	xmlOutput=_bstr_t(calypsoTrade); 
	if (xmlOutput.empty())
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_CalypsoToolkit::GetTrade:"<<calypsoId<<":xml is empty"); 
	c.Stop();
	
	MSXML2::IXMLDOMDocumentPtr xmlDoc; 
	xmlDoc = XMLTools::LoadXML(xmlOutput); 
	#ifdef _DEBUG
	_variant_t tmp("c:\\temp\\Calypso_Trade.xml"); 
	xmlDoc->save(tmp); 
	#endif

	// ICMLOG("ARM_CalypsoToolkit::GetTrade:"<<calypsoId<<":"<<c.GetDurationInMilliseconds()<<" ms"); 
}
//	------------------------------------------------------------------------------------
// static 
void 
ARM_CalypsoToolkit::GetFixing(const std::string& index,
		const std::string& term,
		const std::string& ccy,
		const std::string& page,
		const std::string& forceCurveName,
		const ARM_Date& date,
		double& output)
{
	ICMELAPSER("ARM_CalypsoToolkit::GetFixing:"<<index
		<<","<<term<<","<<ccy<<","<<page<<","<<forceCurveName<<","<<date<<" ...")
	PreciseChrono c; 
	c.Start();
	try 
	{
		_variant_t ret ; 
		JavaExcelCom::get()->MarketData_GetQuoteSet(forceCurveName.c_str(),&ret); 
		long quoteId = ret ;
		if (quoteId==0) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_CalypsoToolkit::ARM_CalypsoToolkitGetFixing: can't get QuoteSet "<<forceCurveName); 
		_variant_t vdate; 
		VariantTools::convert(date,vdate); 	
		std::vector<std::string> resRequest; 
		std::string arg = "MM." ;
		arg += ccy + "." + index + "." + term + "." + page ; 
		JavaExcelCom::get()->MarketData_GetQuoteByName(quoteId,arg.c_str(),vdate,&ret); 
		VariantTools::extractColumnFromRowOfRow(ret,1,resRequest); 
		if (resRequest.size() == 1){
			VariantTools::convertFromRowOfRow(ret,0,7,output) ;
			if (output != 0) {
				output*=100.;
				return;
			}
		}
		if ( resRequest.size() == 0 || output == 0)
			ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_CalypsoToolkit::ARM_CalypsoToolkitGetFixing: can't find  "<<arg); 
	}
	catch(_com_error &e)
	{
		std::string desc (e.Description()); 
		std::string src (e.Source()); 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_CalypsoToolkit::GetTrade:"<<src<<":"<<desc); 
	}
	c.Stop(); 
	// ICMLOG("ARM_CalypsoToolkit::GetFixing:"<<c.GetDurationInMilliseconds()<<" ms"); 
}
//	------------------------------------------------------------------------------------
// 	static 
void 
ARM_CalypsoToolkit::GetFXRate(const std::string& ccy1,
		const std::string& ccy2,
		const std::string& forceCurveName,
		const ARM_Date& date,
		double& output) 
{
	ICMELAPSER("ARM_CalypsoToolkit::GetFXRate: "<<ccy1<<","<<ccy2<<","<<forceCurveName
		<<","<<date); 
	PreciseChrono c; 
	c.Start();
	try 
	{
		// same ccy always returns 1 
		if (ccy1==ccy2) { output=1; return; } 

		_variant_t ret ; 
		JavaExcelCom::get()->MarketData_GetQuoteSet(forceCurveName.c_str(),&ret); 
		long quoteId = ret ;
		if (quoteId==0) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_CalypsoToolkit::GetFXRate: can't get QuoteSet "<<forceCurveName); 
		_variant_t vdate; 
		VariantTools::convert(date,vdate);
		std::vector<std::string> resRequest;
		std::string arg("");
		arg = "FX."+ccy1+"."+ccy2; 
		JavaExcelCom::get()->MarketData_GetQuoteByName(quoteId,arg.c_str(),vdate,&ret); 
		VariantTools::extractColumnFromRowOfRow(ret,1,resRequest); 
		if (resRequest.size() == 1){
			VariantTools::convertFromRowOfRow(ret,0,7,output) ;
			if (output != 0) return;
		}
		// currency 
		arg = "FX."+ccy2+"."+ccy1;
		JavaExcelCom::get()->MarketData_GetQuoteByName(quoteId,arg.c_str(),vdate,&ret);
		VariantTools::extractColumnFromRowOfRow(ret,1,resRequest); 
		if (resRequest.size() == 1){
			VariantTools::convertFromRowOfRow(ret,0,7,output) ;
			if (output != 0) {
				output = 1./output;
				return;
			}
		}
		// Look for USD/ccy1 in any direction
		arg = "FX.USD."+ccy1; 
		double FixCcy1USD=0.;
		JavaExcelCom::get()->MarketData_GetQuoteByName(quoteId,arg.c_str(),vdate,&ret); 
		VariantTools::extractColumnFromRowOfRow(ret,1,resRequest); 
		if (resRequest.size() == 1){
			VariantTools::convertFromRowOfRow(ret,0,7,FixCcy1USD) ;
			if(FixCcy1USD != 0)	FixCcy1USD=1./FixCcy1USD ;
		}
		if (resRequest.size() != 1 || FixCcy1USD == 0 ){
			arg = "FX."+ccy1+".USD"; 
			JavaExcelCom::get()->MarketData_GetQuoteByName(quoteId,arg.c_str(),vdate,&ret);
			VariantTools::extractColumnFromRowOfRow(ret,1,resRequest); 
			if (resRequest.size() == 1){
				VariantTools::convertFromRowOfRow(ret,0,7,FixCcy1USD) ;
				if(FixCcy1USD != 0)	FixCcy1USD=1./FixCcy1USD ;
			} else {
				ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_CalypsoToolkit::GetFXRate: not found USD/"<<ccy1);
			}
		}
		// Look for USD/ccy2 in any direction
		arg = "FX.USD."+ccy2;
		double FixUDSCcy2 =0.;
		JavaExcelCom::get()->MarketData_GetQuoteByName(quoteId,arg.c_str(),vdate,&ret);
		VariantTools::extractColumnFromRowOfRow(ret,1,resRequest); 
		if (resRequest.size() == 1)
			VariantTools::convertFromRowOfRow(ret,0,7,FixUDSCcy2) ;

		if (resRequest.size() != 1 || FixUDSCcy2 == 0 ){
			arg = "FX."+ccy2+".USD"; 
			JavaExcelCom::get()->MarketData_GetQuoteByName(quoteId,arg.c_str(),vdate,&ret);
			VariantTools::extractColumnFromRowOfRow(ret,1,resRequest); 
			if (resRequest.size() == 1){
				VariantTools::convertFromRowOfRow(ret,0,7,FixUDSCcy2) ;
				if(FixUDSCcy2 != 0)	FixUDSCcy2=1./FixUDSCcy2 ;
			} else {
				ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_CalypsoToolkit::GetFXRate: not found USD/"<<ccy2);
			}
		}

		output=FixCcy1USD*FixUDSCcy2;
		return ;
	}
	catch(_com_error &e)
	{
		std::string desc (e.Description()); 
		std::string src (e.Source()); 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_CalypsoToolkit::GetFXRate:"<<src<<":"<<desc);
		
	}
	c.Stop(); 
	ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_CalypsoToolkit::GetFXRate: can't get "<<ccy1<<" "<<ccy2); 
}
//---------------------------------------------------------------------------------------
void 
ARM_CalypsoToolkit::GetVolatilitySurface(const std::string&  volSurfName,
                                        const std::string& cvName, 
                                        const ARM_Date& date,
                                        std::string& xmlOutput){
ICMELAPSER("ARM_CalypsoToolkit::GetVolatilitySurface:"<<volSurfName); 

	PreciseChrono c;
	c.Start(); 

	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetMonth()<<"/"
		<<setw(2)<<setfill('0')<<date.GetDay()<<"/"<<date.GetYear(); 

	_variant_t calypsoVolatility; 
	_variant_t err ; 
	try {
		JavaExcelCom::get()->getVolatilitySurface(volSurfName.c_str(),cvName.c_str(),sstr.str().c_str(),&calypsoVolatility,&err); 
	}
	catch(_com_error &e)
	{
		std::string desc (e.Description()); 
		std::string src (e.Source()); 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_CalypsoToolkit::GetVolatilitySurface:"<<src<<":"<<desc); 
	}
	xmlOutput=_bstr_t(calypsoVolatility); 
	if (xmlOutput.empty())
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_CalypsoToolkit::GetVolatilitySurface:"<<volSurfName<<":xml is empty"); 
	c.Stop();
	
	IXMLDOMDocumentPtr xmlDoc; 
	xmlDoc = XMLTools::LoadXML(xmlOutput); 
	#ifdef _DEBUG
	_variant_t tmp("c:\\temp\\Calypso_Volatility.xml"); 
	xmlDoc->save(tmp); 
	#endif
}

//	------------------------------------------------------------------------------------
	ARM_ZeroLInterpol* ARMLOCAL_ParseXMLForCalypsoZC(const std::string& chaineXML,
												 const ARM_Date & aSdate,
												 const std::string& sCcy,
												 int interp)
{
	VARIANT_BOOL bOK;
	MSXML2::IXMLDOMDocumentPtr xmlDoc; 

	std::vector<double> pdRate ; 
	std::vector<double> pdMatu; 
	long interpMethodId;


		xmlDoc=XMLTools::CreateDOMDocument30(); 
									   
		_bstr_t b;
		VariantTools::convert(chaineXML,b); 
		xmlDoc->loadXML(b,&bOK); 

		#ifdef _DEBUG
		CCString chtmp = "c:\\temp\\Calypso_zerocurve.xml";
		VARIANT v_chtmp;
		CCString2VARIANT (chtmp, &v_chtmp);
		xmlDoc->save(v_chtmp);
		#endif


	long nbNodes;


		MSXML2::IXMLDOMNodeListPtr resultList; 
		
		resultList = XMLTools::selectNodes(xmlDoc,"/ZeroCurveList/CurveZero/ZCPoints/ZCPoint") ;
		{
			resultList->get_length(&nbNodes);

			if ( nbNodes == 0 )
				ICMTHROW(ERR_INVALID_ARGUMENT,"Invalid XML string for ZC"); 

			pdRate.resize( nbNodes) ; 
			pdMatu.resize( nbNodes) ; 

			for (long indexNode = 0 ; indexNode < nbNodes ; indexNode++)
			{

				MSXML2::IXMLDOMNodePtr listItem = XMLTools::get_item(resultList,indexNode); 


				if (  listItem != NULL )
				{
					ARM_Date endDate = GetDateFromXMLNode(listItem,"Date");
					pdMatu[indexNode] = (endDate - aSdate) / 365.;
					pdRate[indexNode] = GetDoubleFromXMLNode(listItem,"Value");

				}
			}
		}
		
        //
		
		MSXML2::IXMLDOMNodePtr theNode  = XMLTools::selectSingleNode(xmlDoc,"/ZeroCurveList/CurveZero/CurveInfo/Interpolator"); 
		// {
			if (theNode!=NULL)
			{
				std::string ff1; 
				XMLTools::convert(theNode,ff1); 
				if (ff1=="InterpolatorContinuousARM") 
					interpMethodId = K_CONTINUOUS;
				else
					interpMethodId = K_LINEAR;

			}



	ARM_ZeroLInterpol* newZcLin = NULL;


		newZcLin = new ARM_ZeroLInterpol(ARM_Date(aSdate), 
			&ARM_Vector(pdMatu), 
			&ARM_Vector(pdRate), 0, 0, interpMethodId);


		newZcLin->SetCurrencyUnit(&ARM_Currency(sCcy.c_str()));

	

		return newZcLin;
}
