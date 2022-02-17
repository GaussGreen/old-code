#include <ARM\libarm_local\firstToBeIncluded.h>


#include <ARM\libarm_frometk\arm_local_parsexml_common.h> 
#include <ARM\libarm_frometk\arm_local_parsexml_util.h> 
#include <ARM\libarm_frometk\arm_local_xgiga.h> 
#include "PaserManagerUtilities.h"

#include <fstream>
#include <libCCatl\CCatl.h>
#include <ARM\libarm_frometk\XMLTools.h>
#include <ARM\libarm_frometk\VariantTools.h>


#import <XGiga.exe>	// 
using namespace MSXML2;

// this is part of ICM Kernel, but should be located on ARM_Date 
extern std::ostream& operator<<(std::ostream&,const ARM_Date&); 

//	------------------------------------------------------------------------------------
//
//		This is the singleton class holding a XGiga.CGetData COM instance
//
class XGigaCGetData
{
public:	
	static XGiga::IGetDataPtr& get()
	{
		if (itsInstance.get()==0) 
		{
			itsInstance=std::auto_ptr<XGigaCGetData>(new XGigaCGetData) ;
		}
		// return itsInstance->addinptr(); 
		if (*(itsInstance->itsAddinPtr)) return *(itsInstance->itsAddinPtr) ; 
		throw std::runtime_error("XGigaCGetData: not valid ptr"); 
	}
	static void setEnv(const std::string& env) 
	{
		itsEnv=env; 
	}
public:
	// XGiga::IGetDataPtr& addinptr() { if !(*itsAddinPtr) throw std::runtime_error("XGigaCGetData: not valid ptr"); return *itsAddinPtr; }
private:
	XGigaCGetData()
	{
//		ICMLOG("XGigaCGetData::XGigaCGetData: Creating the COM instance"); 
		itsAddinPtr = new XGiga::IGetDataPtr; 
		itsAddinPtr->CreateInstance(__uuidof(XGiga::CGetData)); 
		if (*itsAddinPtr) (*itsAddinPtr)->GetInstance(itsEnv.c_str()); 
	}
	~XGigaCGetData()
	{
//		ICMLOG("XGigaCGetData::~XGigaCGetData"); 
		try 
		{ 
			delete itsAddinPtr ; 
		} 
		catch(...)
		{
		} 
	}
	XGigaCGetData(const XGigaCGetData&); //NA
	XGigaCGetData& operator=(const XGigaCGetData&); // NA 
private:
	static std::auto_ptr<XGigaCGetData> itsInstance; 
	static std::string itsEnv; 
	XGiga::IGetDataPtr * itsAddinPtr; 
	friend class std::auto_ptr<XGigaCGetData>; 
}; 
std::auto_ptr<XGigaCGetData> XGigaCGetData::itsInstance(0);
std::string XGigaCGetData::itsEnv("NotDefined");

//	------------------------------------------------------------------------------------
//	static
void 
ARM_XGigaToolKit::init(const std::string& env)
{
//	ICMLOG("ARM_XGigaToolKit::init:"<<env<<" ..."); 
//	PreciseChrono c; 
//	c.Start() ;
	try {
		/** if (XGigaCGetData::get())
			XGigaCGetData::get()->GetInstance(env.c_str()); 
			**/ 
		XGigaCGetData::setEnv(env); 
	}
	catch(_com_error&)
	{
//		std::string desc (e.Description()); 
//		std::string src (e.Source()); 
//		ICMTHROW(ERR_INVALID_ARGUMENT,src<<":"<<desc); 
	}
	catch(...)
	{
	}
//	c.Stop() ;
//	ICMLOG("ARM_XGigaToolKit::init:"<<env<<":"<<c.GetDurationInMilliseconds()<<" ms"); 
}
//	------------------------------------------------------------------------------------
//	static
void 
ARM_XGigaToolKit::doGetYLD(const std::string& index,
										const std::string& ccy,
										const std::string& term,
										const std::string& pricingEnv,
										const ARM_Date& date,
										std::string& xmlOutput) 
{
//	ICMLOG("ARM_XGigaToolKit::doGetYLD:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<" ..."); 
//	PreciseChrono c;	
//	c.Start(); 

	
	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetDay()
		<<setw(2)<<setfill('0')<<date.GetMonth()
		<<setw(2)<<setfill('0')<<date.GetYear()-2000 ; 
	_variant_t probaCurve; 
	
	try {
		// if (XGigaCGetData::get())
			XGigaCGetData::get()->doGetYLD(index.c_str(),ccy.c_str(),pricingEnv.c_str(),sstr.str().c_str(),&probaCurve); 
	}
	catch(_com_error&)
	{
//		std::string desc (e.Description()); 
//		std::string src (e.Source()); 
//		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_XGigaToolKit::doGetYLD:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date
//			<<":"<<src<<":"<<desc); 
	}
	catch(...)
	{
	}
	// xmlOutput=_bstr_t(probaCurve); 
	VariantTools::convert(probaCurve,xmlOutput); 
//	c.Stop(); 
//	ICMLOG("ARM_XGigaToolKit::doGetYLD:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<":"<<c.GetDurationInMilliseconds()<<" ms");  

}
//	------------------------------------------------------------------------------------
//	static
void 
ARM_XGigaToolKit::doGetZC(const std::string& index,
										const std::string& ccy,
										const std::string& term,
										const std::string& pricingEnv,
										const ARM_Date& date,
										std::string& xmlOutput) 
{
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<" ..."); 
//	PreciseChrono c;	
//	c.Start(); 

	
	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetDay()
		<<setw(2)<<setfill('0')<<date.GetMonth()
		<<setw(2)<<setfill('0')<<date.GetYear()-2000 ; 
	_variant_t probaCurve; 
	
	try {
		// if (XGigaCGetData::get())
			XGigaCGetData::get()->doGetZC(index.c_str(),ccy.c_str(),pricingEnv.c_str(),sstr.str().c_str(),&probaCurve); 
	}
	catch(_com_error&)
	{
//		std::string desc (e.Description()); 
//		std::string src (e.Source()); 
//		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date
//			<<":"<<src<<":"<<desc); 
	}
	catch(...)
	{
	}
	// xmlOutput=_bstr_t(probaCurve); 
	VariantTools::convert(probaCurve,xmlOutput); 
//	c.Stop(); 
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<":"<<c.GetDurationInMilliseconds()<<" ms");  

}
//	------------------------------------------------------------------------------------
//	static
void 
ARM_XGigaToolKit::doGetFXCORREL(const std::string& index,
										const std::string& ccy,
										const std::string& ccy1,
										const std::string& ccy2,
										const std::string& term,
										const std::string& pricingEnv,
										const ARM_Date& date,
										std::string& xmlOutput) 
{
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<" ..."); 
//	PreciseChrono c;	
//	c.Start(); 

	
	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetDay()
		<<setw(2)<<setfill('0')<<date.GetMonth()
		<<setw(2)<<setfill('0')<<date.GetYear()-2000 ; 
	_variant_t probaCurve; 
	
	try {
		//if (XGigaCGetData::get())
			XGigaCGetData::get()->doGetFXCORREL(index.c_str(),ccy.c_str(),ccy1.c_str(),ccy2.c_str(),term.c_str(),pricingEnv.c_str(),sstr.str().c_str(),&probaCurve); 
	}
	catch(_com_error&)
	{
//		std::string desc (e.Description()); 
//		std::string src (e.Source()); 
//		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_XGigaToolKit::doGetFXCORREL:"<<index<<","<<ccy<<","<<ccy1<<","<<ccy2<<","<<term<<","<<pricingEnv<<","<<date
//			<<":"<<src<<":"<<desc); 
	}
	catch(...)
	{
	}
	// xmlOutput=_bstr_t(probaCurve); 
	VariantTools::convert(probaCurve,xmlOutput); 
//	c.Stop(); 
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<":"<<c.GetDurationInMilliseconds()<<" ms");  

}
//	------------------------------------------------------------------------------------
//	static
void 
ARM_XGigaToolKit::doGetRESET(const std::string& ccy,
										const std::string& index,
										const std::string& source,
										const std::string& term,
										std::string& xmlOutput) 
{
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<" ..."); 
//	PreciseChrono c;	
//	c.Start(); 

	ARM_Date date;	
	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetDay()
		<<setw(2)<<setfill('0')<<date.GetMonth()
		<<setw(2)<<setfill('0')<<date.GetYear()-2000 ; 
	_variant_t probaCurve; 
	
	try {
		//if (XGigaCGetData::get())
			XGigaCGetData::get()->doGetRESET(ccy.c_str(),index.c_str(),source.c_str(),term.c_str(),sstr.str().c_str(),&probaCurve); 
	}
	catch(_com_error&)
	{
//		std::string desc (e.Description()); 
//		std::string src (e.Source()); 
//		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_XGigaToolKit::doGetRESET:"<<ccy<<","<<index<<","<<source<<","<<term<<","<<date
//			<<":"<<src<<":"<<desc); 
	}
	catch(...)
	{
	}
	// xmlOutput=_bstr_t(probaCurve); 
	VariantTools::convert(probaCurve,xmlOutput); 

//	c.Stop(); 
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<":"<<c.GetDurationInMilliseconds()<<" ms");  

}

//	------------------------------------------------------------------------------------
//	static
void 
ARM_XGigaToolKit::doGetSMILETenorList(const std::string& index,
										const std::string& ccy,
										const std::string& cvname,
										const std::string& type,
										const ARM_Date& date,
										std::string& xmlOutput) 
{
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<" ..."); 
//	PreciseChrono c;	
//	c.Start(); 

	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetDay()
		<<setw(2)<<setfill('0')<<date.GetMonth()
		<<setw(2)<<setfill('0')<<date.GetYear()-2000 ; 
	_variant_t probaCurve; 
	
	try {
		//if (XGigaCGetData::get())
			XGigaCGetData::get()->doGetSMILETenorsList(index.c_str(),ccy.c_str(),cvname.c_str(),type.c_str(),sstr.str().c_str(),&probaCurve); 
	}
	catch(_com_error&)
	{
//		std::string desc (e.Description()); 
//		std::string src (e.Source()); 
//		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_XGigaToolKit::doGetSMILETenorList:"<<index<<","<<ccy<<","<<cvname<<","<<type<<","<<date
//			<<":"<<src<<":"<<desc);
	}
	catch(...)
	{
	}
	// xmlOutput=_bstr_t(probaCurve); 
	VariantTools::convert(probaCurve,xmlOutput); 

//	c.Stop(); 
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<":"<<c.GetDurationInMilliseconds()<<" ms");  
}

//	------------------------------------------------------------------------------------
//	static
void 
ARM_XGigaToolKit::doGetFX(const std::string& ccy1,
										const std::string& ccy2,
										const std::string& cvname,
										const ARM_Date& date,
										std::string& xmlOutput) 
{
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<" ..."); 
//	PreciseChrono c;	
//	c.Start(); 

	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetDay()
		<<setw(2)<<setfill('0')<<date.GetMonth()
		<<setw(2)<<setfill('0')<<date.GetYear()-2000 ; 
	_variant_t probaCurve; 
	
	try {
		//if (XGigaCGetData::get())
			XGigaCGetData::get()->doGetFX(ccy1.c_str(),ccy2.c_str(),cvname.c_str(),sstr.str().c_str(),&probaCurve); 
	}
	catch(_com_error&)
	{
//		std::string desc (e.Description()); 
//		std::string src (e.Source()); 
//		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_XGigaToolKit::doGetFX:"<<ccy1<<","<<ccy2<<","<<cvname<<","<<date
//			<<":"<<src<<":"<<desc);
	}
	catch(...)
	{
	}
	// xmlOutput=_bstr_t(probaCurve); 
	VariantTools::convert(probaCurve,xmlOutput); 

//	c.Stop(); 
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<":"<<c.GetDurationInMilliseconds()<<" ms");  
}

//	------------------------------------------------------------------------------------
//	static
void 
ARM_XGigaToolKit::doGetVOLTenorsList(const std::string& type,
									 const std::string& index,
									 const std::string& ccy,
									 const std::string& cvname,
									 const std::string& rqtype,
									 const ARM_Date& date,
									 std::string& xmlOutput) 
{
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<" ..."); 
//	PreciseChrono c;	
//	c.Start(); 

	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetDay()
		<<setw(2)<<setfill('0')<<date.GetMonth()
		<<setw(2)<<setfill('0')<<date.GetYear()-2000 ; 
	_variant_t probaCurve; 
	
	try {
		//if (XGigaCGetData::get())
			XGigaCGetData::get()->doGetVOLTenorsList(type.c_str(),index.c_str(),ccy.c_str(),cvname.c_str(),rqtype.c_str(),sstr.str().c_str(),&probaCurve); 
	}
	catch(_com_error&)
	{
//		std::string desc (e.Description()); 
//		std::string src (e.Source()); 
//		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_XGigaToolKit::doGetVOLTenorsList:"<<type<<","<<index<<","<<ccy<<","<<cvname<<","<<rqtype<<","<<date
//			<<":"<<src<<":"<<desc);
	}
	catch(...)
	{
	}
	// xmlOutput=_bstr_t(probaCurve); 
	VariantTools::convert(probaCurve,xmlOutput); 

//	c.Stop(); 
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<":"<<c.GetDurationInMilliseconds()<<" ms");  
}

//	------------------------------------------------------------------------------------
//	static
void 
ARM_XGigaToolKit::doGetVOLTenor(const std::string& type,
								const std::string& index,
								const std::string& ccy,
								const std::string& cvname,
								const std::string& rqtype,
								const std::string& tenor,
								const ARM_Date& date,
								std::string& xmlOutput) 
{
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<" ..."); 
//	PreciseChrono c;	
//	c.Start(); 

	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetDay()
		<<setw(2)<<setfill('0')<<date.GetMonth()
		<<setw(2)<<setfill('0')<<date.GetYear()-2000 ; 
	_variant_t probaCurve; 
	
	try {
		//if (XGigaCGetData::get())
			XGigaCGetData::get()->doGetVOLTenor(type.c_str(),index.c_str(),ccy.c_str(),cvname.c_str(),rqtype.c_str(),tenor.c_str(),sstr.str().c_str(),&probaCurve); 
	}
	catch(_com_error&)
	{
//		std::string desc (e.Description()); 
//		std::string src (e.Source()); 
//		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_XGigaToolKit::doGetVOLTenor:"<<type<<","<<index<<","<<ccy<<","<<cvname<<","<<rqtype<<","<<tenor<<","<<date
//			<<":"<<src<<":"<<desc);
	}
	catch(...)
	{
	}
	// xmlOutput=_bstr_t(probaCurve); 
	VariantTools::convert(probaCurve,xmlOutput); 
//	c.Stop(); 
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<":"<<c.GetDurationInMilliseconds()<<" ms");  
}

//	------------------------------------------------------------------------------------
//	static
void 
ARM_XGigaToolKit::doGetSMILETenorStrikeList(
								const std::string& index,
								const std::string& ccy,
								const std::string& cvname,
								const std::string& type,
								const std::string& tenor,
								const ARM_Date& date,
								std::string& xmlOutput) 
{
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<" ..."); 
//	PreciseChrono c;	
//	c.Start(); 

	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetDay()
		<<setw(2)<<setfill('0')<<date.GetMonth()
		<<setw(2)<<setfill('0')<<date.GetYear()-2000 ; 
	_variant_t probaCurve; 
	
	try {
		//if (XGigaCGetData::get())
			XGigaCGetData::get()->doGetSMILETenorStrikeList(index.c_str(),ccy.c_str(),cvname.c_str(),type.c_str(),tenor.c_str(),sstr.str().c_str(),&probaCurve); 
	}
	catch(_com_error&)
	{
//		std::string desc (e.Description()); 
//		std::string src (e.Source()); 
//		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_XGigaToolKit::doGetSMILETenorStrikeList:"<<index<<","<<ccy<<","<<cvname<<","<<type<<","<<tenor<<","<<date
//			<<":"<<src<<":"<<desc);
	}
	catch(...)
	{
	}
	VariantTools::convert(probaCurve,xmlOutput); 
	// xmlOutput=_bstr_t(probaCurve); 

//	c.Stop(); 
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<":"<<c.GetDurationInMilliseconds()<<" ms");  
}

//	------------------------------------------------------------------------------------
//	static
void 
ARM_XGigaToolKit::doGetSMILETenorStrike(
								const std::string& index,
								const std::string& ccy,
								const std::string& cvname,
								const std::string& type,
								const std::string& tenor,
								const std::string& strike,
								const ARM_Date& date,
								std::string& xmlOutput) 
{
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<" ..."); 
//	PreciseChrono c;	
//	c.Start(); 

	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetDay()
		<<setw(2)<<setfill('0')<<date.GetMonth()
		<<setw(2)<<setfill('0')<<date.GetYear()-2000 ; 
	_variant_t probaCurve; 
	
	try {
		//if (XGigaCGetData::get())
			XGigaCGetData::get()->doGetSMILETenorStrike(index.c_str(),ccy.c_str(),cvname.c_str(),type.c_str(),tenor.c_str(),strike.c_str(),sstr.str().c_str(),&probaCurve); 
	}
	catch(_com_error&)
	{
//		std::string desc (e.Description()); 
//		std::string src (e.Source()); 
//		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_XGigaToolKit::doGetSMILETenorStrike:"<<index<<","<<ccy<<","<<cvname<<","<<type<<","<<tenor<<","<<strike<<","<<date
//			<<":"<<src<<":"<<desc);
	}
	catch(...)
	{
	}
	// xmlOutput=_bstr_t(probaCurve); 
	VariantTools::convert(probaCurve,xmlOutput); 

//	c.Stop(); 
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<":"<<c.GetDurationInMilliseconds()<<" ms");  
}

//	------------------------------------------------------------------------------------
//	static
void 
ARM_XGigaToolKit::doGetFXVOL(
								const std::string& ccy1,
								const std::string& ccy2,
								const std::string& cvname,
								const ARM_Date& date,
								std::string& xmlOutput) 
{
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<" ..."); 
//	PreciseChrono c;	
//	c.Start(); 

	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetDay()
		<<setw(2)<<setfill('0')<<date.GetMonth()
		<<setw(2)<<setfill('0')<<date.GetYear()-2000 ; 
	_variant_t probaCurve; 
	
	try {
		//if (XGigaCGetData::get())
			XGigaCGetData::get()->doGetFXVOL(ccy1.c_str(),ccy2.c_str(),cvname.c_str(),sstr.str().c_str(),&probaCurve); 
	}
	catch(_com_error&)
	{
//		std::string desc (e.Description()); 
//		std::string src (e.Source()); 
//		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_XGigaToolKit::doGetFXVOL:"<<ccy1<<","<<ccy2<<","<<cvname<<","<<date
//			<<":"<<src<<":"<<desc);
	}
	catch(...)
	{
	}
	// xmlOutput=_bstr_t(probaCurve); 
	VariantTools::convert(probaCurve,xmlOutput); 

//	c.Stop(); 
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<":"<<c.GetDurationInMilliseconds()<<" ms");  
}

//	------------------------------------------------------------------------------------
//	static
void 
ARM_XGigaToolKit::doGetCORRELTenorsList(
								const std::string& index1,
								const std::string& ccy1,
								const std::string& index2,
								const std::string& ccy2,
								const std::string& cvname,
								const ARM_Date& date,
								std::string& xmlOutput) 
{
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<" ..."); 
//	PreciseChrono c;	
//	c.Start(); 

	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetDay()
		<<setw(2)<<setfill('0')<<date.GetMonth()
		<<setw(2)<<setfill('0')<<date.GetYear()-2000 ; 
	_variant_t probaCurve; 
	
	try {
		//if (XGigaCGetData::get())
			XGigaCGetData::get()->doGetCORRELTenorsList(index1.c_str(),ccy1.c_str(),index2.c_str(),ccy2.c_str(),cvname.c_str(),sstr.str().c_str(),&probaCurve); 
	}
	catch(_com_error&)
	{
//		std::string desc (e.Description()); 
//		std::string src (e.Source()); 
//		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_XGigaToolKit::doGetCORRELTenorsList:"<<index1<<","<<ccy1<<","<<index2<<","<<ccy2<<","<<cvname<<","<<date
//			<<":"<<src<<":"<<desc);
	}
	catch(...)
	{
	}
	// xmlOutput=_bstr_t(probaCurve); 
	VariantTools::convert(probaCurve,xmlOutput); 

//	c.Stop(); 
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<":"<<c.GetDurationInMilliseconds()<<" ms");  
}

//	------------------------------------------------------------------------------------
//	static
void 
ARM_XGigaToolKit::doGetCORRELTenor(
								const std::string& index1,
								const std::string& ccy1,
								const std::string& index2,
								const std::string& ccy2,
								const std::string& tenor,
								const std::string& cvname,
								const ARM_Date& date,
								std::string& xmlOutput) 
{
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<" ..."); 
//	PreciseChrono c;	
//	c.Start(); 

	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetDay()
		<<setw(2)<<setfill('0')<<date.GetMonth()
		<<setw(2)<<setfill('0')<<date.GetYear()-2000 ; 
	_variant_t probaCurve; 
	
	try {
		//if (XGigaCGetData::get())
			XGigaCGetData::get()->doGetCORRELTenor(index1.c_str(),ccy1.c_str(),index2.c_str(),ccy2.c_str(),tenor.c_str(),cvname.c_str(),sstr.str().c_str(),&probaCurve); 
	}
	catch(_com_error&)
	{
//		std::string desc (e.Description()); 
//		std::string src (e.Source()); 
//		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_XGigaToolKit::doGetCORRELTenor:"<<index1<<","<<ccy1<<","<<index2<<","<<ccy2<<","<<tenor<<","<<cvname<<","<<date
//			<<":"<<src<<":"<<desc);
	}
	catch(...)
	{
	}
	// xmlOutput=_bstr_t(probaCurve); 
	VariantTools::convert(probaCurve,xmlOutput); 

//	c.Stop(); 
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<":"<<c.GetDurationInMilliseconds()<<" ms");  
}

//	------------------------------------------------------------------------------------
//	static
void 
ARM_XGigaToolKit::doGetMODELFACTOR(
		const std::string& ccy,
		const std::string& index,
		const std::string& model,
		const std::string& type,
		const std::string& factorname,
		const std::string& cvname,
		const ARM_Date& date,
		std::string& xmlOutput
		)
{
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<" ..."); 
//	PreciseChrono c;	
//	c.Start(); 

	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetDay()
		<<setw(2)<<setfill('0')<<date.GetMonth()
		<<setw(2)<<setfill('0')<<date.GetYear()-2000 ; 
	_variant_t probaCurve; 
	
	try {
		//if (XGigaCGetData::get())
			XGigaCGetData::get()->doGetMODELFACTOR(ccy.c_str(),index.c_str(),model.c_str(),type.c_str(),factorname.c_str(),cvname.c_str(),sstr.str().c_str(),&probaCurve); 
	}
	catch(_com_error&)
	{
//		std::string desc (e.Description()); 
//		std::string src (e.Source()); 
//		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_XGigaToolKit::doGetCORRELTenor:"<<index1<<","<<ccy1<<","<<index2<<","<<ccy2<<","<<tenor<<","<<cvname<<","<<date
//			<<":"<<src<<":"<<desc);
	}
	catch(...)
	{
	}
	// xmlOutput=_bstr_t(probaCurve); 
	VariantTools::convert(probaCurve,xmlOutput); 

//	c.Stop(); 
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<":"<<c.GetDurationInMilliseconds()<<" ms");  
}


//	------------------------------------------------------------------------------------
//	static
void 
ARM_XGigaToolKit::doGetSmileFXTenorList(
		const std::string& ccy1,
		const std::string& ccy2,
		const std::string& callOrPut,
		const std::string& cvname,
		const ARM_Date& date,
		std::string& xmlOutput
		)
{
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<" ..."); 
//	PreciseChrono c;	
//	c.Start(); 

	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetDay()
		<<setw(2)<<setfill('0')<<date.GetMonth()
		<<setw(2)<<setfill('0')<<date.GetYear()-2000 ; 
	_variant_t probaCurve; 
	
	try {
		//if (XGigaCGetData::get())
			XGigaCGetData::get()->doGetSMILEFXTenorsList(ccy1.c_str(),ccy2.c_str(),cvname.c_str(),callOrPut.c_str(),sstr.str().c_str(),&probaCurve); 
	}
	catch(_com_error&)
	{
//		std::string desc (e.Description()); 
//		std::string src (e.Source()); 
//		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_XGigaToolKit::doGetCORRELTenor:"<<index1<<","<<ccy1<<","<<index2<<","<<ccy2<<","<<tenor<<","<<cvname<<","<<date
//			<<":"<<src<<":"<<desc);
	}
	catch(...)
	{
	}
	// xmlOutput=_bstr_t(probaCurve); 
	VariantTools::convert(probaCurve,xmlOutput); 
//	c.Stop(); 
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<":"<<c.GetDurationInMilliseconds()<<" ms");  
}


//	------------------------------------------------------------------------------------
//	static
void 
ARM_XGigaToolKit::doGetSmileFXTenor(
		const std::string& ccy1,
		const std::string& ccy2,
		const std::string& callOrPut,
		const std::string& tenor,
		const std::string& cvname,
		const ARM_Date& date,
		std::string& xmlOutput
		)
{
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<" ..."); 
//	PreciseChrono c;	
//	c.Start(); 

	std::stringstream sstr; 
	sstr<<setw(2)<<setfill('0')<<date.GetDay()
		<<setw(2)<<setfill('0')<<date.GetMonth()
		<<setw(2)<<setfill('0')<<date.GetYear()-2000 ; 
	_variant_t probaCurve; 
	
	try {
		//if (XGigaCGetData::get())
			XGigaCGetData::get()->doGetSMILEFXTenor(ccy1.c_str(),ccy2.c_str(),cvname.c_str(),callOrPut.c_str(),tenor.c_str(),sstr.str().c_str(),&probaCurve); 
	}
	catch(_com_error&)
	{
//		std::string desc (e.Description()); 
//		std::string src (e.Source()); 
//		ICMTHROW(ERR_INVALID_ARGUMENT,"ARM_XGigaToolKit::doGetCORRELTenor:"<<index1<<","<<ccy1<<","<<index2<<","<<ccy2<<","<<tenor<<","<<cvname<<","<<date
//			<<":"<<src<<":"<<desc);
	}
	catch(...)
	{
	}
	// xmlOutput=_bstr_t(probaCurve); 
	VariantTools::convert(probaCurve,xmlOutput); 

//	c.Stop(); 
//	ICMLOG("ARM_XGigaToolKit::doGetZC:"<<index<<","<<ccy<<","<<term<<","<<pricingEnv<<","<<date<<":"<<c.GetDurationInMilliseconds()<<" ms");  
}
