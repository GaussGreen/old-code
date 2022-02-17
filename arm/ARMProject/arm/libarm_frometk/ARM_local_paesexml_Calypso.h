#ifndef ARM_LOCAL_PARSEXML_CALYPSO_H
#define ARM_LOCAL_PARSEXML_CALYPSO_H

#include <string>
#include <dates.h>
#include <zeroint.h>
//	-----------------------------------------------------------------------------------
//
//		This is the API to retrieve XML content from Calypso. 
//
//		All 
//
//			xmlFileName : optional string parameter ;
//			if not empty the file will be read and stored in xmlContent
//
//		MarketData
//
//			forceCurveName : optional string parameter
//			if not empty, the specific instance will be retrieved
//			otherwise, ARM will generate the curve name 
//
//		Note that there is no fallback from XML File to Calypso. 
//		
class ARM_CalypsoToolkit 
{
public:
	static void init(const std::string &cmdline) ; 
	static void GetProbabilityCurve(	
		const std::string& issuer,
		const std::string& ccy,
		const std::string& seniority,
		const std::string& forceCurveName,
		const std::string& pricingEnv,
		const ARM_Date& date,
		const std::string& xmlFileName, 
		std::string& xmlOutput
		); 
	static void GetCurveZero(
		const std::string& index,
		const std::string& ccy,
		const std::string& term,
		const std::string& forceCurveName,
		const std::string& pricingEnv,
		const ARM_Date& date,
		const std::string& xmlFileName,
		std::string& xmlOutput
		); 
	static void GetBasketCorrel(
		const std::string& forceCurveName,
		const std::string& pricingEnv,
		const ARM_Date& date,
		const std::string& xmlFileName,
		std::string& xmlOutput
		); 
	static void GetTrade(
		const std::string& idCalypso,
		const ARM_Date& date,
		std::string& xmlOutput
		); 
	static void GetFixing(const std::string& index,
		const std::string& term,
		const std::string& ccy,
		const std::string& page,
		const std::string& forceCurveName,
		const ARM_Date& date,
		double& output) ;
	static void GetFXRate(const std::string& ccy1,
		const std::string& ccy2,
		const std::string& forceCurveName,
		const ARM_Date& date,
		double& output) ;
    static void GetVolatilitySurface(const std::string&  volSurfName,
                                        const std::string& cvName, 
                                        const ARM_Date& date,
                                        std::string& xmlOutput);
private:
	static void readTextFile(const std::string& xmlFileName,std::string& xmlContent); 
} ;

//	----------------------------------------------------------------------------------
//
//	
ARM_ZeroLInterpol* ARMLOCAL_ParseXMLForCalypsoZC(const std::string& chaineXML,
												 const ARM_Date &aSdate,
												 const std::string& sCcy,
												 int interp);



#endif