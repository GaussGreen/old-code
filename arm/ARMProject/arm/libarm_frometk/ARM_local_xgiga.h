#ifndef ARM_LOCAL_XGIGA_H
#define ARM_LOCAL_XGIGA_H

#include <string>

class ARM_Date;
//	-----------------------------------------------------------------------------------
//
//		This is the API to retrieve XML content from GigaSpaces. 
//
//
//

//		
class ARM_XGigaToolKit 
{
public:
	static void init(const std::string &cmdline) ; 
	static void doGetZC(	
		const std::string& index,
		const std::string& ccy,
		const std::string& term,
		const std::string& pricingEnv,
		const ARM_Date& date,
		std::string& xmlOutput
		);
	static void doGetYLD(
		const std::string& index,
		const std::string& ccy,
		const std::string& term,
		const std::string& pricingEnv,
		const ARM_Date& date,
		std::string& xmlOutput
		);
	static void doGetFXCORREL(
		const std::string& index,
		const std::string& ccy,
		const std::string& ccy1,
		const std::string& ccy2,
		const std::string& term,
		const std::string& pricingEnv,
		const ARM_Date& date,
		std::string& xmlOutput
		);
	static void doGetRESET(
		const std::string& ccy,
		const std::string& index,
		const std::string& source,
		const std::string& term,
		std::string& xmlOutput
		);
	static void doGetSMILETenorList(
		const std::string& index,
		const std::string& ccy,
		const std::string& pricingEnv,
		const std::string& type,
		const ARM_Date& date,
		std::string& xmlOutput
		);
	static void doGetFX(
		const std::string& ccy1,
		const std::string& ccy2,
		const std::string& cvname,
		const ARM_Date& date,
		std::string& xmlOutput
		);
	static void doGetVOLTenorsList(
		const std::string& type,
		const std::string& index,
		const std::string& ccy,
		const std::string& cvname,
		const std::string& rqtype,
		const ARM_Date& date,
		std::string& xmlOutput
		);
	static void doGetVOLTenor(
		const std::string& type,
		const std::string& index,
		const std::string& ccy,
		const std::string& cvname,
		const std::string& rqtype,
		const std::string& tenor,
		const ARM_Date& date,
		std::string& xmlOutput
		);
	static void doGetSMILETenorStrikeList(
		const std::string& index,
		const std::string& ccy,
		const std::string& cvname,
		const std::string& type,
		const std::string& tenor,
		const ARM_Date& date,
		std::string& xmlOutput
		);
	static void doGetSMILETenorStrike(
		const std::string& index,
		const std::string& ccy,
		const std::string& cvname,
		const std::string& type,
		const std::string& tenor,
		const std::string& strike,
		const ARM_Date& date,
		std::string& xmlOutput
		);
	static void doGetFXVOL(
		const std::string& ccy1,
		const std::string& ccy2,
		const std::string& cvname,
		const ARM_Date& date,
		std::string& xmlOutput
		);
	static void doGetCORRELTenorsList(
		const std::string& index1,
		const std::string& ccy1,
		const std::string& index2,
		const std::string& ccy2,
		const std::string& cvname,
		const ARM_Date& date,
		std::string& xmlOutput
		);
	static void doGetCORRELTenor(
		const std::string& index1,
		const std::string& ccy1,
		const std::string& index2,
		const std::string& ccy2,
		const std::string& tenor,
		const std::string& cvname,
		const ARM_Date& date,
		std::string& xmlOutput
		);
	static void doGetMODELFACTOR(
		const std::string& ccy,
		const std::string& index,
		const std::string& model,
		const std::string& type,
		const std::string& factorname,
		const std::string& cvname,
		const ARM_Date& date,
		std::string& xmlOutput
		);
	static void doGetSmileFXTenorList(
		const std::string& ccy1,
		const std::string& ccy2,
		const std::string& callOrPut,
		const std::string& cvname,
		const ARM_Date& date,
		std::string& xmlOutput
		);
	static void doGetSmileFXTenor(
		const std::string& ccy1,
		const std::string& ccy2,
		const std::string& callOrPut,
		const std::string& tenor,
		const std::string& cvname,
		const ARM_Date& date,
		std::string& xmlOutput
		);
};

//	----------------------------------------------------------------------------------
//
//	




#endif // ARM_LOCAL_XGIGA_H