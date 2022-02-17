/**********************************************************************
*
*	C++ %name:		CalypsoTradeDescription.cpp %
*	Instance:		1
*	Description:	
*	%created_by:	hbelefquih %
*	%date_created:	Fri Mar 30 11:04:12 2007 %
*
**********************************************************************/


#include "CalypsoTradeDescription.h"

CalypsoTradeDescription getTradeDescription(IXMLDOMDocumentPtr xmlDoc){
    CalypsoTradeDescription desc;
    
    // Start Date
	XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Description/StartDate"),desc.startDate,"YYYYMMDD");
    
    // End Date
	XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Description/EndDate"),desc.endDate,"YYYYMMDD");
     
    //BuySell
	double quantity; XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Flows/quantity"),quantity);
    if (quantity >=0 ) 
   		desc.buySell = -1;
    else 
  		desc.buySell = 1;
	desc.buySell= calypso2ARMPayOrRec(desc.buySell);


    // Receive or Pay ...
    string payLegType;XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Flows/payLeg/fixFloat"),payLegType);

    desc.fundingLeg= "payLeg";
    desc.couponLeg= "receiveLeg";
	long payRec=  K_RCV;
    if((payLegType == "Fix")||(payLegType == "Form")){
    	desc.payRec = K_PAY;
	   	desc.fundingLeg = "receiveLeg";
       	desc.couponLeg = "payLeg";
    }
    
    XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Description/ProductDefinition"),desc.pruductDefinition);
	if(desc.pruductDefinition !="Swaption")
        desc.payRec*=-desc.buySell;

    return desc;
}