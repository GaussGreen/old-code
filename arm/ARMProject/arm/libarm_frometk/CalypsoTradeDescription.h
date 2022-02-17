/**********************************************************************
*
*	Header %name:	CalypsoTradeDescription.h %
*	Instance:		1
*	Description:	
*	%created_by:	hbelefquih %
*	%date_created:	Fri Mar 30 16:10:45 2007 %
*
**********************************************************************/

#ifndef _CALYPSOTRADEDESCRIPTION_H
#define _CALYPSOTRADEDESCRIPTION_H

#include <util\fromto.h>
#include <ARM\libarm_frometk\arm_local_parsexml_util.h> 
#include "XMLTools.h"


class ARM_Date;


struct CalypsoTradeDescription  
{

    ARM_Date startDate;
    ARM_Date endDate;
    int buySell;
    int payRec;
    string couponLeg;
    string fundingLeg;
    string pruductDefinition;

};

CalypsoTradeDescription getTradeDescription(IXMLDOMDocumentPtr xmlDoc);
#endif 
