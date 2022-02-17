/**********************************************************************
*
*	Header %name:	CalypsoSwapLeg.h %
*	Instance:		1
*	Description:	
*	%created_by:	jpriaudel %
*	%date_created:	Mon Jul 09 15:06:56 2007 %
*
**********************************************************************/
#ifndef _CALYPSOSWAPLEG_H
#define _CALYPSOSWAPLEG_H


/* Everything else goes here */
#include <refvalue.h>
#include "swapleg.h"
#include <util\fromto.h>
#include <ARM\libarm_frometk\arm_local_parsexml_util.h> 
#include "XMLTools.h"


class ARM_Date;
class ARM_ReferenceValue;

#define FIX       0
#define FLOAT       1
#define EXOTIC       2

struct CalypsoLeg {
	int type;
    double nominal;
	double fixRate;
	double spread;
	ARM_Currency currency;
	int stubRule;
	int intRule;
	int fwdRule;
	long notionalExchangeFlag;
	int decompFlag;
	int freq;
	int decompFreq;
	int dayCount;
	string payCal;
	int payTiming;
	int payGap;
	int resetTiming;
	int resetGap;
	int resetFreq;
	string resetCal;
	int indexType; 
	string indexTerm;
	int indexDayCount;
	long flowCount;
	double leverage;
	vector<double> flowResetDatesVector;
	vector<double> flowStartDatesVector;
	vector<double> flowEndDatesVector;
	vector<double> flowNominalVector;
	vector<double> flowSpreadVector;
	vector<double> flowRateVector;

    vector<double> strikeVector;
    vector<double> leverageVector;
    vector<double> cpnMinVector;
    vector<double> cpnMaxVector;

    ARM_Date startDate;
    ARM_Date endDate;
    ARM_Date fixEndDate;
};

ARM_SwapLeg * convertFixLeg(MSXML2::IXMLDOMNodePtr xmlNode, ARM_Date startDate, ARM_Date endDate, int rcvPay, ARM_ReferenceValue*& resetRates);
ARM_SwapLeg * convertVarLeg(MSXML2::IXMLDOMNodePtr xmlNode, ARM_Date startDate, ARM_Date endDate, int rcvPay, ARM_ReferenceValue*& spreads,  ARM_ReferenceValue*& notionals);

CalypsoLeg getSwapLeg(MSXML2::IXMLDOMNodePtr xmlNode);
CalypsoLeg getCRFLeg(MSXML2::IXMLDOMNodePtr xmlNode);

#endif
