/**********************************************************************
*
*	C++ %name:		CalypsoSwapLeg.cpp %
*	Instance:		1
*	Description:	
*	%created_by:	jpriaudel %
*	%date_created:	Mon Jul 09 15:04:33 2007 %
*
**********************************************************************/



#include "CalypsoSwapLeg.h"


ARM_SwapLeg * convertFixLeg(MSXML2::IXMLDOMNodePtr xmlNode, ARM_Date startDate, ARM_Date endDate, int rcvPay, ARM_ReferenceValue*& resetRates) {
	ARM_SwapLeg * fixLeg = NULL;

	CalypsoLeg leg = getSwapLeg(xmlNode);
	char * payCal = new char[ leg.payCal.size() + 1 ];
	strcpy(payCal, leg.payCal.c_str());

	fixLeg = new ARM_SwapLeg(startDate, 
							 endDate, 
							 leg.fixRate, 
							 rcvPay,
							 leg.freq, 
							 leg.dayCount, 
							 leg.decompFreq,
							 leg.payTiming,
							 leg.intRule,
							 leg.stubRule,
							 &leg.currency,
							 payCal,
							 leg.notionalExchangeFlag,
							 NULL, //ref date
							 1);   //adjust start

	// Rates
	resetRates = createRefValue(leg.flowResetDatesVector, leg.flowRateVector, K_STEPUP_LEFT, 1);

	if(payCal) {
		delete[] payCal;
		payCal = NULL;
	}
	return fixLeg;
}

ARM_SwapLeg * convertVarLeg(MSXML2::IXMLDOMNodePtr xmlNode, ARM_Date startDate, ARM_Date endDate, int rcvPay, ARM_ReferenceValue*& spreads,  ARM_ReferenceValue*& notionals) {
	ARM_SwapLeg * varLeg = NULL;

	CalypsoLeg leg = getSwapLeg(xmlNode);
	char * payCal = new char[ leg.payCal.size() + 1 ];
	strcpy(payCal, leg.payCal.c_str());
	char * resetCal = new char[ leg.resetCal.size() + 1 ];
	strcpy(resetCal, leg.resetCal.c_str());

	varLeg = new ARM_SwapLeg(startDate,
							 endDate,
							 (ARM_INDEX_TYPE)leg.indexType,
							 rcvPay,
							 leg.spread, 
							 leg.resetFreq,
							 leg.freq,
							 leg.resetTiming,
							 leg.payTiming,
							 &leg.currency,
							 leg.intRule,
							 leg.resetGap,
							 resetCal,
							 payCal,
							 leg.decompFlag,
							 leg.notionalExchangeFlag,
							 leg.stubRule,
							 NULL, //ref date
							 1,    //adjust start
							 leg.dayCount,
							 leg.fwdRule,
							 leg.payGap);   

	// Spreads
	vector<double> spreadsPct;
	for (int i = 0; i< leg.flowCount; i++) {
		spreadsPct.push_back(100 * leg.flowSpreadVector[i]);
	}

	spreads = createRefValue(leg.flowResetDatesVector, spreadsPct, K_STEPUP_LEFT, 1, true, true);
	// Notionals
	notionals = createRefValue(leg.flowEndDatesVector, leg.flowNominalVector, K_STEPUP_RIGHT, 1);

	if(payCal) {
		delete[] payCal;
		payCal = NULL;
	}
	if(resetCal) {
		delete[] resetCal;
		resetCal = NULL;
	}
	return varLeg;
}


CalypsoLeg getSwapLeg(MSXML2::IXMLDOMNodePtr xmlNode) {

	CalypsoLeg leg;

    //Type
    string stype;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"fixFloat"),stype);
    leg.type=(stype=="Form")?EXOTIC:(stype=="Float")?FLOAT:FIX;


    //startDate
    XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"startDate"),leg.startDate,"YYYYMMDD");


    //endDate
    XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"endDate"),leg.endDate,"YYYYMMDD");


    //fixEndDate
    leg.fixEndDate=leg.endDate;

	// Fix rate
	try{
        XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"fixedRate"),leg.fixRate );
        leg.fixRate *= 100;
    }catch(Exception e){
        leg.fixRate= 0;
    }

	// Spread
	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"spread"),leg.spread);
	
    // Nominal
	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"notional/quantity"),leg.nominal);
	
	// Currency
	string ccy;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"notional/currency"),ccy);
	leg.currency = ARM_Currency(ccy.c_str());
	
	// Payment Frequency
	string sfreq;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"couponFrequency"),sfreq);
    leg.freq = calypso2ARMFreq(sfreq.c_str());


    // Decomp flag
    string cmp;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"Cmp"),cmp);
    leg.decompFlag =(cmp=="true")?1:0;
	
    // TODO remove special conversion
    // Decomp Frequency
    leg.decompFreq = leg.freq;
    if(leg.decompFlag == 1){
        string sCompoundFrequency;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"CompoundFrequency"),sCompoundFrequency);
    	if(sCompoundFrequency   == "Semi-annually")
            sCompoundFrequency="SA";
        if(sCompoundFrequency   == "Per annum")
            sCompoundFrequency="PA";
        leg.decompFreq =calypso2ARMFreq(sCompoundFrequency.c_str());
    }
    if(leg.decompFreq==K_DEF_FREQ){
        leg.decompFreq = leg.freq;
    }
    

	// Day count
	string sdayCount;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"dayCount"),sdayCount);
    leg.dayCount = calypso2ARMDaycount(sdayCount.c_str());

	// Payment calendar
	string spayCal;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"couponDateAdjustement/holidays"),spayCal);
	leg.payCal = calypso2ARMCalendar(spayCal);

	// Payment timing
	//string spayTiming;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"resetTiming2"),spayTiming);
	string spayTiming;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"PmtResetTiming"),spayTiming);
	leg.payTiming = calypso2ARMTiming(spayTiming.c_str());


	// Stub Rule
	string sstubRule;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"couponStubRule"),sstubRule);
	leg.stubRule =  calypso2ARMStubRule(sstubRule.c_str());


	// Payment int rule
	string sintRule;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"periodRule"),sintRule);
	leg.intRule =  calypso2ARMInterestRule(sintRule.c_str());



	// Notional exchange flag
	bool actual =false;
	bool initEx=false;
	bool finalEx=false;
	try {
		string sactual;	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"Actual"),sactual);
		actual= (sactual=="true")?true:false;
		//TODO add Condition XCCySwap instatnce to add in argument
        if(actual) {
            try{
              //  string sinitEx;	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"PrincipalExchange/Initial"),sinitEx);
		        string sinitEx;	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"PrincipalExchangeType/Initial"),sinitEx);
		        initEx = (sinitEx=="true")?true:false;
		    //    string sfinalEx;	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"PrincipalExchange/Final"),sfinalEx);
		        string sfinalEx;	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"PrincipalExchangeType/Final"),sfinalEx);
		        finalEx = (sfinalEx=="true")?true:false;
            }catch(Exception e){
                initEx = false;
	            finalEx = false;
	         }
        }	
       }catch (Exception e){
	    actual = false;
	    initEx = false;
	    finalEx = false;
	   }
    string snotionalEx = calypso2ARMNotionalExchange(initEx, finalEx);
	leg.notionalExchangeFlag = ARM_NotionalExchange(snotionalEx .c_str());


    // Index term
	try{
        XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"rateIndex/tenor"),leg.indexTerm);
    }catch(Exception e){
        leg.indexTerm = "1Y";
    }
    
    // Index type
	try{
        string sindexName;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"rateIndex/indexName"),sindexName);
	    leg.indexType = calypso2ARMIndexType(sindexName,leg.indexTerm);
    }catch(Exception e){}
	

	
    

	// Reset timing
    //string sresetTiming;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"resetTiming1"),sresetTiming);
    string sresetTiming;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"FixingResetTiming"),sresetTiming);
	leg.resetTiming = calypso2ARMTiming(sresetTiming.c_str());


    	// Reset Frequency  // TODO bad path to chamge
	//mid = jni->GetMethodID(clsLeg, "getResetFrequency", "()Ljava/lang/String;");
	string sfundResetFreq;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"couponFrequency"),sfundResetFreq);
	leg.resetFreq = calypso2ARMFreq(sfundResetFreq.c_str());
    
    

    // Forward Rule  
	string sfwdRule;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"couponDateAdjustement/dateRoll"),sfwdRule);
	leg.fwdRule = calypso2ARMDateRoll(sfwdRule);
	
    // Leverage  //TODO Find correct path
	//mid = jni->GetMethodID(clsLeg, "getLeverage", "()D");
	leg.leverage = 1;


	// Payment gap
	//mid = jni->GetMethodID(clsLeg, "getPayGap", "()I");
	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"couponDateAdjustement/days"),leg.payGap);

	string index;
	// Cash flows
    MSXML2::IXMLDOMNodeListPtr xmlNodeList= XMLTools::selectNodes(xmlNode,"Listecashflow/cashflow"); 
	xmlNodeList->get_length(&leg.flowCount);
   	ARM_Date flowResetDate,flowStartDate,flowEndDate;
    double flowNominalValue,flowSpreadValue,flowRateValue;
    int dateInit=0;
    leg.resetGap =0;
    for (int i = 0; i< leg.flowCount; i++) {
        string stype;XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"type"),stype);
        if(stype == "INTEREST"){
            XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"startDate"),flowStartDate,"YYYYMMDD");
            
            try{
                XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"resetDate"),flowResetDate,"YYYYMMDD");
                leg.flowResetDatesVector.push_back(flowResetDate.GetJulian());
            }catch(Exception e){
                //TODO finalize version
                ARM_Date* tmpDate =  new ARM_Date(flowStartDate.GetJulian());

                int gap =(leg.resetGap ==0)?-2:leg.resetGap ;
                
                if(gap >0)
                    tmpDate->NextBusinessDay(ABS(gap),leg.currency.GetCcyName());
                 else
                    tmpDate->PreviousBusinessDay(ABS(gap), (char*)(leg.payCal.c_str()));

                leg.flowResetDatesVector.push_back(tmpDate->GetJulian());
                if(tmpDate) {
		             delete tmpDate;
	        	     tmpDate = NULL;
                }
            }

            XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"endDate"),flowEndDate,"YYYYMMDD");
            XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"notional"),flowNominalValue);
            XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"spread"),flowSpreadValue);
            XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"rate"),flowRateValue);

            leg.flowStartDatesVector.push_back(flowStartDate.GetJulian());
		    leg.flowEndDatesVector.push_back(flowEndDate.GetJulian());
		    leg.flowNominalVector.push_back(ABS(flowNominalValue));
		    leg.flowSpreadVector.push_back(flowSpreadValue);
            leg.flowRateVector.push_back(100 * flowRateValue);

/*            if(leg.type==EXOTIC){
                long nbFormulaParameters=0;
                MSXML2::IXMLDOMNodeListPtr	xmlNodeList2 = XMLTools::selectNodes(XMLTools::get_item(xmlNodeList,i),"ListFormulaParameters/FormulaParameter"); 
				xmlNodeList2->get_length(&nbFormulaParameters); 
				string formulaParameterName;
				
                double dlevrage, dcpnMin, dcpnMax,dtstrike;
                
				for(int j=0;j<nbFormulaParameters;j++) {
					XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList2,j),"Parameter_Name"),formulaParameterName);
					
					if(formulaParameterName =="aLower")
						XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList2,j),"Parameter_Value"),dcpnMin);
					if(formulaParameterName == "aCoef")
						XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList2,j),"Parameter_Value"),dlevrage);
					if(formulaParameterName == "aUpper")
						XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList2,j),"Parameter_Value"),dcpnMax);
					if(formulaParameterName == "aFixedRate")
						XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList2,j),"Parameter_Value"),dtstrike);
					if(formulaParameterName == "aFloatingRate")
						XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList2,j),"Parameter_Value"),index);
				
				
				}
                if(dlevrage==0){
					dlevrage=1;
					index = "";
				}
	            
                if (dcpnMax == 0) dcpnMax= 100;

                if (!(index == "") && (leg.fixEndDate==leg.endDate)&&(dateInit==0)) {
        				leg.fixEndDate = flowStartDate;
						dateInit=1;
        		}

                leg.strikeVector.push_back(dtstrike);
				leg.leverageVector.push_back(dlevrage);
				leg.cpnMinVector.push_back(dcpnMin);
				leg.cpnMaxVector.push_back(dcpnMax);
            }else {
                //leg.strikeVector.push_back(100 * flowRateValue);
				leg.strikeVector.push_back(flowRateValue);
				//leg.leverageVector.push_back(0);
				//leg.cpnMinVector.push_back(0);
				//leg.cpnMaxVector.push_back(10000);
            }
*/        }
    }
    leg.flowCount=leg.flowStartDatesVector.size();
	
    
    if(leg.type==EXOTIC){
        leg.indexDayCount = leg.dayCount;
	
		MSXML2::IXMLDOMNodeListPtr xmlNodeList3 = XMLTools::selectNodes(xmlNode,"ListeExoticVariables/ExoticVariable");
        if(xmlNodeList3==NULL)
                xmlNodeList3 = XMLTools::selectNodes(xmlNode,"ListeExoticVariable/ExoticVariable"); 
        try{
            string sname;XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList3,0),"variableName"),sname);
        }catch(Exception e){
            xmlNodeList3 = XMLTools::selectNodes(xmlNode,"ListeExoticVariable/ExoticVariable"); 
        }
		
		long exoticVariableCount;
		xmlNodeList3->get_length(&exoticVariableCount);

		string sname,sindexDaycount,sresetCal;
		for( i=0;i<exoticVariableCount;i++) 
		{
			XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList3,i),"variableName"),sname);
            
			if(sname==index){
				// Reset Gap
                XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList3,i),"resetDays"),leg.resetGap);
				leg.resetGap =-ABS(leg.resetGap);
				// Reset Calendar
                XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList3,i),"resetHolidays"),sresetCal);
				leg.resetCal = calypso2ARMCalendar(sresetCal);
				// Index day count
                XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList3,i),"dayCount"),sindexDaycount);
				break;
			}
		}

       leg.indexTerm= index.substr(4,index.size()-1);
	   if(leg.indexTerm =="")	 leg.indexTerm = "1Y";	
		
		if(sindexDaycount!="")
			leg.indexDayCount = calypso2ARMDaycount(sindexDaycount.c_str());
    }else{
	        // Reset Gap
            //string sresetGap="";XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"FundingResetGap"),sresetGap);
            string sresetGap="";XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"FixingResetGap"),sresetGap);
            try{
                if((sresetGap!="")&&(sresetGap.size()>4))
                    leg.resetGap = atoi(sresetGap.substr(4,sresetGap.find_first_of(' ' ,4)).c_str());
                else
                    leg.resetGap = 0;
            }catch(Exception e){
                leg.resetGap = 0;
            }


            //TODO verify if its the real path
            // Index day count
            try{
                leg.indexDayCount = leg.dayCount;
            }catch(Exception e){
                leg.indexDayCount = leg.dayCount;
            }


	        // Reset Calendar // TODO bad path to chamge
	        //mid = jni->GetMethodID(clsLeg, "getResetCalendar","()Ljava/lang/String;");	
	        string sresetPayCal;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"couponDateAdjustement/holidays"),sresetPayCal);
	        if(leg.resetCal =="")
		        leg.resetCal =leg.payCal ;
	        else
		        leg.resetCal = calypso2ARMCalendar(sresetPayCal);
    }
    
    
    return leg;
}



CalypsoLeg getCRFLeg(MSXML2::IXMLDOMNodePtr xmlNode) {

	CalypsoLeg leg;

    //Type
    string stype;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"fixFloat"),stype);
    leg.type=(stype=="Form")?EXOTIC:(stype=="Float")?FLOAT:FIX;


    //startDate
    XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"startDate"),leg.startDate,"YYYYMMDD");


    //endDate
    XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"endDate"),leg.endDate,"YYYYMMDD");


    //fixEndDate
    leg.fixEndDate=leg.endDate;

	// Fix rate
	try{
        XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"fixedRate"),leg.fixRate );
        leg.fixRate *= 100;
    }catch(Exception e){
        leg.fixRate= 0;
    }

	// Spread
	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"spread"),leg.spread);
	
    // Nominal
	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"notional/quantity"),leg.nominal);
	
	// Currency
	string ccy;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"notional/currency"),ccy);
	leg.currency = ARM_Currency(ccy.c_str());
	
	// Payment Frequency
	string sfreq;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"couponFrequency"),sfreq);
    leg.freq = calypso2ARMFreq(sfreq.c_str());


    // Decomp flag
    string cmp;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"Cmp"),cmp);
    leg.decompFlag =(cmp=="true")?1:0;
	
    // TODO remove special conversion
    // Decomp Frequency
    leg.decompFreq = leg.freq;
    if(leg.decompFlag == 1){
        string sCompoundFrequency;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"CompoundFrequency"),sCompoundFrequency);
    	if(sCompoundFrequency   == "Semi-annually")
            sCompoundFrequency="SA";
        if(sCompoundFrequency   == "Per annum")
            sCompoundFrequency="PA";
        leg.decompFreq =calypso2ARMFreq(sCompoundFrequency.c_str());
    }
    if(leg.decompFreq==K_DEF_FREQ){
        leg.decompFreq = leg.freq;
    }
    

	// Day count
	string sdayCount;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"dayCount"),sdayCount);
    leg.dayCount = calypso2ARMDaycount(sdayCount.c_str());

	// Payment calendar
	string spayCal;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"couponDateAdjustement/holidays"),spayCal);
	leg.payCal = calypso2ARMCalendar(spayCal);

	// Payment timing
	//string spayTiming;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"resetTiming2"),spayTiming);
	string spayTiming;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"PmtResetTiming"),spayTiming);
	leg.payTiming = calypso2ARMTiming(spayTiming.c_str());


	// Stub Rule
	string sstubRule;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"couponStubRule"),sstubRule);
	leg.stubRule =  calypso2ARMStubRule(sstubRule.c_str());


	// Payment int rule
	string sintRule;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"periodRule"),sintRule);
	leg.intRule =  calypso2ARMInterestRule(sintRule.c_str());



	// Notional exchange flag
	bool actual =false;
	bool initEx=false;
	bool finalEx=false;
	try {
		string sactual;	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"Actual"),sactual);
		actual= (sactual=="true")?true:false;
		//TODO add Condition XCCySwap instatnce to add in argument
        if(actual) {
            try{
              //  string sinitEx;	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"PrincipalExchange/Initial"),sinitEx);
		        string sinitEx;	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"PrincipalExchangeType/Initial"),sinitEx);
		        initEx = (sinitEx=="true")?true:false;
		    //    string sfinalEx;	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"PrincipalExchange/Final"),sfinalEx);
		        string sfinalEx;	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"PrincipalExchangeType/Final"),sfinalEx);
		        finalEx = (sfinalEx=="true")?true:false;
            }catch(Exception e){
                initEx = false;
	            finalEx = false;
	         }
        }	
       }catch (Exception e){
	    actual = false;
	    initEx = false;
	    finalEx = false;
	   }
    string snotionalEx = calypso2ARMNotionalExchange(initEx, finalEx);
	leg.notionalExchangeFlag = ARM_NotionalExchange(snotionalEx .c_str());


    // Index term
	try{
        XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"rateIndex/tenor"),leg.indexTerm);
    }catch(Exception e){
        leg.indexTerm = "1Y";
    }
    
    // Index type
	try{
        string sindexName;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"rateIndex/indexName"),sindexName);
	    leg.indexType = calypso2ARMIndexType(sindexName,leg.indexTerm);
    }catch(Exception e){}
	

	
    

	// Reset timing
    //string sresetTiming;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"resetTiming1"),sresetTiming);
//    string sresetTiming;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"FixingResetTiming"),sresetTiming);
//	leg.resetTiming = calypso2ARMTiming(sresetTiming.c_str());


    	// Reset Frequency  // TODO bad path to chamge
	//mid = jni->GetMethodID(clsLeg, "getResetFrequency", "()Ljava/lang/String;");
	string sfundResetFreq;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"couponFrequency"),sfundResetFreq);
	leg.resetFreq = calypso2ARMFreq(sfundResetFreq.c_str());
    
    

    // Forward Rule  
	string sfwdRule;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"couponDateAdjustement/dateRoll"),sfwdRule);
	leg.fwdRule = calypso2ARMDateRoll(sfwdRule);
	
    // Leverage  //TODO Find correct path
	//mid = jni->GetMethodID(clsLeg, "getLeverage", "()D");
	leg.leverage = 1;


	// Payment gap
	//mid = jni->GetMethodID(clsLeg, "getPayGap", "()I");
	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"couponDateAdjustement/days"),leg.payGap);

	string index;
	// Cash flows
    MSXML2::IXMLDOMNodeListPtr xmlNodeList= XMLTools::selectNodes(xmlNode,"Listecashflow/cashflow"); 
	xmlNodeList->get_length(&leg.flowCount);
   	ARM_Date flowResetDate,flowStartDate,flowEndDate;
    double flowNominalValue,flowSpreadValue,flowRateValue;
    int dateInit=0;
    leg.resetGap =0;

    for (int i = 0; i< leg.flowCount; i++) {
        string stype;XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"type"),stype);
        if(stype == "INTEREST"){
            XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"startDate"),flowStartDate,"YYYYMMDD");
            
            try{
                XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"resetDate"),flowResetDate,"YYYYMMDD");
                leg.flowResetDatesVector.push_back(flowResetDate.GetJulian());
            }catch(Exception e){
                //TODO finalize version
                ARM_Date* tmpDate =  new ARM_Date(flowStartDate.GetJulian());

                int gap =(leg.resetGap ==0)?-2:leg.resetGap ;
                
                if(gap >0)
                    tmpDate->NextBusinessDay(ABS(gap),leg.currency.GetCcyName());
                 else
                    tmpDate->PreviousBusinessDay(ABS(gap), (char*)(leg.payCal.c_str()));

                leg.flowResetDatesVector.push_back(tmpDate->GetJulian());
                if(tmpDate) {
		             delete tmpDate;
	        	     tmpDate = NULL;
                }
            }

            XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"endDate"),flowEndDate,"YYYYMMDD");
            XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"notional"),flowNominalValue);
            XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"spread"),flowSpreadValue);
            XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"rate"),flowRateValue);

            leg.flowStartDatesVector.push_back(flowStartDate.GetJulian());
		    leg.flowEndDatesVector.push_back(flowEndDate.GetJulian());
		    leg.flowNominalVector.push_back(ABS(flowNominalValue));
		    leg.flowSpreadVector.push_back(flowSpreadValue);
            leg.flowRateVector.push_back(100 * flowRateValue);

            if(leg.type==EXOTIC){
                long nbFormulaParameters=0;
                MSXML2::IXMLDOMNodeListPtr	xmlNodeList2 = XMLTools::selectNodes(XMLTools::get_item(xmlNodeList,i),"ListFormulaParameters/FormulaParameter"); 
				xmlNodeList2->get_length(&nbFormulaParameters); 
				string formulaParameterName;
				
                double dlevrage, dcpnMin, dcpnMax,dtstrike;
                
				for(int j=0;j<nbFormulaParameters;j++) {
					XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList2,j),"Parameter_Name"),formulaParameterName);
					
					if(formulaParameterName =="aLower")
						XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList2,j),"Parameter_Value"),dcpnMin);
					if(formulaParameterName == "aCoef")
						XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList2,j),"Parameter_Value"),dlevrage);
					if(formulaParameterName == "aUpper")
						XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList2,j),"Parameter_Value"),dcpnMax);
					if(formulaParameterName == "aFixedRate")
						XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList2,j),"Parameter_Value"),dtstrike);
					if(formulaParameterName == "aFloatingRate")
						XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList2,j),"Parameter_Value"),index);
				}
                if(dlevrage==0){
					dlevrage=1;
					index = "";
				}
	            
                if (dcpnMax == 0) dcpnMax= 100;

                if (!(index == "") && (leg.fixEndDate==leg.endDate)&&(dateInit==0)) {
        				leg.fixEndDate = flowStartDate;
						dateInit=1;
        		}

                leg.strikeVector.push_back(dtstrike);
				leg.leverageVector.push_back(dlevrage);
				leg.cpnMinVector.push_back(dcpnMin);
				leg.cpnMaxVector.push_back(dcpnMax);
            }else {
                //leg.strikeVector.push_back(100 * flowRateValue);
				leg.strikeVector.push_back(flowRateValue);
				//leg.leverageVector.push_back(0);
				//leg.cpnMinVector.push_back(0);
				//leg.cpnMaxVector.push_back(10000);
            }
        }
    }
    leg.flowCount=leg.flowStartDatesVector.size();
	
    
    if(leg.type==EXOTIC){
        leg.indexDayCount = leg.dayCount;
	
		MSXML2::IXMLDOMNodeListPtr xmlNodeList3 = XMLTools::selectNodes(xmlNode,"ListeExoticVariables/ExoticVariable");
        if(xmlNodeList3==NULL)
                xmlNodeList3 = XMLTools::selectNodes(xmlNode,"ListeExoticVariable/ExoticVariable"); 
        try{
            string sname;XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList3,0),"variableName"),sname);
        }catch(Exception e){
            xmlNodeList3 = XMLTools::selectNodes(xmlNode,"ListeExoticVariable/ExoticVariable"); 
        }
		
		long exoticVariableCount;
		xmlNodeList3->get_length(&exoticVariableCount);

		string sname,sindexDaycount,sresetCal,sresetInArrears;
		for( i=0;i<exoticVariableCount;i++) 
		{
			XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList3,i),"variableName"),sname);
     
			if(sname==index){
				// Reset Gap
                XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList3,i),"resetDays"),leg.resetGap);
				leg.resetGap =-ABS(leg.resetGap);
				// Reset Calendar
                XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList3,i),"resetHolidays"),sresetCal);
				leg.resetCal = calypso2ARMCalendar(sresetCal);
				// Index day count
                XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList3,i),"dayCount"),sindexDaycount);
                //resetTiming
				XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList3,i),"resetInArrears"),sresetCal);
				if (sresetInArrears == "true")
					leg.resetTiming = K_ARREARS;
				else
					leg.resetTiming = K_ADVANCE;

				break;
			}
		}

       leg.indexTerm= index.substr(4,index.size()-1);
	   if(leg.indexTerm =="")	 leg.indexTerm = "1Y";	
		
		if(sindexDaycount!="")
			leg.indexDayCount = calypso2ARMDaycount(sindexDaycount.c_str());
    }else{
	        // Reset Gap
            //string sresetGap="";XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"FundingResetGap"),sresetGap);
            string sresetGap="";XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"FixingResetGap"),sresetGap);
            try{
                if((sresetGap!="")&&(sresetGap.size()>4))
                    leg.resetGap = atoi(sresetGap.substr(4,sresetGap.find_first_of(' ' ,4)).c_str());
                else
                    leg.resetGap = 0;
            }catch(Exception e){
                leg.resetGap = 0;
            }


            //TODO verify if its the real path
            // Index day count
            try{
                leg.indexDayCount = leg.dayCount;
            }catch(Exception e){
                leg.indexDayCount = leg.dayCount;
            }


	        // Reset Calendar // TODO bad path to chamge
	        //mid = jni->GetMethodID(clsLeg, "getResetCalendar","()Ljava/lang/String;");	
	        string sresetPayCal;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"couponDateAdjustement/holidays"),sresetPayCal);
	        if(leg.resetCal =="")
		        leg.resetCal =leg.payCal ;
	        else
		        leg.resetCal = calypso2ARMCalendar(sresetPayCal);
    }
    
    
    return leg;
}
