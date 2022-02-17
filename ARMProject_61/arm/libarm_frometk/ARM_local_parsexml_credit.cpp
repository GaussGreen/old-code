
#include <ARM\libarm_local\firstToBeIncluded.h>

#include <ICMKernel\crv\icm_constant_piecewise.h>
#include <ICMKernel\glob\icm_corrmatrix.h>


#include <atlbase.h>
#include <ARM\libarm_frometk\XMLTools.h>
#include <ARM\libarm_frometk\arm_local_parsexml_util.h>
#include "VariantTools.h"

#include <ICMKernel\inst\icm_cds.h>
#include <ICMKernel\inst\icm_mez.h>
#include <ICMKernel\inst\icm_nthtd.h>
#include <ICMKernel\inst\icm_collateral.h>

#include <ARM\libarm_frometk\arm_local_parsexml.h>

#include <ARM\libarm_frometk\arm_local_etoolkit_for_icm.h>
#include <ARM\libarm_frometk\arm_local_parsexml_for_icm.h>


#include <ARM\libarm_frometk\PaserManagerUtilities.h>

#include <ARM\libarm_frometk\arm_local_paesexml_calypso.h> 

#include <libCCatl\CCatl.h>
#include <ARM\libarm_frometk\arm_local_parsexml_credit.h>


ICM_Cds* ARMLOCAL_ParseCalypsoCDS(const string& chaineXML ,
										 const ARM_Date& date,
										 const string modeltype,
										 string& BookName,
										 string& calypsoId)
{
	double spread ;
	ARM_Date EffectiveDate , Maturitydate ; //, ReferenceDate , FirstCpneffDate; not set in xml
	ARM_Date ProtectionStartDate , ProtectionEndDate;

	int Frequency ;
	int DayCount ;
	
	double FixedPayerAmount ; 
	double FloatingPayerAmount ;
	int StubRule , intRule, CreditLag , StartAdj ;
	bool IncludeMaturity ; 
	std::string Currency ;
	double Binary ;
	ICM_Cds* cds = NULL;

	try
	{
		MSXML2::IXMLDOMDocumentPtr XMLDoc;
		XMLDoc = XMLTools::LoadXML(chaineXML);


		MSXML2::IXMLDOMNodePtr theNode ;

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/PremiumLeg/FixedRate") ;
		XMLTools::convert(theNode,spread);

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/Description/StartDate");
		
		XMLTools::convert(theNode,EffectiveDate,"YYYYMMDD");

		
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/Description/EndDate");
		
		XMLTools::convert(theNode,Maturitydate,"YYYYMMDD");

		
		ProtectionStartDate = EffectiveDate ;
		ProtectionEndDate = Maturitydate ;

		std::string tmp ;
		
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/PremiumLeg/PremiumCpnFreq");
		XMLTools::convert(theNode,tmp);
		Frequency = FromSummitFreqToARMFreq(tmp.c_str());
		
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/PremiumLeg/PremiumDayCount");
		XMLTools::convert(theNode,tmp);
		
		DayCount = FromSummitDaycountToARMDaycount(tmp.c_str()) ;


		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/Description/PVCcy");
		XMLTools::convert(theNode,Currency);
		

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/Settlement/SettlementLagDays");
		XMLTools::convert(theNode,CreditLag);
		
		if ( CreditLag  <= 0) CreditLag = 30 ;
		
		
		//Notional = fisrt Notional in FluxPrime

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/PremiumLeg/ListeFluxPrime/FluxPrime[0]/Notional");
		XMLTools::convert(theNode,FixedPayerAmount);

		FloatingPayerAmount = FixedPayerAmount;

		/*theNode = XMLTools::selectSingleNode(XMLDoc,"ListTrade/Trade[0]/ProtectionLeg/Settlement/RecoveryRateFix");
		XMLTools::convert(theNode,Binary);*/
		
		// Not Set parameter in XML ?
		Binary = -999.0;
		
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/Settlement/RecoveryRateFix");
		XMLTools::convert(theNode,Binary);

		if ( Binary  <= 0) Binary = -999 ;

		intRule = K_ADJUSTED;
		StartAdj = K_ADJUSTED;
		
		StubRule = ARM_ConvStubRule("SS");
		IncludeMaturity = true ;
		std::string name ="" ;



		cds = new ICM_Cds(EffectiveDate,
						Maturitydate,
						0, // reference date
						0,// fstCpneffDate
						ProtectionStartDate,
						ProtectionEndDate,
						spread,
						FixedPayerAmount, FloatingPayerAmount,
						Frequency,
						DayCount,
						qACCRUED_SETTLED,
						Currency,
						StubRule,
						CreditLag,
						-1 , // Freqdefleg
						intRule,
						IncludeMaturity,
						StartAdj,
						Currency,
						qRunning_Leg,
						qStandart_Recovery_Leg,
						name,
						Binary);


	}
	catch(...)
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Pb in creating the CDS From Calypso");
	}



	return cds ;
}

 ICM_Nthtd* ARMLOCAL_ParseCalypsoNTD(const string& chaineXML,
										   const ARM_Date& date,
										   const string modeltype,
										   string& BookName,
										   string& calypsoId)
{
	double spread ;
	ARM_Date EffectiveDate , Maturitydate ; //, ReferenceDate , FirstCpneffDate; not set in xml
	int FirstNumDefault, LastNumDefault;
	int Frequency ;
	int DayCount ;
	
	double TradedNotional ; 
	double RcvFee;
	double IssuerNotional ;
	int StubRule , intRule, CreditLag , StartAdj ;
	bool IncludeMaturity ; 
	std::string Currency ;
	double Binary ;
	qPAYMENT_PREMIUM_LEG		AccruedOnDefault;
	ICM_Nthtd* ntd = NULL;

	try
	{
		MSXML2::IXMLDOMDocumentPtr XMLDoc;
		XMLDoc = XMLTools::LoadXML(chaineXML);


		MSXML2::IXMLDOMNodePtr theNode ;

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/PremiumLeg/FixedRate") ;
		XMLTools::convert(theNode,spread);

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/Description/StartDate");
		
		XMLTools::convert(theNode,EffectiveDate,"YYYYMMDD");

		
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/Description/EndDate");
		
		XMLTools::convert(theNode,Maturitydate,"YYYYMMDD");


		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/StartDefault");
		
		XMLTools::convert(theNode,FirstNumDefault);

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/EndDefault");
		
		XMLTools::convert(theNode,LastNumDefault);

		std::string tmp ;
		
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/PremiumLeg/PremiumCpnFreq");
		XMLTools::convert(theNode,tmp);
		Frequency = FromSummitFreqToARMFreq(tmp.c_str());

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/PremiumLeg/PremiumDayCount");
		XMLTools::convert(theNode,tmp);
		
		DayCount = FromSummitDaycountToARMDaycount(tmp.c_str()) ;

		
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/PremiumLeg/PremiumType");
		XMLTools::convert(theNode,tmp);
		if (tmp == "PAY_ACCRUAL") 
			AccruedOnDefault = qACCRUED_SETTLED ;
		//ICM_EnumsCnv::cnv(tmp,AccruedOnDefault); 

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/Description/PVCcy");
		XMLTools::convert(theNode,Currency);
		

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/Settlement/SettlementLagDays");
		XMLTools::convert(theNode,CreditLag);
		
		if ( CreditLag  <= 0) CreditLag = 30 ;
		
		
		//Notional = fisrt Notional in FluxPrime

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/PremiumLeg/ListeFluxPrime/FluxPrime[0]/Notional");
		XMLTools::convert(theNode,TradedNotional);
		
		RcvFee = abs(TradedNotional)/TradedNotional;
		
		TradedNotional = abs(TradedNotional);
		
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/Settlement/RecoveryRateFix");
		XMLTools::convert(theNode,Binary);

		if ( Binary  <= 0) Binary = -999 ;

		intRule = K_ADJUSTED;
		StartAdj = K_ADJUSTED;
		
		StubRule = ARM_ConvStubRule("SS");
		IncludeMaturity = true ;

		// Issuers
		
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/BasketElements/BasketItem[0]/Size");
		XMLTools::convert(theNode,IssuerNotional);
		

		MSXML2::IXMLDOMNodeListPtr Underlyings ;

		Underlyings = XMLTools::selectNodes(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/BasketElements/BasketItem");
		
		long NbNodes(0);
		Underlyings ->get_length(&NbNodes);
		vector<string> labels_(NbNodes); 

		vector<double> notionals(NbNodes);
		
		string Issuer;
		double Notional;
		for ( int i = 0;i<labels_.size();i++)
		{
			XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(Underlyings,i),"Issuer"),Issuer);
			// Correlator ISIN (adding _J when subordinate name)???
			XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(Underlyings,i),"Size"),Notional);
			labels_[i] = Issuer;
			notionals[i] = Notional;
			
		}

		ntd = new ICM_Nthtd((ARM_Date) EffectiveDate,
							  (ARM_Date) Maturitydate,
							  0,
							  0,
							   spread,
							   intRule, //	//K_ADJUSTED, intRule
							   StartAdj , //, K_ADJUSTED,	// startAdj
							   FirstNumDefault,
							   LastNumDefault,
							   // nbissuers,
							   labels_,
							   notionals,
							   Frequency,
							   DayCount,
							   RcvFee*IssuerNotional, 
							   AccruedOnDefault,
							   Currency, 
							   RcvFee*IssuerNotional,
							   StubRule,
							   CreditLag,
							   Frequency,
							   Binary,
							   Currency,
							   IncludeMaturity
							   );



	}
	catch(...)
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Pb in creating the NTD from Calypso");
	}

	return ntd;

}



extern ICM_Mez* ARMLOCAL_ParseCalypsoCDO(const string& chaineXML,
										   const ARM_Date& date,
										   const string modeltype,
										   string& BookName,
										   string& calypsoId)
{

	ICM_Mez* Mezz = NULL;

	double spread ;
	ARM_Date EffectiveDate , Maturitydate ; //, ReferenceDate , FirstCpneffDate; not set in xml
	double MezzAmount;
	double SubAmount ;
	int Frequency ;
	int DayCount ;
	
	double TradedNotional ; 
	double RcvFee;
	double IssuerNotional ;
	int StubRule , intRule, CreditLag , StartAdj ;
	bool IncludeMaturity ; 
	std::string Currency ;
	double Binary ;
	qPAYMENT_PREMIUM_LEG		AccruedOnDefault;


	try
	{

		MSXML2::IXMLDOMDocumentPtr XMLDoc;
		XMLDoc = XMLTools::LoadXML(chaineXML);


		MSXML2::IXMLDOMNodePtr theNode ;

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/PremiumLeg/FixedRate") ;
		XMLTools::convert(theNode,spread);

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/Description/StartDate");
		
		XMLTools::convert(theNode,EffectiveDate,"YYYYMMDD");

		
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/Description/EndDate");
		
		XMLTools::convert(theNode,Maturitydate,"YYYYMMDD");



		std::string tmp ;
		
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/PremiumLeg/PremiumCpnFreq");
		XMLTools::convert(theNode,tmp);
		Frequency = FromSummitFreqToARMFreq(tmp.c_str());

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/PremiumLeg/PremiumDayCount");
		XMLTools::convert(theNode,tmp);
		
		DayCount = FromSummitDaycountToARMDaycount(tmp.c_str()) ;

		
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/PremiumLeg/PremiumType");
		XMLTools::convert(theNode,tmp);
		if (tmp == "PAY_ACCRUAL") 
			AccruedOnDefault = qACCRUED_SETTLED ;
		//ICM_EnumsCnv::cnv(tmp,AccruedOnDefault); 

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/Description/PVCcy");
		XMLTools::convert(theNode,Currency);
		

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/Settlement/SettlementLagDays");
		XMLTools::convert(theNode,CreditLag);
		
		if ( CreditLag  <= 0) CreditLag = 30 ;
		
		
		//Notional = fisrt Notional in FluxPrime

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/PremiumLeg/ListeFluxPrime/FluxPrime[0]/Notional");
		XMLTools::convert(theNode,TradedNotional);
		
		RcvFee = abs(TradedNotional)/TradedNotional;
		
		TradedNotional = abs(TradedNotional);
		
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/Settlement/RecoveryRateFix");
		XMLTools::convert(theNode,Binary);

		if ( Binary  <= 0) Binary = -999 ;

		intRule = K_ADJUSTED;
		StartAdj = K_ADJUSTED;
		
		// default values: TODO get from XML
		StubRule = ARM_ConvStubRule("SS");
		IncludeMaturity = true ;

		vector<string> labels_; 

		ARM_Vector* notionals = new ARM_Vector();

		string Issuer;
		double Notional;
		string seniority;
		double PercentLow ;
		double PercentHigh ;
		double TradedCoeff ;

		/****************      if case of an Index CDO        ****************/
		
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/Description/ProductDefinition");
		XMLTools::convert(theNode,tmp);
		if (tmp == "CDSIndexTranche")
		{

			labels_.resize(125); 
			notionals->Resize(125);

			Notional = 10000000. ;
			theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/Issuer");
			XMLTools::convert(theNode,Issuer);

			for ( int i = 0;i<labels_.size();i++)
			{
				labels_[i] = Issuer;
				notionals->Elt(i) = Notional;
			
			}

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/StartLoss");
		XMLTools::convert(theNode,PercentLow);
		
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/EndLoss");
		XMLTools::convert(theNode,PercentHigh);

		MezzAmount = 125 * Notional * (PercentLow) ;
		SubAmount = 125* Notional * (PercentHigh - PercentLow) ;

		TradedCoeff=TradedNotional/SubAmount;
		}
		else
		{

		// Issuers
		
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/BasketElements/BasketItem[0]/Size");
		XMLTools::convert(theNode,IssuerNotional);
		

		MSXML2::IXMLDOMNodeListPtr Underlyings ;
	

		Underlyings = XMLTools::selectNodes(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/BasketElements/BasketItem");
		
		long NbNodes(0);
		Underlyings ->get_length(&NbNodes);
		labels_.resize(NbNodes); 
		notionals->Resize(NbNodes);
		
		MezzAmount =0. ;
		SubAmount=0. ;
		for ( int i = 0;i<labels_.size();i++)
		{
			XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(Underlyings,i),"Issuer"),Issuer);
			// Correlator ISIN (adding _J when subordinate name) : TO CHANGE
			XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(Underlyings,i),"Seniority"),seniority);
			if ( seniority == "SUBORDINATE")
				Issuer+="_J";
			
			XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(Underlyings,i),"Size"),Notional);
			labels_[i] = Issuer;
			notionals->Elt(i) = Notional;
			MezzAmount += Notional;
		}

		double PercentLow ;
		double PercentHigh ;
		double init_sub = SubAmount;
		double init_mezz = MezzAmount;
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/StartLossPercent");
		XMLTools::convert(theNode,PercentLow);
		
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/EndLossPercent");
		XMLTools::convert(theNode,PercentHigh);

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/StartLoss");
		XMLTools::convert(theNode,init_sub);

		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/ProtectionLeg/EndLoss");
		XMLTools::convert(theNode,init_mezz);


// FIXMEFRED: mig.vc8 (30/05/2007 17:53:21):cast
		ARM_Vector* tmp_Notionals = new ARM_Vector(notionals);

		ConversionAmounts(labels_.size(),tmp_Notionals,init_sub,init_mezz,SubAmount,MezzAmount,notionals);
		TradedCoeff=TradedNotional/SubAmount;
		if (tmp_Notionals) 
			delete tmp_Notionals;
		tmp_Notionals = NULL;
		}

		ICM_Leg* FeeLeg = NULL;
		ICM_Leg* DefLeg = NULL;
		ICM_Collateral* Collat = NULL ;
		ICM_Cds* CdsGen = NULL; 
		ARM_IRIndex* Index = NULL;
		ICM_Credit_Index* CreditIndex = NULL ;
		
		std::string subtype ;
		theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/Description/IXIS_Classification");
		XMLTools::convert(theNode,subtype);
		
		/***********************  Case of a CM Tranche       *******************************/
		
		if (subtype == "Constant Maturity CDS" )
		{
			// TODO

			// Underlying Index


			// FeeLeg

			// defleg




		}

		
		else
		{

			theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/PremiumLeg/FixFloat");
			XMLTools::convert(theNode,subtype);

		/***********************  Case of a RiskyLeg Tranche *******************************/

			if ( subtype=="Float")
			{
				ARM_INDEX_TYPE armIndexType;
				std::string FloatSpread ;
				theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/PremiumLeg/Index");
				XMLTools::convert(theNode,FloatSpread);
				
				theNode = XMLTools::selectSingleNode(XMLDoc,"/Extraction/ListTrade/Trade[0]/PremiumLeg/IndexTerm");
				XMLTools::convert(theNode,tmp);
				if (FloatSpread == "EURIB") FloatSpread = "EURIBOR";
				FloatSpread+=tmp ;

				ARM_Currency* ccy=NULL;
				ccy = new ARM_Currency(Currency.c_str()) ;
				armIndexType =(ARM_INDEX_TYPE) ARM_ConvIrType (FloatSpread.c_str());
				Index = new ARM_IRIndex(armIndexType,
											 Frequency,
											 Frequency,
											 ccy); // amlways a credit index ???

				FeeLeg = new ICM_Leg((ARM_Date)EffectiveDate, 
									 (ARM_Date)Maturitydate, 
									 0,
									 NULL,
									 spread,
									 SubAmount,	
									 NULL,
									 NULL,
									 NULL,
									 Frequency , 
									 DayCount , 
									 K_ARREARS, // "ARR"
									 intRule ,
									 StubRule ,
									 Currency,
									 Currency,
									 qSwapLeg,
									 IncludeMaturity,
									 StartAdj,
									 Index,
									 Binary,
									 ISSUER_UNDEFINE,
									 0,
									 AccruedOnDefault);
				DefLeg = new ICM_Leg((ARM_Date)EffectiveDate, 
									 (ARM_Date)Maturitydate, 
									 0,
									 NULL,
									 1,
									 SubAmount,	
									 NULL,
									 NULL,
									 NULL,
									 Frequency , 
									 DayCount , 
									 K_ARREARS, // "ARR"
									 intRule ,
									 StubRule ,
									 Currency,
									 Currency,
									 qStandart_Recovery_Leg,
									 IncludeMaturity,
									 StartAdj,
									 NULL,
									 Binary,
									 ISSUER_UNDEFINE,
									 0,
									 AccruedOnDefault);

			}
			else
			{
			
		/***********************  Case of a bespoke Tranche	 *******************************/

				FeeLeg = new ICM_Leg((ARM_Date)EffectiveDate, 
									 (ARM_Date)Maturitydate, 
									 0,
									 NULL,
									 spread,
									 SubAmount,	
									 NULL,
									 NULL,
									 NULL,
									 Frequency , 
									 DayCount , 
									 K_ARREARS, // "ARR"
									 intRule ,
									 StubRule ,
									 Currency,
									 Currency,
									 qRunning_Leg,
									 IncludeMaturity,
									 StartAdj,
									 NULL,
									 Binary,
									 ISSUER_UNDEFINE,
									 0,
									 AccruedOnDefault);


				DefLeg = new ICM_Leg((ARM_Date)EffectiveDate, 
									 (ARM_Date)Maturitydate, 
									 0,
									 NULL,
									 1,
									 SubAmount,	
									 NULL,
									 NULL,
									 NULL,
									 Frequency , 
									 DayCount , 
									 K_ARREARS, // "ARR"
									 intRule ,
									 StubRule ,
									 Currency,
									 Currency,
									 qStandart_Recovery_Leg,
									 IncludeMaturity,
									 StartAdj,
									 NULL,
									 Binary,
									 ISSUER_UNDEFINE,
									 0,
									 AccruedOnDefault);

				}


			CdsGen = new ICM_Cds(FeeLeg,DefLeg);
			CdsGen->SetPorS(RcvFee) ;
			CdsGen ->SetTradedCoef(TradedCoeff);

			//Collat = new ICM_Collateral(labels_,ARM_Vector(CreateARMVectorFromVECTOR(notionals,notionals.size())));
			Collat = new ICM_Collateral(labels_,*notionals);

			Mezz = new ICM_Mez(CdsGen,MezzAmount,*Collat,Binary);
			Mezz->SetPorS(1);
		

		}

		if (notionals)
			delete notionals;
		notionals = NULL;
/*

		Mezz = new ICM_Mez((ARM_Date) EffectiveDate,
								(ARM_Date) Maturitydate,
								0,
								0,
							   (spread/10000.),
							   intRule,  //	K_ADJUSTED  intRule
							   StartAdj, //K_ADJUSTED  adjStartDate	
							   MezzAmount,
							   TradedNotional,
							   labels_,
							   notionals,
							   Frequency,
							   DayCount,
							   RcvFee*TradedNotional, 
							   AccruedOnDefault,
							   Currency, 
							   RcvFee*TradedNotional,
							   StubRule,
							   CreditLag,
							   Frequency,
							   Binary,
							   Currency,
							   qRunning_Leg,
							   qStandart_Recovery_Leg,
							   IncludeMaturity);
	
*/

	}
	catch(...)
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Pb in creating the CDO from Calypso");
	}


	return Mezz ;

}
// function using calypso addins dev
void ConversionAmounts(	int nbnames ,
							  ARM_Vector* init_notionals ,
							  double sub_init,
							  double end_init,
							  double& sub_out,
							  double& tranche_out,
							  ARM_Vector* out_notionals)
{
	const double NOTIONEL_UNIT = 10000000. ;
	int i=0;
	double initial_ptf_size =0;
	double ptf_size =0;
	double weight=0 ;
	out_notionals->Resize(nbnames); 
	
	for(i=0;i<nbnames;i++)
		initial_ptf_size+= init_notionals->Elt(i);

	initial_ptf_size = round(initial_ptf_size);
	bool UseMin =true;
	const double EPSILON = 1.e-8 ;
	int indexmin = 0;

	for(i=0;i<nbnames;i++)
	{
		if ( fabs(init_notionals->Elt(i)) <fabs(init_notionals->Elt(indexmin)))
			indexmin=i ;
		if( init_notionals->Elt(i) >0)
			ptf_size+=init_notionals->Elt(i);
	}

	weight = fabs(init_notionals->Elt(indexmin)/ptf_size) ;
	
	if (weight==0)
	{
		UseMin=false;
	}
	else
	{
		double div=0;
		double w =0;
		for(i=0;i<nbnames;i++)
		{
			w= init_notionals->Elt(i)/ptf_size ;
			div = ( w/weight) * 2 ;
			if ( fabs(div - round(div)) > EPSILON)
			{
				UseMin = false;
				break ;
			}
		}
	}


	double startprotect = sub_init/initial_ptf_size;
     double endprotect = end_init/initial_ptf_size;

	 if(UseMin) {
		 double refNot = init_notionals->Elt(indexmin);
		 if(refNot<0) {
			out_notionals->Elt(indexmin)= -NOTIONEL_UNIT;
		 } else {
			out_notionals->Elt(indexmin) = NOTIONEL_UNIT;
		 }
 
		 for (i=0;i<nbnames;i++) {
			 if(i!=indexmin) {
				out_notionals->Elt(i) =round(out_notionals->Elt(indexmin)*init_notionals->Elt(i)/refNot);
			 }
		 }
	 } else {
		 for (i=0;i<nbnames;i++) {
			out_notionals->Elt(i) = NOTIONEL_UNIT*100 * init_notionals->Elt(i)/initial_ptf_size;
		 }
	 }


	  double size_ptf_out = 0.;
     for (i=0;i<nbnames;i++) {
		 size_ptf_out += out_notionals->Elt(i);
	 }
     size_ptf_out = round(size_ptf_out);
     tranche_out = round(size_ptf_out*startprotect);
     sub_out = round(size_ptf_out*(endprotect-startprotect));


}