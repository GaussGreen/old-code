
#include "firsttoinc.h"
#include "icm_enums.h"
#include "ARMKernel/glob/expt.h"
//	-----------------------------------------------------------------------------------------
void 
ICM_EnumsCnv::cnv(const std::string& name,qAccelerationStyle&ref,bool&ret,std::string& list)
{
	ret=true; 
	if (name=="ACC") ref=qACC ;
	else if (name=="1") ref=qACC; 
	else if (name=="NACC") ref=qNACC; 
	else if (name=="0") ref=qNACC; 
	else { ret=false;  list="ACC(1),NACC(0)"; }
	return; 
}
//	-----------------------------------------------------------------------------------------
void 
ICM_EnumsCnv::toString(qAccelerationStyle ref,std::string& name) 
{
	switch (ref) 
	{
	case qACC: name="ACC"; return; 
	case qNACC: name="NACC"; return; 
	} ; 
	name="" ; 
	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_EnumsCnv::toString: undefined "<<ref<<" for qAccelerationStyle"); 
}
//	-----------------------------------------------------------------------------------------
void 
ICM_EnumsCnv::cnv(const std::string& name,qKoStyle&ref,bool&ret,std::string& list)
{
	ret=true; 
	if (name=="KO") ref=qKO;
	else if (name=="1") ref=qKO;
	else if (name=="NKO") ref=qNKO; 
	else if (name=="0") ref=qNKO; 
	else { ret=false;  list="KO(1),NKO(0)"; }
	return; 
}

//	-----------------------------------------------------------------------------------------
void 
ICM_EnumsCnv::toString(qKoStyle ref,std::string& name) 
{
	switch (ref) 
	{
	case qKO: name="KO"; return; 
	case qNKO: name="NKO"; return; 
	} ; 
	name=""; 
	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_EnumsCnv::toString: undefined "<<ref<<" for qKoStyle"); 
}
//	-----------------------------------------------------------------------------------------
void 
ICM_EnumsCnv::cnv(const std::string& name,qUnderlying_Maturity_Style&ref,bool&ret,std::string& list)
{
	ret=true; 
	if (name=="CONSTANT") ref=qConstantMaturity;
	else if (name=="RESIDUAL") ref=qResidualMaturity;
	else { ret=false;  list="CONSTANT,RESIDUAL"; }
	return; 
}
//	-----------------------------------------------------------------------------------------
void 
ICM_EnumsCnv::toString(qUnderlying_Maturity_Style ref,std::string& name) 
{
	switch (ref) 
	{
	case qConstantMaturity: name="CONSTANT"; return; 
	case qResidualMaturity: name="RESIDUAL"; return; 
	} ; 
	name=""; 
	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_EnumsCnv::toString: undefined "<<ref<<" for qUnderlying_Maturity_Style"); 
}
//	-----------------------------------------------------------------------------------------
void 
ICM_EnumsCnv::cnv(const std::string& name,qCDS_ADJ&ref,bool&ret,std::string& list)
{
	ret=true; 
	const ICM_Traits<qCDS_ADJ>::list_t & items = ICM_Traits<qCDS_ADJ>::getList(); 
	std::map<std::string,qCDS_ADJ>::const_iterator it = 
		items.find(name); 
	if (it==items.end()) 
	{ 
		ret=false;  
		list="STDCDS,STDIDH5(TRX),STDIDX(M5),STDIDZ(Z5),STDID[M|Z][6|7|8],STDID[M|U|Z]4,STDINDEX,STDCDSQN,STDCDSSA,NONE";  
		return ;
	}
	ref=it->second; 
}
//	-----------------------------------------------------------------------------------------
void 
ICM_EnumsCnv::cnv(const std::string& name,qDEFCURVE_CALIB_ALGO&ref,bool&ret,std::string& list)
{
	ret=true; 
	const ICM_Traits<qDEFCURVE_CALIB_ALGO>::list_t & items = ICM_Traits<qDEFCURVE_CALIB_ALGO>::getList(); 
	std::map<std::string,qDEFCURVE_CALIB_ALGO>::const_iterator it = 
		items.find(name); 
	if (it==items.end()) 
	{ 
		ret=false;  
		list=""; 
		for(it=items.begin();it!=items.end();++it) 
		{ list += it->first ; list+=" " ; }
		return ;
	}
	ref=it->second; 
}
//	-----------------------------------------------------------------------------------------
void 
ICM_EnumsCnv::cnv(const std::string& name,qSENSITIVITY_TYPE&ref,bool&ret,std::string& list)
{
	ret=true; 
	const ICM_Traits<qSENSITIVITY_TYPE>::list_t & items = ICM_Traits<qSENSITIVITY_TYPE>::getList(); 
	std::map<std::string,qSENSITIVITY_TYPE>::const_iterator it = 
		items.find(name); 
	if (it==items.end()) 
	{ 
		ret=false;  
		list=""; 
		for(it=items.begin();it!=items.end();++it) 
		{ list += it->first ; list+=" " ; }
		return ;
	}
	ref=it->second; 
}
//	-----------------------------------------------------------------------------------------
void 
ICM_EnumsCnv::cnv(const std::string& name,qPAYMENT_PREMIUM_LEG&ref,bool&ret,std::string& list)
{
	ret=true; 
	const ICM_Traits<qPAYMENT_PREMIUM_LEG>::list_t & items = ICM_Traits<qPAYMENT_PREMIUM_LEG>::getList(); 
	std::map<std::string,qPAYMENT_PREMIUM_LEG>::const_iterator it = 
		items.find(name); 
	if (it==items.end()) 
	{ 
		ret=false;  
		list=""; 
		for(it=items.begin();it!=items.end();++it) 
		{ list += it->first ; list+=" " ; }
		return ;
	}
	ref=it->second; 
}
//	---------------------------------------------------------------------------------------------
const std::map<std::string,qCMPMETH>& 
ICM_Traits<qCMPMETH>::getList()
{
	static std::map<std::string,qCMPMETH> list; 
	if (list.empty())
	{
		list["NPV"]=qCMPPRICE; 
		list["FEELEG"]=qCMPFEELEGPV;  
		list["SPREAD"]= qCMPSPREAD; 
		list["DEFLEG"]=qCMPDEFLEGPV;
		list["ACCRUED_PV"]=qCMPACCRUED ; 
		list["ACCRUED"]=qCMPACCRUED;
		list["PREMIUM"]=qCMPPREMIUM ;
		list["FWDSPREAD"]=qCMPFWDSPREAD; 
		list["DURATION"]=qCMPDURATION;
		list["CORREL_UP"]=qCMPCORRELUP;
		list["CORREL_DOWN"]=qCMPCORRELDOWN;
		list["AVGCORRDEF"]=qCMPAVGCORRDEF;	
		list["FLATCORR"]=qCMPFLATCORR;
		list["OPT_DELTA"]=qCMP_OPT_AN_DELTA;
		list["OPT_VEGA"]=qCMP_OPT_AN_VEGA;
		list["OPT_GAMMA"]=qCMP_OPT_AN_GAMMA;
		list["FEES"] = qCMPFEES ; 
		list["EL"] = qCMPEL ; 
	}
	return list; 
}
//	---------------------------------------------------------------------------------------------
const std::map<std::string,qVECTMETH>& 
ICM_Traits<qVECTMETH>::getList()
{
	static std::map<std::string,qVECTMETH> list; 
	if (list.empty())
	{
		list["SOME"]=qSOME; 
		list["INDEXBASE"]=qINDEXBASE;
		list["MARKETBASE"]=qMARKETBASE; 
		list["TRUE_MARKETBASE"]=qTRUEMARKETBASE;
		list["TRUE_MARKETSPREAD"]=qTRUEMARKETSPREAD;
	}
	return list; 
}
//	---------------------------------------------------------------------------------------------
const std::map<std::string,qTWO_FACTORS_CORRELATION_TYPE>& 
ICM_Traits<qTWO_FACTORS_CORRELATION_TYPE>::getList()
{
	static std::map<std::string,qTWO_FACTORS_CORRELATION_TYPE> list; 
	if (list.empty())
	{
		list["FULL"]=TFCT_FULL; 
		list["DIFF_INTER_DIFF_INTRA"]=TFCT_FULL;  
		list["0"]=TFCT_FULL ; 
		//
		list["SAME_INTER_DIFF_INTRA"]= TFCT_SAME_INTER_DIFF_INTRA; 
		list["1"]=TFCT_SAME_INTER_DIFF_INTRA;
		//
		list["SAME_INTER_SAME_INTRA"]=TFCT_SAME_INTER_SAME_INTRA;
		list["2"]=TFCT_SAME_INTER_SAME_INTRA ;
	}
	return list; 
}
//	---------------------------------------------------------------------------------------------
// FIXMEFRED: mig.vc8 (28/05/2007 10:23:20):explicit specialization
template<>
const std::map<std::string,qCDS_ADJ>& 
ICM_Traits<qCDS_ADJ>::getList()
{
	static std::map<std::string,qCDS_ADJ> adjList; 
	if (adjList.empty())
	{
		adjList["NONE"]=qCredit_Default;
		adjList["STDCDS"]=qCredit_Adjust20;
		adjList["STDINDEX"]=qCredit_STDINDEX;
		adjList["STDCDSQN"]=qCredit_Adjust20N;
		adjList["STDCDSSA"]=qCredit_Adjust20SA;
		//	index series 
		adjList["STDIDX"]=qCredit_CDSDIND;
		adjList["STDTRX"]=qCredit_CDSDTRX;
		adjList["STDIDZ"]=qCredit_CDSINDZ;
		adjList["STDIDM6"]=qCredit_CDSINDM6;
		adjList["STDIDZ6"]=qCredit_CDSINDZ6;
		adjList["STDIDM7"]=qCredit_CDSINDM7;
		adjList["STDIDZ7"]=qCredit_CDSINDZ7;
		adjList["STDIDM8"]=qCredit_CDSINDM8;
		adjList["STDIDZ8"]=qCredit_CDSINDZ8;
		//	for standardization 
		adjList["STDCDSQ"]=qCredit_Adjust20;
		adjList["STDIDH5"]=qCredit_CDSDTRX;
		adjList["STDIDM5"]=qCredit_CDSDIND;
		adjList["STDIDZ5"]=qCredit_CDSINDZ;
		// summit codif
		adjList["RAW"]=qCredit_Default;
		adjList["CCDS"]=qCredit_Adjust20;
		adjList["CDSDRULE"]=qCredit_Adjust20;
		adjList["CDSDIND"]=qCredit_CDSDIND;
		adjList["CDSDTRX"]=qCredit_CDSDTRX;
		adjList["CDSDIDZ"]=qCredit_CDSINDZ;
		// new former series 
		adjList["STDIDM4"]=qCredit_CDSINDM4;
		adjList["STDIDU4"]=qCredit_CDSINDU4;
		adjList["STDIDZ4"]=qCredit_CDSINDZ4;
	}
	return adjList; 
}
//	---------------------------------------------------------------------------------------------
// FIXMEFRED: mig.vc8 (28/05/2007 10:23:20):explicit specialization
template<>
const std::map<std::string,qSENSITIVITY_TYPE>& 
ICM_Traits<qSENSITIVITY_TYPE>::getList()
{
	static std::map<std::string,qSENSITIVITY_TYPE> adjList; 
	if (adjList.empty())
	{
		adjList["SPREAD"]=ICMSPREAD_TYPE;
		adjList["IRCURVE"]=ICMIRCURVE_TYPE;
		adjList["RECOVERY"]=ICMRECOVERY_TYPE;
		adjList["CORRELATION"]=ICMCORRELATION_TYPE;
		adjList["BETA"]=ICMBETA_TYPE;
		adjList["CORR_BARY"]=ICMRECOVERY_BAR_TYPE;
		adjList["ISSUER_DEFAULT"]=ICM_ISSUER_DEFAULT;
		adjList["xFAST_SPREAD"]=xICM_FAST_SPREAD_TYPE;
		adjList["SAME_CORRELATION"]=ICM_SAMECORRELATION_TYPE;
		adjList["SAME_BETA"]=ICM_SAMEBETA_TYPE;
		adjList["REL_SHIFTSP"]=ICMSPRELSHIFT_TYPE;
		adjList["IRCURVE_ONLY"]=ICMIRCURVE_WITHOUT_DEFCURVE_TYPE;
		adjList["IRCURVE_WITH_CPN"]=ICM_IRCURVE_WITH_CPN;
		adjList["BETA_WITH_SPREAD_SHIFT"]=ICMBETA_WITH_SPREAD_SHIFT;
		adjList["THETA"]=ICM_THETA_Type;
		adjList["xISSUER_FAST_DEFAULT"]=xICM_FAST_RECOVERY_TYPE;
		adjList["CORREL_STR_DOWN"]=ICMCORREL_STRIKE_DOWN_TYPE;
		adjList["CORREL_STR_UP"]=ICMCORREL_STRIKE_UP_TYPE;
		adjList["DTR"]=ICM_DTR_TYPE;
		adjList["GREEK_DELTA"]=ICM_GREEK_DELTA_TYPE;
		adjList["GREEK_GAMMA"]=ICM_GREEK_GAMMA_TYPE;
		adjList["GREEK_VEGA"]=ICM_GREEK_VEGA_TYPE;
		adjList["GREEK_VEGA_ATM"]=ICM_GREEK_VEGA_ATM_TYPE;
		adjList["GREEK_IRVEGA_ATM"]=ICM_GREEK_IRVEGA_ATM_TYPE;
		adjList["GREEK_RHO"]=ICM_GREEK_RHO_TYPE;
		adjList["CORR_INTER_UP"]=ICM_CORREL_BET_CDO_UP_TYPE;
		adjList["CORR_INTER_DW"]=ICM_CORREL_BET_CDO_DW_TYPE;
		adjList["CPN_INFLATION"]=ICM_INFLATION_CPN_CURVE;
		adjList["CPN_IRCURVE"]=ICM_INTEREST_CPN_CURVE;
		adjList["SPREAD_RESCALING"]=ICM_INDX_SPREAD_RESCALING;
		adjList["SPREAD_PARALLEL"]=ICM_CM_SPREAD_PARALLEL;
		adjList["DEFAULT"]=ICM_CM_DEFAULT;
		adjList["RECOVERY_SENS"]=ICM_CM_RECOVERY_SENS;
		adjList["RECOVERY_LOSS"]=ICM_CM_RECOVERY_LOSS;
		adjList["RECOVERY_GLOBAL"]=ICM_CM_RECOVERY_GLOBAL;
		adjList["MC_CORRELATION"]=ICM_CM_CORRELATION;
		adjList["SUBORDINATION"]=ICM_SUBORDINATION;
	}
	return adjList; 
}
//	---------------------------------------------------------------------------------------------
/** 
const std::map<std::string,qINTERPOL_TYPE>& 
ICM_Traits<qINTERPOL_TYPE>::getList()
{
	static std::map<std::string,qINTERPOL_TYPE> adjList; 
	if (adjList.empty())
	{
		adjList["CONSTANT"]=qINTERPOL_CONSTANT;
		adjList["LINEAR"]=qINTERPOL_LINEAR;
	}
	return adjList; 
}
**/ 
//	---------------------------------------------------------------------------------------------
// FIXMEFRED: mig.vc8 (28/05/2007 10:23:20):explicit specialization
template<>
const std::map<std::string,qPAYMENT_PREMIUM_LEG>& 
ICM_Traits<qPAYMENT_PREMIUM_LEG>::getList()
{
	static std::map<std::string,qPAYMENT_PREMIUM_LEG> adjList; 
	if (adjList.empty())
	{
		adjList["CONT"]=qCONTINUE_TO_MATURITY;
		adjList["ACC_NOT_SET"]=qACCRUED_NOT_SETTLED;
		adjList["ACC"]=qACCRUED_SETTLED;
		adjList["COMP"]=qCOMPLETE_CURRENT_PERIOD;
	}
	return adjList; 
}
//	---------------------------------------------------------------------------------------------
// FIXMEFRED: mig.vc8 (28/05/2007 10:23:20):explicit specialization
template<>
const std::map<std::string,qDEFCURVE_CALIB_ALGO>& 
ICM_Traits<qDEFCURVE_CALIB_ALGO>::getList()
{
	static std::map<std::string,qDEFCURVE_CALIB_ALGO> adjList; 
	if (adjList.empty())
	{
		adjList["DICHO"]=qDEFCURVE_DICHO;
		adjList["NEWTON"]=qDEFCURVE_NEWTON;
		adjList["BRENT"]=qDEFCURVE_BRENT;
	}
	return adjList; 
}

// ----------------------------------------------------------------------------------------------
const std::map<std::string,qRAN_GEN>& 
ICM_Traits<qRAN_GEN>::getList()
{
	static std::map<std::string,qRAN_GEN> list; 
	if (list.empty())
	{
		list["NAG"]		=q_NAG; 
		list["RAN1"]	=q_RAN1;
		list["RAN2"]	=q_RAN2; 
		list["DEF"]		=q_DEF;
		list["RANMAR"]	=q_RANMAR;
		list["RNG_STR"]	=q_RNG_STR;
		list["KISS"]	=q_KISS;
		list["INV_CUM_NORM_ACKLAM"] =q_INV_CUM_NORM_ACKLAM;
		list["INV_CUM_NORM_MORO"]	=q_INV_CUM_NORM_MORO;
	}
	return list; 
}
//	---------------------------------------------------------------------------------------------
void 
ICM_EnumsCnv::cnv(const std::string& name,qTWO_FACTORS_CORRELATION_TYPE& ref,bool& ret,std::string& list)
{
	ret=true; 
	const ICM_Traits<qTWO_FACTORS_CORRELATION_TYPE>::list_t& items = ICM_Traits<qTWO_FACTORS_CORRELATION_TYPE>::getList(); 
	ICM_Traits<qTWO_FACTORS_CORRELATION_TYPE>::list_t::const_iterator it = 
		items.find(name) ; 
	if (it==items.end()) 
	{
		ret=false;  
		list=""; 
		for(it=items.begin();it!=items.end();++it) 
		{ list += it->first ; list+=" " ; }
		return ; 
	}
	ref=it->second; 
}
//	---------------------------------------------------------------------------------------------
void 
ICM_EnumsCnv::cnv(const std::string& name,qVECTMETH& ref,bool& ret,std::string& list)
{
	ret=true; 
	const ICM_Traits<qVECTMETH>::list_t& items = ICM_Traits<qVECTMETH>::getList(); 
	ICM_Traits<qVECTMETH>::list_t::const_iterator it = 
		items.find(name) ; 
	if (it==items.end()) 
	{
		ret=false;  
		list=""; 
		for(it=items.begin();it!=items.end();++it) 
		{ list += it->first ; list+=" " ; }
		return ; 
	}
	ref=it->second; 
}
//	---------------------------------------------------------------------------------------------
void 
ICM_EnumsCnv::cnv(const std::string& name,qCMPMETH& ref,bool& ret,std::string& list)
{
	ret=true; 
	const ICM_Traits<qCMPMETH>::list_t& items = ICM_Traits<qCMPMETH>::getList(); 
	ICM_Traits<qCMPMETH>::list_t::const_iterator it = 
		items.find(name) ; 
	if (it==items.end()) 
	{
		ret=false ;
		list="NPV,SPREAD, FEELEG, DEFLEG,ACCRUED,ACCRUED_PV,PREMIUM,FWDSPREAD,DURATION,CORREL_DOWN,CORREL_UP, AVGCORRDEF, FEES, EL, OPT_DELTA, OPT_GAMMA, OPT_VEGA";
		return ; 
	}
	ref=it->second; 
}

//	---------------------------------------------------------------------------------------------
void 
ICM_EnumsCnv::toString(qCDS_ADJ ref,std::string& name) 
{
	switch (ref) 
	{
	/*case qCredit_Special_None_Date :	name="qCredit_Special_None_Date" ; return ;
	case qCredit_Special_None_YF:		name="qCredit_Special_None_Date" ;	 return ;*/
    case qCredit_Default:				name="NONE" ;	return ;		
	case qCredit_Adjust20:				name="STDCDS" ;		return ;
	case qCredit_CDSDTRX:				name="STDTRX" ;return ;
	case qCredit_CDSDIND:				name="STDIDX" ;return ;
	case qCredit_CDSINDZ:				name="STDIDZ" ;return ;
	case qCredit_Adjust20N:				name="STDCDSQN" ;return ;
	case qCredit_Adjust20SA:			name="STDCDSSA" ;return ;
	case qCredit_STDINDEX:				name="STDINDEX" ;return ;
	case qCredit_CDSINDH5:				name="STDIDH5" ;return ;
	case qCredit_CDSINDM5:				name="STDIDM5" ;return ;
	case qCredit_CDSINDZ5:				name="STDIDZ5" ;return ;
	case qCredit_CDSINDM6:				name="STDIDM6" ;return ;
	case qCredit_CDSINDZ6:				name="STDIDZ6" ;return ;
	case qCredit_CDSINDM7:				name="STDIDM7" ;return ;
	case qCredit_CDSINDZ7:				name="STDIDZ7" ;return ;
	case qCredit_CDSINDM8:				name="STDIDM8" ;return ;
	case qCredit_CDSINDZ8:				name="STDIDZ8" ;return ;
	case qCredit_CDSINDM4:				name="STDIDM4" ;return ;
	case qCredit_CDSINDU4:				name="STDIDU4" ;return ;
	case qCredit_CDSINDZ4:				name="STDIDZ4" ;return ;

	} ;
	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_EnumsCnv::toString: undefined "<<ref<<" for qCDS_ADJ"); 
}
//	---------------------------------------------------------------------------------------------
void ICM_EnumsCnv::toString(qCMPMETH ref ,std::string& name)
{
  switch (ref) 
	{
   case qCMPPRICE  :		name="NPV"; return; 				//Net Present Value
   case	qCMPSPREAD :		name="SPREAD"; return; 				//Break even Spread
   case	qCMPFEELEGPV :		name="FEELEG"; return; 			//Fee Leg Present value
   case	qCMPDEFLEGPV :		name="DEFLEG"; return; 			//Default Leg Present value 
   case qCMPDURATION :		name="DURATION"; return; 			//Duration 
   case qCMPACCRUED :		name="ACCRUED"; return; 			//Standard Accrued
   case	qCMPPREMIUM :		name="PPREMIUM"; return; 			//(JLA) Premium for option products
   case	qCMPFWDSPREAD :		name="FWDSPREAD"; return; 			//(JLA)	FwdSpread for option products
   case	qCMPFEES :			name="FEES"; return; 				//Fees Pricing		
   case	qCMPCORRELUP :		name="CORREL_UP"; return; 			// for tranches
   case	qCMPCORRELDOWN :	name="CORREL_DOWN"; return; 
   case	qCMPAVGCORRDEF :	name="AVGCORRDEF"; return;
   case qCMPFLATCORR :		name="FLATCORR";return;
   case qCMPEL		:		name="EL"; return;
   case qCMP_OPT_AN_DELTA :		name="OPT_DELTA";return;
   case qCMP_OPT_AN_VEGA :		name="OPT_VEGA";return;
   case qCMP_OPT_AN_GAMMA :		name="OPT_GAMMA";return;
	};
	name = "";
	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_EnumsCnv::toString: undefined "<<ref<<" for qCMPMETH"); 
}
//	---------------------------------------------------------------------------------------------
void ICM_EnumsCnv::cnv(const int& value ,qRescalType& type ,bool& ret)
{
	ret=true; 
	switch (value)
	{
		case -10:	{type = qRescal_Eqty_Maturity;break;}
		case -102: 	{type = qRescal_Eqty;break;}
		case -103: 	{type = qRescal_Eqty_NormByEL_Maturity;break;}
		case -104: 	{type = qRescal_Eqty_Mixed_Maturity;break;}
		case -105: 	{type = qRescal_Eqty_Mixed_NormByEL_Maturity;break;}
		case -106:	{type = qRescal_Eqty_Digital_Maturity;break;}
		case -11: 	{type = qRescal_Eloss_Maturity;break;}
		case -112:	{type = qRescal_Eloss;break;}
		case -12:	{type = qRescal_Std;break;}
		case -136:	{type = qRescal_ING_Maturity;break;}
		default :	{type = qRescal_Std_Maturity;break;}
	}
	return; 
}
//	---------------------------------------------------------------------------------------------
/**
void ICM_EnumsCnv::cnv(const std::string& name,qINTERPOL_TYPE&ref,bool&ret,std::string&list)
{
	ret=true; 
	const ICM_Traits<qINTERPOL_TYPE>::list_t& items = ICM_Traits<qINTERPOL_TYPE>::getList(); 
	ICM_Traits<qINTERPOL_TYPE>::list_t::const_iterator it = 
		items.find(name) ; 
	if (it==items.end()) 
	{
		ret=false ;
		list="CONSTANT,LINEAR";
		return ; 
	}
	ref=it->second; 
}
**/ 
//	---------------------------------------------------------------------------------------------
void ICM_EnumsCnv::toString(const qRescalType type ,std::string& restype )
{
	switch (type)
	{
		case qRescal_Eqty_Maturity :				{restype = "Rescal_Eqty_Maturity";  break;}
		case qRescal_Eqty :							{restype = "Rescal_Eqty";  break;}
		case qRescal_Eqty_NormByEL_Maturity :		{restype = "Rescal_Eqty_NormByEL_Maturity";  break;}
		case qRescal_Eqty_Mixed_Maturity :			{restype = "Rescal_Eqty_Mixed_Maturity";  break;}
		case qRescal_Eqty_Mixed_NormByEL_Maturity :	{restype = "Rescal_Eqty_Mixed_NormByEL_Maturity";  break;}
		case qRescal_Eloss_Maturity :				{restype = "Rescal_Eloss_Maturity";  break;}
		case qRescal_Eloss :						{restype = "Rescal_Eloss";  break;}
		case qRescal_Std :							{restype = "Rescal_Std";  break;}
		case qRescal_Std_Maturity :					{restype = "Rescal_Std_Maturity";  break;}
		case qRescal_ING_Maturity :					{restype = "Rescal_ING_Maturity";  break;}
	}
}
//	---------------------------------------------------------------------------------------------
void ICM_EnumsCnv::toInt(const qRescalType type ,int& restype) 
{
	switch (type)
	{
		case qRescal_Eqty_Maturity :				{restype = -10;  break;}
		case qRescal_Eqty :							{restype = -102;  break;}
		case qRescal_Eqty_NormByEL_Maturity :		{restype = -103;  break;}
		case qRescal_Eqty_Mixed_Maturity :			{restype = -104;  break;}
		case qRescal_Eqty_Mixed_NormByEL_Maturity :	{restype = -105;  break;}
		case qRescal_Eloss_Maturity :				{restype = -11;  break;}
		case qRescal_Eloss :						{restype = -112;  break;}
		case qRescal_Std :							{restype = -12;  break;}
		case qRescal_ING_Maturity :					{restype = -136;  break;}
		case qRescal_Std_Maturity :					{restype = 0;  break;}
	}
}
//	---------------------------------------------------------------------------------------------


void ICM_EnumsCnv::cnv(const int& value ,qTERM_STRUCTURE& type ,bool& ret)
{
	ret=true; 
	switch (value)
	{
		case 0:		{type = qNoTermStructure;break;}
		case 1: 	{type = qTermStructure;break;}
		case 2: 	{type = qTermStructureR;break;}
		default :	{type = qNoTermStructure;break;}
	}
	return; 
}

void 
ICM_EnumsCnv::cnv(const std::string& name,qRAN_GEN& ref,bool& ret,std::string& list)
{
	ret=true; 
	const ICM_Traits<qRAN_GEN>::list_t& items = ICM_Traits<qRAN_GEN>::getList(); 
	ICM_Traits<qRAN_GEN>::list_t::const_iterator it = 
		items.find(name) ; 
	if (it==items.end()) 
	{
		ret=false;  
		list=""; 
		for(it=items.begin();it!=items.end();++it) 
		{ list += it->first ; list+=" " ; }
		return ; 
	}
	ref=it->second; 
}

//	---------------------------------------------------------------------------------------------
void ICM_EnumsCnv::toString(const qRAN_GEN& type ,std::string& RandomName )
{
	switch (type)
	{
		case q_NAG :		{RandomName = "NAG";  break;}
		case q_RAN1 :		{RandomName = "RAN1";  break;}
		case q_RAN2 :		{RandomName = "RAN2";  break;}
		case q_DEF :		{RandomName = "DEF";  break;}
		case q_RANMAR :		{RandomName = "RANMAR";  break;}
		case q_RNG_STR :	{RandomName = "RNG_STR";  break;}
		case q_KISS :		{RandomName = "KISS";  break;}
		case q_INV_CUM_NORM_ACKLAM :	{RandomName = "INV_CUM_NORM_ACKLAM";  break;}
		case q_INV_CUM_NORM_MORO :		{RandomName = "INV_CUM_NORM_MORO";  break;}
	}
}