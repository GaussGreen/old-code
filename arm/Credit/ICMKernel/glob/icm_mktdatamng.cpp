#include "ARMKernel\glob\firsttoinc.h"
#include "ICMKernel/util/icm_pch.h"
#include "ICMKernel/glob/icm_mktdatamng.h"
#include "ICMKernel/glob/icm_smile_correlation.h"
#include "ICMKernel/crv/icm_implied_loss_tree.h"
#include "ARMKernel/util/fromto.h"
#include "gpinflation/infcurv.h"
#include "ICMKernel/crv/icm_fixing_curve.h"
#include "ARMKernel/inst/irindex.h"
#include "ICMKernel/util/icm_utils.h"

using namespace std;
using namespace ARM; 

string GetZeroCurveName(const string& ccy,const ARM_Date& AsOf)
{
	string requete="ZRC:";
	char* tmp = unconst(AsOf).GetStrDate(); 
	string asof = tmp; delete [] tmp ;
	requete += ccy;
	requete += "_";
	// char* index = new char[10];
	char index[10]; 
	GetSummitCurveIndexFromCurrency((char*)ccy.c_str(),index);
	requete += (string)index;
	requete += "_";
	requete += asof;
	// delete[] index;

	return (requete);
}

string GetInflaCurveName(const const string& label,const ARM_Date& AsOf)
{
	string requete="INF:";
	char* tmp = unconst(AsOf).GetStrDate(); 
	string asof=tmp; delete [] tmp ;
	requete += (string)label;
	requete += "_";
	requete += asof;

	return (requete);
}

string GetDefaultCurveName(const string& label,const ARM_Date& AsOf)
{
	string requete="DEF:";
	char* tmp = unconst(AsOf).GetStrDate(); 
	string asof=tmp; delete [] tmp ;
	requete += (string)label;
	requete += "_";
	requete += asof;

	return (requete);
}

string GetVolCurveName(const string& name,const string& ccy,const ARM_Date& AsOf)
{
	string requete="VOL:";
	char* tmp=unconst(AsOf).GetStrDate(); 
	string asof =tmp; delete [] tmp ;
	requete += name;
	requete += "_";
	requete += ccy;
	requete += "_";
	requete += asof;

	return (requete);
}

string GetCorrelationName(const string& structName,const string& index1, const string& index2, const ARM_Date& AsOf)
{
	string requete="COR:";
	char* tmp = unconst(AsOf).GetStrDate(); 
	string asof = tmp; delete [] tmp ;
	string label = GetCorrelationIndex(structName,index1, index2);
	requete += label;
	requete += "_";
	requete += asof;

	return (requete);
}



string GetLossTreeName(const string& label,const ARM_Date& AsOf)
{
	string requete="TREE:";
	char* tmp = unconst(AsOf).GetStrDate(); 
	string asof = tmp; delete [] tmp ;
	requete += label;
	requete += "_";
	requete += asof;

	return (requete);
}

std::string GetFixingCurveName(const std::string&indexName,const std::string&ccy)
{
	string requete="FIX:" ;
	requete += indexName; 
	requete += "_"; 
	requete += ccy; 
	return requete; 
}


void deducelabelforobject(const ARM_Object&object_,string& label)
{
	ARM_Object* object= const_cast<ARM_Object*>(&object_); 
	ARM_CLASS_NAME name = object->GetName();
	char TMP[1000];

	switch (name)
	{
	case ICM_FIXING_CURVE: 
	{
		std::string name; 
		const ICM_Fixing_Curve& curve = dynamic_cast<const ICM_Fixing_Curve&>(object_); 
		if (curve.hasIndex()) 
		{
			name = GetIndexName(curve.getIndex()) ;
			name +="_" ;
			name += curve.getIndex().GetCurrencyUnit()->GetCcyName(); 
		}
		else name = curve.GetIndexName(); 
		label="FIX:"+name; 
		return;
	}
	break ;
	case ICM_DEFAULTCURVE:
	case ICM_CST_PIECEWISE:
	case ICM_LINEAR_PIECEWISE:
	case ICM_DEFAULT_CURVE_MODEL:
	{//DEFPROB_MO_20050302
	ICM_DefaultCurve* curve = (ICM_DefaultCurve*) object;
	string plabel = curve->GetLabel();
	char* date = curve->GetAsOfDate().GetStrDate();
	sprintf(TMP,"DEF:%s_%s",plabel.c_str(),date);
	delete[] date;
	}
	break;
	case ARM_ZERO_CURVE:
	case ARM_ZERO_INTERPOLATION:
	case ARM_ZERO_LIN_INTERPOL:
	case ARM_ZERO_CUBDIFF:
	case ARM_ZERO_SPLINES:
	case ARM_ZERO_SPLICUB:
	case ARM_ZERO_VASICEK:
	case ARM_ZERO_FLAT:
	{//EUR_EURIB_MO_20050302
	ARM_ZeroCurve* curve = (ARM_ZeroCurve*) object;
	char* ccy = curve->GetCurrencyUnit()->GetCcyName();
	// char* index = new char[10];
	char index[10]; 
	char* date = curve->GetAsOfDate().GetStrDate();
	GetSummitCurveIndexFromCurrency(ccy,index);
	sprintf(TMP,"ZRC:%s_%s_%s",ccy,index,date);
	// delete[] index;
	delete[] date;
	}
	break;
	case ARM_SMILE_CURVE:
	case ARM_VOL_CURVE:
	case ARM_VOL_FLAT:
	case ARM_VOL_LIN_INTERPOL:
	case ARM_VOL_CUBE:
	{//EUR_MO_20050302
	ARM_VolCurve* curve = (ARM_VolCurve*) object;
	std::string name; 
	if (curve->hasIndex()) name = GetIndexName(curve->GetIndex()) ;
	else name = curve->GetIndexName(); 
	
	char* ccy = curve->GetCurrency()->GetCcyName();
	char* date = curve->GetAsOfDate().GetStrDate();
	sprintf(TMP,"VOL:%s_%s_%s",name.c_str(),ccy,date);
	delete[] date;
	}
	break;
	case ICM_FLAT_CORRELATION:
	case ICM_CORRMATRIX:
	case ICM_SMILE_CORRMATRIX:
	case ICM_BETA_CORRMATRIX:
	{//LABEL_MO_20050302
	// ICM_Smile_Correlation* curve = (ICM_Smile_Correlation*) object;
	ICM_Correlation * curve = dynamic_cast<ICM_Correlation*>(object);
	string structname = curve->GetStructName();
	char* date = curve->GetAsOfDate().GetStrDate();
	std::string name1("") ; if (curve->GetIndex1()) name1=GetIndexName(*curve->GetIndex1()); 
	std::string name2("") ; if (curve->GetIndex2()) name2=GetIndexName(*curve->GetIndex2()); 
	std::string name = GetCorrelationIndex(structname, name1, name2);
	sprintf(TMP,"COR:%s_%s",name.c_str(),date);
	delete[] date;
	}
	break;
	case ICM_IMPLIED_LOSS_TREE :
	{//LABEL_MO_20050302
	ICM_ImpLossTree* tree = (ICM_ImpLossTree*) object;
	ICM_Smile_Correlation* correl = tree->GetCorrelation();
	string structname = correl->GetStructName();
	char* date = correl->GetAsOfDate().GetStrDate();
	sprintf(TMP,"TREE:%s_%s",structname.c_str(),date);
	delete[] date;
	}
	break;
	case ARM_INFCURV :
	{ARM_InfCurv* infla = dynamic_cast<ARM_InfCurv*>(object);
	string structname = infla->GetInfIdxName();
	char * date = infla->GetAsOf().GetStrDate();
	sprintf(TMP,"INF:%s_%s",structname.c_str(), date);
	delete[] date;
	}
	break;
	default:
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_MktDataMng::GenerateShiftMng : bad market data for IR sensi");  	
	}
	
	label = (string) TMP;
}

ARM_ZeroCurve* ICM_MktDataMng::GetZeroCurve(const string& ccy,const ARM_Date& AsOf) const 
{
	string requete= GetZeroCurveName(ccy,AsOf);
	return dynamic_cast<ARM_ZeroCurve*>(find(requete));
}

ARM_InfCurv*   ICM_MktDataMng::GetInflaCurve(const const string& label,const ARM_Date& AsOf)const 
{
	string requete= GetInflaCurveName(label,AsOf);
	return dynamic_cast<ARM_InfCurv*>(find(requete));	
}

ICM_DefaultCurve* ICM_MktDataMng::GetDefaultCurve(const string& label,const ARM_Date& AsOf) const 
{
	string requete=GetDefaultCurveName(label,AsOf);
	return dynamic_cast<ICM_DefaultCurve*>(find(requete));
}

ARM_VolCurve* ICM_MktDataMng::GetVolCurve(const string& name,const string& ccy,const ARM_Date& AsOf) const 
{
	string requete=GetVolCurveName(name,ccy,AsOf);
	return dynamic_cast<ARM_VolCurve*>(find(requete));
}

ICM_Correlation* ICM_MktDataMng::GetCorrelation(const string& label,const std::string& index1, const std::string& index2,const ARM_Date& AsOf) const 
{
	string requete=GetCorrelationName(label,index1,index2,AsOf);
	return dynamic_cast<ICM_Correlation*>(find(requete));
}

ICM_Correlation* ICM_MktDataMng::GetCorrelation(const string& label,const ARM_Date& AsOf) const 
{
	string requete=GetCorrelationName(label,"","",AsOf);
	return dynamic_cast<ICM_Correlation*>(find(requete));
}

ICM_ImpLossTree* ICM_MktDataMng::GetLossTree(const string& label,const ARM_Date& AsOf) const 
{
	string requete=GetLossTreeName(label,AsOf);
	return dynamic_cast<ICM_ImpLossTree*>(find(requete));
}

ICM_Correlation* ICM_MktDataMng::GetUnicCorrelation() const 
{
	string requete(CORREL);
	return dynamic_cast<ICM_Correlation*>(find(requete));
}

ARM_ZeroCurve* ICM_MktDataMng::GetUnicZeroCurve() const 
{
	string requete(IRCURVE);
	return dynamic_cast<ARM_ZeroCurve*>(find(requete));
}

ARM_InfCurv*   ICM_MktDataMng::GetUnicInflaCurve() const
{
	string requete(INFCURV);
	return dynamic_cast<ARM_InfCurv*>(find(requete));
}


/**ARM_ZeroCurve* ICM_MktDataMng::GetZeroCurveFromKey(const string& Key) const {
		return dynamic_cast<ARM_ZeroCurve*>(find(Key));
}
ARM_InfCurv*   ICM_MktDataMng::GetInflaCurveFromKey(const string& Key) const {
		return dynamic_cast<ARM_InfCurv*>(find(Key));	
}
ICM_DefaultCurve* ICM_MktDataMng::GetDefaultCurveFromKey(const string& Key) const {
		return dynamic_cast<ICM_DefaultCurve*>(find(Key));
}
ARM_VolCurve* ICM_MktDataMng::GetVolCurveFromKey(const string& Key) const {
		return dynamic_cast<ARM_VolCurve*>(find(Key));
}
ICM_Correlation* ICM_MktDataMng::GetCorrelationFromKey(const string& Key) const {
		return dynamic_cast<ICM_Correlation*>(find(Key));
}
ICM_ImpLossTree* ICM_MktDataMng::GetLossTreeFromKey(const string& Key) const {
		return dynamic_cast<ICM_ImpLossTree*>(find(Key));
}
**/ 
// -------------------------------------------------------------------------
// Bump a package of Default Curves according a given bump profile
// -------------------------------------------------------------------------
ICM_MktDataMng* ICM_MktDataMng::GenerateShiftMng(ICM_BUMP<ICM_BASIS_BUMP>& parameters)
{
    double sensitivity =0.;
	int NoCurve = -1,i=0,h=0,k=0;
	bool res = false;
	ARM_Date AsOf;

	char TermZC[ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS];
	ARM_Vector epsilon(1, parameters.itsEpsilon) ;

	ICM_MktDataMng* mng = (ICM_MktDataMng*) Clone();
	ARM_Object* Obj=NULL;
	ICM_DefaultCurve* DefCurve = NULL;
	ICM_Correlation* correlation = NULL;
	ARM_ZeroCurve* InitShortCurve=NULL; 
	ARM_VolCurve* VolCurve = NULL;
	
    
		for (k=0;k<parameters.itsBumps.size();k++)
		{
			// ---------------------------------------------------------------------------
			// Bump Recovery / IrCurve / Spread 
			// ---------------------------------------------------------------------------
			switch (parameters.itsSensiType)
			{
				case ICM_DTR_TYPE :
				{
					std::string name = GetDefaultCurveName(parameters.itsBumps[k].itsLabel,AsOf); 
					Obj = mng->find(name);
					parameters.itsBumps[k].GetTenors(TermZC,epsilon);
					vector<string> vTerm;
					for ( int i = 0; i< ARM_NB_TERMS; i++)  {
						if (! strcmp(vTerm[i].c_str(), "X")) vTerm.push_back(TermZC[i]);
					}
					DefCurve = ((ICM_DefaultCurve*)Obj)->GenerateShiftCurve(vTerm,epsilon,ICM_DTR_TYPE);
					DefCurve->SetLabel(((ICM_DefaultCurve*)Obj)->GetLabel());
					mng->erase(name) ;
					mng->adopt(DefCurve,name); 
					// mng->replace(GetDefaultCurveName(parameters.itsBumps[k].itsLabel,AsOf),DefCurve);
				}
				break;
				case ICMRECOVERY_TYPE :
				{
					std::string name = GetDefaultCurveName(parameters.itsBumps[k].itsLabel,AsOf); 
					Obj = mng->find(name);
					if (parameters.itsBumps[k].IsParallelShift())
					{
						vector<string> vTerm;
						vTerm.push_back("ALL"); 
						ARM_Vector epsilon(1); 
						epsilon.Elt(0) =parameters.itsBumps[k].itsEpsilon ;
						DefCurve = ((ICM_DefaultCurve*)Obj)->GenerateShiftCurve(vTerm,epsilon,ICMRECOVERY_TYPE);
					}
					else
					{
						parameters.itsBumps[k].GetTenors(TermZC,epsilon);
						vector<string> vTerm;
						for ( int i = 0; i< ARM_NB_TERMS; i++)  {
							if (! strcmp(vTerm[i].c_str(), "X")) vTerm.push_back(TermZC[i]);
						}
						DefCurve = ((ICM_DefaultCurve*)Obj)->GenerateShiftCurve(vTerm,epsilon,ICMRECOVERY_TYPE);
					}
					DefCurve->SetLabel(((ICM_DefaultCurve*)Obj)->GetLabel());
					mng->erase(name) ;
					mng->adopt(DefCurve,name); 
					// mng->replace(GetDefaultCurveName(parameters.itsBumps[k].itsLabel,AsOf),DefCurve);

				}
				break;
				case ICMIRCURVE_TYPE :
				case ICMIRCURVE_WITHOUT_DEFCURVE_TYPE :
				{
					std::string name(GetZeroCurveName(parameters.itsBumps[k].itsLabel,AsOf)); 
					Obj = mng->find(name);
					ARM_ZeroCurve* ModifShortCurve = (ARM_ZeroCurve*) Obj;
					if (parameters.itsBumps[k].IsParallelShift())
					{
						ModifShortCurve = (ARM_ZeroCurve*)ModifShortCurve->Clone() ; 
						ModifShortCurve->ParallelShift(parameters.itsBumps[k].itsEpsilon);
					}
					else
					{
						InitShortCurve = ModifShortCurve;
						parameters.itsBumps[k].GetTenors(TermZC,epsilon);
						ModifShortCurve = (ARM_ZeroCurve*) InitShortCurve->GenerateShiftCurve(TermZC,&epsilon);
						// if (InitShortCurve) delete InitShortCurve;
					}
					// mng->replace(GetDefaultCurveName(parameters.itsBumps[k].itsLabel,AsOf),ModifShortCurve);
					mng->erase(name) ;
					mng->adopt(ModifShortCurve,name); 
				}
				break;
				case ICMSPREAD_TYPE :
				{
					std::string name(GetDefaultCurveName(parameters.itsBumps[k].itsLabel,AsOf));
					Obj = mng->find(name);
					if (parameters.itsBumps[k].IsParallelShift())
					{
						vector<string> vTerm;
						vTerm.push_back("ALL"); 
						ARM_Vector epsilon(1); 
						epsilon.Elt(0) =parameters.itsBumps[k].itsEpsilon ;
						DefCurve = ((ICM_DefaultCurve*)Obj)->GenerateShiftCurve(vTerm,epsilon,ICMSPREAD_TYPE);
					}
					else
					{
						parameters.itsBumps[k].GetTenors(TermZC,epsilon);
						vector<string> vTerm;
						for ( int i = 0; i< ARM_NB_TERMS; i++)  {
							if (! strcmp(vTerm[i].c_str(), "X")) vTerm.push_back(TermZC[i]);
						}
						DefCurve = ((ICM_DefaultCurve*)Obj)->GenerateShiftCurve(vTerm,epsilon,ICMSPREAD_TYPE);
					}
					DefCurve->SetLabel(((ICM_DefaultCurve*)Obj)->GetLabel());
					mng->erase(name) ;
					mng->adopt(DefCurve,name); 
					// mng->replace(GetDefaultCurveName(parameters.itsBumps[k].itsLabel,AsOf),DefCurve);
				
				}
				break;
				case ICM_GREEK_VEGA_TYPE :
				{

					Obj = mng->find(GetDefaultCurveName(parameters.itsBumps[k].itsLabel,AsOf));
					DefCurve = (ICM_DefaultCurve*)Obj ;

					std::string name = GetVolCurveName(parameters.itsBumps[k].itsLabel,DefCurve->GetCurrency(),AsOf) ;

					Obj = mng->find(name);
					
					VolCurve = ((ARM_VolCurve*)Obj->Clone()) ;
					VolCurve->ParallelShift(parameters.itsBumps[k].itsEpsilon);
					mng->erase(name);
					mng->adopt(VolCurve,name); 
					// mng->replace(GetVolCurveName(parameters.itsBumps[k].itsLabel,
					// 							 DefCurve->GetCurrency()->GetCcyName(),AsOf),
					// 							 VolCurve);
				
				}
				break;
				default :
				;
			}

		// ---------------------------------------------------------------------------
		// Bump Beta / Correlation
		// ---------------------------------------------------------------------------

		switch (parameters.itsSensiType)
		{
			case ICMBETA_TYPE :
			case ICMCORRELATION_TYPE :
			case ICMCORREL_STRIKE_DOWN_TYPE :
			case ICMCORREL_STRIKE_UP_TYPE :
			{
				ICMTHROW(ERR_INVALID_ARGUMENT,"JLA: not implemented"); 
			parameters.itsBumps[k].GetTenors(TermZC,epsilon);
			correlation = ((ICM_Correlation*)Obj)->GenerateShiftBetas(parameters.itsBumps[k].itsLabel,ICMSPREAD_TYPE,parameters.itsBumps[k].itsEpsilon);
			// mng->replace(GetDefaultCurveName(parameters.itsBumps[k].itsLabel,AsOf),DefCurve);
			}
			break;
			default :;
			}
		}
				
	
    

	return (mng);
}

void ICM_MktDataMng::GetMktDataVectotKey(std::vector<std::string>& mktdata) const 
{
	for ( MapMgn::const_iterator it = its_Manager.begin(); it != its_Manager.end(); it++){
		mktdata.push_back(it->first);
	}
}


ICM_MktDataMng* ICM_MktDataMng::GenerateShiftMng(qSENSITIVITY_TYPE typesensi,
										const std::string& plot, 
										const string& label,
										double epsvalue)
{
	double sensitivity =0.;
	int NoCurve = -1;
	int i=0,h=0;
	double result = 0.0;
	bool Parallelshift = false ;
	vector<string> Term;
	char TermZC[ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS];
	memset(TermZC,'\0',sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);
	// ARM_Vector* epsilon = NULL;
	vector<string> labelKey(1);
	labelKey[0] = label;

	double EPSL = 0.01;		//Bump sur les taux fixé à 1bp	
	double EPSL_DEF = 0.01; //Bump sur le spread fixé à 1 bp.
	double EPSL_REC = 0.1;  //Bump fixe sur le recovery de 0.1
	double EPSL_CORR = 0.1;  //Bump fixe sur la correl de 0.1
	double EPSL_BETA = 0.1;  //Bump fixe sur la correl de 0.1
	double EPSL_SPREL = 0.1 ;//Bump relatif des spreads de 0.1
	double EPSL_DTR = CREDIT_DEFAULT_VALUE;//Bump relatif des spreads de 0.1

	if (epsvalue != CREDIT_DEFAULT_VALUE)
		EPSL_DTR = EPSL_BETA = EPSL_CORR = EPSL_REC = EPSL_DEF = EPSL = EPSL_SPREL = epsvalue;

	ICM_DefaultCurve* ModifDefCurve = NULL;
	ARM_ZeroCurve* ModifShortCurve = NULL;
	ARM_ReferenceValue* Recovery = NULL;
	ICM_Correlation* Correlation = NULL;
	ARM_InfCurv* ModifInfCurve = NULL;
	ARM_ZeroCurve* ModifCpnIRCurve = NULL;
	
	ARM_Object* Obj=NULL;
	ICM_MktDataMng* modelMkt = (ICM_MktDataMng*) Clone();

	ARM_Vector epsilon (1,epsvalue);

	// if (!strcmp(plot,"NONE"))  //Detection d'un parallel Shift
	if (plot=="NONE")
	{
		Parallelshift = true;
		Term.push_back("ALL"); 
	}
	else
	{
		Term.push_back(plot);
		strcpy(TermZC[0],plot.c_str());
		switch (typesensi)
		{
			case ICMRECOVERY_TYPE :
			{
				epsilon.InitElt(0,EPSL_REC);
			}
			break;
			case ICMCORRELATION_TYPE :
			{
				epsilon.InitElt(0,EPSL_CORR);
			}
			break;
			case ICMIRCURVE_TYPE :
			{
				epsilon.InitElt(0,EPSL);
			}
			break;
			default:
			case ICM_INDX_SPREAD_RESCALING :
			case ICMSPREAD_TYPE :
			{
				epsilon.InitElt(0,EPSL_DEF);
			}
			break;
			case ICMBETA_TYPE :
			{
				epsilon.InitElt(0,EPSL_BETA);
			}
			break;
			case ICMSPRELSHIFT_TYPE :
			{
				epsilon.InitElt(0,EPSL_SPREL);
			}
			break;
			case ICM_DTR_TYPE :
			{
				epsilon.InitElt(0,EPSL_DTR);
			}
			break;
		}
	}


	
	/** try
    {**/ 
		// ---------------------------------------------------------------------------
		// Bump Recovery / IrCurve / Spread 
		// ---------------------------------------------------------------------------
		switch (typesensi)
		{
			case ICM_DTR_TYPE :
			{			
				for ( int i=0; i< labelKey.size(); i++){
					Obj = dynamic_cast<ICM_DefaultCurve*>(modelMkt->find(labelKey[i]));
					if(Obj) {
						ModifDefCurve = ((ICM_DefaultCurve*)Obj)->GenerateShiftCurve(Term,epsilon,ICM_DTR_TYPE);
						modelMkt->adopt(ModifDefCurve);
					} else {
						ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_MktDataMng::GenerateShiftMng : bad market data for DTR sensi");  	
					}
				}
			}break;
		
			case ICMIRCURVE_WITHOUT_DEFCURVE_TYPE :
			{
				for ( int i=0; i<labelKey.size(); i++){
					Obj = dynamic_cast<ARM_ZeroCurve*>(modelMkt->find(labelKey[i]));
					if (Obj) {
						if (Parallelshift)
						{	
							ModifShortCurve = (ARM_ZeroCurve*)(Obj->Clone());
							ModifShortCurve->ParallelShift(EPSL);
						}else {
							ModifShortCurve = ((ARM_ZeroCurve*)Obj)->GenerateShiftCurve(TermZC,&epsilon);
						}
						modelMkt->adopt(ModifShortCurve);
					}else{
						ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_MktDataMng::GenerateShiftMng : bad market data for IR sensi");  	
					}
				}
			}
			break;
			//case ICMIRCURVE_TYPE : with calibration 
			case ICM_INFLATION_CPN_CURVE :
			{
				vector<double> vepsilon;
				vector<string> mktterms;
				vepsilon.resize(epsilon.GetSize());
				mktterms.resize(epsilon.GetSize());
				for (int il_ = 0; il_<epsilon.GetSize(); il_++)
				{vepsilon[il_]=epsilon.Elt(il_);
				 mktterms[il_]= (string)TermZC[il_];}

				for ( int i=0; i<labelKey.size(); i++){
					Obj = dynamic_cast<ARM_InfCurv*>(modelMkt->find(labelKey[i]));
					if (Obj) {
						ModifInfCurve = ((ARM_InfCurv*)Obj)->GenerateShiftCurve(mktterms,vepsilon);
						modelMkt->adopt(ModifInfCurve);
					}else{
						ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_MktDataMng::GenerateShiftMng : bad market data for Infla sensi");  	
					}
				}	
			}
			break;
			case ICM_IRCURVE_WITH_CPN :
			{
				for ( int i=0; i<labelKey.size(); i++){
					Obj = dynamic_cast<ARM_ZeroCurve*>(modelMkt->find(labelKey[i]));
					if (Obj) {
						if (Parallelshift)
						{	
							ModifShortCurve = (ARM_ZeroCurve*)(Obj->Clone());
							ModifShortCurve->ParallelShift(EPSL);
						}else {
							ModifShortCurve = ((ARM_ZeroCurve*)Obj)->GenerateShiftCurve(TermZC,&epsilon);
						}
						modelMkt->adopt(ModifShortCurve);
					}else{
						ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_MktDataMng::GenerateShiftMng : bad market data for IR sensi");  	
					}
				}
					
				ARM_ZeroCurve* ModifCpnIRCurve = NULL;
				for (  i=0; i<labelKey.size(); i++){
					Obj = dynamic_cast<ARM_ZeroCurve*>(modelMkt->find(labelKey[i]));
					if (Obj) {
						if (Parallelshift)
						{
							ModifCpnIRCurve = (ARM_ZeroCurve*) Obj->Clone();
							ModifCpnIRCurve->ParallelShift(EPSL);
						}
						else {
							ModifCpnIRCurve = ((ARM_ZeroCurve*)Obj)->GenerateShiftCurve(TermZC,&epsilon);
						}
						modelMkt->adopt(ModifCpnIRCurve);
					}else{
						ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_MktDataMng::GenerateShiftMng : bad market data for IR CPN sensi");  	
					}
				}

			}
			break;
			case ICM_INTEREST_CPN_CURVE :
			{
				ARM_ZeroCurve* ModifCpnIRCurve = NULL;
				for ( int i=0; i<labelKey.size(); i++){
					Obj = dynamic_cast<ARM_ZeroCurve*>(modelMkt->find(labelKey[i]));
					if (Obj) {
						if (Parallelshift)
						{
							ModifCpnIRCurve = (ARM_ZeroCurve*) Obj->Clone();
							ModifCpnIRCurve->ParallelShift(EPSL);
						}
						else {
							ModifCpnIRCurve = ((ARM_ZeroCurve*)Obj)->GenerateShiftCurve(TermZC,&epsilon);
						}
						modelMkt->adopt(ModifCpnIRCurve);
					}else{
						ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_MktDataMng::GenerateShiftMng : bad market data for IR CPN sensi");  	
					}
				}
			}
			break;
			case ICMSPRELSHIFT_TYPE :
			{

				for ( int i=0; i<labelKey.size(); i++){
					Obj = dynamic_cast<ICM_DefaultCurve*>(modelMkt->find(labelKey[i]));
					if(Obj) {
						// if (Parallelshift) 
						// 	ModifDefCurve = ((ICM_DefaultCurve*)Obj)->GenerateShiftCurve(EPSL_DEF,ICMSPRELSHIFT_TYPE);
						// else 
						ModifDefCurve = ((ICM_DefaultCurve*)Obj)->GenerateShiftCurve(Term,epsilon,ICMSPRELSHIFT_TYPE);
						modelMkt->adopt(ModifDefCurve);
					} else {
						ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_MktDataMng::GenerateShiftMng : bad market data for SPRELSHIFT sensi");  	
					}
				}		
			}
			break;
			case ICM_INDX_SPREAD_RESCALING :
			{
				for ( int i=0; i<labelKey.size(); i++){
					Obj = dynamic_cast<ICM_DefaultCurve*>(modelMkt->find(labelKey[i]));
					if(Obj) {
						// if (Parallelshift) 
						// 	ModifDefCurve = ((ICM_DefaultCurve*)Obj)->GenerateShiftCurve(EPSL_DEF,ICMSPREAD_TYPE);
						// else 
						ModifDefCurve = ((ICM_DefaultCurve*)Obj)->GenerateShiftCurve(Term,epsilon,ICMSPREAD_TYPE);
						modelMkt->adopt(ModifDefCurve);
					} else {
						ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_MktDataMng::GenerateShiftMng : bad market data for SPRELSHIFT sensi");  	
					}
				}
				// case ICM_INDX_SPREAD_RESCALING :
				Obj = dynamic_cast<ICM_Correlation*>(modelMkt->find(labelKey[0]));
				string label = ((ICM_Correlation*)Obj)->GetStructName();
				Correlation = (ICM_Correlation*)((ICM_Correlation*)Obj)->Clone();
				Correlation->ResetBetas();
				modelMkt->adopt(Correlation);

			}
			break;
			case ICMSPREAD_TYPE :
			{
				for ( int i=0; i<labelKey.size(); i++){
					Obj = dynamic_cast<ICM_DefaultCurve*>(modelMkt->find(labelKey[i]));
					if(Obj) {
						// if (Parallelshift) 
						// 	ModifDefCurve = ((ICM_DefaultCurve*)Obj)->GenerateShiftCurve(EPSL_DEF,ICMSPREAD_TYPE);
						// else 
						ModifDefCurve = ((ICM_DefaultCurve*)Obj)->GenerateShiftCurve(Term,epsilon,ICMSPREAD_TYPE);
						modelMkt->adopt(ModifDefCurve);
					} else {
						ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_MktDataMng::GenerateShiftMng : bad market data for SPRELSHIFT sensi");  	
					}
				}
			}break;
			case ICMBETA_TYPE :
			{
				Obj = dynamic_cast<ICM_Correlation*>(modelMkt->find(labelKey[0]));
				string label = ((ICM_Correlation*)Obj)->GetStructName();
				Correlation = ((ICM_Correlation*)Obj)->GenerateShiftBetas(label,typesensi,EPSL_BETA);
				modelMkt->adopt(Correlation);
			}
			break;

			case ICMCORRELATION_TYPE :
			{
				Obj = dynamic_cast<ICM_Correlation*>(modelMkt->find(labelKey[0]));
				string label = ((ICM_Correlation*)Obj)->GetStructName();
				Correlation = ((ICM_Correlation*)Obj)->GenerateShiftBetas(label,typesensi,EPSL_CORR);
				modelMkt->adopt(Correlation);
			}
			break;
			case ICMCORREL_STRIKE_DOWN_TYPE :
			{
				Obj = dynamic_cast<ICM_Correlation*>(modelMkt->find(labelKey[0]));
				string label = ((ICM_Correlation*)Obj)->GetStructName();
				Correlation = ((ICM_Correlation*)Obj)->GenerateShiftBetas(label,typesensi,-EPSL_CORR*10.);
				modelMkt->adopt(Correlation);
			}
			break;
			case ICMCORREL_STRIKE_UP_TYPE :
			{
				Obj = dynamic_cast<ICM_Correlation*>(modelMkt->find(labelKey[0]));
				string label = ((ICM_Correlation*)Obj)->GetStructName();
				Correlation = ((ICM_Correlation*)Obj)->GenerateShiftBetas(label,typesensi,EPSL_CORR*10.);
				modelMkt->adopt(Correlation);
			}
			break;
			default :
			{
				ICMTHROW(ERR_INVALID_MODEL,"Can't process SensiType "<<typesensi); 
			}
		}
	// ---------------------------------------------------------------------------
	// Bump Beta / Correlation
	// ---------------------------------------------------------------------------
 
	 

	return (modelMkt);


}

string GetCorrelationIndex(const string& structname, const string& index1, const string& index2){

	string name("");
	if (index1<index2) name=index1+"/"+index2 ;else name=index2+"/"+index1;
	if (name=="/") name=structname; 
	return name;
}