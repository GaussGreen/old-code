#include "firsttoinc.h"
#include "ICMKernel/pricer/ICM_Pricer_MktMng.h"
#include "ICMKernel/glob/icm_mktdatamng.h"
#include "ICMKernel/mod/icm_defcurvemodel.h"
#include "ICMKernel/inst/icm_ftd.h"
#include "ICMKernel/inst/icm_mez.h"
#include "ICMKernel/inst/icm_credit_index.h"
#include "ICMKernel/glob/icm_smile_correlation.h"
#include "ICMKernel/pricer/icm_pricer_adviser.h"
#include "ICMKernel/util/icm_utils.h"
#include "ICMKernel/crv/icm_distriblosscalculator.h" 
#include "ARMKernel\crv\volcube.h"

void ICM_Pricer_MktMng::Set(ARM_Security *sec, ARM_Object *mod,const ICM_Parameters& params,const ARM_Date&asof)
{
	if ((mod->GetName()!=ICM_MKTDATAMNG))			
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
            "ERROR : pricer not compatible with this  model ");

	ICM_Pricer::Set(sec, mod,params,&asof);
}


void ICM_Pricer_MktMng::View(char* id, FILE* ficOut)
{	
	FILE* fOut;
	char  fOutName[200];

	if ( ficOut == NULL )
	{
	ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
	
	   (void) unlink(fOutName);

       fOut = fopen(fOutName, "w"); 
    }
	else
	{
	fOut = ficOut;
	} 

	ICM_Pricer::View(id, fOut);

	
	if ( ficOut == NULL )
	{
		fclose(fOut);
	}
}


// *************************************************************
// Compute Sensitivity Method For Baskets
// *************************************************************
double ICM_Pricer_MktMng::ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
											 const std::string& plot, 
											 const std::string& labelKey, 
											 double epsvalue,
											 double epsvalueGamma)
{
	double sensitivity =0.;

	ICM_MktDataMng* MktModel = dynamic_cast<ICM_MktDataMng*> (GetMktDataMng());
	ARM_Date ExecutionDate = GetAsOfDate();
	
	std::auto_ptr<ICM_MktDataMng> MktModelShift ; 

	string formatedDate = ExecutionDate.toString('.'); // for format in ICM_Pricer
	ICM_Pricer_Advisor recupPricer;
	int i = 0;
	double result = 0.,initialprice = 0.,modifiedprice = 0.;

	//Perte subie effectivement par la tranche
	double Actual_Loss = 0. ;
    // try
    // {
	if (typesensi != ICMBETA_WITH_SPREAD_SHIFT)
	{
		initialprice = Price(qCMPPRICE);
	}
	switch (typesensi)
	{
		case ICMRECOVERY_TYPE :
		case ICMIRCURVE_TYPE :
		case ICMIRCURVE_WITHOUT_DEFCURVE_TYPE :
		case ICMSPREAD_TYPE :
		case ICMSPRELSHIFT_TYPE :
		case ICM_DTR_TYPE :
		case ICM_SAMEBETA_TYPE :
		case ICM_SAMECORRELATION_TYPE :
		case ICMCORRELATION_TYPE :
		case ICMBETA_TYPE :
		case ICM_INFLATION_CPN_CURVE :
		case ICM_IRCURVE_WITH_CPN :
		case ICM_INTEREST_CPN_CURVE :
		case ICM_INDX_SPREAD_RESCALING :
 		{				
			MktModelShift = std::auto_ptr<ICM_MktDataMng> (MktModel->GenerateShiftMng(typesensi,plot, labelKey,epsvalue));
		}
		break;
		case ICM_GREEK_VEGA_ATM_TYPE: 
		{
			MktModelShift = std::auto_ptr<ICM_MktDataMng>( (ICM_MktDataMng*)MktModel->Clone() );
			ICM_Credit_Index* index = NULL;
			ICM_Cds* cds = dynamic_cast<ICM_Cds*> (GetSecurity()); //For CM-CDS sensitivity
			if(cds) index = cds->GetFeeLeg()->GetCreditIndex();

			ICM_Leg* pLeg = dynamic_cast<ICM_Leg*>( GetSecurity()); //ex For Corridor leg 
			if (pLeg) index = pLeg->GetCreditIndex();
			// name of the Vol in MKtMng :
			string key = GetVolCurveName(index->GetLabels()[0], index->GetCurrencyUnit()->GetCcyName(), GetAsOfDate());
			ARM_VolCube* pIdxVolatility = dynamic_cast<ARM_VolCube*>(MktModelShift->find(key));
			if (!pIdxVolatility) 
				ICMTHROW(ERR_INVALID_MODEL,"ICM_GREEK_VEGA_ATM_TYPE: not a volCube "<<key); 
			
			std::auto_ptr<ARM_VolCube> cloned ( dynamic_cast<ARM_VolCube*> (pIdxVolatility->Clone()) ) ;
			ARM_VolCurve* atm = cloned->GetATMVol(); 
			if (!atm) 
				ICMTHROW(ERR_INVALID_MODEL,"ICM_GREEK_VEGA_ATM_TYPE: atm vol not defined in volCube "<<key); 
			
			int nthLine,nthCol ;
			SearchCoordForVolatility(*atm,plot,nthLine,nthCol) ;
			atm->BumpVolatility(epsvalue/100.,nthLine+1,nthCol+1); 
			MktModelShift->adopt(cloned.release()); 
		} 
		break; 
		case ICM_GREEK_IRVEGA_ATM_TYPE: 
		{
			MktModelShift = std::auto_ptr<ICM_MktDataMng>( (ICM_MktDataMng*)MktModel->Clone() );
			ARM_IRIndex * index = NULL;
			ICM_Cds* cds = dynamic_cast<ICM_Cds*> (GetSecurity()); //For CM-CDS sensitivity
			if(cds) index = cds->GetFeeLeg()->GetIRIndex();

			ICM_Leg* pLeg = dynamic_cast<ICM_Leg*>( GetSecurity()); //ex For Corridor leg 
			if (pLeg) index = pLeg->GetIRIndex();
			// name of the Vol in MKtMng :
			string key = GetVolCurveName(ARM_ParamView::GetMappingName(S_INDEX_TYPES, index->GetIndexType()),index->GetCurrencyUnit()->GetCcyName(), GetAsOfDate());
			ARM_VolCube* pIdxVolatility = dynamic_cast<ARM_VolCube*>(MktModelShift->find(key));
			if (!pIdxVolatility) 
				ICMTHROW(ERR_INVALID_MODEL,"ICM_GREEK_IRVEGA_ATM_TYPE: not a volCube "<<key); 
			
			std::auto_ptr<ARM_VolCube> cloned ( dynamic_cast<ARM_VolCube*> (pIdxVolatility->Clone()) ) ;
			ARM_VolCurve* atm = cloned->GetATMVol(); 
			if (!atm) 
				ICMTHROW(ERR_INVALID_MODEL,"ICM_GREEK_IRVEGA_ATM_TYPE: atm vol not defined in volCube "<<key); 
			
			int nthLine,nthCol ;
			SearchCoordForVolatility(*atm,plot,nthLine,nthCol) ;
			atm->BumpVolatility(epsvalue/100.,nthLine+1,nthCol+1); 
			MktModelShift->adopt(cloned.release()); 
		} 
		break; 
		case ICM_GREEK_GAMMA_TYPE :
		{
			// Compute Delta
			double Delta = ComputeSensitivity(ICMSPREAD_TYPE, plot, labelKey,  epsvalue);
			// Shift curve

			MktModelShift = std::auto_ptr<ICM_MktDataMng> (MktModel->GenerateShiftMng(ICMSPREAD_TYPE,"NONE", labelKey,epsvalueGamma)); // NONE for // Shift
			auto_ptr<ICM_Pricer_MktMng> newPricerShift(dynamic_cast<ICM_Pricer_MktMng*>(recupPricer.GeneratePricer(GetSecurity(),
																			(ARM_Model*)MktModelShift.get(), // CRAD !! mais obligatoire
																			GetName(),
																			CREDIT_DEFAULT_VALUE, // because nbpaths not used 
																			&GetParameters(),
																			(char*)formatedDate.c_str()) ));
			double DeltaShift = newPricerShift.get()->ComputeSensitivity(ICMSPREAD_TYPE, plot, labelKey,  epsvalue);
			return (DeltaShift - Delta);

		}
		break;
		default :
			{
				ICMTHROW(ERR_INVALID_MODEL,"Sensitivity not implemented: "<<typesensi);
			}
	}
	//
	//	At this stage the model is done. 
	//
	auto_ptr<ICM_Pricer> newPricerShift(recupPricer.GeneratePricer(GetSecurity(),
																			(ARM_Model*)MktModelShift.get(), // CRAD !! mais obligatoire
																			GetName(),
																			CREDIT_DEFAULT_VALUE, // because nbpaths not used 
																			&GetParameters(),
																			(char*)formatedDate.c_str()) );
	modifiedprice = newPricerShift.get()->ComputePrice(qCMPPRICE);

	result=modifiedprice - initialprice;
	return (result);
}

