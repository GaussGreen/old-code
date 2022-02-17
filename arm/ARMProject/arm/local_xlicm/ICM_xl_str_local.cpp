
#pragma warning(disable :4005 4786)

#include <libCCxll\CCxll.h>
#include <stdio.h>
#include <stdlib.h>


#include <ARM\libarm_local\ARM_local_glob.h>
#include <ARM\libarm_local\ARM_local_class.h>

#include <ARM\libicm_local\ICM_local_Str.h>

#include "ARM_local_interface.h"
#include "ARM_local_interglob.h"

#include "ARM_xl_trycatch_local.h"
#include <util\fromto.h>

// Faire correspondre cette valeur avec le Credit_Def_Value qui se trouve dans Credit\ICMKERNEL\glob\enums.h
double Credit_Def_Value = -999.;

#ifndef unix 
#define ITOA _itoa
#else
#define ITOA intToStr
#endif


__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Str_Cppi ( LPXLOPER XL_DefCurve,
														    LPXLOPER XL_Recov,
															LPXLOPER XL_Correl,
															LPXLOPER XL_DataCppi,
															LPXLOPER XL_RiskyFactor,
															LPXLOPER XL_DataPlacement,
															LPXLOPER XL_DataSoult,
															LPXLOPER XL_DataVasicek,
															LPXLOPER XL_DataModel,
															LPXLOPER XL_DataCoupon,
															LPXLOPER XL_DataWidener,
															LPXLOPER XL_DataSimul,
															LPXLOPER XL_DataBarriere,
															LPXLOPER XL_taux,
															LPXLOPER XL_Fonction)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;
	
	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	vector<double> CrbDefaut;
	long nbRowCrbDef;
	long nbColCrbDef;

	double Recov;

	vector<double> correl;
	long nbRowCorrel;
	long nbColCorrel;

	vector<double> DataCppi;

	vector<double> riskyFactor;
	long nbRowRiskyFactor;
	long nbColRiskyFactor;
	
	vector<double> DataPlacement;
	
	vector<double> DataSoult;
	
	vector<double> DataVasicek;
	
	vector<double> DataModel;
	
	vector<double> DataCoupon;

	vector<double> DataWidener;

	vector<double> DataSimul;

	vector<double> DataBarriere;

	vector<double> taux;
	long nbRowtaux;
	long nbColtaux;
	
	vector<double> Resultat;

	CCString Fonction;
	
	// error
	static int error;
	static char* reason = "";

	XL_readNumVectorAndSize(XL_DefCurve,nbRowCrbDef,nbColCrbDef,CrbDefaut," ARM_ERR: CrbDefaut : object matrix expected",C_result);
	int RowCrbDef = (int) nbRowCrbDef;
	int ColCrbDef = (int) nbColCrbDef;

	XL_readNumCell(XL_Recov,Recov,"ARM_ERR : Recov: object numeric expected",C_result);

	XL_readNumVectorAndSize(XL_Correl,nbRowCorrel,nbColCorrel,correl," ARM_ERR: correl : object matrix expected",C_result);
	int RowCorrel = (int) nbRowCorrel;
	int ColCorrel = (int) nbColCorrel;

	XL_readNumVector(XL_DataCppi,DataCppi,"ARM_ERR : Cppi : object vector of numeric expected",C_result);
	
	XL_readNumVectorAndSize(XL_RiskyFactor,nbRowRiskyFactor,nbColRiskyFactor,riskyFactor," ARM_ERR: RiskyFactor : object matrix expected",C_result);
	int RowRiskyFactor = (int) nbRowRiskyFactor;
	int ColRiskyFactor = (int) nbColRiskyFactor;

	XL_readNumVector(XL_DataPlacement,DataPlacement,"ARM_ERR : DataPlacement : object vector of numeric expected",C_result);
	XL_readNumVector(XL_DataSoult,DataSoult,"ARM_ERR : DataSoult : object vector of numeric expected",C_result);
	XL_readNumVector(XL_DataVasicek,DataVasicek,"ARM_ERR : DataVasicek : object vector of numeric expected",C_result);
	XL_readNumVector(XL_DataModel,DataModel,"ARM_ERR : DataModel : object vector of numeric expected",C_result);
	XL_readNumVector(XL_DataCoupon,DataCoupon,"ARM_ERR : DataCoupon : object vector of numeric expected",C_result);
	XL_readNumVector(XL_DataWidener,DataWidener,"ARM_ERR : DataWidener: object vector of numeric expected",C_result);
	XL_readNumVector(XL_DataSimul,DataSimul,"ARM_ERR : DataSimul : object vector of numeric expected",C_result);
	XL_readNumVector(XL_DataBarriere,DataBarriere,"ARM_ERR : DataBarriere : object vector of numeric expected",C_result);

	XL_readNumVectorAndSize(XL_taux,nbRowtaux,nbColtaux,taux,"ARM_ERR : taux : object numeric expected",C_result);
	int Rowtaux = (int) nbRowtaux;
	int Coltaux = (int) nbColtaux;

	XL_readStrCell(XL_Fonction,Fonction,"ARM_ERR : Fonction : String expected",C_result);
	
	long retCode;
	
	retCode = ICMLOCAL_STR_Cppi(CrbDefaut,
								RowCrbDef,
								ColCrbDef,
								Recov,
								correl,
								RowCorrel,
								ColCorrel,
								DataCppi,
								riskyFactor,
								RowRiskyFactor,
								ColRiskyFactor,
								DataPlacement,
								DataSoult,
								DataVasicek,
								DataModel,
								DataCoupon,
								DataWidener,
								DataSimul,
								DataBarriere,
								taux,
								Rowtaux,
								Coltaux,
								Fonction,
								Resultat,
								C_result);


if(retCode == ARM_OK)
{
	int nbcolumns = 20;
	int nbrows = (int)(Resultat.size()/nbcolumns);
	
	FreeCurCellErr ();
	XL_result.xltype = xltypeMulti| xlbitDLLFree;
	XL_result.val.array.columns = nbcolumns;
	XL_result.val.array.rows = nbrows; 
	XL_result.val.array.lparray = new XLOPER [nbrows * nbcolumns] ;

	for(int i = 0; i < nbrows; i++)
	{
		for (int j=0; j<nbcolumns;j++)
		{
			XL_result.val.array.lparray [i*nbcolumns + j].xltype = xltypeNum; 
			XL_result.val.array.lparray [i*nbcolumns + j].val.num= Resultat[i*nbcolumns + j	]; 
		}
	}
}
else
{
	ARM_ERR();
}

}

	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Str_Cppi" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

// -----------------------------------------------------------------------------------
// Generateur d'ALEA
// -----------------------------------------------------------------------------------
/*__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Str_GenAlea (LPXLOPER XL_Seed,
															 LPXLOPER XL_NbTir)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double dSeed;
	double SeedDefaultValue = Credit_Def_Value;
	double dNbTir;
	double NbTirDefaultValue = 1.;
		
	// error
	static int error;
	static char* reason = "";
	
	XL_readNumCellWD(XL_Seed,dSeed,Credit_Def_Value," ARM_ERR: Seed : numeric expected",C_result);
	XL_readNumCellWD(XL_NbTir,dNbTir,NbTirDefaultValue," ARM_ERR: NbTir : numeric expected",C_result);
	
	int iSeed = (int)dSeed;
	int iNbTir = (int)dNbTir;
	vector<double> Resultat;

	long retCode;
	retCode = ICMLOCAL_STR_GenAlea(iSeed,
								   iNbTir,
								   Resultat,
								   C_result);
		

if(retCode == ARM_OK)
{
	int nbrows = Resultat.size();
	int nbcolumns = 1;
	
	FreeCurCellErr ();
	XL_result.xltype = xltypeMulti| xlbitDLLFree;
	XL_result.val.array.columns = nbcolumns;
	XL_result.val.array.rows = nbrows; 
	XL_result.val.array.lparray = new XLOPER [nbrows * nbcolumns] ;

	for(int i = 0; i < nbrows; i++)
	{
		XL_result.val.array.lparray [i].xltype=xltypeNum; 
		XL_result.val.array.lparray [i].val.num=Resultat[i]; 
	}
}
else
{
	ARM_ERR();
}

}


	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Str_GenAlea" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
*/

// CDS
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
// Cette fonction permet le calcule de la FL, DL , NPV, Duration, Imply et Sensi par plot d'un CDS
// Utilise la calibration d'une courbe de défaut
/** JLA. Syn with ActiveX 
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Str_CDS (LPXLOPER XL_SpreadCDS,
														  LPXLOPER XL_Maturity,
														  LPXLOPER XL_PasPaiement,
														  LPXLOPER XL_Recov,
														  LPXLOPER XL_Nominal,
														  LPXLOPER XL_SpreadCourbe,
														  LPXLOPER XL_PlotCourbe,
														  LPXLOPER XL_Taux,
														  LPXLOPER XL_DatetoShift,
														  LPXLOPER XL_ValShift,
														  LPXLOPER XL_Stripping,
														  LPXLOPER XL_Fonction)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;
	LPXLOPER pxArray;

/// to make the infrastructure very robust we put a try catch!
ARM_XL_TRY_BLOCK_BEGIN
{
	ARM_NOCALCIFWIZ();

	// C variable
	double SpreadCds;
	double Maturity;
	double PasPaiement;
	double PasDefaultValue = 0.25;
	vector<double> VRecov;
	double Nominal;
	vector<double> VSpreadCourbe;
	double dSpreadCourbe;
	vector<double> PlotCourbe;
	vector<double> PlotCourbeDefaultValue;
	PlotCourbeDefaultValue.resize(5);
	PlotCourbeDefaultValue.push_back(1.);
	PlotCourbeDefaultValue.push_back(3.);
	PlotCourbeDefaultValue.push_back(5.);
	PlotCourbeDefaultValue.push_back(7.);
	PlotCourbeDefaultValue.push_back(10.);

	double Taux;
	double TauxDefaultValue = 0.;
	CCString C_Stripping;
	CCString C_StrippingDefaultValue = "FAST";
	
	CCString C_Fonction;

	// Date Shifte
	double dDateShifte = Credit_Def_Value;
	double DateShifte;

	// Val Shifte
	double ValShift = 0.01/100.;
	double ValShiftdefaultValue = 0.01/100.;
	// error
	static int error;
	static char* reason = "";
	
	XL_readNumCell(XL_SpreadCDS,SpreadCds," ARM_ERR: SpreadCds : numeric expected",C_result);
	XL_readNumCell(XL_Maturity,Maturity," ARM_ERR: Maturity : numeric expected",C_result);
	XL_readNumCellWD(XL_PasPaiement,PasPaiement,PasDefaultValue," ARM_ERR: PasPaiement : numeric expected",C_result);
	XL_readNumVector(XL_Recov,VRecov," ARM_ERR: Recov : vector of numeric expected",C_result);
	XL_readNumCell(XL_Nominal,Nominal," ARM_ERR: Nominal : numeric expected",C_result);
	XL_readNumVectorWD(XL_PlotCourbe,PlotCourbe,PlotCourbeDefaultValue," ARM_ERR: PlotCourbe : Vector of numeric expected",C_result);
	XL_readNumCellWD(XL_Taux,Taux,TauxDefaultValue," ARM_ERR: Taux : numeric expected",C_result);
	XL_readNumCellWD(XL_ValShift,ValShift,ValShiftdefaultValue," ARM_ERR: ValShift : numeric expected",C_result);
	XL_readStrCellWD(XL_Stripping,C_Stripping,C_StrippingDefaultValue," ARM_ERR: Stripping : String expected",C_result);
	XL_readStrCell(XL_Fonction,C_Fonction," ARM_ERR: Fonction : String expected",C_result);
	XL_readNumCellWD(XL_DatetoShift,DateShifte,dDateShifte," ARM_ERR: DateShifte : numeric expected",C_result);

	// On gere les cas ou l'utilisateur entre une valeur au clavier alors que l'on attend un vecteur
	if (XL_SpreadCourbe->xltype == 0x0001)
	{
		XL_readNumCell(XL_SpreadCourbe,dSpreadCourbe," ARM_ERR: SpreadCourbe : Vector of numeric expected",C_result);
		for (int i=0;i<PlotCourbe.size();i++)
			VSpreadCourbe.push_back(dSpreadCourbe);
	}
	else
		XL_readNumVectorWithHole(XL_SpreadCourbe,VSpreadCourbe," ARM_ERR: SpreadCourbe : Vector of numeric expected",C_result);

		

	long retCode;
	vector<double> Resultat;
	retCode = ICMLOCAL_STR_CDS(SpreadCds,
							   Maturity,
							   PasPaiement,
							   VRecov,
							   Nominal,
							   VSpreadCourbe,
							   PlotCourbe,
							   Taux,
							   DateShifte,
							   ValShift,
							   C_Stripping,
							   C_Fonction,
							   Resultat,
							   C_result);
	
if(retCode == ARM_OK)
	{
		int nbrows = 1;
		int nbcolumns = Resultat.size();
		
		FreeCurCellErr ();
		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		for (int j=0;j<nbcolumns;j++)
		{
			pxArray[XL_Coordonnate2Rank (0, j, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (0, j, nbcolumns)].val.num = Resultat[j]; 
		}
	}
	else
	{
		ARM_ERR();
	}
}


	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Str_CDS" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
***/ 

// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
// Fonction calculant le spread Forward d'une courbe de défaut calibrée
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Str_SpreadFwdCurves(LPXLOPER XL_Spread,
																	 LPXLOPER XL_Plot,
																	 LPXLOPER XL_Recovery,
																	 LPXLOPER XL_Taux,
																	 LPXLOPER XL_TimeInit,
																	 LPXLOPER XL_Duree,
																	 LPXLOPER XL_Stripping)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	vector<double> VSpread;
	double dSpread;
	vector<double> Plot;
	vector<double> PlotDefaultValue;
	PlotDefaultValue.push_back(1);
	PlotDefaultValue.push_back(3);
	PlotDefaultValue.push_back(5);
	PlotDefaultValue.push_back(7);
	PlotDefaultValue.push_back(10);
	double Taux;
	double TauxDefaultValue = 0.;
	double Recovery;
	double TimeInit;
	double Duree;
		
	CCString C_Stripping;
	CCString C_StrippingDefaultValue = "FAST";

	// error
	static int error;
	static char* reason = "";
	
	// On gere le cas ou l'utilisateur entre une valeur au clavier alors que l'on attend un vecteur
	XL_readNumVectorWD(XL_Plot,Plot,PlotDefaultValue," ARM_ERR: Plot : Vector of numeric expected",C_result);
	XL_readNumCell(XL_Recovery,Recovery," ARM_ERR: Recovery : numeric expected",C_result);
	XL_readNumCellWD(XL_Taux,Taux,TauxDefaultValue," ARM_ERR: Recovery : numeric expected",C_result);
	XL_readNumCell(XL_TimeInit,TimeInit," ARM_ERR: TimeInit : numeric expected",C_result);
	XL_readNumCell(XL_Duree,Duree," ARM_ERR: Duree : numeric expected",C_result);
	XL_readStrCellWD(XL_Stripping,C_Stripping,C_StrippingDefaultValue," ARM_ERR: Stripping : string expected",C_result);

	// Spread
	vector<double> VSpreadDefaultValue;
	VSpreadDefaultValue.resize(Plot.size());
	for (int i =0;i<Plot.size();i++)
		VSpreadDefaultValue.push_back(Credit_Def_Value);

	if (XL_Spread->xltype == 0x0001)
	{
		XL_readNumCell(XL_Spread,dSpread," ARM_ERR: Spread : Vector of numeric expected",C_result);
		for (int i=0;i<Plot.size();i++)
			VSpread.push_back(dSpread);
	}
	else
		XL_readNumVector(XL_Spread,VSpread," ARM_ERR: Spread : Vector of numeric expected",C_result);
	
	long retCode;
	retCode = ICMLOCAL_STR_SpreadForwardCurves(VSpread,
											   Plot,
											   Recovery,
											   Taux,
										       TimeInit,
										       Duree,
											   C_Stripping,
										       C_result);
	
	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble ();
	}
	else
	{
		ARM_ERR();
	}
	}

	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Price" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

// ------------------------------------------------------------------------------------

// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
// Renvoie le discount Factor d'une courbe ZeroCoupon Flat. Taux
/*__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Str_DPrice ( LPXLOPER XL_YF_Maturity,
															  LPXLOPER XL_Taux)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double YF_Maturity;
	double Taux;
	vector<double> DP;
	
	// error
	static int error;
	static char* reason = "";
	
	XL_readNumCell(XL_YF_Maturity,YF_Maturity," ARM_ERR: Maturity : numeric expected",C_result);
	XL_readNumCell(XL_Taux,Taux," ARM_ERR: Taux : numeric expected",C_result);
	
	long retCode;
	retCode = ICMLOCAL_STR_DPrice(YF_Maturity,
								  Taux,
								  C_result);
		

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble ();
	}
	else
	{
		ARM_ERR();
	}
	}


	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Price" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
*/
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
/*
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Str_ProbDefCurves (LPXLOPER XL_YF_Maturity,
																	LPXLOPER XL_IssuersRecovery,
																	LPXLOPER XL_Plot,
																	LPXLOPER XL_Spread,
																	LPXLOPER XL_Taux,
																	LPXLOPER XL_Stripping)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double YF_Maturity;
	double IssuersRecovery;
	vector<double> VSpread;
	double dSpread;
	vector<double> Plot;
	vector<double> PlotDefaultValue;
	PlotDefaultValue.push_back(1.);
	PlotDefaultValue.push_back(3.);
	PlotDefaultValue.push_back(5.);
	PlotDefaultValue.push_back(7.);
	PlotDefaultValue.push_back(10.);
	double Taux;
	double TauxDefaultValue = 0.;

	CCString C_Stripping;
	CCString C_StrippingDefaultValue = "FAST";
	// error
	static int error;
	static char* reason = "";
	
	XL_readNumVector(XL_Plot, Plot," ARM_ERR: Plot : Vector of numerics expected",C_result);
	XL_readNumCell(XL_YF_Maturity,YF_Maturity," ARM_ERR: Maturity : numeric expected",C_result);
	XL_readNumCell(XL_IssuersRecovery,IssuersRecovery," ARM_ERR: IssuersRecovery : numeric expected",C_result);
	XL_readNumCellWD(XL_Taux,Taux,TauxDefaultValue," ARM_ERR: Taux : numeric expected",C_result);
	XL_readStrCellWD(XL_Stripping,C_Stripping,C_StrippingDefaultValue," ARM_ERR: Stripping : string expected",C_result);
	// Spread
	if (XL_Spread->xltype == 0x0001)
	{
		XL_readNumCell(XL_Spread,dSpread," ARM_ERR: Spread : Vector of numeric expected",C_result);	
		for (int i=0;i<Plot.size();i++)
			VSpread.push_back(dSpread);
	}
	else
		XL_readNumVector(XL_Spread, VSpread," ARM_ERR: Spread : Vector of numerics expected",C_result);

	// ------------------------------------------------------
	long retCode;

	retCode = ICMLOCAL_STR_PDefaultCurves(YF_Maturity,
										  IssuersRecovery,
										  VSpread,
										  Plot,
										  Taux,
										  C_Stripping,
										  C_result);
	
	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble ();
	}
	else
	{
		ARM_ERR();
	}
	}

	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Price" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
*/
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
/*
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Str_ProbDefautTrancheZC(LPXLOPER XL_YF_Maturity,
																	     LPXLOPER XL_Attachement,
																	     LPXLOPER XL_NbIssuers,
																	     LPXLOPER XL_Recovery,
																	     LPXLOPER XL_MatSpreads,
																	     LPXLOPER XL_Plot,
																	     LPXLOPER XL_BetaPrice,
																	     LPXLOPER XL_Taux,
																		 LPXLOPER XL_Stripping)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double YF_Maturity;
	vector<double> Attachement;
	double dNb_Issuers;
	vector<double> Recovery;
	double dRecovery;
	vector<double> MatSpread;
	double dMatSpread;
	vector<double> Plot;
	vector<double> VBeta;
	double dBeta;
	double Taux;
	double TauxDefaultValue = 0.;

	CCString C_Stripping;
	CCString C_StrippingDefaultValue = "FAST";
		
	// error
	static int error;
	static char* reason = "";
	
	XL_readNumVector(XL_Plot, Plot," ARM_ERR: Plot : Vector of numerics expected",C_result);
	XL_readNumCell(XL_YF_Maturity,YF_Maturity," ARM_ERR: Maturity : numeric expected",C_result);
	XL_readNumVector(XL_Attachement,Attachement," ARM_ERR: Attachement : vector of numerics expected",C_result);
	XL_readNumCell(XL_NbIssuers,dNb_Issuers,"ARM_ERR: NbIssuers : numeric expected",C_result);
	XL_readNumCellWD(XL_Taux,Taux,TauxDefaultValue," ARM_ERR: Taux : numeric expected",C_result);
	XL_readStrCellWD(XL_Stripping,C_Stripping,C_StrippingDefaultValue," ARM_ERR: Stripping : numeric expected",C_result);

	// Recov
	if (XL_Recovery->xltype == 0x0001)
	{
		XL_readNumCell(XL_Recovery,dRecovery,"ARM_ERR: Recovery : vector of numeric expected",C_result);
		Recovery.push_back(dRecovery);
		Recovery.push_back(dRecovery);
	}
	else
		XL_readNumVector(XL_Recovery,Recovery," ARM_ERR: Recovery : vector of numerics expected",C_result);
	
	// Spreads
	if (XL_MatSpreads->xltype == 0x0001)
	{
		XL_readNumCell(XL_MatSpreads,dMatSpread," ARM_ERR: MatSpread : Matrix of numeric expected",C_result);
		for (int i = 0;i<(int)dNb_Issuers*Plot.size();i++)
			MatSpread.push_back(dMatSpread);
	}
	else
		XL_readNumVectorWithHole(XL_MatSpreads, MatSpread," ARM_ERR: MatSpread : Matrix of numerics expected",C_result);
	
	// Beta
	if (XL_BetaPrice->xltype == 0x0001)
	{
		XL_readNumCell(XL_BetaPrice,dBeta,"ARM_ERR: Beta : vector of numeric expected",C_result);
		for (int i = 0 ;i<(int)dNb_Issuers;i++)
			VBeta.push_back(dBeta);
	}
	else
		XL_readNumVector(XL_BetaPrice, VBeta," ARM_ERR: Beta : vector of numerics expected",C_result);


	long retCode;

	retCode = ICMLOCAL_STR_ProbDefTranche(YF_Maturity,
										  Attachement,
										  dNb_Issuers,
										  Recovery,
										  MatSpread,
										  Plot,
										  VBeta,
										  Taux,
										  C_Stripping,
										  C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble ();
	}
	else
	{
		ARM_ERR();
	}
}


	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Price" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
*/
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
/*
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Str_DensityCurve(LPXLOPER XL_YF_Maturity,									
																  LPXLOPER XL_NbIssuers,
																  LPXLOPER XL_IssuersRecovery,
																  LPXLOPER XL_IssuersSpread,
																  LPXLOPER XL_Plot,
																  LPXLOPER XL_IssuersBetas,
																  LPXLOPER XL_Taux)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;
	LPXLOPER pxArray;

	/// to make the infrastructure very robust we put a try catch!
ARM_XL_TRY_BLOCK_BEGIN
{
	ARM_NOCALCIFWIZ();

	// C variable
	double YF_Maturity;
	double dNb_Issuers;
	double IssuersRecovery;
	double Taux;
	double TauxDefaultValue = 0.;
	vector<double> VIssuersSpread;
	double dIssuersSpread;
	vector<double> Plot;
	vector<double> VBeta;
	double dBeta;
	
	vector<double> Density;
	
	// error
	static int error;
	static char* reason = "";
	
	XL_readNumCell(XL_YF_Maturity,YF_Maturity," ARM_ERR: Maturity : numeric expected",C_result);
	XL_readNumCell(XL_NbIssuers,dNb_Issuers,"ARM_ERR: NbIssuers : numeric expected",C_result);
	XL_readNumCell(XL_IssuersRecovery,IssuersRecovery," ARM_ERR: IssuersRecovery : numeric expected",C_result);
	XL_readNumCellWD(XL_Taux,Taux,TauxDefaultValue," ARM_ERR: Taux : numeric expected",C_result);
	XL_readNumVector(XL_Plot, Plot," ARM_ERR: Plot : vector of numerics expected",C_result);
	// Spread
	if (XL_IssuersSpread->xltype == 0x0001)
	{
		XL_readNumCell(XL_IssuersSpread,dIssuersSpread,"ARM_ERR: IssuersSpread : vector of numeric expected",C_result);
		VIssuersSpread.push_back(dIssuersSpread);
	}
	else
		XL_readNumVector(XL_IssuersSpread, VIssuersSpread," ARM_ERR: IssuersSpread : vector of numerics expected",C_result);
	
	// Beta
	if (XL_IssuersBetas->xltype == 0x0001)
	{
		XL_readNumCell(XL_IssuersBetas,dBeta,"ARM_ERR: Beta : vector of numeric expected",C_result);
		VBeta.push_back(dBeta);
	}
	else
		XL_readNumVector(XL_IssuersBetas, VBeta," ARM_ERR: Beta : vector of numerics expected",C_result);

	int NbIssuers = (int) dNb_Issuers;
	char** Label = new char*[NbIssuers];
	for (int i=0;i<NbIssuers;i++)
	{
		Label[i] = new char[60];
		sprintf(Label[i],"ISSUER_%i",i);
	}

	long retCode;
	retCode = ICMLOCAL_STR_DensityCurves (YF_Maturity,
										  dNb_Issuers,
										  IssuersRecovery,
										  VIssuersSpread,
										  Plot,
										  Label,
										  VBeta,
										  Taux,
										  Density,
										  C_result);
	
	if (Label)
	{
		for (int i=0;i<NbIssuers;i++)
		{
			if (Label[i])
			{
				delete[] Label[i];
				Label[i] = NULL;
			}
		}
		delete[] Label;
		Label = NULL;
	}
		

	if(retCode == ARM_OK)
	{
		int nbrows = Density.size ();
		int nbcolumns = 1;
		
		FreeCurCellErr ();
		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		for(int i = 0; i < nbrows; i++)
		{
			pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.num = Density[i]; 
		}
	}
	else
	{
		ARM_ERR();
	}
}



	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Price" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
*/
// -------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------

_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Str_CDO2SMILE(LPXLOPER XL_DataMez,
															  LPXLOPER XL_NotionalsMez,
															  LPXLOPER XL_StrikeCDO,
															  LPXLOPER XL_DataDefaultCurve,
															  LPXLOPER XL_IssuersLabel,
															  LPXLOPER XL_SpreadTraxx,
															  LPXLOPER XL_MatriceOVerlap,
															  LPXLOPER XL_DataTraxx,
															  LPXLOPER XL_Recov,
															  LPXLOPER XL_ForcedCorrel,
															  LPXLOPER XL_CorrelParStrike,
															  LPXLOPER XL_LabelToShift,
															  LPXLOPER XL_DateToShift,
															  LPXLOPER XL_ValShift,
															  LPXLOPER XL_Numerical_Method,
															  LPXLOPER XL_Fonction)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	vector<double> DataMez;
	vector<double> DataMezDefaultValue;
	DataMezDefaultValue.resize(8);
	DataMezDefaultValue[0] = 5.;    // maturity
	DataMezDefaultValue[1] = 0.;    // Strike Down
	DataMezDefaultValue[2] = 0.3;   // Strike Up 
	DataMezDefaultValue[3] = 0.25;  // Pas
	DataMezDefaultValue[4] = 0.01;    // Fixed Rate
	DataMezDefaultValue[5] = 125;   // nb nom
	DataMezDefaultValue[6] = 10;    // Nb tranche fille
	DataMezDefaultValue[7] = 0.;    // Taux
	
	vector<double> NotionalsTranche;
	vector<double> NotionalsTrancheDefaultValue;
	NotionalsTrancheDefaultValue.push_back(10000000);
	double dNotionalsTranche;
	double dNotionalsTrancheDefaultValue = 10000000;
	
	vector<CCString> StrikeCDO;
	long nbColStrikeCdo = 0;
	long nbRowStrikeCdo = 0;
	
	vector<CCString> DataDefaultCurve;
	long nbColDefaultCurve = 0;
	long nbRowDefaultCurve = 0;
		
	vector<CCString> IssuersLabel;

	vector<double> SpreadTraxx;
	vector<double> SpreadTraxxDefaultValue;
	SpreadTraxxDefaultValue.push_back(Credit_Def_Value);

	vector<CCString> MatriceOverlap;
	long nbColOverlap = 0;
	long nbRowOverlap = 0;

	vector<double> DataTraxx;
	
	vector<CCString> Recov;
	long nbColonneRecov = 0;
	long nbLigneRecov = 0;

	vector<CCString> ForcedCorrel;
	long nbColForcedCorrel = 0;
	long nbRowForcedCorrel = 0;
	
	vector<CCString> CorrelParStrike;
	long nbColCorrelParStrike = 0;
	long nbRowCorrelParStrike = 0;
	
	CCString LabelToShift;
	CCString LabelToShiftDefaultValue = "NONE";
	
	double DateToShift;
	double DateToShiftDefaultValue = Credit_Def_Value;
	
	double ValShift = 0.;
	double ValShiftDefaultValue = 0.01/100.;

	vector<CCString> NumMethod;
	vector<CCString> NumMethodDefaultValue;
	NumMethodDefaultValue.resize(3);
	NumMethodDefaultValue[0] = "20";
	NumMethodDefaultValue[1] = "H";
	NumMethodDefaultValue[2] = "FAST";
	
	CCString C_Function;
 
	// error
	static int error;
	static char* reason = "";

	// DataMez
	XL_readNumVectorWD(XL_DataMez,DataMez,DataMezDefaultValue," ARM_ERR: DataMez : Vector of numeric expected",C_result);

	// Notionals
	if (XL_NotionalsMez->xltype == 0x0001)
	{
		XL_readNumCellWD(XL_NotionalsMez,dNotionalsTranche,dNotionalsTrancheDefaultValue,"ARM_ERR: Notional : vector of numeric expected",C_result);
		NotionalsTranche.push_back(dNotionalsTranche);
	}
	else
		XL_readNumVectorWD(XL_NotionalsMez,NotionalsTranche,NotionalsTrancheDefaultValue," ARM_ERR: NotionalsTranche : Vector of numeric expected",C_result);

	// StrikeCDO
	XL_readStrVectorAndSize(XL_StrikeCDO,nbRowStrikeCdo,nbColStrikeCdo,StrikeCDO," ARM_ERR: StrikeCDO : object matrix expected",DOUBLE_TYPE,C_result);
	int nbColStrikCdo = (int) nbColStrikeCdo;
	int nbRowStrikCdo = (int) nbRowStrikeCdo;
	
	// DataDefaultCurve
	XL_readStrVectorAndSize(XL_DataDefaultCurve,nbRowDefaultCurve,nbColDefaultCurve,DataDefaultCurve," ARM_ERR: DataDefaultCurve: object matrix expected",DOUBLE_TYPE,C_result);
	int nbColDefCurve = (int) nbColDefaultCurve;
	int nbRowDefCurve = (int) nbRowDefaultCurve;
	
	// IssuersLabel
	int nbLabel = (int)DataMez[5];
	vector<CCString> LabelDefaultValue;
	LabelDefaultValue.resize(nbLabel);
	char* Lbl= new char[5];
	for (int i = 0;i<nbLabel;i++)
	{
		ITOA(i,Lbl,10);
		LabelDefaultValue[i] = Lbl;
	}
	
	if (Lbl)
		delete[] Lbl;
	Lbl = NULL;

	XL_readStrVectorWD(XL_IssuersLabel,IssuersLabel,LabelDefaultValue," ARM_ERR: LabelToShift : String expected",DOUBLE_TYPE,C_result);
	
	// SpreadTraxx
	XL_readNumVectorWD(XL_SpreadTraxx,SpreadTraxx,SpreadTraxxDefaultValue," ARM_ERR: SpreadTraxx : Vector of numeric expected",C_result);


	// Overlap
	XL_readStrVectorAndSize(XL_MatriceOVerlap,nbRowOverlap,nbColOverlap,MatriceOverlap," ARM_ERR: MatriceOverlap: object matrix expected",DOUBLE_TYPE,C_result);
	int nbColOver = (int) nbColOverlap;
	int nbRowOver = (int) nbRowOverlap;


	
	// DataTraxx
	XL_readNumVector(XL_DataTraxx,DataTraxx," ARM_ERR: DataTraxx : Vector of numeric expected",C_result);

	// Recov
	XL_readStrVectorAndSize(XL_Recov,nbLigneRecov,nbColonneRecov,Recov," ARM_ERR: Recov : object matrix expected",DOUBLE_TYPE,C_result);
	int nbColRecov = (int) nbColonneRecov;
	int nbRowRecov = (int) nbLigneRecov;

	// Forced Correl
	XL_readStrVectorAndSize(XL_ForcedCorrel,nbRowForcedCorrel,nbColForcedCorrel,ForcedCorrel," ARM_ERR: ForcedCorrel : object matrix expected",DOUBLE_TYPE,C_result);
	int nbColCorrelForced= (int) nbColForcedCorrel;
	int nbRowCorrelForced = (int) nbRowForcedCorrel;

	// Correl par strike
	XL_readStrVectorAndSize(XL_CorrelParStrike,nbRowCorrelParStrike,nbColCorrelParStrike,CorrelParStrike," ARM_ERR: CorrelParStrike : object matrix expected",DOUBLE_TYPE,C_result);
	int nbColCorrStrike = (int) nbColCorrelParStrike;
	int nbRowCorrStrike = (int) nbRowCorrelParStrike;

	if(CorrelParStrike.size() == 0)
	{
		C_result.setMsg ("ARM_ERR: check your CorrelParStrike Matrix");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}


	// Label To Shift
	XL_readStrCellWD(XL_LabelToShift,LabelToShift,LabelToShiftDefaultValue," ARM_ERR: LabelToShift : string expected",C_result);

	// Date To Shift
	XL_readNumCellWD(XL_DateToShift,DateToShift,DateToShiftDefaultValue," ARM_ERR: DateToShift : Numeric expected",C_result);

	// ValShift
	XL_readNumCellWD(XL_ValShift,ValShift,ValShiftDefaultValue," ARM_ERR: ValShift : numeric expected", C_result);

	// Data Nnum methode
	XL_readStrVectorWD(XL_Numerical_Method,NumMethod,NumMethodDefaultValue,"ARM_ERR : Numerical Method : Vecteur expcected",DOUBLE_TYPE, C_result);

	// Fonction
	XL_readStrCell(XL_Fonction,C_Function,"ARM_ERR: Fonction : string expected", C_result);




// -------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------
	long retCode;

	retCode = ICMLOCAL_STR_CdO2Smile(DataMez,
								    NotionalsTranche,
								    StrikeCDO,
								    nbColStrikCdo,
								    DataDefaultCurve,
								    nbColDefCurve ,
								    nbRowDefCurve ,								
								    IssuersLabel,
								    SpreadTraxx,
								    MatriceOverlap,
									nbColOver,
									nbRowOver,
								    DataTraxx,
								    Recov,
								    nbColRecov,
								    nbRowRecov,
								    ForcedCorrel,
								    nbColCorrelForced,
								    nbRowCorrelForced,
								    CorrelParStrike,
								    nbRowCorrStrike,
								    nbColCorrStrike,
								    LabelToShift,
								    DateToShift,
								    ValShift,
								    NumMethod,
								    C_Function,
								    C_result);

	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble ();
	}
	else
	{
		ARM_ERR();
	}

}


	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Str_CDO2SMILE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

// --------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Str_Basket(LPXLOPER XL_YF_Maturity,
															LPXLOPER XL_IndiceDefault,
															LPXLOPER XL_NbIssuers,
															LPXLOPER XL_IssuersNotionals,
															LPXLOPER XL_Calib_Recovery,
															LPXLOPER XL_Loss_Recovery,
															LPXLOPER XL_YFPlot,
															LPXLOPER XL_MatIssuersSpread,
															LPXLOPER XL_IssuersLabel,
															LPXLOPER XL_Pas,
															LPXLOPER XL_FixedRate,
															LPXLOPER XL_IssuerBeta,
															LPXLOPER XL_Taux,
															LPXLOPER XL_LabelToShift,
															LPXLOPER XL_DateToShift,
															LPXLOPER XL_ValShift,
															LPXLOPER XL_Stripping,
															LPXLOPER XL_Function)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double YF_Maturity;
	
	vector<double> IndiceDefault;
	double dIndiceDefault;

	double dNb_Issuers;

	double IssuersNotionals;
	
	vector<double> CalibRecovery;
	double dCalibRecovery;

	vector<double> LossRecovery;
	double dLossRecovery;
	
	vector<double> YFPlot;
	double dYFPlot;
	
	vector<double> MatIssuersSpread;
	double dIssuersSpread;
	
	vector<CCString> IssuersLabel;

	double Pas;
	double PasDefaultValue = 0.25;
	
	double FixedRate;
		
	vector<double> IssuersBeta;
	double dIssuersBeta;

	double Taux;
	double TauxDefaultValue = 0.;

	CCString LabelToShift;
	CCString LabelToShiftDefaultValue = "NONE";

	double ValShift = 0.;
	double ValShiftdefaultValue = 0.01/100.;

	CCString C_Stripping;
	CCString C_StrippingDefaultValue = "FAST";
	CCString C_Function;

	// error
	static int error;
	static char* reason = "";
	
	XL_readNumCell(XL_YF_Maturity,YF_Maturity," ARM_ERR: Maturity : numeric expected",C_result);
	XL_readNumCell(XL_IssuersNotionals,IssuersNotionals," ARM_ERR: IssuersNotionals : numeric expected",C_result);
	XL_readNumCell(XL_FixedRate,FixedRate,"ARM_ERR: FixedRate : numeric expected",C_result);
	XL_readNumCell(XL_NbIssuers,dNb_Issuers,"ARM_ERR: NbIssuers : numeric expected",C_result);
	XL_readNumCellWD(XL_Pas,Pas,PasDefaultValue," ARM_ERR: Pas : numeric expected",C_result);
	XL_readNumCellWD(XL_Taux,Taux,TauxDefaultValue," ARM_ERR: Taux : numeric expected",C_result);
	XL_readNumCellWD(XL_ValShift,ValShift,ValShiftdefaultValue," ARM_ERR: ValShift : numeric expected",C_result);
	XL_readStrCellWD(XL_Stripping,C_Stripping,C_StrippingDefaultValue," ARM_ERR: Stripping : String expected",C_result);
	XL_readStrCell(XL_Function,C_Function," ARM_ERR: C_Function : String expected",C_result);
	XL_readStrCellWD(XL_LabelToShift,LabelToShift,LabelToShiftDefaultValue," ARM_ERR: LabelToShift : String expected",C_result);

	// Defaults
	if (XL_IndiceDefault->xltype == 0x0001)
	{
		XL_readNumCell(XL_IndiceDefault,dIndiceDefault,"ARM_ERR : IndiceDefault : Vector of numeric expected",C_result);
		IndiceDefault.push_back(dIndiceDefault);
	}
	else
		XL_readNumVector(XL_IndiceDefault,IndiceDefault," ARM_ERR: IndiceDefault : Vector of numeric expected",C_result);

	// YF plot
	if (XL_YFPlot->xltype == 0x0001)
	{
		XL_readNumCell(XL_YFPlot,dYFPlot,"ARM_ERR : Plot : Vector of numeric expected",C_result);
		YFPlot.push_back(dYFPlot);
	}
	else
		XL_readNumVector(XL_YFPlot,YFPlot," ARM_ERR: YFPlot : Vector numeric expected",C_result);

	// Recov
	if (XL_Calib_Recovery->xltype == 0x0001)
	{
		XL_readNumCell(XL_Calib_Recovery,dCalibRecovery,"ARM_ERR: CalibRecovery : vector of numeric expected",C_result);
		CalibRecovery.push_back(dCalibRecovery);
		CalibRecovery.push_back(dCalibRecovery);
	}
	else
		XL_readNumVector(XL_Calib_Recovery,CalibRecovery," ARM_ERR: CalibRecovery : vector of numerics expected",C_result);

	vector<double> LossRecoveryDefaultValue;
	LossRecoveryDefaultValue.resize(CalibRecovery.size());
	LossRecoveryDefaultValue = CalibRecovery;
	// Recov
	if (XL_Loss_Recovery->xltype == 0x0001)
	{
		XL_readNumCell(XL_Loss_Recovery,dLossRecovery,"ARM_ERR: LossRecovery : vector of numeric expected",C_result);
		LossRecovery.push_back(dLossRecovery);
		LossRecovery.push_back(dLossRecovery);
	}
	else
		XL_readNumVectorWD(XL_Loss_Recovery,LossRecovery,LossRecoveryDefaultValue," ARM_ERR: LossRecovery : vector of numerics expected",C_result);

	// BetaIssuers
	if (XL_IssuerBeta->xltype == 0x0001)
	{
		XL_readNumCell(XL_IssuerBeta,dIssuersBeta,"ARM_ERR: Beta : vector of numeric expected",C_result);
		for (int i = 0 ;i<(int)dNb_Issuers;i++)
			IssuersBeta.push_back(dIssuersBeta);
	}
	else
		XL_readNumVector(XL_IssuerBeta, IssuersBeta," ARM_ERR: Beta : vector of numerics expected",C_result);
	
	
	int taille = (int)(YF_Maturity/Pas) + 5;
	
	// Spread	
	int size = YFPlot.size()*(int)dNb_Issuers;
	vector<double> MatSpreadDefaultValue;
	for (int i = 0;i<size;i++) 
		MatSpreadDefaultValue.push_back(Credit_Def_Value);

	if (XL_MatIssuersSpread->xltype == 0x0001)
	{
		XL_readNumCellWD(XL_MatIssuersSpread,dIssuersSpread,Credit_Def_Value,"ARM_ERR: MatIssuersSpread : Matrix of numerics expected",C_result);
		for (i = 0;i<size;i++) 
			MatIssuersSpread.push_back(dIssuersSpread);
	}
	else
		XL_readNumVectorWithHoleWD(XL_MatIssuersSpread, MatIssuersSpread,MatSpreadDefaultValue," ARM_ERR: MatIssuersSpread : Matrix of numerics expected",C_result);

	// Labels
	int nbLabel = (int)dNb_Issuers;
	vector<CCString> LabelDefaultValue;
	LabelDefaultValue.resize(nbLabel);
	char* Lbl= new char[5];
	for (i = 0;i<nbLabel;i++)
	{
		ITOA(i,Lbl,10);
		LabelDefaultValue[i] = Lbl;
	}
	
	if (Lbl)
		delete[] Lbl;
	Lbl = NULL;
	
	XL_readStrVectorWD(XL_IssuersLabel,IssuersLabel,LabelDefaultValue," ARM_ERR: Label : vector of string expected",DOUBLE_TYPE,C_result);


	// Date Shifte
	double dDateShifte = Credit_Def_Value;
	double DateShifte;
	
	XL_readNumCellWD(XL_DateToShift,DateShifte,dDateShifte," ARM_ERR: DateShifte : numeric expected",C_result);
// -------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------
	long retCode;
	double Freq = 0.;
	if ( (fabs(Pas)>1E-6) && (fabs(Pas-YF_Maturity)>1E-6))
		Freq = 1./Pas;
	
	retCode = ICMLOCAL_STR_Basket(YF_Maturity,
							      IndiceDefault,
								  dNb_Issuers,
								  IssuersNotionals,
								  CalibRecovery,
								  LossRecovery,
								  YFPlot,
								  MatIssuersSpread,								
								  IssuersLabel,
								  Freq,
								  FixedRate,
								  IssuersBeta,
								  Taux,
								  LabelToShift,
								  DateShifte,
								  ValShift,
								  C_Stripping,
								  C_Function,
								  C_result);
	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble ();
	}
	else
	{
		ARM_ERR();
	}
}


	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Price" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


// -----------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------

_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Str_CDO(LPXLOPER XL_MezTyp,
														LPXLOPER XL_DataMez,
														LPXLOPER XL_Notionals,
														LPXLOPER XL_Calib_Recovery,
														LPXLOPER XL_Loss_Recovery,
														LPXLOPER XL_YFPlot,
														LPXLOPER XL_MatIssuersSpread,
														LPXLOPER XL_IssuersLabel,
														LPXLOPER XL_IssuerBeta,
														LPXLOPER XL_Taux,
														LPXLOPER XL_LabelToShift,
														LPXLOPER XL_DateToShift,
														LPXLOPER XL_ValShift,
														LPXLOPER XL_SpreadCrbFwd,
														LPXLOPER XL_VolConvexity,
														LPXLOPER XL_Numerical_Method,
														LPXLOPER XL_Fonction)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_MezType;

	vector<double> DataMez;

	vector<double> Notionals;
	vector<double> NotionalsDefaultValue;
	NotionalsDefaultValue.push_back(10000000);
	double dNotional;
	double NotionalDefaultValue = 10000000;
	
	vector<double> CalibRecovery;
	double dCalibRecovery;

	vector<double> LossRecovery;
	double dLossRecovery;
	
	vector<double> YFPlot;
	double dYFPlot;
	
	vector<double> MatIssuersSpread;
	double dIssuersSpread;
	
	vector<CCString> IssuersLabel;
	
	vector<double> IssuersBeta;
	double dIssuersBeta;

	double Taux;
	double TauxDefaultValue = 0.;

	CCString LabelToShift;
	CCString LabelToShiftDefaultValue = "NONE";

	double ValShift = 0.01/100.;
	double ValShiftDefaultValue = 0.01/100.;

	vector<double> SpreadCrbFwd;
	vector<double> SpreadCrbFwdDefaultValue;
	SpreadCrbFwdDefaultValue.push_back(Credit_Def_Value);
	double dSpreadCrbFwd;

	double VolConvexity;
	double VolConvexityDefaultValue = 0.;

	vector<CCString> NumMethod;
	vector<CCString> NumMethodDefaultValue;
	NumMethodDefaultValue.resize(3);
	NumMethodDefaultValue[0] = "20";
	NumMethodDefaultValue[1] = "H";
	NumMethodDefaultValue[2] = "FAST";
	
	CCString C_Function;

	// error
	static int error;
	static char* reason = "";
	
	XL_readStrCellWD(XL_MezTyp,C_MezType,"RUNNING"," ARM_ERR: MezType : String expected",C_result);
	XL_readStrCell(XL_Fonction,C_Function," ARM_ERR: C_Function : String expected",C_result);
	XL_readStrCellWD(XL_LabelToShift,LabelToShift,LabelToShiftDefaultValue," ARM_ERR: LabelToShift : String expected",C_result);

	XL_readNumCellWD(XL_Taux,Taux,TauxDefaultValue," ARM_ERR: Taux : numeric expected",C_result);
	XL_readNumCellWD(XL_VolConvexity,VolConvexity,VolConvexityDefaultValue," ARM_ERR: VolConvexity : numeric expected",C_result);
	XL_readNumCellWD(XL_ValShift,ValShift,ValShiftDefaultValue," ARM_ERR: NValShift : numeric expected",C_result);
	
	// Data Mez
	XL_readNumVector(XL_DataMez,DataMez," ARM_ERR: DataMez : Vector of numeric expected",C_result);
	
	// Data Nnum methode
	XL_readStrVectorWD(XL_Numerical_Method,NumMethod,NumMethodDefaultValue,"ARM_ERR : Numerical Method : Vecteur expcected",DOUBLE_TYPE, C_result);

	// Notionals
	if (XL_Notionals->xltype == 0x0001)
	{
		XL_readNumCellWD(XL_Notionals,dNotional,NotionalDefaultValue,"ARM_ERR: Notional : vector of numeric expected",C_result);
		Notionals.push_back(dNotional);
	}
	else
		XL_readNumVectorWD(XL_Notionals,Notionals,NotionalsDefaultValue," ARM_ERR: Notional : vector of numerics expected",C_result);


	// YF plot
	if (XL_YFPlot->xltype == 0x0001)
	{
		XL_readNumCell(XL_YFPlot,dYFPlot,"ARM_ERR : Plot : Vector of numeric expected",C_result);
		YFPlot.push_back(dYFPlot);
	}
	else
		XL_readNumVector(XL_YFPlot,YFPlot," ARM_ERR: YFPlot : Vector numeric expected",C_result);

	// Recov
	if (XL_Calib_Recovery->xltype == 0x0001)
	{
		XL_readNumCell(XL_Calib_Recovery,dCalibRecovery,"ARM_ERR: CalibRecovery : vector of numeric expected",C_result);
		CalibRecovery.push_back(dCalibRecovery);
		CalibRecovery.push_back(dCalibRecovery);
	}
	else
		XL_readNumVector(XL_Calib_Recovery,CalibRecovery," ARM_ERR: CalibRecovery : vector of numerics expected",C_result);

	vector<double> LossRecoveryDefaultValue;
	LossRecoveryDefaultValue.resize(CalibRecovery.size());
	LossRecoveryDefaultValue = CalibRecovery;
	// Recov
	if (XL_Loss_Recovery->xltype == 0x0001)
	{
		XL_readNumCell(XL_Loss_Recovery,dLossRecovery,"ARM_ERR: LossRecovery : vector of numeric expected",C_result);
		LossRecovery.push_back(dLossRecovery);
		LossRecovery.push_back(dLossRecovery);
	}
	else
		XL_readNumVectorWD(XL_Loss_Recovery,LossRecovery,LossRecoveryDefaultValue," ARM_ERR: LossRecovery : vector of numerics expected",C_result);

	int Nbissuer = (int)DataMez[5];
	// BetaIssuers
	if (XL_IssuerBeta->xltype == 0x0001)
	{
		XL_readNumCell(XL_IssuerBeta,dIssuersBeta,"ARM_ERR: Beta : vector of numeric expected",C_result);
		for (int i = 0 ;i<Nbissuer;i++)
			IssuersBeta.push_back(dIssuersBeta);
	}
	else
		XL_readNumVector(XL_IssuerBeta, IssuersBeta," ARM_ERR: Beta : vector of numerics expected",C_result);
	
	
	int taille = (int)(DataMez[0]/DataMez[3]) + 5;
	
	// Spread	
	int size = YFPlot.size()*Nbissuer;
	vector<double> MatSpreadDefaultValue;
	for (int i = 0;i<size;i++) 
		MatSpreadDefaultValue.push_back(Credit_Def_Value);

	if (XL_MatIssuersSpread->xltype == 0x0001)
	{
		XL_readNumCellWD(XL_MatIssuersSpread,dIssuersSpread,Credit_Def_Value,"ARM_ERR: MatIssuersSpread : Matrix of numerics expected",C_result);
		for (i = 0;i<size;i++) 
			MatIssuersSpread.push_back(dIssuersSpread);
	}
	else
		XL_readNumVectorWithHoleWD(XL_MatIssuersSpread, MatIssuersSpread,MatSpreadDefaultValue," ARM_ERR: MatIssuersSpread : Matrix of numerics expected",C_result);

	// Labels
	vector<CCString> LabelDefaultValue;
	LabelDefaultValue.resize(Nbissuer);
	char* Lbl= new char[5];
	for (i = 0;i<Nbissuer;i++)
	{
		ITOA(i,Lbl,10);
		LabelDefaultValue[i] = Lbl;
	}
	
	if (Lbl)
		delete[] Lbl;
	Lbl = NULL;
	
	XL_readStrVectorWD(XL_IssuersLabel,IssuersLabel,LabelDefaultValue," ARM_ERR: Label : vector of string expected",DOUBLE_TYPE,C_result);


	// Date Shifte
	double dDateShifte = Credit_Def_Value;
	double DateShifte;
	
	XL_readNumCellWD(XL_DateToShift,DateShifte,dDateShifte," ARM_ERR: DateShifte : numeric expected",C_result);

	// Spread Courbe Forward
	if (XL_SpreadCrbFwd->xltype == 0x0001)
	{
		XL_readNumCellWD(XL_SpreadCrbFwd,dSpreadCrbFwd,Credit_Def_Value,"ARM_ERR: SpreadCrbFwd :  Vector of numerics expected",C_result);
		SpreadCrbFwd.push_back(dSpreadCrbFwd);
	}
	else
		XL_readNumVectorWithHoleWD(XL_SpreadCrbFwd,SpreadCrbFwd,SpreadCrbFwdDefaultValue," ARM_ERR: SpreadCrbFwd : numeric expected",C_result);
// -------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------
	long retCode;
	
	int MezType = 0;
	if ((MezType = ARM_ConvCreditMezzType(C_MezType,C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	retCode = ICMLOCAL_STR_CDO(MezType,
							   DataMez,
							   Notionals,
							   CalibRecovery,
							   LossRecovery,
							   YFPlot,
							   MatIssuersSpread,								
							   IssuersLabel,
							   IssuersBeta,
							   Taux,
							   LabelToShift,
							   DateShifte,
							   ValShift,
							   SpreadCrbFwd,
							   VolConvexity,
							   NumMethod,
							   C_Function,
							   C_result);
	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble ();
	}
	else
	{
		ARM_ERR();
	}
}


	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Price" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}

// -----------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------

_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Str_GetExpectedLoss(LPXLOPER XL_YearTerm,
																	LPXLOPER XL_Attachement,
																	LPXLOPER XL_NbIssuers,
																	LPXLOPER XL_Calib_Recovery,
																	LPXLOPER XL_Loss_Recovery,
																	LPXLOPER XL_YFPlot,
																	LPXLOPER XL_MatIssuersSpread,
																	LPXLOPER XL_IssuerBeta,
																	LPXLOPER XL_Taux,
																	LPXLOPER XL_Numerical_Method)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	double YearTerm;
	
	vector<double> Attachement;
	
	double dNb_Issuers;

	vector<double> CalibRecovery;
	double dCalibRecovery;

	vector<double> LossRecovery;
	double dLossRecovery;
	
	vector<double> YFPlot;
	double dYFPlot;
	
	vector<double> MatIssuersSpread;
	double dIssuersSpread;
			
	vector<double> IssuersBeta;
	double dIssuersBeta;

	double Taux;
	double TauxDefaultValue = 0.;

	vector<CCString> NumMethod;
	vector<CCString> NumMethodDefaultValue;
	NumMethodDefaultValue.resize(3);
	NumMethodDefaultValue[0] = "20";
	NumMethodDefaultValue[1] = "H";
	NumMethodDefaultValue[2] = "FAST";
	
		// error
	static int error;
	static char* reason = "";
	
	XL_readNumCell(XL_YearTerm,YearTerm," ARM_ERR: YearTerm : numeric expected",C_result);
	XL_readNumCell(XL_NbIssuers,dNb_Issuers,"ARM_ERR: NbIssuers : numeric expected",C_result);
	XL_readNumCellWD(XL_Taux,Taux,TauxDefaultValue," ARM_ERR: Taux : numeric expected",C_result);
	
	XL_readNumVector(XL_Attachement,Attachement," ARM_ERR: Attachement : Vector of numeric expected",C_result);

	// Data Nnum methode
	XL_readStrVectorWD(XL_Numerical_Method,NumMethod,NumMethodDefaultValue,"ARM_ERR : Numerical Method : Vecteur expcected",DOUBLE_TYPE, C_result);

	// YF plot
	if (XL_YFPlot->xltype == 0x0001)
	{
		XL_readNumCell(XL_YFPlot,dYFPlot,"ARM_ERR : Plot : Vector of numeric expected",C_result);
		YFPlot.push_back(dYFPlot);
	}
	else
		XL_readNumVector(XL_YFPlot,YFPlot," ARM_ERR: YFPlot : Vector numeric expected",C_result);

	// Recov
	if (XL_Calib_Recovery->xltype == 0x0001)
	{
		XL_readNumCell(XL_Calib_Recovery,dCalibRecovery,"ARM_ERR: CalibRecovery : vector of numeric expected",C_result);
		CalibRecovery.push_back(dCalibRecovery);
		CalibRecovery.push_back(dCalibRecovery);
	}
	else
		XL_readNumVector(XL_Calib_Recovery,CalibRecovery," ARM_ERR: CalibRecovery : vector of numerics expected",C_result);

	vector<double> LossRecoveryDefaultValue;
	LossRecoveryDefaultValue.resize(CalibRecovery.size());
	LossRecoveryDefaultValue = CalibRecovery;
	// Recov
	if (XL_Loss_Recovery->xltype == 0x0001)
	{
		XL_readNumCell(XL_Loss_Recovery,dLossRecovery,"ARM_ERR: LossRecovery : vector of numeric expected",C_result);
		LossRecovery.push_back(dLossRecovery);
		LossRecovery.push_back(dLossRecovery);
	}
	else
		XL_readNumVectorWD(XL_Loss_Recovery,LossRecovery,LossRecoveryDefaultValue," ARM_ERR: LossRecovery : vector of numerics expected",C_result);

	// BetaIssuers
	if (XL_IssuerBeta->xltype == 0x0001)
	{
		XL_readNumCell(XL_IssuerBeta,dIssuersBeta,"ARM_ERR: Beta : vector of numeric expected",C_result);
		for (int i = 0 ;i<(int)dNb_Issuers;i++)
			IssuersBeta.push_back(dIssuersBeta);
	}
	else
		XL_readNumVector(XL_IssuerBeta, IssuersBeta," ARM_ERR: Beta : vector of numerics expected",C_result);
	
	// Spread	
	int size = YFPlot.size()*(int)dNb_Issuers;
	vector<double> MatSpreadDefaultValue;
	for (int i = 0;i<size;i++) 
		MatSpreadDefaultValue.push_back(Credit_Def_Value);

	if (XL_MatIssuersSpread->xltype == 0x0001)
	{
		XL_readNumCellWD(XL_MatIssuersSpread,dIssuersSpread,Credit_Def_Value,"ARM_ERR: MatIssuersSpread : Matrix of numerics expected",C_result);
		for (i = 0;i<size;i++) 
			MatIssuersSpread.push_back(dIssuersSpread);
	}
	else
		XL_readNumVectorWithHoleWD(XL_MatIssuersSpread, MatIssuersSpread,MatSpreadDefaultValue," ARM_ERR: MatIssuersSpread : Matrix of numerics expected",C_result);
// -------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------
	long retCode;
		
	retCode = ICMLOCAL_STR_GetExpectedLoss(YearTerm,
										   Attachement,
										   dNb_Issuers,
										   CalibRecovery,
										   LossRecovery,
										   YFPlot,
										   MatIssuersSpread,								
										   IssuersBeta,
										   Taux,
										   NumMethod,
										   C_result);
	if(retCode == ARM_OK)
	{
		FreeCurCellErr ();
		XL_result.xltype = xltypeNum;
		XL_result.val.num = C_result.getDouble ();
	}
	else
	{
		ARM_ERR();
	}
}


	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Str_GetExpectedLoss" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}


// -----------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------
/**

_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Str_CdoSmile(LPXLOPER XL_MezType,
																 LPXLOPER XL_DataMez,
																 LPXLOPER XL_Recovery,
																 LPXLOPER XL_PlotPtf,
																 LPXLOPER XL_MatIssuersSpread,
																 LPXLOPER XL_IssuersLabel,
																 LPXLOPER XL_PlotTraxx,
																 LPXLOPER XL_SpreadTraxx,
																 LPXLOPER XL_MaturityBaseCorrel,
																 LPXLOPER XL_CorrelParStrike,
																 LPXLOPER XL_DataTraxx,
																 LPXLOPER XL_Taux,
																 LPXLOPER XL_LabelToShift,
																 LPXLOPER XL_DateToShift,
																 LPXLOPER XL_ValShift,
																 LPXLOPER XL_SpreadCrbFwd,
																 LPXLOPER XL_VolConvexity,
																 LPXLOPER XL_Numerical_Method,
																 LPXLOPER XL_Fonction)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;
	LPXLOPER pxArray;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	CCString C_MezType;

	vector<double> DataMez;
	vector<double> DataMezDefaultValue;
	DataMezDefaultValue.resize(7);
	DataMezDefaultValue[0] = 5.;		// maturity
	DataMezDefaultValue[1] = 0.;		// Strike Down
	DataMezDefaultValue[2] = 0.2;		// Strike Up 
	DataMezDefaultValue[3] = 0.25;		// Pas
	DataMezDefaultValue[4] = 0.01;		// Fixed Rate
	DataMezDefaultValue[5] = 125;		// nb nom
	DataMezDefaultValue[6] = 10000000;  // Traded Notional
	
	vector<double> Recovery;
	double dRecovery;

	vector<double> PlotPtf;
	double dPlotPtf;
	
	vector<double> MatIssuersSpread;
	double dIssuersSpread;
	
	vector<CCString> IssuersLabel;

	vector<double> PlotTraxx;
	double dPlotTraxx;

	vector<double> SpreadTraxx;
	vector<double> SpreadTraxxDefaultValue;
	SpreadTraxxDefaultValue.push_back(Credit_Def_Value);

	vector<double> MaturityBaseCorrel;
	vector<double> MaturityBaseCorrelDefaultValue;
	MaturityBaseCorrelDefaultValue.resize(2);
	MaturityBaseCorrelDefaultValue[0] = 5.;
	MaturityBaseCorrelDefaultValue[1] = 10.;
	double dMatBaseCorrel;
	double dMaturityBaseCorrelDefaultValue = 5.;

	long nbRowCorrelation=0;
	long nbColCorrelation=0;
	vector<CCString> CorrelParStrike;
	CorrelParStrike.clear();

	vector<double> DataTraxx ;
	vector<double> DataTraxxDefaultValue ;
	DataTraxxDefaultValue.push_back(1);		// Projection on Traxx Eur
	DataTraxxDefaultValue.push_back(0.4);	// Recovery Traxx

	double Taux;
	double TauxDefaultValue = 0.;

	CCString LabelToShift;
	CCString LabelToShiftDefaultValue = "NONE";

	double ValShift = 0.;
	double ValShiftDefaultValue = 0.01/100.;

	vector<double> SpreadCrbFwd;
	vector<double> SpreadCrbFwdDefaultValue;
	SpreadCrbFwdDefaultValue.push_back(Credit_Def_Value);
	double dSpreadCrbFwd;

	double VolConvexity;
	double VolConvexityDefaultValue = 0.;

	vector<CCString> NumMethod;
	vector<CCString> NumMethodDefaultValue;
	NumMethodDefaultValue.resize(3);
	NumMethodDefaultValue[0] = "20";
	NumMethodDefaultValue[1] = "H";
	NumMethodDefaultValue[2] = "FAST";

	CCString C_Function;

	// error
	static int error;
	static char* reason = "";

	// Mez Type
	XL_readStrCellWD(XL_MezType,C_MezType,"RUNNING"," ARM_ERR: MezType : String expected",C_result);

	// Data Mez
	XL_readNumVectorWD(XL_DataMez,DataMez,DataMezDefaultValue," ARM_ERR: DataMez : Vector of numeric expected",C_result);
	
	
	// Taux
	XL_readNumCellWD(XL_Taux,Taux,TauxDefaultValue," ARM_ERR: Taux : numeric expected",C_result);

	// Fonction
	XL_readStrCell(XL_Fonction,C_Function," ARM_ERR: C_Function : String expected",C_result);

	// Label To shift
	XL_readStrCellWD(XL_LabelToShift,LabelToShift,LabelToShiftDefaultValue," ARM_ERR: LabelToShift : String expected",C_result);
	
	// Val Shifte
	XL_readNumCellWD(XL_ValShift,ValShift,ValShiftDefaultValue,"ARM_ERR: ValShift: numeric expected",C_result);

	// Date Shifte
	double dDateShifte = Credit_Def_Value;
	double DateShifte;
	
	XL_readNumCellWD(XL_DateToShift,DateShifte,dDateShifte," ARM_ERR: DateShifte : numeric expected",C_result);


	// Data Nnum methode
	XL_readStrVectorWD(XL_Numerical_Method,NumMethod,NumMethodDefaultValue,"ARM_ERR : Numerical Method : Vecteur expcected",DOUBLE_TYPE, C_result);
	
	// SpreadTraxx
	XL_readNumVectorWD(XL_SpreadTraxx,SpreadTraxx,SpreadTraxxDefaultValue," ARM_ERR: SpreadTraxx : Vector of numeric expected",C_result);

	// Maturity Base Correl
	if (XL_MaturityBaseCorrel->xltype == 0x0001)
	{
		XL_readNumCellWD(XL_MaturityBaseCorrel,dMatBaseCorrel,dMaturityBaseCorrelDefaultValue,"ARM_ERR: MaturityBaseCorrel : vector of numeric expected",C_result);
		MaturityBaseCorrel.push_back(dMatBaseCorrel);
	}
	else
		XL_readNumVectorWD(XL_MaturityBaseCorrel,MaturityBaseCorrel,MaturityBaseCorrelDefaultValue," ARM_ERR: MaturityBaseCorrel : vector of numerics expected",C_result);

	// Correl Par strike
	XL_readStrVectorAndSize(XL_CorrelParStrike,nbRowCorrelation,nbColCorrelation,CorrelParStrike," ARM_ERR: CorrelParStrike : object matrix expected",DOUBLE_TYPE,C_result);
	int nbRowCorrel = (int) nbRowCorrelation;
	int nbColCorrel = (int) nbColCorrelation;	

	double sizeCorrelParStrike = CorrelParStrike.size();

	if(sizeCorrelParStrike == 0)
	{
		C_result.setMsg ("ARM_ERR: check your CorrelParStrike Matrix");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	// Data Traxx
	XL_readNumVectorWD(XL_DataTraxx,DataTraxx,DataTraxxDefaultValue," ARM_ERR: DataTraxx : Vector of numeric expected",C_result);

	// PlotPtf
	if (XL_PlotPtf->xltype == 0x0001)
	{
		XL_readNumCell(XL_PlotPtf,dPlotPtf,"ARM_ERR : PlotPtf : Vector of numeric expected",C_result);
		PlotPtf.push_back(dPlotPtf);
	}
	else
		XL_readNumVector(XL_PlotPtf,PlotPtf," ARM_ERR: PlotPtf : Vector numeric expected",C_result);

	// PlotTraxx
	if (XL_PlotTraxx->xltype == 0x0001)
	{
		XL_readNumCell(XL_PlotTraxx,dPlotTraxx,"ARM_ERR : PlotTraxx : Vector of numeric expected",C_result);
		PlotTraxx.push_back(dPlotTraxx);
	}
	else
		XL_readNumVector(XL_PlotTraxx,PlotTraxx," ARM_ERR: PlotTraxx : Vector numeric expected",C_result);

	// Recov
	if (XL_Recovery->xltype == 0x0001)
	{
		XL_readNumCell(XL_Recovery,dRecovery,"ARM_ERR: Recovery : vector of numeric expected",C_result);
		Recovery.push_back(dRecovery);
		Recovery.push_back(dRecovery);
	}
	else
		XL_readNumVector(XL_Recovery,Recovery," ARM_ERR: Recovery : vector of numerics expected",C_result);


	int taille = (int)(DataMez[0]/DataMez[4]) + 5;
	
	// Spread	
	int NbIssuers = (int)DataMez[5];
	int size = PlotPtf.size()*NbIssuers;
	vector<double> MatSpreadDefaultValue;
	for (int i = 0;i<size;i++) 
		MatSpreadDefaultValue.push_back(Credit_Def_Value);

	if (XL_MatIssuersSpread->xltype == 0x0001)
	{
		XL_readNumCellWD(XL_MatIssuersSpread,dIssuersSpread,Credit_Def_Value,"ARM_ERR: MatIssuersSpread : Matrix of numerics expected",C_result);
		for (i = 0;i<size;i++) 
			MatIssuersSpread.push_back(dIssuersSpread);
	}
	else
		XL_readNumVectorWithHoleWD(XL_MatIssuersSpread, MatIssuersSpread,MatSpreadDefaultValue," ARM_ERR: MatIssuersSpread : Matrix of numerics expected",C_result);

	// Labels
	int nbLabel = NbIssuers;
	vector<CCString> LabelDefaultValue;
	LabelDefaultValue.resize(nbLabel);
	char* Lbl= new char[5];
	for (i = 0;i<nbLabel;i++)
	{
		ITOA(i,Lbl,10);
		LabelDefaultValue[i] = Lbl;
	}
	
	if (Lbl)
		delete[] Lbl;
	Lbl = NULL;
	
	XL_readStrVectorWD(XL_IssuersLabel,IssuersLabel,LabelDefaultValue," ARM_ERR: Label : vector of string expected",DOUBLE_TYPE,C_result);

	// Spread Courbe Forward
	if (XL_SpreadCrbFwd->xltype == 0x0001)
	{
		XL_readNumCellWD(XL_SpreadCrbFwd,dSpreadCrbFwd,Credit_Def_Value,"ARM_ERR: SpreadCrbFwd :  Vector of numerics expected",C_result);
		SpreadCrbFwd.push_back(dSpreadCrbFwd);
	}
	else
		XL_readNumVectorWithHoleWD(XL_SpreadCrbFwd,SpreadCrbFwd,SpreadCrbFwdDefaultValue," ARM_ERR: SpreadCrbFwd : numeric expected",C_result);

	// Vol convexity
	XL_readNumCellWD(XL_VolConvexity,VolConvexity,VolConvexityDefaultValue," ARM_ERR: VolConvexity : numeric expected",C_result);
// -------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------
	long retCode;

	// On gere la frequency
	double FreqTemp = 0.;
	if ( (fabs(DataMez[3])>1E-6) && (fabs(DataMez[3] - DataMez[0])>1E-6) )
		FreqTemp = 1./DataMez[3];

	DataMez[3] = FreqTemp;
	
	// On gere le type de la tranche
	int MezType = 0;
	if ((MezType = ARM_ConvCreditMezzType(C_MezType,C_result)) == ARM_DEFAULT_ERR)
	{
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}
	
	vector<double> Resultat;

	retCode = ICMLOCAL_STR_CDOSMILE(MezType,
									DataMez,
									Recovery,
									PlotPtf,
									MatIssuersSpread,								
									IssuersLabel,
									PlotTraxx,
									SpreadTraxx,
									MaturityBaseCorrel,
									CorrelParStrike,
									nbRowCorrel,
									nbColCorrel,
									DataTraxx,
									Taux,
									LabelToShift,
									DateShifte,
									ValShift,
									SpreadCrbFwd,
									VolConvexity,
									NumMethod,
									C_Function,
									Resultat,
									C_result);

if(retCode == ARM_OK)
{
	int nbrows = 1;
	int nbcolumns = Resultat.size();
	
	FreeCurCellErr ();
	XL_result.xltype = xltypeMulti;
	XL_result.val.array.columns = nbcolumns;
	XL_result.val.array.rows = nbrows; 
	XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

	for(int i = 0; i < nbcolumns; i++)
	{
		pxArray[XL_Coordonnate2Rank (1, 0, i)].xltype = xltypeNum;
		pxArray[XL_Coordonnate2Rank (1, 0, i)].val.num = Resultat[i]; 
	}
}
else
{
	ARM_ERR();
}

}


	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Str_CDOSMILESHOCKED" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
**/ 

// -----------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------

_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Str_GetExpectedLossSMILE(LPXLOPER XL_DataMez,
																	     LPXLOPER XL_Calib_Recovery,
																	     LPXLOPER XL_Loss_Recovery,
																	     LPXLOPER XL_YFPlot,
																	     LPXLOPER XL_MatIssuersSpread,
																	     LPXLOPER XL_SpreadTraxx,
																		 LPXLOPER XL_MaturityBaseCorrel,
																	     LPXLOPER XL_CorrelParStrike,
																	     LPXLOPER XL_DataTraxx,
																	     LPXLOPER XL_Numerical_Method)
{
	// return
	static XLOPER XL_result;
	ARM_result C_result;
	LPXLOPER pxArray;

	/// to make the infrastructure very robust we put a try catch!
	ARM_XL_TRY_BLOCK_BEGIN
	{
	ARM_NOCALCIFWIZ();

	// C variable
	vector<double> DataMez;
	vector<double> DataMezDefaultValue;
	DataMezDefaultValue.resize(6);
	DataMezDefaultValue[0] = 5.;   // maturity
	DataMezDefaultValue[1] = 0.;   // Strike Down
	DataMezDefaultValue[2] = 0.2;  // Strike Up 
	DataMezDefaultValue[3] = 0.25; // Pas
	DataMezDefaultValue[4] = 0.;   // Taux
	DataMezDefaultValue[5] = 125;  // NbNom
	

	vector<double> CalibRecovery;
	double dCalibRecovery;

	vector<double> LossRecovery;
	double dLossRecovery;
	
	vector<double> YFPlot;
	double dYFPlot;
	
	vector<double> MatIssuersSpread;
	double dIssuersSpread;
	
	vector<double> SpreadTraxx;
	vector<double> SpreadTraxxDefaultValue;
	SpreadTraxxDefaultValue.push_back(Credit_Def_Value);

	vector<double> MaturityBaseCorrel;
	vector<double> MaturityBaseCorrelDefaultValue;
	MaturityBaseCorrelDefaultValue.resize(2);
	MaturityBaseCorrelDefaultValue[0] = 5.;
	MaturityBaseCorrelDefaultValue[1] = 10.;
	double dMatBaseCorrel;
	double dMaturityBaseCorrelDefaultValue = 5.;


	long nbRowCorrelation=0;
	long nbColCorrelation=0;
	vector<CCString> CorrelParStrike;
	CorrelParStrike.clear();

	vector<double> DataTraxx ;
	vector<double> DataTraxxDefaultValue ;
	DataTraxxDefaultValue.push_back(1);		// Projection on Traxx Eur
	DataTraxxDefaultValue.push_back(0.4);	// Recovery Traxx

	vector<CCString> NumMethod;
	vector<CCString> NumMethodDefaultValue;
	NumMethodDefaultValue.resize(3);
	NumMethodDefaultValue[0] = "20";
	NumMethodDefaultValue[1] = "H";
	NumMethodDefaultValue[2] = "FAST";
	
	// error
	static int error;
	static char* reason = "";
	
	XL_readNumVectorWD(XL_DataMez,DataMez,DataMezDefaultValue," ARM_ERR: DataMez : Vector of numeric expected",C_result);
	
	// Data Nnum methode
	XL_readStrVectorWD(XL_Numerical_Method,NumMethod,NumMethodDefaultValue,"ARM_ERR : Numerical Method : Vecteur expcected",DOUBLE_TYPE, C_result);

	// Spread Traxx
	XL_readNumVectorWD(XL_SpreadTraxx,SpreadTraxx,SpreadTraxxDefaultValue," ARM_ERR: SpreadTraxx : Vector of numeric expected",C_result);
	

	// Maturity Base Correl
	if (XL_MaturityBaseCorrel->xltype == 0x0001)
	{
		XL_readNumCellWD(XL_MaturityBaseCorrel,dMatBaseCorrel,dMaturityBaseCorrelDefaultValue,"ARM_ERR: MaturityBaseCorrel : vector of numeric expected",C_result);
		MaturityBaseCorrel.push_back(dMatBaseCorrel);
	}
	else
		XL_readNumVectorWD(XL_MaturityBaseCorrel,MaturityBaseCorrel,MaturityBaseCorrelDefaultValue," ARM_ERR: MaturityBaseCorrel : vector of numerics expected",C_result);

	

	// Base Correl
	XL_readStrVectorAndSize(XL_CorrelParStrike,nbRowCorrelation,nbColCorrelation,CorrelParStrike," ARM_ERR: CorrelParStrike : object matrix expected",DOUBLE_TYPE,C_result);
	int nbRowCorrel = (int) nbRowCorrelation;
	int nbColCorrel = (int) nbColCorrelation;	

	double sizeCorrelParStrike = CorrelParStrike.size();

	if(sizeCorrelParStrike == 0)
	{
		C_result.setMsg ("ARM_ERR: check your CorrelParStrike Matrix");
		ARM_ARG_ERR();
		return (LPXLOPER)&XL_result;
	}

	XL_readNumVectorWD(XL_DataTraxx,DataTraxx,DataTraxxDefaultValue," ARM_ERR: DataTraxx : Vector of numeric expected",C_result);

	// YF plot
	if (XL_YFPlot->xltype == 0x0001)
	{
		XL_readNumCell(XL_YFPlot,dYFPlot,"ARM_ERR : Plot : Vector of numeric expected",C_result);
		YFPlot.push_back(dYFPlot);
	}
	else
		XL_readNumVector(XL_YFPlot,YFPlot," ARM_ERR: YFPlot : Vector numeric expected",C_result);

	// Recov
	if (XL_Calib_Recovery->xltype == 0x0001)
	{
		XL_readNumCell(XL_Calib_Recovery,dCalibRecovery,"ARM_ERR: CalibRecovery : vector of numeric expected",C_result);
		CalibRecovery.push_back(dCalibRecovery);
		CalibRecovery.push_back(dCalibRecovery);
	}
	else
		XL_readNumVector(XL_Calib_Recovery,CalibRecovery," ARM_ERR: CalibRecovery : vector of numerics expected",C_result);

	// Recov
	vector<double> LossRecoveryDefaultValue;
	LossRecoveryDefaultValue.resize(CalibRecovery.size());
	LossRecoveryDefaultValue = CalibRecovery;

	if (XL_Loss_Recovery->xltype == 0x0001)
	{
		XL_readNumCell(XL_Loss_Recovery,dLossRecovery,"ARM_ERR: LossRecovery : vector of numeric expected",C_result);
		LossRecovery.push_back(dLossRecovery);
		LossRecovery.push_back(dLossRecovery);
	}
	else
		XL_readNumVectorWD(XL_Loss_Recovery,LossRecovery,LossRecoveryDefaultValue," ARM_ERR: LossRecovery : vector of numerics expected",C_result);
	
	int taille = (int)(DataMez[0]/DataMez[3]) + 5;
	
	// Spread	
	int NbIssuers = DataMez[5];
	int size = YFPlot.size()*NbIssuers;
	vector<double> MatSpreadDefaultValue;
	for (int i = 0;i<size;i++) 
		MatSpreadDefaultValue.push_back(Credit_Def_Value);

	if (XL_MatIssuersSpread->xltype == 0x0001)
	{
		XL_readNumCellWD(XL_MatIssuersSpread,dIssuersSpread,Credit_Def_Value,"ARM_ERR: MatIssuersSpread : Matrix of numerics expected",C_result);
		for (i = 0;i<size;i++) 
			MatIssuersSpread.push_back(dIssuersSpread);
	}
	else
		XL_readNumVectorWithHoleWD(XL_MatIssuersSpread, MatIssuersSpread,MatSpreadDefaultValue," ARM_ERR: MatIssuersSpread : Matrix of numerics expected",C_result);
// -------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------
	long retCode;

	// On gere la frequency
	double FreqTemp = 0.;
	if ( (fabs(DataMez[3])>1E-6) && (fabs(DataMez[3] - DataMez[0])>1E-6) )
		FreqTemp = 1./DataMez[3];

	DataMez[3] = FreqTemp;
	
	vector<double> Resultat;

	retCode = ICMLOCAL_STR_ExpectedLossSMILE(DataMez,
											 CalibRecovery,
											 LossRecovery,
											 YFPlot,
											 MatIssuersSpread,								
											 SpreadTraxx,
											 MaturityBaseCorrel,
											 CorrelParStrike,
											 nbRowCorrel,
											 nbColCorrel,
											 DataTraxx,
											 NumMethod,
											 Resultat,
											 C_result);

if(retCode == ARM_OK)
{
	int nbrows = Resultat.size();
	int nbcolumns = 1;
	
	FreeCurCellErr ();
	XL_result.xltype = xltypeMulti;
	XL_result.val.array.columns = nbcolumns;
	XL_result.val.array.rows = nbrows; 
	XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

	for(int i = 0; i < nbrows; i++)
	{
		pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].xltype = xltypeNum;
		pxArray[XL_Coordonnate2Rank (i, 0, nbcolumns)].val.num = Resultat[i]; 
	}
}
else
{
	ARM_ERR();
}

}


	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Str_CDOSMILE" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}



// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
/** 
__declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Str_LinearInterpolation (LPXLOPER XL_Xref,
																		  LPXLOPER XL_Yref,
																		  LPXLOPER XL_Xtarget)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;
	LPXLOPER pxArray;

	/// to make the infrastructure very robust we put a try catch!
ARM_XL_TRY_BLOCK_BEGIN
{
	ARM_NOCALCIFWIZ();

	// C variable
	vector<double> Xref;
	vector<double> Yref;
	
	vector<double> Xtarget;
	vector<double> Ytarget;
	
	// error
	static int error;
	static char* reason = "";
	
	XL_readNumVector(XL_Xref, Xref," ARM_ERR: Xref : vector of numerics expected",C_result);
	XL_readNumVector(XL_Yref, Yref," ARM_ERR: Yref : vector of numerics expected",C_result);
	XL_readNumVector(XL_Xtarget, Xtarget," ARM_ERR: Xtarget : vector of numerics expected",C_result);

	Ytarget.resize(Xtarget.size());

	long retCode;
	retCode = ICMLOCAL_STR_LinearInterpol(Xref,
									      Yref,
										  Xtarget,
										  Ytarget,
										  C_result);
	
	if(retCode == ARM_OK)
	{
		int nbrows = 1;
		int nbcolumns = Ytarget.size();
		
		FreeCurCellErr ();
		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		for(int i = 0; i < nbcolumns; i++)
		{
			pxArray[XL_Coordonnate2Rank (nbrows, 0, i)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (nbrows, 0, i)].val.num = Ytarget[i]; 
		}
	}
	else
	{
		ARM_ERR();
	}
}



	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Price" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
**/ 
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
// MATRICE DE CORREL INTER TRANCHE
/*
_declspec(dllexport) LPXLOPER WINAPI ARM_Credit_Str_CorrelInterTranche (LPXLOPER XL_NbSimul,
																		LPXLOPER XL_NbTranche,
																		LPXLOPER XL_NbName,
																		LPXLOPER XL_BetaIssuers,
																		LPXLOPER XL_Calib_Recovery,
																		LPXLOPER XL_Loss_Recovery,
																		LPXLOPER XL_Maturity,
																		LPXLOPER XL_StrikeDownUp,
																		LPXLOPER XL_NameInCDO,
																		LPXLOPER XL_YFPlot,
																		LPXLOPER XL_MatIssuersSpread,
																		LPXLOPER XL_Taux)
{
//	ARM_BEGIN();

	// return
	static XLOPER XL_result;
	ARM_result C_result;
	LPXLOPER pxArray;

	/// to make the infrastructure very robust we put a try catch!
ARM_XL_TRY_BLOCK_BEGIN
{
	ARM_NOCALCIFWIZ();

	// C variable
	double dNbSimul;
	double dNbTranche;
	double dNbName;
	double BetaIssuers;
	vector<double> CalibRecovery;
	double dCalibRecovery;
	double LossRecovery;
	double Maturity;
	
	vector<double> YFPlot;
	double dYFPlot;

	vector<double> MatIssuersSpread;
	double dIssuersSpread;

	double Taux;
	double TauxDefaultValue = 0.;
	
	long lnbRowStrikeDownUp=0;
	long lnbColStrikeDownUp=0;
	vector<CCString> StrikeDownUp;
	StrikeDownUp.clear();

	long nbRowNameInCDO=0;
	long nbColNameInCDO=0;
	vector<CCString> NameInCDO;
	NameInCDO.clear();


	// error
	static int error;
	static char* reason = "";
	
	XL_readNumCell(XL_NbSimul,dNbSimul," ARM_ERR: NbSimul : numeric expected",C_result);
	XL_readNumCell(XL_NbTranche,dNbTranche," ARM_ERR: NbTranche : numeric expected",C_result);
	XL_readNumCell(XL_NbName,dNbName," ARM_ERR: NbName : numeric expected",C_result);
	XL_readNumCell(XL_Maturity,Maturity," ARM_ERR: Maturity : numeric expected",C_result);
	XL_readNumCell(XL_BetaIssuers, BetaIssuers," ARM_ERR: BetaIssuers : numeric expected",C_result);
	XL_readNumCellWD(XL_Taux,Taux,TauxDefaultValue," ARM_ERR: Taux : numeric expected",C_result);
		
	// Spread	
	int size = YFPlot.size()*(int)dNbName;
	vector<double> MatSpreadDefaultValue;
	for (int i = 0;i<size;i++) 
		MatSpreadDefaultValue.push_back(Credit_Def_Value);

	if (XL_MatIssuersSpread->xltype == 0x0001)
	{
		XL_readNumCellWD(XL_MatIssuersSpread,dIssuersSpread,Credit_Def_Value,"ARM_ERR: MatIssuersSpread : Matrix of numerics expected",C_result);
		for (i = 0;i<size;i++) 
			MatIssuersSpread.push_back(dIssuersSpread);
	}
	else
		XL_readNumVectorWithHoleWD(XL_MatIssuersSpread, MatIssuersSpread,MatSpreadDefaultValue," ARM_ERR: MatIssuersSpread : Matrix of numerics expected",C_result);
// -------------------------------------------------------------------------------------------------------------------



	XL_readStrVectorAndSize(XL_StrikeDownUp,lnbRowStrikeDownUp,lnbColStrikeDownUp,StrikeDownUp," ARM_ERR: StrikeDownUp : object matrix expected",DOUBLE_TYPE,C_result);
	int nbColStrikeDownUp = (int) lnbColStrikeDownUp;

	XL_readStrVectorAndSize(XL_NameInCDO,nbRowNameInCDO,nbColNameInCDO,NameInCDO," ARM_ERR: NameInCDO : object matrix expected",DOUBLE_TYPE,C_result);
	int nbColIssuersInCDO = (int) nbColNameInCDO;

	// YF plot
	if (XL_YFPlot->xltype == 0x0001)
	{
		XL_readNumCell(XL_YFPlot,dYFPlot,"ARM_ERR : Plot : Vector of numeric expected",C_result);
		YFPlot.push_back(dYFPlot);
	}
	else
		XL_readNumVector(XL_YFPlot,YFPlot," ARM_ERR: YFPlot : Vector numeric expected",C_result);


	// Recovery de Calibration
	if (XL_Calib_Recovery->xltype == 0x0001)
	{
		XL_readNumCell(XL_Calib_Recovery,dCalibRecovery,"ARM_ERR: CalibRecovery : vector of numeric expected",C_result);
		CalibRecovery.push_back(dCalibRecovery);
		CalibRecovery.push_back(dCalibRecovery);
	}
	else
		XL_readNumVector(XL_Calib_Recovery,CalibRecovery," ARM_ERR: CalibRecovery : vector of numerics expected",C_result);

	// Recov pour calculer la lossunit
	double LossRecoveryDefaultValue = CalibRecovery[0];
	XL_readNumCellWD(XL_Loss_Recovery,LossRecovery,LossRecoveryDefaultValue," ARM_ERR: LossRecovery : numerics expected",C_result);


	vector<double*> MatriceEvtDefaut;
	vector<double*> MatriceTpsDefaut;

	ARM_Vector  ProbaDefautTranche;
	ARM_Vector VarianceTranche;

	long NbSimul  = (long)dNbSimul;
	int NbTranche = (int)dNbTranche;
	int NbName	  = (int)dNbName;

	long retCode;
	retCode = ICMLOCAL_Str_CorrelInterTranche(NbSimul,
											  NbTranche,
											  NbName,
											  BetaIssuers,
											  CalibRecovery,
											  LossRecovery,
											  Maturity,
											  StrikeDownUp,
											  nbColStrikeDownUp,
											  NameInCDO,
											  nbColIssuersInCDO,
											  YFPlot,
											  MatIssuersSpread,
											  Taux,
											  ProbaDefautTranche,
											  VarianceTranche,
											  MatriceEvtDefaut,
											  MatriceTpsDefaut,
											  C_result);
	
	if(retCode == ARM_OK)
	{
		int nbrows =  2*NbTranche + 4;
		int nbcolumns = NbTranche +2;		// On rajoute un vecteur (proba de défaut + variance) + 1 colonne pour séparer les matrices des vecteurs colonnes
		
		FreeCurCellErr ();
		XL_result.xltype = xltypeMulti;
		XL_result.val.array.columns = nbcolumns;
		XL_result.val.array.rows = nbrows; 
		XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbrows * nbcolumns * sizeof (XLOPER));

		// On met du blanc partout
		for (int i = 0;i<nbrows;i++)
		{
			for (int j = 0;j<nbcolumns;j++)
			{
				pxArray[XL_Coordonnate2Rank (i,j, nbcolumns)].xltype = xltypeStr;
				pxArray[XL_Coordonnate2Rank (i,j, nbcolumns)].val.str = XL_StrC2StrPascal ("");
			}
		}
				
		// En tete de la matrice Evenement de défaut
		pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].xltype = xltypeStr;
		pxArray[XL_Coordonnate2Rank (0, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("Correl Inter Tranche Evenement Defaut");
		
		// Matrice Evenement de défaut
		for(i = 1; i < NbTranche + 1; i++)
		{
			for (int j=0;j<nbcolumns-2;j++)
			{
				pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].xltype = xltypeNum;
				pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].val.num = MatriceEvtDefaut[i-1][j]; 
			}
		}
	
		// En tete du vecteur proba de defaut
		pxArray[XL_Coordonnate2Rank (0, nbcolumns-1, nbcolumns)].xltype = xltypeStr;
		pxArray[XL_Coordonnate2Rank (0, nbcolumns-1, nbcolumns)].val.str = XL_StrC2StrPascal ("Proba Defaut Tranche");
		
		// Vecteur proba de défaut
		for(i = 1; i < NbTranche + 1; i++)
		{
			pxArray[XL_Coordonnate2Rank (i, nbcolumns-1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (i, nbcolumns-1, nbcolumns)].val.num = ProbaDefautTranche[i-1]; 
		}
	
		// En tete de la matrice Temps de défaut
		pxArray[XL_Coordonnate2Rank (NbTranche + 2, 0, nbcolumns)].xltype = xltypeStr;
		pxArray[XL_Coordonnate2Rank (NbTranche + 2, 0, nbcolumns)].val.str = XL_StrC2StrPascal ("Correl Inter tranche Temps Defaut");
		
		// Matrice Tps de défaut
		for(i = NbTranche + 3; i < nbrows-1; i++)
		{
			for (int j=0;j<nbcolumns-2;j++)
			{
				pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].xltype = xltypeNum;
				pxArray[XL_Coordonnate2Rank (i, j, nbcolumns)].val.num = MatriceTpsDefaut[i-NbTranche - 3][j]; 
			}
		}

		// En tete du vecteur de Variance
		pxArray[XL_Coordonnate2Rank (NbTranche + 2, nbcolumns-1, nbcolumns)].xltype = xltypeStr;
		pxArray[XL_Coordonnate2Rank (NbTranche + 2, nbcolumns-1, nbcolumns)].val.str = XL_StrC2StrPascal ("Variance");
		
		// Variance
		for(i = NbTranche + 3; i < nbrows-1; i++)
		{
			pxArray[XL_Coordonnate2Rank (i, nbcolumns-1, nbcolumns)].xltype = xltypeNum;
			pxArray[XL_Coordonnate2Rank (i, nbcolumns-1, nbcolumns)].val.num = VarianceTranche[i-NbTranche - 3]; 
		}
	}
	else
	{
		ARM_ERR();
	}
}



	/// end of try block
	ARM_XL_TRY_BLOCK_END
	
	/// to catch arm exception
	ARM_XL_CATCH_ARM_EXPT

	/// to cath all the other exceptions
	ARM_XL_CATCH_GAL_EXPT_WITH_MSSG( "Unrecognized failure in ARM_Credit_Price" )

	/// return the result as an LPXLOPER
	return (LPXLOPER)&XL_result;
}
*/

// EOF %M% 


extern "C" __declspec(dllexport) void WINAPI xlAutoFree (LPXLOPER px)
{
	
	// only for strings
	if(px)
	{
		px->xltype &= ~xlbitDLLFree ;
		if (px->xltype==xltypeStr) 
		{
			if(px->val.str)
			{
				free (px->val.str);
				px->val.str = NULL;
			}
		}
	if (px->xltype==xltypeMulti) 
	{
		unsigned int rows,cols; 
		unsigned int i; 
		//	Deleting an array.
		XLOPER* rgx=(XLOPER*)px->val.array.lparray ;
		if (!rgx) return ;			
		if (px->val.array.rows*px->val.array.columns==0) return; 
		rows = px->val.array.rows; 
		cols = px->val.array.columns;
		//	checking element type and deleteling if string
		for (i=0;i<rows*cols;i++) 
		{
			if (rgx[i].xltype==xltypeStr) 
			{
				char* ptr=rgx[i].val.str ;
				free( ptr) ;
			}
		}
		delete [] rgx ; 
	}

	return;
}
}
