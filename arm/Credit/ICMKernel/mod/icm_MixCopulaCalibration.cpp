#error no longer part of the project
#include "ARMKernel/glob/firsttoinc.h"
#include "icm_MixCopulaCalibration.h"
#include "ICMKernel/mod/modelmulticurves.h"
#include <nag.h>
#include <nage04.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include "icm_MixCopulaCalibration.h"


static int Count ;
static double Value0, V0, V1 ;


/***************************************************************************************************/
// Calibration du modèle de Mix Copula avec Proportions Full et indep
/***************************************************************************************************/

/*---------------------------------------------------------------------------------------*/
// Construction/Destruction
/*---------------------------------------------------------------------------------------*/

void MixCopulaCalibration::Init()
{
	itsImplied_Proportion_independant = 0. ;
	itsImplied_Proportion_FullCorrel = 0. ;
	itsPricer= NULL;
	itsPtf = NULL;
	itsSmiledPrice1 = 0. ;
	itsSmiledPrice2 = 0. ;
	itsFirstPrice = 0.;
}

MixCopulaCalibration::MixCopulaCalibration(ICM_Pricer* pricer, 
											ICM_Portfolio* ptf, 
											int BCType, 
											double BC1, 
											double BC2, 
											double Seed1, 
											double Seed2)
{
	Init();

	itsBCType = BCType ;
	itsBaseCorrelation1 = BC1 ;
	itsBaseCorrelation2 = BC2 ;
	itsSeed1 = Seed1;
	itsSeed2 = Seed2;
	itsPricer = pricer;
	itsPtf = ptf;
}

MixCopulaCalibration::~MixCopulaCalibration() {}

void MixCopulaCalibration::Copy(const ARM_Object* src)
{
	ARM_Object::Copy(src);
     BitwiseCopy(src);
}

/*---------------------------------------------------------------------------------------*/
// Cpt Price pour une tranche étant donnée le smile de correl
/*---------------------------------------------------------------------------------------*/

double MixCopulaCalibration::CptPrice_BC(ARM_Security* Sec,double BC, int BCType)
{
	itsPricer->ResetPricer();
	itsPricer->SetInitialPriceFlg(false);

	itsPricer->SetSecurity(Sec);
	ICM_ModelMultiCurves* DefModel = (ICM_ModelMultiCurves*)itsPricer->GetModel();
	ICM_ModelMultiCurves* DefModel2 = (ICM_ModelMultiCurves*)itsPricer->GetModel();

	ARM_Date ExecutionDate = DefModel->GetStartDate();
	
	qSENSITIVITY_TYPE typesensi ;	
	
	double Result = 0. ;

	double price = itsPricer->ComputePrice(qCMPPRICE);

	if(BCType == 0) // Correlation as Input
	{
		typesensi = ICM_SAMECORRELATION_TYPE;
	}
	else // Beta Smile as input
	{
		typesensi = ICMBETA_TYPE;
	}

	Result = price + itsPricer->Hedge(typesensi,"NONE","NONE",BC);
	
	itsPricer->ResetPricer();
	itsPricer->SetModel(DefModel2) ;

	return Result;
}

/*---------------------------------------------------------------------------------------*/
// Cpt Price pour une tranche étant donnés les paramètres du modèle de Mix Copula
/*---------------------------------------------------------------------------------------*/

double MixCopulaCalibration::CptPrice_MixCopula(ARM_Security* Sec,double proportion_independant,double proportion_full_correl)
{
	itsPricer->ResetPricer();
	itsPricer->SetSecurity(Sec);
	ICM_ModelMultiCurves* DefModel = (ICM_ModelMultiCurves*)itsPricer->GetModel();

	ARM_Date ExecutionDate = DefModel->GetStartDate();

	double aux1 = DefModel->GetIndependantPart();
	double aux2 = DefModel->GetFullCorrelPart();

	DefModel->SetIndependantPart(proportion_independant);
	DefModel->SetFullCorrelPart(proportion_full_correl);

	double Result = itsPricer->ComputePrice(qCMPPRICE);

	//On reset les infos initiales
	DefModel->SetIndependantPart(aux1);
	DefModel->SetFullCorrelPart(aux2);

	itsPricer->ResetPricer();

	return (Result);
}

//---------------------------------------------------------------------------------------//
// Fonctions utilisés dans l'appel à NAG
//---------------------------------------------------------------------------------------//

ARM_Vector* MixCopulaCalibration::FunctionToMinimize(double* x)
{
	// entrer le système d'équations
	ARM_Vector* fvec = new ARM_Vector(1);
	
	ARM_Security* sec1 = itsPtf->GetSecurity(0);
	ARM_Security* sec2 = itsPtf->GetSecurity(1);

	double aux0 = 0. , aux1 = 0.;
	if( (!itsSmiledPrice1) || (!itsSmiledPrice2) )
	{
		aux0 = CptPrice_BC(sec1,itsBaseCorrelation1, itsBCType);
		SetitsSmiledPrice1(aux0) ;

		aux1 = CptPrice_BC(sec2,itsBaseCorrelation2, itsBCType);
		SetitsSmiledPrice2(aux1) ;
	}

	aux0 = CptPrice_MixCopula(sec1, x[0], x[1]);
	aux1 = CptPrice_MixCopula(sec2, x[0], x[1]);

	fvec->InitElt(0, ((itsSmiledPrice1 - aux0) * (itsSmiledPrice1 - aux0))
		           + ((itsSmiledPrice2 - aux1) * (itsSmiledPrice2 - aux1)) );
	
	return fvec;
}

//---------------------------------------------------------------------------------------//
// Fonction Global qui sert dans l'appel à NAG
//---------------------------------------------------------------------------------------//


void NAG_CALL FunctionToMin (Integer m, double* x, double* fvec, double * g, Nag_Comm *comm)
{
	MixCopulaCalibration* MultiDim = (MixCopulaCalibration*)(comm->p);
	ARM_Vector * Res = new ARM_Vector(1);
	
	Res = MultiDim->FunctionToMinimize(x);

	fvec[0] = Res->Elt(0);
		
	if(Res)
		delete Res ;
	Res = NULL ;
	
	if (( fvec[0] >= Value0 - 1000) &&  ( fvec[0] <= Value0 + 1000))
		Count++ ;
	else
	{
		Value0 = fvec[0] ;
		Count = 0 ;
	}

	if (Count == 3)
		comm->flag = -1 ;
}

//---------------------------------------------------------------------------------------//
// Fonction qui encapsule l'appel à Nag
//---------------------------------------------------------------------------------------//

void MixCopulaCalibration::Calibrate()
{
	Count  = 0 ;
	Value0 = 0 ;

	int n = 2;

	double* x = new double[2];
	double *g = new double[2] ;
	double * bl = new double[2] ;
	double * bu = new double[2] ;
	double* fvec = new double[2];
	
	x[0] = itsSeed1;
	x[1] = itsSeed2;

	bl[0] = 0. ;
	bl[1] = 0. ;
	
	bu[0] = 1. ;
	bu[1] = 1. ;

    Nag_Comm comm;
    comm.p = this;

    NagError fail;
    INIT_FAIL(fail); 

	Nag_BoundType bound;
	bound = Nag_Bounds ;

    // Initialise options structure 
	Nag_E04_Opt options ;
	ICM_NAGSILENT(options); 
	options.max_iter = 10 ;
    options.optim_tol = 1.0e-2;
    options.step_max = 1.0e+003;
    options.local_search = FALSE;
	
	e04jbc(n, &FunctionToMin, bound, bl, bu, x, fvec, g, &options, &comm, &fail);

	ICM_NAGFREE(options); 

	SetitsImplied_Proportion_independant (x[0]);
	SetitsImplied_Proportion_FullCorrel  (x[1]);

	delete [] g;
	delete [] bl;
	delete [] bu;
	delete [] x;
	delete [] fvec;
}


/***************************************************************************************************/
// Calibration du modèle de Mix Copula avec Proportion indep uniquement et un Beta unique
/***************************************************************************************************/


/*---------------------------------------------------------------------------------------------------*/
// Cpt Price pour une tranche étant donnée un Beta Unique et uniquement la proportion de copule indep
/*---------------------------------------------------------------------------------------------------*/

double MixCopulaCalibration::CptPrice_SameBeta_Indep(ARM_Security* Sec, double proportion_independant,double Beta)
{
	itsPricer->ResetPricer();
	
	itsPricer->SetSecurity(Sec);
	ICM_ModelMultiCurves* DefModel = (ICM_ModelMultiCurves*)itsPricer->GetModel();
	ARM_Date ExecutionDate = DefModel->GetStartDate();
	
	qSENSITIVITY_TYPE typesensi ;	
	typesensi = ICM_SAMEBETA_TYPE;

	double aux1 = DefModel->GetIndependantPart();
	//DefModel->SetIndependantPart(proportion_independant);
	DefModel->SetFullCorrelPart(proportion_independant);

	if(!itsFirstPrice)
	{
		double price = itsPricer->ComputePrice(qCMPPRICE);
		SetitsFirstPrice(price) ;
	}

	itsPricer->SetInitialPriceFlg(false) ;
	itsPricer->SetPrice(itsFirstPrice) ;

	double Result = itsFirstPrice + itsPricer->Hedge(typesensi,"NONE","NONE",Beta);
	
	//DefModel->SetIndependantPart(aux1);
	DefModel->SetFullCorrelPart(aux1); 

	itsPricer->ResetPricer();

	return Result;
}
//---------------------------------------------------------------------------------------//
// Fonctions utilisés dans l'appel à NAG
//---------------------------------------------------------------------------------------//

ARM_Vector* MixCopulaCalibration::FunctionToMinimize3(double* x)
{
	ARM_Vector* fvec = new ARM_Vector(1);
	
	ARM_Security* sec1 = itsPtf->GetSecurity(0);
	ARM_Security* sec2 = itsPtf->GetSecurity(1);

	double aux0 = 0. , aux1 = 0.;
	if( (!itsSmiledPrice1) || (!itsSmiledPrice2) )
	{
		aux0 = CptPrice_BC(sec1,itsBaseCorrelation1, itsBCType);
		SetitsSmiledPrice1(aux0) ;

		aux1 = CptPrice_BC(sec2,itsBaseCorrelation2, itsBCType);
		SetitsSmiledPrice2(aux1) ;
	}

	aux0 = CptPrice_SameBeta_Indep(sec1, x[0], x[1]);
	aux1 = CptPrice_SameBeta_Indep(sec2, x[0], x[1]);

	fvec->InitElt(0, ((itsSmiledPrice1 - aux0) * (itsSmiledPrice1 - aux0))
		           + ((itsSmiledPrice2 - aux1) * (itsSmiledPrice2 - aux1)) );
	
	return fvec;
}

//---------------------------------------------------------------------------------------//
// Fonction Global qui sert dans l'appel à NAG
//---------------------------------------------------------------------------------------//


void NAG_CALL FunctionToMin3 (Integer m, double* x, double* fvec, double * g, Nag_Comm *comm)
{
	MixCopulaCalibration* MultiDim = (MixCopulaCalibration*)(comm->p);
	ARM_Vector * Res = new ARM_Vector(1);
	
	Res = MultiDim->FunctionToMinimize3(x);

	fvec[0] = Res->Elt(0);
		
	if(Res)
		delete Res ;
	Res = NULL ;
	
	if (( fvec[0] >= Value0 - 100000) &&  ( fvec[0] <= Value0 + 100000))
		Count++ ;
	else
	{
		Value0 = fvec[0] ;
		Count = 0 ;
	}

	if (Count == 3)
	{
		comm->flag = -1 ;
		V0 = x[0] ;
		V1 = x[1] ;
	}		
}

//---------------------------------------------------------------------------------------//
// Fonction qui encapsule l'appel à Nag
//---------------------------------------------------------------------------------------//

void MixCopulaCalibration::Calibrate3()
{
	Count  = 0 ;
	Value0 = 0 ;

	int n = 2;

	double* x = new double[2];
	double *g = new double[2] ;
	double * bl = new double[2] ;
	double * bu = new double[2] ;
	double* fvec = new double[2];
	
	x[0] = itsSeed1; // Part de Correl indep
	x[1] = itsSeed2 ; // Beta unique

	bl[0] = 0. ;
	bl[1] = -0.999 ;
	
	bu[0] = 1.0 ;
	bu[1] = 0.999 ;

    Nag_Comm comm;
    comm.p = this;

    NagError fail;
	INIT_FAIL(fail); 
   

	Nag_BoundType bound;
	bound = Nag_Bounds ;

    // Initialise options structure 
	Nag_E04_Opt options ;
	ICM_NAGSILENT(options); 
	options.max_iter = 10 ;
    options.optim_tol = 1.0e-2;
    options.step_max = 1.0e+003;
    options.local_search = FALSE;
	
	e04jbc(n, &FunctionToMin3, bound, bl, bu, x, fvec, g, &options, &comm, &fail);

	
	ICM_NAGFREE(options); 
	
	if(Count == 3)
	{
		x[0] = V0 ;
		x[1] = V1 ;
	}

	
	SetitsImplied_Proportion_independant (x[0]);
	SetitsImplied_UniqueBeta  (x[1]);
	
	delete [] g;
	delete [] bl;
	delete [] bu;
	delete [] x;
	delete [] fvec;
}
