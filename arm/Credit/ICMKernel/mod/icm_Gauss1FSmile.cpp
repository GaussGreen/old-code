#error no longer part of the project

#include "ARMKernel/glob/firsttoinc.h"
#include <nag.h>
#include <nage04.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include "ICMKernel\mod\icm_Gauss1FSmile.h"
#include "ICMKernel\inst\icm_nthtd.h"

static int Count ;
static double Value0 ;


/***************************************************************************************************/
// Calcul du smile de correlation dans le modèle de Gauss à 1 facteur
/***************************************************************************************************/

/*---------------------------------------------------------------------------------------*/
// Construction/Destruction
/*---------------------------------------------------------------------------------------*/

void Gauss1FSmile::Init()
{
	itsPricer= NULL;
	itsFirstPrice = 0.;
	itsSmile = 0. ;
}

Gauss1FSmile::Gauss1FSmile(ICM_Pricer* pricer, 
						   int DataType, 
						   double MktData,
						   double UpfrontPay,
						   double Seed,
						   int SmileType)
{
	Init();

	itsPricer = pricer;
	itsDataType = DataType ;
	itsMktData = MktData ;
	itsUpfrontPay = UpfrontPay;
	itsSeed = Seed;
	itsSmileType = SmileType;
}

Gauss1FSmile::~Gauss1FSmile() {}

void Gauss1FSmile::Copy(const ARM_Object* src)
{
	ARM_Object::Copy(src);
     BitwiseCopy(src);
}


/*---------------------------------------------------------------------------------------*/
// Fonctions dont on cherche la racine
/*---------------------------------------------------------------------------------------*/
ARM_Vector* Gauss1FSmile::FunctionToMinimize2 (double *x)
{
	ICM_Ftd* ftd = (ICM_Ftd*) itsPricer->GetSecurity();
	ICM_Security* security = ftd->GetFeeLeg()->GetCreditInfos();
	ARM_Date ExecutionDate = itsPricer->GetModel()->GetStartDate();
	
	int DataType = itsDataType;  // Mkt Spread or Mkt DefLeg
	double MktData = itsMktData;
	int SmileType = itsSmileType ; // Correlation or Beta
	
	qSENSITIVITY_TYPE typesensi ;

	itsPricer->ResetPricer();

	double Result = 0. ;
	double price = 0. ;

	ARM_Vector* fvec = new ARM_Vector(1);

	if(!itsFirstPrice)
	{
		if (DataType == 0)  //Mkt Spread as Input
		{
			if(itsUpfrontPay >0.)
				price = itsPricer->ComputePrice(qCMPPRICE) - itsUpfrontPay ;
			else
				price = itsPricer->ComputePrice(qCMPPRICE);
		}
		else //Mkt DefLeg as Input
			price = itsPricer->Price(qCMPDEFLEGPV) - MktData;
				
		SetitsFirstPrice(price) ;
	}

	if(SmileType == 0) // Correlation as Input
		typesensi = ICM_SAMECORRELATION_TYPE;
	else // Beta Smile as input
		typesensi = ICMBETA_TYPE;
		
	if (DataType == 0)  //Mkt Spread as Input
	{
		security->SetCreditSpread(MktData/100);
		itsPricer->ResetPricer();

		itsPricer->SetInitialPriceFlg(false) ;
		itsPricer->SetPrice(itsFirstPrice) ;
		
		if(itsUpfrontPay >0.)
			Result = itsFirstPrice + itsPricer->Hedge(typesensi,"NONE","NONE",x[0]) + itsUpfrontPay ;
		else
			Result = itsFirstPrice + itsPricer->Hedge(typesensi,"NONE","NONE",x[0]);
	}
	
	else //Mkt DefLeg as Input
	{
		security->SetCreditSpread(0.);
		itsPricer->ResetPricer();

		itsPricer->SetInitialPriceFlg(false) ;
		itsPricer->SetPrice(itsFirstPrice) ;

		Result = itsFirstPrice + itsPricer->Hedge(typesensi,"NONE","NONE",x[0]) + MktData;
	}

	fvec->Elt(0)= fabs(Result);

	return fvec;

}

//---------------------------------------------------------------------------------------//
// Fonction Global qui sert dans l'appel à NAG
//---------------------------------------------------------------------------------------//

void NAG_CALL FunctionToMin2 (Integer m, double* x, double* fvec, double * g, Nag_Comm *comm)
{
	Gauss1FSmile* MultiDim = (Gauss1FSmile*)(comm->p);
	ARM_Vector * Res = new ARM_Vector(1);
	
	Res = MultiDim->FunctionToMinimize2(x);

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
	{
		comm->flag = -1 ;
		Value0 = x[0] ;
	}
}

//---------------------------------------------------------------------------------------//
// Fonction qui encapsule l'appel à Nag
//---------------------------------------------------------------------------------------//

void Gauss1FSmile::Calibrate2()
{
	Count  = 0 ;
	Value0 = 0 ;

	int n = 1;
	
	double* x = new double[1];
	double* fvec = new double[1];
	double * bl = new double[1] ;
	double * bu = new double[1] ;
	double *g = new double[1] ;

	x[0] = itsSeed;

	if(itsSmileType == 0) // Correlation as Input
	{
		bl[0] = 0. ;
		bu[0] = 0.999 ;
	}
	else // Beta Smile as input
	{
		bl[0] = - 0.999 ;
		bu[0] =   0.999 ;
	}

    Nag_Comm comm;
    comm.p = this;

    NagError fail;
	INIT_FAIL(fail); 

	Nag_BoundType bound;
	bound = Nag_Bounds ;
	
    Nag_E04_Opt options;
	ICM_NAGSILENT(options); 
	options.max_iter = 10 ;
    options.optim_tol = 1.0e-2 ;
    options.step_max = 1.0e+003 ;
    options.local_search = FALSE ;

	e04jbc(n, &FunctionToMin2, bound, bl, bu, x, fvec, g, &options, &comm, &fail);
	
	if(Count == 3)
		x[0] = Value0 ;
	
	SetitsSmile (x[0]);

	delete [] g;
	delete [] bl;
	delete [] bu;
	delete [] x;
	delete [] fvec;

	ICM_NAGFREE(options); 
}