
#include "ICMKernel\mod\icm_defcurvemodel.h"


//	----------------------------------------------------------------------------------------------
ICM_DefaultCurveModel::ICM_DefaultCurveModel(const ICM_DefaultCurve* DefProbCurve,
						  ARM_ZeroCurve* ZeroCurve,
						  ARM_VolCurve*	VolCurve ,
						  bool FlgClone )
{
	Init();
	Set(DefProbCurve,ZeroCurve,VolCurve,FlgClone);
	ARM_Model::SetZeroCurve(itsZeroCurve);
}
//	----------------------------------------------------------------------------------------------
ICM_DefaultCurveModel::ICM_DefaultCurveModel(const ICM_DefaultCurve* DefProbCurve,
						  ARM_ZeroCurve* ZeroCurve,
						  const ARM_VolCurve* VolCurve  ) 
{
	Init();
	Set(const_cast<ICM_DefaultCurve*>(DefProbCurve),ZeroCurve,const_cast<ARM_VolCurve*>(VolCurve),true);
	ARM_Model::SetZeroCurve(itsZeroCurve);
}
//	----------------------------------------------------------------------------------------------
ICM_DefaultCurveModel::ICM_DefaultCurveModel(const ICM_DefaultCurveModel&ref) : ARM_Model(ref)
{
	Init(); 

	// JLA. my understanding is that a copy will alwaus result in cloning arguments. 
	if (ref.itsDefaultCurve)
		itsDefaultCurve = (ICM_DefaultCurve*) ref.itsDefaultCurve->Clone();
	
	if (ref.itsZeroCurve)
		itsZeroCurve = (ARM_ZeroCurve*) ref.itsZeroCurve->Clone();
	
	if (ref.itsVolCurve)
		itsVolCurve = (ARM_VolCurve*) ref.itsVolCurve->Clone();
	
}
//	----------------------------------------------------------------------------------------------
// virtual 
ICM_DefaultCurveModel::~ICM_DefaultCurveModel(void)
{
	if (itsFlgClone)
	{
		if (itsZeroCurve) delete itsZeroCurve;itsZeroCurve=NULL;
		if (itsDefaultCurve) delete itsDefaultCurve;itsDefaultCurve=NULL;
		if (itsVolCurve) delete itsVolCurve;itsVolCurve=NULL;
	}
}

void ICM_DefaultCurveModel::Init(void)
{
	SetName(ICM_DEFAULT_CURVE_MODEL);

	itsDefaultCurve = NULL;
	itsZeroCurve = NULL;
	itsVolCurve = NULL;
	itsFlgClone = true;
}

/** JLA 
void ICM_DefaultCurveModel::BitwiseCopy(const ARM_Object* src)
{
	ICM_DefaultCurveModel *srcMod = (ICM_DefaultCurveModel *) src;

	if (itsFlgClone)
	{
		if (srcMod->itsDefaultCurve)
		{
			if (srcMod->itsDefaultCurve)
				itsDefaultCurve = (ICM_DefaultCurve*) srcMod->itsDefaultCurve->Clone();
		}

		if (srcMod->itsZeroCurve)
		{		
			if (srcMod->itsZeroCurve)
				itsZeroCurve = (ARM_ZeroCurve*) srcMod->itsZeroCurve->Clone();
		}

		if (srcMod->itsVolCurve)
		{
			if (srcMod->itsVolCurve)
				itsVolCurve = (ARM_VolCurve*) srcMod->itsVolCurve->Clone();
		}
	}
	else
	{
		//JLA: my understanding is we never get there. 
		itsDefaultCurve = srcMod->itsDefaultCurve;
		itsZeroCurve = srcMod->itsZeroCurve;
		itsVolCurve = srcMod->itsVolCurve;
	}
}
**/ 

// virtual 
/** 
void ICM_DefaultCurveModel::Copy(const ARM_Object* src)
{
    ARM_Model::Copy(src);

    BitwiseCopy(src);
} **/ 

// virtual 
ARM_Object* ICM_DefaultCurveModel::Clone(void)
{
/**     ICM_DefaultCurveModel* theClone = new ICM_DefaultCurveModel();

    theClone->Copy(this);

    return(theClone); **/ 
	return new ICM_DefaultCurveModel(*this); 
}

void ICM_DefaultCurveModel::SetDefaultCurve(const ICM_DefaultCurve* notcrv) 
{ 
	if (itsFlgClone)
	{
		if (itsDefaultCurve) delete itsDefaultCurve;
		itsDefaultCurve = (ICM_DefaultCurve*) notcrv->Clone();
	}
	else
		itsDefaultCurve = (ICM_DefaultCurve*) notcrv;
}
// virtual 
void ICM_DefaultCurveModel::SetZeroCurve(ARM_ZeroCurve* crv) 
{ 
	if (itsFlgClone)
	{
		if (itsZeroCurve) delete itsZeroCurve;	
		itsZeroCurve = (ARM_ZeroCurve*) crv->Clone();
	}
	else
		itsZeroCurve = (ARM_ZeroCurve*) crv;
}
void ICM_DefaultCurveModel::SetVolCurve(ARM_VolCurve* crv) 
{ 
	if (itsFlgClone)
	{
		if (itsVolCurve) delete itsVolCurve;	
		itsVolCurve = (ARM_VolCurve*) crv->Clone();
	}
	else 
		itsVolCurve = (ARM_VolCurve*) crv;
}
void ICM_DefaultCurveModel::Set(const ICM_DefaultCurve* DefaultCurve,ARM_ZeroCurve* ZeroCurve, ARM_VolCurve* VolCurve,bool FlgClone)
{
	itsFlgClone = FlgClone;
	SetStartDate(ZeroCurve->GetAsOfDate());

	if (itsFlgClone)
	{
		SetDefaultCurve(DefaultCurve);
		SetZeroCurve(ZeroCurve);
		if (VolCurve) SetVolCurve(VolCurve);
	}
	else
	{
		itsDefaultCurve = DefaultCurve;
		itsZeroCurve = ZeroCurve;
		if (VolCurve) itsVolCurve = VolCurve;
	}
}


void ICM_DefaultCurveModel::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[40];
//    char strDate[20];
 
    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

   fprintf(fOut, "\n\n------------------------------------------------------------------\n");
   fprintf(fOut, "---------------- Default Model -------------------------------------------\n");
   fprintf(fOut, "------------------------------------------------------------------\n\n\n");

   itsDefaultCurve->View(id,fOut);

   fprintf(fOut, "\n\n---------------- Zero Curve -------------------------------------------\n");

   itsZeroCurve->View(id,fOut);

   fprintf(fOut, "\n\n---------------- Volatility Curve -------------------------------------------\n");

   if (itsVolCurve) itsVolCurve->View(id,fOut);

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}


ICM_DefaultCurveModel* ICM_DefaultCurveModel::GenerateShiftModel(qSENSITIVITY_TYPE typesensi,
																 const std::string& plot, 
																 double epsvalue)
{
    double sensitivity =0.;
	int NoCurve = -1;
	int i=0,h=0;

	double result = 0.0;
	bool Parallelshift = false; 

	vector<string> Term;
	char TermZC[ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS];
	memset(TermZC,'\0',sizeof(char)*ARM_NB_TERMS*ARM_NB_MAX_CHAR_TERMS);

	// ARM_Vector* Epsilon = NULL;

	double EPSL = 0.01;		//Bump sur les taux fixé à 1bp	
	double EPSL_DEF = 0.01; //Bump sur le spread fixé à 1 bp.
	double EPSL_REC = 0.1;  //Bump fixe sur le recovery de 0.1
	double EPSL_SPREL = 0.1 ;//Bump relatif des spreads de 0.1
	double EPSL_DTR = CREDIT_DEFAULT_VALUE;//Bump relatif des spreads de 0.1

	if (epsvalue != CREDIT_DEFAULT_VALUE)
		EPSL_DTR = EPSL_REC = EPSL_DEF = EPSL = EPSL_SPREL = epsvalue;

	ICM_DefaultCurve* InitDefCurve = NULL;
	ICM_DefaultCurve* ModifDefCurve = NULL;
	ARM_ZeroCurve* InitShortCurve = NULL;
	ARM_ZeroCurve* ModifShortCurve = NULL;
	ARM_ReferenceValue* Recovery = NULL;
	ARM_Vector* DiscreteValues = NULL,DiscreteDates = NULL;

	ARM_Vector Epsilon (1,epsvalue);

	// if (!strcmp(plot,"NONE"))  //Detection d'un parallel Shift
	if (plot=="NONE") 
	{
		Term.push_back("ALL");
		Parallelshift = true;
	}
	else
	{
		Term.push_back(plot);
		strcpy(TermZC[0],plot.c_str());
		switch (typesensi)
		{
			case ICMRECOVERY_TYPE :
			{
				Epsilon.InitElt(0,EPSL_REC);
			}
			break;
			case ICMIRCURVE_TYPE :
			case ICM_GREEK_RHO_TYPE :
			{
				Epsilon.InitElt(0,EPSL);
			}
			break;
			default:
			case ICMSPREAD_TYPE :
			{
				Epsilon.InitElt(0,EPSL_DEF);
			}
			break;
			case ICMSPRELSHIFT_TYPE :
			{
				Epsilon.InitElt(0,EPSL_SPREL);
			}
			break;
			case ICM_DTR_TYPE :
			{
				Epsilon.InitElt(0,EPSL_DTR);
			}
			break;
		}
	}

	InitDefCurve = (ICM_DefaultCurve*) GetDefaultCurve()->Clone();
	InitShortCurve = (ARM_ZeroCurve*) GetZeroCurve()->Clone();

    // 
    // 
	{// ---------------------------------------------------------------------------
		// Bump Recovery / IrCurve / Spread 
		// ---------------------------------------------------------------------------
		switch (typesensi)
		{
			case ICM_DTR_TYPE :
			{

					ModifDefCurve = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(Term,Epsilon,ICM_DTR_TYPE);
					ModifDefCurve->SetLabel(InitDefCurve->GetLabel());
					ModifShortCurve = (ARM_ZeroCurve*) InitShortCurve->Clone();

			}
			break;
			case ICMRECOVERY_TYPE :
			{


					// if (Parallelshift)
					// 	ModifDefCurve = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(EPSL_REC,ICMRECOVERY_TYPE);
					// else
					ModifDefCurve = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(Term,Epsilon,ICMRECOVERY_TYPE);

					// JLA: useless. ModifDefCurve->SetLabel(InitDefCurve->GetLabel());
					ModifShortCurve = (ARM_ZeroCurve*) InitShortCurve->Clone();
			}
			break;
			case ICMIRCURVE_TYPE :
			case ICM_GREEK_RHO_TYPE :
			{

					if (Parallelshift)
					{	
						ModifShortCurve = (ARM_ZeroCurve*)(InitShortCurve->Clone());
						ModifShortCurve->ParallelShift(EPSL);
					}
					else
						ModifShortCurve = (ARM_ZeroCurve*) InitShortCurve->GenerateShiftCurve(TermZC,&Epsilon);


					// Impact sur la DefProbCurve  ********************************
					ModifDefCurve = (ICM_DefaultCurve*) InitDefCurve->Clone();
					ModifDefCurve->SetZeroCurve((ARM_ZeroCurve*)ModifShortCurve->Clone());
			
			}
			break;
			case ICMIRCURVE_WITHOUT_DEFCURVE_TYPE :
			{
				if (Parallelshift)
				{
					ModifShortCurve = (ARM_ZeroCurve*) InitShortCurve->Clone();
					ModifShortCurve->ParallelShift(EPSL);
				}
				else
					ModifShortCurve = (ARM_ZeroCurve*) InitShortCurve->GenerateShiftCurve(TermZC,&Epsilon);

				ModifDefCurve->SetZeroCurve((ARM_ZeroCurve*)ModifShortCurve->Clone());
			}
			break;
			case ICMSPRELSHIFT_TYPE :
			{
					// if (Parallelshift)
					// 	ModifDefCurve = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(EPSL_DEF,ICMSPRELSHIFT_TYPE);
					// else
					// 	ModifDefCurve = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(Term,Epsilon,ICMSPRELSHIFT_TYPE);
					ModifDefCurve = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(Term,Epsilon,ICMSPRELSHIFT_TYPE);
					// JLA. Useless. ModifDefCurve->SetLabel(InitDefCurve->GetLabel());
					ModifShortCurve = (ARM_ZeroCurve*) InitShortCurve->Clone();
			}
			break;
			case ICMSPREAD_TYPE :
			{

					// if (Parallelshift)
					// 	ModifDefCurve = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(EPSL_DEF,ICMSPREAD_TYPE);
					// else
					ModifDefCurve = (ICM_DefaultCurve*) InitDefCurve->GenerateShiftCurve(Term,Epsilon,ICMSPREAD_TYPE);

					//JLA. Useless . ModifDefCurve->SetLabel(InitDefCurve->GetLabel());
					ModifShortCurve = (ARM_ZeroCurve*) InitShortCurve->Clone();	
			}
			break;
			default :
			;
		}


	// ---------------------------------------------------------------------------
	// Reconstruction Model
	// ---------------------------------------------------------------------------

	ICM_DefaultCurveModel* model = new ICM_DefaultCurveModel(ModifDefCurve,ModifShortCurve,GetVolCurve());

	// if (Epsilon) delete Epsilon;
	if (ModifDefCurve) delete ModifDefCurve;
	if (ModifShortCurve) delete ModifShortCurve;
	if (InitDefCurve) delete InitDefCurve;
	if (InitShortCurve) delete InitShortCurve;

	return (model);
	}
}
