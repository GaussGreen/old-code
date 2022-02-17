/*
 *$Log: armfrmmarkovvol.cpp,v $
 *Revision 1.6  2004/02/24 12:28:51  emezzine
 *only to compile (constructor modified)
 *
 *Revision 1.5  2003/08/06 10:46:48  emezzine
 *Modif de constructeur
 *
 *Revision 1.4  2003/07/31 12:49:48  emezzine
 * Modif du corps de constructeur
 *
 *Revision 1.3  2003/07/15 11:53:48  emezzine
 *add delete
 *
 *Revision 1.2  2003/07/11 11:53:29  emezzine
 *Modif SearchIndex()
 *
 *Revision 1.1  2003/06/13 17:04:23  mab
 *Initial revision
 *
 *Revision 1.8  2003/06/13 15:55:12  nbiala
 *Replaced std !!!! by stdVal
 *
 *Revision 1.7  2003/06/10 17:14:58  emezzine
 *Update new version
 *!l
 *
 *Revision 1.5  2003/05/26 09:39:21  emezzine
 *Update all method of calibration
 *
 *Revision 1.4  2003/04/16 08:51:30  jpriaudel
 *Modif pour la compilation
 *
 *Revision 1.3  2003/04/11 07:53:17  emezzine
 *Enlever qui bloquait la compilation
 *
 *Revision 1.2  2003/04/04 15:40:46  emezzine
 *ajout 
 *ajout 
 *ajout qui bloquait la compilation
 *
 */

/*--------------------------------------------------------------------------*/
/*                                                                          */
/* ARM_FRMMARKOVVOL.cpp: implementation of the ARM_FRMMarkovVol class.       */
/*                                                                          */
/*--------------------------------------------------------------------------*/


#include "armfrmmarkovvol.h"
#include "swaption.h"
#include "capfloor.h"

ARM_FRMMarkovVol::ARM_FRMMarkovVol(void)
{
	Init();
}



ARM_FRMMarkovVol::ARM_FRMMarkovVol(const ARM_FRMMarkovVol& FRMVol):	ARM_FRMVol(FRMVol)
{
    Init();
    
    BitwiseCopy(&FRMVol);
}

ARM_FRMMarkovVol::ARM_FRMMarkovVol(double NbFactors,
                                   ARM_Date AsOfDate,
								   ARM_Vector* ACalibrationSchedule,
								   ARM_Vector* KCalibrationSchedule,
								   ARM_ReferenceValue* AInitialise,
								   ARM_ReferenceValue* KInitialise,
								   ARM_Vector* MeanReversionA,
								   ARM_Vector* MeanReversionK,
                                   ARM_Vector* SchedToPrice,
                                   ARM_Vector* Mdec,
								   ARM_Matrix* Bounds,
								   ARM_Vector* PowerData,
								   ARM_Vector* SmoothData,
                                   ARM_IRIndex* Index)

{
	Init();    
    //SetIndex(Index);
    itsAsOfDate = AsOfDate;

    int i ;
    if(ACalibrationSchedule)
    {
        int sizeA = ACalibrationSchedule->GetSize();
        for(i = 0; i < sizeA; i++)        
            ACalibrationSchedule->Elt(i) -=  itsAsOfDate.GetJulian();         
    }
    if(KCalibrationSchedule)
    {
        int sizeK = KCalibrationSchedule->GetSize();
        for(i = 0; i < sizeK; i++)        
            KCalibrationSchedule->Elt(i) -=  itsAsOfDate.GetJulian();       
    }

    for(i = 0; i < SchedToPrice->GetSize(); i++)    
        SchedToPrice->Elt(i) -=  itsAsOfDate.GetJulian();     

    if(AInitialise)
    { 
        for(i = 0; i < KInitialise->GetDiscreteDates()->GetSize(); i++) 
            AInitialise->GetDiscreteDates()->Elt(i) -=  itsAsOfDate.GetJulian(); 

        itsACurve = AInitialise->CptReferenceValues(SchedToPrice);
        itsACurvePrice = (ARM_ReferenceValue*) itsACurve->Clone();
        delete AInitialise;
    }

    
    if(KInitialise)
    {
        for(i = 0; i < KInitialise->GetDiscreteDates()->GetSize(); i++) 
            KInitialise->GetDiscreteDates()->Elt(i) -=  itsAsOfDate.GetJulian(); 

        itsKCurve = KInitialise->CptReferenceValues(SchedToPrice);
        itsKCurvePrice = (ARM_ReferenceValue*) itsKCurve->Clone();
        delete KInitialise;
    }

    delete SchedToPrice;
	
    if(ACalibrationSchedule)
    {
        ARM_ReferenceValue* ACurveToCalibrate = itsACurve->CptReferenceValues(ACalibrationSchedule);
        itsACurveToCalibrate = ACurveToCalibrate;
    }
	
    if(KCalibrationSchedule)
    {
	    ARM_ReferenceValue* KCurveToCalibrate = itsKCurve->CptReferenceValues(KCalibrationSchedule);
        itsKCurveToCalibrate = KCurveToCalibrate;
    }

    itsMainRevA = MeanReversionA;
	itsMainRevK = MeanReversionK;

    ARM_ReferenceValue* InitMdec = new ARM_ReferenceValue();
	InitMdec->SetDiscreteDates(itsKCurve->GetDiscreteDates());
	InitMdec->SetDiscreteValues(Mdec);

    SetMCurve(InitMdec);

    delete InitMdec;

    UpdateCurves();

    if(Bounds)
    {
        itsBounds = (ARM_Matrix *) Bounds->Clone();
        SetPowers(int(PowerData->Elt(0)), int(PowerData->Elt(1))); 
	    SetSmoothParameters(SmoothData->Elt(0),SmoothData->Elt(1));
    }

}



ARM_FRMMarkovVol::~ARM_FRMMarkovVol(void)
{
	if (itsACurvePrice)
		delete itsACurvePrice;

	if (itsKCurvePrice)
		delete itsKCurvePrice;

    if (itsACurve)
		delete itsACurve;

	if (itsKCurve)
		delete itsKCurve;
	
	
	if (itsACurveToCalibrate)
		delete itsACurveToCalibrate;
	
	if (itsKCurveToCalibrate)
		delete itsKCurveToCalibrate;
	
	if (itsMainRevA)
		delete itsMainRevA;
	
	if (itsMainRevK)
		delete itsMainRevK;

	if (itsBounds)
		delete itsBounds;
	
}



//___________________________ Standards Functions ________________________________

void ARM_FRMMarkovVol::Init(void)
{
	SetName(ARM_FRMMARKOVVOL);

	itsACurvePrice = NULL;
	itsKCurvePrice = NULL;

    itsACurve = NULL;
	itsKCurve = NULL;

	itsACurveToCalibrate = NULL;
	itsKCurveToCalibrate = NULL;
	
	itsMainRevA = NULL;
	itsMainRevK = NULL;	
	itsBounds   = NULL;

	itsNbParms = 0;
}

void ARM_FRMMarkovVol::BitwiseCopy(const ARM_Object *FRMMarkovVol)
{
	ARM_FRMMarkovVol* src = (ARM_FRMMarkovVol *) FRMMarkovVol;
	
	if (itsACurvePrice)
		delete itsACurvePrice;
	itsACurvePrice = NULL;
	if (src->itsACurvePrice)
		itsACurvePrice = (ARM_ReferenceValue *) src->itsACurvePrice->Clone();
	
	if (itsKCurvePrice)
		delete itsKCurvePrice;
	itsKCurvePrice = NULL;
	if (src->itsKCurvePrice)
		itsKCurvePrice = (ARM_ReferenceValue *) src->itsKCurvePrice->Clone();

    if (itsACurve)
		delete itsACurve;
	itsACurve = NULL;
	if (src->itsACurve)
		itsACurve = (ARM_ReferenceValue *) src->itsACurve->Clone();
	
	if (itsKCurve)
		delete itsKCurve;
	itsKCurve = NULL;
	if (src->itsKCurve)
		itsKCurve = (ARM_ReferenceValue *) src->itsKCurve->Clone();

	if (itsACurveToCalibrate)
		delete itsACurveToCalibrate;
	itsACurveToCalibrate = NULL;
	if (src->itsACurveToCalibrate)
		itsACurveToCalibrate = (ARM_ReferenceValue *) src->itsACurveToCalibrate->Clone();

	if (itsKCurveToCalibrate)
		delete itsKCurveToCalibrate;
	itsKCurveToCalibrate = NULL;
	if (src->itsKCurveToCalibrate)
		itsKCurveToCalibrate = (ARM_ReferenceValue *) src->itsKCurveToCalibrate->Clone();

	if (itsMainRevA)
		delete itsMainRevA;
	itsMainRevA = NULL;
	if (src->itsMainRevA)
		itsMainRevA = (ARM_Vector *) src->itsMainRevA->Clone();

	if (itsMainRevK)
		delete itsMainRevK;
	itsMainRevK = NULL;
	if (src->itsMainRevK)
		itsMainRevK = (ARM_Vector *) src->itsMainRevK->Clone();

	if (itsBounds)
		delete itsBounds;
	itsBounds = NULL;
	if (src->itsBounds)
		itsBounds = (ARM_Matrix *) src->itsBounds->Clone();


	itsSmoothACurve = src->itsSmoothACurve;
	itsSmoothKCurve = src->itsSmoothKCurve;

	itsAsOfDate = src->itsAsOfDate;
	itsNbParms  = src->itsNbParms;
}



void ARM_FRMMarkovVol::Copy(const ARM_Object* FRMMarkovVol)
{
	ARM_FRMVol::Copy(FRMMarkovVol);

	BitwiseCopy(FRMMarkovVol);
}

ARM_Object* ARM_FRMMarkovVol::Clone(void)
{
	ARM_FRMMarkovVol* theClone = new ARM_FRMMarkovVol();
	
    theClone->Copy(this);

    return(theClone);
}

//__________________________ InPut/OutPut DATA __________________________

void ARM_FRMMarkovVol::SetSmoothParameters(double SmoothACurve, double SmoothKCurve)
{
	itsSmoothACurve = SmoothACurve;
	itsSmoothKCurve = SmoothKCurve;
}

void ARM_FRMMarkovVol::GetSmoothParameters(double& SmoothACurve, double& SmoothKCurve)
{
	SmoothACurve = itsSmoothACurve;
	SmoothKCurve = itsSmoothKCurve;
}


	
void ARM_FRMMarkovVol::SetCurvesPrice(ARM_ReferenceValue* ACurve,ARM_ReferenceValue* KCurve)
{
    if (itsACurvePrice)
        delete itsACurvePrice;
	itsACurvePrice = (ARM_ReferenceValue *) ACurve->Clone();

    if (itsKCurvePrice)
        delete itsKCurvePrice;
	itsKCurvePrice = (ARM_ReferenceValue *) KCurve->Clone();
}

void ARM_FRMMarkovVol::GetCurvesPrice(ARM_ReferenceValue*& ACurve,ARM_ReferenceValue*& KCurve)
{
	ACurve = itsACurvePrice;
	KCurve = itsKCurvePrice;
}

void ARM_FRMMarkovVol::SetACurve(ARM_ReferenceValue* ACurve)
{
    if (itsACurve)
        delete itsACurve;
	itsACurve = (ARM_ReferenceValue *) ACurve->Clone();
}
void ARM_FRMMarkovVol::SetKCurve(ARM_ReferenceValue* KCurve)
{
    if (itsKCurve)
        delete itsKCurve;
	itsKCurve = (ARM_ReferenceValue *) KCurve->Clone();
}


ARM_ReferenceValue* ARM_FRMMarkovVol::GetACurve(void)
{
	return(itsACurve);

}

ARM_ReferenceValue* ARM_FRMMarkovVol::GetKCurve(void)
{
	return (itsKCurve);
}



void ARM_FRMMarkovVol::SetCurvesToCalibrate(ARM_ReferenceValue* ACurve,ARM_ReferenceValue* KCurve)
{

    if (itsACurveToCalibrate)
        delete itsACurveToCalibrate;
	itsACurveToCalibrate = (ARM_ReferenceValue *) ACurve->Clone();

    if (itsKCurveToCalibrate)
        delete itsKCurveToCalibrate;
	itsKCurveToCalibrate = (ARM_ReferenceValue *) KCurve->Clone();
		
}

void ARM_FRMMarkovVol::GetCurvesToCalibrate(ARM_ReferenceValue*& ACurve,ARM_ReferenceValue*& KCurve)
{
    ACurve = itsACurveToCalibrate;
    KCurve = itsKCurveToCalibrate;
}



void ARM_FRMMarkovVol::SetMeanRevParams(ARM_Vector* MainRevA,ARM_Vector* MainRevK)
{
    if (itsMainRevA)
        delete itsMainRevA;
	itsMainRevA = (ARM_Vector *) MainRevA->Clone();

    if (itsMainRevK)
        delete itsMainRevK;
	itsMainRevK = (ARM_Vector *) MainRevK->Clone();
		
}
void ARM_FRMMarkovVol::GetMeanRevParams(ARM_Vector*& MainRevA,ARM_Vector*& MainRevK)
{

	MainRevA = itsMainRevA;
	MainRevK = itsMainRevK;
}


size_t ARM_FRMMarkovVol::GetNbParams(void)
{

	size_t  nbRowA = itsACurveToCalibrate->GetSize();
	
	size_t  nbRowK = itsKCurveToCalibrate->GetSize();

	size_t  nbMRRowA = itsMainRevA->GetSize();
	size_t  nbMRRowK = itsMainRevK->GetSize();

	size_t nbcurve = itsACurveToCalibrate->NbCurves();
	
	if ( nbcurve !=itsKCurveToCalibrate->NbCurves() )
	{
	   throw Exception(__LINE__, __FILE__,ARM_FRMMARKOVVOL,"NbACurve isn't equal NbKCurve");
    }
	
	if(itsBounds->Elt(0,0) == itsBounds->Elt(0,1))
	{
		nbRowA = 0;
	}
	if(itsBounds->Elt(1,0) == itsBounds->Elt(1,1))
	{
		nbRowK = 0;
			
	}
	if(itsBounds->Elt(2,0) == itsBounds->Elt(2,1))
	{
		nbMRRowA = 0;
	}
	if(itsBounds->Elt(3,0) == itsBounds->Elt(3,1))
	{
		nbMRRowK = 0;
	}

	size_t Nb = nbcurve*(nbRowA+nbRowK)+nbMRRowA+nbMRRowK;

	return(Nb);
}

void ARM_FRMMarkovVol::SetNbParms(size_t nbParams)
{
	itsNbParms = nbParams;
}


//__________________ Initialization/ Up Date and Consturction ______________________

void ARM_FRMMarkovVol::SetCurvesToPrice(ARM_Vector* vectCurve)
{
	size_t i,j;
	size_t nbcurve = itsACurveToCalibrate->NbCurves();

	if ( nbcurve != itsKCurveToCalibrate->NbCurves())
	{
	   throw Exception(__LINE__, __FILE__,ARM_FRMMARKOVVOL,
		               "NbACurve is different than NbKCurve");
	}

	// Set des piliers des courbes A
	size_t Idx = 0;
	size_t Start = 0;

   	size_t Size = itsACurveToCalibrate->GetDiscreteDates()->GetSize();

	for (j = 0; j < nbcurve; j++)
	{
		Idx = 0;

		if(itsBounds->Elt(2,0) != itsBounds->Elt(2,1))
		{
			itsMainRevA->Elt(j) = vectCurve->Elt(Start+Idx);
			Idx++;
		}
	}
		Start+=Idx; 

	for (j = 0; j < nbcurve; j++)
	{
		Idx = 0; 

		if(itsBounds->Elt(3,0) != itsBounds->Elt(3,1))
		{
			itsMainRevK->Elt(j) = vectCurve->Elt(Start+Idx); 
			Idx++;
		}
	}
	Start += Idx;	

	for (j = 0; j < nbcurve; j++)
	{
		Idx = 0;
		for (i = 0; i < Size ; i++)
		{

		if(itsBounds->Elt(0,0) != itsBounds->Elt(0,1))
			{

			itsACurveToCalibrate->GetDiscreteValues(j)->Elt(i) = vectCurve->Elt(Start+Idx);
				Idx++;
			}
		}

        Start +=Idx;
	}

    size_t Size1 = itsKCurveToCalibrate->GetDiscreteDates()->GetSize();

	for (j = 0; j < nbcurve;j++)
	{
		Idx = 0;
		for (i = 0; i < Size1 ; i++)
		{
			if(itsBounds->Elt(1,0) != itsBounds->Elt(1,1))
			{
				if (i>1)
				{
					itsKCurveToCalibrate->GetDiscreteValues(j)->Elt(i) = vectCurve->Elt(Start+Idx);
					Idx++;
				}
			}
		}

	Start += Idx;
	}

}

size_t ARM_FRMMarkovVol::GetNbParamstoNag(void)
{

	size_t  nbRowA = itsACurveToCalibrate->GetSize();
	size_t  nbRowK = itsKCurveToCalibrate->GetSize();
	
	size_t nbcurve = itsACurveToCalibrate->NbCurves();
	
	if ( nbcurve !=itsKCurveToCalibrate->NbCurves() )
	{
	   throw Exception(__LINE__, __FILE__,ARM_FRMMARKOVVOL,"NbACurve isn't equal NbKCurve");
    }

	size_t Nb = nbcurve*(nbRowA+nbRowK);

	return(Nb);
}

ARM_Vector* ARM_FRMMarkovVol::SetCurvesToOptimize(ARM_Vector*& LBounds,ARM_Vector*& UBounds, size_t Int)
{
	size_t  i, j, nbcurve;

	nbcurve = itsACurveToCalibrate->NbCurves();
				
	size_t Idx = 0;
	size_t  Start = 0;
	size_t Size = GetNbParams();
	

	// Allocate Y 

	ARM_Vector* Y = new ARM_Vector(Size);

	for (j = 0; j < nbcurve;j++)
	{
		Idx = 0;	
        
        if(itsBounds->Elt(2,0) != itsBounds->Elt(2,1))
		{
			LBounds->Elt(Start+Idx)= itsBounds->Elt(2,0);
			UBounds->Elt(Start+Idx)= itsBounds->Elt(2,1);

			Y->Elt(Start+Idx)=itsMainRevA->Elt(j);
			Idx++;
		}
			
	}
	Start += Idx;

	for (j = 0; j < nbcurve;j++)
	{
		Idx = 0;

		

		
		if(itsBounds->Elt(3,0) != itsBounds->Elt(3,1))
		{
			LBounds->Elt(Start+Idx)= itsBounds->Elt(3,0);
			UBounds->Elt(Start+Idx)= itsBounds->Elt(3,1);

			Y->Elt(Start+Idx) = itsMainRevK->Elt(j);
			Idx++;
		}
	}

	Start += Idx;
	
	size_t Size1 = itsACurveToCalibrate->GetDiscreteDates()->GetSize();

	for(j = 0; j < nbcurve; j++)
	{
		Idx = 0;		
				
		for (i = 0 ; i < Size1; i++)
		{
			if(itsBounds->Elt(0,0) != itsBounds->Elt(0,1))
			{
				LBounds->Elt(Start+Idx)= itsBounds->Elt(0,0);
				UBounds->Elt(Start+Idx)= itsBounds->Elt(0,1);

				Y->Elt(Start+Idx) = itsACurveToCalibrate->GetDiscreteValues(j)->Elt(i);
				Idx++;
			}
		}
	
		Start+=Idx;
	}

	size_t Size2 = itsKCurveToCalibrate->GetDiscreteDates()->GetSize();

	for (j = 0;j < nbcurve;j++)
	{
		Idx = 0;	

		for (i = 0 ; i < Size2; i++)
		{
			if (i>1)
			{
				if(itsBounds->Elt(1,0) != itsBounds->Elt(1,1))
				{	
					LBounds->Elt(Start+Idx)= itsBounds->Elt(1,0);
					UBounds->Elt(Start+Idx)= itsBounds->Elt(1,1);

					Y->Elt(Start+Idx) = itsKCurveToCalibrate->GetDiscreteValues(j)->Elt(i);
					Idx++;
				}
			}
		}

		Start+=Idx;
	}

	return(Y);
}


// __________________ Compute volatility __________________________________

ARM_Vector* ARM_FRMMarkovVol::VolatilityVector(double t0, double T,
       										   ARM_CapFloor* capOrFloor)
{

	if ( capOrFloor == NULL)
	{
	   throw Exception(__LINE__, __FILE__,ARM_FRMMARKOVVOL, "Security expected");
    }

	ARM_Date StartDate = capOrFloor->GetStartDate();

	double Ti = StartDate.GetJulian()-itsAsOfDate.GetJulian();

	size_t NbFactor = GetNbFactors();
	ARM_Vector* VolVect = NULL;

	VolVect = InterpolateCurves(t0, T, Ti);

	return(VolVect);
}


ARM_Vector* ARM_FRMMarkovVol::VolatilityVector(double t0, double T1,
   										       ARM_Vector* Mu,
										       ARM_Swaption* swaption)
{
	size_t i, j;


	if (( Mu->GetSize() == 0 ) || ( swaption == NULL ))
    {
	   throw Exception(__LINE__, __FILE__, ARM_FRMMARKOVVOL, "Security expected");
	}

	ARM_Vector* JulianStartDates = swaption->GetResetDates();
	size_t NbW = JulianStartDates->GetSize();

	if ( NbW != Mu->GetSize() )
	{
	   throw Exception(__LINE__, __FILE__, ARM_FRMMARKOVVOL, "Size expected");
	}

	double Ti;
	size_t NbFactor=GetNbFactors();
	ARM_Vector* Factor = NULL;

	ARM_Vector* VolVect = new ARM_Vector(NbFactor);

	for (i = 0; i < NbW; i++)
	{
		Ti = JulianStartDates->Elt(i)-itsAsOfDate.GetJulian();
	
		Factor = InterpolateCurves(t0, T1, Ti);
		
		for (j = 0; j < NbFactor; j++)
		{
			VolVect->Elt(j) += (Factor->Elt(j)*Mu->Elt(i));
        }

	    if (Factor)
	       delete Factor;
	    Factor = NULL;

	}
	
	return(VolVect);
}

// Actualisation des courbes de pricing en fonction des courbes de calage

void ARM_FRMMarkovVol::UpdateCurves(size_t Idx)
{
	Integratepower2Curve(itsACurvePrice,itsMainRevA);
}

void ARM_FRMMarkovVol::Integratepower2Curve(ARM_ReferenceValue* Curve, ARM_Vector* MeanRev)
{
	// Integration du carre de la courbe 
	double x, Nextx, a, sqr;
	size_t j, i;

	size_t Size = Curve->GetSize();
	size_t NbY  = Curve->NbCurves();
	size_t in = 0;
	
	ARM_Vector* YSum = new ARM_Vector(NbY);
		
	x = 0.0;

    for (i = 0; i < Size; i++)
	{
		Nextx = Curve->GetDiscreteDates()->Elt(i);

		ARM_Vector* NextY = Curve->InterpolateMulti(Nextx);

		for (j = 0; j < NbY; j++)
		{
			a = 2.0*MeanRev->Elt(j);

			if (fabs(a)<EPS)
			   sqr = pow(NextY->Elt(j),2.0)*(Nextx-x)/K_YEAR_LEN;
			else
			   sqr =(pow(NextY->Elt(j),2.0)/a *(exp(a*Nextx/K_YEAR_LEN)-exp(a*x/K_YEAR_LEN)));

			if ( sqr < 0.0)
			{
			   throw Exception(__LINE__, __FILE__, ARM_FRMMARKOVVOL,
				               "The integral of A curve must be increasing ");
			}

			YSum->Elt(j) += sqr;

			itsACurvePrice->SetXYValue(i, Nextx, YSum->Elt(j),j);		
		}
		
		x = Nextx;
        delete NextY;
	}

    delete YSum;
}

size_t ARM_FRMMarkovVol::GetNbFactors(void)
{
    size_t  NbF = itsACurvePrice->NbCurves();


    if ( NbF != itsKCurvePrice->NbCurves() ) 
	{
       throw Exception(__LINE__, __FILE__,ARM_FRMMARKOVVOL,
			            "ASizeCurve must be equal to KSizeCurve");
    }
	
	return(NbF);
}

ARM_Vector* ARM_FRMMarkovVol::InterpolateCurves(double t0,double T,double U)
{
	size_t i;
	size_t Size = GetNbFactors();

    ARM_Vector* A = itsACurvePrice->InterpolateMulti(T);
	ARM_Vector* K = itsKCurvePrice->InterpolateMulti(U);	

	// Allocate 
	ARM_Vector* Factor = new ARM_Vector(Size);

	for (i = 0; i < Size; i++)
	{
		double sqr = A->Elt(i);
		
		if ( sqr < 0.0 )
		{
		   throw Exception(__LINE__, __FILE__,ARM_FRMMARKOVVOL,
			               "The integral of A curve must be increasing ");
        }

		double ex = -itsMainRevK->Elt(i)*U/K_YEAR_LEN;
		double eltk = K->Elt(i);
		double x = sqrt(sqr)*eltk* exp(ex);
		
		Factor->Elt(i) = x; 
	}
	
	return(Factor);
}

ARM_Vector* ARM_FRMMarkovVol::StdDevVect(double s,double t,double T)
{
    
    ARM_Vector* var = VarianceForMC(s,t);
    ARM_Vector* coef = CorrelCoefVect(T);
    int size = var->GetSize();

    if ( size != coef->GetSize())
	{
	   throw Exception(__LINE__, __FILE__,ARM_FRMMARKOVVOL,"Pb in Nb Factors");
    }
    ARM_Vector* Y = new ARM_Vector(size);
    for (int i = 0; i < size; i++)
	{
        double stdev = var->Elt(i)* pow(coef->Elt(i),2.0);

        Y->Elt(i) = sqrt(stdev);
    }

    delete var;
    delete coef;
    
    return(Y);
}

ARM_Vector* ARM_FRMMarkovVol::CorrelCoefVect(double U)
{
    int size = itsKCurve->NbCurves();
    ARM_Vector* Y = new ARM_Vector(size);

    double ex = 1.0;
    for (int i = 0; i < size; i++)
	{       
        double a = itsMainRevK->Elt(i);
        if (U>EPS)        
        {
            ex = exp(- a * U/K_YEAR_LEN);
        }
    
        double y = itsKCurve->Interpolate(U,i);
        ex *= y ;
        Y->Elt(i) = ex;     
    }
    
    return (Y);
}

ARM_Vector* ARM_FRMMarkovVol::VarianceForMC(double s, double t)
{
    if (s >t)
	{
		double tmp = t;
		t = s;
		s = t;
	}
	double y,Tj,Ti,sqrsTi, sqrtTj, sqrs, sqrtt, sqrTi, sqrTj; 
      
	size_t nbfactor = itsACurve->NbCurves();
    size_t size = itsACurve->GetSize();
    ARM_Vector* Xsqrt = new ARM_Vector(nbfactor);

	int i = SearchIdx(itsACurve->GetDiscreteDates(),s,K_FRM_TOL);
	int j = SearchIdx(itsACurve->GetDiscreteDates(),t,K_FRM_TOL);

    for(int k = 0; k < nbfactor; k++)
    {
        double a = 2.0*itsMainRevK->Elt(k);

        if(j == -1)
	    {
            if (t < EPS)
                sqrtt = 0.0;
            else
            {
	            y = itsACurve->GetDiscreteValues(k)->Elt(0);

	            if (fabs(a)<EPS)
	            {
		            sqrtt = pow(y,2.0)*(t)/K_YEAR_LEN;
	            }
	            else
	            {
		            sqrtt = (1.0/a)*(pow(y,2.0)*( exp(a*t/K_YEAR_LEN)- 1.0));
	            }
	        }

	    }
	    else
	    {
		    sqrTj = itsACurvePrice->GetDiscreteValues(k)->Elt(j);
		    Tj    = itsACurve->GetDiscreteDates()->Elt(j);

		    if( fabs(t-Tj)<EPS )
		    {
			    sqrtt = sqrTj;
		    }
		    else
		    {
			    if(j == size -1)
			    {
				    y = itsACurve->GetDiscreteValues(k)->Elt(size-1);
			    }
			    else
			    {
				    y = itsACurve->GetDiscreteValues(k)->Elt(j+1);
			    }

			    if (fabs(a)<EPS)
			    {
				    sqrtTj = pow(y,2.0)*(t-Tj)/K_YEAR_LEN;
			    }
			    else
			    {
				    sqrtTj =(1.0/a)*(pow(y,2.0)*(exp(a*t/K_YEAR_LEN)-exp(a*Tj/K_YEAR_LEN)));
			    }
			    			    
			    sqrtt = sqrTj + sqrtTj;
		    }
        }

	    if(i == -1)
	    {
            if (s < EPS)
                sqrs = 0.0;
            else
            {

	            y = itsACurve->GetDiscreteValues(k)->Elt(0);

	            if (fabs(a)<EPS)
	            {
		            sqrs = pow(y,2.0)*(s)/K_YEAR_LEN;
	            }
	            else
	            {
		            sqrs = (1.0/a)*(pow(y,2.0)*( exp(a*s/K_YEAR_LEN)- 1.0));
	            }
	            
            }
		    
	    }
	    else
	    {
		    sqrTi = itsACurvePrice->GetDiscreteValues(k)->Elt(i);
		    Ti    = itsACurve->GetDiscreteDates()->Elt(i);

		    if( fabs(s-Ti)<EPS )
		    {
			    sqrs = sqrTi;
		    }
		    else
		    {
			    if(i == size -1)
			    {
				    y = itsACurve->GetDiscreteValues(k)->Elt(size-1);
			    }
			    else
			    {
				    y = itsACurve->GetDiscreteValues(k)->Elt(i+1);
			    }

			    if (fabs(a)<EPS)
			    {
				    sqrsTi = pow(y,2.0)*(s-Ti)/K_YEAR_LEN;
			    }
			    else
			    {
				    sqrsTi = (1.0/a)*(pow(y,2.0)*(exp(a*s/K_YEAR_LEN)-exp(a*Ti/K_YEAR_LEN)));
			    }
			    			    
			    sqrs = sqrTi + sqrsTi;
		    }
        }
        
        Xsqrt->Elt(k) = sqrtt - sqrs;
    }
       
	return(Xsqrt);
}

double ARM_FRMMarkovVol::VarianceForTree(double s, double t)
{
    ARM_Vector* X = VarianceForMC(s,t);

    double var = X->Elt(0);

    delete X;

    return (var);
}

double ARM_FRMMarkovVol::CorrelCoef(double T)
{
    ARM_Vector* X = this->CorrelCoefVect(T);
    double cor = X->Elt(0);

    delete X;

    return (cor);
}

double ARM_FRMMarkovVol::TimeCoef(double s)
{
    throw Exception(__LINE__, __FILE__,ARM_FRMMARKOVVOL,"Unimplemented <TimeCoef> method");
}

double ARM_FRMMarkovVol::StdDev(double s,double t,double T)
{
    ARM_Vector* X = StdDevVect(s,t,T);

    double stdVal = X->StdNorm();

    delete X;

    return(stdVal);
}

double ARM_FRMMarkovVol::VarianceToTime(double var)
{
    if(var < EPS)
    {
        throw Exception(__LINE__, __FILE__,ARM_FRMMARKOVVOL,"Var is null");
    }

    double t=0.0;

	// On suppose une interpolation et extrapolation constante...
	// sinon c'est bcp plus complique

	// Recuperation des donnees discretisees
	ARM_Vector *ACurve=itsACurve->GetDiscreteValues();
	ARM_Vector *ASched=itsACurve->GetDiscreteDates();
	ARM_Vector *A2Curve=itsACurvePrice->GetDiscreteValues();

	
    double a = 1.0;
    double deuxMRS = 2.0*itsMainRevA->Elt(0);

	int nbA=A2Curve->GetSize();
	int i=0;
    int idx=0;

	if(var<A2Curve->Elt(0))
	{
        a=ACurve->Elt(0);
        a*=a;

        if(fabs(deuxMRS)<EPS)
		{
			t = var/a;
		}
		else
		{
			t=log(var*deuxMRS/a + 1.0)/deuxMRS;
		}
	}
	else
	{
		for(i=1;i<nbA;i++)
		{
			if(var<A2Curve->Elt(i))
			{
				idx=i-1;
				break;
			}
		}
		if(i==nbA)
        {
            idx=nbA-1;            
            a=ACurve->Elt(idx);
            a*=a;
            
        }
        else
        {           
            a=ACurve->Elt(idx+1);
            a*=a;
            
        }
		double ta=ASched->Elt(idx);

		if(fabs(deuxMRS)<EPS)
		{
			t = (var-A2Curve->Elt(idx))/a + ta/K_YEAR_LEN;
		}
		else
		{
			t=log((var-A2Curve->Elt(idx))*deuxMRS/a + exp(deuxMRS*ta/K_YEAR_LEN))/deuxMRS;
		}
	}

	return (t*K_YEAR_LEN);
}


ARM_Vector* ARM_FRMMarkovVol::CalculateSmoothError(void)
{
	size_t i,j, Size, NbCurve;


	NbCurve = itsACurveToCalibrate->NbCurves();

	size_t SizeA = itsACurveToCalibrate->GetSize();
	size_t SizeK = itsKCurveToCalibrate->GetSize();



	if ( (NbCurve != itsKCurveToCalibrate->NbCurves())  || (NbCurve !=GetNbFactors()) )
	{
	   throw Exception(__LINE__, __FILE__,ARM_FRMMARKOVVOL,
		               "NbACurve size is different than NbKCurve size  or than NbFactors");
	}

	// Lissage par penalisation par rapport a la moyenne glissante

	double Err = 0.0;
	Size = NbCurve* (SizeA + SizeK);
	
	ARM_Vector* Smooth = new ARM_Vector(Size);

	size_t Idx = 0;
	size_t  Start = 0;

	for (j = 0; j < NbCurve; j++)
	{
		Idx = 0;
		ARM_Vector* ErrA = FuncSmoothOneCurve(itsACurveToCalibrate->GetDiscreteValues(j), 
			                             itsSmoothACurve,GetSmoothPower());
		
		for (i = 0; i < ErrA->GetSize(); i++)
		{
			Smooth->Elt(Start +Idx ) = ErrA->Elt(i);
			Idx++;
		}
				
		Start += Idx;

		delete ErrA;
	}

	for (j = 0; j < NbCurve; j++)
	{
		Idx = 0;
		
		ARM_Vector* ErrK = FuncSmoothOneCurve(itsKCurveToCalibrate->GetDiscreteValues(j), 
		                                 itsSmoothKCurve,GetSmoothPower());
			
		for (i = 0; i < ErrK->GetSize(); i++)
		{
			Smooth->Elt(Start +Idx ) = ErrK->Elt(i);
			Idx++;
		}
		Start += Idx;

		delete ErrK;
	}


	
	return(Smooth);
}


/*-------------------------------------------------------------------------------------------*/
/*---- End Of File ----*/
