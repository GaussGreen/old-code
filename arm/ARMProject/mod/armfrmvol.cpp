/*
 * $Log: armfrmvol.cpp,v $
 * Revision 1.24  2004/05/17 13:14:28  emezzine
 * formating
 *
 * Revision 1.23  2004/04/07 11:27:32  emezzine
 *  delete itsLogProbaCurve if we want to set it
 *
 * Revision 1.22  2004/03/30 11:26:39  emezzine
 * Correct bug in Clone
 *
 * Revision 1.21  2004/03/19 15:48:51  emezzine
 * correct bug (interpolation).
 *
 * Revision 1.20  2004/02/26 14:48:23  emezzine
 * Change interpolation method of shift curve.
 *
 * Revision 1.19  2004/02/24 12:31:03  emezzine
 * change the interpolation methods and function's names
 *
 * Revision 1.18  2004/02/17 10:55:52  emezzine
 *  Modif Constructor
 *
 * Revision 1.17  2004/02/09 18:13:01  emezzine
 * added SetBetaCurve().
 *
 * Revision 1.16  2003/11/13 18:17:32  emezzine
 * added exeption
 *
 * Revision 1.15  2003/10/31 08:13:06  emezzine
 * New version to calibrate LogProba Curve.
 *
 * Revision 1.14  2003/10/24 07:13:07  rguillemot
 * Prise en compte du spread de la mixture dans MC
 *
 * Revision 1.13  2003/10/07 13:25:01  emezzine
 * Added a new functions Set/GetInitSigmaCurve.
 *
 * Revision 1.12  2003/10/06 15:51:51  rguillemot
 * Ajout d'un spread dans le modele FRM (Modele Mixture)
 *
 * Revision 1.11  2003/07/11 11:42:50  emezzine
 * Cut SearchIndex()
 *
 * Revision 1.10  2003/06/23 10:03:30  emezzine
 * added clone to set IRIndex
 *
 * Revision 1.9  2003/06/17 14:16:34  emezzine
 * Add view virtual function.
 *
 * Revision 1.8  2003/06/10 17:17:25  emezzine
 * define a new virtual  methods
 *
 * Revision 1.7  2003/05/26 09:44:21  emezzine
 * added  functions
 *
 * Revision 1.6  2003/04/02 15:40:20  emezzine
 * Chgt de ARM_INDEX_TYPE en ARM_IRIndex
 *
 * Revision 1.5  2003/03/28 13:20:27  emezzine
 * SetMeanToPrice
 *
 * Revision 1.4  2003/03/24 16:03:45  mab
 * Suppress #define finite
 *
 * Revision 1.3  2003/03/24 16:01:01  mab
 * #define finite _finite
 *
 */

/*--------------------------------------------------------------------------*/
/*                                                                          */
/* ARM_FRMVOL.cpp: implementation of the ARM_FRMVol class.                  */ 
/*                                                                          */
/*--------------------------------------------------------------------------*/

#include "refvalue.h"
#include "irindex.h"

#include "armfrmvol.h"

ARM_Vector* FuncSmoothOneCurve(const ARM_Vector* Curve,
                               double SmoothCoef,
                               double SmoothPower)                        
{
    double x;
    double SumValue = 0.0;
    size_t i;

    size_t Size = Curve->GetSize();
    
    ARM_Vector* SmoothError = new ARM_Vector(Size);

    double Value;
    double LastValue = 0.0;
    double Sum = LastValue;

    for (i = 0; i < Size; i++)
    {
        Value = Curve->Elt(i);
    
        if (fabs(Sum) < 0.00000001 )
           x = (Value-LastValue)/SmoothCoef;
        else
           x = (Value-LastValue)/(SmoothCoef*Sum/i);

        if (!finite(x))
        {
            throw Exception(__LINE__, __FILE__, ARM_FRMVOL,
                "ComputePrice() : Returned an NaN or INFINITY, SmoothError"); 
        }

        SmoothError->Elt(i) = x;
           
        Sum += Value;

        LastValue = Value;
    }

    return SmoothError;
}

ARM_FRMVol::ARM_FRMVol(ARM_ReferenceValue* MCurve)
{
    Init();
    if (MCurve)    
        itsInitShift = (ARM_ReferenceValue *) MCurve->Clone();
}


ARM_FRMVol::ARM_FRMVol(void)
{
    Init();
    itsInitShift = new ARM_ReferenceValue();
    itsInitShift->SetCalcMethod(K_LININTERPOL_REF);
    itsInitShift->SetExtrapolMeth(K_CONSTANT_EXTRAPOL);
    itsLogProbaCurve = new ARM_ReferenceValue();
    itsLogProbaCurve->SetCalcMethod(K_LININTERPOL_REF);
    itsLogProbaCurve->SetExtrapolMeth(K_CONSTANT_EXTRAPOL);
}



ARM_FRMVol::ARM_FRMVol(const ARM_FRMVol& FRMVol): ARM_Object(FRMVol)                                    
{
    Init();
    BitwiseCopy(&FRMVol);
}

ARM_FRMVol::~ARM_FRMVol(void)
{
    if (itsInitShift)
       delete itsInitShift;

    if (itsLogProbaCurve)
        delete itsLogProbaCurve;
    
    if (ITSCeoff)
        delete ITSCeoff;

	if (itsSpreadCurve)
		delete itsSpreadCurve;
}



//_________________ Standards Functions ____________________________


void ARM_FRMVol::Init(void)
{
    SetName(ARM_FRMVOL);

    itsNbFactors = 1;
    itsInitShift = NULL;
    itsLogProbaCurve = NULL;
    ITSCeoff = NULL;
	itsSpreadCurve = NULL;
    itsSmoothPower = 1;
    itsMarketPower = 1;
    itsIndexType = IDXFIXED;
}



void ARM_FRMVol::BitwiseCopy(const ARM_Object *FRMVol)
{
    ARM_FRMVol* src = (ARM_FRMVol *) FRMVol;

    delete itsInitShift;
    itsInitShift = src->itsInitShift ? (ARM_ReferenceValue *) src->itsInitShift->Clone(): NULL;

    delete itsLogProbaCurve;
    itsLogProbaCurve = src->itsLogProbaCurve ? (ARM_ReferenceValue *) src->itsLogProbaCurve->Clone(): NULL; 
    
    delete ITSCeoff;
    ITSCeoff = src->ITSCeoff ? (ARM_Vector *) src->ITSCeoff->Clone() : NULL;    

    delete itsSpreadCurve;
    itsSpreadCurve = src->itsSpreadCurve ? (ARM_ReferenceValue *) src->itsSpreadCurve->Clone() : NULL;

    itsNbFactors = src->itsNbFactors;            
    itsSmoothPower = src->itsSmoothPower;    
    itsMarketPower = src->itsMarketPower;
    itsIndexType = src->itsIndexType;
}



void ARM_FRMVol::Copy(const ARM_Object* FRMVol)
{
    ARM_Object::Copy(FRMVol);
    BitwiseCopy(FRMVol);
}



ARM_Object* ARM_FRMVol::Clone(void)
{
    ARM_FRMVol* theClone = new ARM_FRMVol();
    theClone->Copy(this);

    return((ARM_FRMVol*) theClone);
}

void ARM_FRMVol::SetPowers(size_t MarketPower, size_t SmoothPower)
{
    itsSmoothPower = SmoothPower;
    itsMarketPower = MarketPower;
}



size_t ARM_FRMVol::GetMarketPower(void)
{
    return(itsMarketPower);
}


size_t ARM_FRMVol::GetSmoothPower(void)
{
     return(itsSmoothPower);
}

void ARM_FRMVol::SetMCurve(ARM_ReferenceValue* MCurve)
{
    delete itsInitShift;
    itsInitShift= MCurve ? (ARM_ReferenceValue*)MCurve->Clone() : NULL;    
}

void ARM_FRMVol::ComputeShiftCurve(void)
{    
    if (itsLogProbaCurve)
    {        
        if (ITSCeoff)
        {
            size_t size = ITSCeoff->GetSize();
            ARM_Vector* Mvector = new ARM_Vector(size);
            for (int i = 0; i < size; i++)
            {
                double beta  = LogProbaParam(itsInitShift->GetDiscreteDates()->Elt(i));
                if(beta <K_FRM_TOL)
                    beta = K_FRM_TOL;
                    
                (*Mvector)[i] = (*ITSCeoff)[i]*(1-beta)/beta;             
            }
	        itsInitShift->SetDiscreteValues(Mvector);
        } 
        else
            throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
					       "ITSCeoff is NULL, please advise!!");
    }
    else
        throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
					       "The beta curve (probability) is null, please advise!!");
}

void ARM_FRMVol::ComputeBetaCurve(void)
{    
    if (itsInitShift)
    {        
        if (ITSCeoff)
        {
            size_t size = ITSCeoff->GetSize();

            ARM_Vector* Betavector = new ARM_Vector(size);
            for (int i = 0; i < size; i++)
            {
                double shift  = (*itsInitShift->GetDiscreteValues())[i];
                
                (*Betavector)[i] = (*ITSCeoff)[i]/((*ITSCeoff)[i]+shift);             
            }
	        itsLogProbaCurve->SetDiscreteValues(Betavector);
            itsLogProbaCurve->SetDiscreteDates(itsInitShift->GetDiscreteDates());
        } 
        else
            throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
					       "ITSCeoff is NULL, please advise!!");
    } 
    else
        throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
					       "The Shift curve is null, please advise!!");
}


ARM_ReferenceValue* ARM_FRMVol::GetMCurve()
{
    return itsInitShift;
}

void ARM_FRMVol::SetLogProbaCurve (ARM_ReferenceValue* LogProbaCurve)
{
    delete itsLogProbaCurve;
    itsLogProbaCurve  = LogProbaCurve ;
}

ARM_ReferenceValue* ARM_FRMVol::GetLogProbaCurve()
{
    return itsLogProbaCurve;        
}
void ARM_FRMVol::SetSpreadCurve(ARM_ReferenceValue* SpreadCurve)
{
    if (itsSpreadCurve)
       delete itsSpreadCurve;
    itsSpreadCurve = NULL;
    if (SpreadCurve)
       itsSpreadCurve = (ARM_ReferenceValue *) SpreadCurve->Clone();
}

ARM_ReferenceValue* ARM_FRMVol::GetSpreadCurve(void)
{
    return itsSpreadCurve;
}

double ARM_FRMVol::Mdec(double t)
{
    return itsInitShift->Interpolate(t);
}

double ARM_FRMVol::LogProbaParam(double date)
{
    return itsLogProbaCurve->Interpolate(date);
}

double ARM_FRMVol::Spread(double t)
{
	double spread = 0.0;

	if (itsSpreadCurve)
	{
		itsSpreadCurve->SetCalcMethod(K_LININTERPOL_REF);
		spread = itsSpreadCurve->Interpolate(t);
	}
	
    return spread;
}

// Fonctions virtuelles pures d×finies ici en fonction de la forme de vol
void ARM_FRMVol::SetMeanRevParams(double MainRev)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
					       "Unimplemented <SetMeanRevParams> method");
}

double ARM_FRMVol::GetMeanRevParams(void)
{    
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
					       "Unimplemented <GetMeanRevParams> method");    
    return(0.0);
}

void ARM_FRMVol::SetInitSigmaCurve(ARM_ReferenceValue* Curve)
{
     throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
                     "Unimplemented <SetCurrentCurve> method");
}

ARM_ReferenceValue* ARM_FRMVol::GetInitSigmaCurve(void)
{
     throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
                     "Unimplemented <GetCurrentCurve> method");
     return NULL;
}

void ARM_FRMVol::SetCurrentCurve(ARM_ReferenceValue* Curve)
{
     throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
                     "Unimplemented <SetCurrentCurve> method");
}

ARM_ReferenceValue* ARM_FRMVol::GetCurrentCurve(void)
{
     throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
                     "Unimplemented <GetCurrentCurve> method");
     return NULL;
}

void ARM_FRMVol::SetEltToPrice(double x,size_t Idx)
{
     throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
                     "Unimplemented <SetEltToPrice> method");
}
        
int ARM_FRMVol::GetNbFactors(void)
{
   return itsNbFactors;
}

void ARM_FRMVol::SetNbFactors(int NbFactors)
{
   itsNbFactors = NbFactors;
}

     // FRMVol :
void ARM_FRMVol::CurvesToOptimizerParam(ARM_Vector& Data)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
				       "Unimplemented <CurvesToOptimizerParam> method");
}

void ARM_FRMVol::UpdateCurvesToCalibrate(ARM_Vector& Try)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
				       "Unimplemented <UpdateCurvesToCalibrate> method");
}


ARM_Vector*  ARM_FRMVol::VolatilityVector(double t, double T, ARM_CapFloor* Data) 
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
				       "Unimplemented <VolatilityVector> method");
    return NULL;
}

ARM_Vector*  ARM_FRMVol::VolatilityVector(double t, double T, ARM_Vector* Mu, ARM_Swaption* Data) 
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
				       "Unimplemented <VolatilityVector> method");
    return NULL;
}

size_t ARM_FRMVol::GetNbParams(void)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
				       "Unimplemented <GetNbParams> method");
    return 0.0;
}

void ARM_FRMVol::SetNbParams(size_t nbParmas)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
				       "Unimplemented <SetNbParams> method");
}

void ARM_FRMVol::SetCurvesToPrice(ARM_Vector* Try) 
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
				       "Unimplemented <SetCurvesToPrice> method");
}

void ARM_FRMVol::SetCurvesPrice(ARM_ReferenceValue* ACurve,ARM_ReferenceValue* KCurve)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
				       "Unimplemented <SetCurvesPrice> method");
}

void ARM_FRMVol::GetCurvesPrice(ARM_ReferenceValue*& ACurve,ARM_ReferenceValue*& KCurve) 
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
				       "Unimplemented <GetCurvesPrice> method");
}

ARM_Vector* ARM_FRMVol::SetCurvesToOptimize(ARM_Vector*& UBounds,ARM_Vector*& LBounds, size_t Idx )
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
				       "Unimplemented <SetCurvesToOptimize> method");
    return NULL;
}

ARM_Vector*  ARM_FRMVol::CalculateSmoothError(void) 
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
				       "Unimplemented <CalculateSmoothError> method");
    return NULL;
}

void ARM_FRMVol::UpdateCurves(size_t Idx)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
				       "Unimplemented <UpdateCurves> method");
}
void ARM_FRMVol::OutPutCurves(void)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
				       "Unimplemented <OutPutCurves> method");
}

// For FRMHWVol
double ARM_FRMVol::VarianceForTree(double s,double t)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
				       "Unimplemented <VarianceForTree> method");
    return 0.0;
}

double ARM_FRMVol::StdDev(double s,double t,double T)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
				       "Unimplemented <StdDev> method");
    return 0.0;
}

double ARM_FRMVol::CorrelCoef(double T)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
				       "Unimplemented <CorrelCoef> method");
    return 0.0;
}

double ARM_FRMVol::TimeCoef(double s)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
				       "Unimplemented <TimeCoef> method");
    return 0.0;
}

double ARM_FRMVol::VarianceToTime(double var)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMHWVOL,
				       "Unimplemented <VarianceToTime> method");
    return 0.0;
}

// For FRMMarkovVol
ARM_Vector* ARM_FRMVol::StdDevVect(double s,double t,double T)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMVOL,
                   "Unimplemented <StdDevVect> method");
    return NULL;
}

ARM_Vector* ARM_FRMVol::CorrelCoefVect(double T)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMVOL,
                   "Unimplemented <CorrelCoefVect> method");
    return NULL;
}

ARM_Vector* ARM_FRMVol::VarianceForMC(double s, double t)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMVOL,
                   "Unimplemented <VarianceForMC> method");
    return NULL;
}

ARM_Matrix* ARM_FRMVol::MatrixVolatility(void)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMVOL,
                   "Unimplemented <MatrixVolatility> method");
    return NULL;
}

ARM_Vector* ARM_FRMVol::InterpolateCurves(double t0,double T,double U)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMVOL,
                   "Unimplemented <InterpolateCurves> method");
     return NULL;
}

void ARM_FRMVol::View(char* id, FILE* fOut)
{
    throw Exception(__LINE__, __FILE__, ARM_FRMVOL,
                    "Unimplemented <View> method");
}


/*----------------------------------------------------------------------------------*/
/*---- End Of File ----*/
