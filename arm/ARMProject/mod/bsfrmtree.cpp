/*
 $Log: bsfrmtree.cpp,v $
 Revision 1.2  2002/07/02 12:57:07  mab
 Yan LI add-on

 */


#include "bsfrmtree.h"

ARM_BSFrmTree::ARM_BSFrmTree(ARM_ZeroCurve *zc, ARM_ZeroCurve *baZc, 
				ARM_VolCurve *vol, ARM_VolCurve *smile, int mfine,
                int autoMode, ARM_Date &endDate,  
            	double decay, double slope, double asymptote, int NbFactor, 
            	ARM_Vector *CorrelatedIndexes, ARM_Vector *indexes,
            	ARM_Matrix *correlations)
				: ARM_Model(zc),
				ARM_Frm(zc, NULL, NbFactor, false, AUTO_MODE::ValAutoShapeMode(autoMode), decay, slope, asymptote),
				ARM_FrmTree1(zc, vol, smile, mfine, autoMode, endDate, 
				decay, slope, asymptote, NbFactor,
				CorrelatedIndexes, indexes, correlations)
{
    Init();

    itsUnAdjustedCurve = zc;
    itsAdjustedCurve   = baZc;
	itsAdjWorkingCurve   = (ARM_ZeroLInterpol*) baZc->Clone();
}

ARM_BSFrmTree::ARM_BSFrmTree(ARM_ZeroCurve *zc, ARM_ZeroCurve *baZc, 
				ARM_VolCurve *vol, ARM_VolCurve *smile, 
				ARM_VolCurve *irgvol, ARM_VolCurve *irgsmile, int mfine,
                int autoMode, ARM_Date &endDate,  
            	double decay, double slope, double asymptote, int NbFactor, 
            	ARM_Vector *CorrelatedIndexes, ARM_Vector *indexes,
            	ARM_Matrix *correlations)
				: ARM_Model(zc),
				ARM_Frm(zc, NULL, NbFactor, false, AUTO_MODE::ValAutoShapeMode(autoMode), decay, slope, asymptote),
				ARM_FrmTree1(zc, vol, smile, irgvol, irgsmile, mfine, autoMode, endDate, 
				decay, slope, asymptote, NbFactor,
				CorrelatedIndexes, indexes, correlations)
{
    Init();

    itsUnAdjustedCurve = zc;
    itsAdjustedCurve   = baZc;
	itsAdjWorkingCurve   = (ARM_ZeroLInterpol*) baZc->Clone();
}



// Construction par modele anlytique    
ARM_BSFrmTree::ARM_BSFrmTree(ARM_FrmAna* AnaModel, ARM_ZeroCurve* baZc, ARM_Date& endDate,
            int fineSteps,
            ARM_Vector *CorrelatedIndexes, 
            ARM_Vector *indexes, ARM_Matrix *correlations)
            : ARM_FrmTree1(AnaModel, endDate, fineSteps, CorrelatedIndexes, 
                         indexes, correlations)
{
    Init();

    itsUnAdjustedCurve = AnaModel->GetZeroCurve();
    itsAdjustedCurve   = baZc;
	itsAdjWorkingCurve   = (ARM_ZeroLInterpol*) baZc->Clone();
}


/*
// Construction pour autocalibration (1 classe de risque en vol)
ARM_BSFrmTree::ARM_BSFrmTree(ARM_ZeroCurve *zc, ARM_ZeroCurve* baZc, ARM_VolCurve *vol, ARM_VolCurve *smile, 
                        int productType, ARM_Date &endDate,  
						int Ntraj, int  MCGeneratorType,
            		   double decay, double slope, double asymptote, int NbFactor, 
            		   ARM_Vector *CorrelatedIndexes, ARM_Vector *indexes,
            		   ARM_Matrix *correlations, int NoControl, long seed)
  : ARM_Model(zc),
	ARM_MCModel(endDate, Ntraj, 0),
    ARM_Frm(zc, NULL, NbFactor, false, K_ROW, decay, slope, asymptote),        
    ARM_FrmTree1(zc, vol, smile, productType, endDate, Ntraj, MCGeneratorType,
               decay, slope, asymptote, NbFactor, 
               CorrelatedIndexes, indexes, correlations, NoControl, seed)
{
    Init();

    itsUnAdjustedCurve = zc;
    itsAdjustedCurve   = baZc;
	itsAdjWorkingCurve   = (ARM_ZeroLInterpol*) zc->Clone();
}
 
// Construction pour autocalibration (2 classes de risque en vol)    
ARM_BSFrmTree::ARM_BSFrmTree(ARM_ZeroCurve *zc, ARM_ZeroCurve* baZc, ARM_VolCurve *swoptVol, ARM_VolCurve *swoptSmile, 
           ARM_VolCurve *irgVol, ARM_VolCurve *irgSmile, 
           int productType, ARM_Date &endDate,  
           int Ntraj, int  MCGeneratorType,
           double decay, double slope, double asymptote, int NbFactor, 
           ARM_Vector *CorrelatedIndexes, ARM_Vector *indexes,
           ARM_Matrix *correlations, int NoControl, long seed)
  : ARM_Model(zc),
	ARM_MCModel(endDate, Ntraj, 0),
    ARM_Frm(zc, NULL, NbFactor, false, K_ROW, decay, slope, asymptote),        
    ARM_FrmTree1(zc, swoptVol, swoptSmile, irgVol, irgSmile, productType, endDate, Ntraj, MCGeneratorType,
            		    decay, slope, asymptote, NbFactor, 
            		    CorrelatedIndexes, indexes, correlations, NoControl, seed)
{
    Init();

    itsUnAdjustedCurve = zc;
    itsAdjustedCurve   = baZc;
	itsAdjWorkingCurve   = (ARM_ZeroLInterpol*) baZc->Clone();
}

*/

void ARM_BSFrmTree::BeFittedTo(ARM_Security *sec,char *ccyName, int mfine)
{
    ARM_FrmTree1::BeFittedTo(sec, ccyName, mfine);

    cleanit(itsAdjustedFwds);

    itsAdjustedFwds = new ARM_Forwards(this, 0, itsAdjustedCurve);

    ARM_Forwards* fwds = GetFwds();

    int size = itsAdjustedFwds->GetSize();

    cleanit(itsBasisSpreads);

    itsBasisSpreads = new ARM_Vector(size);

/* old code 

    for (int i=0; i<size; i++)
    {
        itsBasisSpreads->Elt(i) = itsAdjustedFwds->Elt(i) - fwds->Elt(i);
    }
*/

	double dti;

    for (int i=0; i<size; i++)
    {
		dti = GetEpochs()->get_dt(i);
        itsBasisSpreads->Elt(i) = log( 1.0 + itsAdjustedFwds->Elt(i) * dti ) - log( 1.0 + fwds->Elt(i) * dti );
    }	
}


double ARM_BSFrmTree::FirstPeriodDiscountFactor(int nPath)
{

	double df;

    int timeStep = GetTimeStep();
/* old code
    df = 1 / (1 + GetEpochs()->get_dt(timeStep) *
              (GetWorkingForwards()->Elt(timeStep) + itsBasisSpreads->Elt(timeStep)));

*/


	 df  = 1.0 / 
		 (1 + GetEpochs()->get_dt(timeStep) * GetWorkingForwards()->Elt(timeStep));

	 if(itsBasisDiscountFlag)
		 df *= exp( -itsBasisSpreads->Elt(timeStep) );
	 return df;

}


double ARM_BSFrmTree::ZeroPriceOnDiscountCurve(double calcDate, double zMaturity, int DomOrFrg) 
 {
	double zPrice = itsAdjWorkingCurve->DiscountPrice(zMaturity - calcDate);

	return zPrice;
}



double ARM_BSFrmTree::ZeroPrice(double calcDate, double zMaturity, int DomOrFrg) 
{
	double zPrice = 0;

	if (DomOrFrg != 0)
	{
		zPrice = GetWorkingCurve()->DiscountPrice(zMaturity - calcDate);
	}
	else
	{
		zPrice = itsAdjWorkingCurve->DiscountPrice(zMaturity - calcDate);
	}

	

	return zPrice;
}

void ARM_BSFrmTree::SetWorkingFor_n_Curve(ARM_Forwards *inFwds)
{
	ARM_FrmTree1::SetWorkingFor_n_Curve(inFwds);
	GetFwds()->ShiftWithArg( (ARM_ZeroLInterpol *)itsAdjWorkingCurve, itsBasisSpreads);
}

void ARM_BSFrmTree::Propagate(
        ARM_Object *OutCurve,
        ARM_Vector *dr1,
        int time1, 
        int time2,
        ARM_Object *InpCurve, 
        bool propFullCurve,
        ARM_Vector *driversTarget,
        ARM_Vector *driversCurrent)
{
    bool killIn = false;
    bool killOut = false;
    ARM_Forwards *InpFwds = NULL;
    ARM_Forwards *OutFwds = NULL;

    if (InpCurve->GetName() == ARM_FORWARDS)
    {
        InpFwds = static_cast<ARM_Forwards *> (InpCurve);
    }
    else
    {
        InpFwds = new ARM_Forwards();
        *InpFwds << (ARM_ZeroCurve*) InpCurve;
        killIn = true;
    }

    if (OutCurve != NULL)
    {
        switch (OutCurve->GetName())
        {
            case ARM_FORWARDS :
                {
                    OutFwds = static_cast<ARM_Forwards *> (OutCurve);
                }
                break;

            case ARM_ZERO_LIN_INTERPOL :
			case ARM_ZERO_FLAT:
                {
                    OutFwds = new ARM_Forwards();
                    killOut = true;
                }
                break;

            default :
            		throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET,
                        "Invalid output type for FRM propagate");
        }

    }
    else 
   		throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET,
             "Output cannot be NULL for FRM propagate");


        
    
       

    BasePropagate(OutFwds, dr1, time1, time2, InpFwds, propFullCurve, driversTarget, driversCurrent);

    if (killIn)
        cleanit(InpFwds);

    if (killOut)
    {
        *OutFwds >> (ARM_ZeroLInterpol*) OutCurve;
		OutFwds->ShiftWithArg( (ARM_ZeroLInterpol *)itsAdjWorkingCurve, itsBasisSpreads);
        cleanit(OutFwds);
    }

}
