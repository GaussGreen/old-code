/*
 * $Log: bmc_frm.cpp,v $
 * Revision 1.5  2002/07/02 10:36:30  mab
 *  Yan LI add-on
 *
 * Revision 1.4  2002/03/01 16:11:19  mab
 * MODIF : YLI Correction
 *
 * Revision 1.3  2001/10/03 16:27:16  smysona
 * *** empty log message ***
 *
 * Revision 1.2  2001/09/25 13:53:22  mab
 * Rajout de Commentaire : RCS
 *
 */

#include "firsttoinc.h"
#include "bmc_frm.h"




ARM_BMCFrm::ARM_BMCFrm(ARM_ZeroCurve* zc, ARM_ZeroCurve* baZc,
                       ARM_FRM_Calibrator* aCalibrator, 
                       ARM_Date& endDate, int Ntraj,  
                       int MCGeneratorType, int NumFact,
                        bool decal, long seed)
           : ARM_SMCFrm(zc, aCalibrator, endDate, Ntraj, 
                        MCGeneratorType, NumFact, decal, seed)
{
    Init();

    itsUnAdjustedCurve = zc;
    itsAdjustedCurve   = baZc;

// Fin YAN 02/2002
    itsAdjWorkingCurve   = (ARM_ZeroLInterpol*) baZc->Clone();
}



// Construction par modele anlytique    
ARM_BMCFrm::ARM_BMCFrm(ARM_FrmAna* AnaModel, ARM_ZeroCurve* baZc,
                       ARM_Date& endDate,
                       int Ntraj, int MCGeneratorType,
                       ARM_Vector* CorrelatedIndexes, 
                       ARM_Vector* indexes, ARM_Matrix* correlations,
                       long seed)
            : ARM_SMCFrm(AnaModel, endDate, Ntraj, 
                         MCGeneratorType, CorrelatedIndexes, 
                         indexes, correlations, seed)
{
    Init();

    itsUnAdjustedCurve = AnaModel->GetZeroCurve();
    itsAdjustedCurve   = baZc;

// Fin YAN 02/2002
    itsAdjWorkingCurve   = (ARM_ZeroLInterpol*) baZc->Clone();
}



// Construction pour autocalibration (1 classe de risque en vol)

ARM_BMCFrm::ARM_BMCFrm(ARM_ZeroCurve* zc, ARM_ZeroCurve* baZc,
                       ARM_VolCurve* vol, ARM_VolCurve* smile, 
                       int productType, ARM_Date& endDate,  
                       int Ntraj, int  MCGeneratorType,
                       double decay, double slope,
                       double asymptote, int NbFactor, 
                       ARM_Vector* CorrelatedIndexes,
                       ARM_Vector* indexes,
                       ARM_Matrix* correlations, int NoControl, long seed)
                   : ARM_Model(zc),
                     ARM_MCModel(endDate, Ntraj, 0),
                     ARM_Frm(zc, NULL, NbFactor, false, K_ROW,
                             decay, slope, asymptote),        
                     ARM_SMCFrm(zc, vol, smile, productType,
                                endDate, Ntraj, MCGeneratorType,
                                decay, slope, asymptote, NbFactor, 
                                CorrelatedIndexes, indexes,
                                correlations, NoControl, seed)
{
    Init();

    itsUnAdjustedCurve = zc;
    itsAdjustedCurve   = baZc;

// Fin YAN 02/2002
    itsAdjWorkingCurve   = (ARM_ZeroLInterpol*) zc->Clone();
}
 
// Construction pour autocalibration (2 classes de risque en vol)    
ARM_BMCFrm::ARM_BMCFrm(ARM_ZeroCurve* zc, ARM_ZeroCurve* baZc,
                       ARM_VolCurve* swoptVol, ARM_VolCurve* swoptSmile, 
                       ARM_VolCurve* irgVol, ARM_VolCurve* irgSmile, 
                       int productType, ARM_Date& endDate,  
                       int Ntraj, int  MCGeneratorType,
                       double decay, double slope, double asymptote,
                       int NbFactor, 
                       ARM_Vector* CorrelatedIndexes, ARM_Vector* indexes,
                       ARM_Matrix* correlations, int NoControl, long seed)
                   : ARM_Model(zc),
                     ARM_MCModel(endDate, Ntraj, 0),
                     ARM_Frm(zc, NULL, NbFactor, false,
                             K_ROW, decay, slope, asymptote),        
                             ARM_SMCFrm(zc, swoptVol, swoptSmile,
                             irgVol, irgSmile, productType, endDate,
                             Ntraj, MCGeneratorType,
                             decay, slope, asymptote, NbFactor, 
                             CorrelatedIndexes, indexes,
                             correlations, NoControl, seed)
{
    Init();

    itsUnAdjustedCurve = zc;
    itsAdjustedCurve   = baZc;

// Fin YAN 02/2002
    itsAdjWorkingCurve   = (ARM_ZeroLInterpol*) baZc->Clone();
}



void ARM_BMCFrm::BeFittedTo(ARM_Security *sec)
{
    ARM_SMCFrm::BeFittedTo(sec);

    cleanit(itsAdjustedFwds);

    itsAdjustedFwds = new ARM_Forwards(this, 0, itsAdjustedCurve);

    ARM_Forwards* fwds = GetFwds();

    int size = itsAdjustedFwds->GetSize();

    cleanit(itsBasisSpreads);

    itsBasisSpreads = new ARM_Vector(size);

/* YAN 02/2002 old code 

    for (int i=0; i<size; i++)
    {
        itsBasisSpreads->Elt(i) = itsAdjustedFwds->Elt(i) - fwds->Elt(i);
    }
*/

// YAN 02/2002
    double dti;

    for (int i = 0; i < size; i++)
    {
        dti = GetEpochs()->get_dt(i);

        itsBasisSpreads->Elt(i) = log(1.0+itsAdjustedFwds->Elt(i)*dti )
                                     -log(1.0+fwds->Elt(i)*dti);
    }    
// Fin YAN

/* OLD MEC  YAN  02/2002
    cleanit(itsSpreadCurve);

    ARM_Forwards basis(this,0);

    basis.SetFwds(itsBasisSpreads);

    basis >> (itsSpreadCurve);
*/
}


double ARM_BMCFrm::FirstPeriodDiscountFactor(int nPath)
{

    double df;

    int timeStep = GetTimeStep();
// YAN 02/2002
/* old code
    df = 1 / (1 + GetEpochs()->get_dt(timeStep) *
              (GetWorkingForwards()->Elt(timeStep) 
               +itsBasisSpreads->Elt(timeStep)));

*/


     df  = 1.0/ 
            (1+GetEpochs()->get_dt(timeStep)
             *GetWorkingForwards()->Elt(timeStep));

     if (itsBasisDiscountFlag)
        df *= exp(-itsBasisSpreads->Elt(timeStep) );

// Fin YAN
    return df;
}



double ARM_BMCFrm::ZeroPriceOnDiscountCurve(double calcDate,
                                            double zMaturity,
                                            int DomOrFrg) 
{
// YAN 02/2002

    double zPrice = itsAdjWorkingCurve->DiscountPrice(zMaturity - calcDate);
    //double zPrice = itsAdjustedCurve->DiscountPrice(zMaturity - calcDate);

    return zPrice;
}



double ARM_BMCFrm::ZeroPrice(double calcDate, double zMaturity, int DomOrFrg) 
{
// YAN 02/2002
    double zPrice = 0;

    if ( DomOrFrg != 0 )
    {
       zPrice = GetWorkingCurve()->DiscountPrice(zMaturity - calcDate);
    }
    else
    {
       zPrice = itsAdjWorkingCurve->DiscountPrice(zMaturity - calcDate);
    }

// Fin YAN

    return zPrice;
}



// YAN 02/2002
void ARM_BMCFrm::PropagateCurve(int nPath, int nTimeStep)
{
    ARM_SMCFrm::PropagateCurve(nPath, nTimeStep);

    GetFwds()->ShiftWithArg((ARM_ZeroLInterpol *) itsAdjWorkingCurve,
                            itsBasisSpreads);
}
// Fin YAN



// YAN 02/2002
void ARM_BMCFrm::InitTrajectory(int nPath)
{

    ARM_SMCFrm::InitTrajectory(nPath);

    cleanit(itsAdjWorkingCurve);

    itsAdjWorkingCurve = (ARM_ZeroLInterpol*) itsAdjustedCurve->Clone();

}
// Fin YAN



/*-----------------------------------------------------------------*/
/*---- End Of File ----*/
