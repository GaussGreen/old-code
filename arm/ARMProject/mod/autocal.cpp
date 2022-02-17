/*
 * $Log: autocal.cpp,v $
 * Revision 1.16  2004/02/18 15:43:10  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.15  2003/06/30 09:07:28  ebenhamou
 * remove unused var
 *
 * Revision 1.14  2002/07/02 10:35:00  mab
 * Yan LI add-on
 *
 * Revision 1.13  2002/03/01 16:08:28  mab
 * MODIF : YLI Corrections
 *
 * Revision 1.12  2001/10/03 16:25:55  smysona
 * *** empty log message ***
 *
 * Revision 1.11  2001/08/03 10:03:56  smysona
 * Generation de la vol avant les output
 * dans EndCalibration
 *
 * Revision 1.10  2001/07/30 08:52:51  smysona
 * Modifs des proto et pour les cas cap
 *
 * Revision 1.9  2001/05/23 17:40:25  smysona
 * modifs des reports
 *
 * Revision 1.8  2001/04/27 09:33:10  smysona
 *  Ajout des dates des instruments de calage
 * dans le modele analytiaue
 *
 * Revision 1.7  2001/04/23 09:10:16  smysona
 * Modif pour les vols a caler et les portfolios
 *
 * Revision 1.6  2001/04/03 10:45:35  nicolasm
 * Major revision
 *
 * Revision 1.5  2001/03/12 19:39:38  smysona
 * Optim des setmodels inutiles
 *
 * Revision 1.4  2001/03/06 17:10:29  smysona
 * Le frmana de calibration est construit avec les epochs du frmmc
 *
 * Revision 1.3  2001/02/22 19:46:28  smysona
 * Ajout de l'autocalage cap
 *
 * Revision 1.2  2001/02/01 19:33:09  smysona
 * Decodage propre des Shapes
 *
 * Revision 1.1  2001/01/30 15:32:24  smysona
 * Initial revision
 *
 */


/*----------------------------------------------------------------------------*/


#include "firsttoinc.h"
#include "calibrator.h"
#include "swaption.h"
#include "portfolio.h"
#include "capfloor.h"
#include "frmana.h"
#include "bootstrapcalibration.h"
#include "volflat.h"
#include "merge.h"
#include "autocal.h"
#include "volcube.h"
#include "fromto.h"


/*!
    Automatic calibration constructor with one risk class.
    The type specifies if we provide cap or swaption
    volatility and smile. <BR>
    The smile can be a cube, but only a cube a smile, without ATM vols.
*/    



ARM_Frm_AutoCalibrator::ARM_Frm_AutoCalibrator(int type, 
                                 ARM_VolCurve* VolCurve, ARM_VolCurve* Smile,
                                 ARM_Vector* CorrelatedIndexes, 
                                 ARM_Vector* indexes,
                                 ARM_Matrix* correlations)
                         : ARM_FRM_InterpCalibrator(NULL,
                           CorrelatedIndexes, indexes, correlations)
{
    Init();

    switch (type)
    {
        case PT_IRG : 
        case PT_IRG_SWOPT :
        {
            itsCapBlackVols = (ARM_VolCurve*) VolCurve->Clone();
            itsCapSmile     = (ARM_VolCurve*) Smile->Clone();
        }
        break;

        case PT_SWOPT :
        case PT_SWOPT_IRG :
        {
            itsSwpBlackVols = (ARM_VolCurve*) VolCurve->Clone();
            itsSwapSmile    = (ARM_VolCurve*) Smile->Clone();
        }
        break;

        default :
            throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET,
                            "Invalid autocal type");
    }
}



/*!
    Automatic calibration constructor with two risk classes
     (caps and swaptions volatility and smile). <BR>
    The smile can be a cube, but only a cube a smile, without ATM vols.
*/    

ARM_Frm_AutoCalibrator::ARM_Frm_AutoCalibrator(ARM_VolCurve* VolIdxCurve,
                                 ARM_VolCurve* IdxSmile,
                                 ARM_VolCurve* VolIrgCurve,
                                 ARM_VolCurve* IrgSmile,
                                 ARM_Vector* CorrelatedIndexes, 
                                 ARM_Vector* indexes,
                                 ARM_Matrix* correlations)
                  : ARM_FRM_InterpCalibrator(NULL, CorrelatedIndexes,
                                 indexes, correlations)
{
    Init();

    itsCapBlackVols = (ARM_VolCurve*) VolIrgCurve->Clone();
    itsCapSmile     = (ARM_VolCurve*) IrgSmile->Clone();
    itsSwpBlackVols = (ARM_VolCurve*) VolIdxCurve->Clone();
    itsSwapSmile    = (ARM_VolCurve*) IdxSmile->Clone();
}


ARM_Frm_AutoCalibrator::~ARM_Frm_AutoCalibrator(void)
{
    cleanit(itsSwpBlackVols);
    cleanit(itsCapBlackVols);

    cleanit(itsSwapSmile);
    cleanit(itsCapSmile);


    cleanit(itsCalibPortfolio);

    cleanit(itsMatCurve);
    cleanit(itsValCurve);
    cleanit(itsMatCorr);

    cleanit(itsSwp2CalVols);
    cleanit(itsCap2CalVols);
    cleanit(itsCalVols);
}



/*!
    Creates the portfolio that will be used during the bootstrap. <BR>
    Analytical prices are priced "at hand" in this function.
    The portfolio can only
    contain vanilla caps or swaptions. <BR>
    We also build the matCurve for the bootstrap. <BR> <BR>

    The dates chosen for the matCurve may be inapropriate for 
    a cap diagonal pricing.
*/    

void ARM_Frm_AutoCalibrator::CreateCalibPortfolio(void)
{
    double amount, discount, optMat, optType, term;
    double fixCoupon, parCoupon, smile, vol, price;
    int i, size, iPort;
    ARM_Date maturity, settle;

    if (GetCalibratingModel())
        delete GetCalibratingModel();

    SetCalibratingModel(NULL);


    switch (itsPricedSecurity->CalibrationPT())
    {
        case PT_IRG : 
            if (( itsCapBlackVols == NULL ) || ( itsCapSmile == NULL ))
                throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET,
                "No available vol/smile for IRG calibration");
            break;

        case PT_SWOPT : 
            if (( itsSwpBlackVols == NULL ) || ( itsSwapSmile == NULL ))
               throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET,
                    "No available vol/smile for SWOP calibration");
            break;
    }

    itsPricedSecurity->SetModel(GetModel());

    itsCalibPortfolio = 
           itsPricedSecurity->CalibrationPF(GetModel()->GetStartDate());

    ARM_Portfolio*  itsNewCalibPortfolio = NULL;

    size = itsCalibPortfolio->GetSize();

    cleanit(itsCalVols);

    itsCalVols = new ARM_Vector(size);

    
    ARM_Vector* matCurve     = new ARM_Vector(size);
    double*     matCurveElt  = matCurve->GetElt();


    ARM_Swaption* currentSwaption;
    ARM_CapFloor* currentCap;


    for (i = iPort = 0; i < size; i++, iPort++)
    {
        switch (itsCalibPortfolio->GetAsset(i)->GetName())
        {
            case ARM_SWAPTION :
            {
                if (( itsSwpBlackVols == NULL ) || ( itsSwapSmile == NULL ))
                {
                   delete itsCalibPortfolio->GetAsset(i);

                   itsCalibPortfolio->SetAsset(NULL, i);
                   itsCalibPortfolio->SetWeight(0,i);
                   itsCalibPortfolio->SetPrice(0,i);

                   iPort--;
                   break;
                }

                currentSwaption = (ARM_Swaption *) 
                                  itsCalibPortfolio->GetAsset(i);    

                term = (currentSwaption->GetEndDate().GetJulian()
                        -currentSwaption->GetStartDate().GetJulian()) 
                          / (double) K_YEAR_LEN;

                maturity = currentSwaption->GetExpiryDate();

                settle = GetModel()->GetStartDate();

                matCurveElt[i] = optMat = (maturity.GetJulian()
                                -settle.GetJulian()) / (double) K_YEAR_LEN;    

                optType = currentSwaption->IsPayer()
                             -currentSwaption->IsReceiver();

                fixCoupon = currentSwaption->GetFixedLeg()->GetFixedRate();

                currentSwaption->SetModel(GetModel());

                parCoupon = currentSwaption->PriceToRate(maturity, 0.0);

                smile = itsSwapSmile->ComputeVolatility(maturity, 
                                       fixCoupon-parCoupon, term);

                vol = itsSwpBlackVols->ComputeVolatility(maturity, term);

                vol += smile;

                int freq0 = currentSwaption->GetFixedLeg()->GetIRIndex()->GetPayFrequency();
                double stk0 = fixCoupon/100.0;

                if( (  strcmp(currentSwaption->GetFixedLeg()->GetCurrencyUnit()->GetCcyName(),"USD") == 0  ||
                       strcmp(currentSwaption->GetFixedLeg()->GetCurrencyUnit()->GetCcyName(), "JPY") == 0 ) &&
                    ( freq0 !=2 ) )
                {
                    vol *= 2.0*(exp((freq0/2.0)*log(1.0+stk0/freq0))-1.0)
                             *exp((1.0-freq0/2.0)*log(1.0+stk0/freq0))/stk0;
                }
                else if ( strcmp(currentSwaption->GetFixedLeg()->GetCurrencyUnit()->GetCcyName(),"EUR") == 0 &&
                    ( freq0 !=1 ) )
                {
                    vol *= ( exp(freq0*log(1.0+stk0/freq0))-1.0)
                            *exp((1-freq0)*log(1.0+stk0/freq0))/stk0;
                }

                itsCalVols->Elt(i) = vol;

                vol /= 100.0;

                discount = currentSwaption->GetFixedLeg()->Compute1BP();

                amount = currentSwaption->GetFixedLeg()->GetRcvOrPay()
                    *currentSwaption->GetAmount()->CptReferenceValue(maturity);

                price = amount*discount*bsOption(parCoupon, fixCoupon, 
                         vol, 0.0, 0.0,  optMat, optType);        

                itsCalibPortfolio->SetPrice(price,i);


                if ( fabs(price) > 1.e-20 )
                   itsCalibPortfolio->SetWeight(1/(price*price), i);
                else
                   itsCalibPortfolio->SetWeight(0.0, i);

            }
            break;

            case ARM_CAPFLOOR :
            {
                if (( itsCapBlackVols == NULL ) || ( itsCapSmile == NULL ))
                {
                    delete itsCalibPortfolio->GetAsset(i);

                    itsCalibPortfolio->SetAsset(NULL, i);
                    itsCalibPortfolio->SetWeight(0,i);
                    itsCalibPortfolio->SetPrice(0,i);

                    iPort--;
                    break;
                }

                currentCap = (ARM_CapFloor*) itsCalibPortfolio->GetAsset(i);    

                fixCoupon = currentCap->GetStrike();

                currentCap->SetModel(GetModel());

                settle = GetModel()->GetStartDate();

                term = (currentCap->GetSwapLeg()->GetEndDate().GetJulian()
                        -currentCap->GetSwapLeg()->GetStartDate().GetJulian()) 
                        / (double) K_YEAR_LEN;

                maturity = currentCap->GetExpiryDate();

                matCurveElt[i] = optMat = (maturity.GetJulian()
                                 -settle.GetJulian()) / (double) K_YEAR_LEN;    


                vol = itsCapBlackVols->ComputeVolatility(maturity, term);

                // CHECK
                parCoupon = currentCap->GetSwapLeg()->ComputePrice()
                              /currentCap->GetSwapLeg()->Compute1BP();

                parCoupon/= 100.0;

                smile = itsCapSmile->ComputeVolatility(maturity, 
                                        fixCoupon-parCoupon, term);

                vol += smile;

                itsCalVols->Elt(i) = vol;                

                ARM_VolFlat capVol(itsCapBlackVols->GetAsOfDate(),vol);

                ARM_BSModel bsmodel(GetModel()->GetZeroCurve(), &capVol,
                                    K_YIELD);

                currentCap->SetModel(&bsmodel);

                price = currentCap->ComputePrice();        

                currentCap->SetModel(GetModel());
                
                itsCalibPortfolio->SetPrice(price,i);

                if ( fabs(price)>1.e-20 )
                   itsCalibPortfolio->SetWeight(1/(price*price), i);
                else
                   itsCalibPortfolio->SetWeight(0.0, i);
            }
            break;

            default :            
                throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET,
                   "Invalid calibration instrument");
        }
    }

    // Compact the mat curve
// YAN 02/2002
/*
    ARM_Vector* tmpVect = matCurve->Sort_Compact();

    cleanit(matCurve);

    itsCurveSize = tmpVect->GetSize();

    itsMatCurve= new double[size];

    itsValCurve= new double[size];

    MEMSET(itsValCurve, 0.0, itsCurveSize * sizeof(double));

    MEMCPY(itsMatCurve, tmpVect->GetElt(), itsCurveSize * sizeof(double));

    cleanit(tmpVect);
*/
    ARM_Vector* tmp0Vect = matCurve->Sort_Compact();

    itsCurveSize = tmp0Vect->GetSize()+1;

    for (i = 0; i < tmp0Vect->GetSize(); i++)
    {
        if ( tmp0Vect->Elt(i) <= 0.0)
           itsCurveSize--;
    }

    ARM_Vector* tmpVect = new ARM_Vector(itsCurveSize);
    
    int j=0;

    for (i = 0; i < tmp0Vect->GetSize(); i++)
    {
        if ( tmp0Vect->Elt(i) > 0.0)
        {
           tmpVect->Elt(j++) = tmp0Vect->Elt(i);
        }
    }

    cleanit(matCurve);

    itsMatCurve= new double[itsCurveSize];

    itsValCurve= new double[itsCurveSize];


    MEMSET(itsValCurve, 0.0, itsCurveSize * sizeof(double));

    MEMCPY(itsMatCurve, tmpVect->GetElt(), (itsCurveSize-1) * sizeof(double));

    itsMatCurve[itsCurveSize-1] = itsMatCurve[itsCurveSize-2] + 0.25;

    cleanit(tmpVect);
    cleanit(tmp0Vect);
// Fin YAN

    //Clean portfolio from instruments we cannot price

    itsNewCalibPortfolio = new ARM_Portfolio(iPort);
    ARM_Vector* itsNewCalVols = new ARM_Vector(iPort);

    for (i = iPort = 0; i < size; i++)
    {
        if (itsCalibPortfolio->GetAsset(i))
        {
           itsNewCalibPortfolio->SetAsset(itsCalibPortfolio->GetAsset(i),iPort);
           itsNewCalibPortfolio->SetPrice(itsCalibPortfolio->GetMktPrices()->Elt(i), iPort);
           itsNewCalibPortfolio->SetWeight(itsCalibPortfolio->GetWeights()->Elt(i), iPort);

           itsNewCalVols->Elt(iPort) =  itsCalVols->Elt(i);

           iPort++;
        }
    }

    cleanit(itsCalVols);
    itsCalVols = itsNewCalVols;

    int nSwopt = 0;
    int nCap   = 0;

    for (i = 0; i < itsNewCalibPortfolio->GetSize(); i++)
    {
        switch (itsCalibPortfolio->GetAsset(i)->GetName())
        {
            case ARM_SWAPTION : nSwopt++; break;
            case ARM_CAPFLOOR : nCap++; break;

            default :
               throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET,
                "Invalid calibration instrument for auto calibration");

        }
    }

    cleanit(itsSwp2CalVols);
    cleanit(itsCap2CalVols);

    itsSwp2CalVols = new ARM_Vector(nSwopt);
    itsCap2CalVols = new ARM_Vector(nCap);

    int iCap, iSwopt;

    for (i = iCap = iSwopt = 0; i < itsNewCalibPortfolio->GetSize(); i++)
    {
        switch (itsCalibPortfolio->GetAsset(i)->GetName())
        {
            case ARM_SWAPTION :
                itsSwp2CalVols->Elt(iSwopt++) = itsCalVols->Elt(i); break;

            case ARM_CAPFLOOR : 
                itsCap2CalVols->Elt(iCap++) = itsCalVols->Elt(i); break;

            default :
               throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET,
                "Invalid calibration instrument for auto calibration");

        }
    }    

    cleanit(itsCalibPortfolio);

    itsCalibPortfolio = itsNewCalibPortfolio;

    if (itsMatCorr)
        cleanit(itsMatCorr);

    if ( GetCorrelations()->size() > 0 )
    {
        int nLig, nCol, i, j;

        nLig = GetCorrelatedIndexes()->GetSize();
        nCol = GetIndexes()->GetSize();

        itsMatCorr = new ARM_Matrix(nCol,nLig,0.0);

        for (i=0; i< nLig; i++)
            for (j=0; j < nCol; j++)
                itsMatCorr->Elt(j,i) = GetCorrelations()->at(i)[j];
    }

}



/*!
    Creates the analytical model that will be bootstrapped, and on which we will
    interpolate variances for the model we want to calibrate.
*/    

void ARM_Frm_AutoCalibrator::CreateAnaModel(void)
{
    double    decay, slope, asymptote;


    int shapeSourceType = GetModel()->GetRowShapeSource()->GetSource();

    if (IsShapeFlat(shapeSourceType))
    {
       decay = 0;
       slope = 0;
       asymptote = 0;
    }
    else if (IsShapeFunc1(shapeSourceType))
    {
       decay = GetModel()->GetRowShapeSource()->GetAlpha1();
       slope = GetModel()->GetRowShapeSource()->GetSlope1();
       asymptote = GetModel()->GetRowShapeSource()->GetAsymptote1();
    }
    else
    {
       throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET,
                "Invalid shape for auto calibration");
    }

    ARM_Frm* frmmod = dynamic_cast<ARM_Frm*> (GetModel());
    

    // Include expriry dates of calibration instruments in the epochs
    // of the analytic model

    int i, size = itsCalibPortfolio->GetSize();
    ARM_Vector* epochsDates = frmmod->GetEpochs()->get_NewJulianVect();
    ARM_Vector* expiryDates = new ARM_Vector(size);
    ARM_Vector* newEpochs   = NULL;

    for (i = 0; i < size; i++)
    {
        expiryDates->Elt(i) = 
                  itsCalibPortfolio->GetAsset(i)->GetExpiryDate().GetJulian();
    }

    MergeDates(&newEpochs, epochsDates, expiryDates);

    SetCalibratingModel(new ARM_FrmAna(GetModel()->GetZeroCurve(), 
                        newEpochs, 0, K_ROW, 
                        decay, slope, asymptote, GetModel()->getNumOfFactors(),
                        GetCorrelatedIndexes(), GetIndexes(), itsMatCorr));

    cleanit(epochsDates);
    cleanit(expiryDates);
    cleanit(newEpochs);
}



/*!
    Launch the bootstrap.
*/    

double  ARM_Frm_AutoCalibrator::BootstrapAnaModel(int mode)
{
//    int size = itsCalibPortfolio->GetSize();

//    double error = ComputeParasCurve(GetCalibratingModel(), itsCalibPortfolio, itsCalibPortfolio->GetSize(),

    int size = GetCurveSize();

    double error;

    if (GetModel()->GetRowShapeSource()->IsRowShape())
    {
       error = ComputeParasCurve(GetCalibratingModel(), itsCalibPortfolio, 
                                 size,
                                 itsMatCurve, itsAccuracy, itsMinVol, 
                                 itsMaxVol, itsNbMaxStep, itsValCurve);
    }
    else
    {
        GetCalibratingModel()->SwitchToDiagShape();

        error = ComputeParasCurveDiag(GetCalibratingModel(), 
                                      itsCalibPortfolio, size,
                                      itsMatCurve,  itsAccuracy, itsMinVol,
                                      itsMaxVol, itsNbMaxStep, itsValCurve);
    }

    return error;                      
}


/*!
    Generates the factot strip, build the portfolio, bootstrap 
    ana model and interpolate volatilities.
*/    
void ARM_Frm_AutoCalibrator::GenerateVolMatrix(int vol)
{
    GenerateStructure();
    
    
    CreateCalibPortfolio();    
    CreateAnaModel();
    PrepareCalibration();

    BootstrapAnaModel();

    ARM_FRM_InterpCalibrator::GenerateVolMatrix(vol);


    EndCalibration();    
}


/*!
    SetModel of the calibrating instruments. <BR>
    Choice of the smile type (sitcky strike or sticky delta).
*/    
void ARM_Frm_AutoCalibrator::PrepareCalibration(void)
{
    int i,size = itsCalibPortfolio->GetSize();

    for (i = 0; i < size; i++)
    {
        itsCalibPortfolio->GetAsset(i)->SetModel(GetCalibratingModel());
        itsCalibPortfolio->GetAsset(i)->InitHedge();
    }

    // Set sticky strike
// YAN 02/2002
    SetSwoptSmileStick(1);
}


/*!
    Generates output and cleans memory.
*/

void ARM_Frm_AutoCalibrator::EndCalibration(void)
{
    int i,size = itsCalibPortfolio->GetSize();

    ARM_Model* model = GetModel();

    for (i = 0; i < size; ++i)
    {
        itsCalibPortfolio->GetAsset(i)->SetModel(model);
    }
    

    ARM_FRM_InterpCalibrator::EndCalibration();

    ARM_Frm_AutoCalibrator::Output();

    itsCalibPortfolio->FreePortfolioAndAssets();

    /*
    if (GetSwoptSmileStick() == 1)
    {
        SwitchSwopt2StickyStrike();
    }
    */
    
    cleanit(itsCalibPortfolio);

    cleanit(itsMatCurve);
    cleanit(itsValCurve);
    cleanit(itsMatCorr);
}


/*!
    Calibration report output (used by Summit pricers).
*/

void ARM_Frm_AutoCalibrator::Output(void)
{
    ARM_Frm* frmmod = NULL;
    SmartCast(GetModel(), frmmod);

    double price, vol, strike;


    if (frmmod)
    {
       if (frmmod->GetCalibOutpoutFilename())
       {
          ofstream fout(frmmod->GetCalibOutpoutFilename(),
                        ios::out | ios::app);
          int i;
          char aDate[10+1];

          fout << endl;
          fout << " AutoMatic calibration" << endl;

          switch (itsPricedSecurity->CalibrationPT())
          {
               case PT_SWOPT_IRG :
               case PT_SWOPT     :
               {
                   fout << " Bootstrapping on swaptions" << endl;
                   ARM_Swaption* swop = NULL;

                   fout<<"[n]   Expiry       Maturity     Strike   MktPrice  CalibPrice   MktVol   ModelVol" << endl;

                   for (i = 0; i < itsCalibPortfolio->GetSize(); i++)
                   {
                       swop = (ARM_Swaption*) itsCalibPortfolio->GetAsset(i);
                       price  = swop->ComputePrice();
                       vol    = swop->ComputeImpliedVol(price);

                       strike = swop->GetStrike();

                       fout << "[" << i << "]   ";

                       swop->GetExpiryDate().JulianToStrDate(aDate);
                       fout << aDate << "   ";
                       swop->GetEndDate().JulianToStrDate(aDate);
                       fout << aDate << "   ";
                       fout << strike << "      ";
                       fout << itsCalibPortfolio->GetMktPrices()->Elt(i) << "   ";
                       fout << price << "      ";
                       fout << itsSwp2CalVols->Elt(i) << "       ";
                       fout << vol << endl;
                   }
               } break;

               case PT_IRG_SWOPT : 
               case PT_IRG :

               default :           ;
            }

            fout.close();
       }
    }
    
}


/*!
    Switch input smile from sticky delta to sticky strike.
    We need to store the ATM level
    for each expiry for interpolation purpose.
*/    

void ARM_Frm_AutoCalibrator::SwitchSwopt2StickyStrike(void)
{
    int i,j;

    ARM_Vector* tenors = NULL;

    int nExpiry = itsSwpBlackVols->GetVolatilities()->GetNumLines();

    ARM_VolCube* newVolCube = NULL;

    switch(itsSwapSmile->GetName())
    {
        case ARM_VOL_LIN_INTERPOL :
            tenors = new ARM_Vector(1, 10.0);
            newVolCube = new ARM_VolCube(itsSwpBlackVols, 
                                         &itsSwapSmile, 1, tenors);
        break;

        case ARM_VOL_CUBE :
            tenors = ((ARM_VolCube*) itsSwapSmile)->GetUnderLyings();
            newVolCube = (ARM_VolCube*) ((ARM_VolCube*) itsSwapSmile)->Clone();
            newVolCube->SetATMref(true);
            newVolCube->SetATMVol((ARM_VolCurve*) itsSwpBlackVols->Clone());
        break;

        default :
            throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET,
                "Invalid volatility type");
        
    }

    // switch to sticky strike
    newVolCube->SetStickyDelta(false);

    double* expiryTermsElt = itsSwapSmile->GetExpiryTerms()->GetElt();
    int nTenor  = tenors->GetSize();
    double* tenorTermsElt  = tenors->GetElt();

    ARM_Date start = itsSwpBlackVols->GetAsOfDate();
    ARM_Currency*  ccy = itsSwpBlackVols->GetCurrency();

    ARM_INDEX_TYPE freq = (ARM_INDEX_TYPE) 
        FromFrequencyToXiborType(ccy->GetLiborTerm(), ccy->GetCcyName());

    start.NextBusinessDay(ccy->GetSpotDays() ,ccy->GetCcyName());


    double swapRate;
    ARM_Matrix* strikeLevels = newVolCube->GetStrikeLevels();

    for (i = 0; i < nExpiry; i++)
    {
        for (j = 0; j < nTenor; j++)
        {
            ARM_Date startSwap = start;
            startSwap.AddDays(expiryTermsElt[i] * K_YEAR_LEN);

            ARM_Date endSwap = startSwap;
            endSwap.AddDays(tenorTermsElt[j] * K_YEAR_LEN);

            ARM_Swap* parswap = new ARM_Swap(startSwap, endSwap, 
                                             freq,  
                                             0.0, K_MARKET_RATE,
                                             K_RCV, K_DEF_FREQ,
                                             K_DEF_FREQ, ccy);

            parswap->SetModel(GetModel());

            swapRate =  parswap->CptMarketSwapRate();

            strikeLevels->Elt(j,i) = swapRate;
        }
    }

    SetSwoptBlackCurveForPricing(newVolCube);
}




/*----------------------------------------------------------------------------*//*---- End Of File ----*/
