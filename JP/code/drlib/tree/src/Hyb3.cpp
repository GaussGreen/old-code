//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : Hyb3.cpp
//
//   Description : hyb3 FD model class implementation
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Hyb3.hpp"
#include "edginc/IRVolPair.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/IrConverter.hpp"
#include "edginc/IndexSpecIR.hpp"
#include "edginc/IndexSpecFX.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Results.hpp"
#include "esl_log.h"
#include "esl_market.h"
#include "edginc/MDFUtil.hpp"
#include "cupslib.h"

DRLIB_BEGIN_NAMESPACE

/**************************************** Hyb3 ************************************/

static
int HYB3CorrelOverwrites(
    FX_DATA            *fx_data,                    // (O) Fx data
    char               OverWriteString[3][MAXBUFF], // (I) Overwrite strings
    CORRELATION_DATA   *correlation_data)
{
    double temp[3];

    int status = FAILURE;             /* Error status = FAILURE initially */

              
    temp[0] = correlation_data->CorrIR;
    temp[1] = correlation_data->CorrDomIRFX;
    temp[2] = correlation_data->CorrForIRFX;


    /*********************************************/
    /* Although the fx_data structure supports a */
    /* term-structure of correlations, we only   */
    /* use the first element in Build Tree       */
    /*********************************************/
    fx_data->Rho[0][0] = temp[0];   /* Flat correlation curves */
    fx_data->Rho[1][0] = temp[2];
    fx_data->Rho[2][0] = temp[1];   /* Reverse here: Rho[1]=temp[2] and Rho[2]=temp[1] */
                
    if (strcmp (OverWriteString[0], "nil"))                      /* Correlation overwrite */
    {                    
         fx_data->Rho[0][0] = atof (OverWriteString[0]);/* same comment as above*/
    }
          
    if (strcmp (OverWriteString[1], "nil"))
    {
        fx_data->Rho[1][0] = atof (OverWriteString[1]);          
    }
          
    if (strcmp (OverWriteString[2], "nil"))
    {                    
        fx_data->Rho[2][0] = atof (OverWriteString[2]);                               
    }
          

    if (Hyb3_Correl_Check_WType3 (fx_data) != SUCCESS)                  /* Check validity of input */
    {          
         goto RETURN;
    }
          

    status = SUCCESS;
          
    RETURN:
          
    return (status);

}


static
int HYB3FXSmileOverwrites(
    FX_DATA            *fx_data,                 // (O) Fx data
    char               OverWriteString[MAXBUFF], // (I) FX Owrite str
    FXVOLATILITY_DATA  *Volatility_data,         // (I) FXVols file name
    FXSMILE_DATA       *Smile_data,              // (I) FX smile file name
    char               NbFXSmileOWS[MAXBUFF],    // (I)
    char               FXSmileParamOWS[MAXNBDATE][MAXBUFF])
{
     long
          Year[2],
          Month[2],
          Day[2];
     int
          i,
          readerror,                 /* Reading error status             */
          status = FAILURE;          /* Error status = FAILURE initially */


     fx_data->Today = Volatility_data->ValueDate;
          
     fx_data->ValueDate = fx_data->Today;
     fx_data->SpotDays  = 0;

     fx_data->Spot = Volatility_data->FXSpotRate;
     // ignore BaseVolFreq entry

     fx_data->NbVol = Volatility_data->NbBaseVols;
     if (fx_data->NbVol > MAXNBDATE)
     {
          DR_Error("Nb of vols exceeds max limit of %d in "
                             "(Fx_Input_W_WithSmile_DRI)", MAXNBDATE);
          
          goto RETURN;
     }
     
     Dsplit(fx_data->ValueDate, /* Split value date into month, day and year */
              &(Month[0]), 
              &(Day[0]), 
              &(Year[0]));

     for (i = 1; i <= fx_data->NbVol; i++)  /* Composite vols -> 1 offset*/
     {
        fx_data->VolDate[i] = Volatility_data->BaseVolDates[i-1];
        fx_data->FxVol[i] = Volatility_data->BaseVols[i-1]/100.0;

     }  /* for i */

     /* Repeat the first values for interpolation of pseudocomposite vols */
     /* in Get_TreeSpotVols                                               */

    if (fx_data->NbVol > 0)
    {
        fx_data->VolDate[0] = fx_data->ValueDate;
        fx_data->FxVol[0]   = fx_data->FxVol[1];
    }
     

    fx_data->NbInpSpotVol = Volatility_data->NbSpotVols;
    if (fx_data->NbInpSpotVol > MAXNBDATE)
    {
         DR_Error("Nb of Inp Spot Vol exceeds max limit of %d "
                   "in ! (Fx_Input_W)", MAXNBDATE);
         
         goto RETURN;
    }

    if (fx_data->NbInpSpotVol > 0)
    {
       for (i = 0; i < fx_data->NbInpSpotVol; i++)   /* Spot vols -> 0 Offset*/
       {
             fx_data->InpSpotVolDate[i] = Volatility_data->SpotVolDates[i];
             fx_data->InpSpotVol[i] = Volatility_data->SpotVols[i]/100.0;
       }
    }


    /* deal with possible overlap between composite and spot vols */
    if ((fx_data->NbInpSpotVol > 0) && (fx_data->NbVol > 0))
    {
        long Idx;
        Idx = GetDLOffset(fx_data->NbVol,
                          &(fx_data->VolDate[1]),
                          fx_data->InpSpotVolDate[0],
                          CbkHIGHER);
        /*******************************************************/
        /* if all composite dates are < spotvoldate[0]         */
        /* then Idx = -999, in which case NbVol does not need  */
        /* to be changed.                                      */
        /*******************************************************/
        if (Idx >= 0L)
        {
            fx_data->NbVol = Idx;
        }        
    }
    
    if (strcmp (OverWriteString, "nil") == 0)    /* 1.Disallow cut off      */
    {

        fx_data->FxCutOffFlag = FALSE;
        fx_data->FxCutOffLast = FALSE;
        fx_data->FxCutOffLevel   = FXSPOT_CUTOFF_R; /* Ratio to fwd vol level*/
        fx_data->FxBootStrapMode = FX_NO_FAILURE_ALLOWED; 

     }
     else if (strcmp(OverWriteString,"last") == 0) /* 2.Cut off at last level */
     {
        /* Cut off is allowed and there */
        /* is volatility cutoff data    */
        fx_data->FxCutOffFlag = TRUE;
        fx_data->FxCutOffLast = TRUE;
        fx_data->FxCutOffLevel = FXSPOT_CUTOFF_R; /* Ratio to fwd vol level  */
        fx_data->FxBootStrapMode = FX_USE_LAST_LEVEL ;
     } 
     else                                         /* 3.Cut off at user level */
     {
          /* Cut off is allowed and there */
          /* is volatility cutoff data    */
          fx_data->FxCutOffFlag = TRUE;
          fx_data->FxCutOffLast = FALSE;
          fx_data->FxCutOffLevel = atof(OverWriteString);
          if (ABS(fx_data->FxCutOffLevel) < TINY) /* The FAILURE return */
          {
              DR_Error("Unable to convert FX cut off value!\n");
              goto RETURN;
          }
          fx_data->FxCutOffLevel /= 100.;
          fx_data->FxBootStrapMode = FX_CONSTANT_SPOT_VOL;
          
     }
     
    // Read FX Smile data
    if (strstr(NbFXSmileOWS, "nil") != NULL)
    {

        // if smile data structure is NULL, turn off the FX smile by setting all
        // the smile parameters to 0.0
        if (Smile_data == NULL)
        {
            fx_data->NbSmilePt = 1;

            /* add arbitary number of days after the value date to ensure the smile date is
               > today as in some cases valueDate = Today. - chose to use 5 days but could have
               use any number */
            fx_data->SmileDate[0] = Nxtday(fx_data->ValueDate, 5);
            fx_data->a1[0] = 0.0;
            fx_data->a2[0] = 0.0;
            fx_data->a3[0] = 0.0;
        }
        else
        {
            fx_data->NbSmilePt = Smile_data->NbSmileDates;
            if ((fx_data->NbSmilePt < 0) ||
                (fx_data->NbSmilePt > MAXNBDATE))
            {
                DR_Error("Fx_Input_W: Nb FX smile param lines out of range!");
                goto RETURN;
            }

            for (i=0; i < fx_data->NbSmilePt; i++)
            {
                fx_data->SmileDate[i] = Smile_data->SmileDates[i];
                fx_data->a1[i] = Smile_data->A1[i];
                fx_data->a2[i] = Smile_data->A2[i];
                fx_data->a3[i] = Smile_data->A3[i];
            }  
        }
    }
    else
    {
        /* read the overwrite strings */
        readerror = sscanf (NbFXSmileOWS, 
                            "%d \n", 
                            &(fx_data->NbSmilePt));

        if (readerror != 1)
        {        
            DR_Error("Fx_Input_W: Cannot read Nb FX Smile Param OWS");
            goto RETURN;
        }

        if ((fx_data->NbSmilePt < 0) ||
            (fx_data->NbSmilePt > MAXNBDATE))
        {
            DR_Error("Fx_Input_W: Nb FX smile param lines OWS out of range!");
            goto RETURN;
        }

        for (i = 0; i < fx_data->NbSmilePt; i++)
        {
            readerror = sscanf (FXSmileParamOWS[i], 
                                "%ld\t%lf\t%lf\t%lf \n",
                                &(fx_data->SmileDate[i]),
                                &(fx_data->a1[i]),
                                &(fx_data->a2[i]),
                                &(fx_data->a3[i]));

            if (readerror != 4)
            {          
                DR_Error("Fx_Input_W: Cannot read FX smile params OWS");
                goto RETURN;
            }
        }
    }


     /* Check validity of input */ 
     if (Hyb3_Fx_Check_W (fx_data) != SUCCESS)
     {
          goto RETURN;
     }

     status = SUCCESS;
          
     RETURN:
          
     return (status);

}

void Hyb3::load(CClassSP& clazz){

    clazz->setPublic();
    REGISTER(Hyb3, clazz);
    SUPERCLASS(RateTree);

    FIELD(cetSmoothing,"");
    FIELD_MAKE_OPTIONAL(cetSmoothing);

    // this fields are transient to allow tweaking to copy 
    // the model after ::getMarket() has been called 
    FIELD(mToday, "");
    FIELD_MAKE_TRANSIENT(mToday);
}

CClassConstSP const Hyb3::TYPE = CClass::registerClassLoadMethod(
    "Hyb3", typeid(Hyb3), Hyb3::load);

// for class loading
bool Hyb3Load(void) { return (Hyb3::TYPE != 0); }


void Hyb3::initModel()
{
try {
    int status = FAILURE;          /* Error status = FAILURE initially */
    mTreeBuilt = true; 

    EslSetZeroInterpolation(zeroInterpStyle == ZeroInterpStyle::FLAT_FWD ?
        ESL_INTERP_FLATFWD : ESL_INTERP_LINEAR);
    //EslSetZeroInterpolationStub(ESL_INTERP_LINEAR);

    mToday = DateTime::fromIrDate(mFxData.ValueDate);

    CRIT_DATE* critDate=NULL;
    int nbCritDate = 0;
    IrConverter::AutoDeleteDRArray<CRIT_DATE>
        critDateDelete(&critDate, CRITDATE, &nbCritDate);

    int i,j;
    int printCET;

    map<string, ZeroBond>::iterator zero;
    map<string, NamedZeroBank>::iterator zb;

    // needed to setup timeline
    critDate = (CRIT_DATE *) DR_Array (CRITDATE, 0, 0);
    if (critDate == NULL)
        goto RETURN;

    // value date always critical date
    if (Add_To_DateList(&nbCritDate, 
                        &critDate,
                        getValueDate().toIrDate(),
                        0, 0, 0, 0, 0, 0, 0, 0, 0) != SUCCESS)
        goto RETURN;

//    arrangeDates();

    for (j = 0; j < critDates.size(); j++)
    {

        DateTime critDateL = critDates[j];
        if (critDateL > getValueDate())
        {
            if (Add_To_DateList(&nbCritDate, 
                                &critDate,
                                critDateL.toIrDate(),
                                0, 0, 0, 0, 0, 0, 0, 0, 0) != SUCCESS)
                goto RETURN;
        }
    }

     // if in CLAIMBANK mode, insert critical dates associated with the claimbank
    if (zeroBankMode == ZeroBankMode::CLAIMBANK)
    {
        DR_Error("Claim bank not yet suppported");
        goto RETURN;
    }

    if (Sort_CritDate (nbCritDate, critDate) != SUCCESS)
        goto RETURN;

    sortedCritDates.resize(1);
    sortedCritDates[0] = critDate[0].CritDate;
    for (i = 1; i < nbCritDate; i++)
    {
        if (sortedCritDates.back() == critDate[i].CritDate)
            continue;
        sortedCritDates.push_back(critDate[i].CritDate);
    }
/*
    if (!treeDataDebugFile.empty())
    {
        FILE* stream = fopen("c:\\hyb3CritDates.dat", "w");
        if (stream)
        {
            // sort dates before printing them
            if (Sort_CritDate(nbCritDate, critDate) != SUCCESS)
                goto RETURN;

            PrintCriticalDates(stream, nbCritDate, critDate);
            fclose(stream);
        }
    }
*/
    // ??? for now assume all critical dates are discrete - may have to be revised
    // if we want to approximate american/continous option/exercise/observation style
    // products as we currently do in the wrappers
    for (i = 0; i < NBCRITDATE; i++)
    {
        mTreeData.CritType[i] = 'D';
        mTreeData.NbZeros[i] = 0;
    }

    // make the timeline 
    if (Hyb3_Time_Line (mFxData.ValueDate, nbCritDate, critDate, 'I', &mTreeData) != SUCCESS)
        goto RETURN;

    // if printing out a debug file, turn on printing the default CET debug file also
    printCET = FALSE;
    if (!treeDataDebugFile.empty())
        printCET = TRUE;
    if (Hyb3_Build_Tree(printCET, mTCurves, mMktVolData, &mFxData, &mEqData, 
                        &mTreeData) != SUCCESS)
    {
        goto RETURN;
    }

    // allocate ranges/tree geometry structure referenced by all tree slices
    switch(mTreeData.TreeType)
    {
        case TTYPE_1IR:
        {
            range.reset(new TreeSliceRates::Range(
                -mTreeData.HalfWidth[0],
                mTreeData.HalfWidth[0]));  // ??? check what this is doing (rates=true)
            break;
        }
        case TTYPE_2IR:
        case TTYPE_EQ1IR:
        {
            range.reset(new TreeSliceRates::Range(
                -mTreeData.HalfWidth[0],
                mTreeData.HalfWidth[0],
                -mTreeData.HalfWidth[1],
                mTreeData.HalfWidth[1]));
            break;
        }
        case TTYPE_FX2IR:
        case TTYPE_EQD2IR:
        case TTYPE_EQF2IR:
        case TTYPE_EQC2IR:
        case TTYPE_2IR2F1D:
        {
            range.reset(new TreeSliceRates::Range(
                -mTreeData.HalfWidth[0],
                mTreeData.HalfWidth[0],
                -mTreeData.HalfWidth[1],
                mTreeData.HalfWidth[1],
                -mTreeData.HalfWidth[2],
                mTreeData.HalfWidth[2]));
            break;
        }
        default:
        {
            DR_Error("Specified hyb3 tree mode (enum = %d) invalide", mTreeData.TreeType);
            goto RETURN;
        }
    }

    mTpIdx = mTreeData.NbTP;

    if (zeroBankMode == ZeroBankMode::ZEROBANK)
    {   
        for (zb = mZeroBank.begin(); zb != mZeroBank.end(); zb++)
        {
            string crvName = zb->second.curveName;  // ??? to help debugging
            int currIdx = getCurrIdx(crvName);
            int crvIdx = getCrvIdx(crvName, currIdx);

            // ??? tidy up who defines currency and members of zeroBanks
            // assume for time being that foreign is always before domestic 
            // asset in tree building and are the two first assets
            // ??? this structure might need clean up in the underlying Hyb3 code (make it more general)
            if (zb->second.currency == mCurrency[FOREIGN])
            {
                // zb->second.zeroBankDim = 1;
                zb->second.zeroBankDim = mMktVolData[FOREIGN].NbFactor;
                zb->second.currencyEnum = FOREIGN;
                zb->second.discCurve = crvIdx;
            }
            else if (zb->second.currency == mCurrency[DOMESTIC])
            {
                // zb->second.zeroBankDim = 2;
                zb->second.zeroBankDim = mMktVolData[FOREIGN].NbFactor+mMktVolData[DOMESTIC].NbFactor;
                zb->second.currencyEnum = DOMESTIC;
                zb->second.discCurve = crvIdx;
            }
            else
            {
                DR_Error("Internal error - currIdx %d is not valid - must be 0(foreign) "
                         "or 1 (domestic)",
                         currIdx);
                goto RETURN;
            }

            switch (mTreeData.TreeType)
            {
                case TTYPE_1IR:
                {
                    zb->second.devMode = DISC_1D_NOCUPS;
                    break;
                }
                case TTYPE_2IR:
                case TTYPE_FX2IR:
                {
                    if (zb->second.currencyEnum == FOREIGN)
                        zb->second.devMode = DISC_1D_NOCUPS;
                    else if (zb->second.currencyEnum == DOMESTIC)
                        zb->second.devMode = DISC_2D_CUPS;
                    break;
                }
                case TTYPE_2IR2F1D: // AK, new dev-modes for the 2d+1 case ???
                    {
                        if (zb->second.currencyEnum == FOREIGN)
                            zb->second.devMode = DISC_2D_1IR2F_NOCUPS;
                        else if (zb->second.currencyEnum == DOMESTIC)
                            zb->second.devMode = DISC_3D_2IR2F1D_CUPS;
                        break;
                    }
                default:
                {
                    DR_Error("Required hyb3 tree mode (enum = %d) not supported", mTreeData.TreeType);
                    goto RETURN;
                }
            }

            int nbZeros = zb->second.nbZeros;
            zb->second.zeroCritDates.resize(nbZeros);
            zb->second.zeroBank.resize(nbZeros);

            for (i = 0; i < nbZeros; i++)
            {
                zb->second.zeroBank[i] = Hyb3_Alloc_Slice(&mTreeData, zb->second.zeroBankDim);
                if (zb->second.zeroBank[i] == NULL)
                {
                    DR_Error("Unable to allocate memory for zero banks");
                    goto RETURN;
                }
            }
        }

        // loop through each stored named zero bond
        for (zero = mZeroBond.begin(); zero != mZeroBond.end(); zero++)
        {
            int currIdx = getCurrIdx(zero->second.curveName);
            if (currIdx == 0)
            {
                zero->second.currencyEnum = FOREIGN;
                //zero->second.sliceDim = 1;
                zero->second.sliceDim = mMktVolData[FOREIGN].NbFactor;
                zero->second.currency = mCurrency[currIdx];
            }
            else if (currIdx == 1)
            {
                zero->second.currencyEnum = DOMESTIC;
                // zero->second.sliceDim = 2;
                zero->second.sliceDim = mMktVolData[FOREIGN].NbFactor + mMktVolData[DOMESTIC].NbFactor;
                zero->second.currency = mCurrency[currIdx];
            }
            else
            {
                DR_Error("Internal error - currIdx %d is not valid - must be 0(foreign) "
                         "or 1 (domestic)",
                         currIdx);
                goto RETURN;
            }

            switch (mTreeData.TreeType)
            {
                case TTYPE_1IR:
                {
                    zero->second.devMode = DISC_1D_NOCUPS;
                    break;
                }

                case TTYPE_2IR:
                case TTYPE_FX2IR:
                {
                    if (zero->second.currencyEnum == FOREIGN)
                        zero->second.devMode = DISC_1D_NOCUPS;
                    else if (zero->second.currencyEnum == DOMESTIC)
                        zero->second.devMode = DISC_2D_CUPS;
                    else
                    {
                        DR_Error("Internal Error - currency Enum is invalid value (%d) in zeroBond %s",
                                 zero->second.currencyEnum,
                                 (zero->first).c_str());
                        goto RETURN;
                    }
                    break;
                }
                case TTYPE_2IR2F1D: // AK, new dev-modes for the 2d+1 case ???
                    {
                        if (zero->second.currencyEnum == FOREIGN)
                            zero->second.devMode = DISC_2D_1IR2F_NOCUPS;
                        else if (zero->second.currencyEnum == DOMESTIC)
                            zero->second.devMode = DISC_3D_2IR2F1D_CUPS;
                        else
                        {
                            DR_Error("Internal Error - currency Enum is invalid value (%d) in zeroBond %s",
                                zero->second.currencyEnum,
                                (zero->first).c_str());
                            goto RETURN;
                        }
                        break;
                    }
                default:
                {
                    DR_Error("Required hyb3 tree mode (enum = %d) not supported", mTreeData.TreeType);
                    goto RETURN;
                }
            }

            zero->second.zeroSlice = Hyb3_Alloc_Slice(&mTreeData, zero->second.sliceDim);
            if (zero->second.zeroSlice == NULL)
            {
                DR_Error("Unable to allocate memory for named zero bond %s",
                         zero->second.getName().c_str());
                goto RETURN;
            }
        }
    }
    
    if (Hyb3_Dev_Alloc(&mDevData, &mTreeData) != SUCCESS)
    {
        DR_Error("Failed to allocate storage for DEV data.");
        goto RETURN;
    }

    print();

    status = SUCCESS;

    RETURN:

    if (status != SUCCESS)
        throw ModelException(IrConverter::slog.pop());
}
catch (exception& e) {
    throw ModelException(e, "Hyb3::initModel");  // catch anything else
}
}


//------------------------------------------------------
// Register the zero bond slices required by the pricer
//------------------------------------------------------
void Hyb3::insertZero(const ZeroBond& zero)
{
    const char* routine = "Hyb3::insertZero";
    int i;

    if (zero.startDates.size() != zero.matDates.size())
        throw ModelException(routine, "Number of start dates " + 
                             Format::toString((int)zero.startDates.size()) +
                             " must equal number of maturity dates " +
                             Format::toString((int)zero.matDates.size()));

    // check all the zero dates are in ascending order
    // don't strictly have to be, but makes things simpler
    for (i = 1; i < zero.matDates.size(); i++)
    {
        if (zero.matDates[i-1] >= zero.matDates[i])
            throw ModelException(routine, "zero mat dates for named zero " +
                                 zero.getName() +
                                 "must be in ascending order");

        if (zero.startDates[i-1] >= zero.startDates[i])
            throw ModelException(routine, "zero start dates for named zero " +
                                 zero.getName() +
                                 "must be in ascending order");
    }

    for (i = 0; i < zero.matDates.size(); i++)
    {
        if (zero.startDates[i] >= zero.matDates[i])
            throw ModelException(routine, "zero start date " + 
                                 zero.startDates[i].toString() + 
                                 " must be before the associated zero maturity date " +
                                 zero.matDates[i].toString());
    }

    // simply insert the zero -map ensures a unique name
    mZeroBond[zero.getName()] = zero;
}

void Hyb3::registerZero(
        DateTime obsDate, 
        DateTime matDate, 
        string curveName)
{
    try {
        if (matDate < obsDate) {
            throw ModelException("matDate < obsDate");
        }
        DateTimeArray startList, endList;
        startList.push_back(obsDate);
        endList.push_back(matDate);
        ZeroBondProdSP zeroProd = DYNAMIC_POINTER_CAST<ZeroBondProd>(createProduct(IProdCreatorSP(new ZeroBond(
            startList, endList, curveName))));

        string key = DateTime::dateFormat(obsDate.getDate())+"_"+DateTime::dateFormat(matDate.getDate());
        ZeroProdMap::iterator it = zeroProdList.find(key);
        if (it == zeroProdList.end())
            zeroProdList.insert(pair<string, ZeroBondProdSP>(key, zeroProd));
    }
    catch (exception& e) { 
        throw ModelException(e, __FUNCTION__); 
    }
}


void Hyb3::getZero(
        DateTime obsDate, 
        DateTime matDate, 
        string curveName, 
        TreeSliceSP &slice)
{
    try {
        if (!dynamic_cast<TreeSliceRates*>(slice.get())) {
            slice = RateTree::createSimpleSlice(curveName);
            slice->name = "Zero "+obsDate.toString()+" to "+matDate.toString();
        }

        if (obsDate==matDate) {
            *slice = 1.;
            return;
        }
        string key = DateTime::dateFormat(obsDate.getDate())+"_"+DateTime::dateFormat(matDate.getDate());
        ZeroProdMap::iterator it = zeroProdList.find(key);
        if (it == zeroProdList.end())
            throw ModelException("ZeroProd "+key+" not found");

        *slice = it->second->getValue(range->treeStep, matDate);
    }
    catch (exception& e) { 
        throw ModelException(e, __FUNCTION__); 
    }
}

void Hyb3::getZero(const ZeroBond& zero, int step, TreeSlice& treeSlice)
{
    const char* routine = "Hyb3::getZero";
    TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

    map<string, ZeroBond>::iterator zeroIter = mZeroBond.find(zero.getName());
    if (zeroIter == mZeroBond.end())
        throw ModelException(routine, "Unable to find named zero bond " + zero.getName());

    if (zeroIter->second.sliceDim != slice.getDim())
        slice.allocDim(zeroIter->second.sliceDim);

    // ?? put logic to find zero slice associated with current zero.  For the 
    // meantime, we only have 1 zero slice so it is assumed it is correctly
    // reset, dev'd etc. - must assumes contiguous coupons for now

    // ?? use this for now until sorting out use of TreeSlice internally
    if (Hyb3_CopySlice(slice.getValuePtr(),
                       zeroIter->second.zeroSlice,
                       zeroIter->second.sliceDim,
                       step,
                       &mTreeData) != SUCCESS)
    {
        throw ModelException(routine, IrConverter::slog.pop());
    }
    slice.treeStep = range->treeStep;
}

//------------------------------------------------------------------------------
// Insert directly into the tree the zero bank dates and support fields 
// that have been calculated by the product when the tree is running in
// zeroBank mode
//------------------------------------------------------------------------------
void Hyb3::insertNamedZeroBank(NamedZeroBank& zeroBank)
{
    const char* routine = "Hyb3::insertNamedZeroBank";
    if (zeroBankMode == ZeroBankMode::CLAIMBANK)
        throw ModelException(routine, "Tree must be in ZEROBANK mode to insert " 
                             "client calculated zero bank");

    // simply insert the zero -map ensures a unique name
    mZeroBank[zeroBank.name] = zeroBank; 
}


///////////////////////////////////
// roll back tree from t+t to t
// maintaining internal zero banks
///////////////////////////////////
void Hyb3::update(int t, FDProduct::UpdateType type)
{
    int status = FAILURE;
    static char const* routine  = "Hyb3::update";

    int T = mTreeData.NbTP;

    // set current range
    range->limits.bot1 = mTreeData.Bottom1[t];
    range->limits.top1 = mTreeData.Top1[t];
    range->limits.bot2 = mTreeData.Bottom2[t];
    range->limits.top2 = mTreeData.Top2[t];
    range->limits.bot3 = mTreeData.Bottom3[t];
    range->limits.top3 = mTreeData.Top3[t];
    range->treeStep = t;

    // Update tree to current time slice
    if (Hyb3_Lattice(&mDevData,
                     t,
                     T,
                     mMktVolData,
                     &mTreeData) != SUCCESS)
    {
        goto RETURN;
    }

    // potential dangerous variable: this needs to be the same than 
    // the value in the range
    mTpIdx = t;

    // update the zerobanks
    if (zeroBankMode == ZeroBankMode::ZEROBANK)
    {
        map<string, NamedZeroBank>::iterator zb;
        // loop through each stored named zero bank
        for (zb = mZeroBank.begin(); zb != mZeroBank.end(); zb++)
        {
            int indexCurve = zb->second.discCurve;

            if (Hyb3_Zero_Bank(&zb->second.zeroBank[0],
                               &zb->second.zeroCritDates[0],
                               &zb->second.nbCurrentZeros,
                               zb->second.nbZeros,
                               FALSE,           // Do not reset now, just DEV
                               getDate(t).toIrDate(),
                               t,
                               T,
                               indexCurve,   // discount on the index curve
                               zb->second.devMode,
                               &mDevData,
                               &mTreeData) != SUCCESS)
            {
               DR_Error("Error DEVing zero bank name %s", (zb->first).c_str());
               goto RETURN;
            }
        }

        // update any of the zero bonds maintained by the tree
        map<string, ZeroBond>::iterator zero;
        for (zero = mZeroBond.begin(); zero != mZeroBond.end(); zero++)
        {
            int devMode = zero->second.devMode;
            int discCurve = mTreeData.CvDisc[0];

            if (Hyb3_Zero_t(zero->second.zeroSlice,
                            FALSE,
                            t,
                            T,
                            discCurve,  // Discount on COF curve
                            devMode,
                            &mDevData,
                            &mTreeData) != SUCCESS)
            {
                DR_Error("Error DEVing zero discounter name %s", (zero->first).c_str());
                goto RETURN;
            }
        }
    }
    else
    {
        throw ModelException(routine, "Claim bank mode not yet supported");
    }

    status = SUCCESS;

    RETURN:

    if (status != SUCCESS)
        throw ModelException(routine, IrConverter::slog.pop());
}


void Hyb3::updateFinalize(int t, FDProduct::UpdateType type)
{
    const char* routine = "Hyb3::updateFinalize";
    int status = FAILURE;

    DateTime currentDate = getDate(t);

    if (!mTreeBuilt)
        throw ModelException(routine, "Not able to update tree as the tree has "
                                      "not been built");

    // if zero bond reset date, reset zero bond (ie. set slice to 1)
    if (zeroBankMode == ZeroBankMode::ZEROBANK)
    {
        int i;
        map<string, NamedZeroBank>::iterator zb;
        map<string, ZeroBond>::iterator zero;

        // loop through each stored named zero bank
        for (zb = mZeroBank.begin(); zb != mZeroBank.end(); zb++)
        {
            bool zeroMatDate = false;
            int T = mTreeData.NbTP;
            long currentDateL = currentDate.toIrDate();

            // check if today is a zero bank maturity date  
            // ?? could use find but would be superceeded by event driver to be done
            for (i = 0; i < zb->second.zeroMatDateList.size(); i++)
            {
                if (currentDate == zb->second.zeroMatDateList[i])
                {
                    zeroMatDate = true;
                    break;
                }
            }

            if (zeroMatDate)
            {
                if (zb->second.zeroBank.empty())
                    throw ModelException(routine, "Named zero bank " + zb->second.name +
                                                  " zero bank slice collection is not defined");

                if (Hyb3_Zero_Bank(&zb->second.zeroBank[0],
                                   &zb->second.zeroCritDates[0],
                                   &zb->second.nbCurrentZeros,
                                   zb->second.nbZeros,
                                   1,       // reset zerobank slice
                                   currentDateL,
                                   t,
                                   T,
                                   mTreeData.CvDisc[zb->second.currencyEnum],
                                   zb->second.devMode,
                                   &mDevData,
                                   &mTreeData) != SUCCESS)
                {
                    goto RETURN;
                }
            }
        }

        for (zero = mZeroBond.begin(); zero != mZeroBond.end(); zero++)
        {
            bool zeroMatDate = false;
            int T = mTreeData.NbTP;

            // check if today is a zero bond maturity date  
            for (i = 0; i < zero->second.matDates.size(); i++)
            {
                if (currentDate == zero->second.matDates[i])
                {
                    zeroMatDate = true;
                    break;
                }
            }

            if (zeroMatDate)
            {
                if (Hyb3_Zero_t(zero->second.zeroSlice,
                                1,  // force zero bond reset
                                t,
                                T,
                                mTreeData.CvDisc[zero->second.currencyEnum],
                                zero->second.devMode,
                                &mDevData,
                                &mTreeData) != SUCCESS)
                {
                    goto RETURN;
                }
            }
        }
    }


    status = SUCCESS;

    RETURN:

    if (status != SUCCESS)
        throw ModelException(routine, IrConverter::slog.pop());

}


///////////////////////////////////
// roll forward tree from t to t+1
// maintaining internal zero banks
///////////////////////////////////
void Hyb3::updateStatePrice(int t, FDProduct::UpdateType type)
{
    int status = FAILURE;
    static char const* routine  = "Hyb3::update";
    TreeSliceRatesSP  tempSlice;

    int statePriceDim = -1; 

    //int T = mTreeData.NbTP;

    // set current range
    range->limits.bot1 = mTreeData.Bottom1[t];
    range->limits.top1 = mTreeData.Top1[t];
    range->limits.bot2 = mTreeData.Bottom2[t];
    range->limits.top2 = mTreeData.Top2[t];
    range->limits.bot3 = mTreeData.Bottom3[t];
    range->limits.top3 = mTreeData.Top3[t];
    range->treeStep = t;


    // allocate the slices on the first step
    if (type == FDProduct::FWD_0) {

        // Set the dimension of the state price (dependent on tree-type)
        switch (mTreeData.TreeType)
        {
        case TTYPE_FX2IR:
        case TTYPE_EQD2IR:
        case TTYPE_EQF2IR:
        case TTYPE_EQC2IR:
        case TTYPE_2IR2F1D:
            statePriceDim = 3;
            break;
        case TTYPE_EQ1IR:
        case TTYPE_2IR:
            statePriceDim = 2;
            break;
        default:
            throw ModelException("forward induction not implemented for "
                "model without EQ or FX dependency (i.e. CUPS2D etc)!");
            break;
        }

        statePriceSlice = DYNAMIC_POINTER_CAST<TreeSliceRates>(createSlice( ));
        statePriceSlice->name = "statePriceSlice";
        tempStatePriceSlice = DYNAMIC_POINTER_CAST<TreeSliceRates>(createSlice( ));
        tempStatePriceSlice->name = "tempStatePriceSlice";
        statePriceSlice->allocDim(statePriceDim);
        tempStatePriceSlice->allocDim(statePriceDim);
    }

    // swap the pointer around: we need to compute t+1 in order
    // to populate the FXRate correctly but need statePrice at t 
    // for express DEV calculation (might need to change in Lattice
    // in order to allow for forward induction as well - to do later)
    tempSlice           = tempStatePriceSlice;
    tempStatePriceSlice = statePriceSlice;
    statePriceSlice     = tempSlice;

    // check if we have the correct time
    statePriceSlice->testTreeStep() ; // should be this, but is protected currently...

    if (Hyb3_UpdateStatePrices( t,
                mMktVolData,
                &mTreeData,
                &mDevData,                
                statePriceSlice->getValuePtr(),
                tempStatePriceSlice->getValuePtr() ) != SUCCESS )
    {
        DR_Error("Error updating the State Prices at time %d", t);
        goto RETURN;
    }

    // update the TimePoint Index
    mTpIdx = t;

    // update the next tree slices time step
    tempStatePriceSlice->treeStep = t+1;

    // StatePriceUpdate- function() : we need a state containing that info

    status = SUCCESS;

    RETURN:

    if (status != SUCCESS)
        throw ModelException(routine, IrConverter::slog.pop());
}


/**  getStatePrices from the underlying model                               *
  *  supports the expressDEV (expressBACK) and the forward inductive scheme */
const TreeSlice & Hyb3::getStatePriceSlice(const int t) const
{
    // if we go forward in the tree, simply return the current state-price slice
    if (getInductionType() == FDModel::InductionType::FWD )
    {
        return *statePriceSlice; 
    }
    else if (getInductionType() == FDModel::InductionType::EXPRESSBACK )
    {
        map<int, TreeSliceSP>::const_iterator currentStatePriceSP = statePrices.find( t );
        if (currentStatePriceSP != statePrices.end() )
            return (*currentStatePriceSP->second);
    }
    else
        throw ModelException("state prices are only stored for forward and express backward"
                "induction schemes" );

    // return value to avoid warning of control path not returning values
    return *statePriceSlice;
};



//------------------------------------------------------------------------------
// Perform one backward induction step by moving the slice from t+1 to t - no discounting
//------------------------------------------------------------------------------
void Hyb3::sliceEv(TreeSlice& treeSlice) const
{
    // nothing to do at the back of the tree 
    if (mTpIdx == mTreeData.NbTP)
    {
        return;
    }

//    Hyb3_Ev (slice.getValuePtr(), mTpIdx, mTreeData.NbTP, &mDevData, &mTreeData);
}

//------------------------------------------------------------------------------
void Hyb3::sliceExpand(TreeSlice& treeSlice, const string& curveName) const {
    try {
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

        int currIdx = getCurrIdx(curveName);

        if (currIdx == FOREIGN) {
            slice.expand( mMktVolData[FOREIGN].NbFactor );
            return;
        }
        if (slice.getDim() == 3) 
            return;

        slice.expand( mMktVolData[FOREIGN].NbFactor + mMktVolData[DOMESTIC].NbFactor );
    }
    catch (exception& e) {
        throw ModelException(e, "Hyb3::sliceExpand");
    }
}

int Hyb3::getCurrIdx(const string& curveName) const
{
    if(curveName.empty())
        throw ModelException(__FUNCTION__, "curve id not supplied");
    
    // loop through all curves, find a name match and double check that there
    // are no other duplicate curve names in the other currency

    if (mTreeData.TreeType == TTYPE_1IR)
        return 0;

    int i, j, index=-1;
    bool match=false;
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 3; j++) {
            if (curveName == mCurveName[i][j]) {
                if (match == true)
                    throw ModelException(__FUNCTION__, "Curve name " + curveName +
                                         "appears as both a domestic and foreign currency "
                                         "in hyb3 tree registered curve names");                
                index = i;
                match=true;
                break;  // if match found check the other currency as well
            }
        }
    }

    if (index == -1)
        throw ModelException(__FUNCTION__, "Unable to find curveName " +
                             curveName + 
                             " in list of hyb3 registered yield curve names");
    
    return index;
}


// ??? currently only works for 2IR mode
// ??? AK: the whole tree-slice setup might need rethinking with view on flexible discounting!
void Hyb3::sliceDev(TreeSlice& treeSlice, int curveIdx) const
{
    TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);
    try {
        string sliceName = slice.name;  // to help with debugging

        // nothing to do at the back of the tree 
        if (mTpIdx == mTreeData.NbTP)
            return;

        // curveIdx == crvIdx * 10 + currIdx, see getCurveIdx
        int currIdx = curveIdx % 10;
        int crvIdx = curveIdx / 10;

        int devMode;
        int sliceDim = slice.getDim();

        if (currIdx == FOREIGN) {
            if ( (sliceDim > 1) && (mTreeData.TreeType == TTYPE_2IR2F1D ) )
                throw ModelException("Hyb3 does not currently support discounting a " +
                                    Format::toString(sliceDim) + " dimension Foreign currency "
                                    "denominated slice unless it's in Hyb2+1 mode.  This usually "
                                    "arises when a foreign denominated component is fixing/reseting "
                                    "on a domestic index/rate (which Hyb3 considers a reverse cups).");

            if (mTreeData.TreeType == TTYPE_2IR2F1D )
                devMode = DISC_2D_1IR2F_NOCUPS;
            else
                devMode = DISC_1D_NOCUPS;
        }
        else if (sliceDim == 3) {
            if (mTreeData.TreeType == TTYPE_2IR2F1D )
                devMode = DISC_3D_2IR2F1D_CUPS;
            else
                devMode = DISC_3D_CUPS;
        }
        else {
            devMode = DISC_2D_CUPS;
        }

        if (Hyb3_Dev(slice.getValuePtr(),
                    mTpIdx, 
                    mTreeData.NbTP, 
                    crvIdx, 
                    devMode, 
                    (/*const_cast*/HYB3_DEV_DATA*)&mDevData, 
                    (/*const_cast*/HYB3_TREE_DATA*)&mTreeData) != SUCCESS)
        {
            throw ModelException(IrConverter::slog.pop());
        }
    }
    catch (exception& e) {
        throw ModelException(e, "Hyb3::sliceExpand, slice=\""+slice.name+"\"");
    }
}

int Hyb3::getCurveIdx( const string & curveName ) const
{
    if( curveName.empty())
        return -1;

    int currIdx = getCurrIdx( curveName );
    int crvIdx = getCrvIdx( curveName, currIdx );
    return crvIdx * 10 + currIdx;
}

//
//
int Hyb3::getCrvIdx(const string& curveName, int currIdx) const
{
    try {
        if(curveName.empty())
            throw ModelException("curveName not supplied");
        
        for (int i = 0; i < 3; i++)
        {
            if (curveName == mCurveName[currIdx][i])
                return i;
        }
        throw ModelException("Unable to find curveName " +
                            curveName + 
                            " in list of hyb3 registered curve names");
    }
    catch (exception& e) {
        throw ModelException(e, "Hyb3::getCrvIdx");
    }
}

// ??? needs to work for all modes, or override 
/*
int Hyb3::getCurrIdx(const string& curveName) const
{
    try {
        if(curveName.empty())
            throw ModelException("curve id not supplied");
        
        // loop through all curves, find a name match and double check that there
        // are no other duplicate curve names in the other currency

        if (mTreeData.TreeType == TTYPE_1IR)
            return 0;

        int i, j, index=-1;
        bool match = false;
        for (i = 0; i < 2; i++)
        {
            for (j = 0; j < 3; j++)
            {
                if (curveName == mCurveName[i][j])
                {
                    if (match)
                        throw ModelException("Curve name " + curveName +
                                            "appears as both a domestic and foreign currency "
                                            "in hyb3 tree registered curve names");                
                    index = i;
                    match = true;
                    break;  // if match found check the other currency as well
                }
            }
        }

        if (!match)
            throw ModelException("Unable to find curveName " +
                                curveName + 
                                " in list of hyb3 registered yield curve names");
        
        return index;
    }
    catch (exception& e) {
        throw ModelException(e, "Hyb3::getCurrIdx");
    }
}*/

//------------------------------------------------------------------------------
// Get date on the timeline pointed by the index
//------------------------------------------------------------------------------
DateTime Hyb3::getDate(int dateIdx) const
{
    if (dateIdx < 0 || dateIdx > mTreeData.NbTP) {
        throw ModelException(__FUNCTION__, 
        "Invalid time point index (" + Format::toString(dateIdx) 
        + "). Must be in the range 0..." + Format::toString(mTreeData.NbTP));
    }
    return DateTime::fromIrDate(mTreeData.TPDate[dateIdx]);
}

DateTimeArray Hyb3::getDates() const
{
    DateTimeArray dates( mTreeData.NbTP + 1 );
    for( int i = 0; i <= mTreeData.NbTP; ++i )
        dates[ i ] = DateTime::fromIrDate( mTreeData.TPDate[ i ] );
    return dates;
}

/** get last step (total num of steps) on time line */
int Hyb3::getLastStep() const
{
    return mTreeData.NbTP;
}

// ?? depends on tree mode - only working for 2IR
DateTime Hyb3::getCurveValueDate(string curveName) const
{
    try {
        if (curveName.empty())
        {
            return DateTime::fromIrDate(mFxData.ValueDate);
        }
        else  // return value date of interest rate
        {
            int currIdx = getCurrIdx(curveName);
            int index = getCrvIdx(curveName, currIdx);
            return DateTime::fromIrDate(mTCurves[currIdx][index].ValueDate);
        }
        return DateTime();
    }
    catch (exception& e) {
        throw ModelException(e, "Hyb3Tree::getValueDate");
    }
}

//------------------------------------------------------------------------------
// Print
//------------------------------------------------------------------------------
void Hyb3::print()
{
    static char const* routine  = "Hyb3Tree::Print";
    int i,j,k;

    map<string, ZeroBond>::iterator zero;
    map<string, NamedZeroBank>::iterator zb;

    char   CurrencyString[2][MAXBUFF];
    strcpy(CurrencyString[0], "foreign");
    strcpy(CurrencyString[1], "domestic");

#if 0
    {
        FILE *stream = fopen("b.txt", "wb");

        fprintf(stream, "\nMKTVOL_DATA structure for Foreign IR\n\n");
        MktVol_PrintStructure(stream, &mMktVolData[FOREIGN]);
        fprintf(stream, "\nMKTVOL_DATA structure for Domestic IR\n\n");
        MktVol_PrintStructure(stream, &mMktVolData[DOMESTIC]);

        fprintf(stream, "\nTREE_DATA structure\n\n");
        PrintTree(stream, mTreeData, 0);
        fclose(stream);

    }
#endif

    if (treeDataDebugFile.empty())
        return;

    FILE *stream = fopen(treeDataDebugFile.c_str(), "w");
    if (stream == NULL)
    {
        throw ModelException(routine, "Unable to open tree debug info file  (" + 
                             treeDataDebugFile + ") for writing");
    }
    int count=0;

    fprintf(stream, "Following critical dates registered with the engine:\n\n");
    for (i = 0; i < (int)sortedCritDates.size(); i++)
    {
        fprintf(stream, "%d   %ld\n", i, sortedCritDates[i]);
    }

    fprintf(stream, "Foreign Currency %s\n", mCurrency[FOREIGN].c_str());
    fprintf(stream, "Domestic Currency %s\n\n", mCurrency[DOMESTIC].c_str());
    fprintf(stream, "Nb  ForeignCurves    DomesticCurves\n");
    for (i = 0; i < 3; i++)
    {
        fprintf(stream, "%d  %s  %s\n", 
                i,
                mCurveName[FOREIGN][i].c_str(),
                mCurveName[DOMESTIC][i].c_str());
    }

    fprintf(stream, "\nFollowing zero banks registered with the engine:\n\n");
    for (zb = mZeroBank.begin(); zb != mZeroBank.end(); zb++)
    {
        fprintf(stream, "%d, name: %s, internalName: %s, nbZeroInBankUse %d, zeroBankDim %d, currency %s\n",
                count,
                zb->first.c_str(),
                zb->second.name.c_str(),
                zb->second.nbZeros,
                zb->second.zeroBankDim,
                zb->second.currency.c_str());

        // print out the zero dates associated with the zerobank
        for (i = 0; i < (int)zb->second.zeroMatDateList.size(); i++)
        {
            fprintf(stream, "%d zeroDate: %ld\n",
                    i,
                    zb->second.zeroMatDateList[i].toIrDate());
        }
        count++;
    }

    count = 0;
    fprintf(stream, "\n\nFollowing zero bond/discounting slices registered with the engine:\n\n");
    for (zero = mZeroBond.begin(); zero != mZeroBond.end(); zero++)
    {
        fprintf(stream, "%d, name: %s, nbZeros(ie. non-overlapping) %d, sliceDim %d, currency %s\n",
                count,
                zero->first.c_str(),
                zero->second.startDates.size(),
                zero->second.sliceDim,
                zero->second.currency.c_str());

        fprintf(stream, "   Non-overlapping zero start/end dates stored in zeroIndexSpec\n");
        for (i = 0; i < zero->second.startDates.size(); i++)
        {
            fprintf(stream, "   %d  %ld  %ld\n", i, 
                    zero->second.startDates[i].toIrDate(),
                    zero->second.matDates[i].toIrDate());
        }
        count++;
    }

    fprintf (stream, "\n\n");

    // Slice sizes used for allocation
    fprintf(stream, "  W1   HW1    W2   HW2    W3   HW3\n");
    fprintf(stream, "%5d %5d %5d %5d %5d %5d\n\n\n",
                    mTreeData.Width[0], mTreeData.HalfWidth[0],
                    mTreeData.Width[1], mTreeData.HalfWidth[1],
                    mTreeData.Width[2], mTreeData.HalfWidth[2]);

    if (mTreeData.TreeType == TTYPE_2IR)
    {
        for (j=0; j<2; j++)
        {
            fprintf (stream,"\nTIMELINE INFORMATION (%s IR)\n",CurrencyString[j]);

            fprintf (stream,
                     "Node    Date       Days  Max  Forward   "
                     "Zero0   Discount0   Zero1   Discount1   "
                     "Zero2   Discount2    IrMidNode    SpotVol \n");
        
            for (i = 0; i <= mTreeData.NbTP; i++)
            {
                double Forward = 100. * mTreeData.FwdRate[j][0][i]/mTreeData.Length[i];

                double days = Daysact (mTreeData.TPDate[0], mTreeData.TPDate[i]);

                fprintf(stream,
                        "[%3d] \t%8ld  %5.0f  %3d  %7.4f  %7.4f   %8.6f  "
                        "%7.4f   %8.6f  %7.4f   %8.6f    %9.6f    %6.2f \n",
                        i,
                        mTreeData.TPDate[i],
                        days,
                        j ? mTreeData.Top2[i][0] : mTreeData.Top1[i],
                        Forward,
                        mTreeData.ZeroRate[j][0][i] * 100.,
                        mTreeData.ZeroCoupon[j][0][i],
                        mTreeData.ZeroRate[j][1][i] * 100.,
                        mTreeData.ZeroCoupon[j][1][i],
                        mTreeData.ZeroRate[j][2][i] * 100.,
                        mTreeData.ZeroCoupon[j][2][i],
                        exp(mTreeData.IrZCenter[j][i]),
                        mTreeData.SpotVol[j][i] * 100.);

            }
        }
    }


    fprintf (stream,"\nDomestic Rate Timeline Information (%s IR)\n",mCurrency[DOMESTIC].c_str());

    //fprintf (stream, "Node      Date     days  forward   zero0  discount0  zero1  discount1  zero2  discount2   Drift spotVol\n");
    const char labelsStr[]="Node      Date     days  forward    zero0  discount0    zero1  discount1    zero2  discount2   Drift  spotVol \n";
    const char formatStr[]="[%3d] \t%8ld  %5.0f  %7.4f  %7.4f   %8.6f  %7.4f   %8.6f  %7.4f   %8.6f  %9.6f  %6.2f\n";
    fprintf (stream, labelsStr);
        
    for (i = 0; i <= mTreeData.NbTP; i++)
    {
        double forward = 100. * mTreeData.FwdRate[DOMESTIC][0][i]/mTreeData.Length[i];
        double days = Daysact (mTreeData.TPDate[0], mTreeData.TPDate[i]);

        fprintf(stream, formatStr,
                i,
                mTreeData.TPDate[i],
                days,
                forward, // mTreeData.FwdRate[DOMESTIC][0][i] * 100.   ,  // Forward,
                mTreeData.ZeroRate[DOMESTIC][0][i] * 100.,
                mTreeData.ZeroCoupon[DOMESTIC][0][i],
                mTreeData.ZeroRate[DOMESTIC][1][i] * 100.,
                mTreeData.ZeroCoupon[DOMESTIC][1][i],
                mTreeData.ZeroRate[DOMESTIC][2][i] * 100.,
                mTreeData.ZeroCoupon[DOMESTIC][2][i],
                exp(mTreeData.IrZCenter[DOMESTIC][i]),
                mTreeData.SpotVol[DOMESTIC][i] * 100.);
    }

    // print foreign rate data if tree mode
    if (mTreeData.TreeType == TTYPE_2IR ||
        mTreeData.TreeType == TTYPE_FX2IR)
    {
        fprintf (stream,"\n\nForeign Rate Timeline Information (%s IR)\n",mCurrency[FOREIGN].c_str());

        //fprintf (stream, "Node      Date     days  forward   zero0  discount0  zero1  discount1  zero2  discount2   Drift spotVol\n");
        fprintf (stream, "Node      Date     days  forward   zero0  discount0  zero1  discount1  zero2  discount2   Drift  spotVol\n");
        
        for (i = 0; i <= mTreeData.NbTP; i++)
        {
            double forward = 100. * mTreeData.FwdRate[FOREIGN][0][i]/mTreeData.Length[i];
            double days = Daysact (mTreeData.TPDate[0], mTreeData.TPDate[i]);

            fprintf(stream, formatStr,
                    i,
                    mTreeData.TPDate[i],
                    days,
                    forward, // mTreeData.FwdRate[FOREIGN][0][i] * 100.   ,  // Forward,
                    mTreeData.ZeroRate[FOREIGN][0][i] * 100.,
                    mTreeData.ZeroCoupon[FOREIGN][0][i],
                    mTreeData.ZeroRate[FOREIGN][1][i] * 100.,
                    mTreeData.ZeroCoupon[FOREIGN][1][i],
                    mTreeData.ZeroRate[FOREIGN][2][i] * 100.,
                    mTreeData.ZeroCoupon[FOREIGN][2][i],
                    exp(mTreeData.IrZCenter[FOREIGN][i]),
                    mTreeData.SpotVol[FOREIGN][i] * 100.);
        }
    }

    if(mTreeData.TreeType == TTYPE_FX2IR)
    {   /* If we have FX info, print it */
        fprintf (stream,"\nTIMELINE INFORMATION (FX)\n");
        
        fprintf (stream, "Node   Date        CWidth     Forward    FxVol   "
                         " Spotvol   Rho's (Irf-Ird, Irf-FX, Ird-FX)\n");
                                        
        for (i = 0; i <= mTreeData.NbTP; i++)
        {
            fprintf (stream, "[%3d]  %8ld  %6d  %12.6f    "
                             "%5.2f    %5.2f       %4.2f       %4.2f      %4.2f\n",
                             i, 
                             mTreeData.TPDate[i],
                             mTreeData.Top3[i][0][0],
                             mTreeData.FwdFx[i],
                             mTreeData.FxVol[i] * 100,
                             mTreeData.SpotFxVol[i] * 100.,
                             mTreeData.Rho[0][i],
                             mTreeData.Rho[1][i],
                             mTreeData.Rho[2][i]);

        }  /* for i */
    }

    // Print out input zero curves
    for (j=0; j<2; j++)
    {
        fprintf (stream, "\n\n");
        fprintf (stream, 
                 "%s Currency Curves (Index, COF, Risk)\n\n",CurrencyString[j]);
        for (k = 0; k < 3; k++)
        {
                
            fprintf (stream, "Maturity      Zero       Discount \n");

            for (i = 0; i < (mTCurves[j][k]).NbZero; i++)
            {

                double days = Daysact((mTCurves[j][k]).ValueDate, 
                                      (mTCurves[j][k]).ZeroDate[i]);
                // Discount factor up to the current date
                double discount = ::pow(1. + (mTCurves[j][k]).Zero[i], -days/365.);

                fprintf (stream,
                         "%ld   %9.6f     %8.6f \n",
                         (mTCurves[j][k]).ZeroDate[i],
                         (mTCurves[j][k]).Zero[i] * 100.,
                         discount);
            }
            fprintf (stream, "\n");
        }
    }

    // print out raw hyb3 structures
/*
    fprintf(stream, "\nMKTVOL_DATA structure for Foreign IR\n\n");
    MktVol_PrintStructure(stream, &mMktVolData[FOREIGN]);
    fprintf(stream, "\nMKTVOL_DATA structure for Domestic IR\n\n");
    MktVol_PrintStructure(stream, &mMktVolData[DOMESTIC]);

    fprintf(stream, "\nTREE_DATA structure\n\n");
    PrintTree(stream, mTreeData, 0);
*/
    fclose(stream);

}

//------------------------------------------------------------------------------
// Get today's date
//------------------------------------------------------------------------------
DateTime Hyb3::getToday() const
{
    return mToday;
}

//------------------------------------------------------------------------------
// Get current date
//------------------------------------------------------------------------------
/*DateTime Hyb3::getCurrentDate() const
{
    return DateTime::fromIrDate(mTreeData.TPDate[mTpIdx]);
}*/

//------------------------------------------------------------------------------
// Null-initialize tree structure
//------------------------------------------------------------------------------
void Hyb3::init()
{

    // initialise structures
    MktVol_Init(&mMktVolData[FOR]);
    MktVol_Init(&mMktVolData[DOM]);
    Hyb3_Tree_Init(&mTreeData);
    Hyb3_Dev_Init(&mDevData);

    mTpIdx = -1;
    mSmoothFlag = false;

    mTreeBuilt = false;

    // set the dimension of the individual IR assets
    mMktVolData[DOMESTIC].NbFactor = nbFactorsDom;
    mMktVolData[FOREIGN].NbFactor = nbFactorsFgn;

}

//------------------------------------------------------------------------------
// Delete parts of the model that have been created by build.
//------------------------------------------------------------------------------
Hyb3::~Hyb3() {
    clear();
}

//------------------------------------------------------------------------------
// Delete parts of the model that have been created by build.
//------------------------------------------------------------------------------
void Hyb3::clear()
{
    if (!mTreeBuilt)
        return;
    mTreeBuilt = false;    

    ASSERT (zeroBankMode == ZeroBankMode::ZEROBANK);

    map<string, ZeroBond>::iterator zero;
    map<string, NamedZeroBank>::iterator zb;

    // ??? delete the zero bank
    for (zb = mZeroBank.begin(); zb != mZeroBank.end(); zb++)
    {
        vector<double*> &bank = zb->second.zeroBank;

        for (int i = bank.size()-1; i >=0; --i) {
            Hyb3_Free_Slice(bank[i], &mTreeData, zb->second.zeroBankDim);
            bank[i] = NULL;
        }
        bank.clear();
    }

    for (zero = mZeroBond.begin(); zero != mZeroBond.end(); zero++)
    {
        Hyb3_Free_Slice (zero->second.zeroSlice, &mTreeData, zero->second.sliceDim);
        zero->second.zeroSlice = NULL;
    }

    Hyb3_Dev_Free (&mDevData, &mTreeData);
    Hyb3_Tree_Free(&mTreeData);

    init();
}

//------------------------------------------------------
//
//------------------------------------------------------
void Hyb3::insertIRIndex(const IndexSpecIR& rate, DateTime date) {}

//----------------------------------------------------
//
//----------------------------------------------------
void Hyb3::getIRIndex(const IndexSpecIR& rateSpec, DateTime currentDate, TreeSlice& treeSlice)
{
    int status = FAILURE;
    const char* routine = "Hyb3::getIRIndex";
    TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

    if (zeroBankMode == ZeroBankMode::ZEROBANK)
    {
        string name = rateSpec.getName();
        char dcc = IrConverter::dccTo035A(*rateSpec.dcc);
        char freq = IrConverter::toEslFreq(rateSpec.frequency.get());

        // get namedZeroBank
        map<string, NamedZeroBank>::iterator zb = mZeroBank.find(name);
        if (zb == mZeroBank.end())
            throw ModelException(routine, "Unable to find named zero bank " + name);
        
        if (zb->second.zeroBankDim != slice.getDim())
            slice.allocDim(zb->second.zeroBankDim);

        if (Hyb3_Par_Yield_t(slice.getValuePtr(),
                             zb->second.zeroBankDim,
                             zb->second.nbCurrentZeros,
                             &zb->second.zeroBank[0],
                             &zb->second.zeroCritDates[0],
                             1,  // definite reset
                             currentDate.toIrDate(),
                             currentDate.toIrDate(),  // spot starting
                             rateSpec.tenor->toMonths(),
                             dcc, 
                             freq,
                             0,  // float coupon spread = 0
                             mTpIdx,
                             &mTreeData) != SUCCESS)
        {
            goto RETURN;
        }
        slice.treeStep = range->treeStep;
    }
    else  // CLAIMBANK
    {
        DR_Error("Claim bank not yet supported");
        goto RETURN;
    }

    status = SUCCESS;

    RETURN:

    if (status != SUCCESS)
        throw ModelException(routine, IrConverter::slog.pop());

}

void Hyb3::getFXIndex(TreeSlice& treeSlice)
{
    const char* routine = "Hyb3::getFXIndex";
    TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

    if (slice.getDim() != 3)
        slice.allocDim(3);

    if (mTpIdx != range->treeStep )
        throw ModelException(routine , "tree step index mTpIdx is different from range "
                        "index: probably forgot to update one or the other");


    // for now simply copy the FX Spot slice into the provided slice
    if (Hyb3_CopySlice(slice.getValuePtr(),
                       mDevData.FxSpot,
                       3,       // FX always 3D slice
                       mTpIdx,  // current step
                       &mTreeData) != SUCCESS)
    {
        throw ModelException(routine, IrConverter::slog.pop());
    }
    slice.treeStep = range->treeStep;
}

void Hyb3::populateIRParams(MKTVOL_DATA* mktvol, HYB3_TREE_DATA* tree, RateTree::IRModelParams* mp)
{
    mktvol->NbFactor = mp->nbFactors;
    mktvol->QLeft = 1.0 - mp->QLeft;
    mktvol->QRight = 1.0 - mp->QRight;
    mktvol->FwdShift = mp->FwdShift;
    mktvol->CetNbIter = mp->CetNbIter;

    for (int i = 0; i < 3; ++i)
    {
        mktvol->Alpha[i] = mp->FactorVol[i];
        mktvol->Beta[i] = mp->FactorMR[i];
        mktvol->Rho[i] = mp->FactorCorr[i];
    }


    // Check validity of input
    if (Hyb3_Param_Check(mp->nbFactors,
                         mktvol,
                         tree) != SUCCESS)              
    {        
        throw ModelException("Hyb3::PopulateIRParams", IrConverter::slog.pop());
    }
}


/*********************************** Hyb3 - FX Mode *******************************/

void Hyb3FX::load(CClassSP& clazz){

    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(Hyb3FX, clazz);
    SUPERCLASS(Hyb3);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(ppy, "minimum number of tree nodes(points)/year");
    FIELD(maxStdDeviations, "number of standard deviations to trim the tree at");
    FIELD(forIRParams, "Collection of foreign interest rate model parameters");
    FIELD(domIRParams, "Collection of domestic interest rate model parameters");
    FIELD(FXSmileParams, "Collection of named FX smile parameters - assumes no FX smile if not supplied");
    FIELD_MAKE_OPTIONAL(FXSmileParams);
    FIELD(FXIRCorrelations, "Collection of user defined correlation overrides")
    FIELD_MAKE_OPTIONAL(FXIRCorrelations);
    FIELD(treeDataDebugFile, "Print tree debug data to named output file");
    FIELD_MAKE_OPTIONAL(treeDataDebugFile);
    // though these are members of hyb3 base class, export them in the derived class to allow
    // each mode of hyb3 to decide if it wants to export them or not
    FIELD(zeroInterpStyle, "zero curve interpolation style. Defaults to "
                                  "LINEAR if not supplied");
    FIELD_MAKE_OPTIONAL(zeroInterpStyle);
    FIELD(zeroBankMode, "internal tree zero bank mode - defaults to"
                               "to ZEROBANK mode if not supplied");
    FIELD_MAKE_OPTIONAL(zeroBankMode);
    FIELD(FXSmileDriftAdjustment, "Turn on or off FX Smile CUPS drift adjustment calculation, "
                                         "FALSE/off by default");
    FIELD_MAKE_OPTIONAL(FXSmileDriftAdjustment);
    FIELD(FXSpotVolOverride,"Overrides market comp/spot vols for FX - sets constant spot vol for FX "
                            "diffusion (units = decimal NOT percentage)");
    FIELD_MAKE_OPTIONAL(FXSpotVolOverride);

    FIELD(nbFactorsDom,"Number of interest rate factors for domestic rate (default = 1)");
    FIELD_MAKE_OPTIONAL(nbFactorsDom);
    FIELD(nbFactorsFgn,"Number of interest rate factors for foreign rate (default = 1)");
    FIELD_MAKE_OPTIONAL(nbFactorsFgn);

    // these fields are transient to allow tweaking to copy 
    // the model after ::getMarket() has been called 
    FIELD(fxIndexSpecSupplied, "")
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(fxIndexSpecSupplied);  // ??? 
    FIELD(corrIR, "");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrIR);
    FIELD(corrForIRFX, "");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrForIRFX);
    FIELD(corrDomIRFX, "");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrDomIRFX);

}

CClassConstSP const Hyb3FX::TYPE = CClass::registerClassLoadMethod(
    "Hyb3FX", typeid(Hyb3FX), Hyb3FX::load);

// for class loading
bool Hyb3FXLoad(void) { return (Hyb3FX::TYPE != 0); }


void Hyb3FX::getMarket(const MarketData*  market, IInstrumentCollectionSP instruments)
{
    static const string routine = "Hyb3FX::getMarket";
try {

    mToday = market->GetReferenceDate();

    forIRParams->getData(this, market);
    domIRParams->getData(this, market);

    // in FX mode, must provide curve to discount
    if (forIRParams->curveToDiscount.isEmpty())
        throw ModelException(routine, "must supply foreign curve to discount for Hyb3FX");
    if (domIRParams->curveToDiscount.isEmpty())
        throw ModelException(routine, "must supply domestic curve to discount for Hyb3FX");

    // top level instrument must supply discount factor curve in order to define the
    // pricing/PnL currency of the trade.  This is retrieved through a complicated tangle
    // of code that eventually populates a marketDataFetcher by KComponent implementing the
    // IInstrumentCollection interface function:
    // virtual string discountYieldCurveName() const = 0
    // which returns the name of the discount curve defined in the instrument
    string discYCName = getDomesticYCName();
    if (discYCName.empty())
        throw ModelException("instrument must supply a discount Curve");

    string modelDiscYCName = domIRParams->curveToDiscount->getName();
    if (discYCName != modelDiscYCName)
        throw ModelException("Name of curveToDiscount field set in model " + modelDiscYCName +
                             " must be same as defined on the instrument " + discYCName);

    discYC = domIRParams->curveToDiscount.getSP();
    
    // Get Engine
    if (!engineTable.isEmpty())
        engineTable.getData(this, market);

    // recurse all the potential components in the instrument and look for the type
    // IMarketFactor::TYPE, which represent marketData type fields (yield curves, 
    // index Specs etc) and store them for retrieval later on
    class RetrieveFactors : public ObjectIteration::IAction
    {
        FDModel * model;
        const MarketData * market;
        IMarketFactorArray & factors;

    public:
        RetrieveFactors( FDModel * model, const MarketData * market, IMarketFactorArray & factors ) :
            model( model ), market( market ), factors( factors ) { factors.clear(); }

        virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj)
        {
            IMarketFactor * factor = dynamic_cast< IMarketFactor * >( state.getObject().get() );

            string name = factor->getName();
            string type = factor->getClass()->getName();

            if (FXAsset::TYPE->isInstance(factor)) return true;

            int i;
            for( i = 0; i < factors.size(); ++i )
            {
                if( name == factors[ i ]->getName() && type == factors[ i ]->getClass()->getName() )
                     break;
            }
            if( i >= factors.size() && model->acceptFactor( factor ) )
            {
                factors.push_back( IMarketFactorSP::attachToRef( factor ) );
                return false;
            }
            else
                return true;
        }
    };
    RetrieveFactors retrieveFactors( this, market, factors );
    ObjectIteration iteration( IMarketFactor::TYPE );
    iteration.recurse( retrieveFactors, instruments );

    fxIndexSpecSupplied = false;

    // recurse all the potential components in the instrument and look for the type
    // IndexSpec::TYPE, and store them for later retrieval.  Store the IndexSpecFX
    // types in a separate array to simplify logic of determining which mode the tree
    // should be configured for
    class RetrieveIndexSpecs : public ObjectIteration::IAction
    {
        IndexSpecArray & indexSpecs;
        IMarketFactorArray & fxFactors;
        bool & fxIndexSpecSupplied;

    public:
        RetrieveIndexSpecs(IndexSpecArray & indexSpecs, IMarketFactorArray & fxFactors, 
            bool & fxIndexSpecSupplied) :
            indexSpecs( indexSpecs ), fxFactors( fxFactors ), fxIndexSpecSupplied(fxIndexSpecSupplied)
            { indexSpecs.clear(); fxFactors.clear(); }

        virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj)
        {
            int i;
            IndexSpec * indexSpec = dynamic_cast< IndexSpec * >( state.getObject().get() );

            string name = indexSpec->getName();
            string type = indexSpec->getClass()->getName();

            for( i = 0; i < indexSpecs.size(); ++i )
            {
                if( name == indexSpecs[ i ]->getName() && type == indexSpecs[ i ]->getClass()->getName() )
                    break;
            }
            if( i >= indexSpecs.size() )
            {
                indexSpecs.push_back( IndexSpecSP::attachToRef( indexSpec ) );
                // if FX, get the FX market factor at the same time
                IndexSpecFX * fxIndexSpec = dynamic_cast< IndexSpecFX * >( indexSpec );
                if (fxIndexSpec && !fxIndexSpec->sameCurve) {
                    fxIndexSpecSupplied = true;

                    IMarketFactor *factor = fxIndexSpec->getFactor().get();
                    if (factor && !factor->getName().empty()) {
                        // test for duplicates
                        int j;
                        for( j = 0; j < fxFactors.size(); ++j )
                        {
                            if( name == fxFactors[ j ]->getName())
                            break;
                        }
                        if( j >= fxFactors.size() )
                            fxFactors.push_back(fxIndexSpec->getFactor());
                    }
                }
                return false;
            }
            else
                return true;
        }
    };
    RetrieveIndexSpecs retrieveIS ( indexSpecs, fxFactors, fxIndexSpecSupplied);
    ObjectIteration retrieveISIter( IndexSpec::TYPE );
    retrieveISIter.recurse( retrieveIS, instruments );

    string domCurrency = discYC->getCcy();

    // FX between foreign and domestic discount curves must be supplied
    // ??? ignore for now any others, but check later if supplied that FX spot is the same
    string fxName = market->getFXName(forIRParams->curveToDiscount->getName(),
                                      domIRParams->curveToDiscount->getName());
    mFXName = fxName;
    MarketObjectSP fxObj = market->GetData(fxName, FXAsset::TYPE);

    FXAsset *fx = dynamic_cast<FXAsset*>(fxObj.get());
    if (!fx)
        throw ModelException(routine, "Internal error - Unable to cast an FXAsset::TYPE MarketFactor "
                             "to a FXAsset* pointer, which should never happen");

    fx->getMarket(this, market);

    IMarketFactor* fxFact = dynamic_cast<IMarketFactor*>(fxObj.get());
    string name = fxFact->getName();
    string type = fxFact->getClass()->getName();
    int j;
    for( j = 0; j < fxFactors.size(); ++j )
    {
        if( name == fxFactors[ j ]->getName())
        break;
    }
    if( j >= fxFactors.size() )
        fxFactors.push_back( IMarketFactorSP::attachToRef( fxFact ) );


    // retrieve correlations, into structure of forIR-domIR, forIR-FX, domIR-FX
    double nan = sqrt(-1.0);
    corrIR = nan;
    corrForIRFX = nan;
    corrDomIRFX = nan;

    // forIR/domIR
    // check all yield curves of different currencies and model curves - must be same correlation
    // value or at least one provided
    string factor1, factor2;

    factor1 = forIRParams->curveToDiscount->getName();
    factor2 = domIRParams->curveToDiscount->getName();
    if (market->hasCorrelationData(factor1, factor2))
    {
        string corrName = market->getCorrelationName(factor1, factor2);

        CorrelationBaseSP corrBase = getCorrelation(corrName,
                                                    forIRParams->curveToDiscount->getClass(),
                                                    domIRParams->curveToDiscount->getClass(),
                                                    Correlation::TYPE,
                                                    market);
        corrIR = CorrelationSP::dynamicCast(corrBase)->getCorrelation(); 
    }

    factor1 = forIRParams->curveToDiffuse->getName();
    factor2 = domIRParams->curveToDiffuse->getName();
    if (market->hasCorrelationData(factor1, factor2))
    {
        string corrName = market->getCorrelationName(factor1, factor2);

        CorrelationBaseSP corrBase = getCorrelation(corrName,
                                                    forIRParams->curveToDiscount->getClass(),
                                                    domIRParams->curveToDiscount->getClass(),
                                                    Correlation::TYPE,
                                                    market);
        double result = CorrelationSP::dynamicCast( corrBase )->getCorrelation();
        if (!Maths::finite(corrIR))
            corrIR = result;
        else
        {
            if (!Maths::equals(corrIR, result))
                throw ModelException(routine, "Correlation defined for curveToDiffuse factors " +
                                     factor1 + " and " + factor2 + "(value = " + 
                                     Format::toString(result) +
                                     ") does not equal previously calculated correlation for " +
                                     "model discount curves");
        }
    }
    
    // ??? more to do here to check all other correlations of yield curve to yield curve
    if (!Maths::finite(corrIR))
        throw ModelException(routine, "??? tmp error - correlation between forIR and domIR not found");


    // foreignIR/FX correlation
    // ??? just check curve to diffuse and curve to discount FXs for now, but have to check all
    // of them later

    string foreignIR = forIRParams->curveToDiscount->getName();
    string domesticIR = domIRParams->curveToDiscount->getName();

    if (market->hasFXData(foreignIR, domesticIR))
    {
        string fxName = market->getFXName(foreignIR, domesticIR);
        if (market->hasCorrelationData(foreignIR, fxName))
        {
            string corrName = market->getCorrelationName(foreignIR, fxName);

            CorrelationBaseSP corrBase = getCorrelation(corrName,
                                                        forIRParams->curveToDiscount->getClass(),
                                                        FXAsset::TYPE,
                                                        Correlation::TYPE,
                                                        market);
            corrForIRFX = CorrelationSP::dynamicCast(corrBase)->getCorrelation();
        }
    }

    if (!Maths::finite(corrForIRFX))
        throw ModelException(routine, "Correlation not found in market for foreignIR ("+foreignIR+"), FX ("+fxName+") on discount curve");

    // domesticIR/FX correlation
    // ??? just check curve to diffuse and curve to discount FXs for now, but have to check all
    // of them later

    foreignIR = forIRParams->curveToDiscount->getName();
    domesticIR = domIRParams->curveToDiscount->getName();

    if (market->hasFXData(foreignIR, domesticIR))
    {
        string fxName = market->getFXName(foreignIR, domesticIR);
        if (market->hasCorrelationData(domesticIR, fxName))
        {
            string corrName = market->getCorrelationName(domesticIR, fxName);

            CorrelationBaseSP corrBase = getCorrelation(corrName,
                                                        domIRParams->curveToDiscount->getClass(),
                                                        FXAsset::TYPE,
                                                        Correlation::TYPE,
                                                        market);
            corrDomIRFX = CorrelationSP::dynamicCast(corrBase)->getCorrelation();
        }
    }

    if (!Maths::finite(corrForIRFX))
        throw ModelException(routine, "Correlation not found in market for domesticIR ("+domesticIR+"), FX ("+fxName+") on discount curve");

}
catch (exception& e) {
    throw ModelException(e, routine);
}
}

inline void Hyb3FX::getForeignDiffusionTCurve(T_CURVE &diffusionCurve) const
{
    IrConverter::to_T_CURVE(diffusionCurve, forIRParams->curveToDiffuse.get(),
                            false);
}

inline void Hyb3FX::getDomesticDiffusionTCurve(T_CURVE &diffusionCurve) const
{
    IrConverter::to_T_CURVE(diffusionCurve, domIRParams->curveToDiffuse.get(),
                            false);
}

inline IRVolRawSP Hyb3FX::getForeignIRVolRaw() const
{
    return getIRVolRaw(forIRParams);
}

inline IRVolRawSP Hyb3FX::getDomesticIRVolRaw() const
{
    return getIRVolRaw(domIRParams);
}

void Hyb3FX::populateTreeIRParams(Currency curr)
{
    CcyIRParamsSP modelIRParams;
    IRVolRawSP volSP;

    if (curr == FOREIGN)
    {
        modelIRParams = forIRParams;
        volSP = getForeignIRVolRaw();
    }
    else
    {
        modelIRParams = domIRParams;
        volSP = getDomesticIRVolRaw();
    }

    IRVolSelector volSelector(volSP, mTCurves[curr][0],
                              modelIRParams->volCalibIndex,
                              cetSmoothing, modelIRParams);
    volSelector.getMktVolData(mMktVolData[curr]);

    // ??? temp "intermediate" structure to reuse logic from hyb3 for now
    RateTree::IRModelParams mp;

    // retrieve model parameters from market cache
    IRExoticParamTable* pSmileTable = modelIRParams->smileTable.get();
    IRExoticParamTable* pModelTable = modelIRParams->modelTable.get();
    if (pSmileTable || pModelTable)
    {
        // Use new MarketTable object
        IrConverter::to_MODELPARAMETERS_DATA(mp, 
                                            modelIRParams->smileSet, 
                                            modelIRParams->modelSet,
                                            engineSet,
                                            *(engineTable.get()),
                                            *pSmileTable,
                                            *pModelTable);
    }
    else
    {
        // Use deprecated IRCalib object 
        IrConverter::to_MODELPARAMETERS_DATA(mp,
                    modelIRParams->smileSet, 
                    modelIRParams->modelSet,
                    *(modelIRParams->irCalib.get()));
    }

    // populate tree structures from the model overwrite parameters or
    // underlying model parameters in the market data */
    populateIRParams(&mMktVolData[curr], &mTreeData, &mp);
}

// Implementation of IRVegaPointwise::ISensitivePoints for smart vega tweaking
// and volatility exposure reporting.
IRGridPointAbsArraySP Hyb3FX::getSensitiveIRVolPoints(
    OutputNameConstSP  outputName,
    const CInstrument* inst) const
{
    IRGridPointAbsArraySP volExposuresSP(new IRGridPointAbsArray);

    IRVolRawSP volSP;
    T_CURVE diffusionCurve;

    // Note that checks below for matching outputName are just for efficiency
    // (e.g. to avoid unnecessary construction of diffusion T_CURVEs, vol
    // selection, and exposure processing.

    volSP = getForeignIRVolRaw();
    if (outputName->equals(volSP->getBaseVol()->getName()) ||
        outputName->equals(volSP->getSwapVol()->getName()))
    {
        getForeignDiffusionTCurve(diffusionCurve);
        IRVolSelector volSelector(volSP, diffusionCurve,
                                  forIRParams->volCalibIndex,
                                  cetSmoothing, forIRParams);
        volSelector.getVolExposures(volExposuresSP, outputName);
    }

    volSP = getDomesticIRVolRaw();
    if (outputName->equals(volSP->getBaseVol()->getName()) ||
        outputName->equals(volSP->getSwapVol()->getName()))
    {
        getDomesticDiffusionTCurve(diffusionCurve);
        IRVolSelector volSelector(volSP, diffusionCurve,
                                  domIRParams->volCalibIndex,
                                  cetSmoothing, domIRParams);
        volSelector.getVolExposures(volExposuresSP, outputName);
    }
   
    return volExposuresSP;
}

void Hyb3FX::initModel()
{
    if (FXSmileDriftAdjustment)
        mTreeData.FxMomentMatching = TRUE;
    else
        mTreeData.FxMomentMatching = FALSE;

    Hyb3::initModel();

}

// ??? rethink this later on - maybe mix in with initModel
void Hyb3FX::retrieveFactor(){
    static const string method = "Hyb3::retrieveFactor";
    try{
        int i,j;

        if (!discYC)
            throw ModelException("FDModel::discYC field is not populated - must contain the top product's "
                                 "discount curve. Internal error in FDModel/product code");

        mCurrency[DOMESTIC] = discYC->getCcy();

        if (domIRParams->curveToDiffuse->getCcy() != mCurrency[DOMESTIC])
            throw ModelException(method, "Currency of supplied domestic curve to diffuse " + 
                                 domIRParams->curveToDiffuse->getCcy() +
                                 " must be the same as domestic currency defined by the product " +
                                 mCurrency[DOMESTIC]);

        if (domIRParams->curveToDiscount->getCcy() != mCurrency[DOMESTIC])
            throw ModelException(method, "Currency of supplied domestic curve to diffuse " + 
                                 domIRParams->curveToDiffuse->getCcy() +
                                 " must be the same as domestic currency defined by the product " +
                                 mCurrency[DOMESTIC]);

        mCurrency[FOREIGN] = forIRParams->curveToDiffuse->getCcy();

        if (forIRParams->curveToDiscount->getCcy() != mCurrency[FOREIGN])
            throw ModelException(method, "Currency of supplied foreign curve to discount " + 
                                 forIRParams->curveToDiscount->getCcy() +
                                 " must be the same as foreign currency defined by curveToDiffuse" +
                                 mCurrency[FOREIGN]);

        Hyb3::init();
        mTreeData.Ppy = ppy;
        mTreeData.NbSigmaMax = maxStdDeviations;

        for (i=0; i<2; ++i) 
        {
            mTreeData.CvDisc[i] = 1;
            mTreeData.CvDiff[i] = 0;
            mTreeData.CvIdx1[i] = 1;
            mTreeData.CvIdx2[i] = 2;
        }

        // by convention set 0=diff curve, 1=discount curve, 2=other
        mCurveName[FOREIGN][0] = forIRParams->curveToDiffuse->getName();
        IrConverter::to_T_CURVE(mTCurves[FOREIGN][0], forIRParams->curveToDiffuse.get(), false);
        mCurveName[FOREIGN][1] = forIRParams->curveToDiscount->getName();
        IrConverter::to_T_CURVE(mTCurves[FOREIGN][1], forIRParams->curveToDiscount.get(), false);

        mCurveName[DOMESTIC][0] = domIRParams->curveToDiffuse->getName();
        IrConverter::to_T_CURVE(mTCurves[DOMESTIC][0], domIRParams->curveToDiffuse.get(), false);
        mCurveName[DOMESTIC][1] = domIRParams->curveToDiscount->getName();
        IrConverter::to_T_CURVE(mTCurves[DOMESTIC][1], domIRParams->curveToDiscount.get(), false);

        // for now only allow 1 FX to be registered for now ???
        if (fxFactors.size() != 1)
            throw ModelException("Currently hyb3 only supports one FX factor - need to expand this later");
 
        // check types and number of factors defined by instrument
        int nbYC=0;

        bool thirdCurve[2];
        thirdCurve[0] = false;
        thirdCurve[1] = false;
        for (i = 0; i < factors.size(); i++)
        {
            IMarketFactorConstSP &mf = factors[i];
            const YieldCurve* yc = dynamic_cast<const YieldCurve*>(mf.get());

            if (yc)
            {
                // for now must be either Foreign or Domestic
                if (yc->getCcy() != mCurrency[FOREIGN] &&
                    yc->getCcy() != mCurrency[DOMESTIC])
                    throw ModelException("Yield curve defined in instrument of currency " +
                                     yc->getCcy() +
                                     ". Hyb3 is configured as foreign = " + mCurrency[DOMESTIC] +
                                     " and domestic = " + mCurrency[FOREIGN]);

                for (j = 0; j < 2; j++)
                {
                    if (yc->getName() == mCurveName[j][0] ||
                        yc->getName() == mCurveName[j][1])
                        continue;
                    else
                    {
                        if (yc->getCcy() == mCurrency[j])
                        {
                            // define third currency or error if already defined
                            if (thirdCurve[j])
                                throw ModelException("Hyb3 engine supports maximum of three yield curves/currency "
                                                     "and supplied yield curve/factor " + yc->getName() +
                                                     "is a fourth for currency " + yc->getCcy());
                            mCurveName[j][2] = yc->getName();
                            IrConverter::to_T_CURVE(mTCurves[j][2], yc, false);
                            thirdCurve[j] = true;
                        }
                    }
                }
                nbYC++;
            }
            else
            {
                throw ModelException("Hyb3FX only supports yield curve market factors - not type " +
                                     mf->getClass()->getName());
            }
        }
        // if not supplied, set the third curve to equal the second/discount
        for (j = 0; j < 2; j++)
        {
            if (thirdCurve[j] == false)
            {
                mCurveName[j][2] = mCurveName[j][1];
                mTCurves[j][2] = mTCurves[j][1];
            }
        }

        // only FX and IR indexSpecs supported
        for (i = 0; i < indexSpecs.size(); i++)
        {
            string type = indexSpecs[i]->getClass()->getName();
            if (type != "IndexSpecIR" && 
                type != "IndexSpecFX")
            {
                throw ModelException(method, "Index spec type " + type +
                                     " is not supported by Hyb3FX model - must be either FX or IR");
            }
        }

        const FXAsset* fx = dynamic_cast<FXAsset*>(fxFactors[0].get());
        if (!fx)
            throw ModelException(method, "Internal error - FX market factor can not be cast to FXAsset");

        /*
        // ??? doesn't work as riskCcy has not got market - still just yield curve wrappers
        // dom ycs in fx assset must be domestic discount ycs
        string fxDomYCName = fx->getBaseCcy()->getName();
        if (fxDomYCName != domIRParams->curveToDiscount->getName())
            throw ModelException("Domestic yield curve " + fxDomYCName + " defined for FX Asset " + 
                                 fx->getName() +
                                 " by convention must be the model's domestic curveToDiscount yield curve " +
                                 domIRParams->curveToDiscount->getName());

        string fxForYCName = fx->getRiskCcy()->getName();
        if (fxForYCName != forIRParams->curveToDiscount->getName())
            throw ModelException("Foreign yield curve " + fxForYCName + " defined for FX Asset " + 
                                 fx->getName() +
                                 " by convention must be the model's foreign curveToDiscount yield curve " +
                                 forIRParams->curveToDiscount->getName());
        */
        bool fxMode = false;

        if (fxIndexSpecSupplied == true)
        {
            fxMode = true;
        }
        else 
        {
            vector<bool> factorInd(factors.size());
            for (i = 0; i < factors.size(); i++)
                factorInd[i] = false;
            for (i = 0; i < indexSpecs.size(); i++)
            {
                IndexSpecSP is = indexSpecs[i];
                const IndexSpecIR* ir = dynamic_cast<const IndexSpecIR*>(is.get());
                if (ir)
                {
                    int j;
                    // find factor associated with indexSpec and mark with indicator
                    for (j = 0; j < factors.size(); j++)
                    {
                       string factorName = factors[j]->getName();
                       string irFactorName = ir->getFactor()->getName();
                       if (factorName == irFactorName)
                       {
                            factorInd[j] = true;  // match found - factor derived from indexSpec
                            break;
                       }
                    }
                }
            }
            // now for all IR factors not marked as true (ie., derived from index spec), check currency
            // and if any is foreign, must be in 3D mode
            for (i = 0; i < factors.size(); i++)
            {
                IMarketFactorConstSP &mf = factors[i];

                const YieldCurve* yc = dynamic_cast<const YieldCurve*>(mf.get());
                if (yc)
                {
                    if (factorInd[i] == false)
                    {
                        string ccy=yc->getCcy();  // help debugging
                        if (ccy != discYC->getCcy())
                        {
                            fxMode = true;
                            break;
                        }
                    }
                }
            }
        }
        
        if (fxMode)
            mTreeData.TreeType = TTYPE_FX2IR;
        else
        {
            if ( mMktVolData[FOREIGN].NbFactor ==1 ) 
                mTreeData.TreeType = TTYPE_2IR;
            else if (  mMktVolData[FOREIGN].NbFactor ==2 )
                mTreeData.TreeType = TTYPE_2IR2F1D;
            else
                throw ModelException(method, "Internal error - more than two factors per IR asset not supported in Hyb3");
        }

        FXVOLATILITY_DATA fxVol;
        CORRELATION_DATA correlations;

        // only support a maximum of three distinct ycs, so check if more were supplied
        // ??? not sure if factors has repeats or these are sorted out before hand
        /*if (factors.size() > 3)
            throw ModelException("A maximum of three distinct " + discYC->getCcy() +
                                 " currency curves are supported by hyb3 - " +
                                 Format::toString(factors.size()) + " supplied"); */

        // populate foreign dimension with domesic values for
        populateTreeIRParams(FOREIGN);
        populateTreeIRParams(DOMESTIC); 

        strcpy(mTreeData.Index[FOREIGN], forIRParams->volCalibIndex.c_str());
        strcpy(mTreeData.Index[DOMESTIC], domIRParams->volCalibIndex.c_str());

        IrConverter::to_FXVOLATILITY_DATA(fxVol, fx);
        fxVol.ValueDate = mToday.toIrDate();

        // Correlations - if in single currency mode, set correlations to 0
        if (!FXIRCorrelations.empty())
            throw ModelException("Model correlation family overrides not yet supported, "
                                 "but will be afuture requirement to be able to easily "
                                 "override correlations");

        correlations.CorrIR = corrIR;
        correlations.CorrForIRFX = corrForIRFX;
        correlations.CorrDomIRFX = corrDomIRFX;

        // AK: do we need to overwrite corr-information here??
        // currently hard coded information, needs to done properly
        mMktVolData[FOREIGN].CorrSwapSt    = Nxtmth( mToday.toIrDate(), 120L, 1L);
        mMktVolData[FOREIGN].CorrSwapMat   = Nxtmth(mMktVolData[FOREIGN].CorrSwapSt,  120L, 1L);
        mMktVolData[FOREIGN].CorrSwapDCC   = '3';
        mMktVolData[FOREIGN].CorrSwapFreq  = 'A';
        mMktVolData[DOMESTIC].CorrSwapSt    = Nxtmth( mToday.toIrDate(), 120L, 1L);
        mMktVolData[DOMESTIC].CorrSwapMat   = Nxtmth(mMktVolData[DOMESTIC].CorrSwapSt,  120L, 1L);
        mMktVolData[DOMESTIC].CorrSwapDCC   = '3';
        mMktVolData[DOMESTIC].CorrSwapFreq  = 'A';


        // ??? for now OwriteCorrel strings are set to nill, but may be used when model
        // supplies correlation overwrites
        char OwriteCorrel[3][MAXBUFF];
        strcpy(OwriteCorrel[0], "nil");
        strcpy(OwriteCorrel[1], "nil");
        strcpy(OwriteCorrel[2], "nil");

        if (HYB3CorrelOverwrites(&mFxData,
                                 OwriteCorrel,
                                 &correlations) != SUCCESS)
        {
            throw ModelException("HYB3CorrelOverwrites falied: "+IrConverter::slog.pop());
        }

        if (!FXSmileParams.empty())
            throw ModelException("Model defined FXSmile Parameter Family name is not yet "
                                 "supported - model reads FXSmile params directly from "
                                 "the SRMFXVol market structure");
        
        FXSMILE_DATA dummyFXSmile;
        IrConverter::to_FXSMILE_DATA(dummyFXSmile, fx);

        if (dummyFXSmile.NbSmileDates == 0)
            throw ModelException("At least one FX smile date and a1,a2,a3 data set must "
                                 "be provided.  To turn smile off set a1=a2=a3 = 0.0. "
                                 "Check that smileExpiry, smileA1... in "
                                 "FXAsset::SRMFX::Vol are not empty.");

        char OwriteFxSpot[MAXBUFF];
        char NbFXSmileOWS[MAXBUFF];
        char FXSmileParamOWS[MAXNBDATE][MAXBUFF];

        if (FXSpotVolOverride.get()) {  //  if supplied
            double value = FXSpotVolOverride->doubleValue();
            sprintf(OwriteFxSpot, "%lf", value * 100);
        }
        else
            strcpy(OwriteFxSpot, "nil");

        strcpy(NbFXSmileOWS, "nil");  // read smile params from dummyFXSmile structure

        // ??? reuse this for now, but rewrite later on
        if (HYB3FXSmileOverwrites(&mFxData,
                                  OwriteFxSpot,
                                  &fxVol,
                                  &dummyFXSmile,
                                  NbFXSmileOWS,
                                  FXSmileParamOWS) != SUCCESS)
        {
            throw ModelException("HYB3FXSmileOverwrites falied: "+IrConverter::slog.pop());
        }

    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

// ??? what to do
void Hyb3FX::initTreeData(void) 
{
try {

    }
    catch (exception& e) {
        throw ModelException(e, "Hyb3::initTreeData");
    }
}

// assuming call option
static double deltaToStrike(double delta, double vol, double fwdFX)
{
    double d1 = Normal_InvH(delta);
    double dummy = - d1 * vol + 0.5 * vol * vol;
    double impliedStrike = fwdFX * exp(dummy);

    return impliedStrike;
}

// ??? figure out how to turn this on/off but preferably not defining a full
// outputRequest for each one - maybe have a "debug info" output request
// that switches on all of this debug packet info
void Hyb3FX::recordOutput(Control* ctrl, Results* results) const
{
    int i, j;

    if (!ctrl->requestsOutput(OutputRequest::DBG))
        return; // no debug info requested

    // store tree critical dates
    DateTimeArraySP treeDates(new DateTimeArray(mTreeData.NbTP+1));
    for (i = 0; i <= mTreeData.NbTP; i++)
        (*treeDates)[i] = DateTime::fromIrDate(mTreeData.TPDate[i]);

    results->storeGreek(treeDates,
                        Results::DEBUG_PACKET, 
                        OutputNameConstSP(new OutputName("TREE_DATES")));

    // store number of days offsets from start of tree for each critical date
    IntArraySP daysOffsets(new IntArray(mTreeData.NbTP+1));
    for (i = 0; i <= mTreeData.NbTP; i++)
    {
        double days = Daysact(mTreeData.TPDate[0], mTreeData.TPDate[i]);
        (*daysOffsets)[i] = (int)days;
    }

    results->storeGreek(daysOffsets,
                        Results::DEBUG_PACKET, 
                        OutputNameConstSP(new OutputName("TREE_DATE_OFFSETS")));

    // store FOREIGN ir rates
    CDoubleMatrixSP fgnZeros(new CDoubleMatrix(3, mTreeData.NbTP+1));
    for (i = 0; i < 3; i++)
        for (j = 0; j <= mTreeData.NbTP; j++)
            (*fgnZeros)[i][j] = mTreeData.ZeroRate[0][i][j];

    results->storeGreek(fgnZeros,
                        Results::DEBUG_PACKET, 
                        OutputNameConstSP(new OutputName("FOREIGN_IR_ZEROS")));

    // store DOMESTIC ir rates
    CDoubleMatrixSP domZeros(new CDoubleMatrix(3, mTreeData.NbTP+1));
    for (i = 0; i < 3; i++)
        for (j = 0; j <= mTreeData.NbTP; j++)
            (*domZeros)[i][j] = mTreeData.ZeroRate[1][i][j];

    results->storeGreek(domZeros,
                        Results::DEBUG_PACKET, 
                        OutputNameConstSP(new OutputName("DOMESTIC_IR_ZEROS")));

    // store foreign IR spotVols
    DoubleArraySP forIRSpotVol(new DoubleArray(mTreeData.NbTP+1));
    for (i = 0; i <= mTreeData.NbTP; i++)
        (*forIRSpotVol)[i] = mTreeData.SpotVol[FOREIGN][i];

    results->storeGreek(forIRSpotVol,
                        Results::DEBUG_PACKET, 
                        OutputNameConstSP(new OutputName("FOREIGN_IR_SPOT_VOLS")));

    // store domestic IR spotVols
    DoubleArraySP domIRSpotVol(new DoubleArray(mTreeData.NbTP+1));
    for (i = 0; i <= mTreeData.NbTP; i++)
        (*domIRSpotVol)[i] = mTreeData.SpotVol[DOMESTIC][i];

    results->storeGreek(domIRSpotVol,
                        Results::DEBUG_PACKET, 
                        OutputNameConstSP(new OutputName("DOMESTIC_IR_SPOT_VOLS")));

    // store FX spot vols
    DoubleArraySP fxSpotVols(new DoubleArray(mTreeData.NbTP+1));
    for (i = 0; i <= mTreeData.NbTP; i++)
        (*fxSpotVols)[i] = mTreeData.SpotFxVol[i];

    results->storeGreek(fxSpotVols,
                        Results::DEBUG_PACKET, 
                        OutputNameConstSP(new OutputName("FX_SPOT_VOLS")));

    // store FX composite vols
    DoubleArraySP fxCompVols(new DoubleArray(mTreeData.NbTP+1));
    for (i = 0; i <= mTreeData.NbTP; i++)
        (*fxCompVols)[i] = mTreeData.FxVol[i];

    results->storeGreek(fxCompVols,
                        Results::DEBUG_PACKET, 
                        OutputNameConstSP(new OutputName("FX_COMP_VOLS")));

    // store forward FX 
    DoubleArraySP fxFwds(new DoubleArray(mTreeData.NbTP+1));
    for (i = 0; i <= mTreeData.NbTP; i++)
        (*fxFwds)[i] = mTreeData.FwdFx[i];

    results->storeGreek(fxFwds,
                        Results::DEBUG_PACKET, 
                        OutputNameConstSP(new OutputName("FX_FWDS")));

    // Store local volatility function for various strikes
    CDoubleMatrixSP smileVols(new CDoubleMatrix(5, mTreeData.NbTP+1));
    DoubleArraySP smileDeltas(new DoubleArray(5));
    (*smileDeltas)[0] = 0.1;
    (*smileDeltas)[1] = 0.25;
    (*smileDeltas)[2] = 0.5;
    (*smileDeltas)[3] = 0.75;
    (*smileDeltas)[4] = 0.9;

    for (i = 0; i <= mTreeData.NbTP;i++)
    {
        double time = (*daysOffsets)[i] / 365.0;
        double vol = sqrt(time)* mTreeData.FxVol[i]; // use comp Vol for estiamting delta
        double fxSpotVol = mTreeData.SpotFxVol[i];   // spot Vol for local vol function

        for (j = 0; j < 5; j++)
        {
            double delta = (*smileDeltas)[j];
            double strikeFwdFx = deltaToStrike(delta, vol, mTreeData.FwdFx[i]);
            double moneyness = strikeFwdFx/mTreeData.FwdFx[i];
            double smile;

            if (Hyb3_Gfunc(&smile,
                           moneyness,
                           mTreeData.A1C[i],
                           mTreeData.A2C[i],
                           mTreeData.A3C[i]) != SUCCESS)
               throw ModelException("Hyb3FX::recordOutput", IrConverter::slog.pop());

            smile *= fxSpotVol/moneyness;

            (*smileVols)[j][i] = smile;
        }
    }

    results->storeGreek(smileVols,
                        Results::DEBUG_PACKET, 
                        OutputNameConstSP(new OutputName("FX_SMILE_VOLS")));

    results->storeGreek(smileDeltas,
                        Results::DEBUG_PACKET, 
                        OutputNameConstSP(new OutputName("FX_SMILE_DELTAS")));


/*
            OutputRequest* request = 
            control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            if (paydates && paydates->size() == 0) {
                results->storeNotApplicable(request); 
            } 
            else {
                DateTimeListSP datelist(new DateTimeList(paydates));
                results->storeRequestResult(request, datelist); 
            }
        }       
*/

    // 
}


/** Create a MarketDataFetcher which will be used for retrieving market data etc */
MarketDataFetcherSP Hyb3::createMDF() const {
    MarketDataFetcherSP mdf(new RateTreeMDF(1));
    MDFUtil::setUseCurrencyBasis(*mdf, true);
    return mdf;
}

Hyb3::Hyb3(const CClassConstSP &type) : RateTree(type), cetSmoothing(false),  
    mSmoothFlag(false), zeroBankMode(ZeroBankMode::ZEROBANK), zeroInterpStyle(ZeroInterpStyle::LINEAR), 
    mTreeBuilt(false), momentMatching(false),nbFactorsDom(1),nbFactorsFgn(1)
{
    mFxData.ValueDate = -1; // to prevent getValueDate() from returning crazy results
    memset(&mTreeData, 0, sizeof(mTreeData));
    memset(&mMktVolData, 0, sizeof(mMktVolData));
    memset(&mDevData, 0, sizeof(mDevData));
}


DRLIB_END_NAMESPACE
