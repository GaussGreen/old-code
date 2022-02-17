/************************************************************************
 * Module:	PenGuin
 * File:    vtlkionew.cxx
 * Function:    
 * Author:  Changhong He
 ************************************************************************/
#include "vtlkionew.h"
#include "vtlindicator.h"
#include "kutilios.h"
#include "kmrntree.h"


extern "C" {
#include "drltime.h"
#include "drlmath.h"
};


//===============================================================
//
//===============================================================


//---------------------------------------------------------------

KVPToolKnockIONew::KVPToolKnockIONew(SharedPointer<KVPKnockIONew> ins, KVTree &vt)
    : KVPToolAtom(vt)
{
static  char    routine[] = "KVPToolKnockIONew::KVPToolKnockIONew";
    int idx, n;

    TDate   today = vt.TPToday();

    try {
    // Check valid
    if (ins->mSettleDates.size() <= 0) {
        throw KFailure("%s: no dates found.\n", routine);
    }
        
    /* Trigger rate index curve name */
    mRateIdxName = ins->mRateIndex[0]->CplxRate().CurveName(); 


    // Only include events with observation date >= today 
    //
    ins->ValidEvents(today);

    mKnockIO = ins;

    setCurveName();

    // Add events in the tree.
    //
    n = mKnockIO->mObsDates.size();

    mModelObsDates.resize(n);

    for (idx=0; idx<n; idx++) 
    {
        vt.Insert(mKnockIO->mObsDates[idx]);
        vt.Insert(mKnockIO->mSettleDates[idx]);

        mModelObsDates[idx] = vt.Insert(*(mKnockIO->mRateIndex[idx]));
    }

    mValue = NULL;
    mValueUS  = NULL;
    mGetValue = NULL;
    mDefValue = NULL;

    mRebateValue   = NULL;
    mRebateValueUS = NULL;

    mXT0 = NULL;
    mXT1 = NULL;
    mXT2 = NULL;

    // Run statistics
    mRunStat = (debugLevel > 0 ? TRUE : FALSE);

    // Determine where is the ko index in the dependency list
    for (idx=0; idx<this->NumDep(); idx++) 
    {
        if (Dep(idx)->IsType("KVPToolIndicator"))
            mDepKoind = idx;
    }

    }
    catch (KFailure) 
    {
        throw KFailure("%s: failed on knock IO schedule:\n", routine);
    }
}

//---------------------------------------------------------------

KVPToolKnockIONew::~KVPToolKnockIONew()
{

    delete mValue;
    delete mValueUS;
    delete mGetValue;

    delete mRebateValue;
    delete mRebateValueUS;
    delete mDefValue;

    delete mXT0;
    delete mXT1;
    delete mXT2;

    for (KMap(TDate,KTSlice*)::iterator itEx = mExer.begin();
        itEx != mExer.end(); ++itEx)
    {
        delete (*itEx).second;
    }
    mExer.clear();

    for (KMap(TDate,KTSlice*)::iterator itRebate = mRebate.begin();
        itRebate != mRebate.end(); ++itRebate)
    {
        delete (*itRebate).second;
    }
    mRebate.clear();

    if(debugLevel)
        dppLog << GetName() << ": deleted." << endl;
}


//---------------------------------------------------------------
const String&
KVPToolKnockIONew::GetCurveName()
{
    if (mCurveName.empty()) {
      try 
      {
        int idx, idxMax, cIdxMax, cIdx;
        int discIdx;

        //
        // Identify largest geometry to allocate slice
        // The geometry is given by the curve index associated to
        // a name by the virtual tree.
        // THE CONVENTION IS THAT INCREASING INDICES REPRESENT
        // INCREASING GEOMETRIES.
        //
        cIdxMax = -100000;
        idxMax = -1;
        for(idx=0; idx<NumDep(); idx++) {   
            cIdx = mVTree->GetCurveIdx(Dep(idx)->GetCurveName());
            if (cIdx > cIdxMax) {
                idxMax = idx;
                cIdxMax = cIdx;
            }
        }

        // include discount curve
        discIdx = mVTree->GetCurveIdx(mKnockIO->GetDiscName());
        if (discIdx > cIdxMax)
            cIdxMax = discIdx;

        mCurveName = mVTree->GetCurveName(cIdxMax);
      }
      catch (KFailure) {
        throw KFailure("KVPToolKnockIO::GetCurveName: "
        "failed on knock IO `%s'.\n", mKnockIO->GetName());
      }
    }
    return mCurveName;
}




//---------------------------------------------------------------


void
KVPToolKnockIONew::Update()
{
static  char    routine[] = "KVPToolKnockIONew::Update";
    int idx, n = mKnockIO->mObsDates.size();
    TDate   curDate = mVTree->TPDateCurrent(),
            settleDate, obsDate, modelObsDate;
    int     idxD;

    int     depCurveIdx, rateCurveIdx;

    bool    isKnockOut = false;

    KCplxRateReset      *indexRate = NULL;

    KTSlice     *koidxslice1 = NULL, // knock index slice
                *koidxslice2 = NULL, // 1-ko index slice
                *tmpslice    = NULL; 
    KTSlice     *koidxslice_unsmooth = NULL, // no_smooth ko index slice
                *unsmoothValue       = NULL,
                *unsmoothRebate      = NULL,
                *lbAdj               = NULL,
                *hbAdj               = NULL;

    KTSliceNode tsNode;         // Slice node

    double      underlyingValue = 0e0;

    double      xt0,            // Probability
                xt1,            // expected ko time
                xt2;            // expected ko time^2 

const   String& discCurveName = mKnockIO->mDiscZcName;

    String      curveName = GetCurveName(); // max dim of underlying slice
    String      rateIndexCV = mRateIdxName; 

    SharedPointer<KVPToolIndicator> koindex;    // ko index


 try 
 {
    if (!NeedsUpdate()) return;

    double  exTime = mVTree->TPTimeCurrent();

    // 
    // The slice dimensions of mValue and mExer should be
    // the maximum of rate slice and underlying slice
    //
    depCurveIdx  = mVTree->GetCurveIdx(curveName);
    rateCurveIdx = mVTree->GetCurveIdx(mRateIdxName);

    if (rateCurveIdx > depCurveIdx)
        curveName = mRateIdxName;


    // Set Knock out flag
    //
    if (mKnockIO->mIOType == CRX_KNOCK_OUT)
        isKnockOut = true;

    //------------------------------------------------------
    //
    // (1) Perform dev on the value (both smoothed and unsmoothed)
    // and on the exercise bank.
    //
    //------------------------------------------------------
    if (mValue != NULL) 
    {
        mValue->Dev(discCurveName);
        if (mRunStat) {
            mXT0->Ev();
            mXT1->Ev();
            mXT2->Ev();
        }
    }


    if (mValueUS != NULL) 
    {
        mValueUS->Dev(discCurveName);
    }

    // Dev the exercise value
    //
    for(KMap(TDate,KTSlice*)::iterator itEx=mExer.begin();
         itEx != mExer.end();
         ++itEx) {
        ((*itEx).second)->Dev(discCurveName);
    }

    // Dev the rebate. Apply only to knock-out
    //
    if (isKnockOut)
    {
        for(KMap(TDate,KTSlice*)::iterator itRebate=mRebate.begin();
            itRebate != mRebate.end(); ++itRebate) 
        {
            ((*itRebate).second)->Dev(discCurveName);
        }


        //
        // Ensure mRebate slice exists
        //
        if (mRebateValue != NULL) 
        {
            mRebateValue->Dev(discCurveName);
        }

        if (mRebateValueUS != NULL) 
        {
            mRebateValueUS->Dev(discCurveName);
        }
    }

    //------------------------------------------------------
    //
    // (2) Add settlements
    //
    //------------------------------------------------------
    for (idx=0; idx<n; idx++) 
    {
        settleDate = mKnockIO->mSettleDates[idx];
    
        if (settleDate == curDate) 
        {
            if(debugLevel)
                dppLog << GetName() << format(": Processing settlement %d "
                " (settlement %s)",
                idx,
                DrlTDatePrint(NULL, mKnockIO->mSettleDates[idx]))
                << endl;

        //
        // Create new settlement value and exercise slice
        //
            mExer[settleDate] = new KTSlice(*mVTree, "Exercise", curveName);
            (*mExer[settleDate]) = 0e0;

        // Calculate value of the underlying
        // by summimg all dependencies.
        //
            for (idxD=0; idxD<this->NumDep(); idxD++) 
            {
                if (!(Dep(idxD)->IsType("KVPToolRate")) &&
                    !(Dep(idxD)->IsType("KVPToolCplxRate")) &&
                    !(Dep(idxD)->IsType("KVPToolIndicator")))
                    (*mExer[settleDate]) += Dep(idxD)->GetValue();
            }


        //
        // Create new rebate value and slice
        //
            if (isKnockOut)
            {
                mRebate[settleDate] 
                    = new KTSlice(*mVTree, "Rebate", curveName);
                (*mRebate[settleDate]) = mKnockIO->mRebates[idx];
            }

            if (debugLevel == DEBUG_LEVEL_KIO)
            {
                mVTree->TSliceSpecialOper(
                        *mExer[settleDate],
                        "TSPRT_MINMAX",
                        (ostream*) &dppLog);
                dppLog << endl;

                if (isKnockOut)
                {
                    dppLog << routine << ": Rebate" << endl;
                    mVTree->TSliceSpecialOper(
                         *mRebate[settleDate],
                         "TSPRT_MINMAX",
                         (ostream*) &dppLog);
                    dppLog << endl;
                }
            }

        }   // if settleDate
    }   // for idx


    //------------------------------------------------------
    //
    // (3) Perform knock-in exercises.
    //
    //------------------------------------------------------

    for (idx=0; idx<n; idx++) 
    {
        obsDate      = mKnockIO->mObsDates[idx];
        modelObsDate = mModelObsDates[idx];

        if (modelObsDate == curDate) 
        {
            if(debugLevel)
                dppLog << GetName() << format(": Processing observation "
                    "%d on date %s (settlement %s)", 
                    idx,
                    DrlTDatePrint(NULL, curDate), 
                    DrlTDatePrint(NULL, mKnockIO->mSettleDates[idx]))
                    << endl;

        //
        // Ensure value slice exists
        //
            if (mValue == NULL) 
            {
                mValue = new KTSlice(*mVTree, "mValue", curveName);
                *mValue = 0e0;
                if (mRunStat) {
                    mXT0 = new KTSlice(*mVTree, "mXT0", curveName);
                    *mXT0 = 0e0;
                    mXT1 = new KTSlice(*mVTree, "mXT1", curveName);
                    *mXT1 = 0e0;
                    mXT2 = new KTSlice(*mVTree, "mXT2", curveName);
                    *mXT2 = 0e0;
                }
            }

            if (mKnockIO->mSmooth == DOUBLE_SMOOTH && mValueUS == NULL) 
            {
                mValueUS = new KTSlice(*mVTree, "mValueUS", curveName);
                *mValueUS = 0e0;
            }

            if (isKnockOut)
            {
                if (mRebateValue == NULL) 
                {
                    mRebateValue 
                        = new KTSlice(*mVTree, "mRebateValue", curveName);
                    *mRebateValue = 0e0;
                }

                if (mRebateValueUS == NULL && 
                    mKnockIO->mSmooth == DOUBLE_SMOOTH)
                {
                    mRebateValueUS = 
                        new KTSlice(*mVTree, "mRebateValueUS", curveName);
                    *mRebateValueUS = 0e0;
                }
            }

            if (debugLevel == DEBUG_LEVEL_KIO)
            {
                dppLog << "Knock IO value slice BEFORE exercise: " 
                       << endl;
                dppLog << *mValue << endl;
                
                if (mKnockIO->mSmooth == DOUBLE_SMOOTH) 
                {
                    dppLog << *mValueUS << endl;
                }

                if (isKnockOut)
                    dppLog << *mRebateValue << endl;
            }


        //
        // Retrieve settlement value time slice
        //
            settleDate = mKnockIO->mSettleDates[idx];

            KMap(TDate,KTSlice*)::iterator itEx = mExer.find(settleDate);

            if (itEx == mExer.end()) 
            {
                throw KFailure("%s: Current observation date is %s. "
                       "Can't find exercise slice on "
                       "settlement date %s.\n",
                       routine, 
                       GtoFormatDate(curDate),
                       GtoFormatDate(settleDate));
            }

            KTSlice& exerTS = *((*itEx).second);


        //
        // Retrieve rebate value time slice
        //
            KTSlice *rebateTS = NULL;
            KMap(TDate,KTSlice*)::iterator itRebate 
                = mRebate.find(settleDate);

            if (isKnockOut)
            {
                if (itRebate == mRebate.end()) 
                {
                    throw KFailure("%s: Current observation date is "
                       "%s.  Can't find rebate slice on "
                       "settlement date %s.\n",
                       routine, 
                       GtoFormatDate(curDate),
                       GtoFormatDate(settleDate));
                }
                else
                    rebateTS = (*itRebate).second;
            }


        //
        // Perform exercise
        //

        //
        // Create new observation rate slice
        //
/*            cout << "idx = " << idx << '\t';
            cout << "lb = " << mKnockIO->mBarrierLo[idx] << '\t'
                 << "hb = " << mKnockIO->mBarrierHi[idx] << endl;
            cout << "IoWindow = " << mKnockIO->mIOWindow << endl;
            cout << "KO smooth = " << mKnockIO->mSmooth << endl;

            cout << "mDepKoind = " << mDepKoind << endl;
*/
            SharedPointerDownCastTo(Dep(mDepKoind), koindex);

//          cout << "original formula = " << koindex->getFormula() << endl;
            
            // Reset the formula
            koindex->setFormula(
                        mKnockIO->mBarrierLo[idx],
                        mKnockIO->mBarrierHi[idx],
                        mKnockIO->mIOWindow,
                        mKnockIO->mSmooth);
            //cout << "updated formula = " << koindex->getFormula() << endl;
            

            rateIndexCV = (mKnockIO->mRateIndex[idx])->CplxRate().CurveName(); 

            koidxslice1 = new KTSlice(*mVTree, "koidxslice1", rateIndexCV);

            *koidxslice1 = Dep(mDepKoind)->GetValue();

            koidxslice2 = new KTSlice(*mVTree, "koidxslice2", rateIndexCV);
            
            *koidxslice2  = *koidxslice1;
            *koidxslice2 -= 1.;
            *koidxslice2 *= -1.;

            if (debugLevel == DEBUG_LEVEL_KIO)
            {
                dppLog << "Knock IO rate index slice" << endl;
                dppLog << *koidxslice1 << endl;

                dppLog << "Underlying slice" << endl;
                dppLog << exerTS << endl;
            }
        

            if (mKnockIO->mSmooth != DOUBLE_SMOOTH)
            {
                // NO_SMOOTH or SINGLE_SMOOTH
            
                // mValue = koidxslice1 * exerTS + (1 - koidxslice1) * mValue
                BinaryExpectation(koidxslice1,
                                  &exerTS,
                                  mValue,
                                  mValue);

                if (isKnockOut)
                {
                    //mRebateValue = koidxslice1 * rebateTS + 
                    //               (1-koidxslice1) * mRebateValue
                    BinaryExpectation( koidxslice1,
                                       rebateTS,
                                       mRebateValue,
                                       mRebateValue);
                }

                if (mRunStat)
                {
                    double  ExProb = 1.;
                    double  ExTimeSquare = exTime * exTime;

                    BinaryExpectation( koidxslice1,
                                      &ExProb,
                                       mXT0,
                                       mXT0);
    
                    BinaryExpectation( koidxslice1,
                                      &exTime,
                                       mXT1,
                                       mXT1);

                    BinaryExpectation( koidxslice1,
                                      &ExTimeSquare,
                                       mXT2,
                                       mXT2);
                }

            }
            else  // double smooth
            {
                // (1) unsmoothValue = 
                //         koidxslice1*exerTS + (1-koidxslice1)*mValueUS
                // (2) tmpslice = 
                //         STEP(1-koidxslice1,1,999.INSIDE_RANGE,SMOOTH_NONE)
                // (3) mValue = unsmoothValue + tmpslice*(mValue - mValueUS)
                
                //(1)------------------------------------
                unsmoothValue = new KTSlice(
                                       *mVTree, 
                                       "unsmoothValue", 
                                       rateIndexCV);
                
                unsmoothRebate = new KTSlice(
                                       *mVTree, 
                                       "unsmoothRebate", 
                                       rateIndexCV);

                BinaryExpectation(koidxslice1,
                                  &exerTS,
                                  mValueUS,
                                  unsmoothValue);

                if (isKnockOut)
                {
                    BinaryExpectation(koidxslice1,
                                      rebateTS,
                                      mRebateValueUS,
                                      unsmoothRebate);
                }

                //(2)-------------------------------------
                lbAdj = new KTSlice(*mVTree, "LBAdj", rateIndexCV);
                hbAdj = new KTSlice(*mVTree, "HBAdj", rateIndexCV);
                tmpslice = new KTSlice(*mVTree, "tmpslice", rateIndexCV);

                KKnockIO    ioWindowAdj = CRX_KNOCK_IN;
                KSmooth     smoothAdj   = NO_SMOOTH;

                *lbAdj = 1.0;
                *hbAdj = 99.;

                mVTree->TSliceSpecialOper(
                                    *tmpslice, 
                                    "STEP",
                                    (void*) koidxslice2,
                                    (void*) lbAdj,
                                    (void*) hbAdj,
                                    (void*) &ioWindowAdj,
                                    (void*) &smoothAdj);

                //(3)------------------------------------
                *mValue -= *mValueUS;
                *mValue *= *tmpslice;
                *mValue += *unsmoothValue;
                
                if (isKnockOut)
                {
                    *mRebateValue -= *mRebateValueUS;
                    *mRebateValue *= *tmpslice;
                    *mRebateValue += *unsmoothRebate;
                }
                
                if (mRunStat)
                {
                    double  ExProb = 1.;
                    double  ExTimeSquare = exTime * exTime;

                    BinaryExpectation( koidxslice1,
                                      &ExProb,
                                       mXT0,
                                       mXT0);
    
                    BinaryExpectation( koidxslice1,
                                      &exTime,
                                       mXT1,
                                       mXT1);

                    BinaryExpectation( koidxslice1,
                                      &ExTimeSquare,
                                       mXT2,
                                       mXT2);
                }


                // Update unsmoothed value slices
                // Change ko index formula to unsmoothed fashion
                koindex->setFormula(
                        mKnockIO->mBarrierLo[idx],
                        mKnockIO->mBarrierHi[idx],
                        mKnockIO->mIOWindow,
                        NO_SMOOTH);
            
                koidxslice_unsmooth = 
                    new KTSlice(*mVTree, "koidxslice_unsmooth", rateIndexCV);

                *koidxslice_unsmooth = Dep(mDepKoind)->GetValue();
        
                BinaryExpectation(koidxslice_unsmooth,
                                  &exerTS,
                                  mValueUS,
                                  mValueUS);

                *mRebateValueUS = *rebateTS;
            }// double smooth
        
            if (debugLevel == DEBUG_LEVEL_KIO)
            {
                dppLog << "Knock IO value slice AFTER exercise: " 
                       << endl;
                dppLog << *mValue << endl;

                if (isKnockOut)
                    dppLog << *mRebateValue << endl;
            }



        //
        // Erase exerTs and rebateTS
        //
            delete (*itEx).second;
            mExer.erase(itEx);

            if (isKnockOut)
            {
                delete (*itRebate).second;
                mRebate.erase(itRebate);
                rebateTS = NULL;
            }

        }   // if modelObsDate

    }     // for idx


    //
    // Last TP: store the result.
    //
    if (mVTree->TPIdxCurrent() == 0) 
    {
        if (mValue == NULL) 
        {
            throw KFailure("%s: no observation event "
                "after today date.\n", routine);
        }


        // Store PV
        mResults.insert(
            String_Double("PV", mValue->GetCenter())
        );

        // Calculate PV value of the underlying
        // by summimg all dependencies.
        //
        for (idxD=0; idxD<this->NumDep(); idxD++) 
        {
            if (!(Dep(idxD)->IsType("KVPToolRate")) &&
                !(Dep(idxD)->IsType("KVPToolCplxRate")) &&
                !(Dep(idxD)->IsType("KVPToolIndicator")))
            underlyingValue += Dep(idxD)->GetResults()["PV"];
        }

        mResults.insert(
            String_Double("UNDER", underlyingValue)
        );

        // Expected time to hit
        //
        if (mRunStat) 
        {
            xt0 = mXT0->GetCenter();
            xt1 = mXT1->GetCenter();
            xt2 = mXT2->GetCenter();

            double  pe,         // prob exercise
                    te,         // expected exer time
                    se,         // var exer time
                    tZ;         // last exer time
            TDate   lastDate = mKnockIO->mObsDates.back();

            tZ = (lastDate - mVTree->TPToday()) / 365e0;

            pe = xt0;
            te = xt1 + (1-pe)*tZ;

            mResults["XTPROB"] = pe;
            if (isKnockOut)
                mResults["XTPROB"] = 1e0 - pe;
            mResults["XTEXP"]  = xt1;
            mResults["XTFUG"]  = te;

            se = (xt2 + (1-pe)*tZ*tZ) - te * te;
            if (se < 0e0) 
            {
                mResults["XTSDV"] = sqrt(-se);
            } 
            else 
            {
                mResults["XTSDV"] = sqrt(se);
            }
        } 
        else 
        {
            mResults["XTPROB"] = -9999.99;
            mResults["XTEXP"]  = -9999.99;
            mResults["XTFUG"]  = -9999.99;
            mResults["XTSDV"]  = -9999.99;
        }

        // Correct for knock-out using parity, plus the rebate.
        //
        if(isKnockOut)
        {
            double koPV   = mResults["UNDER"]
                - mResults["PV"]
                + mRebateValue->GetCenter();

            mResults["PV"] = koPV;
        }

    }


    //
    // Done updating, set flag
    //
    UpdateDone();


    delete indexRate;
    delete koidxslice1;
    delete tmpslice;
    delete koidxslice_unsmooth;
    delete unsmoothValue;
    delete unsmoothRebate;
    delete lbAdj;
    delete hbAdj;
 }
 catch (KFailure) 
 {
    delete indexRate;
    delete koidxslice1;
    delete tmpslice;
    delete koidxslice_unsmooth;
    delete unsmoothValue;
    delete unsmoothRebate;
    delete lbAdj;
    delete hbAdj;

    // Don't need to be here anymore. 
    // Default destructor should take care of it
    //
    for (idx=0; idx<n; idx++) 
    {
        //
        // Retrieve settlement value time slice
        //
        settleDate = mKnockIO->mSettleDates[idx];

        KMap(TDate,KTSlice*)::iterator itEx = mExer.find(settleDate);

        if (itEx != mExer.end()) 
        {
            delete (*itEx).second;
            mExer.erase(itEx);
        }

        KMap(TDate,KTSlice*)::iterator itRebate 
            = mRebate.find(settleDate);

        if (itRebate != mRebate.end()) 
        {
            delete (*itRebate).second;
            mRebate.erase(itRebate);
        }

    }

    throw KFailure("%s: failed.\n", routine);
 }
}


//---------------------------------------------------------------

KTSlice&
KVPToolKnockIONew::GetValue()
{
    bool    isKnockOut = false;
    if (mKnockIO->mIOType == CRX_KNOCK_OUT)
        isKnockOut = true;


    // Calculate value of the underlying
    // by summing all dependencies.
    // Correct for knock-out using parity, plus the rebate.
    //
    if (mGetValue == NULL) 
    {
        mGetValue = new KTSlice(*mVTree, "value", mCurveName);
        ASSERT_OR_THROW(mGetValue != NULL);
    }

    if (isKnockOut) 
    {
        *mGetValue = 0e0;
        for (int idxD=0; idxD<this->NumDep(); idxD++) 
        {
            if (!(Dep(idxD)->IsType("KVPToolRate")) &&
                !(Dep(idxD)->IsType("KVPToolCplxRate")) &&
                !(Dep(idxD)->IsType("KVPToolIndicator")))
                *mGetValue += Dep(idxD)->GetValue();
        }
        if (mValue!=NULL)
            *mGetValue -= *mValue;
        *mGetValue += *mRebateValue;
    } 
    else 
    {
        if (mValue!=NULL)
            *mGetValue = *mValue;
        else
            *mGetValue = 0e0;
    }

    return(*mGetValue);
}


void
KVPToolKnockIONew::setCurveName()
{
static  char    routine[] = "KVPToolKnockIONew::setCurveName";
    String       curveName;
    int          nbinstr, i;
    int          cIdx, cIdxMax, discIdx;
    try
    {
        cIdxMax     = -100000;
        curveName   = mKnockIO->GetDiscName();
        discIdx     = mVTree->GetCurveIdx(curveName);
        nbinstr     = mKnockIO->nbIndexInstr();

        for ( i = 0; i < nbinstr; i++)
        {
           if (mKnockIO->mRateIndex[0]->RateIndex(i).IsFloating())
                curveName = 
                       mKnockIO->mRateIndex[0]->RateIndex(i).CurveName();

           cIdx = mVTree->GetCurveIdx(curveName);
           if (cIdx > cIdxMax)
               cIdxMax = cIdx;

        }
        if (discIdx > cIdxMax)
            cIdxMax = discIdx;

        mRateIdxName = mVTree->GetCurveName(cIdxMax);
        for (i = 0; i < mKnockIO->mRateIndex.size(); i++)
            mKnockIO->mRateIndex[i]->CplxRate().SetCurveName(mRateIdxName);
   }
   catch(KFailure) 
   {
       throw KFailure("%s: failed.\n",routine);
   }
}

template<class T1, class T2>
void    
KVPToolKnockIONew::BinaryExpectation(
                              const KTSlice* probUp, 
                              const T1*      valueUp,
                              const T2*      valueDown,
                              KTSlice*       meanValue)
{
static  char    routine[] = "KVPToolKnockIONew::BinaryExpectation";

    KTSlice *probDown = NULL;
    KTSlice *tmpslice = NULL;
        
    try
    {
        probDown = new KTSlice (*probUp);
        tmpslice = new KTSlice (*probUp);

        *probDown -= 1.;
        *probDown *= -1.;
        
        *tmpslice  = *valueUp;
        *tmpslice *= *probUp;

        *meanValue  = *valueDown;
        *meanValue *= *probDown;

        *meanValue += *tmpslice;

        delete probDown;
        delete tmpslice;

    }
    catch (KFailure)
    {
    	delete probDown;
        delete tmpslice;
       throw KFailure("%s: failed.\n",routine);
    }
}
