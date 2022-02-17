#include "ktmxirtree.h"
#include "kmrntree.h"
#include <math.h>

extern "C"{
#include "drlmem.h"
#include "drlvtype.h"			/* DrlVTypeVectAdd */
};

#include <fstream>
using namespace TMX;


KTMXirTree::KTMXirTree() : KPirTree()
{
    NmrToCcy    = FALSE;
    CcyToNmr    = FALSE;
    TMXFlag     = FALSE;

    NbDailyPts  = Nb_Daily_Pts;

}

KTMXirTree::~KTMXirTree()
{
    if (TMXFlag == TRUE)
    {
        if( mktvol_data != NULL)
        {
            delete mktvol_data;
        }
        if (dev_data != NULL)
        {
            Dev_Free(dev_data, tree_data);
            delete dev_data;
        }
        if (tree_data != NULL)
        { 
            Tree_Free(tree_data);
            delete tree_data;
        }
        KMap(int, double*)::iterator iter;
        for( iter = NmrInv.begin();
        iter != NmrInv.end();
        ++iter)
        {
            sliceDelete( (*iter).second );
        }

        for(iter = NmrInvLag.begin();
        iter != NmrInvLag.end();
        ++iter)
        {
            sliceDelete( (*iter).second );
        }

        for(KMap(int, KTSlice*)::iterator iterTSlice = mCurrIRNmrInv.begin();
        iterTSlice != mCurrIRNmrInv.end();
        ++iterTSlice)
        {
            TSliceDestroy( *(iterTSlice->second) );
            delete iterTSlice->second;
        }

        for (KMap(TDate, double*)::iterator iterDev = mDiffNmrInv.begin();
            iterDev != mDiffNmrInv.end();
            ++iterDev)
        { 
            sliceDelete( (*iterDev).second);
        }
    }
}

void KTMXirTree::DeleteMemory()
{
    KPirTree::DeleteMemory();
    
    if (TMXFlag == TRUE)
    {
        if( mktvol_data != NULL)
        {
            delete mktvol_data;
        }
        if (dev_data != NULL)
        {
           Dev_Free(dev_data, tree_data);
           delete dev_data;
        }
        if (tree_data != NULL)
        { 
            Tree_Free(tree_data);
            delete tree_data;
        }
        KMap(int, double*)::iterator iter;
        for( iter = NmrInv.begin();
        iter != NmrInv.end();
        ++iter)
        {
            sliceDelete( (*iter).second );
        }

        for(iter = NmrInvLag.begin();
        iter != NmrInvLag.end();
        ++iter)
        {
            sliceDelete( (*iter).second );
        }

        for(KMap(int, KTSlice*)::iterator iterTSlice = mCurrIRNmrInv.begin();
        iterTSlice != mCurrIRNmrInv.end();
        ++iterTSlice)
        {
            TSliceDestroy( *(iterTSlice->second) );
            delete iterTSlice->second;
        }

        for (KMap(TDate, double*)::iterator iterDev = mDiffNmrInv.begin();
            iterDev != mDiffNmrInv.end();
            ++iterDev)
        { 
            sliceDelete( (*iterDev).second);
        }

        NmrDates.clear();
        CritDates.clear();
        NmrInv.clear();
        NmrInvLag.clear();
        mDiffNmrInv.clear();
        mCurrIRNmrInv.clear();    
    }
}

void
KTMXirTree::Initialize(
        TDate              todayDt,     // (I) today's date
        KMrParam           &mrPar,      // (I) full dimension mr info
        KSmileParam        &irSmilePar, // (I) ir smile info
        int                nIRDim,      // (I) IR dimension
        KVector(int)       &cvTypes,    // (I) cv types (KV_DIFF...)
        KMap(int, KZCurve) &cv,         // (I) array of curves
        KMap(int, String)  &cvNames,    // (I) array of curve names
        KMap(int, TDate)   &cvValueDates, // (I) array of cv value dates
        KVector(TDate)     &volDates,   // (I) volatility dates
        KVector(KVector(double)) &factVol, // (I) spot vols
        KResetBank         &resetBank)  // (I) rate reset bank
    {
        static char routine[] = "KTMXirTree::Initialize";
        
        try
        {
           dppLog << "===========" << routine << "============="<< endl;
           KPirTree::Initialize(
                            todayDt,
                            mrPar,
                            irSmilePar,
                            nIRDim,
                            cvTypes,
                            cv,
                            cvNames,
                            cvValueDates,
                            volDates,
                            factVol,
                            resetBank);
            if (TMXFlag == TRUE)
            {
                mktvol_data = new(TMX::MKTVOL_DATA);
                tree_data   = new(TMX::TREE_DATA);
                dev_data    = new(TMX::DEV_DATA);

                ASSERT_OR_THROW( 
                    mktvol_data != NULL &&
                    tree_data   != NULL &&
                    dev_data    != NULL);

                MktVol_Init(mktvol_data);
                Dev_Init(dev_data);
                Tree_Init(tree_data);
 
                tree_data->NbSigmaMax = mrPar.mNumStdevCut; 

                tree_data->CvDisc = 1;
                tree_data->CvDiff = 0;
                tree_data->CvIdx1 = 1;
                tree_data->CvIdx2 = 2;   

                tree_data->NbFactor = mIRDim; 
                tree_data->Ppy      = mrPar.mPpy;
                ASSERT_OR_THROW( tree_data->NbFactor == 1);/*IR is one factor */

                outStdDevs = 0;             // Nb std devs in tree outer limits
                spotVolInterp = CEILIDX;    // spot vol interpolation method
     
            }
        }
        catch(KFailure)
        {
            throw KFailure("%s: failed\n", routine);
        }

    }

void
KTMXirTree::WrapperEnvPack( KMarketCurves &mktCurves,
                            KVolDiag      &irVolDiag,
                            KMrParam      &irMrParam,
                            KSmileParam   &irSmileParam,
                            const String        &CvDiscName)
{
        static char routine[] = "KTMXirTree::WrapperEnvPack";

        try
        {
            TDate today       = mktCurves.Today();
            KVector(int)::iterator idx = mktCurves.mCVTypes.begin();
            TDate valueDate   = mktCurves.GetValueDate(*idx);

            // Check whether it is one factor
            ASSERT_OR_THROW(irMrParam.mNumFact == 1);

            if (today != valueDate)
                throw KFailure("%s failed: Today (%s) != ValueDate (%s)",
                                routine,
                                GtoFormatDate(today),
                                GtoFormatDate(valueDate));
            
            // In the current TMX tree, tree starts from the valueDate 
            // instead of today.  

            setToday(valueDate);
        
            // Set TMX tree valueDate.  This is necessary so that the benchmark
            // dates in the CET are correct.

            setValueDate(valueDate); 

            // Pack curves and vols
            PackCurves(mktCurves, CvDiscName);


            PackVolData(irVolDiag, mktCurves, irSmileParam, irMrParam);

            PrintEnv();
 
        }
        catch(KFailure)
        {
            throw KFailure("%s: failed\n", routine);
        }
        

}

//---------------------------------------------------------------
// Compute a zero reset.  The pure interest tree can
// deal with two types of zeros: 
// 1. rate out of diffuse curve.
// 2. rate with deterministic spread over diffuse curve.
//
// Must be used after TSliceDev
void 
KTMXirTree::Get(KTSlice &zeroTS, const KZeroReset &zeroReset)
{
    static char routine[] = "KTMXirTree::Get(zeroReset)";

    if ( TMXFlag == FALSE )
    {
        KPirTree::Get( zeroTS, zeroReset);
        return;
    }
    try
    {
        if ( NmrToCcy == TRUE )
        {
            // zBanks are Ccy
            KPirTree::Get(zeroTS, zeroReset);
        }
        else
        {
            throw KFailure("%s failed: %s is not a critical date",
                           routine,
                           GtoFormatDate(TPIdxCurrent()));
        }
    }
    catch(KFailure)
    {
        throw KFailure("%s: failed\n", routine);
    }

}


void KTMXirTree::ClearNmrCcyFlag()
{
    NmrToCcy = FALSE;
    CcyToNmr = FALSE;
}

void KTMXirTree::SetNmrToCcy(int tpIdx)
{
    static char routine[] = "KTMXirTree::SetNmrToCcy";
    try
    {
        if (tpIdx == 0) NmrToCcy = TRUE;
        
        KVector(TDate)::iterator itDate = 
            find(CritDates.begin(), CritDates.end(), TPDates[tpIdx]);

        if (itDate != CritDates.end())
        {
            NmrToCcy = TRUE;
        }
        else
        {
            NmrToCcy = FALSE;
        }
    }
    catch(KFailure)
    {
        throw KFailure("%s: failed\n", routine);
    }
    
}

void
KTMXirTree::SetCcyToNmr()
{
    CcyToNmr = NmrToCcy;
}

bool
KTMXirTree::IsNmrDate(int tpIdx)
{
    static char routine[] = "KTMXirTree::IsNmrDate";
    try
    {
        if (tpIdx == 0) return true;
        
        KVector(TDate)::iterator itDate = 
            find(NmrDates.begin(), NmrDates.end(), TPDates[tpIdx]);

        if (itDate != NmrDates.end())
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    catch(KFailure)
    {
        throw KFailure("%s: failed\n", routine);
    }

}

bool
KTMXirTree::IsCritDate(int tpIdx)
{
    static char routine[] = "KTMXirTree::IsCritDate";
    try
    {
        if (tpIdx == 0) return true;
        
        KVector(TDate)::iterator itDate = 
            find(CritDates.begin(), CritDates.end(), TPDates[tpIdx]);

        if (itDate != CritDates.end())
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    catch(KFailure)
    {
        throw KFailure("%s: failed\n", routine);
    }

}

// Update NmrInv and NmrInvLag
void 
KTMXirTree::NmrUpdate(int tpIdx)
{
    static char routine[] = "KTMXirTree::NmrUpdate";

    try
    {
        if (tpIdx == NbTP)
            return;

        
        double  *sliceDiff;
        double  ZRatio;
        double  DiffTermZero, IdxTermZero;

        KVector(int)::iterator iter;
        KMap(int, double*)::iterator iterNmrInv;
        KMap(int, double*)::iterator iterNmrInvLag;

        if (CcyToNmr == TRUE)
        {
            tpIdxCurrent = tpIdx + 1;
            for(iter = mCVTypes.begin();
                iter != mCVTypes.end();
                ++iter)
            {
        
                iterNmrInvLag = NmrInvLag.find(*iter);
                iterNmrInv    = NmrInv.find(*iter); 

                KMrNTree::sliceUnaryOper(
                            iterNmrInvLag->second,
                            mIRDim,
                            iterNmrInv->second,
                            COPY);
            }
            tpIdxCurrent = tpIdx;
        }

        if (NmrToCcy == FALSE) return; 

        // update NmrInv
        KMap(int, KTSlice*)::iterator iterSlice;
        KMap(TDate, double*)::iterator itDate = mDiffNmrInv.find(TPDate(tpIdx));

        if(itDate == mDiffNmrInv.end())
        {
            throw KFailure( "%s: NmrInv is not available on %s\n", 
                            routine,
                            GtoFormatDate(TPDate(tpIdx)));
        }
        sliceDiff = itDate->second;
        DiffTermZero = GetZeroPrice(KV_DIFF)[NbTP] / 
                       GetZeroPrice(KV_DIFF)[tpIdx];
        // 
        for(iter = mCVTypes.begin();
            iter != mCVTypes.end();
            ++iter)
        {
        
            //iterNmrInvLag = NmrInvLag.find(*iter);
            iterNmrInv    = NmrInv.find(*iter); 
            iterSlice     = mCurrIRNmrInv.find(*iter);

 

            IdxTermZero = GetZeroPrice(*iter)[TPNum()] / 
                          GetZeroPrice(*iter)[tpIdx];
            ZRatio      = DiffTermZero / IdxTermZero;

            KMrNTree::sliceUnaryOper(
                            iterNmrInv->second,
                            mIRDim,
                            sliceDiff,                            
                            COPY);

            KMrNTree::sliceScalarOper( 
                             iterNmrInv->second,
                             mIRDim,
                             ZRatio,
                             MULT);

                    
            // update mCurrIRNmrInv
            if ( iterSlice != mCurrIRNmrInv.end())
            {
                sliceUnaryOper( (double*) ((*iterSlice).second->mData),
                                mIRDim,
                                iterNmrInv->second,
                                COPY);
            }
            else
            {
                throw KFailure("%s failed: mCurrIRNmrInv is empty\n",routine);
            }
            (*iterSlice).second->SetTpIdx(tpIdx);
            
        }
        //remove current NmrInv from mDiffNmrInv
        sliceDelete((*itDate).second);
        mDiffNmrInv.erase(itDate);
        
        if (debugLevel > DEBUG_LEVEL_TIMELINE)
        {
            dppLog << "--------------" << routine << endl;
            dppLog << "---NmrInvLag-----" << endl;
            dppLog << format("NmrInvLag at TPDate = %s:",
                      GtoFormatDate(TPDate(tpIdx))) 
                   << " Tpidx = " << TPIdxCurrent()
                   << " Bottom = " << mBottom1[tpIdx+1] 
                   << " Top = " << mTop1[tpIdx+1] << endl;
            KVector(int)::iterator iter = mCVTypes.begin();
            slicePrint(
                            (NmrInvLag.find(*iter))->second,
                            mIRDim,
                            tpIdx+1,
                            FALSE,
                            dppLog);
            dppLog << "---NmrInv-----" << endl;
            dppLog << format("NmrInv at TPDate = %s:",
                      GtoFormatDate(TPDate(tpIdx))) 
                   << " Bottom = " << mBottom1[tpIdx]
                   << " Top = " << mTop1[tpIdx] << endl;
            slicePrint(
                            (NmrInv.find(*iter))->second,
                            mIRDim,
                            tpIdx,
                            FALSE,
                            dppLog);
            
        }
    }
    catch(KFailure)
    {
        throw KFailure("%s: failed\n", routine);
    }
}


double*
KTMXirTree::GetNmrInv(int curveIdx)
{
    static char routine[] = "KTMXirTree::GetNmrInv";

    KMap(int, double*)::iterator iter = NmrInv.find(curveIdx);

    if(iter == NmrInv.end())
		throw KFailure("%s: invalid curve index (%d). NmrInv "
			       "factor for curve [%d] does not exist.\n",
				routine, curveIdx, curveIdx);

	return (*iter).second;
    
}

double*
KTMXirTree::GetNmrInvLag(int curveIdx)
{
    static char routine[] = "KTMXirTree::GetNmrInvLag";

    KMap(int, double*)::iterator iter = NmrInvLag.find(curveIdx);

    if(iter == NmrInvLag.end())
		throw KFailure("%s: invalid curve index (%d). NmrInvLag "
			       "factor for curve [%d] does not exist.\n",
				routine, curveIdx, curveIdx);

	return (*iter).second;
    
}



KTSlice& 
KTMXirTree::TSliceDev(KTSlice &ts, const String &NmrCurveName)
{
    static char routine[] = "KTMXirTree::TSliceDev";

    int        NmrCurveIdx;             // Nmr curve index
    int        tpIdxTree,               // current time point of the slice
               tpIdxSlice;              // current time point of the tree

    double     *unitSlice = NULL;       // Discount slice is set to 1.
    double     *NmrInvSlice = NULL;

    if ( TMXFlag == FALSE ) 
    {
        KPirTree::TSliceDev( ts, NmrCurveName );
        return ts;
    }
    try
    {

        tpIdxTree   = this->TPIdxCurrent();
        tpIdxSlice  = ts.GetTpIdx();

        // Dev empty slices does nothing
        if (ts.mData == NULL) return ts;
 
        // Don't do anything at last timept
        if (this->TPIdxCurrent() == TPNum()) return ts;
 
        // Check timepoint consistency
        if (tpIdxSlice != tpIdxTree+1)
                throw KFailure("%s: attempting to Ev slice `%s' "
                        "with TP %d on tree at TP %d.\n",
                         routine, ts.GetSliceName().c_str(),
                        tpIdxSlice, tpIdxTree);
 

	    // Check slice dimension
	    if (ts.GetSliceDim() != mIRDim)
		    throw KFailure("%s: slice dimension (%d) does not equal to "
			       "the IR tree dimension (%d).\n",
			       routine, ts.GetSliceDim(), mIRDim);
        

        // Nmr curve index
         NmrCurveIdx = GetCurveIdx(NmrCurveName);
         
        //
        // discount slices are set to unit vector
         unitSlice = sliceNew(this->mIRDim);
         sliceScalarOper( unitSlice,
                          this->mIRDim,
                          1.,
                          COPY);


        if ( CcyToNmr == TRUE )
        {
            NmrInvSlice = GetNmrInvLag(NmrCurveIdx);
            // Multiply NmrInv
            sliceUnaryOper( (double*) ts.mData,
                            this->mIRDim,
                            NmrInvSlice,
                            MULT);

        }

        KMrNTree::sliceEv(
                (double*) ts.mData,
                mIRDim,
                unitSlice,
                tpIdxTree);
 
        if ( NmrToCcy == TRUE )
        {
            NmrInvSlice = GetNmrInv(NmrCurveIdx);
            // Divide NmrInv
            sliceUnaryOper( (double*) ts.mData,
                            this->mIRDim,
                            NmrInvSlice,
                            DIV);            
        }

        if (debugLevel >= DEBUG_LEVEL_GEOMETRY) 
        {
              dppLog << format("=============%s: CurrStatePr at %s\n", 
                     routine, GtoFormatDate(TPDateCurrent()));
              dppLog << "curve " << ts.GetSliceName() << endl;
              slicePrint(
                        (double*) ts.mData,
                        mIRDim, 
                        tpIdxTree,   
                        FALSE, 
                        dppLog);
        }

        ts.SetTpIdx(tpIdxTree);
        sliceDelete(unitSlice);


        return ts;
    }
    catch(KFailure)
    {
        sliceDelete(unitSlice);
        throw KFailure("%s: failed\n", routine);
    }

}

void 
KTMXirTree::Update(int tpIdx)
{
    static char routine[] = "KTMXirTree::Update";

    if ( TMXFlag == FALSE )
    {
        KPirTree::Update(tpIdx);
        return;
    }
    try
    {
        TDate       currDate = TPDate(tpIdx);

        if (tpIdx != TPNum())
        {
            KMrNTree::Update(tpIdx);
        }

        // Update NmrInv and NmrInvLag
        NmrUpdate(tpIdx);

        for (KMap(int, KZeroBank)::iterator itZB = mZBanks.begin();
             itZB != mZBanks.end(); ++itZB)
             {
                 (*itZB).second.Update(tpIdx);
             }


    }
    catch(KFailure)
    {
        throw KFailure("%s: failed\n", routine);
    }

}

void
KTMXirTree::CetTreeFree()
{
    static char routine[] = "KTMXirTree::CetTreeFree";
    try
    {
        int i;
        int Area;
    
        if (tree_data->NbFactor == 1)
        {
            Area = tree_data->Width[0];
        }
        if (tree_data->NbFactor == 2)
        {
            Area = tree_data->Width[0] * tree_data->Width[1];
        }
        if (tree_data->NbFactor == 3)
        {
            Area = tree_data->Width[0] * tree_data->Width[1] * 
                   tree_data->Width[2];
        }

        if (tree_data->EDevDate != NULL)
        {
            TMX::Free_DR_Array (tree_data->EDevDate,LONG,0, 
                                tree_data->NbEDevDates-1);
        }

        if (tree_data->EDevStPrice != NULL)
        {
            for (i = 0; i < tree_data->NbEDevDates; i++)
            {
                TMX::Free_Slice(tree_data->EDevStPrice[i],tree_data);
            }

            TMX::Free_DR_Array(tree_data->EDevStPrice,DOUBLE_PTR,0,
                               tree_data->NbEDevDates-1);
        }
    
        if (tree_data->NmrInv != NULL)
        {
            for (i = 0; i < tree_data->NbNmr; i++)
            {
                TMX::Free_Slice(tree_data->NmrInv[i],tree_data);
            }

            TMX::Free_DR_Array (tree_data->NmrInv,DOUBLE_PTR,0,tree_data->NbNmr-1);
        }
    

        if (tree_data->LastZero != NULL)
        {
            TMX::Free_Slice(tree_data->LastZero, tree_data);
        }

        if (tree_data->Libor != NULL)
        {
            for (i = 0; i < tree_data->NbNmr; i++)
            {
                TMX::Free_DR_Array(tree_data->Libor[i],DOUBLE,-1,Area);
            }

            TMX::Free_DR_Array(tree_data->Libor,DOUBLE_PTR,0,tree_data->NbNmr-1);
        }
    }
    catch(KFailure)
    {
        throw KFailure("%s: failed\n", routine);
    }
}


void 
KTMXirTree::TMXirTreeCet( KMrParam& irMrParam)
{
    static char routine[] = "KTMXirTree::TMXirTreeCet";

    try
    {
        int i, j;
        KVector(TDate)::iterator itDate;
        KVector(double) factVol;

        IF_FAILED_THROW(Cet_Main(TRUE,
                                 t_curve,
                                 mktvol_data,
                                 tree_data));

        IF_FAILED_THROW( Build_Tree( t_curve,
                                     mktvol_data,
                                     tree_data));

//        CetTreeFree();

        if (debugLevel > DEBUG_LEVEL_TIMELINE) 
        {
            dppLog << "----------------tree_data->NmrInv------------" << endl;
            
            int NmrIdx = tree_data->NbNmr-1;
            long curDate;
            int offset;
            for (int t = tree_data->NbTP; t>=0; t--)
            {
                curDate = tree_data->TPDate[t];
                if (mktvol_data->NmrDate[NmrIdx] != curDate)
                    continue;
                dppLog << "NmrDate " << curDate << endl;
                offset = Node_Offset(1,0,0,t,tree_data);
                for(i= tree_data->Bottom1[t]; i<= tree_data->Top1[t]; i++)
                {
                    dppLog << i << '\t' << setprecision(10) 
                           <<(tree_data->NmrInv[NmrIdx])[i+offset] << endl;
                }
                NmrIdx--;
            
            }

            dppLog << format("----------------%s-----------", routine) <<endl;
            dppLog << "vollogn = 1: TMX always uses lognormal backbone"<< endl;
            for ( i = 0; i <= tree_data->NbTP; i++)
            {
                dppLog << "TPDate[" <<i<<"] = " << tree_data->TPDate[i] << '\t'  
                        << "MRTree Date = " << GtoFormatDate(TPDate(i)) 
                        << '\t' << " Aweight " ;
                for (j = 0; j < 6; j++)
                    dppLog << setprecision(18) << tree_data->Aweight[j][i] << '\t';
                dppLog << "Top = " << tree_data->Top1[i] << " OutTop = " 
                       << tree_data->OutTop1[i] << " Bottom = " 
                       << tree_data->Bottom1[i] << " outbottom = "
                       << tree_data->OutBottom1[i] << '\t';
                dppLog << "Length= " << tree_data->Length[i] << '\t';
                dppLog << "LengthJ= " << tree_data->LengthJ[i];
                dppLog << endl;
            }
            dppLog << "benchmark factor vols" << endl;
            for (i = 0; i < mktvol_data->NbVol; i++)
            {
                dppLog << "VolDate = " << mktvol_data->VolDate[i] << '\t'
                       << "Aweight " << mktvol_data->Aweight[0][i] << endl;
            }
            dppLog << "MR tree benchmar factor vols dates" << endl;
            for (itDate = mVolDates.begin();
                 itDate != mVolDates.end();
                 itDate++)
                 {
                     dppLog << "VolDate = " << GtoFormatDate(*itDate) << endl;
                 }

            dppLog << routine << " OLD MR tree benchmar factor vols " << endl;
            KVector(KVector(double))::iterator itVol = mFactVol.begin();
            for (itDate = mVolDates.begin(), i = 0;
                 itDate != mVolDates.end();
                 itDate++, i++)
                 {
                     dppLog << "VolDate = " << GtoFormatDate(*itDate) << '\t'
                            << "spotVol = "<< (*itVol)[i] << endl;
                 }
        }

        // clear factor vols
        mVolDates.clear();
        mFactVol.clear();

        if (irDim() != 1)
            throw KFailure("%s failed: The TMX IR dimension should be 1\n",
                           routine);

        double norm = irMrParam.AlphaNorm();
 //       norm2 = norm2 * norm2;
        for ( i = 0; i < mktvol_data->NbVol; i++)
        {
            mVolDates.push_back(DrlDrDate2TDate(mktvol_data->VolDate[i]));
            factVol.push_back(mktvol_data->Aweight[0][i] / norm);
        }
        mFactVol.insert(mFactVol.end(), irDim(), factVol);
        factVol.clear();

        if (debugLevel > DEBUG_LEVEL_TIMELINE) 
        {
            dppLog << routine << " NEW MR tree benchmark factor vols dates" 
                   << endl;
            KVector(KVector(double))::iterator itVol = mFactVol.begin();
            for (itDate = mVolDates.begin(), i = 0;
                 itDate != mVolDates.end();
                 itDate++, i++)
                 {
                     dppLog << "VolDate = " << GtoFormatDate(*itDate) << '\t'
                            << "spotVol = "<< (*itVol)[i] <<'\t' 
                            << (*itVol)[i]*norm << endl;
                 }
        }


    }
    catch(KFailure)
    {
        throw KFailure("%s: failed\n", routine);
    }
}

void KTMXirTree::TreeSetUp()
{
    static char routine [] = "KTMXirTree::TreeSetUp";
    
    if ( TMXFlag == FALSE )
    {
        KPirTree::TreeSetUp();
        return;
    }

    try
    {
        KPirTree::TreeSetUp();
        

        for(KVector(int)::iterator iter = mCVTypes.begin();
        iter != mCVTypes.end();
        ++iter)
        {
            NmrInv.insert(KMap(int, double*)::value_type(
                    *iter,
                    sliceNew(this->mIRDim)));
            NmrInvLag.insert(KMap(int, double*)::value_type(
                    *iter,
                    sliceNew(this->mIRDim)));
            mCurrIRNmrInv.insert(KMap(int, KTSlice*)::value_type(
                    *iter,
                    new KTSlice(*this, "Curr IR NmrInv",diffCvName)));
        }
        

    }
    catch(KFailure)
    {
        throw KFailure("%s: failed\n", routine);
    }
}

void KTMXirTree::CheckTreeValid()
{
    static char routine[] = "KTMXirTree::CheckTreeValid";

    if ( TMXFlag == FALSE)
    {
        KPirTree::CheckTreeValid();
        return;
    }
    // Check if IR dim is 1
    if( mIRDim != 1)
        throw KFailure("%s: IR dimension of the TMX tree must be 1.\n", routine);
    
    // Check if lognormal backbone
    if ( !IS_EQUAL(mBbq, 1.))
        throw KFailure("%s: TMX requires lognormal backbone.\n", routine);

    // Check if KV_DIFF curve exists
    if(find(mCVTypes.begin(), mCVTypes.end(), KV_DIFF) == mCVTypes.end())	
        throw KFailure("%s: No diffuse curve specified in the tree.\n",
        routine);


    if (debugLevel > DEBUG_LEVEL_TIMELINE) 
    {
        dppLog << format("%s:\n", routine);
        dppLog << format("mVolNorm = %lf\n", 	mVolNorm);
        dppLog << format("mVolLogn = %lf\n", 	mVolLogn);
    }

}

void KTMXirTree::CpTreeParam()
{
    static char routine[] = "KTMXirTree::CpTreeParam";
    
    try
    {

        int i, t;

        if (debugLevel > DEBUG_LEVEL_TIMELINE) 
        {
            dppLog << format("----------------%s-----------", routine) <<endl;
            for ( t = 0; t <=NbTP; t++)
            {
                dppLog << GtoFormatDate(TPDates[t]) << " Aweight ";
                for (i = 0; i < mNbFactor*(mNbFactor+1)/2; i++)
                {
                    dppLog << mAweight[t][i] << '\t';
                }
                dppLog << "Top = " << mTop1[t] << " OutTop = "
                       << mOutTop1[t] << " Bottom = " << mBottom1[t]
                       << " OutBottom = " << mOutBottom1[t] <<'\t';
                dppLog << "LengthJ = " << LengthJ[t] << '\t' << endl;
            }       
        }
    
    
        for (i = 0; i < mIRDim; i++)
        {
        tree_data->Width[i]      = mWidth[i];
        tree_data->HalfWidth[i]  = mHalfWidth[i];
        }

        for (t = -1; t <= NbTP; t++)
        {
        // copy spotvol
            for (i = 0; i < tree_data->NbFactor*(tree_data->NbFactor+1)/2; i++)
            {
                tree_data->Aweight[i][t] = mAweight[t][i];
            }
        // copy Top1 & Bottom1 
        // Since the current TMX only accepts 1 factor,
        // copy 1D only 
            tree_data->Bottom1[t]    = mBottom1[t];
            tree_data->Top1[t]       = mTop1[t];
            tree_data->OutBottom1[t] = mOutBottom1[t];
            tree_data->OutTop1[t]    = mOutTop1[t];

        }

    }
    catch(KFailure)
    {
        throw KFailure("%s: failed\n", routine);
    }

}

//
//=============================================================== 
// Compute mDevStatePr such that
// sumproduct(mDevStatePr(T), unit(T)) = zero(0,T)
// mDevStatePr = statePr * NmrInv * mktvol_data.TermZero0;
//===============================================================
void KTMXirTree::eDevSetUp()
{
    static char routine[] = "KTMXirTree::eDevSetUp";

    int t;          // Time step index
    int nDim, nIRDim;

    double  *CurrStatePr  = NULL,
            *CurrStatePrL = NULL;

    double  *CurrNmrInvIR = NULL,
            *CurrNmrInv   = NULL;

    double  *devStatePr  = NULL,
            *devStatePrL = NULL;

    double  *unitSlice = NULL;
 
    double  zeroPrice,
            termZero,
            p;

    KVector(TDate)::iterator iterDev;
    KMap(TDate, double*)::iterator iterNmr;
    
    try 
    {

        if ( mDevDates.empty() )
        {
            return;
        }

        if (debugLevel > DEBUG_LEVEL_TIMELINE) 
        {
            dppLog << format("----------------%s-----------", routine) <<endl;
        }

        mDevOn   = true;
        iterDev  = mDevDates.begin();
        termZero = GetZeroPrice(KV_DIFF)[NbTP];  // Z(0, LastTP)
    
        nDim     = mNbFactor;	// Total dimension
        nIRDim   = mIRDim;	// IR dimension


        ASSERT_OR_THROW( termZero > TINY);

        // Check dimension valid
        if (nDim   < 1 || nDim   > 2 ||
            nIRDim < 1 || nIRDim > 1)
                throw KFailure("%s: invalid factor numbers.\n"
                "Total dimension is %d, IR dimension is %d.\n",
                routine, nDim, nIRDim);

        //
        // Allocate tmp working space
        //
        CurrStatePr  = sliceNew(nDim);
        unitSlice    = sliceNew(nDim);
        
        CurrStatePrL = CurrStatePr + NodeOffset(nDim, 0, 0, 0);
        
        CurrNmrInvIR = sliceNew(nIRDim);
        CurrNmrInv   = sliceNew(nDim);
        //
        // Set initial state price to 1 at TPIdx = 0.
        //
        CurrStatePrL[0] = 1e0;
 
        //
        // Set up ExpressDEV (r) tool.
        //
        mDevOn = true;
        iterDev = mDevDates.begin();

        // If the start date is a express dev date
        if ( *iterDev == TPToday())
        {
                devStatePr   = sliceNew(nDim);
                devStatePrL = devStatePr + NodeOffset(nDim, 0, 0, 0);
                devStatePrL[0] = 1e0;

                mDevStatePr.insert(KMap(TDate, double*)::value_type(
                       TPToday(), devStatePr));

                devStatePr = NULL;

                iterDev++;
        }

        for (t = 0; t < NbTP; t++)               
        {                                   
            if (iterDev == mDevDates.end())
                break;
            
            tpIdxCurrent = t;   // set the current time point

            //
            // Compute transition probabilities between t and t+1
            //
            KMrNTree::Update(t);

            sliceScalarOper( unitSlice,
                             nDim,
                             1.,
                             COPY);

            sliceFw(CurrStatePr, nDim, unitSlice, t);

            tpIdxCurrent = t+1;

            //
            // sum of all state prices should equal to one
            //
            sliceSpecialOper(CurrStatePr, nDim, "sum", &p);
#ifndef __NO_CALIB__
            if (fabs (p - 1.) > STATE_PRICE_TOL) 
            {
                        throw KFailure("%s: at time point %d, EDev %d-D state "
                            "prices don't add up to Z (res=%lf).\n",
                            routine, t, nDim, fabs (p - 1.));
            }
#endif
            //
            // Express DEV tool: store state prices in full nDim dimension.
            // mDevStatePr = statePr * NmrInv * TermZero0;
            
            if ( *iterDev == TPDate(tpIdxCurrent))
            {
              /*  if (debugLevel >= DEBUG_LEVEL_GEOMETRY) 
                {
                    dppLog << format("=============%s: CurrStatePr at %s\n", 
                        routine, GtoFormatDate(TPDate(tpIdxCurrent)));
                    slicePrint(
                        CurrStatePr,
                        nDim, 
                        tpIdxCurrent,   
                        FALSE, 
                        dppLog);
                }*/
                iterNmr = mDiffNmrInv.find(TPDate(tpIdxCurrent));

                if (iterNmr == mDiffNmrInv.end())
                    throw KFailure("%s failed: Express dev date %s "
                                   "are not in mDiffNmrInv\n",
                                   routine, GtoFormatDate(TPDate(tpIdxCurrent)));


                zeroPrice = GetZeroPrice(KV_DIFF)[tpIdxCurrent];  // Z(0, LastTP)
                devStatePr   = sliceNew(nDim);

                sliceUnaryOper( devStatePr,
                                nDim,
                                CurrStatePr,
                                COPY);

                sliceUnaryOper( CurrNmrInvIR,
                                nIRDim,
                                iterNmr->second,
                                COPY);

                if (nDim > nIRDim)
                {
                    //
                    // Expand NmrInv in nDim dimensions, 
                    // which only varies in the first nIRDim 
                    // dimensions and is constant in the 
                    // higher nDim-nIRDim dimensions
                    //
  /*                  if (debugLevel >= DEBUG_LEVEL_GEOMETRY) 
                    {
                        dppLog << format("%s: CurrNmrInvIR at %ld\n", 
                                routine, DrlTDate2DrDate(TPDate(tpIdxCurrent)));
                        slicePrint(
                            CurrNmrInvIR,
                            nIRDim, 
                            tpIdxCurrent,   
                            FALSE, 
                            dppLog);
                    }*/
                    sliceExpand(tpIdxCurrent, 
                                nDim, 
                                nIRDim, 
                                CurrNmrInvIR, 
                                CurrNmrInv);
                }
                else
                {
                    sliceUnaryOper( CurrNmrInv,
                                    nDim,
                                    CurrNmrInvIR,
                                    COPY);
                }

                if (debugLevel >= DEBUG_LEVEL_GEOMETRY) 
                {
                    dppLog << format("%s: CurrNmrInv at %ld\n", 
                            routine, DrlTDate2DrDate(TPDate(tpIdxCurrent)));
                    slicePrint(
                        CurrNmrInv,
                        nDim, 
                        tpIdxCurrent,   
                        FALSE, 
                        dppLog);
                }

                sliceUnaryOper( devStatePr,
                                nDim,
                                CurrNmrInv,
                                MULT);

                sliceScalarOper(devStatePr,
                                nDim,
                                termZero,
                                MULT);
                //
                // sum of all dev state prices should equal zero bond
                //
                sliceSpecialOper(devStatePr, nDim, "sum", &p);

#ifndef __NO_CALIB__
                if (debugLevel >= DEBUG_LEVEL_DRIFT) 
                {
                    double TmxZeroPrice = TMX::ZeroPrice( 
                                       DrlTDate2DrDate(*iterDev),
                                       t_curve[tree_data->CvDiff].ValueDate,
                                       t_curve[tree_data->CvDiff].NbZero,
                                       t_curve[tree_data->CvDiff].ZeroDate,
                                       t_curve[tree_data->CvDiff].Zero);
                    dppLog   << routine << " TMXzero = " << setprecision(10) 
                             << TmxZeroPrice << "\tAlib zero = "
                             << setprecision(10) << zeroPrice 
                             << "\ton " << GtoFormatDate(*iterDev) << endl;
                    dppLog   << " dev zero = " << setprecision(10) << p
                             << "\tdiff= " << (TmxZeroPrice - p) << endl;
                }
                /*if (fabs (p - zeroPrice) > STATE_PRICE_TOL) 
                {
                        throw KFailure("%s: at time point %d, EDev %d-D state "
                            "prices don't add up to Z (res=%lf).\n",
                            routine, tpIdxCurrent, nDim, fabs (p - zeroPrice));
                }*/
#endif
                //
                // Store state price at special dev dates
                mDevStatePr.insert(KMap(TDate, double*)::value_type(
                       *iterDev, devStatePr));

                devStatePr      = NULL;

                iterDev++;
            }   //iterDev   


        }  // for t                                         


    // Free tmp memory
    sliceDelete(CurrStatePr);
    sliceDelete(unitSlice);
    sliceDelete(CurrNmrInvIR);
    sliceDelete(CurrNmrInv);
 
    }
    catch (KFailure) {
        // Free tmp memory
        sliceDelete(CurrStatePr);
        sliceDelete(unitSlice);
        sliceDelete(CurrNmrInvIR);
        sliceDelete(CurrNmrInv);
 
        throw KFailure("%s: failed.\n", routine);
    }

} // eDevSetUp

// Calibrate Diffusion index NmrInv at each TP
void
KTMXirTree::CalibrateNmrInv()
{
    static char routine[] = "KTMXirTree::CalibrateNmrInv";
    try
    {
        long    TermDate, CurrDate;
        double  TermZeroDiff;
        double *EDevStPrice;
        double *CurrNmrInv;
        int     idxDev, offset;
        int     i;

        EDevStPrice = Alloc_Slice(tree_data);
        ASSERT_OR_THROW( EDevStPrice != NULL);

        TermDate = tree_data->TPDate[NbTP];
        TermZeroDiff = ZeroPrice(TermDate,
                                 treeToday,
                                 t_curve[tree_data->CvDiff].NbZero,
                                 t_curve[tree_data->CvDiff].ZeroDate,
                                 t_curve[tree_data->CvDiff].Zero);
//        TermZeroDiff = GetZeroPrice(KV_DIFF)[TPNum()];

        ASSERT_OR_THROW(TermZeroDiff > TINY);

//        cout << routine << "tree_data->nbnmr= " << tree_data->NbNmr << endl;

        for (int tpIdx = TPNum(); tpIdx >=0; tpIdx--)
        {
            SetNmrToCcy(tpIdx);
            tpIdxCurrent = tpIdx;
            dev_data->NmrToCcy = NmrToCcy;

            CurrDate = tree_data->TPDate[tpIdx];
            offset   =  Node_Offset(1, 0, 0, tpIdx, tree_data);

            /* Interp Nmr */
            if (tpIdx == 0 &&
                TMX::Daysact(treeToday, treeValueDate) > 0)
            {
                // Due to curve extension, the first tree point is not
                // the NmrDate.
                cout << routine << ": extending curves affect Nmr Inv\n";
                sliceScalarOper( dev_data->NmrInv[tree_data->CvDiff],
                                 this->mIRDim,
                                 1./TermZeroDiff,
                                 COPY);
            }
            else
            {
                IF_FAILED_THROW( Lattice(   dev_data,
                                            tpIdx,
                                            NbTP,
                                            mktvol_data,
                                            tree_data,
                                            t_curve));
            }
            
//            if (debugLevel > DEBUG_LEVEL_GEOMETRY)
            {
                    dppLog << "-----------------" << endl;
                    dppLog << format("%s: Tree prob & NmrInv at TPDate = %s:",
                            routine,
                            GtoFormatDate(TPDate(tpIdx))) << endl;
                    dppLog << "CurrTpIdx " << TPIdxCurrent() << " Bottom = " 
                           << mBottom1[tpIdx] << " Top = "
                           << mTop1[tpIdx] << endl;
                              
                    for (i = tree_data->Bottom1[tpIdx]; 
                        i <= tree_data->Top1[tpIdx]; 
                        i++)
                    {
                        dppLog << i << '\t' << (dev_data->pu)[i + offset]
                               << '\t' << (dev_data->p0)[i+offset] 
                               << '\t' << (dev_data->pd)[i+offset]
                               << '\t' << dev_data->NmrInv[tree_data->CvDiff][i+offset]
                               <<'\t' << tree_data->LastZero[i+offset]
                               << endl;
                    }
            }

            if (IsNmrDate(tpIdx) ||
                IsCritDate(tpIdx))
            {
                CurrNmrInv = sliceNew(mIRDim);
                ASSERT_OR_THROW( CurrNmrInv != NULL);
    
                Copy_Slice ( CurrNmrInv,
                             dev_data->NmrInv[tree_data->CvDiff],
                             tpIdx,
                             tree_data);

                mDiffNmrInv.insert(KMap(TDate, double*)::value_type(
                    TPDate(tpIdx),
                    CurrNmrInv));

                if(debugLevel > DEBUG_LEVEL_TIMELINE &&
                    dev_data->NmrToCcy)
                {
                    dppLog << "-----------------" << endl;
                    dppLog << format("%s: Tree prob, NmrInv, CpNmrInv "
                            "at TPDate = %s:",
                            routine,
                            GtoFormatDate(TPDate(tpIdx))) << endl;
                    dppLog << "CurrTpIdx " << TPIdxCurrent() << " Bottom = " 
                           << mBottom1[tpIdx] << " Top = "
                           << mTop1[tpIdx] << endl;
                              
                    for (i = tree_data->Bottom1[tpIdx]; 
                        i <= tree_data->Top1[tpIdx]; 
                        i++)
                    {
                        dppLog << i << '\t' << (dev_data->pu)[i + offset]
                               << '\t' << (dev_data->p0)[i+offset] 
                               << '\t' << (dev_data->pd)[i+offset]
                               << '\t' << dev_data->NmrInv[tree_data->CvDiff][i+offset]
                               << endl;
                    }
                    dppLog << "cpNmrInv" << endl;
                    slicePrint(
                            (mDiffNmrInv.find(TPDate(tpIdx)))->second,
                            mIRDim,
                            tpIdx,
                            FALSE,
                            dppLog);
                }



                idxDev = TMX::GetDLOffset( tree_data->NbEDevDates,
                                       tree_data->EDevDate,
                                       CurrDate,
                                       CbkEXACT);
                if (idxDev != -999)
                {
                    // update EDevStPrice
                    IF_FAILED_THROW( Copy_Slice( 
                                             EDevStPrice,
                                             tree_data->EDevStPrice[idxDev],
                                             tpIdx,
                                             tree_data));

                    if (debugLevel > DEBUG_LEVEL_GEOMETRY)
                    {
                        dppLog << "-----------------" << endl;
                        dppLog << format(
                            "%s: Tree DevState price at TPDate = %s:",
                            routine,
                            GtoFormatDate(TPDate(tpIdx))) << endl;
                        dppLog <<"tmx tree: \n";
                           
                        for (i = tree_data->Bottom1[tpIdx]; 
                            i <= tree_data->Top1[tpIdx]; 
                            i++)
                        {
                            dppLog << i <<'\t'
                                   << tree_data->EDevStPrice[idxDev][i+offset]
                                   << endl;
                        }
                        dppLog << "copy\n";
                        slicePrint( EDevStPrice,
                                    mIRDim,
                                    tpIdx,
                                    FALSE,
                                    dppLog);
                    }

                    // Check cumlative probability
                    double prob = 0;
                    sliceSpecialOper(EDevStPrice,
                                     mIRDim,
                                     "sum",
                                     &prob);
    
                    if (fabs (prob - 1.) > STATE_PRICE_TOL) 
                    {
                            throw KFailure(
                             "%s: at time point %d, EDev %d-D state "
                             "prices don't add up to 1 (res=%lf).\n",
                              routine, tpIdx, mIRDim, fabs (prob - 1.));
                    }

                    sliceUnaryOper(EDevStPrice,
                                   mIRDim,
                                   CurrNmrInv,
                                   MULT);
                    sliceScalarOper(EDevStPrice,
                                    mIRDim,
                                    TermZeroDiff,
                                    MULT);

                    // check zero price on Critical dates
                    if(IsCritDate(tpIdx)) 
                    {
                        double zbond = 0;
                        sliceSpecialOper(EDevStPrice,
                                         mIRDim,
                                         "sum",
                                         &zbond);
                        double BmkZero = GetZeroPrice(KV_DIFF)[tpIdx];
                        if (fabs (zbond - BmkZero) > STATE_PRICE_TOL) 
                        {
                            throw KFailure(
                             "%s: at time point %s, zero price mismatch "
                             "the benchmark (res=%lf).\n",
                              routine, 
                              GtoFormatDate(TPDate(tpIdx)), 
                              fabs (zbond - BmkZero));
                        }
                        
                    }

                }
                CurrNmrInv = NULL;
                dev_data->CcyToNmr = CcyToNmr;
                CcyToNmr           = NmrToCcy;
            
                if (debugLevel > DEBUG_LEVEL_GEOMETRY)
                {
                    dppLog << "-----------------" << endl;
                    dppLog << format("%s: NmrInv at TPDate = %s:",
                            routine,
                            GtoFormatDate(TPDate(tpIdx))) << endl;
                    slicePrint(
                            (mDiffNmrInv.find(TPDate(tpIdx)))->second,
                            mIRDim,
                            tpIdx,
                            FALSE,
                            dppLog);
                }

            }

        }//for

        if (debugLevel >= DEBUG_LEVEL_DRIFT) 
        {
            dppLog << "nb of dev dates = " << mDevDates.size() << endl;
            for (KVector(TDate)::iterator itDevDate = mDevDates.begin();
                itDevDate != mDevDates.end();
                ++itDevDate)
            {
                dppLog << "DevDate: " << GtoFormatDate(*itDevDate) << endl;
            }
        }
        Free_Slice(EDevStPrice, tree_data);
    }
    catch(KFailure)
    {
        throw KFailure("%s: failed\n", routine);
    }
}


void
KTMXirTree::CalibrateDrift()
{
    static char routine[] = "KTMXirTree::CalibrateDrift";
    
    if ( TMXFlag == FALSE )
    {
        KPirTree::CalibrateDrift();
        return;
    }
    try
    {

        int     i;

        // Initialize dev_data
        IF_FAILED_THROW(Dev_Alloc(  dev_data, 
                                    tree_data));

        // Copy NbNmr
        CalibrateNmrInv();

        eDevSetUp();

        if (debugLevel > DEBUG_LEVEL_TIMELINE) 
        {
            int j;
            KVector(TDate)::iterator itDate;            
            dppLog << format("----------------%s-----------", routine) <<endl;
           // dppLog << "vollogn "<< mktvol_data->VolLogn << endl;
            for ( i = 0; i <= tree_data->NbTP; i++)
            {
                dppLog << "TPDate[" <<i<<"] = " << tree_data->TPDate[i] << '\t'  
                        << "MRTree Date = " << GtoFormatDate(TPDate(i)) 
                        << '\t' << " Aweight " ;
                for (j = 0; j < 6; j++)
                    dppLog << tree_data->Aweight[j][i] << '\t';
                dppLog << "Top = " << tree_data->Top1[i] << " OutTop = " 
                       << tree_data->OutTop1[i] << " Bottom = " 
                       << tree_data->Bottom1[i] << " outbottom = "
                       << tree_data->OutBottom1[i] << '\t';
                dppLog << "LengthJ= " << tree_data->LengthJ[i] << '\t';
                dppLog << endl;
            }
            dppLog << "benchmark factor vols" << endl;
            for (i = 0; i < mktvol_data->NbVol; i++)
            {
                dppLog << "VolDate = " << mktvol_data->VolDate[i] << '\t'
                       << "Aweight " << mktvol_data->Aweight[0][i] << endl;
            }
            dppLog << "MR tree benchmar factor vols dates" << endl;
            for (itDate = mVolDates.begin();
                 itDate != mVolDates.end();
                 itDate++)
                 {
                     dppLog << "VolDate = " << GtoFormatDate(*itDate) << endl;
                 }

            dppLog << "OLD MR tree benchmark factor vols " << endl;
            KVector(KVector(double))::iterator itVol = mFactVol.begin();
            for (itDate = mVolDates.begin(), i = 0;
                 itDate != mVolDates.end();
                 itDate++, i++)
                 {
                     dppLog << "VolDate = " << GtoFormatDate(*itDate) << '\t'
                            << "spotVol = "<< (*itVol)[i] << endl;
                 }
        }

        // set the current index to the last TP
        tpIdxCurrent = NbTP;

    }
    catch(KFailure)
    {
        throw KFailure("%s: failed\n", routine);
    }
}

void 
KTMXirTree::Calibrate()
{
    static char routine[] = "KTMXirTree::Calibrate";
    
    if ( TMXFlag == FALSE )
    {
        KPirTree::Calibrate();
        return;
    }

    try
    {
   
        TreeSetUp();

        CheckTreeValid();

        CalibrateDrift();


    }
    catch(KFailure)
    {
        throw KFailure("%s: failed\n", routine);
    }
}



/* -----------------------------------------------------------------
 * 
 * Setup Nmr dates.  Nmr_Schedule
 */
void 
KTMXirTree::InsertNmrDates_old(KVolDiag& BmkVol)
{

    static char routine[] = "KTMXirTree::InsertNmrDates_old";

    long        SwapSt, SwapMat;
    long        LastProdDate;
    long        LibMatDate;
    long        BmkDate;
    long        TermDate;
    TDate       AlibCurDate;
    TDate      *mNmrDates = NULL;
    int         intvl, NbNmr, NbExpiry;
    int        *NbNmrMax = NULL;
    int         i;
    int         SwapMatMos;
    char        Freq = mktvol_data->Freq;
    
    try
    {
        IF_FAILED_THROW(DrlVTypeVectSort(TPDates, 
                                         &NbTP, 
                                         DRL_TDATE_T, 
                                         TRUE));

        // Last critical date
        mLastDate = TPDate(NbTP-1);
        TermDate  = DrlTDate2DrDate(mLastDate);
        LastProdDate = TermDate;
        cout << mLastDate << "\t" << TermDate <<" NbTP = " << NbTP << endl;

        // Determine number of benchmarks to insert :
        // last one must expire after mLastDate
        i = 0;
        KVector(TDate)::const_iterator itDATE;
        KVector(double)::const_iterator itMAT = BmkVol.mVolMats.begin();
        for (itDATE = BmkVol.mVolDates.begin();
             itDATE != BmkVol.mVolDates.end() && itMAT != BmkVol.mVolMats.end(); 
             ++itDATE, ++itMAT, ++i)
             {
                 SwapSt  = DrlTDate2DrDate(*itDATE);
                 SwapMatMos = static_cast<int> ((*itMAT)*12 + 0.5);
                 SwapMat = Nxtmth(SwapSt, SwapMatMos, 1);
                 TermDate= MAX( LastProdDate, SwapMat );
                 cout << "swapst " << SwapSt << " swapmat " << SwapMat << endl;
                 if (SwapSt >= LastProdDate) break;        
             }
        cout << "TermDate= "<< TermDate << endl;
        NbExpiry = MIN(i+1, BmkVol.mVolDates.size());
        intvl    = 12 / Conv_Freq(Freq);
        TermDate = Nxtmth(TermDate, 2*intvl, 1L);


        NbNmr = 0;

        // add value date
        AlibCurDate = DrlDrDate2TDate(treeValueDate);
        IF_FAILED_THROW( DrlVTypeVectAdd(
                            (void **) &mNmrDates,
                            &NbNmr,
                            NbNmrMax,
                            512,
                            (void *) &AlibCurDate,
                            TRUE,
                            DRL_TDATE_T));
            
        // add Term Date
        AlibCurDate = DrlDrDate2TDate(TermDate);
        IF_FAILED_THROW( DrlVTypeVectAdd(
                            (void **) &mNmrDates,
                            &NbNmr,
                            NbNmrMax,
                            512,
                            (void *) &AlibCurDate,
                            TRUE,
                            DRL_TDATE_T));


        // add benchmark dates
        for (i= 0, itDATE = BmkVol.mVolDates.begin();
             i < NbExpiry ;
             ++itDATE, i++)
             {
                 IF_FAILED_THROW( DrlVTypeVectAdd(
                                    (void **) &mNmrDates,
                                    &NbNmr,
                                    NbNmrMax,
                                    512,
                                    (void *) itDATE,
                                    TRUE,
                                    DRL_TDATE_T));
             }



        // add offset of valueDate
        BmkDate    = DrlTDate2DrDate(*(BmkVol.mVolDates.begin())); // the first benchmark date
        LibMatDate = Nxtmth(treeValueDate,
                            intvl,
                            1L);

        while (LibMatDate < AlibCurDate)
        {
            /* To handle the case where tenor is 1 day - which    */
            /* could potentially result in dccfrac = 0 if DCC = 3 */
            /* & MatDate = 31st. Assumes intvl > 1 day            */
            if (LibMatDate >= BmkDate - 2)
                LibMatDate -= 1L;
            AlibCurDate = DrlDrDate2TDate(LibMatDate);
            IF_FAILED_THROW( DrlVTypeVectAdd(
                            (void **) &mNmrDates,
                            &NbNmr,
                            NbNmrMax,
                            512,
                            (void *) &AlibCurDate,
                            TRUE,
                            DRL_TDATE_T));

            LibMatDate = Nxtmth(LibMatDate, intvl, 1L);
        }


        // add offsets of benchmark dates 
        for (i = 0, itDATE = BmkVol.mVolDates.begin(); 
            i < NbExpiry - 1; 
            i++, itDATE++)
        {
            LibMatDate = DrlTDate2DrDate(*itDATE);
            BmkDate    = DrlTDate2DrDate(*(itDATE+1));

            while( LibMatDate < BmkDate)
            {
                if(LibMatDate >= BmkDate - 2)
                    LibMatDate -= 1L;
                AlibCurDate = DrlDrDate2TDate(LibMatDate);
                IF_FAILED_THROW( DrlVTypeVectAdd(
                                (void **) &mNmrDates,
                                &NbNmr,
                                NbNmrMax,
                                512,
                                (void *) &AlibCurDate,
                                TRUE,
                                DRL_TDATE_T));
                LibMatDate = Nxtmth(LibMatDate, intvl, 1L);
            }
        }   

        // add offset of last benchmark date till the tree last date

        LibMatDate = BmkDate;
        while( LibMatDate < TermDate)
        {
            if (LibMatDate >= TermDate - 2)
                LibMatDate -= 1L;
            AlibCurDate = DrlDrDate2TDate(LibMatDate);
            IF_FAILED_THROW( DrlVTypeVectAdd(
                                (void **) &mNmrDates,
                                &NbNmr,
                                NbNmrMax,
                                512,
                                (void *) &AlibCurDate,
                                TRUE,
                                DRL_TDATE_T));
            LibMatDate = Nxtmth(LibMatDate, intvl, 1L);
        }

	    IF_FAILED_THROW(DrlVTypeVectSort(mNmrDates, 
                                         &NbNmr, 
                                         DRL_TDATE_T, 
                                         TRUE));
        for (i = 0; i < NbNmr; i++)
        {
            NmrDates.push_back(mNmrDates[i]); 
            KMrNTree::Insert(mNmrDates[i]);
        }

	    DrlTDateVectFree(mNmrDates, 0, NbNmr);

             
    }
    catch(KFailure)
    {
        if (mNmrDates != NULL)
            DrlTDateVectFree(mNmrDates, 0, NbNmr);
        throw KFailure("%s: failed\n", routine);
    }
}


void KTMXirTree::SetUpCritDate()
{
    static char routine[] = "KTMXirTree::SetUpCritDate";

    try
    {
        // Insert additional Critical dates
        KVector(TDate)::iterator iterDate;
        for(KVector(int)::iterator iterCV = mCVTypes.begin();
            iterCV != mCVTypes.end(); ++iterCV)
            {
                for (iterDate = GetZeroBank(*iterCV).mEvDates.begin();
                iterDate != GetZeroBank(*iterCV).mEvDates.end();
                ++iterDate)
                {
                    Insert(*iterDate);
                    Insert(GetZeroBank(*iterCV).mErDates[*iterDate]);
                }
            }

	    IF_FAILED_THROW(DrlVTypeVectSort(TPDates, 
                                         &NbTP, 
                                         DRL_TDATE_T, 
                                         TRUE));      

        // Save CritDates
        for (int i = 0; i < NbTP; i++)
        {
            CritDates.push_back(TPDate(i));
//            cout << routine << GtoFormatDate(TPDates[i]) << endl;            
        }
    }
    catch(KFailure)
    {
        throw KFailure("%s: failed\n", routine);
    }
}

void KTMXirTree::PrintTPDate()
{
    static char routine[] = "KTMXirTree::PrintTPDate";

    try
    {
        KVector(TDate)::const_iterator iDATE;
        cout << "TPDate " << endl;
        for (int i = 0; i <= NbTP; i++)
        {
            cout << i << "\t" << DrlTDate2DrDate(TPDates[i]) << endl;
        }
        cout << "LastDate " << DrlTDate2DrDate(mLastDate) << endl;

        cout << "CritDate " << endl;
        for (iDATE = CritDates.begin();
             iDATE != CritDates.end();
             ++iDATE)
             {
                 cout << DrlTDate2DrDate(*iDATE) << endl;
             }
        cout << "NmrDate"<< endl;
        for (iDATE = NmrDates.begin(); iDATE != NmrDates.end(); iDATE++)       
        {
            cout << DrlTDate2DrDate(*iDATE) << endl;
        }

    }
    catch(KFailure)
    {
        throw KFailure("%s: failed\n", routine);
    }
}


/* -----------------------------------------------------------------
 * 
 * Setup Nmr dates.  Nmr_Schedule
 */

void
KTMXirTree::InsertNmrDates(KVolDiag& BmkVol)
{
    static char routine[] = "KTMXirTree::InsertNmrDates";
    
    long        LastProdDate;
    long        nmrDate;
    int         i;

    try
    {
        LastProdDate = DrlTDate2DrDate(TPDate(NbTP-1));
//        cout << routine << " valuedate = " << treeValueDate << '\t'
//             << "lastproddate=" << LastProdDate << endl;
        IF_FAILED_THROW( Nmr_Schedule (
                                treeValueDate,
                                LastProdDate,
                                mktvol_data));

        for (i = 0; i < mktvol_data->NbNmr; i++)
        {
            nmrDate = DrlDrDate2TDate(mktvol_data->NmrDate[i]);
//            cout << routine << "insert nmr date " << GtoFormatDate(nmrDate) << endl;
            NmrDates.push_back(nmrDate); 
            KMrNTree::Insert(nmrDate);
        }
    }
    catch(KFailure)
    {
        throw KFailure("%s: failed\n", routine);
    }
}

void
KTMXirTree::SetUpTimeline()
{
static	char	routine[] = "KTMXirTree::SetUpTimeline";

  try 
  {
    if(TMXFlag == TRUE)
    {
	    // After all critical dates are inserted,
	    // sort and merge the critical dates in ascending order
	    //
	    IF_FAILED_THROW(DrlVTypeVectSort(TPDates, &NbTP, DRL_TDATE_T, TRUE));

	    // Last critical date
	    mLastDate = TPDate(NbTP-1);


	    // Add all the benchmark vol calibration dates which are less than
	    // the last critical dates.
	    for(KVector(TDate)::iterator iterVolDate=mVolDates.begin();
				iterVolDate != mVolDates.end(); ++iterVolDate)
	    {
		    if (*iterVolDate < mLastDate)
			    Insert(*iterVolDate);
	    }

	    // Set up timeline based on all the critical dates and PPY.
	    InitializeTimeline_TMX(mEoI, mPPY);
    }
    else
    {
        KPirTree::SetUpTimeline();
    }

  }
  catch(KFailure) 
  {
	throw KFailure("%s: failed.\n", routine);
  }
    
}

void
KTMXirTree::InitializeTimeline_TMX(
    char EoI,               // (I) Equal or increasing time steps
    int  ppy)               // (I) period per year
{
static	char	routine[] = "KTMXirTree::InitializeTimeline_TMX";

    CRIT_DATE   *TreeCritDate = NULL;
    int         NbCritDate = 0;
    int         i;
    long        CurrentNmrDate;

    try
    {
    	// Sort and merge the critical dates in ascending order
	    //
	    IF_FAILED_THROW(DrlVTypeVectSort(TPDates, &NbTP, DRL_TDATE_T, TRUE));
//        cout << "ppy=" << ppy << " EoI = " << EoI << endl;
        // 
        ASSERT_OR_THROW(NbTP > 0);
     //   tree_data->NbTP = NbTP;
     //   Tree_Alloc(tree_data);

        TreeCritDate = (CRIT_DATE *) DR_Array (CRITDATE,0,0);
        ASSERT_OR_THROW(TreeCritDate !=NULL);

        IF_FAILED_THROW(Add_To_DateList(&NbCritDate,
                                        &TreeCritDate,
                                         treeValueDate,
                                         NBCRITDATE,
                                         0,0,0,0,0,
                                         0,0,0));

        // Insert CritDate
        for (i = 0; i < CritDates.size(); i++)
        {
            tree_data->LastProdDate = DrlTDate2DrDate(CritDate(i));
            IF_FAILED_THROW(Add_To_DateList(&NbCritDate,
                                            &TreeCritDate,
                                            tree_data->LastProdDate,
                                            3, 
                                            0,0,0,0,0,
                                            0,0,0));           
        }

        // Insert Nmr dates
//        cout << "nmrsize=" << NmrDates.size() << endl;
        for (i = 0; i < NmrDates.size(); i++)
        {
            CurrentNmrDate = DrlTDate2DrDate(NmrDate(i));
//            cout << routine << "nmr date " << CurrentNmrDate << endl;
            IF_FAILED_THROW(Add_To_DateList(&NbCritDate,
                                            &TreeCritDate,
                                            CurrentNmrDate,
                                            NMREVENT,
                                            i, 
                                            0,0,0,0,
                                            0,0,0));
        }
 
        for (i = 0; i<NBCRITDATE; i++ )
        {
            tree_data->CritType[i] = 'D';
            tree_data->NbZeros[i] = 0;
        }

        IF_FAILED_THROW(Sort_CritDate( NbCritDate,
                                       TreeCritDate));
//        cout << routine << "last crit date=" << TreeCritDate[NbCritDate-1].CritDate << endl;
        // Call TMX time line function to generate tree time line
        IF_FAILED_THROW( Time_Line(
                         treeValueDate,
                         NbCritDate,
                         TreeCritDate,
                         EoI,
                         tree_data));


        // free initial critical dates array
        // so that we can replace it.
        DrlVTypeVectFree((void*) TPDates, NbTP, DRL_TDATE_T);
        TPDates = NULL;

        // Copy TMX tree time line to Basis tree
        NbTP    = tree_data->NbTP;
        NbTPMax = MAX(NbTP+2, NbTPMax);

        ASSERT_OR_THROW((TPDates = DrlTDateVectAlloc(0, NbTP+1)) != NULL);
        ASSERT_OR_THROW((TPTimes = DrlDoubleVectAlloc(0, NbTP+1)) != NULL);

        ASSERT_OR_THROW((Length = DrlDoubleVectAlloc(-1, NbTP+1)) != NULL);
        ASSERT_OR_THROW((LengthJ = DrlDoubleVectAlloc(-1, NbTP+1)) != NULL);

        for (i = 0; i <=NbTP; i++)
        {
            TPDates[i] = DrlDrDate2TDate(tree_data->TPDate[i]);
            TPTimes[i] = (TPDates[i] - TPDates[0]) / 365.;

            //dppLog << tree_data->TPDate[i] << " type= " << tree_data->TPtype[NMREVENT][i] << endl;
        }

        // Copy the time step size (Length) and
        // time jump size (LengthJ)
        for (i = -1; i < NbTP+1; i++)
        {
            Length[i]  = tree_data->Length[i];
            LengthJ[i] = tree_data->LengthJ[i];
        }

        TMX::Free_DR_Array(TreeCritDate,    CRITDATE,   0,  NbCritDate-1);

    }
    catch(KFailure) 
    { 
        TMX::Free_DR_Array(TreeCritDate,    CRITDATE,   0,  NbCritDate-1);
        throw KFailure("%s: failed.\n", routine);
    }
}

void KTMXirTree::WrapTreeTimeLine (KVolDiag &irVolDiag,
                                   KMrParam &irMrParam)
{
    static char routine[] = "KTMXirTree::WrapTreeTimeLine";

    try{

        if (TMXFlag == TRUE)
        {
            SetUpCritDate();
            InsertNmrDates(irVolDiag);
            KTMXirTree::SetUpTimeline();
    //        PrintTPDate();
        }
        else
        {
            KPirTree::SetUpTimeline();
        }
      }
    catch (KFailure) {
	    throw KFailure ("%s: failed.\n", routine);
    }
}

void KTMXirTree::PrintTreeTimeLine()
{
    static char routine[] = "KTMXirTree::PrintTreeTimeLine";
    int i, TmxNbTP;

    ofstream fstream("TREE.PRN", ios::out);
    
    if (!fstream)
    {
        cerr << "TERM.prn could not be opened" << endl;
        exit(1);
    }

    try{
        TmxNbTP = tree_data->NbTP;
        ASSERT_OR_THROW(TmxNbTP > 0);

        fstream << " Tp " << endl;
        for (i = 0; i <= TmxNbTP; i++)
        {
            fstream << setw(3) << i << setw(10) << tree_data->TPDate[i]
                    << "\tNmrtype " << tree_data->TPtype[NMREVENT][i]  
                    << "\tCritype " << tree_data->TPtype[3][i] << endl;
        }
        fstream << endl;
    }
    catch (KFailure) {
	throw KFailure ("%s: failed.\n", routine);
    }
}

void 
KTMXirTree::PackCurves( const KMarketCurves& MktCurves, 
                        const String& CvDiscName)
{
    static char routine[] = "KTMXirTree::PackCurves";

    try
    {
        KZCurve Zcurve;
   //     long    AlibDcc;
        
        /* Pack the diffusion curve*/
        KMap(int, String)::const_iterator iterName = MktCurves.mCVNames.find(KV_DIFF);
        KMap(int, KZCurve)::const_iterator iterCV = MktCurves.mCV.find(KV_DIFF);
        if( iterName == MktCurves.mCVNames.end() ||
            iterCV == MktCurves.mCV.end() )
        {
            throw KFailure("%s: no diffusion curve\n", routine);
        }
        Zcurve = iterCV->second;
        strcpy(diffCvName, iterName->second.c_str());

        // Assign tree_data->CvDiff
        iterName = MktCurves.mZCFiles.find(KV_DIFF);
 //       cout << iterName->second.c_str() << endl;
        if ( strcmp(iterName->second.c_str(),"zero.dat") == 0 )
        {
            tree_data->CvDiff = 0;
        }
        else if (strcmp(iterName->second.c_str(), "disczero.dat") == 0 )
        {
            tree_data->CvDiff = 1;
        }
        else
        {
            tree_data->CvDiff = 2;
        }
 //       cout << "Cvdiff: " << tree_data->CvDiff << endl;

        if (Zcurve.ZeroInterpType() == GTO_LINEAR_INTERP)
        {
            TMX::ZeroInterpTypeFlag = 0;
        }
        else if (Zcurve.ZeroInterpType() == GTO_FLAT_FORWARDS)
        {
            TMX::ZeroInterpTypeFlag = 1;
        }
        else
        {
            throw KFailure("%s: TMX only allows linear and flat zcurve "
                           "interpolation\n", routine);
        }
//        cout << "zcurve interp type " << TMX::ZeroInterpTypeFlag << endl;
        // In default, all curves are same as the diffusion curve.
        PackSingleDRCurve(Zcurve,0);
        PackSingleDRCurve(Zcurve,1);
        PackSingleDRCurve(Zcurve,2);

        // Overwrite the discounting curve.
        int CvDiscIdx;
        for (iterName = MktCurves.mCVNames.begin();
             iterName != MktCurves.mCVNames.end();
             iterName++)
             {
                 if (iterName->second.c_str() == CvDiscName)
                 {
                     CvDiscIdx = iterName->first;
                     break;
                 }
             }
        if (iterName == MktCurves.mCVNames.end())
            throw KFailure("%s: invalide disc curve name\n");
        iterName = MktCurves.mZCFiles.find(CvDiscIdx);
        iterCV = MktCurves.mCV.find(CvDiscIdx);
        Zcurve = iterCV->second;
//        cout << iterName->second.c_str() << endl;
        if ( strcmp(iterName->second.c_str(),"zero.dat") == 0 )
        {
            tree_data->CvDisc = 0;
        }
        else if (strcmp(iterName->second.c_str(), "disczero.dat") == 0 )
        {
            tree_data->CvDisc = 1;
        }
        else
        {
            tree_data->CvDisc = 2;
        }
//        cout << "CvDisc: " << tree_data->CvDisc << endl;
        PackSingleDRCurve(Zcurve,tree_data->CvDisc);

        
    }
    catch (KFailure) {
	throw KFailure ("%s: failed.\n", routine);
    }
}

void
KTMXirTree::PackSingleDRCurve(const KZCurve& Zcurve, int CvIdx)
{
    static char routine[] = "KTMXirTree::PackSingleDRCurve";
    try
    {
        int     i;

        /* curve base date */
        t_curve[CvIdx].ValueDate = DrlTDate2DrDate(Zcurve.BaseDate());
        t_curve[CvIdx].Today     = t_curve[CvIdx].ValueDate;
        t_curve[CvIdx].SpotDays  = 0;
        /* Frequency */
        switch(Zcurve.Freq())
        {
        case 1:
             t_curve[CvIdx].SwapFreq = 'A';
             break;
        case 2:
             t_curve[CvIdx].SwapFreq = 'S';
             break;
        case 4:
             t_curve[CvIdx].SwapFreq = 'Q';
             break;
        default:
             printf("Error: invalid frequency\n");
             SUCCESS_OR_THROW(FAILURE);
        }

        
         /* DCC */
        // The Dcc is hard coded as GTO_ACT_365F in dpp.
        t_curve[CvIdx].MMB = 360;  
        strcpy(t_curve[CvIdx].SwapDCC, "ACT");

        /* dates and rates  */
        t_curve[CvIdx].NbZero = Zcurve.NumItems();
        for (i = 0; i < t_curve[CvIdx].NbZero; i++)
        {
             (t_curve[CvIdx]).ZeroDate[i] = DrlTDate2DrDate(
                                                Zcurve.Date(i));
             (t_curve[CvIdx]).Zero[i]     = Zcurve.Rate(i);
        }
        long LastDate = TMX::Nxtmth(t_curve[CvIdx].ValueDate, 1200L, 1L);
            
        if (t_curve[CvIdx].ZeroDate[t_curve[CvIdx].NbZero-1] < LastDate)
        {
            if ( t_curve[CvIdx].NbZero < MAXNBDATE)
            {
                t_curve[CvIdx].NbZero++;
                t_curve[CvIdx].ZeroDate[t_curve[CvIdx].NbZero-1] = LastDate;
                t_curve[CvIdx].Zero[t_curve[CvIdx].NbZero-1] = 
                        t_curve[CvIdx].Zero[t_curve[CvIdx].NbZero-2];
            }
            else
            {
                t_curve[CvIdx].ZeroDate[t_curve[CvIdx].NbZero-1] = LastDate;
            }
        }

    }
    catch (KFailure) 
    {
	    throw KFailure ("%s: failed.\n", routine);
    }

}


bool
KTMXirTree::IsCalibIdxBase(const KVolDiag& volDiag)
{
    static char routine[] = "KTMXirTree::IsCalibIdxBase";
    try
    {
        KVector(double)::const_iterator iter = volDiag.mVolMats.begin();
        ASSERT_OR_THROW ( iter != volDiag.mVolMats.end());
        if (*iter < 1.0)
            return true;
        else
            return false;
    }
    catch (...)
    {
        throw KFailure("%s: failed\n", routine);    	
    }
}


void KTMXirTree::PackVolData( const KVolDiag& volDiag,
                              const KMarketCurves& MktCurves,
                              const KSmileParam& SmileParam,
                              const KMrParam& MrParam)
{
    static char routine[] = "KTMXirTree::PackVolData";

    long        baseDate;
    long        SmlExp;
    int         i, j;
    int         ILiqBmk;
    int         nbfactor = tree_data->NbFactor;

    baseDate = DrlTDate2DrDate(MktCurves.mToday);

    try
    {
        mktvol_data->BaseDate = baseDate;
        mktvol_data->NbVol    = volDiag.Size();

        mktvol_data->CalibSmileFlag = TRUE;
        mktvol_data->SkipFlag       = FALSE;
     
        if (IsCalibIdxBase(volDiag))
        {
            if (t_curve[0].MMB == 365)
            {
                mktvol_data->DCC = '5';
            }
            else 
            {
                mktvol_data->DCC = '0';
            }

        }
        else
        {
            if (!strcmp(t_curve[0].SwapDCC, "360"))
            {
                mktvol_data->DCC = '0';
            }
            else if (!strcmp(t_curve[0].SwapDCC, "365"))
            {
                mktvol_data->DCC = '5';
            }
            else
            {
                mktvol_data->DCC = '3';
            }
        }

        KVector(int)::const_iterator itFREQ = volDiag.mVolFreqs.begin();
        switch (*itFREQ)
        {
        case 1:
            mktvol_data->Freq = 'A';
            break;
        case 2:
            mktvol_data->Freq = 'S';
            break;
        case 4:
            mktvol_data->Freq = 'Q';
            break;
        case 12:
            mktvol_data->Freq = 'M';
            break;
        default:
            SUCCESS_OR_THROW(FAILURE);
        }

        i = 0;
        for (KVector(TDate)::const_iterator itDATE = volDiag.mVolDates.begin();
             itDATE != volDiag.mVolDates.end();
             ++itDATE)
             {
                 mktvol_data->VolDate[i] = DrlTDate2DrDate(*itDATE);
                 mktvol_data->SwapSt[i]  = mktvol_data->VolDate[i];
                 i++;
             }
        ASSERT_OR_THROW(i == mktvol_data->NbVol);
        i = 0;
        for (KVector(double)::const_iterator itRATE = volDiag.mVolRates.begin();
             itRATE != volDiag.mVolRates.end();
             ++itRATE)
             {
                 mktvol_data->Vol[0][i]      = (*itRATE);
                 mktvol_data->VolUsed[i]     = TRUE;
                 mktvol_data->SmlLiqDate[i]  = 0;
                 i++;
             }
        ASSERT_OR_THROW(i == mktvol_data->NbVol);

        i = 0;
        for (KVector(double)::const_iterator itMAT = volDiag.mVolMats.begin();
             itMAT != volDiag.mVolMats.end();
             ++itMAT)
             {
                 mktvol_data->SwapMatMos[i] = static_cast<int> ((*itMAT)*12 + 0.5);
                 mktvol_data->SwapMat[i] = Nxtmth(mktvol_data->SwapSt[i],
                                                  mktvol_data->SwapMatMos[i],
                                                  1);
                 i++;
             }
        ASSERT_OR_THROW(i == mktvol_data->NbVol);

        /*Convert bp vol to % vol*/
        /*if (volDiag.mVolType == NORMAL)
        {
            ConvertBpVolToPercVol();
        }*/
        /* CetNbIter*/
        mktvol_data->CetNbIter = SmileParam.mNumIter;

        /* MR parameters*/
        for (i = 0; i < nbfactor; i++)
        {
            mktvol_data->Beta[i]  = MrParam.mBeta[i];
            mktvol_data->Alpha[i] = MrParam.mAlpha[i];
        }
        for (i = 0; i < nbfactor*(nbfactor-1)/2; i++)
        {
            mktvol_data->Rho[i]   = MrParam.mRho[i];
        }

        /* Initialize unused variables */
        for (i = nbfactor; i < 3; i++) 
        {
            mktvol_data->Beta[i]  = -999;
            mktvol_data->Alpha[i] = -999;
        }
        for (i = nbfactor*(nbfactor-1)/2; i < 3; i++)
        {
            mktvol_data->Rho[i] = -999;
        }


        /* Backbone */
        if ( !IS_EQUAL(MrParam.mBackboneQ, 1.))
            throw KFailure("%s failed: TMX only accepts lognormal backbone %f\n",
                            routine, MrParam.mBackboneQ);
        mktvol_data->Bbq = 1.;

        if ( IS_EQUAL(MrParam.mSmoothFact, 1.))
            mktvol_data->SmoothingFlag = 'Y';
        else
            mktvol_data->SmoothingFlag = 'N';
//        cout << "smooth = " << MrParam.mSmoothFact << endl;
//        cout << mktvol_data->SmoothingFlag << endl;

       /* mktvol_data->Bbq = MrParam.mBackboneQ;
        ASSERT_OR_THROW( IS_EQUAL(mktvol_data->Bbq, 1) ||
                         IS_EQUAL(mktvol_data->Bbq, 0));
        double norm = 0.;
        for (i = 0; i < MrParam.mNumFact; i++)
        {
            norm += mktvol_data->Alpha[i] * mktvol_data->Alpha[i];
        }
        norm = sqrt(norm);
        if ( IS_EQUAL(mktvol_data->Bbq, 1))
        {
            mktvol_data->VolNorm = 0.;
            mktvol_data->VolLogn = norm;
        }
        else 
        {
            mktvol_data->VolNorm = norm;
            mktvol_data->VolLogn = 0.;
        }*/

        /* Liquid Benchmark Expiries */
        KVector(TDate)::const_iterator itDate;
        KVector(double*)::const_iterator itSmilePar;
        for (itDate = SmileParam.mLiqBmkDates.begin(), 
                itSmilePar = SmileParam.mMQSmile.begin(); 
             itDate != SmileParam.mLiqBmkDates.end() &&
                itSmilePar != SmileParam.mMQSmile.end(); 
             itDate++, itSmilePar++)
        {
            SmlExp = DrlTDate2DrDate(*itDate) ;
            if (SmlExp <= mktvol_data->SwapSt[mktvol_data->NbVol-1])
            {
                ILiqBmk = TMX::GetDLOffset (mktvol_data->NbVol,
                                       mktvol_data->SwapSt,
                                       SmlExp,
                                       CbkEXACT);
                ASSERT_OR_THROW( ILiqBmk >=0 );
                mktvol_data->SmlLiqDate[ILiqBmk] = 1;
            }
            
        }

        /* Copy MultiQ Smile */
        for (itSmilePar = SmileParam.mMQSmile.begin(), j = 0; 
             itSmilePar != SmileParam.mMQSmile.end() && 
             j < SmileParam.mSmlDates.size(); 
             itSmilePar++, j++)
        {
            for (i = 0; i < NBMQSMLPAR; i++)
            {
                mktvol_data->Vol[i+1][j] = (*itSmilePar)[i];
            }
        }
        mktvol_data->NbSigmaMQ = SmileParam.MQNormT;
        mktvol_data->NckMQ     = SmileParam.MQNck;

        /* convert vov*/
        IF_FAILED_THROW( Convert_VoV(mktvol_data, t_curve) );
        /* check */
        IF_FAILED_THROW( MktVol_Check_W(mktvol_data) );

    }
    catch(...){
        throw KFailure("%s: failed\n", routine);
    }    

    
}

void
KTMXirTree::ConvertBpVolToPercVol()
{
    static char routine[] = "KTMXirTree::ConvertBpVolToPercVol";
    try
    {
        int     i;
        double  fwdRate, annuity;

        for (i = 0; i < mktvol_data->NbVol; i++)
        {
            IF_FAILED_THROW(Par_Yield_From_Dates(
                            &fwdRate,
                            &annuity,
                            mktvol_data->SwapSt[i],
                            mktvol_data->SwapMat[i],
                            mktvol_data->DCC,
                            mktvol_data->Freq,
                            'F',
                            t_curve[tree_data->CvDisc].NbZero,
                            t_curve[tree_data->CvDisc].Zero,
                            t_curve[tree_data->CvDisc].ZeroDate,
                            t_curve[tree_data->CvDisc].ValueDate));           
            ASSERT_OR_THROW(fwdRate > TINY);
            mktvol_data->Vol[0][i] /= fwdRate;
        }
    }
    catch(...)
    {
        throw KFailure("%s: failed\n", routine);
    }    
}

void KTMXirTree::PrintEnv()
{
    static char *CurveName[3] = { "basis", "zero", "risk"};
    ofstream fstream("TERM.prn", ios::out);
    
    if (!fstream)
    {
        cerr << "TERM.prn could not be opened" << endl;
        exit(1);
    }

    int i, j;

    fstream << "today: " << treeToday << endl;
    fstream << "valuedate: " << treeValueDate << endl;
    /* Print curves */
    for (i = 0; i < 3; i++)
    {
        fstream << "curve: " << CurveName[i] << endl;
        fstream << "# Start date " << endl;
        fstream << t_curve[i].ValueDate << endl;
        fstream << "# Money Market basis " << endl;
        fstream << t_curve[i].MMB << endl;
        fstream << "# Annual or semi-annual curve" << endl;
        fstream << t_curve[i].SwapFreq << endl;
        fstream << "# Year basis for benchmark swaps" << endl;
        fstream << t_curve[i].SwapDCC << endl;
        fstream << "# Number of entries" << endl;
        fstream << t_curve[i].NbZero << endl;
        fstream << "#zero maturity and rates" << endl;
        for (j = 0; j < t_curve[i].NbZero; j++)
        {
            fstream << setw(8) << (t_curve[i]).ZeroDate[j] << "\t"
                << setprecision(13) << (t_curve[i]).Zero[j] << endl;
        }
        fstream << endl;
    }

    /* IR vol */
    fstream << "IR vol: "<< endl;
    fstream << "vol base date: " << mktvol_data->BaseDate << endl;
    fstream << "Nb vol: " << mktvol_data->NbVol << endl;
    fstream << "Freq: " << mktvol_data->Freq << endl;
    fstream << "VolDate SwapSt SwapMat SwapMatMos VolUsed SmlLiqDate Vol\t\t\t Smile"<< endl;
    for (i = 0; i < mktvol_data->NbVol; i++ )
    {
        fstream << setw(8) << mktvol_data->VolDate[i] << "\t" 
                << setw(8) << mktvol_data->SwapSt[i] << "\t"
                << setw(8) << mktvol_data->SwapMat[i] << "\t"
                << setw(4) << mktvol_data->SwapMatMos[i] << "\t"
                << setw(2) << mktvol_data->VolUsed[i] << "\t"
                << setw(2) << mktvol_data->SmlLiqDate[i]<< "\t"
                << setw(9) << setprecision(5) 
                << mktvol_data->Vol[0][i] << "\t";
        for (j = 1; j <= NBMQSMLPAR; j++)
        {
            fstream << setw(9) << setprecision(5) 
                    << mktvol_data->Vol[j][i] << "\t";
        }
        fstream << endl;           
    }
    /* Constants */
    fstream << "CetNbIter: " << mktvol_data->CetNbIter << endl;
    fstream << "Alpha: " << mktvol_data->Alpha[0] << endl;
    fstream << "Beta: " << mktvol_data->Beta[0] << '\t' << mktvol_data->Beta[1] << '\t' << mktvol_data->Beta[2]  << endl;
    fstream << "NbSigmaMQ: " << mktvol_data->NbSigmaMQ << endl;
    fstream << "NCKMQ: " << setw(15) << mktvol_data->NckMQ << endl;
    //fstream << "DeltaNMQ: " << mktvol_data->DeltaNMQ << endl;
    //fstream << "TauNMQ: " << mktvol_data->TauNMQ << endl;
    fstream << "Backbone: " << mktvol_data->Bbq << endl;
    //fstream << "VolNorm = " << mktvol_data->VolNorm << " VolLogn = " 
    //        << mktvol_data->VolLogn << endl;
    //fstream << "DeltaK: " << endl;
   /* for (i=0; i < NBSTRIKE; i++ )
    {
        fstream << "mktvol_data.Delta[" << i << "]=" << mktvol_data->Delta[i] << endl;
    }*/
    fstream << endl;

    /* Mr tree mVolDates */
    fstream << "MrTree vol dates" << endl;
    for (KVector(TDate)::iterator itDate = mVolDates.begin();
        itDate != mVolDates.end();
        itDate++)
        {
            fstream << GtoFormatDate(*itDate) << endl;
        }
    fstream << endl;

}


void KTMXirTree::setToday(long date)
{
    treeToday = DrlTDate2DrDate(date);
}

void KTMXirTree::setValueDate(long date)
{
    treeValueDate = DrlTDate2DrDate(date);
}

long KTMXirTree::getToday()
{
    return treeToday;
}

long KTMXirTree::getValueDate()
{
    return treeValueDate;
}

