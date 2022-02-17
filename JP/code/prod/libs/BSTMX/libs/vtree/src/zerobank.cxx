/***************************************************************
 * Module:        BasisTree
 * Submodule:        
 * File:        
 * Function:        
 * Author:        David Liu
 ***************************************************************/
#include "kstdinc.h"    /* Standard definitions & error hadling */
#include "kutilios.h"

#include "kzrbank.h"
#include "kpirtree.h"

//--------------------------------------------------------------
//
KZeroBank::KZeroBank() 
          : KCBank()
{
static        char        routine[] = "KZeroBank::KZeroBank";

	mZeroInterpType = GTO_LINEAR_INTERP;

}


//--------------------------------------------------------------
//
KZeroBank::KZeroBank(KPirTree           &irTree, 
                     const String &zbName,
                     const String &curveName,
		     const int    interpType)
          : KCBank(irTree, zbName)
{
static        char        routine[] = "KZeroBank::KZeroBank";

        SetCurveIdx(curveName);

	mZeroInterpType = interpType;
}


//---------------------------------------------------------------
// Copy constructor.
 
KZeroBank::KZeroBank(const KZeroBank &zb)
          :KCBank(zb)
{
        mCurveIdx = zb.GetCurveIdx();
	mZeroInterpType = zb.ZeroInterpType();
}
 


//--------------------------------------------------------------
//
KZeroBank::~KZeroBank()
{
static  char    routine[] = "KZeroBank::~KZeroBank";

}


//---------------------------------------------------------------^
// Assignment operator.
KZeroBank&
KZeroBank::operator=(const KZeroBank &zb)
{
        *this = zb;
        mCurveIdx = zb.GetCurveIdx();
	mZeroInterpType = zb.ZeroInterpType();

        return(*this);
}



//--------------------------------------------------------------
//
void
KZeroBank::Update(int tpIdx)           // (I) Index of time point 
{
static        char        routine[] = "KZeroBank::Update";

        TDate        currDate = (mVTree->TPDate)(tpIdx);
        
 try{
        // compact and dev the zero bank 
        KCBank::Update(tpIdx);

        // Add a new zero slice if is a zero maturity date
        KVector(TDate)::iterator iterEvDate =
                find(mEvDates.begin(), mEvDates.end(), currDate);

        if (iterEvDate != mEvDates.end()) 
        {
            KTSlice *newSlice = new KTSlice(*mVTree,
                                            format("Zero slice with maturity "
                                                  "on %s", 
                                                  GtoFormatDate(currDate)),
                                            GetCurveName());

        // IR or risky Dev vs. protection Dev 
        if (GetCurveIdx() == KV_CREDIT_PROT_DEFRECOV ||
            GetCurveIdx() == KV_CREDIT_PROT_BINRECOV)
                *newSlice = 0.0;
        else
                *newSlice = 1.0;

            InsertSlice(currDate, newSlice);
        }

    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }


}




//--------------------------------------------------------------
// Only linear and flat interpolations are supported
//
void
KZeroBank::GetZeroSlice(KTSlice &ts, TDate  EvDate)
{
static  char    routine[] = "KZeroBank::GetZeroSlice";
 
        double        dcfLeft, dcfRight, dcfEv;
        double  tFactor;

        TDate        EvDateL, EvDateR, currDate;

 try{

    currDate = TPDateCurrent();

    if (currDate == EvDate)
    {
        // Discount Factor (riskless or risky)
        if (GetCurveIdx() != KV_CREDIT_PROT_DEFRECOV)    
            ts = 1e0;
        else  // Default protection
            ts = 0e0;
    }
    else
    {
        if (ZeroInterpType() != GTO_LINEAR_INTERP &&
            ZeroInterpType() != GTO_FLAT_FORWARDS)
            throw KFailure("%s: unsupported interpolation type (%ld)"
                           " for zero bank %d.\n",
                           routine,
                           ZeroInterpType(),
			   GetCurveIdx());   

        if (mEvDates.empty())
        {
            throw KFailure("%s: zero bank (%s) is empty.\n",
                           routine,
                           GetCurveName().c_str());
        }

        if (mSlices.count(EvDate))    // Exact date match
                GetSlice(ts, EvDate);
        else {          // linear interpolate from two adjacent slices
            // Lower bound
            KVector(TDate)::iterator iterEvRight =
                        lower_bound(mEvDates.begin(), mEvDates.end(), EvDate);
        
            if(iterEvRight == mEvDates.end()) // out of bound 
                  throw KFailure("%s: Failed to interpolate.\n"
                                 "EvDate (%s) is larger than the maximum "
                                 "active date (%s) in the zero bank %s.\n",
                                 routine, 
                                 GtoFormatDate(EvDate),
                                 GtoFormatDate(mEvDates.back()),
                                 GetCurveName().c_str());

            EvDateR = *iterEvRight;
            EvDateL = MAX(*(iterEvRight-1), currDate);

            tFactor = (double)(EvDate  - EvDateL)/
                          (double)(EvDateR - EvDateL);

            KTSlice tmpSliceR(*mVTree,
                              format("Zero slice with maturity "
                                     "on %s",
                                     GtoFormatDate(EvDateR)),
                                     GetCurveName());

            KTSlice tmpSliceL(*mVTree,
                              format("Zero slice with maturity "
                                     "on %s",
                                     GtoFormatDate(EvDateL)),
                                     GetCurveName());


            // Riskless or risky discount factor
            //
            if (GetCurveIdx() != KV_CREDIT_PROT_DEFRECOV)    
            {

                GetSlice(tmpSliceR, EvDateR);

                KPirTree *pirTree = dynamic_cast<KPirTree*>(mVTree);
                KDayCc   dcc = (pirTree->GetKZCurve(GetCurveIdx())).DayCc();
                double   basis = (pirTree->GetKZCurve(GetCurveIdx())).Basis();

                // day count fraction from EvDate date to current date
                dcfEv = DayCountFraction(currDate, 
                                         EvDate, 
                                         dcc);


                // day count fraction from right date to current date
                dcfRight = DayCountFraction(currDate, 
                                            EvDateR, 
                                            dcc);

                // 
                //  Interpolate using EvDateL and EvDateR
                // 
                //  == curr === EvDateL ===== EvDate ===== EvDateR ======
                // 
                if (EvDateL > currDate)
                {
                    GetSlice(tmpSliceL, EvDateL);

                    // day count fraction from left date to current date
                    dcfLeft = DayCountFraction(currDate, 
                                           EvDateL, 
                                           dcc); 

                    switch (ZeroInterpType()) {
                    case GTO_LINEAR_INTERP:

                        // Convert left zero price to "1 + zrate/basis" slice
                        mVTree->TSliceScalarOper(tmpSliceL,
                                                     -1.0/(basis*dcfLeft),
                                                 POW);

                        mVTree->TSliceScalarOper(tmpSliceR,
                                                 -1.0/(basis*dcfRight),
                                                 POW);

                        // Linear interpolate on the rate 
                        // (actually interpolating on "1 + zrate/basis")
                        //
                        tmpSliceL *= (1. - tFactor);

                        tmpSliceR *= tFactor;

                        tmpSliceL += tmpSliceR;

                        // Convert back to zero price at EvDate from 
                        // "1 + zrate/basis"
                        //
                        mVTree->TSliceScalarOper(tmpSliceL,
                                                 -(basis*dcfEv),
                                                 POW);
                        break;

                    case GTO_FLAT_FORWARDS:
                        tmpSliceR /= tmpSliceL;
                        mVTree->TSliceScalarOper(tmpSliceR,
                                                 tFactor,
                                                 POW);
                        tmpSliceL *= tmpSliceR;
                        break;
                    }        

                        // Copy the results
                    ts = tmpSliceL;

                }
                //  If there is no zero slice between current and EvDate,
                //  then use current date to interpolate on the right.
                // 
                //  == EvDateL === curr ===== EvDate ===== EvDateR ======
                else          // EvDateL <= currDate
                {
                    // Extroplate using flat zero rate from EvDateR
                    //
                    switch (ZeroInterpType()) {
                    case GTO_LINEAR_INTERP:
                        // Convert left zero price to "1 + zrate/basis" slice
                        mVTree->TSliceScalarOper(tmpSliceR,
                                                 -1.0/(basis*dcfRight),
                                                 POW);

                        // Convert back to zero price at EvDate from 
                        // "1 + zrate/basis"
                        //
                        mVTree->TSliceScalarOper(tmpSliceR,
                                                 -(basis*dcfEv),
                                                 POW);
                        break;

                    // Extroplate using flat forward rate EvDateR
                    //
                    case GTO_FLAT_FORWARDS:
                        mVTree->TSliceScalarOper(tmpSliceR,
                                                 tFactor,
                                                 POW);
                        
                        break;
                    }

                        // Copy the results
                    ts = tmpSliceR;
                }
            }      // riskless or risky discount factor
            else   // default protection 
            {
                // Risky discount factor
                //
                KPirTree *pirTree = dynamic_cast<KPirTree*>(mVTree);
                KZeroBank& riskyZBank = pirTree->GetZeroBank(KV_CREDIT_RISKY);

                riskyZBank.GetZeroSlice(tmpSliceR,    EvDateR);
                riskyZBank.GetZeroSlice(ts,           EvDate);

                // 
                //  Interpolate using EvDateL and EvDateR
                //  Only FLAT interp is supported
                //
                // Prot(T) = Prot(TL) + {Prot(TR) - Prot(TL)} 
                //         * {Z(TL) - Z(T)} / {Z(TL) - Z(TR)}
                //
                //  == curr === EvDateL ===== EvDate ===== EvDateR ======
                // 
                if (EvDateL > currDate)
                {
                    riskyZBank.GetZeroSlice(tmpSliceL,    EvDateL);

                    tmpSliceR -= tmpSliceL;
                    ts        -= tmpSliceL;
                    ts        /= tmpSliceR;

                    // Protection slices on EvDateL and EvDateR
                    //
                    GetSlice(tmpSliceR, EvDateR);
                    GetSlice(tmpSliceL, EvDateL);
 
                    // Prot(TL,TR)
                    tmpSliceR -= tmpSliceL;

                    //
                    // Flat interpolation formula
                    //
                    ts *= tmpSliceR;
                    ts += tmpSliceL;
                }
                //  If there is no zero slice between current and EvDate,
                //  then use current date to interpolate on the right.
                // 
                // Prot(T) = Prot(TR) * {1 - Z(T)} / {1 - Z(TR)}
                //
                //  == EvDateL === curr ===== EvDate ===== EvDateR ======
                else
                {
                    tmpSliceR -= 1.0;
                    ts        -= 1.0;
                    ts        /= tmpSliceR;

                    // Protection slices on EvDateL and EvDateR
                    //
                    GetSlice(tmpSliceR, EvDateR);
 
                    //
                    // Flat interpolation formula
                    //
                    ts *= tmpSliceR;
                }
            }
        }

    }  // if currDate == EvDate

    }
    catch(KFailure) {
        throw KFailure("%s: failed to get zero slice "
                       "of maturity %s on %s.\n", 
                       routine,
                       GtoFormatDate(EvDate),
                       GtoFormatDate(currDate));
    }   
}





//--------------------------------------------------------------
// Only linear interpolation is supported
//
void
KZeroBank::ComputeParYieldAndAnnuity(
        int tpIdx,                        // (I) tree time point
        const KRateReset &rateReset,        // (I) rate index
        KTSlice        &parYield,                // (O) par yield
        KTSlice        &annuity)                // (O) annuity
{
static  char    routine[] = "KZeroBank::ComputeParYieldAndAnnuity";
 

        TDate        currDate = (mVTree->TPDate)(tpIdx);

        TDate        resetDate = rateReset.ResetDate(),
                effDate          = rateReset.EffDate();

        KVector(KZeroReset)     zeroResetList;
 
        KTSlice                zeroSlice(*mVTree,
                                  format("Temporary zero slice on %s",
                                          GtoFormatDate(currDate)),
                                  GetCurveName());

        double                dcfrn;

 try{
        
        if (resetDate < currDate)
            throw KFailure("%s: resetDate (%s) < current date (%s).\n",
                           routine, 
                           GtoFormatDate(resetDate),
                           GtoFormatDate(currDate));
        else {
            // Create zero dates from KRateReset
            rateReset.GetZeroResetList(zeroResetList);

            if(zeroResetList.size() < 2)  // minimum: eff date + maturity date
            {
                dppLog << endl;
                dppLog << "DEBUG Rate Reset " << rateReset.Rate().CurveName();

                dppLog << endl << endl;
                dppLog << rateReset << endl;
                dppLog << endl;
                for (KVector(KZeroReset)::iterator it=zeroResetList.begin();
                                it != zeroResetList.end(); ++it)
                        dppLog << *it << endl;

                throw KFailure(        "%s: invalid input rateReset. "
                                "Number of zero resets is %d < 2.\n"
                                "Minimum requires effective date and "
                                "maturity date for a simple rate.\n",
                                routine,
                                zeroResetList.size());
            }


            // Copy initial values to parYield and annuity slices
            KVector(KZeroReset)::iterator iterZero=zeroResetList.begin();

            GetZeroSlice(parYield, iterZero->mMaturityDate);
            annuity  = 0.0;
            
            // Start looping from first cash flow date
            ++iterZero;
            for (; iterZero!=zeroResetList.end(); ++iterZero)
            {
                GetZeroSlice(zeroSlice, iterZero->mMaturityDate); 

                if (iterZero==(zeroResetList.end()-1))
                        parYield -= zeroSlice;

                // Day count fraction
                ASSERT_OR_THROW( GtoDayCountFraction(
                        (iterZero-1)->mMaturityDate,
                        iterZero->mMaturityDate,
                        labs((long)(rateReset.Rate().DayCc())),        
                        &dcfrn) == SUCCESS);

                //
                // Risky accrual adjustment
                // Use deterministic default probability
                // This is a quick and dirty HACK!!!!
                //
                if (GetCurveIdx() == KV_CREDIT_RISKY &&
                    rateReset.Rate().DayCc().isNegative())
                {
                    KPirTree *pirTree = dynamic_cast<KPirTree*>(mVTree);
                    double defProb;

                    // Accrual end date
                    defProb = pirTree->GetZeroPrice(KV_CREDIT_RISKY,
                                                    iterZero->mMaturityDate);

                    // Accrual start date
                    defProb /= pirTree->GetZeroPrice(KV_CREDIT_RISKY,
                                                 (iterZero-1)->mMaturityDate);

                    // Use 1/2 approximation for default accrual
                    //
                    defProb = 0.5 * (1.0 - defProb);

                    zeroSlice *= dcfrn * (1. + defProb);
                }
                else
                    zeroSlice *= dcfrn;

                annuity += zeroSlice;
            }
            
            // Compute the par yield.  
            // Apply a cut-off yield level in case annuity is too small.
            // to be implemented later with each slice node operation.
            // Also check for zero annuity node.
            //
            parYield  /= annuity;
            parYield  += rateReset.Rate().Spread();

        }
        
    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }   
}




//--------------------------------------------------------------
//

void
KZeroBank::InitializeZeroDates()
{
static        char        routine[] = "KZeroBank::InitializeZeroDates";

 try{

        // Given all the unprocessed zero maturity dates
        // and corresponding usage dates,  find an "optimal" 
        // set of zero-bank maturities using David Fung's
        // algorithm. The algorithm minimises the number of 
        // zero-bank dates.
/*
        OptimizeDates();

*/

        // Add the final optimized maturity and 
        // early usage dates to the critical dates.
        for(KVector(TDate)::iterator iter =  mEvDates.begin();
                                        iter != mEvDates.end(); ++iter)
        {
                mVTree->Insert(*iter);
                mVTree->Insert(mErDates[*iter]);
        }

    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }

}








/***************************************************************


//--------------------------------------------------------------
//
void
KZeroBank::OptimizeDates()
//        int        NbTP,                // (I) Nb of critical dates
//        TDate   *TPDates,        // (I) critical dates       
//        TDate        today)                // (I) today's date             
{
static        char        routine[] = "KZeroBank::OptimizeDates";

        int     bCritDates = 0;
        long   *CritDates = NULL;

        int     NbMat = 0;
        long   *Mat   = NULL;
        long   *LaMat = NULL;
        long   *ErMat = NULL;

        int     NbZ      = 0;
        long   *Z        = NULL;     // intermediate zerobank mats 
        long   *ErZ      = NULL;     // earliest usage for Z's 
        long   *FirstMat = NULL;

        int     N; // size of unprocessed Mat (Mat[0] .. Mat[N-1]) 
        int     n; // smallest processed idx (Mat[n] .. Mat[N-1] is processed)

 try{
        // basic checks 
        if (mNbMatDates == 0) return;
        if (mNbUseDates != mNbErDates) 
                throw KFailure("%s: mNbMatDates(%ld) != mNbUseDates(%ld).\n", 
                              routine, mNbMatDates, mNbUseDates);

        if ((mZMatDates == NULL) || (mZUseDates == NULL)) 
                throw KFailure("%s: empty mZMatDates or mZUseDates array.\n", 
                              routine);

        if ((mNbEvDate == NULL) ||
            (mNbErDate == NULL) ||
            (mEvDates  == NULL)   ||
            (mErDates  == NULL)) 
                throw KFailure("%s: empty mNbEvDates, mNbErDate, mEvDates,"
                                " mErDates.\n", routine);


        // sort and merge the critical dates 
        if (DrlVTypeVectSort(TPDates,
                             &NbTP,
                             DRL_TDATE_T,
                             TRUE) == FAILURE) 
                throw KFailure("%s: failed sorting critical dates.\n", routine);


        // sort/check/merge mat datelists, and calc latest usage (LaMat) 
        // and earliest usage (ErMat) 
        ProcessDates(mNbMatDates,
                     mZMatDates,
                     mZUseDates,
                     &NbMat,
                     &Mat,
                     &LaMat,
                     &ErMat);
                

        N = NbMat;

        while (N > 0) // size of unprocessed Mat (N) is > 0 
        {
                // First pass: span the maturities 
                n = ZbkG(N, Mat, LaMat, ErMat, today,
                         &NbZ, &Z, &ErZ, &FirstMat);
                if (n < 0) 
                        throw KFailure("%s: number of zero bank maturities"
                                       " is negative (%ld).\n", 
                                        routine, n);

                // now Z has covered Mat[n] .. Mat[N-1] 

                // Second pass: shift the dates 
                ZbkH(NbZ, Z, ErZ, FirstMat,
                         NbCritDates, CritDates,
                         today,
                         mNbEvDate, mEvDates,
                         mNbErDate, mErDates);

                N=n;

                // reset the intermediate zerobank 
                DrlVTypeVectFree((void*) Z, NbZ-1, DRL_LONG_T);
                DrlVTypeVectFree((void*) ErZ, NbZ-1, DRL_LONG_T);
                DrlVTypeVectFree((void*) FirstMatZ, NbZ-1, DRL_LONG_T);
                Z = ErZ = FirstMat = NULL;
                NbZ = 0;

        } // while N > 0 

        &mNbEvDates = &

    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }


}



//--------------------------------------------------------------
//        This routine performs the following to a raw list of zero mat
//      dates and the corresponding usage dates:
//      - checks that mErDates[i] < mEvDates[i]
//      - sort the mEvDates in ascending order
//      - remove duplicate mEvDates
//      - set the earliest usage date as the minimum usage dates
//        of all duplicates
//      - set the latest usage date as the max usage dates of all
//        duplicates
//
void     
KZeroBank::ProcessDates(
        int    mNbEvDates,  // (I) Nb of input mats       
        TDate  *mEvDates,   // (I) zero mat datelist       
        TDate  *mErDates,   // (I) zero usage list          
        int   *NbOutMat,    // (O) Nb of processed zero mats 
        TDate **Mat,        // (O) processed zero mats     
        TDate **LaMat,      // (O) latest usage date list   
        TDate **ErMat)      // (O) earliest usage date list  
{
static        char        routine[] = "KZeroBank::ProcessDates";

        int     i;

        // local datelists for the output 
        int     NbMat = 0;
        int     NbLaMat = 0;
        int     NbErMat = 0;
        long   *MatLocal   = NULL;
        long   *ErMatLocal = NULL;
        long   *LaMatLocal = NULL;

 try{
        // ensure usage date < mat date 
        for (i=0; i<mNbEvDates; i++) {
            if (mErDates[i] >= mEvDates[i])
                throw KFailure("%s: mErDates[%d](%s) >= mEvDates[%d](%s).\n", 
                              routine, i, GtoFormatDate(mErDates[i]), 
                              GtoFormatDate(mEvrDates[i]), i);
        }

        // sort the date lists 
        IF_FAILED_THROW( SortDateList(mNbEvDates, mEvDates, mErDates)); 

        // compact the list to have unique Mat 
        InsertDateToList(&NbMat,
                         &MatLocal,
                         mEvDates[0]);

        InsertDateToList(&NbLaMat,
                         &LaMatLocal,
                         mErDates[0]);

        InsertDateToList(&NbErMat,
                         &ErMatLocal,
                         mErDates[0]);


        for (i=1; i<mNbEvDates; i++)
        {
            if (mEvDates[i] > MatLocal[NbMat-1])
            {
                    // keep this mat date by adding it to output list 
                InsertDateToList(&NbMat,
                                  &MatLocal,
                                  mEvDates[i]);

                InsertDateToList(&NbLaMat,
                                  &LaMatLocal,
                                  mErDates[i]);

                InsertDateToList(&NbErMat,
                                  &ErMatLocal,
                                  mErDates[i]);
           }
           else
           if (mEvDates[i] == MatLocal[NbMat-1])
           {
                    // ignore this mat, only update the latest usage date
                    // the earliest usage is already done due to the sort 
                    LaMatLocal[NbLaMat-1] = mErDates[i];
           }
           else goto RETURN; // should not get here due to the sort 

        } // for i 

        *NbOutMat = NbMat;
        *Mat      = MatLocal;
        *LaMat    = LaMatLocal;
        *ErMat    = ErMatLocal;

    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);

        DrlVTypeVectFree((void*) MatLocal, NbMat-1, DRL_LONG_T);
        DrlVTypeVectFree((void*) LaMatLocal, NbLaMat-1, DRL_LONG_T);
        DrlVTypeVectFree((void*) ErMatLocal, NbErMat-1, DRL_LONG_T);
    }


} // ProcessDates 



//-------------------------------------------------------------------------
//
//      Given the start date and the bounds date (the date the Bnds is
//      calculated), returns a new date satisfying the following:
//      if mode = +1:
//          new date - start date = Bnd(bounds date)
//      if mode = -1:
//          start date - new date = Bnd(bounds date)
//

long    ZbkBnd(TDate  Sdate,     // (I) date the offset is calculated from 
               TDate  BDate,     // (I) date used to determine the bound   
               TDate  today,     // (I) value date                         
               int    mode)      // (I) +1 = Fwd; -1 = Backward            
{
        long   sign;
        double BdateInYrs;
        TDate  newDate;
        TDateInterval        interval;

        sign = ((mode > 0) ? +1L : -1L);
        IF_FAILED_THROW(GtoDayCountFraction(
                                        today,
                                        BDate,
                                        GTO_ACT_365F,
                                        &BDateInYrs));

        if (BdateInYrs <= 2.0)  // 1 month
                IF_FAILED_THROW(GtoMakeDateInterval(sign*1, 'M', &interval));
        else
        if (BdateInYrs <= 5.0)  // 3 months
                IF_FAILED_THROW(GtoMakeDateInterval(sign*3, 'M', &interval));
        else                            // 6 months
                IF_FAILED_THROW(GtoMakeDateInterval(sign*6, 'M', &interval));


        IF_FAILED_THROW(GtoDtFwdAny(Sdate, 
                                    interval,
                                    &newDate));

        return(newDate);

} // ZbkBnd 




//--------------------------------------------------------------------
//
//      Performs the first pass of the algorithm. This routine attempts
//      to span a set of target maturities with the fewest number of
//      zero-bank maturity dates.
//      Returns the zero-bank mats (in descending order) and their
//      corresponding earliest usage dates, together with the index of
//      the smallest processed Mat.
//      Returns -999 on FAILURE
//
//      Only called by OptimizeDates, do NOT call directly
//

int     ZbkG(int          NbMat,     // (I) size of target mats          
             TDate        *Mat,      // (I) target zero mats            
             TDate        *LaMat,    // (I) latest usage for zero mats 
             TDate        *ErMat,    // (I) earliest usage for zero mats 
             TDate        today,     // (I) value date                   
             int         *OutNbZ,    // (O) size of OutZ and ErOutZ      
             TDate       **OutZ,     // (O) zerobank mat array           
             TDate       **OutErZ,   // (O) zerobank earliest usage list 
             TDate       **OutFstMat)// (O) 1st mat in each zero intval  
{
        int    outputN  = -999;     // output: smallest processed Mat[] idx  

        int    k;

        TDate *ZLocal      = NULL;  // local zerobank mat array             
        TDate *ErZLocal    = NULL;  // local zerobank earliest usage array   
        TDate *FstMatLocal = NULL;  // local 1st mat array                   
        int    NbZ      = 0;        // size of ZLocal   
        int    NbErZ    = 0;        // size of ErZLocal 
        int    NbFstMat = 0;        // size of FstMat   

        int    n;                   // curr smallest processed Mat[] index  
        TDate  Zi, ErZi;            // curr zbank mat and earliest use dates 
        TDate  Zdash;               // proposed next zbank mat 

        TDate  minErMat;            // min ErMat[k] with M[k] in (Zdash, Zi) 


        // don't proceed if the output is not empty 
        if (*OutZ   != NULL) goto FREE_MEM_AND_RETURN;
        if (*OutErZ != NULL) goto FREE_MEM_AND_RETURN;

        // 4.4: process the last target mat 
        n    = NbMat-1;
        Zi   = Mat[n];
        ErZi = ErMat[n];

        InsertDateToList(&NbZ,   &ZLocal,   Zi);
        InsertDateToList(&NbErZ, &ErZLocal, ErZi);
        InsertDateToList(&NbFstMat, &FstMatLocal, Mat[n]);

        while (n>0) // there is Mat[0] to Mat[n-1] to process 
        {
        // 4.5: find max backward jump 
        Zdash   = ZbkBnd(Zi,                            // start date 
                         ZbkBnd(Zi, Zi, today, -1),     // bnd date    
                         today,
                         -1);                           // go backwards 

        if (Mat[n-1] < Zdash) break; // 4.6: empty period 

        // 4.7, 4.8: find next n, next Zi, and next ErZi 

        minErMat = 99999999L;
        for (k=n-1; k>=0; k--)
        {
            // if Mat[k] is outside [Zdash, Zi)... 
            if (Mat[k] < Zdash) break;

            Zdash = MAX(Zdash, Nxtday(LaMat[k],1L));
            if (Mat[k] > Zdash) minErMat = MIN(minErMat, ErMat[k]);

        } // for 

        // prepare for next loop 
        n = k+1;
        Zi = Zdash;
        ErZi = MIN(minErMat, ErMat[n]);

        if (Mat[n] == Zi)
            ErZLocal[NbErZ-1] = MIN(minErMat, ErZLocal[NbErZ-1]);
        else
            ErZLocal[NbErZ-1] = MIN(ErZi,     ErZLocal[NbErZ-1]);

        if (AddDateToList(&NbZ,   &ZLocal,   Zi)   == FAILURE)
                goto FREE_MEM_AND_RETURN;
        if (AddDateToList(&NbErZ, &ErZLocal, ErZi) == FAILURE)
                goto FREE_MEM_AND_RETURN;
        if (AddDateToList(&NbFstMat, &FstMatLocal, Mat[n]) == FAILURE)
                goto FREE_MEM_AND_RETURN;

    } // while 

    outputN = n;
    *OutNbZ = NbZ;
    *OutZ = ZLocal;
    *OutErZ = ErZLocal;
    *OutFstMat = FstMatLocal;

FREE_MEM_AND_RETURN:

    if (outputN == -999)
    {
        Free_DR_Array(ZLocal,      LONG, 0, NbZ);
        Free_DR_Array(ErZLocal,    LONG, 0, NbErZ);
        Free_DR_Array(FstMatLocal, LONG, 0, NbFstMat);
    }
    return (outputN);

} // ZbkG 




//****  ZbkH  *****************************************************
//
//      Performs the second pass of the algorithm. This routine attempts
//      to shift the zero-bank mats from the first pass to the right, so
//      that they may coincide with the critical dates
//      Returns the final zero-bank mats (in ascending order) and their
//      corresponding earliest usage dates.
//
//      Only called by ZbkOptDates, do NOT call directly
//

int     ZbkH(int       NbZ,         // (I) size of 1st pass zeros          
             long     *Z,           // (I) 1st pass zero mats (desc order)
             long     *ErZ,         // (I) earliest usage for zero mats  
             long     *FirstMat,    // (I) 1st mat date in [Zi, Zi+1)   
             int       NbC,         // (I) size of critical dates      
             long     *C,           // (I) critical dates             
             long      ValueDate,   // (I) value date                
             int      *NbMatDates,  // (O) Nb of output zerobank mats    
             long    **ZbkMats,     // (O) output zerobank mats         
             int      *NbErDates,   // (O) Nb of output earliest usage 
             long    **ZbkErs)      // (O) zerobank earliest usage list
{

    int  status = FAILURE;

    int  i;
    long UpBound;          // max right shift allowed               
    int  UseMaxC;          // TRUE => shift to a critical date      
    int  MaxT;             // index of the max critdate to shift to


    // 4.11: process the first zerobank mat 
    if (AddDateToList(NbMatDates, ZbkMats, FirstMat[NbZ-1]) == FAILURE)
            goto RETURN;
    if (AddDateToList(NbErDates, ZbkErs, ErZ[NbZ-1]) == FAILURE)
            goto RETURN;

    for (i=NbZ-2; i>=1; i--)
    {
        // 4.12: find upper bound for right shift 

        UpBound = ZbkBnd((*ZbkMats)[*NbMatDates-1],   // start date 
                         (*ZbkMats)[*NbMatDates-1],   // bnd date  
                         ValueDate,
                         +1);                         // go forward 

        UpBound = MIN(UpBound, Z[i-1]);
        UpBound = MIN(UpBound, FirstMat[i]);

        // find max critical date index <= UpBound 
        MaxT = GetDLOffset(NbC, C, UpBound, CbkLOWER);
        UseMaxC = ((MaxT<0) ? FALSE : (C[MaxT] > Z[i]));

        // 4.13: do the shift 
        if (UseMaxC)
        {
            // shift zerobank mat to the critical date 
            if (AddDateToList(NbMatDates, ZbkMats, C[MaxT]) == FAILURE)
                    goto RETURN;
            if (AddDateToList(NbErDates, ZbkErs, ErZ[i]) == FAILURE)
                    goto RETURN;
        }
        else
        {
            // shift to the max allowed date to allow further right shift 
            if (AddDateToList(NbMatDates, ZbkMats, UpBound) == FAILURE)
                    goto RETURN;
            if (AddDateToList(NbErDates, ZbkErs, ErZ[i]) == FAILURE)
                    goto RETURN;
        }
    } // for 

    // 4.15: process the last zerobank mat 
    if (NbZ > 1)
    {
        if (AddDateToList(NbMatDates, ZbkMats, Z[0]) == FAILURE)
                goto RETURN;
        if (AddDateToList(NbErDates, ZbkErs, ErZ[0]) == FAILURE)
                goto RETURN;
    }

    status = SUCCESS;

RETURN:

    return (status);

} // ZbkH 



//****  CbkProcessDL  ********************************************
//
//      This routine performs the following to a raw list of payoff known
//      dates (EvalDates) and the corresponding earliestUsage dates:
//      - sort the dates by EvalDates ascending order
//      - remove duplicate EvalDates and replace the earliest usage date by
//        the minimum earliest usage dates of all duplicates
//      - return the two lists and their size through the same arguments
//      Returns SUCCESS or FAILURE
//

int     CbkProcessDL
            (int       *NbEvDates,
             long     **EvDL,
             int       *NbErDates,
             long     **ErDL)
{
    int     status = FAILURE;   // Error status = FAILURE initially 
    char    ErrorMsg[MAXBUFF];
    int     i;
    int     freeOfs;            // lowest ofs with obsolete data 
    int     newNbDates = 0;     // new number of dates in list   

    long   *EvDL_Local = NULL;
    long   *ErDL_Local = NULL;

    if ((NbEvDates == NULL) ||
        (NbErDates == NULL) ||
        (EvDL      == NULL) ||
        (ErDL      == NULL))
        goto RETURN;

    if (*NbEvDates != *NbErDates) goto RETURN;
    if (*NbEvDates == 0) return(SUCCESS); // nothing to do
    if (((*EvDL) == NULL) || ((*ErDL) == NULL)) goto RETURN;

    // ensure ErDates <= EvDates 
    for (i=0; i<*NbEvDates; i++)
    {
        if ((*ErDL)[i] > (*EvDL)[i])
        {
            sprintf (ErrorMsg,
                     "CbkProcessDL: Earliest usage date #%d (%ld) > "
                     "evaluation date (%ld)!\n",
                     i+1, (*ErDL)[i], (*EvDL)[i]);
            DR_Error (ErrorMsg);
            goto RETURN;
        }
    }

    // sort the date lists 
    if (SortDateList(*NbEvDates, *EvDL, *ErDL) == FAILURE) goto RETURN;

    // compact the list to have unique EvalDates and 
    // smallest EarliestUseDate                     
    freeOfs = 1;
    newNbDates = *NbEvDates;
   
    for (i=1; i<*NbEvDates; i++)
    {
        if ((*EvDL)[i] > (*EvDL)[freeOfs-1])
        {
            // keep this EvDate by moving it to the free slot 
            (*EvDL)[freeOfs] = (*EvDL)[i];
            (*ErDL)[freeOfs] = (*ErDL)[i];
            freeOfs++;
        }
        else
        if ((*EvDL)[i] == (*EvDL)[freeOfs-1])
        {
            // remove this date, the sort already ensures that
            // the smaller earliestUse date is at freeOfs-1    
            newNbDates--;
        }
        else goto RETURN; // should not get here due to the sort 

    } // for i 

    // realloc the date lists if necessary 
    if ((*NbEvDates) != newNbDates)
    {
        EvDL_Local = (long *)realloc(*EvDL,newNbDates*sizeof(long));
        ErDL_Local = (long *)realloc(*ErDL,newNbDates*sizeof(long));

        if ((EvDL_Local == NULL) || (ErDL_Local == NULL))
        {
            DR_Error ("CbkProcessDL: re-allocation failure!");
            goto RETURN;
        }

        *EvDL = EvDL_Local;
        *ErDL = ErDL_Local;
    }

    *NbEvDates = newNbDates;
    *NbErDates = newNbDates;

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("CbkProcessDL: Failed!\n");
        Free_DR_Array(EvDL_Local, LONG, 0, newNbDates-1);
        Free_DR_Array(ErDL_Local, LONG, 0, newNbDates-1);
    }

    return (status);

} // CbkProcessDL 


***************************************************************/

