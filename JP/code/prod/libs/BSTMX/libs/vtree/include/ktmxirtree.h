#ifndef	_TMXIRTREE_H
#define	_TMXIRTREE_H

#include <fstream>
#include "kutilios.h"
#include "kpirtree.h"

#ifdef PI
#undef PI
#endif

#ifdef TINY
#undef TINY
#endif

#ifdef ERROR
#undef ERROR
#endif

#ifndef SHIFT_ZERO
#define SHIFT_ZERO(x) ((x) = (IS_ALMOST_ZERO(x) ? 1e-5 : (x)))
#endif


namespace TMX
{
extern "C"{
#include "q3.h"
#include "tmx123head.h"
#include "l_date.h"
#include "drtypes.h"
}
}

//--------------------------------------------------------------
/**
 * A basis interest rate tree class based
 * on an TMX interest rate tree.
 */

class KTMXirTree : public KPirTree {
public:

    KTMXirTree();

    virtual ~KTMXirTree();
  
    virtual void    DeleteMemory();
    
    virtual void    Calibrate();
    
    virtual void    Update(int tpIdx);
    
    virtual KTSlice& TSliceDev(KTSlice&, const String&);

    // No CalcDiscount in TMX
    virtual void    CalcDiscount(int tpIdx,
                                 int nIRDim,
                                 double Zt,
                                 KMap(int, double*) &mDiscountIR)
    {
        if ( TMXFlag == TRUE)
        {
            throw KFailure(" Failure: TMX does not need CalcDiscount()\n "); 
        }
        else
        {
            KPirTree::CalcDiscount(tpIdx, nIRDim, Zt, mDiscountIR);
        }
        return;

    }

    // No offset computation in TMX
    virtual void    SolveOffset(
            int tpIdx,          // (I) time point index
            int nIRDim,         // (I) IR dimemsion
            double *discount,   // (I) slice disc bet t and t+1
            double *statePr,    // (I) slice state price at t
            double zeroPrice,   // (I) zero price at t+1
            double Zt,          // (I) intial offset
            double *DelZt)      // (O) incremental offset 
    {   
        if ( TMXFlag == TRUE)
        {
            throw KFailure(" Failure: TMX does not need SolveOffset()\n "); 
        }
        else
        {
            KPirTree::SolveOffset(tpIdx, nIRDim, discount, statePr, zeroPrice, Zt, DelZt);
        }
        return;
    }

    
    
    /** Get a zero reset at the current time point
     */
    virtual void    Get(KTSlice &ts, const KZeroReset &zeroReset);
   
    /** Check tree validity */
    virtual void    CheckTreeValid();

    void    Initialize(
        TDate              todayDt,     // (I) today's date
        KMrParam           &mrPar,      // (I) full dimension mr info
        KSmileParam        &irSmilePar, // (I) ir smile info
        int                 nIRDim,     // (I) IR dimension
        KVector(int)       &cvTypes,    // (I) cv types (KV_DIFF...)
        KMap(int, KZCurve) &cv,         // (I) array of curves
        KMap(int, String)  &cvNames,    // (I) array of curve names
        KMap(int, TDate)   &cvValueDates, // (I) array of cv value dates
        KVector(TDate)     &volDates,   // (I) volatility dates
        KVector(KVector(double)) &factVol, // (I) spot vols
        KResetBank         &resetBank); // (I) rate reset bank

    void TMXirTreeCet(KMrParam&);

    void CetTreeFree();

    void CalibrateNmrInv();

    virtual void    SetUpTimeline(); // Use the new algorithm to setup the timeline.
 
	void	InitializeTimeline_TMX(
		char EoI,		// (I) Equal or increasing time steps
		int  ppy);		// (I) period per year
    
protected:

    virtual void    eDevSetUp();

    //------------------------------------------------------
    // Protected methods and data
    //------------------------------------------------------
    //

    virtual void    CalibrateDrift();

    virtual void    TreeSetUp();


    // virtual void     DeleteMemory();

    // -----------------------------------------------------
    // TMX basis tree specific public methods
    // -----------------------------------------------------
public:
    int irDim()
    { return (mIRDim);}

// Critical and Nmr dates functions

    int NBCritDates()
    {   return (CritDates.size());}

    int NBNmrDates()
    {   return (NmrDates.size()); }


    TDate CritDate(int idx)
    {   return (CritDates.at(idx));}

    TDate NmrDate(int idx)
    {   return (NmrDates.at(idx));}


    void    SetNmrToCcy(int);

    void    SetCcyToNmr();

    int     GetNmrToCcy()
    {   return (NmrToCcy); };

    int     GetCcyToNmr()
    {   return (CcyToNmr); };

    void    ClearNmrCcyFlag();

    bool    IsNmrDate(int tpIdx);

    bool    IsCritDate(int tpIdx);

    //--------------------------------------------------------
    // Interface function to set up TMX wrapper objects
    //--------------------------------------------------------
    void    SetTmxFlag(int flag)
    {   TMXFlag = flag;  };

    int     GetTmxFlag()
    {   return (TMXFlag); };


    long    getToday();
    long    getValueDate();

    void    WrapperEnvPack( KMarketCurves&,
                             KVolDiag&,
                             KMrParam&,
                             KSmileParam&,
                             const String&);
    

    void    WrapTreeTimeLine(KVolDiag&,
                              KMrParam&);
    
    void    PrintTreeTimeLine();
    

//---------------------------------------------------------------------



protected:
    KVector(TDate)          NmrDates;          // Nmr dates
    KVector(TDate)          CritDates;         // CritDates to set NmrToCcy and CcyToNmr
    KMap(int, double*)      NmrInv;            // Nmr Inv at current TP
    KMap(int, double*)      NmrInvLag;         // Nmr Inv at previous TP

    KMap(TDate, double*)    mDiffNmrInv;       // Diffusion idx Nmr Inv at each TP

    KMap(int, KTSlice*)     mCurrIRNmrInv;     // Current NmrInv Slices for all 
                                               // indicies

    int     NmrToCcy;
    int     CcyToNmr;

    int     TMXFlag;        // TRUE/FALSE


    TMX::T_CURVE              t_curve[3];       /* Structure of zero curve data   */
    TMX::MKTVOL_DATA*         mktvol_data;      /* Structure of IR vol data       */
    TMX::TREE_DATA*           tree_data;        /* Structure of tree data         */
    TMX::DEV_DATA*            dev_data;


    char    diffCvName[MAXBUFF];

    long    treeToday;      // TMX tree today, in DRDate format
    long    treeValueDate;  // TMX tree valuedate, In DRDate format
    int     NbDailyPts;     // Number of daily points


    void NmrUpdate(int tpIdx);

    virtual double*     GetNmrInv(int curveIdx);

    virtual double*     GetNmrInvLag(int curveIdx);

    // Functions for wrapping the TMX engine
    void setToday(long);
    
    void setValueDate(long);

    // Insert Numeraire dates
    void InsertNmrDates_old(KVolDiag&); 
    void InsertNmrDates(KVolDiag&); 

    // Set up critical dates
    void SetUpCritDate();

    void PrintTPDate();

    void ConvertBpVolToPercVol();

    bool IsCalibIdxBase(const KVolDiag &);

    void PackVolData( const KVolDiag & ,
                      const KMarketCurves &,
                      const KSmileParam &,
                      const KMrParam&);      /* Pack IR vol                 */

    void PackCurves ( const KMarketCurves&,
                      const String& );  /* Pack curves                 */

    void PackSingleDRCurve (const KZCurve&, int);
    
    void PrintEnv();

    void CpTreeParam();    


};

#endif

