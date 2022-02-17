/****************************************************************
 * Module:  VirtualProduct
 * Submodule:
 * File:    vtlkionew.h
 * Function:
 * Author:  Changhong He
 *****************************************************************/
#ifndef _vtlkionew_H
#define _vtlkionew_H

#include "vtlbase.h"
#include "vpkionew.h"


//--------------------------------------------------------------
/**
 * Class for a weighted vector of assets.
 */

class KVPToolKnockIONew : public KVPToolAtom {
public:
    /**
     * Constructor.
     */
    KVPToolKnockIONew(SharedPointer<KVPKnockIONew> ins, KVTree &vt);

    /**
     * Destructor.
     */
virtual ~KVPToolKnockIONew();


    /**
     * Type name.
     */
virtual const char* TypeName() const {return("KVPToolKnockIONew");};

    /**
     * Returns the name for the geometry of the slice
     * (see description of method in KVPToolAtom).
     */
virtual	const String&	GetCurveName();


    /**
     * Returns the KVPAtom to which the KVPToolAtom is associated.
     */
virtual SharedPointer<KVPAtom>	Atom()
    {
        SharedPointer<KVPAtom> vp;
        SharedPointerConvertTo(mKnockIO, vp);
        return (vp);
    }

    /**
     * Update tool.
     */
virtual void    Update();


    /**
     * Returns the current value of the stream.
     */
virtual KTSlice&    GetValue();

    /**
     * Calculates and returns the results at time.
     * In addition to "PV", also contains the fields
     * "XTPROB" (probability of hit), "XTEXP" (expected hit time)
     * "XTSDV" (standard dev of hit time) and "XTFUG" (fugit).
     */
virtual KMap(String,double)  GetResults()
    {
        return mResults;
    }

void setCurveName();

private:
                /** Pointer to inetrument. */
    SharedPointer<KVPKnockIONew>    mKnockIO;
                /** Curve name for slice geometry. */
    String      mCurveName;
                /** Curve name for trigger rate.  */
    String      mRateIdxName;

                /** Default value
                  * This is the default value assigned to
                  * the component, to avoid return error in GetValue.
                  */
    KTSlice     *mDefValue;

                /** Current Value Slice */
    KTSlice     *mValue;
                /** Unsmoothed current Value Slice */
    KTSlice     *mValueUS;
                /** Current rebate slice */
    KTSlice     *mRebateValue;
                /** Current unsmoothed rebate slice */
    KTSlice     *mRebateValueUS;
                /** Used to compute the current KIO value. */
    KTSlice     *mGetValue;


                /** Model observation reset dates */
    KVector(TDate)  mModelObsDates;
                /** Claim bank for underlying settlements */
    KMap(TDate, KTSlice*)   mExer;
                /** Claim bank for rebate at settlement dates */
    KMap(TDate, KTSlice*)   mRebate;

                /** Run the fugit calculations, etc */
    int         mRunStat;
                /** Indicates which dependency is the ko index*/
    int         mDepKoind;
                /** Probability of hit */
    KTSlice     *mXT0;
                /** Expected hit time */
    KTSlice     *mXT1;
                /** Expected hit time^2 */
    KTSlice     *mXT2;

                /** Results (PV,UNDER,FUGIT,ETC) */
    KMap(String,double) mResults;

    template<class T1, class T2> // T1, T2 can only be double or KTSlice
    void    BinaryExpectation(const KTSlice* probUp, 
                              const T1*      valueUp,
                              const T2*      valueDown,
                              KTSlice*       meanValue);
};




#endif




