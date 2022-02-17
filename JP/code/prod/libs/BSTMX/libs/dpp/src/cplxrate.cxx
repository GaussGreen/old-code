/****************************************************************
 * Module:  Complex rate
 * Submodule:	
 * File:    cplxrate.cxx
 * Function:
 * Author:  Changhong He
 * Revision:$Header$
 *****************************************************************/
#include "kstdinc.h"    /* Standard definitions & error hadling */
#include "kutilios.h"

#include "kcplxrate.h"



KCplxRate::KCplxRate(
                const String&          name,
                const KVector(KRate*)& rateInstr,    // (I) rate indexes
                const String&          formula)      // (I) coupon formula
{
    SetName(name);
    mformula = formula;

    for( KVector(KRate*)::const_iterator iterRate = rateInstr.begin();
         iterRate != rateInstr.end();
         ++iterRate)      
            mRates.push_back(new KRate (**iterRate));

    mCurveName = mRates[0]->CurveName();
    cout << "standard constructor " << mCurveName << endl;
}

//**********************************************************
// Construct object from KRate
KCplxRate::KCplxRate( const KRate& Rate)
:KVPAtom(Rate)
{
    mRates.push_back(new KRate(Rate));
    mformula = "x0";
    mCurveName = Rate.CurveName();
}

KCplxRate&
KCplxRate:: operator = (const KRate& Rate)
{
    for(KVector(KRate*)::iterator it= mRates.begin();
        it!=mRates.end(); ++it)
            delete (*it);
    mRates.clear();

    mformula = "x0";
    mCurveName = Rate.CurveName();
    
    mRates.push_back( new KRate (Rate));
    
    return *this;
}

KCplxRate::KCplxRate(
        const KDateInterval &mat,       // (I) rate maturity
        const KDateInterval &freq,      // (I) rate payment freq
        KDayCc dayCc,                   // (I) rate day count conv
        const KDateInterval &spotOffset,// (I) spot offset
        double spread,                  // (I) spread
        double weight)                  // (I) weight
{
    mRates.push_back( new KRate(
                            mat,
                            freq,
                            dayCc,
                            spotOffset,
                            spread,
                            weight));
    mformula = "x0";
    mCurveName = mRates[0]->CurveName();
}

KCplxRate::KCplxRate(
        const String  &curveName,       // (I) curve name
        const KDateInterval &mat,       // (I) rate maturity
        const KDateInterval &freq,      // (I) rate payment freq
        KDayCc dayCc,                   // (I) rate day count conv
        const KDateInterval &spotOffset,// (I) spot offset
        double spread,                  // (I) spread
        double weight)                  // (I) weight
{
    mRates.push_back(new KRate(
                            curveName,
                            mat,
                            freq,
                            dayCc,
                            spotOffset,
                            spread,
                            weight));
    mformula = "x0";
    mCurveName = curveName;
}

KCplxRate::KCplxRate(
        const String  &curveName,       // (I) curve name
        const TDate   startDate,        // (I) start date
        const TDate   endDate,          // (I) end date
        KDayCc dayCc,                   // (I) rate day count conv
        const KDateInterval &spotOffset,// (I) spot offset
        double spread,                  // (I) spread
        double weight)                  // (I) weight
{
    mRates.push_back(new KRate(
                            curveName,
                            startDate,
                            endDate,
                            dayCc,
                            spotOffset,
                            spread,
                            weight));
    mformula = "x0";
    mCurveName = curveName;
}


KCplxRate::KCplxRate(
        const String  &name,            // (I) debug name
        const String  &curveName,       // (I) curve name
        const KDateInterval &mat,       // (I) rate maturity
        const KDateInterval &freq,      // (I) rate payment freq
        KDayCc dayCc,                   // (I) rate day count conv
        const KDateInterval &spotOffset,// (I) spot offset
        double spread,                  // (I) spread
        double weight)                  // (I) weight
{
    mRates.push_back(new KRate(
                            name,
                            curveName,
                            mat,
                            freq,
                            dayCc,
                            spotOffset,
                            spread,
                            weight));
    mformula = "x0";
    mCurveName = curveName;
}

//**********************************************************

KCplxRate::KCplxRate( const KCplxRate& cplxRate)
:KVPAtom(cplxRate)
{
    mformula = cplxRate.mformula;
    mCurveName = cplxRate.CurveName();
    
    for (int i = 0; i < cplxRate.nbInstr(); i++)
    {
        //mRates.push_back(&cplxRate.instrIdx(i));
        mRates.push_back( new KRate (cplxRate.instrIdx(i)));
    }
    cout << "copy constructor" << mCurveName << endl;
    
}


KCplxRate&
KCplxRate:: operator = (const KCplxRate& cplxRate)
{
    if (this != &cplxRate)
    {
    	for(KVector(KRate*)::iterator it= mRates.begin();
            it!=mRates.end(); ++it)
                delete (*it);
        mRates.clear();

        mformula = cplxRate.mformula;
        mCurveName = cplxRate.CurveName();
    
        for (int i = 0; i < cplxRate.nbInstr(); i++)
        {
    //        mRates.push_back(&cplxRate.instrIdx(i));
            mRates.push_back( new KRate (cplxRate.instrIdx(i)));
        }
    }
    return *this;
}



KCplxRate::KCplxRate()
{
    mformula = " ";
    mCurveName = K_DEFAULT_NAME;
}

//---------------------------------------------------------------
//
KCplxRate::~KCplxRate()
{
	for(KVector(KRate*)::iterator it= mRates.begin();
        it!=mRates.end(); ++it)
                delete (*it);
 
        mRates.clear();
}


/*KRate*
KCplxRate::instrIdx(int idx) const
{
    return mRates.at(idx);
}*/


KRate&
KCplxRate::instrIdx(int idx) const
{
    return *(mRates[idx]);
}


const String&
KCplxRate::getFormula() const
{
    return mformula;
}

istream&
KCplxRate::Get(istream& is, int drwFmt)
{

	return(is);
}





//---------------------------------------------------------------

ostream&
KCplxRate::Put(ostream& os, int indent) const
{
    os << "COMPLEX RATE: Instruments\n" ;
    for (int i = 0; i < mRates.size(); i++)
    {
        os << instrIdx(i).Put(os, indent) << endl;
    }
    os << "Formula: " << mformula;
	return(os);
}



//---------------------------------------------------------------

ostream&
KCplxRate::WriteSimple(ostream& os)
{
   /* for( KVector(KRate*)::iterator iterRateIdx = mRates.begin();
    iterRateIdx != mRates.end();
    ++iterRateIdx)
    {
        (*iterRateIdx)->WriteSimple(os);
    }*/

    for (int i = 0; i < mRates.size(); i++)
    {
        instrIdx(i).WriteSimple(os);
    }
	WriteDone();

	return(os);
}




//---------------------------------------------------------------

ostream&
KCplxRate::YacctionWrite(ostream& os, int indent)
{

	if (GetWriteFlag())
	{
	    os  << GetName() << "=";

	    WriteSimple(os); 

	    os << ";" << endl << endl;

	    WriteDone();
	}

	return(os);
}

//---------------------------------------------------------------
// return TRUE/FALSE if all indexes are floating
int
KCplxRate::IsAllFloating() const
{
    int flag;
    int i;
    flag = instrIdx(0).IsFloating();
    
    for (i = 1; i < nbInstr(); i++)
    {
        flag = flag && instrIdx(i).IsFloating();
    }
    return flag;

}
//***************************************************************
//          CplxRateReset functions
//***************************************************************

//---------------------------------------------------------------
KCplxRateReset::KCplxRateReset( 
    const KVector(KRateReset*)& rateResetIndexes,
    const KCplxRate&            cplxRate)
{
    mcplxRate = cplxRate;
    for ( KVector(KRateReset*)::const_iterator iterRateReset 
                    = rateResetIndexes.begin();
        iterRateReset != rateResetIndexes.end();
        ++iterRateReset)
    {
        mRateResets.push_back(
            new KRateReset( **iterRateReset ));
    }

    mResetDate = rateResetIndexes[0]->ResetDate();
    mEffDate   = rateResetIndexes[0]->EffDate();
}

//---------------------------------------------------------------
// All indexes reset at same dates

KCplxRateReset::KCplxRateReset(
    TDate                   resetDate,
    const KCplxRate&        cplxRate)
{

    for ( int i = 0; i < cplxRate.nbInstr(); i++)
    {
        mEffDate   = resetDate + cplxRate.instrIdx(i).SpotOffset();
        mRateResets.push_back(
            new KRateReset( resetDate,
                            mEffDate,
                            cplxRate.instrIdx(i)));   
    }

    mResetDate = resetDate;
    mEffDate   = resetDate + cplxRate.instrIdx(0).SpotOffset();
    mcplxRate  = cplxRate;

}

//---------------------------------------------------------------
// All indexes reset at same dates

KCplxRateReset::KCplxRateReset(
    TDate   resetDate,
    TDate   effDate,
    const KCplxRate &       cplxRate)
{
    for ( int i = 0; i < cplxRate.nbInstr(); i++)
    {
        mRateResets.push_back(
            new KRateReset( resetDate,
                            effDate,
                            cplxRate.instrIdx(i)));   
    }
    mResetDate = resetDate;
    mEffDate   = effDate;
    mcplxRate  = cplxRate;

}

//---------------------------------------------------------------

KCplxRateReset::KCplxRateReset (
        TDate   resetDate,     // (I) observation dates
        TDate   effDate,       // (I) effective dates
        TDate   endDate,       // (I) end dates
        const KCplxRate &cplxRate)        // (I) rate description
{
    for ( int i = 0; i < cplxRate.nbInstr(); i++)
    {
        mRateResets.push_back(
            new KRateReset( resetDate,
                            effDate,
                            endDate,
                            cplxRate.instrIdx(i)));   
    }
    mResetDate = resetDate;
    mEffDate   = effDate;
    mcplxRate  = cplxRate;

}

//---------------------------------------------------------------

KCplxRateReset::KCplxRateReset( KRateReset& rateReset)
{
    mRateResets.push_back(
            new KRateReset( rateReset ));   

    mcplxRate  = KCplxRate( rateReset.Rate());
    mResetDate = rateReset.ResetDate();
    mEffDate   = rateReset.EffDate();

}

//---------------------------------------------------------------

KCplxRateReset::KCplxRateReset( 
                 TDate resetDate,
                 const KRate& rate)
{
    mRateResets.push_back(
            new KRateReset( resetDate,
                            rate));   
    mcplxRate  = KCplxRate( RateResetIndex(0).Rate());
    mResetDate = RateResetIndex(0).ResetDate();
    mEffDate   = RateResetIndex(0).EffDate();    

}

KCplxRateReset::KCplxRateReset(
        TDate resetDate,        // (I) observation date
        TDate effDate,          // (I) effective date
        const KRate& rate)      // (I) rate description
{
    mRateResets.push_back(
            new KRateReset( resetDate, 
                            effDate,
                            rate));
    mcplxRate  = KCplxRate( RateResetIndex(0).Rate());
    mResetDate = RateResetIndex(0).ResetDate();
    mEffDate   = RateResetIndex(0).EffDate();    


}

KCplxRateReset::KCplxRateReset(
        const String& curveName,    // (I) curve name
        TDate resetDate,            // (I) observation date
        TDate effDate,              // (I) effective date
        const KRate& rate)          // (I) rate description
{
    mRateResets.push_back(
            new KRateReset( curveName,
                            resetDate, 
                            effDate,
                            rate));
    mcplxRate  = KCplxRate( RateResetIndex(0).Rate());
    mResetDate = RateResetIndex(0).ResetDate();
    mEffDate   = RateResetIndex(0).EffDate();    


}

KCplxRateReset::KCplxRateReset(
        TDate resetDate,        // (I) observation date
        TDate effDate,          // (I) effective date
        TDate endDate,          // (I) End date
        const KRate& rate)      // (I) rate description
{
    mRateResets.push_back(
            new KRateReset( resetDate, 
                            effDate,
                            endDate,
                            rate));
    mcplxRate  = KCplxRate( RateResetIndex(0).Rate());
    mResetDate = RateResetIndex(0).ResetDate();
    mEffDate   = RateResetIndex(0).EffDate();    
}


    
//---------------------------------------------------------------

KCplxRateReset::KCplxRateReset( const KCplxRateReset& cplxRateReset)
{
    for (int i = 0; i < cplxRateReset.mRateResets.size(); i++)
    {
        mRateResets.push_back(
            new KRateReset(*(cplxRateReset.mRateResets[i])));
    }
    mResetDate = cplxRateReset.mResetDate;
    mEffDate   = cplxRateReset.mEffDate;
    mcplxRate  = KCplxRate(cplxRateReset.CplxRate());
}

//---------------------------------------------------------------
bool
KCplxRateReset::operator < (const KCplxRateReset& cplxRateReset) const
{
    return (RateIndex(0) < cplxRateReset.RateIndex(0));
}

//---------------------------------------------------------------
// Destructor
KCplxRateReset::~KCplxRateReset()
{
	for(KVector(KRateReset*)::iterator it= mRateResets.begin();
        it!=mRateResets.end(); ++it)
                delete (*it);
 
        mRateResets.clear();
}

ostream&
KCplxRateReset::Put(ostream& os, int indent) const
{
    for (int i = 0; i < mRateResets.size(); i++)
    {
        os << "Idx " << setw(2) << i << ": ";
        mRateResets[i]->Put(os, indent);
        os << endl;
    }
    os << "Formula: " << this->CplxRate().getFormula();
    return (os);
}

