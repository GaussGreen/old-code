/****************************************************************
 * Module:  Indicator
 * Submodule:	
 * File:    indicator.cxx
 * Function:
 * Author:  Changhong He
 * Revision:$Header$
 *****************************************************************/
#include "kstdinc.h"    /* Standard definitions & error hadling */
#include "kutilios.h"

#include "kindicator.h"

KIndicator::KIndicator( 
                const String&          name,         // (I) indicator name
                const KVector(KRate*)& rateInstr,    // (I) instruments
                const String&          formula,      // (I) formula
                double                 lbarrier,     // (I) lower barrier
                double                 hbarrier,     // (I) higher barrier
                const KKnockIO        &ioWindow,     // (I) Inside or Outside
                const KSmooth         &smooth)       // (I) Node smoothing flag
{
        mCplxRate =  new KCplxRate (name,
                                    rateInstr,
                                    formula);

        mLbarrier = lbarrier;
        mHbarrier = hbarrier;
        mIoWindow = ioWindow;
        mSmooth   = smooth;
        mName     = name;
}


KIndicator::KIndicator( 
        const String&          name,         // (I) indicator name
        const KCplxRate       *rateIndex,    // (I) complex rate index
        double                 lbarrier,     // (I) lower barrier
        double                 hbarrier,     // (I) higher barrier
        const KKnockIO        &ioWindow,     // (I) Inside or Outside
        const KSmooth         &smooth)       // (I) Node smoothing flag
{
        mCplxRate = new KCplxRate (*rateIndex);
        mLbarrier = lbarrier;
        mHbarrier = hbarrier;
        mIoWindow = ioWindow;
        mSmooth   = smooth;
        mName     = name;
}

KIndicator::KIndicator( 
        const String&          name,         // (I) indicator name
        const KRate           *rateIndex,    // (I) complex rate index
        double                 lbarrier,     // (I) lower barrier
        double                 hbarrier,     // (I) higher barrier
        const KKnockIO        &ioWindow,     // (I) Inside or Outside
        const KSmooth         &smooth)       // (I) Node smoothing flag
{
        mCplxRate = new KCplxRate (*rateIndex);
        mLbarrier = lbarrier;
        mHbarrier = hbarrier;
        mIoWindow = ioWindow;
        mSmooth   = smooth;
        mName     = name;
}


KIndicator::KIndicator(const KIndicator& indicator)
{
        mCplxRate = new KCplxRate(indicator.CplxRate());
        mLbarrier = indicator.mLbarrier;
        mHbarrier = indicator.mHbarrier;
        mIoWindow = indicator.mIoWindow;
        mSmooth   = indicator.mSmooth;
        mName     = indicator.mName;
}

KIndicator&
KIndicator::operator =(const KIndicator& indicator)
{
        mCplxRate = new KCplxRate(indicator.CplxRate());
        mLbarrier = indicator.mLbarrier;
        mHbarrier = indicator.mHbarrier;
        mIoWindow = indicator.mIoWindow;
        mSmooth   = indicator.mSmooth;
        mName     = indicator.mName;

        return *this;
}


KIndicator::KIndicator()
{
    mCplxRate = NULL;
}

KIndicator::~KIndicator()
{
    delete (mCplxRate);
}

//---------------------------------------------------------------

istream&
KIndicator::Get(istream& is, int drwFmt)
{

	return(is);
}

//---------------------------------------------------------------

ostream&
KIndicator::Put(ostream& os, int indent) const
{
    os << "Indicator: Instrument" << endl;
    mCplxRate->Put(os, indent);
    os << endl;
    os << "Lower Barrier = " << mLbarrier << '\t'
       << "Higher Barrier = " << mHbarrier << '\t'
       << "Inside / Outside = " << mIoWindow << '\t'
       << "Smooth = " << mSmooth << endl;
	return(os);
}



//---------------------------------------------------------------

ostream&
KIndicator::WriteSimple(ostream& os)
{

    mCplxRate->WriteSimple(os);
    os << "Lower Barrier = " << mLbarrier << '\t'
       << "Higher Barrier = " << mHbarrier << '\t'
       << "Inside / Outside = " << mIoWindow << '\t'
       << "Smooth = " << mSmooth << endl;        
    
    WriteDone();

	return(os);
}




//---------------------------------------------------------------

ostream&
KIndicator::YacctionWrite(ostream& os, int indent)
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

