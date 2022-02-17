/************************************************************************
 * Module:        PenGuin
 * File:
 * Function:    
 * Author:      David Liu
 ************************************************************************/
#include "vpprotleg.h"

#include "kutilios.h"                // ios utilities

extern "C" {
#include "drltime.h"
};

//---------------------------------------------------------------


KVPProtLeg::KVPProtLeg(
        const char *name,           // (I) object name
        const TDate     stDate,     // (I) start date
        const TDate     endDate,    // (I) end date
        const double    notional,   // (I) notional
        const char     *recovery,   // (I) recovery rate
        KProtPayConv    payType,    // (I) default payment type
        const char     *discZcName) // (I) CDS curve name
        : KVPInstr(name), mNtlStartDts(1, stDate), 
          mNtlEndDts(1, endDate), mNtlAmts(1, notional), mPayType(payType)  
{
static        char        routine[] = "KVPProtLeg::KVPProtLeg";

        //
        // Check recovery passed
        //
        if ((recovery == NULL)         ||
            (!strlen(recovery) != 0)   ||
            (!strcmp(recovery, "NIL")) ||
            (!strcmp(recovery, "nil")) ||
            (!strcmp(recovery, "Nil")) )
        {
                mIsDefRecovery = true;
        }
        else
        {
            mIsDefRecovery = false;
            sscanf(recovery, "%lf", &mRecovery);
        }

        mDiscZcName = String(discZcName);

}

KVPProtLeg::KVPProtLeg(
        const char *name,           // (I) object name
        const KVector(TDate)&     ntlStartDates,     // (I) dates
        const KVector(TDate)&     ntlEndDates,
        const KVector(double)&    ntlAmts,   // (I) protected notionals
        const char     *recovery,   // (I) recovery rate
        KProtPayConv    payType,    // (I) default payment type
        const char     *discZcName) // (I) CDS curve name
        : KVPInstr(name), mNtlStartDts(ntlStartDates), mNtlEndDts(ntlEndDates), mNtlAmts(ntlAmts), mPayType(payType)  
{
static        char        routine[] = "KVPProtLeg::KVPProtLeg";

        //
        // Check recovery passed
        //
        if ((recovery == NULL)         ||
            (!strlen(recovery) != 0)   ||
            (!strcmp(recovery, "NIL")) ||
            (!strcmp(recovery, "nil")) ||
            (!strcmp(recovery, "Nil")) )
        {
                mIsDefRecovery = true;
        }
        else
        {
            mIsDefRecovery = false;
            sscanf(recovery, "%lf", &mRecovery);
        }

        mDiscZcName = String(discZcName);

}



//---------------------------------------------------------------
// Only called after construction

void
KVPProtLeg::CheckValid() const
{
static        char        routine[] = "KVPProtLeg::CheckValid";

    try {

        /* Notional and dates vectors must be same positive size */
        ASSERT_OR_THROW(mNtlAmts.size()>=1);
        ASSERT_OR_THROW(mNtlStartDts.size()==mNtlAmts.size());
        ASSERT_OR_THROW(mNtlEndDts.size()==mNtlAmts.size());

        /* Check that dates are all in strict order */
        int i;
        for (i=0; i<mNtlStartDts.size(); i++) {
            ASSERT_OR_THROW(mNtlStartDts[i]<mNtlEndDts[i]);
            if (i>0) {
                ASSERT_OR_THROW(mNtlStartDts[i]>=mNtlEndDts[i-1]);
            }
        }

        /* Recovery rate must be between 0 and 1 if digital */
        if (!mIsDefRecovery)
        {
            ASSERT_OR_THROW(mRecovery >= 0.);
            ASSERT_OR_THROW(mRecovery <= 1.);
        }
    }
    catch (KFailure) {
        throw KFailure("KVPProtLeg::CheckValid: failed for %s.\n", GetName());
    }

}



//---------------------------------------------------------------

istream&
KVPProtLeg::Get(istream& is, int drw)
{

    try {
        if (drw) {

            mNtlStartDts.insert(mNtlStartDts.end(), getTDate(is, "KVPProtLeg::Get: start date."));
            mNtlEndDts.insert(mNtlStartDts.end(), getTDate(is, "KVPProtLeg::Get: end date."));
            mNtlAmts.insert(mNtlAmts.end(),getDouble(is, "KVPProtLeg::Get: notional."));
            mRecovery = getDouble(is, "KVPProtLeg::Get: recovery.");

        } else {
                throw KFailure("KVPProtLeg::Get: format N/A.\n");
        }


        return(is);
    }
    catch (...) {
        throw KFailure("KVPProtLeg::Get: failed.\n");
    }
}


//---------------------------------------------------------------

ostream&
KVPProtLeg::Put(ostream& os, int indent) const
{

    try {

        CheckValid();

        os << "NAME: `" << GetName() << "'" << endl;

        if (mNtlStartDts.size()==1) {
            os << "Start Date:    " << DrlTDatePrint(NULL, mNtlStartDts[0])  << endl;
            os << "End Date:      " << DrlTDatePrint(NULL, mNtlEndDts[0]) << endl;
            os << "Notional     : " << format("%12.4f",  mNtlAmts[0])  << endl;
        } else {
            os << "Varying Notional Schedule" << endl;
            int i;
            for (i=0; i<mNtlStartDts.size(); i++) {
                os << i << "\t";
                os << DrlTDatePrint(NULL, mNtlStartDts[i]) << "\t";
                os << DrlTDatePrint(NULL, mNtlEndDts[i]) << "\t";
                os << format("%12.4f",  mNtlAmts[i])  << endl;
            }
        }
        if (!mIsDefRecovery)
            os << "Recovery Rate: " << format("%12.4f",  mRecovery)  << endl;
        else
            os << "Recovery Rate: nil" <<  endl;
            
        os << "Pay Conv:      " << mPayType << endl;
        os << "Credit Curve:  " << mDiscZcName << endl;

        return(os);
    }
    catch (KFailure) {
        throw KFailure("KVPProtLeg::Get: failed for %s.\n", GetName());
    }
}




//---------------------------------------------------------------
// Only allowed one rate index (with spread schedule), 
// So only the first dependent rate index is printed.
// Other rate indices are assumed only differ by constant spreads.
// 
ostream&
KVPProtLeg::YacctionWriteRecursive(ostream& os)
{
        YacctionWrite(os);
 
        return(os);
}




//---------------------------------------------------------------
// Only allowed one rate index (with spread schedule), 
// and one payment formula
//

ostream&
KVPProtLeg::YacctionWrite(ostream& os, int indent) 
{

    try {

        if (GetWriteFlag())
        {
            int i;
            CheckValid();

            os << GetName() << "=PROTLEG(" << endl;
            for (i=0; i<mNtlStartDts.size(); i++) {
                os << format(" %10s, %10s, %12.4f",
                            GtoFormatDate(mNtlStartDts[i]),
                            GtoFormatDate(mNtlEndDts[i]),
                            mNtlAmts[i]) << "," << endl;
            }

            if (mIsDefRecovery)
                os << "  \"nil\", ";
            else
                os << "\"" << mRecovery << "\", ";

            if (mPayType == PAY_DEF)
                os << "PAY_DEF, " << endl;
            else
                os << "PAY_MAT," << endl;

            os << "\t\"" << mDiscZcName << "\");"
                << endl << endl;

            WriteDone();

        }

        return(os);
    }
    catch (KFailure) {
        throw KFailure("KVPProtLeg::YacctionWrite: failed for %s.\n", GetName());
    }
}

