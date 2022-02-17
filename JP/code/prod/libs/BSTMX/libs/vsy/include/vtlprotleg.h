/****************************************************************
 * Module:        VirtualProduct
 * Submodule:        
 * File:          vtlprotleg.h
 * Function:        
 * Author:        David Liu
 *****************************************************************/
#ifndef        _vtlprotleg_H
#define        _vtlprotleg_H

#include "vtlbase.h"
#include "vpprotleg.h"


//--------------------------------------------------------------
/**
 * Valuation class for protection leg
 */

class KVPToolProtLeg : public KVPToolAtom {
public:
        /**
         * Default Constructor
         */
        KVPToolProtLeg(KVTree &vt) : KVPToolAtom(vt) {};

        /**
         * Constructor.
         */
        KVPToolProtLeg(
                SharedPointer<KVPProtLeg> ins,        // (I) inetrument to value
                KVTree &vt);                          // (I) virtual tree

        /**
         * Destructor.
         */
virtual        ~KVPToolProtLeg();


        /**
         * Type name.
         */
virtual        const char*        TypeName() const {return("KVPToolProtLeg");};

        /**
         * Returns the name for the geometry of the slice
         * (see description of method in KVPToolAtom).
         */
virtual        const String&        GetCurveName() { return mProtCurveName;}


        /**
         * Returns the KVPAtom to which the KVPToolAtom is associated.
         */
virtual        SharedPointer<KVPAtom>        Atom()
        {
                SharedPointer<KVPAtom> vp;
                SharedPointerConvertTo(mProtLeg, vp);
                return (vp);
        }

        /**
         * Update tool.
         */
virtual        void        Update();


        /**
         * Returns the current value of the stream.
         */
virtual        KTSlice&    GetValue();



protected:
                     /** Pointer to instrument. */
        SharedPointer<KVPProtLeg>        mProtLeg;


                     /** Effective start date, MAX(startDate, today) */
        //TDate        mEffStDate;    

                     /** Default value
                      * This is the default value assigned to
                      * the component, to avoid return error in GetValue.
                      */
        KTSlice                        *mDefValue;
                     /** Value Slice */
        KTSlice                        *mValue;
        /** Temp slice */
        KTSlice                        *mTmpSlice;
                     /** Zero-resets for protection end dates */
        KVector(KZeroReset*)  mProtZero;
                     /** Prot name for geometry of the slice */
        String                         mProtCurveName;

};


#endif




