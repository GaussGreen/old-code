//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : LatticeModelEDR.hpp
//
//   Description : EDR base FDModel class.
//
//   Date        : Apr 12, 2006
//
//----------------------------------------------------------------------------

#ifndef LATTICE_MODEL_EDR_HPP
#define LATTICE_MODEL_EDR_HPP

#include "edginc/LatticeModelEQ.hpp"

DRLIB_BEGIN_NAMESPACE

class LatticeModelEDR : public FDModel
{
public:
    int getSliceCount() const { return range.get() ? range->count() : 1; }
    int getSliceIndex( int step ) const { return step % getSliceCount(); }

    int getDimBit( const string & factorName ) const
    {
        if( factorName == "" )
            return 0;

        for( int i = 0; i < factors.size(); ++i )
        {
            if ( factorName == factors[ i ]->getName() )
                return ( 1 << i );
        }

        throw ModelException( "LatticeModelEDR::getDimBit", factorName + " factor is not supported" );
    }

    /** creates new slice */
    virtual TreeSliceSP createSlice(
        const string & curveToDEV = "",
        const string & factorName1 = "",
        const string & factorName2 = "",
        const string & factorName3 = "" ) const
    {
        int dimBits = getDimBit( factorName1 ) | getDimBit( factorName2 ) | getDimBit( factorName3 );

        // state variable support
        return createLayer( createSlice( dimBits ) );
    }

    /** get date at step on the time line */
    double getTrdYrFrac( int step ) const
    {
        return timeLine->TradeYrFrac[step];
    }

protected:
    TreeSliceGeneral::RangeSP range;

    virtual TreeSliceSP createSlice( int dimBits ) const
    {
        return TreeSliceGeneral::create( *range, dimBits ? dimBits : ( 1 << (*range)->nDim ) - 1 );
    }

    LatticeModelEDR( const CClassConstSP type ) : FDModel( type ) {}
};

DRLIB_END_NAMESPACE

#endif
