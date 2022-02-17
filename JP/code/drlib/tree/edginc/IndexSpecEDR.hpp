//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : IndexSpecEDR.hpp
//
//   Description : specification for EDR FD/tree underlying/payoff index
//
//----------------------------------------------------------------------------

#ifndef INDEX_SPEC_EDR_HPP
#define INDEX_SPEC_EDR_HPP

#include "edginc/IndexSpecEQ.hpp"
#include "edginc/FDModel.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL IndexSpecEDR : virtual public FDModel::IIntoProduct, public CObject
{
public:
    explicit IndexSpecEDR( const string & factorName ) :
        CObject( CClassConstSP() ),
        factorName( factorName )
    {}

    class TREE_DLL Product : public FDProduct
    {
    public:    
        Product( FDModel * model, const string & factorName ) :
            FDProduct( model ),
            factorName( factorName ) {calcOff = true;}

    protected:
        virtual bool isElementary() const { return true; }

        virtual void init( Control*  control ) const {}
        virtual void initProd() { startDEV( slice = model->createSlice( "", factorName ) ); }

        virtual void update( int & step, FDProduct::UpdateType type ) {}
        virtual const TreeSlice & getValue( int step ) const { return *slice; }

        virtual void recordOutput( Control * control, YieldCurveConstSP disc, Results * results ) {}

        virtual DateTime getStartDate() const { return model->getDate( 0 ); }

        string factorName;
        TreeSliceSP slice;
    };

protected:
    string factorName;

private:
    virtual FDProductSP createProduct( FDModel * model ) const
    {
        return FDProductSP( new Product( model, factorName ) );
    }
};

DRLIB_END_NAMESPACE

#endif
