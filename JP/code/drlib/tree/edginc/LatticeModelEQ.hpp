#ifndef LATTICE_MODEL_EQ_HPP
#define LATTICE_MODEL_EQ_HPP

#include "edginc/FDModel.hpp"
#include "edginc/TreeSlice.hpp"

DRLIB_BEGIN_NAMESPACE

// leave these here for now
const double FP_MIN = 1.0e-10;

class LatticeModelEQ : public FDModel
{
public:
    int getDim( const string & factorName ) const
    {
        for( int i = 0; i < factors.size(); ++i )
        {
            if ( factorName == factors[ i ]->getName() )
                return i;
        }

        throw ModelException( "LatticeModelEQ::getDim", "Factor '" + factorName + "' not found" );
    }

    int getDimBit( const string & factorName ) const
    {
        return factorName == "" ? 0 : ( 1 << getDim( factorName ) );
    }

    virtual TreeSliceSP createSlice(
        const string & curveToDEV = "",
        const string & factorName1 = "",
        const string & factorName2 = "",
        const string & factorName3 = "" ) const
    {
        int dimBits =
            getDimBit( factorName1 ) |
            getDimBit( factorName2 ) |
            getDimBit( factorName3 );

        // state variable support
        return createLayer( createSlice( dimBits ) );
    }

protected:
    TreeSliceEQ::RangeSP range;
    TreeSliceEQRef xValue;
    TreeSliceEQRef xVar;

    TreeSliceEQSP createSlice( int dimBits ) const
    {
        return TreeSliceEQ::create( *range, dimBits /*? dimBits : ( 1 << range->nDim ) - 1*/ );
    }

    LatticeModelEQ( const CClassConstSP & type ) : FDModel( type ) {}

    class TREE_DLL Spot : virtual public FDModel::IIntoProduct, public CObject
    {
    public:
        explicit Spot() : CObject( CClassConstSP() ) {}

        class TREE_DLL Product : public FDProduct
        {
        public:
            Product( FDModel * model ) : FDProduct( model ) {}

        protected:
            virtual bool isElementary() const { return true; }

            virtual void init( Control*  control ) const {}
            virtual void initProd() {}

            virtual void update( int & step, FDProduct::UpdateType type ) {}
            virtual const TreeSlice & getValue( int step ) const
            {
                return *static_cast< const LatticeModelEQ & >( *model ).xValue;
            }

            virtual void recordOutput( Control * control, YieldCurveConstSP disc, Results * results ) {}

            virtual DateTime getStartDate() const { return model->getDate( 0 ); }
        };

    protected:
        virtual FDProductSP createProduct( FDModel * model ) const
        {
            return FDProductSP( new Product( model ) );
        }
    };
};

DRLIB_END_NAMESPACE

#endif
