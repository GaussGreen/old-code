//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LatticeProdEDR.hpp
//
//   Description : bridging class from EDR to new FD Product
//                 this is temparory and to be reviewed and removed in future.
//
//   Author      : Ning Shen
//
//----------------------------------------------------------------------------

#ifndef LATTICE_PROD_EDR_HPP
#define LATTICE_PROD_EDR_HPP

#include "edginc/FDProduct.hpp"
#include "edginc/LatticeModelEDR.hpp"
#include "edginc/Tree1f.hpp"
#include "edginc/FD1DRet.hpp"
#include "edginc/FD1DMulti.hpp"

DRLIB_BEGIN_NAMESPACE


// Interface required of product to support LogNormal models, to be reviewed
class TREE_DLL IFDProductLN {
public:
    virtual CVolRequestLNSP getVolInterp(int iAsset) const = 0;
    virtual ~IFDProductLN(){}
};

/////////////////////// for bridging to new FDProduct from EDR //////////////
class TREE_DLL ILatticeProdEDR
{
public:
    /** !!! this is to make EDR legacy trees work */
    //for set up timeline only,  to be reviewed/removed */
    virtual CVolRequestConstSP GetLNRequest() const = 0 ;
    // to be reviewed/removed

    // map inserted node and prices to slice values for legacy EDR tree
    virtual void mapInsertNode(
        TreeSliceGeneral::RangeSP & range,
        TreeSliceGeneralSP & nodes,
        TreeSliceGeneralContSP & prices,
        int ** priority) = 0;

    /** To make it possible to change InsertNode level,*/
    virtual bool moveInsertNode(int currStep, int iPrice) = 0;
    /** !!! this is to make EDR legacy trees work */
};

/////////////////////////////
// basic implementation
class TREE_DLL LatticeProdEDR : public FDProduct, public ILatticeProdEDR
{
protected:
    LatticeProdEDR( FDModel * model ) :
        FDProduct( model ),
        tree1f( dynamic_cast< CTree1f * >( model ) ),
        fd1dRet( dynamic_cast< FD1DRet * >( model ) ),
		fd1dMulti( dynamic_cast< FD1DMulti * >( model ) )
    {}

    /** !!! this is to make EDR legacy trees work */
    //for set up timeline only,  to be reviewed/removed */
    virtual CVolRequestConstSP GetLNRequest() const { return CVolRequestConstSP(); }

    // map inserted node and prices to slice values for legacy EDR tree
    virtual void mapInsertNode(
        TreeSliceGeneral::RangeSP & range,
        TreeSliceGeneralSP & nodes,
        TreeSliceGeneralContSP & prices,
        int ** priority)
    {}

    /** To make it possible to change InsertNode level,*/
    virtual bool moveInsertNode(int currStep, int iPrice) { return false; }

    /** get product value */
    virtual const TreeSlice & getValue( int step ) const { return *slices[0]; }

    /** create slices for multiple payoffs */
    void initSlices( int size, const string & curveToDEV = "" )
    {
        slices.resize( size );
        for( int i = 0; i < size; ++i )
            startDEV( slices[ i ] = model->createSlice( curveToDEV ) );
    }

    /** return sices values as vector of double pointers */
    static const vector< double * > getValues( const vector< TreeSliceSP > & slices )
    {
        int size = slices.size();
        vector< double * > values( size );
        for( int i = 0; i < size; ++i )
            values[ i ] = slices[ i ]->getValues();

        return values;
    }

    // down-casted to Tree1f
    CTree1f* tree1f;
    // down-casted to FD1DRet
    FD1DRet* fd1dRet;
	// down-casted to FD1DMulti
    FD1DMulti* fd1dMulti;

    // one and only payoff index (usually spot)
    FDProductSP payoffIndex;

    vector< TreeSliceSP > slices;
};
/////////////////////////////
// basic implementation

/////////////////////////////
// inserted node for product with barriers, a helper class for Tree only
class TREE_DLL LatticeProdEDRIns : public LatticeProdEDR
{
protected:
    LatticeProdEDRIns( FDModel * model, int nIns, int nPrices );
    ~LatticeProdEDRIns();

    // creates inserted node slices
    void initInsertNode();

    // map inserted node and prices to slice values for legacy EDR tree
    virtual void mapInsertNode(
        TreeSliceGeneral::RangeSP & range,
        TreeSliceGeneralSP & nodes,
        TreeSliceGeneralContSP & prices,
        int ** priority);

    //number of insert nodes
    int numIns;
    //number of price elements
    int numPrices;

    // for node insertion, we now keep these in products
    TreeSliceGeneral::RangeSP insRange;  // inserted nodes range
    TreeSliceGeneralSP        insNodes;  // inserted stock levels
    TreeSliceGeneralContSP    insPrices; // price arrays for inserted node
    vector< int * >           insOrders; // priority of inserted nodes if more than one found within branching range
};

inline LatticeProdEDRIns::LatticeProdEDRIns( FDModel * model, int nIns, int nPrices ) :
    LatticeProdEDR( model ),
    numIns( nIns ),
    numPrices( nPrices )
{}

inline void LatticeProdEDRIns::initInsertNode()
{
    if( numIns )
    {
        const LatticeModelEDR * modelEDR = dynamic_cast< const LatticeModelEDR * >( model );
        if( ! modelEDR )
            return;

        int numSlices = modelEDR->getSliceCount();

        insRange  = TreeSliceGeneral::Range::create( numSlices, 0, numIns - 1 );
        int dimBits = (1<<(*insRange)->nDim)-1;
        insNodes  = TreeSliceGeneral::create( *insRange, dimBits );
        insPrices = TreeSliceGeneralCont::create( *insRange, numPrices, dimBits );

        insOrders.resize( numSlices );
        for( int i = 0; i < (int)numSlices; ++i )
            insOrders[ i ] = new int[ numIns ];
    }
}

inline LatticeProdEDRIns::~LatticeProdEDRIns()
{
    if( numIns )
    {
        for( int i = 0; i < (int)insOrders.size(); ++i )
            delete [] insOrders[ i ];
    }
}

inline void LatticeProdEDRIns::mapInsertNode(
    TreeSliceGeneral::RangeSP & range,
    TreeSliceGeneralSP & nodes,
    TreeSliceGeneralContSP & prices,
    int ** priority)
{
    if( ! model || ! numIns )
        throw ModelException( "LatticeProdEDRIns::mapInsertNode", "no insert nodes to map" );

    range = insRange;
    nodes = insNodes;
    prices = insPrices;

    if( tree1f )
    {
        for( int i = 0; i < (int)insOrders.size(); ++i )
            priority[ i ] = insOrders[ i ];
    }
}
/////////////////////////////
// inserted node for product with barriers, a helper class

DRLIB_END_NAMESPACE

#endif
