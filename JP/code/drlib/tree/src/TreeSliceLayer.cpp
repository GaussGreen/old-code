#include "edginc/config.hpp"
#include "edginc/TreeSlice.hpp"
#include "edginc/FDModel.hpp"

DRLIB_BEGIN_NAMESPACE

/**********************************************************************************************************/
namespace
{

typedef TreeSliceLayer::StateSupport::InterpMethod InterpMethod;

//----------------------------------------------------------------------------
// slices interpolation operator
//----------------------------------------------------------------------------
template< InterpMethod interpMethod >
class InterpolateOper
{
    const TreeSliceLayer::StateSupport & owner;

    const int nbGridIn;
    const int nbGridOut;
    const int nbLayers;

    const vector< double > & transitionIn;
    mutable vector< double > transitionOut;

    const vector< double > & gridIn;
    const vector< double > & gridOut;
    const vector< TreeSliceLayer * > & layers;
    const vector< const TreeSlice * > & srcSlices;

    vector< vector< TreeSliceSP > * > slicesIn;
    vector< vector< TreeSliceSP > * > slicesOut;

public:
    InterpolateOper(
        const TreeSliceLayer::StateSupport & owner,
        bool isFwdTransition,
        const vector< double > & prevGrid,
        const vector< double > & currGrid,
        const vector< TreeSliceLayer * > & layers,
        const vector< const TreeSlice * > & srcSlices )
    :
        owner( owner ),
        nbGridIn( prevGrid.size() ),
        nbGridOut( currGrid.size() ),
        nbLayers( layers.size() ),
        transitionIn( isFwdTransition ? currGrid : prevGrid ),
        transitionOut( transitionIn.size() ),
        gridIn( isFwdTransition ? prevGrid : transitionOut ),
        gridOut( isFwdTransition ? transitionOut : currGrid ),
        layers( layers ),
        srcSlices( srcSlices ),
        sliceCount( nbLayers * ( nbGridIn + nbGridOut ) + srcSlices.size() )
    {
        ASSERT( nbGridIn && nbGridOut );

        slicesIn.resize( nbLayers );
        slicesOut.resize( nbLayers );
        for( int l = 0; l < nbLayers; ++l )
        {
            slicesOut[ l ] = &layers[ l ]->getSlices();
            slicesIn[ l ] = &layers[ l ]->getTempSlices();
            swapT( *slicesIn[ l ], *slicesOut[ l ] );

            TreeSlice::resizeVector(
                nbGridOut,
                *(*slicesIn[ l ])[ 0 ], *slicesOut[ l ], layers[ l ]->name );

            // set correct treeStep on interpolated slices
            int treeStep = (*slicesIn[ l ])[ 0 ]->treeStep;
            for( int i = 0; i < nbGridOut; ++i )
                (*slicesOut[ l ])[ i ]->treeStep = treeStep;
        }
    }

    // TreeSlice "expression template" primitives
    const int sliceCount;
    template< typename S >
    const S** listInputSlices(const S** list) const
    {
        for( int l = 0; l < nbLayers; ++l )
        {
            for( int i = 0; i < nbGridIn; ++i )
               list = (*slicesIn[ l ])[ i ]->listSlices( list );
        }

        for( int i = 0; i < (int)srcSlices.size(); ++i )
            list = srcSlices[ i ]->listSlices( list );

        return list;
    }
    template< typename S >
    S** listOutputSlices(S** list) const
    {
        for( int l = 0; l < nbLayers; ++l )
        {
            for( int i = 0; i < nbGridOut; ++i )
                *list++ = (*slicesOut[ l ])[ i ].get();
        }

        return list;
    }

    void compute1() const;
    void compute2() const;
    void compute3() const;

    void compute() const;

    void printDebug(char *s) const { strcat(s, "InterpolateOper"); }
};
template< InterpMethod interpMethod >
inline
void InterpolateOper< interpMethod >::compute1() const
{
    for( int i = 0; i < nbGridOut; ++i )
    {
        for( int l = 0; l < nbLayers; ++l )
        {
            *(*slicesOut[ l ])[ i ]->iter =
                *(*slicesIn[ l ])[ 0 ]->iter;
        }
    }
}
template< InterpMethod interpMethod >
inline
void InterpolateOper< interpMethod >::compute2() const
{
    if( gridIn.size() < 2 )
        compute1();

    owner.transition( srcSlices, transitionIn, transitionOut );

    double x1 = gridIn[ 0 ];
    double x2 = gridIn[ 1 ];

    int k = 1;
    for( int i = 0; i < nbGridOut; ++i )
    {
        double x = gridOut[ i ];

        while( k < nbGridIn - 1 && x2 < x )
        {
            x1 = x2;
            x2 = gridIn[ ++k ];
        }

        double coeff1;
        double coeff2;

        if( Maths::isZero( x1 - x2 ) )
        {
            coeff1 = 1.;
            coeff2 = 0.;
        }
        else
        {
            coeff1 = ( x - x2 ) / ( x1 - x2 );
            coeff2 = ( x - x1 ) / ( x2 - x1 );
        }

        for( int l = 0; l < nbLayers; ++l )
        {
            *(*slicesOut[ l ])[ i ]->iter =
                coeff1 * *(*slicesIn[ l ])[ k - 1 ]->iter +
                coeff2 * *(*slicesIn[ l ])[ k     ]->iter;
        }
    }
}
template< InterpMethod interpMethod >
inline
void InterpolateOper< interpMethod >::compute3() const
{
    if( gridIn.size() < 3 )
        compute2();

    owner.transition( srcSlices, transitionIn, transitionOut );

    int i = 0;
    double x;

    double x1 = gridIn[ 0 ];
    double x2 = gridIn[ 1 ];
    double x3 = gridIn[ 2 ];

    while( i < nbGridOut && ( x = gridOut[ i ] ) < x1 )
    {
        double coeff1;
        double coeff2;

        if( Maths::isZero( x1 - x2 ) )
        {
            coeff1 = 0.;
            coeff2 = 1.;
        }
        else
        {
            coeff1 = ( x - x2 ) / ( x1 - x2 );
            coeff2 = ( x - x1 ) / ( x2 - x1 );
        }

        for( int l = 0; l < nbLayers; ++l )
        {
            *(*slicesOut[ l ])[ i ]->iter =
                coeff1 * *(*slicesIn[ l ])[ 0 ]->iter +
                coeff2 * *(*slicesIn[ l ])[ 1 ]->iter;
        }

        ++i;
    }

    double xn = gridIn[ nbGridIn - 1 ];

    int k = 2;
    while( i < nbGridOut && ( x = gridOut[ i ] ) < xn )
    {
        while( k < nbGridIn - 1 && x2 < x )
        {
            x1 = x2;
            x2 = x3;
            x3 = gridIn[ ++k ];
        }

        double coeff1;
        double coeff2;
        double coeff3;

        if( Maths::isZero( x1 - x2 ) || Maths::isZero( x2 - x3 ) || Maths::isZero( x3 - x1 ) )
        {
            if( x < x2 )
            {
                if( Maths::isZero( x1 - x2 ) )
                {
                    coeff1 = 0.;
                    coeff2 = 1.;
                    coeff3 = 0.;
                }
                else
                {
                    coeff1 = ( x - x2 ) / ( x1 - x2 );
                    coeff2 = ( x - x1 ) / ( x2 - x1 );
                    coeff3 = 0.;
                }
            }
            else
            {
                if( Maths::isZero( x2 - x3 ) )
                {
                    coeff1 = 0.;
                    coeff2 = 1.;
                    coeff3 = 0.;
                }
                else
                {
                    coeff1 = 0.;
                    coeff2 = ( x - x3 ) / ( x2 - x3 );
                    coeff3 = ( x - x2 ) / ( x3 - x2 );
                }
            }
        }
        else
        {
            coeff1 = ( x - x2 ) * ( x - x3 ) / ( x1 - x2 ) / ( x1 - x3 );
            coeff2 = ( x - x3 ) * ( x - x1 ) / ( x2 - x3 ) / ( x2 - x1 );
            coeff3 = ( x - x1 ) * ( x - x2 ) / ( x3 - x1 ) / ( x3 - x2 );
        }

        for( int l = 0; l < nbLayers; ++l )
        {
            *(*slicesOut[ l ])[ i ]->iter =
                coeff1 * *(*slicesIn[ l ])[ k - 2 ]->iter +
                coeff2 * *(*slicesIn[ l ])[ k - 1 ]->iter +
                coeff3 * *(*slicesIn[ l ])[ k     ]->iter;
        }

        ++i;
    }

    x1 = gridIn[ nbGridIn - 2 ];
    x2 = xn;

    while( i < nbGridOut )
    {
        x = gridOut[ i ];

        double coeff1;
        double coeff2;

        if( Maths::isZero( x1 - x2 ) )
        {
            coeff1 = 1.;
            coeff2 = 0.;
        }
        else
        {
            coeff1 = ( x - x2 ) / ( x1 - x2 );
            coeff2 = ( x - x1 ) / ( x2 - x1 );
        }

        for( int l = 0; l < nbLayers; ++l )
        {
            *(*slicesOut[ l ])[ i ]->iter =
                coeff1 * *(*slicesIn[ l ])[ nbGridIn - 2 ]->iter +
                coeff2 * *(*slicesIn[ l ])[ nbGridIn - 1 ]->iter;
        }

        ++i;
    }
}
template<>
inline
void InterpolateOper< TreeSliceLayer::StateSupport::INTERP_LINEAR >::compute() const
{
    compute2();
}
template<>
inline
void InterpolateOper< TreeSliceLayer::StateSupport::INTERP_QUADRATIC >::compute() const
{
    compute3();
}
//----------------------------------------------------------------------------

}
//****************************************************************************

void TreeSliceLayer::StateSupport::interpolate( int step, const vector< const TreeSlice * > & srcSlices )
{
    vector< TreeSliceLayer * > layers;
    model->collectLayers( this, layers );

    int nbLayers = layers.size();
    for( int l = 0; l < nbLayers; ++l )
        checkInterpStep( *layers[ l ], true );

    switch( interpMethod )
    {
        case INTERP_LINEAR:
            TreeSlice::loopOnSlices( InterpolateOper< INTERP_LINEAR >( 
                *this, isFwdTransition, prevGrid, currGrid, layers, srcSlices ) );
            break;

        case INTERP_QUADRATIC:
            TreeSlice::loopOnSlices( InterpolateOper< INTERP_QUADRATIC >( 
                *this, isFwdTransition, prevGrid, currGrid, layers, srcSlices ) );
            break;

        default:
            throw ModelException( "TreeSliceLayer::StateSupport::interpPrepare",
                "Unsupported interpolation type" );
    }

    prevGrid = currGrid;
    prevInterpStep = currInterpStep;
    currInterpStep = step;

    for( int l = 0; l < nbLayers; ++l )
        layers[ l ]->lastInterpStep = currInterpStep;
}

void TreeSliceLayer::StateSupport::checkInterpStep( const TreeSliceLayer & layer, bool usePrev ) const
{
    int step = usePrev ? prevInterpStep : currInterpStep;
    if( step >= 0 && layer.lastInterpStep >= 0 && layer.lastInterpStep > step )
    {
        throw ModelException( "TreeSliceLayer::StateSupport::checkInterpStep",
            "TreeSlice " + layer.name + " is not sampled at the values for state variable " + name );
    }
}
    
TreeSliceLayer::StateSupport::StateSupport(
    const string & name,
    double todayValue,
    InterpMethod interpMethod,
    bool isFwdTransition )
:
    todayValue( todayValue ),
    interpMethod( interpMethod ),
    isFwdTransition( isFwdTransition ),
    model( 0 ),
    level( -1 ), 
    currInterpStep( -1 ),  
    prevInterpStep( -1 ),
    name( name )
{}
    
void TreeSliceLayer::StateSupport::init( const FDModel * model, int level )
{
    this->model = model;
    this->level = level;

    gridLayer = DYNAMIC_POINTER_CAST< TreeSliceLayer >( model->createSlice() );
    gridLayer->name = name;
}

void TreeSliceLayer::StateSupport::populateGrid( int step )
{
    int nbCurrGrid = currGrid.size();
    ASSERT( nbCurrGrid );

    vector< TreeSliceSP > & slices = gridLayer->getSlices();
    TreeSlice::resizeVector( nbCurrGrid, *slices[ 0 ], slices, gridLayer->name );

    // set the slices of "grid" to the values of "currGrid"
    for( int i = 0; i < nbCurrGrid; ++i )
    {
        *slices[ i ] = currGrid[ i ];
        // avoid treeStep consistency check for grid slices
        slices[ i ]->treeStep = -1;
#ifdef DEBUG
        slices[ i ]->name = Format::toString( i ) + "(" + Format::toString( currGrid[ i ] ) + ")";
#endif
    }

    gridLayer->lastInterpStep = step;
    // avoid treeStep consistency check for grid slices
    gridLayer->treeStep = -1;
}

double TreeSliceLayer::StateSupport::getPrice0( const TreeSliceLayer & layer ) const
{
    try
    {
        checkInterpStep( layer );

        const vector< TreeSliceSP > & slices = layer.getSlices();
        int nbSlices = slices.size();
        
        if( nbSlices == 1 )
            return model->getPrice0( *slices[ 0 ] );

        ASSERT( nbSlices && nbSlices == (int)currGrid.size() );

        int lo = 0;
        int hi = nbSlices - 1;

        // do linear interpolation at todayValue
        double x = todayValue;

        while( hi - lo > 1 )
        {
            int mid = ( hi + lo ) / 2;
            if( x > currGrid[ mid ] )
                lo = mid;
            else
                hi = mid;
        }

        double x1 = currGrid[ lo ];
        double x2 = currGrid[ hi ];
        ASSERT( ! Maths::isZero( x1 - x2 ) );

        double y1 = model->getPrice0( *slices[ lo ] );
        double y2 = model->getPrice0( *slices[ hi ] );

        return ( x - x2 ) / ( x1 - x2 ) * y1 + ( x - x1 ) / ( x2 - x1 ) * y2;
    }
    catch( exception & e )
    { 
        throw ModelException( e, "TreeSliceLayer::StateSupport::getPrice0" );
    }
}

/******************************** TreeSliceLayer *****************************/

TreeSliceLayer::TreeSliceLayer(
    StateSupport * owner,
    const TreeSlice & slice,
    bool copyValues,
    int lastInterpStep )
    :
    lastInterpStep( lastInterpStep ),
    owner( owner )
{
    treeStep = slice.treeStep;
    const TreeSliceLayer * layer = dynamic_cast< const TreeSliceLayer * >( &slice );
    if( layer && layer->owner == owner )
    {
        int size = layer->slices.size();
        slices.resize( size );
        for( int i = 0; i < size; ++i )
            slices[ i ] = layer->slices[ i ]->clone( copyValues );

        return;
    }
    slices.resize( 1, slice.clone( copyValues ) ); // 1 slice layer is created
}

double TreeSliceLayer::getPrice0() const 
{
    return owner->getPrice0( *this );
}

TreeSliceSP TreeSliceLayer::clone( bool copyValues ) const
{
    return TreeSliceLayerSP( new TreeSliceLayer( owner, *this, copyValues, lastInterpStep ) );
}

TreeSlice& TreeSliceLayer::operator =( double v )
{
    slices.resize(1);
    *slices[0] = v;
    lastInterpStep = -1;
    treeStep = -1;
    return *this;
}

double TreeSliceLayer::getCentre() const {
    return slices[slices.size()/2]->getCentre();
}

void TreeSliceLayer::collectLayers(
    const StateSupport * owner,
    const vector< TreeSliceSP > & slices,
    vector< TreeSliceLayer * > & layers,
    bool expandedOnly )
{
    for( size_t j = 0; j < slices.size(); ++j )
    {
        TreeSliceLayer * layer = dynamic_cast< TreeSliceLayer * >( slices[ j ].get() );
        if( ! layer )
            continue;

        if( layer->owner == owner )
        {
            // check whether to skip layers not expanded in state variable dimension
            if( ! expandedOnly || layer->slices.size() > 1 )
                layers.push_back( layer );
        }
        else
            collectLayers( owner, layer->slices, layers, expandedOnly );
    }
}

bool TreeSliceLayer::isZero() const {
    if (slices.size()==0) 
        return true;

    for (int i=0; i<(int)slices.size(); ++i) {
        if (!slices[i]->isZero())
            return false;
    }
    return true;
}

// !!! TO BE REMOVED
void TreeSliceLayer::getCalcRange( int & bot, int & top ) const {
    throw ModelException("TreeSliceLayer::getCalcRange()", "not implemented");
}
// !!! TO BE REMOVED
double * TreeSliceLayer::getValues() const {
    throw ModelException("TreeSliceLayer::operator *()", "not implemented");
}

void TreeSliceLayer::printDetails(char *s) const {
    TreeSlice::printDetails(s);
    strcat(s, " nbSlices=");
    static char buf[30];
    sprintf(buf,"%d", slices.size());
    strcat(s, buf);
    for (int i=0; i<(int)slices.size(); ++i) {
        strcat(s, "\n");
        strcat(s, "[");
        sprintf(buf,"%d", i);
        strcat(s, buf);
        strcat(s, "] ");
        slices[i]->printDetails(s);
    }
}

DRLIB_END_NAMESPACE
