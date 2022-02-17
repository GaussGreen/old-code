#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/FD1D.hpp"
#include "edginc/FD1DSolver.hpp"
#include "edginc/LatticeProdEDR.hpp"
#include "edginc/IInstrumentCollection.hpp"

DRLIB_BEGIN_NAMESPACE

//----------------------------------------------------------------

bool FD1D::acceptFactor( IMarketFactor * factor )
{
    // only support CAsset type factors
    return dynamic_cast< CAsset * >( factor ) != 0;
}

void FD1D::retrieveFactor()
{
    static const string method = "FD1D::retrieveFactor";
    try
    {
        FDModel::retrieveFactor();

        if( factors.size() != 1 )
            throw ModelException( method, "only 1 factor suppported" );

        underlying = CAssetConstSP::dynamicCast( factors[ 0 ] );
        if( ! underlying )
            throw ModelException( method, "only asset underlier suppported" );
    }
    catch( exception & e )
    {
        throw ModelException( e, method );
    }
}

void FD1D::initModel()
{
    static const string method = "FD1D::initModel";
    try
    {
        stepsPerYearFD = stepsPerYear; // set FDModel::stepsPerYearFD
        FDModel::initModel();

        convertInputStrings();

        // create working range
        range = TreeSliceEQ::Range::create( 0, dim1 - 1 );

        // create underlying Value and solver Variable slices
        xValue = STATIC_POINTER_CAST< TreeSliceEQ >( createSlice( "", factors[ 0 ]->getName() ) );
        xVar = STATIC_POINTER_CAST< TreeSliceEQ >( xValue->clone( false ) );

        // allocate stepForward
        stepForward = DoubleArray( timeLine->NumOfStep + 1, 1. );
    }
    catch( exception & e )
    {
        throw ModelException( e, method );
    }
}

void FD1D::finaliseModel( CControl * control )
{
    static const string method = "FD1D::finaliseModel";
    try
    {
        // use same grid for greeks
        if( isRebuilt( control ) )
        {
            initGrid();

            // cache solver
            solver = IFDSolverSP( new FD1DSolver( *this ) );
        }

        // initialize xValue from xVar
        initValue( timeLine->NumOfStep );
    }
    catch( exception & e )
    {
        throw ModelException( e, method );
    }
}

//----------------------------------------------------------------

bool FD1D::insertNode( int step, int index, double x, double x1, double x2 )
{
    // ignore if level is close enough to an already existing one
    if( x < x1 * ( 1. + FP_MIN ) || x2 * ( 1. - FP_MIN ) < x )
        return false;

    // when node is inserted grid becomes variable
    isVariableGrid = true;

    TreeSliceEQ::Range::InsNodeSP insNode( new TreeSliceEQ::Range::InsNode(
        index, ( x - x2 ) / ( x1 - x2 ), ( x - x1 ) / ( x2 - x1 ) ) );
    // insert node in 1st dimension
    range->insertNode( 0, insNode );

    // insert new node into xValue right away
    xValue->insertNodes();

    // insert new node into xVar right away
    xVar->insertNodes();
    // set value of xVar at inserted node
    switch( changeOfVar )
    {
        case X:
            xVar[ index ] = x;
            break;
        case LOG_X:
            xVar[ index ] = ::log( x );
            break;
        case LOG_FWDX:
            xVar[ index ] = ::log( stepForward[ 0 ] / stepForward[ step ] * x );
            break;
    }

    return true;
}

inline
void FD1D::insertNodesBelow( int count, int step, int index, double x, double x1, double x2 )
{
    x2 = x;
    for( int i = 0; i < count; ++i )
    {
        x = ( x1 + x2 ) / 2.;
        insertNode( step, index, x, x1, x2 );
        x1 = x;
        ++index;
    }
}

inline
void FD1D::insertNodesAbove( int count, int step, int index, double x, double x1, double x2 )
{
    x1 = x;
    ++index;
    for( int i = 0; i < count; ++i )
    {
        x = ( x1 + x2 ) / 2.;
        insertNode( step, index, x, x1, x2 );
        x2 = x;
    }
}

//----------------------------------------------------------------

void FD1D::addCriticalLevel(
    int step,
    const TreeSlice & base,
    double level,
    TreeSlice & target,
    LevelKind kind,
    double value )
{
    // return if no nodes should be insterted
    if( ! nodesPerCritLevel )
        return;

    // fail if critical level is not based on xValue
    if( &base != xValue.get() )
        throw ModelException( "FD1D::addCriticalLevel", "critical level must be based on underlying" );

    // assuming values in xValue are rising monotonicaly...

    // get current range
    int bot1 = range->limits.bot[ 0 ];
    int top1 = range->limits.top[ 0 ];

    double x = level;

    // return if critical level is outside the value range
    if( x < xValue[ bot1 ] || xValue[ top1 ] < x )
        return;

    // find nodes adjacent to critical level
    int lo = bot1;
    int hi = top1;
    double xAdj = x * ( 1. - FP_MIN );
    while( hi - lo > 1 )
    {
        int mid = ( hi + lo ) / 2;
        if( xAdj > xValue[ mid ] )
            lo = mid;
        else
            hi = mid;
    }

    int index = hi;
    double x1 = xValue[ lo ];
    double x2 = xValue[ hi ];

    TreeSliceEQ & xTarget = static_cast< TreeSliceEQ & >( target );

    // insert all nodes added prior
    xTarget.insertNodes();

    // set new boundary and value at it
    if( kind == LEVEL_BARRIER_DOWN )
    {
        xTarget.setBotBound( 0, index );
        xTarget.setBotValue( 0, value );
    }
    else if( kind == LEVEL_BARRIER_UP )
    {
        xTarget.setTopBound( 0, index );
        xTarget.setTopValue( 0, value );
    }

    // first insert critical level
    if( insertNode( step, index, x, x1, x2 ) )
    {
        // insert additional nodes
        int nodeCount = nodesPerCritLevel - 1;
        if( kind == LEVEL_BARRIER_DOWN )
        {
            insertNodesAbove( nodeCount, step, index, x, x1, x2 );
        }
        else if( kind == LEVEL_BARRIER_UP )
        {
            insertNodesBelow( nodeCount, step, index, x, x1, x2 );
        }
        else if( kind == LEVEL_GENERIC )
        {
            nodeCount /= 2;
            insertNodesBelow( nodeCount, step, index, x, x1, x2 );
            insertNodesAbove( nodeCount, step, index, x, x1, x2 );
        }
    }

    // insert new nodes into target slice
    xTarget.insertNodes();
}

//----------------------------------------------------------------

FD1D::FD1D( const CClassConstSP & type ) :
    LatticeModelEQ( type ),
    truncation1D( 4. ),
    dim1( 51 ),
    isVariableGrid( false ),
    changeOfVarStr( "LOG_X" ),
    varStepSpacingStr( "NONE" ),
    nodesPerCritLevel( 1 )
{
    sameGridTweak = false; // this must be false here !!!
    isFwdInduction = false; // only backward for now !!!
}

void FD1D::validatePop2Object()
{
    static const string method = "FD1D::validatePop2Object";
    try
    {
        if( truncation1D < 1. || 15. < truncation1D )
            throw ModelException( method, "truncation1D must be between 1. and 15.");

        if( dim1 < 10 )
            throw ModelException( method, "dim1 should be at least 10");

        if( stepsPerYear < 10 )
            throw ModelException( method, "stepsPerYear should be at least 10");

        if( nodesPerCritLevel < 0 )
            throw ModelException( method, "nodesPerCritLevel should be at least 0");
    }
    catch( exception & e )
    {
        throw ModelException( e, method );
    }
}

// convert input string to enum type
void FD1D::convertInputStrings()
{
    static const string method = "FD1D::convertInputStrings";

    if( changeOfVarStr == "X" )
        changeOfVar = X;
    else if( changeOfVarStr == "LOG_X" )
        changeOfVar = LOG_X;
    else if( changeOfVarStr == "LOG_FWDX" )
        changeOfVar = LOG_FWDX;
    else
    {
        throw ModelException( "FD1D::convertInputStrings",
                              "Unknown change of variable '" + changeOfVarStr + "', "
                              "Only 'X', 'LOG_X' and 'LOG_FWDX' are supported" );
    }

    if( varStepSpacingStr == "NONE" )
        varStepSpacing = NONE;
    else if( varStepSpacingStr == "STD" )
        varStepSpacing = STD;
    else
    {
        throw ModelException( "FD1D::convertInputStrings",
                              "Uknown variable step spacing '" + varStepSpacingStr + "', "
                              "Only 'NONE' and 'STD' are supported" );
    }
}

//----------------------------------------------------------------

void FD1D::initGrid()
{
    static const string method = "FD1D::initGrid";
    try
    {
        double botX, topX;
        setFdBounds( truncation1D, botX, topX );

        int bot1 = range->limits.bot[ 0 ];
        int top1 = range->limits.top[ 0 ];

        double dx = ( topX - botX ) / ( top1 - bot1 );

        for( int i = bot1; i <= top1; ++i )
            xVar[ i ] = botX + dx * i;

        switch( varStepSpacing )
        {
            case NONE:
                break;

            case STD:
                if( changeOfVar == LOG_X )
                {
                    isVariableGrid  = true;
                    setVarSpaceSteps( bot1, top1 );
                }
                else
                    throw ModelException( method, "'STD' step spacing only implemented for 'LOG_X' change of variable" );
                break;

            default:
                throw ModelException( method, "Step spacing '" + varStepSpacingStr + "' not implemented");
        }
    }
    catch( exception & e )
    {
        throw ModelException( e, method );
    }
}

void FD1D::setVarSpaceSteps( int bot1, int top1 )
{
    int n = top1 - bot1;

    //hard code for now
    double trunc5 = truncation1D;
    double trunc3 = 3. / 5. * truncation1D;
    double trunc2 = 2. / 5. * truncation1D;

    double botX5, topX5, botX3, topX3, botX2, topX2;

    setFdBounds( trunc5, botX5, topX5 );
    setFdBounds( trunc3, botX3, topX3 );
    setFdBounds( trunc2, botX2, topX2 );

    // 10% of nodes inside trunc5 outside trunc3
    // 20% of nodes inside trunc3 outside trunc2
    // 70% of nodes inside trunc2
    double w5 = .05; // 5% on each side
    double w3 = .1;  // 10% on each side

    int n5 = int( w5 * n );
    int n3 = int( w3 * n );
    int n2 = n - 2 * n5 - 2 * n3;

    int i = bot1 + 1;

    double dx = ( botX3 - botX5 ) / ( n5 + 1 );
    for( ; i <= bot1 + n5; ++i )
        xVar[ i ] = xVar[ i - 1 ] + dx;

    dx = ( botX2 - botX3 ) / n3;
    for( ; i <= bot1 + n3 + n5; ++i )
        xVar[ i ] = xVar[ i - 1 ] + dx;

    dx = ( topX2 - botX2 ) / n2;
    for( ; i <= bot1 + n2 + n3 + n5; ++i )
        xVar[ i ] = xVar[ i - 1 ] + dx;

    dx = ( topX3 - topX2 ) / n3;
    for( ; i <= bot1 + n3 + n2 + n3 + n5; ++i )
        xVar[ i ] = xVar[ i - 1 ] + dx;

    dx = ( topX5 - topX3 ) / n5;
    for( ; i <= bot1 + n5 + n3 + n2 + n3 + n5; ++i )
        xVar[ i ] = xVar[ i - 1 ] + dx;
}

//----------------------------------------------------------------

void FD1D::initValue( int step ) const
{
    switch( changeOfVar )
    {
        case X:
            *xValue = *xVar;
            break;
        case LOG_X:
            *xValue = exp( *xVar );
            break;
        case LOG_FWDX:
            *xValue = stepForward[ step ] / stepForward[ 0 ] * exp( *xVar );
            break;
    }
}

//----------------------------------------------------------------

void FD1D::updateValue( int step ) const
{
    if( changeOfVar == LOG_FWDX )
        *xValue = stepForward[ step ] / stepForward[ step + 1 ] * *xValue;
}

//----------------------------------------------------------------

FDProductSP FD1D::makeProduct(const IProdCreatorSP & creator)
{
    const IndexSpecEQ * spec = dynamic_cast< const IndexSpecEQ * >( creator.get() );
    if( spec )
    {
        if( spec->getFactor()->getName() == factors[0]->getName() )
        {
            // avoid recursion: use FDModel::makeProduct instead of FDModel::createProduct
            // FDModel::createProduct calls makeProduct back
            return FDModel::makeProduct( IProdCreatorSP( new Spot() ) );
        }
        else
        {
            throw ModelException( "FD1D::makeProduct",
                "IndexSpec factor is '" + spec->getFactor()->getName() + "', "
                "while '" + factors[0]->getName() + "' is expected" );
        }
    }

    return FDModel::makeProduct( creator );
}

//----------------------------------------------------------------

FDModel::IFDSolverSP FD1D::createSolver()
{
    return solver;
}

//----------------------------------------------------------------

double FD1D::getPrice0( const TreeSlice & price ) const
{
    const string method = "FD1D::getPrice0";
    try
    {
        const TreeSliceLayer * layer = dynamic_cast< const TreeSliceLayer * >( &price );
        if( layer )
            return layer->getPrice0();

        double s0 = stepForward[ 0 ];
        double * sVal = xValue;
        double * pVal = price.getValues();

        double price0, deltaDummy, gammaDummy;
        int done = FDInterpolationD(
            range->limits.top[ 0 ] - range->limits.bot[ 0 ] + 1,
            sVal, pVal, 1, &s0, &price0, &deltaDummy, &gammaDummy );

        return price0;
    }
    catch( exception & e )
    {
        throw ModelException( e, method );
    }
}

//----------------------------------------------------------------

void FD1D::load( CClassSP & type )
{
    REGISTER(FD1D, type);
    SUPERCLASS(FDModel);

    FIELD(truncation1D, "num of stdev to truncate for 1st dimention");
    FIELD(dim1, "num of point in FD grid for 1st dimension");
    FIELD(stepsPerYear, "num of time points per year");

    FIELD_USING_ALIAS(changeOfVarStr, changeOfVar, "X: underlying, LOG_X: log(underlying), LOG_FWDX: log(fwd underlying)");
    FIELD_MAKE_OPTIONAL(changeOfVar);

    FIELD_USING_ALIAS(varStepSpacingStr, varStepSpacing, "NONE: fixed spacing, STD: variable spacing");
    FIELD_MAKE_OPTIONAL(varStepSpacing);

    FIELD(nodesPerCritLevel, "number of nodes inserted per critical level");
    FIELD_MAKE_OPTIONAL(nodesPerCritLevel);

    FIELD(DEBUG_DumpToFile, "file name to dump slices to");
    FIELD_MAKE_OPTIONAL(DEBUG_DumpToFile);
}

CClassConstSP const FD1D::TYPE =
    CClass::registerClassLoadMethod( "FD1D", typeid( FD1D ), load );

DRLIB_END_NAMESPACE
