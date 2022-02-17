//----------------------------------------------------------------------------
//   Group       : EDG Derivatives Research
//
//   Filename    : FD1DMulti.cpp
//
//   Description : 1-factor finite difference engine to solve multiple PDEs, which have the form of 
//
//					U(i,t,x)_t + a(i)*U(i,t,x)_x + c(i)*U(i,t,x)_xx + f(i)*U(i,t,x) 
//						
//						+ Sum_j{q(i,j)*U(j,t,x + jump(i,j))} = g(i)
//
//				   for i = 1, ..., I
//
//   Author      : Zhijiang Huang
//
//   Date        : Aug 08, 2006
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/FD1DMulti.hpp"
#include "edginc/FD1DMultiSolver.hpp"
#include "edginc/LatticeProdEDR.hpp"
#include "edginc/IInstrumentCollection.hpp"

DRLIB_BEGIN_NAMESPACE

//----------------------------------------------------------------

bool FD1DMulti::acceptFactor( IMarketFactor * factor )
{
    // only support CAsset type factors
    return dynamic_cast< CAsset * >( factor ) != 0;
}

void FD1DMulti::retrieveFactor()
{
    static const string method = "FD1DMulti::retrieveFactor";
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

void FD1DMulti::initModel()
{
    static const string method = "FD1DMulti::initModel";
    try
    {
        stepsPerYearFD = stepsPerYear; // set FDModel::stepsPerYearFD
        FDModel::initModel();

        convertInputStrings();

        // create working range
        range = TreeSliceEQ::Range::create( 0, dim * nPdes - 1 );

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

void FD1DMulti::finaliseModel( CControl * control )
{
    static const string method = "FD1DMulti::finaliseModel";
    try
    {
        if( isRebuilt( control ) )
        {
            // set FD Grids. Use same grids to calc Greeks
            initVarGrid();
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

void FD1DMulti::addCriticalLevel( int step, const TreeSlice & slice, double value )
{
    // ignore if critical level is not based on xValue
    if( &slice != xValue.get() )
        throw ModelException( "FD1DMulti::addCriticalLevel", "critical level must be based on underlying" );

    double level;
    switch( changeOfVar )
    {
        case X:
            level = value;
            break;
        case LOG_X:
            level = log( value );
            break;
        case LOG_FWDX:
            level = log( stepForward[ 0 ] / stepForward[ step ] * value );
            break;
    }

    // values in xVar are rising monotonicaly

    int bot = range->bot[ 0 ];
    int top = range->top[ 0 ];
    double x = level;

    // ignore if critical level is outside the variable range
    if( x < xVar[ bot ] || xVar[ top ] < x )
        throw ModelException( "FD1DMulti::addCriticalLevel", "critical level is ouside the slice range" );

    int k = bot + 1;
    while( k < top && xVar[ k ] < x )
        ++k;

    double x1 = xVar[ k - 1 ];
    double x2 = xVar[ k     ];

    // ignore if level is close enough to an already existing one
    static const double tolerance = 1.e-10;
    if( x - x1 < tolerance || x2 - x < tolerance )
        return;

//!!!    range->insertNode( 0, k, ( x - x2 ) / ( x1 - x2 ), ( x - x1 ) / ( x2 - x1 ) );
}

//----------------------------------------------------------------

FD1DMulti::FD1DMulti( const CClassConstSP & type ) :
    LatticeModelEQ( type ),
    truncation( 4. ),
    dim( 51 ),
    isVariableGrid( false ),
    changeOfVarStr( "LOG_X" ),
    varStepSpacingStr( "NONE" )
{
    sameGridTweak = false; // this must be false here !!!
	isFwdInduction = false;
}

void FD1DMulti::validatePop2Object()
{
    static const string method = "FD1DMulti::validatePop2Object";
    try
    {
        if( truncation < 1. || 15. < truncation )
            throw ModelException( method, "truncation1D must be between 1. and 15.");

        if( dim < 10 )
            throw ModelException( method, "dim1 should be at least 10");

        if( stepsPerYear < 10 )
            throw ModelException( method, "stepsPerYear should be at least 10");
    }
    catch( exception & e )
    {
        throw ModelException( e, method );
    }
}

// convert input string to enum type
void FD1DMulti::convertInputStrings()
{
    static const string method = "FD1DMulti::convertInputStrings";

    if( changeOfVarStr == "X" )
        changeOfVar = X;
    else if( changeOfVarStr == "LOG_X" )
        changeOfVar = LOG_X;
    else if( changeOfVarStr == "LOG_FWDX" )
        changeOfVar = LOG_FWDX;
    else
    {
        throw ModelException( "FD1DMulti::convertInputStrings",
                              "Unknown change of variable '" + changeOfVarStr + "', "
                              "Only 'X', 'LOG_X' and 'LOG_FWDX' are supported" );
    }

    if( varStepSpacingStr == "NONE" )
        varStepSpacing = NONE;
    else if( varStepSpacingStr == "STD" )
        varStepSpacing = STD;
    else
    {
        throw ModelException( "FD1DMulti::convertInputStrings",
                              "Uknown variable step spacing '" + varStepSpacingStr + "', "
                              "Only 'NONE' and 'STD' are supported" );
    }
}

//----------------------------------------------------------------

void FD1DMulti::initVarGrid()
{
    static const string method = "FD1DMulti::initVarGrid";
    try
    {
        double botX, topX;
        setFdBounds( truncation, botX, topX );

        double dx = ( topX - botX ) / ( dim - 1 );

		for( int iPde = 0; iPde < nPdes; iPde++){
			for( int i = 0; i < dim; ++i ){
				xVar[ iPde * dim + i ] = botX + dx * i;
			}
		}

        switch( varStepSpacing )
        {
            case NONE:
                break;

            case STD:
                if( changeOfVar == LOG_X )
                {
                    isVariableGrid  = true;
                    setVarSpaceSteps();
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

void FD1DMulti::setVarSpaceSteps()
{
    int n = dim - 1;

    //hard code for now
    double trunc5 = truncation;
    double trunc3 = 3. / 5. * truncation;
    double trunc2 = 2. / 5. * truncation;

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

    int i = 1;

    double dx = ( botX3 - botX5 ) / ( n5 + 1 );
    for( ; i <= n5; ++i )
        xVar[ i ] = xVar[ i - 1 ] + dx;

    dx = ( botX2 - botX3 ) / n3;
    for( ; i <= n3 + n5; ++i )
        xVar[ i ] = xVar[ i - 1 ] + dx;

    dx = ( topX2 - botX2 ) / n2;
    for( ; i <= n2 + n3 + n5; ++i )
        xVar[ i ] = xVar[ i - 1 ] + dx;

    dx = ( topX3 - topX2 ) / n3;
    for( ; i <= n3 + n2 + n3 + n5; ++i )
        xVar[ i ] = xVar[ i - 1 ] + dx;

    dx = ( topX5 - topX3 ) / n5;
    for( ; i <= n5 + n3 + n2 + n3 + n5; ++i )
        xVar[ i ] = xVar[ i - 1 ] + dx;
}

//----------------------------------------------------------------

void FD1DMulti::initValue( int step ) const
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
            *xValue =
                stepForward[ step ] / stepForward[ 0 ] *
                exp( *xVar );
            break;
    }
}

//----------------------------------------------------------------

void FD1DMulti::updateValue( int step ) const
{
    if( changeOfVar == LOG_FWDX )
        *xValue = stepForward[ step ] / stepForward[ step + 1 ] * *xValue;
}

//----------------------------------------------------------------

FDProductSP FD1DMulti::makeProduct(const IProdCreatorSP & creator)
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
            throw ModelException( "FD1DMulti::makeProduct",
                "IndexSpec factor is '" + spec->getFactor()->getName() + "', "
                "while '" + factors[0]->getName() + "' is expected" );
        }
    }

    return FDModel::makeProduct( creator );
}

//----------------------------------------------------------------

FDModel::IFDSolverSP FD1DMulti::createSolver()
{
    return IFDSolverSP( new FD1DMultiSolver( *this ) );
}

//----------------------------------------------------------------

double FD1DMulti::getPrice0( const TreeSlice & price ) const
{
    const string method = "FD1DMulti::getPrice0";
    try
    {
        const TreeSliceLayer * layer = dynamic_cast<const TreeSliceLayer * >( &price );
        if( layer )
            return layer->getPrice0();

        double s0 = stepForward[ 0 ];

		double * sVal = *xValue;
		TreeSliceEQ * sValEq = dynamic_cast< TreeSliceEQ * >( &(*xValue) );
		if ( sValEq ) {
			sVal = &((*sValEq)(iState*dim));
		}

		double * pVal = price.getValues();
		TreeSliceEQ * pValEq = const_cast<TreeSliceEQ *>(dynamic_cast<const TreeSliceEQ * >( &price ));
		if ( pValEq ) {
			pVal = &((*pValEq)(iState*dim));
		}

        double price0, deltaDummy, gammaDummy;
        int done = FDInterpolationD( dim, sVal, pVal, 1, &s0, &price0, &deltaDummy, &gammaDummy );

        return price0;
    }
    catch( exception & e )
    {
        throw ModelException( e, method );
    }
}

double FD1DMulti::getTruncationStd(){
	return truncation;
}

double FD1DMulti::getTrdYrFrac( int step ) const
{
    return timeLine->TradeYrFrac[step];
}

void FD1DMulti::load( CClassSP & type )
{
    REGISTER(FD1DMulti, type);
    SUPERCLASS(FDModel);

    FIELD(truncation, "num of stdev to truncate for spatial dimention");
    FIELD(dim, "num of point in FD grid for spatial dimension");
    FIELD(stepsPerYear, "num of time points per year");

	FIELD(nPdes, "num of PDEs");
	FIELD_MAKE_TRANSIENT(nPdes);

	FIELD(iState, "the index of the solution required");
	FIELD_MAKE_TRANSIENT(iState);

    FIELD_USING_ALIAS(changeOfVarStr, changeOfVar, "X: underlying, LOG_X: log(underlying), LOG_FWDX: log(fwd underlying)");
    FIELD_MAKE_OPTIONAL(changeOfVar);

    FIELD_USING_ALIAS(varStepSpacingStr, varStepSpacing, "NONE: fixed spacing, STD: variable spacing");
    FIELD_MAKE_OPTIONAL(varStepSpacing);

    FIELD(DEBUG_DumpToFile, "file name to dump slices to")
    FIELD_MAKE_OPTIONAL(DEBUG_DumpToFile);
}

CClassConstSP const FD1DMulti::TYPE =
    CClass::registerClassLoadMethod( "FD1DMulti", typeid( FD1DMulti ), load );

DRLIB_END_NAMESPACE
