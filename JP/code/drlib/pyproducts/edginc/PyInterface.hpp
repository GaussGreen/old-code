#ifndef QLIB_PY_INTERFACE_HPP
#define QLIB_PY_INTERFACE_HPP

#include "edginc/KPython.hpp"
#include "edginc/MCPython.hpp"
#include "edginc/PyTreeSlice.hpp"

DRLIB_BEGIN_NAMESPACE

class DateTimePy
{
public:
    // We need a local copy of DateTime because most QLib methods return it by value
    DateTime obj;

    DateTimePy( const string & date, const string & time ) : obj( date, time ) {}
    DateTimePy( const DateTime & obj ) : obj( obj ) {}

    int getDate() { return obj.getDate(); }
    int getTime() { return obj.getTime(); }
    string toString() { return obj.toString(); }

    bool __lt__( const DateTimePy & other) { return obj < other.obj; }
    bool __le__( const DateTimePy & other) { return obj <= other.obj; }
    bool __eq__( const DateTimePy & other) { return obj == other.obj; }
    bool __ne__( const DateTimePy & other) { return obj != other.obj; }
    bool __gt__( const DateTimePy & other) { return obj > other.obj; }
    bool __ge__( const DateTimePy & other) { return obj >= other.obj; }
};

class DateTimeArrayPy
{
public:
    const DateTimeArray & obj;

    DateTimeArrayPy( const DateTimeArray & obj ) : obj( obj ) {}

    int size() { return obj.size(); }

    DateTimePy __getitem__( int index ) { return DateTimePy( obj[index] ); }
};

class DoubleArrayPy
{
public:
    DoubleArray & obj;

    DoubleArrayPy( DoubleArray & obj ) : obj( obj ) {}

    void assign( const DoubleArrayPy & src ) { obj = src.obj; }

    double __getitem__( int index ) { return obj[index]; }
    void __setitem__( int index, double value ) { obj[index] = value; }
};

class IMCPricesPy
{
public:
    IMCPrices & obj;

    IMCPricesPy( IMCPrices & obj ) : obj( obj ) {}

    double maxWithZero( double value ) { return obj.maxWithZero( value ); }
    void add( double value ) { obj.add( value ); }
};

class IPathGeneratorPy
{
public:
    const IPathGenerator & obj;

    IPathGeneratorPy( const IPathGenerator & obj ) : obj( obj ) {}
};

class IRefLevelIStateVarPy
{
public:
    IRefLevel::IStateVar & obj;

    IRefLevelIStateVarPy( IRefLevel::IStateVar & obj ) : obj( obj ) {}

    double refLevel( int iAsset ) { return obj.refLevel( iAsset ); }
};

class SVPathPy
{
public:
    const SVPath & obj;

    SVPathPy( const SVPath & obj ) : obj( obj ) {}

    int begin() { return obj.begin(); }
    int end() { return obj.end(); }

    double __getitem__( int index ) { return obj[index]; }
};

class SVGenSpotIStateVarPy
{
public:
    const SVGenSpot::IStateVar & obj;

    SVGenSpotIStateVarPy( const SVGenSpot::IStateVar & obj ) : obj( obj ) {}

    bool doingPast() { return obj.doingPast(); }
    SVPathPy path( int index ) { return obj.path( index ); }
};

class FlexDatesPy
{
public:
    const FlexDates & obj;

    FlexDatesPy( const FlexDates & obj ) : obj( obj ) {}

    DateTimeArrayPy getDates() { return obj.getDates(); }
};

class OptionSchedDatesPy
{
public:
    const OptionSchedDates & obj;

    OptionSchedDatesPy( const OptionSchedDates & obj ) : obj( obj ), 
        notifDate( obj.notifDate ),
        exerciseDate( obj.exerciseDate ) {}

    FlexDatesPy notifDate;
    FlexDatesPy exerciseDate;
};

class FDModelPy
{
public:
    FDModel & obj;

    FDModelPy( FDModel & obj ) : obj( obj ) {}

    DateTimePy getDate( int step ) { return DateTimePy( obj.getDate( step ) ); }
    DateTimePy getToday() { return DateTimePy( obj.getToday() ); }

    void getZero(
        const DateTimePy & useDate,
        const DateTimePy & matDate,
        const string & curveName,
        TreeSliceSPPy & slice )
    {
        obj.getZero( useDate.obj, matDate.obj, curveName, slice.obj ); 
    }
};

class FDProductPy
{
public:
    const FDProduct & obj;

    FDProductPy( const FDProduct & obj ) : obj( obj ) {}
    
    TreeSliceConstPy getValue( int step, const DateTimePy & eventDate )
    {
        return obj.getValue( step, eventDate.obj ); 
    }
};

class FDProductArrayPy
{
public:
    const FDProductArray & obj;

    FDProductArrayPy( const FDProductArray & obj ) : obj( obj ) {}

    FDProductPy __getitem__( int index ) { return *obj[index]; }
};

class YieldCurvePy
{
public:
    const YieldCurveWrapper & obj;

    YieldCurvePy( const YieldCurveWrapper & obj ) : obj( obj ) {}

    string getName() { return obj.getName(); }
};

class KPythonTreePy
{
public:
    KPythonTree & obj;

    KPythonTreePy( KPythonTree & obj ) :
        obj( obj ), 
        sched( *obj.inst->sched ),
        strikes( *obj.inst->strikes ),
        discount( obj.inst->discount ), 
        model( *obj.model ),
        undProd( obj.undProd ),
        keepValueIdx( obj.keepValueIdx ), 
        optionPrice( obj.optionPrice ),
        optSkipExerPrice( obj.optSkipExerPrice ),
        undExer( obj.undExer ) {}

    bool startDEV( const TreeSlicePy & slice ) { return obj.startDEV( slice.obj ); }

    OptionSchedDatesPy sched;
    DoubleArrayPy strikes;
    YieldCurvePy discount;
    FDModelPy model;
    FDProductArrayPy undProd;
    int keepValueIdx;
    TreeSlicePy optionPrice;
    TreeSlicePy optSkipExerPrice;
    TreeSlicePy undExer;
};

class MCPythonProdPy
{
public:
    MCPythonProd & obj;

    MCPythonProdPy( MCPythonProd & obj ) :
        obj( obj ),
        isBest( obj.inst->isBest ), 
        strike( obj.inst->strike ),
        averageOutDates( obj.inst->averageOutDates ), 
        notional( obj.inst->notional ),
        sumOut( obj.sumOut ),
        sumSoFar( obj.sumSoFar ), 
        CoP( obj.CoP ),
        nbAssets( obj.nbAssets ),
        spotSV( *obj.spotSV ),
        refLevelSV( *obj.refLevelSV ) {}

    bool isBest;
    double strike;
    DateTimeArrayPy averageOutDates;
    double notional;
    DoubleArrayPy sumOut;
    DoubleArrayPy sumSoFar;
    double CoP;
    int nbAssets;
    SVGenSpotIStateVarPy spotSV;
    IRefLevelIStateVarPy refLevelSV;
};

void pyInterfaceInit();

DRLIB_END_NAMESPACE

#include <boost/python.hpp>

#endif
