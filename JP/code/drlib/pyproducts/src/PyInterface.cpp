#include "edginc/PyInterface.hpp"

DRLIB_BEGIN_NAMESPACE

BOOST_PYTHON_MODULE(qlib)
{
    boost::python::class_< DateTimePy >( "DateTime", boost::python::init< string, string >() )
        .def( "getDate", &DateTimePy::getDate )
        .def( "getTime", &DateTimePy::getTime )
        .def( "toString", &DateTimePy::toString )
        .def( "__lt__", &DateTimePy::__lt__ )
        .def( "__le__", &DateTimePy::__le__ )
        .def( "__eq__", &DateTimePy::__eq__ )
        .def( "__ne__", &DateTimePy::__ne__ )
        .def( "__gt__", &DateTimePy::__gt__ )
        .def( "__ge__", &DateTimePy::__ge__ )
    ;

    boost::python::class_< DateTimeArrayPy >( "DateTimeArray", boost::python::no_init )
        .def( "size", &DateTimeArrayPy::size )
        .def( "__getitem__", &DateTimeArrayPy::__getitem__ )
    ;

    boost::python::class_< DoubleArrayPy >( "DoubleArray", boost::python::no_init )
        .def( "assign", &DoubleArrayPy::assign )
        .def( "__getitem__", &DoubleArrayPy::__getitem__ )
        .def( "__setitem__", &DoubleArrayPy::__setitem__ )
    ;

    boost::python::class_< IMCPricesPy >( "IMCPrices", boost::python::no_init )
        .def( "add", &IMCPricesPy::add )
        .def( "maxWithZero", &IMCPricesPy::maxWithZero )
    ;

    boost::python::class_< IPathGeneratorPy >( "IPathGenerator", boost::python::no_init )
    ;

    boost::python::class_< IRefLevelIStateVarPy >( "IRefLevelIStateVar", boost::python::no_init )
        .def( "refLevel", &IRefLevelIStateVarPy::refLevel )
    ;

    boost::python::class_< SVGenSpotIStateVarPy >( "SVGenSpotIStateVar", boost::python::no_init )
        .def( "doingPast", &SVGenSpotIStateVarPy::doingPast )
        .def( "path", &SVGenSpotIStateVarPy::path )
    ;

    boost::python::class_< SVPathPy >( "SVPath", boost::python::no_init )
        .def( "begin", &SVPathPy::begin )
        .def( "end", &SVPathPy::end )
        .def( "__getitem__", &SVPathPy::__getitem__ )
    ;

    boost::python::enum_< FDProduct::UpdateType >( "UpdateType" )
        .value("BWD", FDProduct::BWD )
        .value("BWD_T", FDProduct::BWD_T )
        .value("FWD_0", FDProduct::FWD_0 )
        .value("FWD", FDProduct::FWD )
        .value("BWD_NODE_INSERTION", FDProduct::BWD_NODE_INSERTION )
    ;

    boost::python::class_< FlexDatesPy >( "FlexDates", boost::python::no_init )
        .def( "getDates", &FlexDatesPy::getDates )
    ;

    boost::python::class_< OptionSchedDatesPy >( "OptionSchedDates", boost::python::no_init )
        .def_readonly( "notifDate", &OptionSchedDatesPy::notifDate )
        .def_readonly( "exerciseDate", &OptionSchedDatesPy::exerciseDate )
    ;

    boost::python::def( "min", &SliceOperandPy::__smin );
    boost::python::def( "min", &SliceOperandPy::__sminL );
    boost::python::def( "min", &SliceOperandPy::__sminR );
    boost::python::def( "max", &SliceOperandPy::__smax );
    boost::python::def( "max", &SliceOperandPy::__smaxL );
    boost::python::def( "max", &SliceOperandPy::__smaxR );
    boost::python::def( "log", &SliceOperandPy::__log );
    boost::python::def( "exp", &SliceOperandPy::__exp );
    boost::python::def( "sqrt", &SliceOperandPy::__sqrt );
    boost::python::def( "cond", &SliceOperandPy::__cond );
    boost::python::def( "cond", &SliceOperandPy::__condL );
    boost::python::def( "cond", &SliceOperandPy::__condR );
    boost::python::def( "cond", &SliceOperandPy::__condB );
    boost::python::def( "cond", &SliceOperandPy::__condBL );
    boost::python::def( "cond", &SliceOperandPy::__condBR );

    boost::python::class_< SliceOperandPy >( "SliceOperandPy", boost::python::no_init )
        .def( "__add__", &SliceOperandPy::__add )
        .def( "__add__", &SliceOperandPy::__addL )
        .def( "__radd__", &SliceOperandPy::__addR )
        .def( "__sub__", &SliceOperandPy::__sub )
        .def( "__sub__", &SliceOperandPy::__subL )
        .def( "__rsub__", &SliceOperandPy::__subR )
        .def( "__mul__", &SliceOperandPy::__mul )
        .def( "__mul__", &SliceOperandPy::__mulL )
        .def( "__rmul__", &SliceOperandPy::__mulR )
        .def( "__div__", &SliceOperandPy::__div )
        .def( "__div__", &SliceOperandPy::__divL )
        .def( "__rdiv__", &SliceOperandPy::__divR )
        .def( "__pow__", &SliceOperandPy::__pow )
        .def( "__pow__", &SliceOperandPy::__powL )
        .def( "__pos__", &SliceOperandPy::__pos )
        .def( "__neg__", &SliceOperandPy::__neg )
        .def( "__abs__", &SliceOperandPy::__abs )
        .def( "__lt__", &SliceOperandPy::__lt )
        .def( "__lt__", &SliceOperandPy::__ltL )
        .def( "__rlt__", &SliceOperandPy::__ltR )
        .def( "__le__", &SliceOperandPy::__le )
        .def( "__le__", &SliceOperandPy::__leL )
        .def( "__rle__", &SliceOperandPy::__leR )
        .def( "__gt__", &SliceOperandPy::__gt )
        .def( "__gt__", &SliceOperandPy::__gtL )
        .def( "__rgt__", &SliceOperandPy::__gtR )
        .def( "__ge__", &SliceOperandPy::__ge )
        .def( "__ge__", &SliceOperandPy::__geL )
        .def( "__rge__", &SliceOperandPy::__geR )
        .def( "__or__", &SliceOperandPy::__or )
        .def( "__or__", &SliceOperandPy::__orL )
        .def( "__ror__", &SliceOperandPy::__orR )
        .def( "__and__", &SliceOperandPy::__and )
        .def( "__and__", &SliceOperandPy::__andL )
        .def( "__rand__", &SliceOperandPy::__andR )
        .def( "__not__", &SliceOperandPy::__not )
    ;

    boost::python::class_< TreeSliceConstPy, boost::python::bases< SliceOperandPy > >( "TreeSliceConst", boost::python::no_init )
    ;

    boost::python::class_< TreeSlicePy, boost::python::bases< SliceOperandPy > >( "TreeSlice", boost::python::no_init )
        .def( "assign", &TreeSlicePy::assignSliceMarker )
        .def( "assign", &TreeSlicePy::assignTreeSliceConst )
        .def( "assign", &TreeSlicePy::assignTreeSlice )
        .def( "assign", &TreeSlicePy::assignDouble )
        .def( "isZero", &TreeSlicePy::isZero )
    ;

    boost::python::class_< TreeSliceSPPy >( "TreeSliceSP" )
        .def( "getObject",  &TreeSliceSPPy::getObject )
    ;

    boost::python::class_< FDModelPy >( "FDModel", boost::python::no_init )
        .def( "getDate", &FDModelPy::getDate )
        .def( "getToday", &FDModelPy::getToday )
        .def( "getZero", &FDModelPy::getZero )
    ;

    boost::python::class_< FDProductPy >( "FDProduct", boost::python::no_init )
        .def( "getValue", &FDProductPy::getValue )
    ;

    boost::python::class_< FDProductArrayPy >( "FDProductArray", boost::python::no_init )
        .def( "__getitem__", &FDProductArrayPy::__getitem__ )
    ;

    boost::python::class_< YieldCurvePy >( "YieldCurvePy", boost::python::no_init )
        .def( "getName", &YieldCurvePy::getName )
    ;

    boost::python::class_< KPythonTreePy >( "KPythonTree", boost::python::no_init )
        .def( "startDEV", &KPythonTreePy::startDEV )
        .def_readonly( "sched", &KPythonTreePy::sched )
        .def_readonly( "strikes", &KPythonTreePy::strikes )
        .def_readonly( "discount", &KPythonTreePy::discount )
        .def_readonly( "model", &KPythonTreePy::model )
        .def_readonly( "undProd", &KPythonTreePy::undProd )
        .def_readonly( "keepValueIdx", &KPythonTreePy::keepValueIdx )
        .def_readonly( "optionPrice", &KPythonTreePy::optionPrice )
        .def_readonly( "optSkipExerPrice", &KPythonTreePy::optSkipExerPrice )
        .def_readonly( "undExer", &KPythonTreePy::undExer )
    ;

    boost::python::class_< MCPythonProdPy >( "MCPythonProd", boost::python::no_init )
        .def_readonly( "isBest", &MCPythonProdPy::isBest )
        .def_readonly( "strike", &MCPythonProdPy::strike )
        .def_readonly( "averageOutDates", &MCPythonProdPy::averageOutDates )
        .def_readonly( "notional", &MCPythonProdPy::notional )
        .def_readonly( "sumOut", &MCPythonProdPy::sumOut )
        .def_readonly( "sumSoFar", &MCPythonProdPy::sumSoFar )
        .def_readonly( "CoP", &MCPythonProdPy::CoP )
        .def_readonly( "nbAssets", &MCPythonProdPy::nbAssets )
        .def_readonly( "spotSV", &MCPythonProdPy::spotSV )
        .def_readonly( "refLevelSV", &MCPythonProdPy::refLevelSV )
    ;
}

void pyInterfaceInit()
{
    static bool done = false;

    if( ! done )
    {
        try
        {
            Py_Initialize();
            initqlib();
            done = true;
        }
        catch( boost::python::error_already_set )
        {
            PyErr_Print();
        }
    }
}

DRLIB_END_NAMESPACE
