//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRUtils.cpp
//
//   Description : Holds utility 'precompiled' FR functions
//
//   Author      : Mark A Robson
//
//   Date        : 1st August 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FRFunction.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/TabulatedFunc.hpp"
#include "edginc/Algorithm.hpp"
#include "edginc/Format.hpp"
#include "edginc/Black.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE
#define FR_GET_VALUE(x) ((x)->func(x))
#define FR_GET_SIZE(x) ((x)->size(x))

/** Holds utility 'precompiled' FR functions. Please note that functions
    don't have to live here - it's just a convenient dumping ground */
class FRUtils{
    /** test function: takes args: (double p1, int p2, bool p3) */
    static double testFunc1(FRFunction::Val* args){
        // args 0 requested to come through as expression
        double d = FR_GET_VALUE(args[0].dExp);
        // get hold of magic extra parameter which is the current index
        int index = args[3].i;
        // args 1 requested to come through as vector of expressions
        int i =  FR_GET_VALUE((*args[1].rValues)[index].iExp);
        return (d + (double)(args[2].b? i:-i));
    }

    /** sums up double array */
    static double sum(FRFunction::Val* args){
        FRIfaces::IRValueDoubleArray::RT* rt = args[0].dExpArray; // for ease
        int size = FR_GET_SIZE(rt);
        FRIfaces::IRValueDoubleArray::TGetValue* arrayGetValue = rt->func;
        double sum = 0.0;
        for (int i = 0; i < size; i++){
            sum += arrayGetValue(rt, i);
        }
        return sum;
    }

    /** Returns number of element in a bool array that are true */
    static int numTrue(FRFunction::Val* args){
        FRIfaces::IRValueBoolArray::RT* rt = args[0].bExpArray; // for ease
        int size = FR_GET_SIZE(rt);
        FRIfaces::IRValueBoolArray::TGetValue* arrayGetValue = rt->func;
        int num = 0;
        for (int i = 0; i < size; i++){
            if (arrayGetValue(rt, i)){
                num++;
            }
        }
        return num;
    }

    /** does product of double array */
    static double product(FRFunction::Val* args){
        FRIfaces::IRValueDoubleArray::RT* rt = args[0].dExpArray; // for ease
        int size = FR_GET_SIZE(rt);
        FRIfaces::IRValueDoubleArray::TGetValue* arrayGetValue = rt->func;
        double total = 1.0;
        for (int i = 0; i < size; i++){
            total *= arrayGetValue(rt, i);
        }
        return total;
    }

    /** sums up weighted double array */
    static double sumProduct(FRFunction::Val* args){
        FRIfaces::IRValueDoubleArray::RT* rtPerf = args[0].dExpArray;
        FRIfaces::IRValueDoubleArray::RT* rtWeights = args[1].dExpArray;
        int size = FR_GET_SIZE(rtPerf);
        if (size != FR_GET_SIZE(rtWeights)){
            throw ModelException("FRUtils::sumProduct", "Both arrays must "
                                 "be the same length");
        } 
        FRIfaces::IRValueDoubleArray::TGetValue* perfGetValue = rtPerf->func;
        FRIfaces::IRValueDoubleArray::TGetValue* weightsGetValue = 
            rtWeights->func;
        double perf = 0.0;
        for (int i = 0; i < size; i++){
            perf += perfGetValue(rtPerf, i) * weightsGetValue(rtWeights, i);
        }
        return perf;
    }

    /** sums up sorted weighted double array */
    static double rainbow(FRFunction::Val* args){
        FRIfaces::IRValueDoubleArray::RT* rtPerf = args[0].dExpArray;
        FRIfaces::IRValueDoubleArray::RT* rtWeights = args[1].dExpArray;
        int size = FR_GET_SIZE(rtPerf);
        if (size != FR_GET_SIZE(rtWeights)){
            throw ModelException("FRUtils::rainbow", "Both arrays must "
                                 "be the same length");
        } 
        FRIfaces::IRValueDoubleArray::TGetValue* perfGetValue = rtPerf->func;
        FRIfaces::IRValueDoubleArray::TGetValue* weightsGetValue = 
            rtWeights->func;
        // note: taking a copy - could be slow
        DoubleArray dbArray;
        dbArray.reserve(size);
        for (int i = 0; i < size; i++){
            dbArray.push_back(perfGetValue(rtPerf, i));
        }
        Algorithm::shellSort(dbArray);
        double perf = 0.0;
        for (int j = 0; j < size; j++){
            perf += dbArray[j] * weightsGetValue(rtWeights, j);
        }
        return perf;
    }

    static double interpSchedule(FRFunction::Val* args){
        const Schedule*       schedule = args[0].sched;
        const DateTime::Date* date = args[1].dt;
        return schedule->interpolate(*date);
    }

    static double interpTabulatedFunc(FRFunction::Val* args){
        const TabulatedFunc*  tabFunc = args[0].tabFunc;
        double x = args[1].d;
        return tabFunc->interpolate(x);
    }

    //// calculates ln(x)
    static double naturalLog(FRFunction::Val* args){
        const static string method("FRUtils::naturalLog");
        double value = args[0].d;
        
        if (value > 0.0) {
            return log(value);
        }
        else {
            throw ModelException(method, "Can not calculate the logarithm of the negative or nul number: "+
                                 Format::toString(value));
        }
    }

    /** log of a double array */
    static void naturalLogArray(FRFunction::Val* args, DoubleArray& dbles){
        const static string method("FRUtils::naturalLogArray");
        FRIfaces::IRValueDoubleArray::RT* rt = args[0].dExpArray; // for ease
        FRIfaces::IRValueDoubleArray::TGetValue* arrayGetValue = rt->func;
        
        // size of the input array
        int size = FR_GET_SIZE(rt);
        // resize the result array
        dbles.resize(size);
        
        for (int i = 0; i < size; i++){
            double value = arrayGetValue(rt, i);
            if (value > 0.0) {
                dbles[i] = log(value);
            }
            else {
                throw ModelException(method, "Can not calculate the logarithm of the negative or nul number: "+
                                     Format::toString(value));
            }
        }
    }

    //// calculates floor(x)
    static int flexFloor(FRFunction::Val* args){
        const static string method("FRUtils::flexFloor");
        double value = args[0].d;
        
        return (int)floor(value);
    }

    /** floor of a double array */
    static void floorArray(FRFunction::Val* args, IntArray& ints){
        const static string method("FRUtils::floor");
        FRIfaces::IRValueDoubleArray::RT* rt = args[0].dExpArray; // for ease
        FRIfaces::IRValueDoubleArray::TGetValue* arrayGetValue = rt->func;
        
        // size of the input array
        int size = FR_GET_SIZE(rt);
        // resize the result array
        ints.resize(size);
        
        for (int i = 0; i < size; i++){
            double value = arrayGetValue(rt, i);
            ints[i] = (int)floor(value);
        }
    }

    //// calculates exp(x)
    static double exponential(FRFunction::Val* args){
        return exp(args[0].d);
    }

    /** sums up double array */
    static double black(FRFunction::Val* args){
        enum {
            eIsCall = 0,
            eFwd, 
            eStrike,
            eTerm,
            eVol
        };
        double vol = args[eVol].d;
        double var = vol * vol * args[eTerm].d;
        return Black::price(args[eIsCall].b,args[eFwd].d, args[eStrike].d,
                            1.0, // no discounting
                            var);
    }

    //// don't build instances of this class
    FRUtils();

    /** Returns a double array which is the average value of each of
        the components in the supplied array across a specified set of
        time points */
    static void average(FRFunction::Val* args, DoubleArray& dbles){
        const static string method("FRUtils::average");
        enum {
            eValues = 0,
            eStartIndex, 
            eLength
        };
        const vector<FRIfaces::RValUnion>& rValues = *(args[eValues].rValues);
        int start = args[eStartIndex].i;
        int length = args[eLength].i;
        int end = start + length;
        if (start < 0 || end > (int) rValues.size()){
            throw ModelException(method, "Specified interval ["+
                                 Format::toString(start)+", "+
                                 Format::toString(end)+") must lie within "
                                 "simulation");
        }
        if (end <= start){
            throw ModelException(method,
                                 "Number of values to average <= 0!");
        }
        // get hold of magic extra parameter which is the current index
        int index = args[3].i;
        if (index < end-1) {
            throw ModelException(method,
                                 "Cannot reference future sample at index " + Format::toString(end-1));
        }
        FRIfaces::IRValueDoubleArray::RT* rt = rValues[start].dArrayExp;
        int numValues = FR_GET_SIZE(rt);
        dbles.resize(numValues);
        // handle first date explicitly
        FRIfaces::IRValueDoubleArray::TGetValue* arrayGetValue = rt->func;
        for (int j = 0; j < numValues; j++){
            dbles[j] = arrayGetValue(rt, j)/length;
        }
        // thenloop across remaining simulation dates
        for (int i = start+1; i < end; i++){
            FRIfaces::IRValueDoubleArray::RT* rt = rValues[i].dArrayExp;
            // validate length of array
            if (FR_GET_SIZE(rt) != numValues){
                throw ModelException(method, "Array to average must have "
                                     "the same number of elements at each "
                                     "averaging date");
            }
            // loop across values
            for (int j = 0; j < numValues; j++){
                arrayGetValue = rt->func;
                dbles[j] += arrayGetValue(rt, j)/length;
            }
        }
    }

    /** For sort algorithm - allows an int array to be sorted based upon
        the value of a corresponding double array */
    class DoubleSortHelper{
        vector<double>::iterator dbles;
    public:
        DoubleSortHelper(vector<double>::iterator dbles): dbles(dbles){}
        //// we want best first (ie largest to smallest)
        bool operator()(int i1, int i2){
            return (dbles[i1] > dbles[i2]);
        }
    };

    /* Order is a function which takes an array of doubles and
       the returns the index of each of the components from best to
       worst. (return type is IntArray) */
    static void order(FRFunction::Val* args, IntArray& ints){
        const static string method("FRUtils::order");
        enum {
            eValues = 0,
        };
        // start by finding length of array
        FRIfaces::IRValueDoubleArray::RT* rt = args[eValues].dExpArray;
        int numValues = FR_GET_SIZE(rt);
        ints.resize(numValues);
        /* Then cache value of performances (better than using getValue at
           every operation in the sort) */
        vector<double>  dbles;
        dbles.reserve(numValues);
        // we have to create an initial array to point to each element before
        // we can sort
        FRIfaces::IRValueDoubleArray::TGetValue* rtFunc = rt->func;
        for (int i = 0; i < numValues; i++){
            ints[i] = i;
            dbles.push_back(rtFunc(rt, i));
        }
        // then create class which holds access to performances
        DoubleSortHelper sortHelper(dbles.begin());
        // finally sort
        sort(ints.begin(), ints.end(), sortHelper);
    }

    /* subOrder is a function which takes a subset of the supplied
       array of doubles and returns the index of each of the
       components from best to worst. If the second parameter is
       lastOrder then the subset chosen to order is {perf[lastOrder[i]],
       perf[lastOrder[i+1]], ..., perf[lastOrder[i+length-1]} 
       (return type is IntArray) */
    static void subOrder(FRFunction::Val* args, IntArray& ints){
        const static string method("FRUtils::subOrder");
        enum {
            eValues = 0,
            eSubSet,
            eStart,
            eLength
        };
        // start by finding/checking length of arrays
        FRIfaces::IRValueDoubleArray::RT* rtValues = args[eValues].dExpArray;
        int numValues = FR_GET_SIZE(rtValues);
        FRIfaces::IRValueIntArray::RT* rtSubSet = args[eSubSet].iExpArray;
        int numSubSet = FR_GET_SIZE(rtSubSet);
        int indexStart = args[eStart].i;
        int length = args[eLength].i;
        if (indexStart < 0 || length < 0 || indexStart+length > numSubSet){
            throw ModelException(method, "Asked to order "+
                                 Format::toString(length)+" items from array"
                                 " using an array of length "+
                                 Format::toString(numSubSet)+ " starting at "
                                 "index "+Format::toString(indexStart));
        }
        ints.resize(length);
        /* Then cache value of performances (better than using getValue at
           every operation in the sort) */
        vector<double>  dbles;
        dbles.reserve(length);
        // Also cache the subset (see final step)
        vector<int>  subSet;
        subSet.reserve(length);
        // we have to create an initial array to point to each element before
        // we can sort
        FRIfaces::IRValueDoubleArray::TGetValue* rtValuesFunc = rtValues->func;
        FRIfaces::IRValueIntArray::TGetValue* rtSubSetFunc = rtSubSet->func;
        for (int i = 0; i < length; i++){
            int index = rtSubSetFunc(rtSubSet, indexStart+i);
            if (index < 0 || index >= numValues){
                throw ModelException(method, "Asked to include index "+
                                     Format::toString(index)+" of array of"
                                     " length "+Format::toString(numValues)+
                                     " when ordering");
            }
            subSet.push_back(index); // cache index
            dbles.push_back(rtValuesFunc(rtValues, index));
            ints[i] = i;
        }
        // then create class which holds access to performances
        DoubleSortHelper sortHelper(dbles.begin());
        // sort our subset
        sort(ints.begin(), ints.end(), sortHelper);
        // and then map back to true id's
        for (int j = 0; j < length; j++){
            ints[j] = subSet[ints[j]];
        }
    }

    /* select is a function which here returns {dbles[subset[start]],
       dbles[subset[1]], ..., dbles[subset[start+l-1]]}.  (return type
       is DoubleArray) */
    static void select(FRFunction::Val* args, DoubleArray& dbles){
        const static string method("FRUtils::select");
        enum {
            eValues = 0,
            eSubSet,
            eStart,
            eLength
        };
        // start by finding/checking length of arrays
        FRIfaces::IRValueDoubleArray::RT* rtValues = args[eValues].dExpArray;
        int numValues = FR_GET_SIZE(rtValues);
        FRIfaces::IRValueIntArray::RT* rtSubSet = args[eSubSet].iExpArray;
        int numSubSet = FR_GET_SIZE(rtSubSet);
        int indexStart = args[eStart].i;
        int length = args[eLength].i;
        if (indexStart < 0 || length < 0 || indexStart+length > numSubSet){
            throw ModelException(method, "Asked to select "+
                                 Format::toString(length)+" items from array"
                                 " using an array of length "+
                                 Format::toString(numSubSet)+ " starting at "
                                 "index "+Format::toString(indexStart));
        }
        dbles.resize(length);
        FRIfaces::IRValueDoubleArray::TGetValue* rtValuesFunc = rtValues->func;
        FRIfaces::IRValueIntArray::TGetValue* rtSubSetFunc = rtSubSet->func;
        for (int i = 0; i < length; i++){
            int index = rtSubSetFunc(rtSubSet, indexStart+i);
            if (index < 0 || index >= numValues){
                throw ModelException(method, "Asked to include index "+
                                     Format::toString(index)+" of array of"
                                     " length "+Format::toString(numValues));
            }
            dbles[i] = rtValuesFunc(rtValues, index);
        }
    }
    class Slice{
    public:
        enum {
            eValues = 0,
            eStart,
            eLength,
            eCurrentIndex
        };
        static void checkTimeLimits(int subStart, int subLength, int length,
                                    int currentIndex){
            static const string method("Slice::checkTimeLimits");
            int subEnd = subStart + subLength;
            if (subStart < 0 || subEnd > length){
                throw ModelException(method, "Specified interval ["+
                                     Format::toString(subStart)+", "+
                                     Format::toString(subEnd)+") must lie "
                                     "within simulation");
            }
            if (subLength < 0){
                throw ModelException(method,
                                     "Number of time points < 0!");
            }
            if (subEnd-1 > currentIndex) {
                throw ModelException(method,
                                     "Slice ends at index " + 
                                     Format::toString(subEnd-1) +
                                     " which is after current index " + 
                                     Format::toString(currentIndex));
            }
        }
        /** Capture a slice through time of a double variable returned as an
            double array */
        static void dSlice(FRFunction::Val* args, DoubleArray& dbles){
            const static string method("FRUtils::dSlice");
            const vector<FRIfaces::RValUnion>& rValues = 
                *(args[eValues].rValues);
            int start = args[eStart].i;
            int length = args[eLength].i;
            checkTimeLimits(start, length, rValues.size(),
                            args[eCurrentIndex].i);
            dbles.resize(length);
            for (int i = 0; i < length; i++){
                FRIfaces::IRValueDouble::RT* rt = rValues[start+i].dExp;
                dbles[i] = FR_GET_VALUE(rt);
            }
        }
        /** Capture a slice through time of a int variable returned as an
            int array */
        static void iSlice(FRFunction::Val* args, IntArray& dbles){
            const static string method("FRUtils::iSlice");
            const vector<FRIfaces::RValUnion>& rValues = 
                *(args[eValues].rValues);
            int start = args[eStart].i;
            int length = args[eLength].i;
            checkTimeLimits(start, length, rValues.size(),
                            args[eCurrentIndex].i);
            dbles.resize(length);
            for (int i = 0; i < length; i++){
                FRIfaces::IRValueInt::RT* rt = rValues[start+i].iExp;
                dbles[i] = FR_GET_VALUE(rt);
            }
        }
        /** Capture a slice through time of a bool variable returned as an
            bool array */
        static void bSlice(FRFunction::Val* args, BoolArray& dbles){
            const static string method("FRUtils::bSlice");
            const vector<FRIfaces::RValUnion>& rValues = 
                *(args[eValues].rValues);
            int start = args[eStart].i;
            int length = args[eLength].i;
            checkTimeLimits(start, length, rValues.size(),
                            args[eCurrentIndex].i);
            dbles.resize(length);
            for (int i = 0; i < length; i++){
                FRIfaces::IRValueBool::RT* rt = rValues[start+i].bExp;
                dbles[i] = FR_GET_VALUE(rt);
            }
        }
    };
public:
    static void registerFuncs(){
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleType, FRIfaces::intType, FRIfaces::boolType};
            static const char* paramNames[] = {"p1", "p2", "p3"};
            static const FRFunction::ValCalcType calcTypes[] = 
            {FRFunction::expression, FRFunction::arrayExpression, 
             FRFunction::native};
            FRFunction::registerGenericDoubleFunction(
                "__test1", 
                0, // description - null means the function is hidden from help
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                calcTypes, // mixture
                testFunc1);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleArrayType
            };
            static const char* paramNames[] = {"dbles"};
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::expression};
            FRFunction::registerGenericDoubleFunction(
                "sum", "Sums an array of doubles",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                calcTypes,
                sum);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {
                FRIfaces::boolArrayType
            };
            static const char* paramNames[] = {"bools"};
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::expression};
            FRFunction::registerGenericIntFunction(
                "numTrue", 
                "Counts the number of bools in an array that are true",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                calcTypes,
                numTrue);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleArrayType
            };
            static const char* paramNames[] = {"dbles"};
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::expression};
            FRFunction::registerGenericDoubleFunction(
                "product", "Calculates product of an array of doubles",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                calcTypes,
                product);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleArrayType, FRIfaces::doubleArrayType
            };
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::expression, FRFunction::expression};
            static const char* paramNames[] = {"dbles", "weights"};
            FRFunction::registerGenericDoubleFunction(
                "sumProduct", "Sums a weighted array of doubles",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                calcTypes, 
                sumProduct);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleArrayType, FRIfaces::doubleArrayType
            };
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::expression, FRFunction::expression};
            static const char* paramNames[] = {"perfs", "weights"};
            FRFunction::registerGenericDoubleFunction(
                "rainbow", "Sums a sorted weighted array of doubles",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                calcTypes,
                rainbow);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {
                FRIfaces::scheduleType, FRIfaces::dateType
            };
            static const char* paramNames[] = {"schedule", "dateForInterp"};
            FRFunction::registerGenericDoubleFunction(
                "interpSchedule", "Interpolates a schedule",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                0, // all params as primitive
                interpSchedule);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {
                FRIfaces::tabulatedFuncType, FRIfaces::doubleType
            };
            static const char* paramNames[] = {"tabulatedFunc", 
                                               "whereToInterp"};
            FRFunction::registerGenericDoubleFunction(
                "interpTabulatedFunc", "Evaluates a tabulated function using "
                "interpolation",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                0, // all params as primitive
                interpTabulatedFunc);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {FRIfaces::doubleType};
            static const char* paramNames[] = {"x"};
            FRFunction::registerGenericDoubleFunction(
                "log", "Calculates natural logarithm of x",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                0, // all params as primitive
                naturalLog);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleArrayType};
            static const char* paramNames[] = {"doubleArray"};
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::expression};
            FRFunction::Func func;
            func.dArrayFunc = naturalLogArray;
            FRFunction::registerGenericFunction(
                "logArray", "Calculates natural logarithm of each element of the supplied array",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                calcTypes,
                FRIfaces::doubleArrayType,
                func);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {FRIfaces::doubleType};
            static const char* paramNames[] = {"x"};
            FRFunction::registerGenericIntFunction(
                "floor", "Calculates the floor of x",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                0, // all params as primitive
                flexFloor);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleArrayType};
            static const char* paramNames[] = {"doubleArray"};
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::expression};
            FRFunction::Func func;
            func.iArrayFunc = floorArray;
            FRFunction::registerGenericFunction(
                "floorArray", "Calculates the floor of each element of the supplied array",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                calcTypes,
                FRIfaces::intArrayType,
                func);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {FRIfaces::doubleType};
            static const char* paramNames[] = {"x"};
            FRFunction::registerGenericDoubleFunction(
                "exp", "Calculates exponential of x",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                0, // all params as primitive
                exponential);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {
                FRIfaces::boolType,
                FRIfaces::doubleType,
                FRIfaces::doubleType,
                FRIfaces::doubleType,
                FRIfaces::doubleType};
            static const char* paramNames[] = {"isCall",
            "forward", "strike", "term", "vol"};
            FRFunction::registerGenericDoubleFunction(
                "black", "Returns the non discounted fair value of the "
                "specified european option",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                0, // all params as primitive
                black);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleArrayType, 
                FRIfaces::intType, 
                FRIfaces::intType};
            static const char* paramNames[] = {"valuesArray", "startIndex",
                                               "length"};
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::arrayExpression, FRFunction::native,
                FRFunction::native};
            FRFunction::Func func;
            func.dArrayFunc = average;
            FRFunction::registerGenericFunction(
                "average", "Calculates the average of each component in"
                " an array across specified interval",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                calcTypes,
                FRIfaces::doubleArrayType,
                func);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleArrayType};
            static const char* paramNames[] = {"valuesToOrder"};
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::expression};
            FRFunction::Func func;
            func.iArrayFunc = order;
            FRFunction::registerGenericFunction(
                "order", "Returns the index of each of the components in "
                "the supplied array from best to worst.",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                calcTypes,
                FRIfaces::intArrayType,
                func);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleArrayType,
                FRIfaces::intArrayType,
                FRIfaces::intType, 
                FRIfaces::intType};
            static const char* paramNames[] = {"valuesToOrder",
            "subSetToUse", "start", "length"};
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::expression, FRFunction::expression,
                FRFunction::native, FRFunction::native};
            FRFunction::Func func;
            func.iArrayFunc = subOrder;
            FRFunction::registerGenericFunction(
                "subOrder", "Returns the index of each of the components in "
                "a subset of the supplied array from best to worst. The "
                "subset used is defined via the last 3 parameters",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                calcTypes,
                FRIfaces::intArrayType,
                func);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleArrayType,
                FRIfaces::intArrayType,
                FRIfaces::intType, 
                FRIfaces::intType};
            static const char* paramNames[] = {"valuesToSelect",
            "subSetToUse", "start", "length"};
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::expression, FRFunction::expression,
                FRFunction::native, FRFunction::native};
            FRFunction::Func func;
            func.dArrayFunc = select;
            FRFunction::registerGenericFunction(
                "select", "Returns a subset of the supplied array of doubles. "
                "The subset used is defined via the last 3 parameters",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                calcTypes,
                FRIfaces::doubleArrayType,
                func);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {
                FRIfaces::doubleType, 
                FRIfaces::intType, 
                FRIfaces::intType};
            static const char* paramNames[] = {"doubleVariable", "startIndex",
                                               "length"};
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::arrayExpression, FRFunction::native,
                FRFunction::native};
            FRFunction::Func func;
            func.dArrayFunc = Slice::dSlice;
            FRFunction::registerGenericFunction(
                "dSlice", "Captures the values of a double variable across"
                " simulation dates as a double array",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                calcTypes,
                FRIfaces::doubleArrayType,
                func);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {
                FRIfaces::intType, 
                FRIfaces::intType, 
                FRIfaces::intType};
            static const char* paramNames[] = {"intVariable", "startIndex",
                                               "length"};
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::arrayExpression, FRFunction::native,
                FRFunction::native};
            FRFunction::Func func;
            func.iArrayFunc = Slice::iSlice;
            FRFunction::registerGenericFunction(
                "iSlice", "Captures the values of a int variable across"
                " simulation dates as a int array",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                calcTypes,
                FRIfaces::intArrayType,
                func);
        }
        {
            // NB Must be static
            static const FRIfaces::VarType types[] = {
                FRIfaces::boolType, 
                FRIfaces::intType, 
                FRIfaces::intType};
            static const char* paramNames[] = {"boolVariable", "startIndex",
                                               "length"};
            static const FRFunction::ValCalcType calcTypes[] = {
                FRFunction::arrayExpression, FRFunction::native,
                FRFunction::native};
            FRFunction::Func func;
            func.bArrayFunc = Slice::bSlice;
            FRFunction::registerGenericFunction(
                "bSlice", "Captures the values of a bool variable across"
                " simulation dates as a bool array",
                sizeof(types)/sizeof(FRIfaces::VarType),
                types,
                paramNames,
                calcTypes,
                FRIfaces::boolArrayType,
                func);
        }
    }
    
};

bool loadFRUtils() {
    FRUtils::registerFuncs();
    return true;
}


DRLIB_END_NAMESPACE
