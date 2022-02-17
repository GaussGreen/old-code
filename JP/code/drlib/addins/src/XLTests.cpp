//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLTests.hpp
//
//   Description : Tests for XL Kit
//
//   Author      : Mark A Robson
//
//   Date        : 20 Feb 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Maths.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/Asset.hpp"

DRLIB_BEGIN_NAMESPACE

class XLMultiply: public CObject{
    double   val1;
    double   val2;
public:
    static CClassConstSP const TYPE;

    double multiply(){
        return (val1 * val2);
    }

    XLMultiply(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLMultiply, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLMultiply);
        FIELD(val1, "first number");
        FIELD(val2, "second number");
        Addin::registerDoubleMethod("MULTIPLY",
                                    Addin::XL_TESTS,
                                    "Multiplies two numbers together",
                                    &XLMultiply::multiply);
    }

    static IObject* defaultXLMultiply(){
        return new XLMultiply();
    }
};

CClassConstSP const XLMultiply::TYPE = CClass::registerClassLoadMethod(
    "XLMultiply", typeid(XLMultiply), load);


class XLTestDouble1: public CObject{
public:
    static CClassConstSP const TYPE;
    CDoubleArray array1;
    double       val1;
    double       val2;
    CDoubleArray array2;

    static IObjectSP run(XLTestDouble1 *params){
        /* return original object - but add 1 to values */
        if (params->array1.size() != params->array2.size()){
            throw ModelException("XLTestDouble1::run", "Arrays are of "
                                 "different length");
        }
        XLTestDouble1  *addin = copy(params);
        int idx;
        for (idx = 0; idx < addin->array1.size(); idx++) {
            addin->array1[idx]++;
            addin->array2[idx]++;
        }
        addin->val1++;
        addin->val2++;
        return IObjectSP(addin);
    }
    

    XLTestDouble1(): CObject(TYPE), val1(-7.1), val2(0){
        array1 = CDoubleArray(3);
        array1[0] = 3.4;
        array1[1] = -2.1;
        array1[2] = 0.0;
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLTestDouble1, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLTestDouble1);
        FIELD(array1, "First array of doubles");
        FIELD_MAKE_OPTIONAL(array1);
        FIELD(val1, "Another double");
        FIELD_MAKE_OPTIONAL(val1);
        FIELD(val2, "Yet another double");
        FIELD(array2, "Array of doubles of equal length to first");
        Addin::registerClassObjectMethod("TEST_DOUBLE1",
                                         Addin::XL_TESTS,
                                         "Adds one to each number",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)run);
    }
    static IObject* defaultXLTestDouble1(){
        return new XLTestDouble1();
    }
};

CClassConstSP const XLTestDouble1::TYPE = CClass::registerClassLoadMethod(
    "XLTestDouble1", typeid(XLTestDouble1), load);

class XLTestDouble2: public CObject{
public:
    static CClassConstSP const TYPE;
    double       count;

    static IObjectSP run(XLTestDouble2 *addin){
        int             count    = Maths::max(0, (int) addin->count);
        IObjectSP       outObj;
        if (count < 10) {
            CDoubleArraySP  dbList(new CDoubleArray(count));
            for (int idx = 0; idx < count; idx++) {
                (*dbList)[idx] = idx;
            }
            outObj = dbList;
        } else {
            XLTestDouble1 *testDb1 = new XLTestDouble1;
            count -= 10;
            outObj = IObjectSP(testDb1);
            if (count > 0) {
                testDb1->array1 = CDoubleArray(count);
                testDb1->array2 = CDoubleArray(count);
            } else {
                testDb1->array2 = testDb1->array1;
                testDb1->array1 = CDoubleArray(0);
            }
            testDb1->val1 = testDb1->val2 = count;
        }
        return outObj;
    }
   

    XLTestDouble2(): CObject(TYPE), count(0) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLTestDouble2, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLTestDouble2);
        FIELD(count, "Length to make list");
        Addin::registerClassObjectMethod("TEST_DOUBLE2",
                                         Addin::XL_TESTS,
                                         "Creates a list of given length",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)run);
    }
    static IObject* defaultXLTestDouble2(){
        return new XLTestDouble2();
    }
};

CClassConstSP const XLTestDouble2::TYPE = CClass::registerClassLoadMethod(
    "XLTestDouble2", typeid(XLTestDouble2), load);

class XLTestInt1: public CObject{
public:
    static CClassConstSP const TYPE;
    CIntArray array1;
    int       val1;
    int       val2;
    CIntArray array2;

    static IObjectSP run(XLTestInt1 *params){
        /* return original object - but add 1 to values */
        if (params->array1.size() != params->array2.size()){
            throw ModelException("XLTestInt1::run", "Arrays are of "
                                 "different length");
        }
        XLTestInt1  *addin = copy(params);
        for (int idx = 0; idx < addin->array1.size(); idx++) {
            addin->array1[idx]++;
            addin->array2[idx]++;
        }
        addin->val1++;
        addin->val2++;
        return IObjectSP(addin);
    }
    

    XLTestInt1(): CObject(TYPE), val1(-7), val2(0){
        array1 = CIntArray(3);
        array1[0] = 3;
        array1[1] = -2;
        array1[2] = 0;
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLTestInt1, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLTestInt1);
        FIELD(array1, "First array of ints");
        FIELD_MAKE_OPTIONAL(array1);
        FIELD(val1, "Another int");
        FIELD_MAKE_OPTIONAL(val1);
        FIELD(val2, "Yet another int");
        FIELD(array2, "Array of ints of equal length to first");
        Addin::registerClassObjectMethod("TEST_INT1",
                                         Addin::XL_TESTS,
                                         "Adds one to each number",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)run);
    }
    static IObject* defaultXLTestInt1(){
        return new XLTestInt1();
    }
};

CClassConstSP const XLTestInt1::TYPE = CClass::registerClassLoadMethod(
    "XLTestInt1", typeid(XLTestInt1), load);

class XLTestInt2: public CObject{
public:
    static CClassConstSP const TYPE;
    int       count;

    static IObjectSP run(XLTestInt2 *addin){
        int             count    = Maths::max(0, addin->count);
        IObjectSP       outObj;
        if (count < 10) {
            CIntArraySP  intList(new CIntArray(count));
            for (int idx = 0; idx < count; idx++) {
                (*intList)[idx] = idx;
            }
            outObj = intList;
        } else {
            XLTestInt1 *testInt1 = new XLTestInt1;
            count -= 10;
            outObj = IObjectSP(testInt1);
            if (count > 0) {
                testInt1->array1 = CIntArray(count);
                testInt1->array2 = CIntArray(count);
            } else {
                testInt1->array2 = testInt1->array1;
                testInt1->array1 = CIntArray(0);
            }
            testInt1->val1 = testInt1->val2 = count;
        }
        return outObj;
    }
   

    XLTestInt2(): CObject(TYPE), count(0) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLTestInt2, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLTestInt2);
        FIELD(count, "Length to make list");
        Addin::registerClassObjectMethod("TEST_INT2",
                                         Addin::XL_TESTS,
                                         "Creates a list of given length",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)run);
    }
    static IObject* defaultXLTestInt2(){
        return new XLTestInt2();
    }
};

CClassConstSP const XLTestInt2::TYPE = CClass::registerClassLoadMethod(
    "XLTestInt2", typeid(XLTestInt2), load);

class XLTestString1: public CObject{
public:
    static CClassConstSP const TYPE;
    CStringArray array1;
    string       val1;
    string       val2;
    CStringArray array2;

    static IObjectSP run(XLTestString1 *params){
        /* return original object - but add 1 to values */
        if (params->array1.size() != params->array2.size()){
            throw ModelException("XLTestString1::run", "Arrays are of "
                                 "different length");
        }
        XLTestString1  *addin = copy(params);
        for (int idx = 0; idx < addin->array1.size(); idx++) {
            if (!addin->array1[idx].empty()){
                addin->array1[idx][0]++;
            }
            if (!addin->array2[idx].empty()){
                addin->array2[idx][0]++;
            }
        }
        if (!addin->val1.empty()){
            addin->val1[0]++;
        }
        if (!addin->val2.empty()){
            addin->val2[0]++;
        }
        return IObjectSP(addin);
    }
    

    XLTestString1(): CObject(TYPE), val1("hi") {
        array1 = CStringArray(3);
        array1[0] = "c";
        array1[1] = "b";
        array1[2] = "a";
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLTestString1, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLTestString1);
        FIELD(array1, "First array of strings");
        FIELD_MAKE_OPTIONAL(array1);
        FIELD(val1, "Another string");
        FIELD_MAKE_OPTIONAL(val1);
        FIELD(val2, "Yet another string");
        FIELD(array2, "Array of strings of equal length to first");
        Addin::registerClassObjectMethod("TEST_STRING1",
                                         Addin::XL_TESTS,
                                         "Adds one to each number",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)run);
    }
    static IObject* defaultXLTestString1(){
        return new XLTestString1();
    }
};

CClassConstSP const XLTestString1::TYPE = CClass::registerClassLoadMethod(
    "XLTestString1", typeid(XLTestString1), load);

class XLTestString2: public CObject{
public:
    static CClassConstSP const TYPE;
    string       count;

    static IObjectSP run(XLTestString2 *addin){
        int  count    = addin->count.empty()?
            0: Maths::max(0, addin->count[0] - 'J');
        IObjectSP       outObj;
        if (count < 10) {
            CStringArraySP  stringList(new CStringArray(count));
            for (int idx = 0; idx < count; idx++) {
                (*stringList)[idx] = 'A' + (char) idx;
            }
            outObj = stringList;
        } else {
            XLTestString1 *testString1 = new XLTestString1;
            count -= 10;
            outObj = IObjectSP(testString1);
            if (count > 0) {
                testString1->array1 = CStringArray(count);
                testString1->array2 = CStringArray(count);
                for (int i = 0; i < count; i++){
                    testString1->array1[i] = 'A' + (char) i;
                    testString1->array2[i] = 'Z' - (char) i;
                }
            } else {
                testString1->array2 = testString1->array1;
                testString1->array1 = CStringArray(0);
            }
            testString1->val1 = testString1->val2 = '0' + (char) count;
        }
        return outObj;
    }
   

    XLTestString2(): CObject(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLTestString2, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLTestString2);
        FIELD(count, "Length to make list");
        Addin::registerClassObjectMethod("TEST_STRING2",
                                         Addin::XL_TESTS,
                                         "Creates a list of given length",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)run);
    }
    static IObject* defaultXLTestString2(){
        return new XLTestString2();
    }
};

CClassConstSP const XLTestString2::TYPE = CClass::registerClassLoadMethod(
    "XLTestString2", typeid(XLTestString2), load);


class XLTestBool1: public CObject{
public:
    static CClassConstSP const TYPE;
    CBoolArray array1;
    bool       val1;
    bool       val2;
    CBoolArray array2;

    static IObjectSP run(XLTestBool1 *params){
        /* return original object - but add 1 to values */
        if (params->array1.size() != params->array2.size()){
            throw ModelException("XLTestBool1::run", "Arrays are of "
                                 "different length");
        }
        XLTestBool1  *addin = copy(params);
        if (addin->val1 && addin->val2){
            addin->array1.pop_back();
            addin->array2.pop_back();
        }

        for (int idx = 0; idx < addin->array1.size(); idx++) {
            addin->array1[idx] = ! addin->array1[idx];
            addin->array2[idx] = ! addin->array2[idx];
        }
        
        addin->val1 = !addin->val1;
        addin->val2 = !addin->val2;
        return IObjectSP(addin);
    }
    

    XLTestBool1(): CObject(TYPE), val1(true), val2(false){
        array1 = CBoolArray(3);
        array1[0] = true;
        array1[1] = false;
        array1[2] = true;
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLTestBool1, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLTestBool1);
        FIELD(array1, "First array of bools");
        FIELD_MAKE_OPTIONAL(array1);
        FIELD(val1, "Another bool");
        FIELD_MAKE_OPTIONAL(val1);
        FIELD(val2, "Yet another bool");
        FIELD(array2, "Array of bools of equal length to first");
        Addin::registerClassObjectMethod("TEST_BOOL1",
                                         Addin::XL_TESTS,
                                         "Adds one to each number",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)run);
    }
    static IObject* defaultXLTestBool1(){
        return new XLTestBool1();
    }
};

CClassConstSP const XLTestBool1::TYPE = CClass::registerClassLoadMethod(
    "XLTestBool1", typeid(XLTestBool1), load);

class XLTestBool2: public CObject{
public:
    static CClassConstSP const TYPE;
    bool       count;

    XLTestBool2(): CObject(TYPE), count(0){}

    static IObjectSP run(XLTestBool2 *params){
        XLTestBool1* test1 = new XLTestBool1();
        IObjectSP  returnObj(test1);
        if (params->count){
            test1->array1 = CBoolArray(0);
            test1->array2 = CBoolArray(2);
        }
        return returnObj;
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLTestBool2, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLTestBool2);
        FIELD(count, "Length to make list");
        Addin::registerClassObjectMethod("TEST_BOOL2",
                                         Addin::XL_TESTS,
                                         "Creates a list of given length",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)run);
    }
    static IObject* defaultXLTestBool2(){
        return new XLTestBool2();
    }
};

CClassConstSP const XLTestBool2::TYPE = CClass::registerClassLoadMethod(
    "XLTestBool2", typeid(XLTestBool2), load);

class XLTestDateTime1: public CObject{
public:
    static CClassConstSP const TYPE;
    DateTimeArray array1;
    DateTime      val1;
    DateTime      val2;
    DateTimeArray array2;

    static IObjectSP run(XLTestDateTime1 *params){
        /* return original object - but add 1 to values */
        if (params->array1.size() != params->array2.size()){
            throw ModelException("XLTestDateTime1::run", "Arrays are of "
                                 "different length");
        }
        XLTestDateTime1  *addin = copy(params);
        for (int idx = 0; idx < addin->array1.size(); idx++) {
            addin->array1[idx] = addin->array1[idx].rollDate(1);
            addin->array2[idx] = addin->array2[idx].rollDate(1);
        }
        DateTime magicDate("26-DEC-2020", "EOD");
        if (magicDate.equals(addin->val1)){
            addin->array2 = DateTimeArray(2);
            addin->array2[0] = magicDate;
            addin->array2[1] = magicDate.rollDate(1);
        }
        addin->val1 = addin->val1.rollDate(1);
        addin->val2 = addin->val2.rollDate(1);
        return IObjectSP(addin);
    }
    

    XLTestDateTime1(): CObject(TYPE), val1("25-Dec-1999",
                                           DateTime::START_OF_DAY){
        array1 = DateTimeArray(3);
        array1[0] = DateTime("25-Dec-2000", DateTime::START_OF_DAY);
        array1[1] = DateTime("25-Dec-2001", DateTime::START_OF_DAY);
        array1[2] = DateTime("25-Dec-2002", DateTime::END_OF_DAY);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLTestDateTime1, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLTestDateTime1);
        FIELD(array1, "First array of dateTimes");
        FIELD_MAKE_OPTIONAL(array1);
        FIELD(val1, "Another dateTime");
        FIELD_MAKE_OPTIONAL(val1);
        FIELD(val2, "Yet another dateTime");
        FIELD(array2, "Array of dateTimes of equal length to first");
        Addin::registerClassObjectMethod("TEST_DATETIME1",
                                         Addin::XL_TESTS,
                                         "Adds one to each number",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)run);
    }
    static IObject* defaultXLTestDateTime1(){
        return new XLTestDateTime1();
    }
};

CClassConstSP const XLTestDateTime1::TYPE = CClass::registerClassLoadMethod(
    "XLTestDateTime1", typeid(XLTestDateTime1), load);


class XLTestDoubleMatrix: public CObject{
public:
    static CClassConstSP const TYPE;
    CDoubleArray      strikes;
    DateTimeArray     dates;
    CDoubleMatrix     matrix;

    static IObjectSP run(XLTestDoubleMatrix *params){
        /* return original object -  - but double strikes */
        XLTestDoubleMatrix  *addin = copy(params);
        if (addin->strikes.size() == 2 &&
            addin->strikes[0] == 0.0){
            // magic case
            addin->matrix = CDoubleMatrix();
        } else {
            for (int idx = 0; idx < addin->strikes.size(); idx++) {
                addin->strikes[idx] *= 1.1;
            }
        }
        return IObjectSP(addin);
    }
    

    XLTestDoubleMatrix(): CObject(TYPE) {};

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLTestDoubleMatrix, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLTestDoubleMatrix);
        FIELD(strikes, "strikes");
        FIELD_MAKE_OPTIONAL(strikes);
        FIELD(dates, "dates");
        FIELD(matrix, "matrix");
;
        Addin::registerClassObjectMethod("TEST_DOUBLE_MAT",
                                         Addin::XL_TESTS,
                                         "Scales strikes",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)run);
    }
    static IObject* defaultXLTestDoubleMatrix(){
        return new XLTestDoubleMatrix();
    }
};

CClassConstSP const XLTestDoubleMatrix::TYPE = CClass::registerClassLoadMethod(
    "XLTestDoubleMatrix", typeid(XLTestDoubleMatrix), load);


class XLTestDoubleMatrix2: public CObject{
public:
    static CClassConstSP const TYPE;
    CDoubleArray      strikes;
    DateTimeArray     dates;
    CIntArray         intList;
    CStringArray      stringList;
    CDoubleMatrix     matrix;


    static IObjectSP run(XLTestDoubleMatrix2 *params){
        /* return empty lists for all components */
        return IObjectSP(new XLTestDoubleMatrix2());
    }
    

    XLTestDoubleMatrix2(): CObject(TYPE) {};

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(XLTestDoubleMatrix2, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLTestDoubleMatrix2);
        FIELD(strikes, "strikes");
        FIELD_MAKE_OPTIONAL(strikes);
        FIELD(dates, "dates");
        FIELD_MAKE_OPTIONAL(dates);
        FIELD(intList, "integers");
        FIELD_MAKE_OPTIONAL(intList);
        FIELD(stringList, "strings");
        FIELD_MAKE_OPTIONAL(stringList);
        FIELD(matrix, "matrix");
        FIELD_MAKE_OPTIONAL(matrix);
        Addin::registerClassObjectMethod("TEST_DOUBLE_MAT2",
                                         Addin::XL_TESTS,
                                         "Tests empty arrays etc",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)run);
    }
    static IObject* defaultXLTestDoubleMatrix2(){
        return new XLTestDoubleMatrix2();
    }
};

CClassConstSP const XLTestDoubleMatrix2::TYPE = 
CClass::registerClassLoadMethod("XLTestDoubleMatrix2", 
                                typeid(XLTestDoubleMatrix2), load);

class XLTestZeroParams: public CObject{
public:
    static CClassConstSP const TYPE;

    XLTestZeroParams(): CObject(TYPE){}

    string zp1(){
        return "Hello";
    }
    
    static double zp2(void *params){
        return 3.1415;
    }
    bool zp3(){
        return true;
    }
    int zp4(){
        return -1;
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLTestZeroParams, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLTestZeroParams);
        Addin::registerStringMethod("TEST_ZERO_PARAM1",
                                    Addin::XL_TESTS,
                                    "Tests funcs taking no inputs",
                                    &XLTestZeroParams::zp1);
        Addin::registerClassDoubleMethod("TEST_ZERO_PARAM2",
                                         Addin::XL_TESTS,
                                         "Tests funcs taking no inputs",
                                         TYPE, zp2);
        Addin::registerBoolMethod("TEST_ZERO_PARAM3",
                                  Addin::XL_TESTS,
                                  "Tests funcs taking no inputs",
                                  &XLTestZeroParams::zp3);
        Addin::registerIntMethod("TEST_ZERO_PARAM4",
                                 Addin::XL_TESTS,
                                 "Tests funcs taking no inputs",
                                 &XLTestZeroParams::zp4);
    }
    static IObject* defaultXLTestZeroParams(){
        return new XLTestZeroParams();
    }
};

CClassConstSP const XLTestZeroParams::TYPE = 
CClass::registerClassLoadMethod("XLTestZeroParams", 
                                typeid(XLTestZeroParams), load);


class XLTestExpiry1: public CObject{
public:
    static CClassConstSP const TYPE;
    ExpiryArray array1;
    ExpirySP    val1;
    ExpirySP    val2;
    ExpiryArray array2;

    static IObjectSP run(XLTestExpiry1 *params){
        /* return original object - but add 1 to values */
        if (params->array1.size() != params->array2.size()){
            throw ModelException("XLTestExpiry1::run", "Arrays are of "
                                 "different length");
        }
        XLTestExpiry1  *addin = copy(params);
        DateTime magicDate("26-DEC-2020", "EOD");
        BenchmarkDate magicExpiry(magicDate);
        if (magicExpiry.equals(addin->val1.get())){
            addin->array2 = ExpiryArray(2);
            addin->array2[0] = ExpirySP(new BenchmarkDate(magicDate));
            addin->array2[1] = 
                ExpirySP(new BenchmarkDate(magicDate.rollDate(1)));
        }
        return IObjectSP(addin);
    }
    

    XLTestExpiry1(): CObject(TYPE), val1(new MaturityPeriod(2, "Q")){
        array1 = ExpiryArray(3);
        array1[0] = ExpirySP(
            new BenchmarkDate(DateTime("25-Dec-2000",DateTime::START_OF_DAY)));
        array1[1] = ExpirySP(
            new BenchmarkDate(DateTime("25-Dec-2001",DateTime::START_OF_DAY)));
        array1[2] = ExpirySP(
            new BenchmarkDate(DateTime("25-Dec-2002",DateTime::END_OF_DAY)));
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLTestExpiry1, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLTestExpiry1);
        FIELD(array1, "First array of Expiries");
        FIELD_MAKE_OPTIONAL(array1);
        FIELD(val1, "Another Expiry");
        FIELD_MAKE_OPTIONAL(val1);
        FIELD(val2, "Yet another Expiry");
        FIELD(array2, "Array of Expiries of equal length to first");
        Addin::registerClassObjectMethod("TEST_EXPIRY",
                                         Addin::XL_TESTS,
                                         "Tests expiries",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)run);
    }
    static IObject* defaultXLTestExpiry1(){
        return new XLTestExpiry1();
    }
};

CClassConstSP const XLTestExpiry1::TYPE = CClass::registerClassLoadMethod(
    "XLTestExpiry1", typeid(XLTestExpiry1), load);

class XLTestMktObjWrapper: public CObject{
public:
    static CClassConstSP const TYPE;
    CAssetWrapperArray   assets;      // array of assets
    YieldCurveWrapper    yc;
 

    XLTestMktObjWrapper(): CObject(TYPE) {};

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLTestMktObjWrapper, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLTestMktObjWrapper);
        FIELD(assets, "array of assets");
        FIELD(yc, "single asset");
        Addin::registerConstructor("TEST_MARKET_WRAPPER",
                                   Addin::XL_TESTS,
                                   "Tests market wrappers",
                                   TYPE);
    }
    static IObject* defaultXLTestMktObjWrapper(){
        return new XLTestMktObjWrapper();
    }
};

CClassConstSP const XLTestMktObjWrapper::TYPE =CClass::registerClassLoadMethod(
    "XLTestMktObjWrapper", typeid(XLTestMktObjWrapper), load);
#if 0
class XLLotsOfParams: public CObject{
    double   val1;
    double   val2;
    double   val3;
    double   val4;
    double   val5;
    double   val6;
    double   val7;
    double   val8;
    double   val9;
    double   val10;
    double   val11;
    double   val12;
    double   val13;
    double   val14;
    double   val15;
    double   val16;
    double   val17;
    double   val18;
    double   val19;
    double   val20;
    double   val21;
    double   val22;
    double   val23;
    double   val24;
    double   val25;
    
public:
    static CClassConstSP const TYPE;

    double notMuch(){
        return (1.0);
    }

    XLLotsOfParams(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(XLLotsOfParams, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLLotsOfParams);
        FIELD_NO_DESC(val1);
        FIELD_NO_DESC(val2);
        FIELD_NO_DESC(val3);
        FIELD_NO_DESC(val4);
        FIELD_NO_DESC(val5);
        FIELD_NO_DESC(val6);
        FIELD_NO_DESC(val7);
        FIELD_NO_DESC(val8);
        FIELD_NO_DESC(val9);
        FIELD_NO_DESC(val10);
        FIELD_NO_DESC(val11);
        FIELD_NO_DESC(val12);
        FIELD_NO_DESC(val13);
        FIELD_NO_DESC(val14);
        FIELD_NO_DESC(val15);
        FIELD_NO_DESC(val16);
        FIELD_NO_DESC(val17);
        FIELD_NO_DESC(val18);
        FIELD_NO_DESC(val19);
        FIELD_NO_DESC(val20);
        FIELD_NO_DESC(val21);
        FIELD_NO_DESC(val22);
        FIELD_NO_DESC(val23);
        FIELD_NO_DESC(val24);
        FIELD_NO_DESC(val25);
        Addin::registerDoubleMethod("TEST_LOTS_PARAMS",
                                    Addin::XL_TESTS,
                                    "Tests addin func with 25 parameters",
                                    &XLLotsOfParams::notMuch);
    }

    static IObject* defaultXLLotsOfParams(){
        return new XLLotsOfParams();
    }
};

CClassConstSP const XLLotsOfParams::TYPE = CClass::registerClassLoadMethod(
    "XLLotsOfParams", typeid(XLLotsOfParams), load);
#endif

/* non static variable to force linkage of this file. Outside of class to
   avoid necessity of header file */
bool XLTestsRegistered = true;


DRLIB_END_NAMESPACE
