//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MatrixResult.cpp
//
//   Description : Captures result of tweaking a matrix
//
//   Author      : Mark A Robson
//
//   Date        : 03 Apr 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MatrixResult.hpp"
#include "edginc/Maths.hpp"
#include "edginc/MatrixShift.hpp"
#include "edginc/ExpiryResult.hpp"
#include "edginc/NotApplicable.hpp"
#include "edginc/Untweakable.hpp"

DRLIB_BEGIN_NAMESPACE

// PUBLIC methods

/** Constructor - note: does not clone inputs */
MatrixResult::MatrixResult(const CDoubleArraySP&     strikes,
                           const ExpiryArrayConstSP& expiries,
                           const CDoubleMatrixSP&    matrix):
    CObject(TYPE), strikesSP(strikes), expiriesSP(expiries), matrix(matrix){
    validatePop2Object();
}

MatrixResult::MatrixResult():
    CObject(TYPE),
    strikesSP(CDoubleArray::SP()),
    expiriesSP(ExpiryArray::SP()),
    matrix(new CDoubleMatrix())
{}

//// ensures that the size of the arrays and matrix are consistent
void MatrixResult::validatePop2Object(){
    if ( matrix->numCols() != strikesSP->size() ||
         matrix->numRows() != expiriesSP->size()) {
        throw ModelException("MatrixResult::validatePop2Object",
                             "Malformed MatrixResult");
    }
}

/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void MatrixResult::outputWrite(const string& linePrefix,
                               const string& prefix, ostream& stream) const
{
    ios::fmtflags oldFlags     = stream.flags();
    long          oldPrecision = stream.precision();
    char buffer[512]; // can cope with biggest possible double
    CDoubleArray&   strikes = *strikesSP;
    const ExpiryArray&    expiries = *expiriesSP;

    for (int x = 0; x < strikes.size(); x++){
        for (int y = 0; y < expiries.size(); y++){
            sprintf(buffer, "%.6f", strikes[x]);
            double db = (*matrix)[x][y];
            if (db < 1 && db > -1){
                // 8 dp
                stream.flags(ios::fixed);
            }
            else {
                // 8 sf
                stream.unsetf(ios::fixed);
            }        
            stream.precision(8);
            
            stream << linePrefix << prefix << "_" << buffer << "_" <<
                expiries[y]->toString() << ": " << db << endl;
        }
    }
    stream.flags(oldFlags);
    stream.precision(oldPrecision);
}

/** scale by factor x */
void MatrixResult::scale(double x) {
    matrix->scale(x);
}

/** add MatrixResult object (scaled by scaleFactor) to this result
    (Implementation of CombinableResult) */
void MatrixResult::add(const CombinableResult& x, double scaleFactor) {
    static const string method = "MatrixResult::add";
    // gcc bug: force to IObject before dynamic cast
    try {
        const MatrixResult& mr =dynamic_cast<const MatrixResult&>(static_cast<const IObject&>(x));
        const ExpiryArray&  expiries1 = *expiriesSP;
        const DoubleArray&  strikes1  = *strikesSP;
        const ExpiryArray&  expiries2 = *mr.expiriesSP;
        const DoubleArray&  strikes2  = *mr.strikesSP;
        const CDoubleMatrix matrix1   = *matrix;
        const CDoubleMatrix matrix2   = *mr.matrix;
        // first check that expiries match (they always should if we've used
        // the same market data)
        bool match = false;
        if (expiries1.size() == expiries2.size()){
            match = true;
            for (int i = 0; i < expiries1.size() && match; i++){
                if (!expiries1[i]->equals(expiries2[i].get())){
                    match = false;
                }
            }
        }
        if (!match){
            throw ModelException(method, "Mismatch between expiry arrays");
        }

        // create new double matrix of maximum size needed (trim later)
        CDoubleMatrixSP newDbMatrix(new CDoubleMatrix(strikes1.size()+
                                                      strikes2.size(),
                                                      expiries1.size()));
        CDoubleMatrix&  newMatrix = *newDbMatrix;
        CDoubleArraySP  newDbXAxis(new CDoubleArray(0));
        CDoubleArray&   newXAxis = *newDbXAxis;
        newXAxis.reserve(strikes1.size()+strikes2.size());
        int i1, i2, j;
        for (i1 = 0, i2 = 0, j = 0; 
             !(i1 == strikes1.size() && i2 == strikes2.size()); j++){
            // NB Be careful with the boundary cases
            bool okayToCompare = i1 < strikes1.size() && i2 < strikes2.size();
            if (i2 == strikes2.size() ||
                (okayToCompare && strikes1[i1] < strikes2[i2] && 
                 strikes2[i2]/strikes1[i1] >= MatrixShift::BETA)){
                     newXAxis.push_back(strikes1[i1]);
                for (int i = 0; i < expiries1.size(); i++){
                    newMatrix[j][i] = matrix1[i1][i];
                }
                i1++;
            } else if (i1 == strikes1.size() ||
                       (okayToCompare && strikes1[i1] > strikes2[i2] && 
                        strikes1[i1]/strikes2[i2] >= MatrixShift::BETA)){
                newXAxis.push_back(strikes2[i2]);
                for (int i = 0; i < expiries1.size(); i++){
                    newMatrix[j][i] = scaleFactor * matrix2[i2][i];
                }
                i2++;
            } else {
                // merge strikes
                newXAxis.push_back(Maths::min(strikes1[i1], strikes2[i2]));
                for (int i = 0; i < expiries1.size(); i++){
                    newMatrix[j][i] = matrix1[i1][i] + 
                        scaleFactor * matrix2[i2][i];
                }
                i1++; i2++;
            }
        }
        // trim matrix
        for (/* where we are */; j < strikes1.size()+strikes2.size(); j++){
            newMatrix.removeLastCol();
        }

        // finally update this object
        strikesSP = newDbXAxis;
        matrix = newDbMatrix;

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** create a new CombinableResult object by adding an object (scaled by
    scaleFactor) to this result. Supports adding ExpiryResultArray,
    NotApplicable and Untweakable */
IObject* MatrixResult::addResult(const  IObject& x, 
                                 double scaleFactor) const{
    static const string method = "MatrixResult::addResult";
    if (ExpiryResultArray::TYPE->isInstance(&x)){
        IObjectSP arrayCopy(x.clone());
        ExpiryResultArray& resultArray = 
            dynamic_cast<ExpiryResultArray&>(*arrayCopy);
        // first check that expiries match (they always should if we've used
        // the same market data)
        bool match = false;
        if (resultArray.size() == expiriesSP->size()){
            match = true;
            for (int i = 0; i < resultArray.size() && match; i++){
                if (!resultArray[i].getExpiry()->equals(
                    (*expiriesSP)[i].get())){
                    match = false;
                }
            }
        }
        if (!match){
            throw ModelException(method, "Mismatch between expiry arrays");
        }
        for (int i = 0; i < resultArray.size(); i++){
            // then sum the sensitivities for the different strikes into 1
            double sens = 0.0;
            for (int j = 0; j < strikesSP->size();j++){
                sens += (*matrix)[j][i];
            }
            resultArray[i].scale(scaleFactor); // scale the original
            resultArray[i].addToResult(sens);  // then add this to it
        }
        return arrayCopy.release();
    } else if (!CombinableResult::TYPE->isInstance(&x)){
        // fall through to error case below
    } else {
        const CombinableResult& cr = dynamic_cast<const CombinableResult&>(x);
        if (TYPE->isInstance(&x)){
            MatrixResultSP copyMS(copy(this));
            copyMS->add(cr, scaleFactor);
            return copyMS.release();
        } else if (NotApplicable::TYPE->isInstance(&x)){
            return copy(this);
        } else if (Untweakable::TYPE->isInstance(&x)){
            return new Untweakable(string("MatrixShift::addResult: Component "
                                          "was untweakable"));
        }
    }
    throw ModelException("MatrixShift::addResult", "Can't add MatrixShift"
                         " to "+x.getClass()->getName());
}

// this is really lame, should be in CArray

static int indexOf(const DoubleArray& xs, double x) {
    int s = xs.size() - 1;
    while (s >= 0 && xs[s] != x) --s;
    return s;
}

static int indexOf(const ExpiryArray& xs, ExpiryConstSP x) {
    int e = xs.size() - 1;
    while (e >= 0 && (!x ? !!xs[e] : !x->equalTo(xs[e].get()))) --e;
    return e;
}

double MatrixResult::getOr0(const ExpiryAndStrike& expiryAndStrike) const {
    int s = indexOf(*strikesSP, expiryAndStrike.strike);
    if (s < 0) return 0.;
    int e = indexOf(*expiriesSP, expiryAndStrike.expiry);
    if (e < 0) return 0.;

    return (*matrix)[s][e];
}

void MatrixResult::set(const ExpiryAndStrike& expiryAndStrike, double value) {
    int s = indexOf(*strikesSP, expiryAndStrike.strike);
    if (s < 0) {
        strikesSP.reset(strikesSP.clone());
        strikesSP->push_back(expiryAndStrike.strike);
        s = matrix->numCols();
        matrix->insertCols(s, 1);
    }

    int e = indexOf(*expiriesSP, expiryAndStrike.expiry);
    if (e < 0) {
        expiriesSP.reset(expiriesSP.clone());
        ExpiryArraySP::constCast(expiriesSP)->push_back(ExpirySP::constCast(expiryAndStrike.expiry));
        e = matrix->numRows();
        matrix->insertRows(e, 1);
    }

    (*matrix)[s][e] = value;
}

class MatrixResultHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(MatrixResult, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(CombinableMixedResult);
        EMPTY_SHELL_METHOD(defaultMatrixResult);
        FIELD(strikesSP,  "strikes");
        FIELD(expiriesSP,  "expiries");
        FIELD(matrix, "values");
    }

    static IObject* defaultMatrixResult(){
        return new MatrixResult();
    }
};

CClassConstSP const MatrixResult::TYPE = CClass::registerClassLoadMethod(
    "MatrixResult", typeid(MatrixResult), MatrixResultHelper::load);


DRLIB_END_NAMESPACE
