//----------------------------------------------------------------------------
//
//   Filename    : MatrixResult.hpp
//
//   Description : Captures result of tweaking a matrix
//
//   Author      : Mark A Robson
//
//   Date        : 15 June 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_MATRIX_RESULT_HPP
#define EDG_MATRIX_RESULT_HPP

#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/CombinableMixedResult.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/ExpiryAndStrike.hpp"


DRLIB_BEGIN_NAMESPACE

/** Captures result of tweaking a [volatility] matrix (You could replace the
    axes with arbitrary arrays but probably better to have a separate class) */
class RISKMGR_DLL MatrixResult: public CObject,
                                public CombinableMixedResult {
public:
    static CClassConstSP const TYPE;
    friend class MatrixResultHelper;

    /** Constructor - note: does not clone inputs */
    MatrixResult(const CDoubleArraySP&     xAxis,
                 const ExpiryArrayConstSP& yAxis,
                 const CDoubleMatrixSP&    matrix);

    /** A 0x0 MatrixResult */
    MatrixResult();

    //// ensures that the size of the arrays and matrix are consistent
    virtual void validatePop2Object();

    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    /** scale by factor x */
    virtual void scale(double x);

    /** add MatrixResult object to this result (Implementation of
        CombinableResult) */
    virtual void add(const CombinableResult& x, double scaleFactor);

    /** create a new CombinableResult object by adding an object (scaled by
        scaleFactor) to this result. Supports adding ExpiryResultArray,
        NotApplicable and Untweakable */
    virtual IObject* addResult(const IObject& x, 
                               double scaleFactor) const;

    /** returns 0 if the entry does not exist */
    virtual double getOr0(const ExpiryAndStrike& expiryAndStrike) const;

    /** set an entry, adding a new row and/or column if necessary */
    virtual void set(const ExpiryAndStrike& expiryAndStrike,
                     double value);

private:
    MatrixResult(const MatrixResult &rhs);
    MatrixResult& operator=(const MatrixResult& rhs);

    CDoubleArraySP     strikesSP;
    ExpiryArrayConstSP expiriesSP;
    CDoubleMatrixSP    matrix;  
};

typedef smartPtr<MatrixResult> MatrixResultSP;
typedef smartConstPtr<MatrixResult> MatrixResultConstSP;

DRLIB_END_NAMESPACE

#endif
