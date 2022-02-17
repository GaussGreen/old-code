
#ifndef EDG_DOUBLEMATRIX_H
#define EDG_DOUBLEMATRIX_H
#include "edginc/AtomicArray.hpp"
#include "edginc/CombinableResult.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/mathlib.hpp"

DRLIB_BEGIN_NAMESPACE
/** For recording eigenvalues */
class EigenVectorAnalysis;

typedef smartConstPtr<EigenVectorAnalysis> EigenVectorAnalysisConstSP;
typedef smartPtr<EigenVectorAnalysis> EigenVectorAnalysisSP;
#ifndef QLIB_DOUBLEMATRIX_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<EigenVectorAnalysis>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<EigenVectorAnalysis>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<EigenVectorAnalysis>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<EigenVectorAnalysis>);
#endif

/** For recording SVD */
class SVDAnalysis;

typedef smartConstPtr<SVDAnalysis> SVDAnalysisConstSP;
typedef smartPtr<SVDAnalysis> SVDAnalysisSP;
#ifndef QLIB_DOUBLEMATRIX_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<SVDAnalysis>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<SVDAnalysis>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<SVDAnalysis>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<SVDAnalysis>);
#endif

/** For recording beta factorization */
class BetaFactorization;

typedef smartConstPtr<BetaFactorization> BetaFactorizationConstSP;
typedef smartPtr<BetaFactorization> BetaFactorizationSP;
#ifndef QLIB_DOUBLEMATRIX_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<BetaFactorization>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<BetaFactorization>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<BetaFactorization>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<BetaFactorization>);
#endif


/** Represents a m x n matrix where the elements are doubles. */
class UTIL_DLL CDoubleMatrix: public CObject,
                     public virtual CombinableResult{
    friend class ShutTheCompilerUp;
public:
    static CClassConstSP const TYPE;

    /** tag for recording number of columns in a matrix */
    static const string NUM_COLUMNS;
    /** tag for recording number of rows in a matrix */
    static const string NUM_ROWS;

    /** Creates an 0x0 matrix */
    CDoubleMatrix();
    
    /** Creates an mxn matrix, the elements of which are set to 0.0 */
    CDoubleMatrix(int numCols, int numRows);
    /** Creates an mxn matrix, the elements of which are itialised using
        the given data */
    CDoubleMatrix(int numCols, int numRows, const double** data);

    /** Creates a 1 x n matrix (i.e a curve) from a double array */
    explicit CDoubleMatrix(const DoubleArray& theArray);

    /** Creates a m x n matrix (i.e a curve) from a double array */
    CDoubleMatrix(int numCols, int numRows, const DoubleArray& theArray);

    /** Creates a n x m matrix from a 2-dimensional [n][m] double array */
    explicit CDoubleMatrix(const DoubleArrayArray& theArray);

    /** Copy constructor */
    CDoubleMatrix(const CDoubleMatrix& rhs);

    /** Hash code function */
    virtual int hashCode() const;

    /** Comparison function */
    virtual bool equalTo(const IObject* obj) const;

    /** assignment operator */
    CDoubleMatrix& operator=(const CDoubleMatrix& rhs);

    virtual IObject* clone() const;    // inherited 

    /** Returns a column of the matrix with 0 <= colIndex < numCols().
        Inline for performance */
    double* operator[](int colIndex);

    /** Returns a column of the matrix with 0 <= colIndex < numCols()
        Inline for performance */
    const double* operator[](int colIndex) const;

    /** Returns the number of columns in the matrix. Inline for performance */
    int numCols() const;

    /** Returns the number if rows in the matrix. Inline for performance */
    int numRows() const;

    /** Returns true if the double matrix is empty ie 0x0 dimensions */
    bool empty() const;

    /** Returns true of the double matrix is the identiy matrix */
    bool isIdentity() const;

    /** write object out in XML format */
    virtual void write(const string& tag, Writer* writer) const;

    /** populate an empty object from XML description */
    virtual void import(Reader::Node* elem, Reader* reader);

    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    /** resize the matrix to a new number of columns and rows.
        do nothing with the values. */
    void resizeOnly(int numCols, int numRows);
    
    /** Resize the matrix to a new number of columns and rows. 
        Fill it with 0s if the size has been changed */
    void resize(int numCols, int numRows);

	/** set all the element of the matrix to a given number */
    void fill(double val);
   
    /** resize the matrix to a new number of columns and rows and
        set all elements to val (by default, val = 0.0)
    **/
    void resizeAndFill(int numCols, int numRows, double val = 0.0);

    /** checks all values are >= 0 (actually, exception if <-DBL_EPSILON)*/
    void checkNonNegative() const;

    /** Checks that all values are >0 (actually, exception if <DBL_EPSILON) */
    void checkPositive() const;

    /** checks matrix is symmetric */
    void checkSymmetric() const;

    /** adds given value to all elements of matrix */
    void scalarAdd(double val);

    /** adds given value to all elements of given row in matrix */
    void rowAdd(int rowIndex, double val);

    /** adds given value to all elements of given column in matrix */
    void colAdd(int colIndex, double val);

    /** Multiplies a specified row by a given value */
    void rowMultiply(int rowIndex, double val);

    /** Multiplies a specified column by a given value */
    void colMultiply(int colIndex, double val);

    /** Removes last column of matrix */
    void removeLastCol();

    void insertCols(int start, int num);

    /** Removes numRowsToRemove rows from matrix starting at startRow */
    void removeRows(int startRow, int numRowsToRemove);

    void insertRows(int start, int num);

    /** transpose (even for m*n matrix) */
    void transpose();

    /** squeeze, ie compute rho + (1-rho)*s, s>0, resp rho(1+s), s<0 */
    void squeeze(double s);

    /** scale by factor x */
    void scale(double x);
    static void scale(const CDoubleMatrix& x, 
                      double gFactor, CDoubleMatrix * gResult);
    static void scale(double gFactor, 
                      const CDoubleMatrix& x, CDoubleMatrix * gResult);

    /** Negate all entries (equivalent to scale(-1.0) */
    void negate();

    /** add a CDoubleMatrix (scaled by scaleFactor) to this
        result.  */
    virtual void add(const CombinableResult& x, double scaleFactor);

    /** Arithmetic operations */
    // Returns this * x. Dimensions must match.
    CDoubleMatrix mult(const CDoubleMatrix& x) const;
    static void mult(const CDoubleMatrix& x, 
                     const CDoubleMatrix& y, CDoubleMatrix * gResult);
    // Returns this * array. Dimensions must match.
    CDoubleArray mult(const CDoubleArray &a) const;
    CDoubleArray mult(const CIntArray& a) const;
    CDoubleArray mult(const vector<int>& a) const;
    // Returns this - x. Dimensions must match.
    CDoubleMatrix minus(const CDoubleMatrix& x) const;
    // Returns this + x. Dimensions must match.
    CDoubleMatrix plus(const CDoubleMatrix& x) const;
    // Operator versions of -, + and *
    CDoubleMatrix operator*(const CDoubleMatrix& x) const;
    CDoubleMatrix operator-(const CDoubleMatrix& x) const;
    CDoubleMatrix operator+(const CDoubleMatrix& x) const;

    /** computes the square root of a matrix aka Choleski / Cholesky etc. */
    CDoubleMatrix computeSquareRoot() const;
    
    /** computes the inverse of a square matrix using LU */
    CDoubleMatrix computeInverse() const;

    /** computes the inverse of a square matrix using LU. No allocation */
    void computeInverse(CDoubleMatrix * yWrapper,  // this is the result
                        CDoubleMatrix * aWrapper,
                        DoublePointerVector * a,
                        DoublePointerVector * yVector,
                        DoubleVector        * colWrapper,
                        IntVector           * indxWrapper) const;

    /** solves for x in Ax=b, with A square matrix and b=rhs */
    CDoubleArray solve(CDoubleArray& rhs) const;

    /** input:  
            symmetric matrix 
        algorithm:
            1) floor the (real) eigenvalues at 0 -> covariance matrix
            2) normalise covariance matrix -> correlation matrix
        output: 
            1) correlation matrix 
            2) average squared difference between symmetric and covariance matrices */
    CDoubleMatrix symmToCorrel(double* squareErr, double eigenValueFloor = 0.0000000001) const;
    
    /** Computes eigenvalues and eigenvectors of a given matrix */
    EigenVectorAnalysisSP computeEigenVectors() const;

    /** Computes the Singular Value Decomposition of a given matrix */
    SVDAnalysisSP computeSVD() const;

    /** Computes the Beta Factorization of a given matrix */
    BetaFactorizationSP computeBetaFactorization(int    nbFactors,
                                                 double precision,
                                                 int    maxIter) const;

    /**Frobenius norm of a matrix is sqrt(Tr(A'A)).*/
    double frobeniusNorm() const;

    /**Spectral norm of a matrix is the square root of the modulus of the biggest eigenvalue of A'A*/
    double spectralNorm() const;

    /**Compute the exponential exp(t M) of this (M) times a constant, t, using a Pade approximation 
    of specified order (numerator and denominator are the same). If |Mt|>=0.5, then t is halved
    repeatedly until the norm is less than 0.5 to ensure that the Pade approximation is
    likely to be in range, and the result is squared the same numer of times for the
    final answer. The argument for the spectral norm is there for efficiency: computing it
    is expensive enough that, for repeated exponentiations with different t values, it
    is better to compute it once only. If you submit a negative spectral norm, a fresh one
    will be calculated. */
    void expPade(     double gT, 
                      double gSpectralNorm, 
                      int gPadeOrder,
                      CDoubleMatrix * gResult,
                      CDoubleMatrix * gAux1,
                      CDoubleMatrix * gAux11,
                      CDoubleMatrix * gAux2,
                      CDoubleMatrix * gAux21,
                      CDoubleMatrix * gAux3) const;

    /** Same as expPade but result is allocated inside the method **/
    CDoubleMatrix exp(  double          gT, 
                        double          gSpectralNorm, 
                        int             gPadeOrder) const;

    ~CDoubleMatrix();
private:
    int      NumCols;        // $unregistered
    int      NumRows;        // $unregistered
    // organised so that data[i] gives a column of values
    double** data;           // $unregistered
    void allocateMem();
    void deallocateMem();
    void populateMatrix(const double** rawData);
    void populateMatrix(const DoubleArray& theArray);
    void populateMatrix(const DoubleArrayArray& theArray);

    static void load(CClassSP& clazz);

    // build matrix from array of double arrays
    static IObjectSP fromArray(const IObjectSP& object, 
                               CClassConstSP    requiredType);
    
    static CoeffsPQ4Pade PADEcoeffs; // could be const ...
};

#if !defined(DEBUG) || defined(QLIB_DOUBLEMATRIX_CPP)
OPT_INLINE double* CDoubleMatrix::operator[](int colIndex){
    #ifdef DEBUG
        if(colIndex >= NumCols)
        {
            throw ModelException("CDoubleMatrix::operator[]", "colIndex >= NumCols");
        }
        if(colIndex < 0)
        {
            throw ModelException("colIndex < 0");
        }
    #endif    
    return data[colIndex];
}

OPT_INLINE const double* CDoubleMatrix::operator[](int colIndex) const{
    #ifdef DEBUG
        if(colIndex >= NumCols)
        {
            throw ModelException("CDoubleMatrix::operator[]", "colIndex >= NumCols");
        }
        if(colIndex < 0)
        {
            throw ModelException("colIndex < 0");
        }
    #endif    
    return data[colIndex];
}

OPT_INLINE int CDoubleMatrix::numCols() const{
    return NumCols;
}

OPT_INLINE int CDoubleMatrix::numRows() const{
    return NumRows;
}
#endif

typedef smartConstPtr<CDoubleMatrix> CDoubleMatrixConstSP;
typedef smartPtr<CDoubleMatrix> CDoubleMatrixSP;
#ifndef QLIB_DOUBLEMATRIX_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<CDoubleMatrix>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<CDoubleMatrix>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<CDoubleMatrix>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<CDoubleMatrix>);
#endif

typedef CDoubleMatrix DoubleMatrix;

typedef array<CDoubleMatrixSP, DoubleMatrix> DoubleMatrixArray;
typedef smartPtr<DoubleMatrixArray> DoubleMatrixArraySP;
typedef smartConstPtr<DoubleMatrixArray> DoubleMatrixArrayConstSP;
#ifndef QLIB_DOUBLEMATRIX_CPP
EXTERN_TEMPLATE(class UTIL_DLL array<CDoubleMatrixSP _COMMA_ DoubleMatrix>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<DoubleMatrixArray>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<DoubleMatrixArray>);
EXTERN_TEMPLATE(IObjectSP UTIL_DLL FieldGetInLine<DoubleMatrix>(DoubleMatrix* t));
EXTERN_TEMPLATE(void UTIL_DLL FieldSetInLine<DoubleMatrix>(DoubleMatrix* t,IObjectSP o));
EXTERN_TEMPLATE(IObjectSP UTIL_DLL FieldGetSmartPtr<CDoubleMatrixSP>(
                    CDoubleMatrixSP* t));
EXTERN_TEMPLATE(void UTIL_DLL FieldSetSmartPtr<CDoubleMatrixSP>(CDoubleMatrixSP* t,
                                                       IObjectSP o));
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL array<CDoubleMatrixSP _COMMA_ DoubleMatrix>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<DoubleMatrixArray>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<DoubleMatrixArray>);
INSTANTIATE_TEMPLATE(IObjectSP UTIL_DLL FieldGetInLine<DoubleMatrix>(DoubleMatrix* t));
INSTANTIATE_TEMPLATE(void UTIL_DLL FieldSetInLine<DoubleMatrix>(DoubleMatrix* t,
                                                       IObjectSP o));
INSTANTIATE_TEMPLATE(IObjectSP UTIL_DLL FieldGetSmartPtr<CDoubleMatrixSP>(
                         CDoubleMatrixSP* t));
INSTANTIATE_TEMPLATE(void UTIL_DLL FieldSetSmartPtr<CDoubleMatrixSP>(CDoubleMatrixSP* t,
                                                            IObjectSP o));
#endif


/** Holds the eigen values and eigen vectors of a matrix */
class UTIL_DLL EigenVectorAnalysis: public CObject {
public:
    static CClassConstSP const TYPE;

    EigenVectorAnalysis();
    EigenVectorAnalysis(int n);
    
    CDoubleMatrixSP eigenVectors;
    DoubleArraySP   eigenValues;
    int             dimension;

private:
    static IObject* defaultEigenVectorAnalysis();
    static void load(CClassSP& clazz);
};

/** Holds the SVD of a matrix */
class UTIL_DLL SVDAnalysis: public CObject {
public:
    static CClassConstSP const TYPE;

    SVDAnalysis();
    SVDAnalysis(int numCols, int numRows);
    
    // A = U diag(sValues) V'
    CDoubleMatrixSP matrixU;
    CDoubleMatrixSP matrixVTranspose;
    DoubleArraySP   sValues;

private:
    static IObject* defaultSVDAnalysis();
    static void load(CClassSP& clazz);
};

/** Holds the Beta Factorization of a matrix */
class UTIL_DLL BetaFactorization: public CObject {
public:
    static CClassConstSP const TYPE;

    BetaFactorization();
    BetaFactorization(double            avgDistance,                      
                      CDoubleMatrixSP   betaFactors, 
                      CDoubleMatrixSP   betaMatrix,
                      double            squeezeLimitLow,
                      double            squeezeLimitHigh);
    
    double          avgDistance;
    CDoubleMatrixSP betaFactors;
    CDoubleMatrixSP betaMatrix; 
    double          squeezeLimitLow;
    double          squeezeLimitHigh;

private:
    static IObject* defaultBetaFactorization();
    static void load(CClassSP& clazz);
};


DRLIB_END_NAMESPACE
#endif
