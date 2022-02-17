//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Optimizer.hpp
//
//   Description : 
//
//   Date        : 20 May 2002
//
//----------------------------------------------------------------------------
#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/Function.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/OneToOneMapping.hpp"
#include "edginc/Atomic.hpp"

#ifndef EDR_OPTIMIZER_HPP
#define EDR_OPTIMIZER_HPP

DRLIB_BEGIN_NAMESPACE

/** Interface for N-dimensional optimizers */
class UTIL_DLL OptimizerND: public CObject{
public:
    static CClassConstSP const TYPE;

    virtual void minimize(const MFunctionND&  func,
                          const CDoubleArray& xguess,        // initial guess
                          CDoubleArray&       x) const = 0;  // result

    virtual void minimize(const MFunctionND&  func,
                          const CDoubleArray& xguess,        // initial guess
                          const CStringArray& ids,           // identifiers
                          CDoubleArray&       x) const = 0;  // result

    virtual void minimize(const MFunctionND&  func,
                          const CDoubleArray& xguess,        // initial guess
                          CDoubleArray&       x,
                          CDoubleArray&       f) const;

    virtual void minimize(const MFunctionND&  func,
                          const CDoubleArray& xguess,        // initial guess
                          const CStringArray& ids,           // identifiers
                          CDoubleArray&       x,
                          CDoubleArray&       f) const;

  

protected:
    OptimizerND(const CClassConstSP& clazz);

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<OptimizerND> OptimizerNDSP;
typedef smartConstPtr<OptimizerND> OptimizerNDConstSP;

#ifndef QLIB_OPTIMIZER_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<OptimizerND>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<OptimizerND>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<OptimizerND>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<OptimizerND>);
#endif

/** Interface for least-square based optimizers */
class UTIL_DLL LeastSquareOptimizer: public OptimizerND{
public:
    static CClassConstSP const TYPE;

    virtual bool statsIsAvailable() const = 0;
    virtual double getVariance() const = 0;
    virtual const DoubleMatrix& getCovarianceMatrix() const = 0;
    virtual const DoubleMatrix& getCorrelationMatrix() const = 0;
    virtual const DoubleMatrix& getHessianInverse() const = 0;
    virtual const DoubleArray& getLowerLimits() const = 0;
    virtual const DoubleArray& getUpperLimits() const = 0;
    virtual const DoubleArray& getStandardErrors() const = 0;

protected:
    LeastSquareOptimizer(const CClassConstSP& clazz);

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<LeastSquareOptimizer> LeastSquareOptimizerSP;
typedef smartConstPtr<LeastSquareOptimizer> LeastSquareOptimizerConstSP;
#ifndef QLIB_OPTIMIZER_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<LeastSquareOptimizer>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<LeastSquareOptimizer>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<LeastSquareOptimizer>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<LeastSquareOptimizer>);
#endif

/** Wrapper around imsl's nonlin_least_squares */
class UTIL_DLL LevenbergMarquardt: public LeastSquareOptimizer{
public:
    static CClassConstSP const TYPE;
    friend class LevenbergMarquardtHelper;

    virtual void minimize(const MFunctionND&  func,
                          const CDoubleArray& guess,        // initial guess
                          CDoubleArray&       x) const;     // result

    virtual void minimize(const MFunctionND&  func,
                          const CDoubleArray& xguess,        // initial guess
                          const CStringArray& ids,           // identifiers
                          CDoubleArray&       x) const;      // result

    virtual bool statsIsAvailable() const;
    virtual double getVariance() const;
    virtual const DoubleMatrix& getCovarianceMatrix() const;
    virtual const DoubleMatrix& getCorrelationMatrix() const;
    virtual const DoubleMatrix& getHessianInverse() const;
    virtual const DoubleArray& getLowerLimits() const;
    virtual const DoubleArray& getUpperLimits() const;
    virtual const DoubleArray& getStandardErrors() const;

    LevenbergMarquardt();
    void validatePop2Object();

private:
    // registered
    CDoubleSP gradTol;
    CDoubleSP stepTol;
    CDoubleSP relFcnTol;
    CDoubleSP absFcnTol;
    CDoubleSP maxStepSize;
    CIntSP    maxItnNb;
    CIntSP    maxFcnEvalNb;
    CIntSP    maxJacobEvalNb;

    // not registered
    const static double DefaultMachinePrecision;
    const static double DEFAULT_gradTol;
    const static double DEFAULT_stepTol;
    const static double DEFAULT_relFcnTol;
    const static double DEFAULT_absFcnTol;
    const static double DEFAULT_maxStepSize;
    const static int    DEFAULT_maxItnNb;
    const static int    DEFAULT_maxFcnEvalNb;
    const static int    DEFAULT_maxJacobEvalNb;


    mutable bool            haveStats; // $unregistered
    mutable double          sqSigma; // $unregistered
    mutable CDoubleMatrixSP jtj_inv; // $unregistered
    mutable CDoubleMatrixSP covar; // $unregistered
    mutable CDoubleMatrixSP correl; // $unregistered
    mutable CDoubleArray    lowerLims; // $unregistered
    mutable CDoubleArray    upperLims; // $unregistered
    mutable CDoubleArray    stdError; // $unregistered
};

typedef smartPtr<LevenbergMarquardt> LevenbergMarquardtSP;
typedef smartConstPtr<LevenbergMarquardt> LevenbergMarquardtConstSP;
#ifndef QLIB_OPTIMIZER_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<LevenbergMarquardt>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<LevenbergMarquardt>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<LevenbergMarquardt>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<LevenbergMarquardt>);
#endif

/** Wrapper around imsl's min_uncon_multivar */
class UTIL_DLL QuasiNewton: public OptimizerND{
public:
    static CClassConstSP const TYPE;
    friend class QuasiNewtonHelper;

    virtual void minimize(const MFunctionND&  func,
                          const CDoubleArray& guess,        // initial guess
                          CDoubleArray&       x) const;     // result

    virtual void minimize(const MFunctionND&  func,
                          const CDoubleArray& xguess,        // initial guess
                          const CStringArray& ids,           // identifiers
                          CDoubleArray&       x) const;      // result

    QuasiNewton();
};

typedef smartPtr<QuasiNewton> QuasiNewtonSP;
typedef smartConstPtr<QuasiNewton> QuasiNewtonConstSP;
#ifndef QLIB_OPTIMIZER_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<QuasiNewton>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<QuasiNewton>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<QuasiNewton>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<QuasiNewton>);
#endif


// Forward declaration of Simplex
class Simplex;
typedef smartPtr<Simplex> SimplexSP;
typedef smartConstPtr<Simplex> SimplexConstSP;
#ifndef QLIB_OPTIMIZER_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartConstPtr<Simplex>);
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<Simplex>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartConstPtr<Simplex>);
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<Simplex>);
#endif

/** Wrapper around nr's simplex method */
class UTIL_DLL Simplex: public OptimizerND{
public:
    static CClassConstSP const TYPE;

    // on windows an insane clash with <wingdi.h> where the next two ids are macroses
#if defined (_WINGDI_) && defined( RELATIVE )
#undef RELATIVE
#undef ABSOLUTE
#endif
	
	static const string RELATIVE;
	static const string ABSOLUTE;
	static const string COMPOSITE;
	static const string USE_AMOEBA2;



    friend class SimplexHelper;

    virtual void validatePop2Object();

    virtual void minimize(const MFunctionND&  func,
                          const CDoubleArray& guess,        // initial guess
                          CDoubleArray&       x) const;     // result

    virtual void minimize(const MFunctionND&  func,
                          const CDoubleArray& xguess,        // initial guess
                          const CStringArray& ids,           // identifiers
                          CDoubleArray&       x) const;      // result


	/** mimics public constructor */
    static SimplexSP create(double _lengthScale = 1.0);

    void setAbsTol(double tol);
    void setLengthScale(double _lengthScale);
    void setStoppingCriterionType(string _stoppingCriterionType);
    void useSRMSimplex(string choose);
    void useMyMapping(string _useNewMapping);
    double getLengthScale();

    /** Class that is equipped with a method that returns the coordinates 
        of a basis in a 'n'-dimensional Euclidean space. 
        Will be invoqued to calculate the n+1 vertices of the simplex */
    class UTIL_DLL IBasisMaker: virtual public IObject{
    public:
        static CClassConstSP const TYPE;
        
        /* Returns a square matrix of size n x n
           The i-th 'col' of the returned matrix is the i-th vector
           of the basis*/
        virtual DoubleMatrix make(int n) const = 0;
    };
    typedef smartPtr<IBasisMaker> IBasisMakerSP; 

    /** Returns the coordinates of a basis in a 'n'-dimensional 
        Euclidean space. Will be invoqued to calculate the n+1 
        vertices of the simplex */
    class UTIL_DLL CanonicalBasisMaker: public CObject,
                               public virtual IBasisMaker{
    public:
        static CClassConstSP const TYPE;
        
        CanonicalBasisMaker();

        /* Returns a square matrix of size n x n
           The i-th 'col' of the returned matrix is the i-th vector
           of the basis*/
        virtual DoubleMatrix make(int n) const;


    private:
        static void load(CClassSP& clazz);

        static IObject* defaultCtor();
    };
    
private:
    static double funk(double []);

    void amoeba(double **p, double y[], int ndim, double ftol,
        double (*funk)(double []), int *nfunk) const;

    void amoeba2(double **p, double y[], int ndim, double ftol,double ftolAbs,
                       double (*funk)(double []), int *nfunk, int idx) const;

  
    static double amotry(double **p, double y[], double psum[], int ndim, 
                         double (*funk)(double []), int ihi, int ilo, double fac);

    static void nrerror(const string& text);

    static const double TINY; // A small number.
    static const int NMAX_default;       // Maximum allowed number of function evaluations.
    static const double EXCEPTIONALLY_LARGE_NUMBER;     // returned by obj func upon failure


    Simplex();

    // registered vars
    bool          exitUponFailure;
    int           maxNbIter;   // Maximum allowed number of function evaluations.
    double        ftol;
    double        ftolAbs;
    string        stoppingCriterionType;
    double        lengthScale;
    IBasisMakerSP basisMaker;
    

    // not registered vars
    static Simplex*        me;
    int                    NMAX;
    int                    ndim; // $unregistered
    DoubleArray            x; // $unregistered
    DoubleArray            f; // $unregistered
    const MFunctionND*     func; // $unregistered
    OneToOneMappingConstArray   mappings; // $unregistered
    const CStringArray*    ids; // $unregistered
    int                    nfunk; // $unregistered
    string                 isSRM;
    string                 useNewMapping;
};



/** Interface for 1-dimensional optimizers */
class UTIL_DLL Optimizer1D: public CObject{
public:
    static CClassConstSP const TYPE;

    virtual double minimize(const Function1DDouble& func,
                            double                  xa,
                            double                  xb,
                            double                  xc) const = 0;

    class UTIL_DLL Bracketer{
    public:
        /** Wrapper around numerical recipe's mnbrak */
        static void bracket(const Function1DDouble& func,
                            double&                 xa,
                            double&                 xb,
                            double&                 xc,
                            double&                 fa,
                            double&                 fb,
                            double&                 fc);
    private:
        static double local_func(double x);

        static const Function1DDouble* func;
        static OneToOneMappingConstSP  mapping;
    };

protected:
    Optimizer1D(const CClassConstSP& clazz);

private:
    static void load(CClassSP& clazz);
};

class UTIL_DLL Optimizer1DBrent: public Optimizer1D{
public:
    static CClassConstSP const TYPE;

    virtual double minimize(const Function1DDouble& func,
                            double                  xa,
                            double                  xb,
                            double                  xc) const;

    Optimizer1DBrent(double tol);

private:
    Optimizer1DBrent();

    static void load(CClassSP& clazz);

    static double local_func(double x);

    // registered var
    double tol;
    // transient
    mutable double fmin;
    // not registered vars
    static Optimizer1DBrent* me;
    const Function1DDouble*  func; // $unregistered
    OneToOneMappingConstSP   mapping; // $unregistered
};

DRLIB_END_NAMESPACE

#endif




