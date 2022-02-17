#include "edginc/config.hpp"
#include "edginc/FunctionOperations.hpp"

DRLIB_BEGIN_NAMESPACE

/** Utility class for FunctionOperations::shift */
class _MFunctionNDShift:
    public MFunctionND
{
public:
    // Constructor
    _MFunctionNDShift(
        const CDoubleArray& shift,
        const MFunctionND& f)
    : MFunctionND(f.getNbVars(), f.getNbFuncs(), f.getIntervals()),
      m_shift(shift),
      m_f(f),
      m_arg(shift.size())
    {
      // empty
    }
    
    void operator()(const CDoubleArray& x, CDoubleArray& result) const
    {
      for (size_t i = 0; i < x.size(); ++i)
        m_arg[i] = x[i] + m_shift[i]; 
      m_f(m_arg, result);
    }

private:
  // data
  CDoubleArray m_shift;
  const MFunctionND& m_f;

  // workspace
  mutable CDoubleArray m_arg;
};

MFunctionNDSP FunctionOperations::shift(const CDoubleArray& shift, const MFunctionND& f)
{
  return MFunctionNDSP(new _MFunctionNDShift(shift, f));
}

/** Utility class for FunctionOperations::multiply */
class _MFunctionNDMultiplier:
    public virtual MFunctionND
{
public:

    virtual ~_MFunctionNDMultiplier(){}
    
    // Constructor
    _MFunctionNDMultiplier(
        const MFunctionND& f1,
        const MFunctionND& f2,
        int N,
        int M,
        const RangeArray& ranges):
            MFunctionND(N, M, ranges),
            f1(f1),
            f2(f2),
            tmpResult(M)
    {}
    
    virtual void operator()(const CDoubleArray& x, CDoubleArray& result) const
    {
        f1(x, tmpResult);
        f2(x, result);
        for (int i = 0; i < result.size(); ++i) {
			result[i] *= tmpResult[i];
		}
    }
    
private:
    const MFunctionND& f1;
    const MFunctionND& f2;
    mutable CDoubleArray tmpResult;
};

/** Build g(x) = f1(x) * f2(x) defined on the input ranges */
MFunctionNDSP FunctionOperations::multiply(const MFunctionND& f1, const MFunctionND& f2, const RangeArray& ranges)
{
    int N = f1.getNbVars();
    int M = f1.getNbFuncs();

    // sanity checks
    if (N != f2.getNbVars() || M != f2.getNbFuncs())
    {
        throw ModelException("FunctionOperations::multiply",
            "Incompatible functions: functions are not defined on the same spaces.");
    }
    if (ranges.size() != N)
    {
        throw ModelException("FunctionOperations::multiply",
            "Incompatible ranges.");
    }
    
    return MFunctionNDSP(new _MFunctionNDMultiplier(f1, f2, N, M, ranges));
}

/** Utility class for FunctionOperations::add */
class _MFunctionNDAdder:
    public virtual MFunctionND
{
public:

    virtual ~_MFunctionNDAdder(){}
    
    // Constructor
    _MFunctionNDAdder(
        const MFunctionND& f1,
        const MFunctionND& f2,
        int N,
        int M,
        const RangeArray& ranges):
            MFunctionND(N, M, ranges),
            f1(f1),
            f2(f2),
            tmpResult(M)
    {}
    
    virtual void operator()(const CDoubleArray& x, CDoubleArray& result) const
    {
        f1(x, tmpResult);
        f2(x, result);
        for (int i = 0; i < result.size(); ++i) {
            result[i] += tmpResult[i];
        }
    }
    
private:
    const MFunctionND& f1;
    const MFunctionND& f2;
    mutable CDoubleArray tmpResult;
};

/** Build g(x) = f1(x) + f2(x) defined on the input ranges */
MFunctionNDSP FunctionOperations::add(const MFunctionND& f1, const MFunctionND& f2, const RangeArray& ranges)
{
    int N = f1.getNbVars();
    int M = f1.getNbFuncs();
    
    // sanity checks
    if (N != f2.getNbVars() || M != f2.getNbFuncs())
    {
        throw ModelException("FunctionOperations::add",
            "Incompatible functions: functions are not defined on the same spaces.");
    }
    if (ranges.size() != N)
    {
        throw ModelException("FunctionOperations::add",
            "Incompatible ranges.");
    }

    return MFunctionNDSP(new _MFunctionNDAdder(f1, f2, N, M, ranges));
}

/** Utility class for FunctionOperations::project */
class _MFunctionNDProjector:
    public virtual FunctionNDDouble
{
public:

    virtual ~_MFunctionNDProjector(){}
    
    // Constructor
    _MFunctionNDProjector(
        const MFunctionND& refFunction,
        int projectionSpace,
        int N,
        int M,
        const RangeArray& ranges):
            FunctionNDDouble(N, ranges),
            refFunction(refFunction),
            projectionSpace(projectionSpace),
            tmpResult(M)
    {}
    
    virtual double operator()(const CDoubleArray&  x) const
    {
        refFunction(x, tmpResult);
        return tmpResult[projectionSpace];
    }
    
private:
    const MFunctionND& refFunction;
    int projectionSpace;
    mutable CDoubleArray tmpResult;
};

/** Utility class for FunctionOperations::project */
class _MFunction1DProjector:
    public virtual Function1DDouble
{
public:

    virtual ~_MFunction1DProjector(){}

    // Constructor
    _MFunction1DProjector(
        const MFunctionND& refFunction,
        int projectionSpace,
        int M,
        const Range& range):
            Function1DDouble(range),
            refFunction(refFunction),
            projectionSpace(projectionSpace),
            tmpResult(M),
            tmpInput(1)
    {}
    
    virtual double operator()(double x) const
    {
        tmpInput[0] = x;
        refFunction(tmpInput, tmpResult);
        return tmpResult[projectionSpace];
    }
    
private:
    const MFunctionND& refFunction;
    int projectionSpace;
    mutable CDoubleArray tmpResult;
    mutable CDoubleArray tmpInput;
};

/**
 * If f:R^N->R^M such that f(x) = (f_0(x), f_1(x), ..., f_M-1(x)),
 * builds g(x) = f_space(x)
 * 
 * NB: "space" belongs to [0, ..., M-1]
 * */
MFunctionNDSP FunctionOperations::project(const MFunctionND& f, int space)
{
    int N = f.getNbVars();
    int M = f.getNbFuncs();
    
    if (space < 0 || space >= M)
    {
        throw ModelException(
            "FunctionOperations::project",
            "Invalid space: " +
            Format::toString(space) +
            " (valid range is [0-" +
            Format::toString(M-1) +
            "]).");
    }
    
    if (N != 1)
    {
        return MFunctionNDSP(new _MFunctionNDProjector(f, space, N, M, f.getIntervals()));
    }
    else
    {
        return MFunctionNDSP(new _MFunction1DProjector(f, space, M, *(f.getIntervals()[0])));
    }
}

/** Utility class for FunctionOperations::partialFunction */
class _PartialMFunctionND :
    public virtual MFunctionND
{
public:
    // destructor
    virtual ~_PartialMFunctionND() {}

    // constructor
    _PartialMFunctionND(
        const MFunctionND& globalFunction,
        int N,
        int M,
        const RangeArray& ranges,
        int space,
        double y) :
            MFunctionND(N, M, ranges),
            globalFunction(globalFunction),
            value(N+1),
            space(space),
            y(y)
    {
        value[space] = y;
    }

    // function
    virtual void operator()(const CDoubleArray& x, CDoubleArray& f) const
    {
        int i;
        for (i = 0; i < x.size() && i < space; ++i) {
			value[i] = x[i];
		}
        for (i = space; i < x.size(); ++i) {
            value[i+1] = x[i];
        }
        globalFunction(value, f);
    }

private:
    const MFunctionND& globalFunction;
    mutable DoubleArray value;
    int space;
    double y;
};

/**
 * If f:R^N->R^M with f(x) = f(x_0, x_1, ..., x_N-1),
 * builds g:R^(N-1)->R^M such that
 * g(x_1, ...,x_space-1, x_space+1, ...,x_N-1) = f(x_1, ...,x_space-1, y, x_space+1, ...,x_N-1)
 * 
 * NB: "space" belongs to [0, ..., M-1]
 * */
MFunctionNDSP FunctionOperations::partialFunction(const MFunctionND& f, int space, double y)
{
    int N = f.getNbVars();
    int M = f.getNbFuncs();
    
    return MFunctionNDSP(new _PartialMFunctionND(f, N-1, M, f.getIntervals(), space, y));
}

DRLIB_END_NAMESPACE
