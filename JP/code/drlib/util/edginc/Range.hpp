//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Range.hpp
//
//   Description : 
//
//   Date        : 21 May 02
//
//
//----------------------------------------------------------------------------

#include "edginc/Object.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Array.hpp"
#include "edginc/DECLARE.hpp"

#ifndef EDR_RANGE_HPP
#define EDR_RANGE_HPP

DRLIB_BEGIN_NAMESPACE

/** Boundary classes */
class UTIL_DLL Boundary : public CObject
{
public:
	static CClassConstSP const TYPE;

	/** Called immediately after object constructed */
	virtual void validatePop2Object();  

    friend class Range;

    Boundary& operator=(const Boundary& rhs);

    bool isClosedBracket() const{
        return isclosed;
    }

    bool isInfinite() const{
        return isinf;
    }

    double getValue() const{
        if(isinf){
            throw ModelException("Boundary::getValue",
                                 "Cannot get value as boundary is infinite");
        }
        return value;
    }
//protected:
    Boundary(const Boundary& rhs);
private:

	// For reflection
	static void load (CClassSP& clazz);
	/** private constructor */
	Boundary();
	/** default 'constructor' */
	static IObject* defaultBoundary();


    bool isGreater(const Boundary& rhs) const;
    bool isLess(const Boundary& rhs) const;
    Boundary subtract(double var) const;
    string toString() const;

protected:
    explicit Boundary(double value, bool isinf, bool isclosed);

    double value;
    bool isinf;     // NB: if +infinity, value = +1.0; if -infinity, value = -1.0
    bool isclosed;  // NB: if +/-infinity, isclosed=false
};

DECLARE(Boundary);

class UTIL_DLL ClosedBoundary: public Boundary{
public:
    explicit ClosedBoundary(double value);
    ClosedBoundary(const ClosedBoundary& rhs);

private:
    ClosedBoundary& operator=(const ClosedBoundary& rhs);
};

class UTIL_DLL OpenBoundary: public Boundary{
public:
    explicit OpenBoundary(double value);
    OpenBoundary(const OpenBoundary& rhs);

protected:
    explicit OpenBoundary(double value, bool isinf);

private:
    OpenBoundary& operator=(const OpenBoundary& rhs);
};

class UTIL_DLL Infinity: public OpenBoundary{
public:
    enum Sign {
        Minus = -1,
        Plus = 1
    };

    explicit Infinity(Sign sign);

private:
    Infinity(const OpenBoundary& rhs);
    Infinity& operator=(const OpenBoundary& rhs);
};

class UTIL_DLL Range : public CObject
{
public:
	// For reflection
	static CClassConstSP const TYPE;  

    Range(const Boundary& lower, const Boundary& upper);

    Range(double lower, bool lowerClosed,
          double upper, bool upperClosed);

    Range(const Range& rhs);

    Range& operator=(const Range& rhs);

    //// needed for RangeArray on VC71 and dlls
    //Range();

    const Boundary& getLower() const{
        return l;
    }

    const Boundary& getUpper() const{
        return r;
    }

    bool isNonEmpty() const;
    static void checkIsNonEmpty(const Range& interval);

    bool isSingleton() const;
    static void checkIsNotSingleton(const Range& interval);

    bool isOpen() const;
    static void checkIsOpen(const Range& interval);

    static bool variableIsInRange(const Range& range, double variable);
    static void checkVariableIsInRange(const Range& range, 
                                       double variable, 
                                       string variableName = "Variable");

    Range subtract(double var) const;
    static Range subtract(const Range& range, double var);

    string toString() const;

    // convert a possibly closed range into an open range
    void makeOpen();

    //// need this for refCountPtr
    ~Range();
private:
	// For reflection
	static void load (CClassSP& clazz);
	/** private constructor */
	Range();
	/** default 'constructor' */
	static IObject* defaultRange();

    Boundary l; // left / lower bound
    Boundary r; // right / upper bound
};


DECLARE(Range);

class UTIL_DLL InfiniteRange: public Range{
public:
    InfiniteRange();
    static RangeArraySP createInfiniteRangeArray(int size);
private:
    InfiniteRange& operator=(const Range& rhs);
    InfiniteRange(const Range& rhs);
};

typedef smartPtr<InfiniteRange> InfiniteRangeSP;

DRLIB_END_NAMESPACE

#endif
