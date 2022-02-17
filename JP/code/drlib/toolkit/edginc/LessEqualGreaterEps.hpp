//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : LessEqualGreaterEps.hpp
//
//   Description : Comparisons with epsilon (un)tolerance as the objects
//
//   Date        : March 2006
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_LESSEQUALGREATEREPS_HPP
#define QLIB_LESSEQUALGREATEREPS_HPP

DRLIB_BEGIN_NAMESPACE
///////////////////////////////////////////////////////////////////////////////
class TOOLKIT_DLL strictLessEps : public binary_function<double, double, bool>
{
public:
    strictLessEps(double gEpsilon);
    ~strictLessEps(){;};
	bool operator()(double g1, double g2) const;
	const double EPSILON;
private:
    void operator=(strictLessEps & g2Clone);     // not defined
    strictLessEps(strictLessEps & g2Clone);      // not defined
};
///////////////////////////////////////////////////////////////////////////////
template <class CType> 
class TOOLKIT_DLL lessOrEqual : public binary_function<CType, CType, bool>
{
public:
    lessOrEqual(){;};
    ~lessOrEqual(){;};
    bool operator()(const CType & g1, const CType & g2) const
    {
        return (!(g2 > g1));
    }
private:
    void operator=(lessOrEqual & g2Clone);     // not defined
    lessOrEqual(lessOrEqual & g2Clone);        // not defined
};
///////////////////////////////////////////////////////////////////////////////
class TOOLKIT_DLL looseLessEps : public binary_function<double, double, bool>
{
public:
	looseLessEps(double gEpsilon);
    ~looseLessEps(){;};
	bool operator()(double g1, double g2) const;
	const double EPSILON;
private:
    void operator=(looseLessEps & g2Clone);     // not defined
    looseLessEps(looseLessEps & g2Clone);       // not defined
};
///////////////////////////////////////////////////////////////////////////////
class TOOLKIT_DLL equalEps : public binary_function<double, double, bool>
{
public:
	equalEps(double gEpsilon);
    ~equalEps(){;};
	bool operator()(double g1, double g2) const;
	const double EPSILON;
private:
    void operator=(equalEps & g2Clone);     // not defined
    equalEps(equalEps & g2Clone);           // not defined
};
///////////////////////////////////////////////////////////////////////////////
class TOOLKIT_DLL looseGreaterEps: public binary_function<double, double, bool>
{
public:
	looseGreaterEps(double gEpsilon);
    ~looseGreaterEps(){;};
    bool operator()(double g1, double g2) const;
	const double EPSILON;
private:
    void operator=(looseGreaterEps & g2Clone);     // not defined
    looseGreaterEps(looseGreaterEps & g2Clone);    // not defined
};
///////////////////////////////////////////////////////////////////////////////
template <class CType> 
class TOOLKIT_DLL greaterOrEqual : public binary_function<CType, CType, bool>
{
public:
	greaterOrEqual(){;};
    greaterOrEqual(greaterOrEqual & g2Clone){;};
    ~greaterOrEqual(){;};
	bool operator()(CType g1, CType g2) const
    {
        return (!(g1<g2));
    };
private:
    void operator=(greaterOrEqual & g2Clone);     // not defined
};
///////////////////////////////////////////////////////////////////////////////
class TOOLKIT_DLL strictGreaterEps:public binary_function<double, double, bool>
{
public:
	strictGreaterEps(double gEpsilon);
    ~strictGreaterEps(){;};
	bool operator()(double g1, double g2) const;
	const double EPSILON;
private:
    void operator=(strictGreaterEps & g2Clone);   // not defined
    strictGreaterEps(strictGreaterEps & g2Clone); // not defined
};
///////////////////////////////////////////////////////////////////////////////
template <class ElemType> class TOOLKIT_DLL CompareNames : 
                               public binary_function<ElemType, ElemType, bool>
{
public:
    CompareNames(){;};
    ~CompareNames(){;};
	bool operator()(const ElemType & g1, const ElemType & g2) const
    {
        return (g1.getName() < g2.getName());
    };
private:
    void operator=(CompareNames & g2Clone);  // not defined
    CompareNames(CompareNames & g2Clone);    // not defined
};
///////////////////////////////////////////////////////////////////////////////
DRLIB_END_NAMESPACE
#endif // QLIB_LESSEQUALGREATEREPS_HPP
