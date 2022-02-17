// Function.h: interface for the Function class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FUNCTION_H__E79F7879_6A23_11D2_865A_000000000000__INCLUDED_)
/**@#-*/
#define AFX_FUNCTION_H__E79F7879_6A23_11D2_865A_000000000000__INCLUDED_
/**@#+*/

#include "kfunc.h"

/** Math Function

  Provides a function from Arg to Results.
  <br>
  This function can be initialized in two different ways.
  <ol>
  <li> Vector of (arg,result) pairs with an interpolation type.  Extrapolation is NOT done on points
  beyond the endpoints, but rather the endpoint values are returned.
  <li> Constant functions.
  </ol>

  After creating a function, the usual binary operations +, -, *, / and composition can be performed.

  <p>
  For efficiency reasons, function (arg, value) pairs are memoized (cached).  The default setting is to only
  memoize the last pair.  However, you can adjust this by using setCacheSize.  
  If there is a cache overflow, the (argument, value) with the least ordered argument is removed.
  @see KArgFunction
  @see KDoubleFunction
  @see compose
  */
template <class Arg, class Result>
class KFunction
{
public:
	typedef std::pair<Arg,Result> point;
	typedef KPtr< BaseFunc<Arg,Result> > BaseFuncPtr;
	
	/** Create an empty function that returns the default value fo Result */
	KFunction() {}

#ifdef K_MEMBER_TEMPLATES
	template <class InIt>
		KFunction(InIt first, InIt last, const KInterp<Arg,Result>& i) {_fp = new	TableFunc<Arg,Result>(first,last, i);}
#else 	
	/** Creates a function from a vector of points and an interpolation type */
	KFunction(KVector(point)::const_iterator first, KVector(point)::const_iterator last, const KInterp<Arg,Result>& i) {_fp = new	TableFunc<Arg,Result>(first,last, i);}
#endif

	/** Creates a constant function */
	KFunction (const Result& r) {_fp = new ConstFunc<Arg,Result> (r);}

	/**@#-*/
	KFunction (const BaseFuncPtr& p):_fp(p){}
	virtual ~KFunction() {}
	/**@#+*/
	
	/** Evaluates the function */
	Result operator()(const Arg& a) {if (_fp.null()) return Result(); else return (*_fp)(a);}
	
	/** Adds two functions */
	KFunction operator+(const KFunction& f) const {return	KFunction(BaseFuncPtr(new AddFunc<Arg,Result>(*this, f)));}
	/** Subtracts two functions */
	KFunction operator-(const KFunction& f) const {return	KFunction(BaseFuncPtr(new SubFunc<Arg,Result>(*this, f)));}
	/** Divides two functions */
	KFunction operator/(const KFunction& f) const {return	KFunction(BaseFuncPtr(new DivFunc<Arg,Result>(*this, f)));}
	/** Multiplies two functions */
	KFunction operator*(const KFunction& f) const {return	KFunction(BaseFuncPtr(new MultFunc<Arg,Result>(*this, f)));}
	
	/** Adds a constant to a function */
	KFunction operator+(const Result& r) const {return (*this) + KFunction(r);}
	/** Subtracts a constant from a function */
	KFunction operator-(const Result& r) const {return (*this) -	KFunction(r);}
	/** Divides a function by a constant */
	KFunction operator/(const Result& r) const {return (*this) /	KFunction(r);}
	/** Multiplies a constant to a function */
	KFunction operator*(const Result& r) const {return (*this) *	KFunction(r);}
	
	/** Adds a constant to a function */
	friend KFunction operator+(const Result& r, const KFunction& f) {return KFunction(r) + f;}
	/** Subtracts a function from a constant*/
	friend KFunction operator-(const Result& r, const KFunction& f) {return KFunction(r) - f;}
	/** Divides a constant by a function */
	friend KFunction operator*(const Result& r, const KFunction& f) {return KFunction(r) * f;}
	/** Multiplies a constant to a function */
	friend KFunction operator/(const Result& r, const KFunction& f) {return KFunction(r) / f;}


	void setCacheSize (int s) {_fp.setMapSize(s);}
	/**@-*/
	operator const BaseFuncPtr() const {return _fp;}
	/**@+*/

protected:
	BaseFuncPtr _fp;
};

/** Function composition */
template <class Arg, class T, class Result>
KFunction<Arg, Result> compose (const KFunction<Arg, T>& f, const KFunction<T, Result>& g)
{
	KFunction<Arg,Result>::BaseFuncPtr temp = new CompFunc<Arg, T, Result> (f, g);
	return KFunction<Arg,Result>::KFunction (temp);
}


/** Math function from Arg to Arg */
template <class Arg>
class KArgFunction : public KFunction<Arg, Arg>
{
public:
	/** Empty function that returns Arg() */
	KArgFunction () {}
	/** Construct from KFunction */
	KArgFunction (const KFunction<Arg,Arg>& f) : KFunction<Arg,Arg>(f) {}

#ifdef K_MEMBER_TEMPLATES
	template <class InIt>
		KArgFunction(InIt first, InIt last, const KInterp<Arg,Arg>& i): Function<Arg,Arg>(first, last, i) {}
#else
	/** Creates a function from a vector of points and an interpolation type */
	KArgFunction(KVector(point)::const_iterator first, KVector(point)::const_iterator last, const KInterp<Arg,Arg>& i): KFunction<Arg,Arg>(first, last, i) {}
#endif

	/** Creates a constant function */
	KArgFunction (const Arg& r) : KFunction<Arg, Arg>(r) {}
	
	/** Creates a function from an expression and unary/binary function pointers */
	KArgFunction (KString s, const KExpression<Arg>::UnaryFunctionMap& uf, 
		const KExpression<Arg>::BinaryFunctionMap& bf) {_fp = new ExpressionFunc<Arg> (s, uf, bf);}
	/** Creates a function from a polynomial */
	KArgFunction (const KPolynomial<Arg>& p) {_fp = new PolyFunc<Arg>(p);}
};

/** Math function from double to double */
class KDoubleFunction : public KArgFunction<double> {
public:
	/** Empty function that returns 0 */
	KDoubleFunction () {}

#ifdef K_MEMBER_TEMPLATES
	template <class InIt>
		KDoubleFunction(InIt first, InIt last, const KInterp<double,double>& i): ArgFunction<double>(first, last, i) {}
#else
	/** Creates a function from a vector of points and an interpolation type */
	KDoubleFunction(KVector(point)::const_iterator first, KVector(point)::const_iterator last, const KInterp<double,double>& i): KArgFunction<double>(first, last, i) {}
#endif

	/** Creates a constant function */
	KDoubleFunction (double r) : KArgFunction<double>(r) {}
	/** Creates a function from a polynomial */
	KDoubleFunction (const KPolynomial<double>& p) : KArgFunction<double>(p) {}
	/** Creates a function from an expression and unary/binary function pointers */
	KDoubleFunction (KString s, const KDoubleExpression::UnaryFunctionMap& uf= KDoubleExpression::UnaryMathFunctions, 
		const KDoubleExpression::BinaryFunctionMap& bf = KDoubleExpression::BinaryMathFunctions) {_fp = new DExpressionFunc (s, uf, bf);}
	/** Construct from KArgFunction */
	KDoubleFunction (const KFunction<double, double>& f) : KArgFunction<double> (f) {}	
};

/** Simple Binary Function for Bing */
template <class Arg1, class Arg2, class Result>
class KBinaryFunction {
public:
	KBinaryFunction () {}
	Result operator() (const Arg1& a1, const Arg2& a2) {return eval(a1,a2);}
protected:
	virtual Result eval(const Arg1&, const Arg2&) = 0;
};

template <class T>
class KPlus:public KBinaryFunction<T, T, T>
{
protected:
	T eval(const T&a, const T&b){return a+b;}
};
template <class T>
class KMinus:public KBinaryFunction<T, T, T>
{
protected:
	T eval(const T&a, const T&b){return a-b;}
};

template <class T>
class KMultiplies:public KBinaryFunction<T, T, T>
{
protected:
	T eval(const T&a, const T&b){return a*b;}
};
template <class T>
class KDivides:public KBinaryFunction<T, T, T>
{
protected:
	T eval(const T&a, const T&b){return a/b;}
};
/** Simple Double Binary Function for Bing */
class KDoubleBinaryFunction : public KBinaryFunction<double, double, double> {
public:
	KDoubleBinaryFunction (): _e("0") {}
	KDoubleBinaryFunction (KString s, const KDoubleExpression::UnaryFunctionMap& uf= KDoubleExpression::UnaryMathFunctions, 
		const KDoubleExpression::BinaryFunctionMap& bf = KDoubleExpression::BinaryMathFunctions) : _e(s,uf,bf)
	{
			{if (_e.getSymbols().size() > 2) 
	throw KException ("Too many symbols in expression");}
	}

protected:
	double eval (const double& a, const double& b)  
	{const KDoubleExpression::SymbolSet& symbols = _e.getSymbols();
	KDoubleExpression::SymbolMap s;
	if (symbols.size() > 0)  s[*(symbols.begin())] = a;
	if (symbols.size() > 1)  s[*(++(symbols.begin()))] = b;
	return _e.eval(s);
	}
private:
	KDoubleExpression _e;
};


#endif //!defined(AFX_FUNCTION_H__E79F7879_6A23_11D2_865A_000000000000__INCLUDED_)


