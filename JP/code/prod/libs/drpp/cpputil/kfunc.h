// Func.h: interface for the Func class.
//
//////////////////////////////////////////////////////////////////////

/**@#-*/
#if !defined(AFX_FUNC_H__C2FDE8CF_69C3_11D2_97E7_00C04FD8EB9A__INCLUDED_)
#define AFX_FUNC_H__C2FDE8CF_69C3_11D2_97E7_00C04FD8EB9A__INCLUDED_

#include "kinterp.h"
#include "kpolynomial.h"
#include "kptr.h"
#include "kexpression.h"

// Support classes for function.h

template <class Arg, class Result>
class BaseFunc {
public:
	typedef KPtr< BaseFunc<Arg,Result> > BaseFuncPtr;
	BaseFunc() : _mapSize(1) {}
	virtual ~BaseFunc() {}

	Result operator() (const Arg& a) 
	{
		FuncMap::iterator iter = _m.find(a); 
		if (iter != _m.end()) return (*iter).second;
		if (_m.size() == _mapSize) _m.erase(_m.begin());
		return(_m[a]=eval(a));
	}

	virtual Result eval(const Arg&) = 0;
	void setMapSize (int s)  {if (s < 0) throw DRException ("Size must be > 0.  Can't use ") << s; _mapSize = s;}
private:
	typedef KMap(Arg, Result) FuncMap;
	FuncMap _m;
	int _mapSize;
};


template <class Arg, class Result>
class ConstFunc : public BaseFunc<Arg,Result> {
public:
	ConstFunc (const Result& r): _r(r){}
	Result eval (const Arg&) {return _r;}
private:
	Result _r;
};

template <class Arg, class Result>
class TableFunc : public BaseFunc<Arg,Result> {
public:
	typedef pair<Arg,Result> point;
	
	
#ifdef DR_MEMBER_TEMPLATES
	template <class InIt>
		TableFunc (InIt first, InIt last, const Interp<Arg,Result>& i):_i(i)
	{
		while (first != last) {
			point& p = *first;
			_m[p.first] = p.second;
			first++;
		}
	}
#else
	TableFunc (KVector(point)::const_iterator first, KVector(point)::const_iterator last, const KInterp<Arg,Result>& i):_i(i)
	{
		while (first != last) {
			const point& p = *first;
			_m[p.first] = p.second;
			first++;
		}
	}
#endif

	Result eval (const Arg& a)
	{
		TMap::const_iterator iter = _m.upper_bound(a);
		if (iter == _m.begin()) return (*iter).second;
		
		TMap::const_iterator prev = iter;
		prev--;
		if (iter == _m.end()) return (*prev).second;
#ifdef _MSC_VER
		return _i(a, *prev, *iter);
#endif
#ifdef __SUNPRO_CC	// Solaris is such a piece of shit
		point p1 ((*prev).first, (*prev).second);
		point p2 ((*iter).first, (*iter).second);
		return _i(a,p1,p2);
#endif	
	}
private:
	typedef KMap (Arg, Result) TMap;
	TMap _m;
	KInterp<Arg,Result> _i;
};




template <class Arg, class Result>
class AddFunc : public BaseFunc<Arg,Result> {
public:
	AddFunc(const BaseFuncPtr& f, const BaseFuncPtr& g) : _f(f), _g(g) {}
	Result eval (const Arg& a) {return (*_f)(a) + (*_g)(a);}
private:
	BaseFuncPtr _f, _g;
};

template <class Arg, class Result>
class MultFunc : public BaseFunc<Arg,Result> {
public:
	MultFunc(const BaseFuncPtr& f, const BaseFuncPtr& g) : _f(f), _g(g) {}
	Result eval(const Arg& a) {return (*_f)(a) * (*_g)(a);}
private:
	BaseFuncPtr _f, _g;
};

template <class Arg, class Result>
class SubFunc : public BaseFunc<Arg,Result> {
public:
	SubFunc(const BaseFuncPtr& f, const BaseFuncPtr& g) : _f(f), _g(g) {}
	Result eval (const Arg& a)  {return (*_f)(a) - (*_g)(a);}
private:
	BaseFuncPtr _f, _g;
};

template <class Arg, class Result>
class DivFunc : public BaseFunc<Arg,Result> {
public:
	DivFunc(const BaseFuncPtr& f, const BaseFuncPtr& g) : _f(f), _g(g) {}
	Result eval (const Arg& a)  {return (*_f)(a) / (*_g)(a);}
private:
	BaseFuncPtr _f, _g;
};

template <class Arg, class T, class Result>
class CompFunc : public BaseFunc<Arg,Result> {
public:
	typedef KPtr< BaseFunc<T,Result> > BaseFuncPtr1;
	typedef KPtr< BaseFunc<Arg,T> > BaseFuncPtr2;
	
	CompFunc(const BaseFuncPtr1& f, const BaseFuncPtr2& g) : _f(f), _g(g) {}
	Result eval (const Arg& a)  {return (*_f)((*_g)(a));}
private:
	BaseFuncPtr1 _f;
	BaseFuncPtr2 _g;
};

template <class Arg>
class PolyFunc : public BaseFunc<Arg, Arg> {
public:
	PolyFunc(const KPolynomial<Arg>& p): _p(p) {}
	Arg eval (const Arg& a)  {return _p(a);}
private:
	KPolynomial<Arg> _p;
};

template <class Arg>
class ExpressionFunc : public BaseFunc<Arg,Arg> {
public:
	ExpressionFunc (const string& s, const KExpression<Arg>::UnaryFunctionMap& uf, 
		const KExpression<Arg>::BinaryFunctionMap& bf): _e(s, uf, bf)
	{if (_e.GetSymbols().size() > 1) 
	throw DRException ("Too many symbols in expression");}
	
	Arg eval (const Arg& a)  
	{const Expression<Arg>::SymbolSet& symbols = _e.GetSymbols();
	Expression<Arg>::SymbolMap s;
	if (symbols.size() != 0)  s[*(symbols.begin())] = a;
	return _e.Eval(s);
	}
private:
	KExpression<Arg> _e;
};

class DExpressionFunc : public BaseFunc<double,double> {
public:
	DExpressionFunc (const string& s, const KDoubleExpression::UnaryFunctionMap& uf, 
		const KDoubleExpression::BinaryFunctionMap& bf): _e(s, uf, bf)
	{if (_e.getSymbols().size() > 1) 
	throw DRException ("Too many symbols in expression");}
	
	double eval (const double& a)  
	{const KDoubleExpression::SymbolSet& symbols = _e.getSymbols();
	KDoubleExpression::SymbolMap s;
	if (symbols.size() != 0)  s[*(symbols.begin())] = a;
	return _e.eval(s);
	}
private:
	KDoubleExpression _e;
};


#endif // !defined(AFX_FUNC_H__C2FDE8CF_69C3_11D2_97E7_00C04FD8EB9A__INCLUDED_)
/**@#+*/
