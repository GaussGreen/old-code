// BaseExp.h: interface for the BaseExp class.
//
//////////////////////////////////////////////////////////////////////

/**@#-*/

#if !defined(AFX_BASEEXP_H__03DCC626_6AF6_11D2_865A_000000000000__INCLUDED_)
#define AFX_BASEEXP_H__03DCC626_6AF6_11D2_865A_000000000000__INCLUDED_

#include "kplatdep.h"
#include <vector>
#include <map>
#include "kptr.h"
#include "kexception.h"
#include "kstring.h"

//  These classes support the expression class.  In a language with packages, they would be private parts of
//  the expression package.


template <class Arg>
class BaseExp {
public:
	typedef KMap(upper_string, Arg) SymbolMap;
	typedef KPtr< BaseExp<Arg> > BaseExpPtr;
	virtual Arg Eval (const SymbolMap&) = 0;
	virtual ~BaseExp() {}
};

template <class Arg>
class Symbol : public BaseExp<Arg> {
public:
	typedef BaseExp<Arg> super;
	Symbol (const upper_string& s) : _s(s) {}
	virtual Arg Eval (const super::SymbolMap& m)
	{
		super::SymbolMap::const_iterator iter = m.find(_s);
		if (iter == m.end()) throw KException ("Undefined symbol ") << _s;
		return (*iter).second;
	}
private:
	upper_string _s;
};

template <class Arg>
class IfExp : public BaseExp<Arg> {
public:
	typedef BaseExp<Arg> super;
	IfExp (const BaseExpPtr& p, const BaseExpPtr& r1, const BaseExpPtr& r2) : _p(p), _r1(r1), _r2(r2) {}
	virtual Arg Eval (const super::SymbolMap& m) {if ((*_p).Eval(m)) return (*_r1).Eval(m); else return (*_r2).Eval(m);}
private:
	BaseExpPtr _p, _r1, _r2;
};

template <class Arg>
class Value : public BaseExp<Arg> {
public:
	typedef BaseExp<Arg> super;
	Value (const Arg& a) : _a(a) {}
	virtual Arg Eval (const super::SymbolMap&) {return _a;}
private:
	Arg _a;
};

template <class Arg>
class ParenExp : public BaseExp<Arg> {
public:
	typedef BaseExp<Arg> super;
	ParenExp (const BaseExpPtr& p) : _p(p) {}
	virtual Arg Eval (const super::SymbolMap& m) {return (*_p).Eval(m);}
private:
	BaseExpPtr _p;
};

template <class Arg>
class UnaryFuncExp : public BaseExp<Arg> {
public:
	typedef BaseExp<Arg> super;
	typedef Arg (*UnaryFunc) (Arg);
	UnaryFuncExp (const UnaryFunc& f, const BaseExpPtr& p) : _f(f), _p(p) {}
	virtual Arg Eval (const super::SymbolMap& m) {return _f((*_p).Eval(m));}
private:
	BaseExpPtr _p;
	UnaryFunc _f;
};

template <class Arg>
class BinaryFuncExp : public BaseExp<Arg> {
public:
	typedef BaseExp<Arg> super;
	typedef Arg (*BinaryFunc) (Arg, Arg);
	BinaryFuncExp (const BinaryFunc& f, const BaseExpPtr& p1, const BaseExpPtr& p2) : _f(f), _p1(p1), _p2(p2) {}
	virtual Arg Eval (const super::SymbolMap& m) {return _f((*_p1).Eval(m), (*_p2).Eval(m));}
private:
	BaseExpPtr _p1, _p2;
	BinaryFunc _f;
};

template <class Arg>
class UnaryOpExp : public BaseExp<Arg> {
public:
	typedef BaseExp<Arg> super;
	enum Type { PLUS, MINUS } ;
	UnaryOpExp (Type t, const BaseExpPtr& p) : _t(t), _p(p) {}
	virtual Arg Eval (const super::SymbolMap& m) {return (_t==PLUS) ? +(*_p).Eval(m) : - (*_p).Eval(m);}
private:
	BaseExpPtr _p;
	Type _t;
};

template <class Arg>
class BinaryOpExp : public BaseExp<Arg> {
public:
	typedef BaseExp<Arg> super;
	enum Type { PLUS, MINUS, MULT, DIV } ;
	BinaryOpExp (const KVector(BaseExpPtr)& le, const KVector(Type)& lo) : _le(le), _lo(lo) {_lv.resize(_le.size());}
	virtual Arg Eval (const super::SymbolMap& m) 
	{
		int i;

		for (i = 0; i < _lv.size(); i++) {
			_lv[i] = _le[i]->Eval(m);
		}
		
		int hold = 0;
		for (i = 0; i < _lo.size(); i++) {
			if (_lo[i] == MULT) {
				_lv[hold] *= _lv[i+1];
			}
			else if (_lo[i] == DIV) {
				_lv[hold] /= _lv[i+1];
			}
			else {
				hold = i+1;
			}
		}
		for (i = 0; i < _lo.size(); i++) {
			if (_lo[i] == PLUS) {
				_lv[0] += _lv[i+1];
			}
			else if (_lo[i] == MINUS) {
				_lv[0] -= _lv[i+1];
			}
		}
		return _lv[0];						
	}

		

	
	/*		DRList(Arg) la;
		DRList(Type) lo = _lo;
		for (list<BaseExpPtr>::const_iterator iter = _le.begin(); iter != _le.end(); iter++) {
			const BaseExpPtr& exp = *iter;
			Arg a = (*exp).Eval(m);
			la.push_back (a);
		}
		DRList(Arg)::iterator argIter = la.begin();
		DRList(Type)::const_iterator opIter = lo.begin();
		while (opIter != lo.end()) {
			if (*opIter == MULT) {
				UseOperator (la, lo, argIter, opIter, Mult);
			}
			else if (*opIter == DIV) {
				UseOperator (la, lo, argIter, opIter, Div);
			}
			else {
				opIter++;
				argIter++;
			}
		}
		argIter = la.begin();
		opIter = lo.begin();
		while (opIter != lo.end()) {
			if (*opIter == PLUS) {
				UseOperator (la, lo, argIter, opIter, Plus);
			}
			else if (*opIter == MINUS) {
				UseOperator (la, lo, argIter, opIter, Minus);
			}
			else {
				throw DRException ("Screw up in Binary Op Exp");
			}
		}
		return la.front();
	}
*/
private:
	static Arg Plus (Arg a, Arg b) {return a+b;}
	static Arg Minus (Arg a, Arg b) {return a-b;}
	static Arg Mult (Arg a, Arg b) {return a*b;}
	static Arg Div (Arg a, Arg b) {return a/b;}

	typedef Arg (*Op) (Arg, Arg);
//	void UseOperator (DRList(Arg)& la, DRList(Type)& lo, DRList(Arg)::iterator& aIter, DRList(Type)::const_iterator& opIter, Op op)
//	{
/*		list <Arg>::iterator aIter2 = aIter;
		aIter2++;
		(*aIter) = op((*aIter), (*aIter2));
		la.erase(aIter2);
		opIter = lo.erase(opIter);
		*/
//	}
	KVector(BaseExpPtr) _le;
	KVector(Type) _lo;
	KVector(double) _lv;
};

#endif // !defined(AFX_BASEEXP_H__03DCC626_6AF6_11D2_865A_000000000000__INCLUDED_)

/**@#+*/
