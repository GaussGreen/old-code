// Expression.h: interface for the Expression class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_EXPRESSION_H__9DB90283_6ACB_11D2_97E7_00C04FD8EB9A__INCLUDED_)
/**@#-*/
#define AFX_EXPRESSION_H__9DB90283_6ACB_11D2_97E7_00C04FD8EB9A__INCLUDED_
/**@#+*/

#include "kplatdep.h"
#include SSTREAM_H
#include <set>
#include <list>
#include <vector>
#include "kbaseexp.h"


/** Math expression using +, -, *, /

	Parses a math expression into an internal represenation which follows the normal
	infix (not post-fix/reverse polish) rules of math.
	<br>
	To evalue this expression, you must provide values for each symbol in the expression.
	<p>
	Notes:
	<ol>
	<li> Symbols are NOT case-sensitive
	<li> This templated class supports: <b> if </b> 
		Usage:
	<pre><dir>
		if (a, b, c)  
	</dir></pre>
	which returns <b>b</b> if <b>a</b> is TRUE else returns <b>c</b>
	<p>
	Therefore, the templated type MUST be able to cast into bool.
	(This decision was difficult and may be incorrect, but I didn't want to write too many classes.)
*/

template <class Arg>
class KExpression  
{
public:
	typedef Arg (*UnaryFunc) (Arg);
	typedef KMap (upper_string, UnaryFunc) UnaryFunctionMap;
	
	typedef Arg (*BinaryFunc) (Arg, Arg);
	typedef KMap (upper_string, BinaryFunc) BinaryFunctionMap;
	
	/** Constructs an expression from a string.
		@param uf Map of strings to unary function pointers
		@param bf Map of strings to binary function pointers
		*/
	KExpression(upper_string s, const UnaryFunctionMap& uf, const BinaryFunctionMap& bf) 
		: _uf(uf), _bf(bf), _text(s) {_e = Parse(s);}
	/**@#-*/
	virtual ~KExpression() {}
	/**@#+*/
	
	typedef upper_string ExpSymbol;
	typedef KSet(ExpSymbol) SymbolSet;
	typedef KMap(ExpSymbol, Arg) SymbolMap;
	
	/** Returns the set of symbols in the expression */
	const SymbolSet getSymbols() const {return _s;}
	/** Evaluates the expression 

		For the given set of symbol values (provided in the SymbolMap), evaluates the expression 
	*/
	Arg eval (const SymbolMap& m) const {if (_e.null()) throw DRException ("Can't eval null exp"); return (*_e).Eval(m);}
	
	friend ostream& operator<< (ostream& s, const KExpression& e) {s << e._text; return s;}
private:
	upper_string _text;
	SymbolSet _s;
	UnaryFunctionMap _uf;
	BinaryFunctionMap _bf;
	typedef KPtr< BaseExp <Arg> > BaseExpPtr;
	BaseExpPtr _e;
	
	BaseExpPtr Parse (upper_string s)
	{
		ISTRINGSTREAM i(s.c_str()); 
		return Parse(i);
	}
	
	BaseExpPtr Parse (istream& s)
	{
		BaseExpPtr first = GetNextExpression(s);
		if (s.eof() || s.peek() < 0) return first;
		else {
			KVector(BaseExpPtr) le;
			KVector(BinaryOpExp<Arg>::Type) lo;
			le.push_back (first);
			while (!s.eof() && s.peek() > 0) {
				BinaryOpExp<Arg>::Type op = GetBinaryOp(s);
				BaseExpPtr exp = GetNextExpression(s);
				le.push_back(exp);
				lo.push_back(op);
			}			
			return BaseExpPtr (new BinaryOpExp<Arg>(le, lo));			
		}
	}
	
	BinaryOpExp<Arg>::Type GetBinaryOp(istream& s)
	{
		char c = s.get();
		BinaryOpExp<Arg>::Type t;
		if (c == '+') t = BinaryOpExp<Arg>::PLUS;
		else if (c == '-') t = BinaryOpExp<Arg>::MINUS;
		else if (c == '*') t = BinaryOpExp<Arg>::MULT;
		else if (c == '/') t = BinaryOpExp<Arg>::DIV;
		else
			throw DRException ("Invalid binary operator ") << c;
		
		CleanWhite(s);
		return t;
	}
	
	BaseExpPtr GetNextExpression (istream& s)
	{
		BaseExpPtr p;
		char c = s.peek();
		if (c == '-') {
			s.get(); CleanWhite(s);
			p = BaseExpPtr (new UnaryOpExp<Arg> (UnaryOpExp<Arg>::MINUS, GetNextExpression(s)));
		}
		else if (c == '+') {
			s.get(); CleanWhite(s);
			p = BaseExpPtr (new UnaryOpExp<Arg> (UnaryOpExp<Arg>::PLUS, GetNextExpression(s)));
		}
		else if (c == '(') {
			s.get(); CleanWhite(s);
			BaseExpPtr exp = Parse (GetTillCharInSameParenScope(s,')'));
			p = BaseExpPtr (new ParenExp<Arg> (exp));
		}
		else if (IsLetter(c)) {
			upper_string text = GetText (s);
			char c2 = s.peek();
			if (c2 == '(') {
				s.get(); CleanWhite(s);
				if (text == "IF") {
						BaseExpPtr pred = Parse (GetTillCharInSameParenScope(s,','));
						BaseExpPtr argExp1 = Parse (GetTillCharInSameParenScope(s,','));
						BaseExpPtr argExp2 = Parse (GetTillCharInSameParenScope(s,')'));
						p = BaseExpPtr (new IfExp<Arg> (pred, argExp1, argExp2));
				}
				else {
					
					UnaryFunctionMap::const_iterator uIter = _uf.find(text);
					if (uIter != _uf.end()) {
						BaseExpPtr argExp = Parse (GetTillCharInSameParenScope(s,')'));
						p= BaseExpPtr (new UnaryFuncExp<Arg> ((*uIter).second, argExp));
					}
					else {
						BinaryFunctionMap::const_iterator bIter = _bf.find(text);
						if (bIter == _bf.end())
							throw DRException ("Bad function name ") << text;
						
						BaseExpPtr argExp1 = Parse (GetTillCharInSameParenScope(s,','));
						BaseExpPtr argExp2 = Parse (GetTillCharInSameParenScope(s,')'));
						p = BaseExpPtr (new BinaryFuncExp<Arg> ((*bIter).second, argExp1, argExp2));
					}
				}
			}
			else {
				_s.insert(text);
				p = BaseExpPtr (new Symbol<Arg> (text));
			}
		}
		else {
			Arg a;
			s >> a;
			p = BaseExpPtr (new Value<Arg> (a));
		}
		
		CleanWhite(s);
		return p;
	}
	
	upper_string GetTillCharInSameParenScope (istream& s, char t)
	{
		int count = 0;
		upper_string a;
		char c;
		while (!s.eof() && s.peek() > 0) {
			c = s.peek();
			if (c == t && count ==0) break;
			a += s.get();
			if (c == '(') count++;
			if (c == ')') count--;
		}
		if (s.eof() || s.peek() < 0)
			throw DRException ("Failed to find ") << c;
		s.get();
		CleanWhite(s);
		return a;
	}
	
	upper_string GetText(istream& s)
	{
		upper_string ans;
		while (!s.eof() && (IsLetter(s.peek()) || IsNumber(s.peek()))) {
			ans += s.get();
		}
		CleanWhite(s);
		return ans;
	}
	static void CleanWhite (istream& s) 
	{char c; while ((c = s.peek()) == ' ' || c == '\n' || c =='\t' || c == 13) s.get();}
	
	static bool IsLetter (char c) {return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');}
	static bool IsNumber (char c) {return (c >= '0' && c <= '9');}
};

/** Math expression using doubles.

  By restricting ourselves to doubles, we provide support for the following standard unary and binary
  functions (by default).
	<ol>
	<li> <b>Unary Functions:</b>  acos, asin, atan, cos, cosh, exp, fabs, log, log10, sin, sinh, tan, 
	tanh, sqrt, ceil, floor
	<li> <b>Binary Functions:</b> atan2, fmod, max, min, pow
	</ol>
	These functions perform as specified in &ltcmath&gt, except for max and min which have
	the obvious meaning.
*/
class KDoubleExpression : public KExpression<double> {
public:
	/** Constructs the expression */
	KDoubleExpression (upper_string s, const UnaryFunctionMap& uf = UnaryMathFunctions,
		const BinaryFunctionMap& bf = BinaryMathFunctions) : KExpression<double>(s,uf,bf) {}

	/** Standard set of unary math functions */
	static UnaryFunctionMap UnaryMathFunctions;
	/** Standard set of binary math functions */
	static BinaryFunctionMap BinaryMathFunctions;
private:
	static double max (double a, double b);
	static double min (double a, double b);
	static UnaryFunctionMap InitUnary();

	static BinaryFunctionMap InitBinary();
};

#endif // !defined(AFX_EXPRESSION_H__9DB90283_6ACB_11D2_97E7_00C04FD8EB9A__INCLUDED_)
