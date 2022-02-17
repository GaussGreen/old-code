// Interp.h: interface for the Interp class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_INTERP_H__C2FDE8CD_69C3_11D2_97E7_00C04FD8EB9A__INCLUDED_)
/**@#-*/
#define AFX_INTERP_H__C2FDE8CD_69C3_11D2_97E7_00C04FD8EB9A__INCLUDED_
/**@#+*/


#include "kplatdep.h"
#include <map>

/** Interpolation class

  This class provides interpolation of values.  It supports two types of interpolation:
  <ol>
  <li> Step.  For arguments between the end points, the right most value is used (ie. if f(0) = 0 
	and f(1) = 1, then f(1/2) = 1)
  <li> Linear.  Fit the two points to a straight line and evaluate the argument along this line.
  </ol>

  The way to get the interpolation types is slightly usual.  To get an interpolation type, call
  one of the static member functions: step or linear.
  <p>
  The reason for this oddity is the lack of template support in SunCC-4.2.
  */
template <class Arg, class Result>
class KInterp {
public:
	KInterp(){}
	/** Gets a left continuous step interpolation */
	static KInterp l_step() {return KInterp(L_STEP);}
	/** Gets a right continuous step interpolation */
	static KInterp r_step() {return KInterp(R_STEP);}
	/** Gets a linear interpolation */
	static KInterp linear () {return KInterp(LINEAR_INTERP);}

	/** Performs interpolation 
		
		  A end point is a pair (arg, value).  For step interpolations, we assume arg of p1 &lt 
		  arg of p2.

		  @param a  Evaluation argument
		  @param p1 First end point.
		  @param p2 Second end point.
	*/
	Result operator()(const Arg& a, const std::pair<Arg,Result>& p1, const std::pair<Arg,Result>& p2) const
	{
		if (_t == L_STEP)
			return p1.second;
		else if(_t == R_STEP)
			return p2.second;
	
		Result ans = p1.second;
		ans += (a - p1.first)/(p2.first-p1.first)*(p2.second-p1.second);
		return ans;
	}
	
private:
	enum Type { L_STEP,R_STEP, LINEAR_INTERP };
	KInterp(Type t) : _t(t) {}
	Type _t;
};



#endif // !defined(AFX_INTERP_H__C2FDE8CD_69C3_11D2_97E7_00C04FD8EB9A__INCLUDED_)
