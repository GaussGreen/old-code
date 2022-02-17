// drfunction.cpp: implementation of the drfunction class.
//
//////////////////////////////////////////////////////////////////////

#include "drutils.h"
#include "drvalarray.h"
#include "drfunction.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

DRFunction::DRFunction(DArray& xvals, DArray& yvals, FuncInterpType interp)
{
	if (xvals.size() != yvals.size()) 
		throw DRException("Table arrays must be the same size");

	for (int i = 0; i < xvals.size(); i++)
		operator[](xvals[i]) = yvals[i];

	if (interp == LINEAR) {
		Func = LinearInterp;
	}
	else {
		Func = GeometricInterp;
	}
}


double DRFunction::operator() (double x)
{
	iterator upper = upper_bound (x);	

	if (upper == begin()) {
		return (*begin()).second;
	}
	else if (upper == end()) {
		iterator temp = end();
		temp--;
		return (*temp).second;
	}
	else {
		iterator lower = upper;
		lower--;
		double ans = Func (x, (*lower).first, (*upper).first,
			(*lower).second, (*upper).second);
		return ans;
	}
}

ostream& operator<<(ostream& s, const DRFunction& a)
{
	for (DRFunction::const_iterator iter = a.begin() ; iter != a.end(); iter++)
		s << (*iter).first << "\t" << (*iter).second << endl;

	return s;
}

