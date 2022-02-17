// drcurvearray.cpp: implementation of the drcurvearray class.
//
//////////////////////////////////////////////////////////////////////

#include "drcurvearray.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

DRCurveArray::operator TCurve**() 
{
	m_tcurves.resize(size());
	for (int i = 0; i < size(); i++) {
		m_tcurves[i] = (TCurve*) operator[](i);
	}
	return (TCurve**) m_tcurves;
}

ostream& operator<< (ostream& s, const DRCurveArray& a)
{
	for (int i = 0 ; i < a.size(); i++) {
		s << "Curve " << i  << endl;
		s << a.operator[](i);
	}

	return s;
}
