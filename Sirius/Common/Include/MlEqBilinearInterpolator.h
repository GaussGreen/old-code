//	MlEqBilinearInterpolator.h:	Stores values of a function f(t_s, t_e).
//								Bilinear interpolation is implemented.
//
//	author:						David Cuin
//
/////////////////////////////////////////////////////////////////////

#ifndef _MLEQBILINEARINTERPOLATOR_H
#define _MLEQBILINEARINTERPOLATOR_H

#pragma once
#include <vector>
#include "smart.h"

template <class ARG_TYPE, class VAL_TYPE>				// e.g. set ARG_TYPE to a long (to represent a date) and VAL_TYPE to double.
class MlEqBilinearInterpolator: public RCObject
{	
public:
	bool								GetValue(ARG_TYPE t_s, ARG_TYPE t_e, bool bInterpolate, bool bExtrapolate, VAL_TYPE* pvt);
	void								AddValue(ARG_TYPE t_s, ARG_TYPE t_e, VAL_TYPE vt);

protected:
	class E
	{
	public:
		ARG_TYPE						at;
		VAL_TYPE						vt;		
	};

	class S
	{
	public:
		ARG_TYPE						at;
		std::vector<E>					v;		
	};

	std::vector<S>						m_v;

	//	Returns pn1 and pn2 (pn1 <= pn2) which are the elements in an input array closest to an input ARG_TYPE.
	//	The parameters bInterpolate and bExtrapolate only affect the return value - pn1 and pn2 are always set.
	//	If no data are found then pn1 and pn2 are set to -1.
	template <class T> bool				GetElement(const std::vector<T>& v, ARG_TYPE at, bool bInterpolate, bool bExtrapolate, long* pn1, long* pn2) const
	{		
		long							nLow, nHigh, nMid;
		nLow = 0;
		nHigh = v.size() - 1;
		while (nLow <= nHigh){
			nMid = (nLow + nHigh) / 2;
			if (at == v[nMid].at){
				// direct hit
				*pn1 = *pn2 = nMid;
				return true;
			} else if (nMid + 1 < v.size() && at > v[nMid].at && at < v[nMid + 1].at){
				// interpolation is required
				*pn1 = nMid;
				*pn2 = nMid + 1;
				return bInterpolate;
			} else if (at > v[nMid].at){
				nLow = nMid + 1;
			} else {
				nHigh = nMid - 1;
			}
		}

		if (!v.size()){
			*pn1 = *pn2 = - 1;
			return false;
		} else if (v[0].at == at){
			// direct hit on the first element (not picked up earlier for single v case)
			*pn1 = *pn2 = 0;
			return true;			
		} else if (v[0].at < at){
			*pn1 = *pn2 = v.size() - 1;
			return bExtrapolate;
		} else {
			*pn1 = *pn2 = 0;
			return bExtrapolate;
		}
	}


	template <class T> 
	std::vector<T>::iterator			GetElementIterator(std::vector<T>& v, long n) const
	{
		std::vector<T>::iterator it = v.begin();
		for (; it != v.end() && n ; it++, n--);
		return it;
	}
};

//	Add a value to the interpolator schedule.
template <class ARG_TYPE, class VAL_TYPE>
void MlEqBilinearInterpolator<ARG_TYPE, VAL_TYPE>::AddValue(ARG_TYPE t_s, ARG_TYPE t_e, VAL_TYPE vt)
{	
	long								n_l, n_h;

	// Find a vector in m_v corresponding to t_s. We then assign a reference variable, s, to this (for clarity).
	if (GetElement(m_v, t_s, false, false, &n_l, &n_h)){
		// Direct hit - m_v[n_l] is the S value that we need to modify.		
	} else if (n_l != n_h){
		// Interpolation case. We need to insert a new S just before n_h.
		n_l = n_h;
		m_v.insert(GetElementIterator(m_v, n_l));		
	} else if (!m_v.size() || t_s < m_v[0].at){
		// Insert a new S at the start of m_v.		
		n_l = n_h = 0;
		m_v.insert(GetElementIterator(m_v, n_l));
	} else {
		// Insert a new S at the end of m_v.
		n_l = n_h = m_v.size();
		m_v.push_back(S());		
	}	
	S& s = m_v[n_l];
	s.at = t_s;

	// Find a vector in s.v corresponding to t_e. We then assign a reference variable, e, to this (for clarity).
	if (GetElement(s.v, t_e, false, false, &n_l, &n_h)){
		// Direct hit - s.v[n_l] is the VAL_TYPE value that we need to modify.
	} else if (n_l != n_h){
		// Interpolation case. We need to insert vt just before n_h.		
		n_l = n_h;
		s.v.insert(GetElementIterator(s.v, n_l));
	} else if (!s.v.size() || t_e < s.v[0].at){
		// Insert a new E at the start of m_v.
		n_l = n_h = 0;
		s.v.insert(GetElementIterator(s.v, n_l));
	} else {
		// Insert a new E at the end of s.v.	
		n_l = n_h = s.v.size();
		s.v.push_back(E());		
	}
	E& e = s.v[n_l];
	
	// Now set the value of e
	e.at = t_e;
	e.vt = vt;
}

template <class ARG_TYPE, class VAL_TYPE>
bool MlEqBilinearInterpolator<ARG_TYPE, VAL_TYPE>::GetValue(ARG_TYPE t_s, ARG_TYPE t_e, bool bInterpolate, bool bExtrapolate, VAL_TYPE* pvt)
{
	ARG_TYPE							sl, sh;
	ARG_TYPE							el, eh;
	VAL_TYPE							l, h;							// low and high values for the return value (corresponding to elements of m_v around t_s)
			
	if (!GetElement(m_v, t_s, bInterpolate, bExtrapolate, &sl, &sh)) return false;

	// Get l
	if (!GetElement(m_v[sl].v, t_e, bInterpolate, bExtrapolate, &el, &eh)) return false;
	if (el == eh){		
		l = m_v[sl].v[el].vt;
	} else {
		l = ((m_v[sl].v[eh].vt - m_v[sl].v[el].vt) * t_e + (m_v[sl].v[el].vt * m_v[sl].v[eh].at - m_v[sl].v[eh].vt * m_v[sl].v[el].at)) / (m_v[sl].v[eh].at - m_v[sl].v[el].at);
	}

	// Get h
	if (sl == sh){		
		*pvt = l;
		return true;
	}
	if (!GetElement(m_v[sh].v, t_e, bInterpolate, bExtrapolate, &el, &eh)) return false;
	if (el == eh){	
		h = m_v[sh].v[el].vt;
	} else {
		h = ((m_v[sh].v[eh].vt - m_v[sh].v[el].vt) * t_e + (m_v[sh].v[el].vt * m_v[sh].v[eh].at - m_v[sh].v[eh].vt * m_v[sh].v[el].at)) / (m_v[sh].v[eh].at - m_v[sh].v[el].at);
	}

	// Interpolate on l and h
	*pvt = ((h - l) * t_s + (l * m_v[sh].at - h * m_v[sl].at)) / (m_v[sh].at - m_v[sl].at);
	return true;
}

#endif