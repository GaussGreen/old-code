// drbasevols.h: interface for the drbasevols class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DRBASEVOLS_H__5DD73653_A7A6_11D1_80A4_00C04FB91C08__INCLUDED_)
#define AFX_DRBASEVOLS_H__5DD73653_A7A6_11D1_80A4_00C04FB91C08__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include "drcurve.h"
#include "drswapvols.h"

//k Represents a base vol curve
//k	In my mind's eye, a base vol curve is:
//k		A set of dates and vols.
//k		A rate description: freq, maturity, and 
class DRBaseVols : public DRCurve {
public:
	DRBaseVols (); //
	DRBaseVols (TCurve* volCurve, int volFreq, int volMat); //

	DRBaseVols (FileType fileType, DRDate& envDate, 
		DRString filename = "Base_vol.usd", int mat = 3);

	DRBaseVols (DRDate&, DRSwapVols&, double matForCalib);

	DRBaseVols (DRDate& envDate, MbsIO&);
	~DRBaseVols();
	
	void ChangeToSimpleBaseVolInLognormalSpace (long newMat, DRDate& envDate, 
		DRCurve& zeroCurve, DArray& calibParams);

	int maturity() const;	// maturity in months
	int freq() const {return m_freq;}		// number of cpns per year
	bool fromCapVol() const;
	bool operator==(const DRBaseVols&) const;
	
protected:
	int m_volMat;	// maturity in months
	int m_freq;		// number of coupons per year
	bool m_fromCapVol;

	DRCurve m_originalCurve;
	int m_originalMat, m_originalFreq;

	void Print(ostream& stream) const;
	void LoadFromKapitalFiles(DRDate&);
	void LoadFromSwapFiles (DRDate&, DRString, int mat);
};

inline int DRBaseVols::maturity() const {return m_volMat;}
inline bool DRBaseVols::fromCapVol() const {return m_fromCapVol;}

#endif // !defined(AFX_DRBASEVOLS_H__5DD73653_A7A6_11D1_80A4_00C04FB91C08__INCLUDED_)

