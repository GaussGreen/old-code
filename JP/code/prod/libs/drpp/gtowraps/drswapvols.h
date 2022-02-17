#ifndef __drswapvol__h
#define __drswapvol__h

#include "drstruct.h"
#include "drgtowrapenum.h"

class DRSwapVols {
public:
	DRSwapVols() {}
	DRSwapVols(FileType, DRString filename = "swo_vol.usd");

	int freq();
	DArray& exps();
	DArray& mats();
	DMatrix& matrix();

	friend ostream& operator<<(ostream&, const DRSwapVols&);

protected:
	
	int	m_freq;
	DArray	m_exps;
	DArray	m_mats;
	DMatrix	m_volMatrix;

	void LoadFromSwapFiles(DRString filename);
	void LoadFromKapitalFiles ();
};

inline int DRSwapVols::freq() {return m_freq;}
inline DArray& DRSwapVols::exps() {return m_exps;}
inline DArray& DRSwapVols::mats() {return m_mats;}
inline DMatrix& DRSwapVols::matrix() {return m_volMatrix;}

#endif
