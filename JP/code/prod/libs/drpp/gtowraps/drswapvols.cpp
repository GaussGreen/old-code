#include "drswapvols.h"
#include "drsymboltable.h"

const DRString KAPITAL_SWAPVOL_FILE = "swapvol.dat";

DRSwapVols::DRSwapVols(FileType type, DRString filename)
{
	if (type == SWAP) LoadFromSwapFiles(filename);
	else LoadFromKapitalFiles();
}

void DRSwapVols::LoadFromSwapFiles (DRString filename)
{
	m_freq		= 2;
	int numExp	= 20;
	int numMat	= 16;
	
	m_exps.resize(numExp);
	m_mats.resize(numMat);
	m_volMatrix.resize(numExp, numMat);
	
	char temp[20];
	
	ifstream inFile(CheckPathDelimiters(theSymbolTable.get(filename)).c_str());
	if (!inFile) 
		throw DRException("Failed to open file: ") << filename;
	
	inFile >> temp ;		// skip "Exp/Und"
	int iMat;
	for (iMat = 0; iMat < numMat; iMat++) inFile >> m_mats[iMat];
	
	for (int iExp = 0; iExp < numExp; iExp++) {
		inFile >> m_exps[iExp];
		for (iMat = 0; iMat < numMat; iMat++) {
			inFile >> temp;
			m_volMatrix [iExp][iMat] = atof(temp) / 100.;
		}
	}
}

void DRSwapVols::LoadFromKapitalFiles ()
{
	int numExp, numMat;
	
	m_freq				= 2;
	
	ifstream inFile(KAPITAL_SWAPVOL_FILE.c_str()); 
	if (!inFile) 
		throw DRException ("Failure to Open Kapital File ") << KAPITAL_SWAPVOL_FILE;
	
	char temp[20], trash[255];
	
	inFile.getline(trash, 255);
	inFile >> numExp;
	
	inFile.getline(trash, 255);
	inFile.getline(trash, 255);
	inFile >> numMat;
	
	m_exps.resize(numExp);
	m_mats.resize(numMat);
	m_volMatrix.resize(numExp, numMat);
	
	inFile.getline(trash, 255);
	inFile.getline(trash, 255);
	int iMat;
	for (iMat = 0; iMat < numMat; iMat++) {
		inFile >> m_mats[iMat];
	}
	
	for (int iExp = 0; iExp < numExp; iExp++) {
		inFile >> m_exps[iExp];
		m_exps [iExp] /= 12.;
		for (iMat = 0; iMat < numMat; iMat++) {
			inFile >> temp;
			m_volMatrix [iExp][iMat] = atof(temp) / 100.;
		}
	}
}

ostream& operator<<(ostream& s, const DRSwapVols& a)
{
	s << "Exps:\t" << a.m_exps;
	s << "Mats:\t" << a.m_mats;
	s << a.m_volMatrix << endl;
	return s;
}
