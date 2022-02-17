#include "ito33/utile.hpp"

DLLEXPORT void ExcelDeleteDouble(long _lAd)
{
	delete [] (double *)_lAd;
}

DLLEXPORT void ExcelDeleteInt(long _lAd)
{
	delete [] (int*)_lAd;
}

DLLEXPORT int ExcelNewDouble(double **_plAddresse, int _iTaille)
{
	double
		*pdV = new double [_iTaille];
	*_plAddresse = pdV;
	return 1;
}

DLLEXPORT int ExcelNewInt(long &_lAddresse, int _iTaille)
{
	int *piV = new int [_iTaille];
	_lAddresse = (long)piV;
	
	return 1;
}

DLLEXPORT void ExcelSetTableauDouble(long _lAd, int _iIndice, double _dValeur)
{
	double *pdV = (double *)_lAd;
	pdV[_iIndice] = _dValeur;
}


DLLEXPORT void ExcelSetTableauInt(long _pi, int _iIndice, int _iValeur)
{
	((int *)_pi)[_iIndice] = _iValeur;
}


DLLEXPORT double ExcelGetTableauDouble(long _lAd, int _iIndice)
{
	return ((double *)_lAd)[_iIndice];
}


DLLEXPORT int ExcelGetTableauInt(long _pi, int _iIndice)
{
	return ((int *)_pi)[_iIndice];
}
