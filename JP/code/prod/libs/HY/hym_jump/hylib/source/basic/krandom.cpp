#include <math.h>
#include "krandom.h"
#include "kexception.h"


KNormRandom::KNormRandom(const KMatrix<double> &corr,int randSeed)
:_corr(corr),_seed2(123456789),_iy(0),_isSecond(false), _seed(randSeed)
{
	int	dim = _corr.row_size();
	int	i, j, k;
	for(j=0; j<dim; j++)
	{
		for(i=0; i<j; i++)
		{
			_corr[j][i] =  corr[j][i];
			for(k=0; k<i; k++)
				_corr[j][i] -= _corr[i][k] * _corr[j][k];
			_corr[j][i] /= _corr[i][i];
		}
		_corr[j][j] = 1.;
		for(i=0; i<j; i++)
			_corr[j][j] -= _corr[j][i] * _corr[j][i];

		if(_corr[j][j] <= 0.)
			throw KException("Correlation matrix has problems!");
		else
			_corr[j][j] = sqrt(_corr[j][j]);
	}
}

KValarray<double>	KNormRandom::get_sample(void)
{
	int	i, j;
	int	dim = _corr.row_size();

	KValarray<double>	rand(dim);
	KValarray<double> &tmpRan = get_random(dim);

	for(i=0; i<dim; i++)
	{
		rand[i] = _corr[i][i] * tmpRan[i];
		for(j=0; j < i; j++)
			rand[i] += _corr[i][j] * tmpRan[j];
	}
	return	rand;
}

KValarray<double>	KNormRandom::get_random(int d)
{
	int	i, k = d/2;
	KValarray<double>	rand(d);
	double	*randPtr = rand;
	double fac,rsq,v1,v2;

	for(i=0; i<k; i++)
	{
		do 
		{
			v1 = 2 * _random_generator() - 1;
			v2 = 2 * _random_generator() - 1;
			rsq = v1 * v1 + v2 * v2;
		} 
		while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0*log(rsq)/rsq);
		*randPtr++ = v1 * fac;
		*randPtr++ = v2 * fac;
	}
	//last one
	if(d%2 != 0)
	{
		do
		{
			v1 = 2 * _random_generator() - 1;
			v2 = 2 * _random_generator() - 1;
			rsq = v1 * v1 + v2 * v2;
		} 
		while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0*log(rsq)/rsq);
		*randPtr = v1 * fac;
	}
	return	rand;
}


#define EPS		1.2e-7
#define RNMX	(1.0-EPS)
#define IM1		2147483563
#define IM2		2147483399
//IMM1 = IM1 - 1
#define IMM1	2147483562
#define IA1		40014
#define IA2		40692
#define IQ1		53668
#define IQ2		52774
#define IR1		12211
#define IR2		3791
//NDIV1 = (1+IMM1/32) 
#define NDIV1	67108862


double	KNormRandom::_random_generator(void)
{
	int j, k;

	if (_seed <= 0) 
	{
		if (_seed > -1) _seed = 1;
		else _seed *= -1;
		_seed2 = _seed;
		for (j = 39; j>=0; j--) 
		{
			k= _seed /IQ1;
			_seed=IA1*(_seed%IQ1)-k*IR1;
			if (_seed < 0)	_seed += IM1;
			if (j < 32) _iv[j] = _seed;
		}
		_iy=_iv[0];
	}

	k = _seed/IQ1;
	_seed = IA1*(_seed%IQ1)-k*IR1;
	if (_seed < 0) _seed += IM1;
	k= _seed2/IQ2;
	_seed2 = IA2*(_seed2%IQ2)-k*IR2;
	if (_seed2 < 0)	_seed2 += IM2;
	j=_iy/NDIV1;
	_iy=_iv[j]-_seed2;
	_iv[j] = _seed;
	if (_iy < 1) _iy += IMM1;
	double temp = _iy;
	if ((temp /= IM1) > RNMX) 
		return RNMX;
	else 
		return temp;
}


#undef IM1
#undef IM2
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NDIV1
#undef EPS
#undef RNMX
