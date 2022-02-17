#ifndef _RANDOM_H_
#define _RANDOM_H_

#include "kmatrix.h"
#include "kvalarray.h"
//Normal random generator 

class KNormRandom
{
public:
	KNormRandom(int randSeed = -1):_seed2(123456789),_iy(0),_isSecond(false), 
							_seed(randSeed){}
	KNormRandom(const KMatrix<double> &corr,int randSeed = -1);
	~KNormRandom(void){}

	KMatrix<double>	get_corr()const{return _corr;}
	//get a sample from normal distribution
	KValarray<double>	get_sample(void); 
	//get a sample form a d-dim random number
	KValarray<double>	get_random(int d);
private:
//NRC: ran2
	int		_seed;
	//ran2 static variables
	int		_seed2;
	int		_iy;
	int		_iv[32];
	//normal random parameters
	bool	_isSecond;
	double	_secondRand;

	double	_random_generator();
//correlation matrix
	KMatrix<double>	_corr;
};

#endif 
