#include "ICMKernel\\random\icm_RandomRanmar.h"
#include "ICMKernel\util\icm_macro.h"


ICM_RandomRanmar::ICM_RandomRanmar() 
{ 
	Init(1802,9373) ;
}

ICM_RandomRanmar::ICM_RandomRanmar(int _ij, int _kl) 
{
	Init(_ij,_kl) ;
}

ICM_RandomRanmar::ICM_RandomRanmar(const ICM_RandomRanmar& ref):ICM_RandomGenerator(ref)
{
	its_ij = ref.its_ij; 
	its_kl = ref.its_kl; 
	
}

ARM_Object * ICM_RandomRanmar::Clone(void)
{
	return new ICM_RandomRanmar(*this);
}

ICM_RandomRanmar& ICM_RandomRanmar::operator=(const ICM_RandomRanmar& ref)
{
	if (this!=&ref)
	{
		this->~ICM_RandomRanmar(); 
		new(this)ICM_RandomRanmar(ref);
	}
	return *this; 
}


ICM_RandomRanmar::~ICM_RandomRanmar() {}


void ICM_RandomRanmar::Init(int _ij, int _kl)
{
		if(_ij <0 || _kl<0)
			ICMTHROW(ERR_INVALID_ARGUMENT,"seeds must be >0") ; 
		its_ij = _ij;
		its_kl = _kl;
		rmarin(its_ij,its_kl); 
		its_testvar = false;

}

void ICM_RandomRanmar::rmarin(int ij, int kl)
{
/*
C This is the initialization routine for the random number generator RANMAR()
C NOTE: The seed variables can have values between:    0 <= IJ <= 31328
C                                                      0 <= KL <= 30081
C The random number sequences created by these two seeds are of sufficient
C length to complete an entire calculation with. For example, if sveral
C different groups are working on different parts of the same calculation,
C each group could be assigned its own IJ seed. This would leave each group
C with 30000 choices for the second seed. That is to say, this random
C number generator can create 900 million different subsequences -- with
C each subsequence having a length of approximately 10^30.
C
C Use IJ = 1802 & KL = 9373 to test the random number generator. The
C subroutine RANMAR should be used to generate 20000 random numbers.
C Then display the next six random numbers generated multiplied by 4096*4096
C If the random number generator is working properly, the random numbers
C should be:
C           6533892.0  14220222.0  7275067.0
C           6172232.0  8354498.0   10633180.0
*/
	int i, j, k, l, ii, jj, m;
	float s, t;

	if (its_ij<0 || its_ij>31328 || its_kl<0 || its_kl>30081) {
	
		throw;
		/// puts("The first random number seed must have a value between 0 and 31328.");
		/// puts("The second seed must have a value between 0 and 30081.");
		/// exit(1);
	}

	i = (its_ij/177)%177 + 2;
	j = its_ij%177 + 2;
	k = (its_kl/169)%178 + 1;
	l = its_kl%169;

	for (ii=1; ii<=97; ii++) {
		s = 0.0;
		t = 0.5;
		for (jj=1; jj<=24; jj++) {
			m = (((i*j)%179)*k) % 179;
			i = j;
			j = k;
			k = m;
			l = (53*l + 1) % 169;
			if ((l*m)%64 >= 32) s += t;
			t *= 0.5;
		}
		its_u[ii] = s;
	}

	its_c1 = 362436.0 / 16777216.0;
	its_cd = 7654321.0 / 16777216.0;
	its_cm = 16777213.0 / 16777216.0;

	its_i97 = 97;
	its_j97 = 33;

	its_testvar = true;
}

void ICM_RandomRanmar::reset(){
	rmarin(its_ij, its_kl);
}


void ICM_RandomRanmar::setParameters(const std::string& sParamName, double dParamValue)
{
	if (sParamName == "ij")
		its_ij = dParamValue; 
	if (sParamName == "kl")
		its_kl = dParamValue; 

}


void ICM_RandomRanmar::View(char* id, FILE* ficOut)
{
	FILE* fOut;
	char  fOutName[200];

	if ( ficOut == NULL )
	{
	ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

	(void) unlink(fOutName);

	fOut = fopen(fOutName, "w"); 
	}
	else
	{
	fOut = ficOut;
	} 

	int size =0;

	fprintf(fOut, "\t\t\t ----------------- Random Generator Ranmar ----------------- \n\n");

	fprintf(fOut, " Seed one : %d \n",its_ij);
	fprintf(fOut, " Seed two: %d \n",its_kl);
//	fprintf(fOut, " Random nb : %f \n",itsRandom);
	fprintf(fOut, "\n");

	ICM_RandomGenerator::View(id, fOut);

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}

}