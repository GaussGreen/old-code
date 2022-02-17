#include "ICMKernel\util\icm_fft.h"
#include "math.h"

//#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

/* About >> and << operators
expr1 << expr2 Left Shift
expr1 >> expr2 Right Shift

(left assoc) 
<< returns expr1 shifted left expr2 bit positions; 0 bits are shifted into the low order bits.
Mathematically equivalent to expr1 * 2expr2 (note that there is no exponentiation operator in C, nor in netscape 1.1n :-)). 

>> returns expr1 shifted right expr2 bit positions; if the type of expr1 is signed, 1 bits are shifted into the high order bits, otherwise 0 bits are shifted in.
Mathematically equivalent to expr1 / 2expr2. 
Example: 
b = 0x12 << 2 

b is assigned the value 0x48.  
*/

//--------------------------------------------------------------------------------------------------------------------------------------------
double
fourier::m_2PI=6.28318530717959;
//--------------------------------------------------------------------------------------------------------------------------------------------
// From Numerical Recipes
void
fourier::fft(double* data, // 1+nn*2 array of complex numbers : [nothing, Re0, Im0, Re1, Im1, ...]
						// arrays in Numerical begins at index 1
			 unsigned long nn, // Numbers of Points
			 int isign) // 1 or -1
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			std::swap(data[j],data[i]);
			std::swap(data[j+1],data[i+1]);
		}
		m = n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		//theta=isign*6.28318530717959/mmax;
		theta=isign*m_2PI/mmax;
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
//--------------------------------------------------------------------------------------------------------------------------------------------
bool
fourier::test(std::string& errStr)
{
	bool ret(true);
	
	errStr = "Testing Fast Fourier Transform";

	double* d = new double [10];
	//fourier::fft(d,  nn, -1);

	delete [] d;

	return ret;
}
//-------------------------------------------------------------------------------
