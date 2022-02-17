#ifndef ICM_FOURIER_TOOLS
#define ICM_FOURIER_TOOLS

#include <string>

class fourier
{
public:
	fourier() {}
	virtual ~fourier() {}

	static void fft(double* data, unsigned long nn, int isign);
	
	static double	m_2PI;
	static bool test(std::string& errStr);
};

#endif