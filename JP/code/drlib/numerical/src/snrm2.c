#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/* INCX is not used, but leave the calling sequence intact. */
#ifdef ANSI
Mfloat imsl_snrm2(Mint n, Mfloat *sx, Mint incx)
#else
Mfloat imsl_snrm2(n, sx, incx)
    Mint n;
    Mfloat *sx;
    Mint incx;
#endif
{
	Mfloat norm = F_ZERO;
	Mfloat l;
	Mfloat u;
	Mfloat eta;
	Mfloat a;
	Mfloat b;
	Mfloat c;
	Mfloat d;
	Mfloat t = F_ZERO;
	Mint i;

	l = imsl_amach(1);
	u = imsl_amach(2);
	eta = imsl_amach(4);

	a = sqrt(l/eta);
	b = F_ONE/(a*eta);
	c = sqrt(u*eta);
#if defined(COMPUTER_ALFAC_IEEE) && defined(DOUBLE)
	d = (sqrt(l)/eta)/1.340781e+154;
#else
	d = (sqrt(l)/eta)/sqrt(u);
#endif

	for (i=0; i<=n-1; ++i)
		t += fabs(*(sx+i));

	if (t < a) {
		for (i=0; i<=n-1; ++i)
			norm += (b**(sx+i))*(b**(sx+i));
		norm = sqrt(norm);
		norm = norm/b;
	}

	else if (t > c) {
		for (i=0; i<=n-1; ++i)
			norm += (d**(sx+i))*(d**(sx+i));
		norm = sqrt(norm);
		norm = norm/d;
	}

	else {
		for (i=0; i<=n-1; ++i)
			norm += *(sx+i)**(sx+i);
		norm = sqrt(norm);
	}

	return (norm);
}
