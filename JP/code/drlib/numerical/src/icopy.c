#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef COMPUTER_VAX
#include <string.h>
#else
#include <memory.h>
#endif

/* in-line BLAS does not handle negative increments */
#ifdef ANSI
void imsl_icopy(Mint n, Mint sx[], Mint incx, Mint sy[], Mint incy)
#else
void imsl_icopy(n, sx, incx, sy, incy)
        Mint             n;
        Mint           sx[];
        Mint             incx;
        Mint           sy[];
        Mint             incy;
#endif
{
        Mint             i, ix, iy;

	if (n <= 0) return;

	if (incx == 1 && incy == 1) {
	    memcpy((char*)sy, (char*)sx, n*sizeof(Mint));
        } else {
                        /* Code for unequal increments */
	    ix = 0;
            iy = 0;
            if (incx < 0) ix = (-n + 1) * incx;
            if (incy < 0) iy = (-n + 1) * incy;
            for (i = 0; i < n; i++) {
		sy[iy] = sx[ix];
                ix += incx;
                iy += incy;
            }
        }
        return;
}
