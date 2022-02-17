/* Inline version of BLAS 1 */

#define INLINE_BLAS

/* convention: prepend  blas__ to name to avoid conflict with names in its arguments */



#define AXPY(PREC, blas__n_, blas__a_, blas__x_, blas__incx_, blas__y_, blas__incy_)	\
    {	register PREC blas__a, *blas__x, *blas__y;			\
	register int blas__k;						\
    	blas__a = blas__a_; blas__x = blas__x_; blas__y = blas__y_;	\
	if (blas__incx_ < 0) blas__x += ((blas__n_)-1)*(blas__incx_);	\
	if (blas__incy_ < 0) blas__y += ((blas__n_)-1)*(blas__incy_);	\
	for (blas__k=0; blas__k< blas__n_; blas__k++) {			\
	    *blas__y += blas__a * (*blas__x);				\
	    blas__x += blas__incx_;   blas__y += blas__incy_; }		}

#define saxpy(N,A,X,INCX,Y,INCY)	AXPY(float,N,A,X,INCX,Y,INCY)

#define daxpy(N,A,X,INCX,Y,INCY)	AXPY(double,N,A,X,INCX,Y,INCY)





#define saxpy2(blas__n_, blas__a_, blas__x_, blas__incx_, blas__y_, blas__incy_) \
    {	register float blas__a, *blas__x, *blas__y;			\
	register int blas__i, blas__j, blas__k;				\
    	blas__a = blas__a_; blas__x = blas__x_; blas__y = blas__y_;	\
	blas__i = (blas__incx_ < 0) ? ((blas__n_)-1)*(blas__incx_) : 0;	\
	blas__j = (blas__incy_ < 0) ? ((blas__n_)-1)*(blas__incy_) : 0;	\
	for (blas__k=0; blas__k< blas__n_; blas__k++) {			\
	    blas__y[blas__j] += blas__a * blas__x[blas__i];		\
	    blas__i += blas__incx_;   blas__j += blas__incy_; }		}



#define COPY(PREC, blas__n_, blas__x_, blas__incx_, blas__y_, blas__incy_) \
    {	register PREC *blas__x, *blas__y;				\
	register int blas__i, blas__j, blas__k;				\
	blas__x = blas__x_;   blas__y = blas__y_;				\
	if (blas__x != blas__y || (blas__incx_) != (blas__incy_)) {	\
	blas__i = (blas__incx_ < 0) ? ((blas__n_)-1)*(blas__incx_) : 0;	\
	blas__j = (blas__incy_ < 0) ? ((blas__n_)-1)*(blas__incy_) : 0;	\
	for (blas__k=0; blas__k< blas__n_; blas__k++) {			\
	    blas__y[blas__j] = blas__x[blas__i];			\
	    blas__i += blas__incx_;   blas__j += blas__incy_; }}}

#if defined(CMATH_MINT_TO_LONG)

#define icopy(N,X,INCX,Y,INCY)		COPY(long,N,X,INCX,Y,INCY)

#else

#define icopy(N,X,INCX,Y,INCY)		COPY(int,N,X,INCX,Y,INCY)

#endif



#define scopy(N,X,INCX,Y,INCY)		COPY(float,N,X,INCX,Y,INCY)

#define dcopy(N,X,INCX,Y,INCY)		COPY(double,N,X,INCX,Y,INCY)





#define SWAP(PREC, blas__n_, blas__x_, blas__incx_, blas__y_, blas__incy_)	\
    {	register PREC *blas__x, *blas__y, blas__t;				\
	register int blas__i, blas__j, blas__k;					\
	blas__x = blas__x_;   blas__y = blas__y_;				\
	blas__i = (blas__incx_ < 0) ? ((blas__n_)-1)*(blas__incx_) : 0;		\
	blas__j = (blas__incy_ < 0) ? ((blas__n_)-1)*(blas__incy_) : 0;		\
	for (blas__k=0; blas__k< blas__n_; blas__k++) {				\
	    blas__t = blas__x[blas__i];  blas__x[blas__i] = blas__y[blas__j];   \
	    blas__y[blas__j] = blas__t;						\
	    blas__i += blas__incx_;   blas__j += blas__incy_;}}

#define sswap(N,X,INCX,Y,INCY)		SWAP(float,N,X,INCX,Y,INCY)

#define dswap(N,X,INCX,Y,INCY)		SWAP(double,N,X,INCX,Y,INCY)





#define SCAL(PREC, blas__n_, blas__a_, blas__x_, blas__incx_)		\
    {	register PREC blas__a, *blas__x;				\
	register int blas__k;						\
    	blas__a = blas__a_; blas__x = blas__x_;				\
	for (blas__k=0; blas__k< blas__n_; blas__k++) {			\
	    *blas__x *= blas__a;					\
	    blas__x += blas__incx_; }}

#define sscal(N,A,X,INCX)		SCAL(float,N,A,X,INCX)

#define dscal(N,A,X,INCX)		SCAL(double,N,A,X,INCX)





#define SET(PREC, blas__n_, blas__a_, blas__x_, blas__incx_)	\
    {	register PREC blas__a, *blas__x;			\
	register int blas__k;					\
    	blas__a = blas__a_; blas__x = blas__x_;			\
	for (blas__k=0; blas__k< blas__n_; blas__k++) {		\
	    *blas__x = blas__a;					\
	    blas__x += blas__incx_; }}



#if defined(CMATH_MINT_TO_LONG)

#define iset(N,A,X,INCX)		SET(long,N,A,X,INCX)

#else

#define iset(N,A,X,INCX)		SET(int,N,A,X,INCX)

#endif

#define sset(N,A,X,INCX)		SET(float,N,A,X,INCX)

#define dset(N,A,X,INCX)		SET(double,N,A,X,INCX)





#define ADD(PREC, blas__n_, blas__a_, blas__x_, blas__incx_)		\
    {	register PREC blas__a, *blas__x;				\
	register int blas__k;						\
    	blas__a = blas__a_;  blas__x = blas__x_;			\
	for (blas__k=0; blas__k< blas__n_; blas__k++) {			\
	    *blas__x += blas__a;					\
	    blas__x += blas__incx_; }}

#define sadd(N,A,X,INCX)		ADD(float,N,A,X,INCX)

#define dadd(N,A,X,INCX)		ADD(double,N,A,X,INCX)

