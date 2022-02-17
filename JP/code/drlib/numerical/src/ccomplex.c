#define CCOMPLEX_C
#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

/*
#define Mac2ms
*/

static Mf_complex PROTO(l_c_proc_cs,(Mf_complex z));

    /* complex unary minus operation */
#ifdef ANSI
Mf_complex imsl_c_neg( Mf_complex z )	
#else
Mf_complex imsl_c_neg( z )	
    Mf_complex z;
#endif
{
        z.re = -z.re;
        z.im = -z.im;	
        return( z );
}

    /* return the sum of two complex numbers */
#ifdef ANSI
Mf_complex imsl_c_add( Mf_complex l, Mf_complex r )	
#else
Mf_complex imsl_c_add( l, r )	
    Mf_complex l, r;
#endif
{
        l.re += r.re;
        l.im += r.im;
        return( l );
}

    /* return the difference of two complex numbers */
#ifdef ANSI
Mf_complex imsl_c_sub( Mf_complex l, Mf_complex r )	
#else
Mf_complex imsl_c_sub( l, r )	
    Mf_complex l, r;
#endif
{
	l.re -= r.re;
	l.im -= r.im;
	return( l );
}
/** Work around for using gcc on intel with -msse2 compiler option.
    The compiler can't cope with the imsl_c_div function below */
#if defined(__i386__) && defined(__linux__) && defined(__GNUC__)
static __attribute__((noinline)) void imsl_c_div1( Mf_complex* l, 
                                                   Mf_complex* r,
                                                   Mf_complex* z)
{
    Mfloat t = r->im/r->re;
    Mfloat den = r->re + r->im*t;
    z->re = (l->re + l->im*t)/den;
    z->im = (l->im - l->re*t)/den;
}

static __attribute__((noinline)) void imsl_c_div2( Mf_complex* l, 
                                                   Mf_complex* r,
                                                   Mf_complex* z)
{
    Mfloat t = r->re/r->im;
    Mfloat den = r->im + r->re*t;
    z->re = (l->im + l->re*t)/den;
    z->im = (l->im*t - l->re)/den;
}
#endif

    /* return the quotient of two complex numbers */
#ifdef ANSI
Mf_complex imsl_c_div( Mf_complex l, Mf_complex r ) 
#else
Mf_complex imsl_c_div( l, r ) 
    Mf_complex l, r;
#endif
{
    /* Smith's Algorithm, ACM number 116 (1962) */

        Mf_complex z;
#if defined(__i386__) && defined(__linux__) && defined(__GNUC__)
#else
        Mfloat    den, t;
#endif
        if( r.re == F_ZERO && r.im == F_ZERO ){
/*	    imsl_ermes(5, 1, "Complex division by zero."); */
            imsl_ermes(IMSL_TERMINAL, IMSL_COMPLEX_DIVIDE_BY_ZERO);
	    return z;
	}

        if (fabs(r.re) > fabs(r.im)){
#if defined(__i386__) && defined(__linux__) && defined(__GNUC__)
            imsl_c_div1(&l, &r, &z);
#else
	    t = r.im/r.re;
            den = r.re + r.im*t;
            z.re = (l.re + l.im*t)/den;
            z.im = (l.im - l.re*t)/den;
#endif
	}
        else{
#if defined(__i386__) && defined(__linux__) && defined(__GNUC__)
            imsl_c_div2(&l, &r, &z);
#else
	    t = r.re/r.im;
            den = r.im + r.re*t;
            z.re = (l.im + l.re*t)/den;
            z.im = (l.im*t - l.re)/den;
#endif
        }
        return( z );
}

    /* return the product of two complex numbers */
#ifdef ANSI
Mf_complex imsl_c_mul( Mf_complex l, Mf_complex r )	
#else
Mf_complex imsl_c_mul( l, r )	
    Mf_complex l, r;
#endif
{
        Mf_complex z;
        z.re = l.re*r.re - l.im*r.im;
        z.im = l.re*r.im + l.im*r.re;
        return( z );
}

#if 0
    /* general complex power */
Mf_complex imsl_c_pow( l, r )	
    Mf_complex l, r;
{
        Mf_complex temp1,ans;
        Mfloat temp2;


        temp2 = imsl_c_abs(l);
        if (temp2 == F_ZERO) {
	    temp1.re = F_ZERO;
            temp1.im = F_ZERO;
            return(temp1);
	}
        temp1.re = F_ZERO;
        temp1.im = r.im*log(temp2)+r.re*imsl_c_arg(l);
        temp1 = imsl_c_exp(temp1);


        ans.re = fabs(pow(temp2,r.re)*exp(-r.im*imsl_c_arg(l)))*temp1.re;
        ans.im = fabs(pow(temp2,r.re)*exp(-r.im*imsl_c_arg(l)))*temp1.im;


        return(ans);
}

    /* general complex raised to Mfloat power */
Mf_complex imsl_c_powf( l, r )	
    Mf_complex l;
    Mdouble r;
{
       Mf_complex temp1, temp3;
       Mfloat temp2;
    
       temp2 = imsl_c_abs(l);
       if (temp2 == F_ZERO)
	   return(l); 
       temp1.re = F_ZERO;
       temp1.im = r*imsl_c_arg(l);
       temp1 = imsl_c_exp(temp1);
       temp3.re = pow(temp2,r)*temp1.re;
       temp3.im = pow(temp2,r)*temp1.im;
       return(temp3);
}
#endif

    /* return TRUE if complex numbers are equal */
#ifdef ANSI
Mint imsl_c_eq( Mf_complex l, Mf_complex r )
#else
Mint imsl_c_eq( l, r )
    Mf_complex l, r;
#endif
{
        if( l.re == r.re  &&  l.im == r.im )
                return( 1 );
        return( 0 );
}

#ifdef DOUBLE

    /* f_complex conversion to d_complex */
#ifdef ANSI
d_complex imsl_zc_convert( f_complex z )        
#else
d_complex imsl_zc_convert( z )        
    f_complex z;
#endif
{
        d_complex dz;

        dz.re = z.re;
        dz.im = z.im;
        return( dz );
}

#else   

    /* d_complex conversion to f_complex */
#ifdef ANSI
f_complex imsl_cz_convert( d_complex dz )
#else
f_complex imsl_cz_convert( dz )
    d_complex dz;
#endif
{
        f_complex z;

        z.re = dz.re;
        z.im = dz.im;
        return( z );
}


#endif  /* ifdef DOUBLE */

    /* convert floats to complex */
#ifdef ANSI
Mf_complex imsl_cf_convert( Mfloat r, Mfloat i )
#else
Mf_complex imsl_cf_convert( r, i )
    Mfloat r, i;
#endif
{
        Mf_complex z;
        z.re = r;
        z.im = i;
	return( z );
}

    /* complex (real component) to float conversion */
#ifdef ANSI
Mfloat imsl_fc_convert( Mf_complex z )
#else
Mfloat imsl_fc_convert( z )
    Mf_complex z;
#endif
{
        return( z.re );
}

    /* return the complex conjugate */
#ifdef ANSI
Mf_complex imsl_c_conjg( Mf_complex z )
#else
Mf_complex imsl_c_conjg( z )
    Mf_complex z;
#endif
{
        z.im = -z.im;
        return( z );
}

    /* return the imaginary component of a complex number */
#ifdef ANSI
Mfloat imsl_c_aimag( Mf_complex z )
#else
Mfloat imsl_c_aimag( z )
    Mf_complex z;
#endif
{
        return( z.im );
}

    /* return square root of complex number */
#ifdef ANSI
Mf_complex imsl_c_sqrt( Mf_complex z )
#else
Mf_complex imsl_c_sqrt( z )
    Mf_complex z;
#endif
{
       Mfloat r;
       Mf_complex cs,z_hat;
       r = imsl_c_abs(z);
       if (r==F_ZERO) 
        return(z);
       if (fabs(z.re)<fabs(z.im)){
	   if (z.im >= F_ZERO){
	       z_hat.re = z.im;
               z_hat.im = -z.re;
               r     = F_HALF*r;
               cs = l_c_proc_cs(z_hat);
               r = sqrt(r);
               z.re = r*(cs.re - cs.im);
               z.im = r*(cs.re + cs.im); 
	   }
           else           {
	       z_hat.re = -z.im;
               z_hat.im = z.re;
               r     = F_HALF*r;
               cs = l_c_proc_cs(z_hat);
               r = sqrt(r);
               z.re = r*(cs.re + cs.im);
               z.im = r*(cs.im - cs.re) ;
	   }
       }
       else                     {
           if (z.re < F_ZERO) {
	       z.re = -z.re;
               z.im = -z.im;
               cs = l_c_proc_cs(z);
               r = sqrt(r);
               if (z.im > F_ZERO) 
               {
                  z.re = r*cs.im;
                  z.im = -r*cs.re; 
	       }
               else
               {
                  z.re = -r*cs.im;
                  z.im = r*cs.re; 
	       }
	   }
           else           {
               cs = l_c_proc_cs(z);
               r = sqrt(r);
               z.re = r*cs.re;
               z.im = r*cs.im; 
	   }
       }
       return(z);
}

    /* Used in support of CSQRT */
#ifdef ANSI
static Mf_complex l_c_proc_cs( Mf_complex z )
#else
static Mf_complex l_c_proc_cs( z )
    Mf_complex z;
#endif
{
    static Mfloat base = 0.0;
    static Mfloat tol  = 0.0;
    Mfloat tau, sigma, t;
    Mf_complex ans;

    if (base == F_ZERO) base = pow(F_TEN,imsl_amach(5));
    if (tol  == F_ZERO) tol = sqrt(imsl_amach(4))/base;
    /* Compute tau = tan(2*theta) = y/x */
    tau = z.im/z.re;
    if (fabs(tau) <= tol){
	t = F_HALF*tau;
        ans.re = F_ONE;
        ans.im = t;
    }
    else{
        sigma = F_ONE/tau;
        t = F_ONE/(sigma+((sigma>=F_ZERO)?sqrt(F_ONE+sigma*sigma):-sqrt(F_ONE+sigma*sigma)));
        ans.re = F_ONE/sqrt(F_ONE+t*t);
        ans.im = ans.re*t;
    }
    return ( ans );
}

    /* complex natural log */
#ifdef ANSI
Mf_complex imsl_c_log( Mf_complex z )  
#else
Mf_complex imsl_c_log( z )  
  Mf_complex z;
#endif
{
Mf_complex l;
        l.re = log( imsl_c_abs(z) );
        l.im = imsl_c_arg(z);     /* the direct angle */
          return( l );
}

    /* complex exponential */
#ifdef ANSI
Mf_complex imsl_c_exp( Mf_complex z )
#else
Mf_complex imsl_c_exp( z )
    Mf_complex z;
#endif
{
    Mfloat exp_re;
        if( z.re == F_ZERO )
                exp_re = F_ONE;
        else
                exp_re = exp(z.re);        
        if( z.im == F_ZERO ){        /* real only complex number */
                z.re = exp_re;
                return( z );
                }
	z.re = exp_re*cos(z.im);
        z.im = exp_re*sin(z.im);
        return( z );
}

    /* sine of complex number */
#ifdef ANSI
Mf_complex imsl_c_sin( Mf_complex z )
#else
Mf_complex imsl_c_sin( z )
    Mf_complex z;
#endif
{
        Mf_complex zs;
        zs.re = sin( z.re )*cosh( z.im );
        zs.im = cos( z.re )*sinh( z.im );
        return( zs );
}

    /* cosine of complex number */
#ifdef ANSI
Mf_complex imsl_c_cos( Mf_complex z )
#else
Mf_complex imsl_c_cos( z )
    Mf_complex z;
#endif
{
        Mf_complex zc;
        zc.re = cos(z.re)*cosh(z.im);
        zc.im = -sin(z.re)*sinh(z.im);
        return( zc );
}

    /* complex absolute value */
#ifdef ANSI
Mfloat imsl_c_abs( Mf_complex z )
#else
Mfloat imsl_c_abs( z )
    Mf_complex z;
#endif
{
/* #ifndef COMPUTER_ALFAC_IEEE */
#ifndef Mac2ms
    return hypot(z.re,z.im);
#else
	static Mfloat big = 0.0;
	static Mfloat small;
	static Mfloat prec;
	static Mfloat one_over_hi;
	static Mfloat one_over_lo;
	static Mfloat one_over_base;
	static Mfloat base;
        static Mfloat hi;
	static Mfloat lo;
	Mfloat t;

        if( z.re == F_ZERO  &&  z.im == F_ZERO )
          return(F_ZERO);

                            /* Get machine dependent constants.*/
	if (big == F_ZERO) {
	    big = imsl_amach(2);
	    small = imsl_amach(1);
	    prec = imsl_amach(4);
	    base = pow(F_TEN,imsl_amach(5));
	    one_over_base = F_ONE/base;
	    hi = sqrt(big*prec)*one_over_base;
	    lo = sqrt(small/prec)*base;
            one_over_hi = F_ONE/hi;
            one_over_lo = F_ONE/lo;
	}                            
        t = fabs(z.re) + fabs(z.im);
                             /*  Check for special cases */
        if (t>hi){
           return(hi*sqrt(  (one_over_hi*z.re)*(one_over_hi*z.re)+
                           (one_over_hi*z.im)*(one_over_hi*z.im)   ));
        }
        if(t<lo){
           return(lo*sqrt(  (one_over_lo*z.re)*(one_over_lo*z.re)+
                           (one_over_lo*z.im)*(one_over_lo*z.im)   ));
        }
                            /* This case will be the one that
                               is usually used for the computation
                               of the returned value. */
        return(sqrt(z.re*z.re + z.im*z.im));
#endif
}

       /* support for complex library functions */

    /* direct angle (in radians) of a complex number */
#ifdef ANSI
Mfloat imsl_c_arg( Mf_complex z )
#else
Mfloat imsl_c_arg( z )
    Mf_complex z;
#endif
{
        if( z.re == F_ZERO && z.im == F_ZERO ) 
           /* arg(z) is arbitrary in this case. */
           return(F_ZERO);
        return ( atan2(z.im,z.re));
}
