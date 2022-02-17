#include "stdafx.h"

#include "MlEqNormal.h"
#include "MlEqMaths.h"

#include <math.h>
//#include <algorithm>



static double normal_error( double  x );
static double normal_error_comp( double  a );
static double drezn( double x,  double y,  double rho );

extern double x20[], w20[], x40[], w40[];

const double  SQRT2OVER2 = 7.07106781186547524401E-1;
const double  MAXLOG     = 8.8029691931113054295988E1;    // log(2**127) 
const double  SQRT2PI = 2.50662827463;


static const double EPS = 1.0e-07;

double normal_density( double x ) 
{
    return( exp(-x*x/2.0)/SQRT2PI );
}

double normal_density( double x, double y, double rho )
{
    double result, a = 1 - rho * rho;
    if( fabs( a ) > EPS ) {
        result  = exp( -( x * x - 2.0 * x * y * rho + y * y ) / ( 2.0 * a ) );
        result /= ( SQRT2PI * SQRT2PI * sqrt( a ) );
    } else {
        result = normal_density( x ) * normal_density( y );
    }
    return( result );
}

double polynomial( double x, double coef[], int N )
{
    double *p = coef;
    double ans = *p++;
    do ans = ans * x + *p++;
    while( --N );
    return( ans );
}

double polynomial_1( double x, double coef[], int N )
{
    double *p = coef;
    double ans = x + *p++;
    int i = N-1;
    do ans = ans * x  + *p++;
    while( --i );
    return( ans );
}

static  double Mtools_nep[] = {
 2.46196981473530512524E-10, 5.64189564831068821977E-1, 7.46321056442269912687E0,
 4.86371970985681366614E1, 1.96520832956077098242E2, 5.26445194995477358631E2,
 9.34528527171957607540E2, 1.02755188689515710272E3, 5.57535335369399327526E2};
static  double Mtools_neq[] = {/* 1.00000000000000000000E0,*/
 1.32281951154744992508E1, 8.67072140885989742329E1, 3.54937778887819891062E2,
 9.75708501743205489753E2, 1.82390916687909736289E3, 2.24633760818710981792E3,
 1.65666309194161350182E3, 5.57535340817727675546E2};
static  double Mtools_ner[] = {
 5.64189583547755073984E-1, 1.27536670759978104416E0, 5.01905042251180477414E0,
 6.16021097993053585195E0, 7.40974269950448939160E0, 2.97886665372100240670E0};
static  double Mtools_nes[] = {/* 1.00000000000000000000E0,*/
 2.26052863220117276590E0, 9.39603524938001434673E0, 1.20489539808096656605E1,
 1.70814450747565897222E1, 9.60896809063285878198E0, 3.36907645100081516050E0};
static  double Mtools_nip0[5] = {
-5.99633501014107895267E1, 9.80010754185999661536E1,-5.66762857469070293439E1,
 1.39312609387279679503E1,-1.23916583867381258016E0,};
static  double Mtools_niq0[8] = {/* 1.00000000000000000000E0,*/
 1.95448858338141759834E0, 4.67627912898881538453E0, 8.63602421390890590575E1,
-2.25462687854119370527E2, 2.00260212380060660359E2,-8.20372256168333339912E1,
 1.59056225126211695515E1,-1.18331621121330003142E0,};
static  double Mtools_nip1[9] = {
 4.05544892305962419923E0, 3.15251094599893866154E1, 5.71628192246421288162E1,
 4.40805073893200834700E1, 1.46849561928858024014E1, 2.18663306850790267539E0,
-1.40256079171354495875E-1,-3.50424626827848203418E-2,-8.57456785154685413611E-4,};
static  double Mtools_nip2[9] = {
  3.23774891776946035970E0,6.91522889068984211695E0,3.93881025292474443415E0,
  1.33303460815807542389E0,2.01485389549179081538E-1,1.23716634817820021358E-2,
  3.01581553508235416007E-4,2.65806974686737550832E-6,6.23974539184983293730E-9,};
static  double Mtools_niq1[8] = {/*  1.00000000000000000000E0,*/
 1.57799883256466749731E1, 4.53907635128879210584E1, 4.13172038254672030440E1,
 1.50425385692907503408E1, 2.50464946208309415979E0,-1.42182922854787788574E-1,
-3.80806407691578277194E-2,-9.33259480895457427372E-4,};
static  double Mtools_niq2[8] = {/*  1.00000000000000000000E0,*/
  6.02427039364742014255E0,3.67983563856160859403E0,1.37702099489081330271E0,
  2.16236993594496635890E-1,1.34204006088543189037E-2,3.28014464682127739104E-4,
  2.89247864745380683936E-6,6.79019408009981274425E-9,};

/*C------------------------------------------------------------------------
      FUNCTION:  double normal( double a);
       PURPOSE:  Normal distribution function
                 Returns the area under the Gaussian probability density
                 function, integrated from Minus infinity to x:
                                            x
                                             -
                                   1        | |          2
                  normal(x)  = ---------    |    exp( - t /2 ) dt
                               sqrt(2pi)  | |
                                           -
                                          -inf.
                             =  ( 1 + normal_error(z) ) / 2
                             =  normal_error_comp(z) / 2
                 where z = x/sqrt(2). Computation is via the functions
                 normal_error and normal_error_comp.
          CALL:  y = normal( x );
     ARGUMENTS:  double x, y, normal();
        EFFECT:
   PORTABILITY:  This function is UNIX-DOS-K&R portable.
     DEVELOPED:  4/87 (GJV)
------------------------------------------------------------------------C*/
double normal( double a )
{
    double x;
    double y;
    double z;

#ifdef WIN32

    x = a * SQRT2OVER2;
    z = fabs( x );
    if( z < SQRT2OVER2 )
        y = 0.5 + 0.5 * normal_error( x );
    else {
        y = 0.5 * normal_error_comp( z );
        if( x > 0.0 )
            y = 1.0 - y;
    }

#else

    a = a/sqrt(2.0);
    y = 0.5*(1.0+erf(a));

#endif

    return( y );
}

double normal_inv( double y0 )
{
    double         x;
    double         y;
    double         z;
    double        y2;
    double        x0;
    double        x1;
    int         code;

    if( y0 <= 0.0 || y0 >= 1.0 )
	{
        throw(   "Argument mut be in [0,1]." );
	}
                
    code = 1;
    y = y0;
    if( y > (1.0 - 0.13533528323661269189) ) { /* 0.135... = exp(-2) */
        y = 1.0 - y;
        code = 0;
    }
    if( y > 0.13533528323661269189 ) {
        y = y - 0.5;
        y2 = y * y;
        x = y + y * ( y2 * polynomial( y2, Mtools_nip0, 4 ) / polynomial_1( y2, Mtools_niq0, 8 ) );
        x = x * SQRT2PI;
        return( x );
    }
    x = sqrt( -2.0 * log( y ) );
    x0 = x - log( x ) / x;
    z = 1.0 / x;
    if( x < 8.0 ) /* y > exp( -32 ) = 1.2664165549e-14 */
        x1 = z * polynomial( z, Mtools_nip1, 8 ) / polynomial_1( z, Mtools_niq1, 8 );
    else
        x1 = z * polynomial( z, Mtools_nip2, 8 ) / polynomial_1( z, Mtools_niq2, 8 );
    x = x0 - x1;
    if( code != 0 )
        x = -x;
    return( x );
}

double normal_inv2( double u )
{
    static double a[4]={2.50662823884,
						-18.61500062529,
 						41.39119773534,
 						-25.44106049637};

    static double b[4]={-8.47351093090,
 						23.08336743743,
 						-21.06224101826,
 						3.13082909833};

    static double c[9]={0.3374754822726147,
 						0.9761690190917186,
 						0.1607979714918209,
 						0.0276438810333863,
 						0.0038405729373609,
 						0.0003951896511919,
 						0.0000321767881768,
 						0.0000002888167364,
 						0.0000003960315187};


    double x = 0;
    double r = 0;

    x = u - 0.5;
    if (fabs(x) < 0.42) 
    {
 	    r = x*x;
 	    r = x*(((a[3]*r+a[2])*r+a[1])*r+a[0])/(((( b[3]*r+b[2])*r+b[1])*r+b[0])*r+1.0);
     return(r);
    }
    r=u;
    if (x > 0.0) 
    {
 	    r=1.0-u;
    }
    r = log(-log(r));
    r=c[0]+r*(c[1]+r*(c[2]+r*(c[3]+r*(c[4]+r*(c[5]+r*(c[6]+r*(c[7]+r*c[8])))))));
    if (x < 0.0) 
 	    r = -r;

    return(r);              
}

static double Mtools_T[] = {
 9.60497373987051638749E0, 9.00260197203842689217E1, 2.23200534594684319226E3,
 7.00332514112805075473E3, 5.55923013010394962768E4};
static double Mtools_U[] = {/* 1.00000000000000000000E0,*/
 3.35617141647503099647E1, 5.21357949780152679795E2, 4.59432382970980127987E3,
 2.26290000613890934246E4, 4.92673942608635921086E4};

double normal_error( double x )
{
    double y;
    double z;

    if( fabs( x ) > 1.0 )
        return( 1.0 - normal_error_comp( x ) );
    z = x * x;
    y = x * polynomial( z, Mtools_T, 4 ) / polynomial_1( z, Mtools_U, 5 );
    return( y );
}

double normal_error_comp( double a )
{	
    double p;
    double q;
    double x;
    double y;
    double z;
	
    if( a < 0.0 )
        x = -a;
    else
        x = a;
    if( x < 1.0 )
        return( 1.0 - normal_error(a) );
    z = -a * a;        
    if( z < - MAXLOG )
        return( 0.0 );        
    z = exp( z );
    if( x < 8.0 ) {
        p = polynomial( x, Mtools_nep, 8 );
        q = polynomial_1( x, Mtools_neq, 8 );
    } else {
        p = polynomial( x, Mtools_ner, 5 );
        q = polynomial_1( x, Mtools_nes, 6 );
    }
    y = ( z * p ) / q;
    if( a < 0.0 )
        y = 2.0 - y;
    return( y );
}	
	
	
//--------------------------------------------------------------------------------------
// The following algorithms for the approximation of the bivariate and trivariate normal
// distribution are from 
// [1] Drezner,Wesolowsky, On the computation of the bivariate normal integral,
//    Journal of Statistical Computation and Simulation,35(1990),101-107
// and
// [2] Drezner, Computation of the trivariate nnormal integral,
//    Mathematics of Compuation,62(1994),289-294.
// with an improved Gauss quadrature using more points.
// The algorithm for the 4-variate normal distribution is developed along the lines of [2].
//--------------------------------------------------------------------------------------

///////////////////////////////////////////////
double normal( double x, double y, double rho )
///////////////////////////////////////////////    
{
    double bi;
    if( MlEqMaths::equal(rho,-1.0) ) 
        return( ( x > -y ? normal(x)-normal(-y) : 0.0 ) );
    if( MlEqMaths::equal(rho,+1.0) )
        return( normal(MlEqMaths::Min(x,y)) );
    if( fabs(rho) > 1.0 )
	{
        throw(   "rho must be in [-1.0,1.0]" );
	}
    
    if( rho * rho < 0.5 )
        bi = drezn( x, y, rho );
    else  {
        double rho1 = sqrt(1.0 - rho * rho);       
        double u = ( y - rho * x ) / rho1;
        if( rho >= 0.0 )
            bi = normal( x ) * normal( u ) + normal( y ) - drezn( u, y, rho1 );
        else
            bi = -1.0 * normal( ( -1.0 * x ) ) * normal( u ) + drezn( u , y , rho1 );
    }      
    return( bi );
} 
///////////////////////////////////////////////////////
inline static  double Z2(  double x,  double y,  double r )
///////////////////////////////////////////////////////
{   
    double a = 1-r*r;
    return( exp((-x*x+2*x*y*r-y*y)/2/a)/sqrt(a) );   
}
/////////////////////////////////////////////////////
static  double drezn(  double x,  double y,  double rho )
/////////////////////////////////////////////////////     
{
    double sum = 0.0;    
    for( int i = 0; i < 20; i++ )           
         sum += w20[i]*Z2(x,y,x20[i]*rho);        
    sum = normal(x)*normal(y)+rho*sum;
    return(sum);
}
/////////////////////////////////////////////////////////////////////////////////////////////////
static  double Z3(  double t,  double h1,  double h2,  double h3,  double r12,  double r13,  double r23 )
/////////////////////////////////////////////////////////////////////////////////////////////////
{
    double result, x, y, sdet, a2, a3;
   sdet = sqrt(1-t*t*r12*r12-t*t*r13*r13-r23*r23+2*t*t*r12*r13*r23);
   a2 = 1-t*t*r12*r12;
   a3 = 1-t*t*r13*r13;
   x = h3*a2+(r23*r12-r13)*t*h1+(t*t*r13*r12-r23)*h2;
   y = h2*a3+(r23*r13-r12)*t*h1+(t*t*r13*r12-r23)*h3;   
   x /= sqrt(a2)*sdet;
   y /= sqrt(a3)*sdet;
   result = r12*Z2(h1,h2,t*r12)*normal(x) + r13*Z2(h1,h3,t*r13)*normal(y);
   return( result );
}

void swap( double &a, double &b )
{
	double z = a;
	a = b;
	b = z;
}
////////////////////////////////////////////////////////////////////////////////////
double normal( double h1, double h2, double h3, double r12, double r13, double r23 )
////////////////////////////////////////////////////////////////////////////////////
{ 
    if( MlEqMaths::equal(r23,-1.0) ) {
        if( !MlEqMaths::equal(r12,-r13) ) 
		{
            throw(   "Invalid correlation matrix." );
		}
        return( h2 > -h3 ? normal(h1,h2,r12) - normal(h1,-h3,r12) : 0.0 );
    }
    if( MlEqMaths::equal(r23,+1.0) ) {
        if( !MlEqMaths::equal(r12,r13) ) 
		{
            throw(   "Invalid correlation matrix." );
		}
        return( normal(h1, MlEqMaths::Min(h2,h3),r12) );
    }
    if( MlEqMaths::equal(r12,-1.0) ) {
        if( !MlEqMaths::equal(r13,-r23) ) 
		{
            throw(   "Invalid correlation matrix." );
		}
        return( h1 > -h2 ? normal(h1,h3,r13) - normal(-h2,h3,r13) : 0.0 );
    }
    if( MlEqMaths::equal(r12,+1.0) ) {
        if( !MlEqMaths::equal(r13,r23) ) 
		{
            throw(   "Invalid correlation matrix." );
		}
        return( normal(MlEqMaths::Min(h1,h2),h3,r13) );
    }
    if( MlEqMaths::equal(r13,-1.0) ) {
        if( !MlEqMaths::equal(r12,-r23) ) 
		{
            throw(   "Invalid correlation matrix." );
		}
        return( h1 > -h3 ? normal(h1,h2,r12) - normal(-h3,h2,r12) : 0.0 );
    }
    if( MlEqMaths::equal(r13,1.0) ) {
        if( !MlEqMaths::equal(r12,r23) ) 
		{
            throw(   "Invalid correlation matrix." );
		}
        return( normal(MlEqMaths::Min(h1,h3),h2,r12) );
    }
    double det = 1-r12*r12-r13*r13-r23*r23+2*r12*r13*r23;
    if( det < 0.0 )
	{
        throw(   "Invalid correlation matrix." );
	}
	
    if( fabs(r13) > fabs(r23) ) {
        swap(h1,h2);
        swap(r13,r23);
    }
    if( fabs(r12) > fabs(r23) ) {
        swap(h1,h3);
        swap(r12,r23);
    }
	
    double sum = 0.0;
    for( int i = 0; i < 20; i++ )          
         sum += w20[i]*Z3(x20[i],h1,h2,h3,r12,r13,r23);    
    sum += normal(h1)*normal(h2,h3,r23);    
    return( sum );
}	
/////////////////////////////////////////////////////////////////////////////////////////////////
static  double Z4(  double t,  double h1,  double h2,  double h3,  double h4, 
                        double r12,  double r13,  double r14,  double r23,  double r24,  double r34 )
/////////////////////////////////////////////////////////////////////////////////////////////////
{	
    double result, x2, x3, x4, y2, y3, y4, det, c33, c34, c44, d33, d34, d44, e33, e34, e44, rho2, rho3, rho4;
   det = 1 - t*t*r12*r12 - t*t*r13*r13 - t*t*r14*r14 - r23*r23 - r24*r24 - r34*r34
           + t*t*r12*r12*r34*r34 + t*t*r13*r13*r24*r24 + t*t*r14*r14*r23*r23
           + 2*t*t*r12*r13*r23 + 2*r23*r24*r34 + 2*t*t*r12*r14*r24 + 2*t*t*r13*r14*r34
           - 2*t*t*r12*r14*r23*r34 - 2*t*t*r12*r13*r24*r34 - 2*t*t*r13*r14*r23*r24 ;
	
   c33 = ( 1 + 2*t*t*r12*r14*r24 - t*t*r12*r12 - t*t*r14*r14 - r24*r24 ) / det;
   c34 = ( -r34 - t*t*r12*r13*r24 - t*t*r12*r14*r23 + t*t*r13*r14 + r23*r24 + t*t*r12*r12*r34 ) / det;
   c44 = ( 1 + 2*t*t*r12*r13*r23 - t*t*r12*r12 - t*t*r13*r13 - r23*r23 ) / det;
   rho2 = -c34/sqrt(c33)/sqrt(c44);
	
   d33 = ( 1 + 2*t*t*r13*r14*r34 - t*t*r13*r13 - t*t*r14*r14 - r34*r34 ) / det;
   d34 = ( -r24 - t*t*r13*r12*r34 - t*t*r13*r14*r23 + t*t*r12*r14 + r23*r34 + t*t*r13*r13*r24 ) / det;
   d44 = c44;
   rho3 = -d34/sqrt(d33)/sqrt(d44);
	
   e33 = c33;
   e34 = ( -r23 - t*t*r14*r13*r24 - t*t*r12*r14*r34 + t*t*r13*r12 + r34*r24 + t*t*r14*r14*r23 ) / det;
   e44 = ( 1 + 2*t*t*r14*r13*r34 - t*t*r14*r14 - t*t*r13*r13 - r34*r34 ) / det;
   rho4 = -e34/sqrt(e33)/sqrt(e44);   
	
   x2 = sqrt(c33)*sqrt(1-rho2*rho2)*( h3*(1-t*t*r12*r12) + (r23*r12-r13)*t*h1 + (t*t*r13*r12-r23)*h2 ) / (1-t*t*r12*r12);
   y2 = sqrt(c44)*sqrt(1-rho2*rho2)*( h4*(1-t*t*r12*r12) + (r24*r12-r14)*t*h1 + (t*t*r14*r12-r24)*h2 ) / (1-t*t*r12*r12);
	
   x3 = sqrt(d33)*sqrt(1-rho3*rho3)*( h2*(1-t*t*r13*r13) + (r23*r13-r12)*t*h1 + (t*t*r13*r12-r23)*h3 ) / (1-t*t*r13*r13);
   y3 = sqrt(d44)*sqrt(1-rho3*rho3)*( h4*(1-t*t*r13*r13) + (r34*r13-r14)*t*h1 + (t*t*r14*r13-r34)*h3 ) / (1-t*t*r13*r13);
	
   x4 = sqrt(e33)*sqrt(1-rho4*rho4)*( h3*(1-t*t*r14*r14) + (r34*r14-r13)*t*h1 + (t*t*r13*r14-r34)*h4 ) / (1-t*t*r14*r14);
   y4 = sqrt(e44)*sqrt(1-rho4*rho4)*( h2*(1-t*t*r14*r14) + (r24*r14-r12)*t*h1 + (t*t*r12*r14-r24)*h4 ) / (1-t*t*r14*r14);
	
   result = r12*Z2(h1,h2,t*r12)*normal(x2,y2,rho2) + r13*Z2(h1,h3,t*r13)*normal(y3,y3,rho3) + r14*Z2(h1,h4,t*r14)*normal(y4,y4,rho4);
   return( result );
}	
////////////////////////////////////////////////////////////////////////////////////
double normal( double h1, double h2, double h3, double h4, double r12, double r13, double r14, double r23, double r24, double r34 )
////////////////////////////////////////////////////////////////////////////////////
{	
    double det = 1 - r12*r12 - r13*r13 - r14*r14 - r23*r23 - r24*r24 - r34*r34
                   + r12*r12*r34*r34 + r13*r13*r24*r24 + r14*r14*r23*r23
                   + 2*r12*r13*r23 + 2*r23*r24*r34 + 2*r12*r14*r24 + 2*r13*r14*r34
                   - 2*r12*r14*r23*r34 - 2*r12*r13*r24*r34 - 2*r13*r14*r23*r24 ;
    if( det < 0.0 )
        throw(   "Invalid correlation matrix." );
	
    double sum = 0.0;
    for( int i = 0; i < 40; i++ )          
         sum += w40[i]*Z4(x40[i],h1,h2,h3,h4,r12,r13,r14,r23,r24,r34);    
    sum += normal(h1)*normal(h2,h3,h4,r23,r24,r34);    
    return( sum );
}	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Global constants for Gauss quadrature, from Abramowitz,Stegun, Pocketbook of Mathematical functions, Table 25.4.
// The weights and abscissas are adapted for the integral \int_0^1 f(s) ds / (2\pi) .
// Notation: xNN[], wNN[] where x are the abscissas, w are the weights and NN is the number of them.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 double x20[] = {   
3.4357004074525580e-003, 
1.8014036361043100e-002, 
4.3882785874337030e-002, 
8.0441514088890610e-002, 
1.2683404676992460e-001,                      
1.8197315963674250e-001, 
2.4456649902458640e-001, 
3.1314695564229020e-001, 
3.8610707442917750e-001, 
4.6173673943325130e-001, 
5.3826326056674870e-001, 
6.1389292557082250e-001, 
6.8685304435770980e-001, 
7.5543350097541360e-001, 
8.1802684036325760e-001,                       
8.7316595323007540e-001, 
9.1955848591110940e-001, 
9.5611721412566300e-001, 
9.8198596363895700e-001, 
9.9656429959254740e-001 };
 double w20[] = {   
1.4016781519259970e-003, 
3.2309591246650840e-003, 
4.9872831430338210e-003, 
6.6269525332594600e-003, 
8.1113412094314890e-003, 
9.4056220040543340e-003, 
1.0479448879113940e-002, 
1.1307649096073420e-002, 
1.1870809086447660e-002, 
1.2155728317942470e-002, 
1.2155728317942470e-002, 
1.1870809086447660e-002, 
1.1307649096073420e-002, 
1.0479448879113940e-002, 
9.4056220040543340e-003, 
8.1113412094314890e-003, 
6.6269525332594600e-003, 
4.9872831430338210e-003, 
3.2309591246650840e-003, 
1.4016781519259970e-003 };
 double x40[] = {
8.811451447203744e-004,
4.636880650271513e-003,
1.137002500811285e-002,
2.104159039310416e-002,
3.359359586066174e-002,
4.895059651556283e-002,
6.702024839387027e-002,
8.769388458334415e-002,
1.108471742867403e-001,
1.363408724050365e-001,
1.640216576929102e-001,
1.937230551660099e-001,
2.252664374524359e-001,
2.584620991569107e-001,
2.931103978141975e-001,
3.290029545871208e-001,
3.659239074963732e-001,
4.036512096493145e-001,
4.419579646623724e-001,
4.806137912469746e-001,
5.193862087530254e-001,
5.580420353376276e-001,
5.963487903506856e-001,
6.340760925036268e-001,
6.709970454128793e-001,
7.068896021858024e-001,
7.415379008430894e-001,
7.747335625475641e-001,
8.062769448339902e-001,
8.359783423070898e-001,
8.636591275949636e-001,
8.891528257132597e-001,
9.123061154166559e-001,
9.329797516061298e-001,
9.510494034844372e-001,
9.664064041393383e-001,
9.789584096068958e-001,
9.886299749918872e-001,
9.953631193497285e-001,
9.991188548552796e-001 };
 double w40[] = {
3.597917996598698e-004,
8.354269385590755e-004,
1.306746306140621e-003,
1.770268431264264e-003,
2.223156378031754e-003,
2.662677738018804e-003,
3.086186868478710e-003,
3.491135947840297e-003,
3.875089246486902e-003,
4.235737478816289e-003,
4.570911591112636e-003,
4.878595767570021e-003,
5.156939536905879e-003,
5.404268893488870e-003,
5.619096361092614e-003,
5.800129936683468e-003,
5.946280859533507e-003,
6.056670158498864e-003,
6.130633937870425e-003,
6.167727369894800e-003,
6.167727369894800e-003,
6.130633937870425e-003,
6.056670158498864e-003,
5.946280859533507e-003,
5.800129936683468e-003,
5.619096361092614e-003,
5.404268893488870e-003,
5.156939536905879e-003,
4.878595767570021e-003,
4.570911591112636e-003,
4.235737478816289e-003,
3.875089246486902e-003,
3.491135947840297e-003,
3.086186868478710e-003,
2.662677738018804e-003,
2.223156378031754e-003,
1.770268431264264e-003,
1.306746306140621e-003,
8.354269385590755e-004,
3.597917996598698e-004 };




