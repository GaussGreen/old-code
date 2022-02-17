#include	"stdafx.h"

#include "MCInt\IntNormal.h"


//----------------------------------------------------------------------------------------------

// Normal density function
double Ndf(double x)
{
  return 0.398942280401433*exp(-0.5*x*x);
}

//----------------------------------------------------------------------------------------------

// Evaluate polynomial C_0 + C_1*x + ... + C_N*x^N
// Note: coefficients are stored in reverse order: coef[i] = C_(N-i)
double polevl(double x, double coef[], int N)
{
  double ans;
  int i;
  double *p;

  p = coef;
  ans = *p++;
  i = N;

  do
    ans = ans*x  +  *p++;
  while (--i);

  return ans;
}

//----------------------------------------------------------------------------------------------

// Evaluate polynomial C_0 + C_1*x + ... + C_(N-1)*x^(N-1) + x^N
// Note: p1evl() assumes that coef[N] = 1.0, which is omitted from the array
double p1evl(double x, double coef[], int N)
{
  double ans;
  double *p;
  int i;

  p = coef;
  ans = x + *p++;
  i = N-1;

  do
    ans = ans*x  + *p++;
  while (--i);

  return ans;
}

//----------------------------------------------------------------------------------------------

// Error function
// From ndtr.c, Cephes Math Library Release 2.1: January, 1989
// http://www.netlib.org/cephes/
// http://scipy.net/cgi-bin/viewcvsx.cgi/scipy/special/cephes/
double erf(double x)
{
  static double T[] = {  9.60497373987051638749E0,   9.00260197203842689217E1, 2.23200534594684319226E3,
			                   7.00332514112805075473E3,   5.55923013010394962768E4 };
  static double U[] = {/*1.00000000000000000000E0,*/ 3.35617141647503099647E1, 5.21357949780152679795E2,
		                     4.59432382970980127987E3,   2.26290000613890934246E4, 4.92673942608635921086E4 };
  double y, z;

  if (fabs(x) > 1.0)
    return 1.0 - erfc(x);
        
  z = x * x;
  y = x * polevl(z, T, 4) / p1evl(z, U, 5);

  return y;
}

//----------------------------------------------------------------------------------------------

// Complementary error function
// From ndtr.c, Cephes Math Library Release 2.1: January, 1989
// http://www.netlib.org/cephes/
// http://scipy.net/cgi-bin/viewcvsx.cgi/scipy/special/cephes/
double erfc(double a)
{
  // log(2**127)
  static double MAXLOG =  8.8029691931113054295988E1;

  static double P[] = {  2.46196981473530512524E-10, 5.64189564831068821977E-1, 7.46321056442269912687E0,
			                   4.86371970985681366614E1,   1.96520832956077098242E2,  5.26445194995477358631E2,
			                   9.34528527171957607540E2,   1.02755188689515710272E3,  5.57535335369399327526E2 };
  static double Q[] = {/*1.00000000000000000000E0,*/ 1.32281951154744992508E1,  8.67072140885989742329E1,
		                     3.54937778887819891062E2,   9.75708501743205489753E2,  1.82390916687909736289E3,
		                     2.24633760818710981792E3,   1.65666309194161350182E3,  5.57535340817727675546E2 };
  static double R[] = {  5.64189583547755073984E-1,  1.27536670759978104416E0,  5.01905042251180477414E0,
			                   6.16021097993053585195E0,   7.40974269950448939160E0,  2.97886665372100240670E0 };
  static double S[] = {/*1.00000000000000000000E0,*/ 2.26052863220117276590E0,  9.39603524938001434673E0,
		                     1.20489539808096656605E1,   1.70814450747565897222E1,  9.60896809063285878198E0,
		                     3.36907645100081516050E0 };
  double p,q,x,y,z;

  if (a < 0.0)
    x = -a;
  else
    x = a;

  if (x < 1.0)
    return 1.0 - erf(a);

  z = - a * a;

  if (z < -MAXLOG)
    {
      if (a < 0)
	return 2.0;
      else
	return 0.0;
    }

  z = exp(z);

  if (x < 8.0)
    {
      p = polevl(x, P, 8);
      q = p1evl(x, Q, 8);
    }
  else
    {
      p = polevl(x, R, 5);
      q = p1evl(x, S, 6);
    }
  y = (z * p)/q;

  if (a < 0)
    y = 2.0 - y;

  if (y == 0.0)
    {
      if (a < 0)
	return 2.0;
      else
	return 0.0;
    }

  return y;
}

//----------------------------------------------------------------------------------------------

// Normal cumulative density function, Implementation 1
// From "nc" in http://www.mathfinance.de/FF/cpplib.html
// Based on a polynomial approximation given in M. Abramowitz and I. Stegun, 
// Handbook of Mathematical Functions, p. 932, eqn. 26.2.17
double Ncdf1(double x)
{
  static double a[5] = {0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429 };

  double res;

  if (x < -7.0)
    {
      res = Ndf(x)/sqrt(1.0 + x*x);
    }
  else if (x > 7.0)
    {
      res = 1.0 - Ndf(x)/sqrt(1.0 + x*x);
    }
  else
    {
      res = 0.2316419;
      res = 1.0/(1.0 + res*fabs(x));
      res = Ndf(x)*(res*(a[0] + res*(a[1] + res*(a[2] + res*(a[3] + res*a[4])))));
                
      if (x > 0) 
	res = 1.0 - res;
    }
        
  return res;
}

//----------------------------------------------------------------------------------------------

// Normal cumulative density function, Implementation 2
// From "ndtr" in ndtr.c, Cephes Math Library Release 2.1: January, 1989
// http://www.netlib.org/cephes/
// http://scipy.net/cgi-bin/viewcvsx.cgi/scipy/special/cephes/
double Ncdf2(double a)
{
  // sqrt(2)/2
  static double SQRTH = 7.07106781186547524401E-1;

  double x, y, z;

  x = a * SQRTH;
  z = fabs(x);

  if (z < SQRTH)
    {
      y = 0.5 + 0.5 * erf(x);
    }
  else
    {
      y = 0.5 * erfc(z);

      if (x > 0)
	y = 1.0 - y;
    }

  return y;
}

//----------------------------------------------------------------------------------------------

// Normal inverse cumulative density function, Implementation 1
// From Boris Moro, The Full Monte, Union Bank of Switzerland, RISK 1995(2)
double Nicdf1(double u)
{
  static double a[4] = {2.50662823884, -18.61500062529, 41.39119773534, -25.44106049637};
    
  static double b[4] = {-8.47351093090, 23.08336743743, -21.06224101826, 3.13082909833};
    
  static double c[9] = {0.3374754822726147, 0.9761690190917186, 0.1607979714918209,
			                  0.0276438810333863, 0.0038405729373609, 0.0003951896511919,
			                  0.0000321767881768, 0.0000002888167364, 0.0000003960315187};

  if ((u <= 0.0) || (u >= 1.0)) 
  {
	  string err = "In Nicdf1: u = ";
	  err += (double)u;
	  err = err + " outside of range (0, 1). ";
	  ERROR(err)
  }

  double x, r;

  x = u - 0.5;    

  if (fabs(x) < 0.42) 
    {
      r = x*x;
      r = x*(((a[3]*r + a[2])*r + a[1])*r + a[0]) / ((((b[3]*r + b[2])*r + b[1])*r + b[0])*r + 1.0);
        
      return r;
    }

  r = u;

  if (x > 0.0) 
    {
      r = 1.0 - u;
    }

  r = log(-log(r));
  r = c[0] + r*(c[1] + r*(c[2] + r*(c[3] + r*(c[4] + r*(c[5] + r*(c[6] + r*(c[7] + r*c[8])))))));
        
  if (x < 0.0) 
    {
      r = -r;
    }

  return r;              
}

//----------------------------------------------------------------------------------------------

// Normal inverse cumulative density function, Implementation 2
// From "ndtri" in ndtri.c, Cephes Math Library Release 2.1: January, 1989
// http://www.netlib.org/cephes/
// http://scipy.net/cgi-bin/viewcvsx.cgi/scipy/special/cephes/
double Nicdf2(double u)
{
  static double s2pi = 2.50662827463100050242E0;

  // approximation for 0 <= |y - 0.5| <= 3/8
  static double P0[5] = { -5.99633501014107895267E1,   9.80010754185999661536E1,  -5.66762857469070293439E1,
			                     1.39312609387279679503E1,  -1.23916583867381258016E0, };
  static double Q0[8] = {/*1.00000000000000000000E0,*/ 1.95448858338141759834E0,   4.67627912898881538453E0,
			                     8.63602421390890590575E1,  -2.25462687854119370527E2,   2.00260212380060660359E2,
			                    -8.20372256168333339912E1,   1.59056225126211695515E1,  -1.18331621121330003142E0, };

  // Approximation for interval z = sqrt(-2 log y ) between 2 and 8 
  // i.e., y between exp(-2) = .135 and exp(-32) = 1.27e-14.
  static double P1[9] = {  4.05544892305962419923E0,   3.15251094599893866154E1,   5.71628192246421288162E1,
			                     4.40805073893200834700E1,   1.46849561928858024014E1,   2.18663306850790267539E0,
			                    -1.40256079171354495875E-1, -3.50424626827848203418E-2, -8.57456785154685413611E-4, };
  static double Q1[8] = {/*1.00000000000000000000E0,*/ 1.57799883256466749731E1,   4.53907635128879210584E1,
		                       4.13172038254672030440E1,   1.50425385692907503408E1,   2.50464946208309415979E0,
		                     - 1.42182922854787788574E-1, -3.80806407691578277194E-2, -9.33259480895457427372E-4, };

  // Approximation for interval z = sqrt(-2 log y ) between 8 and 64
  // i.e., y between exp(-32) = 1.27e-14 and exp(-2048) = 3.67e-890.
  static double P2[9] = {  3.23774891776946035970E0,   6.91522889068984211695E0,   3.93881025292474443415E0,
			                     1.33303460815807542389E0,   2.01485389549179081538E-1,  1.23716634817820021358E-2,
			                     3.01581553508235416007E-4,  2.65806974686737550832E-6,  6.23974539184983293730E-9, };
  static double Q2[8] = {/*1.00000000000000000000E0,*/ 6.02427039364742014255E0,   3.67983563856160859403E0,
			                     1.37702099489081330271E0,   2.16236993594496635890E-1,  1.34204006088543189037E-2,
			                     3.28014464682127739104E-4,  2.89247864745380683936E-6,  6.79019408009981274425E-9, };

  double x, y, z, y2, x0, x1;
  int code;

  if ((u <= 0.0) || (u >= 1.0))
  {
	  string err = "In Nicdf1: u = ";
	  err += u;
	  err += " outside of range (0, 1). ";	
	  ERROR(err)
  }

  code = 1;
  y = u;

  // 0.135... = exp(-2)
  if(y > (1.0 - 0.13533528323661269189)) 
    {
      y = 1.0 - y;
      code = 0;
    }

  if(y > 0.13533528323661269189)
    {
      y = y - 0.5;
      y2 = y * y;
      x = y + y * (y2 * polevl(y2, P0, 4) / p1evl(y2, Q0, 8));
      x = x * s2pi; 
      return x;
    }

  x = sqrt(-2.0 * log(y));
  x0 = x - log(x)/x;

  z = 1.0/x;

  // y > exp(-32) = 1.2664165549e-14
  if(x < 8.0) 
    x1 = z * polevl(z, P1, 8) / p1evl(z, Q1, 8);
  else
    x1 = z * polevl(z, P2, 8) / p1evl(z, Q2, 8);

  x = x0 - x1;

  if(code != 0)
    x = -x;

  return x;
}

//----------------------------------------------------------------------------------------------

