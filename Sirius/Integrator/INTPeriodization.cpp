#include	"stdafx.h"


#include "MCInt\INTPeriodization.h"


static double xMin(double x,double y)
{
	if ( x < y ) return x;
	else return y;
}

static double xMax(double x,double y)
{
	if ( x < y ) return y;
	else return x;
}

	//----------------------------------------------------------------------------------------------

// Plot phi(x), or phi1(x), or (phi(x), phi1(x)), or phi(x), ..., phi1(x), ...
// (according to choice 0, 1, 2 or 3) vs x for given periodizations
void CPeriodization::plot(const vector<CPeriodization*>& p, int choice, const string& filename, int res, int margin)
{
  // Number of rows = evaluations of phi(x)
  int R = res - 2*margin;
  int r;

  // Number of colums = number of phis to plot
  int C = p.size();
  int Ceff = (choice >= 2) ? 2*C : C;
  int c;

  // Row and column labels and table entries
  STLDoubleVector X(R, 0);
  vector<const char*> Y(Ceff, (const char*)NULL);
  STLDoubleVectorVector Z(R, STLDoubleVector(Ceff, 0.0));

  // Validate input
  for (c = 0; c < C; c++) if (!p[c]) ERROR("In CPeriodization::plotphi: periodization invalid. ")
  if ((choice < 0) || (choice > 3)) ERROR("In CPeriodization::plotphi: choice must be in {0, 1, 2, 3}. ")

  // Fill in table

  for (r = 0; r < R; r++)
    {
      X[r] = (double)(r + margin)/(double)(res - 1);

      for (c = 0; c < C; c++)
        {
          switch(choice)
            {
            case 0:
              Z[r][c] = p[c]->phi(X[r]);
              break;
            case 1:
              Z[r][c] = p[c]->phi1(X[r]);
              break;
            case 2:
              Z[r][2*c] = p[c]->phi(X[r]);
              Z[r][2*c+1] = p[c]->phi1(X[r]);
              break;
            case 3:
              Z[r][c] = p[c]->phi(X[r]);
              Z[r][C+c] = p[c]->phi1(X[r]);
            }
        }
    }

  // Store column labels
  vector<string> label(Ceff, "");
  for (c = 0; c < C; c++) 
    {
      switch(choice)
        {
        case 0:
          label[c] = string(p[c]->acronym()) + ".phi(x)";
          Y[c] = label[c].c_str();
          break;
        case 1: 
          label[c] = string(p[c]->acronym()) + ".phi1(x)";
          Y[c] = label[c].c_str();
          break;
        case 2:
          label[2*c] = string(p[c]->acronym()) + ".phi(x)";
          label[2*c+1] = string(p[c]->acronym()) + ".phi1(x)";
          Y[2*c] = label[2*c].c_str();
          Y[2*c+1] = label[2*c+1].c_str();
          break;
        case 3:
          label[c] = string(p[c]->acronym()) + ".phi(x)";
          label[C+c] = string(p[c]->acronym()) + ".phi1(x)";
          Y[c] = label[c].c_str();
          Y[C+c] = label[C+c].c_str();
        }
    }

  // Store table title
  string title;
  switch(choice)
    {
    case 0:
      title = string("x \\ p.phi(x)");
      break;
    case 1:
      title = string("x \\ p.phi1(x)");
      break;
    case 2:
      title = string("x \\ (p.phi, p.phi1(x))");
      break;
    case 3:
      title = string("x \\ p.phi(x) resp. p.phi1(x)");
      break;
    }      

  plotXYZ(X, Y, Z, title.c_str(), filename);
}

// Plot phi(x), or phi1(x), or (phi(x), phi1(x)), or phi(x), ..., phi1(x), ...
// (according to choice 0, 1, 2 or 3) vs x for given periodizations
void CPeriodization::plot(int choice, const string& filename, int res, int margin)
{
  // Generate filename if required
  string fn = (filename.size() == 0) ? INFOPATH + acronym() + "-plot.txt" : filename;

  vector<CPeriodization*> p(1, this);

  plot(p, choice, fn, res, margin);
}

//----------------------------------------------------------------------------------------------

// Change of variables mapping |R^k -> |R^k
STLDoubleVector& CPeriodization::phi(STLDoubleVector& x) const
{
  int k = x.size();

  for (int i = 0; i < k; i++)
    {
      x[i] = xMin( xMax(phi(x[i]), 0.000001), 0.999999);
    }

  return x;
}

//----------------------------------------------------------------------------------------------

double CPeriodization::product(const STLDoubleVector& x) const
{
  int k = x.size();
  double result = 1.0;

  for (int i = 0; i < k; i++)
    {
      result *= phi1(x[i]);
    }

  return result;
}

//----------------------------------------------------------------------------------------------

CPolynomialperiodization::CPolynomialperiodization(int d) : CFinitetypeperiodization(d)
{
  // Note: (2*d+1)!/(d!)^2 = (d+1)*binomial(2*d+1, d)
  _c = (double)(d + 1)*binomial(2*d + 1, d);

  _S = string("p-poly(") + _d + ")";
}

//----------------------------------------------------------------------------------------------

double CPolynomialperiodization::phi(double x) const
{
  double sum = 0.0;

  for (int i = 0; i <= _d; i++)
    {
      sum += A(i)*pow(x, _d + i + 1);
    }

  return _c*sum;
}

//----------------------------------------------------------------------------------------------

CTrigonometricperiodization::CTrigonometricperiodization(int d) : CFinitetypeperiodization(d)
{
  if (d%2 == 1)
    {
     _i_max = (d - 1)/2;
     _c = 0.5*PI*specialbinomial(d);
    }
  else
    {
      _i_max = d/2-1;
      _c = specialbinomial(d);
    }

  _S = string("p-trig(") + _d + ")";
}

//----------------------------------------------------------------------------------------------

/*
double CTrigonometricperiodization::phi(double x) const
{
  double result;
  double sum = 0.0;
  double sinPix2 = pow(sin(PI*x), 2);

  if (d%2 == 1)
    {
      for (int i = 0; i <= _i_max; i++)
        {
          sum += 1.0/specialbinomial(2*i)*pow(sinPix2, i);
        }
       result = 0.5*(1.0 - cos(PI*x)*sum);
    }
  else
    {
      for (int i = 0; i <= _i_max; i++)
        {
          sum += 1.0/specialbinomial(2*i + 1)*pow(sinPix2, i);
        }
      result = x - 1.0/PI*cos(PI*x)*sin(PI*x)*sum;
    }

  return result;
}
*/

//----------------------------------------------------------------------------------------------

double CTrigonometricperiodization::phi(double x) const
{
  double result;
  double sum = 0.0;
  double summand = 1.0;
  double sinPix2 = pow(sin(PI*x), 2);

  if (_d%2 == 1)
    {
      for (int i = 0; i <= _i_max; i++)
        {
          sum += summand;
          summand *= (double)(2*(i + 1) - 1)/(double)(2*(i + 1)) * sinPix2;
        }
      result = 0.5*(1.0 - cos(PI*x)*sum);
    }
  else
    {
      for (int i = 0; i <= _i_max; i++)
        {
          sum += summand;
          summand *= (double)2*(i + 1)/(double)(2*(i + 1) + 1) * sinPix2;
        }
      result = x - 1.0/PI*cos(PI*x)*sin(PI*x)*sum;
    }

  return result;
}

//----------------------------------------------------------------------------------------------

double CKressperiodization::phi(double x) const
{
  double xdplus1 = pow(x, _d + 1);

  return xdplus1 / (xdplus1 + pow(1.0 - x, _d + 1));
}

//----------------------------------------------------------------------------------------------

double CKressperiodization::phi1(double x) const
{
  double xd = pow(x, _d);
  double oneminusxd = pow(1.0 - x, _d);

  return (double)(_d + 1) * xd * oneminusxd / pow(xd*x + oneminusxd*(1.0 - x), 2);
}

//----------------------------------------------------------------------------------------------

double CTANHperiodization::phi(double x) const 
{ 
  return 0.5*(tanh(- 0.5*_a*(oneby(pow(x, _p)) - oneby(pow(1.0 - x, _p)))) + 1.0); 
}

//----------------------------------------------------------------------------------------------

double CTANHperiodization::phi1(double x) const
{ 
  double xp = pow(x, _p);
  double oneminusxp = pow(1.0 - x, _p);

  double num = 0.25*_a*_p*(oneby(xp*x) + oneby(oneminusxp*(1.0 - x)));
  double denom = pow(cosh(0.5*_a*(oneby(xp) - oneby(oneminusxp))), 2);

  return num/denom;
}

//----------------------------------------------------------------------------------------------

double CDEperiodization::phi(double x) const
{
  return 0.5*(tanh(0.5*PI*sinh(_a*(2.0*x - 1.0))) + 1.0);
}

//----------------------------------------------------------------------------------------------

double CDEperiodization::phi1(double x) const
{
  double num = 0.5*PI*_a*cosh(_a*(2.0*x - 1.0));
  double denom = pow(cosh(0.5*PI*sinh(_a*(2.0*x - 1.0))), 2);

  return num/denom;
}

//----------------------------------------------------------------------------------------------
