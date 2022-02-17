#include	"stdafx.h"




#include "MCInt\INTIntegrand.h"
#include "MCInt\INTEnumerator.h"
#include "MCInt\INTPeriodization.h"



//----------------------------------------------------------------------------------------------

// Plot integrand of one or two variables to file after domain transformation
// (except rectilinear domain), manually generating mesh of specified resolution
void CIntegrand::plot(const vector<CIntegrand*>& f, const CDomain& D, 
                      const string& filename, int res, int margin, bool excel)
{  
  int reseff = res - 2*margin;
  
  // Number of rows = function evaluations
  int R = pow(reseff, D._k);
  int r = 0;

  // Number of colums = integrands to plot
  int C = f.size();
  int c;

  // Validate input
  if (excel && (D._k == 2) && (C > 1)) ERROR("In CIntegrand::plot: 3D-plots impossible for multiple functions. ") 
  for (c = 0; c < C; c++) if (!f[c]) ERROR("In CIntegrand::plot: integrand invalid. ")

  // Row and column labels and table entries
  vector<STLDoubleVector> X(R, STLDoubleVector(D._k, 0.0));
  vector<const char*> Y(C);
  STLDoubleVectorVector Z(R, STLDoubleVector(C, 0.0));

  // Fill in table

  // Enumerator for cubic lattice of specified resolution on domain D
  CEnumerator e(D._k, res, margin);

  // Current sample point
  STLDoubleVector x(D._k);

  do
    {
      e.xprop(x);

      // Check if the geometry of the domain D requires to plot f(psi(x))*J(x) vs x in [0, 1)^k 
      // (rather than f(y) vs y in D). True by default except for rectilinear domains.
      if (D.transform_f_to_plot())
        {
          // Abscissa x in [0, 1)^k
          X[r] = x;       

          // Change of variables y = psi(x)
          double Jx = D.J(x);
          STLDoubleVector& y = D.psi(x);

          // Evaluate f(psi(x))*J(x) for each f
          for (c = 0; c < C; c++)
            {           
              Z[r][c] = f[c]->evaluate(y) * Jx;
            }
        }
      else
        {
          // Change of variables y = psi(x)
          STLDoubleVector& y = D.psi(x);

          // Abscissa y in D
          X[r] = y;

          // Evaluate f(y) for each f
          for (c = 0; c < C; c++)
            {           
              Z[r][c] = f[c]->evaluate(y);
            }
        }

      r++;
    } 
  while (++e);

  if (excel && (D._k == 2))
    {
      // 3D-plot of single function f_0 : |R^2 -> |R      

      // The first column gives the x co-ordinate, the first row gives the y co-ordinate.
      // Field (x, y) holds f(x, y).

      // Store column labels y and relabel rows from (x, y) to x
      C = reseff;
      STLDoubleVector X1(C), Y1(C);
      for (c = 0; c < C; c++) 
        {
          X1[c] = Y1[c]= X[c][D._k - 1];
        }

      // Rearrange table entries
      STLDoubleVectorVector Z1(C, STLDoubleVector(C, 0.0));
      int i = 0;
      for (r = 0; r < C; r++)
        {
          for (c = 0; c < C; c++)
            {
              Z1[r][c] = Z[i++][0];
            }
        }

      // Store table title
      string title = string("x \\ y for ") + f[0]->acronym();

      plotXYZ(X1, Y1, Z1, title.c_str(), filename);
    }
  else 
    {
      // Store column labels
      for (c = 0; c < C; c++) 
        {
          if (!f[c]) ERROR("In CIntegrand::plot: integrand invalid. ")
          Y[c] = f[c]->acronym();
        }

      if (D._k == 1)
        {
          // 1D-plot of (multiple) function(s) f(_i) : |R -> |R      

          // Relabel rows from (x) to x
          STLDoubleVector X1(R);
          for (r = 0; r < R; r++) 
            {
              X1[r] = X[r][0];
            }

          // Store table title
          string title("x \\ f(x)");

          plotXYZ(X1, Y, Z, title.c_str(), filename);
        }
      else
        {
          // Tabulate (multiple) function(s) f(_i) : |R^k -> |R, k >= 2

          // Store table title
          string title;
          if (D._k == 2)
            {
              title = string("(x_1, y) \\ f(x, y)");
            }
          else
            {
              title = string("(x_1, ..., x_") + D._k + ") \\ f(X)";
            }

          plotXYZ(X, Y, Z, title.c_str(), filename);
        }          
    }
}

//----------------------------------------------------------------------------------------------

// Plot integrand of one or two variables to file after domain transformation and periodization
void CIntegrand::plot(const vector<CIntegrand*>& f, const CDomain& D, const CPeriodization& p, 
                      const string& filename, int res, int margin, bool excel)
{
  int M = f.size();
  int m;

  // Instantiate transformed integrands
  vector<CIntegrand*> ft;
  for (m = 0; m < M; m++) ft.push_back(new CTransformedintegrand(*(f[m]), D, &p));

  // Plot over the cube [-0.5, 1.5)^k to exhibit effects of periodization
  CCube cube(D._k, -0.5, 1.5);

  // Call static plot function on transformed integrands
  plot(ft, cube, filename, res, margin, excel);

  // Discard transformed integrands
  for (m = 0; m < M; m++) delete ft[0];
}

//----------------------------------------------------------------------------------------------

// Plot this integrand
void CIntegrand::plot(const CDomain& D, const string& filename, int res, int margin, bool excel)
{
  // Generate filename if required
  string fn = (filename.size() == 0) ? INFOPATH + acronym() + "-plot.txt" : filename;

  vector<CIntegrand*> f(1, this);

  plot(f, D, fn, res, margin, excel);
}

void CIntegrand::plot(const CDomain& D, const CPeriodization& p, const string& filename, int res, int margin, bool excel)
{
  // Generate filename if required
  string fn = (filename.size() == 0) ? INFOPATH + acronym() + "-plot.txt" : filename;

  vector<CIntegrand*> f(1, this);

  plot(f, D, p, fn, res, margin, excel);
}

//----------------------------------------------------------------------------------------------

double CIntegrand::evaluate(STLDoubleVector& x, const CDomain& D, const CPeriodization *p) const
{
  double productphi1x = 1.0;
  double Jphix;
  double fpsiphix;

  if (p)
    {
      modto01k(x);
      productphi1x = p->product(x);
      p->phi(x);
    }

  Jphix = D.J(x);

  D.psi(x);
  fpsiphix = evaluate(x);

  return fpsiphix * Jphix * productphi1x;
}

//----------------------------------------------------------------------------------------------

const char* CIntegrand::acronym(void) const 
{ 
  static string _S_other = "f-other"; 
  return (_S.size() == 0) ? _S_other.c_str() : _S.c_str(); 
}

//----------------------------------------------------------------------------------------------

double CIntegrand::I(int k) const
{
	string err = "In CIntegrand::I: analytic integral over unit cube [0, 1)^";
	err = err + k;
	err = err + " not known. ";
  ERROR(err)
}

//----------------------------------------------------------------------------------------------

CTransformedintegrand::CTransformedintegrand(const CIntegrand& f, const CDomain& D, 
                                             const CPeriodization *p) : _f(f), _D(D), _p(p) 
{
  _S = string(_f.acronym()) + string("-t");  
}

//----------------------------------------------------------------------------------------------

