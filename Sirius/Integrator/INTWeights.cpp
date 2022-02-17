#include	"stdafx.h"



#include "MCInt\INTWeights.h"




// Plot weight distributions (m_i, C*A(m_i)) vs i
void CWeights::plotdistribution(const vector<CWeights*>& w, long N, const string& filename, int res)
{
  // Number of rows = evaluations of A_m
  int R = res;
  int r;

  // Number of colums = 2 x number of weight distributions to plot
  int C = w.size();
  int Ceff = 2*C;
  int c;

  // Row and column labels and table entries
  STLIntegerVector X(R, 0);
  vector<const char*> Y(Ceff, (const char*)NULL);
  STLDoubleVectorVector Z(R, STLDoubleVector(Ceff, 0.0));

  // Validate input and obtain Neff, Nrange and C for all weights
  vector<long> Neffw(C, 0);
  vector< EQSP_CC_STL_PAIR(long, long) > Nrangew(C, make_pair((long)0, (long)0));

  STLDoubleVector Cw(C, 0.0);

  for (c = 0; c < C; c++) 
    {
      if (!w[c]) ERROR("In CWeights::plotdistribution: weights invalid. ")

      Neffw[c] = w[c]->Neff(N);
      w[c]->Nrange(Nrangew[c], Neffw[c]);
      Cw[c] = w[c]->C(Neffw[c]);
    }

  // Fill in table 

  long m;
  double mprop;
  double CAm;
  
  for (r = 0; r < R; r++)
    {
      X[r] = r+1;
      mprop = (double)r/(double)(res - 1);

      for (c = 0; c < C; c++)
        {
          m = (long)((double)(Nrangew[c].second - Nrangew[c].first)*mprop) + Nrangew[c].first;
          CAm = Cw[c] * w[c]->A(Neffw[c], m);

          Z[r][2*c] = m;
          Z[r][2*c+1] = CAm;
        }
    }

  // Store column labels
  string empty(" ");
  for (c = 0; c < C; c++) 
    {
      Y[2*c] = empty.c_str();
      Y[2*c+1] = w[c]->acronym();
    }

  // Store table title
  string title("i \\ (m_i, w.C*w.A(m_i))");

  plotXYZ(X, Y, Z, title.c_str(), filename);
}

//----------------------------------------------------------------------------------------------

// Plot this weight distribution
void CWeights::plotdistribution(long N, const string& filename, int res)
{
  // Generate filename if required
  string fn = (filename.size() == 0) ? INFOPATH + acronym() + "-plot.txt" : filename;

  vector<CWeights*> w(1, this);

  plotdistribution(w, N, fn, res);
}

//----------------------------------------------------------------------------------------------

CHaselgroveweights::CHaselgroveweights(int r) : CWeights(r) 
{ 
  if ((r <= 0) || (r >= 5)) 
    {
      ERROR("In CHaselgroveweights::CHaselgroveweights: r must be in {1, 2, 3, 4}. ") 
    }

  _S = string("w-has(") + _r + ")";
}

//----------------------------------------------------------------------------------------------

// Overall constant factor
double CHaselgroveweights::C(long N) const
{
  switch (_r)
    {
    case 1:
      return 1.0/(double)(2*N + 1);
    case 2:
      return 1.0/pow((double)(N + 1), 2);
    case 3:
      return 1.0/(pow((double)(N + 1), 2)*(2.0*(double)N + 3.0));
    case 4:
      return 1.0/pow((double)(N + 1), 4);
    default:
      ERROR("In CHaselgroveweights::C: Haselgrove weights only available for r in {1, 2, 3, 4}. ")
    }

  return 0.0;
}

//----------------------------------------------------------------------------------------------

// Weighting coefficient for f(m*alpha)
double CHaselgroveweights::A(long N, long m) const
{
  switch (_r)
    {
    case 1:
    case 2:
      return Ar(N, m);
    case 3:
      return Ar(2*N+1, m) - (fabs(m) <= N ? 2.0*Ar(N, m) : 0);
    case 4:
      return Ar(2*N, m) - (fabs(m) <= N-1 ? 4.0*Ar(N-1, m) : 0);
    default:
      ERROR("In CHaselgroveweights::C: Haselgrove weights only available for r in {1, 2, 3, 4}. ")
    }

  return 0.0;
}

//----------------------------------------------------------------------------------------------

// Rescaled number of sampling points (to ensure f is evaluated N times independently of weights)
long CHaselgroveweights::Neff(long N) const
{
  switch (_r)
    {
    case 1:
      return (double)(N - 1)/2.0;
    case 2:
      return (double)(N - 1)/2.0;
    case 3:
      return (double)(N - 3)/4.0;
    case 4:
      return (double)(N - 1)/4.0;
    default:
      ERROR("In CHaselgroveweights::C: Haselgrove weights only available for r in {1, 2, 3, 4}. ")
    }

  return 0.0;
}

//----------------------------------------------------------------------------------------------

// Range {N_min, ..., N_max} of m, where Nrange is the pair (N_min, N_max)
EQSP_CC_STL_PAIR(long, long)& CHaselgroveweights::Nrange(EQSP_CC_STL_PAIR(long, long)& Nrange, long N) const
{
  switch (_r)
    {
    case 1:
    case 2:
      Nrange = make_pair(-N, N);
      break;
    case 3:
      Nrange = make_pair(-(2*N + 1), 2*N + 1);
      break;
    case 4:
      Nrange = make_pair(-2*N, 2*N);
      break;
    default:
      ERROR("In CHaselgroveweights::C: Haselgrove weights only available for r in {1, 2, 3, 4}. ")
    }

  return Nrange;
}

//----------------------------------------------------------------------------------------------

CNiederreiterweights::CNiederreiterweights(int r) : CWeights(r) 
{ 
  if (r <= 0) ERROR("In CHaselgroveweights::CHaselgroveweights: r must be positive. ") 

  _S = string("w-nie(") + _r + ")";
}

//----------------------------------------------------------------------------------------------

// Weighting coefficient for f(m*alpha)
double CNiederreiterweights::A(long N, long m) const
{
  int i_max = m/(N+1);
  double sum = 0.0;
  double sign = 1.0;

  for (int i = 0; i <= i_max; i++)
    {
      sum += sign * binomial(_r, i) * binomial(m - (N+1)*(long)i + (long)(_r-1), (long)(_r-1));
      sign *= -1.0;
    }
  
  return sum;
}

//----------------------------------------------------------------------------------------------

// Range {N_min, ..., N_max} of m, where Nrange is the pair (N_min, N_max)
EQSP_CC_STL_PAIR(long, long)& CNiederreiterweights::Nrange(EQSP_CC_STL_PAIR(long, long)& Nrange, long N) const
{ 
  Nrange = make_pair((long)0, (long)_r*N); 
  return Nrange; 
} 

//----------------------------------------------------------------------------------------------

CSugiharaMurotaweights::CSugiharaMurotaweights(int r) : CWeights(r)
{
  if (_r <= 0) ERROR("In CSugiharaMurotaweights::CSugiharaMurotaweights: r must be positive. ")
 
  // Note: (2*r-1)!/((r-1)!)^2 = r*binomial(2*r-1, r-1)
  _c = (double)r*binomial(2*r - 1, r - 1);

  _S = string("w-sm(") + _r + ")"; 
}

