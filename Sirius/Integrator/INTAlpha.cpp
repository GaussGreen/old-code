

#include	"stdafx.h"

#include "MCInt\INTAlpha.h"
#include "MCInt\INTUtilities.h"


//----------------------------------------------------------------------------------------------

// Update buffered alpha (or insert if no alpha buffered for dimension k)
STLDoubleVector& CAlpha::update(const STLDoubleVector& alpha, int k)
{
  if (k == 0) 
    {
      k = alpha.size();
    }
  else 
    {
      if (alpha.size() != k) ERROR("In CAlpha::update: vector alpha has wrong size. ")
    }

  _alphaDB.erase(k);

  EQSP_CC_STL_PAIR( AlphaDB::iterator, bool ) result = _alphaDB.insert(AlphaDB::value_type(k, alpha));

  return result.first->second;
};

//----------------------------------------------------------------------------------------------

// Retrieve buffered alpha for dimension k (compute if not available)
STLDoubleVector& CAlpha::get(STLDoubleVector& alpha, int k)
{
  AlphaDB::iterator result = _alphaDB.find(k);

  if (result == _alphaDB.end()) 
    {
      alpha = compute(k);
    }
  else
    {
      alpha = result->second;
    }

  return alpha;
};

//----------------------------------------------------------------------------------------------

// Plot distribution of sample points {m*alpha} in [0,1]^k (for k in {1, 2, 3})
void CAlpha::plotdistribution(int k, int N, const string& filename)
{
  // Number of rows = sampled points {m*alpha}
  int& R = N;
  int r;
  
  // Number of columns = dimension
  int& C = k;
  int c;

  // Row and column labels and table entries
  STLIntegerVector X(R, 0);
  vector< const char*> Y(C);
  STLDoubleVectorVector Z(R, STLDoubleVector(C, 0.0));

  // Fill in table

  // Retrieve alpha
  STLDoubleVector alpha; 
  get(alpha, k);

  // Current sample point
  STLDoubleVector x(k, 0.0);

  for (r = 0; r < R; r++)
    {
      x += alpha;
      modto01k(x);

      X[r] = r + 1;
      Z[r] = x;
    }

  // Store column labels 
  vector<string> x_i(C);
  for (c = 0; c < C; c++) 
    {
      x_i[c] = string("{n*alpha}_") + (c+1);
      Y[c] = x_i[c].c_str();
    }

  // Store table title
  string title = string("n \\ {n*alpha}_j for ") + acronym();

  // Generate filename if required
  string fn = (filename.size() == 0) ? INFOPATH + acronym() + "-plot.txt" : filename;

  plotXYZ(X, Y, Z, title.c_str(), fn);
}

//----------------------------------------------------------------------------------------------

// Report alpha
void CAlpha::report(int k)
{
  if (k > 0)
    {
      STLDoubleVector alpha;
	  std::cout << " " << get(alpha, k) << " in " << k << " dimensions";
    }
}

//----------------------------------------------------------------------------------------------

// Read in userdefined alpha for new dimension
STLDoubleVector& CUseralpha::compute(int k)
{ 
  std::cout << "Please enter an alpha:" << std::endl;
  STLDoubleVector alpha(k);
  for (int i = 0; i < k; i++)
    {
      std::cout << "Component No. " << i+1 << ": ";
      std::cin >> alpha[i];
      std::cout << std::endl;
    }

  return update(alpha, k);
};

//----------------------------------------------------------------------------------------------

CHaselgrovealpha::CHaselgrovealpha(int s) : _s(s) 
{ 
  if ((s != 2) && (s != 4)) 
    {
      ERROR("In CHaselgrovealpha::CHaselgrovealpha: s must be in {2, 4}. ") 
    }
  
  _S = string("a-has(") + _s + ")"; 
}

//----------------------------------------------------------------------------------------------

// Initialise and buffer alpha for new dimension
STLDoubleVector& CHaselgrovealpha::compute(int k)
{ 
  STLDoubleVector alpha;
  //CreateVectorClass<double> helperClass;

  switch(_s)
    {
    case 2:
      switch(k)
	      {
	      case 1: 
            CreateVectorClass<double>::createvector(alpha, 1, 0.73258893); 
	        break;
	      case 2:  
            CreateVectorClass<double>::createvector(alpha, 2, 0.62055505, 0.22610245);
	        break;
          case 3:
            CreateVectorClass<double>::createvector(alpha, 3, 0.96498949, 0.81091316, 0.46960090);
	        break;
          case 4:
            CreateVectorClass<double>::createvector(alpha, 4, 0.62366851, 0.04150108, 0.48574769, 0.27210703);
	        break;
          case 5:
            CreateVectorClass<double>::createvector(alpha, 5, 0.95734608, 0.86730270, 0.09724025, 0.31301950, 0.48476582);
	        break;
          case 6:
            CreateVectorClass<double>::createvector(alpha, 6, 0.43657951, 0.59185199, 0.05024400, 0.84373919, 0.38104000, 0.75808683);
	        break;
          case 7:
            CreateVectorClass<double>::createvector(alpha, 7, 0.80638723, 0.22584927, 0.72510075, 0.51310685, 0.11080509, 0.60161858, 0.92715171);
	        break;
          case 8:
            CreateVectorClass<double>::createvector(alpha, 8, 0.73750248, 0.08314415, 0.84753682, 0.88989711, 0.80254484, 0.27951501, 0.67340402, 0.53040927);
	        break;
          default:
            ERROR("In CHaselgrovealpha::compute: Haselgrove alphas only available for dimensions 1 <= k <= 8. ")
              }       
      break;
      
    case 4:
      switch(k)
        {
        case 1: 
          CreateVectorClass<double>::createvector(alpha, 1, 0.83969144);
          break;
        case 2: 
          CreateVectorClass<double>::createvector(alpha, 2, 0.59734470, 0.92828094);
          break;
        case 3: 
          CreateVectorClass<double>::createvector(alpha, 3, 0.74235492, 0.57387033, 0.32279917);
          break;
        case 4: 
          CreateVectorClass<double>::createvector(alpha, 4, 0.17665781, 0.71327190, 0.98875216, 0.60299793);
          break;
        case 5: 
          CreateVectorClass<double>::createvector(alpha, 5, 0.44810200, 0.53589831, 0.56039410, 0.83630131, 0.22148205);
          break;
        case 6: 
          CreateVectorClass<double>::createvector(alpha, 6, 0.10613747, 0.40278232, 0.88772556, 0.43554826, 0.17219381, 0.63794472);
          break;
        case 7: 
          CreateVectorClass<double>::createvector(alpha, 7, 0.58505729, 0.50196855, 0.77797734, 0.60504620, 0.62193588, 0.84244165, 0.64543976);
          break;
        case 8: 
          CreateVectorClass<double>::createvector(alpha, 8, 0.23975940, 0.01544979, 0.57794809, 0.81182909, 0.78068912, 0.62319488, 0.70710061, 0.60389317);
          break;
        default:
          ERROR("In CHaselgrovealpha::compute: Haselgrove alphas only available for dimensions 1 <= k <= 8. ")
	      }       
        break;

      default:
        ERROR("In CHaselgrovealpha::compute: Haselgrove alphas only available for s in {2, 4}. ")
    }

  return update(alpha, k);
}

//----------------------------------------------------------------------------------------------

// Compute and buffer cyclotomic alpha for new dimension
STLDoubleVector& CCyclotomicalpha::compute(int k)
{
  if (!isprime(2*k + 3))
    {
      std::cout << "In CCyclotomicalpha::compute: " << std::endl;
      std::cout << "Use of cyclotomic alpha is not recommended for k = " << k << std::endl;
      std::cout << "as 2*k + 3 = " << 2*k + 3 << " is not prime. Proceed anyway? (1/0)" << std::endl;
      
      bool proceedanyway;
      std::cin >> proceedanyway;

      std::cout << std::endl << std::endl;        
      if (!proceedanyway) ERROR("In CCyclotomicalpha::compute: ")
    }

  STLDoubleVector alpha(k);
  for (int i = 0; i < k; i++)
    {
      alpha[i] = 2.0*cos((2.0*PI*(double)(i + 1))/(double)(2*k + 3));
    }

  return update(alpha, k);
};

//----------------------------------------------------------------------------------------------

// Compute and buffer Baker alpha for new dimension
STLDoubleVector& CBakeralpha::compute(int k)
{
  STLDoubleVector alpha(k);
  for (int i = 0; i < k; i++)
    {
      alpha[i] = exp(_m*_r.get(i+1));
    }
  
  return update(alpha, k);
};

//----------------------------------------------------------------------------------------------

// Compute and buffer Niederreiter alpha for new dimension
STLDoubleVector& CNiederreiteralpha::compute(int k) 
{
  double log2bykplus1 = log(2.0)/(double)(k+1);

  STLDoubleVector alpha(k);
  for (int i = 0; i < k; i++)
    {
      alpha[i] = exp(log2bykplus1*(double)(i+1));
    }

  return update(alpha, k);
};

//----------------------------------------------------------------------------------------------

// Compute and buffer prime root alpha for new dimension
STLDoubleVector& CPrimerootalpha::compute(int k) 
{
  STLDoubleVector alpha(k);
  for (int i = 0; i < k; i++)
    {
      alpha[i] = sqrt(_p.get(i+1));
    }

  return update(alpha, k);
};

//----------------------------------------------------------------------------------------------

