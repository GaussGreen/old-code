#include	"stdafx.h"



#include "MCInt\INTEnumerator.h"
#include "INTIntegrator.h"

#include "MCInt\INTSystematic.h"


//----------------------------------------------------------------------------------------------

void systematic(const vector<CIntegrand*>& f,		const CDomain& D, 
				STLDoubleVector& I,							const STLLongVector& N, 
                const vector<CAlpha*>& a,			const vector<CWeights*>& w,  
				const vector<CPeriodization*>& p,	const STLIntegerVectorVector& cat, 
				const STLIntegerVectorVector& levels
				)
{
  // Number of integrands with domains to be included in figure of merit
  int Q = f.size();
  if (Q == 0) ERROR("In systematic: must provide at least one integrand. ")
  bool computeI = (I.size() == 0);
  if ((!computeI) && (I.size() != Q)) ERROR("In systematic: wrong number of analytic integrals. ")

  // Obtain analytic integrals
  int q;  
  for (q = 0; q < Q; q++) 
    {
      if (!f[q]) ERROR("In systematic: integrand invalid. ")
      if (computeI) I[q] = f[q]->I(D._k);
    }

  // Number of alphas, weights, and periodizations to be tested
  int i;
  STLIntegerVector n(3);
  n[0] = a.size();
  for (i = 0; i < n[0]; i++) if (!a[i]) ERROR("In systematic: alpha invalid. ")
  n[1] = w.size();
  for (i = 0; i < n[1]; i++) if (!w[i]) ERROR("In systematic: weights invalid. ")
  n[2] = p.size();
  for (i = 0; i < n[2]; i++) if (!p[i]) ERROR("In systematic: periodization invalid. ")

  // Number of categorizations requested
  STLIntegerVectorVector _cat(cat);
  int G = cat.size();
  int g;

  if (G == 0)
    {
      // Unless otherwise specified all categorizations will be produced
      CPermutator p(3);
      do 
      { 
        _cat.push_back(p._x); 
        G++;
      } 
      while(++p);
    }
  else
    {
      // Validate input
      for (g = 0; g < G; g++)
        {
          if (!CPermutator::ispermutation(3, cat[g])) 
            {
              ERROR("In systematic: cat[g] is not a permutation of {0,1,2}. ")
            }
        }
    }

  // Number of print levels requested per categorization
  STLIntegerVectorVector _levels(levels);
  int H;
  int h;

  if (levels.size() == 0)
    {
      // Unless otherwise specified all plot levels will be produced
      _levels.resize(G, STLIntegerVector(3, 0));
      for (g = 0; g < G; g++) 
        {
          for (h = 0; h < 3; h++)
            {
              _levels[g][h] = h + 1;
            }
        }
    }
  else
    {
      // Validate input
      if (levels.size() != G) ERROR("In systematic: Must request levels for all categorizations. ")
      for (g = 0; g < G; g++) 
        {        
          H = levels[g].size();
          for (h = 0; h < H; h++)
            {
              if ((levels[g][h] <= 0) || (levels[g][h] > 3)) 
                {
                  ERROR("In systematic: levels[g] must be a subset of {1, 2, 3}. ")
                }
            }
        }    
    }

  // Number of rows R = number of supplied N
  int R = N.size();
  int r;

  // Number of columns C = number of NT-parameter configurations
  int C = n[0]*n[1]*n[2];
  int c;

  // Row and column labels and table entries
  // Note: (C + 1)st column is for the results of Monte Carlo Integration
  STLLongVector X(N);
  STLDoubleVectorVector Z(R, STLDoubleVector(C + 1, 0.0));

  // Fill in tables

  // Enumerator for lattice of NT-parameter configurations
  CEnumerator e(3, n);

  // Current NT-parameter configuration
  CNTparameters NTpar;

  // MC- and NT-integrator
  CMCintegrator MCint;
  CNTintegrator NTint;

  // Numerical approximation of integral, relative error and RMS relative error
  double numI, relerrI;
  double RMS;

  for (r = 0; r < R; r++)
    {
      std::cout << "N = " << N[r] << ": " << std::endl;

      for (c = 0; c < C; c++)
        {              
          // Pick next NT-parameter configuration
          NTpar = CNTparameters(a[e._x[0]], w[e._x[1]], p[e._x[2]]);

          // Compute Relative Root Mean Squared Error (RMS) using NT-integration
          RMS = 0.0;
          for (q = 0; q < Q; q++)
            {
              numI = NTint.integrate(*(f[q]), D, N[r], NTpar);
              relerrI = fabs(I[q] - numI)/I[q];
              RMS += relerrI*relerrI;
            }
          Z[r][c] = sqrt(RMS/(double)Q);

          // Output NT result
          std::cout << Z[r][c] << " ";
          if (c % 5 == 4) std::cout << std::endl;
         
          // Increase enumerator
          ++e;

          std::cout << ".";
        }
      
      // Compute Root Mean Square (RMS) relative error using MC-integration
      RMS = 0.0;
      for (q = 0; q < Q; q++)
        {
          numI = MCint.integrate(*(f[q]), D, N[r]);
          relerrI = fabs(I[q] - numI)/I[q];
          RMS += relerrI*relerrI;
        }
      Z[r][C] = sqrt(RMS/(double)Q);

      // Output MC result
      std::cout << Z[r][C] << std::endl;

      // Reset enumerator
      e.reset();

      std::cout << std::endl;
    }

  // For each categorization
  for (g = 0; g < G; g++)
    {
      // Inverse of category permutation
      STLIntegerVector cat_inv(3);
      for (i = 0; i < 3; i++) cat_inv[_cat[g][i]] = i;

      // Anacronyms for categories
      vector<string> catanacronym(3);
      catanacronym[cat_inv[0]] = "a";
      catanacronym[cat_inv[1]] = "w";
      catanacronym[cat_inv[2]] = "p";

      // Size of categories
      STLIntegerVector n_cat;
      CreateVectorClass<int>::createvector(n_cat, 3, n[_cat[g][0]], n[_cat[g][1]], n[_cat[g][2]]);

      // For each print level per categorization
      H = _levels[g].size();
      for (h = 0; h < H; h++)
        {         
          // Number of tables to be produced
          int T;
          int t;

          switch(_levels[g][h])
            {
            case 1: 
                // One table overall
                T = 1; 
                break;
            case 2: 
                // One table per entry x in category 0
                T = n_cat[0]; 
                break;
            case 3:
                // One table per pair of entries (x, y) in category 0 x category 1
                T = n_cat[0] * n_cat[1]; 
                break;
            }
          
          // Number of columns per table
          int C_t = n[0]*n[1]*n[2]/T;

          // Row and column labels and table entries
          vector<const char*> Y_t(C_t + 1);
          STLDoubleVectorVector Z_t(R, STLDoubleVector(C_t + 1, 0.0));

          // Enumerator for lattice of NT-parameter configurations
          CEnumerator e_g(3, n_cat);

          for (t = 0; t < T; t++)
            {
              // Generate table title
              vector<string> titleNTpar(3);
              titleNTpar[cat_inv[0]] = a[e_g._x[cat_inv[0]]]->acronym();
              titleNTpar[cat_inv[1]] = w[e_g._x[cat_inv[1]]]->acronym();
              titleNTpar[cat_inv[2]] = p[e_g._x[cat_inv[2]]]->acronym();
              titleNTpar[2] = catanacronym[2];
              if (_levels[g][h] <= 2) titleNTpar[1] = catanacronym[1];
              if (_levels[g][h] == 1) titleNTpar[0] = catanacronym[0];
              string title = string("N vs RMS(") + Q + " x " + f[0]->acronym() + ", " + D.acronym() + ", N, ";
              title = title + "(" + titleNTpar[0] + ", " + titleNTpar[1] + ", " + titleNTpar[2] + "))";

              // Transcribe (required part of the) table of results and store column labels 
              vector<string> label(C_t + 1);
              vector<string> NTpar(3);
              for (c = 0; c < C_t; c++)
                { 
                  // Retrieve correct column index
                  int c_source = e_g._x[cat_inv[0]]*n[1]*n[2] + e_g._x[cat_inv[1]]*n[2] + e_g._x[cat_inv[2]];
                  
                  // Transcribe all entries down that column
                  for (r = 0; r < R; r++) Z_t[r][c] = Z[r][c_source];

                  // Generate column label
                  NTpar[cat_inv[0]] = a[e_g._x[cat_inv[0]]]->acronym();
                  NTpar[cat_inv[1]] = w[e_g._x[cat_inv[1]]]->acronym();
                  NTpar[cat_inv[2]] = p[e_g._x[cat_inv[2]]]->acronym();
                  if (_levels[g][h] >= 2) NTpar[0] = catanacronym[0];
                  if (_levels[g][h] == 3) NTpar[1] = catanacronym[1];
                  label[c] = string("(") + NTpar[0] + ", " + NTpar[1] + ", " + NTpar[2] + ")";
                  Y_t[c] = label[c].c_str();

                  // Increase enumerator
                  ++e_g;
                }

              // Transcribe MC column with label
              for (r = 0; r < R; r++) Z_t[r][C_t] = Z[r][C];
              label[C_t] = string("MC");
              Y_t[C_t] = label[C_t].c_str();

              // Generate filename
              string fn = INFOPATH + _levels[g][h] + " - " + title + ".txt";

              // Produce plot
              plotXYZ(X, Y_t, Z_t, title.c_str(), fn);  
            }                  

            // Reset enumerator
            e_g.reset();
        }
    }
}

//----------------------------------------------------------------------------------------------

