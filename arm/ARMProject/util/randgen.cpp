/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE         : randgen.cpp                                                 */
/*                                                                            */
/* DESCRIPTION  : Random generators for Monte carlo Kernels                   */
/*                                                                            */
/*                Sep. 10th 1997                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/



#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <memory.h>

#include "armglob.h"
#include "randgen.h"

int* G_PREMIER = NULL;
int PREM_PREV_SIZE = 0;
double** G_HAMMERSLEY = NULL;
long G_NB_TIRAGES = 0;
int G_NB_DIM = 0;
static double MINI_RAND = 0.0000000001;

int first_random = 1; /* 1 : oui, 0 : non */
int rand_du_C = 0;    /* 1 : oui, 0 : non */

double* xn_rand = NULL;
int xn_rand_degre;

double* a_rand;
int a_rand_degre;

double* c_rand;
int c_rand_degre;

double m_rand;
int m_rand_degre;

double baseDefault;




void FreeHammersley(void)
{
    long i;


    for ( i = 0; i < G_NB_TIRAGES; i++)
    {
        free(G_HAMMERSLEY[i]);
    }

    if (G_HAMMERSLEY)
       free(G_HAMMERSLEY);
}


void reductionR(double* x /*entree*/, int degrex/*entree*/, double base/*entree*/, 
                double** y/*sortie*/, int* degrey/*sortie*/)
                /*simplifie x*/
{
    double* yt;
    double* ytbis;
    int degreyt;
    int degreytbis;
    double retenueplus;
    double quotient;
    int i;

    yt = new double[2*degrex+1];

    quotient = 0.0;

    for(i = 0; i <= degrex; i++)
    {
        retenueplus = x[i]+quotient;
        quotient    = floor(retenueplus/base);
        yt[i]       = retenueplus - quotient * base;

        if ( yt[i] > 0 )
        {
           degreyt=i;
        }
    }

    while ( quotient > 0 )
    {
        retenueplus = quotient;
        quotient    = floor(retenueplus/base);
        yt[i]       = retenueplus - quotient * base;

        if ( yt[i] > 0 )
        {
           degreyt=i;
        }

        i++;
    }

    degreytbis = degreyt;
    ytbis = new double[degreytbis+1];

    for (i = 0; i <= degreyt; i++)
    {
        ytbis[i]=yt[i];
    }

    *y = ytbis;
    *degrey = degreytbis;

    delete [] yt;
}


void additionR(double* x1/*entree*/, int degrex1/*entree*/,
               double* x2/*entree*/, int degrex2/*entree*/,
               double base/*entree*/, 
                double** y/*sortie*/, int* degrey/*sortie*/)
                /*retourne x1 + x2 */
{
    double* yt;
    double* ytbis;
    int degreyt;
    int degreytbis;
    int i;
    int del1 = 0;
    int del2 = 0;
    int degremin;



    if ( degrex1 > degrex2 )
    {
       del1=1;
    }
    else
    {
       del2=1;
    }

    degreyt = del1 * degrex1 + del2 * degrex2;
    degremin = degrex1 + degrex2 -degreyt;

    yt = new double[degreyt+1];

    for (i=0;i <= degremin; i++)
    {
        yt[i] = x1[i]+x2[i];
    }

    if ( del1 == 1 )
    {
       for (i=degremin+1;i<=degreyt;i++)
       {
           yt[i] = x1[i];
       }
    }
    else
    {
       for (i=degremin+1;i <= degreyt; i++)
       {
           yt[i] = x2[i];
       }

       if ( degrex1 == degrex2 )
       {
          yt[degreyt]+=x1[degreyt];
       }
    }

    reductionR(yt,degreyt,base, &ytbis,&degreytbis);

    *y = ytbis;
    *degrey = degreytbis;

    delete [] yt;
}



void multiplicationR(double* x1/*entree*/, int degrex1/*entree*/,
               double* x2/*entree*/, int degrex2/*entree*/,
               double base/*entree*/, 
                double** y/*sortie*/, int* degrey/*sortie*/)
                /*retourne x1 * x2 */
{
    double* yt;
    double* ytbis;
    int degreyt;
    int degreytbis;
    int i,j;
    double sommep;

    degreyt = degrex1+degrex2;

    yt = new double[degreyt+1];

    for (i = 0; i <= degreyt;i++)
    {
        sommep = 0.0;

        for(j = MAX(0, i-degrex1) ; j<=MIN(i, degrex2); j++)
        {    
            sommep += x1[(i-j)]*x2[(j)];
        }

        yt[i] = sommep;
    }

    reductionR(yt,degreyt,base, &ytbis,&degreytbis);

    *y = ytbis;
    *degrey = degreytbis;

    delete [] yt;
}



void multsimpR(double* x1/*entree*/, int degrex1/*entree*/,
               double x2/*entree*/, int degrex2/*entree*/,
               double base/*entree*/, 
                double** y/*sortie*/, int* degrey/*sortie*/)
                /*
                retourne x1 * x2*base^degrex2, x2<base.
                */
{
    double* yt;
    double* ytbis;
    int degreyt;
    int degreytbis;
    int i;

    degreyt = degrex1;
    yt = new double [degreyt+1];
    
    for (i=0;i<=degrex1;i++)
    {
        yt[i] = x2 * x1[i];
    }

    reductionR(yt, degreyt, base, &ytbis, &degreytbis);

    delete [] yt;
    
    degreyt = degrex1 + degrex2;

    yt = new double [degreyt+1];

    for(i = 0; i < degrex2; i++)
    {
        yt[i]=0.0;
    }

    for(i=degrex2;i<degreyt;i++)
    {
        yt[i]=ytbis[i-degrex2];
    }

    *y = yt;
    *degrey = degreyt;

    delete [] ytbis;    
}



void quotientR(double* x1/*entree*/, int degrex1/*entree*/,
               double x2/*entree*/, int degrex2/*entree*/,
               double base/*entree*/, 
                double** y/*sortie*/, int* degrey/*sortie*/)
            /*
              retourne le quotient de la division de x1 par x2*base^degrex2, x2<base.
            */
{
    double* yt;
    double* ytbis;
    int degreyt;
    int degreytbis;
    int i;
    double reste;
    double inter;

    if(degrex1<degrex2)
    {
        degreytbis = 0;
        ytbis = new double[degreytbis+1];
        ytbis[0] = 0.0;
        *y = ytbis;
        *degrey = degreytbis;
    }
    else
    {
        degreytbis = degrex1-degrex2;
        ytbis = new double [degreytbis+1];

        reste = 0.0;

        for (i=degreytbis; i>=0; i--)
        {
            inter = x1[degrex2+i]+reste*base;
            ytbis[i] = floor((inter)/x2); 
            reste = inter - ytbis[i]*x2;
        }
        reductionR(ytbis,degreytbis,base, &yt,&degreyt);        
        *y = yt;
        *degrey = degreyt;
        delete [] ytbis;
    }
}



void resteR(double* x1/*entree*/, int degrex1/*entree*/,
               double x2/*entree*/, int degrex2/*entree*/,
               double base/*entree*/, 
               double** y/*sortie*/, int* degrey/*sortie*/)
            /*
               retourne le reste de la division de x1 par x2*base^degrex2, x2<base.
            */
{
    double* yt;
    double* ytbis;
    int degreyt;
    int degreytbis;
    int i;
    double reste;
    double inter;


    if ( degrex1 < degrex2 )
    {
       degreytbis = degrex1;
       ytbis = new double[degreytbis+1];

       for (i=0; i<=degreytbis; i++)
       {
           ytbis[i] = x1[i];
       }
       *y = ytbis;
       *degrey = degreytbis;
    }
    else
    {
       reste = 0.0;

       for(i=degrex1-degrex2; i>=0; i--)
       {
          inter = x1[degrex2+i]+reste*base;
          reste = inter - (floor((inter)/x2))*x2;
       }
    
       degreytbis = degrex2;
       ytbis = new double[degreytbis+1];
       for(i=0; i<degrex2; i++)
       {
          ytbis[i] = x1[i];
       }

       ytbis[degreytbis] = reste;        

       reductionR(ytbis,degreytbis,base,&yt,&degreyt);
        
       *y = yt;
       *degrey = degreyt;
       delete [] ytbis;
    }
}


double engrosR(double* x, int degrex, double base)
{
    double res;
    int i;

    res = 0.0;

    for(i=degrex; i>=0;i--)
    {
        res = res * base + x[i]; 
    }

    return(res);
}



void init_Rand(int deterministe, int meth)
/*
    si deterministe==1, toujours la meme serie, sinon tirage aleatoire

    meth 2 : rand du C, meth 0 ou 1 : congruence de grands nombres
                                    (plus lent mais meilleur)
*/
{
    double initval;

    baseDefault = 16777216.0;

    if (deterministe == 1)
    {
        initval = 1.0;
        srand(1);
    }
    else
    {
        srand( (unsigned)time( NULL ) );
        initval = (double) rand();
    }

    if (meth == 2)
    {
       rand_du_C = 1;
    }
    else
    {
       if (first_random != 1)
       {
          delete [] xn_rand;
       }

       first_random = 0;
   
       reductionR(&initval,0, baseDefault, &xn_rand, &xn_rand_degre);

       if (meth==0)
       {
          /*13^13*/
          a_rand_degre = 2;
          a_rand = new double [a_rand_degre+1];
          a_rand[0] = 2344445.0;
          a_rand[1] = 1275547.0;
          a_rand[2] = 1.0;
    
          c_rand_degre = 0;
          c_rand = new double [c_rand_degre+1];
          c_rand[0] = 0.0;
    
          /*2^59*/
          m_rand_degre = 2;
          m_rand = 2048.0;
       }
       else
       {
          /*5^17*/
          a_rand_degre = 1;
          a_rand = new double [a_rand_degre+1];
          a_rand[0] = 12332741.0;
          a_rand[1] = 45474.0;
   
          c_rand_degre = 0;
          c_rand = new double [c_rand_degre+1];
          c_rand[0] = 1.0;
    
          /*2^48*/
          m_rand_degre = 2;
          m_rand = 1.0;
       }
    }
}



double myRipleyRandom(void)
{
    double res;


    if ( rand_du_C == 1 )
    {
        res = ((double) rand())/RAND_MAX;
    }
    else
    {
        double* xloc;
        int xlocdegre;
        double* yloc;
        int ylocdegre;
        int i;

        multiplicationR(xn_rand,xn_rand_degre,a_rand,a_rand_degre,
                        baseDefault,&xloc,&xlocdegre);

        delete [] xn_rand;

        additionR(xloc,xlocdegre,c_rand,c_rand_degre,baseDefault,&yloc,&ylocdegre);

        resteR(yloc,ylocdegre,m_rand,m_rand_degre, baseDefault,
               &xn_rand,&xn_rand_degre);

        delete[] xloc;
        delete[] yloc;

        res = engrosR(xn_rand,xn_rand_degre,baseDefault);

        for (i=1; i<=m_rand_degre; i++)
        {
            res = res/baseDefault;
        }
        res = res/m_rand;
    }

    res = MINI_RAND + (1-2*MINI_RAND)*res;

    return(res);
}



/******************************************************************************/


/*methode de quasi Monte-Carlo*/


/*----------------------------------------------------------------------------*/
/* Utilities functions                                                        */
/*----------------------------------------------------------------------------*/

long IntValue(double v)/*utile ?*/
{
    char buf[50];
    long  intVal;


    sprintf(buf, "%lf", floor(v));
    intVal = atol(buf);

    return(intVal);
}

 
long DivInt(long p1, int p2)/*utile ?*/
{
    char buf[50];


// FIXMEFRED: mig.vc8 (22/05/2007 15:52:55): long/int != double
	sprintf(buf, "%lf", p1/p2);
    return(atol(buf));
}



/*----------------------------------------------------------------------------*/
/* Generates "nombres premiers" in Global Variable G_PREMIER                  */
/*----------------------------------------------------------------------------*/

void GetPrimeNumber(int nb)
{
    int i,n,test,j,quotient,div;


    if ( nb > PREM_PREV_SIZE )
    {

       if ( PREM_PREV_SIZE == 0 )
       {
          G_PREMIER = (int *) malloc(nb*sizeof(int));

          G_PREMIER[0] = 2;
       }
       else
       {
          G_PREMIER = (int *) realloc(G_PREMIER, nb-PREM_PREV_SIZE);
       }

       if ( (PREM_PREV_SIZE <= 1) && (nb > 1) )
       {
          G_PREMIER[1] = 3;
       }

       for ( i = MAX(PREM_PREV_SIZE, 2); i < nb; i++)
       {
           n = G_PREMIER[i-1];

           test=0;

           while ( test == 0 )
           {
               n = n+2;
               j=0;
               div=0;
               quotient = n;

               while ((( G_PREMIER[MAX(j-1, 0)] < quotient ) 
                         || ( div == 0 )) && (j<i))
               {
                   quotient = DivInt(n, G_PREMIER[j]);

                   if ( n == quotient*G_PREMIER[j] )
                   {
                      div = 1;
                   }

                  j++;
               }

               if ( div == 0 )
               {
                   test = 1;
               }
           }

           G_PREMIER[i]=n;
       }

       PREM_PREV_SIZE = nb; 
    }
}



double PhiDiscrep(long nb, int p)
{
    int k,i;
    long a, quotient, puisp;
    int* decomp;
    double result;



    k = (int) IntValue((double) log((double) nb) / (double) log((double) p)) + 1;

    decomp = new int[k];

    a = nb;

    for ( i = 0 ; i < k ; i++ )
    {
        quotient = DivInt(a,p);

        decomp[i] =(int) (a-quotient*p);

        a = quotient;
    }

    result = 0.0;
    puisp = (long) p; 

    for ( i = 0 ; i < k ; i++ )
    {
        result = result + (double) decomp[i] / (double) puisp;
        puisp = puisp * p;
    }

    delete [] decomp;

    return(result);
}



double Hammersley(long numtirage, int dimInt)
{
    return(PhiDiscrep(numtirage+1, G_PREMIER[dimInt]));
}



void GenerateHammersley(int maxDim, long maxTirages)
{
    long i;
    int  j;



    if ( maxTirages > G_NB_TIRAGES )
    {
       if ( G_NB_TIRAGES == 0 )
       {
          G_HAMMERSLEY = (double **) malloc(sizeof(double*)*maxTirages);

          for ( i = 0; i < maxTirages; i++)
          {
              G_HAMMERSLEY[i] = (double *) malloc(sizeof(double)*G_NB_DIM); 
          }
       }
       else
       {
          G_HAMMERSLEY = (double **) realloc(G_HAMMERSLEY, 
                               sizeof(double*)*maxTirages);

          for ( i = G_NB_TIRAGES; i < maxTirages; i++)
          {
              G_HAMMERSLEY[i] = (double *) malloc(sizeof(double)*G_NB_DIM); 
          }
       }

       for ( i = G_NB_TIRAGES; i < maxTirages; i++)
       {
           for ( j = 0; j < G_NB_DIM; j++)
           {
               G_HAMMERSLEY[i][j] = Hammersley(i, j); 
           } 
       }
    } 

    if ( maxDim > G_NB_DIM )
    {
       for ( i = 0 ; i < maxTirages ; i++)
       {
           G_HAMMERSLEY[i] = 
                (double *) realloc(G_HAMMERSLEY[i], sizeof(double)*maxDim);
        
           for ( j = G_NB_DIM ; j < maxDim ; j++ )
           {
               G_HAMMERSLEY[i][j] = Hammersley(i,j);
           }
       }
    }

    G_NB_TIRAGES = MAX(G_NB_TIRAGES, maxTirages);

    G_NB_DIM = MAX(G_NB_DIM, maxDim);
}



double DinfiStar(int dimInteg, long nbTir)
{
    int i;
    double product;


    if ( dimInteg == 1 )
    {
       return(1.0/nbTir);
    }

    product = 1.0;

    for ( i = 0 ; i < dimInteg-1 ; i++)
    {
        product = product*G_PREMIER[i] 
                   *log((double)(nbTir*G_PREMIER[i]))/log(static_cast<double>(G_PREMIER[i]));
    }

    product = product / (double) nbTir;

    return(product);    
}



double Dinfi(int dimInteg, long nbTir)
{
    int i;
    long puiss;
    
    puiss = 1;

    for (i = 1 ; i <= dimInteg; i++)
    {
         puiss = puiss*2;
    }

    return((double) puiss * DinfiStar(dimInteg, nbTir));
}



double Jinfi(int dimInteg, long nbTir)
{
    double res;


    if ( dimInteg == 1 )
    {
       return(2.0/nbTir);
    }

    res = pow(Dinfi(dimInteg, nbTir), 1.0/dimInteg);

    res = (4.0 * dimInteg * sqrt((double)dimInteg)+1.0)*res;

    return(res); 
}



void HammersleyAndPremiers(int dimInteg, long nbTir)
{
    GetPrimeNumber(dimInteg);

    GenerateHammersley(dimInteg, nbTir);
}


/*----------------------------------------------------------------------------*/
/*---- End of File ----*/
