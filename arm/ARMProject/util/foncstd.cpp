/*
 * $Log: foncstd.cpp,v $
 * Revision 1.10  2003/10/30 20:37:12  ebenhamou
 * change for DotNet compilation
 *
 * Revision 1.9  2003/10/21 12:19:27  jpriaudel
 * modif du define de Precision en ARM_Precision
 *
 * Revision 1.8  2003/01/07 17:46:32  sgasquet
 * Correction bug Correlation a -1
 *
 * Revision 1.7  2003/01/07 16:37:13  mab
 * Formatting
 *
 * Revision 1.6  2002/11/25 16:23:52  mab
 * Formatting
 *
 */

/*---------------------------------------------------------------------------*/


#include "foncstd.h"
#include "math.h"
#include "gaussian.h"
#include "armglob.h"
#include "interpol.h"


#ifdef WIN32
//#include "gpclosedForms/spreadoption_lognormal_interface.h"
//#include "gpclosedForms/spreadoption_lognormal_formula.h"
//#include "gpclosedForms/bisabr_spreadoption.h"
//#include "gpclosedForms/spreadoption_bisabr_interface.h"
//#include "gpclosedForms/spreadoption_sabr_interface.h"
//#include "gpclosedForms/bisabr_interface.h"

//using namespace ARM;
#endif


/*---------------------------------------------------------------------------*/
/* Calcul du prix d'une option classique sur spread                          */
/*---------------------------------------------------------------------------*/

double fonct_std_d0(double S1, double S2, double vol1, double vol2,
                double d1, double d2, double rho,
                double r, double K, double t, double x)
{
    // coefficient associe a K;
    // calcul de f(Xsqrt(2))
    // pc vaut 1 pour un call et - 1 pour un put

    double y;

    double a1,
        a2,
        b1,
        b2,
        c1,
        c2;


    if (K > 0)
    {
        a2 = S2 * exp(t * (r - d2 - vol2 * vol2 * 0.5));
        b2 = vol2 * rho * sqrt(t);
        c2 = vol2 * sqrt(1 - rho * rho) * sqrt(t);
        a1 = S1 * exp(t * (r - d1 - vol1 * vol1 * 0.5));
        b1 = vol1 * sqrt(t);

        y = log(a1 / a2 * exp(b1 * x) + K / a2);
        y = y - b2 * x;
        y = y / c2;

        return (-y);
    }
    else
        //K <= 0
    {
        a1 = S1 * exp(t * (r - d1 - vol1 * vol1 * 0.5));
        b1 = vol1 * rho * sqrt(t);
        c1 = vol1 * sqrt(1 - rho * rho) * sqrt(t);
        a2 = S2 * exp(t * (r - d2 - vol2 * vol2 * 0.5));
        b2 = vol2 * sqrt(t);

        y = log(a2 / a1 * exp(b2 * x) - K / a1);
        y = y - b1 * x;
        y = y / c1;
        return (y);
    }
}



double fonct_std_d2(double S1, double S2, double vol1, double vol2,
                double d1, double d2, double rho, double r,
                double K, double t, double x)
{
    //coefficient associ ë ÇS2
    // introduction proba Q2
    double y;
    double a1,
        a2,
        b1,
        b2,
        c1,
        c2;

    if (K > 0)
    {
        a2 = S2 * exp(t * (r - d2 + vol2 * vol2 * 0.5));
        b2 = vol2 * rho * sqrt(t);
        c2 = vol2 * sqrt(1 - rho * rho) * sqrt(t);
        a1 = S1 * exp(t * (r - d1 - vol1 * vol1 * 0.5 + rho * vol1 * vol2));
        b1 = vol1 * sqrt(t);

        y = log(a1 / a2 * exp(b1 * x) + K / a2);
        y = y - b2 * x;
        y = y / c2;

        return (-y);
    }
    else
        //K <= 0
    {
        a1 = S1 * exp(t * (r - d1 - vol1 * vol1 * 0.5 + rho * vol1 * vol2));
        b1 = vol1 * rho * sqrt(t);
        c1 = vol1 * sqrt(1 - rho * rho) * sqrt(t);
        a2 = S2 * exp(t * (r - d2 + vol2 * vol2 * 0.5));
        b2 = vol2 * sqrt(t);

        y = log(a2 / a1 * exp(b2 * x) - K / a1);
        y = y - b1 * x;
        y = y / c1;

        return (y);
    }
}



double fonct_std_d1(double S1, double S2, double vol1, double vol2,
                double d1, double d2, double rho, double r,
                double K, double t, double x)
{
    //coefficient associe a S1
    // introduction proba Q1

    double y;
    double a1,
        a2,
        b1,
        b2,
        c1,
        c2;

    if (K > 0)
    {
        a2 = S2 * exp(t * (r - d2 - vol2 * vol2 * 0.5 + rho * vol1 * vol2));
        b2 = vol2 * rho * sqrt(t);
        c2 = vol2 * sqrt(1 - rho * rho) * sqrt(t);
        a1 = S1 * exp(t * (r - d1 + vol1 * vol1 * 0.5));
        b1 = vol1 * sqrt(t);

        y = log(a1 / a2 * exp(b1 * x) + K / a2);
        y = y - b2 * x;
        y = y / c2;

        return (-y);
    }
    else
        //K <= 0
    {
        a1 = S1 * exp(t * (r - d1 + vol1 * vol1 * 0.5));
        b1 = vol1 * rho * sqrt(t);
        c1 = vol1 * sqrt(1 - rho * rho) * sqrt(t);
        a2 = S2 * exp(t * (r - d2 - vol2 * vol2 * 0.5 + rho * vol1 * vol2));
        b2 = vol2 * sqrt(t);

        y = log(a2 / a1 * exp(b2 * x) - K / a1);
        y = y - b1 * x;
        y = y / c1;

        return (y);
    }
}



//************************************************
// Calcul du prix d 'une option std sur spread   *
//************************************************

double Proba0(double S1,
              double S2,
              double vol1,
              double vol2,
              double d1,
              double d2,
              double rho,
              double r,
              double K,
              double t)
{
    int n = 49;
    double N_x = 0.0;
    double integrK = 0.0;
    double coor;
    int i;


    // Calcul du coefficient associe a K

    for (i = 0; i < n; i++)
    {
        coor = fonct_std_d0(S1, S2, vol1, vol2, d1, d2, rho,
                    r, K, t, Abs[i] * sqrt(2.0));

        N_x = cdfNormal(coor);

        N_x = N_x * Pds[i];
        integrK += N_x;
    }

    integrK = integrK / (sqrt(M_PI));

    return integrK;
}



double Proba1(double S1,
              double S2,
              double vol1,
              double vol2,
              double d1,
              double d2,
              double rho,
              double r,
              double K,
              double t)
{
    int n = 49;
    double N_x = 0.0;
    double integrS1 = 0.0;
    double coor;
    int i;


    // Calcul du coefficient associe a S1
    for (i = 0; i < n; i++)
    {
        coor = fonct_std_d1(S1, S2, vol1, vol2, d1, d2, rho, r, 
                            K, t, Abs[i] * sqrt(2.0));

        N_x = cdfNormal(coor);

        N_x = N_x * Pds[i];
        integrS1 += N_x;

    }

    integrS1 = integrS1 / (sqrt(M_PI));

    return integrS1;
}



double Proba2(double S1,
              double S2,
              double vol1,
              double vol2,
              double d1,
              double d2,
              double rho,
              double r,
              double K,
              double t)

{

    int n = 49;
    double N_x = 0.0;
    double integrS2 = 0.0;
    double coor;
    int i;

    // Calcul du coefficient associ ë ÇS2
    for (i = 0; i < n; i++)
    {
        coor = fonct_std_d2(S1, S2, vol1, vol2, d1, d2, rho, r,
                            K, t, Abs[i] * sqrt(2.0));

        N_x = cdfNormal(coor);

        N_x = N_x * Pds[i];
        integrS2 += N_x;

    }
    integrS2 = integrS2 / (sqrt(M_PI));

    return integrS2;
}






/*----------------------------------------------------------------------------*
  Formules analytiques pour le calcul d 'une option standard sur spread
  (Cas particulier de correlation egale a 1 ou - 1)
  Pay - Off = MAX((S2 - S1) - K; 0) pour call pc = 1
  Pay - Off = MAX(K - (S2 - S1); 0) pour put pc = -1

  Condition d 'exercice de l' option Spread
  {a1 * exp(x * b1) - a2 * exp(x * b2) - K < 0}
  A1 = S2 * exp(-vol2/2 * racine(Mat))
  A2 = S1 * exp(-vol1/2 * racine(Mat))
  B1 = vol2 * racine(Mat)
  B2 = vol1 * racine(Mat)
*----------------------------------------------------------------------------*/


double condition(double Var, double A1, double A2, double B1, 
                 double B2, double K)
{
    double x;

    x = A1 * exp(Var * B1) - A2 * exp(Var * B2) - K;

    return(x);
}


// Derivee de la fonction Condition d 'exercice de l' option Spread
// exp(b1 * x) *[a1 * b1 * exp{x*(b1 - b2) }-a2*b2]

double derivee_condition(double Var, double A1, double A2,
                         double B1, double B2)
{
    double x;

    x = exp(B1 * Var) * (A1 * B1 - A2 * B2 * exp(Var * (B2 - B1)));

    return(x);
}


// Fonction de calcul de la racine de la fonction condition a partir 
// d'un point donne et d' une precision donne
// Methode de Newton

double Newton(double Pt_depart, double epsilon, double A1, double A2,
              double B1, double B2, double B)
{
    double x;
    double y;

    y = Pt_depart;
    x = y - (condition(y, A1, A2, B1, B2, B)) 
        / derivee_condition(y, A1, A2, B1, B2);

    while (fabs(x - y) > epsilon)
    {
        y = x;
        x = y - (condition(y, A1, A2, B1, B2, B)) 
                 / derivee_condition(y, A1, A2, B1, B2);

    }

    return (x);
}


// Fonction de calcul de la racine de la fonction condition a partir
// de la valeur MAXimale de la racine;
// cas fonction croissante
// Methode de Dichotomie

double Dicho_croissante(double xmax, double epsilon, double A1,
                        double A2, double B1, double B2, double B)
{
    double xmin;
    double x;

    xmin = xmax - ECART_RACINE;
    x = (xmin + xmax) / 2;

    while (condition(x, A1, A2, B1, B2, B) > epsilon)
    {
        if (condition(x, A1, A2, B1, B2, B) > 0)
        {
            xmax = x;
            x = (xmin + xmax) / 2;
        }
        else
        {
            xmin = x;
            x = (xmin + xmax) / 2;
        }
    }
    return (x);
}


//Fonction de calcul de la racine de la fonction condition a partir
// de la valeur maximale de la racine;
// cas fonction decroissante
// Methode de Dichotomie

double Dicho_decroissante(double xmax, double epsilon, double A1,
                          double A2, double B1, double B2, double B)
{
    double xmin;
    double x;

    xmin = xmax - ECART_RACINE;
    x = (xmin + xmax) / 2;

    while (condition(x, A1, A2, B1, B2, B) > epsilon)
    {
        if (condition(x, A1, A2, B1, B2, B) > 0)
        {
            xmin = x;
            x = (xmin + xmax) / 2;
        }
        else
        {
            xmax = x;
            x = (xmin + xmax) / 2;
        }
    }
    return (x);
}




// Fonction Depart

double depart(double var, double A1, double A2, double B1, double B2, double B)
{
    double expo;

    expo = var - condition(var, A1, A2, B1, B2, B) 
                           / derivee_condition(var, A1, A2, B1, B2);

    expo = expo * MAX(B1, B2);
    while (expo > MAX_EXP)
    {
        var = var + alpha_start * MAX(fabs(var),0.0001);
        expo = var - condition(var, A1, A2, B1, B2, B) 
                           / derivee_condition(var, A1, A2, B1, B2);
        expo = expo * MAX(B1, B2);
    }

    return(var);
}




/*----------------------------------------------------------------------------*
  Calcul de l 'option Spread dans le cas ou la correlation est egale a -1
  on calcule sous les proba Q, Q1 et Q2 la proba que S2 - S1 > K
*----------------------------------------------------------------------------*/

double Probamin(double A1, double A2, double B1, double B2, double B)

  // Definition de B1, B2, A1, A2
  // La correlation est egale a -1 et B2 est toujours strictement
  // inferieur a 0
  // La fonction condition est strictement croissante;
  // il  n 'y a qu' une racine
  // les valeurs de A1,
  //    A2,
  //    B1 et B2 varient en fonction de la probabilite sous laquelle 
  // on travaille
{
        double x;
        double y;
        double prix;

        y = 0.0;
        x = condition(y, A1, A2, B1, B2, B) 
                      / derivee_condition(y, A1, A2, B1, B2);

        while (fabs(x - y) > ARM_Precision)
        {
            y = x;
            x = y - (condition(y, A1, A2, B1, B2, B)) 
                               / derivee_condition(y, A1, A2, B1, B2);
        }

        prix = cdfNormal(-x);

        return (prix);
}



/*----------------------------------------------------------------------------*
  Calcul de l'option Spread dans le cas ou la correlation est egale a 1
  on calcule sous les proba Q, Q1 et Q2 la proba que S2 - S1 > K
*----------------------------------------------------------------------------*/

double Probamax(double A1, double A2, double B1, double B2, double B)

//  Definition de B1, B2, A1, A2
//  La correlation est e gale a 1 et B2 est toujours positif
//  Il faudra distinguer les cas ou B2 > B1 et B2 < B1

{
    double x;
    double u,
        v;
    double y;
    double cond_u,
        cond_v;
    double r1;
    double r2;
    double prix;


    // u est le MAX de la fonction condition definie si B2 different de B1

    if (fabs(B2 - B1) > ARM_Precision)
    {
        u = log(A1 * B1 / (A2 * B2));
        u = u / (B2 - B1);

        // on utilise le fait que A1 et A2 > 0
        if (u > MAX_EXP)
        {
            if (B1 > B2)
                cond_u = INFINIT;
            else
                cond_u = -INFINIT;
        }
        else
            cond_u = condition(u, A1, A2, B1, B2, B);
    }


    // u est le MAX de la fonction derive_condition definie si B2 different de B1

    if ((fabs(B2 - B1) > ARM_Precision) && (u < MAX_EXP))
    {
        v = log(A1 * B1 * B1 / (A2 * B2 * B2));
        v = v / (B2 - B1);
        cond_v = condition(v, A1, A2, B1, B2, B);
    }


    if (B2 > B1)
    {
       // Premier cas:    B2 > B1
       // u est positif
       // Il peut y avoir 0, 1 ou 2 racines
       // La fonction condition est alors croissante puis
       // de croissante avec un MAX atteint en u
       // Lim(-inf) = -B et Lim(+inf) = -inf

       if (cond_u <= 0)
       {
          // Il n 'y a pas de racine, la fonction Condition est 
          // toujours negative
          prix = 0;
       }
       else
       {
          if (B <= 0)
            //Il n 'y a qu' une seule racine, Çdroite de u
            // Recherche de la racine
            // Le prix correspondra a la probabilite d etre 
            // inferieur a cette racine
            // Calcul de cette racine a partir de la fonction depart qui
            // assure la non explosion de la valeur passee
            // dans l 'exponentielle 
          {
               x = u + fabs(u) * alpha_start;
               x = depart(x, A1, A2, B1, B2, B);
               y = Newton(x, ARM_Precision, A1, A2, B1, B2, B);
               prix = cdfNormal(y);
          }
          else
             // Il y a deux racines
             // Le prix correspond a la probabilite d'etre compris 
             // entre ces deux racines
             // Recherche de la premiere racine a gauche de u,
             // a partir de v, point d 'inflexion de la fonction derive_condition
          {
             x = v;

             if (((v - condition(v, A1, A2, B1, B2, B) 
                  / derivee_condition(v, A1, A2, B1, B2)) * MAX(B1, B2)) 
                 < MAX_EXP
                )
                  y = Newton(x, ARM_Precision, A1, A2, B1, B2, B);
             else
                y = Dicho_croissante(u, ARM_Precision, A1, A2, B1, B2, B);

             r1 = y;

             // Recherche de la deuxieme racine a droite de u
             // Calcul de cette racine Ç partir de la fonction 
             // depart qui assure la non explosion de la valeur passee
             // dans l 'exponentielle 
             x = u + fabs(u) * alpha_start;
             x = depart(x, A1, A2, B1, B2, B);
             y = Newton(x, ARM_Precision, A1, A2, B1, B2, B);
             r2 = y;

             prix = cdfNormal(r2) - cdfNormal(r1);
          }
       }
    }
    else if (B1 > B2)
       // Deuxieme cas:B1 > B2
       // Il peut y avoir 0, 1 ou 2 racines
       // La fonction condition est alors decroissante puis croissante 
       // avec un MAX atteint en u < 0
       // Lim(-inf) = -B et Lim(+inf) = +inf
    {
        if (cond_u >= 0)
           // Il n 'y a pas de racine
           prix = 1;

        else
        {
            if (B >= 0)
               // Il n 'y a qu' une seule racine
               // Recherche de la racine
               // Le prix correspondra a la probabilite d'etre superieur a
               // cette racine 
               // Calcul de cette racine a partir de la fonction depart 
               // qui assure la non explosion de la valeur pass ë e
               // dans l 'exponentielle 
            {
                x = u + fabs(u) * alpha_start;
                x = depart(x, A1, A2, B1, B2, B);
                y = Newton(x, ARM_Precision, A1, A2, B1, B2, B);

                prix = cdfNormal(-y);
            }
            else
               // Il y a deux racines
               // Le prix correspond a la probabilite detre inferieur 
               // la premiere racine et d'etre superieur a la seconde
               // Recherche de la premiere racine a gauche de u, 
               // a partir de v, point d 'inflexion de 
               // la fonction derivee_condition
            {
                x = v;
                x = depart(x, A1, A2, B1, B2, B);
                // rajoute
                if (((v - condition(v, A1, A2, B1, B2, B) 
                    / derivee_condition(v, A1, A2, B1, B2)) * MAX(B1, B2)) 
                    < MAX_EXP)
                    y = Newton(x, ARM_Precision, A1, A2, B1, B2, B);
                else
                    y = Dicho_decroissante(u, ARM_Precision, A1, A2, B1, B2, B);

                r1 = y;

                // Recherche de la deuxieme racine a droite de u
                // Calcul de cette racine a partir de la fonction 
                // depart qui assure la non explosion de la valeur passee
                // dans l 'exponentielle 

                x = u + fabs(u) * alpha_start;
                x = depart(x, A1, A2, B1, B2, B);
                y = Newton(x, ARM_Precision, A1, A2, B1, B2, B);

                r2 = y;

                prix = cdfNormal(r1) + cdfNormal(-r2);
            }
        }
    }
    else
       // Troisieme cas:B1 = B2
       // Cas ou les volatilites sont egales
    {
        if (A1 > A2)
        {
            if (B > 0)
               // Il y a une racine et la fonction Condition est 
               // strictement croissante
            {
                y = log(B / (A1 - A2));
                y = y / B1;
                prix = cdfNormal(-y);
            }
            else
                // Il n 'y a pas de racine et la fonction Condition 
                // est strictement croissante et positive
                prix = 1;
        }
        else
        {
            if (A1 < A2)
            {
                if (B < 0)
                {
                    // Il y a une racine et la fonction Condition 
                    // est strictement d ë croissante
                    y = log(B / (A1 - A2));
                    y = y / B1;
                    prix = cdfNormal(y);
                }
                else
                    // Il n 'y a pas de racine et la fonction Condition 
                    // est strictement decroissante et negative
                        prix = 0;
            }

            else
               // Si S2 = S1
               // La condition Condition est toujours egale a -B
            {
                if (B < 0)
                    prix = 1;
                else
                    prix = 0;
            }
        }
    }

    return(prix);
}



double Pay_Off_Spread(double S1, double S2, double K, double d1, 
                      double d2, double Rate, double Mat, 
                      double ProbaS1, double ProbaS2, double ProbaK, int pc)
{
    switch (pc)
    {
        case 1:
            return(S2 * exp(-d2 * Mat) * ProbaS2 - S1 * exp(-d1 * Mat) 
                    * ProbaS1 - K * exp(-Rate * Mat) * ProbaK);
        case -1:
            return(-S2 * exp(-d2 * Mat) * (1 - ProbaS2) + S1 * exp(-d1 * Mat) 
               * (1 - ProbaS1) + K * exp(-Rate * Mat) * (1 - ProbaK));

        default:
            return -1;
    }
}


double Spread_opt_norm(double S1,
                       double S2,
                       double vol1,
                       double vol2,
                       double d1,
                       double d2,
                       double rho,
                       double r,
                       double K,
                       double t,
                       int pc)
{

    int n = 49;
    double N_x = 0.0;
    double ProbaK = 0.0;
    double ProbaS1 = 0.0;
    double ProbaS2 = 0.0;

    // Calcul du coefficient associ ë ÇK
    ProbaK = Proba0(S1, S2, vol1, vol2, d1, d2, rho, r, K, t);


    // Calcul du coefficient associ ë ÇS2
    ProbaS2 = Proba2(S1, S2, vol1, vol2, d1, d2, rho, r, K, t);


    // Calcul du coefficient associ ë ÇS1
    ProbaS1 = Proba1(S1, S2, vol1, vol2, d1, d2, rho, r, K, t);


    return Pay_Off_Spread(S1, S2, K, d1, d2, r, t, ProbaS1, ProbaS2, ProbaK, pc);
}



double Spread_opt_min(double S1,
                      double S2,
                      double vol1,
                      double vol2,
                      double d1,
                      double d2,
                      double Rate,
                      double K,
                      double Mat,
                      int pc)
{
    double B1,
        B2,
        A1,
        A2;

    double ProbaK,
        ProbaS1,
        ProbaS2;

    B1 = vol2 * sqrt(Mat);
    B2 = -vol1 * sqrt(Mat);
    A1 = S2 * exp((Rate - d2) * Mat - vol2 * vol2 * 0.5 * Mat);
    A2 = S1 * exp((Rate - d1) * Mat - vol1 * vol1 * 0.5 * Mat);
    ProbaK = Probamin(A1, A2, B1, B2, K);

    A1 = S2 * exp((Rate - d2) * Mat - vol2 * vol2 * 0.5 
         * Mat - vol1 * vol2 * Mat);

    A2 = S1 * exp((Rate - d1) * Mat + vol1 * vol1 * 0.5 * Mat);

    ProbaS1 = Probamin(A1, A2, B1, B2, K);


    A1 = S2 * exp((Rate - d2) * Mat + vol2 * vol2 * 0.5 * Mat);
    A2 = S1 * exp((Rate - d1) * Mat - vol1 * vol1 * 0.5 * Mat 
           - vol1 * vol2 * Mat);

    ProbaS2 = Probamin(A1, A2, B1, B2, K);

    return Pay_Off_Spread(S1, S2, K, d1, d2, Rate, Mat, 
                          ProbaS1, ProbaS2, ProbaK, pc);
}



double Spread_opt_max(double S1,
                      double S2,
                      double vol1,
                      double vol2,
                      double d1,
                      double d2,
                      double Rate,
                      double K,
                      double Mat,
                      int pc)
{
    double B1,
        B2,
        A1,
        A2;
    double ProbaK,
        ProbaS1,
        ProbaS2;

    B1 = vol2 * sqrt(Mat);
    B2 = vol1 * sqrt(Mat);
    A1 = S2 * exp((Rate - d2) * Mat - vol2 * vol2 * 0.5 * Mat);
    A2 = S1 * exp((Rate - d1) * Mat - vol1 * vol1 * 0.5 * Mat);
    ProbaK = Probamax(A1, A2, B1, B2, K);

    A1 = S2 * exp((Rate - d2) * Mat - vol2 * vol2 * 0.5 
                * Mat + vol1 * vol2 * Mat);

    A2 = S1 * exp((Rate - d1) * Mat + vol1 * vol1 * 0.5 * Mat);
    ProbaS1 = Probamax(A1, A2, B1, B2, K);

    A1 = S2 * exp((Rate - d2) * Mat + vol2 * vol2 * 0.5 * Mat);
    A2 = S1 * exp((Rate - d1) * Mat - vol1 * vol1 * 0.5 
                   * Mat + vol1 * vol2 * Mat);

    ProbaS2 = Probamax(A1, A2, B1, B2, K);

    return Pay_Off_Spread(S1, S2, K, d1, d2, Rate, Mat, ProbaS1, 
                          ProbaS2, ProbaK, pc);
}



double MargrabeSpreadOption(double S1,
                            double S2,
                            double vol1,
                            double vol2,
                            double d1,
                            double d2,
                            double Rho,
                            double Rate,
                            double Mat,
                            int pc)
{
    double Sigma,
        d,
        value;

    Sigma = sqrt(Mat * (vol1 * vol1 + vol2 * vol2 - 2 * Rho * vol1 * vol2));
    d = log(S2 * exp(-d2 * Mat) / (S1 * exp(-d1 * Mat)));
    d = d / Sigma + 0.5 * Sigma;

    switch (pc)
    {
        case 1:
            value = S2 * exp(-d2 * Mat) * cdfNormal(d) 
                      - S1 * exp(-d1 * Mat) * cdfNormal(d - Sigma);
        break;

        case -1:
            value = -S2 * exp(-d2 * Mat) * cdfNormal(-d) 
                       + S1 * exp(-d1 * Mat) * cdfNormal(-d + Sigma);
        break;

        default:
            value = -1.0;
        break;
    }

    return value;
}


#ifdef WIN32
double SpreadOption_ClosedForms(double S1, double S2,
								double vol1, double vol2,
								double Rho, double K,
								double Mat, int pc)
{
	return 0;//tbellaj Export_LogNormal_SpreadOption(S1,S2,vol1,vol2,Rho,K,Mat,pc,ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION ,175);
}

#endif


double SpreadOption(double S1,
                    double S2,
                    double vol1,
                    double vol2,
                    double d1,
                    double d2,
                    double Rho,
                    double Rate,
                    double K,
                    double Mat,
                    int pc,
					int ComputedFormula /* 0 formule d'Olivier, 1 Formule Classique*/)
{
#ifdef WIN32
	if (ComputedFormula == 0)
		return SpreadOption_ClosedForms(S2,S1,vol2,vol1,Rho,K,Mat,pc);
#endif

    double ProbaK, ProbaS1, ProbaS2;

    double borne = 0.99;
    double EPS = 0.001;

    double prix;

    // on utilise la formule Margrabe pour le cas particulier ou le
    // strike est nul

    if ( fabs(K) < K_DOUBLE_TOL )
    {
       return MargrabeSpreadOption(S1, S2, vol1, vol2, d1, d2, Rho,
                                   Rate, Mat, pc);
    }

    if ( fabs(Rho) < borne )
    {
       double value = Spread_opt_norm(S1, S2, vol1, vol2, d1, d2, Rho, Rate, K, Mat, pc);

       return(value); 
    }
    else
    {
        double B1,
               B2,
               A1,
               A2;

        if (( Rho >= -1.0 ) && ( Rho <= -borne ))
        {
            // on price avec une correlation de rho=- 1
            double prix_borne, prix_borne_pente, prix_borne_pente2;
/*

            B1 = vol2*sqrt(Mat);
            B2 = -vol1*sqrt(Mat);

            A1 = S2*exp((Rate-d2)*Mat-vol2*vol2*0.5*Mat);
            A2 = S1*exp((Rate-d1)*Mat-vol1*vol1*0.5*Mat);
            ProbaK = Probamin(A1, A2, B1, B2, K);

            A1 = S2*exp((Rate-d2)*Mat-vol2*vol2*0.5*Mat+vol1*vol2*Mat);
            A2 = S1 * exp((Rate-d1)*Mat+vol1*vol1*0.5*Mat);
            ProbaS1 = Probamin(A1, A2, B1, B2, K);

            A1 = S2*exp((Rate-d2)*Mat+vol2*vol2*0.5*Mat);
            A2 = S1*exp((Rate-d1)*Mat-vol1*vol1*0.5*Mat+vol1*vol2*Mat);

            ProbaS2 = Probamin(A1, A2, B1, B2, K);

*/
            //prix_rhoMoins1   = Pay_Off_Spread(S1, S2, K, d1, d2, Rate, Mat, ProbaS1, ProbaS2, ProbaK, pc);

            prix_borne        = Spread_opt_norm(S1, S2, vol1, vol2, d1, d2, -borne, Rate, K, Mat, pc);
         
            prix_borne_pente  = Spread_opt_norm(S1, S2, vol1, vol2, d1, d2, -borne+EPS, Rate, K, Mat, pc);

            prix_borne_pente2 = Spread_opt_norm(S1, S2, vol1, vol2, d1, d2, -borne+2.0*EPS, Rate, K, Mat, pc);


            prix = QuadraticLagrangeInterpol(Rho,
                                             -borne, prix_borne,
                                             -borne+EPS, prix_borne_pente,
                                             -borne+2.0*EPS, prix_borne_pente2);

            return(prix);
        }
        else
        {
            //(Rho <= 1.0) & (Rho >= 0.99)

            //rho=1

            B1 = vol2*sqrt(Mat);
            B2 = vol1*sqrt(Mat);

            A1 = S2*exp((Rate-d2)*Mat-vol2*vol2*0.5*Mat);
            A2 = S1*exp((Rate - d1)*Mat-vol1*vol1*0.5*Mat);
            ProbaK = Probamax(A1, A2, B1, B2, K);
            //LimProbaK = Proba0(S1, S2, vol1, vol2, d1, d2, borne, Rate, K, Mat);

            //ProbaK = LimProbaK+(ProbaK-LimProbaK) 
            //                   /(1.0-borne)*(Rho-borne);

            A1 = S2 * exp((Rate - d2) * Mat - vol2 * vol2 
                           * 0.5 * Mat + vol1 * vol2 * Mat);
            A2 = S1 * exp((Rate - d1) * Mat + vol1 * vol1 * 0.5 * Mat);
            ProbaS1 = Probamax(A1, A2, B1, B2, K);
            //LimProba1 = Proba1(S1, S2, vol1, vol2, d1, d2, borne, Rate, K, Mat);
            //ProbaS1 = LimProba1 + (ProbaS1 - LimProba1) 
            //                       / (1.0 - borne) * (Rho - borne);

            A1 = S2 * exp((Rate - d2) * Mat + vol2 * vol2 * 0.5 * Mat);
            A2 = S1 * exp((Rate - d1) * Mat - vol1 * vol1 * 0.5 
                             * Mat + vol1 * vol2 * Mat);
            ProbaS2 = Probamax(A1, A2, B1, B2, K);
            //LimProba2 = Proba2(S1, S2, vol1, vol2, d1, d2, borne, Rate, K, Mat);
            //ProbaS2 = LimProba2 + (ProbaS2 - LimProba2) 
            //               / (1.0 - borne) * (Rho - borne);

      

            double prix_borne_pente = Spread_opt_norm(S1, S2, vol1, vol2, d1, d2, 
                                                      borne-EPS, Rate, K, Mat, pc);

            double prix_borne       = Spread_opt_norm(S1, S2, vol1, vol2, d1, d2, 
                                                      borne, Rate, K, Mat, pc);

            double prix_rho1 = Pay_Off_Spread(S1, S2, K, d1, d2, Rate, Mat, 
                                              ProbaS1, ProbaS2, ProbaK, pc);

            prix = QuadraticLagrangeInterpol(Rho,
                                             borne-EPS, prix_borne_pente,
                                             borne, prix_borne,
                                             1.0, prix_rho1);
            return(prix);
        }
    }
}




double GaussianSpreadOption(double S1,
                            double S2,
                            double SpreadVol,
                            double d1,
                            double d2,
                            double Rho,
                            double Rate,
                            double K,
                            double Mat,
                            int pc)
{
    //Pricing en faisant l 'hypothese que le spread est gaussien
    double price = 0.0,
        d = 0.0;

    d = pc * (S2 * exp(-d2 * Mat) - S1 * exp(-d1 * Mat) 
              - K * exp(-Rate * Mat)) / (SpreadVol * sqrt(Mat));

    price = pc * (S2 * exp(-d2 * Mat) - S1 * exp(-d1 * Mat) 
              - K * exp(-Rate * Mat)) * cdfNormal(d);
    price = price + SpreadVol * sqrt(Mat) * exp(-d * d * 0.5) / sqrt(2 * M_PI);

    return price;
}



double generic_spreadoption_price(double S1, double S2,
                    marginal_descr *m1, marginal_descr *m2, jointdistrib_descr *jd,
					double Rate, double K, double Mat, int pc, double w1, double w2, 
					int ModelType, int ComputedFormula, int NumSteps)
{
	double price = 0.0;

	//if (ModelType == K_2SABR)
	//{
	//	struct sabr_descr* F1 =   (struct sabr_descr*) m1;
	//	struct sabr_descr* F2 =   (struct sabr_descr*) m2;
	//	struct bi_sabr* JD = (struct bi_sabr*) jd;

	//	double NewS1 = S1 * w1;
	//	double NewS2 = S2 * w2;
	//	double NewAlpha1 = F1->alpha*pow(w1,1.0 - F1->beta);
	//	double NewAlpha2 = F2->alpha*pow(w2,1.0 - F2->beta);
	//	double NewStrike = K;

	//	// Cap / Floor traite par export

	//	price = 0 ;/*Packaged_BiSABR_SpreadOption(NewS2, NewAlpha2, F2->beta, F2->rho, F2->nu, NewS1, NewAlpha1, F1->beta, F1->rho, F1->nu,
	//	 			JD->rho_s1_s2, JD->rho_vol1_vol2, JD->rho_s2_vol1, JD->rho_s1_vol2, NewStrike, Mat, pc, ComputedFormula);*/

	//}
	//
	//if (ModelType == K_GCOP) // In fact Sabr's marginal mixed in a gaussian Copula
	//{
	//	struct sabr_descr* F1 =   (struct sabr_descr*) m1;
	//	struct sabr_descr* F2 =   (struct sabr_descr*) m2;
	//	struct bi_sabr* JD = (struct bi_sabr*) jd;

	//	if (ComputedFormula == 0)
	//	{
	//		if (pc == K_CAP) // Cap Case
	//		{
	//			price = Export_Gaussian_SABR_Power_SpreadOption(S2, S1,
	//							   F2->alpha, F2->beta, F2->rho, F2->nu,
	//							   F1->alpha, F1->beta, F1->rho, F1->nu,
	//							   JD->rho_s1_s2, Mat, K_SABR_IMPLNVOL,
	//							   w2, -w1, K, w2, -w1, K, NumSteps, 0.01, 1.5, 0.02);
	//		}
	//		else // Floor Case
	//		{
	//			price = Export_Gaussian_SABR_Power_SpreadOption(S2, S1,
	//							   F2->alpha, F2->beta, F2->rho, F2->nu,
	//							   F1->alpha, F1->beta, F1->rho, F1->nu,
	//							   JD->rho_s1_s2, Mat, K_SABR_IMPLNVOL,
	//							   -w2, w1, -K, -w2, w1, -K, NumSteps,0.01, 1.5, 0.02);
	//							   // -w2, w1, K, -w2, w1, -K, NumSteps);
	//		}
	//	}
	//	else // nouvelle formule
	//	{
	//		double norm_S2, norm_alpha2, norm_S1, norm_alpha1, norm_K;

	//		norm_S1 = w1*S1;
	//		norm_alpha1 = F1->alpha * pow(w1,1.- F1->beta);
	//		norm_S2 = w2*S2;
	//		norm_alpha2 = F2->alpha * pow(w2,1.- F2->beta);
	//		norm_K = K;

	//		
	//		price = Export_GaussianSABRDigitalCallPayingS1(norm_S2, norm_alpha2, F2->beta, F2->rho, F2->nu, K_SABR_IMPLNVOL,
	//			                                       norm_S1, norm_alpha1, F1->beta, F1->rho, F1->nu, K_SABR_IMPLNVOL,
	//												   JD->rho_s1_s2, K, Mat, NumSteps,0.01, 1.5, 0.02);
	//												   /// JD->rho_s1_s2, norm_K, Mat, NumSteps);
	//		price -= Export_GaussianSABRDigitalCallPayingS2(norm_S2, norm_alpha2, F2->beta, F2->rho, F2->nu, K_SABR_IMPLNVOL,
	//			                                       norm_S1, norm_alpha1, F1->beta, F1->rho, F1->nu, K_SABR_IMPLNVOL,
	//												   JD->rho_s1_s2, K, Mat, NumSteps,0.01, 1.5, 0.02);
	//												   /// JD->rho_s1_s2, norm_K, Mat, NumSteps);
	//		price -= K * Export_GaussianSABRDigitalCall(norm_S2, norm_alpha2, F2->beta, F2->rho, F2->nu, K_SABR_IMPLNVOL,
	//			                                       norm_S1, norm_alpha1, F1->beta, F1->rho, F1->nu, K_SABR_IMPLNVOL,
	//												   JD->rho_s1_s2, K, Mat, NumSteps,0.01, 1.5, 0.02);
	//												   /// JD->rho_s1_s2, norm_K, Mat, NumSteps);

	//		if (pc == K_FLOOR) // Floor dealed with Call-Put Parity
	//		{
	//			price = price - (norm_S2 - norm_S1 - K);
	//		}
	//	}

	//}	

	//if (ModelType == K_2LOG)
	//{
	//	
	//	struct one_parameter_bivariate* JD = (struct one_parameter_bivariate*) jd;

	//	double vol1 = m1->volatility;
	//	double vol2 = m2->volatility;
	//	double Rho = JD->correl;

	//	price = SpreadOption(S1,
	//						S2,
	//						vol1,
	//						vol2,
	//						0., //d1,
	//						0., //d2,
	//						Rho,
	//						Rate,
	//						K,
	//						Mat,
	//						pc,
	//						ComputedFormula /* 0 formule d'Olivier, 1 Formule Classique*/);

	//}

	return price;
}


double generic_DigitalSO_price(double S1, double S2,
                    marginal_descr *m1, marginal_descr *m2, jointdistrib_descr *jd,
					double K, double Mat, int pc, double w1, double w2, 
					int ModelType, int ComputedFormula, int NumSteps)
{
	double price = 0.0;

	/*if (ModelType == K_2SABR)
	{
		struct sabr_descr* F1 =   (struct sabr_descr*) m1;
		struct sabr_descr* F2 =   (struct sabr_descr*) m2;
		struct bi_sabr* JD = (struct bi_sabr*) jd;

		double NewS1 = 0.01 * S1 * w1;
		double NewS2 = 0.01 * S2 * w2;
		double NewAlpha1 = F1->alpha*pow(w1,1.0 - F1->beta);
		double NewAlpha2 = F2->alpha*pow(w2,1.0 - F2->beta);
		double NewStrike = 0.01 * K;


		price = Export_BiSABR_Digital_SpreadOption(NewS2, NewAlpha2, F2->beta, F2->rho, F2->nu, NewS1, NewAlpha1, F1->beta, F1->rho, F1->nu,
		 			NewStrike, Mat, pc, JD->rho_s1_s2, JD->rho_vol1_vol2, JD->rho_s2_vol1, JD->rho_s1_vol2, ComputedFormula);
	}*/

	return price;

}


double generic_DigitalSO_PayingIndex_price(double S1, double S2, double S3,
                    marginal_descr *m1, marginal_descr *m2, marginal_descr *m3, jointdistrib_descr *jd,
					double K, double Mat, int pc, double w1, double w2, 
					int ModelType, int ComputedFormula, int NumSteps)
{
	double price = 0.0;

	//if (ModelType == K_2SABR)
	//{
	//	struct sabr_descr* F1 =   (struct sabr_descr*) m1;
	//	struct sabr_descr* F2 =   (struct sabr_descr*) m2;
	//	struct bi_sabr_payIndex* JD = (struct bi_sabr_payIndex*) jd;

	//	double NewS1 = 0.01 * S1 * w1;
	//	double NewS2 = 0.01 * S2 * w2;
	//	double NewAlpha1 = F1->alpha*pow(w1,1.0 - F1->beta);
	//	double NewAlpha2 = F2->alpha*pow(w2,1.0 - F2->beta);
	//	double NewStrike = 0.01 * K;
	//	double NewS3 = 0.01 * S3;

	//	// attention les index 1 et 2 sont intervertis entre Kernel et ClosedForms
	//	// c-à-d Payoff Kernel = Index2 - Index1 alors que Payoff ClosedForms = Index1 - Index2

	//	// struct bi_sabr* joint_s1_s2 = (struct bi_sabr*) JD->joint_s1_s2;
	//	double rho_s1_s2 = JD->joint_s1_s2->rho_s1_s2;
	//	double rho_vol1_vol2 = JD->joint_s1_s2->rho_vol1_vol2;
	//	double rho_s1_vol2 = JD->joint_s1_s2->rho_s1_vol2;
	//	double rho_s2_vol1 = JD->joint_s1_s2->rho_s2_vol1;
	//	double sigma30 = m3->volatility;
	//	double rho13 = JD->rho_s2_s3;
	//	double rho23 = JD->rho_s1_s3;

	//	
	//	if (fabs(S2 - S3) < 1e-06)
	//	{
	//		price = Export_BiSABR_Digital_SpreadOption_PayS1(
	//								   NewS2, NewAlpha2, F2->beta, F2->rho, F2->nu,
	//								   NewS1, NewAlpha1, F1->beta, F1->rho, F1->nu,
	//								   NewStrike, Mat, pc,
	//								   rho_s1_s2, rho_vol1_vol2, rho_s2_vol1, rho_s1_vol2, ComputedFormula);
	//	}

	//	if (fabs(S1 - S3) < 1e-06)
	//	{
	//		price = Export_BiSABR_Digital_SpreadOption_PayS2(
	//								   NewS2, NewAlpha2, F2->beta, F2->rho, F2->nu,
	//								   NewS1, NewAlpha1, F1->beta, F1->rho, F1->nu,
	//								   NewStrike, Mat, pc,
	//								   rho_s1_s2, rho_vol1_vol2, rho_s2_vol1, rho_s1_vol2, ComputedFormula);
	//	}

	//	if ((S3 != S1) && (S3 != S2))
	//	{
	//		
	//		price = Export_BiSABR_Digital_SpreadOption_PayS3(
	//								   NewS2, NewAlpha2, F2->beta, F2->rho, F2->nu,
	//								   NewS1, NewAlpha1, F1->beta, F1->rho, F1->nu,
	//								   rho_s1_s2, rho_vol1_vol2, rho_s2_vol1, rho_s1_vol2,
	//								   NewS3, sigma30, rho13, rho23,
	//								   NewStrike, Mat, pc, ComputedFormula);
	//	}

	//	price = price * 100.;

	//}
	//

	return price;
}

/*---------------------------------------------------------------------*/
/*---- End Of File ----*/
