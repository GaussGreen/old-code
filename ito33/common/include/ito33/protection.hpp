// ITO 33
// PROTECTION.HPP
// Création : 6 avril 1999
// Mise à jour : 20 septembre 1999
//-------------------------------------------------------------------------------------------------

#ifndef PROTECTION_HPP
#define PROTECTION_HPP

int SetProtection(char *NomFichier, int Perturbation, char *_pcCheminRacine,
  char *_pcCheminRegistry);
int GetProtection(char *NomFichier, int Perturbation, char *_pcCheminRacine);
int GetFonctionCryptee(char *SignatureCryptee, char *DateHeureCryptee, char *Perturbation,
  char *Result, int &_Signature, char *_pch);

void GetFonctionCryptee(char *SignatureCryptee, char *DateHeureCryptee, char *Perturbation,
  char *Result);

int GenerationCle(char *CleSignal, char *CleTemps, int Annee, int Mois, int Jour, char *Perturbation, 
        char *CleITO, char* CleITO2, int &Signature, char* MyDate);

int ClientSideFile(char *NomFichier, int Perturbation, char *_pcCheminRacine);

int OurSideFile(char *NomFichier, char* Perturbation);

int HostOurSideFile(int Signature, char *NomFichier, char* Perturbation);

#endif
