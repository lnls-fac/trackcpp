/***************************************************************************
                          naffutils.h  -  description
                             -------------------
    begin                : Wed May 25 20:24:44 CET 2005
    copyright            : (C) 2002 by Laurent Nadolski
    email                : nadolski@synchrotron-soleil.fr

    modified by:         : Accelerator Physics Group - LNLS !!!

 ***************************************************************************/

#include "pos.h"
#include <vector>

/* Frequency Map Analysis */
void naff_run(const std::vector<Pos<double>>& data, double& tunex, double& tuney);

//void Get_NAFF(int nterm, long ndata, const std::vector<Pos<double>>& Tab, double *fx, double *fz, int nb_freq[2]);


/***************************************************************************
                          modnaff.h  -  description
                             -------------------
    begin                : Fri Dec 22 20:24:44 CET 2002
    copyright            : (C) 2002 by Laurent Nadolski
    email                : nadolski@synchrotron-soleil.fr

 Origine : Bureau des Longitudes (IMCCE)
           J. Laskar, M. Gastineau

 ***************************************************************************/

#ifndef MODNAFF_H
#define MODNAFF_H

#include <cfloat>
#include <stdio.h>
#include <cstdlib>
#include <cmath>

/*--------*/
/* define */
/*--------*/
#define BOOL int
#define FALSE 0
#define TRUE  1
#define Myyerror printf
#define MAX(x,y) ((x)<(y)?(y):(x))
#define INLINE_EXTERN static inline


/*---------------------*/
/*routine d'allocation */
/*---------------------*/
#define SYSCHECKMALLOCSIZE(variable, type, size) \
  if ((variable=(type*)malloc(sizeof(type)*(size)))==NULL)\
  { printf("error : malloc failed!\n"); exit(1);}
#define SYSCHECKMALLOC(variable, type) \
  if ((variable=(type*)malloc(sizeof(type)))==NULL)\
  { printf("error : malloc failed!\n"); exit(1);}
#define DIM2(prow,row,col,type,l) \
  {\
  register type *pdata;\
  int I;\
  SYSCHECKMALLOCSIZE(pdata, type,(row) * (col));\
  SYSCHECKMALLOCSIZE(prow, type *,(row));\
  for (I = 0; I < (row); I++)\
     {\
     prow[I] = pdata;\
     pdata += col;\
     }\
  }
#define SYSFREE(variable) free(variable)
#define HFREE2(variable) {SYSFREE(*variable); SYSFREE(variable);}

struct complexe {
  double reel,imag;
};
typedef struct complexe t_complexe;

/* v0.96 M. GASTINEAU 18/12/98 : ajout */
/*liste des frequences pour NAF */
struct list_fenetre_naf {
 double                   dFreqMin; /* frequence minimale */
 double                   dFreqMax; /* frequence maximale */
 int                      iNbTerme; /* nombre de termes a rechercher */
 struct list_fenetre_naf *suivant; /*fenetre suivante */
};
/* v0.96 M. GASTINEAU 18/12/98 : fin ajout */
typedef struct list_fenetre_naf t_list_fenetre_naf; /* v0.96 M. GASTINEAU 18/12/98 : ajout */

/*v0.96 M. GASTINEAU 04/09/98 : ajout pour la gestion de naf */
/* pour le role de ces champs, cf. modnaff.c */
struct stnaf {
 /*champ utilise par modnaff.c */
 FILE *NFPRT;
 double EPSM;
 int NTERM,KTABS,NFS;
 int ICPLX,IW,ISEC;
 int NERROR,IPRT;

 double *TFS;
 t_complexe *ZAMP;
 t_complexe **ZALP; /*tableau a deux dimensions*/
 t_complexe *ZTABS;

 double DTOUR,UNIANG,FREFON;
 double XH,T0;

 /*autre champ utilise en tant que flag */
 double m_dneps; /*equivaut a DNEPS */
 int m_iNbLineToIgnore; /*nombre de lignes a ignorer en debut du fichier des solutions */
 BOOL m_bFSTAB; /*=TRUE => sauve le tableau ZTABS.*/
 /* v0.96 M. GASTINEAU 06/01/99 : ajout */
 t_list_fenetre_naf *m_pListFen; /*liste des fenetres */
 /* v0.96 M. GASTINEAU 06/01/99 : fin ajout */

 /* Ximenes XRR 2015-08-20, trying to get rid of all global variables */
 double* TWIN;
 double  AF,BF;

};
typedef struct stnaf t_naf;
/*v0.96 M. GASTINEAU 04/09/98 : fin ajout */


/*-----------------*/
/* public functions*/
/*-----------------*/
void naf_initnaf_notab();
void naf_cleannaf_notab();
void naf_initnaf(t_naf& g_NAFVariable);
void naf_cleannaf(t_naf& g_NAFVariable);
BOOL naf_mftnaf(t_naf& g_NAFVariable, int NBTERM, double EPS);
void naf_prtabs(t_naf& g_NAFVariable, int KTABS, t_complexe *ZTABS, int IPAS);
void naf_smoy(t_naf& g_NAFVariable, t_complexe *ZM);

#endif



/***************************************************************************
                          complexe.h  -  description
                             -------------------
    begin                : Fri Dec 22 20:24:44 CET 2002
    copyright            : (C) 2002 by Laurent Nadolski
    email                : nadolski@synchrotron-soleil.fr

 Origine : Bureau des Longitudes (IMCCE)
           J. Laskar, M. Gastineau

 ***************************************************************************/


/* Astronomie et Systemes dynamiques                            */
/* M. GASTINEAU  Bureau des Longitudes, Paris 15/07/98          */
/* ensemble de fonctions inlines traitant uniquement les complexes */
/* v0.97 M. GASTINEAU 15/02/99: ajout des fonctions trigonometrique*/
/* v0.97 M. GASTINEAU 22/03/99: ajout de i_compl_powreel        */

/* Laurent Nadolski : compilation SOLEIL */
/* 8 decembre 2002 */
/* BUG avec gcc, module et phase d'un nombre complexe fausse */
/* car prototype de la fonction inconnu */
/* Prototype des fonctions mis dans un fichier .h propre */
/** Copier coller mis dans complexeheader **/

/* Fin ajout */

/*----------------IMPLEMENTATION--------------------------------*/

/*----------------i_compl_cmplx---------------------------------*/
/*Construit un complexe a partir de deux doubles(retourne c=a+i*b)*/
/*--------------------------------------------------------------*/
 void i_compl_cmplx(t_complexe *c, double a,double b);
/*----------------module----------------------------------------*/
/* retourne le module du nombre complexe c                      */
/*--------------------------------------------------------------*/
 double    i_compl_module(t_complexe c);
 double    i_compl_angle(t_complexe c);
/*----------------i_compl_add-----------------------------------*/
/* addition de deux nombres complexes c1 et c2 : c1+c2          */
/*--------------------------------------------------------------*/
 t_complexe i_compl_add(const t_complexe c1,const t_complexe c2);
/*----------------i_compl_padd----------------------------------*/
/* addition de deux nombres complexes c1 et c2 : c1+=c2         */
/*--------------------------------------------------------------*/
 void i_compl_padd(t_complexe *c1,t_complexe *c2);
/*----------------i_compl_paddconst-----------------------------*/
/* addition de deux nombres complexes c1 et c2 : c1+=c2         */
/*--------------------------------------------------------------*/
 void i_compl_paddconst(t_complexe *c1,t_complexe c2);
/*----------------i_compl_mul----------------------------------*/
/* multiplication de deux nombres complexes c1 et c2 : c1*c2   */
/*-------------------------------------------------------------*/
 t_complexe i_compl_mul(const t_complexe c1,const t_complexe c2);
/*----------------i_compl_pmul---------------------------------*/
/* multiplication de deux nombres complexes c1 et c2 : c1*=c2  */
/*-------------------------------------------------------------*/
 void i_compl_pmul(t_complexe *c1,t_complexe *c2);
/*----------------i_compl_muldoubl-----------------------------*/
/* multiplication d'un double c1 d'un nombre complexe c2: c1*c2*/
/*-------------------------------------------------------------*/
 t_complexe i_compl_muldoubl(const double c1, const t_complexe c2);
/*----------------i_compl_pmuldoubl----------------------------*/
/* multiplication d'un double c1 d'un nombre complexe c2: c2*=c1*/
/*-------------------------------------------------------------*/
 void i_compl_pmuldoubl(t_complexe *c2, double* const c1);
/*----------------i_compl_pdivdoubl----------------------------*/
/* division d'un double c1 d'un nombre complexe c2: c2/=c1     */
/*-------------------------------------------------------------*/
 void i_compl_pdivdoubl(t_complexe *c2, double *c1);
/*----------------i_compl_conj---------------------------------*/
/*retourne le conjugue de c1                                   */
/*-------------------------------------------------------------*/
 t_complexe i_compl_conj(t_complexe *c1);
/*----------------i_compl_pdiv---------------------------------*/
/* division de deux complexes : c1/=c2                         */
/*-------------------------------------------------------------*/
 void i_compl_pdiv(t_complexe * const c1,const t_complexe * const c2);
/*----------------i_compl_div2d--------------------------------*/
/* division d'un reel par un complexe c1/(c2_r+i*c2_i)         */
/*-------------------------------------------------------------*/
t_complexe i_compl_div2d(const double c1,const double c2_r,const double c2_i);
/*----------------i_compl_div4d--------------------------------*/
/* division de deux complexes : (c1_r+i*c1_i)/(c2_r+i*c2_i)    */
/*-------------------------------------------------------------*/
t_complexe i_compl_div4d(register const double c1_r,register const double c1_i,
                         register const double  c2_r,register const double c2_i);
/*----------------i_compl_div----------------------------------*/
/* division de deux complexes : c1/c2                          */
/*-------------------------------------------------------------*/
t_complexe i_compl_div(const t_complexe c1,const t_complexe c2);
/*----------------i_compl_pow2d--------------------------------*/
/*retourne la puissance (p_r+i*p_i)**b                         */
/*-------------------------------------------------------------*/
t_complexe i_compl_pow2d(const double p_r, const double p_i,int n);
/*----------------i_compl_powreel------------------------------*/
/*retourne la puissance a**b avec b reel                       */
/*-------------------------------------------------------------*/
/* v0.97 M. GASTINEAU 22/03/99: ajout de i_compl_powreel       */
t_complexe i_compl_powreel(const t_complexe z,double p_dK);
/*----------------i_compl_pow----------------------------------*/
/*retourne la puissance a**b                                   */
/*-------------------------------------------------------------*/
 t_complexe i_compl_pow(const t_complexe a,int n);
/*----------------expcomplexe-----------------------------------*/
/* calcule et retourne exp(c1)                                  */
/*--------------------------------------------------------------*/
/* v0.96 M. GASTINEAU 04/09/98 : ajout */
t_complexe i_compl_exp(t_complexe c1);
/*----------------i_compl_psub----------------------------------*/
/* soustrait de deux nombres complexes c1 et c2 : c1-=c2        */
/*--------------------------------------------------------------*/
/*v0.96 M. GASTINEAU 01/12/98 : ajout */
void i_compl_psub(t_complexe *c1,t_complexe *c2);
/*----------------i_compl_sub----------------------------------*/
/* soustrait de deux nombres complexes c1 et c2 : c1-c2        */
/*-------------------------------------------------------------*/
/*v0.96 M. GASTINEAU 01/12/98 : ajout */
t_complexe i_compl_sub(t_complexe c1,t_complexe c2);
/*----------------i_compl_cos----------------------------------*/
/* retourne le cosinus de c1 (cos x  cosh y  -  i sin x sinh y)*/
/*-------------------------------------------------------------*/
/* v0.97 M. GASTINEAU 15/02/99: ajout */
t_complexe i_compl_cos(t_complexe c1);
/*----------------i_compl_sin----------------------------------*/
/* retourne le sinus de c1 (= sin x  cosh y  +  i cos x sinh y)*/
/*-------------------------------------------------------------*/
/* v0.97 M. GASTINEAU 15/02/99: ajout */
t_complexe i_compl_sin(t_complexe c1);
/*----------------i_compl_cosh---------------------------------*/
/* retourne le cosinus hyperbolique de c1                      */
/* (= cosh x  cos y  +  i sinh x sin y)                        */
/*-------------------------------------------------------------*/
/* v0.97 M. GASTINEAU 15/02/99: ajout */
t_complexe i_compl_cosh(t_complexe c1);
/*----------------i_compl_sinh---------------------------------*/
/* retourne le sinus hyperbolique de c1                        */
/* (= sinh x  cos y  +  i cosh x sin y)                        */
/*-------------------------------------------------------------*/
/* v0.97 M. GASTINEAU 15/02/99: ajout */
t_complexe i_compl_sinh(t_complexe c1);
/*----------------i_compl_tan----------------------------------*/
/* retourne la tangente de c1                                  */
/*-------------------------------------------------------------*/
/* v0.97 M. GASTINEAU 15/02/99: ajout */
t_complexe i_compl_tan(t_complexe c1);
/*----------------i_compl_tanh---------------------------------*/
/* retourne la tangente hyperbolique de c1                     */
/*-------------------------------------------------------------*/
/* v0.97 M. GASTINEAU 15/02/99: ajout */
t_complexe i_compl_tanh(t_complexe c1);
/*----------------i_compl_log----------------------------------*/
/* retourne le logarithme de c1                                */
/*-------------------------------------------------------------*/
/* v0.97 M. GASTINEAU 15/02/99: ajout */
t_complexe i_compl_log(t_complexe c1);
/*----------------i_compl_log10--------------------------------*/
/* retourne le logarithme  base 10 de c1                       */
/*-------------------------------------------------------------*/
/* v0.97 M. GASTINEAU 15/02/99: ajout */
t_complexe i_compl_log10(t_complexe c1);
