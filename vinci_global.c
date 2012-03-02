/****************************************************************************************/
/*                                                                                      */
/*                                      vinci_global.c                                  */
/*                                                                                      */
/****************************************************************************************/
/*                                                                                      */
/* Authors: Benno Bueeler (bueeler@ifor.math.ethz.ch)                                   */
/*          and                                                                         */
/*          Andreas Enge (enge@ifor.math.ethz.ch)                                       */
/*          Institute for Operations Research                                           */
/*	    Swiss Federal Institute of Technology Zurich                                */
/*	    Switzerland                                                                 */
/*                                                                                      */
/* Last Changes: February 4, 2001                                                       */
/*                                                                                      */
/****************************************************************************************/
/*                                                                                      */
/* definitions of the global variables                                                  */
/*                                                                                      */
/****************************************************************************************/

#include "vinci.h"

/****************************************************************************************/

int G_d = 0;
int G_m = 0;
int G_n = 0;

real **G_Hyperplanes = NULL;
boolean **G_Incidence = NULL;
T_VertexSet G_Vertices;

int G_Storage = -1;
int G_RandomSeed = 4;

rational G_Minus1 = -1;

#ifdef STATISTICS
   unsigned int Stat_Count;
   real Stat_Smallest, Stat_Biggest;
   unsigned int Stat_CountPos [STAT_BIGGEST_EXP - STAT_SMALLEST_EXP + 3];
   unsigned int Stat_CountNeg [STAT_BIGGEST_EXP - STAT_SMALLEST_EXP + 3];
   unsigned int *Stat_CountStored = NULL, *Stat_CountRetrieved = NULL;
   unsigned int Stat_CountShifts;
   long int Stat_ActualMem = 0;
   long int Stat_MaxMem = 0;
#endif

/****************************************************************************************/
