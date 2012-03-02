/****************************************************************************************/
/*                                                                                      */
/*                                      vinci_screen.c                                  */
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
/* some functions for outputting data on the scrren, mainly for debugging information   */
/* The functions demand a parameter designating the file where the output has to be     */
/* directed. In general this will be stdout or stderr.                                  */
/*                                                                                      */
/****************************************************************************************/

#include "vinci.h"

/****************************************************************************************/

void print_hyperplanes (FILE *f)
{
   int i,j;

   fprintf (f, "\n****** Hyperplanes are: ******\n");
   for (j = 0; j < G_m; j++)
   {  fprintf (f, "  Hyperplane [%i]:  ", j);
      for (i=0; i < G_d; i++) fprintf (f, "%10.3f", G_Hyperplanes [j] [i]);
      fprintf (f, " :%10.3f", G_Hyperplanes [j] [G_d]);
      fprintf (f, "       Contains the vertices no. ");
      for (i = 0; i < G_n; i++)
         if (G_Incidence [i][j])
            fprintf (f, "%i ", i+1);
      fprintf (f, "\n");
   }
}

/****************************************************************************************/

void print_incidence (FILE *f)
{
   int i,j;

   fprintf (f, "\n****** The incidence structure is: ******\n");
   for (j = 0; j < G_m; j++)
   {  fprintf (f, "%i:  ", j);
      for (i = 0; i < G_n; i++)
         if (G_Incidence [i][j])
            fprintf (f, "%i ", i+1);
      fprintf (f, "\n");
   }
}

/****************************************************************************************/

void print_coords (FILE *f, T_Vertex *v)
{  int j;
   fprintf (f, "\n\t(%d)", v -> no);
   for (j=0; j < G_d; j++)
#ifdef RATIONAL
   cout << v -> coords [j];
#else
   fprintf (f, "%10.3f", v -> coords [j]);
#endif
}

/****************************************************************************************/

void print_matrix (FILE *f, rational **A, int m, int n)

{  int i, j;

   for (i = 0; i < m; i++)
   {  fprintf (f, "\n");
      for (j = 0; j < n; j++) fprintf (f, "\t%f", (real) A [i] [j]);
   }
}

/****************************************************************************************/
/****************************************************************************************/

#ifdef STATISTICS

void print_statistics (FILE *f, int method)
   /* prints the partial volume distribution to the screen */

{  int exponent, i, sum_stored = 0, sum_retrieved = 0;

   switch (method)
   {
   case LAWD:
   case LAWND:
      fprintf (f, "\n\nStatistics of the volume computation:");
      fprintf (f, "\nNumber of evaluations of Lawrence's formula:       %i", Stat_Count);
      fprintf (f, "\nSmallest partial volume:    %8.2e", Stat_Smallest);
      fprintf (f, "\nBiggest partial volume:     %8.2e", Stat_Biggest);
      fprintf (f, "\n");

      if (Stat_CountPos [0] > 0)
         fprintf (f, "\nNumber of positive partial volumes smaller than 10^%i: %i",
                  STAT_SMALLEST_EXP, Stat_CountPos [0]);
      fprintf (f, "\nNumber of positive partial volumes between 10^e and 10^(e+1):");
      for (exponent = STAT_SMALLEST_EXP; exponent < STAT_BIGGEST_EXP; exponent ++)
         if (Stat_CountPos [exponent - STAT_SMALLEST_EXP + 1] > 0)
            fprintf (f, "\ne = %4i : %7i", exponent,
                    Stat_CountPos [exponent - STAT_SMALLEST_EXP + 1]);
      if (Stat_CountPos [STAT_BIGGEST_EXP - STAT_SMALLEST_EXP + 2] > 0)
         fprintf (f, "\nNumber of positive partial volumes equal to or bigger than 10^%i: %i", STAT_BIGGEST_EXP + 1,
                  Stat_CountPos [STAT_BIGGEST_EXP - STAT_SMALLEST_EXP + 2]);
      fprintf (f, "\n");

      if (Stat_CountNeg [0] > 0)
         fprintf (f, "\n\nNumber of negative partial volumes smaller than 10^%i: %i",
                  STAT_SMALLEST_EXP, Stat_CountNeg [0]);
      fprintf (f, "\nNumber of negative partial volumes between 10^e and 10^(e+1):");
      for (exponent = STAT_SMALLEST_EXP; exponent < STAT_BIGGEST_EXP; exponent ++)
         if (Stat_CountNeg [exponent - STAT_SMALLEST_EXP + 1] > 0)
            fprintf (f, "\ne = %4i : %7i", exponent,
                    Stat_CountNeg [exponent - STAT_SMALLEST_EXP + 1]);
      if (Stat_CountNeg [STAT_BIGGEST_EXP - STAT_SMALLEST_EXP + 2] > 0)
         fprintf (f, "\nNumber of negative partial volumes equal to or bigger than 10^%i: %i", STAT_BIGGEST_EXP + 1,
                  Stat_CountNeg [STAT_BIGGEST_EXP - STAT_SMALLEST_EXP + 2]);
      break;
   case RCH:
      fprintf (f, "\n\nStatistics of this triangulation:");
      fprintf (f, "\nNumber of simplices:       %i", Stat_Count);
      fprintf (f, "\nSmallest simplex:    %8.2e", Stat_Smallest);
      fprintf (f, "\nBiggest simplex:     %8.2e", Stat_Biggest);
      fprintf (f, "\n");

      if (Stat_CountPos [0] > 0)
         fprintf (f, "\nNumber of simplices smaller than 10^%i: %i", STAT_SMALLEST_EXP,
                  Stat_CountPos [0]);
      fprintf (f, "\nNumber of simplices the volumes of which are between 10^e and 10^(e+1):");
      for (exponent = STAT_SMALLEST_EXP; exponent < STAT_BIGGEST_EXP; exponent ++)
         if (Stat_CountPos [exponent - STAT_SMALLEST_EXP + 1] > 0)
            fprintf (f, "\ne = %4i : %7i", exponent,
                    Stat_CountPos [exponent - STAT_SMALLEST_EXP + 1]);
      if (Stat_CountPos [STAT_BIGGEST_EXP - STAT_SMALLEST_EXP + 2] > 0)
         fprintf (f, "\nNumber of simplices equal to or bigger than 10^%i: %i", STAT_BIGGEST_EXP + 1,
                  Stat_CountPos [STAT_BIGGEST_EXP - STAT_SMALLEST_EXP + 2]);
      break;
   case HOT:
   case RLASS:
      fprintf(f, "\n\nStatistics of storing intermediate results:");
      fprintf (f, "\n\ndimension   volumes stored   volumes retrieved\n\n");
      for (i = 2; i < G_d - 1; i++)
      {
         fprintf(f, "%5i", i);
         fprintf(f, "%17i", Stat_CountStored [i]);
         fprintf(f, "%18i", Stat_CountRetrieved [i]);
         fprintf(f, "\n");
         sum_stored    += Stat_CountStored [i];
         sum_retrieved += Stat_CountRetrieved [i];
      }
      fprintf(f, "\n  all");
      fprintf(f, "%17i", sum_stored);
      fprintf(f, "%18i", sum_retrieved);
      if (method == RLASS)
         fprintf(f, "\n\n%i shifts performed.", Stat_CountShifts);
      break;
   }

   if (method != LRS && method != LAWD)
   {
      fprintf (f, "\n\nTotal amount of memory used on the heap:  %6ld KBytes", Stat_MaxMem / 1024);
   }

}

#endif

/****************************************************************************************/

