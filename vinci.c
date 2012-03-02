/****************************************************************************************/
/*                                                                                      */
/*                                      vinci                                           */
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
/* Last Changes: February 9, 2001                                                       */
/*                                                                                      */
/****************************************************************************************/

#include "vinci.h"

/****************************************************************************************/

boolean existing_files (char *basename, int *ext, int  *ine, int *d, int *m, int *n)

   /* Checks if the .ext- and .icd-files exist for the specified filename, in which     */
   /* case their data types are stored in the variables ext and ine. The dimension, the */
   /* number of inequalities and the number of vertices are stored in the variables d,  */
   /* m and n, respectively. If an inconsistency in the dimensions is detected the      */
   /* return value is FALSE.                                                            */

{  FILE    *fp;
   char    filename [255], line [255], data_type [255];
   boolean begin_found, return_value = TRUE;
   int     local_d;

   *ext = *ine = NONE;
   *d = *m = *n = 0;

   /* try to open the .ext-file */
   sprintf (filename, "%s.ext", basename);
   if ((fp = fopen (filename, "r")))
   {  begin_found = FALSE;
      while (!feof (fp) && !begin_found)
      {  fgets (line, 255, fp);
         if (!strncmp (line, "begin", 5)) begin_found = TRUE;
      }
      if (begin_found)
      {  fscanf (fp, "%i %i %s ", n, d, data_type);
         *ext = determine_data_type (data_type);
      }
      fclose (fp);
   }

   /* try to open the .ine-file */
   sprintf (filename, "%s.ine", basename);
   if ((fp = fopen (filename, "r")))
   {  begin_found = FALSE;
      while (!feof (fp) && !begin_found)
      {  fgets (line, 255, fp);
         if (!strncmp (line, "begin", 5)) begin_found = TRUE;
      }
      if (begin_found)
      {  fscanf (fp, "%i %i %s ", m, &local_d, data_type);
         *ine = determine_data_type (data_type);
         if ((*d > 0) && (*d != local_d))
         {  *d = 0;
            printf ("\nThe information about the dimension is contradictory in your .ext- and .ine-");
            printf ("\nfiles. Please check the files and restart the programme.\n");
            return_value = FALSE;
         }
      }
      fclose (fp);
   }

   return return_value;

}

/****************************************************************************************/

void existing_programmes (boolean *lrs)
   /* determines if the specified executables exist */

{  FILE *fp;
   char filename [255];

   *lrs = FALSE;

   /* search "lrs" */
   sprintf (filename, "%s", LRS_EXEC);
   if ((fp = fopen (filename, "rb")))
   {  *lrs = TRUE;
      fclose (fp);
   }

}

/****************************************************************************************/

int method_proposal (int ext, int ine, boolean lrs, int d, int m, int n)
   /* The function tries to determine the optimal method depending on the existing      */
   /* files and their data types, the installed programmes and the values of d, m and   */
   /* n; the corresponding method code is returned (or NONE if no method can be used).  */
   /* The parameters are the results of 'existing_files' and 'existing_programmes'.     */

{  int files_count = 0; /* number of files for the polytope */

   if (ext != NONE) files_count++;
   if (ine != NONE) files_count++;

   printf ("\nThe volume computation methods which can be applied depend on the types of");
   printf ("\nthe existing files as well as on the additional codes you have installed.");
   printf ("\n\nFor the polytope you specified ");

   switch (files_count)
   {
   case 0:
      printf ("no file exists.");
      break;
   case 1:
      printf ("only the .");
      if   (ext != NONE) printf ("ext");
      else               printf ("ine");
      printf ("-file exists.");
      break;
   case 2:
      printf ("the .ext- and .ine-files exist.");
      break;
   }

   printf ("\n(Please note that files which do not contain a line 'begin' or which use a");
   printf ("\ndata type other than 'integer', 'real' or 'rational' are considered as");
   printf ("\nnon-existing.)");

   if (lrs)
      printf ("\nYou installed lrs.");
   else
      printf ("\nNo additional programme is installed.");

   printf ("\n");


   if (files_count == 0)
   {  printf ("\nNo method can be applied; please make sure that at least one of the .ext- and");
      printf ("\n.ine-files is present.\n");
      return NONE;
   }

      if (files_count == 1 && ext != NONE && !(lrs && ext != REAL_T))
   {  printf ("\nActually no method can be applied; please install lrs or create the .ine");
      printf ("\nfile.\n");
      return NONE;
   }

   printf ("\nThus the following methods can be applied (for a description of each method");
   printf ("\nplease read the manual):");
   if (files_count == 1)
      if (ext == NONE)
      {
         printf ("\n- rlass");
         if (lrs)
            printf ("\n- lawd");
      }
      else
      {
         if (lrs && ext != REAL_T)
            printf ("\n- lrs, if the origin is in the interior of the polytope");
      }
   else
   {  printf ("\n- hot");
      printf ("\n- rlass");
      printf ("\n- rch");
      if (lrs)
         printf ("\n- lawd");
      if (lrs && ext != REAL_T)
         printf ("\n- lrs, if the origin is in the interior of the polytope");
   }

   printf ("\n");

   /* determine the method; if possible choose hot or rlass, otherwise lrs */
   if (ext != NONE && ine != NONE)
   {  printf ("\nI recommend 'hot'.\n");
      return HOT;
   }
   else if (ine != NONE)
   {  printf ("\nI recommend 'rlass'.\n");
      return RLASS;
   }
   else /* only ext present */
   {  printf ("\nThe only possible method is 'lrs'; note that it will fail if the origin does");
      printf ("\nnot lie in the polytope interior!\n");
      return LRS;
   }

}

/****************************************************************************************/

boolean method_test (int method, int ext, int ine, boolean lrs, int d, int m, int n)
   /* The function tests if the chosen method can be run in the situation characterised */
   /* by the parameters which are the same as in 'method_proposal'.                     */

{  int     files_count = 0;
   boolean ok = TRUE;

   if (ext != NONE) files_count++;
   if (ine != NONE) files_count++;

   switch (method)
   {
   case LRS:
      if (lrs && ext != NONE)
         if (ext != REAL_T)
         {  printf ("\nWARNING: The method 'lrs' only works if the polytope interior contains the");
            printf ("\norigin; otherwise the computed volume will be wrong!");
         }
         else
         {  printf ("\n'lrs' will only compute the correct volume if the .ext-file uses integer");
            printf ("\nor rational coordinates.");
            ok = FALSE;
         }
      else
      {  printf ("\nTo use the method 'lrs' you need to install lrs and to create the .ext-file.");
         if (!lrs && ext == NONE)
            printf ("Do so and restart the programme.");
         else if (ext == NONE)
            printf ("Please create the vertex file and restart the programme.");
         else
            printf ("\nPlease install 'lrs' and run the programme again.");
         ok = FALSE;
      }
      break;
   case RCH:
      if (ext == NONE || ine == NONE)
      {  printf ("\nTo use 'rch' you need the vertex and the hyperplane file. However, ");
         if (ext == NONE && ine == NONE)
         {  printf ("none of them");
            printf ("\ncould be found. Please create them and restart the programme.");
         }
         else if (ext == NONE)
         {  printf ("the");
            printf ("\n.ext-file could not be found. Please create it and rerun the programme.");
         }
         else
         {  printf ("the");
            printf ("\n.ine-file could not be found. Please create it and rerun the programme.");
         }
         ok = FALSE;
      }
      break;
   case HOT:
      if (ext == NONE || ine == NONE)
      {  printf ("\nTo use 'hot' you need the vertex and the hyperplane file. However, ");
         if (ext == NONE && ine == NONE)
         {  printf ("none of them could be found. Please create them and restart the");
            printf ("\nprogramme.");
         }
         else if (ext == NONE)
         {  printf ("the");
            printf ("\n.ext-file could not be found. Please create it and rerun the programme.");
         }
         else
         {  printf ("the");
            printf ("\n.ine-file could not be found. Please create it and rerun the programme.");
         }
         ok = FALSE;
      }
      break;
   case LAWND:
      if (files_count == 2)
         printf ("\nWARNING: Be aware that 'lawnd' only works for non degenerate input data!");
      else
      {  printf ("\nTo use 'lawnd' you need the vertex and the hyperplane file. However, ");
         if (files_count == 0)
            printf ("\nnone of them could be found. Please create them and restart the programme.");
         else
         {
            if (ext == NONE)
            {  printf ("\nthe .ext-file could not be found. Please create it and rerun");
               printf ("\nthe programme.");
            }
            else
            {  printf ("\nthe .ine-file could not be found. Please create it and rerun");
               printf ("\nthe programme.");
            }
            ok = FALSE;
         }
      }
      break;
   case LAWD:
      if (lrs && ine != NONE)
      {  if (ine == REAL_T)
         {  printf ("\n'lawd' can only be used for integer or rational data, not for floating");
            printf ("\npoint data.");
            ok = FALSE;
         }
      }
      else
      {  printf ("\nTo use the method 'lawd' you need to install lrs and to create the .ine-file.");
         if (!lrs && ine == NONE)
            printf ("\nDo so and restart the programme.");
         else if (ine == NONE)
            printf ("\nPlease create the inequality file and restart the programme.");
         else
            printf ("\nPlease install 'lrs' and run the programme again.");
         ok = FALSE;
      }
      break;
   case RLASS:
      if (ine == NONE)
      {  printf ("\nYou need the inequality file to run 'rlass'. Please make it available and");
         printf ("\nstart the programme once again.");
         ok = FALSE;
      }
      break;
   }
   return ok;
}

/****************************************************************************************/

void print_methods (FILE *f)
   /* prints the existing volume computation methods to the specified file              */

{
   fprintf (f, "\nTo specify a volume computation method use the option '-m' followed by a");
   fprintf (f, "\nblank and one of the following codes:");
   fprintf (f, "\n\n- hot   for the hybrid orthonormalisation technique");
   fprintf (f, "\n- rlass for Lasserre's revised recursive scheme");
   fprintf (f, "\n- rch   for revised Cohen-Hickey-Triangulation");
   fprintf (f, "\n- lawd  for Lawrence's formula in the general case using 'lrs'");
   fprintf (f, "\n- lawnd for Lawrence's formula in the non degenerate case");
   fprintf (f, "\n- lrs   for boundary triangulation via 'lrs'");
   fprintf (f, "\n\nIf no method is specified, the programme tries to find an optimal algorithm.");
}

/****************************************************************************************/

void print_help (FILE *f)
   /* prints a short help text to the specified file */

{
   fprintf (f, "\n%s", T01);
   fprintf (f, "\n%s", T02);
   fprintf (f, "\n%s", T03);
   fprintf (f, "\n%s", T04);
   fprintf (f, "\n%s", T05);
   fprintf (f, "\n%s", T06);
   fprintf (f, "\n%s", T07);
   fprintf (f, "\n%s", T08);
   fprintf (f, "\n%s", T09);
   fprintf (f, "\n%s", T10);
   fprintf (f, "\n%s", T11);
   fprintf (f, "\n%s", T12);
   fprintf (f, "\n%s", T13);
   fprintf (f, "\n%s", T14);
   fprintf (f, "\n%s", T15);
   fprintf (f, "\n%s", T16);
   fprintf (f, "\n%s", T17);
   fprintf (f, "\n%s", T18);
   fprintf (f, "\n%s", T19);
   fprintf (f, "\n%s", T20);
}

/****************************************************************************************/

void print_pivoting (FILE *f, int method)
   /* prints the pivoting strategy to the specified file */

{
   if (method != RLASS)
   {
      if (PIVOTING == 1)
         fprintf (f, "\nThe pivoting strategy is partial pivoting.");
      else if (PIVOTING == 2)
      {
         fprintf (f, "\nThe pivoting strategy is global pivoting where possible or partial pivoting");
         fprintf (f, "\notherwise.");
      }
      else
         fprintf (f, "\nThe pivoting strategy uses a MIN_PIVOT of %g.", MIN_PIVOT);
   }
   else
   {
      if (PIVOTING_LASS == 1)
         fprintf (f, "\nThe pivoting strategy is partial pivoting.");
      else if (PIVOTING_LASS == 2)
      {
         fprintf (f, "\nThe pivoting strategy is global pivoting where possible or partial pivoting");
         fprintf (f, "\notherwise.");
      }
      else
         fprintf (f, "\nThe pivoting strategy uses a MIN_PIVOT of %g.",
         MIN_PIVOT_LASS);
   }
}

/****************************************************************************************/

boolean determine_method (char *choice, int *method)
   /* tries to determine the desired volume computation method from the contents of     */
   /* choice. If this is successfully done, TRUE is returned and FALSE otherwise.       */

{  boolean ok = TRUE;

   if      (!strcmp (choice, "lrs"))
      *method = LRS;
   else if (!strcmp (choice, "rch"))
      *method = RCH;
   else if (!strcmp (choice, "hot"))
      *method = HOT;
   else if (!strcmp (choice, "lawd"))
      *method = LAWD;
   else if (!strcmp (choice, "lawnd"))
      *method = LAWND;
   else if (!strcmp (choice, "rlass"))
      *method = RLASS;
   else
   {  *method = 0;
      ok = FALSE;
      printf ("\nThe desired volume computation method does not exist.");
      print_methods (stdout);
   }
   return ok;
}

/****************************************************************************************/

boolean evaluate_parameters (int argc, char *argv [], char *filename, int *method)
   /* The function tries to determine the parameter values and consequently sets the    */
   /* filename, the desired method and the global variable "precomp". If an error oc-   */
   /* curs which forces the programme to break the return value is FALSE.               */

{  boolean filename_chosen = FALSE, ok = TRUE;
   int     index = 1;
      /* points to the actually considered entry of the parameter list */

   *method = NONE;
   G_Storage = -1;

   while (index < argc && ok)
   {  /* analyse entry "index" of argv */

      if (argv [index] [0] != '-')
      {  /* this must be the file name */
         if (filename_chosen)
         {  printf ("\nThere seem to be two file names in your command line, '%s'", filename);
            printf ("\nand '%s'. Please check this again.", argv [index]);
            ok = FALSE;
         }
         else
         {  filename_chosen = TRUE;
            strcpy (filename, argv [index]);
            index++;
         }
      }

      else /* entry is an option, check for type */
      if (!strcmp (argv [index], "-m"))
      {  if (*method != NONE)
         {  printf ("\nYou specified the option '-m' twice. Please decide for one of them.");
            ok = FALSE;
         }
         else if (index + 1 >= argc)
         {  printf ("\nYou used the option '-m' without anything following.");
            print_methods (stdout);
            ok = FALSE;
         }
         else
         {  ok = determine_method (argv [index + 1], method);
            index += 2;
         }
      }

      else if (strlen (argv [index]) >= 2 && argv [index] [1] == 's')
      {
         if (strlen (argv [index]) == 2)
         {  printf ("\nYou specified the option '-s' without any integer following; use '-s0' to");
            printf ("\nprevent storing of intermediate results or any higher value, e. g. '-s5', for");
            printf ("\ndetermining the level down to which storage is performed.");
            ok = FALSE;
         }
         else if (G_Storage != -1)
         {  printf ("\nYou specified both the options '-s%d' and '%s'; please decide for one of them.",
                    G_Storage, argv [index]);
            ok = FALSE;
         }
         else
         {  G_Storage = atoi (argv [index] + 2);
            index++;
         }
      }

      else if (strlen (argv [index]) >= 2 && argv [index] [1] == 'r')
      {
         if (strlen (argv [index]) == 2)
         {  printf ("\nYou specified the option '-r' without any integer following; use '-r' directly");
            printf ("\nfollowed by an integer to set the random seed to this value, e. g. '-r10'.");
            ok = FALSE;
         }
         else
         {  G_RandomSeed = atoi (argv [index] + 2);
            index++;
         }
      }

      else
      {  printf ("\nYou specified the option '%s' which does not exist. The following text provides", argv [index]);
         printf ("\nsome help on how to use the programme.\n");
         print_help (stdout);
         ok = FALSE;
      }

   }

   if (ok)
   {  if (G_Storage == -1)
         G_Storage = DEFAULT_STORAGE;

      if (!filename_chosen)
      {  printf ("\nYou did not specify any file name, so I suppose that you are not familiar with");
         printf ("\nthe programme. Please read the short help text below, and if you still have");
         printf ("\nopen questions, consult the manual.\n");
         print_help (stdout);
         ok = FALSE;
      }
   }

   return ok;

}

/****************************************************************************************/

int main (int argc, char *argv [])

{  char       filename [255];
   int        method = NONE;
   boolean    ok;
   char       vertexfile [255], planesfile [255];
   rational   volume;
   char       rational_volume [255];
   int        ext, ine;
   boolean    lrs;
   int        d, m, n;
   int        i;
   struct tms time_info;


   printf ("\n                       VINCI - Version %s as of %s\n", VERSION, VERSION_DATE);
   printf ("\n%s", COPYRIGHT1);
   printf ("\n%s", COPYRIGHT2);
   printf ("\n%s\n", COPYRIGHT3);

   ok = evaluate_parameters (argc, argv, filename, &method);

   if (!ok)
      printf ("\n\n");
   else
   {
      sprintf (vertexfile, "%s.ext", filename);
      sprintf (planesfile, "%s.ine", filename);

      ok = existing_files (filename, &ext, &ine, &d, &m, &n);
      existing_programmes (&lrs);

      if (ok)
      {
         if (method == NONE)
         {  printf ("\nYou did not specify any method; let us analyse the situation.\n");
            method = method_proposal (ext, ine, lrs, d, m, n);
         }
         else
         {  ok = method_test (method, ext, ine, lrs, d, m, n);
            if (!ok)
            {  method = NONE;
               printf ("\nAlternatively let us analyse what can be done in the present situation.\n");
               method_proposal (ext, ine, lrs, d, m, n);
            }
         }

         if (method == NONE)
            printf ("\n\n");
         else
         {

            printf ("\n_______________________________________________________________________________\n");

            switch (method)
            {
            case RCH:
               printf ("\nUsing revised Cohen-Hickey-triangulation for computing the volume.");
               print_pivoting (stdout, method);
               printf ("\n");
               volume_ch_file (&volume, vertexfile, planesfile);
               break;
            case HOT:
               printf ("\nUsing the hybrid orthonormalisation technique.");
               printf ("\nThe storage level is set to %i.", G_Storage);
               print_pivoting (stdout, method);
               printf ("\n");
               volume_ortho_file (&volume, vertexfile, planesfile);
               break;
            case LAWND:
               printf ("\nUsing Lawrence's formula in the non-degenerate case for computing ");
               printf ("\nthe volume.");
               printf ("\nThe random seed is set to %i.", G_RandomSeed);
               print_pivoting (stdout, method);
               printf ("\n");
               volume_lawrence_file (&volume, vertexfile, planesfile);
               break;
            case LAWD:
               printf ("\nUsing 'lrs' and Lawrence's formula for computing the volume.");
               printf ("\nThe random seed is set to %i.", G_RandomSeed);
               print_pivoting (stdout, method);
               printf ("\n");
               volume_lawrence_lrs_file (&volume, planesfile);
               break;
            case RLASS:
               printf ("\nUsing Lasserre's revised recursive scheme for computing the volume");
               printf ("\nThe storage level is set to %i.", G_Storage);
               print_pivoting (stdout, method);
               printf ("\n");
               volume_lasserre_file (&volume, planesfile);
               break;
            case LRS:
               printf ("\nUsing 'lrs' for computing a boundary triangulation.\n");
               volume_lrs_file (&volume, rational_volume, vertexfile);
               break;
            }

            times (&time_info);

            printf ("\n_______________________________________________________________________________");
#ifdef STATISTICS
            print_statistics (stdout, method);
            printf ("\n_______________________________________________________________________________");
#endif

#ifdef RATIONAL
            printf ("\n\nVolume: ");
            cout << volume;
#else
            if (method != LRS)
               printf ("\n\nVolume: %20.12e", volume);
            else
            {  printf ("\n\nVolume: ");
               for (i = 0; i < strlen (rational_volume) - 1; i++)
                  printf ("%c", rational_volume [i]);
               if (strpbrk (rational_volume, "/") != NULL)
                  printf (" = %20.12e", volume);
            }

#endif

            printf ("\n\nTime passed with computation: %8.1f s\n",
                    ((double) time_info.tms_utime) / ((double) CLOCKS_PER_SEC));
            printf ("_______________________________________________________________________________\n\n");

         }
      }
   }

   return TRUE;
}

/****************************************************************************************/
