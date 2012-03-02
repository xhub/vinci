/****************************************************************************************/
/*                                                                                      */
/*                                     vinci_volume.c                                   */
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
/*                                                                                      */
/* functions computing the volume of a polytope                                         */
/*                                                                                      */
/****************************************************************************************/

#include "vinci.h"

/* global variables for C&H-triangulation and orthonormalisation */

static T_VertexSet *face;
   /* face considered at each recursion level */
static rational ***ortho_basis;
   /* orthonormal basis considered at each recursion level */
static T_VertexSet S;
   /* simplex constructed so far (C&H) */

static T_Key   key;               /* key for storing the actually considered face */
static T_Tree  *tree_volumes;     /* tree for storing intermediate volumes */

/****************************************************************************************/

static void tri (int d, rational *V)
   /* after recursion, contains the d-dimensional volume of face [d] in V               */
   /* The actually considered face is stored in the global variable face [d], the sim-  */
   /* plex constructed so far in the global variable S.                                 */
   /* In this function we work with the ascending order given by the numbers of the     */
   /* vertices.                                                                         */
   /* (Here all volumes must still be divided by dimension!.)                           */

{  int              k;
   T_VertexSuperset *L = create_empty_superset ();
      /* L contains the faces already examined in this recursion step */
   rational         volume;

   if (face [d].lastel > d)
   {  /* cut face [d] with all hyperplanes. If the result is appropriate, start        */
      /* recursion.                                                                    */
      for (k = 0; k < G_m; k++)
      {  /* determine if the smallest element of face [d] is not contained in plane k */
         if (! is_in_hyperplane ((face [d]).loe [0], k))
         {  /* let face [d-1] be the face (face [d] intersected with hyperplane k) */

            if (d == G_d)
               printf ("\nStarting hyperplane %i", k+1);

            intersect_with_hyperplane (face [d], k, &(face [d-1]));

            /* determine whether face [d-1] is a new face and possibly of dimension d-1 */
            if ((face [d-1]).lastel >= d - 1 && !is_in_superset (face [d-1], L))
            {  add_superelement (&L, duplicate_set (face [d-1]));
               add_element (&S, (face [d-1]).loe [0]);

               tri (d - 1, V);

               delete_element (&S, (face [d-1]).loe [0]);
            }
         }
      }
   }
   else /* face [d] is a simplex, which happens at the latest for d = 1 */
   {  /* Insert the elements of face [d] in S to get a simplex of the triangulation. */
      for (k = 1; k <= face [d].lastel; k++)
         add_element (&S, (face [d]).loe [k]);
      /* compute the volume of the simplex defined by S; S may be of empty interior! */
      simplex_volume (S, &volume, FALSE);
      *V += volume;
      for (k = 1; k <= face [d].lastel; k++)
         delete_element (&S, (face [d]).loe [k]);
   }
   free_superset (&L);
}

/****************************************************************************************/

void volume_ch_file (rational *volume, char *vertexfile, char *planesfile)
   /* The function computes the volume of a polytope using the method of Cohen&Hickey.  */
   /* Most of the work is done by the procedure "tri", see there for further            */
   /* information.                                                                      */

{
   rational local_volume = 0;

   read_vertices (vertexfile);
   read_hyperplanes (planesfile);
   compute_incidence ();

   /* renumber vertices such that highly degenerate ones get a lower number */
   renumber_vertices ();

   /* preparing the global variables */
   face = create_faces ();
   copy_set (G_Vertices, &(face [G_d]));
   S = create_empty_set ();
   add_element (&S, G_Vertices.loe [0]);

#ifdef STATISTICS
   init_statistics ();
#endif

   printf ("\nTotal number of hyperplanes: %i\n", G_m);

   tri (G_d, &local_volume);

   /* The real volume is local_volume / (factorial of dimension) */
   (*volume) = local_volume / factorial (G_d);

   free_incidence ();
   free_set (S);
   free_faces (face);
   free_set_and_vertices (G_Vertices);

}

/****************************************************************************************/

static void tri_ortho (int d, rational *V)
   /* After recursion V contains the volume of the actually considered face (stored in  */
   /* the global variable face [d]), and the global variable ortho_basis [d] contains   */
   /* an orthonormal basis of the corresponding linear subspace.                        */
   /* All volumes must still be divided by dimension!                                   */

{  int              i, j, k, dimdiff;
   T_VertexSuperset *L = create_empty_superset ();
      /* L contains the faces already examined in this recursion step */
   rational         volume, *stored_volume;
   rational         distance, maxdistance = 0;
   boolean          i_balance = FALSE, store_volume = FALSE, compute_volume = TRUE;
   T_Key            *dummy;

   *V = 0;

   if (face [d].lastel > d)
   {
      /* if possible try to retrieve the volume */
      dimdiff = G_d - d;
      if ((G_Storage > (dimdiff-2)) && (dimdiff >= 2))
      {
         copy_set (face [d], &(key.vertices.set));
         key.vertices.d = d;
         tree_out (&tree_volumes, &i_balance, key, &stored_volume, &dummy, KEY_VERTICES);
         if (*stored_volume < -0.5)  /* volume has not yet been computed and is -1 */
            /* stored_volume points to a tree element where the volume has to be stored */
            store_volume = TRUE;
         else if (*stored_volume < EPSILON) /* volume has been computed but is 0 */
         {
#ifdef STATISTICS
            Stat_CountRetrieved [d] ++;
#endif
            compute_volume = FALSE;
         }
         else /* volume has been computed and is not 0 */
         {
#ifdef STATISTICS
            Stat_CountRetrieved [d] ++;
#endif
            *V = *stored_volume;
            compute_volume = FALSE;

            /* compute orthonormal basis of face [d] */
            orthonormal (d, face [d], ortho_basis [d]);

         } /* else retrieved volume bigger 0 */
      } /* if G_Storage big enough */

      if (compute_volume) /* do so */
      {
         /* cut face [d] with all hyperplanes. If the result is appropriate, start      */
         /* recursion. */
         for (k = 0; k < G_m; k++)
         {  /* determine furthermore if the smallest element of face [d] is not */
            /* contained in plane k */
            if (! is_in_hyperplane ((face [d]).loe [0], k))
            {  /* let face [d-1] be the face (face [d] intersected with hyperplane k) */

               if (d == G_d)
                  printf ("\nStarting hyperplane %i", k+1);
               intersect_with_hyperplane (face [d], k, &(face [d-1]));

               /* determine whether face [d-1] is a new face and possibly of dimension  */
               /* d-1 */
               if ((face [d-1]).lastel + 1 >= d && !is_in_superset (face [d-1], L))
               {
                  add_superelement (&L, duplicate_set (face [d-1]));

                  tri_ortho (d - 1, &volume);

                  if (fabs (volume) > EPSILON)
                  {
                     /* orthonormalise with the additional point and determine the dis- */
                     /* tance of this point to face [d-1]                               */
                     distance = add_orthonormal (d, face [d-1],
                                                 ortho_basis [d-1], face [d].loe [0]);

                     *V += volume * distance;
                     if (distance > maxdistance)
                     {
                        maxdistance = distance;
                        for (i = 0; i < d; i++)
                           for (j = 0; j < G_d; j++)
                              ortho_basis [d] [i] [j] = ortho_basis [d-1] [i] [j];
                     }
                  }
               }
            }
         } /* for k */
         if (store_volume) /* do so */
         {
               *stored_volume = *V;
#ifdef STATISTICS
               Stat_CountStored [d] ++;
#endif
         }
      }
#ifdef STATISTICS
      Stat_Count ++;
      if (Stat_Count % 100000 == 0)
      {  printf ("\n%10i partial volumes computed.", Stat_Count);
         if (Stat_Count % 1000000 == 0) print_statistics (stdout, HOT);
      }
#endif
   }

   else /* face [d] is a simplex, which happens at the latest for d = 1 */
   {  /* compute the normal basis and the volume      */
      *V = orthonormal (d, face [d], ortho_basis [d]);
   }

   free_superset (&L);
}

/****************************************************************************************/

void volume_ortho_file (rational *volume, char *vertexfile, char *planesfile)
   /* The function computes the volume of a polytope using the face enumeration scheme  */
   /* of Cohen and Hickey and orthonormalisation of Schmidt or Householder. Most of the */
   /* work is done by the procedure "tri_ortho", see there for further information.     */

{
   T_VertexSet vertices = create_empty_set ();
   rational local_volume = 0, scaling_factor;

   read_vertices (vertexfile);
   read_hyperplanes (planesfile);
   compute_incidence ();

   /* renumber vertices such that highly degenerate ones get a lower number */
   renumber_vertices ();

   /* normalise the vertices so that their coefficients are between -1 and 1; this is   */
   /* necessary because otherwise a too big polytope has caused numerical errors        */
   scaling_factor = normalise_vertices ();

   /* preparing the global variables */
   face = create_faces ();
   face [G_d] = duplicate_set (G_Vertices);
   ortho_basis = create_basis ();
   tree_volumes = NULL;

#ifdef STATISTICS
   init_statistics ();
#endif

   key.vertices.set = create_empty_set ();
   key.vertices.d   = G_d;

   printf ("\nTotal number of hyperplanes: %i\n", G_m);

   tri_ortho (G_d, &local_volume);

   /* The real volume of the scaled polytope is local_volume / (factorial of dimension) */
   (*volume) = scaling_factor * local_volume / factorial (G_d);

   free_incidence ();
   free_faces (face);
   free_basis (ortho_basis);
   free_set_and_vertices (vertices);
}

/****************************************************************************************/
/****************************************************************************************/

static void determine_c_and_d (rational *c, rational *d)

{  int i;

   srand (G_RandomSeed); /* sets a seed for the random number generator */

   printf ("\nCoordinates of the objective function c:");
   for (i = 0; i < G_d; i++)
   {  c [i] = (rational) rand ();
      printf ("\nc [%i]: ", i);
#ifdef RATIONAL
      cout << c [i];
#else
      printf ("%6.0f", c [i]);
#endif
   }
   *d = 0;
   printf ("\nd: ");
#ifdef RATIONAL
   cout << *d;
#else
   printf ("%10.0f", *d);
#endif
   printf ("\n");
}

/****************************************************************************************/

static void volume_lawrence_set (T_VertexSet vertices, rational *volume)

{  int      v, i, j, k;
   rational **A = create_matrix (G_d, G_d + 1);
   rational *c = create_vector (), d;
   rational Nv, fv;
   
   (*volume) = 0;
   determine_c_and_d (c, &d);
   
   for (v = 0; v <= vertices.lastel; v++)
   {  /* find the binding constraints for vertex v and write them (transposed!) in A */
      j = k = 0;
      while (k < G_m && j <= G_d)
      {  /* check if v is contained in hyperplane k */
         if (is_in_hyperplane (vertices.loe [v], k))
         {  for (i = 0; i < G_d; i++)
               A [i] [j] = G_Hyperplanes [k] [i];
            j++;
         }
         k++;
      }

      if (j > G_d)
      {  /* v is contained in at least dimension + 1 hyperplanes */
         fprintf (stderr, "\n***** ERROR: Degenerated vertex in 'volume_lawrence':");
         print_coords (stderr, vertices.loe [v]);
         fprintf (stderr, "\n*** Planes containing the vertex:");
         for (i = 0; i < G_m; i++)
            if (is_in_hyperplane (vertices.loe [v], i))
               fprintf (stderr, "\t%i", i + 1);
         exit (0);
      }
      else
      {
         for (i = 0; i < G_d; i++) A [i] [G_d] = c [i];

         Nv = 1 / det_and_invert (A, G_d, G_d + 1, TRUE);
         if (Nv < 0) Nv = - Nv;
         for (i = 0; i < G_d; i++)
           if (fabs (A [i] [G_d]) < EPSILON)
           {  fprintf (stderr, "\n***** ERROR: Division by zero (");
              fprintf (stderr, "%e) in 'volume_lawrence'.", (real) A [i] [G_d]);
              exit (0);
           }
           else Nv /= A [i] [G_d];

         /* compute f(v) */
         fv = d;
         for (i = 0; i < G_d; i++) fv += c [i] * vertices.loe [v] -> coords [i];
      
         for (i = 0; i < G_d; i++) Nv *= fv;
      }
#ifdef STATISTICS
      update_statistics (Nv);
#endif
      (*volume) += Nv;
   }
   (*volume) /= factorial (G_d);

   free_matrix (A, G_d, G_d + 1);
   free_vector (c);
}

/****************************************************************************************/

void volume_lawrence_file (rational *volume, char *vertexfile, char *planesfile)

{
   read_vertices (vertexfile);
   read_hyperplanes (planesfile);
   compute_incidence ();
   volume_lawrence_set (G_Vertices, volume);
   free_hyperplanes ();
   free_incidence ();
   free_set_and_vertices (G_Vertices);
}

/****************************************************************************************/
/****************************************************************************************/

void volume_lawrence_lrs_file (rational *volume, char *planesfile)

{  FILE     *f;
   int      i, j;
   int      *cobasis;
   rational **A, *c, d;
   rational Nv, fv;
   char     line [255], system_command [255], *pos;
   char     *tmp_file, tmp_file_1 [255], tmp_file_2 [255];

   /* Temporary files are stored in the same directory as the executable with the same  */
   /* base name as the planesfile.                                                      */
   tmp_file = strrchr (planesfile, '/');
   if (tmp_file == NULL)
      tmp_file = planesfile;
   else
      tmp_file ++;

   /* copy the planes file to the temporary file planesfile+".li.tmp" */
   sprintf (tmp_file_1, "%s.li.tmp", tmp_file);
   sprintf (tmp_file_2, "%s.lo.tmp", tmp_file);
   sprintf (system_command, "cp %s %s", planesfile, tmp_file_1);
   system (system_command);
   sprintf (system_command, "chmod a+w %s", tmp_file_1);
   system (system_command);
   sprintf (system_command, "echo >> %s", tmp_file_1);
   system (system_command);
   sprintf (system_command, "echo printcobasis >> %s", tmp_file_1);
   system (system_command);

   /* call lrs */
   printf ("\nRunning 'lrs'.");
   sprintf (system_command, "%s < %s > %s", LRS_EXEC, tmp_file_1, tmp_file_2);
   system (system_command);

   /* prepare variables */
   read_hyperplanes (tmp_file_1);
   A = create_matrix (G_d, G_d + 1);
   c = create_vector ();
   cobasis = create_int_vector (G_d);
   (*volume) = 0;
   determine_c_and_d (c, &d);

   /* determine all feasible cobases from the output of lrs */
   f = open_read (tmp_file_2);
   fscanf (f, "***** %*i %*s ");
   /* read lines until EOF */
   while (fgets (line, 255, f) != NULL)
   {
      /* look for string 'facets' in line */
      pos = strstr (line, "facets  ");
      if (pos != NULL)
      {  pos += 8;
         for (i = 0; i < G_d; i++)
         {  sscanf (pos, "%i", & (cobasis [i]));
            cobasis [i] -= 1;
            pos = strpbrk (pos, " ") + 1;
         }
         /* cobasis has been read in */

         /* copy binding restrictions (transposed) and c to A */
         for (i = 0; i < G_d; i++)
         {  for (j = 0; j < G_d; j++) A [i] [j] = G_Hyperplanes [cobasis [j]] [i];
            A [i] [G_d] = c [i];
         }

         /* compute denominator of Nv */
         Nv = 1 / det_and_invert (A, G_d, G_d + 1, TRUE);
         if (Nv < 0) Nv = -Nv;
         for (i = 0; i < G_d; i++)
           if (fabs (A [i] [G_d]) < EPSILON)
           {  fprintf (stderr, "\n***** ERROR: Division by zero (");
              fprintf (stderr, "%e) in 'volume_lawrence_lrs_file'", (real) A [i] [G_d]);
              exit (0);
           }
           else Nv /= A [i] [G_d];

         /* copy binding restrictions and b to A */
         for (i = 0; i < G_d; i++)
           for (j = 0; j <= G_d; j++) A [i] [j] = G_Hyperplanes [cobasis [i]] [j];
         det_and_invert (A, G_d, G_d + 1, TRUE);

         /* compute f(v) */
         fv = d;
         for (i = 0; i < G_d; i++) fv += c [i] * A [i] [G_d];
            /* column "dimension" of A contains vertex coordinates */

         for (i = 0; i < G_d; i++) Nv *= fv;

#ifdef STATISTICS
         update_statistics (Nv);
#endif

         (*volume) += Nv;
      }
   }
   fclose (f);
   (*volume) /= factorial (G_d);

   free_hyperplanes ();
   free_int_vector (cobasis, G_d);
   free_matrix (A, G_d, G_d + 1);
   free_vector (c);
}

/****************************************************************************************/
/****************************************************************************************/

void volume_lrs_file (rational *volume, char *rational_volume, char *vertexfile)
   /* The routine calls lrs on the vertex file to obtain a boundary triangulation. No   */
   /* test is included if the file is of integer or rational type and if the polytope   */
   /* interior contains the origin; thus the returned volume is wrong in these cases.   */
   /* rational_volume contains the volume output of lrs in form of a string.            */

{  FILE    *f1, *f2;
   char    line [255], system_command [255], *pos = NULL;
   char    *tmp_file, tmp_file_1 [255], tmp_file_2 [255];
   boolean volume_found = FALSE;

   /* Temporary files are stored in the same directory as the executable with the same  */
   /* base name as the vertexfile.                                                      */
   tmp_file = strrchr (vertexfile, '/');
   if (tmp_file == NULL)
      tmp_file = vertexfile;
   else
      tmp_file ++;

   /* copy the vertex file to vertexfile+".li.tmp"; at the same time, add needed lines  */
   sprintf (tmp_file_1, "%s.li.tmp", tmp_file);
   sprintf (tmp_file_2, "%s.lo.tmp", tmp_file);
   f1 = open_read (vertexfile);
   f2 = fopen (tmp_file_1, "w");
   fprintf (f2, "*\nV-representation\nbegin\n");
   while (fgets (line, 255, f1) != NULL && strncmp (line, "end", 3))
      fprintf (f2, "%s", line);
   fprintf (f2, "end\nvolume\n");
   fclose (f1);
   fclose (f2);

   /* call lrs */
   printf ("\nRunning 'lrs'.");
   sprintf (system_command, "%s < %s > %s", LRS_EXEC, tmp_file_1, tmp_file_2);
   system (system_command);

   /* search for the volume in the output */
   f1 = open_read (tmp_file_2);
   while (fgets (line, 255, f1) != NULL)
   {  pos = strstr (line, "*Volume=");
      if (pos != NULL)
      {  volume_found = TRUE;
         break;
      }
   }
   fclose (f1);

   if (volume_found)
   {  pos += 8; /* pos now points to the volume */
      strcpy (rational_volume, pos);
      sread_rational_value (pos, volume);
   }
   else
   {  *volume = 0;
      sprintf (rational_volume, "0\n");
   }
}

/****************************************************************************************/
/****************************************************************************************/
