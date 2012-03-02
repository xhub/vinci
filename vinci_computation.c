/****************************************************************************************/
/*                                                                                      */
/*                                  vinci_computation.c                                 */
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
/* Last Changes: February 2, 2001                                                       */
/*                                                                                      */
/****************************************************************************************/
/*                                                                                      */
/* the main computational routines of the package                                       */
/*                                                                                      */
/****************************************************************************************/

#include "vinci.h"

/****************************************************************************************/

rational factorial (int n)
   /* calculates n! by calculating all i! for 0 <= i <= G_d and storing them in an      */
   /* array */

{  static rational *fact;
   static boolean  first_call = TRUE;
   int i;

   if (first_call)
   {
      fact = create_fact ();
      fact [0] = 1;
      for (i = 1; i <= G_d; i++) fact [i] = i * fact [i - 1];
      first_call = FALSE;
   }

   return (fact [n]);
}

/****************************************************************************************/
/*                computations with simplices and determinants                          */
/****************************************************************************************/

rational det_and_invert (rational **A, int rows, int columns, boolean verbose)
   /* Let B be the quadratic submatrix of A defined by the first "rows" rows and        */
   /* columns of A. Then this function turns B into an upper triangular matrix by means */
   /* of Gaussian pivot operations. At the same time, the colums "rows + 1" up to       */
   /* "columns" of A are replaced by B^-1 times these columns. As a byproduct, the de-  */
   /* terminant of B (modulo its sign) is returned.                                     */
   /* If verbose is TRUE an error message is output for zero volume.                    */

{  rational dummy, det = 1, *row_pointer;
   int      i, j, k, r;
#if PIVOTING == 2
   int      s;
#elif PIVOTING == 0
   boolean  pivot_found = FALSE;
#endif

   for (k = 0; k < rows; k++)
   {  /* look for a pivot element */
      r = k;

#if PIVOTING == 0
      /* choose the first row r with r >= k and A (r, k) > MIN_PIVOT */
      for (i = k + 1; i < rows && !pivot_found; i++)
      {  dummy = fabs (A [i] [k]);
         if (dummy > MIN_PIVOT)
         {  r = i;
            pivot_found = TRUE;
         }
         else if (dummy > fabs (A [r] [k])) r = i;
      }
#elif PIVOTING == 2
      /* total pivoting; choose row r and column s with maximal |A (r, s)|, but only if */
      /* columns is equal to rows (no equations have to be solved); otherwise use par-  */
      /* tial pivoting                                                                  */
      if (columns == rows)
         for (i = k; i < G_d; i++)
            for (j = k; j < G_d; j++)
               if (fabs (A [i] [j]) > fabs (A [r] [s]))
               {  r = i;
                  s = j;}
      else
         for (i = k + 1; i < rows; i++)
            if (fabs (A [i] [k]) > fabs (A [r] [k])) r = i;
#else
      /* partial pivoting; choose row r with maximal |A (r, k)| */
      for (i = k + 1; i < rows; i++)
         if (fabs (A [i] [k]) > fabs (A [r] [k])) r = i;
#endif

      if (r != k) /* exchange the rows r and k */
      {  row_pointer = A [k];
         A [k] = A [r];
         A [r] = row_pointer;
      }

      /* check for singularity */
      if (fabs (A [k] [k]) < EPSILON)
      {  if (verbose)
         {  fprintf (stderr, "\n***** ERROR: (Almost) singular matrix in 'det_and_");
            fprintf (stderr, "invert';\n      One pivot is %20.18e.\n", (real) A [k] [k]);
         }
         return 0;
      }

      /* now make an elimination step */
      for (i = k + 1; i < rows; i++)
      {  dummy = A [i] [k] / A [k] [k]; /* store the factor once for the row */
         for (j = k + 1; j < columns; j++) A [i] [j] -= dummy * A [k] [j];
      };
      det *= A [k] [k];
   }

   /* now B is an upper triangular matrix; calculate B^-1 times the last columns */
   for (k = rows; k < columns; k++)
      for (i = rows - 1; i >= 0; i--)
      {  /* compute (B^-1 * A_k)_i */
         for (j = i + 1; j < rows; j++) A [i] [k] -= A [i] [j] * A [j] [k];
         A [i] [k] /= A [i] [i];
      }

   return det;
}

/****************************************************************************************/

void simplex_volume (T_VertexSet S, rational *volume, boolean verbose)
   /* From a set of vertices which are presumed of number G_d + 1 and affinely indepen- */
   /* dent, the function computes G_d! times the volume of the corresponding simplex    */
   /* and stores it in volume.                                                          */
   /* If verbose is TRUE an error message is output for zero volume.                    */

{  static rational **A;
   static boolean  first_call = TRUE;
   rational dummy;
   int      i, j;

   if (first_call)
   {
      A = create_matrix (G_d, G_d);
      first_call = FALSE;
   }

   /* copy the relevant information into A */
   for (j = 0; j < G_d; j++)
   {  dummy = S.loe [G_d] -> coords [j];
      for (i = 0; i < G_d; i++)
         A [i] [j] = (S.loe [i] -> coords [j]) - dummy;
   }

   *volume = det_and_invert (A, G_d, G_d, verbose);

   if (*volume < 0) *volume = - (*volume);

#ifdef STATISTICS
   if (*volume > EPSILON) update_statistics ((real) *volume / factorial (G_d));
#endif

}

/****************************************************************************************/
/*                      routines for computing orthonormal bases                        */
/****************************************************************************************/

rational add_orthonormal (int d, T_VertexSet face, rational **H, T_Vertex *vertex)
   /* The function assumes Householder vectors of the linear subspace associated with   */
   /* the variable face in the first d-1 rows of H. It computes the Householder vector  */
   /* for the new vertex minus a vertex of face and stores it into the dth row of H.    */
   /* The first k-1 zero entries of the Householder vector k are neglected, its essen-  */
   /* tial part is stored in the columns k to dimension. All Householder vectors are    */
   /* normalised. This may be changed, but then do not forget to divide by its norm in  */
   /* later steps!                                                                      */
   /* Return value is the distance of the vertex to face.                               */

{  int      i, j;
   rational scalar_product, alpha_squared, alpha, divisor, distance, mindistance = 1e99;

   /* look for the vertex in face with the smallest distance to the given vertex  */
   for (i = 0; i <= face.lastel; i++)
   {
      distance = 0;
      for (j = 0; j < G_d; j++)
         distance +=   (vertex -> coords [j] - face.loe [i] -> coords [j])
                     * (vertex -> coords [j] - face.loe [i] -> coords [j]);
      if (distance < mindistance)
      {  mindistance = distance;
         for (j = 0; j < G_d; j++)
            H [d-1] [j] = vertex -> coords [j] - face.loe [i] -> coords [j];
      }
   }

   /* multiply from the left with the previous Householder matrices */
   for (i = 0; i < d-1; i++)
   {  scalar_product = 0;
      for (j = i; j < G_d; j++)
         scalar_product += H [i] [j] * H [d-1] [j];
      scalar_product *= 2;
      for (j = i; j < G_d; j++)
         H [d-1] [j] -= scalar_product * H [i] [j];
   }

   /* compute the new Householder vector and the distance */
   alpha_squared = 0;
   for (j = d-1; j < G_d; j++)
      alpha_squared += H [d-1] [j] * H [d-1] [j];
   alpha = sqrt (alpha_squared);
   if (H [d-1] [d-1] < 0)
      alpha = - alpha;

   divisor = sqrt (2 * (alpha_squared + alpha * H [d-1] [d-1]));
   H [d-1] [d-1] += alpha;
   for (j = d-1; j < G_d; j++)
      H [d-1] [j] /= divisor;

   return fabs (alpha);

}

/****************************************************************************************/

rational orthonormal (int d, T_VertexSet face, rational **H)
   /* The function computes Householder vectors for the linear subspace associated with */
   /* the variable face. The vectors are returned via the first d rows of H.            */
   /* This routine uses pivoting techniques and should be numerically stable.           */
   /* The return value is d! times the volume of the orthonormalised simplex.           */

{  static rational **local_H;
   static int      m = 0;
      /* the number of rows in local_H */
   rational *dummy_row;
   int      i, j, k, maxindex = 0;
   rational scalar_product, alpha_squared = -1, alpha, divisor;
   rational volume = 1;

   /* create local_H in the correct dimension */
   if (m == 0)
   {
      m = face.lastel;
      local_H = create_matrix (m, G_d + 1);
         /* The last component of each row will contain the norm of the essential part  */
         /* of the vector  */
   }
   else if (face.lastel > m)
   {
      redim_matrix (&local_H, m, face.lastel, G_d + 1);
      m = face.lastel;
   }

   /* copy the spanning vectors into local_H; the zeroth vertex of face [d] is trans-   */
   /* lated into the origin                                                             */
   for (i = 0; i < face.lastel; i++)
   {  local_H [i] [G_d] = 0;
      for (j = 0; j < G_d; j++)
      {  local_H [i] [j] =   face.loe [i+1] -> coords [j]
                           - face.loe [0  ] -> coords [j];
         local_H [i] [G_d] += local_H [i] [j] * local_H [i] [j];
      }
      if (local_H [i] [G_d] > alpha_squared)
      {  alpha_squared = local_H [i] [G_d];
         maxindex = i;
      }
   }

   for (k = 0; k < d; k++)
   {  /* compute k-th Householder vector; choose the element whose relevant part has    */
      /* the biggest norm as next element                                               */

      /* If alpha_squared is 0 matrix has not got the full rank. Stop orthonormalising. */
      if (alpha_squared / volume < EPSILON)
      {  return 0;
      }

      /* otherwise exchange rows */
      dummy_row = local_H [k];
      local_H [k] = local_H [maxindex];
      local_H [maxindex] = dummy_row;

      /* compute the new Householder vector */
      alpha = sqrt (alpha_squared);
      volume *= alpha;
      if (local_H [k] [k] < 0)
         alpha = - alpha;
      divisor = sqrt (2 * (alpha_squared + alpha * local_H [k] [k]));
      local_H [k] [k] += alpha;
      for (j = k; j < G_d; j++)
         local_H [k] [j] /= divisor;

      /* apply the Householder matrix to the resting rows of local_H; at the same time  */
      /* compute the norms of the relevant parts                                        */
      alpha_squared = -1;
      for (i = k + 1; i < face.lastel; i++)
      {
         scalar_product = 0;
         for (j = k; j < G_d; j++)
            scalar_product += local_H [k] [j] * local_H [i] [j];
         scalar_product *= 2;
         for (j = k; j < G_d; j++)
            local_H [i] [j] -= scalar_product * local_H [k] [j];

         local_H [i] [G_d] -= local_H [i] [k] * local_H [i] [k];
         if (local_H [i] [G_d] > alpha_squared)
         {  alpha_squared = local_H [i] [G_d];
            maxindex = i;
         }
      }

   } /* for */

   /* copy the results into H */
   for (i = 0; i < d; i++)
      for (j = 0; j < G_d; j++)
         H [i] [j] = local_H [i] [j];

   return volume;
}

/****************************************************************************************/
