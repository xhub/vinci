/****************************************************************************************/
/*                                                                                      */
/*                                   vinci_set.c                                        */
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
/*  functions on the data structures "set of pointers to vertices" and sets of them     */
/*                                                                                      */
/****************************************************************************************/

#include "vinci.h"

/****************************************************************************************/
/*                           functions on sets of vertices                              */
/****************************************************************************************/

static void redimloe (T_VertexSet *s)
   /* adds some space to loe */

{
   s -> maxel += ARRAYSIZESTEP;

   s -> loe = (T_Vertex **)
              my_realloc (s -> loe, (s -> maxel + 1) * sizeof (T_Vertex *),
                          ARRAYSIZESTEP * sizeof (T_Vertex *));
}

/****************************************************************************************/

static int position_of_element (T_Vertex *e, T_VertexSet s)
   /* If e is in s, the function returns the index i with s.loe [i] = e. Otherwise, -1  */
   /* is returned. The function uses a binary search algorithm. */

{  int first = 0, last = s.lastel; /* e is searched between loe [first] and loe [last]. */
   int middle, e_no = e -> no, middle_no;
   
   if (s.lastel == -1) return -1;
   while (last - first > 1)
   {  middle = (first + last) / 2;
      middle_no = s.loe [middle] -> no;
      if (middle_no == e_no)
         return middle;
      else if (middle_no > e_no)
         last = middle;
      else
         first = middle;
   }
   if (s.loe [first] -> no == e_no)
      return first;
   else if (s.loe [last] -> no == e_no)
      return last;
   else
      return -1;
}
 
/****************************************************************************************/
/****************************************************************************************/

T_VertexSet create_empty_set ()
   /* returns an empty set */
   
{  T_VertexSet s;
   
   s.maxel = G_d + ARRAYSIZESTEP;
   s.loe = NULL;
   s.lastel = -1;
   s.loe = (T_Vertex **) my_malloc ((s.maxel + 1) * sizeof (T_Vertex *));
   return s;
}

/****************************************************************************************/

void free_set (T_VertexSet s)
   /* frees the memory space used by s. Attention: You will not be able to use s after! */
   
{  
   my_free (s.loe, (s.maxel + 1) * sizeof (T_Vertex *));
}

/****************************************************************************************/

void free_set_and_vertices (T_VertexSet s)
   /* frees the memory space used by s and all its elements. Attention: You will nei-   */
   /* ther be able to use s nor any of its elements after! */
   
{  int i;
   
   for (i = 0; i <= s.lastel; i++) free_vertex (s.loe [i]);
   free_set (s);
}

/****************************************************************************************/

void clear_set (T_VertexSet *s)
   /* empty the set s */
   
{  (*s).lastel = -1;
}

/****************************************************************************************/

T_VertexSet duplicate_set (T_VertexSet s)
   /* creates a new set with the same elements as s */
   
{  T_VertexSet newset;

   newset.maxel = s.lastel;
   newset.lastel = s.lastel;
   newset.loe = (T_Vertex **) my_malloc ((s.lastel + 1) * sizeof (T_Vertex *));
   memcpy (newset.loe, s.loe, (s.lastel + 1) * sizeof (T_Vertex *));
   return newset;
}
      
/****************************************************************************************/

void copy_set (T_VertexSet s1, T_VertexSet *s2)
   /* copies the set s1 elementwise to s2 */
   
{  int	           i;

   clear_set (s2);
   for (i = 0; i <= s1.lastel; i++) add_element (s2, s1.loe [i]);
}

/****************************************************************************************/

void renumber_vertices ()
   /* The vertices in the global variable G_Vertices are renumbered in the following    */
   /* way: If the i-th vertex is contained in k hyperplanes, it gets the number - k so  */
   /* that highly degenerate vertices get lower numbers. Then the set is reordered and  */
   /* numbered from 0 on. If it is not void, the global variable G_Incidence is         */
   /* changed accordingly.                                                              */


{
   int i, j, k;
   T_Vertex *tmp_vertex;
   boolean *tmp_incidence;

   for (i = 0; i < G_n; i++)
      G_Vertices.loe [i] -> no = 0;

   for (i = 0; i < G_n; i++)
      for (j = 0; j < G_m; j++)
         if (G_Incidence [i][j])
            G_Vertices.loe [i] -> no --;

   /* sort the set of vertices */
   for (i = 0; i < G_n - 1; i++)
   {
      /* look for the smallest vertex with index j >= i and swap with i */
      j = i;
      for (k = i + 1; k < G_n; k++)
         if (G_Vertices.loe [k] -> no < G_Vertices.loe [j] -> no)
            j = k;
      if (i != j)
      {
         tmp_vertex = G_Vertices.loe [i];
         G_Vertices.loe [i] = G_Vertices.loe [j];
         G_Vertices.loe [j] = tmp_vertex;
         if (G_Incidence != NULL)
         {
            tmp_incidence = G_Incidence [i];
            G_Incidence [i] = G_Incidence [j];
            G_Incidence [j] = tmp_incidence;
         }
      }
   }

   for (i = 0; i < G_n; i++)
      G_Vertices.loe [i] -> no = i;

   printf ("Vertices reordered.\n");
}

/****************************************************************************************/

rational normalise_vertices ()
   /* The vertices are expected in the global variable G_Vertices. In each dimension    */
   /* they are multiplied by a factor so that the entries are between -1 and 1. The     */
   /* return value is the product over the scaling factors in all dimensions; it is the */
   /* scaling factor for the volume. G_d must be known!                                 */
   
{
   rational *scaling = create_vector ();
   rational scaling_volume = 1;
   rational absolute;
   int      i, j;
   
   for (j = 0; j < G_d; j ++)
      scaling [j] = -1;

   for (j = 0; j < G_d; j++)
      for (i = 0; i <= G_Vertices.lastel; i++)
      {  absolute = fabs (G_Vertices.loe [i] -> coords [j]);
         if (absolute > scaling [j])
            scaling [j] = absolute;
      }

   for (j = 0; j < G_d; j++)
   {  scaling_volume *= scaling [j];
      for (i = 0; i <= G_Vertices.lastel; i++)
         G_Vertices.loe [i] -> coords [j] /= scaling [j];
   }
   
   return scaling_volume;
}

/****************************************************************************************/

void print_set (FILE *f, T_VertexSet s)
   /* prints the numbers of the elements of s to the (already opened) file f */

{  int i;

   if (empty (s)) fprintf (f, "\nThe set is empty.");
   else
      for (i=0; i <= (s.lastel); i++)
         fprintf (f, "%i\t", s.loe [i] -> no);
}

/****************************************************************************************/

boolean empty (T_VertexSet s)

{
   return (s.lastel == -1);
}

/****************************************************************************************/

boolean is_in_hyperplane (T_Vertex *e, int j)

{
   return G_Incidence [e -> no][j];
}

/****************************************************************************************/

boolean are_equal_sets (T_VertexSet s1, T_VertexSet s2)

{
   if (s1.lastel != s2.lastel) return FALSE;
   else if (empty (s1)) return TRUE;
   else
      return (!(memcmp (s1.loe, s2.loe, (s1.lastel + 1) * sizeof (T_Vertex *))));

}

/****************************************************************************************/

void add_element (T_VertexSet *s, T_Vertex *e)
   /* adds the element e to the set s. */

{  int i, pos, first = 0, last = s -> lastel, e_no = e -> no;
   int middle, middle_no;

   if (last == s -> maxel) redimloe (s);
   if (empty (*s)) 
   {  s -> lastel = 0;
      s -> loe [0] = e;
   }
   else if (e_no > s -> loe [last] -> no)
      /* insert element at the end of the list */
      s -> loe [++ s -> lastel] = e;
   else if (e_no < s -> loe [0] -> no)
   {  /* insert element at the beginning of the list */
      for (i = last; i >= 0; i--) s -> loe [i+1] = s -> loe [i];
      s -> loe [0] = e;
      s -> lastel++;
   }
   else 
   {  /* The element will be inserted after position pos. */
      /* look for pos between first and last */
      while (last - first > 1)
      {  middle = (first + last) / 2;
         middle_no = (*s).loe [middle] -> no;
         if (middle_no >= e_no)
            last = middle;
         else
            first = middle;
      }
      pos = first;
      for (i = s -> lastel; i > pos; i--) s -> loe [i+1] = s -> loe [i];
      s -> loe [pos + 1] = e;
      s -> lastel++;
   }
}

/****************************************************************************************/

boolean delete_element (T_VertexSet *s, T_Vertex *e)
   /* removes the element e from the set s. Return value is TRUE if the element actual- */
   /* ly was in the set, FALSE otherwise. */

{  int j, position;

   /* look for element in the list */
   position = position_of_element (e, *s);
   if (position == -1) return FALSE;
   else
   {  for (j = position; j < s -> lastel; j++)
         s -> loe [j] = s -> loe [j+1];
      s -> lastel --;
      return TRUE;
   }
}

/****************************************************************************************/

void intersect_with_hyperplane (T_VertexSet s, int j, T_VertexSet *inter)
   /* stores the intersection of s with hyperplane j in inter */

{
   int      i;

   inter -> lastel = -1;
   if (!(empty (s)))
      for (i = 0; i <= s.lastel; i++)
         if (is_in_hyperplane (s.loe [i], j))
            {  /* add_element (inter, s.loe [i]); */
               if (inter -> lastel == inter -> maxel) redimloe (inter);
               inter -> loe [++ inter -> lastel] = s.loe [i];
            }
}

/****************************************************************************************/
/*                 functions on sets of sets, called 'supersets'                        */
/****************************************************************************************/

T_VertexSuperset *create_empty_superset ()
   /* returns an empty set of sets */
   
{  return NULL;
}
  
/****************************************************************************************/

void free_superset (T_VertexSuperset **S)
   /* frees the memory space used by S. */
{
   T_VertexSuperset *L;
	
   while (*S != NULL)
   {  
      L = *S;
      (*S) = (*S) -> next;
      free_set (L -> content);
      my_free (L, sizeof (T_VertexSuperset));
   }
}
   
/****************************************************************************************/

void print_superset (FILE *f, T_VertexSuperset *S)
   /* prints the elements of S to the (already opened) file f */
   
{
   if (S == NULL) printf ("The set of sets is empty.\n");
   else
   {  printf ("The set of sets contains the following sets:\n");
      while (S != NULL)
      {  print_set (f, S -> content);
         S = S -> next;}
   }
}

/****************************************************************************************/

boolean is_in_superset (T_VertexSet s, T_VertexSuperset *S)
   /* tests whether s is contained in S */

{  while (S != NULL && !are_equal_sets (S -> content, s)) S = S -> next;
   if (S == NULL) return FALSE;
   else return TRUE;
}
      
/****************************************************************************************/

void add_superelement (T_VertexSuperset **S, T_VertexSet s)
   /* adds the element s at the beginning to the set S. */
   
{  T_VertexSuperset *L;
   
   L = (T_VertexSuperset *) my_malloc (sizeof (T_VertexSuperset));
   L -> content = s;
   L -> next = *S;
   *S = L;
}
   
/****************************************************************************************/
/****************************************************************************************/
