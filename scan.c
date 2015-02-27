#include <Python.h>
#include <structmember.h>

#include <signal.h>

typedef struct {
  FILE *fd;

  /* Set of poly_id's to examine. */
  int start;
  int end;

  /* pid for class_x */
  pid_t class_x_pid;

  /* poly.x command to run */
  char *poly_x_cmd;

  /* Python function to be called for each face, should return a tuple
     of points. */
  PyObject *face_cb;

  /* Python function to be called for each polytope, should return the
     lattice points on the polytope. */
  PyObject *poly_cb;

  /* A callback for checking whether the face is the target model. */
  PyObject *three_sector_cb;

  /* A set of reference 2d models (vertices in the 2d lattice) */
  PyObject *dP0_cube;
  PyObject *dP0_dP0_dP1;

  /* Id for the polytope currently under consideration. */
  int poly_id;
} threadinfo;

typedef struct {
  int vertices[2];
  int length;
} edge_t;

typedef struct {
  int nvertices;

  int *vertices;

  /* Total length of the boundary in lattice units */
  int boundary_length;

  /* Number of interior points in each edge, in traversal order */
  int *edge_ninterior;

  /* Classication of the face, if known. */
  int is_conifold;

  int is_Laba;
  int Laba_a;
  int Laba_b;

  int is_Z2xZ2_orbifold;

  int is_N2_orbifold;
  int N2_n;
  
  int is_Ypq;
  int Y_p;
  int Y_q;

  int is_dP_n;
  int dP_n;

  int is_3dP0s;

  int is_5dP0s;

  int is_dP0dP0dP1;
} face_t;

/* Object keeping statistics of the faces found. */
typedef struct {
  long long *face_count; /* The entries of this array measure the
			    number of faces having the given number of
			    interior points. */
  int max_interior;

  long long num_conifold;

  long long num_Z2xZ2;

  int max_Laba_b;
  long long **num_Laba; /* num_Laba[b][a], (note the ordering) a goes
			   from 0 to b. */

  int max_N2_n;
  long long *num_N2_n;

  long long num_dP_n[4];

  int max_p;
  long long **num_Ypq; /* num_Ypq[p][q], q goes from 0 to p-1. */

  /* Number of faces containing the model with 3 dP0s */
  long long num_3dP0s;

  /* Number of faces containing the model with 5 dP0s */
  long long num_5dP0s;

  /* Number of faces containing the model with 2 dP0s and one dP1 */
  long long num_dP0dP0dP1;

  /* Examples of polytopes */
  PyObject **very_singular_polys;
  int *very_singular_npoints;
  int nvery_singular;

  PyObject *Ypq_examples;

  PyObject *conifold_example;

  PyObject *Laba_example;

  PyObject *dP0dP0dP1_examples;

  PyObject *dP0_cube_examples;

  PyObject *dP0_5_examples;

  PyObject *dP_n_example[4];

  PyObject *N2_examples;
} stats_t;

#if !defined(__linux__)
#define getline _getline
static int _print_debug_info = 1;
#endif

static ssize_t __attribute__ ((unused)) _getline(char **lineptr, size_t *n, FILE *stream)
{
  char *bufptr = NULL;
  char *p = bufptr;
  size_t size;
  int c;

#if !defined(__linux__)
  if (_print_debug_info) {
    printf("Non-linux version detected, using alternative getline\n");
    _print_debug_info = 0;
  }
#endif

  if (lineptr == NULL) {
    return -1;
  }
  if (stream == NULL) {
    return -1;
  }
  if (n == NULL) {
    return -1;
  }

  bufptr = *lineptr;
  size = *n;

  if (bufptr == NULL) {
    bufptr = malloc(128);
    if (bufptr == NULL) {
      return -1;
    }
    size = 128;
  }

  p = bufptr;

  while (1) {
    c = fgetc(stream);
    if (c == EOF)
      break;

    if ((p - bufptr) >= (size - 1)) {
      long int offset = p-bufptr;

      size = size + 128;
      bufptr = realloc(bufptr, size);
      if (bufptr == NULL) {
	return -1;
      }

      p = bufptr+offset;
    }

    *p = c;
    p++;

    if (c == '\n') {
      break;
    }
  }

  *p++ = '\0';
  *lineptr = bufptr;
  *n = size;

  if ((c == EOF) && (p - bufptr - 1 == 0))
    return -1;

  return p - bufptr - 1;
}

static void print_vector(int *v, threadinfo *info)
{
  int i;

  printf("[");
  for (i=0; i<4; i++) {
    printf("%d", v[i]);
    if (i != 3)
      printf(", ");
  }
  printf("]");
}

static int dot(int *a, int *b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
}

/* Compute the gcd of a and b using Euclid's algorithm. Notice that in
   the degenrate case that both a and b are 0 this returns 0, and it
   returns a positive integer otherwise. */
static unsigned int gcd(int a, int b)
{
  int remainder;

  a = abs(a);
  b = abs(b);

  if (a < b) {
    int c = a;
    a = b;
    b = c;
  }

  if (b == 0)
    return a;

  do {
    remainder = a%b;

    a = b;
    b = remainder;
    
  } while (remainder != 0);

  return a;
}

/* Compute the gcd of the components of a 4d vector given by the
   difference of a and v */
static unsigned int gcd_delta(int *a, int *b)
{
  return gcd(a[0] - b[0], gcd(a[1] - b[1], gcd(a[2]-b[2], a[3]-b[3])));
}

static int matrix_minor(int *a, int *b, int *c, int i, int j, int k)
{
  return a[i]*(b[j]*c[k] - b[k]*c[j]) - a[j]*(b[i]*c[k] - b[k]*c[i])
    + a[k]*(b[i]*c[j] - b[j]*c[i]);
}

/* Construct the 4x3 matrix given by the 3 4-dimensional vectors, and
   return 1 if one of the dimension 3 minors is non-vanishing, 0
   otherwise. */
static int is_rank_3(int *a, int *b, int *c)
{
  if (matrix_minor(a,b,c,1,2,3) != 0 ||
      matrix_minor(a,b,c,0,2,3) != 0 ||
      matrix_minor(a,b,c,0,1,3) != 0 ||
      matrix_minor(a,b,c,0,1,2) != 0)
    return 1;

  return 0;
}

/* Given a list of vertices, construct a wrapping Python object */
static PyObject*
list_from_vertices(int nvertices, int vertices[][4])
{
  int i, j;
  PyObject *list = PyTuple_New(nvertices);

  for (i=0; i<nvertices; i++) {
    PyObject *v = PyTuple_New(4);
    for (j=0; j<4; j++)
      PyTuple_SetItem(v, j, PyInt_FromLong(vertices[i][j]));
    PyTuple_SetItem(list, i, v);
  }

  return list;
}

static PyObject*
new_poly(int nvertices, int vertices[][4], int poly_id)
{
  PyObject *p = PyDict_New();
  PyObject *id = PyInt_FromLong(poly_id);
  PyObject *v = list_from_vertices(nvertices, vertices);

  PyDict_SetItemString(p, "poly_id", id);
  Py_DECREF(id);

  PyDict_SetItemString(p, "vertices", v);
  Py_DECREF(v);

  return p;
}

/* Create a two-way pipe for communication with poly.x */
static void poly_x(FILE **in, FILE **out, const char *poly_x_cmd,
		   int compute_points)
{
  int input[2], output[2];
  int child;

  if (pipe(input)<0 || pipe(output) < 0) {
    perror("poly_x");
    
    assert(0);
  }

  child = fork();
  if (child<0) {
    perror("poly_x (fork)");
    assert(0);
  } else if (child == 0) { /* in the child */
    close(input[1]);
    close(output[0]);

    dup2(input[0], 0);
    dup2(output[1], 1);

    close(input[0]);
    close(output[1]);

    if (execlp(poly_x_cmd, poly_x_cmd, compute_points?"-fip":"-fi", NULL) < 0) {
      perror("poly_x execl");
      fprintf(stderr, "%s\n", poly_x_cmd);
      assert(0);
    }
  } else { /* in the parent */
    close(input[0]);
    close(output[1]);
  }

  *in = fdopen(input[1], "w");
  *out = fdopen(output[0], "r");
}

/* Create a one-way (read) pipe for communicating with class_x. */
static FILE* class_x(const char *class_x_cmd, const char *zzdb, int *pid)
{
  int output[2];
  int child;

  if (pipe(output) < 0) {
    perror("class_x");
    
    assert(0);
  }

  child = fork();
  if (child<0) {
    perror("class_x (fork)");
    assert(0);
  } else if (child == 0) { /* in the child */
    close(output[0]);

    dup2(output[1], 1);

    close(output[1]);
    
    if (execlp(class_x_cmd, class_x_cmd, "-b2a", "-di", zzdb, NULL) < 0) {
      perror("class_x execl");
      fprintf(stderr, "%s -b2a -di %s\n", class_x_cmd, zzdb);
      assert(0);
    }
  } else { /* in the parent */
    close(output[1]);
  }

  *pid = child;

  return fdopen(output[0], "r");
}

/* Returns 1 if the given vertex is in the boundary of the given
   face, 0 otherwise. */
static int vertex_in_face(int vertex, int *face_vertices, int nvertices) {
  int i;

  for (i=0; i<nvertices; i++) {
    if (face_vertices[i] == vertex)
      return 1;
  }

  return 0;
}

/* Order the given set of vertices in boundary order, and store the
   result in face.vertices. This also computes boundary_length. */
static void traverse_boundary(edge_t *edges, int nedges,
			      face_t *face, int *face_vertices)
{
  int prev2;
  int j;

  face->vertices[0] = face_vertices[0];
  face->boundary_length = 0;

  /* This the vertex the _previous_ vertex in the sequence was
     connected to (i.e., the previous of the previous). For the
     first vertex we don't have that information yet, so set to some
     invalid value. The routine below will then do the right
     thing. */
  prev2 = -1;

  /* previous vertex in the sequence */
  for (j=1; j<=face->nvertices; j++) {
    int k;

    /* Find an edge in the polytope which connects to the previous
       vertex, and is not prev2. */
    for (k=0; k<nedges; k++) {
      if (edges[k].vertices[0] == face->vertices[j-1]
	  && edges[k].vertices[1] != prev2
	  && vertex_in_face(edges[k].vertices[1],
			    face_vertices, face->nvertices)) {
	if (j < face->nvertices)
	  face->vertices[j] = edges[k].vertices[1];
	else
	  assert(edges[k].vertices[1] == face->vertices[0]);
	break;
      } else if (edges[k].vertices[1] == face->vertices[j-1]
	  && edges[k].vertices[0] != prev2
	  && vertex_in_face(edges[k].vertices[0],
			    face_vertices, face->nvertices)) {
	if (j < face->nvertices)
	  face->vertices[j] = edges[k].vertices[0];
	else
	  assert(edges[k].vertices[0] == face->vertices[0]);
	break;
      }
    }
    
    /* We should have found something */
    assert(k<nedges);

    face->edge_ninterior[j-1] = edges[k].length-1;

    face->boundary_length += edges[k].length;

    prev2 = face->vertices[j-1];
  }
}

/* Allow python code to check whether the object is our target. */
static int python_three_sector_cb(int vertices[][4], face_t *face,
				  PyObject *ref_model,
				  PyObject *three_sector_cb)
{
  int i;
  PyObject *args, *verts, *result;
  int is_3dP0s;

  verts = PyTuple_New(face->nvertices);

  for (i=0; i<face->nvertices; i++) {
    int *coords = vertices[face->vertices[i]];
    PyObject *v = Py_BuildValue("(iiii)",
				coords[0], coords[1], coords[2], coords[3]);

    PyTuple_SetItem(verts, i, v);
  }

  args = Py_BuildValue("(OO)", ref_model, verts);

  result = PyObject_CallObject(three_sector_cb, args);

  Py_DECREF(verts);
  Py_DECREF(args);

  if (result == NULL) {
    printf("Got an exception from Python\n");
    PyErr_Print();
    return -1;
  }

  if (!PyBool_Check(result)) {
    printf("The object three_sector_cb returned is not a boolean!\n");
    return -1;
  }

  is_3dP0s = (result == Py_True);

  Py_DECREF(result);
  
  return is_3dP0s;
}

/* Process a single face, once we have read its structure. Note that
   here vertices refers to the polytope, the set of vertices for the
   face is kept in face.  If the face can be identified as a dP_n or a
   Y^{p,q} the corresponding values will be saved in the passed
   pointers, and they will be otherwise untouched.
*/
static int npoints_in_face(int vertices[][4], int nvertices,
			   face_t *face,
			   int points[][4], int npoints, int debug,
			   threadinfo *info)
{
  int i, count;
  int v0[4], v1[4], vn[4], delta[4];
  int nboundary_interior = face->boundary_length - face->nvertices;
  int ninterior;

  /* Points in the face are those coplanar with the face. */
  for (i=0; i<4; i++) {
    v0[i] = vertices[face->vertices[0]][i];
    v1[i] = vertices[face->vertices[1]][i] - v0[i];
    vn[i] = vertices[face->vertices[face->nvertices-1]][i] - v0[i];
  }

  count = 0;
  for (i=0; i<npoints; i++) {
    int j;

    for (j=0; j<4; j++)
      delta[j] = points[i][j] - v0[j];

    if (!is_rank_3(v1,vn,delta))
      count++;
  }

  ninterior = count - face->boundary_length;

  /* Check for the conifold */
  if (face -> nvertices == 4 && nboundary_interior == 0 &&
      ninterior == 0) {
    face -> is_conifold = 1;
    face -> is_Laba = 1;
    face -> Laba_a = 1;
    face -> Laba_b = 1;
  }

  /* L^{aba} manifolds: hep-th/0505211 */
  if (face -> nvertices == 4 && ninterior == 0 && nboundary_interior > 0) {
    for (i=0; i<4; i++) {
      if (face->edge_ninterior[i] > 0) {
	int a, b;
	a = face->edge_ninterior[i]+1;
	b = face->edge_ninterior[(i+2)%4]+1;

	if (a>b) {
	  int tmp=a;
	  a = b;
	  b = tmp;
	}
	
	face->is_Laba = 1;
	face->Laba_a = a;
	face->Laba_b = b;

	/* Double check, should always succeed */
	if (a+b+2 != face->boundary_length) {
	  int j;

	  face->is_Laba = 0;

	  printf("Error [%d]: misidentified face as L^{%d,%d,%d}.\n",
		 info->poly_id, a,b,a);

	  for (j=0; j<face->nvertices;j++) {
	    print_vector(vertices[face->vertices[j]], info);
	    if (j!=face->nvertices-1)
	      printf(", ");
	  }
	  printf("\n");

	  printf("Polytope:\n");
	  for (i=0; i<nvertices;i++) {
	    print_vector(vertices[i], info);
	    if (i!=nvertices-1)
	      printf(", ");
	  }
	  printf("\n");
	}

	break;
      }
    }
  }

  /* C^3/(Z_2xZ_2) */
  if (face->nvertices == 3 && ninterior == 0 && face->boundary_length == 6) {
    if (face->edge_ninterior[0] == 1) {
      face->is_Z2xZ2_orbifold = 1;
    }

    /* Double check, should always suceed */
    if (face->is_Z2xZ2_orbifold && face->edge_ninterior[1] != 1) {
      face->is_Z2xZ2_orbifold = 0;
	  
      printf("Error [%d]: misidentified face as C^3/(Z_2xZ_2)\n",
	     info->poly_id);

      for (i=0; i<face->nvertices;i++) {
	print_vector(vertices[face->vertices[i]], info);
	if (i!=face->nvertices-1)
	  printf(", ");
      }
      printf("\n");
      printf("Polytope:\n");
      for (i=0; i<nvertices;i++) {
	print_vector(vertices[i], info);
	if (i!=nvertices-1)
	  printf(", ");
      }
      printf("\n");
    }
  }

  /* C x C^2/Z_n */
  if (face->nvertices == 3 && ninterior == 0 &&
      !face->is_Z2xZ2_orbifold) {
    for (i=0; i<3; i++) {
      if (face->edge_ninterior[i] != 0)
	break;
    }

    if (i<3) {
      /* Found an edge with interior points, this is the orbifold edge. */
      face->is_N2_orbifold = 1;
      face->N2_n = face->edge_ninterior[i]+1;

      face->is_Laba = 1;
      face->Laba_a = 0;
      face->Laba_b = face->edge_ninterior[i]+1;

      /* Double check, should always succeed */
      if (face->edge_ninterior[i]+3 != face->boundary_length) {
	  face->is_Laba = 0;
	  face->is_N2_orbifold = 0;

	  printf("Error [%d]: misidentified face as C x C^2/Z_%d.\n",
		 info->poly_id, face->N2_n);

	  for (i=0; i<face->nvertices;i++) {
	    print_vector(vertices[face->vertices[i]], info);
	    if (i!=face->nvertices-1)
	      printf(", ");
	  }
	  printf("\n");

	  printf("Polytope:\n");
	  for (i=0; i<nvertices;i++) {
	    print_vector(vertices[i], info);
	    if (i!=nvertices-1)
	      printf(", ");
	  }
	  printf("\n");
      }
    }
  }

  /* Check for Y^{p,q}: hep-th/0411238 */
  if (face -> nvertices == 4 && nboundary_interior == 0
      && ninterior > 0) {
      int *va, *vb, *vc, *v0 = NULL; 

      if (gcd_delta(vertices[face->vertices[0]],
		    vertices[face->vertices[2]]) == ninterior+1) {
	v0 = vertices[face->vertices[0]];
	va = vertices[face->vertices[2]];
	vb = vertices[face->vertices[1]];
	vc = vertices[face->vertices[3]];
      } else if (gcd_delta(vertices[face->vertices[1]],
			   vertices[face->vertices[3]]) == ninterior+1) {
	v0 = vertices[face->vertices[1]];
	va = vertices[face->vertices[3]];
	vb = vertices[face->vertices[0]];
	vc = vertices[face->vertices[2]];
      }
      
      if (v0 != NULL) {
	/* All the interior points are along a diagonal, so we found a
	   Y^{p,q}, with p = ninterior+1. The following code
	   determines l = p-q. */
	int a[4], b[4], c[4];
	int num, den;
	int a2, ab, b2;
	int l;
	int p = ninterior+1;

	for (i=0; i<4; i++) {
	  a[i] = (va[i] - v0[i])/(ninterior+1);
	  b[i] = vb[i] - v0[i];
	  c[i] = vc[i] - v0[i];
	}

	a2 = dot(a,a);
	b2 = dot(b,b);
	ab = dot(a,b);

	num = dot(a,c)*b2 - dot(b,c)*ab;
	den = a2*b2 - ab*ab;

	l = num/den;

	if (l > p)
	  l = 2*p - l;

	face->is_Ypq = 1;
	face->Y_p = p;
	face->Y_q = p - l;

	/* Double checks, should always succeed */
	if (gcd(num, den) != den) {
	  face->is_Ypq = 0;

	  printf("Error: the denominator %d does not divide %d.\n", den, num);

	  for (i=0; i<4;i++) {
	    print_vector(vertices[face->vertices[i]], info);
	    if (i!=3)
	      printf(", ");
	  }
	  printf("\n");
	}

	if (l < 1 || l > 2*p) {
	  face->is_Ypq = 0;

	  printf("Value for l out of range: %d [p=%d]\n", l, p);

	  for (i=0; i<4;i++) {
	    print_vector(vertices[face->vertices[i]], info);
	    if (i!=3)
	      printf(", ");
	  }
	  printf("\n");
	}
      }
  }
 
  /* Check for the del Pezzo's */
  if ((nboundary_interior == 0) && (ninterior == 1)) {
    int n = face->nvertices - 3;

    /* Distinguish F_0 from dP_1. These are Y^{2,0} and Y^{2,1}
       respectively, so the hard work has already been done by the
       code above. */
    if (n != 1 || face->Y_q == 1) {
      face->is_dP_n = 1;
      face->dP_n = n;

      if (n > 3) {
	int j;

	printf("n is too big for toric del Pezzo: %d\n", n);
	face->is_dP_n = 0;

	for (j=0; j<face->nvertices;j++) {
	  print_vector(vertices[face->vertices[j]], info);
	  if (j!=face->nvertices-1)
	    printf(", ");
	}
	printf("\n");
      }
    }
  }

  /* Check for the model with 3 dP0's */
  if ((nboundary_interior == 1) && (ninterior == 3) &&
      (face->nvertices == 4)) {
    int *va, *vb, *vc, *vd;

    /* Find the edge with the interior point */
    for (i=0; i<face->nvertices; i++) {
      if (face->edge_ninterior[i] == 1)
	break;
    }

    va = vertices[face->vertices[i]];
    vb = vertices[face->vertices[(i+1)%face->nvertices]];
    vc = vertices[face->vertices[(i+2)%face->nvertices]];
    vd = vertices[face->vertices[(i+3)%face->nvertices]];

    if ((gcd_delta(va, vc) == 3) && (gcd_delta(vb, vd) == 3)) {
      /* Found our geometry. To see this, notice that 3 interior
	 points arranged in two lines (the diagonals between ac, bd
	 that we just computed) can always be put in the coners of an
	 elementary triangle of the lattice. Once they are in such a
	 triangle, we can reconstruct the whole polytope in a unique
	 way, so there is a single lattice polytope satisfying the
	 constraints so far.  */
      face->is_3dP0s = 1;

      if (info->three_sector_cb != NULL &&
	  !python_three_sector_cb(vertices, face, info->dP0_cube,
				  info->three_sector_cb)) {
	printf("Found a 3-sector model, but Python disagrees\n");
	for (i=0; i<4;i++) {
	  print_vector(vertices[face->vertices[i]], info);
	  if (i!=3)
	    printf(", ");
	}
	printf("\n");
      }
    } else {
      if (info->three_sector_cb != NULL &&
	  python_three_sector_cb(vertices, face, info->dP0_cube,
				 info->three_sector_cb)) {
	printf("Discarded a face, but Python disagrees\n");
	for (i=0; i<4;i++) {
	  print_vector(vertices[face->vertices[i]], info);
	  if (i!=3)
	    printf(", ");
	}
	printf("\n");
      }
    }
  } 

  /* Check for the model with 5 dP0's*/
  if ((nboundary_interior == 3) && (ninterior == 5) &&
      (face->nvertices == 4)) {
    int *va, *vb, *vc, *vd;

    /* Find the long edge */
    for (i=0; i<face->nvertices; i++) {
      va = vertices[face->vertices[i]];
      vb = vertices[face->vertices[(i+1)%face->nvertices]];
      vc = vertices[face->vertices[(i+2)%face->nvertices]];
      vd = vertices[face->vertices[(i+3)%face->nvertices]];

      if ((gcd_delta(va,vb) == 3) &&
	  (gcd_delta(vc,vd) == 2))
	break;
    }

    if (i != face->nvertices) {
      /* Found a candidate, check the middle points */
      int vab1[4], vab2[4], vcd1[4];

      for (i=0; i<4; i++) {
	vab1[i] = va[i] + ((vb[i]-va[i])/3);
	vab2[i] = va[i] + ((vb[i]-va[i])/3)*2;
	vcd1[i] = vc[i] + ((vd[i]-vc[i])/2);
      }

      if ((gcd_delta(va, vcd1) == 3) &&
	  (gcd_delta(vb, vcd1) == 3) &&
	  (gcd_delta(vab1, vc) == 3) &&
	  (gcd_delta(vab2, vd) == 3)) {
	face->is_5dP0s = 1;
      }
    }
  }

  /* Check for the model with 2 dP0's and one dP1 */
  if ((nboundary_interior == 1) && (ninterior == 3) &&
      (face->nvertices == 5)) {
    int *va, *vb, *vc, *vd, *ve;

    /* Find the edge with the interior point */
    for (i=0; i<face->nvertices; i++) {
      if (face->edge_ninterior[i] == 1)
	break;
    }

    va = vertices[face->vertices[i]];
    vb = vertices[face->vertices[(i+1)%face->nvertices]];
    vc = vertices[face->vertices[(i+2)%face->nvertices]];
    vd = vertices[face->vertices[(i+3)%face->nvertices]];
    ve = vertices[face->vertices[(i+4)%face->nvertices]];

    if ((gcd_delta(va, vc) == 3) && (gcd_delta(vb, vd) == 3)) {
      int mp[4];

      face->is_dP0dP0dP1 = 1;

      /* Double check, should always succeed */
      for (i=0; i<4; i++) {
	mp[i] = va[i] + (vb[i]-va[i])/2;
      }

      if (gcd_delta(mp, ve) != 2) {
	int j;

	face->is_dP0dP0dP1 = 0;

	printf("Misidentified a face as dP0^2 x dP1\n");

	for (j=0; j<face->nvertices;j++) {
	  print_vector(vertices[face->vertices[j]], info);
	  if (j!=face->nvertices-1)
	    printf(", ");
	}
	printf("\n");

      }
    } else if ((gcd_delta(va, vd)==3) && (gcd_delta(vb, ve) == 3)) {
      int mp[4];

      face->is_dP0dP0dP1 = 1;

      /* Double check, should always succeed */
      for (i=0; i<4; i++) {
	mp[i] = va[i] + (vb[i]-va[i])/2;
      }

      if (gcd_delta(mp, vc) != 2) {
	int j;

	face->is_dP0dP0dP1 = 0;

	printf("Misidentified a face as dP0^2 x dP1\n");

	for (j=0; j<face->nvertices;j++) {
	  print_vector(vertices[face->vertices[j]], info);
	  if (j!=face->nvertices-1)
	    printf(", ");
	}
	printf("\n");

      }
    } else {
      /* These are not the models you are looking for */
    }

    if (info->three_sector_cb != NULL &&
	(!!python_three_sector_cb(vertices, face, info->dP0_dP0_dP1,
				  info->three_sector_cb) != !!face->is_dP0dP0dP1)) {
      printf("3-sector dP0^2xdP1 model is %d, but Python disagrees\n",
	     face->is_dP0dP0dP1);
      for (i=0; i<face->nvertices;i++) {
	print_vector(vertices[face->vertices[i]], info);
	if (i!=face->nvertices-1)
	  printf(", ");
      }
      printf("\n");

    }
  }

  return count;
}

/* Do the same thing through the supplied callback. */
static int python_npoints_in_face(int vertices[][4], face_t *face,
				  PyObject *face_cb)
{
  int i;
  PyObject *args, *verts, *result;
  int npoints;

  verts = PyTuple_New(face->nvertices);

  for (i=0; i<face->nvertices; i++) {
    int *coords = vertices[face->vertices[i]];
    PyObject *v = Py_BuildValue("(iiii)",
				coords[0], coords[1], coords[2], coords[3]);

    PyTuple_SetItem(verts, i, v);
  }

  args = Py_BuildValue("(O)", verts);

  result = PyObject_CallObject(face_cb, args);

  Py_DECREF(verts);
  Py_DECREF(args);

  if (result == NULL)
    return -1;

  if (!PySequence_Check(result)) {
    printf("The object face_cb returned is not a sequence!\n");
    return -1;
  }

  npoints = (int)PySequence_Length(result);

  if (npoints < 0) {
    printf("Getting the number of points in python failed!\n");
    return -1;
  }

  Py_DECREF(result);
  
  return npoints;
}

/*
  Read the list of the integral points in the convex hull of the
  vertices by calling poly_cb, which should return a sequence of
  coordinates for the points, each coordinate given by a sequence of
  4 numbers.
  
  On success this function returns the number of points read, or -1 in
  case of error. points is allocated for the caller, who is
  responsible for freeing the returned memory.
*/
static int
read_polytope_points_from_python(int vertices[][4], int nvertices,
				 int (**points)[4], threadinfo *info)
{
  int i;
  PyObject *args, *verts, *result;
  int npoints;

  verts = PyList_New(nvertices);
  for (i=0; i<nvertices; i++) {
    int *coords = vertices[i];
    PyObject *v = Py_BuildValue("(iiii)",
				coords[0], coords[1], coords[2], coords[3]);
    
    PyList_SetItem(verts, i, v);
  }
  args = Py_BuildValue("(O)", verts);

  result = PyObject_CallObject(info->poly_cb, args);

  Py_DECREF(verts);
  Py_DECREF(args);

  if (result == NULL)
    return -1;

  if (!PySequence_Check(result)) {
    printf("The object poly_cb returned is not a sequence!\n");
    return -1;
  }

  npoints = (int)PySequence_Length(result);

  if (npoints < 0) {
    printf("Getting the number of points in python failed!\n");
    return -1;
  }

  *points = malloc(npoints*sizeof(int[4]));

  for (i=0; i<npoints; i++) {
    PyObject *coords = PySequence_GetItem(result, i);
    int j;

    if (!PySequence_Check(coords)) {
      printf("%d-th point returned by poly_cb is not a sequence\n", i);
      Py_DECREF(coords);
      continue;
    }

    if ((int)PySequence_Length(coords) != 4) {
      printf("%d-th point returned by poly_cb has length different from 4\n",
	     i);
      Py_DECREF(coords);
      continue;
    }

    for (j=0; j<4; j++) {
      PyObject *pj = PySequence_GetItem(coords, j);

      if (!PyNumber_Check(pj)) {
	printf("coordinate is not a number\n");
	Py_DECREF(pj);
	continue;
      }

      (*points)[i][j] = (int)PyNumber_AsSsize_t(pj, NULL);

      Py_DECREF(pj);
    }

    Py_DECREF(coords);
  }

  Py_DECREF(result);

  return npoints;
}

static void
read_one_polytope(int rows, int cols, FILE* poly_x_in, FILE *poly_x_out,
		  threadinfo *info, int transposed,
		  int *max_face_points, stats_t *stats)
{
  char *line = NULL;
  size_t nline = 0;
  int i, j;
  edge_t *edges = NULL; int nedges = 0;
  face_t *faces = NULL; int nfaces = 0;
  int vertices[cols][rows];
  int input_rows = rows, input_cols = cols;
  int (*points)[4] = NULL;
  int npoints;
  PyObject *poly = NULL;
  int max_p=-1, max_q=-1;
  int has_conifold = 0;
  int has_L121212 = 0;
  int has_dP0dP0dP1 = 0;
  int has_dP0_cube = 0;
  int has_dP0_5 = 0;
  int has_dP_n[4] = {0,0,0,0};
  int *has_N2_n = NULL;
  int max_N2_n = -1;

  /* Palp (class.x) sometimes outputs polytopes transposed (that is, n
    rows x 4 cols, as opposed to its usual 4 rows x n cols
    format). That is the reason why in this routine we distinguish
    between rows, cols (which always refer to the canonical 4xn form)
    and input_rows, input_cols. */
  if (transposed) {
    input_rows = cols;
    input_cols = rows;
  }

  if (info->poly_id >= info->start)
    fprintf(poly_x_in, "%d %d\n", input_rows, input_cols);

  for (i=0; i<input_rows; i++) {
    if (getline(&line, &nline, info->fd) < 0) {
      printf("#! Assertion failed!\n");
      assert(0);
    }

    if (info->poly_id >= info->start) {
      char *ptr;

      fputs(line, poly_x_in);

      ptr = line;
      
      for (j=0; j<input_cols; j++) {
	char *next;
	long int n = strtol(ptr, &next, 10);
	
	if (next == ptr || n == LONG_MIN || n == LONG_MAX)
	  assert(0);
	
	if (!transposed)
	  vertices[j][i] = (int)n;
	else
	  vertices[i][j] = (int)n;
	
	ptr = next;
      }
    }
  }

  if (info->poly_id < info->start) {
    free(line);

    return;
  }

  /* poly_x_in is a block buffered file, so the previous code does
     not necessarily send the data down the pipe. Flush the buffer
     explicitly so the data is actually sent, and we can read the
     reply. We could also use setvbuf to set the file to unbuffered,
     but this seems to give some performance penalty. */
  fflush(poly_x_in);

  /* Read the list of lattice points in the polytope */
  if (info->poly_cb == NULL) {
    char *next;
    int is_transposed = 0;
    long int n;

    if (getline(&line, &nline, poly_x_out) < 0) {
      //printf("Error reading reply from poly_x\n");
      return;
    }

    n = strtol(line, &next, 10);

    if (next == line || n == LONG_MIN || n == LONG_MAX) {
      printf("Reading number of points failed!\n%s", line);
      return;
    }
    npoints = (int)n;

    /* PALP... */
    if (npoints == 4) {
      char *ptr;

      n = strtol(next, &ptr, 10);
      if (next == ptr || n == LONG_MIN || n == LONG_MAX || n < 6) {
	printf("Cannot make sense of the number of points!\n");
	return;
      }
      npoints = (int)n;

      is_transposed = 1;
    }

    points = malloc(npoints*sizeof(int[4]));

    if (!is_transposed) {
      for (i=0; i<npoints; i++) {
	char *ptr, *next;

	if (getline(&line, &nline, poly_x_out) < 0) {
	  printf("failed to read the list of points\n");
	  return;
	}

	ptr = line;
	for (j=0; j<4; j++) {
	  long int n = strtol(ptr, &next, 10);

	  if (next == ptr || n == LONG_MIN || n == LONG_MAX) {
	    printf("Cannot read point coordinates.\n");
	    return;
	  }

	  points[i][j] = (int)n;

	  ptr = next;
	}
      }
    } else {
      for (i=0; i<4; i++) {
	char *ptr, *next;

	if (getline(&line, &nline, poly_x_out) < 0) {
	  printf("failed to read the transposed list of points\n");
	  return;
	}

	ptr = line;
	for (j=0; j<npoints; j++) {
	  long int n = strtol(ptr, &next, 10);

	  if (next == ptr || n == LONG_MIN || n == LONG_MAX) {
	    printf("Cannot read permuted point coordinates.\n");
	    return;
	  }

	  points[j][i] = (int)n;

	  ptr = next;
	}
      }
    }
  } else {
    /* A callback was provided, get the points of the polytope from
       the Python callback instead. */
    npoints = read_polytope_points_from_python(vertices, cols, &points,
					       info);
    if (npoints <= cols) {
      printf("Invalid number of points from poly_cb: %d\n", npoints);
      return;
    }
  }

  /* Read the list of faces */
  /* The algorithm used here is a C translation of the one in the
     LatticePolytope package in sage. In particular the
     _read_poly_x_incidences and all_faces routines in that file
     describe in some detail the i/o format poly.x -fi used for
     communication. */

  /* We can ignore the first line, it is just some information */
  if (getline(&line, &nline, poly_x_out) < 0) {
    printf("Error here!\n");
    return;
  }

  /* This gives the list of faces in terms of vertices. */
  /* The first line can also be ignored. */
  if (getline(&line, &nline, poly_x_out) < 0) {
    fprintf(stderr, "Ough! %d\n", info->poly_id);
    return;
  }

  /* Process the list of vertices. It is organized by dimension, we
     just need the faces of dimension 2 and 1. Dimension 1 are the
     edges, we need this for being able to traverse the boundaries in
     order. */
  /* This information for each face is given in the shape of an
     incidence matrix: a structure of the form 01001..001 indicates
     that the n-0th vertex is not part of the face (0), the n-1th is
     part of it, the n-2th is not part, etc. (Notice that for some
     reason the order in which this information apears is opposite to
     the order one would expect.) */
  for (i = 0; i < rows; i++) {
    if (getline(&line, &nline, poly_x_out) < 0) {
      printf("Oups!!!<<<\n");
      return;
    }

    /* Process the edges */
    if (i==1 || i == 2) {
      char *p = line;
      int n = 0;
      int face_vertices[cols];
      int nvertex = 0;
      int done = 0;

      while (*p != ' ')
	p++;

      p++;

      while (!done) {
	switch (*p) {
	case '0':
	  n++;
	  break;
	case '1':
	  /* Add the vertex to the list */
	  face_vertices[nvertex++] = (cols - 1) - n;
	  n++;
	  break;
	case '\n':
	  done = 1;
	  if (n == 0)
	    break;
	  /* note the fall-through, it is on purpose. We want to
	     process the last set face/edge, in case the list didn't
	     end in a space */
	case ' ':
	  /* Completed a face/edge. Use the information that we read for
	     constructing the face/edge and add it to the list. */
	  assert(n == cols);

	  if (i==1) {
	    /* add an edge */
	    assert(nvertex == 2);
	    
	    edges = realloc(edges, sizeof(edge_t)*(nedges+1));
	    memcpy(edges[nedges].vertices, face_vertices, sizeof(int)*2);
	    edges[nedges].length = gcd_delta(vertices[face_vertices[0]],
					     vertices[face_vertices[1]]);
	    nedges++;
	  } else {
	    /* add a 2d face */
	    assert(nvertex > 2);

	    faces = realloc(faces, sizeof(face_t)*(nfaces+1));
	    faces[nfaces].nvertices = nvertex;
	    faces[nfaces].vertices = malloc(sizeof(int)*nvertex);
	    faces[nfaces].is_conifold = 0;
	    faces[nfaces].is_Laba = 0;
	    faces[nfaces].is_Z2xZ2_orbifold = 0;
	    faces[nfaces].is_N2_orbifold = 0;
	    faces[nfaces].is_dP_n = 0;
	    faces[nfaces].is_Ypq = 0;
	    faces[nfaces].is_3dP0s = 0;
	    faces[nfaces].is_5dP0s = 0;
	    faces[nfaces].is_dP0dP0dP1 = 0;
	    faces[nfaces].edge_ninterior = malloc(sizeof(int)*nvertex);
	    traverse_boundary(edges, nedges, &faces[nfaces], face_vertices);
	    nfaces++;
	  }

	  nvertex = 0;
	  n = 0;
	  break;
	default:
	  printf("This shouldn't happen!! (unexpected character in poly.x output)\n");
	  return;
	}

	p++;
      }
    }
  }

  /* This gives the list of faces in terms of facets. We do not use
     this, so we can ignore it. */
  if (getline(&line, &nline, poly_x_out) < 0) {
    printf("Too short output!\n");
    return;
  }

  for (i = 0; i < rows; i++) {
    if (getline(&line, &nline, poly_x_out) < 0) {
      printf("What!?\n");
      return;
    }
  }

  free(line);

  /* Process the faces */
  for (i=0; i<nfaces; i++) {
    /* Got a match */
    int n_interior;
    int n_total_points, n_total_points_py;

    n_total_points = npoints_in_face(vertices, cols,
				     &faces[i], points, npoints, 0, info);

    if (info->face_cb != NULL) {
      n_total_points_py = python_npoints_in_face(vertices,
						 &faces[i], info->face_cb);
   
      if (n_total_points != n_total_points_py)
	printf("Python says %d, C says %d\n", n_total_points_py, n_total_points);
    }

    n_interior = n_total_points - faces[i].boundary_length;

    if (n_interior > stats->max_interior) {
      stats->face_count = realloc(stats->face_count,
				  (n_interior+1)*sizeof(long long));
      memset(stats->face_count+stats->max_interior+1, 0,
	     (n_interior-stats->max_interior)*sizeof(long long));

      stats->max_interior = n_interior;
    }

    if (n_interior > *max_face_points)
      *max_face_points = n_interior;

    stats->face_count[n_interior]++;

    if (faces[i].is_conifold) {
      has_conifold = 1;
      stats->num_conifold++;
    }

    if (faces[i].is_dP_n) {
	stats->num_dP_n[faces[i].dP_n]++;
	has_dP_n[faces[i].dP_n] = 1;
    }

    if (faces[i].is_Ypq) {
      int p = faces[i].Y_p;
      int q = faces[i].Y_q;

      if (stats->max_p < p) {
	stats->num_Ypq = realloc(stats->num_Ypq, (p+1)*sizeof(long long*));
	for (j=stats->max_p+1; j<=p; j++) {
	  stats->num_Ypq[j] = calloc(j, sizeof(long long));
	}
	stats->max_p = p;
      }

      if (p > max_p)
	max_p = p;
      if (q > max_q)
	max_q = q;

      stats->num_Ypq[p][q]++;
    }

    if (faces[i].is_3dP0s) {
      has_dP0_cube = 1;
      stats->num_3dP0s++;
    }

    if (faces[i].is_5dP0s) {
      has_dP0_5 = 1;
      stats->num_5dP0s++;
    }

    if (faces[i].is_Z2xZ2_orbifold) {
      stats->num_Z2xZ2++;
    }

    if (faces[i].is_Laba) {
      int a = faces[i].Laba_a;
      int b = faces[i].Laba_b;

      if (b > stats->max_Laba_b) {
	stats->num_Laba = realloc(stats->num_Laba, (b+1)*sizeof(long long*));
	for (j=stats->max_Laba_b+1; j<=b; j++) {
	  stats->num_Laba[j] = calloc(j+1, sizeof(long long));
	}
	stats->max_Laba_b = b;
      }

      if (a == 12 && b == 12)
	has_L121212 = 1;

      stats->num_Laba[b][a]++;
    }

    if (faces[i].is_N2_orbifold) {
      int n = faces[i].N2_n;

      if (n > stats->max_N2_n) {
	stats->num_N2_n = realloc(stats->num_N2_n, (n+1)*sizeof(long long));
	memset(stats->num_N2_n+(stats->max_N2_n+1), 0,
	       (n - stats->max_N2_n)*sizeof(long long));
	stats->max_N2_n = n;
      }

      stats->num_N2_n[n]++;

      if (n > 30) {
	if (n > max_N2_n) {
	  has_N2_n = realloc(has_N2_n, (n+1)*sizeof(int));
	  memset(has_N2_n+(max_N2_n+1), 0, (n - max_N2_n)*sizeof(int));
	  max_N2_n = n;
	}
	has_N2_n[n] = 1;
      }
    }

    if (faces[i].is_dP0dP0dP1) {
      has_dP0dP0dP1 = 1;
      stats -> num_dP0dP0dP1++;
    }
  }

  /* Save examples */
  if (*max_face_points > 200) {
    int npolys = stats->nvery_singular;

    if (!poly)
      poly = new_poly(cols, vertices, info->poly_id);
    Py_INCREF(poly);

    stats->very_singular_polys = realloc(stats->very_singular_polys,
					 (npolys+1)*sizeof(stats->very_singular_polys[0]));
    stats->very_singular_npoints = realloc(stats->very_singular_npoints,
					   (npolys+1)*sizeof(stats->very_singular_npoints[0]));

    stats->very_singular_polys[npolys] = poly;
    stats->very_singular_npoints[npolys] = *max_face_points;
    stats->nvery_singular++;
  }

  if (max_p >= 8) {
    PyObject *key;
    PyObject *list;

    if (!poly)
      poly = new_poly(cols, vertices, info->poly_id);

    key = Py_BuildValue("(ii)", max_p, max_q);
    
    if (stats->Ypq_examples == NULL)
      stats->Ypq_examples = PyDict_New();
    
    list = PyDict_GetItem(stats->Ypq_examples, key);
    if (!list)
      list = PyList_New(0);
    else
      Py_INCREF(list);

    PyList_Append(list, poly);

    PyDict_SetItem(stats->Ypq_examples, key, list);
    Py_DECREF(list);

    Py_DECREF(key);
  }

  if (has_conifold && stats->conifold_example == NULL) {
    if (!poly)
      poly = new_poly(cols, vertices, info->poly_id);
    Py_INCREF(poly);

    stats->conifold_example = poly;
  }

  if (has_L121212 && stats->Laba_example == NULL) {
    if (!poly)
      poly = new_poly(cols, vertices, info->poly_id);
    Py_INCREF(poly);

    stats->Laba_example = poly;
  }

  if (has_dP0dP0dP1) {
    if (!poly)
      poly = new_poly(cols, vertices, info->poly_id);

    if (stats->dP0dP0dP1_examples == NULL)
      stats->dP0dP0dP1_examples = PyList_New(0);

    if (PyList_Size(stats->dP0dP0dP1_examples) < 1000) {
      PyList_Append(stats->dP0dP0dP1_examples, poly);
    }
  }

  if (has_dP0_cube) {
    if (!poly)
      poly = new_poly(cols, vertices, info->poly_id);

    if (stats->dP0_cube_examples == NULL)
      stats->dP0_cube_examples = PyList_New(0);

    if (PyList_Size(stats->dP0_cube_examples) < 1000) {
      PyList_Append(stats->dP0_cube_examples, poly);
    }
  }

  if (has_dP0_5) {
    if (!poly)
      poly = new_poly(cols, vertices, info->poly_id);

    if (stats->dP0_5_examples == NULL)
      stats->dP0_5_examples = PyList_New(0);

    PyList_Append(stats->dP0_5_examples, poly);
  }

  for (i=0; i<4; i++) {
    if (has_dP_n[i] && stats->dP_n_example[i] == NULL) {
      if (!poly)
	poly = new_poly(cols, vertices, info->poly_id);
      
      Py_INCREF(poly);
      stats->dP_n_example[i] = poly;
    }
  }

  for (i=0; i<=max_N2_n; i++) {
    if (has_N2_n[i]) {
      PyObject *key = PyInt_FromLong(i);
      PyObject *list;

      if (!poly)
	poly = new_poly(cols, vertices, info->poly_id);

      if (stats->N2_examples == NULL) {
	stats->N2_examples = PyDict_New();
      }

      list = PyDict_GetItem(stats->N2_examples, key);
      if (list == NULL)
	list = PyList_New(0);
      else
	Py_INCREF(list);

      PyList_Append(list, poly);

      PyDict_SetItem(stats->N2_examples, key, list);

      Py_DECREF(list);

      Py_DECREF(key);
    }
  }

  Py_XDECREF(poly);

  free(has_N2_n);

  /* Free the memory */
  for (i = 0; i<nfaces; i++) {
    free(faces[i].vertices);
    free(faces[i].edge_ninterior);
  }

  free(faces);
  free(edges);
  free(points);
}

static PyObject *pyint_fromlonglong(long long n)
{
  if (n <= LONG_MAX  && n >= LONG_MIN)
    return PyInt_FromLong(n);

  return PyLong_FromLongLong(n);
}

static PyObject *wrap_results(stats_t *stats,
			      int extra_singular_polys[][2],
			      int nextra, threadinfo *info)
{
  PyObject *dict = PyDict_New();
  PyObject *histogram = PyTuple_New(stats->max_interior+1);
  PyObject *extra_singus = PyTuple_New(nextra);
  PyObject *num_dP_n = PyTuple_New(4);
  PyObject *num_Ypq = PyTuple_New(stats->max_p+1);
  PyObject *num_3dP0s = pyint_fromlonglong(stats->num_3dP0s);
  PyObject *num_5dP0s = pyint_fromlonglong(stats->num_5dP0s);
  PyObject *num_dP0dP0dP1 = pyint_fromlonglong(stats->num_dP0dP0dP1);
  PyObject *num_conifold = pyint_fromlonglong(stats->num_conifold);
  PyObject *num_Z2xZ2 = pyint_fromlonglong(stats->num_Z2xZ2);
  PyObject *num_Laba = PyDict_New();
  PyObject *num_N2_n = PyTuple_New(stats->max_N2_n+1);
  PyObject *scan_range = Py_BuildValue("(ii)", info->start,
				       info->poly_id-1);
  PyObject *very_singular = PyTuple_New(stats->nvery_singular);
  PyObject *dP_n_example = PyTuple_New(4);
  int i;

  for (i=0; i<=stats->max_interior; i++) {
    PyTuple_SetItem(histogram, i, pyint_fromlonglong(stats->face_count[i]));
  }

  PyDict_SetItemString(dict, "Interior Point Histogram", histogram);
  Py_DECREF(histogram);

  for (i=0; i<nextra; i++) {
    PyTuple_SetItem(extra_singus, i,
		    Py_BuildValue("(ii)",
				  extra_singular_polys[i][0],
				  extra_singular_polys[i][1]));
  }
  PyDict_SetItemString(dict, "Polytopes with very singular faces",
		       extra_singus);
  Py_DECREF(extra_singus);

  for (i=0; i<4; i++) {
    PyTuple_SetItem(num_dP_n, i,
		    pyint_fromlonglong(stats->num_dP_n[i]));
  }
  PyDict_SetItemString(dict, "Number of del Pezzos", num_dP_n);
  Py_DECREF(num_dP_n);

  for (i=0; i<=stats->max_p; i++) {
    PyObject *nYp = PyTuple_New(i);
    int j;
    for (j=0; j<i; j++) {
      PyTuple_SetItem(nYp, j, pyint_fromlonglong(stats->num_Ypq[i][j]));
    }

    PyTuple_SetItem(num_Ypq, i, nYp);
  }
  PyDict_SetItemString(dict, "Number of Ypq", num_Ypq);
  Py_DECREF(num_Ypq);

  PyDict_SetItemString(dict, "Number of dP0^3 models", num_3dP0s);
  Py_DECREF(num_3dP0s);

  PyDict_SetItemString(dict, "Number of dP0^5 models", num_5dP0s);
  Py_DECREF(num_5dP0s);

  PyDict_SetItemString(dict, "Number of dP0^2 x dP1 models", num_dP0dP0dP1);
  Py_DECREF(num_dP0dP0dP1);

  PyDict_SetItemString(dict, "Number of conifolds", num_conifold);
  Py_DECREF(num_conifold);

  PyDict_SetItemString(dict, "Number of C^3/(Z2xZ2) orbifolds", num_Z2xZ2);
  Py_DECREF(num_Z2xZ2);

  for (i=0; i<=stats->max_N2_n; i++)
    PyTuple_SetItem(num_N2_n, i, pyint_fromlonglong(stats->num_N2_n[i]));
  PyDict_SetItemString(dict, "Number of CxC^2/Z_n orbifolds", num_N2_n);
  Py_DECREF(num_N2_n);

  for (i=0; i<=stats->max_Laba_b; i++) {
    int j;

    for (j=0; j<=i; j++) {
      if (stats->num_Laba[i][j] != 0) {
	PyObject *key = Py_BuildValue("(ii)", j, i);
	PyObject *value = pyint_fromlonglong(stats->num_Laba[i][j]);

	PyDict_SetItem(num_Laba, key, value);

	Py_DECREF(key);
	Py_DECREF(value);
      }
    }
  }
  PyDict_SetItemString(dict, "Number of L^{aba} cones", num_Laba);
  Py_DECREF(num_Laba);

  PyDict_SetItemString(dict, "Scan range", scan_range);
  Py_DECREF(scan_range);

  for (i=0; i<stats->nvery_singular; i++) {
    /* Here we steal the reference for the poly */
    PyTuple_SetItem(very_singular, i,
		   Py_BuildValue("(iO)",
				 stats->very_singular_npoints[i],
				 stats->very_singular_polys[i]));
  }
  PyDict_SetItemString(dict, "Examples of very singular polytopes",
		       very_singular);
  Py_DECREF(very_singular);

  if (stats->Ypq_examples) {
    PyDict_SetItemString(dict, "Examples of Ypq with p>7",
			 stats->Ypq_examples);
    Py_DECREF(stats->Ypq_examples);
  }

  if (stats->conifold_example) {
    PyDict_SetItemString(dict, "Conifold example", stats->conifold_example);
    Py_DECREF(stats->conifold_example);
  }

  if (stats->Laba_example) {
    PyDict_SetItemString(dict, "L^{12,12,12} example", stats->Laba_example);
    Py_DECREF(stats->Laba_example);
  }

  if (stats->dP0dP0dP1_examples) {
    PyDict_SetItemString(dict, "Examples of dP0^2 x dP1 models",
			 stats->dP0dP0dP1_examples);
    Py_DECREF(stats->dP0dP0dP1_examples);
  }

  if (stats->dP0_cube_examples) {
    PyDict_SetItemString(dict, "Examples of dP0^3 models",
			 stats->dP0_cube_examples);
    Py_DECREF(stats->dP0_cube_examples);
  }

  if (stats->dP0_5_examples) {
    PyDict_SetItemString(dict, "dP0^5 models", stats->dP0_5_examples);
    Py_DECREF(stats->dP0_5_examples);
  }

  for (i=0; i<4; i++) {
    if (stats->dP_n_example[i] != NULL)
      PyTuple_SetItem(dP_n_example, i, stats->dP_n_example[i]);
    else {
      Py_INCREF(Py_None);
      PyTuple_SetItem(dP_n_example, i, Py_None);
    }
  }
  PyDict_SetItemString(dict, "dP_n examples", dP_n_example);
  Py_DECREF(dP_n_example);

  if (stats->N2_examples) {
    PyDict_SetItemString(dict, "Examples of CxC^2/Z_n with n>30",
			 stats->N2_examples);
    Py_DECREF(stats->N2_examples);
  }

  return dict;
}

/* Analyze the given sequence of polytopes, and return a summary dictionary. */
static PyObject* read_polytopes(threadinfo *info)
{
  char *line = NULL;
  size_t nline = 0;

  PyObject *result;

  int (*extra_singular_polys)[2] = NULL;
  int nextra = 0;
  int max_face_points;

  int i;

  FILE *poly_x_in;
  FILE *poly_x_out;

  stats_t stats = {
    .max_interior = -1,
    .face_count = NULL,
    .num_conifold = 0,
    .num_dP_n = {0,0,0,0},
    .max_p = -1,
    .num_Ypq = NULL,
    .num_3dP0s = 0,
    .num_5dP0s = 0,
    .num_dP0dP0dP1 = 0,
    .max_Laba_b = -1,
    .num_Laba = NULL,
    .num_Z2xZ2 = 0,
    .max_N2_n = -1,
    .num_N2_n = NULL,
    .very_singular_polys = NULL,
    .very_singular_npoints = NULL,
    .nvery_singular = 0,
    .Ypq_examples = NULL,
    .conifold_example = NULL,
    .Laba_example = NULL,
    .dP0dP0dP1_examples = NULL,
    .dP0_cube_examples = NULL,
    .dP0_5_examples = NULL,
    .dP_n_example = {NULL, NULL, NULL, NULL},
    .N2_examples = NULL
  };

  info->poly_id = 0;

  poly_x(&poly_x_in, &poly_x_out, info->poly_x_cmd, info->poly_cb == NULL);

  while (1) {
    int rows, cols;
    int transposed = 0;

    if (getline(&line, &nline, info->fd) < 0) {
      break;
    }

    sscanf(line, "%d %d", &rows, &cols);

    if (rows > cols) {
      int c = rows;
      rows = cols;
      cols = c;

      transposed = 1;
    }

    assert(rows == 4);

    max_face_points = 0;

    read_one_polytope(rows, cols, poly_x_in, poly_x_out,
		      info, transposed, &max_face_points, &stats);

    if (info->poly_id >= info->start
	&& max_face_points > 200) {
      extra_singular_polys = realloc(extra_singular_polys,
				     (nextra+1)*sizeof(int[2]));
      extra_singular_polys[nextra][0] = info->poly_id;
      extra_singular_polys[nextra][1] = max_face_points;

      nextra++;
    }

    info->poly_id++;

    if (info->poly_id % 10000000 == 0)
      printf("%d\n", info->poly_id);

    if (info->end >= 0 && info->poly_id > info->end)
      break;
  }

  free(line);

  fclose(poly_x_in);
  fclose(poly_x_out);

  result = wrap_results(&stats, extra_singular_polys, nextra, info);

  free(stats.face_count);
  free(extra_singular_polys);

  for (i=0; i<=stats.max_p; i++)
    free(stats.num_Ypq[i]);
  free(stats.num_Ypq);
  free(stats.num_N2_n);
  for (i=0; i<=stats.max_Laba_b; i++)
    free(stats.num_Laba[i]);
  free(stats.num_Laba);

  free(stats.very_singular_polys);
  free(stats.very_singular_npoints);

  return result;
}

/* Construct the reference faces */
static void
construct_reference_models(threadinfo *info)
{
  info -> dP0_cube = PyTuple_New(4);
  PyTuple_SetItem(info->dP0_cube, 0, Py_BuildValue("(ii)", 0, 0));
  PyTuple_SetItem(info->dP0_cube, 1, Py_BuildValue("(ii)", -2, 2));
  PyTuple_SetItem(info->dP0_cube, 2, Py_BuildValue("(ii)", -2, -1));
  PyTuple_SetItem(info->dP0_cube, 3, Py_BuildValue("(ii)", -3, 0));

  info -> dP0_dP0_dP1 = PyTuple_New(5);
  PyTuple_SetItem(info->dP0_dP0_dP1, 0, Py_BuildValue("(ii)", 0, 0));
  PyTuple_SetItem(info->dP0_dP0_dP1, 1, Py_BuildValue("(ii)", -2, 2));
  PyTuple_SetItem(info->dP0_dP0_dP1, 2, Py_BuildValue("(ii)", -2, -1));
  PyTuple_SetItem(info->dP0_dP0_dP1, 3, Py_BuildValue("(ii)", -3, 0));
  PyTuple_SetItem(info->dP0_dP0_dP1, 4, Py_BuildValue("(ii)", -3, 1));
}

static PyObject *
scan_begin_scan(PyObject *self, PyObject *args, PyObject *keywords)
{
  const char *class_x_cmd;
  const char *zzdb;
  threadinfo info = {
    .poly_cb = NULL,
    .face_cb = NULL,
    .three_sector_cb = NULL
  };
  PyObject *result;
  static char *kwlist[] = {"start", "end",
			   "class_x_cmd", "zzdb", "poly_x_cmd",
			   "three_sector_cb", "face_cb", "poly_cb", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, keywords,
				   "iisss|OOO:begin_scan", kwlist,
				   &info.start, &info.end,
				   &class_x_cmd, &zzdb, &info.poly_x_cmd,
				   &info.three_sector_cb,
				   &info.face_cb, &info.poly_cb))
    return NULL;

  if (info.three_sector_cb != NULL &&
      !PyCallable_Check(info.three_sector_cb)) {
    PyErr_SetString(PyExc_TypeError,
		    "three_sector_cb must be callable.");
    return NULL;
  }

  if (info.face_cb != NULL &&
      !PyCallable_Check(info.face_cb)) {
    PyErr_SetString(PyExc_TypeError,
		    "face_cb must be callable.");
    return NULL;
  }

  if (info.poly_cb != NULL &&
      !PyCallable_Check(info.poly_cb)) {
    PyErr_SetString(PyExc_TypeError,
		    "poly_cb must be callable.");
    return NULL;
  }

  info.fd = class_x(class_x_cmd, zzdb, &info.class_x_pid);

  Py_XINCREF(info.face_cb);
  Py_XINCREF(info.poly_cb);

  if (info.three_sector_cb)
    construct_reference_models(&info);

  result = read_polytopes(&info);

  if (kill(info.class_x_pid, SIGTERM) < 0) {
    perror("kill");
  }
  fclose(info.fd);

  Py_XDECREF(info.poly_cb);
  Py_XDECREF(info.face_cb);

  if (info.three_sector_cb) {
    Py_DECREF(info.dP0_cube);
    Py_DECREF(info.dP0_dP0_dP1);
  }

  return result;
}

static PyMethodDef ScanMethods[] = {
  {"begin_scan", (PyCFunction)scan_begin_scan,
   METH_VARARGS | METH_KEYWORDS,
   "begin_scan(start, end, class_x_cmd, zzdb, poly_x_cmd\n"
   "\t\t[, three_sector_cb, face_cb, poly_cb])\n\n"
   "Starts the scan of 4d reflexive polytopes with id in [start,end].\n"
   "Returns the histogram of interior points in the faces of the scanned\n"
   "polytopes.\n\n"
   "start:\tid of the first polytope to scan in zzdb.\n"
   "end:\tid of the last polytope to scan in zzdb, can be negative\n"
   "\tto indicate scan until the end.\n"
   "class_x_cmd:\tPath to the appropriate polytope generator command.\n"
   "zzdb:\tLocation of the polytope database to pass to class_x_cmd.\n"
   "poly_x_cmd:\tPath to the poly.x executable.\n\n"
   "If given, face_cb and poly_cb will be called for each face/polytope,\n"
   "and they should return a list of points in the polytope.\n\n"
   "three_sector_cb(2d_model, 4d_vertices) will be called for each\n"
   "candidate in the list of target models. It should return whether\n"
   "the given 4d vertices form an embedding of the given 2d_model."
  },
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initscan(void)
{
  PyObject *m;

  if (sizeof(long long) < 8) {
    printf("WARNING: sizeof(long long) is not big enough,\n"
	   "some results may wrap around!\n");
  }

  m = Py_InitModule3("scan", ScanMethods,
		     "A module for scanning polytopes in PALP output");
}
