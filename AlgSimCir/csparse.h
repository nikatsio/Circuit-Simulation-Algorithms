#ifndef SPARSE_MATRIX_H_
#define SPARSE_MATRIX_H_

#include <stdlib.h>
#include <stdio.h>

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(w,j) (w [j] < 0)
#define CS_MARK(w,j) { w [j] = CS_FLIP (w [j]) ; }
#define CS_CSC(A) (A && (A->nz == -1))
#define CS_TRIPLET(A) (A && (A->nz >= 0))
#define HEAD(k,j) (ata ? head [k] : j)
#define NEXT(J)   (ata ? next [J] : -1)

/********************************************************************************
 *                                                                              *
 *                       DATA STRUCTURES DEFINITIONS                            *
 *                                                                              *
 ********************************************************************************/

typedef struct cs_sparse /* matrix in compressed-column or triplet form */ {
    int nzmax; /* maximum number of entries */
    int m; /* number of rows */
    int n; /* number of columns */
    int *p; /* column pointers (size n+1) or col indices (size nzmax) */
    int *i; /* row indices, size nzmax */
    double *x; /* numerical values, size nzmax */
    int nz; /* # of entries in triplet matrix, -1 for compressed-col */
} cs;

typedef struct cs_symbolic /* symbolic Cholesky, LU, or QR analysis */ {
    int *pinv; /* inverse row perm. for QR, fill red. perm for Chol */
    int *q; /* fill-reducing column permutation for LU and QR */
    int *parent; /* elimination tree for Cholesky and QR */
    int *cp; /* column pointers for Cholesky, row counts for QR */
    int *leftmost; /* leftmost[i] = min(find(A(i,:))), for QR */
    int m2; /* # of rows for QR, after adding fictitious rows */
    double lnz; /* # entries in L for LU or Cholesky; in V for QR */
    double unz; /* # entries in U for LU; in R for QR */
} css;

typedef struct cs_numeric /* numeric Cholesky, LU, or QR factorization */ {
    cs *L; /* L for LU and Cholesky, V for QR */
    cs *U; /* U for LU, R for QR, not used for Cholesky */
    int *pinv; /* partial pivoting for LU */
    double *B; /* beta [0..n-1] for QR */
} csn;


/********************************************************************************
 *                                                                              *
 *                            FUNCTION DECLARATIONS                             *
 *                                                                              *
 ********************************************************************************/


/**
 *  Wrapper for malloc() function. It is used to allocate at least memory space equal to size.
 *  @param n The number of objects.
 *  @param size The size of each object.
 *  @return Pointer to the allocated space or NULL in case of failure.
 */
void *cs_malloc(int n, size_t size);


/**
 *  Wrapper for calloc() function. It is used to allocate and clear at least memory space equal to size.
 *  @param n The number of objects.
 *  @param size The size of each object.
 *  @return Pointer to the allocated space or NULL in case of failure.
 */
void *cs_calloc(int n, size_t size);


/**
 *  Wrapper for free() function. It is used to deallocate a previously allocated memory space.
 *  @param p Pointer to the allocated memory.
 *  @return NULL in order to simplify the use of cs_free().
 */
void *cs_free(void *p);


/**
 *  Wrapper for realloc() function.
 *  @param p Pointer to a previously allocated memory space.
 *  @param size The new size of the memory space.
 *  @param ok Pointer to a integer used to denote success or failure.
 *  @return Pointer to the newly allocate memory space in case of success or pointer to the original memory space otherwise.
 */
void *cs_realloc(void *p, int n, size_t size, int *ok);


/**
 *  Function for deallocating the allocated memory space for a sparse matrix in the Compressed Column format.
 *  @param A Pointer to the matrix.
 *  @return NULL.
 */
cs *cs_spfree(cs *A);


/**
 *  Function for deallocating the allocated memory space for a matrix numerical factorization.
 *  @param N Pointer to the struct describing the matrix factorization.
 *  @return NULL.
 */
csn *cs_nfree(csn *N);


/**
 *  Function for deallocating the allocated memory space for a matrix symbolic factorization.
 *  @param S Pointer to the struct describing the matrix factorization.
 *  @return NULL.
 */
css *cs_sfree(css *S);


/**
 *  Function for deallocating the internally allocated workspace and returning a sparse matrix result.
 *  @param C Sparse matrix result.
 *  @param w Workspace to free.
 *  @param x Workspace to free.
 *  @param ok Integer denoting whether to free (ok = 0) or keep sparse matrix (ok = 1).
 *  @return C in case of success or NULL otherwise.
 */
cs *cs_done(cs *C, void *w, void *x, int ok);


/**
 *  Function for deallocating the internally allocated workspace and returning a int matrix..
 *  @param p Int array.
 *  @param C Temporary sparse matrix to free.
 *  @param w Workspace to free.
 *  @param ok Integer denoting whether to free (ok = 0) or keep int matrix (ok = 1).
 *  @return p in case of success or NULL otherwise.
 */
int *cs_idone(int *p, cs *C, void *w, int ok);


/**
 *  Function for deallocating the internally allocated workspace and returning a numeric factorization result.
 *  @param N Numeric factorization result.
 *  @param C Temporary sparse matrix to free.
 *  @param w Workspace to free.
 *  @param x Workspace to free.
 *  @param ok Integer denoting whether to free (ok = 0) or keep numeric factorization (ok = 1).
 *  @return N in case of success or NULL otherwise.
 */
csn *cs_ndone(csn *N, cs *C, void *w, void *x, int ok);


/**
 *  Function for allocating the appropriate memory space for a sparse matrix in triplet or compressed-column format.
 *  @param m Number of rows.
 *  @param n Number of columns.
 *  @param nzmax Number of maximum number of non-zero elements.
 *  @param values Flag that is used to denote whether only pattern (values = 0) or both pattern and values (value = 1) will be allocated.
 *  @param triplet Flag that denotes whether the matrix will be stored in the compressed-column (triplet = 0) or triplet format (triplet = 1).
 *  @return Pointer to the struct describing the compressed matrix in case of success and NULL otherwise.
 */
cs *cs_spalloc(int m, int n, int nzmax, int values, int triplet);


/**
 *  Function for changing the maximun number of entries a sparse matrix can store.
 *  @param A Pointer to the struct describing the sparse matrix.
 *  @param nzmax New number of maximum entries.
 *  @return 1 if modification is successful and 0 in case of failure.
 */
int cs_sprealloc(cs *A, int nzmax);


/**
 *  Function for converting a matrix from triplet to compressed-column format. The columns of new matrix
 *  are not sorted and duplicate entries may be present.
 *  @param T Sparse matrix in triplet format.
 *  @return The sparse matrix in compressed-column format or NULL on error.
 */
cs *cs_compress(const cs *T);


/**
 *  Function for computing the cumulative sum of an integer vector.
 *  @param p The cumulative sum of the integer vector.
 *  @param c The input integer vector. It is overwritten with the elements p[0 ... n-1] when function returns.
 *  @param n The length of vector c. Vector p has size n+1.
 *  @return Function returns sum(c) or 0 in case of an error.
 */
double cs_cumsum(int *p, int *c, int n);


/**
 *  Function that is used for computing the tranpose of a sparse matrix.
 *  @param A The sparse matrix.
 *  @param values Flag that denotes whether only pattern (values = 0) or both pattern and values (value = 1) will be transposed.
 *  @return The tranpose matrix or NULL in case of error.
 */
cs *cs_transpose(const cs *A, int values);


/**
 *  Function that removes and sums duplicate entries in a sparse matrix.
 *  @param A The sparse matrix.
 *  @return 1 if successful and 0 in case of failure.
 */
int cs_dupl(cs *A);


/**
 *  Function that computes the permutation x = Pb of a vector.
 *  @param p The permutation vector. If p==NULL then the permutation vector is the identity vector.
 *  @param b Input vector.
 *  @param x Output vector.
 *  @param n Vector length.
 *  @return 1 if successful and 0 in case of error.
 */
int cs_pvec(const int *p, const double *b, double *x, int n);


/**
 *  Function that computes the permutation x = P'b of a vector.
 *  @param p The permutation vector. If p==NULL then the permutation vector is the identity vector.
 *  @param b Input vector.
 *  @param x Output vector.
 *  @param n Vector length.
 *  @return 1 if successful and 0 in case of error.
 */
int cs_ipvec(const int *p, const double *b, double *x, int n);


/**
 *  Function that inverts a permutation vector.
 *  @param p The permutation vector.
 *  @param n Vector length.
 *  @return The inverted permutation or NULL on error.
 */
int *cs_pinv(int const *p, int n);


/**
 *  Function that computes the symmetric permutation C = PAP' for a symmetric matrix A.
 *  @param A Compressed matrix to permute. Only the upper triangular part is used.
 *  @param pinv Inverse permutation vector.
 *  @param values Allocate pattern only if values = 0 and values and pattern otherwise.
 *  @return The symmetric permutation or NULL on error.
 */
cs *cs_symperm(const cs *A, const int *pinv, int values);


/**
 *  Function that scatters and sums a sparse vector A(:,j) into a dense vector, x = x + beta*A(:,j).
 *  @param A The sparse vector A(:,j).
 *  @param j The column of A to use.
 *  @param beta Scalar multiplied by A(:,j).
 *  @param w Auxiliary vector that stores the marked elements of A.
 *  @param x The final vector. It is ignored if it is NULL.
 *  @param mark Mark value for vector w.
 *  @param C Pattern of x accumulated in C->i.
 *  @param nz Pattern of x placed in C starting at C->i[nz].
 *  @return New value of nz or -1 on error.
 */
int cs_scatter(const cs *A, int j, double beta, int *w, double *x, int mark, cs *C, int nz);

/*
 *  Function that performs a sparse matrix permutation.
 *  @param A Sparse matrix to permute.
 *  @param pinv permutation vector of length m.
 *  @param q permutation vector of length n.
 *  @param values allocate pattern only if 0, values and pattern otherwise.
 *  @return C=A(p,q) or NULL on error.
 */
cs *cs_permute(const cs *A, const int *pinv, const int *q, int values);


/**
 *  Function for sparse matrix addition C = alpha * A + beta * B.
 *  @param A The first matrix.
 *  @param B The second matrix.
 *  @param alpha Multiplication factor for matrix A.
 *  @param beta Multiplication factor for matrix B.
 *  @return New sparse matrix or NULL on error.
 */
cs *cs_add(const cs *A, const cs *B, double alpha, double beta);


/**
 *  Function for sparse matrix addition C = alpha * A + beta * B.
 *  @param A Multiplicand matrix.
 *  @param B Multiplier matrix.
 *  @return New sparse matrix or NULL on error.
 */
cs *cs_multiply(const cs *A, const cs *B);


/**
 * Function that implements the following equation: y = A*x + y.
 * @param A Multiplicand matrix.
 * @param x Multiplier Vector.
 * @param y Addition vector on input and the solution on output.
 * @return 1 if successful and 0 in case of error.
 */
int cs_gaxpy(const cs *A, const double *x, double *y);
double *cs_gaxp(const cs *A, const double *x);

/*
 *  1st norm of a matrix.
 *  @param A matrix.
 *  @return norm.
 */
double cs_norm(const cs *A);


/**
 *  Function for dropping entries from a sparse matrix. It drops element a[i][j] if fkeep(i, j, a[i][j], other) is zero.
 *  @param A Sparse matrix.
 *  @param fkeep Pointer to the fkeep function used for testing.
 *  @param other Optional parameter for fkeep function.
 *  @return The new number of entries in matrix A or NULL on error.
 */
int cs_fkeep(cs *A, int(*fkeep)(int, int, double, void *), void *other);

/*
 *  Function for performing LU decomposition with partial (row) pivoting of a sparse matrix.
 *  @param A Sparse matrix.
 *  @param S The symbolic analysis of matrix A, as it is computed from cs_sqr() function.
 *  @param tol partial pivoting threshold (1 for partial pivoting).
 *  @return The numerical analysis of matrix A or NULL on error.
 */
csn *cs_lu(const cs *A, const css *S, double tol);

/**
 *  Function for solving a sparse lower triangular system Lx = b.
 *  @param L The lower triangular matrix. Matrix must have a zero-free diagonal.
 *  @param x The right-hand side vector on input and the solution on output.
 *  @return 1 if successful and 0 in case of error.
 */
int cs_lsolve(const cs *L, double *x);


/*
 *  Function for solving a sparse upper triangular system Ux = b.
 *  @param U The upper triangular matrix. Matrix must have a zero-free diagonal.
 *  @param x The right-hand side vector on input and the solution on output.
 *  @return 1 if successful and 0 in case of error.
 */
int cs_usolve(const cs *U, double *x);


/**
 *  Function for solving a sparse upper triangular system L'x = b.
 *  @param L The lower triangular matrix. Matrix must have a zero-free diagonal.
 *  @param x The right-hand side vector on input and the solution on output.
 *  @return 1 if successful and 0 in case of error.
 */
int cs_ltsolve(const cs *L, double *x);


/*
 *  Function that solves x=A\b. where A is unsymmetric; b overwritten with solution
 *  @param order The ordering method that will be used (0:natural, 1:Chol, 2:LU, 3:QR).
 *  @param A Matrix to analyze.
 *  @param b The right-hand side vector on input and the solution on output.
 *  @param tol partial pivoting threshold (1 for partial pivoting).
 *  @return 1 if successful and 0 in case of error.
 */
int cs_lusol(int order, const cs *A, double *b, double tol);


/*
 *  Column counts for sparse LU factorization of A.
 * 	@param A Matrix to analyse.
 * 	@param S The symbolic analysis of matrix A.
 * 	@return 0 on success and 1 otherwise.
 */
static int cs_vcount(const cs *A, css *S);


/*
 *  Function that performs symbolic ordering and analysis for QR or LU.
 *  @param order The ordering method that will be used (0:natural, 1:Chol, 2:LU, 3:QR).
 *  @param A Matrix to analyze.
 *  @param qr Flag that denotes whether will be performed QR symbolic analysis (qr = 1) or not (qr = 0).
 *  @return The symbolic analysis of matrix A.
 */
css *cs_sqr(int order, const cs *A, int qr);


/**
 *  Function for computing the elimination tree of A or A'A, without forming A'A.
 *  @param A Matrix to analyze.
 *  @param ata Flag that denotes whether we need to analyze A (ata = 0) or A'A (ata = 1).
 *  @return Vector of size n with the elimination pattern of matrix (parent) or NULL on error.
 */
int *cs_etree(const cs *A, int ata);

/**
 *  Finds the nonzero pattern of kth row of Cholesky factor, L(k,1:k-1).
 *  @param G Graph to search.
 *  @param B Right-Hand side, b=B(:,k).
 *  @param k kth column of B.
 *  @param xi Output in x[top ... n-1], size 2*n.
 *  @param pinv Mapping of rows to columns of G, ignored if NULL.
 *  @return top or -1 on error.
 */
int cs_reach(cs *G, const cs *B, int k, int *xi, const int *pinv);

/**
 *  Finds the nonzero pattern of kth row of Cholesky factor, L(k,1:k-1).
 *  @param A L is the Cholesky factor of A.
 *  @param k The number of the row.
 *  @param parent The elimination tree of A.
 *  @param s Vector with the nonzero pattern of L(k,1:k-1).
 *  @param w Temporary vector that holds the mark value of each node.
 *  @return The position in vector s where the nonzero pattern L(k,:) starts and -1 on error.
 */
int cs_ereach(const cs *A, int k, const int *parent, int *s, int *w);

/*
 *  Depth-First-Search of the graph of a matrix.
 *  @param j The starting node of the traversal.
 *  @param G graph to search.
 *  @param top stack xi[top ... n-1] in use on input.
 *  @param xi stack containing nodes traversed, size n.
 *  @param pstack work vector, size n.
 *  @param pinv mapping of rows to columns of G, ignored in NULL.
 *  @return
 */
int cs_dfs(int j, cs *G, int top, int *xi, int *pstack, const int *pinv);


/**
 *  Postorder traversal of a tree.
 *  @param j The starting node of the traversal.
 *  @param k The number of nodes ordered so far.
 *  @param head On input head[i] stores the first child of node i and -1 on output.
 *  @param post Postordering.
 *  @param stack Temporary vector of size n.
 *  @return New value of k and -1 on error.
 */
int cs_tdfs(int j, int k, int *head, const int *next, int *post, int *stack);


/**
 *  Postorder traversal of a tree or forest.
 *  @param parent Defines the tree of n nodes.
 *  @param n Length of parent vector.
 *  @return Int array post where post[k] = i or NULL on error.
 */
int *cs_post(const int *parent, int n);


/**
 *  Function that determines whether j is a leaf and find least common ancestor.
 *  @param i The ith row subtree.
 *  @param j The number of leaf to be checked.
 *  @param first Vector with the first ancestor of each node.
 *  @param maxfirst The maximum ancestor seen so far.
 *  @param prevleaf Vector with the previous leaf of ith subtree.
 *  @param ancestor Stores the ancestors of the ith root subtree.
 *  @param jleaf Pointer to integer that stores whether this is the first or a subsequent leaf.
 *  @return The least common ancestor.
 */
int cs_leaf(int i, int j, const int *first, int *maxfirst, int *prevleaf, int *ancestor, int *jleaf);


/**
 *  Function that initializes the appropriate data structures that column_counts() needs in order to compute
 *  column counts for A'A matrix.
 *  @param AT The transpose of matrix A.
 *  @param post Postordering vector of parent.
 *  @param head The head of the linked list that is formed of the rows of A.
 *  @param next The next pointer of each node in the linked list.
 *  @return Nothing.
 */
void init_ata(cs *AT, const int *post, int *w, int **head, int **next);


/**
 *  Column counts for Cholesky factorization of A or A'A.
 *  @param A Matrix to analyze.
 *  @param parent Elimination tree of A.
 *  @param post Postordering of parent.
 *  @param ata Flag that denotes whether we need to analyze A (ata = 0) or A'A (ata = 1).
 *  @return A vector of length n with the column counts if operation is successful and NULL on error.
 */
int *cs_counts(const cs *A, const int *parent, const int *post, int ata);


/**
 *  Function that clears matrix w.
 *  @param mark Integer that denotes whether matrix w will be cleared (mark < 2). If matrix is cleared, mark is set to 2.
 *  @param lemax Integer that controls whether matrix w will be cleared.
 *  @param w Matrix w to be cleared.
 *  @param n The length of the matrix.
 *  @return The value of mark.
 */
int cs_wclear(int mark, int lemax, int *w, int n);


/**
 *  Function that drops diagonal entries. It is used as the fkeep parameter in cs_fkeep function.
 *  @param i The row of the element.
 *  @param j The column of the element.
 *  @param aij UNUSED.
 *  @param other UNUSED.
 *  @return 1 if i == j and 0 otherwise.
 */
int cs_diag(int i, int j, double aij, void *other);


/**
 *  Function that computes the approximate minimum degree ordering of A+A' or A'A.
 *  @param order The ordering method that will be used (0:natural, 1:Chol, 2:LU, 3:QR).
 *  @param A Matrix to order.
 *  @return The permutation of size n or NULL on error or if natural ordering is used.
 */
int *cs_amd(int order, const cs *A);


/**
 *  Function that computes a symbolic ordering and analysis for a Cholesky factorization.
 *  @param order The ordering option that will be subsequently used in cs_amd function.
 *  @param A Matrix to factorize.
 *  @return The symbolic analysis for cs_chol() function or NULL on error.
 */
css *cs_schol(int order, const cs *A);


/**
 *  Function that computes the sparse Cholesky factorization of a matrix.
 *  @param A Matrix to factorize.
 *  @param S The symbolic analysis of matrix A, as it is computed from cs_schol() function.
 *  @return The numerical analysis of matrix A or NULL on error.
 */
csn *cs_chol(const cs *A, const css *S);


/**
 *  Function that computes the refactorization of a matrix.
 *  @param A Matrix to factorize.
 *  @param N The numerical factorization of A.
 *  @param pinv The permutation vector.
 *  @param c Vector that stores the column pointers of A.
 *  @param x Permuted vector (the result of the invocation to cs_ipvec() function).
 *  @return The numerical analysis of matrix A or NULL on error.
 */
int cs_rechol(const cs *A, const csn *N, int *pinv, int *c, double *x);


/*
 *  Function that solves x=A\b. where A is symmetric and positive definite; b overwritten with solution
 *  @param order The ordering method that will be used (0:natural, 1:Chol, 2:LU, 3:QR).
 *  @param A Matrix to analyze.
 *  @param b The right-hand side vector on input and the solution on output.
 *  @return 1 if successful and 0 in case of error.
 */
int cs_cholsol(int order, const cs *A, double *b);


/**
 *  Function that drops matrix elements below a tolerance value. It is assumed that the matrix is created in such
 *  way that the first element Ax(Ap[j]) of each column is the diagonal element.
 *  @param A Matrix to analyze.
 *  @param tol Tolerance value.
 *  @return The new number of non-zero elements of A.
 */
int cs_reltol(cs *A, double tol);

/**
 *  Function that solves lower or upper triangular system.
 *  @param G is either upper U (lo=0) or lower L (lo=1) triangular.
 *  @param B Right-Hand side, b=B(:,k).
 *  @param k kth column of B.
 *  @param xi Output in xi[top ... n-1], size 2*n.
 *  @param x Output in x[xi[top ... n-1]], size n.
 *  @param pinv Mapping of rows to columns of L, ignored if NULL.
 *  @param lo is 0 for upper triangular and 1 for lower triangular.
 *  @return top or -1 on error.
 */
int cs_spsolve(cs *G, const cs *B, int k, int *xi, double *x, const int *pinv, int lo);


/********************************************************************************
 *                                                                              *
 *                            UTILITY FUNCTIONS                                 *
 *                                                                              *
 ********************************************************************************/

/**
 *  Utility function that inserts a new element in a compressed matrix.
 *  @param T Input matrix.
 *  @param i The row of the new element.
 *  @param j The column of the new element.
 *  @param x The value of the new element.
 *  @return 0 on success and 1 otherwise.
 */
int cs_entry(cs *T, int i, int j, double x);


/**
 *  Utility function that is used to print a matrix in sparse format.
 *  @param A Matrix to print.
 *  @param outputFilePtr Output file name.
 *  @param brief If brief is equal to 1, only the first 20 non-zero elements of each column are printed.
 *  @return 0 on error and 1 otherwise.
 */
int cs_print(const cs *A, const char *outputFilename, int brief);


#endif /* SPARSE_MATRIX_H_ */