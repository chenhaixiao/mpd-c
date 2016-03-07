#ifndef _mpp_h
#define _mpp_h
#define _GNU_SOURCE // for asprintf from stdio.h
#define TRUE 1
#define FALSE 0
#define minim(a,b) ((a<b)?a:b)
#define maxim(a,b) ((a>b)?a:b)
#endif

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "dbg.h"

typedef struct primer_node
{
  int start;
  int end;
  char *sequence;
  double tm;
  double gc;
} PRIMER;

typedef struct primer_pair
{
  int id;
  PRIMER *forward;
  PRIMER *reverse;
  char *sequence;
  int length;
  double gc;
  double tm;
  int chrom;
} PNODE;

void die ( char *message );
int check_hairpin (char *ss, int n);
int check_hairpin_min (char *ss, char *flip, int n, int min);
int check_dimer (char *p1, char *p2, int n);
int check_uneven_dimer (char *p1, char *p2, int n1, int n2);
int check_watson_crick (char a, char b);
int test_dimer (char *p1, char *p2, int n);
int read_primer_pools (const char *filename, int max_ppairs, int max_primer_count,
  int *pool_count, PNODE ***primer_pool);
void Print_isPcr (const char *filename, int max_pools, int *primers_in_pool, PNODE ***primer_pool);
void reverse_string (char *contig, char *s, int n);
int is_poolable_primer (PNODE * p1, PNODE * p2);
void check_poolability (PNODE ** primer_pool, int N, int pool_number);
void Check_all_pools ( int max_pools, int *pool_count, PNODE ***primer_pool );
void print_cmat (char **cmat, int N);
PNODE ***primer_pool_create ( int max_ppairs, int max_ppairs_count);
PRIMER *create_primer (int max_primer_length);
PNODE *create_ppair (int max_primer_length, int max_amplicon_length);

int *create_ivec (int row);
void free_ivec (int *v, int row);
char *create_cvec (int row);
void free_cvec (int *v, int row);
int **create_imat (int row, int col);
void free_imat (int **m, int row, int col);
char **create_cmat (int row, int col);
void free_cmat (char **m, int row, int col);

char *cvector (int nl, int nh);
int *ivector (int nl, int nh);
char **cmatrix (int nrl, int nrh, int ncl, int nch);
void free_ivector (int *v, int nl, int nh);
void free_cmatrix (char **m, int nrl, int nrh, int ncl, int nch);
void read_var (char *line, char *result);
