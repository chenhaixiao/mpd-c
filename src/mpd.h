#ifndef _mpp_h
#define _mpp_h
#define _GNU_SOURCE // for asprintf from stdio.h
#define TRUE 1
#define FALSE 0
#define MAX_PAIRS 10
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
#include "mem.h"

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
    PRIMER *forward;
    PRIMER *reverse;
    char *sequence;
    int length;
    double gc;
    double tm;
    int chrom;
} PNODE;

typedef struct snp_node
{
    char baseA;
    char baseB;
    char name[32];
    unsigned int pos;
    unsigned int chrom;
    int no_pairs;
    int no_disc;
    double het;
    PNODE **pair;
} SNODE;

typedef struct amp_node
{
    char name[32];
    unsigned int start_pos;
    unsigned int stop_pos;
    unsigned int chrom;
    int no_pairs;
    PNODE **pair;
} AMPNODE;

void die ( char *message );
int read_primer_pools (const char *filename, int max_ppairs, int max_primer_count,
                       int *pool_count, PNODE ***primer_pool);
void Print_isPcr (const char *filename, int max_pools, int *primers_in_pool, PNODE ***primer_pool);
void reverse_string (char *contig, char *s, int n);
void check_poolability (PNODE ** primer_pool, int N, int pool_number);
void Check_all_pools ( int max_pools, int *pool_count, PNODE ***primer_pool );

PNODE ***primer_pool_create ( int max_ppairs, int max_ppairs_count);
PRIMER *create_primer (int max_primer_length);
PNODE *create_ppair (int max_primer_length, int max_amplicon_length);
void flat_index_contig_high (int *index, char *contig, int L, int depth, int high_depth, unsigned char *double_high);
void flat_index_contig (int *index, char *contig, int L, int depth);
void high_index (unsigned char *dhigh, char *s, int n);
int
find_primers (SNODE ** snp_list, AMPNODE * tn, int no_snp, int *flat, char *contig,
              int L, int min_primer, int max_primer, int amp_max, int amp_min,
              double min_gc, double max_gc, double min_tm, double max_tm,
              int depth, int local_depth, int target_base, int start_pos,
              unsigned char *highmer, int **repeats, int no_repeats, int end_region, int chrom);
void count_copys (int *flat, char *s, int n, int *reps);
void simple_count_copys (int *flat, char *s, int n, int *reps, int forward);
void reverse_transcribe (char *contig, char *s, int n);
void transcribe (char *contig, char *s, int n);
void reverse_string (char *contig, char *s, int n);
int sort_compare_index (const void *a, const void *b);
int check_hairpin (char *ss, int n);
int check_hairpin_min (char *ss, char *flip, int n, int min);
int check_dimer (char *p1, char *p2, int n);
int check_uneven_dimer (char *p1, char *p2, int n1, int n2);
int check_watson_crick (char a, char b);
int test_dimer (char *p1, char *p2, int n);
int find_seg (int x, int y, int **list, int start, int end, int guess);
unsigned int encode_basepairs (char *ss, int n);
void decode_basepairs (unsigned char *s, char *dest, int n);
void index_string (int *index, char *s, int n);
void convert_int_basepairs (int i, char *s, int k);
void
fill_quality_scores (int *flat, int *local, char *contig, int L, int window, int depth, int ld, int min_primer,
                     int max_primer, double *fq_left, double *gc_left, int *index_left, int *plen_l,
                     double min_gc, double max_gc, double min_tm, double max_tm, int **repeats,
                     int no_repeats, unsigned char *highmer, int start_pos);
int find_frag (int x, int **list, int start, int end, int guess);
void fill_hs (char a, char b, double *h, double *s);
double calc_tm (char *s, int n);
int sort_compare_struct (const void *a, const void *b);
int check_15mer (char *string, unsigned char *map, int n);
int is_not_repeat (int x, int y, int **list, int no_frags, double max);
SNODE **fill_snp_list (FILE * sfile, int *n, unsigned int chrom);
AMPNODE **fill_amp_list (FILE * sfile, int n);
SNODE *snp_alloc (void);
AMPNODE *amp_alloc (void);
int isbase (char c);
double calc_gc (char *s, int n);
PRIMER *primer_alloc (void);
PNODE *product_alloc (void);
int poly_under_primer (unsigned int p_start, unsigned int p_end, SNODE ** list, int start, int stop, int which);
int *select_snps (SNODE ** list, int n, int target, int start, int stop);
int pick_random (double *x, int n);
int select_subset (SNODE ** list, int n, int pick, int *priority, int *selected, int which, int start, int stop);
double fill_dist (SNODE ** list, int n, int *priority, int *selected, int which, int start, int stop, double *dist);
int line_count (FILE * sfile);
int count_compatable_primers (int **poolable, int total_primers, int k);
int **zero_matrix (int **poolable, int total_primers, int k);
int is_poolable_amp (AMPNODE * a1, AMPNODE * a2, int i1, int i2);
void find_min_pools (int n, int *need_pooling, int **same_amplicon, int **poolable, PNODE ** primer_list);
int is_poolable_primer (PNODE * p1, PNODE * p2, int size_diff_threshold, int tm_diff_threshold);
int count_poolable (int t, int tot_primer, int *need_pooling, int **same_amplicon, int **poolable);
void make_greedy_pools (FILE *outfile, PNODE **plist, char **cmat, int *pc, int **redund, int *bstart, int Nregs, int N, int still, int *current_pool, int this_pool, int max_pool);
void make_less_greedy_pools (FILE *outfile, PNODE ** plist, char **cmat, int *pc, int **redund, int *bstart, int Nregs, int N, int still, int *current_pool, int this_pool, int max_pool);
void read_var (char *line, char *result);
int *create_ivec (int row);
char *create_cvec (int row);
int **create_imat (int row, int col);
char **create_cmat (int row, int col);
