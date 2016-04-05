#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <ctype.h>
#include "dbg.h"

#define TRUE 1
#define FALSE 0
#define MAX_PAIRS 10

#define minim(a,b) ((a<b)?a:b)
#define maxim(a,b) ((a>b)?a:b)

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

void read_var (char *line, char *result);
int *ivector (int nl, int nh);
char *cvector (int nl, int nh);
unsigned char *ucvector (int nl, int nh);
double *dvector (int nl, int nh);
double **dmatrix (int nrl, int nrh, int ncl, int nch);
char **cmatrix (int nrl, int nrh, int ncl, int nch);
unsigned char **ucmatrix (int nrl, int nrh, int ncl, int nch);
double ran1 (int *idum);
double gamdev (int ia, int *idum);
double gammln (double xx);
double poidev (double xm, int *idum);
double gs (double alpha);
double gcf (double alpha);
double rgamma (double alpha);
double rbeta (double u, double v);
int pick (int n);
void free_cvector (char *v, int nl, int nh);
void free_ucvector (unsigned char *v, int nl, int nh);
void free_ivector (int *v, int nl, int nh);
void free_dvector (double *v, int nl, int nh);
void free_dmatrix (double **m, int nrl, int nrh, int ncl, int nch);
void free_cmatrix (char **m, int nrl, int nrh, int ncl, int nch);
void free_ucmatrix (unsigned char **m, int nrl, int nrh, int ncl, int nch);
int **imatrix (int nrl, int nrh, int ncl, int nch);
void free_imatrix (int **m, int nrl, int nrh, int ncl, int nch);
void set_next (char *sss, int *i, int *j);
void flat_index_contig (int *index, char *contig, int L, int depth);
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
int is_poolable_primer (PNODE * p1, PNODE * p2);
int count_poolable (int t, int tot_primer, int *need_pooling, int **same_amplicon, int **poolable);
void make_pools (PNODE ** plist, char **cmat, int *pc, int **redund, int *bstart, int Nregs, int N, int still,
		  int *current_pool, int this_pool, int max_pool);

static FILE *outfile;
static int idum;



int
main ()
{
  char ss[256], sss[4196], **filename;
  char *scratch_pad, **contig_descript;
  unsigned char **compressed_map, *highmer;
  int i, j, k, N, N_snpfiles, *contig_snp_count, *contig_length;
  int pad_size, fasta, idepth, target_contig, total_index;
  int *flat_index, old_index, ***repeat_list, *no_repeats;
  int amp_min, amp_max, pool_size, N_targets;
  int min_primer, max_primer, genome_start, genome_stop;
  FILE *sfile, *cfile, *idfile, *rfile, *highfile, *snpfile_idx, *snpfile, *target_ampfile;
  double max_gc, min_gc, min_tm, max_tm, tm_inc;
  SNODE ***snp_list;
  AMPNODE **target_amp_list;

  outfile = stdout;
  read_var ("\nSend Output to Screen or Disk? [S,D]\n", ss);

  if ((strchr (ss, 'D')) || (strchr (ss, 'd')))
  {
    read_var ("Please Enter File Name for Output\n", ss);
    if ((outfile = fopen (ss, "w")) == (FILE *) NULL)
    {
      printf ("\n Can not open file %s\n", ss);
      exit (1);
    }
  }
  else
    outfile = stdout;

  old_index = TRUE;

  read_var ("Primer Picker Summary Filename (e.g., index summary like hg19.sdx)\n", ss);
  if ((sfile = fopen (ss, "r")) == (FILE *) NULL)
  {
    printf ("\nCould Not Open file %s\n", ss);
    exit (1);
  }

  // sdx file: read the 1st line that contains an int of contigs that are indexed
  fgets (sss, 4195, sfile);
  N = atoi (sss);
  printf ("\n There are N chromosomes %d \n\n", N);

  compressed_map = (unsigned char **) malloc ((unsigned) (N + 1) * sizeof (unsigned char *));
  if (!compressed_map)
    log_err ("allocation failure for compressed_map");

  contig_descript = cmatrix (0, N, 0, 4196);
  contig_length = ivector (0, N);
  no_repeats = ivector (0, N);

  repeat_list = (int ***) malloc ((unsigned) (N + 1) * sizeof (int **));
  if (!repeat_list)
    log_err ("allocation failure for repeat_list");


  // sdx file: read in contig information that is in the format: contig_lenght number_repeats contig_description
  for (i = 0; i < N; i++)
  {
    fgets (sss, 4195, sfile);
    sscanf (sss, "%d\t%d\t%s", &contig_length[i], &no_repeats[i], contig_descript[i]);
    repeat_list[i] = imatrix (0, no_repeats[i], 0, 1);
    printf ("\n Contig %d is named %s and is length %d\n", i, contig_descript[i], contig_length[i]);
  }

  // sdx file: read in int representing depth of coverage
  fgets (sss, 4195, sfile);
  idepth = atoi (sss);

  // printf("\n Indexing to a depth of %d \n",idepth);

  // sdx file: read next line - should be hg19.cdx
  fgets (sss, 4195, sfile);
  for (i = 0; i < 4194; i++)
    if (isspace (sss[i]))
    {
      sss[i] = '\0';
      i = 4195;
    }

  // cdx file: set file
  if ((cfile = fopen (sss, "r")) == (FILE *) NULL)
  {
    printf ("\nCould Not Open 1 file \"%s\"\n", sss);
    exit (1);
  }

  // sdx file: read next line - should be hg19.idx
  fgets (sss, 4195, sfile);
  for (i = 0; i < 4194; i++)
    if (isspace (sss[i]))
    {
      sss[i] = '\0';
      i = 4195;
    }

  // idx file: set file
  if ((idfile = fopen (sss, "r")) == (FILE *) NULL)
  {
    printf ("\nCould Not Open 2 file \"%s\"\n", sss);
    exit (1);
  }

  // sdx file: read next line - should be hg19.rdx
  fgets (sss, 4195, sfile);
  for (i = 0; i < 4194; i++)
    if (isspace (sss[i]))
    {
      sss[i] = '\0';
      i = 4195;
    }

  // rdx file: set file
  if ((rfile = fopen (sss, "r")) == (FILE *) NULL)
  {
    printf ("\nCould Not Open 3 file \"%s\"\n", sss);
    exit (1);
  }

  //sdx file: read next line - should be hg19.15x
  fgets (sss, 4195, sfile);
  for (i = 0; i < 4194; i++)
    if (isspace (sss[i]))
    {
      sss[i] = '\0';
      i = 4195;
    }

  // 'highfile' or 'hg19.15x': set file
  if ((highfile = fopen (sss, "r")) == (FILE *) NULL)
  {
    printf ("\nCould Not Open 4 file \"%s\"\n", sss);
    exit (1);
  }
  fclose (sfile);

  printf ("\nAbout to read repeat file\n\n");

  // index  Description i number of contigs: N j number of repeats in a contig: no_repeats[i] k not sure location of start of repeat until end of repeat
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < no_repeats[i]; j++)
      for (k = 0; k < 2; k++)
	      if ((fread (&repeat_list[i][j][k], sizeof (int), 1, rfile)) < 1)
      	{
      	  printf ("\n Expected to read %d repeats for contig %d, but got %d\n\n", no_repeats[i], i, j);
    	    exit (1);
      	}
    // for(j=0;j<1;j++) printf("\nFor contig %d repeat %d goes from %d %d\n\n", i,j,repeat_list[i][j][0],repeat_list[i][j][1]);
  }
  fclose (rfile);
  printf ("\nFinished reading repeat file");

  printf ("\nAbout to read 15mer file\n\n");
  highmer = ucvector (0, 134217728);
  if ((j = fread (highmer, sizeof (unsigned char), 134217728, highfile)) != 134217728)
  {
    printf ("\n Expected to read %d 15mers but got %d\n", 134217728, j);
    exit (1);
  }
  fclose (highfile);
  printf ("\nFinished reading 15mer file\n");

  //read dbSNP information -> summary file => actual files fill ***snp_list
  read_var ("Name of dbSNP summary file\n", ss);
  if ((snpfile_idx = fopen (ss, "r")) == (FILE *) NULL)
  {
    printf ("\nCould Not Open dbsnp file %s\n", ss);
    exit (1);
  }

  // dbsnp summary file line 1 = number of contigs/chr
  fgets (sss, 4195, snpfile_idx);
  N_snpfiles = atoi (sss);

  // read dbsnp files
  contig_snp_count = ivector (0, N_snpfiles);
  snp_list = (SNODE ***) malloc (N_snpfiles * sizeof (SNODE **));
  filename = cmatrix (0, N_snpfiles, 0, 256);
  for (i = 0; i < N_snpfiles; i++)
  {
    fgets (sss, 4195, snpfile_idx);
    sscanf (sss, "%s", filename[i]);
    if ((snpfile = fopen (filename[i], "r")) == (FILE *) NULL)
    {
      printf ("\n Can not open file %s\n", filename[i]);
      exit (1);
    }
    snp_list[i] = fill_snp_list (snpfile, &contig_snp_count[i], i + 1);
    fclose (snpfile);
    printf ("\nFinished reading dbSNP file: %s. Found %d SNPs.\n", filename[i], contig_snp_count[i]);
  }
  fclose(snpfile_idx);

  // read target file
  read_var ("\nName of file with amplicon target coordinates\n", ss);
  if ((target_ampfile = fopen (ss, "r")) == (FILE *) NULL)
  {
    printf ("\nCould Not Open file with amplicon target coordinates %s\n", ss);
    exit (1);
  }
  N_targets = line_count (target_ampfile);
  fseek (target_ampfile, 0, SEEK_SET);
  target_amp_list = fill_amp_list (target_ampfile, N_targets);
  unsigned int total_size_toamp = 0;
  fclose (target_ampfile);

  for (i = 0; i < N_targets; i++)
  {
    printf ("\nSearching for Primers target #%d: %s chr%d:%d-%d\n",
	    i + 1,
	    target_amp_list[i]->name,
	    target_amp_list[i]->chrom, target_amp_list[i]->start_pos + 1, target_amp_list[i]->stop_pos + 1);
    total_size_toamp += 1 + target_amp_list[i]->stop_pos - target_amp_list[i]->start_pos;
  }
  printf ("\nFinished reading amplicon target coordinates\n\tFound %d targets with a total_size of %u\n", N_targets,total_size_toamp);

  //PCR primer parameters
  read_var ("Minimum Primer Length\n", ss);
  min_primer = atoi (ss);

  read_var ("Maximum Primer Length\n", ss);
  max_primer = atoi (ss);

  read_var ("Minimum Amplicon Length\n", ss);
  amp_min = atoi (ss);

  read_var ("Maximum Amplicon Length\n", ss);
  amp_max = atoi (ss);

  read_var ("Minimum GC content [0..1.0]\n", ss);
  min_gc = (double) atof (ss);

  read_var ("Maximum GC content [0..1.0]\n", ss);
  max_gc = (double) atof (ss);

  read_var ("Minimum tm_primer in degrees C\n", ss);
  min_tm = (double) atof (ss);

  read_var ("Maximum tm_primer in degrees C\n", ss);
  max_tm = (double) atof (ss);

  read_var ("Maximum number of primer pairs to pool together\n", ss);
  pool_size = atoi (ss);

  read_var ("Pad size\n", ss);
  pad_size = atoi (ss);

  // sanity check that pad size
  if (pad_size <= max_primer)
  {
    printf("\nERROR: pad size must be larger than the maximum size of a primer\n");
    exit(1);
  }
  else if (pad_size > 3 * amp_max)
  {
    printf("\nERROR: pad size is unrealistically large, i.e. 3 times the size of the maximum amplicon length\n");
    exit(1);
  }

  read_var ("Tm increment (0.5 to 4.0)\n", ss);
  tm_inc = (double) atof (ss);

  if (tm_inc > 4 || tm_inc <0.5)
  {
    printf("\nError please choose an increment between 0.5 and 4 C\n");
    exit(1);
  }

  // print header for outfile
  fprintf (outfile, "Primer_number\tForward_primer\tForward_Tm\tForward_GC\tReverse_primer\tReverse_Tm\tReverse_GC\tChr\t");
  fprintf (outfile, "Forward_start_position\tForward_stop_position\tReverse_start_position\tReverse_stop_position\t");
  fprintf (outfile, "Product_length\tProduct_GC\tProduct_tm\tProduct\n");

  // allocate some memory for indexed genome
  j = 4;
  for (i = 0; i < idepth; i++)
    j *= 4;
  j = (j - 4) / 3;

  total_index = j;
  printf ("\n Determined the size of flat index to be %d * %lu = %lu  bytes\n",
	  total_index, sizeof (int), total_index * sizeof (int));

  flat_index = ivector (0, total_index);

  for (i = 0; i <= total_index; i++)
    flat_index[i] = 0;

  for (fasta = 0; fasta < N; fasta++)
  {
    j = contig_length[fasta];
    if (j % 4 == 0)
      j /= 4;
    else
      j = j / 4 + 1;

    compressed_map[fasta] = ucvector (0, j + 1);
    printf ("\n For chromosome %d we are going to read %d bytes\n", fasta, j);
    k = fread (compressed_map[fasta], sizeof (unsigned char), j, cfile);

    if (k != j)
    {
      printf ("\nCompressed Sequence %d %s should have been length %d but was %d\n",
	      fasta, contig_descript[fasta], j, k);
      exit (1);
    }
  }
  fclose (cfile);
  k = fread (flat_index, sizeof (int), total_index, idfile);

  if (k < total_index)
  {
    printf ("\nIndexed of N'mers should have been length %d but was %d\n", total_index, k);
    exit (1);
  }
  fclose (idfile);

  // allocated memory for amp_pool => array of regions captured.
  PNODE **all_primer_pairs;
  int max_regions = (1 + maxim(5*total_size_toamp / amp_min,N_targets*2));
  int max_ppairs = MAX_PAIRS * max_regions;
  printf ("\n We have max_regions = %d  max_pairs = %d \n", max_regions, max_ppairs);
  all_primer_pairs = malloc ((unsigned) max_ppairs * sizeof (PNODE *));
  if (!all_primer_pairs)
  {
    printf ("\n Could not allocate space for %u primer pairs \n", max_ppairs);
    exit (1);
  }
  int **redundant_list, *best_start;
  redundant_list = imatrix (0, max_ppairs, 0, MAX_PAIRS);
  best_start = ivector (0, max_regions);
  for(i=0;i<=max_regions;i++)
	best_start[i] = -1;
  for (i = 0; i <= max_ppairs; i++)
  {
    for (j = 0; j <= MAX_PAIRS; j++)
      redundant_list[i][j] = -1;
  }
  int amp_pool_count = 0;
  int primer_count = 0;

  // start finding primers for targets
  AMPNODE *loop_amp;
  loop_amp = amp_alloc ();
  printf("\n loop_amp is located at position %ld in memory \n\n",(long)loop_amp);
  for (k = 0; k < N_targets; k++)
  {
    // set chromosome
    target_contig = target_amp_list[k]->chrom - 1;

    // set target start and stop
    int target_start = target_amp_list[k]->start_pos + 1;
    int target_stop = target_amp_list[k]->stop_pos + 1;

    // set temporary start and stop
    int this_start, this_stop, this_midpoint;
    if (target_stop - target_start > amp_max)
    {
      this_start = target_start;
      this_stop = target_start + amp_max;
    }
    else
    {
      this_midpoint = ((target_stop - target_start) / 2) + target_start;
      this_start    = this_midpoint - (amp_min / 2);
      this_stop     = this_midpoint + (amp_min / 2);
    }

    // while loop variable: call it the "region" loop
    int not_covered = 1;

    // print out target info
    printf ("\nSearching for Primers target #%d: %s chr%d:%d-%d\n\n",
	    k + 1,
	    target_amp_list[k]->name,
	    target_amp_list[k]->chrom, target_amp_list[k]->start_pos + 1, target_amp_list[k]->stop_pos + 1);
    printf ("\nStarting params this_start = %d, this_stop = %d\n\n", this_start, this_stop);

    // variables:
    //      target_start, target_stop => region to cover with amplicons
    //      this_start, this_stop => region to cover with the specific iteration of the while loop
    //      genome_start, genome_stop => coordinates used to extract actual genomic region to target
    //


    while (not_covered)
    {
      printf("\n 2 loop_amp is located at position %ld in memory \n\n",(long)loop_amp);

      // start site - make sure it's a multiple of 4
      genome_start = maxim (1, this_start - pad_size);
      j = genome_start / 4;
      genome_start = j * 4;

      // stop site - make sure it's a multiple of 4
      genome_stop = (this_stop + pad_size) / 4;
      genome_stop = minim (genome_stop * 4, contig_length[target_contig]);

      // length of region
      j = (genome_stop - genome_start);
      j = minim (j, contig_length[target_contig]);

      // set loop_amp attributes
      loop_amp->start_pos = this_start + 1;
      loop_amp->stop_pos = this_stop + 1;
      loop_amp->chrom = target_amp_list[k]->chrom;
      loop_amp->no_pairs = 0;
      sprintf (loop_amp->name, "%s_%02d", target_amp_list[k]->name, amp_pool_count);

      // tm range
      int this_min_tm = min_tm;
      int this_max_tm = max_tm;

      // while loop variables: call it the "specific target" loop
      int found_count = 0;
      int this_trial = 0;

      // scratch pad
      if (j > 2000)
      {
				printf ("\nj (length) is too big %d\n\n", j);
				exit (1);
      }
      // allocate memory for scratch pad
      scratch_pad = cvector (0, j + 5);
      scratch_pad[0] = '\0';
      decode_basepairs (&compressed_map[target_contig][genome_start / 4], scratch_pad, j / 4);

      do
      {
        found_count = find_primers (snp_list[target_contig], loop_amp, contig_snp_count[target_contig],	// number of snps in the contig
				  flat_index, scratch_pad,	// copy of the target
				  j,	// length of region
				  min_primer, max_primer, amp_max, amp_min, min_gc, max_gc, this_min_tm, this_max_tm, idepth,	// index depth
				  10,	// local depth
				  j / 2,	// target base (from old primer_snp.c program)
				  genome_start + 1, highmer,	// not sure
				  repeat_list[target_contig], no_repeats[target_contig], pad_size,	// size on either end of the contig to look for primer
				  loop_amp->chrom);
        this_min_tm -= tm_inc;
        this_max_tm += tm_inc;
        this_trial++;
        if (this_min_tm < min_tm || this_max_tm > max_tm) {
          break;
        }
      } while (this_trial < 10 && found_count < 1);
      printf("\n 3 loop_amp is located at position %ld in memory \n\n",(long)loop_amp);


      free_cvector (scratch_pad, 0, j + 5);
      int nearest_stop = 0;
      if (found_count > 0)
      {
				if (loop_amp->no_pairs > 0)
				{
				  if (amp_pool_count == max_ppairs)
				  {
				    printf ("ERROR: exceeded the maximum number of primer pairs allocated: %d\n\n", max_ppairs);
				    exit (1);
				  }

				  printf ("\nBED: chr%d\t%d\t%d\n", loop_amp->chrom, this_start, this_stop);

				  int start_primer_count = primer_count;

				  printf("\n\nabout to fill the array of primers start = %d with %d pairs coming\n\n",
				     start_primer_count,loop_amp->no_pairs);

				  for (j = 0; j < loop_amp->no_pairs; j++)
				    all_primer_pairs[primer_count++] = loop_amp->pair[j];
				  int jj;
				  for (j = start_primer_count; j < primer_count; j++)
          {
				    for (jj = 0; jj < loop_amp->no_pairs; jj++)
            {
              printf("\n\nj = %d, jj = %d, start_primer_count = %d, primer_count = %d\n\n", j, jj, start_primer_count, primer_count);
				      redundant_list[j][jj] = start_primer_count + jj;
      			printf("\n 4 loop_amp is located at position %ld in memory \n\n",(long)loop_amp);
            }
          }
				printf("\n Made it here \n\n");
				  best_start[amp_pool_count] = start_primer_count;
				  amp_pool_count++;

				  for (j = 0; j < loop_amp->no_pairs; j++)
				  {
				    if (nearest_stop == 0)
				      nearest_stop = loop_amp->pair[j]->reverse->start;
				    else if (loop_amp->pair[j]->reverse->start < nearest_stop)
				      nearest_stop = loop_amp->pair[j]->reverse->start;
				  }
				}
      }
      printf("\n\nnearest stop is %d\n\n", nearest_stop);

      if (nearest_stop > 0)
      {
				this_start = nearest_stop;
				this_stop = this_start + amp_min;
      }
      else
      {
				this_start += pad_size;
				this_stop = this_start + amp_min;
      }
      if (nearest_stop >= target_stop || this_start > target_stop)
				not_covered = 0;
    }
  }

  char **poolable_matrix;
  int *poolable_count;
  poolable_count = ivector (0, primer_count);
  for (i = 0; i < primer_count; i++)
    poolable_count[i] = 0;
  poolable_matrix = cmatrix (0, primer_count, 0, primer_count);
  for (i = 0; i < primer_count; i++)
    for (j = i + 1; j < primer_count; j++)
    {
      poolable_matrix[i][j] = is_poolable_primer (all_primer_pairs[i], all_primer_pairs[j]);
      poolable_matrix[j][i] = poolable_matrix[i][j];
      poolable_count[i] += poolable_matrix[i][j];
      poolable_count[j] += poolable_matrix[i][j];
    }
  for (i = 0; i < amp_pool_count; i++)
    if (best_start[i] >= 0)
    {
      int k = best_start[i];
      for (j = 0; j < MAX_PAIRS; j++)
	if (redundant_list[k][j] >= 0)
	  if (poolable_count[redundant_list[k][j]] > poolable_count[best_start[i]])
	    best_start[i] = redundant_list[k][j];
    }
  //
  // print cmat
  //
  //printf ("\n");
	//for (i = 0; i < primer_count; i++)
	//printf ("\t%s_%s_%03d", all_primer_pairs[i]->forward->sequence, all_primer_pairs[i]->reverse->sequence, i);
	//for (i = 0; i < primer_count; i++)
	//{
	//printf ("\n%s_%s_%03d", all_primer_pairs[i]->forward->sequence, all_primer_pairs[i]->reverse->sequence, i);
	//for (j = 0; j < primer_count; j++)
	//if (i != j)
	//printf ("\t%d", (int) poolable_matrix[i][j]);
	//else
	//printf ("\t.");
	//}
	//printf ("\n");
  //
  // print best start matrix
  //
  //  for (i = 0; i < amp_pool_count; i++)
  //  {
  //    int k = best_start[i];
  //    for (j = 0; j < MAX_PAIRS; j++)
  //      if (redundant_list[k][j] >= 0)
  //        printf (" %03d", redundant_list[k][j]);
  //      else
  //        printf ("   .");
  //    printf ("\t| Region: %03d\t| Primer: %03d\t| Poolable Count: %03d\n", i, best_start[i], poolable_count[best_start[i]]);
  //  }
  //
  // redundant matrix
  //
  //  printf("\n\n");
  //  for (i=0; i<primer_count; i++)
  //  {
  //    for (j=0;j<MAX_PAIRS; j++)
  //      if (redundant_list[i][j] >= 0)
  //        printf (" %03d", redundant_list[i][j]);
  //      else
  //        printf ("   .");
  //    printf("\n");
  //  }

  printf ("\n\n going to make_pools with amp_pools = %d, primer_count = %d\n", amp_pool_count, primer_count);
  int *current_pool;
  current_pool = ivector (0, 20);
  make_pools (all_primer_pairs, poolable_matrix, poolable_count, redundant_list, best_start, amp_pool_count,
	       primer_count, primer_count, current_pool, 0, pool_size);
  return 0;
}

/*---------------------------------------------------------------------*/
void
make_pools (PNODE ** plist, char **cmat, int *pc, int **redund, int *bstart, int Nregs, int N, int still,
	     int *current_pool, int this_pool, int max_pool)
{
  int i, j;
  if (still < 1)
    return;

  int this_p = -1;
  if (this_pool == 0)
  {
    // fprintf (outfile, "\nStarting a New Pool\n");
    int best = -1;
    for (i = 0; i < Nregs; i++)
      if (bstart[i] >= 0)
				if (pc[bstart[i]] > best)
				{
				  best = pc[bstart[i]];
				  this_p = bstart[i];
				}
    if (this_p < 1)
    {
      printf ("\n This is impossible \n");
      for (i = 0; i < Nregs; i++)
      	printf ("\n Best start region %d is %d", i, bstart[i]);
      exit (1);
    }
  }
  else
  {
    int best = -1;
    for (i = 0; i < N; i++)
      if (plist[i] != NULL)
	      if (pc[i] > best)
      	{
      	  int it_fits = TRUE;
      	  for (j = 0; j < this_pool; j++)
      	    if (!cmat[i][current_pool[j]])
	            it_fits = FALSE;
      	  if (it_fits)
      	  {
      	    best = pc[i];
	          this_p = i;
          }
	      }
  }

  printf ("\n Got here with this_p = %d \n", this_p);
  if (this_p >= 0)
  {
    fprintf (outfile, "%d\t%s\t%g\t%g\t%s\t%g\t%g\t%d\t%d\t%d\t%d\t%d\t%d\t%g\t%g\t%s\n",
	     this_pool,
	     plist[this_p]->forward->sequence, plist[this_p]->forward->tm, plist[this_p]->forward->gc,
	     plist[this_p]->reverse->sequence, plist[this_p]->reverse->tm, plist[this_p]->reverse->gc,
	     plist[this_p]->chrom,
	     plist[this_p]->forward->start, plist[this_p]->forward->end,
	     plist[this_p]->reverse->start, plist[this_p]->reverse->end,
	     plist[this_p]->length, plist[this_p]->gc, plist[this_p]->tm, plist[this_p]->sequence);
    current_pool[this_pool] = this_p;
    for (i = 0; i < MAX_PAIRS; i++)
      if (redund[this_p][i] >= 0)
      {
				int ii = redund[this_p][i];
				for (j = 0; j < Nregs; j++)
				  if (bstart[j] == ii)
				    bstart[j] = -1;
				plist[ii] = NULL;
				pc[ii] = 0;
				still--;
      }
    this_pool++;
    this_pool %= max_pool;
  }
  else
    this_pool = 0;

  printf ("\ngoing to make_pools with N = %d, still = %d, this_pool = %d, max_pool = %d\n\n", N, still, this_pool,
	  max_pool);

  make_pools (plist, cmat, pc, redund, bstart, Nregs, N, still, current_pool, this_pool, max_pool);
}
/*---------------------------------------------------------------------*/

int
is_poolable_primer (PNODE * p1, PNODE * p2)
{
  int i, f1, r1, f2, r2;
  char flipf1[80], flipf2[80], flipr1[80], flipr2[80];

  for (i = p1->forward->start; i <= p1->reverse->end; i++)
    if ((i >= p2->forward->start) && (i <= p2->reverse->end))
      return FALSE;

  for (i = p2->forward->start; i <= p2->reverse->end; i++)
    if ((i >= p1->forward->start) && (i <= p1->reverse->end))
      return FALSE;

  f1 = abs (p1->forward->end - p1->forward->start) + 1;
  r1 = abs (p1->reverse->end - p1->reverse->start) + 1;
  f2 = abs (p2->forward->end - p2->forward->start) + 1;
  r2 = abs (p2->reverse->end - p2->reverse->start) + 1;

  reverse_string (p1->reverse->sequence, flipr1, r1);
  reverse_string (p2->reverse->sequence, flipr2, r2);
  reverse_string (p1->forward->sequence, flipf1, f1);
  reverse_string (p2->forward->sequence, flipf2, f2);

  if (check_uneven_dimer (p1->forward->sequence, flipf2, f1, f2))
    return FALSE;

  if (check_uneven_dimer (p1->reverse->sequence, flipr2, r1, r2))
    return FALSE;

  if (check_uneven_dimer (p1->forward->sequence, flipr2, f1, r2))
    return FALSE;

  if (check_uneven_dimer (p1->reverse->sequence, flipf2, r1, f2))
    return FALSE;

  if (fabs (p1->forward->tm - p2->forward->tm) > 2)
    return FALSE;

  if (fabs (p1->reverse->tm - p2->reverse->tm) > 2)
    return FALSE;

  if (abs (p1->length - p2->length) > 50)
    return FALSE;

  return TRUE;
}
/*---------------------------------------------------------------------*/

double
fill_dist (SNODE ** list, int n, int *priority, int *selected, int which, int start, int stop, double *dist)
{
  int i, j, flag, m;
  double totd, dtemp;

  totd = 0;
    j = m = 0;

  for (i = 0; i < n; i++)
  {
    if (priority[i] == which)
    {
      flag = TRUE;
      j = i - 1;
      dtemp = 0;
      while ((j >= 0) && (flag))
      {
	j--;
	if (j >= 0)
	  if (selected[j] == 1)
	    flag = FALSE;
      }
      if (j >= 0)
          dtemp = list[i]->pos - list[j]->pos;
      else
          dtemp = list[i]->pos - start;

      dist[i] = dtemp * dtemp;

      j = i + 1;

      flag = TRUE;
      dtemp = 0;
      while ((j < n) && (flag))
      {
          j++;
          if (j < n)
          if (selected[j] == 1)
              flag = FALSE;
      }
      if (j < n)
          dtemp = list[j]->pos - list[i]->pos;
      else
          dtemp = stop - list[i]->pos;

      dist[i] += dtemp * dtemp;
    }
    else
      dist[i] = 0;

    totd += dist[i];
    if (selected[i] == 1)
      m++;
  }

  if (totd > 0)
    for (i = 0; i < n; i++)
      dist[i] /= totd;

  return totd;

}

/*---------------------------------------------------------------------*/
int
pick_random (double *x, int n)
{
  int i;
  double y, last;


  y = ran1 (&idum);

  i = 0;
  last = 0;

  while ((i < n) && (y > last))
  {
    if (x[i] > 0)
      last += x[i];
    i++;
  }

  if (i > 0)
    return i - 1;
  return i;

}

/*---------------------------------------------------------------------*/
SNODE *
snp_alloc (void)
{
  SNODE *tn;

  tn = (SNODE *) malloc ((unsigned) sizeof (struct snp_node));
  if (!tn)
    log_err ("allocation failure in snp_alloc()");

  sprintf (tn->name, "tempname");
  tn->no_disc = 0;
  tn->chrom = 0;
  tn->pos = 0;
  tn->het = 0.0;
  tn->baseA = 'N';
  tn->baseB = 'N';
  tn->no_pairs = 0;
  tn->pair = (PNODE **) malloc ((unsigned) (MAX_PAIRS + 1) * sizeof (PNODE *));
  if (!tn->pair)
    log_err ("allocation failure 2 in snp_alloc()");

  return tn;
}

/*---------------------------------------------------------------------*/
AMPNODE *
amp_alloc (void)
{
  AMPNODE *tn;

  tn = (AMPNODE *) malloc ((unsigned) sizeof (struct amp_node));
  if (!tn)
    log_err ("allocation failure in amp_alloc()");

  sprintf (tn->name, "tempname");
  tn->chrom = 0;
  tn->start_pos = 0;
  tn->stop_pos = 0;
  tn->no_pairs = 0;
  tn->pair = (PNODE **) malloc ((unsigned) (MAX_PAIRS + 1) * sizeof (PNODE *));
  if (!tn->pair)
    log_err ("allocation failure 2 in amp_alloc()");

  return tn;
}

/*---------------------------------------------------------------------*/
int
isbase (char c)
{
  char C;

  C = toupper (c);

  if ((C == 'A') || (C == 'C') || (C == 'G') || (C == 'T'))
    return TRUE;

  return FALSE;
}

/*---------------------------------------------------------------------*/
SNODE **
fill_snp_list (FILE * sfile, int *n, unsigned int chrom)
{
  int flag, i, temp_no, final_no;
  char sss[4096], s[256];
  SNODE **temp_snp;

  i = 0;
  flag = TRUE;

  printf ("reading snps for chr %d\n", chrom);
  temp_snp = (SNODE **) malloc ((unsigned) (2000000) * sizeof (SNODE *));
  if (!temp_snp)
    log_err ("Allocation failure in Temporary SNP storage");

  fgets (sss, 4094, sfile);

  while (flag)
  {
    temp_snp[i] = snp_alloc ();
    sscanf (sss, "%s\t%d\t%d\t%d\t%s\t%c/%c", temp_snp[i]->name, &temp_snp[i]->no_disc, &temp_snp[i]->chrom, &temp_snp[i]->pos, s,
	    &temp_snp[i]->baseA, &temp_snp[i]->baseB);

    temp_snp[i]->pos--;
    temp_snp[i]->het = (double) atof (s);
    temp_snp[i]->no_pairs = 0;

    if ((isbase (temp_snp[i]->baseA)) && (isbase (temp_snp[i]->baseB))
	&& (temp_snp[i]->chrom == chrom) && (temp_snp[i]->het > 0.05))
    {
      i++;
    }
    else
    {
      free (temp_snp[i]);
    }

    sss[0] = '\0';
    if ((!feof (sfile)) && (i < 2000000))
    {
      fgets (sss, 4094, sfile);
      if (strlen (sss) < 3)
      	flag = FALSE;
    }
    else
      flag = FALSE;
  }

  temp_no = i;
  qsort (temp_snp, temp_no, sizeof (SNODE *), sort_compare_struct);
  final_no = i;
  (*n) = final_no;
  printf ("\n Total number of SNPs found = %d\n", final_no);

  return temp_snp;
}

/*---------------------------------------------------------------------*/
AMPNODE **
fill_amp_list (FILE * sfile, int n)
{
  int flag, i, temp_no;
  char sss[4096];
  AMPNODE **temp_amp;

  i = 0;
  temp_no = n;
  flag = TRUE;

  temp_amp = (AMPNODE **) malloc ((unsigned) (temp_no + 1) * sizeof (AMPNODE *));
  if (!temp_amp)
    log_err ("Allocation failure in Temporary SNP storage");

  fgets (sss, 4094, sfile);

  while (flag)
  {
    temp_amp[i] = amp_alloc ();
    sscanf (sss, "chr%d\t%d\t%d\t%s",
	    &temp_amp[i]->chrom, &temp_amp[i]->start_pos, &temp_amp[i]->stop_pos, temp_amp[i]->name);

    temp_amp[i]->start_pos--;
    temp_amp[i]->stop_pos--;
    temp_amp[i]->no_pairs = 0;
    i++;
    sss[0] = '\0';
    if ((!feof (sfile)) && (i < n))
    {
      fgets (sss, 4094, sfile);
      if (strlen (sss) < 3)
	flag = FALSE;
    }
    else
      flag = FALSE;
  }
  printf ("\n Total number of regions found = %d\n", i);
  return temp_amp;
}

/*---------------------------------------------------------------------*/
int
sort_compare_struct (const void *a, const void *b)
{
  //printf("\n Comparing %d with %d ", (*(SNODE **) a)->pos, (*(SNODE **) b)->pos);
  if ((*(SNODE **) a)->pos < (*(SNODE **) b)->pos)
    return -1;
  else if ((*(SNODE **) a)->pos > (*(SNODE **) b)->pos)
    return 1;
  else
    return 0;
}

/*---------------------------------------------------------------------*/
void
fill_hs (char a, char b, double *h, double *s)
{
  if (a == 'A')
  {
    if (b == 'A')
    {
      (*h) += 9100;
      (*s) += 24;
    }
    else if (b == 'C')
    {
      (*h) += 6500;
      (*s) += 17.3;
    }
    else if (b == 'G')
    {
      (*h) += 7800;
      (*s) += 20.8;
    }
    else if (b == 'T')
    {
      (*h) += 8600;
      (*s) += 23.9;
    }
  }
  else if (a == 'C')
  {
    if (b == 'A')
    {
      (*h) += 5800;
      (*s) += 12.9;
    }
    else if (b == 'C')
    {
      (*h) += 11000;
      (*s) += 26.6;
    }
    else if (b == 'G')
    {
      (*h) += 11900;
      (*s) += 27.8;
    }
    else if (b == 'T')
    {
      (*h) += 7800;
      (*s) += 20.8;
    }
  }
  else if (a == 'G')
  {
    if (b == 'A')
    {
      (*h) += 5600;
      (*s) += 13.5;
    }
    else if (b == 'C')
    {
      (*h) += 11100;
      (*s) += 26.7;
    }
    else if (b == 'G')
    {
      (*h) += 11000;
      (*s) += 26.6;
    }
    else if (b == 'T')
    {
      (*h) += 6500;
      (*s) += 17.3;
    }
  }
  else if (a == 'T')
  {
    if (b == 'A')
    {
      (*h) += 6000;
      (*s) += 16.9;
    }
    else if (b == 'C')
    {
      (*h) += 5600;
      (*s) += 13.5;
    }
    else if (b == 'G')
    {
      (*h) += 5800;
      (*s) += 12.9;
    }
    else if (b == 'T')
    {
      (*h) += 9100;
      (*s) += 24.0;
    }
  }
  /* printf("\n For %c %c we have H = %g  S = %g",a,b,(*h),(*s)); */

}

/*---------------------------------------------------------------------*/
double
calc_tm (char *ss, int n)
{
  int i;
  double h, s, tm;

  h = s = 0.0;

  for (i = 0; i < n - 1; i++)
    fill_hs (ss[i], ss[i + 1], &h, &s);

  /* tm = h/(s + 57.6945289) - 21.4624334 - 273.15; */
  /* tm = h/(s + 57.6945289) - 294.6124334; */
  tm = h / (s + 47.16510465) - 294.6124334;

  return tm;
}

/*---------------------------------------------------------------------*/

void
convert_int_basepairs (int i, char *s, int k)
{
  int j, l;
  char ss[256];

  l = k - 1;
  for (j = 0; j < 256; j++)
    ss[j] = 'A';

  while (l >= 0)
  {
    j = i % 4;
    if (j == 0)
      ss[l] = 'A';
    else if (j == 1)
      ss[l] = 'C';
    else if (j == 2)
      ss[l] = 'G';
    else
      ss[l] = 'T';
    i -= j;
    i /= 4;
    l--;
  }

  ss[k] = '\0';
  strcpy (s, ss);

}

/*---------------------------------------------------------------------*/
void
decode_basepairs (unsigned char *s, char *dest, int n)
{
  int i, m, l;
  unsigned char j;
  char a, ss[5];

  ss[4] = '\0';
  for (i = 0; i < n; i++)
  {
    j = s[i];

    //printf("\n Decoding %d ", j);
    for (l = 3; l >= 0; l--)
    {
      m = j % 4;
      if (m == 0)
	a = 'A';
      else if (m == 1)
	a = 'C';
      else if (m == 2)
	a = 'G';
      else
	a = 'T';

      ss[l] = a;
      j = j >> 2;;

    }

    for (m = 0; m < 4; m++)
      *dest++ = ss[m];

    //printf("as %s", ss);
  }
  *dest = '\0';
}

/*---------------------------------------------------------------------*/
unsigned int
encode_basepairs (char *ss, int n)
{
  unsigned int k;
  int i;

  k = 0;
  for (i = 0; i < n; i++)
  {
    k = k << 2;
    if (ss[i] == 'A');
    else if (ss[i] == 'C')
      k++;
    else if (ss[i] == 'G')
      k += 2;
    else
      k += 3;
  }
  // printf("\nn = %d Encoding %c%c%c%c as %d ",n,ss[0],ss[1],ss[2],ss[3],k);

  return k;
}

/*---------------------------------------------------------------------*/
static double *FQ_LIST;

int
find_primers (SNODE ** snp_list, AMPNODE * tn, int no_snp, int *flat, char *contig, int L, int min_primer, int max_primer, int amp_max, int amp_min, double min_gc, double max_gc, double min_tm, double max_tm, int depth, int local_depth, int target_base,	// originally the bp of the SNP that dave's original primer_snp program was looking for
	      int start_pos, unsigned char *highmer, int **repeats, int no_repeats, int end_region, int chrom)
{
    int i, j, k, *index_left, *index_right, fl, fr, this_amp, pairs_todump;
    int *amp_size, *right_side, *left_side, *local_index, local_size, temp_length;
    int *plen_l, *plen_r;
    double *gc_left, *gc_right, *fq_left, *fq_right;
    double *total_fq, *best_gc_left, *best_gc_right, tm_l, tm_r;
    char ss[256], sss[256], flip[256], *rt_contig, **best_left, **best_right;
    PNODE *product;
    PRIMER *p_left, *p_right;

    // printf("\nForward sequence\n\n");
    // printf("\n In tile contig with L = %d; amp_max = %d;  amp_min = %d; depth = %d\n\n", L, amp_max, amp_min, depth);
    // for (i = 0; i < L; i++)
    // {
    //     printf("%c", contig[i]);
    //     //if ((i + 1 == target_base) || (i == target_base)) printf(" * ");
    //     if (i % 80 == 79) printf("\n");
    // }
    // printf("\n");

    fq_left = dvector (0, L);
    fq_right = dvector (0, L);
    gc_left = dvector (0, L);
    gc_right = dvector (0, L);
    index_left = ivector (0, L);
    index_right = ivector (0, L);
    plen_l = ivector (0, L);
    plen_r = ivector (0, L);
    printf("\n Local Depth is %d \n\n",local_depth);

    j = 4;
    for (i = 0; i < local_depth; i++)
        j *= 4;
    j = (j - 4) / 3;

    local_size = j;
    local_index = ivector (0, local_size);
    for (i = 0; i <= local_size; i++)
        local_index[i] = 0;

    flat_index_contig (local_index, contig, L, local_depth);

    printf("\n About to fill quality in forward direction\n\n");

    fill_quality_scores (flat, local_index, contig, L, minim (target_base - 30, end_region), depth, local_depth,
                         min_primer, max_primer, fq_left, gc_left, index_left, plen_l,
                         min_gc, max_gc, min_tm, max_tm, repeats, no_repeats, highmer, start_pos);
    rt_contig = cvector (0, L);
    reverse_transcribe (contig, rt_contig, L);

    // printf("\nReverse Transcribe\n\n");
    // printf("\n");for (i = 0; i < L; i++) {printf("%c", rt_contig[i]);
    // if ((i + 1 == L - target_base) || (i + 2 == L - target_base)) printf(" * ");
    // if (i % 80 == 79)
    //     printf("\n");
    // }
    // printf("\n");
    // printf ("\n About to fill quality in reverse direction\n");

    fill_quality_scores (flat, local_index, rt_contig, L, minim (target_base - 30, end_region), depth, local_depth,
                         min_primer, max_primer, fq_right, gc_right, index_right, plen_r,
                         min_gc, max_gc, min_tm, max_tm, repeats, no_repeats, highmer, start_pos);

    i = 0;
    fl = 0;
    pairs_todump = MAX_PAIRS;
    amp_size = ivector (0, pairs_todump);
    best_left = cmatrix (0, pairs_todump, 0, 256);
    best_right = cmatrix (0, pairs_todump, 0, 256);
    total_fq = dvector (0, pairs_todump);
    right_side = ivector (0, pairs_todump);
    left_side = ivector (0, pairs_todump);
    best_gc_left = dvector (0, pairs_todump);
    best_gc_right = dvector (0, pairs_todump);

    printf("\n About to start finding primers\n");

    while ((fq_left[index_left[i]] < 1e7) && (fl < pairs_todump))
    {
        k = start_pos + index_left[i];
        //printf("\n pos (k) = %d\n\n", k);
        if ((!poly_under_primer (k, k + plen_l[index_left[i]], snp_list, 0, no_snp - 1, (no_snp) / 2)) &&
            (gc_left[index_left[i]] >= min_gc) && (gc_left[index_left[i]] <= max_gc) &&
            (is_not_repeat (k, k + plen_l[index_left[i]], repeats, no_repeats, 0.001)))
        {
            j = 0;
            fr = 0;
            strncpy (ss, contig + index_left[i], plen_l[index_left[i]]);
            ss[plen_l[index_left[i]]] = '\0';
            tm_l = calc_tm (ss, plen_l[index_left[i]]);

            if (fq_right[index_right[j]] >= 1e7)
            {
                printf("\n No right side primers \n\n");
                break;
            }
            // printf("\nMatching 5' %s 3' (fq = %g, gc = %g, len = %d  tm = %g) with", ss, fq_left[index_left[i]], gc_left[index_left[i]], plen_l[index_left[i]], tm_l);

            while ((fq_right[index_right[j]] < 1e7) && (fr < 1))
            {
                temp_length = L - index_right[j] - index_left[i];

        				/*
        				 * strncpy(sss,rt_contig+index_right[j],primer
        				 * ); sss[primer] = '\0'; printf("\n This
        				 * primer pair %s appears to %d
        				 * \n",sss,temp_length);
        				 */

                k = start_pos + L - index_right[j];
                if ((!poly_under_primer (k - plen_r[index_right[j]], k, snp_list, 0, no_snp - 1, (no_snp) / 2)) &&
                    (temp_length >= amp_min) && (temp_length <= amp_max) &&
                    (is_not_repeat (k - plen_r[index_right[j]], k, repeats, no_repeats, 0.001)))
                {
                    strncpy (sss, rt_contig + index_right[j], plen_r[index_right[j]]);
                    sss[plen_r[index_right[j]]] = '\0';
                    reverse_string (sss, flip, plen_r[index_right[j]]);
                    if (check_uneven_dimer (ss, flip, plen_l[index_left[i]], plen_r[index_right[j]]))
                    {
                        printf("\n\t\t\tNot 5' %s 3' because of a dimer", sss);
                    }
                    else
                    {
                        tm_r = calc_tm (sss, plen_r[index_right[j]]);
                        printf("\n In here with tm_l = %g and tm_r  = %g\n\n",tm_r,tm_l);
                        if (fabs (tm_l - tm_r) < 5.0)
                        {
                            if (fr == 0)
                            {
                                this_amp = L - index_right[j];
                                right_side[fl] = this_amp;
                                left_side[fl] = index_left[i];
                                amp_size[fl] = temp_length;
                                sprintf (best_left[fl], "%s", ss);
                                sprintf (best_right[fl], "%s", sss);
                                total_fq[fl] = fq_right[index_right[j]] + fq_left[index_left[i]];
                                best_gc_left[fl] = gc_left[index_left[i]];
                                best_gc_right[fl] = gc_right[index_right[j]];
                            }
                            product = product_alloc ();
                            p_left = primer_alloc ();
                            p_right = primer_alloc ();
                            p_left->sequence = cvector (0, plen_l[index_left[i]] + 1);
                            sprintf (p_left->sequence, "%s", ss);
                            p_right->sequence = cvector (0, plen_r[index_right[j]] + 1);
                            sprintf (p_right->sequence, "%s", sss);
                            p_left->tm = tm_l;
                            p_right->tm = tm_r;
                            p_left->gc = gc_left[index_left[i]];
                            p_right->gc = gc_right[index_right[j]];
                            p_left->start = start_pos + index_left[i];
                            p_left->end = p_left->start + plen_l[index_left[i]];
                            p_right->end = start_pos + L - (index_right[j]);
                            p_right->start = p_right->end - plen_r[index_right[j]];
                            p_left->end--;
                            p_right->end--;
                            product->forward = p_left;
                            product->reverse = p_right;
                            product->sequence = cvector (0, temp_length + 1);
                            strncpy (product->sequence, contig + index_left[i], temp_length);
                            product->sequence[temp_length] = '\0';
                            product->length = temp_length;
                            product->gc = calc_gc (product->sequence, temp_length);
                            product->tm = 41.0 * product->gc - 675.0 / (double) product->length - 21.4624334;
                            product->chrom = chrom;
                            tn->pair[tn->no_pairs++] = product;

                            printf("\n\tSuccess with 5' %s 3' (fq = %g, gc = %g, len = %d tm = %g)", sss, fq_right[index_right[j]], gc_right[index_right[j]], plen_r[index_right[j]], tm_r);
                            fr++;
                        }
                    }
                }
                j++;
            }
            if (fr > 0)
                fl++;
        }
        i++;
    }


    if (fl > 0)
    {
        for (i = 0; i <= pairs_todump; i++)
            index_left[i] = i;


        FQ_LIST = total_fq;
        qsort ((void *) index_left, fl, sizeof (int), sort_compare_index);


        this_amp = right_side[index_left[0]];

        //for(i=0;i<fl;i++) fprintf(outfile,"\n5' %s 3'\twith\t5' %s 3'\tsize\t%d\t%d\tto\t%d\t%g\t%g", best_left[index_left[i]],best_right[index_left[i]],amp_size[index_left[i]], start_pos+left_side[index_left[i]],start_pos+right_side[index_left[i]]-1, best_gc_left[index_left[i]],best_gc_right[index_left[i]]); fprintf(outfile,"\n\n");

    }
    else
    {
        this_amp = 0;
        //fprintf(outfile,"\n\n **** Have completely failed to find a primer pair *******\n\n");
        //printf("\n\n **** Have completely failed to find a primer pair *******\n\n");

    }


    free_ivector (local_index, 0, local_size);
    free_ivector (amp_size, 0, pairs_todump);
    free_cmatrix (best_left, 0, pairs_todump, 0, 256);
    free_cmatrix (best_right, 0, pairs_todump, 0, 256);
    free_dvector (total_fq, 0, pairs_todump);
    free_ivector (right_side, 0, pairs_todump);
    free_ivector (left_side, 0, pairs_todump);
    free_dvector (best_gc_left, 0, pairs_todump);
    free_dvector (best_gc_right, 0, pairs_todump);

    free_ivector (index_left, 0, L);
    free_ivector (index_right, 0, L);
    free_ivector (plen_l, 0, L);
    free_ivector (plen_r, 0, L);

    free_cvector (rt_contig, 0, L);

    free_dvector (gc_left, 0, L);
    free_dvector (fq_left, 0, L);
    free_dvector (gc_right, 0, L);
    free_dvector (fq_right, 0, L);


    return this_amp;

}

/*---------------------------------------------------------------------*/
// !poly_under_primer(k, k + plen_l[index_left[i]], snp_list, 0, no_snp - 1, (no_snp) / 2)

int
poly_under_primer (unsigned int p_start, unsigned int p_end, SNODE ** list, int start, int stop, int which)
{
  //printf("\n In poly_under with primer from %d to %d and the list start = %d stop = %d which = %d, list[which] = %d\n", p_start,p_end,start,stop,which,list[which]->pos);

  if ((list[which]->pos >= p_start) && (list[which]->pos <= p_end))
  {
    //printf("\nA primer which goes from %d to %d appears to have %s under it at pos %d\n\n", p_start,p_end,list[which]->name,list[which]->pos);
    return TRUE;
  }
  if (stop - start <= 1)
    return FALSE;

  if (list[which]->pos > p_end)
  {
    if (start >= which)
      return FALSE;
    else
      return poly_under_primer (p_start, p_end, list, start, which, (which + start) / 2);
  }
  if (list[which]->pos < p_start)
  {
    if (which >= stop)
      return FALSE;
    else
      return poly_under_primer (p_start, p_end, list, which, stop, (stop + which) / 2);
  }
  return FALSE;


}

/*---------------------------------------------------------------------*/
double
calc_gc (char *s, int n)
{
  char c;
  int i, gc;

  if (n <= 0)
    return 0;

  gc = 0;
  for (i = 0; i < n; i++)
  {
    c = toupper (s[i]);
    if ((c == 'G') || (c == 'C'))
      gc++;
  }

  return (double) gc / (double) n;
}

/*---------------------------------------------------------------------*/
PRIMER *
primer_alloc (void)
{
  PRIMER *tn;

  tn = (PRIMER *) malloc ((unsigned) sizeof (struct primer_node));
  if (!tn)
    log_err ("allocation failure in primer_alloc");

  tn->start = 0;
  tn->end = 0;
  tn->sequence = NULL;
  tn->tm = 0;
  tn->gc = 0;

  return tn;

}

/*---------------------------------------------------------------------*/
PNODE *
product_alloc (void)
{
  PNODE *tn;

  tn = (PNODE *) malloc ((unsigned) sizeof (struct primer_pair));
  if (!tn)
    log_err ("allocation failure in primer_alloc()");

  tn->forward = NULL;
  tn->reverse = NULL;
  tn->sequence = NULL;
  tn->length = 0;
  tn->gc = 0;
  tn->tm = 0;
  tn->chrom = 0;

  return tn;
}

/*---------------------------------------------------------------------*/
void
fill_quality_scores (int *flat, int *local, char *contig, int L, int window, int depth, int ld, int min_primer,
		     int max_primer, double *fq_left, double *gc_left, int *index_left, int *plen, double min_gc,
		     double max_gc, double min_tm, double max_tm, int **repeats, int no_repeats, unsigned char *highmer,
		     int start_pos)
{
  int **reps_left, **local_reps, tail, primer;
  int i, j, k, m, rm, min_match, offset, *uflag;
  char flip[256];
  double discount, thisd, tm, fq_temp, gc_temp;

  reps_left = imatrix (0, window, 0, depth);
  local_reps = imatrix (0, window, 0, depth);

  for (j = 0; j <= window; j++)
    for (i = 0; i <= depth; i++)
      local_reps[j][i] = reps_left[j][i] = 0;

  uflag = ivector (0, window);
  k = no_repeats / 2;

  for (i = 0; i < window; i++)
  {
    count_copys (local, contig + i, ld, local_reps[i]);
    count_copys (flat, contig + i, depth, reps_left[i]);

    if (check_15mer (contig + i, highmer, 15) > 0)
      uflag[i] = TRUE;
    else
      uflag[i] = FALSE;

    //for(j=1;j<=depth;j++) printf("\nAt window position %d, depth %d Found %d copys to the left %c", i,j,reps_left[i][j],*(contig+i+j-1));
    //for(j=1;j<=ld;j++) printf("\nAt window position %d, local depth %d Found %d local copys to the left %c", i,j,local_reps[i][j],*(contig+i+j-1));
  }

  // discount = 0.25;   One base less lowers score by 1/4
  // discount = 0.125;   One base less lowers score by 1/8
  // discount = 0.0625;   One base less lowers score by 1/16
  discount = 0.0625;

  for (i = 0; i <= window; i++)
  {
    fq_left[i] = 1e9;
    gc_left[i] = 2.0;
  }
  min_match = 0;
  for (i = 0; i < window - max_primer; i++)
  {
    for (primer = min_primer; primer <= max_primer; primer++)
    {
      offset = primer - ld;
      tail = primer - 15;
      tm = calc_tm (contig + i, primer);
      fq_temp = 1e8;

      if ((tm >= min_tm) && (tm <= max_tm))
      {
	fq_temp = gc_temp = 0.0;
	for (j = 0; j < primer; j++)
	{
	  rm = minim (depth, primer - j);
	  thisd = 1.0;
	  for (m = depth; m > rm; m--)
	    thisd *= discount;
	  k = i + j;
	  for (m = rm; m > min_match; m--)
	  {
	    fq_temp += thisd * reps_left[k][m];
	    thisd *= discount;
	  }

	  k = i + j;
	  if ((contig[k] == 'G') || (contig[k] == 'C'))
	    gc_temp += 1.0;

	  if (j + tail < primer)
      {
	    if (uflag[k])
        {
	      if (j + tail + 1 < primer)
          {
              fq_temp += 100;
          }
	      else
          {
              fq_temp += 2000;
          }
        }
      }
	}
	gc_temp /= (double) primer;
	fq_temp /= (double) primer;

	//printf("\n Primer = %d i = %d fq = %g  gc = %g \n\n",primer,i,fq_left[i],gc_left[i]);

	if ((gc_left[i] >= min_gc) && (gc_left[i] <= max_gc))
	{
	  if (check_hairpin (contig + i, primer))
	  {
	    /*
	     * strncpy(ss,contig+i,primer); ss[primer] = '\0'; printf("\n\t\tDetermined %s to be a hairpin %g",ss,fq_left[i]);
	     */
	    fq_temp += 1e7;
	  }
	  else
	  {
	    reverse_string (contig + i, flip, primer);
	    if (check_dimer (contig + i, flip, primer))
	    {
	      // strncpy(ss,contig+i,primer); ss[primer] = '\0'; printf("\n\t\tDetermined %s to be a self-dimer %g",ss,fq_left[i]);
	      fq_temp += 1e7;
	    }
	  }
	}
	else
	{
	  // strncpy(ss,contig+i,primer); ss[primer] = '\0'; printf("\n\t\tDetermined %s to be outside the gc window %g",ss,fq_left[i]);
	  fq_temp += 1e7;
	}

	// if(uflag[i+tail]) fq_left[i] += 1e8;

	// printf("\nLeft i = %d fq=%g  gc=%g ",i,fq_left[i],gc_left[i]);
      }
      if (fq_temp < fq_left[i])
      {
      	plen[i] = primer;
      	fq_left[i] = fq_temp;
      	gc_left[i] = gc_temp;
      }
    }
  }

  for (i = 0; i <= window; i++)
    index_left[i] = i;

  FQ_LIST = fq_left;
  qsort ((void *) index_left, window - max_primer, sizeof (int), sort_compare_index);

  free_ivector (uflag, 0, window);
  free_imatrix (reps_left, 0, window, 0, depth);
  free_imatrix (local_reps, 0, window, 0, depth);

  // for(i=0;i<window-primer;i++) { printf("\nBest Primer Left %d is site %d with FQ = %g and GC = %g", i,index_left[i],fq_left[index_left[i]],gc_left[index_left[i]]);}
}

/*---------------------------------------------------------------------*/
int
check_15mer (char *string, unsigned char *map, int n)
{
  unsigned int which;
  unsigned char mask, final;
  int i, j;

  which = encode_basepairs (string, n);

  i = which / 8;

  j = which % 8;

  mask = 1 << j;

  final = map[i] & mask;


  /*
   * printf("\nChecking %d which maps to position %d bit %d  with word %d yielding %d ", which,i,j,map[i],final);
   *
   * for(i=0;i<n;i++) printf("%c",string[i]);
   */

  return final;

}

/*---------------------------------------------------------------------*/

int
is_not_repeat (int x, int y, int **list, int no_frags, double max)
{
  int i, j, k, count;

  /* printf("\n In is_not_repeat x = %d y = %d no_frages = %d max = %g\n\n",x,y,no_frags,max); */

  count = 0;

  i = find_frag (x, list, 0, no_frags - 1, no_frags / 2);
  if (i >= 0)
    count += minim (list[i][1], y) - maxim (x, list[i][0]) + 1;

  j = find_frag (y, list, 0, no_frags - 1, no_frags / 2);
  if (j >= 0)
    count += minim (list[j][1], y) - maxim (x, list[j][0]) + 1;

  if ((i < 0) && (j < 0))
  {
    k = find_seg (x, y, list, 0, no_frags - 1, no_frags / 2);
    if (k > 0)
      count += minim (list[k][1], y) - maxim (x, list[k][0]) + 1;
  }
  if ((double) count / (double) (y - x + 1) >= max)
    return FALSE;
  else
    return TRUE;
}

/*---------------------------------------------------------------------*/
int
find_frag (int x, int **list, int start, int end, int guess)
{
  /*
   * if(PRINT_ME) printf("\n In find_frag x = %lu start = %d  end = %d  guess = %d  (%lu,%lu)", x,start,nd,guess,list[guess][0],list[guess][1]);
   */

  if (x < 0)
    return -1;

  if ((x >= list[guess][0]) && (x <= list[guess][1]))
    return guess;

  /*
   * if( (start == guess) && (end == guess) ) return -1;
   */

  if (x < list[guess][0])
  {
    if (guess <= start)
      return -1;

    end = guess - 1;
  }
  else if (x > list[guess][1])
  {
    if (guess >= end)
      return -1;

    start = guess + 1;
  }
  guess = (start + end) / 2;
  return find_frag (x, list, start, end, guess);

}

/*---------------------------------------------------------------------*/
int
find_seg (int x, int y, int **list, int start, int end, int guess)
{
  /*
   * if(PRINT_ME) printf("\n In find_frag x = %lu start = %d  end = %d  guess = %d  (%lu,%lu)", x,start,end,guess,list[guess][0],list[guess][1]);
   */

  if ((x <= list[guess][0]) && (y >= list[guess][1]))
    return guess;

  /*
   * if( (start == guess) && (end == guess) ) return -1;
   */

  if (x < list[guess][0])
  {
    if (guess <= start)
      return -1;

    end = guess - 1;
  }
  else if (x > list[guess][1])
  {
    if (guess >= end)
      return -1;

    start = guess + 1;
  }
  guess = (start + end) / 2;
  return find_seg (x, y, list, start, end, guess);

}

/*---------------------------------------------------------------------*/

int
check_hairpin (char *ss, int n)
{
  int min;
  char flip[1024];

  if (n > 1000)
  {
    printf ("\n Dude what's with a %d base primer .... PPPPLEEASE \n", n);
    exit (1);
  }
  reverse_string (ss, flip, n);

  for (min = 3; min <= n / 2 - 4; min++)
    if (check_hairpin_min (ss, flip, n, min))
      return TRUE;

  return FALSE;
}

/*---------------------------------------------------------------------*/

int
check_hairpin_min (char *ss, char *flip, int n, int min)
{
  int half;

  half = n / 2 - min;

  if (half < 4)
    return FALSE;

  return check_dimer (ss, flip, half);

}

/*---------------------------------------------------------------------*/
int
test_dimer (char *p1, char *p2, int n)
{
  int i, matches;

  if (n < 3)
    return FALSE;

  matches = 0;

  for (i = 0; i < n; i++)
    if (check_watson_crick (p1[i], p2[i]))
      matches++;

  if ((double) matches / (double) n >= 0.75)
    return TRUE;

  return FALSE;

}

/*---------------------------------------------------------------------*/
int
check_uneven_dimer (char *p1, char *p2, int n1, int n2)
{
  int i, n;

  if (n1 == n2)
    return check_dimer (p1, p2, n1);


  if ((n1 < 3) || (n2 < 3))
    return FALSE;

  n = minim (n1, n2);
  if (test_dimer (p1, p2, n))
    return TRUE;

  for (i = 1; i < n1; i++)
  {
    if (test_dimer (p1 + i, p2, minim (n1 - i, n)))
      return TRUE;
  }

  for (i = 1; i < n2; i++)
  {
    if (test_dimer (p1, p2 + i, minim (n2 - i, n)))
      return TRUE;
  }

  return FALSE;
}

/*---------------------------------------------------------------------*/
int
check_dimer (char *p1, char *p2, int n)
{
  int i;


  if (n < 3)
    return FALSE;

  if (test_dimer (p1, p2, n))
    return TRUE;


  for (i = 1; i < n; i++)
  {
    if (test_dimer (p1 + i, p2, n - i))
      return TRUE;
    if (test_dimer (p1, p2 + i, n - i))
      return TRUE;
  }

  return FALSE;

}

/*---------------------------------------------------------------------*/
int
check_watson_crick (char a, char b)
{
  if (a == 'A')
  {
    if (b == 'T')
    {
      return TRUE;
    }
    else
    {
      return FALSE;
    }
  }
  if (a == 'T')
  {
    if (b == 'A')
    {
      return TRUE;
    }
    else
    {
      return FALSE;
    }
  }
  if (a == 'G')
  {
    if (b == 'C')
    {
      return TRUE;
    }
    else
    {
      return FALSE;
    }
  }
  if (a == 'C')
  {
    if (b == 'G')
    {
      return TRUE;
    }
    else
    {
      return FALSE;
    }
  }
  // added on 06-09-2014 to silence the compiler warning that this
  // function might not return anything
  return FALSE;
}

/*---------------------------------------------------------------------*/
void
count_copys (int *flat, char *s, int n, int *reps)
{
  char ss[256];

  simple_count_copys (flat, s, n, reps, TRUE);
  reverse_transcribe (s, ss, n);
  simple_count_copys (flat, ss, n, reps, FALSE);
}

/*---------------------------------------------------------------------*/

void
simple_count_copys (int *flat, char *s, int n, int *reps, int forward)
{
  int i, j;
  unsigned int offset, k;

  offset = 0;
  for (i = 1; i <= n; i++)
  {
    if (s[i - 1] == 'N')
    {
      for (j = 1; j <= n; j++)
	reps[j] += 1e8;
      return;
    }
    if (forward)
      k = encode_basepairs (s, i);
    else
      k = encode_basepairs (s + n - i, i);

    reps[i] += flat[k + offset];

    /*
     * if(i == 1) printf("\n n = %d forward k = %d",n,k+offset);
     */

    offset++;
    offset = offset << 2;
  }

}

/*---------------------------------------------------------------------*/
void
reverse_string (char *contig, char *s, int n)
{
  int i;

  for (i = n - 1; i > -1; i--)
  {
    *s++ = *(contig + i);
  }

}

/*---------------------------------------------------------------------*/
void
reverse_transcribe (char *contig, char *s, int n)
{
  int i;
  char c;

  for (i = n - 1; i > -1; i--)
  {
    if (*(contig + i) == 'A')
      c = 'T';
    else if (*(contig + i) == 'C')
      c = 'G';
    else if (*(contig + i) == 'G')
      c = 'C';
    else if (*(contig + i) == 'T')
      c = 'A';
    else
      c = 'N';

    *s++ = c;
  }
  /*
   * *s = '\0'; printf("\n%s");
   */

}

/*---------------------------------------------------------------------*/
void
transcribe (char *contig, char *s, int n)
{
  int i;
  char c;

  for (i = 0; i < n; i++)
  {
    if (*(contig + i) == 'A')
      c = 'T';
    else if (*(contig + i) == 'C')
      c = 'G';
    else if (*(contig + i) == 'G')
      c = 'C';
    else if (*(contig + i) == 'T')
      c = 'A';
    else
      c = 'N';

    s[i] = c;
  }
}

/*---------------------------------------------------------------------*/

void
flat_index_contig (int *index, char *contig, int L, int depth)
{
  int i, stop;

  for (i = 0; i < L; i++)
  {
    stop = minim (depth, L - i);

    // if(i%100000 == 0) printf("\nIndex %d bases out of %d at depth %d",i,L,stop);

    if ((*(contig + i)) != 'N')
      index_string (index, contig + i, stop);
  }
}

/*---------------------------------------------------------------------*/
void
index_string (int *index, char *s, int n)
{
  unsigned int i, offset;
  int j;

  offset = 0;
  for (j = 1; j <= n; j++)
  {
    i = encode_basepairs (s, j);
    index[i + offset]++;
    offset++;
    offset = offset << 2;
  }
}

/*---------------------------------------------------------------------*/

int
sort_compare (const void *a, const void *b)
{
  if (*((double *) a) < *((double *) b))
    return -1;
  else if (*((double *) a) > *((double *) b))
    return 1;
  else
    return 0;
}

/*---------------------------------------------------------------------*/

int
sort_compare_index (const void *a, const void *b)
{
  double ad, bd;

  ad = FQ_LIST[*((int *) a)];
  bd = FQ_LIST[*((int *) b)];

  if (ad < bd)
    return -1;
  else if (ad > bd)
    return 1;
  else
    return 0;
}

/*---------------------------------------------------------------------*/

void
set_next (char *sss, int *i, int *j)
{

  while (!(isalnum (sss[*i]) || (sss[*i] == '(') || (sss[*i] == '*') || (sss[*i] == '-')))
    (*i)++;
  *j = *i;
  while (!isspace (sss[*j]) && (sss[*j] != '"'))
    (*j)++;
  sss[*j] = '\0';
  /* printf("\n %s",&sss[*i]); */
}

/*---------------------------------------------------------------------*/
int
pick (int n)
{
  return floor ((double) n * ran1 (&idum));
}

/*---------------------------------------------------------------------*/
#define M1 259200l
#define IA1 7141l
#define IC1 54773l
#define RM1 (1.0/M1)
#define M2 134456l
#define IA2 8121l
#define IC2 28411l
#define RM2 (1.0/M2)
#define M3 243000l
#define IA3 4561l
#define IC3 51349l

double
ran1 (int *idum)
{
  static int ix1, ix2, ix3;
  static double r[98];
  double temp;
  static int iff = 0;
  int j;

  if (*idum < 0 || iff == 0)
  {
    iff = 1;
    ix1 = (IC1 - (*idum)) % M1;
    ix1 = (IA1 * ix1 + IC1) % M1;
    ix2 = ix1 % M2;
    ix1 = (IA1 * ix1 + IC1) % M1;
    ix3 = ix1 % M3;
    for (j = 1; j <= 97; j++)
    {
      ix1 = (IA1 * ix1 + IC1) % M1;
      ix2 = (IA2 * ix2 + IC2) % M2;
      r[j] = (ix1 + ix2 * RM2) * RM1;
    }
    *idum = 1;
  }
  ix1 = (IA1 * ix1 + IC1) % M1;
  ix2 = (IA2 * ix2 + IC2) % M2;
  ix3 = (IA3 * ix3 + IC3) % M3;
  j = (int) (1 + ((97 * ix3) / M3));
  if (j > 97 || j < 1)
    log_err ("RAN1: This cannot happen.");
  temp = r[j];
  r[j] = (ix1 + ix2 * RM2) * RM1;
  return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3

/*---------------------------------------------------------------------*/
#define E 2.71828182

double
gs (double alpha)
{
  int flag = 1;
  double b, p, x;

  b = (alpha + E) / E;
  do
  {
    p = b * ran1 (&idum);
    if (p > 1.0)
    {
      x = -log ((b - p) / alpha);
      if (pow (x, alpha - 1.0) >= ran1 (&idum))
	flag = 0;
    }
    else
    {
      x = pow (p, 1.0 / alpha);
      if (x <= -log (ran1 (&idum)))
	flag = 0;
    }
  }
  while (flag);
  return x;
}

/*---------------------------------------------------------------------*/
double
gcf (double alpha)
{
  int flag = 1;
  static double c1, c2, c3, c4, c5;
  static double aprev = 0.0;
  double tmp, u1, u2, w;

  if (alpha != aprev)
  {
    c1 = alpha - 1.0;
    tmp = 1.0 / c1;
    c2 = tmp * (alpha - 1.0 / (6.0 * alpha));
    c3 = 2.0 * tmp;
    c4 = c3 + 2.0;
    if (alpha > 2.5)
      c5 = 1.0 / sqrt (alpha);
    aprev = alpha;
  }
  do
  {
    u1 = ran1 (&idum);
    u2 = ran1 (&idum);
    if (alpha > 2.5)
    {
      u1 = u2 + c5 * (1.0 - 1.86 * u1);
      if (u1 <= 0.0 || u1 >= 1.0)
	continue;
    }
    w = c2 * u2 / u1;
    if ((c3 * u1 + w + 1.0 / w < c4) || (c3 * log (u1) - log (w) + w < 1.0))
      flag = 0;
  }
  while (flag);
  return c1 * w;
}


/*---------------------------------------------------------------------*/
double
rgamma (double alpha)
{
  if (alpha == 1.0)
    return -log (ran1 (&idum));
  if (alpha < 1.0)
    return gs (alpha);
  return gcf (alpha);
}

/*---------------------------------------------------------------------*/

double
gamdev (int ia, int *idum)
{
  int j;
  double am, e, s, v1, v2, x, y;

  if (ia < 1)
    log_err ("Error in routine GAMDEV");
  if (ia < 6)
  {
    x = 1.0;
    for (j = 1; j <= ia; j++)
      x *= ran1 (idum);
    x = -log (x);
  }
  else
  {
    do
    {
      do
      {
	do
	{
	  v1 = 2.0 * ran1 (idum) - 1.0;
	  v2 = 2.0 * ran1 (idum) - 1.0;
	}
	while (v1 * v1 + v2 * v2 > 1.0);
	y = v2 / v1;
	am = ia - 1;
	s = sqrt (2.0 * am + 1.0);
	x = s * y + am;
      }
      while (x <= 0.0);
      e = (1.0 + y * y) * exp (am * log (x / am) - s * y);
    }
    while (ran1 (idum) > e);
  }
  return x;
}

/*---------------------------------------------------------------------*/

#define PI 3.141592654

double
poidev (double xm, int *idum)
{
  static double sq, alxm, g, oldm = (-1.0);
  double em, t, y;

  if (xm < 12.0)
  {
    if (xm != oldm)
    {
      oldm = xm;
      g = exp (-xm);
    }
    em = -1;
    t = 1.0;
    do
    {
      em += 1.0;
      t *= ran1 (idum);
    }
    while (t > g);
  }
  else
  {
    if (xm != oldm)
    {
      oldm = xm;
      sq = sqrt (2.0 * xm);
      alxm = log (xm);
      g = xm * alxm - gammln (xm + 1.0);
    }
    do
    {
      do
      {
	y = tan (PI * ran1 (idum));
	em = sq * y + xm;
      }
      while (em < 0.0);
      em = floor (em);
      t = 0.9 * (1.0 + y * y) * exp (em * alxm - gammln (em + 1.0) - g);
    }
    while (ran1 (idum) > t);
  }
  return em;
}

/* #undef PI */
/*---------------------------------------------------------------------*/

void
read_var (char *line, char *result)
{

  char line1[256];
  unsigned int i;

  sprintf (line1, "%s", line);
  printf ("%s", line1);
  fgets (result, 250, stdin);
  result[strlen (result) - 1] = '\0';
  /* printf("\n You entered %s which is %d characters long\n",result,strlen(result)); */
  if (outfile != stdout)
  {
    for (i = 0; i < minim (strlen (line1), 255); i++)
      if (line1[i] == '\n')
	line1[i] = '\0';
    // fprintf (outfile, "\"%s\",%s\n", line1, result);
  }
}

/*---------------------------------------------------------------------*/

double
gammln (double xx)
{
  double x, tmp, ser;
  static double cof[6] = { 76.18009173, -86.50532033, 24.01409822,
    -1.231739516, 0.120858003e-2, -0.536382e-5
  };
  int j;

  x = xx - 1.0;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log (tmp);
  ser = 1.0;
  for (j = 0; j <= 5; j++)
  {
    x += 1.0;
    ser += cof[j] / x;
  }
  return -tmp + log (2.50662827465 * ser);
}

/*---------------------------------------------------------------------*/
double
rbeta (double u, double v)
{
  double r;

  r = rgamma (u);
  return (r / (r + rgamma (v)));
}

/*---------------------------------------------------------------------*/

int
line_count (FILE * sfile)
{
  int c, nl;
  c = nl = 0;
  while ((c = getc (sfile)) != EOF)
  {
    if (c == '\n')
      nl++;
  }
  return nl;
}

/*---------------------------------------------------------------------*/
