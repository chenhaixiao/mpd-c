#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include "dbg.h"
#include "mem.h"

#define TRUE 1
#define FALSE 0

#define minim(a,b) ((a<b)?a:b)
#define maxim(a,b) ((a>b)?a:b)

void read_var (char *line, char *result);
double ran1 (int *idum);
int pick (int n);
void set_next (char *sss, int *i, int *j);
void flat_index_contig (int *index, char *contig, int L, int depth, int high_depth, unsigned char *double_high);
void high_index (unsigned char *dhigh, char *s, int n);
void reverse_transcribe (char *contig, char *s, int n);
void reverse_string (char *contig, char *s, int n);
int sort_compare_index (const void *a, const void *b);
int check_watson_crick (char a, char b);
unsigned int encode_basepairs (char *ss, int n);
void decode_basepairs (unsigned char *s, char *dest, int n);
void index_string (int *index, char *s, int n);
void convert_int_basepairs (int i, char *s, int k);

static FILE *outfile, **innfile;
static int idum;

int
main ()
{
  char **filename, ss[256], sss[4196], basename[1024];
  char *scratch_pad, **contig_descript;
  unsigned char **compressed_map, *double_high;
  int i, j, k, N, max_sites, *contig_length, not_done, newpos, high_size, high_depth;
  int fasta, idepth, total_index, *repeat_count;
  int *flat_index, in_repeat;
  FILE *sfile, *cfile, *idfile, *rpfile, *highfile;


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

  read_var ("Number of Contig Fasta Files to Process\n", ss);
  N = atoi (ss);
  filename = cmatrix (0, N, 0, 256);
  repeat_count = ivector (0, N);

  innfile = (FILE **) malloc ((unsigned) (N + 1) * sizeof (FILE *));
  if (!innfile)
    log_err ("allocation failure for innfile");

  compressed_map = (unsigned char **) malloc ((unsigned) (N + 1) * sizeof (unsigned char *));
  if (!compressed_map)
    log_err ("allocation failure for compressed_map");

  contig_descript = cmatrix (0, N, 0, 4196);
  contig_length = ivector (0, N);

  for (i = 0; i < N; i++)
  {
    sprintf (sss, "Please Enter Name For Contig Fasta File %d\n", i + 1);
    read_var (sss, filename[i]);
    if ((innfile[i] = fopen (filename[i], "r")) == (FILE *) NULL)
    {
      printf ("\n Can not open file %s\n", filename[i]);
      exit (1);
    }
  }
  read_var ("Basename to save compressed Genome and Indexes\n", basename);

  read_var ("Maximum Number Of Sites in a contig\n", ss);
  max_sites = atoi (ss);

  read_var ("Index Depth\n", ss);
  idepth = atoi (ss);
  sprintf (sss, "%s.sdx", basename);
  if ((sfile = fopen (sss, "w")) == (FILE *) NULL)
  {
    printf ("\nCould Not Open file %s\n", sss);
    exit (1);
  }
  fprintf (sfile, "%d\n", N);

  j = 4;
  for (i = 0; i < idepth; i++)
    j *= 4;
  j = (j - 4) / 3;

  total_index = j;
  printf ("\n Determined the size of flat index to be %d * %ld = %ld bytes\n",
	  total_index, sizeof (int), total_index * sizeof (int));

  flat_index = ivector (0, total_index);

  for (i = 0; i <= total_index; i++)
    flat_index[i] = 0;

  scratch_pad = cvector (0, max_sites);
  not_done = TRUE;

  sprintf (sss, "%s.cdx", basename);
  if ((cfile = fopen (sss, "w")) == (FILE *) NULL)
  {
    printf ("\nCould Not Open file %s\n", sss);
    exit (1);
  }
  sprintf (sss, "%s.rdx", basename);
  if ((rpfile = fopen (sss, "w")) == (FILE *) NULL)
  {
    printf ("\nCould Not Open file %s\n", sss);
    exit (1);
  }
  sprintf (sss, "%s.15x", basename);
  if ((highfile = fopen (sss, "w")) == (FILE *) NULL)
  {
    printf ("\nCould Not Open file %s\n", sss);
    exit (1);
  }
  high_depth = 15;
  high_size = 1;
  for (i = 0; i < high_depth; i++)
    high_size *= 4;
  high_size /= 8;

  double_high = ucvector (0, high_size);
  for (i = 0; i <= high_size; i++)
    double_high[i] = 0;

  for (fasta = 0; fasta < N; fasta++)
  {
    not_done = TRUE;

    in_repeat = FALSE;
    newpos = 0;

    fgets (sss, 4195, innfile[fasta]);
    j = strlen (sss);
    for (i = 0; i <= j; i++)
      if ((sss[i] == '\n') || (sss[i] == 10) || (sss[i] == 13))
	contig_descript[fasta][i] = ' ';
      else
	contig_descript[fasta][i] = sss[i];


    contig_descript[fasta][i] = '\0';

    fgets (sss, 256, innfile[fasta]);
    repeat_count[fasta] = 0;

    while (not_done)
    {
      /* printf("%s",sss);   */
      j = strlen (sss);

      for (i = 0; i < j; i++)
      {
	if (isalpha (sss[i]))
	{
	  scratch_pad[newpos] = toupper (sss[i]);
	  if (islower (sss[i]) || (sss[i] == 'N'))
	  {
	    if (!in_repeat)
	    {
	      in_repeat = TRUE;
	      repeat_count[fasta]++;
	      fwrite (&newpos, sizeof (int), 1, rpfile);
	    }
	  }
	  else if (in_repeat)
	  {
	    in_repeat = FALSE;
	    k = newpos - 1;
	    fwrite (&k, sizeof (int), 1, rpfile);
	  }
	  newpos++;
	}
      }

      if (feof (innfile[fasta]) != 0)
	not_done = FALSE;
      else
	not_done = TRUE;

      if (not_done)
      {
	/* sprintf(sss, "");                   what is this doing? */
	fgets (sss, 256, innfile[fasta]);
	if (strlen (sss) < 1)
	  not_done = FALSE;
      }
    }
    contig_length[fasta] = newpos;
    if (in_repeat)
      fwrite (&newpos, sizeof (int), 1, rpfile);

    printf ("\n Finished reading fasta %d \n\n", fasta);

    for (i = 0; i < 4; i++)
      scratch_pad[newpos + i] = 'A';	/* Pad with A's */
    printf ("\n\t\tBeginning to index contig %s which has length %d\n\n", contig_descript[fasta], contig_length[fasta]);
    flat_index_contig (flat_index, scratch_pad, contig_length[fasta], idepth, high_depth, double_high);
    printf ("\nFinished indexing contig %d \n\n", fasta);

    i = contig_length[fasta] % 4;
    int ts;
    if (i == 0)
      ts = newpos / 4;
    else
      ts = newpos / 4 + 1;

    compressed_map[fasta] = ucvector (0, ts);
    j = 0;
    for (i = 0; i < contig_length[fasta]; i += 4, j++)
      compressed_map[fasta][j] = (unsigned char) encode_basepairs (&scratch_pad[i], 4);
    fwrite (compressed_map[fasta], sizeof (unsigned char), j, cfile);
    free_ucvector (compressed_map[fasta], 0, ts);
  }
  printf ("\n Finishing up now \n\n");
  fclose (cfile);
  fwrite (double_high, sizeof (unsigned char), high_size, highfile);
  fclose (highfile);
  for (i = 0; i < N; i++)
    fprintf (sfile, "%d\t%d\t%s\n", contig_length[i], repeat_count[i], contig_descript[i]);
  fprintf (sfile, "%d\n", idepth);

  sprintf (sss, "%s.idx", basename);
  if ((idfile = fopen (sss, "w")) == (FILE *) NULL)
  {
    printf ("\nCould Not Open file %s\n", sss);
    exit (1);
  }
  fwrite (flat_index, sizeof (int), total_index, idfile);
  fclose (idfile);

  fprintf (sfile, "%s.cdx\n", basename);
  fprintf (sfile, "%s.idx\n", basename);
  fprintf (sfile, "%s.rdx\n", basename);
  fprintf (sfile, "%s.15x\n", basename);
  fclose (sfile);
  printf ("\n Finished reading fasta got here.\n");
  return (0);
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

    /* printf("\n Decoding %d ",j); */
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

    /* printf("as %s",ss); */
  }
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
    if (ss[i] == 'C')
      k++;
    else if (ss[i] == 'G')
      k += 2;
    else if (ss[i] == 'T')
      k += 3;
    else if (ss[i] != 'A')
      return -1;

  }

  /*
   * printf("\nn = %d Encoding %c%c%c%c as %d ",n,ss[0],ss[1],ss[2],ss[3],k);
   */


  return k;
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
  return FALSE;			
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

}

/*---------------------------------------------------------------------*/

void
flat_index_contig (int *index, char *contig, int L, int depth, int high_depth, unsigned char *double_high)
{
  int i, stop;
  char flip[256];


  for (i = 0; i < L; i++)
  {
    stop = minim (depth, L - i);
    if (i % 1000000 == 0)
      printf ("\nIndex %d bases out of %d at depth %d", i, L, stop);

    if ((*(contig + i)) != 'N')
    {
      index_string (index, contig + i, stop);
      if (i + high_depth < L)
      {
	high_index (double_high, contig + i, high_depth);
	reverse_transcribe (contig + i, flip, high_depth);
	high_index (double_high, flip, high_depth);
      }
    }
  }


}

/*---------------------------------------------------------------------*/
void
high_index (unsigned char *dhigh, char *s, int n)
{
  int i, j;
  unsigned char k, bit;

  i = encode_basepairs (s, n);
  /* printf("\nEntered High Index n = %d  i = %u\n\n",n,i); */
  if (i >= 0)
  {
    j = i / 8;
    k = i % 8;
    bit = 1 << k;
    dhigh[j] = dhigh[j] | bit;
  }
}

/*---------------------------------------------------------------------*/
void
index_string (int *index, char *s, int n)
{
  unsigned int offset;
  int i, j;

  offset = 0;
  for (j = 1; j <= n; j++)
  {
    i = encode_basepairs (s, j);
    if (i >= 0)
    {
      index[i + offset]++;
      offset++;
      offset = offset << 2;
    }
    else
      return;
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

void
read_var (char *line, char *result)
{

  char line1[256];
  int i;

  sprintf (line1, "%s", line);
  printf("%s",line1);
  fgets (result, 250, stdin);
  result[strlen (result) - 1] = '\0';
  /*
   * printf("\n You entered %s which is %d characters long\n",result,strlen(result));
   */
  if (outfile != stdout)
  {
    for (i = 0; i < minim (strlen (line1), 255); i++)
      if (line1[i] == '\n')
	line1[i] = '\0';
    fprintf (outfile, "\"%s\",%s\n", line1, result);
  }
}

/*---------------------------------------------------------------------*/
