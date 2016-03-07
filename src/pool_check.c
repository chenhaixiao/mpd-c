#include "mpd.h"

// TS Wingo
// created: 2013-10-28
// updated: 2015-08-02
// checks that the pools created by primer4amplicons.c
// are actually poolable. It needs the input stripped of
// lines that don't contain primers and this is done
// using p4a_2_pool_check.pl

int main (int argc, char **argv)
{
  usage(argc == 3, "pool_check <primer_file> <isPcr_filename>" );

  // initialize
  int max_primers_in_pool = 20;
  int max_primer_pairs = 1000;
  int *primers_in_pool = create_ivec (max_primer_pairs);
  PNODE ***primer_pool = primer_pool_create (max_primers_in_pool, max_primer_pairs);

  // read data
  int max_pools = read_primer_pools( argv[1], max_primer_pairs, max_primers_in_pool,
    primers_in_pool, primer_pool);

  // check pools
  Check_all_pools( max_pools, primers_in_pool, primer_pool);

  // print isPcr
  Print_isPcr( argv[2], max_pools, primers_in_pool, primer_pool );

  return 0;

error:
  return 1;
}

/*---------------------------------------------------------------------*/

int
read_primer_pools (const char *filename, int max_primer_pairs, int max_primers_in_pool,
  int *primers_in_pool, PNODE ***primer_pool)
{
  // open primer file
  FILE *primer_file;
  primer_file = fopen (filename, "r");
  check(primer_file, "cannot open primer file, '%s'", filename);

  int ppairs_count, pool_number;
  ppairs_count = pool_number = 0;

  while (!feof(primer_file))
  {
    debug("creating new temp_ppair");
    PNODE *temp_ppair = create_ppair( 1000, 50 );

    int read_line =
      fscanf (primer_file,
        "%d\t%s\t%lg\t%lg\t%s\t%lg\t%lg\t%d\t%d\t%d\t%d\t%d\t%d\t%lg\t%lg\t%s\n",
        &temp_ppair->id,
        temp_ppair->forward->sequence,
        &temp_ppair->forward->tm,
        &temp_ppair->forward->gc,
        temp_ppair->reverse->sequence,
        &temp_ppair->reverse->tm,
        &temp_ppair->reverse->gc,
        &temp_ppair->chrom,
        &temp_ppair->forward->start,
        &temp_ppair->forward->end,
        &temp_ppair->reverse->start,
        &temp_ppair->reverse->end,
        &temp_ppair->length,
        &temp_ppair->gc,
        &temp_ppair->tm,
        temp_ppair->sequence);
    check((read_line == 16), "Error processing file: %s", filename);
    check((ppairs_count < max_primer_pairs), "Out of Memory. Increase primer_pool size.");
    check((temp_ppair->id <= max_primers_in_pool), "Too many primer pairs in pool %d", pool_number);

    // temp_ppair->id => primer number within the pool
    if ( temp_ppair->id == 0)
    {
      if (ppairs_count != 0) // the first primer => pool 0
      {
        pool_number++;
      }
    }

    primer_pool[pool_number][temp_ppair->id] = temp_ppair;
    debug("assigned primer pair id '%d' to pool number '%d'", temp_ppair->id, pool_number);

    primers_in_pool[pool_number]++;
    ppairs_count++;
    debug("\n\n\tpool = %d, pool primer pairs count = %d, total primer pairs count  = %d\n",
      pool_number, primers_in_pool[pool_number], ppairs_count);
  }

  debug("returning pool_number: %d", pool_number);
  return pool_number;

error:
  exit(1);
}

/*---------------------------------------------------------------------*/

void
Print_isPcr (const char *filename, int max_pools, int *primers_in_pool, PNODE ***primer_pool)
{
  FILE *isPcr_File;
  isPcr_File = fopen (filename, "w");
  check(isPcr_File, "cannot open primer file, '%s'", filename);
  for (int i = 0; i <= max_pools; i++ )
  {
    for (int j = 0; j < primers_in_pool[i]; j++)
    {
      fprintf(isPcr_File, "pool_%d_%02d\t%s\t%s\n", (i+1), (j+1),
        primer_pool[i][j]->forward->sequence,
        primer_pool[i][j]->reverse->sequence );
    }
  }

error:
  exit(1);
}

void
Check_all_pools ( int max_pools, int *primers_in_pool, PNODE ***primer_pool )
{
  printf("%s\t%s\t%s\t%s\n",
    "Pool Number", "Primer Pair Count", "Compatable", "Comparisons");
  for (int i = 0; i <= max_pools; i++ )
  {
    if (primers_in_pool[i] > 1 )
    {
      check_poolability( primer_pool[i], primers_in_pool[i], i);
    }
    else
    {
      printf("%d\t1\tYes\t%d\n", (i + 1), 0);
    }
  }
}

/*---------------------------------------------------------------------*/

void
die ( char *message )
{
  if (errno)
  {
    perror(message);
  }
  else
  {
    printf("ERROR: %s\n", message);
  }
  exit(1);
}

/*---------------------------------------------------------------------*/

void
check_poolability (PNODE ** primer_pool, int primers_in_pool, int pool_number)
{
  int i, j, k, *poolable_count;
  char **cmat;

  poolable_count = create_ivec (primers_in_pool);
  cmat = create_cmat (primers_in_pool, primers_in_pool);

  for (i = 0; i < primers_in_pool; i++)
    for (j = i + 1; j < primers_in_pool; j++)
    {
      // printf("checking, fwd: %s rev: %s with fwd: %s rev: %s\n",
      // primer_pool[i]->forward->sequence,
      // primer_pool[i]->reverse->sequence,
      // primer_pool[j]->forward->sequence,
      // primer_pool[j]->reverse->sequence
      // );
      cmat[i][j] = is_poolable_primer (primer_pool[i], primer_pool[j]);
      cmat[j][i] = cmat[i][j];
      poolable_count[i] += cmat[i][j];
      // debug("%d", poolable_count[i]);
      poolable_count[j] += cmat[i][j];
      // debug("%d", poolable_count[j]);
    }

  k = 0;
  for (i = 0; i < primers_in_pool; i++)
    k += poolable_count[i];

  if ((primers_in_pool - 1) * primers_in_pool == k)
    printf ("%d\t%d\tYes\t%d\n", (pool_number + 1), (primers_in_pool + 1), k);
  else
    printf ("%d\t%d\tNo\t%d\n", (pool_number + 1), (primers_in_pool + 1), k);

  // print_cmat (cmat, primers_in_pool);
  // printf ("\n");

}

/*---------------------------------------------------------------------*/

void
print_cmat (char **cmat, int N)
{
  int i, j;
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
      if (i == j)
	printf (".");
      else
	printf ("%d", cmat[i][j]);

    printf ("\n");
  }
}

/*---------------------------------------------------------------------*/

int is_poolable_primer(PNODE *p1, PNODE *p2)
{
  int i, f1, r1, f2, r2;
  char flipf1[80], flipf2[80], flipr1[80], flipr2[80];

  for (i = p1->forward->start; i <= p1->reverse->end; i++)
    if ((i >= p2->forward->start) && (i <= p2->reverse->end))
    {
      log_err ("Overlapping primers");
      return FALSE;
    }

  for (i = p2->forward->start; i <= p2->reverse->end; i++)
    if ((i >= p1->forward->start) && (i <= p1->reverse->end))
    {
      log_err ("Overlapping primers");
      return FALSE;
    }

  f1 = abs (p1->forward->end - p1->forward->start) + 1;
  r1 = abs (p1->reverse->end - p1->reverse->start) + 1;
  f2 = abs (p2->forward->end - p2->forward->start) + 1;
  r2 = abs (p2->reverse->end - p2->reverse->start) + 1;

  reverse_string (p1->reverse->sequence, flipr1, r1);
  reverse_string (p2->reverse->sequence, flipr2, r2);
  reverse_string (p1->forward->sequence, flipf1, f1);
  reverse_string (p2->forward->sequence, flipf2, f2);

  if (check_uneven_dimer (p1->forward->sequence, flipf2, f1, f2))
  {
    // log_err ("fwd 1 makes dimer with fwd 2\n");
    log_err ("\t%s makes dimer with %s\n", p1->forward->sequence, p2->forward->sequence);
    return FALSE;
  }

  if (check_uneven_dimer (p1->reverse->sequence, flipr2, r1, r2))
  {
    // log_err ("rev 1 makes dimer with rev 2\n");
    log_err ("\t%s makes dimer with %s\n", p1->reverse->sequence, p2->reverse->sequence);
    return FALSE;
  }

  if (check_uneven_dimer (p1->forward->sequence, flipr2, f1, r2))
  {
    // log_err ("fwd 1 makes dimer with rev 2\n");
    log_err ("\t%s makes dimer with %s\n", p1->forward->sequence, p2->reverse->sequence);
    return FALSE;
  }

  if (check_uneven_dimer (p1->reverse->sequence, flipf2, r1, f2))
  {
    // log_err ("rev 1 makes dimer with fwd 2\n");
    log_err ("\t%s makes dimer with %s\n", p1->reverse->sequence, p2->forward->sequence);
    return FALSE;
  }

  if (abs (p1->forward->tm - p2->forward->tm) > 4)
  {
    log_err ("fwd 1 and fwd 2 have wide abs diff in Tm > 4\n");
    return FALSE;
  }

  if (abs (p1->reverse->tm - p2->reverse->tm) > 4)
  {
    log_err ("rev 1 and rev 2 have wide abs diff in Tm > 4\n");
    return FALSE;
  }

  if (abs (p1->length - p2->length) > 50)
  {
    log_err ("row difference of product > 50bp\n");
    return FALSE;
  }

  return TRUE;
}

/*---------------------------------------------------------------------*/

PNODE ***
primer_pool_create ( int max_primers_in_pool, int max_ppairs_count)
{
  PNODE ***primer_pool = (PNODE ***) malloc ((unsigned) max_ppairs_count * sizeof (PNODE **));
  check_mem(primer_pool);

  for (int i = 0; i < max_ppairs_count; i++)
  {
    primer_pool[i] = (PNODE **) malloc ((unsigned) max_primers_in_pool * sizeof (PNODE *));
    check_mem(primer_pool[i]);
  }
  return primer_pool;

error:
  exit(1);
}

/*---------------------------------------------------------------------*/

PRIMER *
create_primer (int max_primer_length)
{
  PRIMER *tn;

  tn = (PRIMER *) malloc( (unsigned) sizeof(PRIMER) );
  assert( tn != NULL );
  check_mem(tn);

  tn->start = 0;
  tn->end = 0;
  tn->sequence = create_cvec(max_primer_length);
  tn->tm = 0;
  tn->gc = 0;

  return tn;

error:
  exit(1);
}

/*---------------------------------------------------------------------*/

PNODE *
create_ppair (int max_primer_length, int max_amplicon_length)
{
  PNODE *tn;

  tn = (PNODE *) malloc( (unsigned) sizeof(PNODE) );
  assert( tn != NULL );
  check_mem(tn);

  tn->id = 0;
  tn->forward = create_primer(max_primer_length);
  tn->reverse = create_primer(max_primer_length);
  tn->sequence = create_cvec(max_amplicon_length);
  tn->length = 0;
  tn->gc = 0;
  tn->tm = 0;
  tn->chrom = 0;

  return tn;

error:
  exit(1);
}

/*---------------------------------------------------------------------*/


int
check_hairpin(char *ss, int n)
{
  int min;
  char flip[1024];

  if (n > 1000) {
  	printf("\n Dude what's with a %d base primer .... PPPPLEEASE \n", n);
  	exit(1);
  }
  reverse_string(ss, flip, n);

  for (min = 3; min <= n / 2 - 4; min++)
  	if (check_hairpin_min(ss, flip, n, min))
  		return TRUE;

  return FALSE;
}

/*---------------------------------------------------------------------*/

int
check_hairpin_min(char *ss, char *flip, int n, int min)
{
  int half;

  half = n / 2 - min;

  if (half < 4)
  	return FALSE;

  return check_dimer(ss, flip, half);

}
/*---------------------------------------------------------------------*/

int
test_dimer(char *p1, char *p2, int n)
{
  int i, matches;

  if (n < 3)
  	return FALSE;

  matches = 0;

  for (i = 0; i < n; i++)
  	if (check_watson_crick(p1[i], p2[i]))
  		matches++;

  if ((double)matches / (double)n >= 0.75)
  	return TRUE;

  return FALSE;

}
/*---------------------------------------------------------------------*/

int
check_uneven_dimer(char *p1, char *p2, int n1, int n2)
{
  int i , n;

  if (n1 == n2)
  	return check_dimer(p1, p2, n1);


  if ((n1 < 3) || (n2 < 3))
  	return FALSE;

  n = minim(n1, n2);
  if (test_dimer(p1, p2, n))
  	return TRUE;

  for (i = 1; i < n1; i++) {
  	if (test_dimer(p1 + i, p2, minim(n1 - i, n)))
  		return TRUE;
  }

  for (i = 1; i < n2; i++) {
  	if (test_dimer(p1, p2 + i, minim(n2 - i, n)))
  		return TRUE;
  }

  return FALSE;
}

/*---------------------------------------------------------------------*/

int
check_dimer(char *p1, char *p2, int n)
{
  int i;


  if (n < 3)
  	return FALSE;

  if (test_dimer(p1, p2, n))
  	return TRUE;


  for (i = 1; i < n; i++) {
  	if (test_dimer(p1 + i, p2, n - i))
  		return TRUE;
  	if (test_dimer(p1, p2 + i, n - i))
  		return TRUE;
  }

  return FALSE;
}

/*---------------------------------------------------------------------*/

int
check_watson_crick(char a, char b)
{
  if (a == 'A') {
  	if (b == 'T') {
  		return TRUE;
  	} else {
  		return FALSE;
  	}
  }
  if (a == 'T') {
  	if (b == 'A') {
  		return TRUE;
  	} else {
  		return FALSE;
  	}
  }
  if (a == 'G') {
  	if (b == 'C') {
  		return TRUE;
  	} else {
  		return FALSE;
  	}
  }
  if (a == 'C') {
  	if (b == 'G') {
  		return TRUE;
  	} else {
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

int *
create_ivec (int row)
{
  int *v = (int *) malloc ((unsigned) row * sizeof (int));
  check_mem(v);
  for (int i = 0; i < row; i++ )
  {
    v[i] = 0;
  }
  return v;

error:
  exit(1);
}

/*---------------------------------------------------------------------*/

void
free_ivec (int *v, int row)
{
  free (v + row);
}

/*---------------------------------------------------------------------*/

char *
create_cvec (int row)
{
  char *v = (char *) malloc ((unsigned) row * sizeof (char));
  check_mem(v);
  for (int i = 0; i < row; i++ )
  {
    v[i] = 0;
  }
  return v;

error:
  exit(1);
}

/*---------------------------------------------------------------------*/

void
free_cvec (int *v, int row)
{
  free (v + row);
}

/*---------------------------------------------------------------------*/

int **
create_imat (int row, int col)
{
  int **m = (int **) malloc ((unsigned) row * sizeof (int *));
  check_mem(m);

  for (int i = 0; i < col; i++)
  {
    m[i] = (int *) malloc ((unsigned) col * sizeof (int));
    check_mem(m[i]);
    for (int j = 0; j < col; j++ )
    {
      m[i][j] = 0;
    }
  }
  return m;

error:
  exit(1);
}
/*---------------------------------------------------------------------*/

void
free_imat (int **m, int row, int col)
{
  for (int i = 0; i <= col; i++)
    free ((void *) (m[i] + row));
  free ((void *) (m + col));
}


/*---------------------------------------------------------------------*/

char **
create_cmat (int row, int col)
{
  char **m = (char **) malloc ((unsigned) row * sizeof (char *));
  check_mem(m);

  for (int i = 0; i < col; i++)
  {
    m[i] = (char *) malloc ((unsigned) col * sizeof (char));
    check_mem(m[i]);
    for (int j = 0; j < col; j++)
    {
      m[i][j] = 0;
    }
  }
  return m;

error:
  exit(1);
}

/*---------------------------------------------------------------------*/

void
free_cmat (char **m, int row, int col)
{
  for (int i = 0; i <= col; i++)
    free ((void *) (m[i] + row));
  free ((void *) (m + col));
}


/*---------------------------------------------------------------------*/

int *
ivector (int nl, int nh)
{
  int *v;

  v = (int *) malloc ((unsigned) (nh - nl + 1) * sizeof (int));
  if (!v) die ("allocation failure in ivector()");
  return v - nl;
}

/*---------------------------------------------------------------------*/

char **
cmatrix (int nrl, int nrh, int ncl, int nch)
{
  int i;
  char **m;

  m = (char **) malloc ((unsigned) (nrh - nrl + 1) * sizeof (char *));
  if (!m) die ("allocation failure 1 in cmatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)
  {
    m[i] = (char *) malloc ((unsigned) (nch - ncl + 1) * sizeof (char));
    if (!m[i]) die ("allocation failure 2 in cmatrix()");
    m[i] -= ncl;
  }
  return m;
}

/*---------------------------------------------------------------------*/

void
free_cmatrix (char **m, int nrl, int nrh, int ncl, int nch)
{
  int i;

  for (i = nrh; i >= nrl; i--)
    free ((void *) (m[i] + ncl));
  free ((void *) (m + nrl));
}

/*---------------------------------------------------------------------*/

void
free_ivector (int *v, int nl, int nh)
{
  free ((void *) (v + nl));
}

/*---------------------------------------------------------------------*/

char *
cvector (int nl, int nh)
{
  char *v;
  v = (char *) malloc ((unsigned) (nh - nl + 1) * sizeof (char));
  if (!v) die ("allocation failure in cvector()");
  return v - nl;
}
