#include "mpd.h"

/*---------------------------------------------------------------------*/
void
make_greedy_pools (FILE *outfile, PNODE **plist, char **cmat, int *pc, int **redund, int *bstart, int Nregs, int N, int still, int *current_pool, int this_pool, int max_pool)
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

  make_greedy_pools (outfile, plist, cmat, pc, redund, bstart, Nregs, N, still, current_pool, this_pool, max_pool);
}

/*---------------------------------------------------------------------*/
void
make_less_greedy_pools (FILE *outfile, PNODE **plist, char **cmat, int *pc, int **redund, int *bstart, int Nregs, int N, int still,
	     int *current_pool, int this_pool, int max_pool)
{
  int i, j;
  if (still < 1)
    return;

  int this_p = -1;
  if (this_pool == 0)
  {
    // fprintf (outfile, "\nStarting a New Pool\n");
    int best = 100000000;
    for (i = 0; i < Nregs; i++)
      if (bstart[i] >= 0)
				if (pc[bstart[i]] < best)
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

  make_less_greedy_pools (outfile, plist, cmat, pc, redund, bstart, Nregs, N, still, current_pool, this_pool, max_pool);
}
/*---------------------------------------------------------------------*/

int
is_poolable_primer (PNODE * p1, PNODE * p2, int size_diff_threshold, int tm_diff_threshold)
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

  if (fabs (p1->forward->tm - p2->forward->tm) > tm_diff_threshold)
    return FALSE;

  if (fabs (p1->reverse->tm - p2->reverse->tm) > tm_diff_threshold)
    return FALSE;

  if (abs (p1->length - p2->length) > size_diff_threshold)
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
    sscanf (sss, "%s\t%d\t%d\t%d\t%s\t%c/%c", temp_snp[i]->name, &temp_snp[i]->no_disc, &temp_snp[i]->chrom, &temp_snp[i]->pos, s,	/* het frequency */
	    &temp_snp[i]->baseA, &temp_snp[i]->baseB);

    temp_snp[i]->pos--;
    temp_snp[i]->het = (double) atof (s);
    temp_snp[i]->no_pairs = 0;

    if ((isbase (temp_snp[i]->baseA)) && (isbase (temp_snp[i]->baseB))
	&& (temp_snp[i]->chrom == chrom) && (temp_snp[i]->het > 0.01))
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
    
    /*
    printf("\nForward sequence\n\n");
    printf("\n In tile contig with L = %d; amp_max = %d;  amp_min = %d; depth = %d\n\n", L, amp_max, amp_min, depth);
    for (i = 0; i < L; i++)
    {
        printf("%c", contig[i]);
        //if ((i + 1 == target_base) || (i == target_base)) printf(" * ");
        if (i % 80 == 79) printf("\n");
    }
    printf("\n");
    */
    
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
    
    // printf("\n About to fill quality in forward direction\n\n");
    
    fill_quality_scores (flat, local_index, contig, L, minim (target_base - 30, end_region), depth, local_depth,
                         min_primer, max_primer, fq_left, gc_left, index_left, plen_l,
                         min_gc, max_gc, min_tm, max_tm, repeats, no_repeats, highmer, start_pos);
    rt_contig = cvector (0, L);
    reverse_transcribe (contig, rt_contig, L);
    
    /*
    printf("\nReverse Transcribe\n\n");
    printf("\n");for (i = 0; i < L; i++) {printf("%c", rt_contig[i]);
    if ((i + 1 == L - target_base) || (i + 2 == L - target_base)) printf(" * ");
    if (i % 80 == 79)
        printf("\n");
    }
    printf("\n");
    */
    
    printf ("\n About to fill quality in reverse direction\n");
    
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
            printf("\nMatching 5' %s 3' (fq = %g, gc = %g, len = %d  tm = %g) with", ss, fq_left[index_left[i]], gc_left[index_left[i]], plen_l[index_left[i]], tm_l);
            
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
        printf("\n\n **** Have completely failed to find a primer pair *******\n\n");
        
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
    printf("\nA primer which goes from %d to %d appears to have %s under it at pos %d\n\n", p_start,p_end,list[which]->name,list[which]->pos);
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
read_var (char *line, char *result)
{

  char line1[256];
  unsigned int i;

  sprintf (line1, "%s", line);
  printf ("%s", line1);
  fgets (result, 250, stdin);
  result[strlen (result) - 1] = '\0';
  for (i = 0; i < minim (strlen (line1), 255); i++)
    if (line1[i] == '\n')
	    line1[i] = '\0';
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
    int this_primer_num;

    int read_line =
      fscanf (primer_file,
        "%d\t%s\t%lg\t%lg\t%s\t%lg\t%lg\t%d\t%d\t%d\t%d\t%d\t%d\t%lg\t%lg\t%s\n",
        &this_primer_num,
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
    check((this_primer_num < max_primers_in_pool), "Too many primer pairs in pool %d", pool_number);

    // temp_ppair->id => primer number within the pool
    if ( this_primer_num == 0)
    {
      if (ppairs_count != 0) // the first primer => pool 0
      {
        pool_number++;
      }
    }

    primer_pool[pool_number][this_primer_num] = temp_ppair;
    debug("assigned primer pair id '%d' to pool number '%d'", this_primer_num, pool_number);

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
  printf("%s\t%s\t%s\t%s\n", "Pool Number", "Primer Pair Count", "Compatable", "Comparisons");
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
      // debug("checking, fwd: %s rev: %s with fwd: %s rev: %s\n",
      // primer_pool[i]->forward->sequence,
      // primer_pool[i]->reverse->sequence,
      // primer_pool[j]->forward->sequence,
      // primer_pool[j]->reverse->sequence
      // );
      
      // this is not ideal, but using 50bp and 2C as for compatibility
      cmat[i][j] = is_poolable_primer (primer_pool[i], primer_pool[j], 50, 2);
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

void
flat_index_contig_high (int *index, char *contig, int L, int depth, int high_depth, unsigned char *double_high)
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
