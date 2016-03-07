#define _GNU_SOURCE
#define TRUE 1
#define FALSE 0
#define minim(a,b) ((a<b)?a:b)
#define maxim(a,b) ((a>b)?a:b)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <ctype.h>
#include "dbg.h"


// TS Wingo
// 10-30-2013
// checks a primer is compat with another primer

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

char *cvector (int nl, int nh);
int *ivector (int nl, int nh);
char **cmatrix (int nrl, int nrh, int ncl, int nch);
void free_ivector (int *v, int nl, int nh);
void free_cmatrix (char **m, int nrl, int nrh, int ncl, int nch);
void read_var (char *line, char *result);
int check_hairpin (char *ss, int n);
int check_hairpin_min (char *ss, char *flip, int n, int min);
int check_dimer (char *p1, char *p2, int n);
int check_uneven_dimer (char *p1, char *p2, int n1, int n2);
int check_watson_crick (char a, char b);
int test_dimer (char *p1, char *p2, int n);
void reverse_string (char *contig, char *s, int n);
int is_poolable_primer (PNODE * p1, PNODE * p2);
void print_cmat (char **cmat, int N);
PRIMER *primer_alloc (void);
PNODE *product_alloc (void);
static FILE *outfile;


int
main (int argc, char **argv)
{
  char fp1[80], fp2[80], rp1[80], rp2[80];
  int i, f1, r1, f2, r2;
  char flipf1[80], flipf2[80], flipr1[80], flipr2[80];

  if (argc != 5)
  {
    printf("\nUsage: primer_compat forward_primer_1 reverse_primer_1 forward_primer_2 reverse_primer_2\n");
    exit(1);
  }

  sprintf(fp1, "%s", argv[1]);
  sprintf(rp1, "%s", argv[2]);
  sprintf(fp2, "%s", argv[3]);
  sprintf(rp2, "%s", argv[4]);

  printf("\nprimer 1 fwd: %s\trev: %s\nprimer 2 fwd: %s\trev: %s\n", fp1, rp1, fp2, rp2);

  f1 = strlen(fp1);
  r1 = strlen(rp1);
  f2 = strlen(fp2);
  r2 = strlen(rp2);

  for(i=0;i<80;i++)
	flipf1[i] = flipf2[i] = flipr1[i] = flipr2[i] = '\0';

  printf("\n\nlen primer 1 fwd: %d\trev: %d\nlen primer 2 fwd: %d\trev: %d\n", f1, r1, f2, r2);
  reverse_string (rp1, flipr1, r1);
  reverse_string (rp2, flipr2, r2);
  reverse_string (fp1, flipf1, f1);
  reverse_string (fp2, flipf2, f2);

  printf("\n\ncomplement of primer 1 fwd: %s\trev: %s\ncomplement of primer 2 fwd: %s\trev: %s\n", flipf1, flipr1, flipf2, flipr2);

  if (check_uneven_dimer (fp1, flipf2, f1, f2))
    printf ("\nprimer 1 forward (%s) makes dimer with primer 2 forward (%s)\n", fp1, fp2);
  else if (check_uneven_dimer (rp1, flipr2, r1, r2))
    printf ("\nprimer 1 reverse (%s) makes dimer with primer 2 reverse (%s)\n", rp1, rp2);
  else if (check_uneven_dimer (fp1, flipr2, f1, r2))
    printf ("\nprimer 1 forward (%s) makes dimer with primer 2 reverse (%s)\n", fp1, rp2);
  else if (check_uneven_dimer (rp1, flipf2, r1, f2))
    printf ("\nprimer 1 reverse (%s) makes dimer with primer 2 forward (%s)\n", rp1, fp2);
  else
    printf("\nprimer pair 1 and 2 seem compatable.\n");

}

/*---------------------------------------------------------------------*/

int is_poolable_primer(PNODE *p1, PNODE *p2)
{
  int i, f1, r1, f2, r2;
  char flipf1[80], flipf2[80], flipr1[80], flipr2[80];

  for (i = p1->forward->start; i <= p1->reverse->end; i++)
    if ((i >= p2->forward->start) && (i <= p2->reverse->end))
    {
      printf ("Overlapping primers\n");
      return FALSE;
    }

  for (i = p2->forward->start; i <= p2->reverse->end; i++)
    if ((i >= p1->forward->start) && (i <= p1->reverse->end))
    {
      printf ("Overlapping primers\n");
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
    printf ("fwd 1 makes dimer with fwd 2\n");
    printf ("\t%s makes dimer with %s\n", p1->forward->sequence, p2->forward->sequence);
    return FALSE;
  }

  if (check_uneven_dimer (p1->reverse->sequence, flipr2, r1, r2))
  {
    printf ("rev 1 makes dimer with rev 2\n");
    printf ("\t%s makes dimer with %s\n", p1->reverse->sequence, p2->reverse->sequence);
    return FALSE;
  }

  if (check_uneven_dimer (p1->forward->sequence, flipr2, f1, r2))
  {
    printf ("fwd 1 makes dimer with rev 2\n");
    printf ("\t%s makes dimer with %s\n", p1->forward->sequence, p2->reverse->sequence);
    return FALSE;
  }

  if (check_uneven_dimer (p1->reverse->sequence, flipf2, r1, f2))
  {
    printf ("rev 1 makes dimer with fwd 2\n");
    printf ("\t%s makes dimer with %s\n", p1->reverse->sequence, p2->forward->sequence);
    return FALSE;
  }

  if (fabs (p1->forward->tm - p2->forward->tm) > 4)
  {
    printf ("fwd 1 and fwd 2 have wide abs diff in Tm > 4\n");
    return FALSE;
  }

  if (fabs (p1->reverse->tm - p2->reverse->tm) > 4)
  {
    printf ("rev 1 and rev 2 have wide abs diff in Tm > 4\n");
    return FALSE;
  }

  if (abs (p1->length - p2->length) > 50)
  {
    printf ("length difference of product > 50bp\n");
    return FALSE;
  }

  return TRUE;
}
/*---------------------------------------------------------------------*/
/*
int is_poolable_primer(PNODE *p1, PNODE *p2)
{
    int i, f1, r1, f2, r2;
    char flipf1[80],flipf2[80],flipr1[80],flipr2[80];

    for(i=p1->forward->start;i<=p1->reverse->end;i++)
        if( (i >= p2->forward->start) && (i <= p2->reverse->end))
            return FALSE;

    for(i=p2->forward->start;i<=p2->reverse->end;i++)
        if( (i >= p1->forward->start) && (i <= p1->reverse->end))
            return FALSE;

    f1 = abs(p1->forward->end - p1->forward->start) + 1;
    r1 = abs(p1->reverse->end - p1->reverse->start) + 1;
    f2 = abs(p2->forward->end - p2->forward->start) + 1;
    r2 = abs(p2->reverse->end - p2->reverse->start) + 1;

  reverse_string(p1->reverse->sequence,flipr1,r1);
  reverse_string(p2->reverse->sequence,flipr2,r2);
  reverse_string(p1->forward->sequence,flipf1,f1);
  reverse_string(p2->forward->sequence,flipf2,f2);
  
  if(check_uneven_dimer(p1->forward->sequence,flipf2,f1,f2))
      return FALSE;
  
  if(check_uneven_dimer(p1->reverse->sequence,flipr2,r1,r2))
      return FALSE;
  
  if(check_uneven_dimer(p1->forward->sequence,flipr2,f1,r2))
      return FALSE;
  
  if(check_uneven_dimer(p1->reverse->sequence,flipf2,r1,f2))
      return FALSE;

    if (abs(p1->forward->tm - p2->forward->tm) > 4)
        return FALSE;

    if (abs(p1->reverse->tm - p2->reverse->tm) > 4)
        return FALSE;

    if (abs(p1->length - p2->length) > 50)
        return FALSE;

  return TRUE;
} 
*/
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


int
check_hairpin(char *ss, int n)
{
  int		min;
  char		flip      [1024];

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
  int		half;

  half = n / 2 - min;

  if (half < 4)
  	return FALSE;

  return check_dimer(ss, flip, half);

}
/*---------------------------------------------------------------------*/

int
test_dimer(char *p1, char *p2, int n)
{
  int		i         , matches;

  if (n < 3)
  	return FALSE;

  matches = 0;

  for (i = 0; i < n; i++)
  	if (check_watson_crick(p1[i], p2[i]))
  		matches++;

  if ((double)matches / (double)n >= 0.75)
  {
	char s1[80];
	char s2[80];
	s1[n] = s2[n] = '\0';
	strncpy(s1,p1,n);
	strncpy(s2,p2,n);
	printf("\n We are calling a %d out %d match between %s and %s \n",matches,n,s1,s2);
  	return TRUE;
  }

  return FALSE;

}
/*---------------------------------------------------------------------*/

int
check_uneven_dimer(char *p1, char *p2, int n1, int n2)
{
  int		i         , n;

  if (n1 == n2)
  	return check_dimer(p1, p2, n1);


  if ((n1 < 8) || (n2 < 8))
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
  int		i;


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
  // printf("\n In here with %s going to %s for %d bytes\n",contig,s,n);

  for (i = n - 1; i > -1; i--)
  { 
    *(s++) = *(contig + i);
  }
  // printf("\n Finished with %s going from %s for %d bytes\n",contig,s,n);

}

/*---------------------------------------------------------------------*/

void
read_var (char *line, char *result)
{

  char line1[256];
  int i;

  sprintf (line1, "%s", line);
  fgets (result, 250, stdin);
  result[strlen (result) - 1] = '\0';
  /* printf("\n You entered %s which is %d characters long\n",result,strlen(result)); */
  if (outfile != stdout)
  {
    for (i = 0; i < minim (strlen (line1), 255); i++)
      if (line1[i] == '\n')
	line1[i] = '\0';
    fprintf (outfile, "\"%s\",%s\n", line1, result);
  }
}

/*---------------------------------------------------------------------*/

int *
ivector (int nl, int nh)
{
  int *v;

  v = (int *) malloc ((unsigned) (nh - nl + 1) * sizeof (int));
  if (!v)
    log_err ("allocation failure in ivector()");
  return v - nl;
}

/*---------------------------------------------------------------------*/

char **
cmatrix (int nrl, int nrh, int ncl, int nch)
{
  int i;
  char **m;

  m = (char **) malloc ((unsigned) (nrh - nrl + 1) * sizeof (char *));
  if (!m)
    log_err ("allocation failure 1 in cmatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)
  {
    m[i] = (char *) malloc ((unsigned) (nch - ncl + 1) * sizeof (char));
    if (!m[i])
      log_err ("allocation failure 2 in cmatrix()");
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
  if (!v)
    log_err ("allocation failure in cvector()");
  return v - nl;
}
