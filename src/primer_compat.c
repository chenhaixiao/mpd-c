#include "mpd.h"

// TS Wingo
// 10-30-2013
// checks a primer is compat with another primer

int main (int argc, char **argv)
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
