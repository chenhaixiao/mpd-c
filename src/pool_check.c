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

