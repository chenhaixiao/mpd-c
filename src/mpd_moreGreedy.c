#include "mpd.h"
#include "mem.h"

static FILE *outfile;

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
                found_count = find_primers (snp_list[target_contig],
                                            loop_amp,
                                            contig_snp_count[target_contig],	// number of snps in the contig
                                            flat_index,
                                            scratch_pad,	// copy of the target
                                            j,	// length of region
                                            min_primer, max_primer, amp_max, amp_min, min_gc, max_gc, this_min_tm, this_max_tm, 
                                            idepth,	// index depth
                                            10,	// local depth
                                            j / 2,	// target base (from old primer_snp.c program)
                                            genome_start + 1,
                                            highmer,
                                            repeat_list[target_contig],
                                            no_repeats[target_contig],
                                            pad_size,	// size on either end of the contig to look for primer
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
            double max_amp_diff = (double) (amp_max * 0.15) + 1;
            poolable_matrix[i][j] = is_poolable_primer (all_primer_pairs[i], all_primer_pairs[j], (int) max_amp_diff, 2);
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
    make_greedy_pools ( outfile, all_primer_pairs, poolable_matrix, poolable_count, redundant_list, best_start, amp_pool_count, primer_count, primer_count, current_pool, 0, pool_size);
    return 0;
}
