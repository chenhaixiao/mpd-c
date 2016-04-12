#include "mpd.h"

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
        flat_index_contig_high(flat_index, scratch_pad, contig_length[fasta], idepth, high_depth, double_high);
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
