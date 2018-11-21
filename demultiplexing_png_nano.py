import os
import subprocess


def demultiplex():
    seq_dir = '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20181120_tara_nanopore_png'

    fastq_list = [file for file in os.listdir(seq_dir) if 'fastq.gz' in file]
    site_file_name_lists = list(set([file.split('_')[0] for file in fastq_list]))

    index_list = ['BC01','BC02', 'BC03', 'BC04', 'BC05', 'BC06', 'BC07', 'BC08', 'BC09', 'BC10', 'BC11', 'BC12']

    index_dict = {'BC01' : 'AAGAAAGTTGTCGGTGTCTTTGTG', 'BC02' : 'TCGATTCCGTTTGTAGTCGTCTGT',
                  'BC03' : 'GAGTCTTGTGTCCCAGTTACCAGG', 'BC04' : 'TTCGGATTCTATCGTGTTTCCCTA',
                  'BC05' : 'CTTGTCCAGGGTTTGTGTAACCTT', 'BC06' : 'TTCTCGCAAAGGCAGAAAGTAGTC',
                  'BC07' : 'GTGTTACCGTGGGAATGAATCCTT', 'BC08' : 'TTCAGGGAACAAACCAAGTTACGT',
                  'BC09' : 'AACTAGGCACAGCGAGTCTTGGTT', 'BC10' : 'AAGCGTTGAAACCTTTGTCCTCTC',
                  'BC11' : 'GTTTCATCTATCGGAGGGAATGGA', 'BC12' : 'CAGGTAGAAAGAAGCAGAATCGGA'}



    # for each of the fite_file_name pairs
    # we want to produce a file that has been made by using a single index in that file
    # to do this using Sabre to demultiplex we need to produce a dab delimited file
    # INDEX\tforward_name\trev_name

    # we can then write that file out and use it to run the sabre command for the given site_file_name

    for sfn in site_file_name_lists:
        # we need to uncompress the files
        r1_sfn_file = '{}_r1.fastq.gz'.format(sfn)
        r2_sfn_file = '{}_r2.fastq.gz'.format(sfn)

        subprocess.run(['gunzip', r1_sfn_file])
        subprocess.run(['gunzip', r2_sfn_file])

        r1_sfn_file_ungz = '{}/{}_r1.fastq'.format(seq_dir, sfn)
        r2_sfn_file_ungz = '{}/{}_r2.fastq'.format(seq_dir,sfn)

        sfn_dir = '{}/{}'.format(seq_dir, sfn)
        os.makedirs(sfn_dir, exist_ok=True)

        tab_delim_file = []
        for index in index_list:
            r1_demultiplexed_file = '{}/{}_{}_r1.fastq'.format(sfn_dir, sfn, index)
            r2_demultiplexed_file = '{}/{}_{}_r2.fastq'.format(sfn_dir, sfn, index)
            tab_delim_file.append('\t'.join([index_dict[index], r1_demultiplexed_file, r2_demultiplexed_file]))

        # here we have the demultiplexed_tab_delim file populated
        # now write out
        sfn_dir = '{}/{}'.format(seq_dir, sfn)
        os.makedirs(sfn_dir, exist_ok=True)

        tab_delim_path = '{}/tab_delim.txt'.format(sfn_dir)
        with open(tab_delim_path, 'w') as f:
            for line in tab_delim_file:
                f.write('{}\n'.format(line))
        r1_reads_no_match_path = '{}/no_bc_match_r1.fastq'.format(sfn_dir)
        r2_reads_no_match_path = '{}/no_bc_match_r2.fastq'.format(sfn_dir)



        # now run sabre
        subprocess.run(['sabre', 'pe', '-f', r1_sfn_file_ungz, '-r', r2_sfn_file_ungz, '-b', tab_delim_path,
                        '-u', r1_reads_no_match_path, '-w', r2_reads_no_match_path])
    # now recompress all of the fastq files
    fastq_list = [file for file in os.listdir(seq_dir) if 'fastq.gz' in file]
    for fastq in fastq_list:
        subprocess.run(['gzip', fastq])

demultiplex()