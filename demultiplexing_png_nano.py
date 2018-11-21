import os
import subprocess
import regex
import pickle
from multiprocessing import Queue, Process

def demultiplex_from_scratch():
    seq_dir = '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20181120_tara_nanopore_png'

    site_file_name_list = ['KimbeS1', 'KimbeS2', 'KimbeS3', 'KimbeS4', 'KimbeS4bis']

    index_list = ['BC01', 'BC02', 'BC03', 'BC04', 'BC05', 'BC06', 'BC07', 'BC08', 'BC09', 'BC10', 'BC11', 'BC12']

    index_dict = {'BC01': 'AAGAAAGTTGTCGGTGTCTTTGTG', 'BC02': 'TCGATTCCGTTTGTAGTCGTCTGT',
                  'BC03': 'GAGTCTTGTGTCCCAGTTACCAGG', 'BC04': 'TTCGGATTCTATCGTGTTTCCCTA',
                  'BC05': 'CTTGTCCAGGGTTTGTGTAACCTT', 'BC06': 'TTCTCGCAAAGGCAGAAAGTAGTC',
                  'BC07': 'GTGTTACCGTGGGAATGAATCCTT', 'BC08': 'TTCAGGGAACAAACCAAGTTACGT',
                  'BC09': 'AACTAGGCACAGCGAGTCTTGGTT', 'BC10': 'AAGCGTTGAAACCTTTGTCCTCTC',
                  'BC11': 'GTTTCATCTATCGGAGGGAATGGA', 'BC12': 'CAGGTAGAAAGAAGCAGAATCGGA'}

    re_list = []
    for index in index_list:
        # compile a regular expression and add it to the re_list
        c_re = regex.compile('('+index_dict[index]+'){e<=2}')
        re_list.append(c_re)

    input_q = Queue()


    for snf in site_file_name_list:
        input_q.put(snf)

    numProc = 10
    for n in range(numProc):
        input_q.put('STOP')

    allProcesses = []

    # http://stackoverflow.com/questions/8242837/django-multiprocessing-and-database-connections
    for n in range(numProc):
        p = Process(target=demul_worker, args=(input_q, index_list, seq_dir, re_list))
        allProcesses.append(p)
        p.start()



    for p in allProcesses:
        p.join()

    apples = 'asdf'





def demul_worker(input_q, index_list, seq_dir, re_list):

    for sfn in iter(input_q.get, 'STOP'):

        # check to see if we have already created the fastq files in question
        done = True
        for i, index in enumerate(index_list):
            for j, read in enumerate(['r1', 'r2']):
                if not os.path.isfile('{}/{}_{}_{}.fastq'.format(seq_dir, sfn, index, read)):
                    done = False

        if done:
            continue


        new_fastq_dict = {i : ([],[]) for i in range(len(index_list))}

        r1_sfn_file_path = '{}/{}_r1.fastq'.format(seq_dir, sfn)
        r2_sfn_file_path = '{}/{}_r2.fastq'.format(seq_dir, sfn)

        # this is going to take some time to read in a pickle
        if os.path.isfile('{}/{}_r1_dict.pickle'.format(seq_dir, sfn)):
            r1_sfn_dict = pickle.load(open('{}/{}_r1_dict.pickle'.format(seq_dir, sfn), 'rb'))
        else:
            with open(r1_sfn_file_path, 'r') as f:
                r1_sfn_file = [line.rstrip() for line in f]

            r1_sfn_dict = {}
            for i in range(len(r1_sfn_file)):
                if r1_sfn_file[i].startswith('@M'):
                    r1_sfn_dict[r1_sfn_file[i][1:-2]] = [r1_sfn_file[i + 1], r1_sfn_file[i + 2], r1_sfn_file[i + 3]]
            pickle.dump(r1_sfn_dict, open('{}/{}_r1_dict.pickle'.format(seq_dir, sfn), 'wb'))


        if os.path.isfile('{}/{}_r2_dict.pickle'.format(seq_dir, sfn)):
            r2_sfn_dict = pickle.load(open('{}/{}_r2_dict.pickle'.format(seq_dir, sfn), 'rb'))
        else:
            with open(r2_sfn_file_path, 'r') as f:
                r2_sfn_file = [line.rstrip() for line in f]

            r2_sfn_dict = {}
            for i in range(len(r2_sfn_file)):
                if r2_sfn_file[i].startswith('@M'):
                    r2_sfn_dict[r2_sfn_file[i][1:-2]] = [r2_sfn_file[i + 1], r2_sfn_file[i + 2], r2_sfn_file[i + 3]]

            pickle.dump(r2_sfn_dict, open('{}/{}_r2_dict.pickle'.format(seq_dir, sfn), 'wb'))

        # then for each name in the r1 dictionary, barcodes using a compiled regular expression
        # best to have a compiled regular expression list
        num_seqs = len(r1_sfn_dict.items())
        count = -1
        for seq_name, seq_lines in r1_sfn_dict.items():
            count += 1
            print('checking seq {}; {}/{}'.format(seq_name, count, num_seqs))
            for i, rx in enumerate(re_list):
                result = rx.search(seq_lines[0])
                if result:
                    apples = 'asdf'
                    # then we have found a match for one of the indices and we should check the corresponding seq
                    # in the other fastq
                    if seq_name in r2_sfn_dict.keys():
                        seq_lines_2 = r2_sfn_dict[seq_name]
                        result_two = rx.search(seq_lines_2[0])
                        if result_two:
                            # then both seqs have the idex in them and we should add the seq to the respective
                            # new fastq
                            # add the r1 seq
                            new_fastq_dict[i][0].append('@{}/1'.format(seq_name))
                            new_fastq_dict[i][0].extend(seq_lines)
                            new_fastq_dict[i][1].append('@{}/2'.format(seq_name))
                            new_fastq_dict[i][1].extend(seq_lines_2)

        # now write out the new fastq files
        for i, index in enumerate(index_list):
            for j, read in enumerate(['r1', 'r2']):
                with open('{}/{}_{}_{}.fastq'.format(seq_dir, sfn, index, read), 'w') as f:
                    for line in new_fastq_dict[i][j]:
                        f.write('{}\n'.format(line))
        # here we should have written out the

def demultiplex():
    seq_dir = '/home/humebc/phylogeneticSoftware/SymPortal_Data/rawData/20181120_tara_nanopore_png'

    fastq_list = [file for file in os.listdir(seq_dir) if 'fastq.gz' in file]

    site_file_name_lists = ['KimbeS1', 'KimbeS2', 'KimbeS3', 'KimbeS4', 'KimbeS4bis']

    # site_file_name_lists = list(set([file.split('_')[0] for file in fastq_list]))

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
        subprocess.run(['sabre', 'pe', '-f', r1_sfn_file, '-r', r2_sfn_file, '-b', tab_delim_path,
                        '-u', r1_reads_no_match_path, '-w', r2_reads_no_match_path, '-m', '4'])

    # now recompress all of the fastq files
    # easiest to do this in bash directly on the command line
    apples = 'pwer'

demultiplex_from_scratch()