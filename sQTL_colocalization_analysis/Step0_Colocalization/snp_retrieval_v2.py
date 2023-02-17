# script to retrieve gtex snps that are within 400kb of the as_event
# the as_event (i.e., gene) was determined to be within 400kb of ebmd gwas loci

import gzip
from collections import defaultdict
import glob
import os
import sys

python_script, slurm_array_id, gtex_sqtl_file_path, tiss_dir_path, tiss_out_path = sys.argv


# get list of lead_snps and as_events within 400kb
# as_events will be grouped and processed by chromosome
# these events are those for which to retrieve snp info from gtex sqtl file
as_events_to_process = defaultdict(lambda: set()) # chrom -> <set of as_events to process>
for fname in glob.glob(os.path.join(tiss_dir_path, '*genes.400kbtxt')):
    for as_event in open(fname):
        as_event = as_event.strip()
        chrom = as_event.split(':')[0]
        as_events_to_process[chrom].add(as_event)

# read in gtex sqtl lines by chromosome
sqtl = defaultdict(lambda: list()) # this will be renewed for every chromosome, initially will be chr1
with gzip.open(gtex_sqtl_file_path) as f:
    f.readline() # skip the header
    #i = 0 # line counter
    current_chrom = 'chr1'
    for line in f:
        line = line.decode('utf8') # bytes to string
        #i += 1
        #if i % 10000000 == 0: print(i) # track the read in of data every 10M lines       
        phenotype_id = line.split('\t')[0] # phenotype_id is as_event
        chrom = phenotype_id.split(':')[0]
        if chrom == current_chrom:   
            sqtl[phenotype_id].append(line)
        else:
            # at this point, all sqtls for a chromosome are in sqtl dictionary
            # retrieve gtex data and output *snps file
            for as_event in as_events_to_process[current_chrom]:
                gtex_lines = sqtl[as_event]
                if len(gtex_lines) > 0:
                    ofile_name = os.path.join(tiss_out_path, as_event + '.snps')
                    with open(ofile_name, 'w') as ofile:
                        ofile.write(''.join(gtex_lines))
            # renew sqtl dictionary with the next chromosome
            sqtl.clear()
            sqtl = defaultdict(lambda: list())
            current_chrom = chrom
            sqtl[phenotype_id].append(line)
    # output data for chr22 (last chromosome)
    for as_event in as_events_to_process[current_chrom]:
        gtex_lines = sqtl[as_event]
        if len(gtex_lines) > 0:
            ofile_name = os.path.join(tiss_dir_path, as_event + '.snps')
            with open(ofile_name, 'w') as ofile:
                ofile.write(''.join(gtex_lines))
