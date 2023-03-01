%
transcripts = {}


# read in the GTF file
with open('./CDS_corrected_3.gtf', 'r') as infile:
    for line in infile:
        # skip comment lines
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        feature_type = fields[2]
        # only process transcript, exon, and CDS features
        if feature_type not in ['transcript', 'exon', 'CDS']:
            continue
        attributes = fields[8].split(';')
        # get the transcript_id attribute
        transcript_id = None
        for attr in attributes:
            if attr.startswith('transcript_id'):
                transcript_id = attr.split('"')[1]
                break
        # skip lines without a transcript_id attribute
        if transcript_id is None:
            continue
        # add the line to the dictionary of transcripts
        if transcript_id not in transcripts:
            transcripts[transcript_id] = [line]
        else:
            transcripts[transcript_id].append(line)

#%%
# write out the unique transcripts to a new GTF file
with open('output.gtf', 'w') as outfile:
    for transcript_id, lines in transcripts.items():
        for line in lines:
            outfile.write(line)


