import gffutils

# Create gff data base with gffutils and write only features annotated as genes to a new bed file.
# This only needs to run once.
# Takes a long time to run and there are probably faster ways, but the generated files are small and can be stored easily.



# mouse
db = gffutils.create_db('gencode.vM27.annotation.gff3', dbfn='mmus.db', force=False, verbose=True, merge_strategy='merge')

all_genes = [i.id for i in db.features_of_type('gene')]

with open('only_genes_mmus.bed', 'w') as f:
    f.write('chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n')
    for gene in all_genes:
        f.write(db[gene].seqid.strip('chr') + '\t')
        f.write(str(db[gene].start) + '\t')
        f.write(str(db[gene].end) + '\t')
        f.write(gene + '\t')
        f.write('0\t')
        f.write(db[gene].strand + '\n')


# human
db = gffutils.create_db('gencode.v38.annotation.gff3', dbfn='hsap.db', force=False, verbose=True, merge_strategy='merge')

all_genes = [i.id for i in db.features_of_type('gene')]

with open('only_genes_hsap.bed', 'w') as f:
    f.write('chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n')
    for gene in all_genes:
        f.write(db[gene].seqid.strip('chr') + '\t')
        f.write(str(db[gene].start) + '\t')
        f.write(str(db[gene].end) + '\t')
        f.write(gene + '\t')
        f.write('0\t')
        f.write(db[gene].strand + '\n')