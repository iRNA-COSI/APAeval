
import sys
import os
import pandas as pd
from apybiomart import query

# Define input/output files
infile = sys.argv[1]
out01  = sys.argv[2]
out03  = sys.argv[3]

# Read in the APA-Scan .tsv file
apa = pd.read_csv(infile, sep="\t")

# Format the APA-Scan df so that the gene names are uppercase
apa["Gene Name"] = apa["Gene Name"].str.upper() # Need to do this as APA-Scan outputs genenames lowercase. Not actually sure what naming convention that is but uppercasing means they should work for the apybiomart queries

# Run apybiomart query
bmart = query(attributes=["ensembl_gene_id", "hgnc_symbol"],
      filters={},dataset="hsapiens_gene_ensembl")

# Merge dataframes to get common
outdf = pd.merge(bmart,apa,left_on="HGNC symbol",right_on="Gene Name")

# Rearrange to BED format (Format 01)/Format 03 and remove any duplicated rows (PolyA sites)
outbed = outdf[["Chrom","Start","End","Gene stable ID","p-value","Strand"]].drop_duplicates()
outbed["p-value"] = "." # Set p-value filed to "." so that it matches the BED format (Format 01) specification from the Execution Workflows README: https://github.com/iRNA-COSI/APAeval/tree/main/execution_workflows
out03 = outdf[["Gene stable ID","p-value"]].drop_duplicates()

# Remove existing files (else pandas will append when saving)
if os.path.exists(out01):
    os.remove(out01)
if os.path.exists(out03):
    os.remove(out03)

# And save to file
outbed.to_csv(out01,sep="\t",index=False,header=False)
out03.to_csv(out03,sep="\t",index=False,header=False)
