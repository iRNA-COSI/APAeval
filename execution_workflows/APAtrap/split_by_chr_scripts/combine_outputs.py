import pandas as pd

# initialize the list of chromosomes to be combined
chromosomes = [str(x) for x  in range(1,23)] + ["X","Y","M"]
df = pd.read_csv("samplesheet_chr1.csv")
# get the name of samples to be processed from the sample sheet
samples = [x.split('_chr')[0] for x in df['sample']]
# assert number of samples to be combined
assert len(samples) == 4

# for each sample, combine all the chromosomes into one output
for sample in samples:
    # for the current sample, create the final output for identification and quantification challenges
    sample_identification = open(sample + "_APAtrap_01.bed", "a+")
    sample_quantification = open(sample + "_APAtrap_02.bed", "a+")
    # loop through the list of chromosomes to be added to the final output files above
    for chromosome in chromosomes:
        data_dir = "/home/fzhuang/apatrap/results/apatrap/challenges_outputs"
        # current chromosome files to be added to the final files
        curr_sample_identification = open(data_dir + "/" + sample + "_chr" + chromosome + "_APAtrap_01.bed")
        curr_sample_quantification = open(data_dir + "/" + sample + "_chr" + chromosome + "_APAtrap_02.bed")
        # write to the final output files for the current sample
        sample_identification.write(curr_sample_identification.read())
        sample_quantification.write(curr_sample_quantification.read())


    sample_identification.close()
    sample_quantification.close()
