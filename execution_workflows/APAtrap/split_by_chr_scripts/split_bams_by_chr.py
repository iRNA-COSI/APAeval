samp_list_h = '''
SRR10109438
SRR10109439
SRR10109440
SRR10109441
SRR10109442
SRR10109443
GTEXsim_cerebellum_R1
GTEXsim_cerebellum_R2
GTEXsim_cerebellum_R3
GTEXsim_cerebellum_R4
GTEXsim_cerebellum_R5
GTEXsim_cerebellum_R6
GTEXsim_cerebellum_R7
GTEXsim_cerebellum_R8
GTEXsim_cerebellum_R9
GTEXsim_cerebellum_R10
GTEXsim_muscle_R1
GTEXsim_muscle_R2
GTEXsim_muscle_R3
GTEXsim_muscle_R4
GTEXsim_muscle_R5
GTEXsim_muscle_R6
GTEXsim_muscle_R7
GTEXsim_muscle_R8
GTEXsim_muscle_R9
GTEXsim_muscle_R10'''.split()

samp_list_m = '''SRR11918577
SRR11918578
SRR11918579
SRR11918580
SRR11918581
SRR11918582
SRR1811005
SRR3067958
SRR3067957
SRR3067959'''.split()

chr_list_m = [x.split()[1].split(':')[1] for x in '''@SQ     SN:chr1 LN:195471971
@SQ     SN:chr2 LN:182113224
@SQ     SN:chr3 LN:160039680
@SQ     SN:chr4 LN:156508116
@SQ     SN:chr5 LN:151834684
@SQ     SN:chr6 LN:149736546
@SQ     SN:chr7 LN:145441459
@SQ     SN:chr8 LN:129401213
@SQ     SN:chr9 LN:124595110
@SQ     SN:chr10        LN:130694993
@SQ     SN:chr11        LN:122082543
@SQ     SN:chr12        LN:120129022
@SQ     SN:chr13        LN:120421639
@SQ     SN:chr14        LN:124902244
@SQ     SN:chr15        LN:104043685
@SQ     SN:chr16        LN:98207768
@SQ     SN:chr17        LN:94987271
@SQ     SN:chr18        LN:90702639
@SQ     SN:chr19        LN:61431566
@SQ     SN:chrX LN:171031299
@SQ     SN:chrY LN:91744698
@SQ     SN:chrM LN:16299'''.split('\n')]

chr_list = [x.split()[1].split(':')[1] for x in '''@SQ     SN:chr1 LN:248956422
@SQ     SN:chr2 LN:242193529
@SQ     SN:chr3 LN:198295559
@SQ     SN:chr4 LN:190214555
@SQ     SN:chr5 LN:181538259
@SQ     SN:chr6 LN:170805979
@SQ     SN:chr7 LN:159345973
@SQ     SN:chr8 LN:145138636
@SQ     SN:chr9 LN:138394717
@SQ     SN:chr10        LN:133797422
@SQ     SN:chr11        LN:135086622
@SQ     SN:chr12        LN:133275309
@SQ     SN:chr13 s       LN:114364328
@SQ     SN:chr14        LN:107043718
@SQ     SN:chr15        LN:101991189
@SQ     SN:chr16        LN:90338345
@SQ     SN:chr17        LN:83257441
@SQ     SN:chr18        LN:80373285
@SQ     SN:chr19        LN:58617616
@SQ     SN:chr20        LN:64444167
@SQ     SN:chr21        LN:46709983
@SQ     SN:chr22        LN:50818468
@SQ     SN:chrX LN:156040895
@SQ     SN:chrY LN:57227415
@SQ     SN:chrM LN:16569'''.split('\n')]

count = 0

samp_list_h = '''GTEXsim_cerebellum_R1
GTEXsim_cerebellum_R2
GTEXsim_cerebellum_R3
GTEXsim_cerebellum_R4
GTEXsim_cerebellum_R5
GTEXsim_cerebellum_R6
GTEXsim_cerebellum_R7
GTEXsim_cerebellum_R8
GTEXsim_cerebellum_R9
GTEXsim_cerebellum_R10'''.split()

fw = open('split_gtex_cer.sh','w')

for s in samp_list_h:
    for c in chr_list:
        orig_bam = '/data/apaeval/nf_rnaseq/bams_sorted/%s.bam'%s
        chr_bam = '/data/apaeval/nf_rnaseq/bams_sorted/bams_by_chr/%s_%s.bam'%(s, c)
        chr_sorted = '/data/apaeval/nf_rnaseq/bams_sorted/bams_by_chr/%s_%s.sorted.bam'%(s, c)
        fw.write('samtools view -b -o %s %s %s\n'%(chr_bam, orig_bam, c))
        fw.write('samtools sort %s --no-PG -O bam -o %s -@ 12\n'%(chr_bam, chr_sorted))
        fw.write('rm %s\nmv %s %s\n'%(chr_bam, chr_sorted, chr_bam))
        fw.write('samtools index %s -@ 12\n'%(chr_bam))
    count += 1
    #if count % 5 == 0:
        #fw.close()
        #fw = open('split_human_%s.sh'%count, 'w')
fw.close()

fw = open('split_mouse.sh','w')
for s in samp_list_m:
    for c in chr_list_m:
        orig_bam = '/data/apaeval/nf_rnaseq/bams_sorted/%s.bam'%s
        chr_bam = '/data/apaeval/nf_rnaseq/bams_sorted/bams_by_chr/%s_%s.bam'%(s, c)
        chr_sorted = '/data/apaeval/nf_rnaseq/bams_sorted/bams_by_chr/%s_%s.sorted.bam'%(s, c)
        fw.write('samtools view -b -o %s %s %s\n'%(chr_bam, orig_bam, c))
        fw.write('samtools sort %s --no-PG -O bam -o %s -@ 12\n'%(chr_bam, chr_sorted))
        fw.write('rm %s\nmv %s %s\n'%(chr_bam, chr_sorted, chr_bam))
        fw.write('samtools index %s -@ 12\n'%(chr_bam))
fw.close()