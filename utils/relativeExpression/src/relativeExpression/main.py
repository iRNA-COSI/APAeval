##metrics to calculate
#The fraction, compared to the # of true terminal exons in our list, of reported terminal exons in those genes that are reported as having APA but are not according to the ground truth (i.e. putative FP)
#The fraction from # of true terminal exons in our list that were not reported as APA (i.e. FN)
#The Jaccard index for 1-2 above i.e. the intersection of identified and true divided(TP) by the union of identified(TP+FP) and true(FN). This Jaccard can be done at the level of terminal exons with APA reported vs true and at the PA site level for the total sites in those.
def get_FPFN_fraction(f_GT_path,matched,only_PD, only_GT):
    f_GT=pd.read_csv(f_GT_path,sep='\t')
    FP=len(only_PD) #sites found in prediction only
    FN=len(only_GT) #sites found in the filtered groud truth only
    TP=len(matched) #sites found in both the prediction and filtered ground truth
    allGT=len(f_GT) #all sites in the filtered ground truth
    frac_FP=FP/(FP+TP) ##based on cells 46-48 from Sam's relative_quant_poc.ipynb
    frac_FN=FN/(TP+FN) ##based on cells 51-53 from Sam's relative_quant_poc.ipynb
    jaccard=TP/(TP+FP+FN) ##based on cell 54 from Sam's relative_quant_poc.ipynb
    return (frac_FP,frac_FN,jaccard)
#frac_FP,frac_FN,jaccard=get_FPFN_fraction(gold_standard,matched,only_PD, only_GT)


#The fraction from the matched terminal exons with APA where both sites where correctly identified according to window size definition
#The fraction from the matched terminal exons with APA where only the proximal site was correctly located
#The fraction from the matched terminal exons with APA where only the distal site was correctly located.
#The fraction from the matched terminal exons with APA where both distal and proximal sites were not correctly located (i.e. the terminal exon was called correctly as having APA but both site locations were “off”).
def match_APA_locations(matched):  ##idk how the window size thing works
    countd={}
    both,proximal,distal=0,0,0
    matched=matched.sort_values('name_g')
    for n in matched['name_g']:
        n=n.split('|')
        te_id,type='|'.join(n[:-1]),n[-1]
        if te_id not in countd:
            countd[te_id]=[type]
        else:
            countd[te_id].append(type)
    for te_id in countd:
         if len(countd[te_id]) > 1:
             both+=1
         else:
             if countd[te_id][0] == '1.0':
                 proximal+=1
             if countd[te_id][0] == '2.0':
                 distal+=1
    totTEs=len(countd)
    frac_both,frac_proximal_only,frac_distal_only=both/totTEs,proximal/totTEs,distal/totTEs
    return (frac_both,frac_proximal_only,frac_distal_only)
#frac_both,frac_proximal_only,frac_distal_only=match_APA_locations(matched)
