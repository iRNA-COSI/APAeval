/*
 * IsoSCM module containing all three steps in IsoSCM program:
 *    - Assemble
 *    - Enumerate (optional)
 *    - Compare
 */


// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)


process ISOSCM_ALL {
        publishDir "${params.outdir}", mode: params.publish_dir_mode
        
        input:
        path samplesheet
        path bamdir

        output:
        path "*"
        
        
        shell:
        '''
        sampleA=$(awk -F',' 'NR==2{print $1}' !{samplesheet})
        sampleB=$(awk -F',' 'NR==3{print $1}' !{samplesheet})
        bamA=!{bamdir}/$sampleA*.bam
        bamB=!{bamdir}/$sampleB*.bam
        strandA=$(awk -F',' 'NR==2{print $3}' !{samplesheet})
        strandB=$(awk -F',' 'NR==3{print $3}' !{samplesheet})
        
        # Assemble
        ## must specify the output dir as 'isoscm' because there's a bug in the enumerate step..
        java -Xmx2048m -jar /IsoSCM-2.0.12.jar assemble -bam $bamA -base $sampleA -s $strandA -dir isoscm
        java -Xmx2048m -jar /IsoSCM-2.0.12.jar assemble -bam $bamB -base $sampleB -s $strandB -dir isoscm
        
        for s in $sampleA $sampleB; do
            ## generate bed files here:
            get_bed_files.sh $s ./isoscm
        
            # Enumerate
            java -Xmx2048m -jar /IsoSCM-2.0.12.jar enumerate -max_isoforms 20 -x isoscm/${s}.assembly_parameters.xml
        done
        
        # Compare
        java -Xmx2048m -jar /IsoSCM-2.0.12.jar compare -base ${sampleA}.vs.${sampleB} -x1 isoscm/${sampleA}.assembly_parameters.xml -x2 isoscm/${sampleB}.assembly_parameters.xml
        
        cp compare/*.txt ./IsoSCM_03.txt
        rm tmp.gtf isoscm*.log
        '''
}

