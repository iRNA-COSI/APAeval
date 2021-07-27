/*
 * RUN ISOSCM
 * Better run the three processes in one working dir (one module)
 * as the program will search for the absolute output file paths
 */

include { ISOSCM_ALL } from '../modules/isoscm_processes' addParams( options: [:] )


workflow RUN_ISOSCM {
    take:
    ch_samplesheet
    ch_bamdir
    
    main:
    
    /*
     * Assemble, Enumerate, Compare
     */
     
     ISOSCM_ALL ( ch_samplesheet, ch_bamdir )
}

