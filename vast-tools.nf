/*
 * Copyright (c) 2016, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'VAST-TOOLS-NF'.
 *
 *   VAST-TOOLS-NF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   VAST-TOOLS-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with VAST-TOOLS-NF.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * Main VAST-TOOLS-NF pipeline script
 *
 * @authors
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Evan Floden <evanfloden@gmail.com> 
 * Veronica Venturi <veronica.venturi@crg.eu>
 * Chris Wyatt <chris.wyatt@crg.eu>
 *
 * VAST-TOOLS is developed by several contributors: 
 * https://github.com/vastgroup/vast-tools#contributions
 *
 */

params.name          = "A Nextflow Implementation of VAST-TOOLS"
params.reads         = "$baseDir/tutorial/reads/*.fastq"
params.groups        = "$baseDir/tutorial/groups/groups.txt"
params.species       = "human"
params.readLen       = 100 
params.ir_version    = 2
params.expr          = false
params.output        = "results/"


log.info "V A S T - T O O L S - N F  ~  version 0.1"
log.info "====================================="
log.info "name                               : ${params.name}"
log.info "reads (FASTQ files)                : ${params.reads}"
log.info "groups (text file)                 : ${params.groups}"
log.info "species (human/mouse/chicken)      : ${params.species}"
log.info "read length (integer)              : ${params.readLen}"
log.info "intron retention version (1 or 2)  : ${params.ir_version}"
log.info "calculate expression (cRPKMs)      : ${params.expr}"
log.info "output (directory)                 : ${params.output}"
log.info "\n"


/*
 *  Validate the species input
 */
species = params.species
if ( species == 'human' ) { vast_sp = 'Hsa' }


/*
 *  Determine Intron Retention Version to be Used (1 or 2)
 */
ir_version = params.ir_version as int


/*
 *  Specify the Read Length of the FASTQ files
 */
readLen = params.readLen as int


/*
 * Validate the groups file
 /
groups = file(params.groups)
if ( !groups.exists() ) exit 1, "Missing groups design file: ${groups}"


/*
 * Check if to calculate expression
 */
expr = params.expr ? '--expr' : ''


/*
 * Create a channel for read files 
 */
 
Channel
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { path -> 
       def prefix = readPrefix(path, params.reads)
       tuple(prefix, path) 
    }
    .groupTuple(sort: true)
    .set { read_files } 


process vast_tools_align {
    tag "align: $name"

    input:
    set val(name), file(reads:'*') from read_files

    output:
    file("vast_out/to_combine/${name}.exskX") into exskX
    file("vast_out/to_combine/${name}.IR2") into IR2
    file("vast_out/to_combine/${name}.IR.summary_v2.txt") into IR_summary
    file("vast_out/to_combine/${name}.micX") into micX
    file("vast_out/to_combine/${name}.eej2") into eej2
    file("vast_out/to_combine/${name}.MULTI3X") into MULTI3X

    script:
    //
    // VAST-TOOLS Align
    //
    def single = reads instanceof Path
    if( !single ) {
        """
        vast-tools align --output vast_out ${reads} --IR_version ${ir_version} ${expr}
        """
    }  
    else {
        """
        vast-tools align --output vast_out ${reads} --IR_version ${ir_version} ${expr}
        """
    }

}


process combine {

    input:
    file exskX from exskX.toSortedList()
    file IR2 from IR2.toSortedList()
    file IR_summary from IR_summary.toSortedList()
    file micX from micX.toSortedList()
    file eej2 from eej2.toSortedList()
    file MULTI3X from MULTI3X.toSortedList()

    output: 
    file 'outdir/INCLUSION_LEVELS*' into combine_tables

    script:
    //
    // VAST-Tools Combine
    //
 
    """
    mkdir -p outdir/to_combine
    mv -t outdir/to_combine/. *.eej2 *.IR2 *.txt *.micX *.MULTI3X *.exskX
    vast-tools combine -o outdir -sp ${vast_sp} --IR_version ${ir_version}
    """
}



process compare {

    input:
    file inclusion_table from combine_tables
    file groups

    output:
    file 'diff_spliced.tab' into diff_spliced

    script:
    //
    // VAST-Tools Compare
    //

    """
    a=`grep 'Sample_A' ${groups} | cut -f 1 | paste -d, -s` 
    b=`grep 'Sample_B' ${groups} | cut -f 1 | paste -d, -s`
    vast-tools compare ${inclusion_table} -a \$a -b \$b --min_dPSI 25 --min_range 5 --outFile diff_spliced.tab --no_plot
    """

}


// ===================== UTILITY FUNCTIONS ============================


/* 
 * Helper function, given a file Path 
 * returns the file name region matching a specified glob pattern
 * starting from the beginning of the name up to last matching group.
 * 
 * For example: 
 *   readPrefix('/some/data/file_alpha_1.fa', 'file*_1.fa' )
 * 
 * Returns: 
 *   'file_alpha'
 */
 
def readPrefix( Path actual, template ) {

    final fileName = actual.getFileName().toString()

    def filePattern = template.toString()
    int p = filePattern.lastIndexOf('/')
    if( p != -1 ) filePattern = filePattern.substring(p+1)
    if( !filePattern.contains('*') && !filePattern.contains('?') ) 
        filePattern = '*' + filePattern 
  
    def regex = filePattern
                    .replace('.','\\.')
                    .replace('*','(.*)')
                    .replace('?','(.?)')
                    .replace('{','(?:')
                    .replace('}',')')
                    .replace(',','|')

    def matcher = (fileName =~ /$regex/)
    if( matcher.matches() ) {  
        def end = matcher.end(matcher.groupCount() )      
        def prefix = fileName.substring(0,end)
        while(prefix.endsWith('-') || prefix.endsWith('_') || prefix.endsWith('.') ) 
          prefix=prefix[0..-2]
          
        return prefix
    }
    
    return fileName
}


