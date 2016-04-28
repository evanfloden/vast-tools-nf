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
params.readLen       = 101 
params.ir_version    = 2
params.expr          = false
params.groupA        = "Sample_A"
params.groupB        = "Sample_B"
params.output        = "$baseDir/results/"


log.info "V A S T - T O O L S - N F  ~  version 0.1"
log.info "====================================="
log.info "name                               : ${params.name}"
log.info "reads (FASTQ files)                : ${params.reads}"
log.info "groups (text file)                 : ${params.groups}"
log.info "group ids (group A : group B)      : ${params.groupA} : ${params.groupB}"
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
 * Validate the groups and group file
 */
groups = file(params.groups)
if ( !groups.exists() ) exit 1, "Missing groups design file: ${groups}"
groupA = params.groupA
groupB = params.groupB


/*
 * Check if to calculate expression
 */
expr = params.expr ? '--expr' : ''

/*
 *  Specify results directory
 */
results_path = file(params.output)


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
    publishDir "$results_path/align", mode: 'copy', overwrite: 'true'

    input:
    set val(name), file(reads:'*') from read_files

    output:
    file("vast_out/to_combine/${name}.exskX") into exskX
    file("vast_out/to_combine/${name}.IR2") into IR2
    file("vast_out/to_combine/${name}.IR.summary_v2.txt") into IR_summary
    file("vast_out/to_combine/${name}.micX") into micX
    file("vast_out/to_combine/${name}.eej2") into eej2
    file("vast_out/to_combine/${name}.MULTI3X") into MULTI3X
    file("vast_out/expr_out") into expression_values

    script:
    //
    // VAST-TOOLS Align
    //
    """
    vast-tools align ${reads} --stepSize 25 --output vast_out --sp ${vast_sp} --IR_version ${ir_version} --cores ${task.cpus} --readLen ${readLen} ${expr}

    #
    # It could be nice to replace this with a condtional output if possible
    #
    if [ ! -f vast_out/expr_out/${name}.cRPKM ] && [ ! -f vast_out/expr_out/${name}.3bias ];
    then 
        mkdir -p vast_out/expr_out
        echo "Please use the '--expr' flag to get expression values (cRPKM)" > vast_out/expr_out/info.txt
    fi
    """
}


process combine {
    publishDir "$results_path/combine", mode: 'copy', overwrite: 'true'

    input:
    file exskX from exskX.toSortedList()
    file IR2 from IR2.toSortedList()
    file IR_summary from IR_summary.toSortedList()
    file micX from micX.toSortedList()
    file eej2 from eej2.toSortedList()
    file MULTI3X from MULTI3X.toSortedList()

    output: 
    file 'vast_out/INCLUSION_LEVELS*' into combine_tables

    script:
    //
    // VAST-Tools Combine
    //
 
    """
    mkdir -p vast_out/to_combine
    mv -t vast_out/to_combine/. *.eej2 *.IR2 *.txt *.micX *.MULTI3X *.exskX
    vast-tools combine -o ./vast_out -sp ${vast_sp} --IR_version ${ir_version}
    """
}



process compare {
    publishDir "$results_path/compare", mode: 'copy', overwrite: 'true'

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
    a=`grep ${groupA} ${groups} | cut -f 1 | paste -d, -s` 
    b=`grep ${groupB} ${groups} | cut -f 1 | paste -d, -s`
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


