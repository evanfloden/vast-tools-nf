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
params.species       = "human"
params.groups        = "$baseDir/tutorial/groups/group_ids.txt"
params.ir_version    = 2
params.output        = "results/"


log.info "V A S T - T O O L S - N F  ~  version 0.1"
log.info "====================================="
log.info "name                       : ${params.name}"
log.info "reads                      : ${params.reads}"
log.info "species                    : ${params.species}"
log.info "groups                     : ${params.groups}"
log.info "intron retention version   : ${params.ir_version}"
log.info "output                     : ${params.output}"
log.info "\n"


/*
 *  Validate the species input
 */
species = params.species

if ( species == 'human' ) { vast_sp = 'Hsa' }

ir_version = params.ir_version as int

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
    set val(name), file('vast_out') into vast_align_out_dirs 

    script:
    //
    // VAST-TOOLS Align
    //
    def single = reads instanceof Path
    if( !single ) {
        """
        vast-tools align --output ${name} ${reads} --IR_version ${ir_version}
        """
    }  
    else {
        """
        vast-tools align --output ${name} ${reads} --IR_version ${ir_version}
        """
    }

}


process combine {
    input:
    file 'to_combine/*' from vast_align_out_dirs.toSortedList()   

    output: 
    file 'combine/INCLUSION_TABLE.tab' into combine_tables

    script:
    //
    // VAST-Tools Combine
    //
 
    """
    vast-tools combine -o combine -sp ${vast_sp} --IR_version ir_version
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


