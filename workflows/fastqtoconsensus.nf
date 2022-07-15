/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowFastqtoconsensus.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, 
                           params.multiqc_config, 
						   params.fasta, 
						   params.fasta_fai, 
						   params.dict,
						   params.bwa]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Create input channel
ch_input_sample = extract_csv(file(params.input, checkIfExists: true), params.test_run)
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REFERENCES AND INDICES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
fasta       = params.fasta     ? Channel.fromPath(params.fasta).collect()     : Channel.empty()
fasta_fai   = params.fasta_fai ? Channel.fromPath(params.fasta_fai).collect() : Channel.empty()
dict        = params.dict      ? Channel.fromPath(params.dict).collect()      : Chaneel.empty()
bwa         = params.bwa       ? Channel.fromPath(params.bwa).collect()       : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FGBIO parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
read_structure       = params.read_structure                ?: Channel.empty()
gr_strategy          = params.group_reads_strategy          ?: Channel.empty()
gr_edits             = params.group_reads_edits             ?: Channel.empty()
con_min_reads        = params.consensus_min_reads           ?: Channel.empty()
con_min_base_quality = params.consensus_min_base_quality    ?: Channel.empty()
fl_min_reads         = params.filter_min_reads              ?: Channel.empty()
fl_max_error_rate    = params.filter_max_error_rate         ?: Channel.empty()
fl_min_base_quality  = params.filter_reads_min_base_quality ?: Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// SUBWORKFLOW: installed from subworkflows
//
//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                            } from '../modules/nf-core/modules/fastqc/main'
include { FGBIO_FASTQTOBAM                  } from '../modules/nf-core/modules/fgbio/fastqtobam/main'
include { PICARD_MERGESAMFILES              } from '../modules/nf-core/modules/picard/mergesamfiles/main'
include { FGBIO_ZIPPER                      } from '../modules/local/fgbio/zipper/main'
include { FGBIO_GROUPREADSBYUMI             } from '../modules/nf-core/modules/fgbio/groupreadsbyumi/main'
include { FGBIO_CALLMOLECULARCONSENSUSREADS } from '../modules/nf-core/modules/fgbio/callmolecularconsensusreads/main'
include { FGBIO_FILTERCONSENSUSREADS        } from '../modules/nf-core/modules/fgbio/filterconsensusreads/main'
include { SAMTOOLS_FASTQ                    } from '../modules/nf-core/modules/samtools/fastq/main'
include { MULTIQC                           } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow FASTQTOCONSENSUS {

    ch_versions = Channel.empty()

    //
    // MODULE: Run FastQC
    //

    FASTQC (
        ch_input_sample
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run FGBIO fastqtobam
    //

    FGBIO_FASTQTOBAM (
        ch_input_sample,read_structure
    )
    ch_versions = ch_versions.mix(FGBIO_FASTQTOBAM.out.versions.first())

    //
    // SUBWORKFLOW: merge runs if necessary
    //

    FGBIO_FASTQTOBAM.out.umibam
        .map{ meta, bam ->
            new_meta = [patient:meta.patient, sample:meta.sample, gender:meta.gender, status:meta.status, id:meta.sample, numLanes:meta.numLanes, data_type:meta.data_type,size:meta.size]
            [new_meta, bam]
        }.groupTuple()
        .set{group_umi_bam}


    PICARD_MERGESAMFILES(group_umi_bam)
    ch_versions = ch_versions.mix(PICARD_MERGESAMFILES.out.versions.first())

    //
    // MODULE: fgbio modules
    //
	
	FGBIO_ZIPPER (
        PICARD_MERGESAMFILES.out.bam,bwa,fasta,fasta_fai,dict
    )
    ch_versions = ch_versions.mix(FGBIO_ZIPPER.out.versions.first())
    FGBIO_ZIPPER.out.zipperbam.view()
    FGBIO_GROUPREADSBYUMI (
        FGBIO_ZIPPER.out.zipperbam,gr_strategy,gr_edits
    )
    ch_versions = ch_versions.mix(FGBIO_GROUPREADSBYUMI.out.versions.first())

    FGBIO_CALLMOLECULARCONSENSUSREADS (
        FGBIO_GROUPREADSBYUMI.out.groupbam,con_min_reads,con_min_base_quality
    )
    ch_versions = ch_versions.mix(FGBIO_CALLMOLECULARCONSENSUSREADS.out.versions.first())

    FGBIO_FILTERCONSENSUSREADS (
        FGBIO_CALLMOLECULARCONSENSUSREADS.out.consensusbam,fasta,fasta_fai,dict,fl_min_reads,fl_min_base_quality,fl_max_error_rate
    )
    ch_versions = ch_versions.mix(FGBIO_FILTERCONSENSUSREADS.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowFastqtoconsensus.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FGBIO_GROUPREADSBYUMI.out.histogram.collect{it[1]}.ifEmpty([]))
    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Function to extract information (meta data + file(s)) from csv file(s)
def extract_csv(csv_file, test_run) {

    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    new File(csv_file.toString()).withReader('UTF-8') { reader ->
        def line, numberOfLinesInSampleSheet = 0;
            while ((line = reader.readLine()) != null) {
            numberOfLinesInSampleSheet++
        }
        if( numberOfLinesInSampleSheet < 2){
            log.error "Sample sheet had less than two lines. The sample sheet must be a csv file with a header, so at least two lines."
            System.exit(1)
        }
    }



    Channel.from(csv_file).splitCsv(header: true)
        //Retrieves number of lanes by grouping together by patient and sample and counting how many entries there are for this combination
        .map{ row ->
            if (!(row.patient && row.sample)){
                log.error "Missing field in csv file header. The csv file must have fields named 'patient' and 'sample'."
                System.exit(1)
            }
            [[row.patient.toString(), row.sample.toString()], row]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [rows, size]
        }.transpose()
        .map{ row, numLanes -> //from here do the usual thing for csv parsing
        def meta = [:]
        def fastq_meta = []
        // Meta data to identify samplesheet
        // Both patient and sample are mandatory
        // Several sample can belong to the same patient
        // Sample should be unique for the patient
        if (row.patient) meta.patient = row.patient.toString()
        if (row.sample)  meta.sample  = row.sample.toString()

        // If no gender specified, gender is not considered
        // gender is only mandatory for somatic CNV
        if (row.gender) meta.gender = row.gender.toString()
        else meta.gender = "NA"

        // If no status specified, sample is assumed normal
        if (row.status) meta.status = row.status.toInteger()
        else meta.status = 0

        // mapping with fastq
        if (row.lane && row.fastq_2) {
            meta.id         = "${row.sample}-${row.lane}".toString()
            def fastq_1     = file(row.fastq_1, checkIfExists: true)
            fastq_meta.add(fastq_1)
            def fastq_2     = file(row.fastq_2, checkIfExists: true)
            fastq_meta.add(fastq_2)
            if(row.fastq_3) {
                def fastq_3  = file(row.fastq_3, checkIfExists: true)
                fastq_meta.add(fastq_3)
            }
            def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''

            def flowcell = flowcellLaneFromFastq(fastq_1,test_run)

            //Don't use a random element for ID, it breaks resuming
            def read_group  = "\"@RG\\tID:${flowcell}.${row.sample}.${row.lane}\\t${CN}PU:${row.lane}\\tSM:${row.patient}_${row.sample}\\tLB:${row.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""

            meta.numLanes   = numLanes.toInteger()
            meta.read_group = read_group.toString()
            meta.data_type  = "fastq"

            meta.size       = 1 // default number of splitted fastq
            return [meta, fastq_meta]
        // start from BAM
        } else if (row.lane && row.bam) {
            meta.id         = "${row.sample}-${row.lane}".toString()
            def bam         = file(row.bam,   checkIfExists: true)
            def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''
            def read_group  = "\"@RG\\tID:${row_sample}_${row.lane}\\t${CN}PU:${row.lane}\\tSM:${row.sample}\\tLB:${row.sample}\\tPL:${params.seq_platform}\""
            meta.numLanes   = numLanes.toInteger()
            meta.read_group = read_group.toString()
            meta.data_type  = "bam"
            meta.size       = 1 // default number of splitted fastq
            return [meta, bam]
        } else {
            log.warn "Missing or unknown field in csv file header"
        }
    }
}

// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path,test_run) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    if(!test_run){
        def line
        path.withInputStream {
            InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
            Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
            BufferedReader buffered = new BufferedReader(decoder)
            line = buffered.readLine()
        }
        assert line.startsWith('@')
        line = line.substring(1)
        def fields = line.split(':')
        String fcid

        if (fields.size() >= 7) {
            // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
            // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
            fcid = fields[2]
        } else if (fields.size() == 5) {
            fcid = fields[0]
        }
        return fcid
    } else {
        fcid = "test"
        return fcid
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
