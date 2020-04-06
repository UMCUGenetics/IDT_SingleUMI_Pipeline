# IDT SingleUMI Pipeline
Pipeline for creating FastQ files with single UMI and first n nucleotides of the read. By using the first n nucleotides of a reads, the pipeline effictively skips the mapping step.

## Installing & Setup
1. [Install Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
2. [Install Singularity](https://sylabs.io/guides/3.5/admin-guide/)
3. [Pull/Clone IDT_SingleUMI_Pipeline](#pull-or-clone)
4. [Get & configure resources](docs/resources.md)
5. [Configure nextflow](docs/nextflow.md)
6. [Configure processes](docs/processes.md)

## Pull or Clone
Be sure to add the --recursive option to also included the neccessary modules.

```
git clone git@github.com:UMCUGenetics/IDT_SingleUMI_Pipeline.git --recursive
```

## Running the workflow
In this section we'll provide you with a few different ways to run the workflow.

### Change the run-template.config to start your analysis

Always keep these lines in your run.config file:
```
includeConfig 'nextflow.config'
includeConfig 'process.config'
includeConfig 'resources.config'
```

All of the parameters in the params section can also be supplied on the commandline or can be pre-filled in the run.config file.

```
params {
  fastq_path = ''
  out_dir = ''
  singleEnd = false <- Set to true for single-end reads. Be sure to also adapt the read_structures appropriately.
  filter = false <- Set to true for QC filtering of UMI-fied reads
  sample_sheet = '/path/to/SampleSheet.csv' <- Path to samplesheet, see below for a template
  read_structures = '8B9M 8B 9M61T 70T' <- define which bases to use for reads,barcodes and UMI. See also: https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures
}
```

### Create a SampleSheet
An example SampleSheet.csv file:
```
Sample_ID,Sample_Name,Sample_Barcode
S1,100ng-5cyc-S1_HYKVFBGX9,TGAGCTAGGAACGGTT
S2,100ng-5cyc-S2_HYKVFBGX9,AGAGTAGCAAGGCGTA

```
The samplesheet is formatted as follows:
<Sample_ID><Samplename_FlowcellID><I1barcodeI2Barcode>
In the case of dual barcoding the I1 & I2 barcode sequences should be appended together.

### Starting the workflow
```
nextflow run idt_singleumi.nf -c run.config --out_dir /raw_data/runX -profile sge -resume
```
