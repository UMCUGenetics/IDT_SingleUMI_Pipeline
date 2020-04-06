# Get & configure resources (per genome)

## Steps

1. [Reference genome](#1-set-up-reference-genome)
2. [Create resources config](#2-create-resources-config)

### 1 Set up reference genome

Download your reference genome of choice from:
https://support.illumina.com/sequencing/sequencing_software/igenome.html

```
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/Ensembl/GRCh37/Homo_sapiens_Ensembl_GRCh37.tar.gz
tar -xzvf Homo_sapiens_Ensembl_GRCh37.tar.gz
```

Copy the files in the /Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/* directory to a resources directory:

```
cp -R /Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/* /IDT_SingleUMI_Pipeline/resources/GRCh37/Sequence
```

### 2 Create resources config
Adapt the configs/resources.config file to include the resources you just gathered.
Also don't forget to set the resource_dir. 
```
params {
  resource_dir = '/full/path/to/resources/dir/'
  genome_fasta = "${params.resource_dir}/GRCh37/Sequence/genome.fa"
}
```
