# SV_detection_pipeline
Structural variation (SV) detection pipeline which uses DELLY, dysgu and Manta to detect SVs and report SVs which are supported by the three tools via 0.8 minimum reciprocal overlap.
## Dependencies
* Nextflow
* DELLY
* dysgu
* Manta
* Samtools
* BCFTools
* BEDTools
## Inputs
A reference file, indexed (using Samtools) BAM files for the samples, list of chromosomes.

List of chromosomes should be a TXT file like this:
```
Chr01
Chr02
Chr03
Chr04
Chr05
```
## Execution
