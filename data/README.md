## FASTQ File Naming Convention for Cell Ranger

Proper naming of FASTQ files is essential for seamless data processing with Cell Ranger. 
This guide outlines the standard naming conventions to ensure compatibility and ease of use when using tools like 'bcl-convert', 'bcl2fastq', and 'mkfastq'.

### Naming Formats
[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz or 
[Sample Name]_S1_[Read Type]_001.fastq.gz

Where Read Type is one of:

I1: Sample index read1 (optional)
I2: Sample index read2 (optional)
R1: Read 1 (forward reads in paired-end sequencing)
R2: Read 2 (reverse reads in paired-end sequencing)

e.g. SampleA_S1_L001_R1_001.fastq.gz
     SampleB_S1_R1_001.fastq.gz


## Sample Description
The data included 12 patient samples representing different histologies, specifically CRPC adenocarcinoma and transformed neuroendocrine prostate cancer (NEPC).

| BioSample     | Run(SRR)                                            | 
|---------------|-----------------------------------------------------|
| SAMN30111353  | SRR20761377, SRR20761378                            | 
| SAMN30111352  | SRR20761375, SRR20761376                            | 
| SAMN30111351  | SRR20761373, SRR20761374                            |
| SAMN30111350  | SRR20761355, SRR20761356                            |
| SAMN30111349  | SRR20761353, SRR20761354                            |
| SAMN30111348  | SRR20761351, SRR20761352                            |
| SAMN30111347  | SRR20761349, SRR20761350                            |
| SAMN30111346  |                                                     |
| SAMN30111345  | SRR20761346, SRR20761347                            |
| SAMN30111344  | SRR20761344, SRR20761345                            |
| SAMN30111343  | SRR20761342, SRR20761343                            |
| SAMN30111342  | SRR20761338, SRR20761339, SRR20761340, SRR20761341  |
| SAMN30111341  | SRR20761332, SRR20761332, SRR20761334               |
| SAMN30111340  | SRR20761328, SRR20761329                            |