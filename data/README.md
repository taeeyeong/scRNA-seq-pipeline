## FASTQ File Naming Convention for Cell Ranger

Proper naming of FASTQ files is essential for seamless data processing with Cell Ranger. 
This guide outlines the standard naming conventions to ensure compatibility and ease of use when using tools like 'bcl-convert', 'bcl2fastq', and 'mkfastq'.

### Naming Formats
- [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz or 
- [Sample Name]_S1_[Read Type]_001.fastq.gz

Where Read Type is one of:

-  I1: Sample index read1 (optional)
- I2: Sample index read2 (optional)
- R1: Read 1 (forward reads in paired-end sequencing)
- R2: Read 2 (reverse reads in paired-end sequencing)

e.g. 
- SampleA_S1_L001_R1_001.fastq.gz
- SampleB_S1_R1_001.fastq.gz


## Sample Description
The data included 12 patient samples representing different histologies, specifically CRPC adenocarcinoma and transformed neuroendocrine prostate cancer (NEPC).

|Sample ID| BioSample    | Run (SRR)                                 | Group    | Pathology (Clinical) | Site of Biopsy | PSA  | Prior to Treatment | # of prior regimen(s)
|---------|--------------|-------------------------------------------|----------|-----------|-----------------|------|--------------------|------------|
|HMP04    | SAMN30111353 | SRR20761377, SRR20761378                  |CRPC-NEPC |Metastatic poorly differentiated neuroendocrine carcinoma, small cell type|Liver metastasis|0.10|Bicalutamide, Lupron, Enzalutamide, Abiraterone/Cabazitaxel, Carboplatin/Etoposide|5|
|HMP05    | SAMN30111352 | SRR20761375, SRR20761376                  |CRPC-adeno|Poorly Differentiated Adenocarcinoma|Liver metastasis|0.15|Bicalutamide, Degarelix, Abiraterone|3|
|HMP08    | SAMN30111351 | SRR20761373, SRR20761374                  |CRPC-adeno|Poorly Differentiated Adenocarcinoma (weakly positive for Nkx3.1)|Mediastinal mass|<0.05|Degarelix, Enzalutamide, Darolutamide|3|
|HMP11_1  | SAMN30111350 | SRR20761355, SRR20761356                  |CRPC-adeno|Metastatic Prostate Adenocarcinoma|Brain metastasis|729|Lupron, Abiraterone, Enzalutamide, Radium-223, Docetaxel, Cabazitaxel, Carboplatin|7|
|HMP11_2  | SAMN30111349 | SRR20761353, SRR20761354                  |CRPC-adeno|Metastatic Prostate Adenocarcinoma|Epidural metastasis|2306|Lupron, Abiraterone, Enzalutamide, Radium-223, Docetaxel, Cabazitaxel, Carboplatin|7|
|HMP13    | SAMN30111348 | SRR20761351, SRR20761352                  |CRPC-adeno|Metastatic Prostate Adenocarcinoma|Supraclavicular lymph node|0.71|Lupron, Bicalutamide, Enzalutamide, Abiraterone/Cabazitaxel, Ipilumimab/Nivolumab|5|
|HMP14    | SAMN30111347 | SRR20761349, SRR20761350                  |CRPC-adeno|Metastatic Prostate Adenocarcinoma|Left posteriormedial external lymph node|0.98|Lupron, Docetaxel|2|
|HMP15    | SAMN30111346 |                                           |          |           |                 |      |                    |
|HMP16    | SAMN30111345 | SRR20761346, SRR20761347                  |CRPC-NEPC |Metastatic poorly differentiated neuroendocrine carcinoma, small cell type|Liver metastasis|34.8|Bicalutamide, Lupron, Abiraterone, Radium-223, Enzalutamide, Docetaxel/Carboplatin|6|
|HMP17    | SAMN30111344 | SRR20761344, SRR20761345                  |CRPC-NEPC |Metastatic poorly differentiated neuroendocrine carcinoma, small cell type|Left retrocrural abdomen mass|399.24|Degarelix, Bicalutamide, Lu-177 PSMA, Enzalutamide|4|
|HMP19    | SAMN30111343 | SRR20761342, SRR20761343                  |CRPC-adeno|Metastatic Prostate Adenocarcinoma|Para-aortic lymph node|1419.66|Nilandron, Ketoconazole, Lupron, Abiraterone/Cabazitaxel, Enzalutamide|5|
|HMP20    | SAMN30111342 | SRR20761338, SRR20761339, SRR20761340, SRR20761341 |CRPC-adeno|Malignant cells, consistent with prostate origin|Pelvic LN|1.35|Lupron, Abiraterone|2|
|HMP25    | SAMN30111341 | SRR20761332, SRR20761332, SRR20761334     |CRPC-adeno|Metastatic Prostate Adenocarcinoma|Liver|<0.05|Degarelix, Casodex, Carboplatin, Docetaxel|4|
|HMP26    | SAMN30111340 | SRR20761328, SRR20761329                  |CRPC-adeno|Poorly Differentiated Adenocarcinoma (weakly positive for Nkx3.1)|Liver|0.10|Lupron, Exanlutamide, Docetaxel, Carboplatin|4|
