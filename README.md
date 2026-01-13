## CNVkit Pipeline for Copy Number Analysis

This pipeline performs a full CNVkit workflow from BAM and VCF files, generating CNV calls, plots, and gene-level metrics. It runs **all steps in order**, waiting for all jobs in each step to finish before moving to the next, and logs job submission, running, completion, and failure with a unique pipeline log.

**Steps:**
1. Coverage – Compute on-target (`_target.cnn`) and off-target (`_antitarget.cnn`) coverage using BAM files.  
2. Reference – Build reference `.cnn` using normal samples.  
3. Fix – Combine target and antitarget coverage with the reference to produce `.cnr`.  
4. Segment – Segment `.cnr` files to produce `.cns`.  
5. Scatter & Diagram – Generate CNV plots (`.pdf`) using `.cnr` and `.cns`.  
6. Call – Assign absolute integer copy numbers to segments using `.cns` and matched VCF files.  
7. Genemetrics – Compute gene-level statistics from `.cnr` and `.cns`.  

**Directory Structure:**
```text
PROJECT_DIR/
├── bam/                  # BAM files
├── vcf/                  # VCF files
├── targets/              # Target and antitarget BED files
├── output/               # CNVkit output (.cnn, .cnr, .cns, .call.cns)
├── plots/                # Scatter and diagram PDFs
├── genemetrics/          # Gene-level metrics
├── job_files/            # Generated SLURM job scripts
├── job_output/           # SLURM stdout
├── job_errors/           # SLURM stderr
├── pipeline_<RAND>.log   # Pipeline log with timestamps
└── Sample-Info.txt       # Sample information (Sequencing ID, Patient ID, IP_status)
```
**Usage:**
```
python cnvkit_pipeline.py <PROJECT_DIR> <Sample-Info.txt>
```
