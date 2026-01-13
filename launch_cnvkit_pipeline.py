#!/usr/bin/env python3
"""
Full CNVkit pipeline launcher

Usage:
    python cnvkit_pipeline.py <PROJECT_DIR> <Sample-Info.txt>
"""

import os
import sys
import subprocess
import pandas as pd
import datetime
import random
import time

if len(sys.argv) != 3:
    print("Usage: python cnvkit_pipeline.py <PROJECT_DIR> <Sample-Info.txt>")
    sys.exit(1)

PROJECT_DIR = os.path.abspath(sys.argv[1])
SAMPLE_INFO = os.path.abspath(sys.argv[2])

# Paths
BAM_DIR = os.path.join(PROJECT_DIR, "bam")
VCF_DIR = os.path.join(PROJECT_DIR, "vcf")
TARGETS = os.path.join(PROJECT_DIR, "targets/hg38_exome_v2.0.2.bed")
ANTITARGET = os.path.join(PROJECT_DIR, "targets/hg38_exome_v2.0.2.antitarget.bed")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "output")
PLOT_DIR = os.path.join(PROJECT_DIR, "plots")
JOB_DIR = os.path.join(PROJECT_DIR, "job_files")
JOB_OUT = os.path.join(PROJECT_DIR, "job_output")
JOB_ERR = os.path.join(PROJECT_DIR, "job_errors")

FASTA = "/home/software/Clonal-Neoantigen-Large/database/hg38_minimal/bwa2/hg38.fa"

for d in [OUTPUT_DIR, PLOT_DIR, JOB_DIR, JOB_OUT, JOB_ERR]:
    os.makedirs(d, exist_ok=True)

# Pipeline log
log_number = random.randint(1000, 9999)
PIPE_LOG = os.path.join(PROJECT_DIR, f"pipeline_{log_number}.log")

def log(msg):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(PIPE_LOG, "a") as f:
        f.write(f"[{timestamp}] {msg}\n")
    print(f"[{timestamp}] {msg}")

# Read sample info
df = pd.read_csv(SAMPLE_INFO, sep="\t")
df.columns = df.columns.str.strip()
assert "Sequencing ID" in df.columns, "Sequencing ID column missing"
assert "Patient ID" in df.columns, "Patient ID column missing"

samples = df["Sequencing ID"].tolist()
patient_map = dict(zip(df["Sequencing ID"], df["Patient ID"]))

# Helper: write and submit SLURM job
def write_submit_job(job_name, script_text):
    job_script = os.path.join(JOB_DIR, f"{job_name}.sh")
    with open(job_script, "w") as f:
        f.write(script_text)
    os.chmod(job_script, 0o755)
    result = subprocess.run(["sbatch", job_script], capture_output=True, text=True)
    if result.returncode != 0:
        log(f"FAILED submitting {job_name}: {result.stderr.strip()}")
        return None
    job_id = result.stdout.strip().split()[-1]
    log(f"JOB WRITTEN & SUBMITTED: {job_name}, JobID: {job_id}")
    return job_id

# Helper: wait for list of job IDs to finish
def wait_for_jobs(job_ids):
    while True:
        if not job_ids:
            return
        running_jobs = []
        for jid in job_ids:
            squeue = subprocess.run(["squeue", "-j", jid], capture_output=True, text=True)
            if jid in squeue.stdout:
                running_jobs.append(jid)
        if not running_jobs:
            break
        time.sleep(30)

# 1. Coverage (target + antitarget)
log("STEP 1: Coverage")
coverage_job_ids = []
for sid in samples:
    bam = os.path.join(BAM_DIR, f"{sid}.consensus.mapped.filtered.recalibrated.exome.dupes-removed.bam")
    if not os.path.exists(bam):
        log(f"Missing BAM: {bam}, skipping")
        continue

    job_name = f"coverage_{sid}"
    script_text = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --output={JOB_OUT}/{job_name}_%j.out
#SBATCH --error={JOB_ERR}/{job_name}_%j.err

cd {PROJECT_DIR}

docker run --rm -v {PROJECT_DIR}:/data etal/cnvkit cnvkit.py coverage /data/bam/{sid}.consensus.mapped.filtered.recalibrated.exome.dupes-removed.bam /data/targets/hg38_exome_v2.0.2.bed --fasta /data/{os.path.relpath(FASTA, PROJECT_DIR)} -o /data/output/{sid}_target.cnn -p 2

docker run --rm -v {PROJECT_DIR}:/data etal/cnvkit cnvkit.py coverage /data/bam/{sid}.consensus.mapped.filtered.recalibrated.exome.dupes-removed.bam /data/targets/hg38_exome_v2.0.2.antitarget.bed --fasta /data/{os.path.relpath(FASTA, PROJECT_DIR)} -o /data/output/{sid}_antitarget.cnn -p 2
"""
    jid = write_submit_job(job_name, script_text)
    if jid: coverage_job_ids.append(jid)

wait_for_jobs(coverage_job_ids)

# 2. Reference using normal samples
log("STEP 2: Reference")
normals = df[df["IP_status"]=="NAT"]["Sequencing ID"].tolist()
ref_job_ids = []
if normals:
    job_name = "reference_build"
    normal_targets = " ".join([f"/data/output/{sid}_target.cnn" for sid in normals])
    script_text = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output={JOB_OUT}/{job_name}_%j.out
#SBATCH --error={JOB_ERR}/{job_name}_%j.err

cd {PROJECT_DIR}

docker run --rm -v {PROJECT_DIR}:/data -v /home/software/Clonal-Neoantigen-Large/database/hg38_minimal/bwa2:/ref etal/cnvkit cnvkit.py reference {normal_targets} --fasta /ref/hg38.fa -x female -o /data/output/reference.cnn > {OUTPUT_DIR}/reference_build.log 2>&1
"""
    jid = write_submit_job(job_name, script_text)
    if jid: ref_job_ids.append(jid)
    wait_for_jobs(ref_job_ids)
else:
    log("No normal samples found for reference build, skipping STEP 2")

# 3. Fix using both target and antitarget
log("STEP 3: Fix")
fix_job_ids = []
for sid in samples:
    job_name = f"fix_{sid}"
    script_text = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output={JOB_OUT}/{job_name}_%j.out
#SBATCH --error={JOB_ERR}/{job_name}_%j.err

cd {PROJECT_DIR}

docker run --rm -v {PROJECT_DIR}:/data etal/cnvkit cnvkit.py fix /data/output/{sid}_target.cnn /data/output/{sid}_antitarget.cnn /data/output/reference.cnn -o /data/output/{sid}.cnr
"""
    jid = write_submit_job(job_name, script_text)
    if jid: fix_job_ids.append(jid)

wait_for_jobs(fix_job_ids)

# 4. Segment
log("STEP 4: Segment")
segment_job_ids = []
for sid in samples:
    job_name = f"segment_{sid}"
    script_text = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output={JOB_OUT}/{job_name}_%j.out
#SBATCH --error={JOB_ERR}/{job_name}_%j.err

cd {PROJECT_DIR}

docker run --rm -v {PROJECT_DIR}:/data etal/cnvkit cnvkit.py segment /data/output/{sid}.cnr -o /data/output/{sid}.cns
"""
    jid = write_submit_job(job_name, script_text)
    if jid: segment_job_ids.append(jid)

wait_for_jobs(segment_job_ids)

# 5. Scatter + Diagram
log("STEP 5: Scatter & Diagram")
plot_job_ids = []
for sid in samples:
    job_name = f"plot_{sid}"
    script_text = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output={JOB_OUT}/{job_name}_%j.out
#SBATCH --error={JOB_ERR}/{job_name}_%j.err

cd {PROJECT_DIR}

docker run --rm -v {PROJECT_DIR}:/data etal/cnvkit bash -c "
cnvkit.py scatter /data/output/{sid}.cnr -s /data/output/{sid}.cns -o /data/plots/{sid}-scatter.pdf &&
cnvkit.py diagram /data/output/{sid}.cnr -s /data/output/{sid}.cns -o /data/plots/{sid}-diagram.pdf
"
"""
    jid = write_submit_job(job_name, script_text)
    if jid: plot_job_ids.append(jid)

wait_for_jobs(plot_job_ids)

# 6. Call using cns + VCF (tumor)
log("STEP 6: Call")
call_job_ids = []
for sid in samples:
    patient = patient_map[sid]
    vcf_file = os.path.join(VCF_DIR, f"{patient}.vcf")
    if not os.path.exists(vcf_file):
        log(f"Missing VCF for {sid} ({patient}), skipping call")
        continue
    job_name = f"call_{sid}"
    script_text = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output={JOB_OUT}/{job_name}_%j.out
#SBATCH --error={JOB_ERR}/{job_name}_%j.err

cd {PROJECT_DIR}

docker run --rm -v {PROJECT_DIR}:/data etal/cnvkit cnvkit.py call /data/output/{sid}.cns -y -v /data/vcf/{patient}.vcf -m clonal --purity 0.7 -o /data/output/{sid}.call.cns
"""
    jid = write_submit_job(job_name, script_text)
    if jid: call_job_ids.append(jid)

wait_for_jobs(call_job_ids)

# 7. Genemetrics
log("STEP 7: Genemetrics")
genemetrics_job_ids = []
for sid in samples:
    seg_file = os.path.join(OUTPUT_DIR, f"{sid}.cns")
    if not os.path.exists(seg_file):
        log(f"Missing segment for {sid}, skipping genemetrics")
        continue
    job_name = f"genemetrics_{sid}"
    script_text = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output={JOB_OUT}/{job_name}_%j.out
#SBATCH --error={JOB_ERR}/{job_name}_%j.err

cd {PROJECT_DIR}

docker run --rm -v {PROJECT_DIR}:/data etal/cnvkit cnvkit.py genemetrics /data/output/{sid}.cnr -s /data/output/{sid}.cns -t 0.2 -m 3 -x female --median --stdev -o /data/genemetrics/{sid}_genemetrics.tsv
"""
    jid = write_submit_job(job_name, script_text)
    if jid: genemetrics_job_ids.append(jid)

wait_for_jobs(genemetrics_job_ids)

log("Pipeline complete!")
