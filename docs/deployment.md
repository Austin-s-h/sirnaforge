# Deployment Guide

Production deployment options for siRNAforge.

## Docker (Recommended)

Complete environment with all dependencies included.

### Dependencies Included
- Python packages (siRNA design, scoring)
- Nextflow ≥25.04.0 (workflow management)
- BWA-MEM2 (sequence alignment)
- SAMtools (BAM/SAM processing)
- ViennaRNA (RNA folding)
- AWS CLI (S3 downloads)
- Java 17+ (Nextflow requirement)

## Deployment Options

### 🐳 Docker (Recommended)

The comprehensive Docker image includes everything you need:

```bash
# Build the comprehensive image
docker build -f docker/Dockerfile -t sirnaforge:latest .

# Quick start - complete workflow
docker run -v $(pwd):/workspace -w /workspace sirnaforge:latest \
    sirnaforge workflow TP53 --output-dir results

# Interactive usage
docker run -it -v $(pwd):/workspace -w /workspace sirnaforge:latest bash

# Run Nextflow pipeline directly
docker run -v $(pwd):/workspace -w /workspace sirnaforge:latest \
    nextflow run nextflow_pipeline/main.nf \
    --input candidates.fasta \
    --outdir results \
    --download_indexes true
```

**Advantages:**
- ✅ **Everything included** - No additional installation needed
- ✅ **Consistent environment** - Same results everywhere
- ✅ **Easy CI/CD integration** - Perfect for automated pipelines
- ✅ **Cloud-ready** - Works with AWS Batch, Kubernetes, etc.

### 🖥️ Local Installation (Alternative)

For local development without Docker, see [Getting Started](getting_started.md) for basic installation, then add bioinformatics tools:

```bash
# 1. Basic installation (see Getting Started guide)
# 2. Add bioinformatics tools via micromamba
micromamba env create -f environment.yml
micromamba activate nextflow
```

## Production Deployment

### ☁️ AWS Batch

The Docker image works seamlessly with AWS Batch:

```bash
# Create Nextflow config for AWS
cat > nextflow.config << 'EOF'
profiles {
    aws {
        process.executor = 'awsbatch'
        process.queue = 'sirnaforge-queue'
        process.container = 'your-account.dkr.ecr.region.amazonaws.com/sirnaforge:latest'
        aws.region = 'us-east-1'
        aws.batch.cliPath = '/usr/local/bin/aws'
    }
}
EOF

# Run on AWS Batch
docker run -v $(pwd):/workspace -w /workspace sirnaforge:latest \
    nextflow run nextflow_pipeline/main.nf \
    -profile aws \
    --input s3://your-bucket/candidates.fasta \
    --outdir s3://your-bucket/results
```

### 🖥️ HPC/SLURM

For HPC environments, use Apptainer (formerly Singularity):

```bash
# Convert Docker image to Apptainer
apptainer build sirnaforge.sif docker://your-registry/sirnaforge:latest

# Submit SLURM job
cat > sirnaforge_job.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=sirnaforge
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=4:00:00

apptainer run sirnaforge.sif \
    sirnaforge workflow $1 --output-dir results
EOF

sbatch sirnaforge_job.sh TP53
```

### 🚢 Kubernetes

Deploy to Kubernetes for scalable analysis:

```yaml
# sirnaforge-job.yaml
apiVersion: batch/v1
kind: Job
metadata:
  name: sirnaforge-analysis
spec:
  template:
    spec:
      containers:
      - name: sirnaforge
        image: your-registry/sirnaforge:latest
        command: ["sirnaforge"]
        args: ["workflow", "TP53", "--output-dir", "/results"]
        volumeMounts:
        - name: results-volume
          mountPath: /results
        resources:
          requests:
            memory: "16Gi"
            cpu: "4"
      restartPolicy: Never
      volumes:
      - name: results-volume
        persistentVolumeClaim:
          claimName: sirnaforge-pvc
```

## Genome Index Setup

### Option 1: S3 Auto-download (Recommended for Testing)

```bash
# Uses public iGenomes repository
nextflow run nextflow_pipeline/main.nf \
    --input candidates.fasta \
    --genome_species human \
    --download_indexes true
```

### Option 2: Local Index Building

```bash
# Download reference genome
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Build BWA index
micromamba activate nextflow
bwa-mem2 index -p /data/genomes/human/GRCh38 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Configure in genomes.yaml
human:
  index_prefix: /data/genomes/human/GRCh38
```

## Environment Variables

Key environment variables for configuration:

```bash
# Nextflow configuration
export NXF_WORK=/tmp/nextflow-work
export NXF_CONDA_CACHEDIR=$HOME/.nextflow-conda-cache

# AWS configuration (if using S3 features)
export AWS_DEFAULT_REGION=us-east-1
export DOWNLOAD_INDEXES=true

# Container configuration
export NXF_DOCKER_SUDO=false
```

## Troubleshooting

### Common Issues

1. **Java not found**: Nextflow requires Java 11+
   ```bash
   sudo apt-get install openjdk-11-jre-headless
   ```

2. **BWA-MEM2 not found**: Check PATH or use full path
   ```bash
   which bwa-mem2
   micromamba activate nextflow  # If using micromamba
   ```

3. **Docker permission denied**:
   ```bash
   sudo usermod -aG docker $USER
   # Logout and login again
   ```

4. **S3 download fails**: Check AWS CLI configuration
   ```bash
   aws configure  # Or use --no-sign-request for public buckets
   ```

### Performance Tuning

- **Memory allocation**: Adjust `process.memory` in `nextflow.config`
- **CPU cores**: Modify `process.cpus` based on available resources
- **Temporary directory**: Set `NXF_WORK` to fast storage (SSD)
- **Index caching**: Use local indices for repeated analyses

## Integration Examples

### Jupyter Notebook Integration

```python
import subprocess
import os

# Ensure Nextflow environment is available
subprocess.run(['micromamba', 'activate', 'nextflow'], shell=True)

# Run siRNAforge workflow
result = subprocess.run([
    'uv', 'run', 'sirnaforge', 'workflow', 'TP53',
    '--output-dir', 'analysis_results'
], capture_output=True, text=True)

print(result.stdout)
```

### CI/CD Integration

```yaml
# .github/workflows/test.yml
name: Test Pipeline
on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v3

    - name: Setup Micromamba
      uses: mamba-org/provision-with-micromamba@main
      with:
        environment-file: environment.yml

    - name: Install Python dependencies
      run: uv sync

    - name: Test pipeline
      run: |
        micromamba activate nextflow
        uv run sirnaforge workflow TP53 --output-dir test_results
```

This comprehensive approach ensures reproducible deployments across different environments while maintaining flexibility for various use cases.
