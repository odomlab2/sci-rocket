# LSF profile

These commands were used to set up the LSF profile on the DKFZ LSF cluster.

## Run cookiecutter template for LSF profile

```shell
# Set to your username
profile_dir="/home/<username>/.config/snakemake"
mkdir -p "$profile_dir"

python3 -m cookiecutter $profile_dir gh:Snakemake-Profiles/lsf
```

When prompted, set the following parameters:

```text
LSF_UNIT_FOR_LIMITS=MB
UNKWN_behaviour=wait
ZOMBI_behaviour=ignore
latency-wait=30
use_conda=n
use-singularity=n
restart_times=2
print_shell_commands=y
jobs=50 (default)
default_mem_mb=1024 (default)
default_cluster_logdir=logs/cluster
default_queue=long
default_project=<leave empty>
max_status_checks_per_second=10 (default)
max_jobs_per_second=10 (default)
max_status_checks=1 (default)
wait_between_tries=0.001 (default)
jobscript_timeout=10 (default)
profile=lsf_dkfz
```
