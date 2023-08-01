# LSF profile

These commands were used to set up the LSF profile on the DKFZ LSF cluster:

```shell
pip install --user cookiecutter

# create configuration directory that snakemake searches for profiles
profile_dir="/home/<username>/.config/snakemake"
mkdir -p "$profile_dir"
# use cookiecutter to create the profile in the config directory
template="gh:Snakemake-Profiles/lsf"
cookiecutter --output-dir "$profile_dir" "$template"

# parameters to set
LSF_UNIT_FOR_LIMITS=MB
UNKWN_behaviour=wait
ZOMBI_behaviour=ignore
use_conda=False
use-singularity=False
latency-wait=30
printshellcmds=True
restart-times=2
jobs=50
max-jobs-per-second=10
max-status-checks-per-second=10
profile=lsf_dkfz
```
