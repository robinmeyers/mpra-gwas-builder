# MPRA GWAS Builder: snakemake workflow

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.12.0-brightgreen.svg)](https://snakemake.bitbucket.io)

## Authors

* Robin Meyers (@robinmeyers)

## Usage

### Simple

#### Step 1: Install workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/robinmeyers/mpra-gwas-builder/releases).
If you intend to modify and further extend this workflow or want to work under version control, fork this repository as outlined in [Advanced](#advanced). The latter way is recommended.

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, if available, its DOI (see above).

#### Step 2: Configure workflow

Install and activate the conda environment

```
conda env create -n mpra-gwas-builder -f envs/conda.yaml
conda activate mpra-gwas-builder
```

Configure the workflow according to your needs via editing the file `config.yaml`.

Request an [API token for LDlink](https://ldlink.nci.nih.gov/?tab=apiaccess) and paste it into a file named `.Renviron` in this directory

```
LDLINK_TOKEN=YourTokenHere123
```

#### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N


See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.


### Advanced

The following recipe provides established best practices for running and extending this workflow in a reproducible way.

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to the desired working directory for the concrete project/run on your machine.
3. [Create a new branch](https://git-scm.com/docs/gittutorial#_managing_branches) (the project-branch) within the clone and switch to it. The branch will contain any project-specific modifications (e.g. to configuration, but also to code).
4. Modify the config, and any necessary sheets (and probably the workflow) as needed.
5. Commit any changes and push the project-branch to your fork on github.
6. Run the analysis.
7. Optional: Merge back any valuable and generalizable changes to the [upstream repo](https://github.com/robinmeyers/mpra-gwas-builder) via a [**pull request**](https://help.github.com/en/articles/creating-a-pull-request). This would be **greatly appreciated**.
8. Optional: Push results (plots/tables) to the remote branch on your fork.
9. Optional: Create a self-contained workflow archive for publication along with the paper (snakemake --archive).
10. Optional: Delete the local clone/workdir to free space.

