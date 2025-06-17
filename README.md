# CWL workflow fo IMPACT trial project
This is a cwl wrapper for streamling all the processes in [impact-trial](https://github.com/childrens-bti/impact-trial)



## To be modified ....

## To reproduce the code in this repository:
This repository contains a docker image and code used to conduct analyses for the project.

1. Clone the repository
```
git clone git@github.com:childrens-bti/impact-trial.git
```

2. Set up AWS SSO to download the data (Make sure you have AWS CLI >= V2 installed.)
```bash
aws configure sso
```

| Prompt           | What to enter                                                        |
| ---------------- | -------------------------------------------------------------------- |
| SSO session name | `impact-session`                                                     |
| SSO start URL    | `https://d-9067576cea.awsapps.com/start/`                            |
| SSO region       | `us-east-1`                                                          |
| SSO registration scopes | [sso:account:access]: (press Enter for default)               |
| AWS account ID   |  <Get from AWS slack channel>                                        |
| Role name        | *(select from list – e.g., `ReadOnly`, `ReadWrite`, or `PowerUser`)* |
| CLI default client Region [None]: |  (press Enter for default)                          |
| CLI default output format [None] |  (press Enter for default)                           |
| CLI profile name | `impact-sso`                                                         |


3. Run the `download_data.sh` shell script to obtain latest data files: 
```
bash download_data.sh
```

4. Pull the docker container:
```
docker pull pgc-images.sbgenomics.com/rokita-lab/impact-trial:latest
```

5. Start the docker container

`cd` into the `impact-trial` folder and run the following code:

```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 8787:8787 -v $PWD:/home/rstudio/impact-trial pgc-images.sbgenomics.com/rokita-lab/impact-trial:latest
```

NOTE: replace `<CONTAINER_NAME>` with an informative name for the docker image (e.g., `impact-trial`)
 
Users can also run Rstudio in the project docker container from a web browser using the instructions below:

__Local Development in Rstudio__ (Max OS X and Linux users only)

Execute the above `docker run` command, then navigate to `localhost:8787` in your web browser. The username for login is `rstudio` and the password will be whatever password is set in the `docker run` command above (default: `pass`).

__Development using Amazon EC2__

Execute the following `docker run` command:

```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=pass -p 80:8787 -v $PWD:/home/rstudio/impact-trial pgc-images.sbgenomics.com/rokita-lab/impact-trial:latest
```

Then, paste the instance IP address into your browser to start Rstudio. 

6. To execute shell within the docker image (locally or on EC2): from the `impact-trial` folder, run:
```
docker exec -ti <CONTAINER_NAME> bash
```

Making sure to update `<CONTAINER_NAME>` with the name you have given the docker image

### Below is the main directory structure listing the analyses and data files used in this repository

```
.
├── Dockerfile
├── LICENSE
├── README.md
├── analyses
├── docs
├── data
│   ├── fusion-arriba-test.tsv.gz -> v1/fusion-arriba-test.tsv.gz
│   ├── fusion-starfusion-test.tsv.gz -> v1/fusion-starfusion-test.tsv.gz
│   ├── histologies-test.tsv -> v1/histologies-test.tsv
│   ├── snv-consensus-plus-hotspots.maf-test.tsv -> v1/snv-consensus-plus-hotspots.maf-test.tsv
│   ├── gencode.v39.primary_assembly.annotation.gtf.gz
│   ├── GCF_000001405.40_GRCh38.p14_protein.faa.gz
│   ├── Homo_sapiens.GRCh38.pep.all.fa.gz
│   ├── release-notes.md -> v1/release-notes.md
│   └── v1
├── download_data.sh
├── figures
└── scripts
```


