```
             
 _ __ ___  __ _  _ ___  _ __  ___ _   ___,  ____,
| '_ ` _ \|  _ \| _' _|| '__|/  _` |/ _   \/    _\
| | | | | | |_ | '/ /_ | |   \ (_| || (_) || |^^_ | 
|_| |_| |_| __,||____| |_|    \__,_|\___, |\_____/ 
          | |                        ___| |
          |_|                       |_____|
A Nextflow pipeline for processing MP2RAGE data in BIDS format using qMRLab.
```
## Use with Docker in 3 steps

1. Install nextflow as decribed in [here](http://nextflow.io)


2. Pull the following Docker images:
```
docker pull qmrlab/minimal:v2.3.1
```
```
docker pull qmrlab/antsfsl:latest
```
3. Run the pipeline: 

Data must be organized in [MP2RAGE BIDS format](https://github.com/qMRLab/qMRFlow/blob/master/mp2rage
/USAGE).

```
cd /qMRFlow/mp2rage
nextflow run mp2rageflow_BIDS.nf --root $path_to_your_dataset -with-report report.html
```
