```
 _      _         __ _   __ _               
| |    /_|       / _|_| / _| |              
| |__ /__| __ _ | |_ _ | |_| | _____      __
|  _ \ | |/ _` ||  _| ||  _| |/ _ \ \ /\ / /
| |_ | | | (_| || | | || | | | (_) \ V  V /  
|_,__| |_|\__,_||_| |_||_| |_|\___/ \_/\_/  

A Nextflow pipeline for processing b1afi data in BIDS format using qMRLab.
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

Data must be organized in [TB1map BIDS format](https://github.com/qMRLab/qMRFlow/blob/master/b1_afi/USAGE).

```
cd /qMRFlow/b1_afi
nextflow run b1afiflow_BIDS.nf --root $path_to_your_dataset -with-report report.html
```
