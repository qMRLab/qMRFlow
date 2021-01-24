```
          __       _   __         
         / _|     | | /, |         
 _    _ | |_  __ _| |_ | |
\ \  / /|  _|/ _` | __|| |
 \ \/ / | |  ,(_| | |  | |
  \__/  |_|  \__ ,|,|  |_|

A Nextflow pipeline for processing VFA T1 data in BIDS format using qMRLab.
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

Data must be organized in [VFA BIDS format](https://github.com/qMRLab/qMRFlow/blob/master/vfa_t1/USAGE).

```
cd /qMRFlow/vfa_t1 
nextflow run vfat1flow_BIDS.nf --root $path_to_your_dataset -with-report report.html
```
