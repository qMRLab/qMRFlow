```
           _             _    __ _               
          | |           | |  / _| |              
 _ __ ___ | |_ ___  __ _| |_| |_| | _____      __
| '_ ` _ \| __/ __|/ _` | __|  _| |/ _ \ \ /\ / /
| | | | | | |_\__ \ (_| | |_| | | | (_) \ V  V / 
|_| |_| |_|\__|___/\__,_|\__|_| |_|\___/ \_/\_/  

A Nextflow pipeline for processing MTSAT data in BIDS format using qMRLab.
```
## Use with Docker

1. Install nextflow as decribed in [here](http://nextflow.io)


2. Pull the following Docker images:
```
docker pull qmrlab/minimal:v2.3.1
```

```
docker pull qmrlab/antsfsl:latest
```

3. In the `/mt_sat/nextflow.config` file, some paramaters must be set to null and
the platform to `"octave"` as shown below: 

```
platform="octave"
    
matlab_path = null

octave_path = null

qmrlab_path = null

wrapper_version = "v1.0.0" 

matlab_path_exception = null

qmrlab_path_exception = null

```

4. Set process-specific parameters in the `nextflow.config` file.
5. Ensure that you have the Docker images by running the following command in 
your terminal: 
```
docker images
```
6. Run the pipeline: 

Data must be organized in [MTsat BIDS format](https://github.com/qMRLab/qMRFlow/blob/master/mt_sat/USAGE).

```
cd /qMRFlow/mt_sat 
nextflow run mtsatflow_BIDS.nf --root $path_to_your_dataset -with-report report.html
```
