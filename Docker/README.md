### Docker collection for qMRflows 

Built images are are available at [qMRLab's public Docker image registry](https://hub.docker.com/u/qmrlab). 

Docker images of `qmrlab/mcrgui` and `qmrlab/octjn` are built at each qMRLab release and respectively tagged. `Dockerfiles` to these images can be found [here](https://github.com/qMRLab/qMRLab/tree/master/Deploy/Docker).  

Below section keeps a log of which Docker image is associated with which `process` of which workflow. 

#### qMRflow: [`mt_sat` pipeline](https://github.com/qMRLab/qMRflow/tree/master/mt_sat)
* `process`: Align_And_Extract
     * `image`: qmrlab/antsfsl:latest 
