### Docker collection for qMRFlows 

Built images are are available at [qMRLab's public Docker image registry](https://hub.docker.com/u/qmrlab). 

Docker images of `qmrlab/mcrgui`, `qmrlab/octjn` and  `qmrlab/minimal` are built at each qMRLab release and respectively tagged. `Dockerfiles` to these images can be found [here](https://github.com/qMRLab/qMRLab/tree/master/Deploy/Docker).  

Below section keeps a log of which Docker image is associated with which `process` of which workflow. 

#### qMRFlow: [`mt_sat` pipeline](https://github.com/qMRLab/qMRFlow/tree/master/mt_sat)
 * `image`: qmrlab/antsfsl:latest 
      * `process`: Extract_Brain
      * `process`: Align_Inpt_Volumes 
 * `image`: qmrlab/minimal:v2.3.1
 	  * `process`: Fit_MTsat_*
