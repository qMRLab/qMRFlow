FROM neurodebian:xenial
MAINTAINER qmrlab

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -y install ants git curl unzip fsl-5.0-core

RUN /bin/bash -c "source /etc/fsl/5.0/fsl.sh"

ENV PATH=$PATH:/usr/lib/ants:/usr/lib/fsl/5.0 \
	FSLDIR=/usr/share/fsl/5.0 \
	LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/fsl/5.0 \
	FSLOUTPUTTYPE=NIFTI_GZ \
	FSLTCLSH=/usr/bin/tclsh