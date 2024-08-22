#!/bin/bash

# DRAP is at http://www.sigenae.org/drap/index.html

docker pull sigenae/drap

docker create --name drap -v /DATA/rchavez/pparvum/drap:/docker/drap -v /DATA/rchavez/pparvum/PRJNA300303_fastq:/docker/fastq -i -t sigenae/drap:latest /bin/bash
docker start drap
docker exec -i -t drap /bin/bash

# a typical drap command line looks like this:
/usr/local/src/drap/runDrap -o /docker/drap/SRR1294411_oases -1 /docker/fastq/SRR1294411_1.fastq -2 /docker/fastq/SRR1294411_2.fastq -l 100 -t contig -v 1 --no-rate --dbg-mem 90 --norm-mem 90 --alignTrim --cfg-file /docker/drap/drap.cfg

# after running drap on all samples run meta:
/usr/local/src/drap/runMeta --outdir PRJNA248394_meta --drap-dirs /docker/drap/SRR1294411_oases,/docker/drap/SRR1294412_oases,/docker/drap/SRR1296710_oases,/docker/drap/SRR1296769_oases,/docker/drap/SRR1296917_oases,/docker/drap/SRR1296973_oases,/docker/drap/SRR1300223_oases -y -l 100 -t contig -v 1 --no-rate --cfg-file /docker/drap/drap.cfg

# when done stop and delete the drap instance
docker stop drap
docker rm drap
