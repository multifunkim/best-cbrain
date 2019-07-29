[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/3316)


# BEst computational pipeline for CBRAIN


------------

### Table of Contents
- [BEst computational pipeline for CBRAIN](#best-computational-pipeline-for-cbrain)
  * [How the BEst container was built](#how-the-best-container-was-built)
    + [1. Compile BEst](#1-compile-best)
      - [Details](#details)
      - [Main requirements](#main-requirements)
    + [2. Build the Docker container and push](#2-build-the-docker-container-and-push)
      - [Main requirements](#main-requirements-1)
    + [3. Build the Singularity container](#3-build-the-singularity-container)
    + [4. Test the image](#4-test-the-image)

------------


## How the BEst container was built

### 1. Compile BEst

Using [best_compile.bash](for_build/best_compile.bash):

```bash
~/programs/best-cbrain/for_build/best_compile.bash -d ~/wkdir/best-190729 -n best-app -c matlab18b
```

#### Details

```bash
~/programs/best-cbrain/for_build/best_compile.bash -h
```

#### Main requirements

- Linux OS
- MATLAB compiler


### 2. Build the Docker container and push

Using [Dockerfile](for_build/containers/Dockerfile):

```bash
cp ~/programs/best-cbrain/for_build/containers/Dockerfile ~/wkdir/best-190729
docker build --no-cache -f ~/wkdir/best-190729/Dockerfile -t multifunkim/best:latest ~/wkdir/best-190729
docker push multifunkim/best:latest
```

#### Main requirements

- Compiled BEst archive
- BEst test dataset


### 3. Build the Singularity container

Choose one of these commands:

- Using [Singularity](for_build/containers/Singularity):
```bash
sudo singularity build ~/programs/best-cbrain.simg ~/programs/best-cbrain/for_build/containers/Singularity
```

- From Docker Hub:
```bash
sudo singularity build ~/programs/best-cbrain.simg docker://multifunkim/best:latest
```


### 4. Test

```bash
singularity exec ~/programs/best-cbrain.simg /bin/bash -c 'BEst inputData $BEST_DATA_DIR/test-data.mat outputName ~/wkdir/best-190729/test-1 memMethod cMEM sensorsType EEG+MEG reconstructionWindow "0.7 0.71" baselineWindow "0 0.5"'
```
Check for the output file `test-1.mat`.
