[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/3316)


# BEst computational pipeline for CBRAIN

## How we built the BEst container?

### 1. Compile BEst
Using [best_compile.bash](for_build/best_compile.bash):

```bash
~/programs/best-cbrain/for_build/best_compile.bash -d ~/wkdir/best-190729 -n best-app -c matlab18b
```

#### Help

```bash
~/programs/best-cbrain/for_build/best_compile.bash -h
```

#### Main requirements

- Linux OS
- MATLAB compiler


### 2. Build the Docker container and push
Using [Dockerfile](for_build/containers/Dockerfile)

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


### 4. Test the image
```bash
singularity exec ~/programs/best-cbrain.simg /bin/bash -c 'BEst inputData $BEST_DATA_DIR/test-data.mat outputName ~/wkdir/best-190729/test-1 memMethod cMEM sensorsType EEG+MEG reconstructionWindow "0.7 0.71" baselineWindow "0 0.5"'
```
Check for the output file `test-1.mat`.
