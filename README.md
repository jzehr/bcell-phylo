# Build and visualize trees from antibody JSON data #

## Installation

### Requirements

- SnakeMake
- Anaconda  
- Python 
- [mafft](https://mafft.cbrc.jp/alignment/software/)
- [FastTree](http://www.microbesonline.org/fasttree/)
- [yarn](https://yarnpkg.com/en/)

### Install instructions

```
conda env create -f environment.yml
source activate bcell
yarn
```

## Pipeline

Obtain a copy of the compressed input data, `bcell-phylo_Ver4.tar.gz`, and place in the `data` directory via:

```
mv /path/to/bcell-phylo_Ver4.tar.gz data/
```

After installing requirements, run the pipeline from the bcell-phylo directory:

```
snakemake -j $NUMBER_OF_CONCURRENT_JOBS $TARGET
```

Or, to distribute jobs with Sun Grid Engine:


```
snakemake $TARGET --cluster "qsub -V" -j $NUMBER_OF_CONCURRENT_JOBS
```

Make sure that the `python` executable for the `conda` environment is on your `$PATH`.

## Visualization (development)

After running the pipeline:

```
yarn start
```

To visualize the steps of the pipeline:

```
yarn dag
```

## Deployment

### Sync JSON files for viewing on dashboard
From location where pipeline has been run (usually head node of a cluster):

```
rsync -arvm --include "*/V*.json" --include="*/" --exclude="*" data/ $USER@$HOST:$BCELLPHYLO/dist/data
```

### Webserver

```
yarn build
npx http-server -p $PORT dist
```