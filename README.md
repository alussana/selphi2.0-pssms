# SELPHI<sub>2.0</sub> PSSM features and prediction analysis

Compute kinase-phosphosite PSSM scores to be used as features in the SELPHI<sub>2.0</sub> model.

Analyze SELPHI<sub>2.0</sub> predictions to generate performance plots against an independent, experimentally supported test set, and to perform hierarchical clustering on the kinases based on the similarity of their predicted substrates. 

> Research paper: [https://doi.org/10.1101/2022.01.15.476449](https://doi.org/10.1101/2022.01.15.476449)

## Build the container image

Requires [Docker](https://www.docker.com) and [Apptainer](https://apptainer.org).

```bash
docker build -t selphi_2.0 - < env/Dockerfile
docker save -o env/selphi_2.0.tar.gz selphi_2.0
singularity build env/selphi_2.0.sif docker-archive://env/selphi_2.0.tar.gz
```

## Customise `nextflow.config`

Modify `process.executor`, `process.queue`, `workDir`, and `env` variables according to the infrastructure where the workflow will be executed. Find the Nextflow [configuration file](https://www.nextflow.io/docs/latest/config.html) documentation.


## PSSMs

This workflow generates PSSMs scores associated to kinase-substrate pairs in

```bash
${env.out_dir}/pssm_score/k_p_ser_thr_pssm_scores.tsv
${env.out_dir}/pssm_score/k_p_tyr_pssm_scores.tsv
```

All other features of kinase-substrate pairs are taken from

```bash
${env.selphi_2_features_table}
```

All features are merged together in

```bash
${env.out_dir}/features_table/k_p_features.tsv
```

## Performance analysis

The workflow computes ROC and PR curves on sets of experientally supported kinase-phosphosite associations ([Sugiyama et al. 2019](https://doi.org/10.1038/s41598-019-46385-4)), obtaining kinase family-specific performance evaluations.

SELPHI<sub>2.0</sub> predictions are taken from

```bash
${env.selphi2_prediction_matrix_dir}/prediction_matrix.csv
```

The results will be saved in 

```bash
${env.out_dir}/selphi2_sugiyama_100_rand_neg_sets/
```

## Kinase clustering

Radial dendrograms are used to display the predicted hierarchical similarity between kinases based on SELPHI<sub>2.0</sub>-assigned substrates.

Substrates identity-based similarity requires

```bash
${env.selphi2_prediction_matrix_dir}/prediction_matrix_high_conf
```

The results are exported at

```bash
${env.out_dir}/selphi2_kinase_dendrogram/
```

## Run the workflow

```bash
nextflow run main.nf -c nextflow.config -resume -with-dag misc/flowchart.svg
```

## Misc

Optionally, linear classifiers of kinase-phosphosite interaction using uniquely the PSSM score predictor are built and evaluated computing ROC and PR curves a hundred times, each time using a different randomly sampled negative set of examples 10 times larger than the positive set. ROC and PR curves can be found in:

```bash
${env.out_dir}/pssm_model_100_rand_neg_sets/*.pdf
```

The same is done with the Phosphormer model:

```bash
${env.out_dir}/phosformer_model_100_rand_neg_sets/*.pdf
```
