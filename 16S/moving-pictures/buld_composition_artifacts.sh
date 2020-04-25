#!/usr/env bash

qiime composition add-pseudocount \
      --i-table gut-table.qza \
      --o-composition-table comp-gut-table.qza

#To build the composition artifact, a FeatureTable[Frequency] artifact must be provided to add-pseudocount (an imputation method), which will produce the FeatureTable[Composition] artifact.
