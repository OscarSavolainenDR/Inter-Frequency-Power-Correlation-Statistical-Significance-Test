# Inter-Frequency-Power-Correlation-Statistical-Significance-Test

Author: Oscar Wiljam Savolainen – NGNI Lab, Imperial College London
Date: 09/08/2021

The directory “Reproducing_results_from_article” contains the code associated with the paper 'The Significance of Neural Inter-Frequency Power Correlations', doi: 10.21203/rs.3.rs-329644/v1
It looks at applying a statistical significance test for the multiple testing of inter-frequency power correlations on a large, publicly available dataset of intracortical broadband neural recordings (Macaque M1).

The directory “Analyzing_an_arbitrary_signal” is a small directory for the application of the same statistical test for an arbitrary signal. It uses functions from “Reproducing_results_from_article”, and performs best for signals of length 1e6 or smaller, given the realistic memory constraints of most desktop computers. If a longer signal needs to be analysed, it is recommended to run the script on a High Performance Computing cluster, e.g. with the elements of the null distributions computed in parallel.

