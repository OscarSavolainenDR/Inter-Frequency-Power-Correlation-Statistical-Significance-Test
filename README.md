# Inter-Frequency-Power-Correlation-Statistical-Significance-Test

Author: Oscar Wiljam Savolainen – NGNI Lab, Imperial College London

Date: 09/08/2021


The directory “Reproducing_results_from_article” contains the code associated with the paper 'The Significance of Neural Inter-Frequency Power Correlations', doi: 10.21203/rs.3.rs-329644/v1
It looks at applying a statistical significance test for the multiple testing of inter-frequency power correlations on a large, publicly available dataset of intracortical broadband neural recordings (Macaque M1). It produces all of the visualization plots from the paper, e.g. the statistically significant inter-frequency power correlations for various recordings.  

![image](https://github.com/OscarSavolainenDR/Inter-Frequency-Power-Correlation-Statistical-Significance-Test/assets/119876479/c4e344bf-5772-4e95-8ae3-3a5889064bd0)

![image](https://github.com/OscarSavolainenDR/Inter-Frequency-Power-Correlation-Statistical-Significance-Test/assets/119876479/61a71041-4c1d-4f16-97a3-b52cdb25b584)



The directory “Analyzing_an_arbitrary_signal” is a small directory for the application of the same statistical test for an arbitrary signal, e.g. for a simple signal with two time-overlapping sinusoidal bursts at 2 different frequencies: (a) time-series signal; (b) time-frequency distribution (Continuous Wavelet Transform); (c) inter-frequency power correlation matrix, with non-statistically significant elements blacked out.

![image](https://github.com/OscarSavolainenDR/Inter-Frequency-Power-Correlation-Statistical-Significance-Test/assets/119876479/20be26d1-12ff-4b05-b74a-41bade5946d3)


It uses functions from “Reproducing_results_from_article”, and performs best for signals of length 1e6 or smaller, given the realistic memory constraints of most desktop computers. If a longer signal needs to be analysed, it is recommended to run the script on a High Performance Computing cluster, e.g. with the elements of the null distributions computed in parallel. 

Please don't hesitate to get in contact if you have any questions.

