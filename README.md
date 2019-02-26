# Bigram-PGK
Phosphoglycerylation prediction using evolutionary information of amino acids

Run the file Bigram_PGK.m to obtain statistical measures of the Bigram_PGK predictor. The algorithm runs on the train and test data to obtain the values. 

To obtain the results of existing predictors on the same train and test sets, run the codes iPGK_PseAAC.m, CKSAAP_PhoglySite.m, and Phogly_PseAAC.m  

Preprocessing_BigramPGK.m file was used to carry out preprocessing for Bigram_PGK predictor, which includes feature extraction for segment size of ±32, filtering, and contruction of train and test sets for 10-fold cross-validation. Everytime the algorithm is executed, it will have a slightly different combination of samples for train and test sets. Nevertheless, the performance will be similar.    
