# Data from Norway

Email from 6/1/16

I hope it is not too late to submit more Crassphage sequences; this time these are from Norway, from a metagenomic virome sequencing.

Please find:
 (a) a table with basic sample data
 (b) sam files and their reference

This is from a limited set of viromes of children longitudinally followed for their high genetic risk of diabetes (our colleagues will have to be included among authors if this is used in papers). These positivities are 3 children: one with two samples (3 months apart), one with 3 samples (3 and 6 months apart, and please see how young the infants are), and finally the last one contributed with two samples (two months apart). Most samples were sequenced from two different stool preparations (just filtered supernatant, or ultracentrifuged virus-like particle prep). I guess this short distance in time between samples from the same individual may be interesting with respect to longitudinal evolvement of the sequence.

And - thanks to the efforts of our kind collaborators from Nigeria and Jordan - we have some 50 samples from these non-Europoid populations sitting in our freezers. Will try soon whether there is any CrAssphage.


individual ID code | sex | date of birth | date of sample | age at sample yrs | bin_id | crassph reads per 100,000 | sample | sample_processing
--- | --- | --- | --- | --- | --- | --- | --- | ---
2848 | F | 10/10/02 | 10/04/03 | 0.50 | 110 | 25643 | 683 | filtering + pelleting + enzymes
2848 | F | 10/10/02 | 10/07/03 | 0.75 | 108 | 359 | 965 | only filtering
2848 | F | 10/10/02 | 10/07/03 | 0.75 | 112 | 3026 | 965 | filtering + pelleting + enzymes
4525 | F | 18/04/03 | 18/07/03 | 0.25 | 133 | 5436 | 998 | only filtering
4525 | F | 18/04/03 | 18/07/03 | 0.25 | 137 | 6633 | 998 | filtering + pelleting + enzymes
4525 | F | 18/04/03 | 18/10/03 | 0.50 | 134 | 2301 | 1327 | only filtering
4525 | F | 18/04/03 | 18/10/03 | 0.50 | 138 | 10020 | 1327 | filtering + pelleting + enzymes
4525 | F | 18/04/03 | 18/01/04 | 0.75 | 135 | 245 | 1686 | only filtering
4525 | F | 18/04/03 | 18/01/04 | 0.75 | 136 | 96 | 1686 | filtering + pelleting + enzymes
20685 | F | 19/01/06 | 29/07/07 | 1.52 | 163 | 721 | 17467 | only filtering
20685 | F | 19/01/06 | 29/07/07 | 1.52 | 183 | 3497 | 17467 | filtering + pelleting + enzymes
20685 | F | 19/01/06 | 30/09/07 | 1.69 | 184 | 207 | 18603 | filtering + pelleting + enzymes


## Analysis

We used our [strain identification pipeline](https://github.com/linsalrob/crAssphage/tree/master/Metagenomes/IdentifyingStrains) to identify the strains in these samples and they were named according to the metadata.
