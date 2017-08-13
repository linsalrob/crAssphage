# Metagenome Metadata

## Getting the data

Recall that the SRA is set up as projects (SRP) have one or more samples. Each sample has one or more experiments (SRX), and each experiment has one or more runs (SRR).

The metadata with a project is asssociated with the *sample*, and not with the run. Therefore, we need to get the sample list, parse the XML about the metadata and then connect it to the run.

A list of [runs that have crAssphage](sample_ids.txt) and the sample IDs associated with those runs. The third column has some alternate IDs for some of the runs.

## PArsing the XML

We have parsed the BioSample [XML](ftp://ftp.ncbi.nlm.nih.gov/biosample/) from NCBI to extract all the metadata associated with the runs that have hits to crAssphage. The parsing was done with two custom Python scripts. The [first script](https://github.com/linsalrob/EdwardsLab/blob/master/sra/sra_xml_print_all_attributes.py) just prints the attributes in the XML file. You can print all tags, or just those you are interested in (this is redundant, so be sure to pass it through sort -u!). The [second script](https://github.com/linsalrob/EdwardsLab/blob/master/sra/sra_xml.py) is essentially a custom written XML parser just for this data, although you can easily adapt it to print data from any biosample set.


## Combining the output 

