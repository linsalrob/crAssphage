# Metagenome Metadata

## Getting the data

Recall that the SRA is set up as projects (SRP) have one or more samples. Each sample has one or more experiments (SRX), and each experiment has one or more runs (SRR).

The metadata with a project is asssociated with the *sample*, and not with the run. Therefore, we need to get the sample list, parse the XML about the metadata and then connect it to the run.

A list of [runs that have crAssphage](sample_ids.txt) and the sample IDs associated with those runs. The third column has some alternate IDs for some of the runs.

## Getting the metadata

The BioSample [XML](ftp://ftp.ncbi.nlm.nih.gov/biosample/) is available from the NCBI ftp site, where you can download a single XML file that has all the data. 

We used the list of sample IDs above to create a small version fo the BioSample XML called [crAssphage.xml](crAssphage.xml) that just has the XML records for the samples that have crAssphage in them.

## Parsing the XML

We parsed the [crAssphage XML](crAssphage.xml) to extract all the metadata associated with the runs that have hits to crAssphage. The parsing was done with two custom Python scripts. The [first script](https://github.com/linsalrob/EdwardsLab/blob/master/sra/sra_xml_print_all_attributes.py) just prints the attributes in the XML file. You can print all tags, or just those you are interested in (this is redundant, so be sure to pass it through sort -u!). The [second script](https://github.com/linsalrob/EdwardsLab/blob/master/sra/sra_xml.py) is essentially a custom written XML parser just for this data, although you can easily adapt it to print data from any biosample set.

The parsing scripts generate a [crAssphage metadata](crAssphage_metadata.tsv) file that is tab separated values. You can import that into excel. There are 137 columns in that file, but many of them are empty of course. A 0-indexed list of columns is shown below. 

A couple of key columns are (these are 0-indexed):

* collection date is column 29
* country is in geo_loc_name column 51
* latitude/longitude is in lat_lon column 80.


## Metadata files

There are two metadata files specifically for country and date:
1. [crAssphage.country.date.tsv](crAssphage.country.date.tsv) has all 137 fields listed below, but only has the sample IDs
2. [runs_country_date.tsv](runs_country_date.tsv) has only four fields: Run ID, Country, Date of Collection, and Lat Lon. At least country and date are populated, but lat lon may not be



# All metadata attributes

0. Sample Accession
1. Ids
2. Description - Title
3. Description - Comment
4. Owner - Name
5. Owner - Email
6. Release Date
7. Links
8. abs_air_humidity
9. age
10. air_temp
11. alkalinity
12. altitude
13. analyte_type
14. biochem_oxygen_dem
15. biomaterial_provider
16. biospecimen_repository
17. biospecimen_repository_sample_id
18. body_habitat
19. body_mass_index
20. body_product
21. breed
22. building_setting
23. build_occup_type
24. carb_dioxide
25. chem_administration
26. chem_oxygen_dem
27. clone
28. collected_by
29. collection_date
30. cultivar
31. depth
32. description
33. dev_stage
34. diet
35. disease
36. elev
37. env_biome
38. env_feature
39. env_material
40. env_package
41. ethnicity
42. family_relationship
43. filter_type
44. gap_accession
45. gap_consent_code
46. gap_consent_short_name
47. gap_sample_id
48. gap_subject_id
49. gastrointest_disord
50. genotype
51. geo_loc_name
52. health_state
53. heat_cool_type
54. host
55. host_age
56. host_body_mass_index
57. host_body_product
58. host_body_temp
59. host_diet
60. host_disease
61. host_family_relationship
62. host_genotype
63. host_height
64. host_last_meal
65. host_occupation
66. host_phenotype
67. host_pulse
68. host_sex
69. host_subject_id
70. host_taxid
71. host_tissue_sampled
72. host_tot_mass
73. ihmc_medication_code
74. indoor_space
75. investigation_type
76. isolate
77. isolation_source
78. isol_growth_condt
79. label
80. lat_lon
81. light_type
82. liver_disord
83. medic_hist_perform
84. misc_param
85. molecular_data_type
86. nitrate
87. occupant_dens_samp
88. occup_samp
89. organism_count
90. oxy_stat_samp
91. perturbation
92. ph
93. phosphate
94. pre_treatment
95. project_name
96. propagation
97. race
98. reactor_type
99. ref_biomaterial
100. rel_air_humidity
101. rel_to_oxygen
102. salinity
103. samp_collect_device
104. sample_name
105. sample_type
106. samp_mat_process
107. samp_salinity
108. samp_size
109. samp_store_dur
110. samp_store_loc
111. samp_store_temp
112. samp_vol_we_dna_ext
113. sewage_type
114. sex
115. sludge_retent_time
116. smoker
117. source_material_id
118. space_typ_state
119. special_diet
120. store_cond
121. strain
122. study_design
123. study_disease
124. study_name
125. subject_is_affected
126. submitted_sample_id
127. submitted_subject_id
128. submitter_handle
129. suspend_solids
130. temp
131. tissue
132. tot_phosphate
133. treatment
134. typ_occupant_dens
135. ventilation_type
136. wastewater_type

