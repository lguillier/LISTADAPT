## Introduction
Sampling is crucial for the pertinence/performance of the genomic analysis carried out. Several strategies of sampling are available and the analyst must adapt to the global objective they have. If the objective is to explore the diversity of strains circulating in a country, the analyst can use the metadata associated to strain to reach that objective. When there is more than two categories of metadata describing the strains, the selection requires an algorithm of selection. 
An objective and reproducible selection procedure for selecting strains to include in panel for source attribution, diversity characterization or outbreak investigation is a high importance.   

## Method os selection
A three-step process for selecting strains based on metadata information (e.g. region, type of animal…) is proposed in this training. This method relies on the Gower's coefficient (1971), which is a metric expressing a dissimilarity: the “distance” between two units is the sum of all the variable-specific distances (associated to metadata categories). GC metric is capable of combining numeric and categorical data. 
GC offers the opportunity to the analyst to select weights for each individual variable, effectively altering the importance of each metadata categories (region more important than year).      

The training is related to metadata2select.r script. It takes a csv file that includes strains ID and metadata information. It provides as output a csv file of selected strains. 
