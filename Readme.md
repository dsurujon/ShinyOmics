![GitHub](https://img.shields.io/github/license/dsurujon/ShinyOmics?color=Green&style=plastic)    

ShinyOmics User Guide
-----------------------------------------

This app allows exploratory data analysis on any 'Omics' dataset, e.g. Transcriptomics, proteomics, phenotypic screens. It comes with pre-loaded datasets from the human pathogens *Streptococcus pneumoniae* and *Mycobacterium tuberculosis*. The *S. pneumoniae* set includes RNA-Seq and Tn-seq data from Zhu et al 2019 *(under review)* where two strains (19F and T4) were treated with 5 different antibiotics (KAN, LVX, VNC, RIF, PEN). The *M. tuberculosis* datasets are microarray and proteomic screens under hypoxic conditions from [Galagan et al., 2013](https://www.nature.com/articles/nature12337) and [Schubert et al., 2015](https://www.sciencedirect.com/science/article/pii/S193131281500222X) respectively. To access the app with this example dataset, visit http://bioinformatics.bc.edu/shiny/ShinyOmics/

There are 4 panels that allow for different types of data exploration: 
* **Single Experiment:** plot experimental value (e.g. differential expression (DE), change in fitness (dW)) of all genes against any other metadata associated with genes (e.g. to answer whether essential genes are more downregulated, we can select Essentiality as the x-axis metadata variable)
* **Compare 2 Experiments:** plot gene value (DE, dW, etc.) from one experiment against the value in another experiment (e.g. to answer whether two antibiotics trigger similar responses, plot T4_PEN against T4_VNC)
* **Compare All Experiments:** Use a heatmap and PCA to see if there are groups of highly similar experiments 
* **Network:** Overlay significant changes on a network, compare network characteristics (e.g. degree) to gene data (e.g. DE or other metadata)    
    
    
Before the app attempts any data processing, it will check whether all files are correctly formatted. If there are issues, it will display error messages informing the user what the problem is. Otherwise, it will display a message saying the validation step was passed without issues    
![Validation Example](https://contattafiles.s3.us-west-1.amazonaws.com/tnt8877/EFREHPJ1WpI0TXF/Pasted%20Image%3A%20Sep%2030%2C%202019%20-%205%3A36%3A03pm)    

### Single Experiment
In this panel, we can select one RNA-Seq experiment (which can include multiple timepoints) using the dropdown menu "Select experiment". The DE data for that experiment and corresponding metadata for the strain will be loaded, and the scatter plot will update accordingly. In this plot, the y-axis is always DE, and the x-axis is the metadata variable we may select from the dropdown menu "Variable". Differentailly Expressed Genes (DEGs) are defined as genes with significant DE. More formally, |log2FoldChange(experiment/control)|>1 and Bonferroni-adjusted p-value < 0.05. DEGs are colored in green, whereas non-DEGs are black.     
There are several visualization options. For categorical x-axis variables, checking "Jitter x axis" prevents overlap of points by adding a small amount of noise to the x-coordinate of each data point, and might make the visualization easier to interpret. It is also possible to overlay a violin plot with mean +- 95% confidence intervals (in red).   
![Jitter Example](https://contattafiles.s3.us-west-1.amazonaws.com/tnt8877/qK2wZmkn9549tMt/Pasted%20Image%3A%20Apr%204%2C%202019%20-%2011%3A58%3A36am)    
Some metadata variables may follow a log-normal distribution. In order to accommodate for this, there is also a checkbox to log-transform the x-axis ("log-scale x axis"). Finally, the slider "Transparency" controls the transparency of the points. For instance, for experiments where a large number of points are plotted simultaneously, a value of ~0.5 is recommended.     
    
If there are a large number of timepoints, the user can select which ones to display with the checkboxes on the right. The scatter plot itself is brushable, meaning it is possble to select a subset of points by dragging a selection window directly on the plot. The selected genes will then be displayed on the table at the bottom of the screen. For example, in the Figure below we selected the genes that are highly downregulated in the functional tag "METABOLISM" at 120 minutes: 
![Brush example](https://contattafiles.s3.us-west-1.amazonaws.com/tnt8877/stTGLkl9nUGpVf1/Pasted%20Image%3A%20Oct%208%2C%202019%20-%206%3A39%3A49pm)    
    
The table generated can be sorted by a specific column, filtered using a search term that applies to the entire table (using the text box on the upper right corner of the table), or using a search term that applies to a specific column (using the search boxes at the bottom of each individual column). For columns with numeric values, we recommend sorting the table by these columns rather than searching by the column. This is because searching through numeric columns will treat the search term as a character string instead of an exact numeric value, and retrieve anything that has the specified character as part of the string (e.g. typing "3" under SequencePrevalence will retrieve "371" as well as "3").      
    
Brushed points will be displayed as a gene list in the text box on the right (Under the heading "Brished Genes"), so that the user can easily copy this gene set and paste it into the gene selector of another tab. The brushing can be reset by clicking anywhere on the plot.     
	
It is also possible to only display a set of genes of interest. To do this, type or paste a gene list in the "Paste gene list" text box. This will filter the dataset such that only these genes will be displayed
![Geneset example](https://contattafiles.s3.us-west-1.amazonaws.com/tnt8877/6UujTZ8qitxJzIO/Pasted%20Image%3A%20Oct%208%2C%202019%20-%206%3A41%3A17pm)    
    
Another way of subsetting the genes to display is to use the "Select genes by metadata variable" checkbox. When checked, this will generate a new set of selectors that allow the user to only display genes that have certain properties. The user can first select which variable to subset on (for example, sequence prevalence, which indicates how many strains share this gene), and then use the slider to subset only very common genes (present in > 335 strains)     
![Subset by prevalence example](https://contattafiles.s3.us-west-1.amazonaws.com/tnt8877/HG8Pd6KvBrluDdi/Pasted%20Image%3A%20Oct%208%2C%202019%20-%206%3A46%3A46pm)
    
The table, and the figure can be downloaded using the download buttons on the app.     

### Compare 2 Experiments
This panel allows us to compare two experiments **from the same organism**. Both experiments need to have the same metadata file associated with them to be able to cross-reference the genes. The comparisons will be made for each timepoint, so another requirement here is that there is at least one overlapping time point in the two experiments.    
Use the dropdown lists to select the two experiments to compare (Experiment 1 will appear as the x-axis, and Experiment 2 will appear as the y-axis). A third dropdown list ("Variable") determines the color variable.     
In this plot we compare the Penicillin (PEN) response to Vancomycin (VNC) in the strain T4. The overall expression changes seem to correlate in the 90-minute timepoint. There are also a number of METABOLISM genes that are downregulated in both. Similar to panel 1, it is possible to brush this plot and update the table below (brushed at 90 minutes). It is also possible to only display a subset of genes, download the table, or the figure.       
![Panel2 Example](https://contattafiles.s3.us-west-1.amazonaws.com/tnt8877/1CdZB4aoPWtnNkC/Pasted%20Image%3A%20Oct%208%2C%202019%20-%206%3A51%3A40pm)

### Compare All Experiments
This panel allows a more global comparison across all conditions (for a given strain/organism). Select the strain we would like to visualize using the dropdown menu. The heatmap on the left will show the experimental value (e.g. expression, fitness, protein abundance) from all experiments (columns) for all genes (rows). If we want, we can subset this heatmap by pasting a gene list in the text box below, and a second, interactive heatmap will be displayed. The plot on the right displays a PCA where we can select which 2 components to plot - here, each point is an experiment/timepoint. Another dropdown menu allows us to color the points by one variable on the experiment sheet. The plot below the PCA plot shows the % variance explained by each of the principal components.     
![Panel3 example](https://contattafiles.s3.us-west-1.amazonaws.com/tnt8877/3agNbO0czFDuJDW/Pasted%20Image%3A%20Sep%2030%2C%202019%20-%205%3A27%3A30pm)

### Network
Here, we can visualize any network (the preloaded data has a metabolic network and transcription regulatory network), and overlay DEGs (up/down-regulated genes are red/blue respectively). Use the corresponding dropdown selectors to selet the appropriate network, experiment and timepoint to visualize. The package ```visNetwork``` is used to generate the interactive network plot, where we can zoom in/out, move nodes around, and select one node to highlight its neighbors.     
The scatter plot on the right side can be used to explore how network characteristics relate to gene expression or metadata. The x-axis selector gives us the option of 3 network characteristics for each gene (Degree, Betweenness Centrality and Eigencentrality, see Table below), and the y-axis selector uses DE and any other metadata columns included in the metadata file. This figure is also brushable, and the brushed genes will appear in the table below.      
![Network example](https://contattafiles.s3.us-west-1.amazonaws.com/tnt8877/AAet6PhfywTydse/Pasted%20Image%3A%20Oct%208%2C%202019%20-%206%3A55%3A00pm)   
    
| Centrality measure | Definition                                                                                                                                                                                                                               |
|--------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Degree             | Number of neighbors a node has.                                                                                                                                                                                                           |
| Betweenness        | Proportion of shortest paths between nodes that pass through this node. This can be interpreted as the amount of information flow that depends on this node.                                                                             |
| Eigencentrality    | Importance in the network, depending on how many important neighbors a node has. Two nodes of equal degree might differ in eigencentrality, if one has more high-degree neighbors (resulting in higher eigencentrality) than the other.  |


## Running the app locally
To run this app locally, we need [R](https://www.r-project.org/), [RStudio](https://www.rstudio.com/products/rstudio/download/) and a few package dependencies. The required packages are ```ggplot2```, ```visNetwork```, ```igraph```, ```RColorBrewer```, ```heatmaply```, ```shinyHeatmaply```, and ```shiny```. To install any of these packages, use the following command in RStudio (once for each package):    
```install.packages('shiny')```    
With ```server.R``` open in RStudio, hit the "Run App" button on the upper right corner.     
    
## Using Shiny Server	
It is also possible to set up an instance of the application on a webserver. Please refer to the [Shiny Server Documentation](https://shiny.rstudio.com/articles/shiny-server.html) for details on how to host shiny applications. Setting up such a webserver instance *will* make the data publicly available. Therefore, if there are concerns about external access to data, ShinyOmics and its accompanying data can be distributed as zipped files, hosted behind firewalls, or the server can be set up with user authentication (which is only possible for Shiny Server Professional subscribers).      
        
## Directory Structure for ShinyOmics
![DirectoryStructure](https://contattafiles.s3.us-west-1.amazonaws.com/tnt8877/ifFyurDlypPkfrk/Pasted%20Image%3A%20Oct%208%2C%202019%20-%207%3A01%3A13pm)    
* All experiments need to be specified in the main experiment sheet, ```exptsheet.csv``` 
* For each organism/strain, there needs to be a metadata file – This file is a table of locus tags and any metadata associated with each locus tag, under ```metadata/```
* There can be as many network files as desired, under ```network/```
* The Omics data could be organized in any way, as long as the file locations are specified in ```exptsheet.csv```. In this case, since the data comes from microarray, proteomics, RNA-Seq, and Tn-seq experiments and from multiple strains, the directory is split first by experiment type, and in the case of RNAseq,  by strain, then by drug. There are several csv files (similar to DESeq2 or limma output) in each of the sub-directories. 


## Using custom data
To use custom data from other experiments in this app, all datasets need to be specified in the experiment sheet ```data/exptsheet.csv```. We can add new columns to specify metadata pertaining to each experiment. The minimum required columns are ```Name``` (which should be a unique identifier, we recommend using ```Experiment-Time``` as the format), ```Experiment``` (this should be the same for all timepoints of the same experiment),```DataFile``` (the csv file corresponding to that experimental datapoint) , ```Strain``` (strain or organism name), ```MetadataFile``` (the csv file corresponding to the strain), ```Time``` (timepoint of the RNAseq experiment).     
The app then refers to the ```DataFile``` and ```MetadataFile``` columns to find the appropriate data. Make sure all data and metadata files are in the file locations specified.     
#### Metadata Files
There should be one metadata file per strain. The metadata files should have at least a ```Gene``` column, and the gene labels under this column should match those in the DESeq data files and network files. We can add any number of additional columns, and the selectors in the app will update accordingly.     
#### Data Files
The data files used can be the standard output of [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). The requirement for the app is to have at least 3 columns: ```Gene```, ```Value``` (DE, fitness, protein abundance...), and ```padj``` (adjusted p-value) OR ```Sig``` (significant value) in each of the data files. If ```padj``` is provided, the app will generate a new ```Sig``` column based on ```Value``` and ```padj```.      
#### Network Files
Network files need to be specified in the ```data/network``` subdirectory, and should contain edge tables. The file name should end with ```_Edges.csv```, and the table should have two columns ```source``` and ```target```. Make sure the gene labels match those in the metadata files and RNAseq data files.     
