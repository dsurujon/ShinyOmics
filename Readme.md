ShinyEDGE: Exploration of Differential Gene Expression
-----------------------------------------

This app allows exploratory data analysis on RNAseq data. It comes pre-loaded with data from Zhu et al 2019 *(in preparation)* where two strains of *Streptococcus pneumoniae* (19F and T4) are treated with 4 different antibiotics (LVX, KAN, VNC, RIF). To access the app with this dataset, visit http://bioinformatics.bc.edu/shiny/ShinyEDGE/

There are 4 panels that allow for different types of data exploration: 
* **Single Experiment:** plot differential expression (DE) against any other metadata associated with genes (e.g. to answer whether essential genes are more downregulated, you can select Essentiality as the x-axis metadata variable)
* **Compare 2 Experiments:** plot DE from one experiment against the DE in another experiment (e.g. to answer whether two antibiotics trigger similar responses, plot T4_LVX against T4_VNC)
* **Compare All Experiments:** Use a heatmap and PCA to see if there are groups of highly similar experiments 
* **Network:** Overlay significant expression changes on a network, compare network characteristics (e.g. degree) to gene data (e.g. DE or other metadata)

### Single Experiment
In this panel, you can select one RNA-Seq experiment (which can include multiple timepoints) using the dropdown menu "Select experiment". The DE data for that experiment and corresponding metadata for the strain will be loaded, and the scatter plot will update accordingly. In this plot, the y-axis is always DE, and the x-axis is the metadata variable you may select from the dropdown menu "Variable". Transcriptionally Important Genes (TIGs) are defined as genes with significant DE. More formally, |log2FoldChange(experiment/control)|>1 and Bonferroni-adjusted p-value < 0.05. TIGs are colored in green, whereas non-TIGs are black.     
There are several visualization options. For categorical x-axis variables, checking "Jitter x axis" prevents overlap of points by adding a small amount of noise to the x-coordinate of each data point, and might make the visualization easier to interpret. It is also possible to overlay a violin plot with mean +- 95% confidence intervals (in red).   
![Jitter Example](https://contattafiles.s3.us-west-1.amazonaws.com/tnt8877/qK2wZmkn9549tMt/Pasted%20Image%3A%20Apr%204%2C%202019%20-%2011%3A58%3A36am)    
Some metadata variables may follow a log-normal distribution. In order to accommodate for this, there is also a checkbox to log-transform the x-axis ("log-scale x axis"). Finally, the slider "Transparency" controls the transparency of the points. For instance, for experiments where a large number of points are plotted simultaneously, a value of ~0.5 is recommended.     
    
The scatter plot itself is brushable, meaning it is possble to select a subset of points by dragging a selection window directly on the plot. The selected genes will then be displayed on the table at the bottom of the screen. For example, in the Figure below we selected the genes that are highly upregulated in the functional tag "GENETIC INFORMATION PROCESSING" at 240 minutes: 
![Brush example](https://contattafiles.s3.us-west-1.amazonaws.com/tnt8877/CuF1dl5AyADSTrB/Pasted%20Image%3A%20Apr%204%2C%202019%20-%2012%3A56%3A25pm)    
    
    
It is also possible to only display a set of previously identified genes. To do this, type or paste a gene list in the "Paste gene list" text box. This will filter the dataset such that only these genes will be displayed
![Geneset example](https://contattafiles.s3.us-west-1.amazonaws.com/tnt8877/Q3Uoiu2tQGEBCfO/Pasted%20Image%3A%20Apr%204%2C%202019%20-%2012%3A58%3A29pm)    

### Compare 2 Experiments
This panel allows you to compare two experiments **from the same organism**. Both experiments need to have the same metadata file associated with them to be able to cross-reference the genes. The comparisons will be made for each timepoint, so another requirement here is that there is at least one overlapping time point in the two experiments.    
Use the dropdown lists to select the two experiments to compare (Experiment 1 will appear as the x-axis, and Experiment 2 will appear as the y-axis). A third dropdown list ("Variable") determines the color variable.     
In this plot we compare the Rifampicin (RIF) response to Kanamycin (KAN) in the strain T4. There are some Genetic Information Processing genes that are uniquely upregulated in KAN, and some Metabolism genes uniquely downregulated in RIF. Similar to panel 1, it is possible to brush this plot and update the table below (brushed at 120 minutes).      
![Panel2 Example](https://contattafiles.s3.us-west-1.amazonaws.com/tnt8877/vwJuQxV5TmEgvfv/Pasted%20Image%3A%20Apr%204%2C%202019%20-%201%3A06%3A38pm)

### Compare All Experiments
This panel allows a more global comparison across all conditions (for a given strain/organism). Select the strain you would like to visualize using the dropdown menu. The heatmap on the left will show DE from all experiments (columns) for all genes (rows). The plot on the right displays the first two principal components, where each point is an experiment/timepoint. The dropdown menu allows you to color the points by one variable on the experiment sheet. The plot below the PCA plot shows the % variance explained by each of the principal components.     
![Panel3 example](https://contattafiles.s3.us-west-1.amazonaws.com/tnt8877/DF2WzmD35c0RxSk/Pasted%20Image%3A%20Apr%209%2C%202019%20-%2011%3A08%3A12am)

### Network
Here, you can visualize any network (the preloaded data has a metabolic network and transcription regulatory network), and overlay TIGs. Use the corresponding dropdown selectors to selet the appropriate network, experiment and timepoint to visualize. The package ```visNetwork``` is used to generate the interactive network plot, where you can zoom in/out, move nodes around, and select one node to highlight its neighbors.     
The scatter plot on the right side can be used to explore how network characteristics relate to gene expression or metadata. The x-axis selector gives you the option of 3 network characteristics for each gene (Degree, Betweenness Centrality and Eigencentrality, see Table below), and the y-axis selector uses DE and any other metadata columns included in the metadata file. This figure is also brushable, and the brushed genes will appear in the table below.      
![Network example](https://contattafiles.s3.us-west-1.amazonaws.com/tnt8877/45I2SQM8byE9YP6/Pasted%20Image%3A%20Apr%204%2C%202019%20-%203%3A29%3A35pm)   
    
| Centrality measure | Definition                                                                                                                                                                                                                               |
|--------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Degree             | Number of neighbors a node has                                                                                                                                                                                                           |
| Betweenness        | Proportion of shortest paths between nodes that pass through this node. This can be interpreted as the amount of information flow that depends on this node.                                                                             |
| Eigencentrality    | Importance in the network, depending on how many important neighbors a node has. Two nodes of equal degree might differ in eigencentrality, if one has more high-degree neighbors (resulting in higher eigencentrality) than the other.  |


## Running the app locally
To run this app locally, you will need [R](https://www.r-project.org/), [RStudio](https://www.rstudio.com/products/rstudio/download/) and a few package dependencies. The required packages are ```ggplot2```, ```visNetwork```, ```igraph```, and ```shiny```. To install any of these packages, use the following command in RStudio (once for each package):    
```install.packages('shiny')```    
With ```server.R``` open in RStudio, hit the "Run App" button on the upper right corner.     

## Using custom data
To use custom data from other experiments in this app, all datasets need to be specified in the experiment sheet ```data/exptsheet.csv```. You can add new columns to specify metadata pertaining to each experiment. The minimum required columns are ```Name``` (which should be a unique identifier, we recommend using ```Experiment-Time``` as the format), ```Experiment``` (this should be the same for all timepoints of the same experiment),```DEseqFile``` (the csv file corresponding to that experimental datapoint) , ```Strain``` (strain or organism name), ```MetadataFile``` (the csv file corresponding to the strain), ```Time``` (timepoint of the RNAseq experiment).     
The app then refers to the ```DEseqFile``` and ```MetadataFile``` columns to find the appropriate data. Make sure all data and metadata files are in the file locations specified.     
#### Metadata Files
There should be one metadata file per strain. The metadata files should have at least a ```Gene``` column, and the gene labels under this column should match those in the DESeq data files and network files. You can add any number of additional columns, and the selectors in the app will update accordingly.     
#### DEseq Data Files
The data files used are the standard output of [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). The requirement for the app is to have at least 3 columns: ```Gene```, ```log2FoldChange``` (DE), and ```padj``` (adjusted p-value) in each of the data files.     
#### Network Files
Network files need to be specified in the ```data/network``` subdirectory, and should contain edge tables. The file name should end with ```_Edges.csv```, and the table should have two columns ```source``` and ```target```. Make sure the gene labels match those in the metadata files and RNAseq data files.     
