# Biorefinery Network Analysis
This repository contains the code for the analysis of the biorefineries' reaction space.

## Table of Contents
1. [How to run the code](#how-to-run-the-code)
2. [Functionality](#functionality)
3. [Data](#data)
4. [Querying REAXYS through browser](#querying-reaxys-through-browser)
5. [Testing](#testing)
6. [Authors](#authors)

## How to run the code
1. Create and activate [_venv_](https://docs.python.org/3/library/venv.html) with _Python 3.11_ and install dependencies from `requirements.txt` by running `pip install -r requirements.txt`.
2. Put your data in the `resources` folder. For more information about data have a look [here](#data).
3. Open `main.py` file, choose one of the below options and change the `option` variable to use the option of your choice.
   - "merge_and_graph_fragmentation" - merges data from `my_file` with the `file_merge` data. Then, performs
           "graph_fragmentation" analysis on such combined data.
   - "graph_fragmentation" - analyses the network fragmentation. Provides info about connected components like
       their size distribution. Creates HTML visualisation of the whole network.
   - "important_molecules" - analyses important molecules in the network according to the highest degree and
       betweenness centrality. Plots the results. Creates the HTML visualisation where the most important molecules
       are highlighted.
   - "degree_distribution" - analyses the degree distribution of the network. Creates plot specified as the `plot_type`
       and finds the exponent for the power law best fitting the data.
   - "graph_property" - calculates network metrics (average shortest path length, clustering coefficient, omega metric)
       and plots the average shortest path length as probability, compared to the Network of Organic Chemistry.
   - "graph_clusters" - performs k-core clustering and girvan-neuman clustering. Finds the best girvan-neuman
       clustering according to created metric and creates its HTML visualisation.
   - "degree_correlation" - calculates and plots the degree correlation for both in and out-degree.
   - "central_point_dominance" - calculates central point dominance
   - "dataset_comparison" - compares the dataset with the data from `file_compare` by calculating Jaccard similarity
       and checking whether important molecules from Biorefinery Network appear in the other dataset.
4. Run the `main.py` file and wait for the results.

## Functionality

`analyse("graph_fragmentation")` from `main.py` - analyses the network fragmentation. Provides info about connected components (and plots it) like their size distribution and creates interactive HTML visualisation of the whole network (`src/visualisations/graph_fragmentation.html`).

![image](https://github.com/user-attachments/assets/2df02eb1-5e9b-444f-8c97-7d17381d6527)


`analyse("graph_property")` from `main.py` - creates the comparison between Biorefinery Network and Network of Organic Chemistry, seen below. Note, it takes a few minutes to run since it also calculates omega metric.

<img src="https://github.com/user-attachments/assets/66006fb0-66f0-4c93-a00c-3001a1185f34" width="450">

`analyse("graph_clusters")` from `main.py` - finds the best partition of the network, using girvan-neuman algorithm and betweenness centrality as the metric. It also produces the plot below and the HTML visualisation, which can be find in `src/visualisations/graph_clusters_betweenness3.html`. Note, it also runs a few minutes.

<img src="https://github.com/user-attachments/assets/937d62a8-6f6a-49ef-afb7-8e6a5427a4ed" width="450">


Screenshot of the visualisation:

<img src="https://github.com/user-attachments/assets/e4e20879-4777-46ee-8c53-398ef1d3a9d3" width="450">


`analyse("dataset_comparison")` from `main.py` - calculates Jaccard similarity to the other dataset. Creates interactive HTML visualisation. Nodes from the first dataset are colored teal, nodes from the second dataset are colored pink, and the intersection nodes are colored lavender.
Screenshot of the visualisation:

<img src="https://github.com/user-attachments/assets/99746053-9b3b-48a2-842f-99c456fd6005" width="450">

`analyse("degree_distribution")` from `main.py` - calculates the best exponent to fit the power law curve with the data. Then plots on the loglog plot the fitted curve (straight line since it's a loglog plot), the original distribution for Biorefinery Network and Network of Organic Chemistry. It is done for both in and out-degree. `degree_distribution.py` allows for creation of more plots representing the degree distribution like bar and scatter plots. 

<img src="https://github.com/user-attachments/assets/34a468cc-3077-4515-8c8e-2c0463346400" width="450">

`analyse("important_molecules")` from `main.py` - finds the molecules with the highes betweenness centrality and the highest degree. Creates plots below.

<img src="https://github.com/user-attachments/assets/a36f0a06-1f90-4d6b-bf4e-7b50f5bdf051" width="450">
<img src="https://github.com/user-attachments/assets/0f580930-6d9f-46ca-ba4a-dd3d90de5c4b" width="450">

`analyse("degree_correlation")` from `main.py` - plots the in and out-degree correlation.

<img src="https://github.com/user-attachments/assets/48c42118-f951-43a2-a655-8fd8b13bf599" width="450">
<img src="https://github.com/user-attachments/assets/ed800b0c-165e-4b9c-a34f-176155f36895" width="450">

## Data
- `resources/reaction_data.tsv` is the main file that should contain all the reactions to analyse. To reproduce the results the data should be taken from REAXYS using keyword search. The _title, abstract and keywords_ parameter has to be set to **biorefiner**. A detailed instruction how to download the data can be found [here](#querying-reaxys-through-browser). The columns from the file used for the analysis are **Reactant**, **Product** and in some of the cases **Number of Reaction Steps** (capitalisation matters). The data **can** contain additional columns but they won't be used. 
- An example datafile can be found in `resources/mock_data.tsv`. It looks as follows: <br/> <img src="README imgs/img_3.png" width="500"> 
- The data in `resources/NOC data` is the data for recreating plots from the paper of [P. M. Jacob and A. Lapkin, 2018](https://pubs.rsc.org/en/content/articlelanding/2018/re/c7re00129k). It was obtained using an online tool for extracting data from plots - [WebPlotDigitizer](https://automeris.io/WebPlotDigitizer.html).

## Querying REAXYS through browser
1. Go to REAXYS [Query builder](https://www.reaxys.com/#/search/advanced) and make sure you are logged in.
2. From the menu on the right open `Topcis and Keywords` tab.
3. A drop-down menu will appear. Choose `Title, Abstracts & Keywords` from it by pressing it.
4. This added a filter represented by the box in the middle of the screen. In that box change `is` option to `contains` by pressing on the arrow pointing downwards in that box and selecting `contains` from the drop-down menu.
5. Next press on the text bar, right of the arrow, and type **biorefiner**. <br/> <img src="README imgs/img.png" width="500"> 
6. The press `Reactions >` button at the top of the page to search for these reactions. Wait until the reactions are found and you are forwarded to the next page.
7. Now press `Export` button at the top of the page. <br/> <img src="README imgs/img_1.png" width="500"> 
8. Change format to `Tab-delimited text` and untick `Additional options` to speed up the export. <br/> <img src="README imgs/img_2.png" width="250"> 
9. After the export is done press download. 
10. Rename the file to `reaction_data.tsv` and put it in the `resources` folder. 
11. The data is ready to use. 


## Testing
A small test suite exists for testing `utils` methods. Other parts of the code were tested manually and assessed qualitatively.

## Authors
This repository was created by [Jakub Kontak](https://github.com/jkontak13) during his Computer Science Honours Programme at TU Delft under the supervision of Jana M. Weber.
