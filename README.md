# top_down_islets_cytokine
This repository is for hosting the scripts and files to reproduce the analyses described in the listed manuscript.

A description of each item in the repository: 
- **limmaFit.R** ### Implemented in "1_DEA.R".
- **limmaDEA.R** ### Implemented in "1_DEA.R".
- **0a_loading_intensity_data.R** Converts the output of TopPIC searches to an MSnSet object that stores label-free quantification data. Final data is on a relative log2-scale.
- **0b_loading_spectralcount_data.R** Converts the output of TopPIC searches to an MSnSet object that stores spectral count data. 
- **1_DEA.R** Takes output of script "0a" and performs differential abundance analysis using the functions "limmaFit.R" and "limmaDEA.R". 
- **Figure1_overview_plots.R** Takes output of script "0a" and generates Figure 1, including quantifying the number of distinct and average proteoforms/genes and the completeness of those identifications. 
- **Figure2_ins_gcg_rect_lines.R** Takes output of script "0b" and generates Figure 2: Rectangle plots of most frequently observed proteoforms on insulin and glucagon.
- **Figure3_volcano_plots_rawp.R** Takes output of script "1_DEA" and generates Figure 3: Volcano plots of all and curated proteoforms.
- **Figure4_CXCL.R** Takes output of script "0a" and generates Figure 4: Heatmap and table summary of observed CXCL proteoforms. Performs hypergeometric probability distribution test to determine which unique observations are statistically significant.
- 



