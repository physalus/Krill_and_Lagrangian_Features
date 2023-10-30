# Krill and Lagrangian Features
 This repository contains the code and processed data to generate the analysis and figures in the manuscript *Submesoscale coupling of krill and whales revealed by aggregative Lagrangian coherent structures*. Authors: James A. Fahlbusch, David E. Cade, Elliott L. Hazen, Meredith L. Elliott, Benjamin T. Saenz, Jeremy A. Goldbogen and Jaime Jahncke

### How to cite

Following publication, please cite this compendium as:

> Fahlbusch, J. A., Cade, D. E., Hazen, E. L.,
> Elliott, M.L., Saenz, B.T., Goldbogen, J. A. and Jahncke, J. (2023). 
> *Compendium of R code and data for “Submesoscale coupling of krill and whales revealed by aggregative Lagrangian coherent structures”*. Accessed 28 Jan 2024.

A knitted R-Markdown with all of the Methods and Results (including the code to produce them) can be found at:
https://physalus.github.io/Krill_and_Lagrangian_Features/

All main text figures and supplemental figures are generated using the scripts found in the main directory and data found in the folders in this repository, however most are then brought into Adobe Illustrator for final figure composition. 

The repository is organized as follows:
* [:file\_folder: Main directory](https://github.com/physalus/Krill_and_Lagrangian_Features): all scripts 
  * Krill_FTLE_Analysis_PDF.Rmd contains all of the analyses for the the manuscript and produces a PDF file when knitted
  * Krill_FTLE_Analysis.Rmd is identical to Krill_FTLE_Analysis_PDF.Rmd but produces the HTML file found at the link above 
  * Appendix1_Data_Processing.R contains the raw data processing steps for the Krill, CTD and Cetacean datasets and including extracting FTLE (Requires the sub-directory "dataRaw")
  * Figure1_StudyAreaMap.R produces the panels for figure 1 in the manuscript
  * Figure2_4PanelFTLE.R produces the panels for figure 2 in the manuscript     
* [:file\_folder: dataProcessed](https://github.com/physalus/Krill_and_Lagrangian_Features/tree/main/dataProcessed): all processed data used in the analyses of this manuscript
  * Files include the outputs from Appendix1_Data_Processing.R and are used as inputs to Krill_FTLE_Analysis.Rmd
* [:file\_folder: Figure_files](https://github.com/physalus/Krill_and_Lagrangian_Features/tree/main/Figure_files): all additional data used to produce figures 1 and 2 of this manuscript
* [:file\_folder: Output](https://github.com/physalus/Krill_and_Lagrangian_Features/tree/main/Output): Plots that are exported during the knitting process and data precessing steps
* [:file\_folder: functions](https://github.com/physalus/Krill_and_Lagrangian_Features/tree/main/functions): Custom functions used this analysis

*Note: Due to space limitations, the following directories are incomplete on github. The files to populate these folders can be found at: https://purl.stanford.edu/kn066mv7984*
* [:file\_folder: dataRaw](https://github.com/physalus/Krill_and_Lagrangian_Features/tree/main/dataRaw): raw krill, CTD and cetacean data and surface current data (including restored and processed FTLE)
  * *These files are only necessary if you want to recreate the processed data. Processed data can be found in the folder "dataProcessed"*
  * FTLEMetadata.xlsx contains the settings used to calculate FTLE for each ACCESS cruise
* [:file\_folder: Models](https://github.com/physalus/Krill_and_Lagrangian_Features/tree/main/Models): contains the processed GLMM models in .rds format.  
  * *These files are the Generalized Linear Mixed Effects Models produced in the analysis and take a considerable amount of time and processing power to produce.*

### Licenses

**Code, Text and Figures :** [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/)

**Data :** [ODbL-1.0 Open Database License](https://opendatacommons.org/licenses/dbcl/1-0/)

