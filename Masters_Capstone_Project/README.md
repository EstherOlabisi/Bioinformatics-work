**Abstract** 

Mass Spectrometry (MS)-based proteomics is a powerful tool for studying biological samples on a large scale. Visualization is an important data analysis step for addressing quality concerns and interpreting MS proteomics data. This project was conducted within the Jennifer Geddes-McAlister Lab, where researchers primarily utilize Perseus for proteomics data analysis but encounter some limitations in generating multidimensional visualizations. To address the limitations, we introduce a web-based proteomics visualization software designed within the R-Shiny GUI platform, enabling users to generate custom visualizations.

The proteomics visualization tool provides a user-friendly interface for users to upload input files namely _protein groups_, _volcano threshold lines_, and _1D annotation enrichment_ files from Perseus. Distinct 1D annotation heatmaps are generated for molecular function, biological process, cellular component, and UniProt Keywords Gene Ontology (GO) terms. The features of the generated volcano plot include the ability to distinguish abundance patterns among compared groups, the display of significant GO terms, and protein ID selection. Additionally, s-curves can be generated for the proteomic states presented in the volcano plot. The PCA plot offers the option to construct data ellipses for proteomic states. Users can download all plots in a high-quality PNG format. The tool's interactive GUI makes it accessible to proteomic researchers irrespective of their experience with the R language while proficient users can leverage its code for customizations. Considerations for future research include making the code open-source and integrating this tool as a plugin within Perseus.


**Using the application**
* To run the Shiny R script, open the script in RStudio and click the "Run App" button at the top right corner of the editor.
* It is preferable to select the external (browser) option from the "Run App" drop-down list but you could also try the app in an RStudio Window.
* For test input files: 
  * prot_frm_perseus.txt is the "Protein Groups" file
  * 1d annot.txt is the "1D annotation" file
  * curve_matrix.txt is the "Volcano threshold lines" file