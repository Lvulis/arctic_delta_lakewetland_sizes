# arctic_delta_lakewetland_sizes
Scripts implementing the waterbody extraction and basic statistical distribution fitting described in Vulis et al. 2021 (under review in GRL)
A quick overview of how to use these scripts: First download Global Surface Water Masks from Google Earth Engine and outline the subaerial delta into a shapefile in a locally relevant UTM.
Then, call `occurrence_and_watermaskplots` to generate an occurrence mask. 
Then, use waterbody_extraction to extract waterbodies.
Then, use plot_together_multipleyr to make multi-delta/year plots
To compare with the PeRL dataset (see final figure in supplementary material), use the PeRL_analysis script after running plot_together_multipleyr.
