ExpandingCell
=============

This code builds an application for multi-scale model validation. The method was developed by Andy and Nick (see Malleson et al. 2009) based on a method by Costanza (1989). The method works by placing a square grid over the two point data sets and counting the number of points within each square. After this, it is possible to use traditional goodness-of-fit statistics such as the Standardised Root Mean Square Error (SRMSE) and R2 to quantitatively estimate the similarity of the two datasets (and hence the validity of the model). Also, the differences in counts in each cell can be mapped to show where the two datasets are similar and where they are different. A nice feature of this type of approach is that it is possible to change the size of each square (and the number of cells in the grid), thereby comparing the datasets at a various spatial resolutions. 

References
Costanza, R. (1989) Model goodness of fit: A multiple resolution procedure. Ecological Modelling 47, 199-215.
Malleson, N., Heppenstall, A., See, L. and Evans, A. (2009) Evaluating an Agent-Based Model of Burglary. Working Paper 10/1, School of Geography, University of Leeds. Available at http://www.geog.leeds.ac.uk/fileadmin/downloads/school/research/wpapers/10_1.pdf
