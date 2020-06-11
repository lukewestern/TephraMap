# TephraMap
These functions aim to use measurements of tephra thickness to create a continuous tephra field, which can be used to draw isopach thicknesses and estimate the total volume. This method is very heavily reliant on the package r-inla (Lindgren & Rue, 2015).
The main function is tephra.map(), which calls the other functions internally. 
The function tephra.map() takes a csv file containing sampled tephra thickness to estimate appropriate isopachs and total volume of tephra. The input data is a csv file with X and Y coordinates in: UTM units, as distance in kilometres from some origin, or in WGS84 latitude/longitude, which will be converted to UTM coordinates internally. If you are estimating for a super eruption (super==TRUE) the coordinates must be in latitude/longitude for the calculation. The third column contains the measured thickness of the tephra deposit, with the column header as the units of the measurement in either mm, cm or m. For example the first 3 columns may appear as:  
 X,Y,mm <br/>
 522100, 1378210, 15 <br/>
 521140, 1377890, 120 <br/>
 Estimates are done on a 2D plain except for super eruptions (super=TRUE) when it is calculated on the globe, i.e. curvature is taken into account.  
 Install it by doing: 
 install.packages(''INLA'', repos=c(getOption(''repos''), INLA=''https://inla.r-inla-download.org/R/stable''), dep=TRUE) 
 install.packages(''<path-to>/TephraMap'')
 References:
 Lindgren, F. and H. Rue, 2015. Bayesian spatial modelling with r-inla. Journal of Statistical Software, 63(19).