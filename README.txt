#This file provides useful information on the operation steps and importing data. For more detailed information, please refer to the code in the corresponding script. 

##Step 1.Univariate cox analysis
##Import Data===================================================================

##File 1:Gene expression file===================================================
##Each row represents one gene and Each column represents one sample============

#TCGA.KO.8406.01 TCGA.KO.8417.01 TCGA.KO.8415.01 TCGA.KO.8413.01 TCGA.KM.8476.01 TCGA.KM.8477.01
#KIF15          7.469677        6.800114        5.980725        5.437914        5.241341        5.771213
#FOXM1          9.972000        9.240809        9.800237        8.050680        8.279642        7.729475
#CENPI          6.707639        5.522011        5.437832        4.535103        2.337723        4.349263

##File 2:Clinical file==========================================================
##This file must include the three columns of Sample, OS and OS.time============
            #Sample OS OS.time Age Gender  M  N T Stage
#TCGA.KO.8406.01 TCGA.KO.8406.01  0    3145  17      0  0  0 1     1
#TCGA.KO.8417.01 TCGA.KO.8417.01  0    2815  30      0  0  0 1     1
#TCGA.KO.8415.01 TCGA.KO.8415.01  0    2939  44      0  0  0 1     1

##File3£ºsignificant genes======================================================

#A1CF
#A2M
#A4GNT


##Step 2.stepwise regression analysis
##Import Data===================================================================

##File 1:Gene expression file===================================================
##Each row represents one gene and Each column represents one sample============

#TCGA.KO.8406.01 TCGA.KO.8417.01 TCGA.KO.8415.01 TCGA.KO.8413.01 TCGA.KM.8476.01 TCGA.KM.8477.01
#KIF15          7.469677        6.800114        5.980725        5.437914        5.241341        5.771213
#FOXM1          9.972000        9.240809        9.800237        8.050680        8.279642        7.729475
#CENPI          6.707639        5.522011        5.437832        4.535103        2.337723        4.349263

##File 2:Clinical file==========================================================
##This file must include the three columns of Sample, OS and OS.time============

#Sample OS OS.time Age Gender  M  N T Stage
#TCGA.KO.8406.01 TCGA.KO.8406.01  0    3145  17      0  0  0 1     1
#TCGA.KO.8417.01 TCGA.KO.8417.01  0    2815  30      0  0  0 1     1
#TCGA.KO.8415.01 TCGA.KO.8415.01  0    2939  44      0  0  0 1     1

##File3£ºsignificant genes======================================================

#gene      P_value          HR  Low.95.CI  High.95.CI          fdr
#KIF15 9.869823e-06  2.30078562 1.59004475   3.3292236 0.0005058284
#FOXM1 9.483495e-05  2.33063261 1.52391529   3.5644031 0.0012960777
#CENPI 8.293118e-05  1.95951886 1.40172432   2.7392791 0.0011724753


##Step 3.K-times LASSO iterative regression analysis 
##Import Data===================================================================

##File 1:Gene expression file===================================================
##Each row represents one gene and Each column represents one sample============

#TCGA.KO.8406.01 TCGA.KO.8417.01 TCGA.KO.8415.01 TCGA.KO.8413.01 TCGA.KM.8476.01 TCGA.KM.8477.01
#KIF15          7.469677        6.800114        5.980725        5.437914        5.241341        5.771213
#FOXM1          9.972000        9.240809        9.800237        8.050680        8.279642        7.729475
#CENPI          6.707639        5.522011        5.437832        4.535103        2.337723        4.349263

##File 2:Clinical file==========================================================
##This file must include the three columns of Sample, OS and OS.time============

#Sample OS OS.time Age Gender  M  N T Stage
#TCGA.KO.8406.01 TCGA.KO.8406.01  0    3145  17      0  0  0 1     1
#TCGA.KO.8417.01 TCGA.KO.8417.01  0    2815  30      0  0  0 1     1
#TCGA.KO.8415.01 TCGA.KO.8415.01  0    2939  44      0  0  0 1     1

##File3£ºsignificant genes======================================================

#gene      P_value          HR  Low.95.CI  High.95.CI          
#KIF15 9.869823e-06  2.30078562 1.59004475   3.3292236 
#FOXM1 9.483495e-05  2.33063261 1.52391529   3.5644031 
#CENPI 8.293118e-05  1.95951886 1.40172432   2.7392791 

