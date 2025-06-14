# Load necessary libraries for PCA and visualization
library(factoextra)     # For PCA visualization
library(FactoMineR)     # For PCA computation
library(gridExtra)      # For arranging multiple ggplots
library(ggbiplot)       # For enhanced biplots
library(corrplot)       # For correlation plots
library(grid)           # For grid text elements
library(writexl)        # To export results to Excel

# Subset relevant numeric columns (traits) from WHEAT dataset
data <- WHEAT[c(2:19)]

# Perform PCA using prcomp (standard PCA function in R)
pca <- prcomp(data, scale. = TRUE)  # Scale standardizes the data
summary(pca)                        # Summary: importance of components
get_eigenvalue(pca)                # Extract eigenvalues (variance explained)

# Perform PCA using FactoMineR with genotypes as row names
TazData <- PCA(data.frame(WHEAT[1:19], row.names = 1), graph = FALSE)
summary(TazData)                   # Detailed PCA results
eig.val <- get_eigenvalue(TazData) # Get eigenvalues from FactoMineR PCA

# Scree plot: shows % variance explained by each component
fviz_eig(TazData, addlabels = TRUE, ylim = c(0, 25)) +
  theme(text = element_text(size = 16, family = "Times New Roman"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 18, hjust = 0.5))

# Individual PCA plot: shows genotype representation and quality (cos2)
fviz_pca_ind(TazData, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) +
  theme(text = element_text(size = 16, family = "Times New Roman"))

# PCA biplot: combines individuals and variables
fviz_pca_biplot(TazData, repel = TRUE, 
                col.var = "red", col.ind = "black",
                labelsize = 6, pointsize = 3, arrowsize = 1.2) +
  theme(text = element_text(size = 16, family = "Times New Roman"))

# Cos2 plot: quality of representation of variables on PC1 and PC2
fviz_cos2(pca, choice = "var", axes = 1:2) +
  theme(text = element_text(size = 16, family = "Times New Roman"))

# Bar plot: contribution of each trait to PC1 and PC2
a <- fviz_contrib(pca, choice = "var", axes = 1) + 
  theme(text = element_text(size = 16, family = "Times New Roman"))
b <- fviz_contrib(pca, choice = "var", axes = 2) + 
  theme(text = element_text(size = 16, family = "Times New Roman"))
grid.arrange(a, b, ncol = 2, 
             top = textGrob("Contribution of the variables to the first two PCs", 
                            gp = gpar(fontsize = 16, family = "Times New Roman")))

# Correlation matrix of squared cosine values
corrplot(pca$rotation^2, is.corr = FALSE)  # Cos2 matrix (squared loadings)

# Contribution of variables to the first 7 principal components
fviz_cos2(pca, choice = "var", axes = 1:7) +
  theme(text = element_text(size = 16, family = "Times New Roman"))

fviz_contrib(pca, choice = "var", axes = 1:7) +
  theme(text = element_text(size = 18, family = "Times New Roman"))

# Top 5 contributing traits to PC1 and PC2
a <- fviz_contrib(pca, choice = "var", axes = 1, top = 5) +
  theme(text = element_text(size = 18, family = "Times New Roman"))
b <- fviz_contrib(pca, choice = "var", axes = 2, top = 5) +
  theme(text = element_text(size = 18, family = "Times New Roman"))
grid.arrange(a, b, ncol = 2,
             top = textGrob("Contribution of the variables to the first two PCs", 
                            gp = gpar(fontsize = 16, family = "Times New Roman")))

# Correlation heatmap of original variables
corrplot(cor(data), diag = FALSE, tl.pos = 'n')  # Visualize correlation without labels

# Another PCA run using FactoMineR for full result object
res.pca <- PCA(data)
res.pca$eig                       # Eigenvalues table
fviz_eig(res.pca)                # Scree plot

# Variable correlation plots for multiple axis combinations
fviz_pca_var(res.pca, axes = 1:2, repel = TRUE) +
  theme(text = element_text(size = 16, family = "Times New Roman"))
fviz_pca_var(res.pca, axes = 2:3, repel = TRUE) +
  theme(text = element_text(size = 16, family = "Times New Roman"))
fviz_pca_var(res.pca, axes = 3:4, repel = TRUE) +
  theme(text = element_text(size = 16, family = "Times New Roman"))

# Contribution plots for each axis (PC1 to PC5)
fviz_contrib(res.pca, choice = 'var')
fviz_contrib(res.pca, choice = 'var', axes = 2)
fviz_contrib(res.pca, choice = 'var', axes = 3)
fviz_contrib(res.pca, choice = 'var', axes = 4)
fviz_contrib(res.pca, choice = 'var', axes = 5)
fviz_contrib(res.pca, choice = 'var', axes = 1:5)

# Extract PCA result components from prcomp
eig.val <- get_eigenvalue(pca)
sdev <- as.data.frame(pca$sdev)           # Standard deviations of PCs
loadings <- as.data.frame(pca$rotation)   # Variable loadings (correlations)
scores <- as.data.frame(pca$x)            # Coordinates of individuals (samples)
center <- as.data.frame(t(pca$center))    # Variable means used for scaling
scale <- as.data.frame(t(pca$scale))      # Standard deviations used for scaling

# Store all extracted components into a list
pca_results <- list(
  StandardDeviations = sdev,
  Loadings = loadings,
  Scores = scores,
  Center = center,
  Scale = scale,
  Eigenvalues = eig.val
)

# Export PCA results to Excel
write_xlsx(pca_results, "PCA_Results.xlsx")

###Results
library(readxl)
> WHEAT <- read_excel("D:/khem data analysis/WHEAT.xlsx", 
+     sheet = "Quantative")
> View(WHEAT)
> library(factoextra)
Loading required package: ggplot2
Welcome! Want to learn more? See two factoextra-related books at https://goo.gl/ve3WBa
Warning messages:
1: package ‘factoextra’ was built under R version 4.3.3 
2: package ‘ggplot2’ was built under R version 4.3.3 
> library(FactoMineR)
Warning message:
package ‘FactoMineR’ was built under R version 4.3.3 
> library(gridExtra)
> library(ggbiplot)
Warning message:
package ‘ggbiplot’ was built under R version 4.3.3 
> library(corrplot)
corrplot 0.95 loaded
Warning message:
package ‘corrplot’ was built under R version 4.3.3 
> library(factoextra)
> library(grid)
> data <- WHEAT[c(2:19)]
> prcomp(data, scale. = TRUE)
Standard deviations (1, .., p=18):
 [1] 1.9777185 1.7639145 1.4988523 1.2795466 1.0998498 1.0041925 0.9467576
 [8] 0.9231898 0.7857069 0.7530864 0.6848779 0.5981610 0.5392674 0.4681611
[15] 0.4422064 0.4323082 0.3626920 0.3023916

Rotation (n x k) = (18 x 18):
             PC1          PC2         PC3         PC4          PC5
DTE   0.09681961 -0.248915171  0.20321573  0.18291811 -0.310866010
DTH   0.12584540 -0.472261847 -0.11757530 -0.02089494  0.227670605
DTF   0.11786668 -0.478033440 -0.10823993  0.07442274  0.260492218
DOM   0.15382854 -0.396641906 -0.01175646  0.12506024  0.015418484
PH   -0.29493475  0.090749471 -0.02647058  0.51611879 -0.002855896
SE   -0.32100654  0.007948727  0.12512940  0.45119013  0.166639957
SPL   0.30414096  0.078260779 -0.29915207  0.09559377  0.213019423
FLL   0.13775089 -0.278355180  0.20422376  0.28393203  0.002210773
FLW   0.15420953 -0.090093908  0.08076752  0.13558631 -0.596214629
NOT  -0.29329146 -0.006130158  0.10613950  0.32507416  0.069803980
SPPS  0.07545754  0.123449623 -0.44166714  0.22780137  0.148851566
SPSP  0.11781520  0.071170214 -0.38644921  0.04356170 -0.363938049
SPS   0.04275478  0.042129811 -0.52071076  0.30040555 -0.166998395
SL    0.38830250  0.049028996  0.08834541  0.15312397  0.215878702
SW    0.30476548  0.313058677  0.22228522  0.20678397  0.045037553
ST    0.26710421  0.178960909  0.28393480  0.17440302  0.114303027
THSW  0.39603450  0.244846434  0.07610327  0.09084347  0.089742048
TY   -0.18075639  0.121010030 -0.09219605 -0.11068834  0.312544160
             PC6         PC7         PC8          PC9         PC10
DTE  -0.46928028  0.03413701  0.31851741 -0.066785870 -0.264969113
DTH  -0.01015222 -0.06168867 -0.02676874  0.024031187  0.162253523
DTF   0.05126315 -0.06852656 -0.04548086  0.090297007  0.184258997
DOM  -0.36270717  0.11344718  0.11912222 -0.161359504 -0.154628867
PH    0.05325463  0.11854808  0.07589190 -0.042899261 -0.191238894
SE    0.03038205 -0.13375927 -0.09453175  0.193883977 -0.314708627
SPL   0.09068765  0.37416328 -0.11007372 -0.048557079 -0.367103688
FLL   0.31952218 -0.29080512 -0.36855450  0.272813546 -0.073589696
FLW   0.11311009  0.47299455 -0.11386904  0.440569701  0.267332957
NOT  -0.10521778  0.23889818 -0.23530316 -0.464909115  0.544125414
SPPS  0.15444705  0.04835073  0.46513892  0.227329589  0.209280346
SPSP -0.16655343 -0.36789168 -0.46278922 -0.221054484 -0.041288165
SPS  -0.14597583 -0.18390496  0.04965776  0.021733506  0.131841701
SL    0.08034421  0.31389703 -0.15313911 -0.261063138  0.017250604
SW   -0.05347546 -0.06602247 -0.01011333  0.041153041  0.009994436
ST   -0.19549966 -0.36145516  0.23831084  0.101654000  0.362607952
THSW -0.14132321 -0.03508035 -0.15258663  0.001709919 -0.052076759
TY   -0.60591750  0.18407058 -0.33759547  0.503103402  0.069215100
            PC11        PC12        PC13         PC14        PC15
DTE  -0.54061634 -0.13466782 -0.01809289 -0.200801407 -0.06070830
DTH   0.11776145 -0.40034255 -0.17735298 -0.248233928  0.23279292
DTF  -0.02166232 -0.25882582  0.10351209 -0.021097785 -0.16208938
DOM   0.44346384  0.45918996  0.03096427  0.353472103  0.19722344
PH    0.32078620 -0.35885942 -0.06279575 -0.043144610  0.21744866
SE   -0.02548608 -0.10539174  0.16873101  0.219696735 -0.17380393
SPL  -0.09167535  0.02434861 -0.60174877  0.048310172 -0.26248742
FLL  -0.23381775  0.42259926 -0.08327764 -0.061354424  0.05571099
FLW   0.14254921 -0.12478704 -0.01671474  0.120988174 -0.04192751
NOT  -0.21590977  0.15462400 -0.20342993 -0.026820748 -0.07140553
SPPS -0.37841976  0.11231962  0.08686916  0.248641626  0.31531345
SPSP -0.13310067 -0.22978324 -0.05149386  0.325331691  0.27925676
SPS   0.21983005  0.24306438  0.18639975 -0.448139673 -0.38006271
SL   -0.10901393 -0.06951377  0.56518582  0.075341388  0.07476128
SW    0.11672507  0.09179200 -0.14458578 -0.451584481  0.50185821
ST    0.12919574 -0.11569746 -0.29245374  0.349808261 -0.22578792
THSW  0.06035859 -0.17803891  0.20356784 -0.008600316 -0.27772092
TY   -0.08991665  0.04853142  0.03434898 -0.023452926  0.11482848
             PC16        PC17        PC18
DTE   0.036329812 -0.05619228 -0.02343687
DTH  -0.111570022 -0.12840123  0.55961150
DTF  -0.001036696  0.33656060 -0.63467819
DOM   0.024898619  0.13902105  0.04659765
PH    0.463881412 -0.21268007 -0.17730568
SE   -0.525630874  0.21035409  0.21510749
SPL  -0.106348685 -0.05289111 -0.05886201
FLL   0.307194402 -0.18813141  0.04528104
FLW  -0.121654683  0.05941029  0.04212071
NOT   0.017527234  0.17161550  0.09445082
SPPS  0.098590122  0.14989011  0.12299422
SPSP -0.078004461  0.01313261 -0.07927746
SPS  -0.110569802 -0.17792081  0.01913831
SL   -0.161370832 -0.44242891 -0.05358645
SW   -0.304798202  0.27766113 -0.19249886
ST   -0.065021938 -0.32655664 -0.07230473
THSW  0.463593550  0.47882984  0.34018789
TY    0.085042689 -0.15508068 -0.06868274
> pca <- prcomp(data, scale. = TRUE)
> summary(pca)
Importance of components:
                          PC1    PC2    PC3     PC4    PC5     PC6    PC7
Standard deviation     1.9777 1.7639 1.4989 1.27955 1.0998 1.00419 0.9468
Proportion of Variance 0.2173 0.1729 0.1248 0.09096 0.0672 0.05602 0.0498
Cumulative Proportion  0.2173 0.3901 0.5150 0.60592 0.6731 0.72915 0.7789
                           PC8    PC9    PC10    PC11    PC12    PC13
Standard deviation     0.92319 0.7857 0.75309 0.68488 0.59816 0.53927
Proportion of Variance 0.04735 0.0343 0.03151 0.02606 0.01988 0.01616
Cumulative Proportion  0.82629 0.8606 0.89210 0.91816 0.93803 0.95419
                          PC14    PC15    PC16    PC17    PC18
Standard deviation     0.46816 0.44221 0.43231 0.36269 0.30239
Proportion of Variance 0.01218 0.01086 0.01038 0.00731 0.00508
Cumulative Proportion  0.96637 0.97723 0.98761 0.99492 1.00000
> # PC values for traits
> get_eigenvalue(pca)
       eigenvalue variance.percent cumulative.variance.percent
Dim.1  3.91137057       21.7298365                    21.72984
Dim.2  3.11139421       17.2855234                    39.01536
Dim.3  2.24655821       12.4808789                    51.49624
Dim.4  1.63723952        9.0957751                    60.59201
Dim.5  1.20966966        6.7203870                    67.31240
Dim.6  1.00840264        5.6022369                    72.91464
Dim.7  0.89635004        4.9797225                    77.89436
Dim.8  0.85227938        4.7348855                    82.62925
Dim.9  0.61733531        3.4296406                    86.05889
Dim.10 0.56713908        3.1507727                    89.20966
Dim.11 0.46905769        2.6058760                    91.81554
Dim.12 0.35779654        1.9877585                    93.80329
Dim.13 0.29080936        1.6156076                    95.41890
Dim.14 0.21917479        1.2176377                    96.63654
Dim.15 0.19554648        1.0863693                    97.72291
Dim.16 0.18689036        1.0382798                    98.76119
Dim.17 0.13154550        0.7308083                    99.49200
Dim.18 0.09144065        0.5080036                   100.00000
> TazData <- PCA(data.frame(WHEAT[1:19], row.names = 1), graph = FALSE)
> summary(TazData)

Call:
PCA(X = data.frame(WHEAT[1:19], row.names = 1), graph = FALSE) 


Eigenvalues
                       Dim.1   Dim.2   Dim.3   Dim.4   Dim.5   Dim.6
Variance               3.911   3.111   2.247   1.637   1.210   1.008
% of var.             21.730  17.286  12.481   9.096   6.720   5.602
Cumulative % of var.  21.730  39.015  51.496  60.592  67.312  72.915
                       Dim.7   Dim.8   Dim.9  Dim.10  Dim.11  Dim.12
Variance               0.896   0.852   0.617   0.567   0.469   0.358
% of var.              4.980   4.735   3.430   3.151   2.606   1.988
Cumulative % of var.  77.894  82.629  86.059  89.210  91.816  93.803
                      Dim.13  Dim.14  Dim.15  Dim.16  Dim.17  Dim.18
Variance               0.291   0.219   0.196   0.187   0.132   0.091
% of var.              1.616   1.218   1.086   1.038   0.731   0.508
Cumulative % of var.  95.419  96.637  97.723  98.761  99.492 100.000

Individuals (the 10 first)
               Dist    Dim.1    ctr   cos2    Dim.2    ctr   cos2    Dim.3
NGRCO 5102 |  4.034 | -2.381  1.250  0.348 | -1.033  0.296  0.066 | -1.294
NGRCO5103  |  9.397 | -0.895  0.176  0.009 |  7.484 15.517  0.634 | -1.176
NGRCO5104  |  4.858 | -3.091  2.106  0.405 | -1.987  1.093  0.167 |  1.815
NGRCO5188  |  4.683 | -1.465  0.473  0.098 | -2.072  1.189  0.196 |  0.727
NGRCO5189  |  3.488 | -2.144  1.013  0.378 | -1.900  1.000  0.297 |  0.553
NGRCO5190  |  3.559 | -2.161  1.029  0.369 | -1.979  1.085  0.309 |  0.366
NGRCO5191  |  3.762 | -0.975  0.209  0.067 | -0.286  0.023  0.006 |  1.544
NGRCO5192  |  3.493 | -0.862  0.164  0.061 | -2.374  1.561  0.462 |  0.457
NGRCO5193  |  2.257 |  0.757  0.126  0.112 | -1.125  0.351  0.248 |  0.862
NGRCO5194  |  3.347 | -0.223  0.011  0.004 | -1.575  0.687  0.221 | -1.809
              ctr   cos2  
NGRCO 5102  0.642  0.103 |
NGRCO5103   0.530  0.016 |
NGRCO5104   1.265  0.140 |
NGRCO5188   0.203  0.024 |
NGRCO5189   0.117  0.025 |
NGRCO5190   0.051  0.011 |
NGRCO5191   0.914  0.168 |
NGRCO5192   0.080  0.017 |
NGRCO5193   0.285  0.146 |
NGRCO5194   1.256  0.292 |

Variables (the 10 first)
              Dim.1    ctr   cos2    Dim.2    ctr   cos2    Dim.3    ctr
DTE        |  0.191  0.937  0.037 |  0.439  6.196  0.193 | -0.305  4.130
DTH        |  0.249  1.584  0.062 |  0.833 22.303  0.694 |  0.176  1.382
DTF        |  0.233  1.389  0.054 |  0.843 22.852  0.711 |  0.162  1.172
DOM        |  0.304  2.366  0.093 |  0.700 15.732  0.489 |  0.018  0.014
PH         | -0.583  8.699  0.340 | -0.160  0.824  0.026 |  0.040  0.070
SE         | -0.635 10.305  0.403 | -0.014  0.006  0.000 | -0.188  1.566
SPL        |  0.602  9.250  0.362 | -0.138  0.612  0.019 |  0.448  8.949
FLL        |  0.272  1.898  0.074 |  0.491  7.748  0.241 | -0.306  4.171
FLW        |  0.305  2.378  0.093 |  0.159  0.812  0.025 | -0.121  0.652
NOT        | -0.580  8.602  0.336 |  0.011  0.004  0.000 | -0.159  1.127
             cos2  
DTE         0.093 |
DTH         0.031 |
DTF         0.026 |
DOM         0.000 |
PH          0.002 |
SE          0.035 |
SPL         0.201 |
FLL         0.094 |
FLW         0.015 |
NOT         0.025 |
> library("factoextra")
> eig.val <- get_eigenvalue(TazData)
> # Scree Plot
> fviz_eig(TazData, addlabels = TRUE, ylim = c(0, 25)) +
+   theme(
+     text = element_text(size = 16, family = "Times New Roman"),  # General text size
+     axis.title = element_text(size = 14, family = "Times New Roman"),  # Axis title size
+     axis.text = element_text(size = 12, family = "Times New Roman"),   # Axis tick label size
+     plot.title = element_text(size = 18, family = "Times New Roman", hjust = 0.5, face = "plain")  # Centered and non-bold plot title
+   )
Warning messages:
1: In grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family not found in Windows font database
2: In grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family not found in Windows font database
3: In grid.Call(C_stringMetric, as.graphicsAnnot(x$label)) :
  font family not found in Windows font database
4: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
5: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
6: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
7: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
8: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
9: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
10: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
> # PCA individual
There were 14 warnings (use warnings() to see them)
> fviz_pca_ind(
+   TazData, col.ind = "cos2",
+   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
+   repel = TRUE
+ ) + 
+   theme(
+     text = element_text(size = 16, family = "Times New Roman"),
+     axis.title = element_text(size = 14, family = "Times New Roman"),
+     axis.text = element_text(size = 12, family = "Times New Roman")
+   )
There were 14 warnings (use warnings() to see them)
> # Biplot
> fviz_pca_biplot(TazData, repel = TRUE, 
+                 col.var = "red", col.ind = "black"
+ ) + 
+   theme(
+     text = element_text(size = 16, family = "Times New Roman"),
+     axis.title = element_text(size = 14, family = "Times New Roman"),
+     axis.text = element_text(size = 12, family = "Times New Roman")
+   )
There were 26 warnings (use warnings() to see them)
> fviz_pca_biplot(TazData, repel = TRUE, 
+                 col.var = "red",
+                 col.ind = "black",
+                 labelsize = 6,
+                 pointsize = 3,
+                 arrowsize = 1.2
+ ) + 
+   theme(
+     text = element_text(size = 16, family = "Times New Roman"),
+     axis.title = element_text(size = 14, family = "Times New Roman"),
+     axis.text = element_text(size = 12, family = "Times New Roman")
+   )
Warning messages:
1: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
2: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
3: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
4: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
5: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
6: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
7: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
8: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
> # To see the most contributing variables for both dimensions
> fviz_cos2(pca, choice = "var", axes = 1:2) + 
+   theme(text = element_text(size = 16, family = "Times New Roman"))
There were 42 warnings (use warnings() to see them)
> # Bar plot of variable contributions
> # Contributions of variables to PC1
> a <- fviz_contrib(pca, choice = "var", axes = 1) + 
+   theme(text = element_text(size = 16, family = "Times New Roman"))
> # Contributions of variables to PC2
There were 50 or more warnings (use warnings() to see the first 50)
> b <- fviz_contrib(pca, choice = "var", axes = 2) + 
+   theme(text = element_text(size = 16, family = "Times New Roman"))
> grid.arrange(a, b, ncol = 2, top = textGrob("Contribution of the variables to the first two PCs", 
+                                             gp = gpar(fontsize = 16, family = "Times New Roman")))
There were 50 or more warnings (use warnings() to see the first 50)
> # Most contributing variables for each dimension
There were 50 or more warnings (use warnings() to see the first 50)
> corrplot(var$cos2, is.corr = FALSE)
Error in var$cos2 : object of type 'closure' is not subsettable
> # Contribution of variables for dimensions 1 to 7
> fviz_cos2(pca, choice = "var", axes = 1:7) + 
+   theme(text = element_text(size = 16, family = "Times New Roman"))
There were 42 warnings (use warnings() to see them)
> # Total contribution on PC1 and PC2
There were 50 or more warnings (use warnings() to see the first 50)
> fviz_contrib(pca, choice = "var", axes = 1:7) + 
+   theme(text = element_text(size = 18, family = "Times New Roman"))
There were 43 warnings (use warnings() to see them)
> # Bar plots of variable contributions
There were 50 or more warnings (use warnings() to see the first 50)
> a <- fviz_contrib(pca, choice = "var", axes = 1, top = 5) + 
+   theme(text = element_text(size = 18, family = "Times New Roman"))
> b <- fviz_contrib(pca, choice = "var", axes = 2, top = 5) + 
+   theme(text = element_text(size = 18, family = "Times New Roman"))
> grid.arrange(a, b, ncol = 2, top = textGrob("Contribution of the variables to the first two PCs", 
+                                             gp = gpar(fontsize = 16, family = "Times New Roman")))
There were 32 warnings (use warnings() to see them)
> # Redundancy
There were 50 or more warnings (use warnings() to see the first 50)
> library(corrplot)
> corrplot(cor(data), diag = FALSE, tl.pos = 'n')
> res.pca <- PCA(data)
Warning message:
ggrepel: 8 unlabeled data points (too many overlaps). Consider increasing max.overlaps 
> # Inspect eigenvalues
> res.pca$eig
        eigenvalue percentage of variance cumulative percentage of variance
comp 1  3.91137057             21.7298365                          21.72984
comp 2  3.11139421             17.2855234                          39.01536
comp 3  2.24655821             12.4808789                          51.49624
comp 4  1.63723952              9.0957751                          60.59201
comp 5  1.20966966              6.7203870                          67.31240
comp 6  1.00840264              5.6022369                          72.91464
comp 7  0.89635004              4.9797225                          77.89436
comp 8  0.85227938              4.7348855                          82.62925
comp 9  0.61733531              3.4296406                          86.05889
comp 10 0.56713908              3.1507727                          89.20966
comp 11 0.46905769              2.6058760                          91.81554
comp 12 0.35779654              1.9877585                          93.80329
comp 13 0.29080936              1.6156076                          95.41890
comp 14 0.21917479              1.2176377                          96.63654
comp 15 0.19554648              1.0863693                          97.72291
comp 16 0.18689036              1.0382798                          98.76119
comp 17 0.13154550              0.7308083                          99.49200
comp 18 0.09144065              0.5080036                         100.00000
> fviz_eig(res.pca)
> # Correlations
> fviz_pca_var(res.pca, axes = 1:2, repel = TRUE) + 
+   theme(text = element_text(size = 16, family = "Times New Roman"))
Warning messages:
1: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
2: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
3: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
4: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
5: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
6: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
7: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
8: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
> fviz_pca_var(res.pca, axes = 2:3, repel = TRUE) + 
There were 14 warnings (use warnings() to see them)
+   theme(text = element_text(size = 16, family = "Times New Roman"))
Warning messages:
1: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
2: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
3: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
4: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
5: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
6: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
7: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
8: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
> fviz_pca_var(res.pca, axes = 3:4, repel = TRUE) + 
There were 14 warnings (use warnings() to see them)
+   theme(text = element_text(size = 16, family = "Times New Roman"))
Warning messages:
1: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
2: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
3: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
4: In grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
5: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
6: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
7: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
8: In grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  :
  font family not found in Windows font database
> corrplot(res.pca$var$coord)
There were 14 warnings (use warnings() to see them)
> fviz_contrib(res.pca, choice = 'var')
> fviz_contrib(res.pca, choice = 'var', axes = 2)
> fviz_contrib(res.pca, choice = 'var', axes = 3)
> fviz_contrib(res.pca, choice = 'var', axes = 4)
> fviz_contrib(res.pca, choice = 'var', axes = 5)
> fviz_contrib(res.pca, choice = 'var', axes = 1:5)
> # Create a list of data frames
> # Get eigenvalues
> eig.val <- get_eigenvalue(pca)
> print(eig.val)
       eigenvalue variance.percent cumulative.variance.percent
Dim.1  3.91137057       21.7298365                    21.72984
Dim.2  3.11139421       17.2855234                    39.01536
Dim.3  2.24655821       12.4808789                    51.49624
Dim.4  1.63723952        9.0957751                    60.59201
Dim.5  1.20966966        6.7203870                    67.31240
Dim.6  1.00840264        5.6022369                    72.91464
Dim.7  0.89635004        4.9797225                    77.89436
Dim.8  0.85227938        4.7348855                    82.62925
Dim.9  0.61733531        3.4296406                    86.05889
Dim.10 0.56713908        3.1507727                    89.20966
Dim.11 0.46905769        2.6058760                    91.81554
Dim.12 0.35779654        1.9877585                    93.80329
Dim.13 0.29080936        1.6156076                    95.41890
Dim.14 0.21917479        1.2176377                    96.63654
Dim.15 0.19554648        1.0863693                    97.72291
Dim.16 0.18689036        1.0382798                    98.76119
Dim.17 0.13154550        0.7308083                    99.49200
Dim.18 0.09144065        0.5080036                   100.00000
> # Extract PCA components
> sdev <- as.data.frame(pca$sdev)  # Standard deviations
> loadings <- as.data.frame(pca$rotation)  # Loadings
> scores <- as.data.frame(pca$x)  # Scores
> center <- as.data.frame(t(pca$center))  # Center
> scale <- as.data.frame(t(pca$scale))  # Scale
> # Create a list of data frames
> pca_results <- list(
+   StandardDeviations = sdev,
+   Loadings = loadings,
+   Scores = scores,
+   Center = center,
+   Scale = scale,
+   Eigenvalues = eig.val
+ )
> # Write to Excel file
library(writexl)
write_xlsx(pca_results, "PCA_Results.xlsx")
