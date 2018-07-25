## This script will open the results from CompareTree algorithm ran in the cluster
## These were 600 output files, each has the results of 10 simulations (repetitions)
## of a specific combination of 2 parameters : Mutation rate (from 0.01 to 0.3) and
## Number of targets (from 10 to 200). 
## This script will take the mean of the 10 repetitions and plot a heatmap with the 
## accuracy of the Mutation Rate - Number of targets parameter space
## -------------------------------------------------------------------------------------
require(graphics); require(grDevices)
require(gplots);
require(ggplot2);
require(reshape2);
require(RColorBrewer);
require(lattice)

## open accuracy values

AccMatrix <- read.csv("./Acc_matrix_Parameter_Space.csv",check.names = F)
rownames(AccMatrix) = AccMatrix[,1]
AccMatrix = AccMatrix[,-1]

## PLOT ------------------------------------

## reshape the data
df.m <- melt(as.data.frame(AccMatrix),factorsAsStrings = T)
colnames(df.m) = c("Mu", "Accuracy")
head(df.m)
df.m$Targets = rep(seq(from = 300,to = 10,by = -10),30)

df.m$Mu = as.numeric(levels(df.m$Mu))[df.m$Mu]

# Define palette
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

## new plot
tiff(filename = "./Parameter_Space_csCompareTree_heatmap_10-300.tiff",width = 600,height = 400,res = 150)

hm1 <- ggplot(df.m,
              aes(x = Mu, y = Targets, fill = Accuracy))
hm1 <- hm1 + geom_tile()
hm1 <- hm1 + ggtitle("Accuracy (R-F in %)",subtitle = "4 states; 16 cell divisions") 
hm1 <- hm1 + scale_fill_gradientn(colours = myPalette(100))
hm1 <- hm1 + scale_x_continuous(expand = c(0,0), breaks = seq(0.05,0.25,0.05))
hm1 <- hm1 + scale_y_continuous(expand = c(0, 0))
hm1 <- hm1 + theme_bw()

hm1

dev.off()

##### PLOT WITH ELEVATION PLOT AND LOESS


# We make use of the loess function to fit a local polynomial trend surface 
# (using weighted least squares) to approximate the elevation across the whole region.
# The function call for a local quadratic surface is shown below:
par(mfrow=c(2,1))

# create data.frame
elevation.df = data.frame(x = 1000 *df.m$Mu,
                          y =  df.m$Targets, 
                          z = 100 * df.m$Accuracy)


elevation.loess = loess(z ~ x*y, data = elevation.df, degree = 2, span = 0.25)

# The expand.grid function creates an array of all combinations of the x and y values 
# that we specify in a list. We choose a range every foot from 10 to 300 feet to create a fine grid:

elevation.fit = expand.grid(list(x = seq(10, 300,1), y = seq(10, 300, 1)))

# The predict function is then used to estimate the surface height at all of these combinations
# covering our grid region. This is saved as an object z which will be used by the base graphics function:

z = predict(elevation.loess, newdata = elevation.fit)

#Lattice Graphics

# The lattice pckg xpect the data in a different format so we make use of the as.numeric function 
# to convert from a table of heights to a single column and append to the object we create based 
# on all combinations of x and y coordinates:

elevation.fit$Height = as.numeric(z)

elevation.fit$Height[elevation.fit$Height<0]=0

# function levelplot for this type of graphical dispaly. We use the data stored in the object elevation.fit
# to create the graph with lattice graphics.

# tiff(filename = "./Parameter_Space_csCompareTree_heatmap.tiff",width = 900,height = 700,res = 200)
# 
# levelplot(Height ~ x*y, data = elevation.fit,
#           xlab = "Mutation Rate (x10^3)", ylab = "Targets",
#           main = "Reconstruction Accuracy",
#           col.regions = colorRampPalette(rev(brewer.pal(11, "Spectral")))(20),
# #          col.regions = colorRampPalette(terrain.colors(10))(100),
# #          col.regions = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(20),
#           at = seq(0,100,5) )
# 
# dev.off()
# ##


tiff(filename = "./Parameter_Space_csCompareTree_LOESS_heatmap_10-300.tiff",width = 1100,height = 1000,res = 300)

contourplot(Height ~ x*y, data = elevation.fit, add =T,
            xlab = list("Mutation Rate",cex=0.8),
            ylab = list("Targets",cex=0.8),
            main = list(label="Accuracy (R-F in %)\n 4 states; 16 cell div (~65k cells)",cex=0.85),
            colorkey = F,
            #            labels = TRUE,
            #            contour = TRUE,
            region = TRUE, alpha.regions=0.7,
            at = c(0,10,25,50,75,90,95,99,100,101),
            #at = c(seq(0,100,20),101),
            labels = list(cex=0.7),
            col.regions = colorRampPalette(rev(brewer.pal(10, "Spectral")))(100),
            lwd=1.5,lty=3,cex=0.8,
            scales = list(x = list(at = c(10,50,100,150,200,250,300),
                                   labels = c(0.01,"" ,0.1,"",0.2,"",0.3)), 
                          y = list(at = c(10,50,100,150,200,250,300),
                                   labels = c(10,"" ,100,"",200,"",300))
                          )
            )

dev.off()
