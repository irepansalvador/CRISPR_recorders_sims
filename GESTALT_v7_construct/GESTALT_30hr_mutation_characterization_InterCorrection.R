# This script analyse the output files from the perl script:
# "GESTALT_30hr_mutation_characterization_InterCorrection.pl"
# which are lists of mutations and the frequencies for each target
# of the v7 construct of GESTALT
#
# The ultimate objective of this analysis is to obtain experimental
# parameters to perform simulations and quantify the accuracy of the GESTALT
# project.
# (DATA from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81713)
#
###################	
sample = 1:6

ALL_MUTS = matrix(data = 0,nrow = 10,ncol = 6)
ALL_MUTS_REL = matrix(data = 0,nrow = 10,ncol = 6)
ALL_INTER = matrix(data = 0,nrow = 10,ncol = 6)
ALL_INTRA = matrix(data = 0,nrow = 10,ncol = 6)
ALL_INTRA_sum1 = matrix(data = 0,nrow = 10,ncol = 6)
ALL_SAT = matrix(data = 0,nrow = 10,ncol = 6)
ALL_NORM = matrix(data = 0,nrow = 10,ncol = 6)
ALL_MUT_FREQ = matrix(data = 0,nrow = 6,ncol = 60)


for (s in sample)	
{	
     # s = 6
    print(sample[s])	
    my.file1= paste("30hr_",sample[s],"_1x.s_INTER_corrected_results.txt",sep = "")	
    my.file2= paste("30hr_",sample[s],"_1x.s_INTRA_corrected_results.txt",sep = "")
    my.file3= paste("30hr_",sample[s],"_1x.s_OVERALL_corrected_results.txt",sep = "")
    
    my.out= paste("30hr_",sample[s],"_1x.s_corrected_results.png",sep = "")	
    mydata <- read.table(file = my.file1,header = T,	
                         colClasses = c("character","numeric","numeric","numeric",
                                        "numeric","character"))	
    mydata$Type = factor(mydata$Type, 
                         levels = c("1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10",
                                    "1-3","2-4","3-5","4-6","5-7","6-8","7-9","8-10",
                                    "1-4","2-5","3-6","4-7","5-8","6-9","7-10",
                                    "1-5","2-6","3-7","4-8","5-9","6-10",
                                    "1-6","2-7","3-8","4-9","5-10",
                                    "1-7","2-8","3-9","4-10",
                                    "1-8","2-9","3-10",
                                    "1-9","2-10",
                                    "1-10"))
    head(mydata)	
    
    
    mydata2 <- read.table(file = my.file2,header = T,	
                          colClasses = c("character","numeric","numeric"))
    #levels(mydata2$Target) = as.factor(1:10)
    mydata2$Target
    mydata2$Target = factor(mydata2$Target,levels=c(1:10))
    head(mydata2)	
    #mydata2$Target
    
    mydata3 = read.table(file = my.file3,header = T, sep = "\t",
                         colClasses = c("character","numeric","numeric","numeric"))
    c <- colnames(mydata3)
    c[2] <- "Unmutated"
    colnames(mydata3) <- c

    head(mydata3)
    ################# Analyse Inter-target mutations ##################	
    Inter.m <- matrix(data = 0,nrow = 10, ncol = 4,	
                      dimnames = list(1:10,c("start","end","both","size")))
    for (i in 1:length(mydata[,1]))
    {	
        a = mydata$Start[i]	
        b = mydata$End[i]	
        sz = mydata$Size[i]	
        Inter.m[a,1] = Inter.m[a,1] + mydata$Count[i];	
        Inter.m[b,2] = Inter.m[b,2] + mydata$Count[i];	
        Inter.m[sz,4] = Inter.m[sz,4] + mydata$Count[i];	
        # print (paste("line ",i," start ",a, " times ", mydata$Count[i],sep = ""))	
    }	
    Inter.m	
    Inter.m[,3] = Inter.m[,1] + Inter.m[,2]	
    Inter.m[Inter.m[,1] == 0,1] = NA	
    Inter.m[Inter.m[,2] == 0,2] = NA	
    Inter.m[Inter.m[,4] == 0,4] = NA	
    #Inter.m = Inter.m/sum(Inter.m[,2],na.rm = T)	
    #### TYPE OF INTERTARGETS	
    Inter_types = matrix(data = 0,nrow = length(levels(mydata$Type)), ncol = 2)	
    Inter_types = as.data.frame(Inter_types)	
    colnames(Inter_types) = c("Type", "Count")	
    Inter_types$Type = levels(mydata$Type)	
    l = levels(mydata$Type)	
    t = length(levels(mydata$Type))	
    for (i in 1:length(mydata[,1]))
    {	
        t =  mydata$Type[i]	
        a <- which(l==t)	
        Inter_types$Count[a] =  Inter_types$Count[a] + mydata$Count[i]	
    }	
    ################# Analyse Intra-target mutations ##################	
    Intra.m <- matrix(data = 0,nrow = 10, ncol = 2,	
                      dimnames = list(1:10,c("count","type")))
    for (i in 1:length(mydata2[,2]))
    {	
        a = mydata2$Count[i]	
        t = as.numeric(mydata2$Target[i])	
        Intra.m[t,1] = Intra.m[t,1] + a	
        Intra.m[t,2] = Intra.m[t,2] + 1	
    }	

    ###### Combine inter and intra	
    all = Intra.m[,1] + Inter.m[,3]	
    mydata3$Mutated <- all
    mydata3$Normalized <- mydata3$Mutated / (mydata3$Unmutated + mydata3$Mutated)
    N.muts = sum(all)	
    ALL_MUTS[,s]  = all/N.muts
    ALL_INTER[,s] = Inter.m[,3]/N.muts
    ALL_INTRA_sum1[,s] = Intra.m[,1]/sum(Intra.m[,1])
    ALL_INTRA[,s] = Intra.m[,1]/N.muts
    ALL_SAT[,s] = mydata3$Proportion
    ALL_NORM[,s] = mydata3$Normalized

    ALL_MUTS_REL[,s] = all/ (mydata3$Reads)
    
    ##################### PLOT ###################
    png(filename = my.out,width = 900,height = 500)	
    par(mfrow=c(2,3))	
    
    ## plot 1	
    plot(0,0, xlim = c(1,10), ylim = c(0,0.26),xlab = "Target", ylab = "Mut Proportion per target",	
         main = paste("Mutations in 30hr sample # ", s, sep = ""),cex.main=2)	
    points(Intra.m[,1]/N.muts, pch = 6, lwd = 2, col = "blue")	
    lines(Intra.m[,1]/N.muts, lty = 6,lwd= 2, col = "blue")	
    points(Inter.m[,3]/N.muts, pch = 2, lwd = 2, col = "orange")	
    lines(Inter.m[,3]/N.muts, lty = 6,lwd= 2, col = "orange")	
    points(all/N.muts, pch = 1, lwd = 2, col = "black")	
    lines(all/N.muts, lty = 1,lwd= 2, col = "black")	
    legend(x = 2, y = 0.26,legend = c("Intra", "Inter", "All"),col = c("blue","orange","black"),	
           pch = c(6,2,1), lwd = c(2,2,2), horiz = T, cex = 1)	
    grid(lwd = 1.5)

    #Plot2
    plot(mydata3$Normalized, pch = 16, xlab = "Targets", ylab = "AllMuts / Unmutated", 
         ylim = c(0,1),main = "Normalized Saturation")
    grid(lwd = 1.5)
    lines(mydata3$Normalized, lwd = 2)
    
    
     #Plot3	
     plot(x = 1:length(Inter_types$Type), y = Inter_types$Count, las = 3, xlab= "Type", ylab= "Count",	
          main = "Total number of mutations", axes=F, pch=8,
          col= c(rep(x = 1,9),rep(x = 2,8),rep(x = 3,7),rep(x = 4,6),rep(x = 5,5),
                 rep(x = 6,4),rep(x = 7,3),rep(x = 8,2),rep(x = 9,1)))	
     axis(side = 1,labels =Inter_types$Type, at= 1:length(Inter_types$Type),las =3 )
     axis(side=2)
     box()
     grid(lwd = 1.5)
     # Plot4
     plot(mydata3$Proportion, ylim = c(0,1), xlab = "Targets", ylab = "Saturation", 
          col="red",pch=20,lwd=2)
     lines(mydata3$Proportion, col="red",lwd=2)
     grid(lwd = 1.5)
     # plot(Intra.m[,2]/sum(Intra.m[,2],na.rm = T), main = "Intra Mut types",	
     #      ylab = "Freq",ylim = c(0,0.25), xlab = "Target",col="blue",lwd=2)	
     # lines(Intra.m[,2]/sum(Intra.m[,2],na.rm = T),lwd=2,col="blue",lty=2)	
     #Plot5	
     #Plot5
     plot(as.factor(mydata$Size),main = "Mutation types by size", xlab= "Size", ylab= "Count",
          col=c(1:9))
    #Plot6
      plot(Inter.m[,4]/sum(Inter.m[,4],na.rm = T),	ylim = c(0,0.5),
          xlab = "Size", ylab = "Freq", main = "Total number of Mutations per size",
          col=c(0:9),pch=22,lwd=5)	
     lines(Inter.m[,4]/sum(Inter.m[,4],na.rm = T),lty=2)	
     legend(x =8 ,y = 0.5,legend = c(1:9),cex = 1,pch=22,horiz = F,
            col = c(1:9),fill = c(1:9),title = "SIZE")
     grid(lwd = 1.5)
     
    # #Plot6	
    #  plot(mydata$Type,main = "Different types of mutations", xlab= "Size", ylab= "Count", las= 3,
    #       col= c(rep(x = 1,9),rep(x = 2,8),rep(x = 3,7),rep(x = 4,6),rep(x = 5,5),
    #              rep(x = 6,4),rep(x = 7,3),rep(x = 8,2),rep(x = 9,1)))	
    dev.off()	
    
    ##################### PLOT Freqs of mutations per target
    
    svg(filename = paste("30hr_",s,"_1x_Mutations_Freqs_by_Target_corrected.svg",sep=""),width = 9,height = 6)	
    par(mfrow=c(3,4))	
    for (i in 1:10)
    {
      plot(sort(mydata2$Count[mydata2$Target==i]/sum(mydata2$Count[mydata2$Target==i])), las = 1,
           pch = 1, ylab = "Count", main = paste("Target #",i,sep=""),lwd= 0.6,xlab= "")
      lines(sort(mydata2$Count[mydata2$Target==i]/sum(mydata2$Count[mydata2$Target==i])), lwd= 1)
    }
    dev.off()
    
    ##################### DATAFRAME with Freqs of mutations per target
    Freqs_df <- matrix(data = NA,nrow = 10,ncol = 60)
    for (i in 1:10)
    {
      ss = sort(mydata2$Count[mydata2$Target==i]/sum(mydata2$Count[mydata2$Target==i]),decreasing = T)
      Freqs_df[i,] = ss[1:60]
    }
    
    svg(filename = paste("30hr_",s,"_1x_Mutations_Freqs_ALL_Targets.svg",sep=""),width = 3.5,height = 4)
    par(mfrow=c(1,1))	
    plot(NA,xlim = c(1,60), ylim = c(0,0.5),ylab = "Freqs", xlab = "Mutated states")
    boxplot(Freqs_df,outline = F,axes=F,main = paste("Replicate #",s,sep=""), add = T)
    lines(x = 1:60, y =apply(Freqs_df, 2, function(x) mean(x,na.rm=T)), col= "red", lwd=3)
    axis(side = 1,at = c(10,20,30,40,50,60))    
    axis(side = 2, at= c(0,0.1,0.2,0.30,0.4))
    box()
    dev.off()
    
    ALL_MUT_FREQ[s,] = apply(Freqs_df, 2, function(x) mean(x,na.rm=T))
    
}	

write.table(x = ALL_MUT_FREQ, file = "All_muts_FREQ.txt",
            quote = F, col.names = F )

### add gamma dist line
SIM_FREQS = matrix(data = 0,nrow = 1000,ncol = 60)
for (g in 1:1000)
{
  x = rgamma(60, shape = 0.1, scale = 2)
  SIM_FREQS[g,] <- sort(x,decreasing = T)/sum(x)
}
boxplot(SIM_FREQS, outline=F)
lines(apply(SIM_FREQS, 2, mean), col="blue", lwd=3)


#####
svg(filename = "Mutational_complexity_Obs_vs_Sim.svg",width = 4.5,height = 3.5)

plot(NA,xlim = c(1,60), ylim = c(0,0.3),ylab = "Frequency", xlab = "Mutated states",  main = "Mutational complexity",axes =F)
#boxplot(ALL_MUT_FREQ, outline= F, axes=F, add= T)
lines(x = 1:60, y =apply(ALL_MUT_FREQ, 2, mean), col= "red", lwd=5)
#lines(x = 1:60, y =ALL_MUT_FREQ[s,], col= "red", lwd=4)
axis(side = 1,at = c(10,20,30,40,50,60))    
axis(side = 2, at= c(0,0.1,0.2,0.30,0.4))
lines(apply(SIM_FREQS, 2, mean), col="blue", lwd=5, lty = 2)
box()

legend(x = 6, y = 0.25,c("GESTALT v7\n(replicate/target mean)", expression(paste("Gamma dist (",kappa,"= 0.1,",theta,"= 2)"))), 
       horiz=F, lty= c(1,2), lwd = 3, col=c("red", "blue"), cex = 1.1,bty = "n")

dev.off()
###
svg(filename = "INTRA+INTER_mutations_30hr_corrected.svg",width = 12,height = 3)

par(mfrow=c(1,3))	
#### PLOT 1
#boxplot(ALL_MUTS,use.cols = F,boxwex=0.2,ylim=c(0,0.24), at = c((1:10)+0))
boxplot(ALL_INTER,use.cols = F,boxwex=0.1,ylim=c(0,0.2), at = c((1:10)-0),
        col="orange",outline=F, xlab="Target",ylab="Proportion of mutations")
boxplot(ALL_INTRA,use.cols = F,boxwex=0.1,ylim=c(0,0.25), at = c((1:10)+0.15),
         add=T,axes=F,col="blue",outline=F)

lines(apply(X = ALL_MUTS,MARGIN = 1,FUN = median),lwd= 2,col= "black")
points(apply(X = ALL_MUTS,MARGIN = 1,FUN = median),lwd= 2,col= "black")
sum(ALL_MUTS)

lines(apply(X = ALL_INTER,MARGIN = 1,FUN = median),lwd= 2,col= "orange",lty=6)
points(apply(X = ALL_INTER,MARGIN = 1,FUN = median),lwd= 2,col= "orange",lty=2,pch=6)

lines(apply(X = ALL_INTRA,MARGIN = 1,FUN = median),lwd= 2,col= "blue",lty=6)
points(apply(X = ALL_INTRA,MARGIN = 1,FUN = median),lwd= 2,col= "blue",lty=2,pch=2)

legend(x = 1.5, y = 0.2,legend = c("Intra     ", "Inter     ", "All"),
col = c("blue","orange","black"),lwd = c(2,2,2),
horiz = T, cex = 0.9,lty=c(6,6,1),pch=c(2,6,1))	
    
###### NORMALIZED SATURATION OF ALL TARGETS
boxplot(ALL_NORM,use.cols = F,boxwex=0.3,ylim=c(0,1), at = c(1:10),
        col="white",outline=F, xlab="Target",ylab="AllMuts / Unmutated", 
        main = "Normalized Saturation")
lines(apply(X = ALL_NORM,MARGIN = 1,FUN = median),lwd= 2,col= "purple")
points(apply(X = ALL_NORM,MARGIN = 1,FUN = median),lwd= 2,col= "purple")
grid(lwd = 1.5)
###### SATURATION OF ALL TARGETS
boxplot(ALL_SAT,use.cols = F,boxwex=0.3,ylim=c(0,1), at = c(1:10),
        col="white",outline=F, xlab="Target",ylab="Unmutated / Reads", main= "Saturation")
lines(apply(X = ALL_SAT,MARGIN = 1,FUN = median),lwd= 2,col= "red")
points(apply(X = ALL_SAT,MARGIN = 1,FUN = median),lwd= 2,col= "red")
grid(lwd = 1.5)

#######

dev.off()


######## EXPORT TO TABLES
# Proportion of mutations in each target
write.table(x = apply(X = ALL_MUTS,MARGIN = 1,FUN = median), file = "All_muts_median_corrected.txt",
            quote = F,row.names = c(1:10), col.names = F )

# Normalized Saturation of the targets
write.table(x = apply(X = ALL_NORM,MARGIN = 1,FUN = median), file = "All_Normalized_saturation_median_corrected.txt",
            quote = F,row.names = c(1:10), col.names = F )

# Saturation of the targets

write.table(x = apply(X = ALL_SAT,MARGIN = 1,FUN = median), file = "All_saturation_median_corrected.txt",
            quote = F,row.names = c(1:10), col.names = F )

