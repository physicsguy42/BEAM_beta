# Script to generate 3d plot of particle field
rm(list=ls()) # clear workspace
library(rgl)
myframe = read.table("/NASA2013/Radiative_transfer/Beam07292014/BEAM_beta4a/nowakes_D01_tau10_nbins40.partdat",skip=3)
mymatrix<-as.matrix(myframe,rownames.force = NA)
cat("number of rows is ", nrow(mymatrix))
npart = nrow(mymatrix)
Rvals1=mymatrix[1:npart,1]
Xvals1=mymatrix[1:npart,2]
Yvals1=mymatrix[1:npart,3]
Zvals1=mymatrix[1:npart,4]
L2 = max(Xvals1)
zmax = max(Zvals1)
# zmax = 200
cat("maximum value of Z is ", max(Zvals1))
cat("maximum value of X is ", max(Xvals1))
# light3d(theta = 0,phi=15,diffuse="gray75")
material3d(lit=FALSE,shininess=0,alpha=0.9)
#plot3d(x=Xvals1,y=Yvals1,z=Zvals1,type="s", col="grey", xlab="X", ylab="Y", zlab="Z", size=1,radius=Rvals1, box=T,axes=T,xlim=c(-L2,L2),ylim=c(-L2,L2),zlim=c(-zmax,zmax))
plot3d(x=Xvals1,y=Yvals1,z=Zvals1,type="s", col="blue", xlab="X", ylab="Y", zlab="Z", size=1,radius=Rvals1, box=T,axes=T,xlim=c(-L2,L2),ylim=c(-L2,L2),zlim=c(-zmax,zmax))

 
