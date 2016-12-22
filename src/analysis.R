#PCA analysis of results. 
library( ANTsR )
library( geomorph )
library( abind )
library( Morpho )
# set the base directory manually. This should point out to the data folder of the cloned repository.
# calculations require large amount of RAM (64+ GB) to be available.
bd=path.expand( "/scratch/mouse_CT_atlas/data/")
reference=paste( bd, "/templates/low_res_skull_only.nii.gz",sep="")

if ( ! dir.exists( bd ) )
  stop("set base directory to point to the data folder of the cloned repository.")
sims = Sys.glob( paste(bd,'/targets/masked/registered/Skull/step1/*Similarity*.mat',sep='') )
mats = Sys.glob( paste(bd,'/targets/masked/registered/Skull/step1/step2/*Warped1Affine.mat',sep='') )
defs = Sys.glob( paste(bd,'/targets/masked/registered/Skull/step1/step2/*[0-9]Warp.nii.gz',sep='') )
total_warps = Sys.glob(paste(bd,"/targets/masked/registered/Skull/step1/step2/aff_def/*.nii.gz",sep=''))
imgs = Sys.glob( paste(bd,'/targets/masked/registered/Skull/step1/step2/*WarpedWarped.nii.gz',sep='') )
template = antsImageRead(reference)
msk = thresholdImage(template,2,255, 1,0)     
#msk = iMath(msk,"MD",1)      #to dilate the mask 
samples=defs
samples=sub("_Warped2Warp.nii.gz","",samples)
samples=sub(paste(bd,'/targets/masked/registered/Skull/step1/step2/',sep=''),"",samples)
samples=sub("Skull_","",samples)

# read the landmark collected from the skull template and map them back to samples.
templatePointsFile = paste(bd,'/templates/low_res_skull_only_LM_LPS.csv',sep='')
templatePoints = read.csv( templatePointsFile )
nlandmarks = nrow( templatePoints )
lmlist = list( )
for ( i in 1:length( sims ) ) {
  lotx = c( defs[i], mats[i], sims[i] )
  wpts <- antsApplyTransformsToPoints( dim=3, points=templatePoints,
                                       transformlist=lotx )
  lmlist[[ i ]] = wpts
}
mydims = c( nrow( lmlist[[1]] ), ncol( lmlist[[1]]), length(lmlist) )
pseudoLM = array( unlist( lmlist ), dim = mydims )
#radius = mean( antsGetSpacing(msk) ) * 5     #these steps are there to confirm positioning of pseudoLMs on the skull template
#mmm = makePointsImage( templatePoints, msk, radius )
#print( range(mmm )  )  # should have max value equal to nlandmarks
#####

#load LM data
load(paste(bd,"/LMs.RData",sep="/")) #LMs manually annotated  from the individual skull.
GPA=gpagen(LM)
GPA.pca=plotTangentSpace(GPA$coords)
exclude.lm=which(! as.character(1:51) %in% rownames(LM[,,1]))   #Total of six LMs were missing from various mice. Exlcuding those from the pseudoLMs.
pseudoGPA=gpagen(pseudoLM[-exclude.lm,,])
pseudoGPA.pca=plotTangentSpace(pseudoGPA$coords)
joint=abind(LM,pseudoLM[-exclude.lm,,])
joint.GPA=gpagen(joint)

par(mfrow=c(1,2))
plot(GPA$Csize, pseudoGPA$Csize,xlab="Manual LMs Centroid Size",ylab="Pseudo LM Centroid Size",pch=20)
plot(GPA.pca$pc.scores[,1],pseudoGPA.pca$pc.scores[,1],xlab="Manual LMs PC1",ylab="Pseudo LM PC1",pch=20)
for (i in 1:5) print(cor.test(GPA.pca$pc.scores[,i],pseudoGPA.pca$pc.scores[,i]))

#MEASURES OF SIZE
Skulls=dir(path=paste(bd,'/targets/masked/registered/Skull/',sep=''),patt="gz",full.names = T)
voxel_count=function(X){
  my=antsImageRead(X)
  my_thresh=thresholdImage(my,25,"Inf",1,0)
  stat=labelStats(my,my_thresh)
  return(stat$Count[2])
}
skull_sizes=NULL
for (i in Skulls) skull_sizes=c(skull_sizes, voxel_count(i))


#Calculate deformation fields
odr = paste(bd,"/targets/masked/registered/Skull/step1/step2/",sep='')
jmean = rep( NA, length( defs ) )
jlist = list()
ylist = list( )
for ( i in 1:length( defs ) ) {
  # build this from composite maps
  myscale = getAntsrTransformParameters( readAntsrTransform( sims[ i ] ) )[7]
  print( paste( i, myscale, skull_sizes[i] ) )
  mymat = readAntsrTransform( mats[i] )
  fixedparms = getAntsrTransformFixedParameters( mymat )
  mymatparms = getAntsrTransformParameters( mymat )
  mymatparms[ c(1,5,9) ] = mymatparms[ c(1,5,9) ] * myscale
  newtx = createAntsrTransform( type = "AffineTransform", precision = "float",
                                dimension = 3, parameters=mymatparms, fixed.parameters=fixedparms )
  newmatfn = paste( odr, samples[i], "_scaledaff.mat",sep='')
  writeAntsrTransform( newtx, newmatfn )
  mytx = c( defs[i], newmatfn )
  ofn = paste( odr, samples[i], "_",sep='')
  compfield = antsApplyTransforms( template, template, transformlist=mytx, compose=ofn )
  jlist[[ i ]] = createJacobianDeterminantImage( template, compfield )
  jmean[ i ] = mean( jlist[[ i ]][ msk == 1 ] )
  ylist[[ i ]] = antsImageRead( compfield )
}
print(  cor.test( jmean , skull_sizes ) )  # should be high

wlist = list( )
for ( i in 1:length( defs ) ) wlist[[ i ]] = antsImageRead( defs[i] )


#Compute the PCA basis on SyN portion only:

pca_type='randPCA'  #randomized PCA requires rsvd package to be installed

if ( !exists("diGPA_SyN_PCA") )
  ####
  diGPA_SyN_PCA = multichannelPCA( wlist, msk, pcaOption=pca_type, verbose=FALSE )

plot(diGPA_SyN_PCA$pca$u[,1:2],pch=20,cex=2,main="SyN Only PCA",xlab="PC1", ylab="PC2")
text(diGPA_SyN_PCA$pca$u[,1:2],labels=samples,cex=.6,pos=2)


#Compute the PCA basis on total warp :
if ( !exists("diGPA_total_PCA") )
  ####
  #diGPA_total_PCA= multichannelPCA( ylist, msk, k=myk, pcaOption=pca_type, verbose=FALSE )
  diGPA_total_PCA= multichannelPCA( ylist, msk, pcaOption=pca_type, verbose=FALSE )

plot(diGPA_total_PCA$pca$u[,1:2],pch=20,cex=2,main="Total Warps PCA",xlab="PC1", ylab="PC2")
text(diGPA_total_PCA$pca$u[,1:2],labels=samples,cex=.6,pos=2)

#Now visualize the first two PC (positive component) and save the as volumes
#these steps require rgl, misc3d and pixmap libraries installed in R

mxk=1 #number of PCs to visualize
for ( ww in 1:mxk )
{
  print( ww )
  myw = diGPA_SyN_PCA$pcaWarps[[ww]] * 50  #arbitrarily scaled
  #myw = smoothImage( myw, antsGetSpacing( template ) * 5.0 )
  warpTx = antsrTransformFromDisplacementField(  myw  )
  # compose several times to get a visual effect
  wtxlong=warpTx
  warped = applyAntsrTransform( wtxlong, data = template,
                                reference = template )
  plot(template)
  title('template')
  plot(warped)    #plot the deformed axis as overlay on template
  title(paste("Warp for +PC",ww))
  # look at the magnitude ...
  splt = splitChannels( diGPA_SyN_PCA$pcaWarps[[ww]]  )
  mgntd = iMath( abs( splt[[1]] ) + abs( splt[[2]]) + abs( splt[[3]]  ), "Normalize" )
  par(mfrow=c(3,1))
  for (k in 1:3) plot( template, mgntd, window.overlay=c(0.1,1),axis=k )
  title(paste("Magnitude for +PC",ww))
  
}

#R's default save function will not save the values stored in pointers 
#please uncomments the line below, if you want to save your analysis. 
#save.ANTsR("~/PCA_results")      

#to reload your session
#load.ANTsR("~/PCA_results")
#and rerun Lines 6-22 to correct the filenames. 


#plots and tables in the paper
#Figure 3.
plot(joint.GPA$consensus[,1:2],pch=32,col="black",cex=2)
for (i in 1:51) points(joint.GPA$coords[,1:2,i],pch=20,col="red")
for (i in 52:102) points(joint.GPA$coords[,1:2,i],pch=20,col="cyan")
#points(joint.GPA$consensus[,1:2],pch=3,col="black",cex=2,lwd=4)
pseudo_Mean=apply(joint.GPA$coords[,,52:102],c(1,2),mean)
GPA_Mean=apply(joint.GPA$coords[,,1:51],c(1,2),mean)
points(GPA_Mean[,1:2],col="black",pch=3,cex=2,lwd=2)
points(pseudo_Mean[,1:2],col="black",pch=3,cex=2,lwd=2)

groups=as.factor(c(rep("GPA",51),rep("Pseudo",51)))
proc.anova=procD.lm(joint.GPA$coords~groups,iter = 1000,RRPP = T)



#Figure 4.
par(mfrow=c(1,2))
plot(GPA$Csize,pseudoGPA$Csize,pch=20,cex=2,cex.axis=1.5,ylab="pseudoLM CS",xlab="LM CS")
plot(log(skull_sizes)~log(GPA$Csize),pch=20,cex=2,cex.axis=1.5,ylab="log Skull Voxel Count",xlab="log LM Centroid Size")


#Figure 5.
par(mfrow=c(2,2))
plot(round((GPA.pca$sdev^2)/sum(GPA.pca$sdev^2),3),pch=20,cex=2,cex.axis=2,main="GPA PCA",
     ylab="Variance Explained",xlab="Principal Components")
plot(round((GPA.pca$sdev^2)/sum(GPA.pca$sdev^2),3),pch=20,cex=2,cex.axis=2,main="Pseudo LM PCA",
     ylab="Variance Explained",xlab="Principal Components")
plot(round((diGPA_SyN_PCA$pca$d^2)/sum(diGPA_SyN_PCA$pca$d^2),3),pch=20,cex=2,cex.axis=2,main="diGPA SyN PCA",
     ylab="Variance Explained",xlab="Principal Components")
plot(round((diGPA_total_PCA$pca$d^2)/sum(diGPA_total_PCA$pca$d^2),3),pch=20,cex=2,cex.axis=2,main="diGPA Total Warp PCA",
     ylab="Variance Explained",xlab="Principal Components")


#Figure 6.
par(mfrow=c(2,2))
plot(GPA.pca$pc.scores[,1],pseudoGPA.pca$pc.scores[,1],pch=20, xlab="GPA PC1",ylab="Pseudo LM PC1")
plot(GPA.pca$pc.scores[,1],diGPA_SyN_PCA$pca$u[,1],pch=20, xlab="GPA PC1",ylab="SyN PC1")
plot(GPA.pca$pc.scores[,1],diGPA_total_PCA$pca$u[,1],pch=20, xlab="GPA PC1",ylab="Total Warp PC1")
plot(diGPA_SyN_PCA$pca$u[,1],diGPA_total_PCA$pca$u[,1],pch=20, xlab="SyN PC1",ylab="Total Warp PC1")


#correlations between corresponding PC scores.
for (i in 1:1) {
  print(cor.test(GPA.pca$pc.scores[,i],pseudoGPA.pca$pc.scores[,i]))
  print(cor.test(GPA.pca$pc.scores[,i],diGPA_SyN_PCA$pca$u[,i]))
  print(cor.test(GPA.pca$pc.scores[,i],diGPA_total_PCA$pca$u[,i]))
  print(cor.test(diGPA_SyN_PCA$pca$u[,i],diGPA_total_PCA$pca$u[,i]))
  
}
