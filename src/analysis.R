MM=proc.time()
library( ANTsR )
library( geomorph )
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
load(paste(bd,"mouse/LMs.RData",sep="/")) #LMs manually annotated  from the individual skull.
GPA=gpagen(LM)
GPA.pca=plotTangentSpace(GPA$coords)
exclude.lm=which(! as.character(1:51) %in% rownames(LM[,,1]))   #Total of six LMs were missing from various mice. Exlcuding those from the pseudoLMs.
diGPA=gpagen(pseudoLM[-exclude.lm,,])
diGPA.pca=plotTangentSpace(diGPA$coords)
plot(GPA$Csize,diGPA$Csize,main="Centroid Size")
plot(GPA.pca$pc.scores[,1],diGPA.pca$pc.scores[,1])
for (i in 1:10) print(cor.test(GPA.pca$pc.scores[,i],diGPA.pca$pc.scores[,i])$estimate) #correlations between PC scores.

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

odr = paste(bd,"/targets/masked/registered/Skull/step1/step2/",sep='')
jmean = rep( NA, length( defs ) )
jlist = list()
ylist = list( )
for ( i in 1:length( defs ) ) {
  # build this from composite maps
  myscale = getAntsrTransformParameters( readAntsrTransform( sims[ i ] ) )[7]
  print( paste( i, myscale, counts[i] ) )
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
print(  cor.test( jmean , counts ) )  # should be high

wlist = list( )
for ( i in 1:length( defs ) ) wlist[[ i ]] = antsImageRead( defs[i] )


#Compute the PCA basis on SyN portion only:

pca_type='randPCA'  #randomized PCA requires rsvd package to be installed
myk = 5
if ( !exists("diGPA_SyN_PCA") )
  ####
diGPA_SyN_PCA = multichannelPCA( wlist, msk, k=myk, pcaOption=pca_type, verbose=FALSE )
plot(diGPA_SyN_PCA$pca$u[,1:2],pch=20,cex=2,main="SyN Only PCA")
text(diGPA_SyN_PCA$pca$u[,1:2],labels=samples,cex=.6,pos=2)

#Compute the PCA basis on total warp :
if ( !exists("diGPA_total_PCA") )
  ####
diGPA_total_PCA= multichannelPCA( ylist, msk, k=myk, pcaOption=pca_type, verbose=FALSE )

summary( lm( log(counts) ~ diGPA_total_PCA$pca$u[,1:myk] ) )

plot(diGPA_total_PCA$pca$u[,1:2],pch=20,cex=2,main="Total Warps PCA")
text(diGPA_total_PCA$pca$u[,1:2],labels=samples,cex=.6,pos=2)

print(proc.time()-MM)


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

#allometry code from geomorph
allometry <- procD.lm(GPA$coords~GPA$Csize,iter = 1)
GPA.resid <- arrayspecs(allometry$residuals,
                          p=dim(GPA$coords)[1], k=dim(GPA$coords)[2]) # size-adjusted residuals
adj.shape <- GPA.resid + array(GPA$consensus, dim(GPA.resid)) # allometry-free shapes
GPA.resi.pca=plotTangentSpace(adj.shape) # PCA of allometry-free shape

par(mfrow=c(2,2))
plot(GPA.pca$pc.scores[,1],diGPA_SyN_PCA$pca$u[,1],xlab="GPA PC1",ylab="diGPA SyN PC1",pch=20,cex=2,cex.axis=2)
plot(GPA.resi.pca$pc.scores[,1],diGPA_SyN_PCA$pca$u[,1],xlab="GPA resi PC1",ylab="diGPA SyN PC1",pch=20,cex=2,cex.axis=2)
plot(GPA.pca$pc.scores[,1],diGPA_total_PCA$pca$u[,1],xlab="GPA PC1",ylab="diGPA totalWarp PC1",pch=20,cex=2,cex.axis=2)
plot(GPA.resi.pca$pc.scores[,1],diGPA_total_PCA$pca$u[,1],xlab="GPA resi PC1",ylab="diGPA totalWarp PC1",pch=20,cex=2,cex.axis=2)

#correlations between corresponding PC scores.
for (i in 1:10) {
  cor.test(GPA.pca$pc.scores[,i],diGPA_SyN_PCA$pca$u[,i])
  cor.test(GPA.resi.pca$pc.scores[,i],diGPA_SyN_PCA$pca$u[,i])
  cor.test(GPA.pca$pc.scores[,i],diGPA_total_PCA$pca$u[,i])
  cor.test(GPA.resi.pca$pc.scores[,i],diGPA_total_PCA$pca$u[,i])
}

par(mfrow=c(2,2))
plot(round((diGPA_SyN_PCA$pca$d^2)/sum(diGPA_SyN_PCA$pca$d^2),3),pch=20,cex=2,cex.axis=2,main="diGPA SyN PCA",
      ylab="Variance Explained",xlab="Principal Components")
plot(round((diGPA_total_PCA$pca$d^2)/sum(diGPA_total_PCA$pca$d^2),3),pch=20,cex=2,cex.axis=2,main="diGPA Total Warp PCA",
     ylab="Variance Explained",xlab="Principal Components")
plot(round((GPA.pca$sdev^2)/sum(GPA.pca$sdev^2),3),pch=20,cex=2,cex.axis=2,main="GPA PCA",
     ylab="Variance Explained",xlab="Principal Components")
plot(round((GPA.resi.pca$sdev^2)/sum(GPA.resi.pca$sdev^2),3),pch=20,cex=2,cex.axis=2,main="Allometry free GPA PCA",
     ylab="Variance Explained",xlab="Principal Components")

plot(log(counts)~log(GPA$Csize),pch=20,cex=2,cex.axis=1.5,ylab="Skull Voxel Count",xlab="LM Centroid Size")

#how does skull size impact the PC1 scores
regress_diGPA_total=lm(diGPA_total_PCA$pca$u[,1]~log(counts))
regress_diGPA_SyN=lm(diGPA_SyN_PCA$pca$u[,1]~log(counts))
regress_GPA_pca=lm(GPA.pca$pc.scores[,1]~log(GPA$Csize))
regress_GPA_resi.pca=lm(GPA.resi.pca$pc.scores[,1]~log(GPA$Csize))


#R's default save function will not save the values stored in pointers 
#please uncomments the line below, if you want to save your analysis. 
#save.ANTsR("~/PCA_results")      

#to reload your session
#load.ANTsR("~/PCA_results")
#and rerun Lines 6-22 to correct the filenames. 


