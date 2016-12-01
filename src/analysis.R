MM=proc.time()
library( ANTsR )
library( geomorph )
# set the base directory manually. This should be the output of data folder from the ISA checkout
# calculations require large amount of RAM (64+ GB) to be available.
bd=path.expand( "/scratch/isa/data/")
reference=paste( bd, "/mouse/templates/low_res_skull_only.nii.gz",sep="")

if ( ! dir.exists( bd ) )
  stop("set base directory to point to the base of the ISA data folder")
sims = Sys.glob( paste(bd,'/mouse/targets/masked/registered/Skull/step1/*Similarity*.mat',sep='') )
mats = Sys.glob( paste(bd,'/mouse/targets/masked/registered/Skull/step1/step2/*Warped1Affine.mat',sep='') )
defs = Sys.glob( paste(bd,'/mouse/targets/masked/registered/Skull/step1/step2/*[0-9]Warp.nii.gz',sep='') )
total_warps = Sys.glob(paste(bd,"/mouse/targets/masked/registered/Skull/step1/step2/aff_def/*.nii.gz",sep=''))
imgs = Sys.glob( paste(bd,'/mouse/targets/masked/registered/Skull/step1/step2/*WarpedWarped.nii.gz',sep='') )
template = antsImageRead(reference)
msk = thresholdImage(template,2,255, 1,0)     #define a mask for the bone only areas by thresholding the template. deformation fields will be calculated only for this mask.
#msk = iMath(msk,"MD",1)
samples=defs
samples=sub("_Warped2Warp.nii.gz","",samples)
samples=sub(paste(bd,'/mouse/targets/masked/registered/Skull/step1/step2/',sep=''),"",samples)
samples=sub("Skull_","",samples)

# get points and map them
templatePointsFile = paste(bd,'/mouse/templates/low_res_skull_only_LM_LPS.csv',sep='')
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
# not sure if the above is the correct format for GPA but this is the idea
pseudoLM = array( unlist( lmlist ), dim = mydims )
# obviously, this is not correct until the landmarks are in the right space
# we need to verify that this following code produces the "right" output
# by creating the mmm image and rendering it in itk-snap on the template
radius = mean( antsGetSpacing(msk) ) * 5
mmm = makePointsImage( templatePoints, msk, radius )
print( range(mmm )  )  # should have max value equal to nlandmarks
#DIGPA = gpagen( pseudoLM )
#####

#load LM data
load(paste(bd,"mouse/LMs.RData",sep="/"))
GPA=gpagen(LM)
GPA.pca=plotTangentSpace(GPA$coords)
exclude.lm=which(! as.character(1:51) %in% rownames(LM[,,1]))
diGPA=gpagen(pseudoLM[-exclude.lm,,])
diGPA.pca=plotTangentSpace(diGPA$coords)
plot(GPA$Csize,diGPA$Csize,main="Centroid Size")
plot(GPA.pca$pc.scores[,1],diGPA.pca$pc.scores[,1])
for (i in 1:10) print(cor.test(GPA.pca$pc.scores[,i],diGPA.pca$pc.scores[,i])$estimate)



#MEASURES OF SIZE
Skulls=dir(path=paste(bd,'/mouse/targets/masked/registered/Skull/',sep=''),patt="gz",full.names = T)
fun=function(X){
  my=antsImageRead(X)
  my_thresh=thresholdImage(my,25,"Inf",1,0)
  stat=labelStats(my,my_thresh)
  return(stat$Count[2])
}
counts=NULL
for (i in Skulls) counts=c(counts, fun(i))

odr = paste(bd,"/mouse/targets/masked/registered/Skull/step1/step2/",sep='')
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
if ( !exists("IBSA_SyN_PCA") )
  ####
IBSA_SyN_PCA = multichannelPCA( wlist, msk, k=myk, pcaOption=pca_type, verbose=FALSE )
plot(IBSA_SyN_PCA$pca$u[,1:2],pch=20,cex=2,main="SyN Only PCA")
text(IBSA_SyN_PCA$pca$u[,1:2],labels=samples,cex=.6,pos=2)

#Compute the PCA basis on total warp :
if ( !exists("IBSA_total_PCA") )
  ####
IBSA_total_PCA= multichannelPCA( ylist, msk, k=myk, pcaOption=pca_type, verbose=FALSE )

summary( lm( log(counts) ~ IBSA_total_PCA$pca$u[,1:myk] ) )

plot(IBSA_total_PCA$pca$u[,1:2],pch=20,cex=2,main="Total Warps PCA")
text(IBSA_total_PCA$pca$u[,1:2],labels=samples,cex=.6,pos=2)

print(proc.time()-MM)

#save.ANTsR("/home/maga/mask2_PCA")

#Now visualize the first two PC (positive component) and save the as volumes
#these require rgl, misc3d and pixmap libraries installed in R

mxk=1 #number of PCs to visualize
for ( ww in 1:mxk )
{
  print( ww )
  myw = IBSA_SyN_PCA$pcaWarps[[ww]] * 50  #arbitrarily scaled
  #myw = smoothImage( myw, antsGetSpacing( template ) * 5.0 )
  warpTx = antsrTransformFromDisplacementField(  myw  )
  # compose several times to get a visual effect
  #wtxlong = list( ) ; for ( i in 1:20 ) wtxlong[[i]]=warpTx
  wtxlong=warpTx
  warped = applyAntsrTransform( wtxlong, data = template,
                                reference = template )
  plot(template)
  title('template')
  plot(warped)    #plot the deformed axis as overlay on template
  title(paste("Warp for +PC",ww))
  # look at the magnitude ...
  splt = splitChannels( IBSA_SyN_PCA$pcaWarps[[ww]]  )
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
plot(GPA.pca$pc.scores[,1],IBSA_SyN_PCA$pca$u[,1],xlab="GPA PC1",ylab="IBSA SyN PC1",pch=20,cex=2,cex.axis=2)
plot(GPA.resi.pca$pc.scores[,1],IBSA_SyN_PCA$pca$u[,1],xlab="GPA resi PC1",ylab="IBSA SyN PC1",pch=20,cex=2,cex.axis=2)
plot(GPA.pca$pc.scores[,1],IBSA_total_PCA$pca$u[,1],xlab="GPA PC1",ylab="IBSA totalWarp PC1",pch=20,cex=2,cex.axis=2)
plot(GPA.resi.pca$pc.scores[,1],IBSA_total_PCA$pca$u[,1],xlab="GPA resi PC1",ylab="IBSA totalWarp PC1",pch=20,cex=2,cex.axis=2)

#correlations
cor(GPA.pca$pc.scores[,1],IBSA_SyN_PCA$pca$u[,1])
cor(GPA.resi.pca$pc.scores[,1],IBSA_SyN_PCA$pca$u[,1])
cor(GPA.pca$pc.scores[,1],IBSA_total_PCA$pca$u[,1])
cor(GPA.resi.pca$pc.scores[,1],IBSA_total_PCA$pca$u[,1])

par(mfrow=c(2,2))
plot(round((IBSA_SyN_PCA$pca$d^2)/sum(IBSA_SyN_PCA$pca$d^2),3),pch=20,cex=2,cex.axis=2,main="IBSA SyN PCA",
      ylab="Variance Explained",xlab="Principal Components")
plot(round((IBSA_total_PCA$pca$d^2)/sum(IBSA_total_PCA$pca$d^2),3),pch=20,cex=2,cex.axis=2,main="IBSA Total Warp PCA",
     ylab="Variance Explained",xlab="Principal Components")
plot(round((GPA.pca$sdev^2)/sum(GPA.pca$sdev^2),3),pch=20,cex=2,cex.axis=2,main="GPA PCA",
     ylab="Variance Explained",xlab="Principal Components")
plot(round((GPA.resi.pca$sdev^2)/sum(GPA.resi.pca$sdev^2),3),pch=20,cex=2,cex.axis=2,main="Allometry free GPA PCA",
     ylab="Variance Explained",xlab="Principal Components")

plot(log(counts)~log(GPA$Csize),pch=20,cex=2,cex.axis=1.5,ylab="Skull Voxel Count",xlab="LM Centroid Size")

#how does skull size impact the PC1 scores
regress_IBSA_total=lm(IBSA_total_PCA$pca$u[,1]~log(counts))
regress_IBSA_SyN=lm(IBSA_SyN_PCA$pca$u[,1]~log(counts))
regress_GPA_pca=lm(GPA.pca$pc.scores[,1]~log(GPA$Csize))
regress_GPA_resi.pca=lm(GPA.resi.pca$pc.scores[,1]~log(GPA$Csize))
