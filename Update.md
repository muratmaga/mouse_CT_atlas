Since the publication, R implementation of the ANTs (ANTsR) began to support more registration options, which wasn't available at the time of publication, hence we used a combination of R, and shell scripts to call ANTs tools. The whole workflow can now be completely done in R/ANTsR, simplfying it greatly. Below is an example where fixed anatomical landmarks in the template space are transfered to individual skulls using deformable image registration:

```
require(ANTsR)
template="path.to.the.template.volume" #use this template to work with the provided data.
ref.LMs=read.csv(file="path.to.the.template.LMs.fcsv", skip=2)[,2:4]

#used to output landmark directly in Slicer compatible format.
source("https://raw.githubusercontent.com/muratmaga/SlicerMorph_Rexamples/main/write.markups.fcsv.R")

low.res=c(0.1, 0.1, 0.1)
lowres.templ = resampleImage(antsImageRead(template), low.res)

setwd("path.to.the.folder.with.volumes.tobe.landmarked")
f=dir(patt='nii.gz')

for (X in 1:length(f)) {
  img=antsImageRead(f[X])
  
  #simple connected component to remove non-cranial elements
  temp=thresholdImage(img,25, "Inf",1,0)
  temp=iMath(temp,"GetLargestComponent")
  img=temp*img
  
  
  img2 = resampleImage( img, low.res, interpType=0 )
  
  nout = 5
  thetas = seq( from = 0, to = 360, length.out = nout )[1:(nout-1)]
  
# find an close enough starting position for samples with very different orientations
  
  mival<-invariantImageSimilarity( lowres.templ, img2,
                                   thetas=thetas, thetas2=thetas, thetas3=thetas,metric="GC",
                                   localSearchIterations=20, transform="Rigid" )
    
  # deformable image registration
  reg = antsRegistration(fixed=lowres.templ, 
                              moving = img2, 
                              initialTransform=mival[[2]], 
                              "SyN")
                              
  # apply the transforms to transfer template LMS to subject space. 
  pts = antsApplyTransformsToPoints( 3, ref.LMs, reg$fwdtransforms) 
  
  #write them in Slicer format.
  write.markups.fcsv(pts=pts, outfile=paste0("/mnt/c/temp/", f[X], ".fcsv"))

}

```
