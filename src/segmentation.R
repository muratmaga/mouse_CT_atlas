segmentation<-function(X){
  print(X)
  mi=antsImageRead(X)
  #mi=resampleImage(mi,c(antsGetSpacing(ffull)),0,1)      #not necessary for these downsampled volumes
  mytx=antsRegistration(fixed = ffull,moving = mi,"SyN")
  mywarped_template=antsApplyTransforms(fixed = mi ,moving=segs,transformlist = mytx$invtransforms,interpolator = "nearestNeighbor")
 
   for (s in 1:length(segments)) {
    if (s!=4) {
      mask=thresholdImage(mywarped_template,lothresh = segments[s], hithresh = segments[s],inval = 1,outval = 0)
      mask=iMath(mask,operation = "MD",1)
      masked=maskImage(mi,mask)
      masked=cropImage(masked,mask)
      antsImageWrite(masked,paste(segment_labels[s],"/",segment_labels[s],"_",X,sep=""))
    } else {
      mask=thresholdImage(mywarped_template,lothresh = segments[s], hithresh = segments[s],inval = 1,outval = 0)
      mask=iMath(mask,operation = "MD",1)
      masked=cropImage(mask,mask)
      antsImageWrite(masked,paste(segment_labels[s],"/",segment_labels[s],"_",X,sep=""))
      
    }
  }
}
