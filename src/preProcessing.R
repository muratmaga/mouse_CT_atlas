preProcessing=function(X,bone_low_thresh=20,bone_hi_thresh=255,dilate_radius=1,trim_padding=5) {
  img=antsImageRead(X)
  temp=thresholdImage(img,bone_low_thresh,bone_hi_thresh,1,0)
  temp=iMath(temp,"GetLargestComponent")
  temp=iMath(temp,"MD",dilate_radius)
  img=temp*img
  temp=cropImage(img,temp,label=1)
  temp=iMath(temp,"PadImage",trim_padding)
  antsSetOrigin(temp,c(0,0,0))
  antsImageWrite(antsImageClone(temp,out_pixeltype="unsigned char"),paste(getwd(),"masked",X,sep="/"))
  print(paste("wrote masked/",X,sep=""))
  gc()
}