inVariant=function(moving,f,n.ang=5) {
  
  require(ANTsR)
  temp=tempdir()
  if (!dir.exists(temp)) dir.create(temp,recursive=T)
  print(moving)
  mfull=antsImageRead(moving)
  m = resampleImage( mfull %>% smoothImage(resample_resolution), rep(resample_resolution,3), interpType=0 )
  gc()
  print( f )
  nout = n.ang
  thetas = seq( from = 0, to = 360, length.out = nout )[1:(nout-1)]
  mytx = paste(temp,"/",moving,"_best.mat",sep="")
  mival<-invariantImageSimilarity( f, m,
                                   thetas=thetas, thetas2=thetas, thetas3=thetas,metric="GC",
                                   localSearchIterations=20, txfn=mytx,transform="Rigid" )
  print( mival )
  
  m2f = antsApplyTransforms( f, m, transformlist=mytx )
  www = antsApplyTransforms( ffull, mfull, transformlist=mytx )
  antsImageWrite(antsImageClone(www,out_pixeltype = "unsigned char"),paste("registered",moving,sep="/"))
  plot.antsImage(ffull,www,alpha=.5,axis=3,nslices=20,ncolumn=10,out=paste("cross_sections/cross_",moving,".jpg",sep=""))

}
