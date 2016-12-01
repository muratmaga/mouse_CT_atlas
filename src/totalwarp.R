total_warp_template=paste(ANTSPATH,"antsApplyTransforms -d 3 -r TEMPLATE -o [aff_def/TOTALWARP.nii.gz,1] -t WARPED2WARP -t AFFINE",sep="/")
f=dir(patt="*Warped2Warp.nii.gz")
aff=dir(patt="*Affine.mat")
samples=sub("_Warped2Warp.nii.gz","",f)
dir.create("aff_def")

for (i in 1:length(f)) {
  job=sub("TEMPLATE", template,total_warp_template)
  job=sub("TOTALWARP",samples[i],job)
  job=sub("WARPED2WARP",f[i],job)
  job=sub("AFFINE",aff[i],job)
  system(job)
}

