#deformable registration
path=getwd()
f=dir(patt="*InverseWarped.nii.gz")
f=sub("WarpedInverse","",f)
inputs=paste(path,f,sep="/")

outputs1=paste(path,"step2",f,sep="/")
outputs1=sub(".nii.gz","Warped.nii.gz",outputs1)
samples=sub("Warped.nii.gz","", outputs1)
outputs2=sub("Warped.nii.gz","InverseWarped.nii.gz",outputs1)

job_template=readLines(paste(root,"src/step2_template.txt",sep="/"))[1]

jobs=NULL
for (i in 1:length(f)){
  tmp=job_template
  tmp=sub('antsRegistration',paste(ANTSPATH,'antsRegistration',sep="/"),tmp)
  tmp=gsub('TEMPLATE',template,tmp)
  
  tmp=sub("OUTPUT",samples[i],tmp)
  tmp=sub("OUTPUTWarped.nii.gz",outputs1[i],tmp)
  tmp=sub("OUTPUTInverseWarped.nii.gz",outputs2[i],tmp)
  tmp=gsub("INPUT",inputs[i],tmp)
  jobs=c(jobs,tmp)
}
dir.create('step2')
cl2=makeCluster(4)
registerDoParallel(cl2)
foreach (i=1:length(jobs)) %dopar%  system(jobs[i])
stopCluster(cl2)
