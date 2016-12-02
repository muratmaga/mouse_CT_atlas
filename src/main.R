MM=proc.time()
#prior to invoking R, make sure TEMP environmental variable points to a folder with a large amount of space available.
#Parts of the workflow use parallel processing (doParallel library) to speed up the execution 
# image registrations require large amount of RAM (64+ GB) to be available and will output 17GB of data.

require(ANTsR)
require(doParallel)

root="/scratch/mouse_CT_atlas/" #point it to the base of the checked out directory
ANTSPATH='/home/apps/ants/bin'  #specify where the ANTs image registration library executables are located. This will be used during the shape analysis. 

template=paste(root,"/data/templates/low_res_template_UCHAR.nii.gz",sep="") #use this template to work with the provided data.
segs=paste(root,"/data/templates/low_res_template_UCHAR-label.nii.gz",sep="")

source(paste(root,"src/preProcessing.R",sep="/")) #to remove non-craniofacial elements
source(paste(root,"src/inVariant.R",sep='/'))  #to rigidly register cleaned NIFTIs to the template
source(paste(root,"src/segmentation.R",sep="/")) #to isolate skull and hemi-mandibles

###########
#PART I. preprocessing to remove non-cranial elements from the field of view. 
###########
out_folder=paste(root,"/data/targets/masked/registered",sep="/")
if (!dir.exists(out_folder)) dir.create(out_folder,recursive=T)

#image processing settings
bone_low=20      #this works with sample scans provided. You may have to work this out by trial and error for your own datasets. 
bone_hi=255
dilate_radius=1 #How much to dilate the mask, in voxels. 
padding=5 #padding to provide in voxels. 

setwd(paste(root,"/data/targets",sep=''))
targets=dir(patt="nii.gz")
for (i in targets) preProcessing(i,bone_low_thresh=bone_low,bone_hi_thresh=bone_hi,dilate_radius=dilate_radius,trim_padding=padding)
#outputs are saved in the masked folder

##########
#PART II. Rigidly registering to the template to remove positional differences in scans by using the multistart function. 
##########
setwd(paste(root,"/data/targets/masked",sep=""))

if (!dir.exists("registered")) dir.create("registered")     #this will where the rigidly registered output volumes be saved
if (!dir.exists("cross-sections")) dir.create("cross_sections") #cross-sectional views superimposed on template for quick QC checks

targets=dir(patt="gz")
print(targets)

ffull=antsImageRead(template)
resample_resolution=0.5 #resample to very coarse resolution to align the sample 

f = resampleImage( ffull %>% smoothImage(resample_resolution), rep(resample_resolution,3), interpType=0 )   

#  start the search with 90 degree intervals (n.ang=90)
n.ang=5  

for (i in targets) inVariant(i,f,n.ang)
#output NIFTIs are in the registered  folder.
#you can check the slice views to evaluate the outcome of the rigid registrations in cross-sections folder
#for the failed ones increase the n.ang to 9. Note the number of iterations are to the cube of n.angle (i.e. 9^3=729 registrations per sample)   
#e.g. sample C57L_J needs more iterations
inVariant(targets[grep("C57L_J",targets)],f,9)
   
  
############
#PART III.    Segmenting individual CF structures.
############
#template label maps are 1: SKull, 2: Right Mandible, 21: Left Mandible, 15: Endocranial space
segments=c(1,2,21,15)
segment_labels=c("Skull","R_Mand","L_Mand","Endo")
setwd(paste(root,"/data/targets/masked/registered",sep=""))
targets=dir(patt="*nii.gz")
segs=antsImageRead(segs)

for (i in 1:4) dir.create(segment_labels[i])
for (i in targets) segmentation(i)
#segmented elements are in following folders
#Skull
#L_Mand
#R_Mand
#Endo


#####END OF IMAGE PROCESSING


#PART IV. Registrations for Image Based Shape Analysis (diGPA)

#PART IVA. Remove the nuisance parameters through rigid + uniform scaling
#need to switch to skull only template to calculate the deformations 
template=paste(root,"/data/templates/low_res_skull_only.nii.gz",sep="")
setwd(paste(root,"/data/targets/masked/registered/Skull",sep=""))
targets=dir(patt="*nii.gz")

source(paste(root,"src/job_builder_for_step1.R",sep=""))    #this will invoke the antsRegistration from the shell.

#PART IVB. Now register the output of step 1 to the template via affine + SyN 
setwd(paste(root,"/data/targets/masked/registered/Skull/step1",sep=""))

source(paste(root,"src/job_builder_for_step2.R",sep=""))      #this will invoke the antsRegistration from the shell.

#PART V. Calculate the combined AFF+SyN deformation field (aka total warp)
setwd(paste(root,"/data/targets/masked/registered/Skull/step1/step2",sep=""))     #this will invoke the antsApplyTransforms from the shell.
source(paste(root,"/src/totalwarp.R",sep=''))

#IMAGE PROCESSING FINISHED!
print(proc.time()-MM)




