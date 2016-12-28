# Mus musculus craniofacial atlas and image based shape analysis repository
Clone this repository somewhere with sufficiently large space. 
When executed, sample workflow produces 17GB of output. 
You need to have ANTs and ANTsR installed in your system:

  Instructions for installing ANTS: https://brianavants.wordpress.com/2012/04/13/updated-ants-compile-instructions-april-12-2012/
  Instruction for installing ANTsR: https://github.com/stnava/ANTsR#installation-from-source (method #2 works well).
  
After both of them are installed, edit the first two lines of src/main.R to specify the root directory where you cloned the repository and where ants/bin is installed on the system. You can then execute the image processing pipeline as

nohup Rscript src/main.R > processing.log&

Image segmentation pipeline (Part I-III) takes approximately 20-30 minutes depending on the availability of computing power. Registrations (Part IV) will take 5-10h, again depending on the available computational power. 

Reported results and graphs in Maga et al. can be reproduced by running the src/analysis.R script interactively. 






