# BATCH-JD
Python script that allows for Batch Processing and Calculation of Dice Overlap and Jaccard Index.

**FolderDirectories and File Types:** 


BATCH JD requires files that are being compared to be in these following formats: 
  (.nii.gz)
  (.mgz)
  (.nii)
  (.nrrd)
  (.nii.seg.nrrd)
As of current state of the program (.mrb) files are not supported, therefore files that are in (.mrb) formats can be modified and exported via programs such as 3D SLICER. 

Files that are being compared must be in separate folders containing files of the same type. For example, if Dice and Jaccard calculations are executed on files such as 100307_aseg.hires.mgz and 100307_label.nii.gz, both files should be contained in separate folders, along with other file of the same type that require calculations. 

Although the script reorganizes file names in numerical order, it is recommended to make sure that files in comparing folders are organized in such a way that file pairs are in the same position in their corresponding folder. 

**Segmentation label Entries**


The order in which the folder directories are entered is arbitrary, however the entries for segmentation labels must correspond to the entry of the folder directory. For instance, if Directory 1 contains files that have a segment label that is '26' for the L - Accumbens, then the user must enter '26' in the entry for segmentation label 1. 

**Script output**


The script will write 2 different files for each file read with the extension of  <>_pmod.nii, <>_pmod.nii.seg.nrrd or <>_smod.nii, <>_smod.nii.seg.nrrd, these files are mainly created so that the program can read the segmentation labels, however, the nii.seg.nrrd files can be imported into 3D slicer in the case where the user does not know the segmentation label. 
The 3D graphs displayed only refer to the first iteration of pairwise files, and displays the binary mask of the selected segmentation label in both files. Therefore, the binary masks of other iterations will not be displayed. The graphs are solely used for visualization purposes to make sure the user has selected the correct segmentation label. 

**#File prep**


If segmentation is done using 3D slicer, it is recommended to export segmentations as segmentation files rather than volumes to accomodate for computation speed. 
