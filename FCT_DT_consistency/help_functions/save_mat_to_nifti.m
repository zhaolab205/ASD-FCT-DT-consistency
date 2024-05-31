function save_mat_to_nifti(original_filename,mat,output_filename)
% save_mat_to_nifti(original_filename,mat,output_filename)
% 
% original filename is a nifti file with the same dimensions (needed to get
% all the nifti parameters instead of specifying them here) 

new_nifti_file=spm_vol(original_filename);
new_nifti_file.fname=output_filename; 
new_nifti_file.private.dat.fname =  output_filename;
spm_write_vol(new_nifti_file,mat);