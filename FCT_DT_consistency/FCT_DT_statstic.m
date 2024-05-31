% Statstic FCT_DT consistency FA > 0.2 

% dti_fa:           the dti fa data path
% fct_dt_c:         the mask file for the fmri data
% mask_label:       the JHU ICBM-81-DTI WM atlas 
% fct_dt_w:         the FCT_DT consistency was normalized to the MNI space
% output_dir_FA2:   The FCT_DT data output path whose FA value is greater than 0.2 
    


% Written by Qing Peng
% /2023/10/30



%% FCT_DT consistency FA > 0.2 
dti_fa='\dti_fa_data';
dti_fa_images = dir(fullfile(dti_fa, '*.nii')); 

fct_dt_c='\FCT_DT_consistency';
fct_dt_c_images = dir(fullfile(fct_dt_c, '*.nii')); 

output_dir_FA2 = '\FCT_DT_c_FA2\';

addpath(genpath('help_functions'));
for j=1:1
     % Load the  fct_dt data and DT FA data
     dti_fa_image = dti_fa_images(j).name; 
     dti_fa_image_path = fullfile(dti_fa, dti_fa_image); 
     dti_FA=dti_fa_image_path;
     disp(dti_FA)
     dti_FA_image=spm_read_vols(spm_vol(dti_FA));

     fct_dt_c_images_name = fct_dt_c_images(j).name;
     fct_dt_c_images_path = fullfile(fct_dt_c, fct_dt_c_images_name); 
     dti_FA2=fct_dt_c_images_path;
     disp(dti_FA2)
     dti_MD_image=spm_read_vols(spm_vol(dti_FA2));

    % Define in FCT_DT the range in which the FA value in dt is greater than 0.2
     dti_FA_2=zeros(size(dti_FA_image));
     FA2 = find(dti_FA_image>0.2);
     dti_FA_2(FA2)=1;

     FCT_DT_C_FA2 =dti_MD_image.*dti_FA_2;
     sgw=[output_dir_FA2  dti_fa_image '_FA2.nii'];
     disp(sgw)
     y_Write(FCT_DT_C_FA2,spm_vol([dti_FA ',1']),sgw);   % saving the new functional matrix
end


%% Statstic FCT_DT consistency
mask_label='ICBM-81-DTI_WM_atlas.nii';
img_mask_label=spm_read_vols(spm_vol(mask_label));

fct_dt_w='\FCT_DT_w';
fct_dt_w_path = dir(fullfile(fct_dt_w, '*.nii')); 

addpath(genpath('help_functions'));

JBH_WM_FCT_DT_consistency = [];

for i=1:length(fct_dt_w_path)
    % Load the normalized fct_dt data 
     fct_dt_w_name = fct_dt_w_path(i).name; 
     fct_dt_w_name_path = fullfile(fct_dt_w, fct_dt_w_name); 
     img_FCT_DTI=spm_read_vols(spm_vol(fct_dt_w_name_path));
     JBH_WM_FCT_DT_c=[];
    % Calculate the FCT_DT_C mean in each WM tract
     for j=1:50
        WM_mask=find(img_mask_label==j);
        WM_FCT_DT=img_FCT_DT(WM_mask);
        WM_FCT_DT=WM_FCT_DT(~isnan(WM_FCT_DT));

        JBH_WM_FCT_DT_c=mean(WM_FCT_DT);
        JBH_WM_FCT_DT_c=[JBH_WM_FCT_DT_c,JBH_WM_FCT_DT_c];
     end
    JBH_WM_FCT_DT_consistency=[JBH_WM_FCT_DT_consistency;JBH_WM_FCT_DT_c];
end


