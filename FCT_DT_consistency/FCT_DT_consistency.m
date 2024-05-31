% computing FCT_DT consistency 

% dti_dir:          the dti data path
% fct_dir :         the FCT data path
% mask_dir:         the mask file for the fmri data
% output_dir:       the output path



% Written by Qing Peng
% /2023/10/30


%% Read tensor and FA data of FCT and DTI
output_dir='\FCT_DT_consistency\';

dti_dir = '\dti_data';
subjects = dir(dti_dir); subjects = subjects(3:end); 

fct_dir = '\FCT_FA';
fct_images = dir(fct_dir); fct_images=fct_images(3:end);


mask_dir='\dti_mask';
files_mask = dir(mask_dir); 
addpath(genpath('help_functions'));

% Calculate the consistency of FCT and DTI within the mask
for j=1:1
    % Load data
    subjects_files = dir(fullfile(dti_dir,subjects(j).name));  subjects_files=subjects_files(3:end);
    disp(subjects_files(3).name)
    dti_tensor=fullfile(fullfile(dti_dir,subjects(j).name), subjects_files(3).name);
    dti_FA=fullfile(fullfile(dti_dir,subjects(j).name), subjects_files(2).name);
    fct_tensor=fullfile(fct_dir, fct_images(3*j).name);
    fct_FA=fullfile(fct_dir, fct_images(2+(j-1)*3).name);
    file_name2 = files_mask(j+2).name; % 获取文件名
    file_path_mask2 = fullfile(mask_dir, file_name2); % 构建完整的文件路径
    file_mask2 = dir(file_path_mask2);
    mask_DTI = fullfile(file_path_mask2, file_mask2(3).name); 
    mask_DTI_image=spm_read_vols(spm_vol(mask_DTI));
    dti_tensor_image=spm_read_vols(spm_vol(dti_tensor));
    dti_FA_image=spm_read_vols(spm_vol(dti_FA));
    fct_tensor_image=spm_read_vols(spm_vol(fct_tensor));
    fct_FA_image=spm_read_vols(spm_vol(fct_FA));
    mask_DTI_WM=find(fct_FA_image>0);

    
    % Normalization processing
    C=zeros(size(fct_FA_image));
    fct_wm_mask = fct_tensor_image;
    fct_sum = sum(fct_wm_mask,'all');
    fct_wm_mask2=zeros(size(fct_wm_mask));

    for k=1:6
        curr_volume_data = fct_wm_mask(:,:,:,k);          % choosing the current time-point (volume)
        curr_volume_data2 = fct_wm_mask2(:,:,:,k);  
        curr_volume_data2(mask_DTI_WM)=curr_volume_data(mask_DTI_WM)./fct_sum;
        fct_wm_mask2(:,:,:,k)=curr_volume_data2;
        clear current_volume_data curr_volume_data2;
    end

    max_fct = max(fct_wm_mask2(mask_DTI_WM));
    min_fct = min(fct_wm_mask2(mask_DTI_WM));
    scaling_factor2=max_fct-min_fct;

    dti_wm_mask = dti_tensor_image;
    dti_sum = sum(dti_wm_mask,'all');
    dti_wm_mask2=zeros(size(dti_wm_mask));

    for k=1:6
        curr_volume_data = dti_wm_mask(:,:,:,k);          % choosing the current time-point (volume)
        curr_volume_data2 = dti_wm_mask2(:,:,:,k);  
        curr_volume_data2(mask_DTI_WM)=curr_volume_data(mask_DTI_WM)./dti_sum;
        dti_wm_mask2(:,:,:,k)=curr_volume_data2;
        clear current_volume_data curr_volume_data2;
    end

    max_dti = max(dti_wm_mask2(mask_DTI_WM));
    min_dti = min(dti_wm_mask2(mask_DTI_WM));
    scaling_factor=max_dti-min_dti;

    fct_wm_mask21=zeros(size(fct_wm_mask2(:,:,:,1)));
    fct_wm_mask3=fct_wm_mask2(:,:,:,1);
    fct_wm_mask21(mask_DTI_WM)=(fct_wm_mask3(mask_DTI_WM)-min_fct)./scaling_factor2;

    dti_wm_mask21=zeros(size(dti_wm_mask2(:,:,:,1)));
    dti_wm_mask3=dti_wm_mask2(:,:,:,1);
    dti_wm_mask21(mask_DTI_WM)=(dti_wm_mask3(mask_DTI_WM)-min_dti)./scaling_factor;


    % FCT and DTI consistency were calculated in 6 vectors
    C_1 = (fct_wm_mask21 - dti_wm_mask21).^2; 
    C_all = zeros(size(C_1));
    C_all = C_all + C_1;

    for i=2:6
        if i==2 || i==3 || i==5
            fct_wm_mask21=zeros(size(fct_wm_mask(:,:,:,i)));
            fct_wm_mask3=fct_wm_mask2(:,:,:,i);
            fct_wm_mask21(mask_DTI_WM)=(fct_wm_mask3(mask_DTI_WM)-min_fct)/scaling_factor2;
    
            dti_wm_mask21=zeros(size(dti_wm_mask(:,:,:,i)));
            dti_wm_mask3=dti_wm_mask2(:,:,:,i);
            dti_wm_mask21(mask_DTI_WM)=(dti_wm_mask3(mask_DTI_WM)-min_dti)/scaling_factor;
    
            C_2 = (((fct_wm_mask21 - dti_wm_mask21)).^2); 
            C_all = C_all + C_2;

        else
            fct_wm_mask21=zeros(size(fct_wm_mask(:,:,:,i)));
            fct_wm_mask3=fct_wm_mask2(:,:,:,i);
            fct_wm_mask21(mask_DTI_WM)=(fct_wm_mask3(mask_DTI_WM)-min_fct)/scaling_factor2;
    
    
            dti_wm_mask21=zeros(size(dti_wm_mask(:,:,:,i)));
            dti_wm_mask3=dti_wm_mask2(:,:,:,i);
            dti_wm_mask21(mask_DTI_WM)=(dti_wm_mask3(mask_DTI_WM)-min_dti)/scaling_factor;
    
            C_2 = (fct_wm_mask21 - dti_wm_mask21).^2; 
            C_all = C_all + C_2;
        end
    end
    dti_FA_mask = zeros(size(dti_FA_image));
    dti_FA_mask(mask_DTI_WM)=dti_FA_image(mask_DTI_WM);

    fct_FA_mask = zeros(size(fct_FA_image));
    fct_FA_mask(mask_DTI_WM)=fct_FA_image(mask_DTI_WM);

    c_all=(1./((C_all).^(1/2)));           % FCT_DT C values 


    C =c_all.*fct_FA_mask.*dti_FA_mask;    % Weighting the FA values of the two modes

    sgw=[output_dir  subjects(j).name '_FCT_DTI.nii'];
    disp(sgw)
    y_Write(C,spm_vol([dti_tensor ',1']),sgw);   % saving the new functional matrix

end

