% computing Functional Correlation Tensor(FCT) and FCT_FA

% dti_dir:          the dti data path
% func_dir:         the rs-fMRI data path
% mask_dir:         the mask file for the fmri data
% output_dir:       the output path
% out_name:         the prefix of output name,which will be added suffix(_tensor.nii)


% Written by Peng Qing
% /2023/10/30

% The detailed interpretation of method to calculate the functional tensor can be found in references as follows 

% Ding, Z., Newton, A. T., Xu, R., Anderson, A. W., Morgan, V. L., & Gore, J. C. (2013). Spatio-temporal correlation tensors reveal functional structure in human brain. PloS one, 8(12), e82107.
% Ding, Z., Xu, R., Bailey, S. K., Wu, T. L., Morgan, V. L., Cutting, L. E., ... & Gore, J. C. (2016). Visualizing functional pathways in the human brain using correlation tensors and magnetic resonance imaging. Magnetic resonance imaging, 34(1), 8-17.
% Zhao J, Huang C-C, Zhang Y, Liu Y, Tsai S-J, Lin C-P, et al. Structure-function coupling in white matter uncovers the abnormal brain connectivity in Schizophrenia. Transl Psychiatry. 2023;13:214


dti_dir='\dti_data';

mask_dir='\dti_mask\';
files_mask = dir(mask_dir); 

output_dir='\dti_data';

func_dir = '\fMRI_data';
func_images = dir(func_dir);
subjects = dir(func_dir); subjects = subjects(3:end); 

addpath(genpath('help_functions'));

%% FCT calculate
% Construct the unit 3X3 neighborhood matrix
for i=1:length(subjects)
    Vec = [];
    for z=-1:1
        for j = -1:1
            for k = -1:1
                if z==0&&j==0&&k==0
                    continue;
                end
                Vec = [Vec [z;j;k]];
            end
        end
    end
    
    dist = sqrt(sum(Vec.*Vec));
    uVec = Vec./repmat(dist,3,1);
    
    
    % design M for unit vector
    nneib=size(Vec,2);
    M=zeros(nneib,6);
    for t=1:nneib
        tmpx=uVec(1,t);
        tmpy=uVec(2,t);
        tmpz=uVec(3,t);
        tmpM=[tmpx^2,2*tmpx*tmpy,2*tmpx*tmpz,tmpy^2,2*tmpy*tmpz,tmpz^2];
        M(t,:)=tmpM;
    end

    M_transp=M';
    clear tmpx tmpy tmpz tmpM
    
     % Read DTI data path
     sub = func_images(i+2).name;
     file_path_rest = fullfile(dti_dir, sub); 
     file_rest = dir(fullfile(file_path_rest, '*.nii')); 
     rest_addr =[file_path_rest '\' file_rest(1).name];
     disp(rest_addr)

     % Load rs-fMRI data
     file_path=rest_addr;
     mr_nii=load_untouch_nii(file_path);
     mrData=mr_nii.img;
     mrData=single(mrData);
    
    % Load DTI mask data
    file_name2 = files_mask(i+2).name; 
    file_path_mask2 = fullfile(mask_dir, file_name2); 
    file_mask2 = dir(file_path_mask2);
    file_path_mask2 = fullfile(file_path_mask2, file_mask2(3).name); 
    disp(file_path_mask2)
    mask_nii2=load_untouch_nii(file_path_mask2);
    maskData2=mask_nii2.img;
        
    %Load the time series of each voxel in fMRI
    for j=1:size(mrData,4)                              % go over all timepoints (volumes)
        curr_volume_data = mrData(:,:,:,j);             % choosing the current time-point (volume)
        curr_volume_data(maskData2==0)=0;               % voxels are contained in the mask
        mask_rest_mat(:,:,:,j) = curr_volume_data;      % taking only the grey-matter from the current functional image
        clear current_volume_data;
    end
    mrData = mask_rest_mat;
    mrData(isnan(mask_rest_mat)==1)=0;mrData=squeeze(mrData);
    [nx,ny,nz,tc]=size(mrData);
    mr_tc=reshape(mrData,[],tc);
    clear mrData
    vox_ind=find(maskData2~=0);
    [vox_xyz(:,1),vox_xyz(:,2),vox_xyz(:,3)]=ind2sub(size(maskData2),vox_ind);  % get masked range
    

    nvox=size(vox_xyz,1);
    tensorMat=zeros(nvox,6);
    eigenVecMat1=zeros(nvox,3);
    fa=zeros(nvox,1);
%   md=zeros(nvox,1);
%   rd=zeros(nvox,1);
    ad=zeros(nvox,1);
    for x=1:nvox  %  parfor : parallel computing to decrease computation time
        r=vox_xyz(x,1);
        c=vox_xyz(x,2);
        d=vox_xyz(x,3);
        i_ind=vox_ind(x);
        if  sum(mr_tc(i_ind,:)==0)~=tc %make sure timeseries is not empty
            
            neibA_tc=zeros(tc,nneib); % stored the timeseries of vox and its 26 neibs (the 1st row was this voxel)
            A_tc=mr_tc(i_ind,:)';
            A=repmat([r;c;d],1,26);
            neibA=A+Vec;
            ind_zero=neibA==0;
            ind_out=(neibA(1,:)>nx | neibA(2,:)>ny | neibA(3,:)>nz);
            valid_ind=(sum(ind_zero)==0 & ind_out==0);
            
            Ind=sub2ind([nx,ny,nz],neibA(1,valid_ind),neibA(2,valid_ind),neibA(3,valid_ind));
            
            neibA_tc(:,valid_ind)=mr_tc(Ind,:)';
            neib_C=corr(A_tc,neibA_tc);
             
            
            neib_C(isnan(neib_C)) = 0;
            % neib_C=fisherz(neib_C);
            neib_C=atanh(neib_C);
            neib_C=neib_C.^2;
            tensortp=(inv(M_transp*M))*M_transp*neib_C'; 
            % remove the super large values
            [tmpx,~]=find(abs(tensortp)>20);

            tensortp(tmpx,:)=0;
            tensortp(isnan(tensortp)==1)=0;
            
            tensorMat(x,:)=tensortp;
            currMat=[tensortp(1),tensortp(2),tensortp(3);tensortp(2),tensortp(4),tensortp(5);tensortp(3),tensortp(5),tensortp(6)];  % 6 tensor values
            % FA value calculation(MD,AD,RD)
            [coeff,latent]=pcacov(currMat);
            eigenVecMat1(x,:)=coeff(:,1);
            meanEigen=mean(latent);
            fa(x)=sqrt((3*(((latent(1)-meanEigen)^2+(latent(2)-meanEigen)^2+(latent(3)-meanEigen)^2)/(latent(1)^2+latent(2)^2+latent(3)^2)))/2);
%           md(x)=(latent(1)+latent(2)+latent(3))/3;
%           ad(x)=latent(1);
%           rd(x)=(latent(2)+latent(3))/2;
        end
    end
    clear vox_xyz tc;
    
    % save data
    subjects_files = dir(fullfile(func_dir,subjects(i).name));  subjects_files=subjects_files(3:end);
    dti_FA=fullfile(fullfile(func_dir,subjects(i).name), subjects_files(1).name);
    disp(dti_FA)
    mask_path=dti_FA;
    mask_nii=load_untouch_nii(mask_path);

    Tensor_tensor=zeros(nx*ny*nz,6);
    Tensor_tensor(vox_ind,:)=tensorMat;
    Tensor_tensor=reshape(Tensor_tensor,nx,ny,nz,6);
    
%     Tensor_v1=zeros(nx*ny*nz,3);
%     Tensor_v1(vox_ind,:)=eigenVecMat1;
%     Tensor_v1=reshape(Tensor_v1,nx,ny,nz,3);
    
    Tensor_fa=zeros(nx*ny*nz,1);
    Tensor_fa(vox_ind)=fa;
    
%     Tensor_md=zeros(nx*ny*nz,1);
%     Tensor_md(vox_ind)=md;
% 
%     Tensor_ad=zeros(nx*ny*nz,1);
%     Tensor_ad(vox_ind)=ad;
% 
%     Tensor_rd=zeros(nx*ny*nz,1);
%     Tensor_rd(vox_ind)=rd;

             
            
%             filename=[out_name sub '_FCT_MD.nii'];
%             mask_nii.img=Tensor_md;
%             mask_nii.hdr.dime.datatype=16; % float 
%             mask_nii.hdr.dime.bitpix=32;
%             save_untouch_nii(mask_nii,filename)
% 
% 
%             filename=[out_name sub '_FCT_RD.nii'];
%             mask_nii.img=Tensor_rd;
%             mask_nii.hdr.dime.datatype=16; % float 
%             mask_nii.hdr.dime.bitpix=32;
%             save_untouch_nii(mask_nii,filename)
% 
% 
%             filename=[out_name sub '_FCT_AD.nii'];
%             mask_nii.img=Tensor_ad;
%             mask_nii.hdr.dime.datatype=16; % float 
%             mask_nii.hdr.dime.bitpix=32;
%             save_untouch_nii(mask_nii,filename)

            out_name='\FCT_DT_consistency\FCT_FA';
            filename=[out_name sub '_FCT_FA_rest.nii'];
            mask_nii.img=Tensor_fa;
            mask_nii.hdr.dime.datatype=16; % float 
            mask_nii.hdr.dime.bitpix=32;
            save_untouch_nii(mask_nii,filename)

            out_name='\FCT_DT_consistency\FCT_tensor';
            filename=[out_name sub '_FCT_tensor.nii'];
            mask_nii.img=Tensor_tensor;
            mask_nii.hdr.dime.dim(1)=4;
            mask_nii.hdr.dime.dim(5)=6;
            mask_nii.hdr.dime.datatype=16; % float 
            mask_nii.hdr.dime.bitpix=32;
            save_untouch_nii(mask_nii,filename)        
            
end   

