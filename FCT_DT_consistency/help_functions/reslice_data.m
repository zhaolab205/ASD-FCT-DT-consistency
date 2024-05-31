function [output_image_matrix] = reslice_data(InputFile, TargetSpaceFile, interpolation_type)
% [output_image_matrix] = reslice_data(InputFile, TargetSpaceFile, interpolation_type)
%
% Adapted by Michael Peer, from Yan Chao-Gan 2015 y_Reslice.m function
% (part of DPARSFA)
%
% This function reslices InputFile to the dimensions of TargetSpaceFile,
% e.g. an anatomical image to the dimensions of a functional image.
%
% interpolation_type: 0 - Nearest Neighbour. 1: Trilinear. 
% (nearest neighbour should be used for labeled images where values must
% not change after interpolation)
%
%
%__________________________________________________________________________
% Written by YAN Chao-Gan 090302 for DPARSF. Referenced from spm_reslice.m in SPM5.
% State Key Laboratory of Cognitive Neuroscience and Learning, Beijing Normal University, China, 100875
% ycg.yan@gmail.com
%__________________________________________________________________________
% Revised by YAN Chao-Gan 100401. Fixed a bug while calculating the new dimension.
% Revised by YAN Chao-Gan 120229. Simplified the processing.
% Revised by YAN Chao-Gan 121214. Fixed the brain edge artifact when reslice to a bigger FOV. Apply a mask from the source image: don't extend values to outside brain.
% Last Revised by YAN Chao-Gan 151117. Fixed a bug when reslicing .nii.gz files. 


% read the target-space file and input file
RefHead = spm_vol([TargetSpaceFile ',1']);
mat=RefHead.mat;
dim=RefHead.dim;

SourceHead = spm_vol([InputFile ',1']);



%Handle .nii.gz. Referenced from y_Read.m. YAN Chao-Gan, 151117
[pathstr, name, ext] = fileparts(InputFile);
if strcmpi(ext,'.gz')
    gunzip(InputFile);
    FileNameWithoutGZ = fullfile(pathstr,name);
    SourceHead.fname = FileNameWithoutGZ;
end

[x1,x2,x3] = ndgrid(1:dim(1),1:dim(2),1:dim(3));
d     = [interpolation_type*[1 1 1]' [1 1 0]'];
C = spm_bsplinc(SourceHead, d);
v = zeros(dim);

M = inv(SourceHead.mat)*mat; % M = inv(mat\SourceHead.mat) in spm_reslice.m
y1   = M(1,1)*x1+M(1,2)*x2+(M(1,3)*x3+M(1,4));
y2   = M(2,1)*x1+M(2,2)*x2+(M(2,3)*x3+M(2,4));
y3   = M(3,1)*x1+M(3,2)*x2+(M(3,3)*x3+M(3,4));


output_image_matrix    = spm_bsplins(C, y1,y2,y3, d);

%Revised by YAN Chao-Gan 121214. Apply a mask from the source image: don't extend values to outside brain.
tiny = 5e-2; % From spm_vol_utils.c
Mask = true(size(y1));
Mask = Mask & (y1 >= (1-tiny) & y1 <= (SourceHead.dim(1)+tiny));
Mask = Mask & (y2 >= (1-tiny) & y2 <= (SourceHead.dim(2)+tiny));
Mask = Mask & (y3 >= (1-tiny) & y3 <= (SourceHead.dim(3)+tiny));

output_image_matrix(~Mask) = 0;


if strcmpi(ext,'.gz')
    delete(FileNameWithoutGZ);
end
