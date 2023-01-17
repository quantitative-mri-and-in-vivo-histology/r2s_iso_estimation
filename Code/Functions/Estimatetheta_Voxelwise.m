function Estimatetheta_Voxelwise(varargin)
% Varargin:
Porg = varargin{1}; % Original images
Prot = varargin{2}; % Coregister or aligned images
Vvector = varargin{3}; % V1 or main diffusion vector (depending on the model you are using)

Vvector_Vol = load_nii(Vvector);
VvectorImg = double(Vvector_Vol.img());
[pthname_theta, ~, ext_map] = fileparts(Vvector);

if nargin < 4 % In case you are using a mask, so it makes it easy.
    Mask_Applied = ones(size(VvectorImg,1:3));
else
    MaskAux = load_nii(varargin{4});
    Mask_Applied = double(MaskAux.img());
end

VvectorImg = VvectorImg .* Mask_Applied;

% Some variables to "speed up" the calculations:
RotatedB0 = zeros(3,size(Porg,1));
RotatedB0_inDiffSpace = zeros(3,size(Porg,1));
Theta_projection = zeros(2,size(Porg,1),...
                         size(VvectorImg,1),size(VvectorImg,2),size(VvectorImg,3));

% Check that the diffusion vectors are normalised. And for calculation, I just "reduced" dimensions to 2D.                    
VvectorR = reshape(double(VvectorImg),prod(size(VvectorImg,1:3)),size(VvectorImg,4));
VvectorR = VvectorR./sqrt(sum(VvectorR.*VvectorR,2));

% From DiffusionVector_Visualization -> x -> -x and z -> -z for the
% diffusion vectors, otherwise they are badly oriented.
% This needs to be checked manually! It depends on which scanner and which
% nii-convertor is used! In my case, it was necessary to flip the x and z
% axis, maybe you don't need it.
VvectorR(:,1) = -VvectorR(:,1);
VvectorR(:,3) = -VvectorR(:,3);

% Original B0 vector - from the original images (not rotated).
N1             = [0 0 1];

% Fixed transformation matrices:
% Now, these transformations are arbritary and estimated from the
% experimental registration from diff -> R2*. Therefore, here is done
% "inversely".
% This transformation matrix was obtained from the step 1 in Figure 2 of
% the paper.
Trans_Slic = [0.9897    -0.1400     -0.0300     0.0;...
             -0.0300    -0.0200     -0.9991     0.0;...
             0.1400     0.9900      -0.0300     0.0;...
             0.0        0.0         0.0         1.0];   
Swap_mX = [-1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];

param = spm_imatrix(inv(Trans_Slic));
Mrot_inv_TS  = spm_matrix([0 0 0 param(4) param(5) param(6)]);
warning off;
%figure; 
for inx = 1:size(Porg,1)       
    % compute angle theta
    % original volume
    Vorg = spm_vol(Porg{inx});       
    % rotated volume
    Vrot = spm_vol(Prot{inx});
    
    load(fullfile(fileparts(pthname_theta),'Chiasm1_allecho',num2str(inx),'TransformationMatrix.mat'));
    iMall(inx,:,:) = TM;
    
    param = spm_imatrix(squeeze(iMall(inx,:,:)));
    Mrot  = spm_matrix([0 0 0 param(4) param(5) param(6)]);
     
    %VvectorR = permute(double(Vvector_Vol.img()), [4,1,2,3]);
    % It seems that is the B0 field rotated while keeping the diffusion vectors!
    RotatedB0(1:3,inx) = Mrot(1:3,1:3) * N1(:);
    RotatedB0_inDiffSpace(:,inx) = (Swap_mX(1:3,1:3)\(Mrot_inv_TS(1:3,1:3) * RotatedB0(:,inx)));
    %RotatedB0_inDiffSpace(:,inx) = (Swap_mX\(Trans_Slic\RotatedB0(:,inx)));    
    RotatedB0_inDiffSpace(:,inx) = RotatedB0_inDiffSpace(:,inx)/norm(RotatedB0_inDiffSpace(:,inx));
    
    dot_prod = squeeze(VvectorR*RotatedB0_inDiffSpace(:,inx));
    Theta_projection(1,inx,:,:,:) = reshape(acos(dot_prod),size(Vvector_Vol.img(),1:3));

    % From visualisation, it seems that B0^z in diffusion space is -B0^z.
    %RotatedB0_inDiffSpace(3,inx) = -RotatedB0_inDiffSpace(3,inx);
    %dot_prod = squeeze(VvectorR*RotatedB0_inDiffSpace(:,inx));
    %Theta_projection(2,inx,:,:,:) = reshape(acos(dot_prod),size(Vvector_Vol.img(),1:3));
     
    figure(3);
    subplot(2,8,inx), imagesc(squeeze(Theta_projection(1,inx,:,:,26))*180/pi);
    subplot(2,8,inx), colormap('jet');
    subplot(2,8,inx), axis('off');
    %subplot(2,8,inx),caxis([0,180]);
    subplot(2,8,inx), caxis([0,pi]);
    subplot(2,8,inx), title({['B0 rot in R2*: ' num2str(squeeze(acos(RotatedB0(1:3,inx)'*[0 0 1]'))*180/pi)],...
                             ['B0 rot in diff: ' num2str(squeeze(acos(RotatedB0_inDiffSpace(:,inx)'*[0 1 0]'))*180/pi)]});
    disp(['B0 rot in R2* meas ' num2str(inx) ': ' num2str(squeeze(acos(RotatedB0(1:3,inx)'*[0 0 1]'))*180/pi)]);
%     ThetainB0atDiffSpace = squeeze(Theta_projection(1,inx,:,:,:));
%     Vvector_Vol.hdr.dime.dim = [3,83,64,50,1,1,1,1];
%     Vvector_Vol.hdr.dime.glmax = pi;
%     Vvector_Vol.hdr.dime.glmin = 0;
%     Vvector_Vol.img = ThetainB0atDiffSpace;
%     save_nii(Vvector_Vol,fullfile(fileparts(pthname_theta),'AngularOrientation_DiffSpace',['ThetatoB0Meas' num2str(inx) '_DiffSpace_Zdir' ext_map]));
%     
%     ThetaProjection_0toPiHalf = Theta_projection(1,inx,:,:,:);
%     ThetaProjection_0toPiHalf = abs(ThetaProjection_0toPiHalf);
%     ThetaProjection_0toPiHalf(ThetaProjection_0toPiHalf > pi/2) = pi - ThetaProjection_0toPiHalf(ThetaProjection_0toPiHalf > pi/2);
%     ThetaProjection_0toPiHalf = abs(ThetaProjection_0toPiHalf);
%     ThetainB0atDiffSpace = squeeze(ThetaProjection_0toPiHalf);
%     Vvector_Vol.hdr.dime.glmax = pi/2;
%     Vvector_Vol.hdr.dime.glmin = 0;
%     Vvector_Vol.img = ThetainB0atDiffSpace;
%     save_nii(Vvector_Vol,fullfile(fileparts(pthname_theta),'AngularOrientation_DiffSpace',['ThetatoB0Meas' num2str(inx) '_DiffSpace_Zdir_0toPiHalf' ext_map]));
end
    %ThetainB0atDiffSpace = squeeze(Theta_projection(2,inx,:,:,:));
    %Vvector_Vol.img = ThetainB0atDiffSpace;
    %save_nii(Vvector_Vol, [pthname_theta '\ThetatoB0Meas' num2str(inx) '_DiffSpace_minusZdir' ext_map]);
warning on
end

