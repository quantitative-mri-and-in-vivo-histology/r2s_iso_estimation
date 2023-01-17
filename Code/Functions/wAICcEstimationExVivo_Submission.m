function [mwAICc_for_table,sdwAICc_for_table] = ...
         wAICcEstimationExVivo_Submission(AICcValues,Angle_Meas,Kappa_Meas,ICVF_Meas,Mask,params)

     % Start of the code
    % Masking all the input parameters:
    % Mask the OC not only with the manual masking but also for unrealistic
    % values in \kappa and ICVF.
    FullOC_indx = (Mask == 1 & Kappa_Meas > 0 & ICVF_Meas > 0 & ICVF_Meas < 1);
    
    % Reshape the angular maps to leave it in n-total voxels x m-measurements.
    theta_Z = reshape(Angle_Meas,[prod(size(Angle_Meas,1:3)),size(Angle_Meas,4)]);
    theta_Z_ROI = theta_Z(FullOC_indx(:),:);
    clear 'theta_Z'
    
    % Check the angular values and modify accordingly: negative values are
    % changed, angles with > pi/2 are reflected (pi - angle value) and
    % unrealistic values (angle < 0.5° or pi/360 rads) are converted to 0.
    theta_Z_ROI(theta_Z_ROI < 0) = abs(theta_Z_ROI(theta_Z_ROI < 0));
    theta_Z_ROI(theta_Z_ROI > pi/2) = pi - theta_Z_ROI(theta_Z_ROI > pi/2);
    theta_Z_ROI(theta_Z_ROI < 0) = abs(theta_Z_ROI(theta_Z_ROI < 0));
    theta_Z_ROI(theta_Z_ROI < pi/360) = 0;
    % If an specific voxel has the same angular value across ALL the
    % angular measurements, then it is expected that these voxels are
    % irrelevant for this study. Those voxels are also used to filter the
    % other maps.
    remove_indx = all(theta_Z_ROI(:,:) == theta_Z_ROI(:,1),2);
    theta_Z_ROI(remove_indx,:) = [];   

    Kvalues = Kappa_Meas(FullOC_indx(:));
    Kvalues(remove_indx) = [];

    Ficvfvalues = ICVF_Meas(FullOC_indx(:));
    Ficvfvalues(remove_indx) = [];

    AICc_space = reshape(AICcValues,[prod(size(AICcValues,1:3)),size(AICcValues,4:6)]);
    % The AICc_space is converted to i-total voxels (x-y-z axis), j-angular
    % measurements (16), k-max. TE's (5) and l-models (4).
    AICc_space = AICc_space(FullOC_indx,:,:,:);
    AICc_space(remove_indx,:,:,:) = []; 
    Lin_AICc = reshape(AICc_space,[prod(size(AICc_space,1:3)),size(AICc_space,4)]);
    
    Diff_AICc = Lin_AICc - min(Lin_AICc')';
    clear 'Lin_AICc'
    wAICc_estimation = exp(-0.5.*Diff_AICc)./sum(exp(-0.5.*Diff_AICc),2);
    wAICc_space = reshape(wAICc_estimation,[prod(size(AICc_space,1:2)),size(AICc_space,3:4)]);
    
    % To generate the cumulative sampling, it is necessary to repeat the
    % ICVF and \kappa maps equal to the number of angular measurements.
    ficvf_fullsampling_Z = repmat(Ficvfvalues,16,1);
    k_values_Z = repmat(Kvalues,16,1);

    for k_indx = 1:3
        switch k_indx
            case 1
                Range_k = 0 < k_values_Z & k_values_Z < 1 & ficvf_fullsampling_Z > 0.8;
            case 2
                Range_k = 1 <= k_values_Z & k_values_Z < 2.5 & ficvf_fullsampling_Z > 0.8;
            case 3
                Range_k = 2.5 <= k_values_Z & ficvf_fullsampling_Z > 0.8;
        end
        
        mean_wAICc_space(1,k_indx,:,:) = mean(wAICc_space(Range_k,:,:));
        sd_wAICc_space(1,k_indx,:,:) = std(wAICc_space(Range_k,:,:));
    end
    
    % Resize it to fit into the ex vivo table
    mwAICc_for_table = repmat(mean_wAICc_space,[params.nbins,1,1,1]);
    sdwAICc_for_table = repmat(sd_wAICc_space,[params.nbins,1,1,1]);
end