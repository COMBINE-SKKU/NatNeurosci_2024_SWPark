
%% -- Preparation ------------------------------------------
    %% 01) load toolboxes and define paths
    for pathtoolbox = 1

        %download the follwoing tools
        % surfstat
        % brainspace
        % cifti-matlab
        % BCT 
        % HCPPipelines
        
        WB_COMMAND = [ 'wb_command' ];
        addpath(genpath('/local_raid1/01_software/toolboxes/surfstat/'));   %SurfStatView
        addpath(genpath('/local_raid1/01_software/toolboxes/npy-matlab/')); %upload hcp_colormap
        addpath(genpath('/local_raid1/01_software/toolboxes/matlab_util/')); %BoSurfStatView
        addpath(genpath('/local_raid1/01_software/toolboxes/cifti-matlab/')); %used for gifti (uploading surfaces)
        addpath(genpath('./files'));
        addpath(genpath('/local_raid2/03_user/shinwon/03_software/BrainSpace/matlab'));


    end
    %% 02) read phenotypic info
    for for_readdata=1
        data = importfile_NKI_subses(['./files/NKI_sub_ses.xlsx']);        
    end    
    for for_demo_vars = 1

        subid = data.subject;
        session = data.session;

        idx_ses1 = find( (session == 'ses-v1') );
        idx_ses2 = find( (session == 'ses-v2') );

        subid=cellstr(subid)
        session=cellstr(session)
        
    end    
    %% 03) Load surfaces & parcellations & colormaps

    for for_load_surfaces = 1

        %upload 10k surfaces: fsLR space         
        temp = gifti(['./files/S900.R.midthickness_MSMAll.10k_fs_LR.surf.gii']);
        surfR.coord = temp.vertices';
        surfR.tri   = temp.faces;   

        temp = gifti(['./files/S900.L.midthickness_MSMAll.10k_fs_LR.surf.gii']);
        surfL.coord = temp.vertices';
        surfL.tri   = temp.faces;

        surf.coord = [ surfL.coord surfR.coord ];
        surf.tri   = [ surfL.tri; surfR.tri+10242; ];

        temp = gifti(['./files/L.atlasroi.10k_fs_LR.shape.gii']);
        surfL_roi = temp.cdata';     
        temp = gifti(['./files/R.atlasroi.10k_fs_LR.shape.gii']);
        surfR_roi = temp.cdata'; 

        surf_roi = cat(2,surfL_roi, surfR_roi);
        figure; SurfStatView(surf_roi,surf);        

        for parcel_10k = 1
            for load_schaefers_400_parcel = 1
                schaefer_400 = ciftiopen(['./files/Schaefer2018_400Parcels_7Networks_order_10k.dlabel.nii'],  WB_COMMAND); 
                schaefer_400_label = schaefer_400.cdata;  
                schaefer_400_full = zeros(20482,1); 
                schaefer_400_full(logical(surf_roi)) = schaefer_400_label;
               
                figure; SurfStatView(schaefer_400_full,surf);
                schaefer_400_label = schaefer_400_full;                      
            end      
            for load_yeo_7networks_cortex = 1
               cii_yeo = ciftiopen( ['./files/RSN-networks.10k_fs_LR.dlabel.nii'],  WB_COMMAND);
               cii_yeo_label = cii_yeo.cdata;       

                yeo_7_label = cii_yeo_label(:,1);
                yeo_7_label(yeo_7_label==37) = 0;
                yeo_7_label(yeo_7_label==38) = 3;
                yeo_7_label(yeo_7_label==39) = 6;
                yeo_7_label(yeo_7_label==40) = 7;
                yeo_7_label(yeo_7_label==41) = 1;
                yeo_7_label(yeo_7_label==42) = 5;
                yeo_7_label(yeo_7_label==43) = 2;
                yeo_7_label(yeo_7_label==44) = 4;      

                yeo_17_label = cii_yeo_label(:,2);
                yeo_17_label(yeo_17_label==37) = 0;   
                yeo_17_label(yeo_17_label==45) = 2;   
                yeo_17_label(yeo_17_label==46) = 16;   
                yeo_17_label(yeo_17_label==47) = 6;   
                yeo_17_label(yeo_17_label==48) = 12;   
                yeo_17_label(yeo_17_label==49) = 10;   
                yeo_17_label(yeo_17_label==50) = 5;   
                yeo_17_label(yeo_17_label==51) = 13;   
                yeo_17_label(yeo_17_label==52) = 17;   
                yeo_17_label(yeo_17_label==53) = 15;   
                yeo_17_label(yeo_17_label==54) = 1;   
                yeo_17_label(yeo_17_label==55) = 4;   
                yeo_17_label(yeo_17_label==56) = 9;   
                yeo_17_label(yeo_17_label==57) = 11;   
                yeo_17_label(yeo_17_label==58) = 3;   
                yeo_17_label(yeo_17_label==59) = 8;   
                yeo_17_label(yeo_17_label==60) = 7;   
                yeo_17_label(yeo_17_label==61) = 14;   

                yeo_7_label_full = zeros(1,20484); 
                yeo_7_label_full(logical(surf_roi)) = yeo_7_label'
                figure; SurfStatView(yeo_7_label_full,surf);

            end    
        end   

    end    
    for for_colormap_load=1
    videenmap = videen(20); videenmap(1:19,:)
    videenmap = [ videenmap; 0.7 0.7 0.7 ];
    hcp_colormap = readNPY([ './files/hcp_colormap.npy' ]);
% 
%     mycol.blackblue =   flipud( [zeros(1,3)*0.8; zeros(127,1) ...
%         (0:126)'/127 ones(127,1)]);
%     mycol.blue      =   flipud( [ones(1,3)*0.8; ...
%         zeros(127,1) (0:126)'/127 ones(127,1)]);
%     mycol.blue2      =   flipud( [ones(1,3); ...
%         zeros(127,1) (0:126)'/127 ones(127,1)]);    
%     mycol.red       =   [ones(1,3)*0.8; ...
%         ones(500,1) linspace(0,253,500)'/254 zeros(500,1);...
%         ones(64,1) ones(64,1) linspace(0,200,64)'/254];
%     mycol.red2       =   [ones(1,3); ...
%         ones(500,1) linspace(0,253,500)'/254 zeros(500,1);...
%         ones(64,1) ones(64,1) linspace(0,200,64)'/254];
    yeo_colormap = [ 200 200 200;    
                 120 18 134;
                 70 130 180;
                 0 118 14;
                 196 58 250;
                 220 248 164;
                 230 148 34;
                 205 62 78 ]/255;
    yeo_colormap_1=yeo_colormap(2:end,:)
    % yeo_colormap= yeo_colormap(2:8, :)
   

%
    values = [...
        'CC'; '10'; '33'; %_rbgyr20_01
        '99'; '20'; '66'; %_rbgyr20_02
        '66'; '31'; '99'; %_rbgyr20_03
        '34'; '41'; 'CC'; %_rbgyr20_04
        '00'; '51'; 'FF'; %_rbgyr20_05
        '00'; '74'; 'CC'; %_rbgyr20_06
        '00'; '97'; '99'; %_rbgyr20_07
        '00'; 'B9'; '66'; %_rbgyr20_08
        '00'; 'DC'; '33'; %_rbgyr20_09
        '00'; 'FF'; '00'; %_rbgyr20_10
        '33'; 'FF'; '00'; %_rbgyr20_11
        '66'; 'FF'; '00'; %_rbgyr20_12
        '99'; 'FF'; '00'; %_rbgyr20_13
        'CC'; 'FF'; '00'; %_rbgyr20_14
        'FF'; 'FF'; '00'; %_rbgyr20_15
        'FF'; 'CC'; '00'; %_rbgyr20_16
        'FF'; '99'; '00'; %_rbgyr20_17
        'FF'; '66'; '00'; %_rbgyr20_18
        'FF'; '33'; '00'; %_rbgyr20_19
        'FF'; '00'; '00';]; %_rbgyr20_20

    values = reshape(hex2dec(values), [3 numel(values)/6])' ./ 255;

    P = size(values,1); 
    rbgyr20_colormap = interp1(1:size(values,1), values, linspace(1,P,200), 'linear');    
    
    
     values2 = [ 158,1,66;    
                 213,62,79;
                 244,109,67;
                 253,174,97;
                 254,224,139; 
                 255,255,191;
                 230,245,152;
                 171,221,164;
                 102,194,165;
                 50,136,189;
                 94,79,162;
                 ]/255;   
 
    P = size(values2,1); 
    spectral_colormap = interp1(1:size(values2,1), values2, linspace(1,P,200), 'linear');      
    
    
    %addpath '/local_raid2/03_user/shinwon/01_project/01_thalamusgradient/05_analysis_cmap_pmap/customcolormap'

    %pasteljet=customcolormap_preset('pasteljet')
    %rdylbl=customcolormap_preset('red-yellow-blue')
    %rdwhbl=customcolormap_preset('red-white-blue')


    end

%% Extended Fig 1A: correlation between NKI and HCPD: CMAP
    %% upload HCPD group CMAP
     for for_grouptemplate_HCPD_10k =1

        cii = ciftiopen([ './files/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.cmap.dscalar.nii' ], WB_COMMAND);   
        cii_L_thal = cii.cdata(cii.diminfo{1, 1}.models{20, 1}.start:(cii.diminfo{1,1}.models{20,1}.start-1) + cii.diminfo{1,1}.models{20,1}.count,:);
            %timeseries_L_thal(:,:,i) = cii_L_thal;        
        groupGradient_L_hcpd(:, :) = cii_L_thal;

        cii = ciftiopen([ './files/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.49.cmap.dscalar.nii' ], WB_COMMAND);   
        cii_R_thal = cii.cdata(cii.diminfo{1,1}.models{21,1}.start:(cii.diminfo{1,1}.models{21,1}.start-1) + cii.diminfo{1,1}.models{21,1}.count,:);
            %timeseries_L_thal(:,:,i) = cii_L_thal;        
        groupGradient_R_hcpd(:, :) = cii_R_thal;
        
    end
    %% upload NKI group CMAP 
    cii = ciftiopen([ './files/NKI_atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.cmap.dscalar.nii' ], WB_COMMAND);   
    cii_L_thal = cii.cdata(cii.diminfo{1, 1}.models{20, 1}.start:(cii.diminfo{1,1}.models{20,1}.start-1) + cii.diminfo{1,1}.models{20,1}.count,:);
    groupGradient_L_nki(:, :) = cii_L_thal;

    cii = ciftiopen([ './files/NKI_atlasroi.Neo-Subcortical-2mm.10k_fs_LR.49.cmap.dscalar.nii' ], WB_COMMAND);   
    cii_R_thal = cii.cdata(cii.diminfo{1,1}.models{21,1}.start:(cii.diminfo{1,1}.models{21,1}.start-1) + cii.diminfo{1,1}.models{21,1}.count,:);
    groupGradient_R_nki(:, :) = cii_R_thal;
    %% correlation bt NKI & HCPD CMAPs
    groupGradient2_L_nki_flip=groupGradient_L_nki(:,2)*-1
    [R P]= corr(groupGradient2_L_nki_flip(:,1), groupGradient_L_hcpd(:,1)) % R=0.83
    [R P]= corr(groupGradient_L_nki(:,3), groupGradient_L_hcpd(:,2)) %R=0.71

    figure; scatter(groupGradient2_L_nki_flip, groupGradient_L_hcpd(:,1), 5, 'filled', 'MarkerFaceColor', [0.7, 0.7, 0.7]); 
    h=lsline; h.LineWidth=1; h.Color='red'; set(gcf, 'color', 'w');set(gca,'FontSize',14);

    figure; scatter(groupGradient_L_nki(:,3), groupGradient_L_hcpd(:,2), 5, 'filled', 'MarkerFaceColor', [0.7, 0.7, 0.7]); 
    h=lsline; h.LineWidth=1; h.Color='red'; set(gcf, 'color', 'w');set(gca,'FontSize',14);

%% Extended Fig 1B: correlation between NKI and HCPD: NEOMAP 
    %% read subjects
        filename = './files/NKI_subjects_any_session_QC_cpac.txt';
        data1 = readtable(filename, 'Delimiter', ',', 'Format', '%s%s%s', 'ReadVariableNames', false );
        data1=data1(:,1:2)
        data1.subject = strcat('sub-', data1.Var1);
        data1.session = strcat('ses-', data1.Var2);
        data1=data1(:,3:4)
        [~, index] = ismember(data1, data);
        % Display indices
        disp('Indices:');
        disp(index);
    %% NKI group NEOMAPs
    %% upload source data for NKI 
    
    load('sourceData_ExtFig1_NKI_NEOMAPs.mat')
   
        for for_preparation_source_data = 1 %do not run, upload source data 
            %% NEOMAP template: NKI  
            for read_individual_pmaps_from_groupCmap_to_create_template = 1
    %             % extracted individual pmaps using the group cmap template with python code
    %             PATH2               = '/local_raid1/03_user/shinwon/01_project/00_thalamus/01_nki_replication/';
    %             IND_PMAP_NKI_GROUPCMAP = [PATH2 '05_aligned_pmaps/03_groupCmap_qc_cpac/' ];

                for for_nki_ses12=1
                    for i = 1:length(subid)
    % 
    %                     pmap_L = ciftiopen([ IND_PMAP_NKI_GROUPCMAP subid{i} '/' session{i} '/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.10.pmap.dscalar.nii' ], WB_COMMAND);   
    %                     pmap_L_thal = pmap_L.cdata(1:(pmap_L.diminfo{1,1}.models{1,1}.count) + pmap_L.diminfo{1,1}.models{2,1}.count,:);
    % 
    %                     pmap_L_full = zeros(pmap_L.diminfo{1,1}.models{1,1}.numvert+pmap_L.diminfo{1,1}.models{2,1}.numvert, pmap_L.diminfo{1,2}.length);
    %                     pmap_L_full(logical( surf_roi ),:) = pmap_L_thal;
    % 
    %                     pmap_L_thal_ind_groupCmap_nki(i, :, :) = pmap_L_full;   
    % 
    %                     pmap_R = ciftiopen([ IND_PMAP_NKI_GROUPCMAP subid{i} '/' session{i} '/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.49.49.pmap.dscalar.nii' ], WB_COMMAND);   
    %                     pmap_R_thal = pmap_R.cdata(1:(pmap_R.diminfo{1,1}.models{1,1}.count) + pmap_R.diminfo{1,1}.models{2,1}.count,:);
    % 
    %                     pmap_R_full = zeros(pmap_R.diminfo{1,1}.models{1,1}.numvert+pmap_R.diminfo{1,1}.models{2,1}.numvert, pmap_R.diminfo{1,2}.length);
    %                     pmap_R_full(logical( surf_roi ),:) = pmap_R_thal;
    % 
    %                     pmap_R_thal_ind_groupCmap_nki(i, :, :) = pmap_R_full;    

                    end    
                end

            end
            for for_create_group_PMAP_ses1_qc_cpac = 1
    %             %use mean pmap of all subjects as template 
    %             group_pmap_R_thal_ind_groupCmap_nki = nanmean(pmap_R_thal_ind_groupCmap_nki(index,:,:),1); 
    %             group_pmap_R_thal_ind_groupCmap_nki_rs = reshape(group_pmap_R_thal_ind_groupCmap_nki, [20484,10]);
    % 
    %             group_pmap_R_thal_ind_groupCmap_nki_rs_nz = group_pmap_R_thal_ind_groupCmap_nki_rs;
    %             temp_idx=find(group_pmap_R_thal_ind_groupCmap_nki_rs_nz(:,1)==0);
    %             group_pmap_R_thal_ind_groupCmap_nki_rs_nz(temp_idx,:)=[];
    %             group_pmap_R_thal_ind_groupCmap_nki_rs_nz=normalize(group_pmap_R_thal_ind_groupCmap_nki_rs_nz);
    % 
    %             group_pmap_L_thal_ind_groupCmap_nki = nanmean(pmap_L_thal_ind_groupCmap_nki(index,:,:),1);
    %             group_pmap_L_thal_ind_groupCmap_nki_rs = reshape(group_pmap_L_thal_ind_groupCmap_nki, [20484,10]);
    % 
    %             group_pmap_L_thal_ind_groupCmap_nki_rs_nz = group_pmap_L_thal_ind_groupCmap_nki_rs;
    %             temp_idx=find(group_pmap_L_thal_ind_groupCmap_nki_rs_nz(:,1)==0);
    %             group_pmap_L_thal_ind_groupCmap_nki_rs_nz(temp_idx,:)=[];
    %             group_pmap_L_thal_ind_groupCmap_nki_rs_nz=normalize(group_pmap_L_thal_ind_groupCmap_nki_rs_nz);

            end    
            %% upload individual NEOMAPs - Session 1 + align to template 
    % 
    %             %align individual pmap to groupPmap: ses1
    % 
    %             pmap_L_thal_ind_groupCmap_nki_ses1_align = zeros(length(idx_ses1), 20484, 10);
    %             pmap_R_thal_ind_groupCmap_nki_ses1_align = zeros(length(idx_ses1), 20484, 10);
    % 
    %             for i = 1:length(idx_ses1)
    % 
    %                 pmap_L = ciftiopen([ IND_PMAP_NKI_GROUPCMAP subid{idx_ses1(i)} '/' session{idx_ses1(i)} '/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.10.pmap.dscalar.nii' ], WB_COMMAND);   
    %                 pmap_L_thal = pmap_L.cdata(1:(pmap_L.diminfo{1,1}.models{1,1}.count) + pmap_L.diminfo{1,1}.models{2,1}.count,:);
    % 
    %                 pmap_R = ciftiopen([ IND_PMAP_NKI_GROUPCMAP subid{idx_ses1(i)} '/' session{idx_ses1(i)} '/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.49.49.pmap.dscalar.nii' ], WB_COMMAND);   
    %                 pmap_R_thal = pmap_R.cdata(1:(pmap_R.diminfo{1,1}.models{1,1}.count) + pmap_R.diminfo{1,1}.models{2,1}.count,:);
    % 
    %                 %procrustes alignment 
    %                 [d_L,Z_L,transform_L] = procrustes(group_pmap_L_thal_ind_groupCmap_nki_rs_nz,pmap_L_thal); %18722
    %                 [d_R,Z_R,transform_R] = procrustes(group_pmap_L_thal_ind_groupCmap_nki_rs_nz,pmap_R_thal);
    % 
    %                 temp_L = zeros(20484,10);
    %                 temp_L(logical(surf_roi),:) = Z_L;
    %                 temp_R = zeros(20484,10);
    %                 temp_R(logical(surf_roi),:) = Z_R;
    % 
    %                 pmap_L_thal_ind_groupCmap_nki_ses1_align(i, :, :) = temp_L;
    %                 pmap_R_thal_ind_groupCmap_nki_ses1_align(i, :, :) = temp_R;
    % 
    %                 pmap_L_thal_ind_groupCmap_nki_ses1_align_nz(i, :, :) = Z_L;
    %                 pmap_R_thal_ind_groupCmap_nki_ses1_align_nz(i, :, :) = Z_R;
    %             end   
    % 
    %             %normalized gradients
    %             pmap1_L_thal_ind_groupCmap_nki_ses1_align = zeros(length(idx_ses1),20484);
    %             pmap2_L_thal_ind_groupCmap_nki_ses1_align = zeros(length(idx_ses1),20484);
    %             pmap3_L_thal_ind_groupCmap_nki_ses1_align = zeros(length(idx_ses1),20484);
    % 
    %             pmap1_R_thal_ind_groupCmap_nki_ses1_align = zeros(length(idx_ses1),20484);
    %             pmap2_R_thal_ind_groupCmap_nki_ses1_align = zeros(length(idx_ses1),20484);
    %             pmap3_R_thal_ind_groupCmap_nki_ses1_align = zeros(length(idx_ses1),20484);
    % 
    %             for i= 1:length(idx_ses1)
    %                 temp_pmap1 =  squeeze(pmap_L_thal_ind_groupCmap_nki_ses1_align_nz(i,:,1));  
    %                 temp1a = ((temp_pmap1 - min(temp_pmap1)) / (max(temp_pmap1) - min(temp_pmap1))) * (1 - (-1)) + (-1);
    %                 temp1b = normalize(temp_pmap1);        
    %                 pmap1_L_thal_ind_groupCmap_nki_ses1_align(i,logical(surf_roi)) = temp1b;
    %                 pmap1_L_thal_ind_groupCmap_nki_ses1_align_scale(i,logical(surf_roi)) = temp1a;
    % 
    %                 temp_pmap2 =  squeeze(pmap_L_thal_ind_groupCmap_nki_ses1_align_nz(i,:,2));  
    %                 temp2a = ((temp_pmap2 - min(temp_pmap2)) / (max(temp_pmap2) - min(temp_pmap2))) * (1 - (-1)) + (-1);
    %                 temp2b = normalize(temp_pmap2);
    %                 pmap2_L_thal_ind_groupCmap_nki_ses1_align(i,logical(surf_roi)) = temp2b;
    %                 pmap2_L_thal_ind_groupCmap_nki_ses1_align_scale(i,logical(surf_roi)) = temp2a;
    % 
    %                 temp_pmap3 =  squeeze(pmap_L_thal_ind_groupCmap_nki_ses1_align_nz(i,:,3));            
    %                 temp3a = ((temp_pmap3 - min(temp_pmap3)) / (max(temp_pmap3) - min(temp_pmap3))) * (1 - (-1)) + (-1);
    %                 temp3b = normalize(temp_pmap3);
    %                 pmap3_L_thal_ind_groupCmap_nki_ses1_align(i,logical(surf_roi)) = temp3b;
    %                 pmap3_L_thal_ind_groupCmap_nki_ses1_align_scale(i,logical(surf_roi)) = temp3a;
    % 
    %             end
            %% upload individual NEOMAPs - Session 2 + align to template 
% 
%         %align individual pmap to groupPmap: ses2
%         pmap_L_thal_ind_groupCmap_nki_ses2_align = zeros(length(idx_ses2), 20484, 10);
%         pmap_R_thal_ind_groupCmap_nki_ses2_align = zeros(length(idx_ses2), 20484, 10);
% 
%         parfor i = 1:length(idx_ses2)
% 
%             pmap_L = ciftiopen([ IND_PMAP_NKI_GROUPCMAP subid{idx_ses2(i)} '/' session{idx_ses2(i)} '/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.10.pmap.dscalar.nii' ], WB_COMMAND);   
%             pmap_L_thal = pmap_L.cdata(1:(pmap_L.diminfo{1,1}.models{1,1}.count) + pmap_L.diminfo{1,1}.models{2,1}.count,:);
% 
%             pmap_R = ciftiopen([ IND_PMAP_NKI_GROUPCMAP subid{idx_ses2(i)} '/' session{idx_ses2(i)} '/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.49.49.pmap.dscalar.nii' ], WB_COMMAND);   
%             pmap_R_thal = pmap_R.cdata(1:(pmap_R.diminfo{1,1}.models{1,1}.count) + pmap_R.diminfo{1,1}.models{2,1}.count,:);
% 
%             %procrustes alignment 
%             [d_L,Z_L,transform_L] = procrustes(group_pmap_L_thal_ind_groupCmap_nki_rs_nz,pmap_L_thal); %18722
%             [d_R,Z_R,transform_R] = procrustes(group_pmap_L_thal_ind_groupCmap_nki_rs_nz,pmap_R_thal);
% 
%             temp_L = zeros(20484,10);
%             temp_L(logical(surf_roi),:) = Z_L;
%             temp_R = zeros(20484,10);
%             temp_R(logical(surf_roi),:) = Z_R;
% 
%             pmap_L_thal_ind_groupCmap_nki_ses2_align(i, :, :) = temp_L;
%             pmap_R_thal_ind_groupCmap_nki_ses2_align(i, :, :) = temp_R;
% 
%             pmap_L_thal_ind_groupCmap_nki_ses2_align_nz(i, :, :) = Z_L;
%             pmap_R_thal_ind_groupCmap_nki_ses2_align_nz(i, :, :) = Z_R;
%         end   
% 
%         %normalized gradients
%         pmap1_L_thal_ind_groupCmap_nki_ses2_align = zeros(length(idx_ses2),20484);
%         pmap2_L_thal_ind_groupCmap_nki_ses2_align = zeros(length(idx_ses2),20484);
%         pmap3_L_thal_ind_groupCmap_nki_ses2_align = zeros(length(idx_ses2),20484);
% 
%         pmap1_R_thal_ind_groupCmap_nki_ses2_align = zeros(length(idx_ses2),20484);
%         pmap2_R_thal_ind_groupCmap_nki_ses2_align = zeros(length(idx_ses2),20484);
%         pmap3_R_thal_ind_groupCmap_nki_ses2_align = zeros(length(idx_ses2),20484);
% 
%         for i= 1:length(idx_ses2)
%             temp_pmap1 =  squeeze(pmap_L_thal_ind_groupCmap_nki_ses2_align_nz(i,:,1));  
%             temp1a = ((temp_pmap1 - min(temp_pmap1)) / (max(temp_pmap1) - min(temp_pmap1))) * (1 - (-1)) + (-1);
%             temp1b = normalize(temp_pmap1);        
%             pmap1_L_thal_ind_groupCmap_nki_ses2_align(i,logical(surf_roi)) = temp1b;
%             pmap1_L_thal_ind_groupCmap_nki_ses2_align_scale(i,logical(surf_roi)) = temp1a;
% 
%             temp_pmap2 =  squeeze(pmap_L_thal_ind_groupCmap_nki_ses2_align_nz(i,:,2));  
%             temp2a = ((temp_pmap2 - min(temp_pmap2)) / (max(temp_pmap2) - min(temp_pmap2))) * (1 - (-1)) + (-1);
%             temp2b = normalize(temp_pmap2);
%             pmap2_L_thal_ind_groupCmap_nki_ses2_align(i,logical(surf_roi)) = temp2b;
%             pmap2_L_thal_ind_groupCmap_nki_ses2_align_scale(i,logical(surf_roi)) = temp2a;
% 
%             temp_pmap3 =  squeeze(pmap_L_thal_ind_groupCmap_nki_ses2_align_nz(i,:,3));            
%             temp3a = ((temp_pmap3 - min(temp_pmap3)) / (max(temp_pmap3) - min(temp_pmap3))) * (1 - (-1)) + (-1);
%             temp3b = normalize(temp_pmap3);
%             pmap3_L_thal_ind_groupCmap_nki_ses2_align(i,logical(surf_roi)) = temp3b;
%             pmap3_L_thal_ind_groupCmap_nki_ses2_align_scale(i,logical(surf_roi)) = temp3a;
% 
%         end
% 
%         figure; histogram(pmap2_L_thal_ind_groupCmap_nki_ses1_align)
%         figure; histogram(pmap2_L_thal_ind_groupCmap_nki_ses2_align)
% 
%         figure; histogram(pmap3_L_thal_ind_groupCmap_nki_ses1_align)
%         figure; histogram(pmap3_L_thal_ind_groupCmap_nki_ses2_align)
% 
%         figure; BoSurfStatView(mean(pmap2_L_thal_ind_groupCmap_nki_ses1_align,1),surf)
%         figure; BoSurfStatView(mean(pmap2_L_thal_ind_groupCmap_nki_ses2_align,1),surf)
% 
%         figure; BoSurfStatView(mean(pmap3_L_thal_ind_groupCmap_nki_ses1_align,1),surf)
%         figure; BoSurfStatView(mean(pmap3_L_thal_ind_groupCmap_nki_ses2_align,1),surf)
        end
        
    %% upload source data for HCPD NEOMAPs
    
    load('sourceData_Fig2,3_ExtFig3,4_HCPD_NEOMAPs.mat')
    
    pmap1_R_thal_ind_groupCmap_hcpd_align = pmap_R_thal_ind_groupCmap_hcpd_align(:,:,1);
    pmap2_R_thal_ind_groupCmap_hcpd_align = pmap_R_thal_ind_groupCmap_hcpd_align(:,:,2);

    pmap1_L_thal_ind_groupCmap_hcpd_align = pmap_L_thal_ind_groupCmap_hcpd_align(:,:,1);
    pmap2_L_thal_ind_groupCmap_hcpd_align = pmap_L_thal_ind_groupCmap_hcpd_align(:,:,2);

    %% scatter pmap1 - pamp2: NEOMAP space - Extended Figure 1B (SuppFig3)
  
    for allAges=1
        % PMAP data                 
        x = -mean(pmap2_L_thal_ind_groupCmap_nki_ses1_align,1); 
        y = mean(pmap3_L_thal_ind_groupCmap_nki_ses1_align,1); 
        
        groups = yeo_7_label_full;

        x(x==0)=NaN;
        y(y==0)=NaN;

        % Scatter plot with different colors for each group
        figure('Position', [100, 100, 800, 800], 'Color', 'white'); % Adjust the values as needed
        scatter(x, y, 50, groups, 'filled', 'MarkerFaceAlpha', 0.1);colormap(yeo_colormap);
        set(gcf,'color','w');
        set(gca, 'FontSize', 16); 
        xlabel('NEOMAP 1', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('NEOMAP 2', 'FontSize', 18, 'FontWeight', 'bold');

        % Marking centroids
        uniqueGroups = unique(groups); 
        hold on;
        for i = 1:numel(uniqueGroups)
            group = uniqueGroups(i);
            groupPoints = [x(groups == group)', y(groups == group)'];
            centroid = mean(groupPoints);        
            plot(centroid(1), centroid(2), 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', [yeo_colormap(i,:)], 'MarkerEdgeColor', 'w' );
            centroid_networks_age1{i} = centroid;
            groupPoints_networks_age1{i} = groupPoints;

        end

        % visualize external networks' centroids   
        % vis som dan sal lim fpn dmn 
        centroid_networks_age1=centroid_networks_age1(:,2:8);
        external_3_age1 =cell2mat([centroid_networks_age1(1,1);centroid_networks_age1(1,2);centroid_networks_age1(1,3)])
        centroid_external_age1 = mean(external_3_age1);
        plot(centroid_external_age1(1), centroid_external_age1(2), 'Marker', '*', 'MarkerSize', 10, 'LineWidth', 2, 'Color', 'w');
        hold off;

    end
    %% correlation bt HCPD & NKI: NEOMAP 

    %neomap 1
    pmap1_L_thal_ind_groupCmap_nki_align_mean = nanmean(pmap1_L_thal_ind_groupCmap_nki_ses1_align,1);
    pmap2_L_thal_ind_groupCmap_nki_align_mean = nanmean(pmap2_L_thal_ind_groupCmap_nki_ses1_align,1);
    pmap1_L_thal_ind_groupCmap_hcpd_align_mean = nanmean(pmap1_L_thal_ind_groupCmap_hcpd_align,1);

    figure; BoSurfStatView(-pmap2_L_thal_ind_groupCmap_nki_align_mean,surf);colormap(spectral_colormap);%BoSurfStatColLim([-200,200]);
    figure; BoSurfStatView(pmap1_L_thal_ind_groupCmap_hcpd_align_mean,surf);colormap(spectral_colormap);%BoSurfStatColLim([-200,200]);
    
    [R P] = corr(-pmap2_L_thal_ind_groupCmap_nki_align_mean', pmap1_L_thal_ind_groupCmap_hcpd_align_mean')

    % Extended Figure 1B - top scatter plot
    x=-pmap2_L_thal_ind_groupCmap_nki_align_mean;
    y=pmap1_L_thal_ind_groupCmap_hcpd_align_mean;
    figure; scatter(x, y, 10, 'k', 'filled','MarkerEdgeColor','w', 'MarkerFaceAlpha', 0.5);
    xlabel('NKI','FontSize',16); ylabel('HCPD','FontSize',16); 
    hold on; %h= lsline; set(h,'Color','r'); 
    p=polyfit(x,y,1);
    yfit=polyval(p,x);
    plot(x,yfit,'-r','LineWidth',1);
    set(gcf,'Color','w');
    ax=gca; ax.FontSize=14;

    % neomap 2
   
    pmap3_L_thal_ind_groupCmap_nki_align_mean = nanmean(pmap3_L_thal_ind_groupCmap_nki_ses1_align,1);
    pmap2_L_thal_ind_groupCmap_hcpd_align_mean = nanmean(pmap2_L_thal_ind_groupCmap_hcpd_align,1);

    figure; BoSurfStatView(-pmap3_L_thal_ind_groupCmap_nki_align_mean,surf);colormap(spectral_colormap);%BoSurfStatColLim([-200,200]);
    figure; BoSurfStatView(pmap2_L_thal_ind_groupCmap_hcpd_align_mean,surf);colormap(spectral_colormap);%BoSurfStatColLim([-200,200]);

    [R P] = corr(pmap3_L_thal_ind_groupCmap_nki_align_mean', pmap2_L_thal_ind_groupCmap_hcpd_align_mean')

    % Extended Figure 1B - bottom scatter plot
    x=pmap3_L_thal_ind_groupCmap_nki_align_mean;
    y=pmap2_L_thal_ind_groupCmap_hcpd_align_mean;
    figure; scatter(x, y, 20, 'k', 'filled','MarkerEdgeColor','w', 'MarkerFaceAlpha', 0.5);
    xlabel('NKI','FontSize',12); ylabel('HCPD','FontSize',12); 
    hold on; %h= lsline; set(h,'Color','r'); 
    p=polyfit(x,y,1);
    yfit=polyval(p,x);
    plot(x,yfit,'-r','LineWidth',1);
    set(gcf,'Color','w');
    ax=gca; ax.FontSize=14;
   
%% Extended Figure 6: NKI Longitudinal analysis 
    %% read demo file with age information 
    data_nki = importfile_demo_NKI_long(['./files/demo_nki_long.xlsx']);

    age_v1 = data_nki.age_bas1
    age_v2 = data_nki.age_flu1
    sex = data_nki.sex  

%% 1) Mixed effects model 

    for for_prepare_y_segidx = 1
       %% segregation index
        idx_yeo_vis = find(yeo_7_label_full == 1);
        idx_yeo_som = find(yeo_7_label_full == 2);
        idx_yeo_dan = find(yeo_7_label_full == 3);
        idx_yeo_sal = find(yeo_7_label_full == 4);
        idx_yeo_lim = find(yeo_7_label_full == 5);
        idx_yeo_fpn = find(yeo_7_label_full == 6);
        idx_yeo_dmn = find(yeo_7_label_full == 7);

       % baseline salience-external
       pmap2_L_thal_ses1_sal = pmap2_L_thal_ind_groupCmap_nki_ses1_align(:, idx_yeo_sal);
       pmap2_L_thal_ses1_sal_mean = mean(pmap2_L_thal_ses1_sal,2);
       pmap2_L_thal_ses1_ext = pmap2_L_thal_ind_groupCmap_nki_ses1_align(:, [ idx_yeo_vis idx_yeo_som idx_yeo_dan ] );
       pmap2_L_thal_ses1_ext_mean = mean(pmap2_L_thal_ses1_ext,2);
       pmap2_L_thal_ses1_sal_ext = abs(pmap2_L_thal_ses1_sal_mean - pmap2_L_thal_ses1_ext_mean);

       % baseline salience-internal; 
       pmap3_L_thal_ses1_sal = pmap3_L_thal_ind_groupCmap_nki_ses1_align(:, idx_yeo_sal);
       pmap3_L_thal_ses1_sal_mean = mean(pmap3_L_thal_ses1_sal,2);
       pmap3_L_thal_ses1_int = pmap3_L_thal_ind_groupCmap_nki_ses1_align(:, idx_yeo_dmn );
       pmap3_L_thal_ses1_int_mean = mean(pmap3_L_thal_ses1_int,2);
       pmap3_L_thal_ses1_sal_int = abs(pmap3_L_thal_ses1_sal_mean - pmap3_L_thal_ses1_int_mean);
       
       % follow-up salience-external
       pmap2_L_thal_ses2_sal = pmap2_L_thal_ind_groupCmap_nki_ses2_align(:, idx_yeo_sal);
       pmap2_L_thal_ses2_sal_mean = mean(pmap2_L_thal_ses2_sal,2);
       pmap2_L_thal_ses2_ext = pmap2_L_thal_ind_groupCmap_nki_ses2_align(:, [ idx_yeo_vis idx_yeo_som idx_yeo_dan ] );
       pmap2_L_thal_ses2_ext_mean = mean(pmap2_L_thal_ses2_ext,2);
       pmap2_L_thal_ses2_sal_ext = abs(pmap2_L_thal_ses2_sal_mean - pmap2_L_thal_ses2_ext_mean);
       
       % follow-up salience-internal
       pmap3_L_thal_ses2_sal = pmap3_L_thal_ind_groupCmap_nki_ses2_align(:, idx_yeo_sal);
       pmap3_L_thal_ses2_sal_mean = mean(pmap3_L_thal_ses2_sal,2);
       pmap3_L_thal_ses2_int = pmap3_L_thal_ind_groupCmap_nki_ses2_align(:, idx_yeo_dmn );
       pmap3_L_thal_ses2_int_mean = mean(pmap3_L_thal_ses2_int,2);
       pmap3_L_thal_ses2_sal_int = abs(pmap3_L_thal_ses2_sal_mean - pmap3_L_thal_ses2_int_mean);
       
       % baseline external - internal ; 
       
       pmap2_L_thal_ses1_ext = pmap2_L_thal_ind_groupCmap_nki_ses1_align(:,  [ idx_yeo_vis idx_yeo_som idx_yeo_dan ]  );
       pmap2_L_thal_ses1_ext_mean = mean(pmap2_L_thal_ses1_ext,2);
       
       pmap3_L_thal_ses1_int = pmap3_L_thal_ind_groupCmap_nki_ses1_align(:, idx_yeo_dmn );
       pmap3_L_thal_ses1_int_mean = mean(pmap3_L_thal_ses1_int,2);
       
       ses1_ext_int = abs(pmap2_L_thal_ses1_ext_mean - pmap3_L_thal_ses1_int_mean)
       
       % follow-up external - internal ; 
       
       pmap2_L_thal_ses2_ext = pmap2_L_thal_ind_groupCmap_nki_ses2_align(:,  [ idx_yeo_vis idx_yeo_som idx_yeo_dan ]  );
       pmap2_L_thal_ses2_ext_mean = mean(pmap2_L_thal_ses2_ext,2);
       
       pmap3_L_thal_ses2_int = pmap3_L_thal_ind_groupCmap_nki_ses2_align(:, idx_yeo_dmn );
       pmap3_L_thal_ses2_int_mean = mean(pmap3_L_thal_ses2_int,2);
       
       ses2_ext_int = abs(pmap2_L_thal_ses2_ext_mean - pmap3_L_thal_ses2_int_mean)
       %% prepare longitudinal graph 

       data_graph = [data_nki.id data_nki.age_bas1 data_nki.age_flu1 pmap2_L_thal_ses1_sal_ext  pmap2_L_thal_ses2_sal_ext pmap3_L_thal_ses1_sal_int pmap3_L_thal_ses2_sal_int]
       id_long = [data_nki.id; data_nki.id]
       age_long = [data_nki.age_bas1; data_nki.age_flu1]
       sal_ext = [pmap2_L_thal_ses1_sal_ext; pmap2_L_thal_ses2_sal_ext]
       sal_int = [pmap3_L_thal_ses1_sal_int; pmap3_L_thal_ses2_sal_int]
       ext_int = [ses1_ext_int; ses2_ext_int]
       timepoint = [ones(1,70), ones(1,70)*2]'
       data_graph_long = table(id_long,age_long, sal_ext, sal_int, ext_int, timepoint)

    end    
    for for_test_segidx = 1 
        % Define data table here
        SubjectID=id_long;
        SubjectAge=age_long;
        SalienceExternal=sal_ext;
        SalienceInternal=sal_int;
        ExternalInternal=ext_int;
        TimePoint=timepoint;

        data_fig = table(SubjectID, SubjectAge, SalienceExternal, SalienceInternal, ExternalInternal, TimePoint);

        %Supp FIGURE 10A
        for for_fig_Ext_Int = 1

            % Initialize figure
            figure; hold on

            % --- Plot for salience-external ---
            title('External-Internal');

            % Find unique subject IDs
            subjects = unique(data_fig.SubjectID);

            % Loop through each subject to plot their data for salience-external
            for i = 1:length(subjects)
                subjData = data_fig(data_fig.SubjectID == subjects(i), :);
                % Plot the points for each visit for salience-external
                plot(subjData.SubjectAge, subjData.ExternalInternal, 'o-'); % Connects the two points with a line
            end

            % Fit a global linear model for all data points for salience-external
            lm_extint = fitlm(data_fig.SubjectAge, data_fig.ExternalInternal);

            % Get the predicted values and confidence intervals for a range of ages
            ageRange = linspace(min(data_fig.SubjectAge), max(data_fig.SubjectAge), 100)';
            [predExtInt ciExtInt] = predict(lm_extint, ageRange);

            % Plot the global fitted line for salience-external
            plot(ageRange, predExtInt, 'k-', 'LineWidth', 2);

            % Plot the confidence interval for salience-external
            fill([ageRange; flipud(ageRange)], [ciExtInt(:,1); flipud(ciExtInt(:,2))], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

            % Labels for salience-external
            xlabel('Age');
            ylabel('Segregation Index');

            ax=gca; 
            ax.FontSize=14;
            set(gcf,'color','w')
            ylim([0 1.6]);
            set(gcf,'Position',[500 500 500 500])


        end
        for for_fig_Sal_Ext_Int_separate = 1

            % Initialize figure
            figure; hold on;

            % --- Plot for salience-external ---
            title('Salience-External');

            % Find unique subject IDs
            subjects = unique(data_fig.SubjectID);

            % Loop through each subject to plot their data for salience-external
            for i = 1:length(subjects)
                subjData = data_fig(data_fig.SubjectID == subjects(i), :);
                % Plot the points for each visit for salience-external
                plot(subjData.SubjectAge, subjData.SalienceExternal, 'o-'); % Connects the two points with a line
            end

            % Fit a global linear model for all data points for salience-external
            lm_external = fitlm(data_fig.SubjectAge, data_fig.SalienceExternal)

            % Get the predicted values and confidence intervals for a range of ages
            ageRange = linspace(min(data_fig.SubjectAge), max(data_fig.SubjectAge), 100)';
            [predExternal, ciExternal] = predict(lm_external, ageRange, 'Alpha',0.05);

            % Plot the global fitted line for salience-external
            plot(ageRange, predExternal, 'k-', 'LineWidth', 2);

            % Plot the confidence interval for salience-external
            fill([ageRange; flipud(ageRange)], [ciExternal(:,1); flipud(ciExternal(:,2))], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

            % Labels for salience-external
            xlabel('Age','FontSize',14);
            ylabel('Segregation Index','FontSize',26);

            ax=gca; 
            ax.FontSize=14;
            ylim([0 1.6]);
            hold off;
            set(gcf,'color','w')
            set(gcf,'Position',[500 500 500 500])

            % --- Plot for salience-internal ---
            figure; hold on;
            title('Salience-Internal');

            % Loop through each subject to plot their data for salience-internal
            for i = 1:length(subjects)
                subjData = data_fig(data_fig.SubjectID == subjects(i), :);

                % Plot the points for each visit for salience-internal
                plot(subjData.SubjectAge, subjData.SalienceInternal, 'o-'); % Connects the two points with a line
            end

            % Fit a global linear model for all data points for salience-internal
            lm_internal = fitlm(data_fig.SubjectAge, data_fig.SalienceInternal)

            % Get the predicted values and confidence intervals for a range of ages
            [predInternal, ciInternal] = predict(lm_internal, ageRange);

            % Plot the global fitted line for salience-internal
            plot(ageRange, predInternal, 'k-', 'LineWidth', 2);

            % Plot the confidence interval for salience-internal
            fill([ageRange; flipud(ageRange)], [ciInternal(:,1); flipud(ciInternal(:,2))], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

            % Labels for salience-internal
            xlabel('Age','FontSize',26);
            ylabel('Segregation Index','FontSize',26);

            ax=gca; 
            ax.FontSize=14;
            ylim([0 1.6]);

            hold off; 
            set(gcf,'color','w')
            set(gcf,'Position',[500 500 500 500])


        end 
       
        %Supp FIGURE 10A - statistics
        for for_test = 1
         
            data_fig.SubjectID = categorical(data_fig.SubjectID);

            % Fit the mixed-effects model for Salience-External
            lmeExternal = fitlme(data_fig, 'SalienceExternal ~ SubjectAge + (1|SubjectID)');
            disp(lmeExternal.Coefficients);

            % Test the significance of the age effect on Salience-External
            [pValueExternal,~,statsExternal] = coefTest(lmeExternal);
            fprintf('P-value for age effect on Salience-External: %3.3f\n', pValueExternal);

            % Fit the mixed-effects model for Salience-Internal
            lmeInternal = fitlme(data_fig, 'SalienceInternal ~ SubjectAge + (1|SubjectID)');
            disp(lmeInternal.Coefficients);

            % Test the significance of the age effect on External-Internal
            [pValueInternal,~,statsInternal] = coefTest(lmeInternal);
            fprintf('P-value for age effect on Salience-Internal: %3.3f\n', pValueInternal);

            % Fit the mixed-effects model for External-Internal
            lmeExtInt = fitlme(data_fig, 'ExternalInternal ~ SubjectAge + (1|SubjectID)');
            disp(lmeExtInt.Coefficients);

            % Test the significance of the age effect on External-Internal
            [pValueExtInt,~,statsExtInt] = coefTest(lmeExtInt);
            fprintf('P-value for age effect on External-Internal: %3.3f\n', pValueExtInt);


        end 
           
    end
    
%% 2) Predictive modeling
    %% upload timeseries

    load('sourceData_ExtFig1_NKI_FC.mat')
    
    for for_prepare_source_data = 1 % do not run, upload source data
        % NKI_TS = '/local_raid1/03_user/shinwon/01_project/00_thalamus/01_nki_replication/02_downsampled10k/'
        %session 2
        for i = 1 : length(idx_ses2)
%             %upload dtseries
%             cii = ciftiopen( [ NKI_TS subid{idx_ses2(i)} '/ses-v2/ses-v2_task-rest_acq-645_Atlas_s0_10k.dtseries.nii'], WB_COMMAND );
%             ts_temp = cii.cdata';
% 
%             surf_roi_idx = find(surf_roi);
% 
%             ts = zeros(size(ts_temp, 1), size(surf_roi, 2));
%             ts(:, surf_roi_idx) = ts_temp(:, 1:length(surf_roi_idx));
%             ts = ts(11:end,:)'; %take out first 10 TR
%             TS_ses2{i} = ts'; %take out first 10 TR
% 
%             TS_Lthal_ses2{i} = cii.cdata(cii.diminfo{1, 1}.models{20, 1}.start:cii.diminfo{1, 1}.models{20, 1}.start+cii.diminfo{1, 1}.models{20, 1}.count - 1 , 11:end); %left thalamus
%             TS_Rthal_ses2{i} = cii.cdata(cii.diminfo{1, 1}.models{21, 1}.start:cii.diminfo{1, 1}.models{21, 1}.start+cii.diminfo{1, 1}.models{21, 1}.count - 1 , 11:end); %left thalamus


            for for_fc_matrix = 1 

                for for_parcelbyparcel = 1
%                     %schaefer400    
%                     ts_clus1 = zeros(size(ts,2), max(schaefer_400_label)); 
%                     for j = 1: max(schaefer_400_label)
%                         ts_clus1(:,j) = mean(ts(schaefer_400_label == j,:),1); 
%                     end
% 
%                     TS_clus.sch400_ses2{i} = ts_clus1;
% 
%                     r2 = corr(double(ts_clus1));
%                     z2 = 1/2*log((1+r2)./(1-r2));    
%                     fc_sch400_nki_ses2_int{i} = int16(z2*10000);  
%                     fc_sch400_nki_ses2{i} = z2;      

                end

            end

        end
        %session 1
        for i = 1 : length(idx_ses1)
            %upload dtseries
%             cii = ciftiopen( [ NKI_TS subid{idx_ses1(i)} '/ses-v1/ses-v1_task-rest_acq-645_Atlas_s0_10k.dtseries.nii'], WB_COMMAND );
%             ts_temp = cii.cdata';
%             %TR{i} = cii.diminfo{1,2}.seriesStep;
% 
%             surf_roi_idx = find(surf_roi);
% 
%             ts = zeros(size(ts_temp, 1), size(surf_roi, 2));
%             ts(:, surf_roi_idx) = ts_temp(:, 1:length(surf_roi_idx));
%             ts = ts(11:end,:)'; %take out first 10 TR
%             TS_ses1{i} = ts'; %take out first 10 TR
% 
%             %ts_sc = ts_temp(11:end,cii.diminfo{1, 1}.models{3, 1}.start:end)'; %take out first 10 TR
%             %TS_SC{i} = ts_temp(11:end,cii.diminfo{1, 1}.models{3, 1}.start:end)'; %take out first 10 TR
% 
%             TS_Lthal_ses1{i} = cii.cdata(cii.diminfo{1, 1}.models{20, 1}.start:cii.diminfo{1, 1}.models{20, 1}.start+cii.diminfo{1, 1}.models{20, 1}.count - 1 , 11:end); %left thalamus
%             TS_Rthal_ses1{i} = cii.cdata(cii.diminfo{1, 1}.models{21, 1}.start:cii.diminfo{1, 1}.models{21, 1}.start+cii.diminfo{1, 1}.models{21, 1}.count - 1 , 11:end); %left thalamus


            for for_fc_matrix = 1 
                for for_parcelbyparcel = 1
%                 %schaefer400    
%                 ts_clus1 = zeros(size(ts,2), max(schaefer_400_label)); 
%                 for j = 1: max(schaefer_400_label)
%                     ts_clus1(:,j) = mean(ts(schaefer_400_label == j,:),1); 
%                 end
% 
%                 TS_clus.sch400_ses1{i} = ts_clus1;
% 
%                 r2 = corr(double(ts_clus1));
%                 z2 = 1/2*log((1+r2)./(1-r2));    
%                 fc_sch400_nki_ses1_int{i} = int16(z2*10000);  
%                 fc_sch400_nki_ses1{i} = z2;      
                end
            end

        end    
    
    end
    
    %% derive individual gradients aligned to group gradient 
        %functional connectivity

        for i = 1 : length(idx_ses2)
            i
            fc_sch400_nki_ses2_a(i,:,:) = fc_sch400_nki_ses2{1, i}';

        end        

        conn_matrix1 = mean(fc_sch400_nki_ses2_a,1);
        conn_matrix1=squeeze(conn_matrix1);
        conn_matrix1(conn_matrix1==Inf) = 1;

        %Gref = GradientMaps('kernel', 'na', 'approach', 'le');
        Gref = GradientMaps('kernel', 'cs', 'approach', 'dm');
        %Gref = GradientMaps('kernel', 'cs', 'approach', 'dm');
        Gref = Gref.fit(conn_matrix1,'sparsity', 90 );
         plot_hemispheres([Gref.gradients{1}(:,1), Gref.gradients{1}(:,2)],...
                    {surfL,surfR}, ...
                    'parcellation', schaefer_400_label);    

        for i = 1: length(idx_ses2)

            conn_matrix_ind = squeeze(fc_sch400_nki_ses2{1, i}');
            conn_matrix_ind(conn_matrix_ind==Inf) = 1;

            gm_ind = GradientMaps('kernel', 'cs', 'approach', 'dm', 'alignment', 'pa');
                %gm_ind = gm_ind.fit({conn_matrix1, conn_matrix_ind}, 'sparsity', 90 );   
            gm_ind = gm_ind.fit(conn_matrix_ind, 'reference', Gref.gradients{1}, 'sparsity', 90)

            orig_grad1_noAlign = gm_ind.gradients{1}(:,1);
            orig_grad2_noAlign = gm_ind.gradients{1}(:,2);   
            orig_grad3_noAlign = gm_ind.gradients{1}(:,3);   

            orig_grad1_align = gm_ind.aligned{1}(:,1);
            orig_grad2_align = gm_ind.aligned{1}(:,2);           
            orig_grad3_align = gm_ind.aligned{1}(:,3);           

            norm_grad1_noAlign = normalize(gm_ind.gradients{1}(:,1));
            norm_grad2_noAlign = normalize(gm_ind.gradients{1}(:,2));   
            norm_grad3_noAlign = normalize(gm_ind.gradients{1}(:,3));   

            norm_grad1_align = normalize(gm_ind.aligned{1}(:,1));
            norm_grad2_align = normalize(gm_ind.aligned{1}(:,2));   
            norm_grad3_align = normalize(gm_ind.aligned{1}(:,3));   

            %store as struct
            ind_grad.orig_grad1_noAlign{i} = orig_grad1_noAlign;
            ind_grad.orig_grad2_noAlign{i} = orig_grad2_noAlign;
            ind_grad.orig_grad3_noAlign{i} = orig_grad3_noAlign;

            ind_grad.orig_grad1_align{i} = orig_grad1_align;
            ind_grad.orig_grad2_align{i} = orig_grad2_align;
            ind_grad.orig_grad3_align{i} = orig_grad3_align;

            ind_grad.norm_grad1_noAlign{i} = norm_grad1_noAlign;
            ind_grad.norm_grad2_noAlign{i} = norm_grad2_noAlign;
            ind_grad.norm_grad3_noAlign{i} = norm_grad3_noAlign;

            ind_grad.norm_grad1_align{i} = norm_grad1_align;
            ind_grad.norm_grad2_align{i} = norm_grad2_align;
            ind_grad.norm_grad3_align{i} = norm_grad3_align;

        end

        %change struct to array
        for i = 1:length(idx_ses2)

                ind_orig_grad1_align(i,:,:) = ind_grad.orig_grad1_align{1, i};
                ind_orig_grad2_align(i,:,:) = ind_grad.orig_grad2_align{1, i};
                ind_orig_grad3_align(i,:,:) = ind_grad.orig_grad3_align{1, i};

                ind_norm_grad1_align(i,:,:) = ind_grad.norm_grad1_align{1, i};
                ind_norm_grad2_align(i,:,:) = ind_grad.norm_grad2_align{1, i};
                ind_norm_grad3_align(i,:,:) = ind_grad.norm_grad3_align{1, i};

            end


    %% prepare Y: Follow-up Ext-Int segregation in CC gradients

    for for_prepare_extInt = 1
        load('./files/yeo7_in_sch400_label.mat')
        yeo7_in_sch400_vis = find(yeo7_in_sch400==1)
        yeo7_in_sch400_som = find(yeo7_in_sch400==2)
        yeo7_in_sch400_dan = find(yeo7_in_sch400==3)
        yeo7_in_sch400_sal = find(yeo7_in_sch400==4)
        yeo7_in_sch400_lim = find(yeo7_in_sch400==5)
        yeo7_in_sch400_fpn = find(yeo7_in_sch400==6)
        yeo7_in_sch400_dmn = find(yeo7_in_sch400==7)

        ind_orig_grad2_align_int = ind_orig_grad2_align(:,[yeo7_in_sch400_dmn] )
        ind_orig_grad2_align_int_mean = mean(ind_orig_grad2_align_int,2)
        ind_orig_grad2_align_ext = ind_orig_grad2_align(:,[yeo7_in_sch400_vis;yeo7_in_sch400_som;yeo7_in_sch400_dan ] )
        ind_orig_grad2_align_ext_mean = mean(ind_orig_grad2_align_ext,2)
        seg_int_ext = abs(ind_orig_grad2_align_ext_mean - ind_orig_grad2_align_int_mean) 

        Y = seg_int_ext;   

    end
    %% prepare X: Baseline NEOMAP 1 & 2 of each participant

    %Delta NEOMAP
    pmap2_delta = pmap2_L_thal_ind_groupCmap_nki_ses2_align_scale - pmap2_L_thal_ind_groupCmap_nki_ses1_align_scale;
    pmap3_delta = pmap3_L_thal_ind_groupCmap_nki_ses2_align_scale - pmap3_L_thal_ind_groupCmap_nki_ses1_align_scale;
    X=[pmap2_delta, pmap3_delta]; 
    %% run prediction  
        i=1 

        for i=1:100 
           CVMD1_ccgrad_int_ext{i} = fitrlinear(X, Y,'KFold',5, ...
               'Learner','leastsquares','Solver','sparsa','Regularization','lasso');
           MSE_ccgrad_int_ext(:, i) = kfoldLoss(CVMD1_ccgrad_int_ext{i});
           yHat_ccgrad_int_ext(:, i) = kfoldPredict(CVMD1_ccgrad_int_ext{i});
           MAE_ccgrad_int_ext(:, i) = mae(Y - yHat_ccgrad_int_ext);
        end

        [yHat_ccgrad_int_ext_Y_r yHat_ccgrad_int_ext_Y_p] = corr(Y, yHat_ccgrad_int_ext)
        yHat_ccgrad_int_ext_r_SEM = std(yHat_ccgrad_int_ext_Y_r)/sqrt(length(yHat_ccgrad_int_ext_Y_r)); 
        temps = tinv([0.025 0.975], length(yHat_ccgrad_int_ext_Y_r)-1); 
        yHat_ccgrad_int_ext_CI = mean(yHat_ccgrad_int_ext_Y_r) + temps*yHat_ccgrad_int_ext_r_SEM

        yHat_ccgrad_int_ext_median = median(yHat_ccgrad_int_ext, 2)
       
        figure; scatter(Y, yHat_ccgrad_int_ext_median)
        [r p] = corr(Y, yHat_ccgrad_int_ext_median)
        
        MAE_ccgrad_int_ext_mean = mean(MAE_ccgrad_int_ext)
        MAE_ccgrad_int_ext_SEM = std(MAE_ccgrad_int_ext)/sqrt(length(MAE_ccgrad_int_ext))
        MAE_ccgrad_int_ext_STD = std(MAE_ccgrad_int_ext)
        
        Obs_Y = Y
        Pred_Y = yHat_ccgrad_int_ext_median
        sz=70
        
        % Regression using matrix notation
        B = [ones(size(Obs_Y)), Obs_Y]\Pred_Y;
        Slope = B(2);
        Intercept = B(1);

        % Calculate predictions and residuals
        Pred_Y_fit = Intercept + Slope * Obs_Y;
        residuals = Pred_Y - Pred_Y_fit;

        % Standard error of the slope
        SE = sqrt(sum(residuals.^2) / (length(Obs_Y) - 2)) / sqrt(sum((Obs_Y - mean(Obs_Y)).^2));
        t_value = tinv(0.975, length(Obs_Y) - 2);  % 95% CI, two-tailed

        % Confidence interval calculation
        CI_upper = (Slope + t_value*SE) * Obs_Y + Intercept;
        CI_lower = (Slope - t_value*SE) * Obs_Y + Intercept;

        % Plotting
        figure;
        scatter(Obs_Y, Pred_Y, sz, 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [0 0 0], 'LineWidth', 1.5);
        hold on;
        plot(Obs_Y, Pred_Y_fit, 'Color', [0 0 0], 'LineWidth', 2);  % regression line
        plot(Obs_Y, CI_upper, 'k:', 'LineWidth', 0.001);  % Upper CI
        plot(Obs_Y, CI_lower, 'k:', 'LineWidth', 0.001);  % Lower CI
        
        xlabel('Observed values');
        ylabel('Predicted values');

        ax=gca; 
        ax.FontSize=14;
        set(gcf,'color','w'); grid on;
        set(gcf,'Position',[500 500 500 500])



