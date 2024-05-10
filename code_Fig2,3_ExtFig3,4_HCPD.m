
%% -- Preparation ------------------------------------------
    %% load toolboxes and define paths
    for pathtoolbox = 1

        %download the following tools
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
        

    end
    %% read phenotypic info
    for for_readdata=1
        %space_flag = '_10k' %'_10k' %_32k
        data = importfile_demo(['./files/demo_dHCP,HCPD_allAges_qcgrad.xlsx']);
        
    end    
    for for_demo_vars = 1
        age_w = data.age_in_weeks;
        age_y = data.age_in_years;
        sex = data.sex;
        subid = data.sub_id;
        study = data.dataset;    
        session = data.session;
        taskrest = data.task_rest;
        meanFD = data.meanFD;
        meanFD04 = data.MeanFD04;
        qc_grad = data.qc_grad;

        idx_hcpd = find( (study == 'HCPD') & (session == 'session1') & (qc_grad == 1) & (meanFD04 == 1));

        age_w_hcpd = age_w(idx_hcpd);
        age_y_hcpd = age_y(idx_hcpd);
        sex_hcpd = sex(idx_hcpd);
        sex_hcpd = cellstr(sex_hcpd);
        subid_hcpd = subid(idx_hcpd);
        meanFD_hcpd = meanFD(idx_hcpd);        

        idx_hcpd_8_10y = find (age_y_hcpd < 10) ;  
        idx_hcpd_10_13y = find (age_y_hcpd >= 10 & age_y_hcpd < 13) ;  
        idx_hcpd_13_16y = find (age_y_hcpd >= 13 & age_y_hcpd < 16) ;  
        idx_hcpd_16_20y = find (age_y_hcpd >= 16 & age_y_hcpd < 20) ;  
        idx_hcpd_20_22y = find (age_y_hcpd >= 20) ;  


        idx_age_8_12y =find(age_y_hcpd<12)
        idx_age_12_18y =find(age_y_hcpd>=12 & age_y_hcpd<18)
        idx_age_18_22y =find(age_y_hcpd>18)     

    end
    %% Load surfaces & parcellations & colormaps

    for for_load_surfaces = 1

        %upload 10k surfaces: fsLR space         
        temp = gifti(['./files/S900.R.midthickness_MSMAll.10k_fs_LR.surf.gii']);
        surfR.coord = temp.vertices';
        surfR.tri   = temp.faces;   

        temp = gifti(['./files/S900.L.midthickness_MSMAll.10k_fs_LR.surf.gii']);
        %temp = gifti([ PATH '01_analysis/03_surface/MNI152_T1_1mm.L.very_inflated.10k_fs_LR.surf.gii']);
        surfL.coord = temp.vertices';
        surfL.tri   = temp.faces;

        surf.coord = [ surfL.coord surfR.coord ];
        surf.tri   = [ surfL.tri; surfR.tri+10242; ];

        %temp = gifti('/local_raid1/01_software/HCPpipelines/global/templates/standard_mesh_atlases/L.atlasroi.10k_fs_LR.shape.gii');
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
        yeo_colormap = [ 200 200 200;    
                     120 18 134;
                     70 130 180;
                     0 118 14;
                     196 58 250;
                     220 248 164;
                     230 148 34;
                     205 62 78 ]/255;
        yeo_colormap_1=yeo_colormap(2:end,:)

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
    
    end

    
%% -- Main Figure 2 ----------------------------------------------
    %% upload source data for NEOMAPs
    
    %includes individual NEOMAPs (603 individuals x 20484 vertices x 10 NEOMAPs)
    load('sourceData_Fig2,3_ExtFig3,4_HCPD_NEOMAPs.mat')
    
    pmap_R_thal_ind_groupCmap_hcpd_p1 = pmap_R_thal_ind_groupCmap_hcpd_align(:,:,1);
    pmap_R_thal_ind_groupCmap_hcpd_p2 = pmap_R_thal_ind_groupCmap_hcpd_align(:,:,2);

    pmap_L_thal_ind_groupCmap_hcpd_p1 = pmap_L_thal_ind_groupCmap_hcpd_align(:,:,1);
    pmap_L_thal_ind_groupCmap_hcpd_p2 = pmap_L_thal_ind_groupCmap_hcpd_align(:,:,2);
    
    %using normalized pmap 
    for i= 1:length(idx_hcpd)
        temp_pmap1 =  squeeze(pmap_L_thal_ind_groupCmap_hcpd_align(i,:,1));     
        temp1 = normalize(temp_pmap1);
        pmap1_L_thal_ind_groupCmap_hcpd_align(i,:) = temp1;

        temp_pmap2 =  squeeze(pmap_L_thal_ind_groupCmap_hcpd_align(i,:,2));     
        temp2 = normalize(temp_pmap2);
        pmap2_L_thal_ind_groupCmap_hcpd_align(i,:) = temp2;

        temp_pmap3 =  squeeze(pmap_L_thal_ind_groupCmap_hcpd_align(i,:,3));     
        temp3 = normalize(temp_pmap3);
        pmap3_L_thal_ind_groupCmap_hcpd_align(i,:) = temp3;


        temp_pmap1 =  squeeze(pmap_R_thal_ind_groupCmap_hcpd_align(i,:,1));     
        temp1 = normalize(temp_pmap1);
        pmap1_R_thal_ind_groupCmap_hcpd_align(i,:) = temp1;

        temp_pmap2 =  squeeze(pmap_R_thal_ind_groupCmap_hcpd_align(i,:,2));     
        temp2 = normalize(temp_pmap2);
        pmap2_R_thal_ind_groupCmap_hcpd_align(i,:) = temp2;

        temp_pmap3 =  squeeze(pmap_R_thal_ind_groupCmap_hcpd_align(i,:,3));     
        temp3 = normalize(temp_pmap3);
        pmap3_R_thal_ind_groupCmap_hcpd_align(i,:) = temp3;
    end
 
            
   
    for for_create_source_data = 1 %do not run
        %% read NEOMAPS from groupCmap 
        for read_individual_pmaps_from_groupCmap = 1
%             
%             % extracted individual pmaps using the group cmap template with python code
%             PATH               = '/local_raid2/03_user/shinwon/01_project/01_thalamusgradient/';
%             IND_PMAP_HCPD_GROUPCMAP = [PATH '03_individual_gradient/06_aligned_pmaps/02_HCPD_groupCmap_rev10vols_cmap10000/' ];

            for hcpd=1
            for i = 1:length(idx_hcpd)

%                 pmap_L = ciftiopen([ IND_PMAP_HCPD_GROUPCMAP 'sub-' subid_hcpd{i} '/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.10.pmap.dscalar.nii' ], WB_COMMAND);   
%                 pmap_L_thal = pmap_L.cdata(1:(pmap_L.diminfo{1,1}.models{1,1}.count) + pmap_L.diminfo{1,1}.models{2,1}.count,:);
% 
%                 pmap_L_full = zeros(pmap_L.diminfo{1,1}.models{1,1}.numvert+pmap_L.diminfo{1,1}.models{2,1}.numvert, pmap_L.diminfo{1,2}.length);
%                 pmap_L_full(logical( surf_roi ),:) = pmap_L_thal;
% 
%                 pmap_L_thal_ind_groupCmap_hcpd(i, :, :) = pmap_L_full;   
% 
%                 pmap_R = ciftiopen([ IND_PMAP_HCPD_GROUPCMAP 'sub-' subid_hcpd{i} '/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.49.49.pmap.dscalar.nii' ], WB_COMMAND);   
%                 pmap_R_thal = pmap_R.cdata(1:(pmap_R.diminfo{1,1}.models{1,1}.count) + pmap_R.diminfo{1,1}.models{2,1}.count,:);
% 
%                 pmap_R_full = zeros(pmap_R.diminfo{1,1}.models{1,1}.numvert+pmap_R.diminfo{1,1}.models{2,1}.numvert, pmap_R.diminfo{1,2}.length);
%                 pmap_R_full(logical( surf_roi ),:) = pmap_R_thal;
% 
%                 pmap_R_thal_ind_groupCmap_hcpd(i, :, :) = pmap_R_full;    

            end    
        end

        end  
    
        %% align individual NEOMAPS to Group NEOMAPS
        for for_align_ind_pmaps_to_groupPmap = 1

            for make_group_pmap = 1
%                 %use mean pmap of term hcpds as template 
%                 group_pmap_R_thal_ind_groupCmap_hcpd = nanmean(pmap_R_thal_ind_groupCmap_hcpd(:,:,:),1);
%                 group_pmap_R_thal_ind_groupCmap_hcpd_rs = reshape(group_pmap_R_thal_ind_groupCmap_hcpd, [20484,10]);
% 
%                 group_pmap_R_thal_ind_groupCmap_hcpd_rs_nz = group_pmap_R_thal_ind_groupCmap_hcpd_rs;
%                 temp_idx=find(group_pmap_R_thal_ind_groupCmap_hcpd_rs_nz(:,1)==0);
%                 group_pmap_R_thal_ind_groupCmap_hcpd_rs_nz(temp_idx,:)=[];
% 
%                 group_pmap_R_thal_ind_groupCmap_hcpd_rs_nz_norm=normalize(group_pmap_R_thal_ind_groupCmap_hcpd_rs_nz);
% 
%                 group_pmap_L_thal_ind_groupCmap_hcpd = nanmean(pmap_L_thal_ind_groupCmap_hcpd(:,:,:),1);
%                 group_pmap_L_thal_ind_groupCmap_hcpd_rs = reshape(group_pmap_L_thal_ind_groupCmap_hcpd, [20484,10]);
% 
%                 group_pmap_L_thal_ind_groupCmap_hcpd_rs_nz = group_pmap_L_thal_ind_groupCmap_hcpd_rs;
%                 temp_idx=find(group_pmap_L_thal_ind_groupCmap_hcpd_rs_nz(:,1)==0);
%                 group_pmap_L_thal_ind_groupCmap_hcpd_rs_nz(temp_idx,:)=[];
% 
%                 group_pmap_L_thal_ind_groupCmap_hcpd_rs_nz_norm=normalize(group_pmap_L_thal_ind_groupCmap_hcpd_rs_nz);
% 
%                 figure; BoSurfStatView (group_pmap_L_thal_ind_groupCmap_hcpd_rs(:,1),surf);%BoSurfStatColLim([-0.05 0.05]);colormap(rbgyr20_colormap)
%                 figure; BoSurfStatView (group_pmap_L_thal_ind_groupCmap_hcpd_rs(:,2),surf);%BoSurfStatColLim([-0.05 0.05]);colormap(rbgyr20_colormap)
%                 figure; BoSurfStatView (group_pmap_L_thal_ind_groupCmap_hcpd_rs(:,3),surf);%BoSurfStatColLim([-0.05 0.05]);colormap(rbgyr20_colormap)
% 
%                 %change normalized to 20484 for visualization before alignment
%                 group_pmap_L_thal_ind_groupCmap_hcpd_rs_norm_full = zeros(20484, 10);
%                 group_pmap_L_thal_ind_groupCmap_hcpd_rs_norm_full(logical(surf_roi),:) = group_pmap_L_thal_ind_groupCmap_hcpd_rs_nz_norm;
% 
%                 figure; BoSurfStatView (group_pmap_L_thal_ind_groupCmap_hcpd_rs_norm_full(:,1),surf);%BoSurfStatColLim([-3 3]);colormap(rbgyr20_colormap)
%                 figure; BoSurfStatView (group_pmap_L_thal_ind_groupCmap_hcpd_rs_norm_full(:,2),surf);%BoSurfStatColLim([-3 3]);colormap(rbgyr20_colormap)

            end

            for align_indPmap_to_groupPMAP = 1
% 
%                 pmap_L_thal_ind_groupCmap_hcpd_align = zeros(length(idx_hcpd), 20484, 10);
%                 pmap_R_thal_ind_groupCmap_hcpd_align = zeros(length(idx_hcpd), 20484, 10);

                parfor i = 1:length(idx_hcpd)
% 
%                     pmap_L = ciftiopen([ IND_PMAP_HCPD_GROUPCMAP 'sub-' subid_hcpd{i} '/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.10.pmap.dscalar.nii' ], WB_COMMAND);   
%                     pmap_L_thal = pmap_L.cdata(1:(pmap_L.diminfo{1,1}.models{1,1}.count) + pmap_L.diminfo{1,1}.models{2,1}.count,:);
% 
% 
%                     pmap_R = ciftiopen([ IND_PMAP_HCPD_GROUPCMAP 'sub-' subid_hcpd{i} '/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.49.49.pmap.dscalar.nii' ], WB_COMMAND);   
%                     pmap_R_thal = pmap_R.cdata(1:(pmap_R.diminfo{1,1}.models{1,1}.count) + pmap_R.diminfo{1,1}.models{2,1}.count,:);
% 
%                     %procrustes alignment 
%                     [d_L,Z_L,transform_L] = procrustes(group_pmap_L_thal_ind_groupCmap_hcpd_rs_nz_norm,pmap_L_thal); %18722
%                     [d_R,Z_R,transform_R] = procrustes(group_pmap_R_thal_ind_groupCmap_hcpd_rs_nz_norm,pmap_R_thal);
% 
%                     temp_L = zeros(20484,10);
%                     temp_L(logical(surf_roi),:) = Z_L;
%                     temp_R = zeros(20484,10);
%                     temp_R(logical(surf_roi),:) = Z_R;
%                     
%                     pmap_L_thal_ind_groupCmap_hcpd_align(i, :, :) = temp_L;
%                     pmap_R_thal_ind_groupCmap_hcpd_align(i, :, :) = temp_R;
% 
%                     pmap_L_thal_ind_groupCmap_hcpd_align_nz(i, :, :) = Z_L;
%                     pmap_R_thal_ind_groupCmap_hcpd_align_nz(i, :, :) = Z_R;
                end   
            end
    
    
    
        end

    %%using aligned pmaps 
% 
%     pmap_R_thal_ind_groupCmap_hcpd_p1 = pmap_R_thal_ind_groupCmap_hcpd_align(:,:,1);
%     pmap_R_thal_ind_groupCmap_hcpd_p2 = pmap_R_thal_ind_groupCmap_hcpd_align(:,:,2);
% 
%     pmap_L_thal_ind_groupCmap_hcpd_p1 = pmap_L_thal_ind_groupCmap_hcpd_align(:,:,1);
%     pmap_L_thal_ind_groupCmap_hcpd_p2 = pmap_L_thal_ind_groupCmap_hcpd_align(:,:,2);
% 
%     pmap_R_thal_ind_groupCmap_hcpd_p1_nz = pmap_R_thal_ind_groupCmap_hcpd_align_nz(:,:,1);
%     pmap_R_thal_ind_groupCmap_hcpd_p2_nz = pmap_R_thal_ind_groupCmap_hcpd_align_nz(:,:,2);
% 
%     pmap_L_thal_ind_groupCmap_hcpd_p1_nz = pmap_L_thal_ind_groupCmap_hcpd_align_nz(:,:,1);
%     pmap_L_thal_ind_groupCmap_hcpd_p2_nz = pmap_L_thal_ind_groupCmap_hcpd_align_nz(:,:,2);


    end
    
    
    %% boxplot for neomap values
    for for_yeo7_networks_neomap1 = 1
        %%pmap1
        pmap_t1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(:,:,1),1);
        figure; SurfStatView (mean(pmap_t1,1),surf)
        SurfStatColLim([-0.5 0.5]);colormap(flipud(spectral_colormap))

        labels = yeo_7_label_full      
        idx=find((labels>0));
        pmap = pmap_t1(1,idx)
        labels=labels(idx)

        uniqueLabels=unique(labels)
        groupMeans = arrayfun(@(x) mean(pmap(labels == x)), uniqueLabels(1:end,:));

        % Sort the groups by their means
        [sortedMeans, sortOrder] = sort(groupMeans);
        sortedLabels = uniqueLabels(sortOrder);
        
        % Sort the colors according to the sorted labels
        my_colormap_0=yeo_colormap_1
        sortedColors = my_colormap_0(sortOrder, :);

        % Create a figure and boxplot with sorted groups
        figure;
        boxplot(pmap, labels, 'GroupOrder', string(sortOrder), 'Colors', sortedColors, 'Symbol', '');
        set(gca,'XTickLabel', get(gca, 'XTickLabel'), 'FontSize',14);

        uniqueLabels=unique(labels);
        stats=struct();
        for i=1:length(uniqueLabels)
            groupData = pmap(labels==uniqueLabels(i));
            stats(i).min=min(groupData);
            stats(i).Q1=prctile(groupData,25); 
            stats(i).median = median(groupData);
            stats(i).Q3 = prctile(groupData, 75); 
            stats(i).max = max(groupData);
            stats(i).IQR = stats(i).Q3 - stats(i).Q1;         
        end


        hold on;

        % Plot a colored scatter plot for sorted groups
        for i = 1:length(sortedLabels)
            label = sortedLabels(i);
            mask = labels == label;
            scatter(i * ones(1, sum(mask)), pmap(mask), 1, sortedColors(i, :), 'filled', 'MarkerFaceAlpha', 0.6, 'jitter', 'on', 'jitterAmount', 0.1);
        end
        set (gcf,'color','w');
        ylim([-1 1]);

        hold off;

    end
    for for_yeo7_networks_neomap2 = 1
        %%pmap2
        pmap_t1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(:,:,2),1);
        labels = yeo_7_label_full      
        idx=find((labels>0));
        pmap = pmap_t1(1,idx)
        labels=labels(idx)

        uniqueLabels=unique(labels)
        groupMeans = arrayfun(@(x) mean(pmap(labels == x)), uniqueLabels(1:end,:));

        % Sort the groups by their means
        [sortedMeans, sortOrder] = sort(groupMeans);
        sortedLabels = uniqueLabels(sortOrder);
        % Sort the colors according to the sorted labels
        my_colormap_0=yeo_colormap_1
        sortedColors = my_colormap_0(sortOrder, :);

        % Create a figure and boxplot with sorted groups
        figure;
        boxplot(pmap, labels, 'GroupOrder', string(sortOrder), 'Colors', sortedColors, 'Symbol', '');
        set(gca,'XTickLabel', get(gca, 'XTickLabel'), 'FontSize',14);
        hold on;

        % Plot a colored scatter plot for sorted groups
        for i = 1:length(sortedLabels)
            label = sortedLabels(i);
            mask = labels == label;
            scatter(i * ones(1, sum(mask)), pmap(mask), 1, sortedColors(i, :), 'filled', 'MarkerFaceAlpha', 0.6, 'jitter', 'on', 'jitterAmount', 0.1);
        end
        set (gcf,'color','w');
        ylim([-1 1]);
        hold off;

    end
    
%% -- Main Figure 3A ----------------------------------------------
    %% Calculate segregation bt Salience & External / Internal

    idx_da = [find(yeo_7_label_full==3), find(yeo_7_label_full==1), find(yeo_7_label_full==2)]';
    idx_sal = [find(yeo_7_label_full==4) ]';
    idx_dmn = [find(yeo_7_label_full==7) ]';
    
    for for_neomap1=1 
        %compute mean of each networks of interest
        mean_da_hcpd_8_10y_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_8_10y,idx_da,1),2); 
        mean_sal_hcpd_8_10y_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_8_10y,idx_sal,1),2);
        mean_dmn_hcpd_8_10y_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_8_10y,idx_dmn,1),2);

        mean_da_hcpd_10_13y_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_10_13y,idx_da,1),2);
        mean_sal_hcpd_10_13y_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_10_13y,idx_sal,1),2);
        mean_dmn_hcpd_10_13y_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_10_13y,idx_dmn,1),2);

        mean_da_hcpd_13_16y_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_13_16y,idx_da,1),2);
        mean_sal_hcpd_13_16y_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_13_16y,idx_sal,1),2);
        mean_dmn_hcpd_13_16y_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_13_16y,idx_dmn,1),2);
      
        mean_da_hcpd_16_20y_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_16_20y,idx_da,1),2);
        mean_sal_hcpd_16_20y_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_16_20y,idx_sal,1),2);
        mean_dmn_hcpd_16_20y_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_16_20y,idx_dmn,1),2);
      
        mean_da_hcpd_20_22y_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_20_22y,idx_da,1),2);
        mean_sal_hcpd_20_22y_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_20_22y,idx_sal,1),2);
        mean_dmn_hcpd_20_22y_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_20_22y,idx_dmn,1),2);

        mean_da_hcpd_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(:,idx_da,1),2);
        mean_sal_hcpd_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(:,idx_sal,1),2);
        mean_dmn_hcpd_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(:,idx_dmn,1),2);
        
        
        mean_da_hcpd_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(:,idx_da,1),2);
        mean_sal_hcpd_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(:,idx_sal,1),2);
        mean_dmn_hcpd_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(:,idx_dmn,1),2);
        
        for for_difference = 1
            %diff between control & da
            diff_con_da_hcpd_8_10y_pmap1 = abs(mean_sal_hcpd_8_10y_pmap1 - mean_da_hcpd_8_10y_pmap1);

            diff_con_da_hcpd_10_13y_pmap1 = abs(mean_sal_hcpd_10_13y_pmap1 - mean_da_hcpd_10_13y_pmap1);

            diff_con_da_hcpd_13_16y_pmap1 = abs(mean_sal_hcpd_13_16y_pmap1 - mean_da_hcpd_13_16y_pmap1);

            diff_con_da_hcpd_16_20y_pmap1 = abs(mean_sal_hcpd_16_20y_pmap1 - mean_da_hcpd_16_20y_pmap1);

            diff_con_da_hcpd_20_22y_pmap1 = abs(mean_sal_hcpd_20_22y_pmap1 - mean_da_hcpd_20_22y_pmap1);

            diff_con_da_hcpd_pmap1 = abs(mean_sal_hcpd_pmap1 - mean_da_hcpd_pmap1);



            %diff between control & dmn
            diff_con_dmn_hcpd_8_10y_pmap1 = abs(mean_sal_hcpd_8_10y_pmap1 - mean_dmn_hcpd_8_10y_pmap1);

            diff_con_dmn_hcpd_10_13y_pmap1 = abs(mean_sal_hcpd_10_13y_pmap1 - mean_dmn_hcpd_10_13y_pmap1);

            diff_con_dmn_hcpd_13_16y_pmap1 = abs(mean_sal_hcpd_13_16y_pmap1 - mean_dmn_hcpd_13_16y_pmap1);

            diff_con_dmn_hcpd_16_20y_pmap1 = abs(mean_sal_hcpd_16_20y_pmap1 - mean_dmn_hcpd_16_20y_pmap1);

            diff_con_dmn_hcpd_20_22y_pmap1 = abs(mean_sal_hcpd_20_22y_pmap1 - mean_dmn_hcpd_20_22y_pmap1);

            diff_con_dmn_hcpd_pmap1 = abs(mean_sal_hcpd_pmap1 - mean_dmn_hcpd_pmap1);
        
        end
    
    end    
    for for_neomap2=1 
        %compute mean of each networks of interest
        mean_da_hcpd_8_10y_pmap2 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_8_10y,idx_da,2),2); 
        mean_sal_hcpd_8_10y_pmap2 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_8_10y,idx_sal,2),2);
        mean_dmn_hcpd_8_10y_pmap2 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_8_10y,idx_dmn,2),2);

        mean_da_hcpd_10_13y_pmap2 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_10_13y,idx_da,2),2);
        mean_sal_hcpd_10_13y_pmap2 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_10_13y,idx_sal,2),2);
        mean_dmn_hcpd_10_13y_pmap2 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_10_13y,idx_dmn,2),2);

        mean_da_hcpd_13_16y_pmap2 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_13_16y,idx_da,2),2);
        mean_sal_hcpd_13_16y_pmap2 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_13_16y,idx_sal,2),2);
        mean_dmn_hcpd_13_16y_pmap2 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_13_16y,idx_dmn,2),2);
        
        mean_da_hcpd_16_20y_pmap2 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_16_20y,idx_da,2),2);
        mean_sal_hcpd_16_20y_pmap2 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_16_20y,idx_sal,2),2);
        mean_dmn_hcpd_16_20y_pmap2 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_16_20y,idx_dmn,2),2);
     
        mean_da_hcpd_20_22y_pmap2 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_20_22y,idx_da,2),2);
        mean_sal_hcpd_20_22y_pmap2 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_20_22y,idx_sal,2),2);
        mean_dmn_hcpd_20_22y_pmap2 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(idx_hcpd_20_22y,idx_dmn,2),2);

        mean_da_hcpd_pmap2 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(:,idx_da,2),2); 
        mean_sal_hcpd_pmap2 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(:,idx_sal,2),2);
        mean_dmn_hcpd_pmap2 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(:,idx_dmn,2),2);
                
        for for_difference = 1
            %diff between control & da
            diff_con_da_hcpd_8_10y_pmap2 = abs(mean_sal_hcpd_8_10y_pmap2 - mean_da_hcpd_8_10y_pmap2);

            diff_con_da_hcpd_10_13y_pmap2 = abs(mean_sal_hcpd_10_13y_pmap2 - mean_da_hcpd_10_13y_pmap2);

            diff_con_da_hcpd_13_16y_pmap2 = abs(mean_sal_hcpd_13_16y_pmap2 - mean_da_hcpd_13_16y_pmap2);

            diff_con_da_hcpd_16_20y_pmap2 = abs(mean_sal_hcpd_16_20y_pmap2 - mean_da_hcpd_16_20y_pmap2);

            diff_con_da_hcpd_20_22y_pmap2 = abs(mean_sal_hcpd_20_22y_pmap2 - mean_da_hcpd_20_22y_pmap2);

            diff_con_da_hcpd_pmap2 = abs(mean_sal_hcpd_pmap2 - mean_da_hcpd_pmap2);

            %diff between control & dmn
            diff_con_dmn_hcpd_8_10y_pmap2 = abs(mean_sal_hcpd_8_10y_pmap2 - mean_dmn_hcpd_8_10y_pmap2);

            diff_con_dmn_hcpd_10_13y_pmap2 = abs(mean_sal_hcpd_10_13y_pmap2 - mean_dmn_hcpd_10_13y_pmap2);

            diff_con_dmn_hcpd_13_16y_pmap2 = abs(mean_sal_hcpd_13_16y_pmap2 - mean_dmn_hcpd_13_16y_pmap2);

            diff_con_dmn_hcpd_16_20y_pmap2 = abs(mean_sal_hcpd_16_20y_pmap2 - mean_dmn_hcpd_16_20y_pmap2);

            diff_con_dmn_hcpd_20_22y_pmap2 = abs(mean_sal_hcpd_20_22y_pmap2 - mean_dmn_hcpd_20_22y_pmap2);

            diff_con_dmn_hcpd_pmap2 = abs(mean_sal_hcpd_pmap2 - mean_dmn_hcpd_pmap2);
        
        end

    end
    %% boxplots
    
    for for_neomap1 = 1

        %combined 
        
        test = [diff_con_da_hcpd_8_10y_pmap1; diff_con_dmn_hcpd_8_10y_pmap1; diff_con_da_hcpd_10_13y_pmap1; diff_con_dmn_hcpd_10_13y_pmap1; ...
            diff_con_da_hcpd_13_16y_pmap1; diff_con_dmn_hcpd_13_16y_pmap1; diff_con_da_hcpd_16_20y_pmap1; diff_con_dmn_hcpd_16_20y_pmap1; ...
            diff_con_da_hcpd_20_22y_pmap1; diff_con_dmn_hcpd_20_22y_pmap1];
        test_g = [repmat({'1A'}, length(diff_con_da_hcpd_8_10y_pmap1),1); repmat({'1B'}, length(diff_con_dmn_hcpd_8_10y_pmap1),1); ... 
            repmat({'2A'}, length(diff_con_da_hcpd_10_13y_pmap1),1); repmat({'2B'}, length(diff_con_dmn_hcpd_10_13y_pmap1),1); ... 
            repmat({'3A'}, length(diff_con_da_hcpd_13_16y_pmap1),1); repmat({'3B'}, length(diff_con_dmn_hcpd_13_16y_pmap1),1); ... 
            repmat({'4A'}, length(diff_con_da_hcpd_16_20y_pmap1),1); repmat({'4B'}, length(diff_con_dmn_hcpd_16_20y_pmap1),1); ... 
            repmat({'5A'}, length(diff_con_da_hcpd_20_22y_pmap1),1); repmat({'5B'}, length(diff_con_dmn_hcpd_20_22y_pmap1),1)];
        %1=8-10y; 2=10-13y; 3=13-16y; 4=16-20y; 5=20-22y
        %A=SAL-EXT; B=SAL-INT
   
        figure;
        h = boxplot(test, test_g, 'Symbol', '');
        ylim([0 1.6])
        set(gca,'YTickLabel', get(gca, 'YTickLabel'), 'FontSize',14);
        set(gca,'XTickLabel',{});

        numBoxes = length(unique(test_g));

        colorsA = yeo_colormap(2,:); 
        colorsB = yeo_colormap(8,:); 
        
        hBox = findobj(gca, 'Tag', 'Box');

         %check hBox order
        xPositions = zeros(length(hBox), 1);

        % Loop through each handle and get the XData
        for iBox = 1:length(hBox)
            xData = get(hBox(iBox), 'XData');
            xPositions(iBox) = median(xData(:));
        end

        % Now xPositions contains the x-positions of the centers of the boxes
        disp(xPositions);

        % Check if the xPositions are in the expected order
        isInOrder = all(diff(xPositions) > 0); % This should be true if the order is ascending

        % Sort xPositions in ascending order and get the indices
        [sortedX, sortOrder] = sort(xPositions);

        % Use sortOrder to reorder hBox
        hBox = hBox(sortOrder);

        % Now hBox is reordered according to the ascending x-positions of the boxes

        
        for iBox = 1:length(hBox)
            patchData = get(hBox(iBox), 'XData');
            if mod(iBox, 2)
                % Odd indices are Group A
                patch(patchData, get(hBox(iBox), 'YData'), colorsA, 'FaceAlpha', .5);
            else
                patch(patchData, get(hBox(iBox), 'YData'), colorsB, 'FaceAlpha', .5);
            end
        end
        
        set(gcf,'color','w');
        
    
end
    for for_neomap2 = 1

        %combined 
        
        test = [diff_con_da_hcpd_8_10y_pmap2; diff_con_dmn_hcpd_8_10y_pmap2; diff_con_da_hcpd_10_13y_pmap2; diff_con_dmn_hcpd_10_13y_pmap2; ...
            diff_con_da_hcpd_13_16y_pmap2; diff_con_dmn_hcpd_13_16y_pmap2; diff_con_da_hcpd_16_20y_pmap2; diff_con_dmn_hcpd_16_20y_pmap2; ...
            diff_con_da_hcpd_20_22y_pmap2; diff_con_dmn_hcpd_20_22y_pmap2]
        %1=8-10y; 2=10-13y; 3=13-16y; 4=16-20y; 5=20-22y
        %A=SAL-EXT; B=SAL-INT
        test_g = [repmat({'1A'}, length(diff_con_da_hcpd_8_10y_pmap2),1); repmat({'1B'}, length(diff_con_dmn_hcpd_8_10y_pmap2),1); ... 
            repmat({'2A'}, length(diff_con_da_hcpd_10_13y_pmap2),1); repmat({'2B'}, length(diff_con_dmn_hcpd_10_13y_pmap2),1); ... 
            repmat({'3A'}, length(diff_con_da_hcpd_13_16y_pmap2),1); repmat({'3B'}, length(diff_con_dmn_hcpd_13_16y_pmap2),1); ... 
            repmat({'4A'}, length(diff_con_da_hcpd_16_20y_pmap2),1); repmat({'4B'}, length(diff_con_dmn_hcpd_16_20y_pmap2),1); ... 
            repmat({'5A'}, length(diff_con_da_hcpd_20_22y_pmap2),1); repmat({'5B'}, length(diff_con_dmn_hcpd_20_22y_pmap2),1)];
        figure;
        h = boxplot(test, test_g, 'Symbol', '');
        ylim([0 1.6])
        set(gca,'YTickLabel', get(gca, 'YTickLabel'), 'FontSize',14);
        set(gca,'XTickLabel',{});

        numBoxes = length(unique(test_g));

        colorsA = yeo_colormap(2,:); 
        colorsB = yeo_colormap(8,:); 

        hBox = findobj(gca, 'Tag', 'Box');

            %check hBox order
            % Assuming hBox contains handles to the box objects
            xPositions = zeros(length(hBox), 1);

            % Loop through each handle and get the XData
            for iBox = 1:length(hBox)
                xData = get(hBox(iBox), 'XData');
                xPositions(iBox) = median(xData(:)); % Get the median as a representative value
            end

            % Now xPositions contains the x-positions of the centers of the boxes
            % You can print them out or examine them in the workspace
            disp(xPositions);

            % Check if the xPositions are in the expected order
            isInOrder = all(diff(xPositions) > 0); % This should be true if the order is ascending

            % Sort xPositions in ascending order and get the indices
            [sortedX, sortOrder] = sort(xPositions);

            % Use sortOrder to reorder hBox
            hBox = hBox(sortOrder);

            % Now hBox is reordered according to the ascending x-positions of the boxes
            
        for iBox = 1:length(hBox)
            patchData = get(hBox(iBox), 'XData');
            if mod(iBox, 2)
                patch(patchData, get(hBox(iBox), 'YData'), colorsA, 'FaceAlpha', .5);
            else
                patch(patchData, get(hBox(iBox), 'YData'), colorsB, 'FaceAlpha', .5);
            end
        end
        set(gcf,'color','w');

    
end

%% -- Main Figure 3B ----------------------------------------------

    %% upload genes & coordinate files
    % genetic expression derived using abagen (https://abagen.readthedocs.io/en/stable/)
    % and thalamic parcellation from Saranathan et al (2021)
    THAL_CALB1 = readtable( ['./files/CALB1_lh_thal_nuclei.csv' ]);
    THAL_PVALB = readtable( ['./files/PVALB_lh_thal_nuclei.csv' ]);

    %% upload thalamic nuclei parcellation file 
    atlas_lthal_nuclei = ciftiopen([ './files/Atlas-Thalamus_space-MNI_hemi-left_label-AllNuclei_desc-MaxProb_MNI_2x2x2.dlabel.nii' ], WB_COMMAND);       
    atlas_lthal_nuclei_label = atlas_lthal_nuclei.cdata;

    %% upload group cmap
        cii = ciftiopen([ './files/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.cmap.dscalar.nii' ], WB_COMMAND);   
        cii_L_thal = cii.cdata(cii.diminfo{1, 1}.models{20, 1}.start:(cii.diminfo{1,1}.models{20,1}.start-1) + cii.diminfo{1,1}.models{20,1}.count,:);
        groupGradient_L_hcpd(:, :) = cii_L_thal;

    %% upload thal-cortex correlation for projections
    % originally we need to upload timeseries data to derive thal-cortex correlation (very large file >100 GB) 
    % instead upload the saved thal-cortex correlations 
    
    load('./Figure3B_correlation_thal_ctx.mat')
    
    %% projection of CALB, PVALB onto cortex 
    % construct average thal x cortex correlation 
    for i=1:length(idx_hcpd)
        a = corrR.Lthal_cortsch400{i}([2,4:14],:);
        a(isnan(a)) = 0; 
        corr_Lthal_cortsch400_rev(:,:,i) = a ;
    end

    mean_corr_Lthal_cortsch400 = mean(corr_Lthal_cortsch400_rev,3);

    %normalize thalamic projection
    CORTEX_CALB1_N2 = zscore(THAL_CALB1.CALB1)' * mean_corr_Lthal_cortsch400;
    project_detection_community_boview(CORTEX_CALB1_N2,'', schaefer_400_label,surf,1);
    colormap(parula);BoSurfStatColLim([-0.5 0.35]); 

    CORTEX_PVALB_N2 = zscore(THAL_PVALB.PVALB)' * mean_corr_Lthal_cortsch400;
    project_detection_community_boview(CORTEX_PVALB_N2,'', schaefer_400_label,surf,1);
    colormap(parula);BoSurfStatColLim([-0.5 0.35]); 

    project_detection_community_boview(CORTEX_CALB1_N2-CORTEX_PVALB_N2,'', schaefer_400_label,surf,1);
    colormap(parula);BoSurfStatColLim([-0.5 0.5]); 

    %% PLS Analysis   
 
        %% read null spins;
        % null spins derived using Brainsmash
        % (https://brainsmash.readthedocs.io/en/latest/)
        null_lh_10000 = readtable('./files/surrogates_lh_10000.csv');
        null_lh_10000 =table2array(null_lh_10000);
        null_lh_10000 =null_lh_10000';
        
        %% yeo colormap for the 400x1 

        for j = 1:max(schaefer_400_label)    
            color_400(j,1) = mode(yeo_7_label_full(1,schaefer_400_label==j))
        end
        project_detection_community_boview(color_400','', schaefer_400_label',surf,1); colormap((yeo_colormap));
        
        %% parcellate pmap into schaefer 400 

        for i = 1:size(pmap_R_thal_ind_groupCmap_hcpd_align,1)
            for j = 1:max(schaefer_400_label)
            pmap_sch400_17nt.pmap_R_thal_ind_groupCmap_hcpd_p1{1,i}(j,:) = nanmean(pmap_R_thal_ind_groupCmap_hcpd_p1(i,schaefer_400_label == j),2 );    
            pmap_sch400_17nt.pmap_R_thal_ind_groupCmap_hcpd_p2{1,i}(j,:) = nanmean(pmap_R_thal_ind_groupCmap_hcpd_p2(i,schaefer_400_label == j),2 );    

            pmap_sch400_17nt.pmap_L_thal_ind_groupCmap_hcpd_p1{1,i}(j,:) = nanmean(pmap_L_thal_ind_groupCmap_hcpd_p1(i,schaefer_400_label == j),2 );    
            pmap_sch400_17nt.pmap_L_thal_ind_groupCmap_hcpd_p2{1,i}(j,:) = nanmean(pmap_L_thal_ind_groupCmap_hcpd_p2(i,schaefer_400_label == j),2 );    

        end        
        end
        for i=1:size(pmap_R_thal_ind_groupCmap_hcpd_align,1) %subject number

            %shaefer400
            sch400_pmap_R_thal_ind_groupCmap_hcpd_p1(i,:) = pmap_sch400_17nt.pmap_R_thal_ind_groupCmap_hcpd_p1{1,i};
            sch400_pmap_R_thal_ind_groupCmap_hcpd_p2(i,:) = pmap_sch400_17nt.pmap_R_thal_ind_groupCmap_hcpd_p2{1,i};

            sch400_pmap_L_thal_ind_groupCmap_hcpd_p1(i,:) = pmap_sch400_17nt.pmap_L_thal_ind_groupCmap_hcpd_p1{1,i};
            sch400_pmap_L_thal_ind_groupCmap_hcpd_p2(i,:) = pmap_sch400_17nt.pmap_L_thal_ind_groupCmap_hcpd_p2{1,i};

        end

        sch400_pmap_L_thal_ind_groupCmap_hcpd_p1_mean = zscore(nanmean(sch400_pmap_L_thal_ind_groupCmap_hcpd_p1,1))';
        project_detection_community_boview(sch400_pmap_L_thal_ind_groupCmap_hcpd_p1_mean,'', schaefer_400_label,surf,1);

        sch400_pmap_L_thal_ind_groupCmap_hcpd_p2_mean = zscore(nanmean(sch400_pmap_L_thal_ind_groupCmap_hcpd_p2,1))';
        project_detection_community_boview(sch400_pmap_L_thal_ind_groupCmap_hcpd_p2_mean,'', schaefer_400_label,surf,1);

        %% run PLS with permutation

        X = [CORTEX_CALB1_N2; CORTEX_PVALB_N2]';
        X_lh = X(1:200,:);
        X_rh = X(201:400,:);

        Y1 = [sch400_pmap_L_thal_ind_groupCmap_hcpd_p1_mean sch400_pmap_L_thal_ind_groupCmap_hcpd_p2_mean];
        Y1_lh=Y1(1:200,:);
        Y1_rh=Y1(201:400,:);

        % pls regression on original data
        %[Xl_pmap1,Yl_pmap1,Xs_pmap1,Ys_pmap1,beta_pmap1,pctVar_pmap1,mse_pmap1,stats_pmap1] = plsregress(X_lh,Y1_lh,1); 
        [Xl_pmap1,Yl_pmap1,Xs_pmap1,Ys_pmap1,beta_pmap1,pctVar_pmap1,mse_pmap1,stats_pmap1] = plsregress(X,Y1,1); 

        [corr_orig_pmap1 P]=corr(Xs_pmap1,Ys_pmap1)

        figure; plot(Xs_pmap1(:,1),Ys_pmap1(:,1),'.');
        %figure; gscatter(Xs_pmap1(:,1),Ys_pmap1(:,1),color_400(1:200,1), yeo_colormap(2:end, :), '.', 40); grid
        
        figure; gscatter(Xs_pmap1(:,1),Ys_pmap1(:,1),color_400(:,1), yeo_colormap(2:end, :), '.', 40); grid
        set (gcf, 'color', 'white', 'Position', [100 100 1400 900]);
        xlim([-0.15 0.15]); ylim([-30 30]); legend('off')    

        c = zeros(200, 3) 
        for i = 1:200
            for j = 1:7 
                if color_400(i,1)==j
                   c(i,:) = yeo_colormap(j+1,:)
                end

            end
        end

        %test significance using permutations   
        % controlling for spatial autocorrelation (permutation) 
        spins = null_lh_10000;          % spatial autocorrelation-preserving permutation assignments
        nspins = 10000;                 % number of permutations ("spins")

        % spin pmap %%%%%%%%
        for k = 1:nspins    
            Y1_lh_spin = Y1_lh(spins(:,k),:);  % permute pmap1

            [Xl_pmap1_sY,Yl_pmap1_sY,Xs_pmap1_sY,Ys_pmap1_sY,beta_pmap1_sY,pctVar_pmap1_sY,mse_pmap1_sY,stats_pmap1_sY] = plsregress(X_lh,Y1_lh_spin,1); %pls regression on permuted pmap
            null_corr_spinY_pmap1(:,:,k) = corr(Xs_pmap1_sY,Ys_pmap1_sY); % save correlation between LV scores
            null_Xl_pmap1_sY(:,:,k) = Xl_pmap1_sY; % save X loadings
            null_pctVar_pmap1_sY(:,:,k)=pctVar_pmap1_sY; % save percentage of variance explained in LV by each PLS component

        end

        figure; histogram(null_corr_spinY_pmap1);

        %calculate the probability that the observed correlation coeffcient is above the null distribution
        pval=length(find(abs(null_corr_spinY_pmap1)>abs(corr_orig_pmap1)))/nspins

%% -- Extended Figure 3 ----------------------------------------------
    for allAges=1
        % PMAP data 
        
        x_1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(:,:,1),1);       
        y_1 = nanmean(pmap_L_thal_ind_groupCmap_hcpd_align(:,:,2),1);       
        
        
        groups = yeo_7_label_full;
        
        % Scatter plot with different colors for each group
        figure('Position', [100, 100, 800, 800], 'Color', 'white'); % Adjust the values as needed
        scatter(x_1, y_1, 50, groups, 'filled', 'MarkerFaceAlpha', 0.1); colormap(yeo_colormap(:,:));
        set(gcf,'color','w');
        xlim([-1 1]); ylim([-1 1]);
        
        % Increase font size for tick labels on both axes
        set(gca, 'FontSize', 16); 
        % Labels and Title
        xlabel('NEOMAP 1', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('NEOMAP 2', 'FontSize', 18, 'FontWeight', 'bold');
        
        % Marking centroids
        uniqueGroups = unique(groups); % Get unique group values
        hold on;
        for i = 2:8
            group = uniqueGroups(i);
            groupPoints = [x_1(groups == group)', y_1(groups == group)'];
            centroid = mean(groupPoints);        
            plot(centroid(1), centroid(2), 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', [yeo_colormap(i,:)], 'MarkerEdgeColor', 'w' );
            centroid_networks_age1{i} = centroid;
            groupPoints_networks_age1{i} = groupPoints;

        end
        %hold off;
        %colorbar;

        % Add lines for top and bottom 10%
        x_1_sorted = sort(x_1);
        y_1_sorted = sort(y_1);
        n = sum(~isnan(x_1));
        top_10_percent_idx = round(n * 0.05);
        bottom_10_percent_idx = round(n * 0.95);

        top_10_x = x_1_sorted(top_10_percent_idx);
        bottom_10_x = x_1_sorted(bottom_10_percent_idx);

        top_10_y = y_1_sorted(top_10_percent_idx);
        bottom_10_y = y_1_sorted(bottom_10_percent_idx);

        % Vertical lines
        plot([top_10_x, top_10_x], ylim, 'k--', 'LineWidth', 2);
        plot([bottom_10_x, bottom_10_x], ylim, 'k--', 'LineWidth', 2);

        % Horizontal lines
        plot(xlim, [top_10_y, top_10_y], 'k--', 'LineWidth', 2);
        plot(xlim, [bottom_10_y, bottom_10_y], 'k--', 'LineWidth', 2);

        % colorbar;

        % visualize external networks' centroids   
        % vis som dan sal lim fpn dmn 
        centroid_networks_age1=centroid_networks_age1(:,2:end);
        external_3_age1 =cell2mat([centroid_networks_age1(1,1);centroid_networks_age1(1,2);centroid_networks_age1(1,3)])
        centroid_external_age1 = mean(external_3_age1);
        plot(centroid_external_age1(1), centroid_external_age1(2), 'Marker', '*', 'MarkerSize', 10, 'LineWidth', 2, 'Color', 'w');
        hold off;
        %plot(centroid_networks{7}(1,1), centroid_networks{7}(1,2), 'Marker', '*', 'MarkerSize', 5, 'LineWidth', 1, 'Color', 'b');

    end

    % Description of average pmap [top, bottom 10%] according to yeo 7 network
    for for_bar_graph = 1
        yeo_7_label_full_1 = length(find(yeo_7_label_full == 1)); 
        yeo_7_label_full_2 = length(find(yeo_7_label_full == 2)); 
        yeo_7_label_full_3 = length(find(yeo_7_label_full == 3)); 
        yeo_7_label_full_4 = length(find(yeo_7_label_full == 4)); 
        yeo_7_label_full_5 = length(find(yeo_7_label_full == 5)); 
        yeo_7_label_full_6 = length(find(yeo_7_label_full == 6)); 
        yeo_7_label_full_7 = length(find(yeo_7_label_full == 7)); 
        
        for for_pmap_L_thal_allAges = 1 
        %using normalized group map
        for for_pmap1_top10 = 1
            %pmap1 - top10%
            pmap_t1 = nanmean(pmap1_L_thal_ind_groupCmap_hcpd_align,1);
            pmap_t2=pmap_t1;
            pmap_t2(logical(~surf_roi))=[];

            [m n] = size(pmap_t2); 
            TN = m*n; TTP = TN*.10; TTP = ceil(TTP); 
            [val, ind] = sort(pmap_t2,'descend'); 
            mask_idx_pmap_1_L_thal_groupCmap_hcpd_top10=(ind(1:TTP));
            mask_val_pmap_1_L_thal_groupCmap_hcpd_top10=(val(1:TTP));

            mask_pmap_1_L_thal_groupCmap_hcpd_top10 = zeros(size(pmap_t2));
            mask_pmap_1_L_thal_groupCmap_hcpd_top10(1,mask_idx_pmap_1_L_thal_groupCmap_hcpd_top10) = 1;

            mask_pmap_1_L_thal_groupCmap_hcpd_top10_full = zeros(1,length(pmap_t1));
            mask_pmap_1_L_thal_groupCmap_hcpd_top10_full(logical(surf_roi))= mask_pmap_1_L_thal_groupCmap_hcpd_top10;
            figure; BoSurfStatView(mask_pmap_1_L_thal_groupCmap_hcpd_top10_full, surf); 

            mask_val_pmap_1_L_thal_groupCmap_hcpd_top10_1 =  zeros(size(pmap_t2));
            mask_val_pmap_1_L_thal_groupCmap_hcpd_top10_1(1,mask_idx_pmap_1_L_thal_groupCmap_hcpd_top10) = mask_val_pmap_1_L_thal_groupCmap_hcpd_top10;      

            mask_val_pmap_1_L_thal_groupCmap_hcpd_top10_full =  zeros(size(pmap_t1));
            mask_val_pmap_1_L_thal_groupCmap_hcpd_top10_full(logical(surf_roi)) = mask_val_pmap_1_L_thal_groupCmap_hcpd_top10_1;      
            figure; BoSurfStatView(mask_val_pmap_1_L_thal_groupCmap_hcpd_top10_full, surf); 

            yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full = yeo_7_label_full'.*mask_pmap_1_L_thal_groupCmap_hcpd_top10_full'
            figure; BoSurfStatView(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full, surf); colormap(yeo_colormap)
            yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full==0)=nan;
           % figure; histogram(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full); xlim([0, 8]); ylim([0, 2000]);

                    count_1= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full==1))/yeo_7_label_full_1
                    count_2= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full==2))/yeo_7_label_full_2
                    count_3= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full==3))/yeo_7_label_full_3
                    count_4= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full==4))/yeo_7_label_full_4
                    count_5= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full==5))/yeo_7_label_full_5
                    count_6= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full==6))/yeo_7_label_full_6
                    count_7= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full==7))/yeo_7_label_full_7
            y=[count_1 count_2 count_3 count_4 count_5 count_6 count_7]
            y_perc=y/sum(y)

            %calculate mean/standard error of actual gradient values
            means = zeros(1, 7);
            SEs = zeros(1, 7);
            for i = 1:7
                area = mask_val_pmap_1_L_thal_groupCmap_hcpd_top10_full(1,find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full==i));
                means(i) = mean(area);
                SEs(i) = std(area) / sqrt(length(area)); % Calculate SE
            end   

            %correct for number of vertices included for each network        
            means_corrected = means.*y_perc

            %sort according to values
            [sorted_means, sorted_indices] = sort(means_corrected);
            sorted_SEs = SEs(sorted_indices);

            % X-axis values (categories)
            x = 1:7;  % Assuming you have 7 categories
            sorted_x = x(sorted_indices);

            figure; set (gcf, 'color','w');
            b = bar(x, sorted_means,'FaceColor','flat');

            b.CData(1,:) = [yeo_colormap(1+sorted_x(1),:)];
            b.CData(2,:) = [yeo_colormap(1+sorted_x(2),:)];
            b.CData(3,:) = [yeo_colormap(1+sorted_x(3),:)];
            b.CData(4,:) = [yeo_colormap(1+sorted_x(4),:)];
            b.CData(5,:) = [yeo_colormap(1+sorted_x(5),:)];
            b.CData(6,:) = [yeo_colormap(1+sorted_x(6),:)];
            b.CData(7,:) = [yeo_colormap(1+sorted_x(7),:)];

            hold on;
            for i = 1:7
                errorbar(x(i), sorted_means(i), sorted_SEs(i), 'k', 'linestyle', 'none', 'linewidth', 1.5);
            end
            ylim([0 1]);

            hold off;
            means_corrected_pmap1_top10=means_corrected

        end
        for for_pmap1_bot10 = 1
            %pmap1 - bottom 10%
            [m n] = size(pmap_t2); 
            TN = m*n; TTP = TN*.10; TTP = ceil(TTP); 
            [val, ind] = sort(pmap_t2,'ascend'); 
            mask_idx_pmap_1_L_thal_groupCmap_hcpd_bot10=(ind(1:TTP));
            mask_val_pmap_1_L_thal_groupCmap_hcpd_bot10=(val(1:TTP));

            mask_pmap_1_L_thal_groupCmap_hcpd_bot10 = zeros(size(pmap_t2));
            mask_pmap_1_L_thal_groupCmap_hcpd_bot10(1,mask_idx_pmap_1_L_thal_groupCmap_hcpd_bot10) = 1;

            mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full = zeros(1,length(pmap_t1));
            mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full(logical(surf_roi))= mask_pmap_1_L_thal_groupCmap_hcpd_bot10;
            %figure; BoSurfStatView(mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full, surf); 

            mask_val_pmap_1_L_thal_groupCmap_hcpd_bot10_1 =  zeros(size(pmap_t2));
            mask_val_pmap_1_L_thal_groupCmap_hcpd_bot10_1(1,mask_idx_pmap_1_L_thal_groupCmap_hcpd_bot10) = mask_val_pmap_1_L_thal_groupCmap_hcpd_bot10;      

            mask_val_pmap_1_L_thal_groupCmap_hcpd_bot10_full =  zeros(size(pmap_t1));
            mask_val_pmap_1_L_thal_groupCmap_hcpd_bot10_full(logical(surf_roi)) = mask_val_pmap_1_L_thal_groupCmap_hcpd_bot10_1;      
            %figure; BoSurfStatView(mask_val_pmap_1_L_thal_groupCmap_hcpd_bot10_full, surf); 

            yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full = yeo_7_label_full'.*mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full'
            %figure; BoSurfStatView(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full, surf); colormap(yeo_colormap)
            yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full==0)=nan;
           % figure; histogram(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full); xlim([0, 8]); ylim([0, 2000]);

                    count_1= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full==1))/yeo_7_label_full_1
                    count_2= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full==2))/yeo_7_label_full_2
                    count_3= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full==3))/yeo_7_label_full_3
                    count_4= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full==4))/yeo_7_label_full_4
                    count_5= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full==5))/yeo_7_label_full_5
                    count_6= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full==6))/yeo_7_label_full_6
                    count_7= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full==7))/yeo_7_label_full_7
            y=[count_1 count_2 count_3 count_4 count_5 count_6 count_7]
            y_perc=y/sum(y)

            %calculate mean/standard error of actual gradient values
            means = zeros(1, 7);
            SEs = zeros(1, 7);
            for i = 1:7
                area = mask_val_pmap_1_L_thal_groupCmap_hcpd_bot10_full(1,find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full==i));
                means(i) = mean(area);
                SEs(i) = std(area) / sqrt(length(area)); % Calculate SE
            end   

            %correct for number of vertices included for each network        
            means_corrected = abs(means.*y_perc)

            %sort according to values
            [sorted_means, sorted_indices] = sort(means_corrected);
            sorted_SEs = SEs(sorted_indices);

            % X-axis values (categories)
            x = 1:7;  % Assuming you have 7 categories
            sorted_x = x(sorted_indices);

            figure; set (gcf, 'color','w');

            b = bar(x, sorted_means,'FaceColor','flat');

            b.CData(1,:) = [yeo_colormap(1+sorted_x(1),:)];
            b.CData(2,:) = [yeo_colormap(1+sorted_x(2),:)];
            b.CData(3,:) = [yeo_colormap(1+sorted_x(3),:)];
            b.CData(4,:) = [yeo_colormap(1+sorted_x(4),:)];
            b.CData(5,:) = [yeo_colormap(1+sorted_x(5),:)];
            b.CData(6,:) = [yeo_colormap(1+sorted_x(6),:)];
            b.CData(7,:) = [yeo_colormap(sorted_x(7),:)];

            hold on;
            for i = 1:7
                errorbar(x(i), sorted_means(i), sorted_SEs(i), 'k', 'linestyle', 'none', 'linewidth', 1.5);
            end
            ylim([0 1]);
            hold off;
            means_corrected_pmap1_bot10=means_corrected

        end
        for for_pmap2_top10 = 1
            %pmap2 - top10%
            pmap_t1 = nanmean(pmap2_L_thal_ind_groupCmap_hcpd_align,1);
            pmap_t2=pmap_t1;
            pmap_t2(logical(~surf_roi))=[];

            [m n] = size(pmap_t2); 
            TN = m*n; TTP = TN*.10; TTP = ceil(TTP); 
            [val, ind] = sort(pmap_t2,'descend'); 
            mask_idx_pmap_1_L_thal_groupCmap_hcpd_top10=(ind(1:TTP));
            mask_val_pmap_1_L_thal_groupCmap_hcpd_top10=(val(1:TTP));

            mask_pmap_1_L_thal_groupCmap_hcpd_top10 = zeros(size(pmap_t2));
            mask_pmap_1_L_thal_groupCmap_hcpd_top10(1,mask_idx_pmap_1_L_thal_groupCmap_hcpd_top10) = 1;

            mask_pmap_1_L_thal_groupCmap_hcpd_top10_full = zeros(1,length(pmap_t1));
            mask_pmap_1_L_thal_groupCmap_hcpd_top10_full(logical(surf_roi))= mask_pmap_1_L_thal_groupCmap_hcpd_top10;
            figure; BoSurfStatView(mask_pmap_1_L_thal_groupCmap_hcpd_top10_full, surf); 

            mask_val_pmap_1_L_thal_groupCmap_hcpd_top10_1 =  zeros(size(pmap_t2));
            mask_val_pmap_1_L_thal_groupCmap_hcpd_top10_1(1,mask_idx_pmap_1_L_thal_groupCmap_hcpd_top10) = mask_val_pmap_1_L_thal_groupCmap_hcpd_top10;      

            mask_val_pmap_1_L_thal_groupCmap_hcpd_top10_full =  zeros(size(pmap_t1));
            mask_val_pmap_1_L_thal_groupCmap_hcpd_top10_full(logical(surf_roi)) = mask_val_pmap_1_L_thal_groupCmap_hcpd_top10_1;      
            figure; BoSurfStatView(mask_val_pmap_1_L_thal_groupCmap_hcpd_top10_full, surf); 


            yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full = yeo_7_label_full'.*mask_pmap_1_L_thal_groupCmap_hcpd_top10_full'
            figure; BoSurfStatView(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full, surf); colormap(yeo_colormap)
            yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full==0)=nan;
           % figure; histogram(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full); xlim([0, 8]); ylim([0, 2000]);
                    count_1= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full==1))/yeo_7_label_full_1
                    count_2= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full==2))/yeo_7_label_full_2
                    count_3= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full==3))/yeo_7_label_full_3
                    count_4= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full==4))/yeo_7_label_full_4
                    count_5= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full==5))/yeo_7_label_full_5
                    count_6= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full==6))/yeo_7_label_full_6
                    count_7= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full==7))/yeo_7_label_full_7
            y=[count_1 count_2 count_3 count_4 count_5 count_6 count_7]
            y_perc=y/sum(y)

            %calculate mean/standard error of actual gradient values
            means = zeros(1, 7);
            SEs = zeros(1, 7);
            for i = 1:7
                area = mask_val_pmap_1_L_thal_groupCmap_hcpd_top10_full(1,find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_top10_full==i));
                means(i) = mean(area);
                SEs(i) = std(area) / sqrt(length(area)); % Calculate SE
            end   

            %correct for number of vertices included for each network        
            means_corrected = means.*y_perc

            %sort according to values
            [sorted_means, sorted_indices] = sort(means_corrected);
            sorted_SEs = SEs(sorted_indices);

            % X-axis values (categories)
            x = 1:7;  % Assuming you have 7 categories
            sorted_x = x(sorted_indices);

            figure; set (gcf, 'color','w');
            b = bar(x, sorted_means,'FaceColor','flat');

            b.CData(1,:) = [yeo_colormap(1+sorted_x(1),:)];
            b.CData(2,:) = [yeo_colormap(1+sorted_x(2),:)];
            b.CData(3,:) = [yeo_colormap(1+sorted_x(3),:)];
            b.CData(4,:) = [yeo_colormap(1+sorted_x(4),:)];
            b.CData(5,:) = [yeo_colormap(1+sorted_x(5),:)];
            b.CData(6,:) = [yeo_colormap(1+sorted_x(6),:)];
            b.CData(7,:) = [yeo_colormap(sorted_x(7),:)];

            hold on;
            for i = 1:7
                errorbar(x(i), sorted_means(i), sorted_SEs(i), 'k', 'linestyle', 'none', 'linewidth', 1.5);
            end
            ylim([0 1]);
            hold off;
            means_corrected_pmap2_top10=means_corrected

        end
        for for_pmap2_bot10 = 1
        %pmap1 - bottom 10%
        [m n] = size(pmap_t2); 
        TN = m*n; TTP = TN*.10; TTP = ceil(TTP); 
        [val, ind] = sort(pmap_t2,'ascend'); 
        mask_idx_pmap_1_L_thal_groupCmap_hcpd_bot10=(ind(1:TTP));
        mask_val_pmap_1_L_thal_groupCmap_hcpd_bot10=(val(1:TTP));

        mask_pmap_1_L_thal_groupCmap_hcpd_bot10 = zeros(size(pmap_t2));
        mask_pmap_1_L_thal_groupCmap_hcpd_bot10(1,mask_idx_pmap_1_L_thal_groupCmap_hcpd_bot10) = 1;

        mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full = zeros(1,length(pmap_t1));
        mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full(logical(surf_roi))= mask_pmap_1_L_thal_groupCmap_hcpd_bot10;
        figure; BoSurfStatView(mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full, surf); 

        mask_val_pmap_1_L_thal_groupCmap_hcpd_bot10_1 =  zeros(size(pmap_t2));
        mask_val_pmap_1_L_thal_groupCmap_hcpd_bot10_1(1,mask_idx_pmap_1_L_thal_groupCmap_hcpd_bot10) = mask_val_pmap_1_L_thal_groupCmap_hcpd_bot10;      

        mask_val_pmap_1_L_thal_groupCmap_hcpd_bot10_full =  zeros(size(pmap_t1));
        mask_val_pmap_1_L_thal_groupCmap_hcpd_bot10_full(logical(surf_roi)) = mask_val_pmap_1_L_thal_groupCmap_hcpd_bot10_1;      
        figure; BoSurfStatView(mask_val_pmap_1_L_thal_groupCmap_hcpd_bot10_full, surf); 

        yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full = yeo_7_label_full'.*mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full'
        figure; BoSurfStatView(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full, surf); colormap(yeo_colormap)
        yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full==0)=nan;
       % figure; histogram(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full); xlim([0, 8]); ylim([0, 2000]);

                count_1= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full==1))/yeo_7_label_full_1
                count_2= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full==2))/yeo_7_label_full_2
                count_3= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full==3))/yeo_7_label_full_3
                count_4= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full==4))/yeo_7_label_full_4
                count_5= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full==5))/yeo_7_label_full_5
                count_6= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full==6))/yeo_7_label_full_6
                count_7= length(find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full==7))/yeo_7_label_full_7
        y=[count_1 count_2 count_3 count_4 count_5 count_6 count_7]
        y_perc=y/sum(y)

        %calculate mean/standard error of actual gradient values
        means = zeros(1, 7);
        SEs = zeros(1, 7);
        for i = 1:7
            area = mask_val_pmap_1_L_thal_groupCmap_hcpd_bot10_full(1,find(yeo7_mask_pmap_1_L_thal_groupCmap_hcpd_bot10_full==i));
            means(i) = mean(area);
            SEs(i) = std(area) / sqrt(length(area)); % Calculate SE
        end   

        %correct for number of vertices included for each network        
        means_corrected = abs(means.*y_perc)

        %sort according to values
        [sorted_means, sorted_indices] = sort(means_corrected);
        sorted_SEs = SEs(sorted_indices);
        
        % X-axis values (categories)
        x = 1:7;  % Assuming you have 7 categories
        sorted_x = x(sorted_indices);

        figure; set (gcf, 'color','w');
        b = bar(x, sorted_means,'FaceColor','flat');
 
        b.CData(1,:) = [yeo_colormap(1+sorted_x(1),:)];
        b.CData(2,:) = [yeo_colormap(1+sorted_x(2),:)];
        b.CData(3,:) = [yeo_colormap(1+sorted_x(3),:)];
        b.CData(4,:) = [yeo_colormap(1+sorted_x(4),:)];
        b.CData(5,:) = [yeo_colormap(1+sorted_x(5),:)];
        b.CData(6,:) = [yeo_colormap(1+sorted_x(6),:)];
        %b.CData(7,:) = [yeo_colormap(sorted_x(7),:)];
        
        hold on;
        for i = 1:7
            errorbar(x(i), sorted_means(i), sorted_SEs(i), 'k', 'linestyle', 'none', 'linewidth', 1.5);
        end
        ylim([0 1]);
        hold off;
        means_corrected_pmap2_bot10=means_corrected;
        %         % X-axis values (categories)
%         x = 1:7;
% 
%         % Create the bar graph
%         figure;
%         b=bar(x, means_corrected,'FaceColor','flat');
%         b.CData(1,:) = [yeo_colormap(2,:)];
%         b.CData(2,:) = [yeo_colormap(3,:)];
%         b.CData(3,:) = [yeo_colormap(4,:)];
%         b.CData(4,:) = [yeo_colormap(5,:)];
%         b.CData(5,:) = [yeo_colormap(6,:)];
%         b.CData(6,:) = [yeo_colormap(7,:)];
%         b.CData(7,:) = [yeo_colormap(8,:)];
%         % Add error bars using standard errors
%         hold on;
%         for i = 1:7
%             errorbar(x(i), means_corrected(i), SEs(i), 'k', 'linestyle', 'none', 'linewidth', 1.5);
%         end
%         hold off;
    end
end  
    end
    
    % combination of bottom and top 10%: Pie Graph
    for for_pie_graph = 1
        means_corrected_pmap12_bot10 = nanmean([means_corrected_pmap1_bot10;means_corrected_pmap2_bot10])
        means_corrected_pmap12_top10 = nanmean([means_corrected_pmap1_top10;means_corrected_pmap2_top10])

        %bar graph
        figure; set (gcf, 'color','w');
        b=bar(x, means_corrected_pmap12_bot10,'FaceColor','flat');
        b.CData(1,:) = [yeo_colormap(2,:)];
        b.CData(2,:) = [yeo_colormap(3,:)];
        b.CData(3,:) = [yeo_colormap(4,:)];
        b.CData(4,:) = [yeo_colormap(5,:)];
        b.CData(5,:) = [yeo_colormap(6,:)];
        b.CData(6,:) = [yeo_colormap(7,:)];
        b.CData(7,:) = [yeo_colormap(8,:)];
        ylim([0 1.5]);

        figure; set (gcf, 'color','w');
        b=bar(x, means_corrected_pmap12_top10,'FaceColor','flat');
        b.CData(1,:) = [yeo_colormap(2,:)];
        b.CData(2,:) = [yeo_colormap(3,:)];
        b.CData(3,:) = [yeo_colormap(4,:)];
        b.CData(4,:) = [yeo_colormap(5,:)];
        b.CData(5,:) = [yeo_colormap(6,:)];
        b.CData(6,:) = [yeo_colormap(7,:)];
        b.CData(7,:) = [yeo_colormap(8,:)];
        ylim([0 1.5]);

        %pie graph
        data=means_corrected_pmap12_bot10;
        data(isnan(data))=0;
        labels={'VN', 'SN', 'DAN', 'SN', 'LN', 'FPN', 'DMN'};

        figure;set (gcf, 'color','w');
        h=pie(data)%, ones(size(data)));
        legend(labels, 'Location', 'BestOutside');
        title('Circle Graph with Percentages');
        % Set colors
        for k = 1:length(data)
            h(k*2-1).FaceColor = yeo_colormap(k+1,:);
        end


        %pie graph
        data=means_corrected_pmap12_top10;
        data(isnan(data))=0;
        labels={'VN', 'SN', 'DAN', 'SN', 'LN', 'FPN', 'DMN'};

        figure;set (gcf, 'color','w');
        h=pie(data)%, ones(size(data)));
        legend(labels, 'Location', 'BestOutside');
        title('Circle Graph with Percentages');
        % Set colors
        for k = 1:length(data)
        h(k*2-1).FaceColor = yeo_colormap(k+1,:);
    end
    end
    
%% --Extended Figure 4A ----------------------------------------------

yeo_parc = yeo_7_label_full;
sublist = idx_hcpd;

%rearrange gradients 
G1 = zeros(20484, 10, 603);
for i=1:length(idx_hcpd) 
    for k=1:10
        temp = pmap_L_thal_ind_groupCmap_hcpd_align(i,:,k);
        G1(:, k, i)=temp';
    end
end

for dispersion_silhouette_calculation = 1

    %for individuals: scatter plots     
    
    % Internal-external overlap of networks 
    nodes_ext = find(yeo_parc == 1 | yeo_parc == 2 | yeo_parc == 3);
    nodes_int = find(yeo_parc == 7 );
    nodes_sal = find(yeo_parc == 4);
    nodes_1 = find(yeo_parc == 1 );
    nodes_2 = find(yeo_parc == 2 );
    nodes_3 = find(yeo_parc == 3 );
    nodes_4 = find(yeo_parc == 4 );
    nodes_5 = find(yeo_parc == 5 );
    nodes_6 = find(yeo_parc == 6 );
    nodes_7 = find(yeo_parc == 7 );

    % Create a grouping vector for silhouette analysis
    yeo_idx_int_ext = [ones(length(nodes_int), 1); ones(length(nodes_ext), 1) * 2];
    yeo_idx_sal_ext = [ones(length(nodes_sal), 1); ones(length(nodes_ext), 1) * 2];
    yeo_idx_sal_int = [ones(length(nodes_sal), 1); ones(length(nodes_int), 1) * 2];
    yeo_idx_all = [ones(length(nodes_1), 1); ones(length(nodes_2), 1) * 2; ...
                ones(length(nodes_3), 1) * 3; ones(length(nodes_4), 1) * 4; ...
                ones(length(nodes_5), 1) * 5; ones(length(nodes_6), 1) * 6; ...
                ones(length(nodes_7), 1) * 7;];
                
    % Silhouette scores for Internal - External: NEOMAP 1& 2 separately              
    %takes a long time    
    for s = 1:length(sublist)

        % Calculate silhouette scores for the current subject: NEOMAP 1 & 2
        silhouette_scores_ind_int_ext = silhouette(G1([nodes_int'; nodes_ext'], [1:2], s), yeo_idx_int_ext);
        % Store the individual silhouette scores
        average_silhouette_score_ind_int_ext(:, s) = nanmean(silhouette_scores_ind_int_ext, 1);

        % Calculate silhouette scores for the current subject: NEOMAP1
        silhouette_scores_ind_int_ext = silhouette(G1([nodes_int'; nodes_ext'], 1, s), yeo_idx_int_ext);
        silhouette_scores_ind_sal_ext = silhouette(G1([nodes_sal'; nodes_ext'], 1, s), yeo_idx_sal_ext);
        silhouette_scores_ind_sal_int = silhouette(G1([nodes_sal'; nodes_int'], 1, s), yeo_idx_sal_int);
        silhouette_scores_ind_all = silhouette(G1([nodes_1'; nodes_2'; nodes_3'; nodes_4'; nodes_5'; nodes_6'; nodes_7'], 1, s), yeo_idx_all);

        % Store the individual silhouette scores
        average_silhouette_score_ind_int_ext_CMAP1(:, s) = nanmean(silhouette_scores_ind_int_ext, 1);
        average_silhouette_score_ind_sal_ext_CMAP1(:, s) = nanmean(silhouette_scores_ind_sal_ext, 1);
        average_silhouette_score_ind_sal_int_CMAP1(:, s) = nanmean(silhouette_scores_ind_sal_int, 1);
        average_silhouette_score_ind_all_CMAP1(:, s) = nanmean(silhouette_scores_ind_all, 1);
        
        % Calculate silhouette scores for the current subject: NEOMAP2
        silhouette_scores_ind_int_ext = silhouette(G1([nodes_int'; nodes_ext'], 2, s), yeo_idx_int_ext);
        silhouette_scores_ind_sal_ext = silhouette(G1([nodes_sal'; nodes_ext'], 2, s), yeo_idx_sal_ext);
        silhouette_scores_ind_sal_int = silhouette(G1([nodes_sal'; nodes_int'], 2, s), yeo_idx_sal_int);
        silhouette_scores_ind_all = silhouette(G1([nodes_1'; nodes_2'; nodes_3'; nodes_4'; nodes_5'; nodes_6'; nodes_7'], 2, s), yeo_idx_all);

        % Store the individual silhouette scores
        average_silhouette_score_ind_int_ext_CMAP2(:, s) = nanmean(silhouette_scores_ind_int_ext, 1);
        average_silhouette_score_ind_sal_ext_CMAP2(:, s) = nanmean(silhouette_scores_ind_sal_ext, 1);
        average_silhouette_score_ind_sal_int_CMAP2(:, s) = nanmean(silhouette_scores_ind_sal_int, 1);
        average_silhouette_score_ind_all_CMAP2(:, s) = nanmean(silhouette_scores_ind_all, 1);
        
    end    
    
    % figure visualization 
    for for_fig_int_ext=1
        x=age_y_hcpd;
        y=average_silhouette_score_ind_int_ext;
        [R P] = corr(x,y')
        % Create a new figure with custom dimensions (width x height)
        figure('Position', [100, 100, 800, 600], 'Color', 'white'); 
        
        % Create a scatter plot
        scatter(x, y, 15, 'filled', 'MarkerFaceColor', 'k');
        ylim([-0.2, 0.6]); xlim([7.5, 22.5]);

        hold on;
        % Fit a regression line (e.g., linear regression)
        coefficients = polyfit(x, y, 1); 

        % Generate points for the regression line
        x_fit = min(x):0.1:max(x); 
        y_fit = polyval(coefficients, x_fit);

        % Create a thick regression line
        plot(x_fit, y_fit, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');

        % Increase font size for tick labels on both axes
        set(gca, 'FontSize', 16);

        % Labels and Title
        xlabel('Age (Year)', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('Silhouette score: External-Internal', 'FontSize', 18, 'FontWeight', 'bold');

        
    end    
    for for_fig_sal_ext_NEOMAP1=1
        x=age_y_hcpd;
        y=average_silhouette_score_ind_sal_ext_CMAP1;
        [R P] = corr(x,y')
        % Create a new figure with custom dimensions (width x height)
        figure('Position', [100, 100, 800, 600], 'Color', 'white');
        
        % Create a scatter plot
        scatter(x, y, 15, 'filled', 'MarkerFaceColor', 'k');
        ylim([-0.2, 0.6]); xlim([7.5, 22.5]);

        hold on;
        % Fit a regression line (e.g., linear regression)
        coefficients = polyfit(x, y, 1); % 1 for linear regression

        % Generate points for the regression line
        x_fit = min(x):0.1:max(x); 
        y_fit = polyval(coefficients, x_fit);

        % Create a thick regression line
        plot(x_fit, y_fit, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');

        % Increase font size for tick labels on both axes
        set(gca, 'FontSize', 16);

        % Labels and Title
        xlabel('Age (Year)', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('Silhouette score: Salience-External', 'FontSize', 18, 'FontWeight', 'bold');
        
        
    end
    for for_fig_sal_int_NEOMAP2=1
        x=age_y_hcpd;
        y=average_silhouette_score_ind_sal_int_CMAP2;
        [R P] = corr(x,y')
        % Create a new figure with custom dimensions (width x height)
        figure('Position', [100, 100, 800, 600], 'Color', 'white'); 
        
        % Create a scatter plot
        scatter(x, y, 15, 'filled', 'MarkerFaceColor', 'k');
        ylim([-0.2, 0.6]); xlim([7.5, 22.5]);

        hold on;
        % Fit a regression line (e.g., linear regression)
        coefficients = polyfit(x, y, 1);

        % Generate points for the regression line
        x_fit = min(x):0.1:max(x); 
        y_fit = polyval(coefficients, x_fit);

        % Create a thick regression line
        plot(x_fit, y_fit, 'LineWidth', 2, 'Color', 'r', 'LineStyle', '--');

        % Increase font size for tick labels on both axes
        set(gca, 'FontSize', 16); 

        % Labels and Title
        xlabel('Age (Year)', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('Silhouette score: Salience-Internal', 'FontSize', 18, 'FontWeight', 'bold');
        
        
    end
    
end

%% --Extended Figure 4B ----------------------------------------------
    for ageGroup1=1
        % PMAP data         
        x_1 = mean(pmap_L_thal_ind_groupCmap_hcpd_align(find(age_y_hcpd<12),:,1));
        y_1 = mean(pmap_L_thal_ind_groupCmap_hcpd_align(find(age_y_hcpd<12),:,2));
        groups = yeo_7_label_full;
        
        % Scatter plot with different colors for each group
        figure('Position', [100, 100, 800, 600], 'Color', 'white'); 
        scatter(x_1, y_1, 50, groups, 'filled', 'MarkerFaceAlpha', 0.1); colormap(yeo_colormap(:,:));
        set(gcf,'color','w');
        xlim([-1 1]); ylim([-1 1]);
        
        % Increase font size for tick labels on both axes
        set(gca, 'FontSize', 16); 
        % Labels and Title
        xlabel('NEOMAP 1', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('NEOMAP 2', 'FontSize', 18, 'FontWeight', 'bold');
        
        % Marking centroids
        uniqueGroups = unique(groups); 
        hold on;
        for i = 2:8
            group = uniqueGroups(i);
            groupPoints = [x_1(groups == group)', y_1(groups == group)'];
            centroid = mean(groupPoints);        
            plot(centroid(1), centroid(2), 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', [yeo_colormap(i,:)], 'MarkerEdgeColor', 'w' );
            centroid_networks_age1{i} = centroid;
            groupPoints_networks_age1{i} = groupPoints;

        end

        % visualize external networks' centroids   
        % vis som dan sal lim fpn dmn 
        centroid_networks_age1=centroid_networks_age1(:,2:end);
        external_3_age1 =cell2mat([centroid_networks_age1(1,1);centroid_networks_age1(1,2);centroid_networks_age1(1,3)])
        centroid_external_age1 = mean(external_3_age1);
        plot(centroid_external_age1(1), centroid_external_age1(2), 'Marker', '*', 'MarkerSize', 10, 'LineWidth', 2, 'Color', 'w');
        hold off;

        % calculate between-centroids distance
        centroid_salience_age1= [centroid_networks_age1{4}(1,1), centroid_networks_age1{4}(1,2)]
        centroid_dmn_age1= [centroid_networks_age1{7}(1,1), centroid_networks_age1{7}(1,2)]

        net_dist_salint_age1 = pdist([centroid_salience_age1; centroid_dmn_age1])
        net_dist_salext_age1 = pdist([centroid_salience_age1; centroid_external_age1])

        % calculate intersection between external=internal 
        groupPoints_networks_age1=groupPoints_networks_age1(:,2:end);

        data_external_age1= cell2mat([groupPoints_networks_age1(1,1);groupPoints_networks_age1(1,2);groupPoints_networks_age1(1,3)]);
        data_internal_age1= cell2mat([groupPoints_networks_age1(1,7)]);

        in_poly_1 = inpolygon(data_external_age1(:,1), data_external_age1(:,2), data_internal_age1(:,1), data_internal_age1(:,2)); % Check if points from data1 fall within scatter plot 2
        in_poly_2 = inpolygon(data_internal_age1(:,1), data_internal_age1(:,2), data_external_age1(:,1), data_external_age1(:,2)); % Check if points from data2 fall within scatter plot 1

        intersection_age1 = sum(in_poly_1) + sum(in_poly_2); % Calculate the total number of points in the intersection

        % Calculate overlap percentage
        overlap_percentage_age1 = (intersection_age1 / (length(data_external_age1) + length(data_internal_age1))) * 100; % Convert intersection to percentage

        % Display intersection and overlap percentage
        disp(['Intersection: ', num2str(intersection_age1)]);
        disp(['Overlap percentage: ', num2str(overlap_percentage_age1), '%']);




    end
    for ageGroup2=1
        % PMAP data        
        x_1 = mean(pmap_L_thal_ind_groupCmap_hcpd_align(find(age_y_hcpd>=12 & age_y_hcpd<18),:,1));
        y_1 = mean(pmap_L_thal_ind_groupCmap_hcpd_align(find(age_y_hcpd>=12 & age_y_hcpd<18),:,2));
        groups = yeo_7_label_full;
        
        % Scatter plot with different colors for each group
        figure('Position', [100, 100, 800, 600], 'Color', 'white'); 
        scatter(x_1, y_1, 50, groups, 'filled', 'MarkerFaceAlpha', 0.1); colormap(yeo_colormap(:,:));
        set(gcf,'color','w');
        xlim([-1 1]); ylim([-1 1]);

        % Increase font size for tick labels on both axes
        set(gca, 'FontSize', 16); 
        % Labels and Title
        xlabel('NEOMAP 1', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('NEOMAP 2', 'FontSize', 18, 'FontWeight', 'bold');
        
        % Marking centroids
        uniqueGroups = unique(groups); 
        hold on;
        for i = 2:8
            group = uniqueGroups(i);
            groupPoints = [x_1(groups == group)', y_1(groups == group)'];
            centroid = mean(groupPoints);        
            plot(centroid(1), centroid(2), 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', [yeo_colormap(i,:)], 'MarkerEdgeColor', 'w' );
            centroid_networks_age1{i} = centroid;
            groupPoints_networks_age1{i} = groupPoints;

        end

        % visualize external networks' centroids   
        % vis som dan sal lim fpn dmn 
        centroid_networks_age1=centroid_networks_age1(:,2:end);
        external_3_age1 =cell2mat([centroid_networks_age1(1,1);centroid_networks_age1(1,2);centroid_networks_age1(1,3)])
        centroid_external_age1 = mean(external_3_age1);
        plot(centroid_external_age1(1), centroid_external_age1(2), 'Marker', '*', 'MarkerSize', 10, 'LineWidth', 2, 'Color', 'w');
        hold off;

        % calculate between-centroids distance
        centroid_salience_age1= [centroid_networks_age1{4}(1,1), centroid_networks_age1{4}(1,2)]
        centroid_dmn_age1= [centroid_networks_age1{7}(1,1), centroid_networks_age1{7}(1,2)]

        net_dist_salint_age1 = pdist([centroid_salience_age1; centroid_dmn_age1])
        net_dist_salext_age1 = pdist([centroid_salience_age1; centroid_external_age1])

        % calculate intersection between external=internal 
        groupPoints_networks_age1=groupPoints_networks_age1(:,2:end);

        data_external_age1= cell2mat([groupPoints_networks_age1(1,1);groupPoints_networks_age1(1,2);groupPoints_networks_age1(1,3)]);
        data_internal_age1= cell2mat([groupPoints_networks_age1(1,7)]);

        in_poly_1 = inpolygon(data_external_age1(:,1), data_external_age1(:,2), data_internal_age1(:,1), data_internal_age1(:,2)); % Check if points from data1 fall within scatter plot 2
        in_poly_2 = inpolygon(data_internal_age1(:,1), data_internal_age1(:,2), data_external_age1(:,1), data_external_age1(:,2)); % Check if points from data2 fall within scatter plot 1

        intersection_age1 = sum(in_poly_1) + sum(in_poly_2); % Calculate the total number of points in the intersection

        % Calculate overlap percentage
        overlap_percentage_age1 = (intersection_age1 / (length(data_external_age1) + length(data_internal_age1))) * 100; % Convert intersection to percentage

        % Display intersection and overlap percentage
        disp(['Intersection: ', num2str(intersection_age1)]);
        disp(['Overlap percentage: ', num2str(overlap_percentage_age1), '%']);




    end
    for ageGroup3=1
        % PMAP data   
        x_1 = mean(pmap_L_thal_ind_groupCmap_hcpd_align(find(age_y_hcpd>=18),:,1));
        y_1 = mean(pmap_L_thal_ind_groupCmap_hcpd_align(find(age_y_hcpd>=18),:,2));

        groups = yeo_7_label_full;
        
        % Scatter plot with different colors for each group
        figure('Position', [100, 100, 800, 600], 'Color', 'white'); 
        scatter(x_1, y_1, 50, groups, 'filled', 'MarkerFaceAlpha', 0.1); colormap(yeo_colormap(:,:));
        set(gcf,'color','w');
        xlim([-1 1]); ylim([-1 1]);

        % Increase font size for tick labels on both axes
        set(gca, 'FontSize', 16); 
        % Labels and Title
        xlabel('NEOMAP 1', 'FontSize', 18, 'FontWeight', 'bold');
        ylabel('NEOMAP 2', 'FontSize', 18, 'FontWeight', 'bold');
        
        % Marking centroids
        uniqueGroups = unique(groups); % Get unique group values
        hold on;
        for i = 2:8
            group = uniqueGroups(i);
            groupPoints = [x_1(groups == group)', y_1(groups == group)'];
            centroid = mean(groupPoints);        
            plot(centroid(1), centroid(2), 'Marker', 'o', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', [yeo_colormap(i,:)], 'MarkerEdgeColor', 'w' );
            centroid_networks_age1{i} = centroid;
            groupPoints_networks_age1{i} = groupPoints;

        end

        % visualize external networks' centroids   
        % vis som dan sal lim fpn dmn 
        centroid_networks_age1=centroid_networks_age1(:,2:end);
        external_3_age1 =cell2mat([centroid_networks_age1(1,1);centroid_networks_age1(1,2);centroid_networks_age1(1,3)])
        centroid_external_age1 = mean(external_3_age1);
        plot(centroid_external_age1(1), centroid_external_age1(2), 'Marker', '*', 'MarkerSize', 10, 'LineWidth', 2, 'Color', 'w');
        hold off;

        % calculate between-centroids distance
        centroid_salience_age1= [centroid_networks_age1{4}(1,1), centroid_networks_age1{4}(1,2)]
        centroid_dmn_age1= [centroid_networks_age1{7}(1,1), centroid_networks_age1{7}(1,2)]

        net_dist_salint_age1 = pdist([centroid_salience_age1; centroid_dmn_age1])
        net_dist_salext_age1 = pdist([centroid_salience_age1; centroid_external_age1])

        % calculate intersection between external=internal 
        groupPoints_networks_age1=groupPoints_networks_age1(:,2:end);

        data_external_age1= cell2mat([groupPoints_networks_age1(1,1);groupPoints_networks_age1(1,2);groupPoints_networks_age1(1,3)]);
        data_internal_age1= cell2mat([groupPoints_networks_age1(1,7)]);

        in_poly_1 = inpolygon(data_external_age1(:,1), data_external_age1(:,2), data_internal_age1(:,1), data_internal_age1(:,2)); % Check if points from data1 fall within scatter plot 2
        in_poly_2 = inpolygon(data_internal_age1(:,1), data_internal_age1(:,2), data_external_age1(:,1), data_external_age1(:,2)); % Check if points from data2 fall within scatter plot 1

        intersection_age1 = sum(in_poly_1) + sum(in_poly_2); % Calculate the total number of points in the intersection

        % Calculate overlap percentage
        overlap_percentage_age1 = (intersection_age1 / (length(data_external_age1) + length(data_internal_age1))) * 100; % Convert intersection to percentage

        % Display intersection and overlap percentage
        disp(['Intersection: ', num2str(intersection_age1)]);
        disp(['Overlap percentage: ', num2str(overlap_percentage_age1), '%']);


    end

%% --Extended Figure 4C ----------------------------------------------

    for for_compare_with_Infant = 1
        load('./files/infant_average_silhouette_score_int_ext_wang10.mat')
        group1= average_silhouette_score_ind_int_ext_infant;
        group2= average_silhouette_score_ind_int_ext;

        [h, p, ci, stats] = ttest2(group1, group2);

        figure;
        histogram(group1, 'Normalization', 'probability', 'FaceColor', 'b', 'FaceAlpha', 0.5);
        hold on;
        histogram(group2, 'Normalization', 'probability', 'FaceColor', 'r', 'FaceAlpha', 0.5);
        xlabel('Value');
        ylabel('Probability');
        legend('dHCP', 'HCPD');%grid on;
        %title('Histogram of individual silhouette scores: Internal - External');
        set(gcf,'color','w')
        hold off;                                           
    end

