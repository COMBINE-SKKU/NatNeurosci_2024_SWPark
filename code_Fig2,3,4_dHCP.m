
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
        
    end
    %% 02) read phenotypic info
    for readdata=1
        data = importfile_demo_dHCP(['./files/demo_dHCP_rev.xlsx']);

    end    
    for demo_vars = 1

        scanage = data.scanage
        birthage = data.birthage
        sex = data.sex
        subid = data.sub_id
        sesid = data.ses_id
        meanFD = data.meanFD
        meanFD04 = data.MeanFD04
        preproc = data.preproc
        study = data.dataset
        qc_grad = data.qc_grad

        idx_infant = find( (study == 'dhcp') & (preproc == 1) & (qc_grad == 1) & (meanFD04 == 1));

        age_infant = scanage(idx_infant);
        sex_infant = sex(idx_infant);
        sex_infant = cellstr(sex_infant);        
        subid_infant = subid(idx_infant);
        sesid_infant = sesid(idx_infant);
        meanFD_infant = meanFD(idx_infant);
        age_birth_infant = birthage(idx_infant);

        idx_infantpt = find (age_infant < 37) ;  
        idx_infant373839 = find (age_infant >= 37 & age_infant < 40) ;  
        idx_infant40 = find (age_infant >= 40 & age_infant < 41) ;  
        idx_infant41 = find (age_infant >= 41 & age_infant < 42) ;  
        idx_infant424344 = find (age_infant >= 42) ;  

    end    
    %% 03) Load left and right surfaces & parcellations

    for load_surfaces = 1

    %upload 10k surfaces: dhcpSym40 space
    temp = gifti(['./files/week-40_hemi-left_space-dhcpSym_dens-10k_midthickness.surf.gii']);
    surfL_infant.coord = temp.vertices';     
    surfL_infant.tri   = temp.faces;   

    temp = gifti(['./files/week-40_hemi-right_space-dhcpSym_dens-10k_midthickness.surf.gii']);
    surfR_infant.coord = temp.vertices';     
    surfR_infant.tri   = temp.faces;   

    surf_infant.coord = [ surfL_infant.coord surfR_infant.coord ];
    surf_infant.tri   = [ surfL_infant.tri; surfR_infant.tri+10242; ];       

    temp = gifti([ './files/week-40_hemi-left_space-dhcpSym_dens-10k_vinflated.surf.gii']);
    surfL_infant_vi.coord = temp.vertices';     
    surfL_infant_vi.tri   = temp.faces;   

    temp = gifti([ './files/week-40_hemi-right_space-dhcpSym_dens-10k_vinflated.surf.gii']);
    surfR_infant_vi.coord = temp.vertices';     
    surfR_infant_vi.tri   = temp.faces;   

    surf_infant_vi.coord = [ surfL_infant_vi.coord surfR_infant_vi.coord ];
    surf_infant_vi.tri   = [ surfL_infant_vi.tri; surfR_infant_vi.tri+10242; ];       


    temp = gifti([ './files/week-40_hemi-left_space-dhcpSym_dens-10k_desc-medialwall_mask.shape.gii']);
    surfL_roi_infant = temp.cdata';  
    surfL_roi_infant_t = surfL_roi_infant;   
    surfL_roi_infant_t(:, surfL_roi_infant_t ~= 0 )=1; 

    temp = gifti([ './files/week-40_hemi-right_space-dhcpSym_dens-10k_desc-medialwall_mask.shape.gii']);
    surfR_roi_infant = temp.cdata'; 
    surfR_roi_infant_t = surfR_roi_infant;   
    surfR_roi_infant_t(:, surfR_roi_infant_t ~= 0 )=1; 

    surf_roi_infant = cat(2, surfL_roi_infant_t, surfR_roi_infant_t);
    figure; SurfStatView(surf_roi_infant,surf_infant);      

    temp = gifti(['./files/L.atlasroi.10k_fs_LR.shape.gii']);
    surfL_roi = temp.cdata';     
    temp = gifti(['./files/R.atlasroi.10k_fs_LR.shape.gii']);
    surfR_roi = temp.cdata'; 

    surf_roi = cat(2,surfL_roi, surfR_roi);
    figure; SurfStatView(surf_roi,surf_infant);        

    for parcel_10k = 1
        for load_schaefers_400_parcel = 1
            schaefer_400_10k = ciftiopen([ './files/Schaefer2018_400Parcels_7Networks_order_10k.dlabel.nii'],  WB_COMMAND); 
            schaefer_400_label = schaefer_400_10k.cdata;   

            schaefer_400_full = zeros(20482,1); 
            schaefer_400_full(logical(surf_roi)) = schaefer_400_label;

            figure; SurfStatView(schaefer_400_full,surf_infant);
            schaefer_400_label = schaefer_400_full;                    


        end      
        for load_yeo_7networks_cortex = 1
           cii_yeo = ciftiopen( './files/RSN-networks.10k_fs_LR.dlabel.nii',  WB_COMMAND);
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
            figure; SurfStatView(yeo_7_label_full,surf_infant);

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

%% -- Figure 2 ------------------------------------------

    %% upload source data for NEOMAPs
    load('sourceData_Fig2,3,4_dHCP_NEOMAPs.mat') 
    
    pmap_R_thal_ind_groupCmap_infant_p1 = pmap_R_thal_ind_groupCmap_infant_align(:,:,1);
    pmap_R_thal_ind_groupCmap_infant_p2 = pmap_R_thal_ind_groupCmap_infant_align(:,:,2);

    pmap_L_thal_ind_groupCmap_infant_p1 = pmap_L_thal_ind_groupCmap_infant_align(:,:,1);
    pmap_L_thal_ind_groupCmap_infant_p2 = pmap_L_thal_ind_groupCmap_infant_align(:,:,2);
    
    %% boxplot for pmap values: all values
    for for_yeo7_networks_pmap1 = 1
        yeo_colormap_1=yeo_colormap(2:end,:)
        %%pmap1
        pmap_t1 = nanmean(pmap_L_thal_ind_groupCmap_infant_align(:,:,1),1);
        figure; SurfStatView (mean(pmap_t1,1),surf_infant)
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
        %figure; SurfStatView(wang_label_LR,surf_infant);colormap(my_colormap)

        % Create a figure and boxplot with sorted groups
        figure;
        boxplot(pmap, labels, 'GroupOrder', string(sortOrder), 'Colors', sortedColors, 'Symbol', '');
        hold on;

        % Plot a colored scatter plot for sorted groups
        for i = 1:length(sortedLabels)
            label = sortedLabels(i);
            mask = labels == label;
            scatter(i * ones(1, sum(mask)), pmap(mask), 1, sortedColors(i, :), 'filled', 'MarkerFaceAlpha', 0.6, 'jitter', 'on', 'jitterAmount', 0.1);
        end
        set (gcf,'color','w');
        ylim([-1 1])

        hold off;

    end
    for for_yeo7_networks_pmap2 = 1
        %%pmap2
        pmap_t2 = nanmean(pmap_L_thal_ind_groupCmap_infant_align(:,:,2),1);
        labels = yeo_7_label_full      
        idx=find((labels>0));
        pmap = pmap_t2(1,idx)
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
        hold on;

        % Plot a colored scatter plot for sorted groups
        for i = 1:length(sortedLabels)
            label = sortedLabels(i);
            mask = labels == label;
            scatter(i * ones(1, sum(mask)), pmap(mask), 1, sortedColors(i, :), 'filled', 'MarkerFaceAlpha', 0.6, 'jitter', 'on', 'jitterAmount', 0.1);
        end
        set (gcf,'color','w');
        ylim([-1 1])

        hold off;

    end

    %% upload group CMAP -> DONT' NEED??!?!?!
    
    % you need to change the P variable to your path
    PATH               = '/local_raid2/03_user/shinwon/01_project/01_thalamusgradient/';
    GROUP_GRADIENTS           = [ PATH '04_group_gradient/' ];

    cii = ciftiopen([ GROUP_GRADIENTS '08_dHCP_rev/func_clean_10k_dhcpSym40/standard_scmask_LRthal_1.5mm_label.2.cmap.dscalar.nii' ], WB_COMMAND);   
    cii_L_thal = cii.cdata(cii.diminfo{1,1}.models{4,1}.start:(cii.diminfo{1,1}.models{4,1}.start-1) + cii.diminfo{1,1}.models{4,1}.count,:);
        %timeseries_L_thal(:,:,i) = cii_L_thal;        
    groupGradient_L_infant(:, :) = cii_L_thal;

    cii = ciftiopen([ GROUP_GRADIENTS '08_dHCP_rev/func_clean_10k_dhcpSym40/standard_scmask_LRthal_1.5mm_label.3.cmap.dscalar.nii' ], WB_COMMAND);   
    cii_R_thal = cii.cdata(cii.diminfo{1,1}.models{5,1}.start:(cii.diminfo{1,1}.models{5,1}.start-1) + cii.diminfo{1,1}.models{5,1}.count,:);
        %timeseries_L_thal(:,:,i) = cii_L_thal;        
    groupGradient_R_infant(:, :) = cii_R_thal;
    
%% -- Figure 3A ------------------------------------------
    %% difference bt SALIENCE(VAN) & DAN/DMN areas 
    idx_sensory = [find(yeo_7_label_full==1) , find(yeo_7_label_full==2)]';
    
    idx_dan = [find(yeo_7_label_full==3) , find(yeo_7_label_full==1) , find(yeo_7_label_full==2)]';
    idx_sal = [find(yeo_7_label_full==4)]';
    idx_dmn = [find(yeo_7_label_full==7)]';
    
    %aligned pmaps: group mean
    for for_pmap1=1 

        %compute mean of each networks of interest
        mean_dan_infantpt_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_infant_align(idx_infantpt,idx_dan,1),2);
        mean_sal_infantpt_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_infant_align(idx_infantpt,idx_sal,1),2);
        mean_dmn_infantpt_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_infant_align(idx_infantpt,idx_dmn,1),2);

        mean_dan_infant373839_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant373839,idx_dan,1),2);
        mean_sal_infant373839_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant373839,idx_sal,1),2);
        mean_dmn_infant373839_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant373839,idx_dmn,1),2);

        mean_dan_infant40_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant40,idx_dan,1),2);
        mean_sal_infant40_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant40,idx_sal,1),2);
        mean_dmn_infant40_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant40,idx_dmn,1),2);

        mean_dan_infant41_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant41,idx_dan,1),2);
        mean_sal_infant41_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant41,idx_sal,1),2);
        mean_dmn_infant41_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant41,idx_dmn,1),2);

        mean_dan_infant424344_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant424344,idx_dan,1),2);
        mean_sal_infant424344_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant424344,idx_sal,1),2);
        mean_dmn_infant424344_pmap1 = nanmean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant424344,idx_dmn,1),2);

        %diff between sal & dan
        diff_sal_dan_infantpt_pmap1 = abs(mean_sal_infantpt_pmap1 - mean_dan_infantpt_pmap1);
        diff_sal_dan_infant373839_pmap1 = abs(mean_sal_infant373839_pmap1 - mean_dan_infant373839_pmap1);
        diff_sal_dan_infant40_pmap1 = abs(mean_sal_infant40_pmap1 - mean_dan_infant40_pmap1);
        diff_sal_dan_infant41_pmap1 = abs(mean_sal_infant41_pmap1 - mean_dan_infant41_pmap1);
        diff_sal_dan_infant424344_pmap1 = abs(mean_sal_infant424344_pmap1 - mean_dan_infant424344_pmap1);
        
        
        %diff between sal & dmn
        diff_sal_dmn_infantpt_pmap1 = abs(mean_sal_infantpt_pmap1 - mean_dmn_infantpt_pmap1);
        diff_sal_dmn_infant373839_pmap1 = abs(mean_sal_infant373839_pmap1 - mean_dmn_infant373839_pmap1);
        diff_sal_dmn_infant40_pmap1 = abs(mean_sal_infant40_pmap1 - mean_dmn_infant40_pmap1);
        diff_sal_dmn_infant41_pmap1 = abs(mean_sal_infant41_pmap1 - mean_dmn_infant41_pmap1);
        diff_sal_dmn_infant424344_pmap1 = abs(mean_sal_infant424344_pmap1 - mean_dmn_infant424344_pmap1);

        
    end
    for for_pmap2=1 

        %compute mean of each networks of interest
        mean_dan_infantpt_pmap2 = mean(pmap_L_thal_ind_groupCmap_infant_align(idx_infantpt,idx_dan,2),2);
        mean_sal_infantpt_pmap2 = mean(pmap_L_thal_ind_groupCmap_infant_align(idx_infantpt,idx_sal,2),2);
        mean_dmn_infantpt_pmap2 = mean(pmap_L_thal_ind_groupCmap_infant_align(idx_infantpt,idx_dmn,2),2);

        mean_dan_infant373839_pmap2 = mean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant373839,idx_dan,2),2);
        mean_sal_infant373839_pmap2 = mean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant373839,idx_sal,2),2);
        mean_dmn_infant373839_pmap2 = mean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant373839,idx_dmn,2),2);

        mean_dan_infant40_pmap2 = mean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant40,idx_dan,2),2);
        mean_sal_infant40_pmap2 = mean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant40,idx_sal,2),2);
        mean_dmn_infant40_pmap2 = mean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant40,idx_dmn,2),2);

        mean_dan_infant41_pmap2 = mean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant41,idx_dan,2),2);
        mean_sal_infant41_pmap2 = mean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant41,idx_sal,2),2);
        mean_dmn_infant41_pmap2 = mean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant41,idx_dmn,2),2);

        mean_dan_infant424344_pmap2 = mean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant424344,idx_dan,2),2);
        mean_sal_infant424344_pmap2 = mean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant424344,idx_sal,2),2);
        mean_dmn_infant424344_pmap2 = mean(pmap_L_thal_ind_groupCmap_infant_align(idx_infant424344,idx_dmn,2),2);

        %diff between sal & dan
        diff_sal_dan_infantpt_pmap2 = abs(mean_sal_infantpt_pmap2 - mean_dan_infantpt_pmap2);
        diff_sal_dan_infant373839_pmap2 = abs(mean_sal_infant373839_pmap2 - mean_dan_infant373839_pmap2);
        diff_sal_dan_infant40_pmap2 = abs(mean_sal_infant40_pmap2 - mean_dan_infant40_pmap2);
        diff_sal_dan_infant41_pmap2 = abs(mean_sal_infant41_pmap2 - mean_dan_infant41_pmap2);
        diff_sal_dan_infant424344_pmap2 = abs(mean_sal_infant424344_pmap2 - mean_dan_infant424344_pmap2);

        %diff between sal & dmn
        diff_sal_dmn_infantpt_pmap2 = abs(mean_sal_infantpt_pmap2 - mean_dmn_infantpt_pmap2);
        diff_sal_dmn_infant373839_pmap2 = abs(mean_sal_infant373839_pmap2 - mean_dmn_infant373839_pmap2);
        diff_sal_dmn_infant40_pmap2 = abs(mean_sal_infant40_pmap2 - mean_dmn_infant40_pmap2);
        diff_sal_dmn_infant41_pmap2 = abs(mean_sal_infant41_pmap2 - mean_dmn_infant41_pmap2);
        diff_sal_dmn_infant424344_pmap2 = abs(mean_sal_infant424344_pmap2 - mean_dmn_infant424344_pmap2);

        
    end
    %% boxplots 
    
    for for_pmap1 = 1

        %combined         
        test = [diff_sal_dan_infantpt_pmap1; diff_sal_dmn_infantpt_pmap1; diff_sal_dan_infant373839_pmap1; diff_sal_dmn_infant373839_pmap1; ...
            diff_sal_dan_infant40_pmap1; diff_sal_dmn_infant40_pmap1; diff_sal_dan_infant41_pmap1; diff_sal_dmn_infant41_pmap1; ...
            diff_sal_dan_infant424344_pmap1; diff_sal_dmn_infant424344_pmap1]
        test_g = [repmat({'1A'}, length(diff_sal_dan_infantpt_pmap1),1); repmat({'1B'}, length(diff_sal_dmn_infantpt_pmap1),1); ... 
            repmat({'2A'}, length(diff_sal_dan_infant373839_pmap1),1); repmat({'2B'}, length(diff_sal_dmn_infant373839_pmap1),1); ... 
            repmat({'3A'}, length(diff_sal_dan_infant40_pmap1),1); repmat({'3B'}, length(diff_sal_dmn_infant40_pmap1),1); ... 
            repmat({'4A'}, length(diff_sal_dan_infant41_pmap1),1); repmat({'4B'}, length(diff_sal_dmn_infant41_pmap1),1); ... 
            repmat({'5A'}, length(diff_sal_dan_infant424344_pmap1),1); repmat({'5B'}, length(diff_sal_dmn_infant424344_pmap1),1)];
   
        figure;
        h = boxplot(test, test_g, 'Symbol', '');
        ylim([0 1.6]) %matched with Child-Adult for comparison
        set(gca,'YTickLabel', get(gca, 'YTickLabel'), 'FontSize',14);
        set(gca,'XTickLabel',{});
        set(gcf,'color','w');


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
            % Odd 
            patch(patchData, get(hBox(iBox), 'YData'), colorsA, 'FaceAlpha', .5);
        else
            % Even 
            patch(patchData, get(hBox(iBox), 'YData'), colorsB, 'FaceAlpha', .5);
        end
    end
        
    
end
    for for_pmap2 = 1

        %combined         
        test = [diff_sal_dan_infantpt_pmap2; diff_sal_dmn_infantpt_pmap2; diff_sal_dan_infant373839_pmap2; diff_sal_dmn_infant373839_pmap2; ...
            diff_sal_dan_infant40_pmap2; diff_sal_dmn_infant40_pmap2; diff_sal_dan_infant41_pmap2; diff_sal_dmn_infant41_pmap2; ...
            diff_sal_dan_infant424344_pmap2; diff_sal_dmn_infant424344_pmap2]
        test_g = [repmat({'1A'}, length(diff_sal_dan_infantpt_pmap2),1); repmat({'1B'}, length(diff_sal_dmn_infantpt_pmap2),1); ... 
            repmat({'2A'}, length(diff_sal_dan_infant373839_pmap2),1); repmat({'2B'}, length(diff_sal_dmn_infant373839_pmap2),1); ... 
            repmat({'3A'}, length(diff_sal_dan_infant40_pmap2),1); repmat({'3B'}, length(diff_sal_dmn_infant40_pmap2),1); ... 
            repmat({'4A'}, length(diff_sal_dan_infant41_pmap2),1); repmat({'4B'}, length(diff_sal_dmn_infant41_pmap2),1); ... 
            repmat({'5A'}, length(diff_sal_dan_infant424344_pmap2),1); repmat({'5B'}, length(diff_sal_dmn_infant424344_pmap2),1)];
   
        figure;
        h = boxplot(test, test_g, 'Symbol', '');
        ylim([0 1.6])
        set(gca,'YTickLabel', get(gca, 'YTickLabel'), 'FontSize',14);
        set(gca,'XTickLabel',{});        
        set(gcf,'color','w');


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
                % Odd Number 
                patch(patchData, get(hBox(iBox), 'YData'), colorsA, 'FaceAlpha', .5);
            else
                % Even Number
                patch(patchData, get(hBox(iBox), 'YData'), colorsB, 'FaceAlpha', .5);
            end
        end
        
    
end

%% -- Figure 4 ------------------------------------------ 

    %% read null spins
    % null spins derived using Brainsmash
    % (https://brainsmash.readthedocs.io/en/latest/)
    null_lh_10000 = readtable('./files/surrogates_lh_10000.csv');
    null_lh_10000 =table2array(null_lh_10000);
    null_lh_10000 =null_lh_10000';

    %% upload gene expression 
    % genetic expression derived using abagen (https://abagen.readthedocs.io/en/stable/)
    for for_allGenes = 1
        fid=fopen([ './files/sch400_expression.csv'  ])
        fmt=['%s', repmat('%f',1,15631)]
        sch400_expression=textscan(fid, fmt, 'collectoutput',true,'headerlines',1,'delimiter',',');
        sch400_expression_name = readtable(['./files/sch400_expression_names.txt'  ]);
        sch400_expression_name_a = table2array(sch400_expression_name);

        sch400_expression = sch400_expression{1, 2};

    end
    
    %% parcellate pmap into schaefer 400

    for i = 1:size(pmap_L_thal_ind_groupCmap_infant_align,1)
        for j = 1:max(schaefer_400_label)
            pmap_sch400_17nt.pmap_R_thal_ind_groupCmap_infant_p1{1,i}(j,:) = nanmean(pmap_R_thal_ind_groupCmap_infant_p1(i,schaefer_400_label == j),2 );    
            pmap_sch400_17nt.pmap_R_thal_ind_groupCmap_infant_p2{1,i}(j,:) = nanmean(pmap_R_thal_ind_groupCmap_infant_p2(i,schaefer_400_label == j),2 );    

            pmap_sch400_17nt.pmap_L_thal_ind_groupCmap_infant_p1{1,i}(j,:) = nanmean(pmap_L_thal_ind_groupCmap_infant_p1(i,schaefer_400_label == j),2 );    
            pmap_sch400_17nt.pmap_L_thal_ind_groupCmap_infant_p2{1,i}(j,:) = nanmean(pmap_L_thal_ind_groupCmap_infant_p2(i,schaefer_400_label == j),2 );    

        end        
    end
    for i = 1:size(pmap_L_thal_ind_groupCmap_infant_align,1) %subject number

            %shaefer400
            sch400_pmap_R_thal_ind_groupCmap_infant_p1(i,:) = pmap_sch400_17nt.pmap_R_thal_ind_groupCmap_infant_p1{1,i};
            sch400_pmap_R_thal_ind_groupCmap_infant_p2(i,:) = pmap_sch400_17nt.pmap_R_thal_ind_groupCmap_infant_p2{1,i};

            sch400_pmap_L_thal_ind_groupCmap_infant_p1(i,:) = pmap_sch400_17nt.pmap_L_thal_ind_groupCmap_infant_p1{1,i};
            sch400_pmap_L_thal_ind_groupCmap_infant_p2(i,:) = pmap_sch400_17nt.pmap_L_thal_ind_groupCmap_infant_p2{1,i};

        end

    %%visualise sch400 pmap 
    sch400_pmap_L_thal_ind_groupCmap_infant_p1_mean = zscore(nanmean(sch400_pmap_L_thal_ind_groupCmap_infant_p1,1))';
    project_detection_community_boview(sch400_pmap_L_thal_ind_groupCmap_infant_p1_mean,'', schaefer_400_label,surf_infant,1);

    sch400_pmap_L_thal_ind_groupCmap_infant_p2_mean = zscore(nanmean(sch400_pmap_L_thal_ind_groupCmap_infant_p2,1))';
    project_detection_community_boview(sch400_pmap_L_thal_ind_groupCmap_infant_p2_mean,'', schaefer_400_label,surf_infant,1);

    %% PLS analysis 

    X = zscore(sch400_expression);
    %X = zscore(sch400_brainGenes_a);
    X_lh = X(1:200,:);
    X_rh = X(201:400,:);

    for pmap1_test=1

        %% original data
        Y1 = zscore(nanmean(sch400_pmap_L_thal_ind_groupCmap_infant_p1,1))';
        Y1_lh=Y1(1:200,:);
        Y1_rh=Y1(201:400,:);

        % pls regression on original data
        [Xl_pmap1,Yl_pmap1,Xs_pmap1,Ys_pmap1,beta_pmap1,pctVar_pmap1,mse_pmap1,stats_pmap1] = plsregress(X_lh,Y1_lh,1); 
        [corr_orig_pmap1 p]=corr(Xs_pmap1,Ys_pmap1)

        % controlling for spatial autocorrelation (permutation) 
        spins = null_lh_10000;          % spatial autocorrelation-preserving permutation assignments
        nspins = 10000;                 % number of permutations ("spins")

        %% spin pmap %%%%%%%%
        for k = 1:nspins    
            Y1_lh_spin = Y1_lh(spins(:,k),:);  % permute pmap1

            [Xl_pmap1_sY,Yl_pmap1_sY,Xs_pmap1_sY,Ys_pmap1_sY,beta_pmap1_sY,pctVar_pmap1_sY,mse_pmap1_sY,stats_pmap1_sY] = plsregress(X_lh,Y1_lh_spin,1); %pls regression on permuted pmap
            null_corr_spinY_pmap1(:,:,k) = corr(Xs_pmap1_sY,Ys_pmap1_sY); % save correlation between LV scores
            null_Xl_pmap1_sY(:,:,k) = Xl_pmap1_sY; % save X loadings
            null_pctVar_pmap1_sY(:,:,k)=pctVar_pmap1_sY; % save percentage of variance explained in LV by each PLS component

        end

        figure; histogram(null_corr_spinY_pmap1);

        %% calculate the probability that the observed correlation coeffcient is above the null distribution
        pval=length(find(abs(null_corr_spinY_pmap1)>abs(corr_orig_pmap1)))/nspins

        %% null distribution of X loadings 
        for ii=1:length(Xl_pmap1)
            pval=length(find(abs(null_Xl_pmap1_sY(ii,:))>abs(Xl_pmap1(ii,1))))/nspins;
            pval_Xl_pmap1_sY(ii) = pval; % save pvalues of each gene

        end
        idx_sig_pval_Xl_pmap1_sY = find(pval_Xl_pmap1_sY<0.05);  % find idx of genes with sig contribution

    end

    for pmap2_test=1

        %% original data
        Y2 = zscore(nanmean(sch400_pmap_L_thal_ind_groupCmap_infant_p2,1))';
        Y2_lh=Y2(1:200,:);
        Y2_rh=Y2(201:400,:);

        %project_detection_community_boview(Y2_lh,'', schaefer_400_label(1:10242,1),surfL,1);

        % pls regression on original data
        [Xl_pmap2,Yl_pmap2,Xs_pmap2,Ys_pmap2,beta_pmap2,pctVar_pmap2,mse_pmap2,stats_pmap2] = plsregress(X_lh,Y2_lh,1);
        [corr_orig_pmap2 p]=corr(Xs_pmap2,Ys_pmap2)

        % visualize latent variables (scores) on cortical surface
    %     project_detection_community_boview(Xs_pmap2,'', schaefer_400_label(1:10242,1),surfL_infant,1);
    %     BoSurfStatColLim([-0.15 0.15]);  
    %     project_detection_community_boview(Ys_pmap2,'', schaefer_400_label(1:10242,1),surfL_infant,1);
    %     BoSurfStatColLim([-20 20]);  

        % controlling for spatial autocorrelation (permutation) 
        spins = null_lh_10000;          % spatial autocorrelation-preserving permutation assignments
        nspins = 10000;                 % number of permutations ("spins")

        %% spin pmap %%%%%%%%
        for k = 1:nspins    
            Y2_lh_spin = Y2_lh(spins(:,k),:);  % permute pmap2

            [Xl_pmap2_sY,Yl_pmap2_sY,Xs_pmap2_sY,Ys_pmap2_sY,beta_pmap2_sY,pctVar_pmap2_sY,mse_pmap2_sY,stats_pmap2_sY] = plsregress(X_lh,Y2_lh_spin,1); %pls regression on permuted pmap
            null_corr_spinY_pmap2(k) = corr(Xs_pmap2_sY,Ys_pmap2_sY); % save correlation between LV scores
            null_Xl_pmap2_sY(:,k) = Xl_pmap2_sY; % save X loadings
            null_pctVar_pmap2_sY(:,k)=pctVar_pmap2_sY; % save percentage of variance explained in LV by each PLS component

        end

        figure; histogram(null_corr_spinY_pmap2);

        %% calculate the probability that the observed correlation coeffcient is above the null distribution
        pval=length(find(abs(null_corr_spinY_pmap2)>abs(corr_orig_pmap2)))/nspins

        %% null distribution of X loadings 
        for ii=1:length(Xl_pmap2)
            pval=length(find(abs(null_Xl_pmap2_sY(ii,:))>abs(Xl_pmap2(ii,1))))/nspins;
            pval_Xl_pmap2_sY(ii) = pval; % save pvalues of each gene

        end
        idx_sig_pval_Xl_pmap2_sY = find(pval_Xl_pmap2_sY<0.05);  % find idx of genes with sig contribution

    end

% pmap2 & significant genes 

pmap2_sig_genes = mean(X(:,idx_sig_pval_Xl_pmap2_sY),2);
project_detection_community_boview(pmap2_sig_genes,'', schaefer_400_label,surf_infant,1);    
BoSurfStatColLim([-0.4 0.4]);  colormap(hcp_colormap(:,1:3));

    %% Gene list for GO analysis 
    gene_list_allGenes = table2array(sch400_expression_name);

    %make gene name + loadings list 
    geneloading_pmap1 = [gene_list_allGenes num2cell(Xl_pmap1)]
    geneloading_pmap2 = [gene_list_allGenes num2cell(Xl_pmap2)]
    
    %% using significant values 
    
    Xl_pmap1_sY_sig=geneloading_pmap1(idx_sig_pval_Xl_pmap1_sY,:)
    Xl_pmap1_sY_sig_sort = sortrows(Xl_pmap1_sY_sig,2,'descend')    

    Xl_pmap2_sY_sig=geneloading_pmap2(idx_sig_pval_Xl_pmap2_sY,:)
    Xl_pmap2_sY_sig_sort = sortrows(Xl_pmap2_sY_sig,2,'descend')    
