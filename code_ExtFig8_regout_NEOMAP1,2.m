
%% -- Preparation ------------------------------------------
    %% 01) load toolboxes and define paths
    for pathtoolbox = 1

        %download the follwoing tools
        % surfstat
        % brainspace
        % cifti-matlab
        % BCT 
        % HCPPipelines
        % Stuart Oldham's code - github: https://github.com/StuartJO/GenerativeNetworkModel
        
        WB_COMMAND = [ 'wb_command' ];
        addpath(genpath('/local_raid1/01_software/toolboxes/surfstat/'));   %SurfStatView
        addpath(genpath('/local_raid1/01_software/toolboxes/npy-matlab/')); %upload hcp_colormap
        addpath(genpath('/local_raid1/01_software/toolboxes/matlab_util/')); %BoSurfStatView
        addpath(genpath('/local_raid1/01_software/toolboxes/cifti-matlab/')); %used for gifti (uploading surfaces)
        addpath(genpath('./files'));
        addpath(genpath('/local_raid2/03_user/shinwon/03_software/BrainSpace/matlab'));
        
        addpath('/local_raid1/01_software/toolboxes/BCT_20190303/')
        addpath('/local_raid3/03_user/shinwon/generative_network_model_analysis/01_analysis/0_reference/GenerativeNetworkModel_Oldham/GenerativeNetworkModel/data/Networks/')
        addpath '/local_raid3/03_user/shinwon/generative_network_model_analysis/01_analysis/0_reference/GenerativeNetworkModel_Oldham/GenerativeNetworkModel/code/analysis'


    end

    %addpath '/local_raid2/03_user/shinwon/01_project/01_thalamusgradient/05_analysis_cmap_pmap'
    load('./files/yeo7_in_sch200_label.mat', 'yeo7_in_sch200')
    %figure; imagesc(yeo7_in_sch200); colormap(yeo_colormap)    
    %% 02) read phenotypic info
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
    %% 03) Load surfaces & parcellations & colormaps

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
            for load_schaefers_200_parcel = 1
                schaefer_200 = ciftiopen([ './files/Schaefer2018_200Parcels_7Networks_order.10k.dlabel.nii'],  WB_COMMAND); 
                %figure; plot_hemispheres(schaefer_400.cdata, {surf_lh,surf_rh});
                %figure; SurfStatView(schaefer_400.cdata,surf);
                schaefer_200_label = schaefer_200.cdata;  


                schaefer_200_full = zeros(20482,1); 
                schaefer_200_full(logical(surf_roi)) = schaefer_200_label;

                figure; SurfStatView(schaefer_200_full,surf);
                schaefer_200_label = schaefer_200_full;                    


        end      
        end   

    end    
    for for_colormap_load=1
        defaultCmap = [102,103,171;...
        204,185,126;...
        210,147,128;]./255;
    
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


        addpath(genpath('./files/customcolormap'))

        pasteljet=customcolormap_preset('pasteljet')
        rdylbl=customcolormap_preset('red-yellow-blue')
        rdwhbl=customcolormap_preset('red-white-blue')


    end
  
%% Generative network modeling - analysis
% run GNM scripts from https://github.com/StuartJO/GenerativeNetworkModel
 
%% ------------------------------------------------------------------
%% Figure code
    %% Ext Figure 8B: upload data- regout NEOMAP 1,2
    load('sourceData_ExtFig8_regout.mat')

    %using OptimMdl
    for i=1:336
        temp2(:,:,i) = a_reg.OptimMdl{1, 1}.min_maxKS.adjmat{1, i};
        temp3(:,:,i) = a_reg.OptimMdl{1, 2}.min_maxKS.adjmat{1, i};

    end

    %'Spatial','Spatial+ThalCort','ThalCort','Spatial+uCGE','uCGE','Thalcort+uCGE','Matching'
    adjmat_mdl2_mean = mean(temp2, 3);
    adjmat_mdl3_mean = mean(temp3, 3);
    %figure; imagesc(adjmat_mdl2_mean);colorbar; caxis([-1 1]); colormap(flipud(rdwhbl))%colormap(hcp_colormap(:,1:3))
    %figure; imagesc(adjmat_mdl3_mean);colorbar; caxis([-1 1]); colormap(flipud(rdwhbl))%colormap(hcp_colormap(:,1:3))
    %% Ext Figure 8B: upload data- ORIG
    load('sourceData_ExtFig8_orig.mat')

    %using OptimMdl
    for i=1:603 

        temp1(:,:,i) = a.OptimMdl{1, 1}.min_maxKS.adjmat{1, i};
        temp2(:,:,i) = a.OptimMdl{1, 2}.min_maxKS.adjmat{1, i};
        temp3(:,:,i) = a.OptimMdl{1, 3}.min_maxKS.adjmat{1, i};
        temp4(:,:,i) = a.OptimMdl{1, 4}.min_maxKS.adjmat{1, i};
        temp5(:,:,i) = a.OptimMdl{1, 5}.min_maxKS.adjmat{1, i};
        temp6(:,:,i) = a.OptimMdl{1, 6}.min_maxKS.adjmat{1, i};
        temp7(:,:,i) = a.OptimMdl{1, 7}.min_maxKS.adjmat{1, i};

    end

    %'Spatial','Spatial+ThalCort','ThalCort','Spatial+uCGE','uCGE','Thalcort+uCGE','Matching'
    adjmat_mdl2_mean_orig = mean(temp2, 3);
    adjmat_mdl3_mean_orig = mean(temp3, 3);
    
    %% Ext Figure 8B: make a surface gradient - regout NEOMAP 1,2 

    %'Spatial','Spatial+ThalCort','ThalCort','Spatial+uCGE','uCGE','Thalcort+uCGE','Matching'

    %using spatial+ThalCort
    conn_matrix2=adjmat_mdl2_mean

    gm_typ_regout = GradientMaps('kernel', 'na', 'approach', 'dm');
    gm_typ_regout = gm_typ_regout.fit(conn_matrix2,'sparsity', 90 );
    scree_plot(gm_typ_regout.lambda{1});

    project_detection_community_boview(gm_typ_regout.gradients{1}(:,1)*-1,'', schaefer_200_label(1:10242,1),surfL,1);
    BoSurfStatColLim([-0.09 0.09]); colormap((hcp_colormap(:,1:3)));
    %% Ext Figure 8B: make a surface gradient - ORIG

    %using spatial+ThalCort
    conn_matrix2=adjmat_mdl2_mean_orig
    gm_typ_orig = GradientMaps('kernel', 'na', 'approach', 'dm');
    gm_typ_orig = gm_typ_orig.fit(conn_matrix2,'sparsity', 90 );
    scree_plot(gm_typ_orig.lambda{1});

    project_detection_community_boview(gm_typ_orig.gradients{1}(:,1)*-1,'', schaefer_200_label(1:10242,1),surfL,1);
    BoSurfStatColLim([-0.09 0.09]); colormap((hcp_colormap(:,1:3)));

    %% Ext Figure 8B: calculate correlation 

    reg_grad=gm_typ_regout.gradients{1}(:,1); 
    orig_grad=gm_typ_orig.gradients{1}(:,1); 

    [R P] = corr(reg_grad, orig_grad)
    figure; scatter(reg_grad, orig_grad, 30, 'k', 'filled'); lsline;
    xlim([-0.2 0.2]); ylim([-0.2 0.2]); set(gcf,'color','w');
    ax=gca;
    ax.FontSize=15;

%% ------------------------------------------------------------------
%% Extended Figure 8A: Correlation between cortico-cortical map & neocortical projection map
    
    %% upload source data for NEOMAPs
    load('sourceData_Fig2,3_ExtFig3,4_HCPD_NEOMAPs.mat')
    
    pmap_R_thal_ind_groupCmap_hcpd_p1 = pmap_R_thal_ind_groupCmap_hcpd_align(:,:,1);
    pmap_R_thal_ind_groupCmap_hcpd_p2 = pmap_R_thal_ind_groupCmap_hcpd_align(:,:,2);

    pmap_L_thal_ind_groupCmap_hcpd_p1 = pmap_L_thal_ind_groupCmap_hcpd_align(:,:,1);
    pmap_L_thal_ind_groupCmap_hcpd_p2 = pmap_L_thal_ind_groupCmap_hcpd_align(:,:,2);        
    %% parcellate pmap into schaefer 400 

    for i = 1:size(pmap_R_thal_ind_groupCmap_hcpd,1)
        for j = 1:max(schaefer_400_label)
            pmap_sch400_17nt.pmap_L_thal_ind_groupCmap_hcpd_p1{1,i}(j,:) = nanmean(pmap_L_thal_ind_groupCmap_hcpd_p1(i,schaefer_400_label == j),2 );    
            pmap_sch400_17nt.pmap_L_thal_ind_groupCmap_hcpd_p2{1,i}(j,:) = nanmean(pmap_L_thal_ind_groupCmap_hcpd_p2(i,schaefer_400_label == j),2 );    

        end        
    end
    for i=1:603 %subject number
        %shaefer400
        sch400_pmap_L_thal_ind_groupCmap_hcpd_p1(i,:) = pmap_sch400_17nt.pmap_L_thal_ind_groupCmap_hcpd_p1{1,i};
        sch400_pmap_L_thal_ind_groupCmap_hcpd_p2(i,:) = pmap_sch400_17nt.pmap_L_thal_ind_groupCmap_hcpd_p2{1,i};
    end

    sch400_pmap_L_thal_ind_groupCmap_hcpd_p1_mean = zscore(nanmean(sch400_pmap_L_thal_ind_groupCmap_hcpd_p1,1))';
    project_detection_community_boview(sch400_pmap_L_thal_ind_groupCmap_hcpd_p1_mean,'', schaefer_400_label,surf,1);
    colormap(rbgyr20_colormap); BoSurfStatColLim([-2 2]); 

    sch400_pmap_L_thal_ind_groupCmap_hcpd_p2_mean = zscore(nanmean(sch400_pmap_L_thal_ind_groupCmap_hcpd_p2,1))';
    project_detection_community_boview(sch400_pmap_L_thal_ind_groupCmap_hcpd_p2_mean,'', schaefer_400_label,surf,1);
    colormap(rbgyr20_colormap); BoSurfStatColLim([-2 2]); 

    %% upload source data for cortico-cortical maps
    
    load('./files/ExtFig8_fc_sch400_hcpd.mat');

    for for_prepare_source_data = 1 %upload source data, do not run
        
        %do not rerun (upload timeseries) ============================================
%         
%         DTSERIES_DIR = '/local_raid2/03_user/shinwon/01_project/01_thalamusgradient/02_preprocessed_data/04_HCPD/fmriresults01/';
%         TS_Lthal = cell(1, length(idx_hcpd));
%         TS_Rthal = cell(1, length(idx_hcpd));
% 
%         fc_sch400_hcpd = cell(1, length(idx_hcpd));
%         
%         TR = cell(1, length(idx_hcpd));

        for i = 1 : length(idx_hcpd)
            %i
% 
%             %upload dtseries
%             cii = ciftiopen( [ DTSERIES_DIR subid_hcpd{i} '_V1_MR/MNINonLinear/Results/rfMRI_REST/rfMRI_REST_Atlas_MSMAll_hp0_clean_subcortical_2mm_10k.dtseries.nii'], WB_COMMAND );
%             ts_temp = cii.cdata';
%             %TR{i} = cii.diminfo{1,2}.seriesStep;
% 
%             surf_roi_idx = find(surf_roi);
% 
%             ts = zeros(size(ts_temp, 1), size(surf_roi, 2));
%             ts(:, surf_roi_idx) = ts_temp(:, 1:length(surf_roi_idx));
%             ts = ts(11:end,:)'; %take out first 10 TR
%             %TS{i} = ts(11:end,:)'; %take out first 10 TR
% 
%             ts_sc = ts_temp(11:end,cii.diminfo{1, 1}.models{3, 1}.start:end)'; %take out first 10 TR
%             %TS_SC{i} = ts_temp(11:end,cii.diminfo{1, 1}.models{3, 1}.start:end)'; %take out first 10 TR
% 
%             TS_Lthal{i} = cii.cdata(cii.diminfo{1, 1}.models{20, 1}.start:cii.diminfo{1, 1}.models{20, 1}.start+cii.diminfo{1, 1}.models{20, 1}.count - 1 , 11:end); %left thalamus
%             TS_Rthal{i} = cii.cdata(cii.diminfo{1, 1}.models{21, 1}.start:cii.diminfo{1, 1}.models{21, 1}.start+cii.diminfo{1, 1}.models{21, 1}.count - 1 , 11:end); %left thalamus
% 
%             %GSR 
%             if(strcmp(GSR_flag, '_GSR'))
%                 rg = mean(ts, 2);
%                 beta=(rg'*rg)\rg'*ts; 
%                 TS{i}=(ts-rg*beta)';
% 
%                 rg_Lthal = mean(TS_Lthal{1,i}', 2);
%                 beta_Lthal=(rg_Lthal'*rg_Lthal)\rg_Lthal'*TS_Lthal{1,i}'; 
%                 TS_Lthal{i}=(TS_Lthal{1,i}' - rg_Lthal*beta_Lthal)';        
% 
%                 rg_Rthal = mean(TS_Rthal{1,i}', 2);
%                 beta_Rthal=(rg_Rthal'*rg_Rthal)\rg_Rthal'*TS_Rthal{1,i}'; 
%                 TS_Rthal{i}=(TS_Rthal{1,i}' - rg_Rthal*beta_Rthal)';              
%             end


            for for_fc_matrix = 1 

                for for_parcelbyparcel = 1
%                     %schaefer400    
%                     ts_clus1 = zeros(size(ts,2), max(schaefer_400_label)); 
%                     for j = 1: max(schaefer_400_label)
%                         ts_clus1(:,j) = mean(ts(schaefer_400_label == j,:),1); 
%                     end
% 
%                     TS_clus.sch400{i} = ts_clus1;
% 
%                     r2 = corr(double(ts_clus1));
%                     z2 = 1/2*log((1+r2)./(1-r2));    
%                     fc_sch400_hcpd_int{i} = int16(z2*10000);  
%                     fc_sch400_hcpd{i} = z2;      

                end

            end

        end

    end
    
    %save('./files/ExtFig8_fc_sch400_hcpd.mat','fc_sch400_hcpd', 'fc_sch400_hcpd_int')

    %% construct cortico-cortical gradient (brainspace)-Schaefer400
    % First load mean connectivity matrix and Schaefer parcellation
    % conn_matrix = load_group_fc('schaefer',400);
    % labeling = load_parcellation('schaefer',400);

    for i = 1 : length(idx_hcpd)
        %i
        fc_cortex_sch400(i,:,:) = fc_sch400_hcpd_int{1, i}';

    end
    conn_matrix1 = mean(fc_cortex_sch400,1);
    conn_matrix1=squeeze(conn_matrix1);
    conn_matrix1(conn_matrix1==Inf) = 1;

    %figure; imagesc(conn_matrix1)%(1:200,1:200))
    %colorbar
    %caxis([0 0.5])

    %h = plot_hemispheres(schaefer_400_label, {surfL,surfR});
    %colormap(h.handles.figure,lines(401))

    % Construct the gradients

    gm_typ = GradientMaps('kernel', 'cs', 'approach', 'le');
    gm_typ = gm_typ.fit(conn_matrix1,'sparsity', 90 );
    scree_plot(gm_typ.lambda{1});

    per = 1./gm_typ.lambda{1};
    per = per/sum(per)*100
    scree_plot(per);


    project_detection_community_boview(gm_typ.gradients{1}(:,1),'', schaefer_400_label,surf,1);
    BoSurfStatColLim([-0.02 0.02]); colormap(hcp_colormap(:,1:3));


    project_detection_community_boview(gm_typ.gradients{1}(:,2),'', schaefer_400_label,surf,1);
    BoSurfStatColLim([-0.02 0.02]); colormap(hcp_colormap(:,1:3));

    project_detection_community_boview(gm_typ.gradients{1}(:,3),'', schaefer_400_label,surf,1);
    BoSurfStatColLim([-0.02 0.02]); colormap(hcp_colormap(:,1:3));

        norm_grad1=normalize(gm_typ.gradients{1}(:,1));
        norm_grad2=normalize(gm_typ.gradients{1}(:,2));
        norm_grad3=normalize(gm_typ.gradients{1}(:,3));

        project_detection_community_boview(norm_grad1*-1,'', schaefer_400_label,surf,1);
        colormap(jet);BoSurfStatColLim([-1.7 1.7]); 
        project_detection_community_boview(norm_grad2*-1,'', schaefer_400_label,surf,1);
        colormap(jet);BoSurfStatColLim([-1.7 1.7]); 
        project_detection_community_boview(norm_grad3*-1,'', schaefer_400_label,surf,1);
        colormap(jet);BoSurfStatColLim([-1.7 1.7]); 


    h=gradient_in_euclidean(gm_typ.gradients{1}(:,1:2), {surfL,surfR},schaefer_400_label);
    h.axes_scatter.XLim = [-0.02 0.02]
    h.axes_scatter.YLim = [-0.02 0.02]

    plot_hemispheres(gm_typ.gradients{1}(:,1:3), {surfL,surfR}, ...
                 'parcellation', schaefer_400_label, ...
                 'labeltext',{'Gradient 1','Gradient 2','Gradient 3'});

    %% calculate correlation & scatter plots

    %upload surrogates
    null_lh_10000 = readtable('./files/surrogates_lh_10000.csv');
    null_lh_10000 =table2array(null_lh_10000);
    null_lh_10000 =null_lh_10000';

   for CCMAP1_and_PMAPS =1
       
       %% CCMAP1 - NEOMAP 1
        project_detection_community_boview(norm_grad1,'', schaefer_400_label,surf,1);
        colormap(rdwhbl); BoSurfStatColLim([-2 2]); 
        project_detection_community_boview(sch400_pmap_L_thal_ind_groupCmap_hcpd_p1_mean,'', schaefer_400_label,surf,1);
        colormap(rdwhbl); BoSurfStatColLim([-2 2]); 

        %pmap 1: original correlation
        [R_real_pmap1_ccmap1 P_real_pmap1_ccmap1] = corr(sch400_pmap_L_thal_ind_groupCmap_hcpd_p1_mean, norm_grad1)
        figure; scatter(sch400_pmap_L_thal_ind_groupCmap_hcpd_p1_mean, norm_grad1, ...
            11, 'filled', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [0 0 0] ); lsline
        xlabel('PMAP1'); 
        ylabel('CCmap1'); 
        xlim([-3 3]); ylim([-3 3])

        %spin test
        X = sch400_pmap_L_thal_ind_groupCmap_hcpd_p1_mean;
        X_lh = X(1:200,:);
        X_rh = X(201:400,:);

        % original data
        Y1 = norm_grad1; %zscored cortico-cortical gradient 1
        Y1_lh=Y1(1:200,:);
        Y1_rh=Y1(201:400,:);
        
        % controlling for spatial autocorrelation (permutation) 
        spins = null_lh_10000;          % spatial autocorrelation-preserving permutation assignments
        nspins = 10000;                 % number of permutations ("spins")

        % spin pmap  %%%%%%%%
        for k = 1:nspins    
            X_lh_spin = X_lh(spins(:,k),:);  % permute pmap1

            [R_null_pmap1_ccmap1 P_null_pmap1_ccmap1] = corr(X_lh_spin,Y1_lh);
            null_R_pmap1_ccmap1(:,k) = R_null_pmap1_ccmap1; % save correlation between LV scores

        end
        
        project_detection_community_boview(Y1_lh,'', schaefer_400_label(1:10242,1),surfL,1);
        colormap(rdwhbl); BoSurfStatColLim([-2 2]); 
        project_detection_community_boview(X_lh_spin,'', schaefer_400_label(1:10242,1),surfL,1);
        colormap(rdwhbl); BoSurfStatColLim([-2 2]); 
        figure; histogram(null_R_pmap1_ccmap1);

        %calculate the probability that the observed correlation coeffcient is above the null distribution
        pval_spin_pmap1_ccmap1=length(find(abs(null_R_pmap1_ccmap1)>abs(R_real_pmap1_ccmap1)))/nspins
        Rsquared_pmap1_ccmap1 = R_real_pmap1_ccmap1^2
        
       %% CCMAP1 - NEOMAP 2

        %pmap 2
        [R_real_pmap2_ccmap1 P_real_pmap2_ccmap1] = corr(sch400_pmap_L_thal_ind_groupCmap_hcpd_p2_mean, norm_grad1)
        figure; scatter(sch400_pmap_L_thal_ind_groupCmap_hcpd_p2_mean, norm_grad1, ...
            11, 'filled', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [0 0 0] ); lsline %90
        xlabel('PMAP2')
        ylabel('CCmap1')        
        xlim([-3 3]); ylim([-3 3])
        
        %spin test
        X = sch400_pmap_L_thal_ind_groupCmap_hcpd_p2_mean;
        X_lh = X(1:200,:);
        X_rh = X(201:400,:);

        % original data
        Y1 = norm_grad1; %zscored cortico-cortical gradient 1
        Y1_lh=Y1(1:200,:);
        Y1_rh=Y1(201:400,:);
        
        % controlling for spatial autocorrelation (permutation) 
        spins = null_lh_10000;          % spatial autocorrelation-preserving permutation assignments
        nspins = 10000;                 % number of permutations ("spins")

        % spin pmap  %%%%%%%%
        for k = 1:nspins    
            X_lh_spin = X_lh(spins(:,k),:);  % permute pmap1

            [R_null_pmap2_ccmap1 P_null_pmap2_ccmap1] = corr(X_lh_spin,Y1_lh);
            null_R_pmap2_ccmap1(:,k) = R_null_pmap2_ccmap1; % save correlation between LV scores

        end
        
        project_detection_community_boview(Y1_lh,'', schaefer_400_label(1:10242,1),surfL,1);
        colormap(rdwhbl); BoSurfStatColLim([-2 2]); 
        project_detection_community_boview(X_lh_spin,'', schaefer_400_label(1:10242,1),surfL,1);
        colormap(rdwhbl); BoSurfStatColLim([-2 2]); 
        figure; histogram(null_R_pmap1_ccmap1);

        %calculate the probability that the observed correlation coeffcient is above the null distribution
        pval_spin_pmap2_ccmap1=length(find(abs(null_R_pmap2_ccmap1)>abs(R_real_pmap2_ccmap1)))/nspins
        Rsquared_pmap2_ccmap1 = R_real_pmap2_ccmap1^2
        

    end

