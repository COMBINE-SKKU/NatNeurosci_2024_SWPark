
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
        
            
            load('./files/yeo7_in_sch200_label.mat')

            project_detection_community_boview(yeo7_in_sch200,'', schaefer_200_label(1:10242,1)',surfL,1);

            yeo7_in_sch200_vis = find(yeo7_in_sch200 == 1); 
            yeo7_in_sch200_som = find(yeo7_in_sch200 == 2); 
            yeo7_in_sch200_dan = find(yeo7_in_sch200 == 3); 
            yeo7_in_sch200_sal = find(yeo7_in_sch200 == 4); 
            yeo7_in_sch200_lim = find(yeo7_in_sch200 == 5); 
            yeo7_in_sch200_fpn = find(yeo7_in_sch200 == 6); 
            yeo7_in_sch200_dmn = find(yeo7_in_sch200 == 7); 
            
            
            
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

%% calculate Euclidean distance 
    %% get centroids of sch200 parcellation 

    %get list of coordinates for each region (100x3)

    for i=1:100
        sch_label{i} = find(schaefer_200_label==i);    
        sch_label_coord{i} = surf.coord(:, sch_label{i});

        sch_label_coord_centroid {i} = mean( sch_label_coord{1,i} , 2 );
        sch_label_coord_centroid_a (i,:) = mean( sch_label_coord{1,i} , 2 );

    end

    % 
     figure; project_detection_community_boview(sch_label_coord_centroid_a(:,1),'', schaefer_200_label(1:10242,1),surfL,1);
     figure; project_detection_community_boview(sch_label_coord_centroid_a(:,2),'', schaefer_200_label(1:10242,1),surfL,1);
     figure; project_detection_community_boview(sch_label_coord_centroid_a(:,3),'', schaefer_200_label(1:10242,1),surfL,1);
    %% calculate euclidean distance between each node (100 x 100) %%%%%%%%%%%%%%%%% A_DIST MATRIX 
    sch_label_coord_centroid_eucdis = squareform(pdist(sch_label_coord_centroid_a,'euclidean')); %%%%%%% Euclidean distance: A_dist 
    
%% calculate binary matrix of FC %%%%%%%%%%%%%%%%% A MATRIX
    %% upload FC matrix for cortex 
    load('./sourceData_Fig5_FCmatrices_HCPD.mat')

    %% calculate affinity matrix (individual & group mean) 
    
    for i = 1: length(idx_hcpd)
        i
        corrR_fc_sch200_hcpd_z(i,:,:) = fc_sch200_hcpd_z{1,i};        
        corrR_fc_sch200_hcpd_z_1sub = fc_sch200_hcpd_z{1,i};
        
        corrR_fc_sch200_hcpd_z_1sub(corrR_fc_sch200_hcpd_z_1sub==Inf)= 0; %change diagonal to 0        
        aff_mat_fc_sch200_hcpd_z{i} = squareform(1-pdist(corrR_fc_sch200_hcpd_z_1sub,'cosine')) + eye(size(corrR_fc_sch200_hcpd_z_1sub,1));     
    
    end   
      
     %% threshold connectomes %TO MAKE A (=cell with binary matrices for each subject) 
    % fc_sch200_hcpd_z -> this needs to be thresholded into binary matrix

    % number of subjects
    nsub = 603;
    % set threshold
    thr  = 0.1;
    % initialise
    binarised_connectomes_sch200_hcpd = [];
    % loop through subjects
    W_thr = [];
    for sub = 1:nsub;
        W     = squeeze(corrR_fc_sch200_hcpd_z(sub,:,:));          % take unthresholded connectome
        W_thr = threshold_proportional(W,thr);                              % absolute threshold
        W_bin = weight_conversion(W_thr, 'binarize');
        binarised_connectomes_sch200_hcpd(sub,:,:) = W_bin;
    end

    %array to cell
    for sub=1:nsub;
        A{sub}=squeeze(binarised_connectomes_sch200_hcpd(sub,:,:)); %A (cell of binarized connecome for each subject)
    end

    %make a mean FC 
    
    A1 = squeeze(mean(binarised_connectomes_sch200_hcpd,1));
    W_thr1 = threshold_proportional(A1,thr);  % absolute threshold
    W_bin1 = weight_conversion(W_thr1, 'binarize');
    A1b = W_bin1;
    
    figure; imagesc(A1b)
    
    temp = cat(3,A1b,zeros(100,100))
    
    % cell to array
    for n=1:2;
        A1b_c{n}=temp(:,:,n); %A
    end
    
%% Calculate cortex-thalamus similarity matrix %%%%%%%%%%%%%%%%% PD3 matrix or A matrix
    %% use code from python to calculate similarity matrix - eta2 

    %for individual participants (100x100)
    %left thalamus x cortex 
    for h=1:length(idx_hcpd) 
        h %to keep track of which subject is being processed
        %X = fc_sch200_lthal_hcpd{1,i}; 
        X = fc_sch200_lthal_hcpd_z{1,h}; 
        %S = zeros(size(fc_sch200_lthal_hcpd{1,1},2), size(fc_sch200_lthal_hcpd{1,1},2));
        S = zeros(100,100);

        for i=1:100
            for j=i:100
                mi = mean([X(:,i), X(:,j)], 2);
                mm = mean(mi);
                ssw = sum( power(double(X(:,i))-mi, 2) + power( double(X(:,j))-mi, 2) );
                sst = sum( power(double(X(:,i))-mm, 2) + power( double(X(:,j))-mm, 2) );
                S(i,j) = 1-ssw/sst;      
            end
        end
        S = triu(S)+triu(S,1)';
        S_set{h} = S; %%this is the similarity matrix (for thalamus) 
    end

    figure; imagesc(S_set{3}); colorbar; %caxis([0.4 0.6])
    test=load('/local_raid2/03_user/shinwon/01_project/01_thalamusgradient/03_individual_gradient/04_HCPD/01_HCPD_editThal3_allAges_meanFD04_sub2mm_remove10vols_notNorm_10000/sub-HCD0021614/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.dlabel.10.eta2.mat')
    figure; imagesc(test.S); colorbar; caxis([0.4 0.6])
% 
% def eta2(X):
%     
%     S = np.zeros((X.shape[0],X.shape[0]))
%     for i in range(0,X.shape[0]):
%         for j in range(i,X.shape[0]):
%             mi = np.mean([X[i,:],X[j,:]],0) 
%             mm = np.mean(mi)
%             ssw = np.sum(np.square(X[i,:]-mi) + np.square(X[j,:]-mi))
%             sst = np.sum(np.square(X[i,:]-mm) + np.square(X[j,:]-mm))
%             S[i,j] = 1-ssw/sst
%     
%     S += S.T 
%     S -= np.eye(S.shape[0])
%     
%     return S    

    %% for age groups (100x100) - calculate mean of individuals' similarity matrix
    for h=1:length(idx_hcpd) 
        S_set_a(h,:,:) = S_set{1,h};

    end

        S_set_hcpd_mean = squeeze(nanmean(S_set_a(:,:,:),1)); %%PD1 for age groups?? 
        %figure; imagesc(S_set_hcpd_mean); caxis([0 1]); colorbar
        S_set_hcpd_8_10y = squeeze(nanmean(S_set_a(idx_hcpd_8_10y,:,:),1));
        %figure; imagesc(S_set_hcpd_8_10y); caxis([0 1]); colorbar
        S_set_hcpd_10_13y = squeeze(nanmean(S_set_a(idx_hcpd_10_13y,:,:),1));
        %figure; imagesc(S_set_hcpd_10_13y); caxis([0 1]); colorbar
        S_set_hcpd_13_16y = squeeze(nanmean(S_set_a(idx_hcpd_13_16y,:,:),1));
        %figure; imagesc(S_set_hcpd_13_16y); caxis([0 1]); colorbar
        S_set_hcpd_16_20y = squeeze(nanmean(S_set_a(idx_hcpd_16_20y,:,:),1));
        %figure; imagesc(S_set_hcpd_16_20y); caxis([0 1]); colorbar
        S_set_hcpd_20_22y = squeeze(nanmean(S_set_a(idx_hcpd_20_22y,:,:),1));   
        %figure; imagesc(S_set_hcpd_20_22y); caxis([0 1]); colorbar
        
    
        S_set_hcpd = {S_set_hcpd_8_10y, S_set_hcpd_10_13y, S_set_hcpd_13_16y, S_set_hcpd_16_20y, S_set_hcpd_20_22y}

    %% make mean age group thal-cort similarity matrix (100x100)

    idx_hcpd_8_12y = find (age_y_hcpd < 12) ;  
    idx_hcpd_12_18y = find (age_y_hcpd >= 12 & age_y_hcpd < 18) ;  
    idx_hcpd_18_22y = find (age_y_hcpd >= 18) ;  
    
        
        A_set_hcpd_mean = squeeze(nanmean(S_set_a(:,:,:),1)); %%PD1 for age groups?? 
        figure; imagesc(A_set_hcpd_mean); caxis([0 1]); colorbar
        
        A_set_hcpd_8_12y = squeeze(nanmean(S_set_a(idx_hcpd_8_12y,:,:),1));
        figure; imagesc(A_set_hcpd_8_12y); caxis([0 1]); colorbar
        A_set_hcpd_12_18y = squeeze(nanmean(S_set_a(idx_hcpd_12_18y,:,:),1));
        figure; imagesc(A_set_hcpd_12_18y); caxis([0 1]); colorbar
        A_set_hcpd_18_22y = squeeze(nanmean(S_set_a(idx_hcpd_18_22y,:,:),1));
        figure; imagesc(A_set_hcpd_18_22y); caxis([0 1]); colorbar
              
%% Make perturbations to networks
    %% make perturbations to salience networks 
   
     random_perm = randperm(size(A_set_hcpd_8_12y,1));
     random_perm_sal = randperm(size(yeo7_in_sch200_sal,1));
     i=1;
     for i=1:length(yeo7_in_sch200_sal)
         rand_yeo7_sal(i, 1) = yeo7_in_sch200_sal(random_perm_sal(i), 1);

     end

    A_set_hcpd_8_12y_p_sal = A_set_hcpd_8_12y;    
    A_set_hcpd_8_12y_p_sal(yeo7_in_sch200_sal, :) = A_set_hcpd_8_12y(rand_yeo7_sal, random_perm);  
    A_set_hcpd_8_12y_p_sal(:, yeo7_in_sch200_sal) = A_set_hcpd_8_12y(random_perm, rand_yeo7_sal);  
    
    
    figure; imagesc(A_set_hcpd_8_12y);
    figure; imagesc(A_set_hcpd_8_12y_p_sal);
    
     random_perm = randperm(size(A_set_hcpd_12_18y,1));
     random_perm_sal = randperm(size(yeo7_in_sch200_sal,1));
     i=1;
     for i=1:length(yeo7_in_sch200_sal)
         rand_yeo7_sal(i, 1) = yeo7_in_sch200_sal(random_perm_sal(i), 1);

     end

    A_set_hcpd_12_18y_p_sal = A_set_hcpd_12_18y;    
    A_set_hcpd_12_18y_p_sal(yeo7_in_sch200_sal, :) = A_set_hcpd_12_18y(rand_yeo7_sal, random_perm);   
    A_set_hcpd_12_18y_p_sal(:, yeo7_in_sch200_sal) = A_set_hcpd_12_18y(random_perm, rand_yeo7_sal);
    
    
    figure; imagesc(A_set_hcpd_12_18y);
    figure; imagesc(A_set_hcpd_12_18y_p_sal);
    

     random_perm = randperm(size(A_set_hcpd_18_22y,1));
     random_perm_sal = randperm(size(yeo7_in_sch200_sal,1));
     i=1;
     for i=1:length(yeo7_in_sch200_sal)
         rand_yeo7_sal(i, 1) = yeo7_in_sch200_sal(random_perm_sal(i), 1);

     end

    A_set_hcpd_18_22y_p_sal = A_set_hcpd_18_22y;    
    A_set_hcpd_18_22y_p_sal(yeo7_in_sch200_sal, :) = A_set_hcpd_18_22y(rand_yeo7_sal, random_perm);  
    A_set_hcpd_18_22y_p_sal(:, yeo7_in_sch200_sal) = A_set_hcpd_18_22y(random_perm, rand_yeo7_sal);
    
    
    figure; imagesc(A_set_hcpd_18_22y);
    figure; imagesc(A_set_hcpd_18_22y_p_sal) ;       
    %% make perturbations to random-11 nodes (because salience=11 nodes) 

     yeo7_in_sch200_rand11=randi([1, 100], 11, 1);
     random_perm = randperm(size(A_set_hcpd_8_12y,1));
     random_perm_rand11 = randperm(size(yeo7_in_sch200_rand11,1));
     i=1;
     for i=1:length(yeo7_in_sch200_rand11)
         rand_yeo7_rand11(i, 1) = yeo7_in_sch200_rand11(random_perm_rand11(i), 1);

     end

    A_set_hcpd_8_12y_p_rand11 = A_set_hcpd_8_12y;    
    A_set_hcpd_8_12y_p_rand11(yeo7_in_sch200_rand11, :) = A_set_hcpd_8_12y(rand_yeo7_rand11, random_perm);    
    A_set_hcpd_8_12y_p_rand11(:, yeo7_in_sch200_rand11) = A_set_hcpd_8_12y(random_perm, rand_yeo7_rand11);  
    
    
    figure; imagesc(A_set_hcpd_8_12y);
    figure; imagesc(A_set_hcpd_8_12y_p_rand11);
    
    
    
     random_perm = randperm(size(A_set_hcpd_12_18y,1));
     random_perm_rand11 = randperm(size(yeo7_in_sch200_rand11,1));
     i=1;
     for i=1:length(yeo7_in_sch200_rand11)
         rand_yeo7_rand11(i, 1) = yeo7_in_sch200_rand11(random_perm_rand11(i), 1);

     end

    A_set_hcpd_12_18y_p_rand11 = A_set_hcpd_12_18y;    
    A_set_hcpd_12_18y_p_rand11(yeo7_in_sch200_rand11, :) = A_set_hcpd_12_18y(rand_yeo7_rand11, random_perm);    
    A_set_hcpd_12_18y_p_rand11(:, yeo7_in_sch200_rand11) = A_set_hcpd_12_18y(random_perm, rand_yeo7_rand11); 
    
    
    figure; imagesc(A_set_hcpd_12_18y);
    figure; imagesc(A_set_hcpd_12_18y_p_rand11);
    

     random_perm = randperm(size(A_set_hcpd_18_22y,1));
     random_perm_rand11 = randperm(size(yeo7_in_sch200_rand11,1));
     i=1;
     for i=1:length(yeo7_in_sch200_rand11)
         rand_yeo7_rand11(i, 1) = yeo7_in_sch200_rand11(random_perm_rand11(i), 1);

     end

    A_set_hcpd_18_22y_p_rand11 = A_set_hcpd_18_22y;    
    A_set_hcpd_18_22y_p_rand11(yeo7_in_sch200_rand11, :) = A_set_hcpd_18_22y(rand_yeo7_rand11, random_perm);    
    A_set_hcpd_18_22y_p_rand11(:, yeo7_in_sch200_rand11) = A_set_hcpd_18_22y(random_perm, rand_yeo7_rand11);  
    
    
    figure; imagesc(A_set_hcpd_18_22y);
    figure; imagesc(A_set_hcpd_18_22y_p_rand11);   
    %% make perturbations to random-11 nodes in primary sensory (because salience=11 nodes) 

     yeo7_in_sch200_rand11_ps=randi([1, 30], 11, 1);
     random_perm = randperm(size(A_set_hcpd_8_12y,1));
     random_perm_rand11_ps = randperm(size(yeo7_in_sch200_rand11_ps,1));
     i=1;
     for i=1:length(yeo7_in_sch200_rand11_ps)
         rand_yeo7_rand11_ps(i, 1) = yeo7_in_sch200_rand11_ps(random_perm_rand11_ps(i), 1);

     end

    A_set_hcpd_8_12y_p_rand11_ps = A_set_hcpd_8_12y;    
    A_set_hcpd_8_12y_p_rand11_ps(yeo7_in_sch200_rand11_ps, :) = A_set_hcpd_8_12y(rand_yeo7_rand11_ps, random_perm);  
    A_set_hcpd_8_12y_p_rand11_ps(:, yeo7_in_sch200_rand11_ps) = A_set_hcpd_8_12y(random_perm, rand_yeo7_rand11_ps);
    
    
    figure; imagesc(A_set_hcpd_8_12y);
    figure; imagesc(A_set_hcpd_8_12y_p_rand11_ps);
    
    
    
     random_perm = randperm(size(A_set_hcpd_12_18y,1));
     random_perm_rand11_ps = randperm(size(yeo7_in_sch200_rand11_ps,1));
     i=1;
     for i=1:length(yeo7_in_sch200_rand11_ps)
         rand_yeo7_rand11_ps(i, 1) = yeo7_in_sch200_rand11_ps(random_perm_rand11_ps(i), 1);

     end

    A_set_hcpd_12_18y_p_rand11_ps = A_set_hcpd_12_18y;    
    A_set_hcpd_12_18y_p_rand11_ps(yeo7_in_sch200_rand11_ps, :) = A_set_hcpd_12_18y(rand_yeo7_rand11_ps, random_perm);    
    A_set_hcpd_12_18y_p_rand11_ps(:, yeo7_in_sch200_rand11_ps) = A_set_hcpd_12_18y(random_perm, rand_yeo7_rand11_ps);  
    
    
    figure; imagesc(A_set_hcpd_12_18y);
    figure; imagesc(A_set_hcpd_12_18y_p_rand11_ps);
    

     random_perm = randperm(size(A_set_hcpd_18_22y,1));
     random_perm_rand11_ps = randperm(size(yeo7_in_sch200_rand11_ps,1));
     i=1;
     for i=1:length(yeo7_in_sch200_rand11_ps)
         rand_yeo7_rand11_ps(i, 1) = yeo7_in_sch200_rand11_ps(random_perm_rand11_ps(i), 1);

     end

    A_set_hcpd_18_22y_p_rand11_ps = A_set_hcpd_18_22y;    
    A_set_hcpd_18_22y_p_rand11_ps(yeo7_in_sch200_rand11_ps, :) = A_set_hcpd_18_22y(rand_yeo7_rand11_ps, random_perm);    
    A_set_hcpd_18_22y_p_rand11_ps(:, yeo7_in_sch200_rand11_ps) = A_set_hcpd_18_22y(random_perm, rand_yeo7_rand11_ps);  
    
    
    figure; imagesc(A_set_hcpd_18_22y);
    figure; imagesc(A_set_hcpd_18_22y_p_rand11_ps);   

    %% group level perturbed network sets: salience network
    
    nnet=3;
    set_static= A_set_hcpd_mean;
    set_noPerturb=cat(3, A_set_hcpd_8_12y, A_set_hcpd_12_18y, A_set_hcpd_18_22y); %perturb no groups
    set_1Perturb_sal=cat(3, A_set_hcpd_8_12y_p_sal, A_set_hcpd_12_18y, A_set_hcpd_18_22y); %perturb only 8-12 yrs 
    set_2Perturb_sal=cat(3, A_set_hcpd_8_12y, A_set_hcpd_12_18y_p_sal, A_set_hcpd_18_22y); %perturb only 12-18 yrs 
    set_3Perturb_sal=cat(3, A_set_hcpd_8_12y, A_set_hcpd_12_18y, A_set_hcpd_18_22y_p_sal); %perturb only 18-22 yrs 
    set_4Perturb_sal=cat(3, A_set_hcpd_8_12y_p_sal, A_set_hcpd_12_18y_p_sal, A_set_hcpd_18_22y_p_sal); %perturb all age groups 

    % cell to array
    for net=1:nnet

        set_noPerturb_c{net}=set_noPerturb(:,:,net); 
        set_1Perturb_sal_c{net}=set_1Perturb_sal(:,:,net); 
        set_2Perturb_sal_c{net}=set_2Perturb_sal(:,:,net); 
        set_3Perturb_sal_c{net}=set_3Perturb_sal(:,:,net); 
        set_4Perturb_sal_c{net}=set_4Perturb_sal(:,:,net); 
    
    end
    
    figure;imagesc(set_noPerturb(:,:,1));colorbar;caxis([0.4 1]);
    figure;imagesc(set_noPerturb(:,:,2));colorbar; caxis([0.4 1]);
    figure;imagesc(set_noPerturb(:,:,3));colorbar; caxis([0.4 1]);
    
    figure;imagesc(set_1Perturb_sal(:,:,1));colorbar; caxis([0.4 1]);
    figure;imagesc(set_1Perturb_sal(:,:,2));colorbar; caxis([0.4 1]);
    figure;imagesc(set_1Perturb_sal(:,:,3));colorbar; caxis([0.4 1]);
    
    figure;imagesc(set_2Perturb_sal(:,:,1));colorbar; caxis([0.4 1]);
    figure;imagesc(set_2Perturb_sal(:,:,2));colorbar; caxis([0.4 1]);
    figure;imagesc(set_2Perturb_sal(:,:,3));colorbar; caxis([0.4 1]);
    
    figure;imagesc(set_3Perturb_sal(:,:,1));colorbar; caxis([0.4 1]);
    figure;imagesc(set_3Perturb_sal(:,:,2));colorbar; caxis([0.4 1]);
    figure;imagesc(set_3Perturb_sal(:,:,3));colorbar; caxis([0.4 1]);
    
    figure;imagesc(set_4Perturb_sal(:,:,1));colorbar; caxis([0.4 1]);
    figure;imagesc(set_4Perturb_sal(:,:,2));colorbar; caxis([0.4 1]);
    figure;imagesc(set_4Perturb_sal(:,:,3));colorbar; caxis([0.4 1]);
    
    figure;imagesc(set_static);      
    %% group level perturbed network sets: random-11
    
    nnet=3;
    set_static= A_set_hcpd_mean;
    set_noPerturb=cat(3, A_set_hcpd_8_12y, A_set_hcpd_12_18y, A_set_hcpd_18_22y);
    set_1Perturb_rand11=cat(3, A_set_hcpd_8_12y_p_rand11, A_set_hcpd_12_18y, A_set_hcpd_18_22y);
    set_2Perturb_rand11=cat(3, A_set_hcpd_8_12y, A_set_hcpd_12_18y_p_rand11, A_set_hcpd_18_22y);
    set_3Perturb_rand11=cat(3, A_set_hcpd_8_12y, A_set_hcpd_12_18y, A_set_hcpd_18_22y_p_rand11);
    set_4Perturb_rand11=cat(3, A_set_hcpd_8_12y_p_rand11, A_set_hcpd_12_18y_p_rand11, A_set_hcpd_18_22y_p_rand11);

    % cell to array
    for net=1:nnet

        set_noPerturb_c{net}=set_noPerturb(:,:,net); 
        set_1Perturb_rand11_c{net}=set_1Perturb_rand11(:,:,net); 
        set_2Perturb_rand11_c{net}=set_2Perturb_rand11(:,:,net); 
        set_3Perturb_rand11_c{net}=set_3Perturb_rand11(:,:,net); 
        set_4Perturb_rand11_c{net}=set_4Perturb_rand11(:,:,net); 
    
    end
    
    figure;imagesc(set_noPerturb(:,:,1));colorbar;caxis([0.4 1]);
    figure;imagesc(set_noPerturb(:,:,2));colorbar; caxis([0.4 1]);
    figure;imagesc(set_noPerturb(:,:,3));colorbar; caxis([0.4 1]);
    
    figure;imagesc(set_1Perturb_rand11(:,:,1));colorbar; caxis([0.4 1]);
    figure;imagesc(set_1Perturb_rand11(:,:,2));colorbar; caxis([0.4 1]);
    figure;imagesc(set_1Perturb_rand11(:,:,3));colorbar; caxis([0.4 1]);
    
    figure;imagesc(set_2Perturb_rand11(:,:,1));colorbar; caxis([0.4 1]);
    figure;imagesc(set_2Perturb_rand11(:,:,2));colorbar; caxis([0.4 1]);
    figure;imagesc(set_2Perturb_rand11(:,:,3));colorbar; caxis([0.4 1]);
    
    figure;imagesc(set_3Perturb_rand11(:,:,1));colorbar; caxis([0.4 1]);
    figure;imagesc(set_3Perturb_rand11(:,:,2));colorbar; caxis([0.4 1]);
    figure;imagesc(set_3Perturb_rand11(:,:,3));colorbar; caxis([0.4 1]);
    
    figure;imagesc(set_4Perturb_rand11(:,:,1));colorbar; caxis([0.4 1]);
    figure;imagesc(set_4Perturb_rand11(:,:,2));colorbar; caxis([0.4 1]);
    figure;imagesc(set_4Perturb_rand11(:,:,3));colorbar; caxis([0.4 1]);
    
    figure;imagesc(set_static);
    %% group level perturbed network sets: random-11 in primary sensory
    
    nnet=3;
    set_static= A_set_hcpd_mean;
    set_noPerturb=cat(3, A_set_hcpd_8_12y, A_set_hcpd_12_18y, A_set_hcpd_18_22y);
    set_1Perturb_rand11_ps=cat(3, A_set_hcpd_8_12y_p_rand11_ps, A_set_hcpd_12_18y, A_set_hcpd_18_22y);
    set_2Perturb_rand11_ps=cat(3, A_set_hcpd_8_12y, A_set_hcpd_12_18y_p_rand11_ps, A_set_hcpd_18_22y);
    set_3Perturb_rand11_ps=cat(3, A_set_hcpd_8_12y, A_set_hcpd_12_18y, A_set_hcpd_18_22y_p_rand11_ps);
    set_4Perturb_rand11_ps=cat(3, A_set_hcpd_8_12y_p_rand11_ps, A_set_hcpd_12_18y_p_rand11_ps, A_set_hcpd_18_22y_p_rand11_ps);

    % cell to array
    for net=1:nnet

        set_noPerturb_c{net}=set_noPerturb(:,:,net); 
        set_1Perturb_rand11_ps_c{net}=set_1Perturb_rand11_ps(:,:,net); 
        set_2Perturb_rand11_ps_c{net}=set_2Perturb_rand11_ps(:,:,net); 
        set_3Perturb_rand11_ps_c{net}=set_3Perturb_rand11_ps(:,:,net); 
        set_4Perturb_rand11_ps_c{net}=set_4Perturb_rand11_ps(:,:,net); 
    
    end
    
    figure;imagesc(set_noPerturb(:,:,1));colorbar;caxis([0.4 1]);
    figure;imagesc(set_noPerturb(:,:,2));colorbar; caxis([0.4 1]);
    figure;imagesc(set_noPerturb(:,:,3));colorbar; caxis([0.4 1]);
    
    figure;imagesc(set_1Perturb_rand11_ps(:,:,1));colorbar; caxis([0.4 1]);
    figure;imagesc(set_1Perturb_rand11_ps(:,:,2));colorbar; caxis([0.4 1]);
    figure;imagesc(set_1Perturb_rand11_ps(:,:,3));colorbar; caxis([0.4 1]);
    
    figure;imagesc(set_2Perturb_rand11_ps(:,:,1));colorbar; caxis([0.4 1]);
    figure;imagesc(set_2Perturb_rand11_ps(:,:,2));colorbar; caxis([0.4 1]);
    figure;imagesc(set_2Perturb_rand11_ps(:,:,3));colorbar; caxis([0.4 1]);
    
    figure;imagesc(set_3Perturb_rand11_ps(:,:,1));colorbar; caxis([0.4 1]);
    figure;imagesc(set_3Perturb_rand11_ps(:,:,2));colorbar; caxis([0.4 1]);
    figure;imagesc(set_3Perturb_rand11_ps(:,:,3));colorbar; caxis([0.4 1]);
    
    figure;imagesc(set_4Perturb_rand11_ps(:,:,1));colorbar; caxis([0.4 1]);
    figure;imagesc(set_4Perturb_rand11_ps(:,:,2));colorbar; caxis([0.4 1]);
    figure;imagesc(set_4Perturb_rand11_ps(:,:,3));colorbar; caxis([0.4 1]);
    
    figure;imagesc(set_static);        
    
%% save as struct for analysis - salience perturbation

mdldata.A_dist = sch_label_coord_centroid_eucdis; 
mdldata.adjs = A1b_c; 
mdldata.sub_id = subid_hcpd;
mdldata.age = age_y_hcpd;
mdldata.thalcort = S_set_hcpd_mean; 
mdldata.thalcort_ageGroups_p0 = set_noPerturb_c; 
mdldata.thalcort_ageGroups_p1 = set_1Perturb_sal_c; 
mdldata.thalcort_ageGroups_p2 = set_2Perturb_sal_c; 
mdldata.thalcort_ageGroups_p3 = set_3Perturb_sal_c; 
mdldata.thalcort_ageGroups_p4 = set_4Perturb_sal_c;    

%% save as struct for analysis - random-11 nodes  perturbation

mdldata.A_dist = sch_label_coord_centroid_eucdis; 
mdldata.adjs = A1b_c; 
mdldata.sub_id = subid_hcpd;
mdldata.age = age_y_hcpd;
mdldata.thalcort = S_set_hcpd_mean; 
mdldata.thalcort_ageGroups_p0 = set_noPerturb_c; 
mdldata.thalcort_ageGroups_p1 = set_1Perturb_rand11_c; 
mdldata.thalcort_ageGroups_p2 = set_2Perturb_rand11_c; 
mdldata.thalcort_ageGroups_p3 = set_3Perturb_rand11_c; 
mdldata.thalcort_ageGroups_p4 = set_4Perturb_rand11_c; 

%% save as struct for analysis - random-11 nodes in primary sensory perturbation

mdldata.A_dist = sch_label_coord_centroid_eucdis; 
mdldata.adjs = A1b_c; 
mdldata.sub_id = subid_hcpd;
mdldata.age = age_y_hcpd;
mdldata.thalcort = S_set_hcpd_mean; 
mdldata.thalcort_ageGroups_p0 = set_noPerturb_c; 
mdldata.thalcort_ageGroups_p1 = set_1Perturb_rand11_ps_c; 
mdldata.thalcort_ageGroups_p2 = set_2Perturb_rand11_ps_c; 
mdldata.thalcort_ageGroups_p3 = set_3Perturb_rand11_ps_c; 
mdldata.thalcort_ageGroups_p4 = set_4Perturb_rand11_ps_c; 

