
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

    load('./files/yeo7_in_sch200_label.mat', 'yeo7_in_sch200')
    %figure; imagesc(yeo7_in_sch200); colormap(yeo_colormap)    
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
    %surfL_roi_infant_t(:, surfL_roi_infant_t < 0.9)=0; 
    surfL_roi_infant_t(:, surfL_roi_infant_t ~= 0 )=1; 

    temp = gifti([ './files/week-40_hemi-right_space-dhcpSym_dens-10k_desc-medialwall_mask.shape.gii']);
    surfR_roi_infant = temp.cdata'; 
    surfR_roi_infant_t = surfR_roi_infant;   
    %surfR_roi_infant_t(:, surfR_roi_infant_t < 0.9)=0; 
    surfR_roi_infant_t(:, surfR_roi_infant_t ~= 0 )=1; 

    surf_roi_infant = cat(2, surfL_roi_infant_t, surfR_roi_infant_t);
    figure; SurfStatView(surf_roi_infant,surf_infant);      

    %temp = gifti('/local_raid1/01_software/HCPpipelines/global/templates/standard_mesh_atlases/L.atlasroi.10k_fs_LR.shape.gii');
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
        for load_schaefers_200_parcel = 1
            schaefer_200 = ciftiopen([ './files/Schaefer2018_200Parcels_7Networks_order.10k.dlabel.nii'],  WB_COMMAND); 
            schaefer_200_label = schaefer_200.cdata;  

            schaefer_200_full = zeros(20482,1); 
            schaefer_200_full(logical(surf_roi)) = schaefer_200_label;

            figure; SurfStatView(schaefer_200_full,surf_infant);
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
    close all

%% Calculate Euclidean distance 
    %% get centroids of sch200 parcellation 

%get list of coordinates for each region (100x3)

for i=1:100
    sch_label{i} = find(schaefer_200_label==i);    
    sch_label_coord{i} = surf_infant.coord(:, sch_label{i});
    
    sch_label_coord_centroid {i} = mean( sch_label_coord{1,i} , 2 );
    sch_label_coord_centroid_a (i,:) = mean( sch_label_coord{1,i} , 2 );

end
  

%figure; project_detection_community_boview(sch_label_coord_centroid_a(:,1),'', schaefer_200_label(1:10242,1),surfL_infant,1);
figure; project_detection_community_boview(sch_label_coord_centroid_a(:,2),'', schaefer_200_label(1:10242,1),surfL_infant,1);
figure; project_detection_community_boview(sch_label_coord_centroid_a(:,3),'', schaefer_200_label(1:10242,1),surfL_infant,1);
    %% calculate euclidean distance between each node (100 x 100) %%%%%%%%%%%%%%%%% A_DIST MATRIX 
    sch_label_coord_centroid_eucdis = squareform(pdist(sch_label_coord_centroid_a,'euclidean')); %%%%%%% Euclidean distance: A_dist 
    %figure; imagesc(sch_label_coord_centroid_eucdis)
    
%% Calculate binary matrix of FC %%%%%%%%%%%%%%%%% A MATRIX
    %% upload FC matrix for cortex 
    load('./sourceData_Fig5_FCmatrices_dHCP.mat')
    
    %% prepare FC matrix source data by uploading timeseries: do not run and load source data
    % 
    % DTSERIES_DIR = '/local_raid2/03_user/shinwon/01_project/01_thalamusgradient/02_preprocessed_data/06_dHCP_rev/';
    % 
    % TS = cell(1, length(idx_infant));
    % TS_Lthal = cell(1, length(idx_infant));
    % TS_Rthal = cell(1, length(idx_infant));
    % 
    % fc_conte69_infant = cell(1, length(idx_infant));
    % TR = cell(1, length(idx_infant));

    for i = 1 : length(idx_infant)
    %i
%     
%     %upload dtseries
%     cii =  ciftiopen( [ DTSERIES_DIR 'sub-' subid_infant{i} '/ses-' num2str(sesid_infant(i)) '/surface/func_clean/func_mesh-fsLR10k_space-extdhcp40wk_dt_flt_3dtp_volume.dtseries.nii'], WB_COMMAND );
%     
%     ts_temp = cii.cdata';
%     
%     surf_roi_idx = find(surf_roi);
%     
%     ts = zeros(size(ts_temp, 1), size(surf_roi, 2));
%     ts(:, surf_roi_idx) = ts_temp(:, 1:length(surf_roi_idx));
% 
%     %save vertexwise timeseries of cortex
%     TS_cortex{i} = ts(11:end,:)'; %take out first 10 TR
%     TS_all{i} = ts_temp(11:end,:)';
%     
%     %save voxelwise timeseries of thalamus
%     TS_Lthal{i} = cii.cdata(cii.diminfo{1, 1}.models{4, 1}.start:cii.diminfo{1, 1}.models{4, 1}.start+cii.diminfo{1, 1}.models{4, 1}.count - 1 , 11:end); %left thalamus
%     TS_Rthal{i} = cii.cdata(cii.diminfo{1, 1}.models{5, 1}.start:cii.diminfo{1, 1}.models{5, 1}.start+cii.diminfo{1, 1}.models{5, 1}.count - 1 , 11:end); %left thalamus
%     
%     %GSR 
%     if(strcmp(GSR_flag, '_GSR'))
%         rg = mean(ts, 2);
%         beta=(rg'*rg)\rg'*ts; 
%         TS{i}=(ts-rg*beta)';
%         
%         rg_Lthal = mean(TS_Lthal{1,i}', 2);
%         beta_Lthal=(rg_Lthal'*rg_Lthal)\rg_Lthal'*TS_Lthal{1,i}'; 
%         TS_Lthal{i}=(TS_Lthal{1,i}' - rg_Lthal*beta_Lthal)';        
%         
%         rg_Rthal = mean(TS_Rthal{1,i}', 2);
%         beta_Rthal=(rg_Rthal'*rg_Rthal)\rg_Rthal'*TS_Rthal{1,i}'; 
%         TS_Rthal{i}=(TS_Rthal{1,i}' - rg_Rthal*beta_Rthal)';              
%     end
%         
%     %schaefer100    
%     ts_clus1 = zeros(size(TS_cortex{1,i},2), 100); %max(schaefer_200_label)); 
%     for j = 1:100 %max(schaefer_200_label)
%         ts_clus1(:,j) = mean(TS_cortex{i}(schaefer_200_label == j,:),1); 
%     end
%     
%     TS_clus.sch200{i} = ts_clus1;
%     
%     %calculate FC matrix for cortex (schaefer 200 parcellation)
%     r1 = corr(double(ts_clus1(:,1:100)));
%     z1 = 1/2*log((1+r1)./(1-r1));    
%     fc_sch200_infant{i} = int16(z1*10000); %%this is the input matrix   
%     fc_sch200_infant_z{i} = z1; %%this is the input matrix   
% 
%     %calculate FC matrix for thal x cortex (schaefer 100)
%     r2 = corr(TS_Lthal{1,i}', TS_clus.sch200{i});   
%     z2 = 1/2*log((1+r2)./(1-r2));    
%     fc_sch200_lthal_infant{i} = int16(z2*10000);  
%     fc_sch200_lthal_infant_z{i} = z2;  
%     
%     r3 = corr(TS_Rthal{1,i}', TS_clus.sch200{i});  
%     z3 = 1/2*log((1+r3)./(1-r3));    
%     fc_sch200_rthal_infant{i} = int16(z3*10000);       
%     fc_sch200_rthal_infant_z{i} = z3;       

end
    % save('sourceData_Fig5_FCmatrices_dHCP.mat','fc_sch200_infant_z','fc_sch200_infant','fc_sch200_lthal_infant_z','fc_sch200_lthal_infant','fc_sch200_rthal_infant_z','fc_sch200_rthal_infant')

    %% calculate affinity matrix (individual & group mean) 
    
    for i = 1: length(idx_infant)
        i
        corrR_fc_sch200_infant_z(i,:,:) = fc_sch200_infant_z{1,i};        
        corrR_fc_sch200_infant_z_1sub = fc_sch200_infant_z{1,i};
        
        corrR_fc_sch200_infant_z_1sub(corrR_fc_sch200_infant_z_1sub==Inf)= 0; %change diagonal to 0        
        aff_mat_fc_sch200_infant_z{i} = squareform(1-pdist(corrR_fc_sch200_infant_z_1sub,'cosine')) + eye(size(corrR_fc_sch200_infant_z_1sub,1));     
    
    end   
     
    figure; imagesc(aff_mat_fc_sch200_infant_z{1})    
   
    %% threshold connectomes %TO MAKE A (=cell with binary matrices for each subject) 
    % fc_sch200_infant_z -> this needs to be thresholded into binary matrix

    % number of subjects
    nsub = 255;
    % set threshold
    thr  = 0.1;
    % initialise
    binarised_connectomes_sch200_infant = [];
    % loop through subjects
    W_thr = [];
    for sub = 1:nsub;
        W     = squeeze(corrR_fc_sch200_infant_z(sub,:,:));          % take unthresholded connectome
        W_thr = threshold_proportional(W,thr);                              % absolute threshold
        W_bin = weight_conversion(W_thr, 'binarize');
        binarised_connectomes_sch200_infant(sub,:,:) = W_bin;
    end

    %array to cell
    for sub=1:nsub;
        A{sub}=squeeze(binarised_connectomes_sch200_infant(sub,:,:)); %A (cell of binarized connecome for each subject)
    end

    figure; imagesc(A{100}); colorbar
    
%% Calculate cortex-thalamus similarity matrix %%%%%%%%%%%%%%%%% PD3 matrix or A matrix
    %% use code from python to calculate similarity matrix - eta2 

    %for individual participants (100x100)
    %left thalamus x cortex 
    S_set = cell(1,255)
    for h=1:length(idx_infant) 
        h %to keep track of which subject is being processed
        X = fc_sch200_lthal_infant_z{1,h}; 
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
    
    %figure; imagesc(S_set{10}); colorbar

    %% for age groups (100x100) - calculate mean of individuals' similarity matrix
    for h=1:length(idx_infant) 
        S_set_a(h,:,:) = S_set{1,h};

    end

        S_set_infant_mean = squeeze(nanmean(S_set_a(:,:,:),1)); %%PD1 for age groups??         
        %figure; imagesc(S_set_infant_mean); colorbar; caxis([0 1])        
        S_set_infantpt = squeeze(nanmean(S_set_a(idx_infantpt,:,:),1));
        %figure; imagesc(S_set_infantpt); colorbar; caxis([0 1])        
        S_set_infant373839 = squeeze(nanmean(S_set_a(idx_infant373839,:,:),1));
        %figure; imagesc(S_set_infant373839); colorbar; caxis([0 1])
        S_set_infant40 = squeeze(nanmean(S_set_a(idx_infant40,:,:),1));
        %figure; imagesc(S_set_infant40); colorbar; caxis([0 1])        
        S_set_infant41 = squeeze(nanmean(S_set_a(idx_infant41,:,:),1));
        %figure; imagesc(S_set_infant41); colorbar; caxis([0 1])        
        S_set_infant424344 = squeeze(nanmean(S_set_a(idx_infant424344,:,:),1));  
        %figure; imagesc(S_set_infant424344); colorbar; caxis([0 1])        
        S_set_infant_term = squeeze(nanmean(S_set_a([ idx_infant373839; idx_infant40; idx_infant41; idx_infant424344],:,:),1));
        %figure; imagesc(S_set_infant_term); colorbar; caxis([0 1])
            
        S_set_infant = {S_set_infantpt, S_set_infant373839, S_set_infant40, S_set_infant41, S_set_infant424344};
   
    
%% create CGE matrix for each pmap (correlated gene expression) %%%%%%%%%%%%%%%%% PD2 MATRIX

for for_allGenes = 1
    fid=fopen('./files/sch200_expression.csv');
    fmt=['%s', repmat('%f',1,15632)];
    sch200_expression=textscan(fid, fmt, 'collectoutput',true,'headerlines',1,'delimiter',',');
    sch200_expression_name = readtable('./files/sch200_expression_names.txt');
    sch200_expression_name_a = table2array(sch200_expression_name);
    
    sch200_expression = sch200_expression{1, 2}(1:100,:);
    
    %calculate 100 x 100 matrix (correlated gene expression) 
    r = corr(double(sch200_expression'));
    z = 1/2*log((1+r)./(1-r));
    cge_sch200_all = int16(z*10000);
    cge_sch200_all_z = z;    
end
for for_brainSpecificGenes=1
    
    brain_genes = readtable( './files/brain-specific_gene_list.txt' );
    brain_genes_a = table2array(brain_genes);
    
    index_brain_genes = cellfun(@(a) strmatch(a,sch200_expression_name_a,'exact'),brain_genes_a,'uniform',false);
    index_brain_genes=cell2mat(index_brain_genes);

    sch200_expression_brain_genes = sch200_expression(:, index_brain_genes);


    %calculate 100 x 100 matrix (correlated gene expression)
    r = corr(double(sch200_expression_brain_genes'));
    z = 1/2*log((1+r)./(1-r));
    cge_sch200_brain_genes = int16(z*10000);
    cge_sch200_brain_genes_z = z;

end


%% make into struct for analysis 

mdldata.A_dist = sch_label_coord_centroid_eucdis; 
mdldata.adjs = A; 
mdldata.uCGE_all = cge_sch200_all_z; 
mdldata.uCGE_brainGenes = cge_sch200_brain_genes_z; 
mdldata.sub_id = subid_infant;
mdldata.thalcort = S_set_infant_mean; 
mdldata.thalcort_term = S_set_infant_term; 
mdldata.thalcort_preterm = S_set_infantpt; 
mdldata.thalcort_ageGroups = S_set_infant; 

%% save mdldata into .mat for input to generative network modeling (Oldham et al., 2022) 
% code for GNM can be found at: https://github.com/StuartJO/GenerativeNetworkModel
%save('mdldata_spark.mat','mdldata');





