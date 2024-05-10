
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

    %no_lim = [2,3,4,5,7,8]
    %yeo_colormap_r = yeo_colormap(no_lim',: )

%% load HCPD - CMAPs
load ('sourceData_ExtFig5_HCPD_CMAPs.mat')

for for_prepare_source_data = 1 %do not run
    %% Group cmap
       for for_grouptemplate_HCPD_10k =1
% 
%             cii = ciftiopen([ './files/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.cmap.dscalar.nii' ], WB_COMMAND);   
%             cii_L_thal = cii.cdata(cii.diminfo{1, 1}.models{20, 1}.start:(cii.diminfo{1,1}.models{20,1}.start-1) + cii.diminfo{1,1}.models{20,1}.count,:);
%             groupGradient_L_hcpd(:, :) = cii_L_thal;
% 
%             cii = ciftiopen([ './files/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.49.cmap.dscalar.nii' ], WB_COMMAND);   
%             cii_R_thal = cii.cdata(cii.diminfo{1,1}.models{21,1}.start:(cii.diminfo{1,1}.models{21,1}.start-1) + cii.diminfo{1,1}.models{21,1}.count,:);
%             groupGradient_R_hcpd(:, :) = cii_R_thal;

        end

    %% read (unaligned) individual gradients (+ do alignment)
    for for_gradient_directory = 1

%         PATH               = '/local_raid2/03_user/shinwon/01_project/01_thalamusgradient/';
%         IND_GRAD_HCPD_ALLAGES         = [ PATH '03_individual_gradient/04_HCPD/' ];

    end
    for for_read_unaligned_ind_gradients =1

        %% HCPD individual gradients (in atlasroi) aligned to HCPD-Group gradient (5-22y)

        %10k - thalamus gradients
        parfor i = 1:length(idx_hcpd)
%             i 
%             cii_L = ciftiopen([ IND_GRAD_HCPD_ALLAGES '01_HCPD_editThal3_allAges_meanFD04_sub2mm_remove10vols_notNorm/sub-' subid_hcpd{i} '/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.cmap.dscalar.nii' ], WB_COMMAND);   
%             cii_L_thal = cii_L.cdata(cii_L.diminfo{1,1}.models{20,1}.start:(cii_L.diminfo{1,1}.models{20,1}.start-1) + cii_L.diminfo{1,1}.models{20,1}.count,:);
%             grad_L_thal_ind_hcpd_noalign(i, :, :) = cii_L_thal;    
% 
%             cii_R = ciftiopen([ IND_GRAD_HCPD_ALLAGES '01_HCPD_editThal3_allAges_meanFD04_sub2mm_remove10vols_notNorm/sub-' subid_hcpd{i} '/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.49.cmap.dscalar.nii' ], WB_COMMAND);   
%             cii_R_thal = cii_R.cdata(cii_R.diminfo{1,1}.models{21,1}.start:(cii_R.diminfo{1,1}.models{21,1}.start-1) + cii_R.diminfo{1,1}.models{21,1}.count,:);
%             grad_R_thal_ind_hcpd_noalign(i, :, :) = cii_R_thal;    
% 
%             %procrustes alignment    
%             [d_L,Z_L,transform_L] = procrustes(groupGradient_L_hcpd,cii_L_thal)%, 'scaling', false, 'reflection', false)
%             [d_R,Z_R,transform_R] = procrustes(groupGradient_R_hcpd,cii_R_thal)%, 'scaling', false, 'reflection', false)
% 
%             cii_L.cdata=zeros(50592,10);
%             cii_L.cdata(cii_L.diminfo{1,1}.models{20,1}.start:(cii_L.diminfo{1,1}.models{20,1}.start-1) + cii_L.diminfo{1,1}.models{20,1}.count,:)=Z_L
%             cii_aligned_L=cii_L  
% 
%             cii_R.cdata=zeros(50592,10);
%             cii_R.cdata(cii_R.diminfo{1,1}.models{21,1}.start:(cii_R.diminfo{1,1}.models{21,1}.start-1) + cii_R.diminfo{1,1}.models{21,1}.count,:)=Z_R
%             cii_aligned_R=cii_R 
% 
%             %ciftisave(cii_aligned_L, [ PATH '05_analysis_cmap_pmap/' results_dir subid_hcpd{i} '_L_thal.cmap.dscalar.nii' ], WB_COMMAND);  
%             %ciftisave(cii_aligned_R, [ PATH '05_analysis_cmap_pmap/' results_dir subid_hcpd{i} '_R_thal.cmap.dscalar.nii' ], WB_COMMAND);  
% 
%             grad_L_thal_ind_hcpd_align(i, :, :) = Z_L;
%             grad_R_thal_ind_hcpd_align(i, :, :) = Z_R;

        end

    end
end

%% load HCPD - NEOMAPs
load('sourceData_Fig2,3_ExtFig3,4_HCPD_NEOMAPs.mat')

pmap_R_thal_ind_groupCmap_hcpd_p1 = pmap_R_thal_ind_groupCmap_hcpd_align(:,:,1);
pmap_R_thal_ind_groupCmap_hcpd_p2 = pmap_R_thal_ind_groupCmap_hcpd_align(:,:,2);

pmap_L_thal_ind_groupCmap_hcpd_p1 = pmap_L_thal_ind_groupCmap_hcpd_align(:,:,1);
pmap_L_thal_ind_groupCmap_hcpd_p2 = pmap_L_thal_ind_groupCmap_hcpd_align(:,:,2);

%% test effects of age: CMAP ===========================
    %% Vertexwise univariate analysis
    % set model 
    %%SurfStat model - set model & contrasts
    % 1. SurfStatLinMod (beta estimation)
    % 2. SurfStatT (contrast)
    % 3. SurfStatP (random field theory) or SurfStatQ (false discovery rate)

    % gradient_R = beta0 + beta1*age + beta2*group
    %GROUP = infant(group); %%female, male
    %M = 1 + AGE + GROUP

    for for_set_model = 1

        for for_set_model=1        

            AGE_HCPD = term(age_w_hcpd);
            SEX_HCPD = term(sex_hcpd);
            MEANFD_HCPD = term(meanFD_hcpd);


        end
        for for_set_model_age_correlations=1

            M_AGE_HCPD = 1 + AGE_HCPD; 
            M_AGE_HCPD_COV_SEXMEANFD = 1 + AGE_HCPD + SEX_HCPD + MEANFD_HCPD; 

        end 



    end
    %% test age effects

    for for_upload_groupGradients_in_struct_format = 1
        groupGradient_L_hcpd_struct = ciftiopen([ './files/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.cmap.dscalar.nii' ], WB_COMMAND);        
        groupGradient_R_hcpd_struct = ciftiopen([ './files/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.49.cmap.dscalar.nii' ], WB_COMMAND);   
        
    end

    for for_test_model_posAge_LThal_hcpd = 1

        sig_cii_1=groupGradient_L_hcpd_struct;
        sig_cii_1.cdata=zeros(length(groupGradient_L_hcpd_struct.cdata),10);
        sig_cii_2=groupGradient_L_hcpd_struct;
        sig_cii_2.cdata=zeros(length(groupGradient_L_hcpd_struct.cdata),10);        
        group_temp=groupGradient_L_hcpd_struct;
        group_temp.cdata(groupGradient_L_hcpd_struct.diminfo{1,1}.models{20,1}.start:(groupGradient_L_hcpd_struct.diminfo{1,1}.models{20,1}.start-1) + groupGradient_L_hcpd_struct.diminfo{1,1}.models{20,1}.count,:)=groupGradient_L_hcpd;
        %ciftisave(group_temp, [PATH '05_analysis_cmap_pmap/' results_dir '/group_gradient_L_thal_hcpd.cmap.dscalar.nii' ], WB_COMMAND);
        
        for comp = 1:2
            slm_posAge_LThal_hcpd = SurfStatLinMod(grad_L_thal_ind_hcpd_align(:,:,comp), M_AGE_HCPD_COV_SEXMEANFD); 
            slm_posAge_LThal_hcpd = SurfStatT(slm_posAge_LThal_hcpd, age_y(idx_hcpd));%behavioral correlation    
            figure; hist(slm_posAge_LThal_hcpd.t);


            %regress out covariates%%%%
            M_regout_temp = 1 + SEX_HCPD + MEANFD_HCPD;
            slm_regout_temp = SurfStatLinMod(grad_L_thal_ind_hcpd_align(:,:,comp), M_regout_temp); 
            grad_L_thal_ind_hcpd_align_regout = grad_L_thal_ind_hcpd_align(:,:,comp) - slm_regout_temp.X*slm_regout_temp.coef;
           %slm_posAge_LThal_hcpd = SurfStatLinMod(grad_L_thal_ind_hcpd_align_regout, M_AGE_hcpd); 

            %FDR-correction
            qval= SurfStatQ(slm_posAge_LThal_hcpd); %
            sig=sum(qval.Q<0.05);
            sig_idx_thal=find(qval.Q<0.05);
            sig_idx_L_thal_hcpd_pos{comp}=[sig_idx_thal];

                    
                    
            if sig~=0 
                grad_values_sigRegion = double(grad_L_thal_ind_hcpd_align_regout(:,sig_idx_thal));     
                figure; SurfStatPlot(age_y_hcpd, mean(grad_values_sigRegion,2)); %group effect
                %ylim([-60 60]); xlim([7 23]);
                [R P] = corr(age_w_hcpd, mean(grad_values_sigRegion,2)) %group effect
                %saveas(gcf, [PATH '/05_analysis_cmap_pmap/' results_dir '/scatter_posAge_L_thal_hcpd_comp' num2str(comp) '.png']) 


                tempq=zeros(1, size(qval.Q,2));             
                tempq(1, qval.Q<0.05)=1;
                
                sig_cii_1.cdata(groupGradient_L_hcpd_struct.diminfo{1,1}.models{20,1}.start:(groupGradient_L_hcpd_struct.diminfo{1,1}.models{20,1}.start-1) + groupGradient_L_hcpd_struct.diminfo{1,1}.models{20,1}.count,comp)=tempq';
                %sig_cii_1.cdata(sig_cii_1.cdata>0.01)=-999;
                sig_cii_qval=sig_cii_1;
                
                sig_idx_wb=find(0<sig_cii_qval.cdata(:,comp) & sig_cii_qval.cdata(:,comp)<0.01);
                sig_idx_wb_L_thal_pos{comp}=[sig_idx_wb];

                sig_cii_2.cdata(groupGradient_L_hcpd_struct.diminfo{1,1}.models{20,1}.start:(groupGradient_L_hcpd_struct.diminfo{1,1}.models{20,1}.start-1) + groupGradient_L_hcpd_struct.diminfo{1,1}.models{20,1}.count,comp)=slm_posAge_LThal_hcpd.t';
                sig_cii_t=sig_cii_2 ;      

                %ciftisave(sig_cii_qval, [PATH '05_analysis_cmap_pmap/' results_dir '/posAge_L_thal_hcpd_qval.cmap.dscalar.nii' ], WB_COMMAND);
                %ciftisave(sig_cii_t, [ PATH '05_analysis_cmap_pmap/' results_dir '/posAge_L_thal_hcpd_tval.cmap.dscalar.nii' ], WB_COMMAND);

            elseif sig == 0 

                %uncorrected p
                p=1-tcdf(slm_posAge_LThal_hcpd.t, slm_posAge_LThal_hcpd.df);
                p = min([ p; 1-p ], [], 1);
                figure; histogram(p); 
                unique(p);
            end

            
        end
        %close all
        clear sig_cii_1 sig_cii_2 group_temp

    end    
    for for_test_model_negAge_LThal_hcpd = 1

        sig_cii_1=groupGradient_L_hcpd_struct;
        sig_cii_1.cdata=zeros(length(groupGradient_L_hcpd_struct.cdata),10);
        sig_cii_2=groupGradient_L_hcpd_struct;
        sig_cii_2.cdata=zeros(length(groupGradient_L_hcpd_struct.cdata),10);        
        group_temp=groupGradient_L_hcpd_struct;
        group_temp.cdata(groupGradient_L_hcpd_struct.diminfo{1,1}.models{20,1}.start:(groupGradient_L_hcpd_struct.diminfo{1,1}.models{20,1}.start-1) + groupGradient_L_hcpd_struct.diminfo{1,1}.models{20,1}.count,:)=groupGradient_L_hcpd;
        %ciftisave(group_temp, [PATH '05_analysis_cmap_pmap/' results_dir '/group_gradient_L_thal_hcpd.cmap.dscalar.nii' ], WB_COMMAND);
                
        for comp = 1:2
            slm_negAge_LThal_hcpd = SurfStatLinMod(grad_L_thal_ind_hcpd_align(:,:,comp), M_AGE_HCPD_COV_SEXMEANFD); 
            slm_negAge_LThal_hcpd = SurfStatT(slm_negAge_LThal_hcpd, -age_y(idx_hcpd));%behavioral correlation    
            %figure; hist(slm_negAge_LThal_hcpd.t);


            %regress out covariates%%%%
            M_regout_temp = 1 + SEX_HCPD + MEANFD_HCPD;
            slm_regout_temp = SurfStatLinMod(grad_L_thal_ind_hcpd_align(:,:,comp), M_regout_temp); 
            grad_L_thal_ind_hcpd_align_regout = grad_L_thal_ind_hcpd_align(:,:,comp) - slm_regout_temp.X*slm_regout_temp.coef;
           %slm_negAge_LThal_hcpd = SurfStatLinMod(grad_L_thal_ind_hcpd_align_regout, M_AGE_hcpd); 

            %FDR-correction
            qval= SurfStatQ(slm_negAge_LThal_hcpd); %
            sig=sum(qval.Q<0.05);
            sig_idx_thal=find(qval.Q<0.05);
            sig_idx_L_thal_hcpd_neg{comp}=[sig_idx_thal];

            if sig~=0 
                grad_values_sigRegion = double(grad_L_thal_ind_hcpd_align_regout(:,sig_idx_thal));     
                figure; SurfStatPlot(age_y_hcpd, mean(grad_values_sigRegion,2)); %group effect
                %ylim([-60 60]); xlim([7 23]);
                [R P] = corr(age_w_hcpd, mean(grad_values_sigRegion,2)) %group effect
               % saveas(gcf, [PATH '/05_analysis_cmap_pmap/' results_dir '/scatter_negAge_L_thal_hcpd_comp' num2str(comp) '.png'])         

                tempq=zeros(1, size(qval.Q,2));             
                tempq(1, qval.Q<0.05)=1;
                
                sig_cii_1.cdata(groupGradient_L_hcpd_struct.diminfo{1,1}.models{20,1}.start:(groupGradient_L_hcpd_struct.diminfo{1,1}.models{20,1}.start-1) + groupGradient_L_hcpd_struct.diminfo{1,1}.models{20,1}.count,comp)=tempq';
                %sig_cii_1.cdata(sig_cii_1.cdata>0.01)=-999;
                sig_cii_qval=sig_cii_1;

                sig_idx_wb=find(0<sig_cii_qval.cdata(:,comp) & sig_cii_qval.cdata(:,comp)<0.01);
                sig_idx_wb_L_thal_neg{comp}=[sig_idx_wb];

                sig_cii_2.cdata(groupGradient_L_hcpd_struct.diminfo{1,1}.models{20,1}.start:(groupGradient_L_hcpd_struct.diminfo{1,1}.models{20,1}.start-1) + groupGradient_L_hcpd_struct.diminfo{1,1}.models{20,1}.count,comp)=slm_negAge_LThal_hcpd.t';
                sig_cii_t=sig_cii_2;

                %ciftisave(sig_cii_qval, [PATH '05_analysis_cmap_pmap/' results_dir '/negAge_L_thal_hcpd_qval.cmap.dscalar.nii' ], WB_COMMAND);
                %ciftisave(sig_cii_t, [ PATH '05_analysis_cmap_pmap/' results_dir '/negAge_L_thal_hcpd_tval.cmap.dscalar.nii' ], WB_COMMAND);
           
            elseif sig == 0 

                %uncorrected p
                p=1-tcdf(slm_negAge_LThal_hcpd.t, slm_negAge_LThal_hcpd.df);
                p = min([ p; 1-p ], [], 1);
                figure; histogram(p); 
                unique(p);
            end

            
        end
        %close all
        clear sig_cii_1 sig_cii_2 group_temp

    end    
    
close all 

%% test age effects: NEOMAP ==d=========================
    %% Vertexwise univariate analysis 
    % set model hcpd
    %%SurfStat model - set model & contrasts
    % 1. SurfStatLinMod (beta estimation)
    % 2. SurfStatT (contrast)
    % 3. SurfStatP (random field theory) or SurfStatQ (false discovery rate)

    % gradient_R = beta0 + beta1*age + beta2*group
    %GROUP = hcpd(group); %%female, male
    %M = 1 + AGE + GROUP

    for for_set_model = 1

        for for_set_model=1        

            AGE_HCPD = term(age_y_hcpd);
            SEX_HCPD = term(sex_hcpd);
            MEANFD_HCPD = term(meanFD_hcpd);


        end
        for for_set_model_age_correlations=1

            M_age_y_hcpd = 1 + AGE_HCPD; 
            M_age_y_hcpd_COV_SEXMEANFD = 1 + AGE_HCPD + SEX_HCPD + MEANFD_HCPD; 

        end 

    end
    %% test age effects
    
    for test_model_age_effect = 1
        for for_SurfStat_LinMod = 1
            %pmap1
            slm_L_thal_pmap1 = SurfStatLinMod(pmap_L_thal_ind_groupCmap_hcpd_align(:,:,1), M_age_y_hcpd_COV_SEXMEANFD, surf); 
            slm_L_thal_pmap1_posAge = SurfStatT(slm_L_thal_pmap1, age_y(idx_hcpd)); 
            slm_L_thal_pmap1_negAge = SurfStatT(slm_L_thal_pmap1, -age_y(idx_hcpd)); 

            %pmap2
            slm_L_thal_pmap2 = SurfStatLinMod(pmap_L_thal_ind_groupCmap_hcpd_align(:,:,2), M_age_y_hcpd_COV_SEXMEANFD, surf); 
            slm_L_thal_pmap2_posAge = SurfStatT(slm_L_thal_pmap2, age_y(idx_hcpd)); 
            slm_L_thal_pmap2_negAge = SurfStatT(slm_L_thal_pmap2, -age_y(idx_hcpd)); 

        end
        for RFT_CORRECTION =1
            %SurfStatP
            [pval_L_thal_pmap1_posAge, peak_L_thal_pmap1_posAge, clus_L_thal_pmap1_posAge, clusid_L_thal_pmap1_posAge] = SurfStatP(slm_L_thal_pmap1_posAge, logical(surf_roi), 0.025);

            figure; BoSurfStatView(pval_L_thal_pmap1_posAge.C, surf);  BoSurfStatColLim([0 0.01]); 

                %threshold by cluster size (>500 vertices)
                %check clus_R_thal_pmap1_posAge.nverts 
                a=sum(clus_L_thal_pmap1_posAge.nverts>500)
                temp_c = clusid_L_thal_pmap1_posAge
                temp_c(temp_c>a) = 0 %change number according to how many clusters to include 
                figure; BoSurfStatView(temp_c, surf);  BoSurfStatColLim([0 0.5]); 

                temp_c(1,temp_c>1) = 1;

                figure; BoSurfStatView(slm_L_thal_pmap1_posAge.t, surf);  BoSurfStatColLim([-5 5]); colormap(parula)
                sig_pmap =pmap_L_thal_ind_groupCmap_hcpd_align(:,:,1).*temp_c;
                sig_pmap(:,sig_pmap(1,:)==0)=[];

            [pval_L_thal_pmap1_negAge, peak_L_thal_pmap1_negAge, clus_L_thal_pmap1_negAge, clusid_L_thal_pmap1_negAge] = SurfStatP(slm_L_thal_pmap1_negAge, logical(surf_roi), 0.025);
            %figure; SurfStatView(pval_L_thal_pmap1_negAge.C, surf);  SurfStatColLim([0 0.05]); 
            figure; BoSurfStatView(pval_L_thal_pmap1_negAge.C, surf);  BoSurfStatColLim([0 0.05]); 

                %threshold by cluster size (>500 vertices)
                %check clus_R_thal_pmap1_posAge.nverts 
                a=sum(clus_L_thal_pmap1_negAge.nverts>500)
                temp_c = clusid_L_thal_pmap1_negAge
                temp_c(temp_c>a) = 0 %change number according to how many clusters to include 
                figure; BoSurfStatView(temp_c, surf);  BoSurfStatColLim([0 0.5]); 


                %figure; BoSurfStatView(group_pmap1_L_thal_ind_groupCmap_hcpd_align, surf); % BoSurfStatColLim([0 0.5]); 
                %figure; BoSurfStatView(group_pmap2_L_thal_ind_groupCmap_hcpd_align, surf); % BoSurfStatColLim([0 0.5]); 

                temp_c(1,temp_c>1) = 1;

                figure; BoSurfStatView(slm_L_thal_pmap1_negAge.t, surf);  BoSurfStatColLim([-5 5]); colormap(parula)
                sig_pmap =pmap_L_thal_ind_groupCmap_hcpd_align(:,:,1).*temp_c;
                sig_pmap(:,sig_pmap(1,:)==0)=[];

            [pval_L_thal_pmap2_posAge, peak_L_thal_pmap2_posAge, clus_L_thal_pmap2_posAge, clusid_L_thal_pmap2_posAge] = SurfStatP(slm_L_thal_pmap2_posAge, logical(surf_roi), 0.025);
            %figure; SurfStatView(pval_L_thal_pmap2_posAge.C, surf);  SurfStatColLim([0 0.05]); 
            figure; BoSurfStatView(pval_L_thal_pmap2_posAge.C, surf);  BoSurfStatColLim([0 0.05]); 

                %threshold by cluster size (>500 vertices)
                %check clus_R_thal_pmap1_posAge.nverts 
                a=sum(clus_L_thal_pmap2_posAge.nverts>500)
                temp_c = clusid_L_thal_pmap2_posAge
                temp_c(temp_c>a) = 0 %change number according to how many clusters to include 
                figure; BoSurfStatView(temp_c, surf);  BoSurfStatColLim([0 0.5]); 

                temp_c(1,temp_c>1) = 1;

                figure; BoSurfStatView(slm_L_thal_pmap2_posAge.t, surf);  BoSurfStatColLim([-5 5]); colormap(parula)
                sig_pmap =pmap_L_thal_ind_groupCmap_hcpd_align(:,:,2).*temp_c;
                sig_pmap(:,sig_pmap(1,:)==0)=[];


            [pval_L_thal_pmap2_negAge, peak_L_thal_pmap2_negAge, clus_L_thal_pmap2_negAge, clusid_L_thal_pmap2_negAge] = SurfStatP(slm_L_thal_pmap2_negAge, logical(surf_roi), 0.025);
            %figure; SurfStatView(pval_L_thal_pmap2_negAge.C, surf);  SurfStatColLim([0 0.05]); 
            figure; BoSurfStatView(pval_L_thal_pmap2_negAge.C, surf);  BoSurfStatColLim([0 0.05]); 

                %threshold by cluster size (>500 vertices)
                %check clus_R_thal_pmap1_posAge.nverts 
                a=sum(clus_L_thal_pmap2_negAge.nverts>500)
                temp_c = clusid_L_thal_pmap2_negAge
                temp_c(temp_c>a) = 0 %change number according to how many clusters to include 
                figure; BoSurfStatView(temp_c, surf);  BoSurfStatColLim([0 0.5]); 

                temp_c(1,temp_c>1) = 1;

                figure; BoSurfStatView(slm_L_thal_pmap2_negAge.t, surf);  BoSurfStatColLim([-5 5]); colormap(parula)
                sig_pmap =pmap_L_thal_ind_groupCmap_hcpd_align(:,:,2).*temp_c;
                sig_pmap(:,sig_pmap(1,:)==0)=[];

        end

        close all
        clear sig_cii_1 sig_cii_2 group_temp

    end     
