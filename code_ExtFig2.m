
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
    %% 02) upload subject id 
    
    t=readtable('./files/Sleep-Wake_subject_session.txt')
    idx_asd = t.subject %58
    ses_asd = t.session
    idx_unique = unique(idx_asd) %29
    
    % sleep / wake demographics     
    data = importfile_demo_sleepwake(['./files/Sleep-Wake_demo.xlsx']);
    age_awakesleep = [data.AgeAwake data.AgeAsleep]
    fdj_awakesleep = [data.FDJ_Awake data.FDJ_Sleep] 
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

%% Align CMAPs to group CMAP - save individual aligned CMAP
    load('sourceData_ExtFig2_Sleep-Wake.mat')
    for for_prepare_source_data= 1 %do not run, upload source data
        %% Group cmap
        for for_grouptemplate_10k =1
    % 
    %         PATH2               = '/local_raid2/03_user/shinwon/02_data/asd_sleep-wake/';
    %         GROUP_GRADIENTS           = [ PATH2 '03_group_gradient/' ];
    %         
    %         cii_awake = ciftiopen([ GROUP_GRADIENTS 'functimeseries_sampler_10k_awake/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.cmap.dscalar.nii' ], WB_COMMAND);   
    %         cii_awake_L_thal = cii_awake.cdata(cii_awake.diminfo{1, 1}.models{20, 1}.start:(cii_awake.diminfo{1,1}.models{20,1}.start-1) + cii_awake.diminfo{1,1}.models{20,1}.count,:);
    %         groupGradient_L_thal_awake(:, :) = cii_awake_L_thal;
    % 
    %         cii_awake = ciftiopen([ GROUP_GRADIENTS 'functimeseries_sampler_10k_awake/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.cmap.dscalar.nii' ], WB_COMMAND);   
    %         cii_awake_R_thal = cii_awake.cdata(cii_awake.diminfo{1, 1}.models{21, 1}.start:(cii_awake.diminfo{1,1}.models{21,1}.start-1) + cii_awake.diminfo{1,1}.models{21,1}.count,:);
    %         groupGradient_R_thal_awake(:, :) = cii_awake_R_thal;
    % 
    %         cii_sleep = ciftiopen([ GROUP_GRADIENTS 'functimeseries_sampler_10k_sleep/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.cmap.dscalar.nii' ], WB_COMMAND);   
    %         cii_sleep_L_thal = cii_sleep.cdata(cii_sleep.diminfo{1, 1}.models{20, 1}.start:(cii_sleep.diminfo{1,1}.models{20,1}.start-1) + cii_sleep.diminfo{1,1}.models{20,1}.count,:);
    %         groupGradient_L_thal_sleep(:, :) = cii_sleep_L_thal;
    % 
    %         cii_sleep = ciftiopen([ GROUP_GRADIENTS 'functimeseries_sampler_10k_awake/atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.cmap.dscalar.nii' ], WB_COMMAND);   
    %         cii_sleep_R_thal = cii_sleep.cdata(cii_sleep.diminfo{1, 1}.models{21, 1}.start:(cii_sleep.diminfo{1,1}.models{21,1}.start-1) + cii_sleep.diminfo{1,1}.models{21,1}.count,:);
    %         groupGradient_R_thal_sleep(:, :) = cii_sleep_R_thal;


        end

        %% Read (unaligned) individual gradients (+ do alignment)

    for for_gradient_directory = 1
%         IND_GRAD         = [ PATH2 '02_individual_gradient/functimeseries_sampler_10k/' ];
    end
    
    for for_read_unaligned_ind_gradients =1

        %% individual gradients (in atlasroi) aligned to Group gradient (5-22y)

        %10k - thalamus gradients - awake
        parfor i = 1:length(idx_unique)
%             i 
%             cii_L_awake = ciftiopen([ IND_GRAD idx_unique{i} '_ses-awake/' 'atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.cmap.dscalar.nii' ], WB_COMMAND);   
%             cii_L_thal_awake = cii_L_awake.cdata(cii_L_awake.diminfo{1,1}.models{20,1}.start:(cii_L_awake.diminfo{1,1}.models{20,1}.start-1) + cii_L_awake.diminfo{1,1}.models{20,1}.count,:);
%             grad_L_thal_ind_noalign_awake(i, :, :) = cii_L_thal_awake;    
%            
%             cii_R_awake = ciftiopen([ IND_GRAD idx_unique{i} '_ses-awake/' 'atlasroi.Neo-Subcortical-2mm.10k_fs_LR.49.cmap.dscalar.nii' ], WB_COMMAND);   
%             cii_R_thal_awake = cii_R_awake.cdata(cii_R_awake.diminfo{1,1}.models{21,1}.start:(cii_R_awake.diminfo{1,1}.models{21,1}.start-1) + cii_R_awake.diminfo{1,1}.models{21,1}.count,:);
%             grad_R_thal_ind_noalign_awake(i, :, :) = cii_R_thal_awake;    
%             
%             %procrustes alignment    
%             [d_L,Z_L,transform_L] = procrustes(groupGradient_L_thal_awake,cii_L_thal_awake, 'scaling', false);
%             [d_R,Z_R,transform_R] = procrustes(groupGradient_R_thal_awake,cii_R_thal_awake, 'scaling', false);
% 
%             cii_L_awake.cdata=zeros(50592,10);
%             cii_L_awake.cdata(cii_L_awake.diminfo{1,1}.models{20,1}.start:(cii_L_awake.diminfo{1,1}.models{20,1}.start-1) + cii_L_awake.diminfo{1,1}.models{20,1}.count,:)=Z_L;
%             cii_aligned_L_awake=cii_L_awake;
%             
%             cii_R_awake.cdata=zeros(50592,10);
%             cii_R_awake.cdata(cii_R_awake.diminfo{1,1}.models{21,1}.start:(cii_R_awake.diminfo{1,1}.models{21,1}.start-1) + cii_R_awake.diminfo{1,1}.models{21,1}.count,:)=Z_R;
%             cii_aligned_R_awake=cii_R_awake;  
% 
%             %ciftisave(cii_aligned_L_awake, [ IND_GRAD idx_unique{i} '_ses-awake/' 'atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.aligned_cmap.dscalar.nii' ], WB_COMMAND);  
%             grad_L_thal_ind_align_awake(i, :, :) = Z_L;
% 
%             %ciftisave(cii_aligned_R_awake, [ IND_GRAD idx_unique{i} '_ses-awake/' 'atlasroi.Neo-Subcortical-2mm.10k_fs_LR.49.aligned_cmap.dscalar.nii' ], WB_COMMAND);  
%             grad_R_thal_ind_align_awake(i, :, :) = Z_R;

        end
        
          %10k - thalamus gradients - sleep
        parfor i = 1:length(idx_unique)
%             i 
%             cii_L_sleep = ciftiopen([ IND_GRAD idx_unique{i} '_ses-sleep/' 'atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.cmap.dscalar.nii' ], WB_COMMAND);   
%             cii_L_thal_sleep = cii_L_sleep.cdata(cii_L_sleep.diminfo{1,1}.models{20,1}.start:(cii_L_sleep.diminfo{1,1}.models{20,1}.start-1) + cii_L_sleep.diminfo{1,1}.models{20,1}.count,:);
%             grad_L_thal_ind_noalign_sleep(i, :, :) = cii_L_thal_sleep;    
%            
%             cii_R_sleep = ciftiopen([ IND_GRAD idx_unique{i} '_ses-sleep/' 'atlasroi.Neo-Subcortical-2mm.10k_fs_LR.49.cmap.dscalar.nii' ], WB_COMMAND);   
%             cii_R_thal_sleep = cii_R_sleep.cdata(cii_R_sleep.diminfo{1,1}.models{21,1}.start:(cii_R_sleep.diminfo{1,1}.models{21,1}.start-1) + cii_R_sleep.diminfo{1,1}.models{21,1}.count,:);
%             grad_R_thal_ind_noalign_sleep(i, :, :) = cii_R_thal_sleep;    
%             
%             %procrustes alignment    
%             [d_L,Z_L,transform_L] = procrustes(groupGradient_L_thal_sleep,cii_L_thal_sleep, 'scaling', false);
%             [d_R,Z_R,transform_R] = procrustes(groupGradient_R_thal_sleep,cii_R_thal_sleep, 'scaling', false);
% 
%             cii_L_sleep.cdata=zeros(50592,10);
%             cii_L_sleep.cdata(cii_L_sleep.diminfo{1,1}.models{20,1}.start:(cii_L_sleep.diminfo{1,1}.models{20,1}.start-1) + cii_L_sleep.diminfo{1,1}.models{20,1}.count,:)=Z_L;
%             cii_aligned_L_sleep=cii_L_sleep;
%             
%             cii_R_sleep.cdata=zeros(50592,10);
%             cii_R_sleep.cdata(cii_R_sleep.diminfo{1,1}.models{21,1}.start:(cii_R_sleep.diminfo{1,1}.models{21,1}.start-1) + cii_R_sleep.diminfo{1,1}.models{21,1}.count,:)=Z_R;
%             cii_aligned_R_sleep=cii_R_sleep;  
% 
%             %ciftisave(cii_aligned_L_sleep, [ IND_GRAD idx_unique{i} '_ses-sleep/' 'atlasroi.Neo-Subcortical-2mm.10k_fs_LR.10.aligned_cmap.dscalar.nii' ], WB_COMMAND);  
%             grad_L_thal_ind_align_sleep(i, :, :) = Z_L;
% 
%             %ciftisave(cii_aligned_R_sleep, [ IND_GRAD idx_unique{i} '_ses-sleep/' 'atlasroi.Neo-Subcortical-2mm.10k_fs_LR.49.aligned_cmap.dscalar.nii' ], WB_COMMAND);  
%             grad_R_thal_ind_align_sleep(i, :, :) = Z_R;

        end
        
    end
    end
    
    %% Extended FIGURE 2A: check spatial correlation of group CMAPs

    [r p] = corr(groupGradient_L_thal_awake(:,1), groupGradient_L_thal_sleep(:,1)) 
    [r p] = corr(groupGradient_L_thal_awake(:,2), groupGradient_L_thal_sleep(:,2)) 
    
    %% Extended FIGURE 2B: check spatial correlation of individual aligned cmaps 
        grad1_L_thal_ind_align_awake = grad_L_thal_ind_align_awake(:,:,1);
        grad2_L_thal_ind_align_awake = grad_L_thal_ind_align_awake(:,:,2);

        grad1_L_thal_ind_align_sleep = grad_L_thal_ind_align_sleep(:,:,1);
        grad2_L_thal_ind_align_sleep = grad_L_thal_ind_align_sleep(:,:,2);

        for i = 1:length(idx_unique)
            i 
            [r p] = corr(grad1_L_thal_ind_align_awake(i,:)', grad1_L_thal_ind_align_sleep(i,:)') ;

            r_ind_aligned_cmap(i)= r;
            p_ind_aligned_cmap(i)= p;

            [r2 p2] = corr(grad2_L_thal_ind_align_awake(i,:)', grad2_L_thal_ind_align_sleep(i,:)') ;

            r2_ind_aligned_cmap(i)= r2;
            p2_ind_aligned_cmap(i)= p2;

        end

        mean(r_ind_aligned_cmap)
        std(r_ind_aligned_cmap)
        mean(r2_ind_aligned_cmap)
        std(r2_ind_aligned_cmap)

        %% Gradient 1
        r_ind_aligned_cmap_s = sort(r_ind_aligned_cmap)

        figure1 = figure('Color',[1 1 1]);

        axes1 = axes('Parent',figure1);
        hold(axes1,'on');

        scatter([1:29] , r_ind_aligned_cmap_s,...
            'MarkerFaceColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
            'MarkerEdgeColor',[0.149019607843137 0.149019607843137 0.149019607843137]);       

        plot(r_ind_aligned_cmap_s,'LineStyle','--','Color',[0.8 0.8 0.8]);
        ylim(axes1,[0.7 1]);
        hold(axes1,'off');
        set(axes1,'XMinorGrid','on','YGrid','on')
        ylim([0.7 1]);set(gcf, 'color', 'w');


        figure1 = figure('Color',[1 1 1]);

        axes1 = axes('Parent',figure1);
        hold(axes1,'on');

        histogram(r_ind_aligned_cmap,'Parent',axes1,...
            'FaceColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
            'BinMethod','auto');

         xlim(axes1,[0.7 1]);
         ylim(axes1,[0 25]);
        box(axes1,'on');
        hold(axes1,'off');
        set(axes1,'YMinorGrid','on');

        %% Gradient 2
        r2_ind_aligned_cmap_s = sort(r2_ind_aligned_cmap)

        figure1 = figure('Color',[1 1 1]);

        axes1 = axes('Parent',figure1);
        hold(axes1,'on');

        scatter([1:29] , r2_ind_aligned_cmap_s,...
            'MarkerFaceColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
            'MarkerEdgeColor',[0.149019607843137 0.149019607843137 0.149019607843137]);       

        plot(r2_ind_aligned_cmap_s,'LineStyle','--','Color',[0.8 0.8 0.8]);
        ylim(axes1,[0.7 1]);
        hold(axes1,'off');
        set(axes1,'XMinorGrid','on','YGrid','on')
        ylim([0.7 1]);set(gcf, 'color', 'w');


        figure1 = figure('Color',[1 1 1]);

        axes1 = axes('Parent',figure1);
        hold(axes1,'on');

        histogram(r2_ind_aligned_cmap,'Parent',axes1,...
            'FaceColor',[0.501960784313725 0.501960784313725 0.501960784313725],...
            'BinMethod','auto');

         xlim(axes1,[0.7 1]);
         ylim(axes1,[0 25]);
        box(axes1,'on');
        hold(axes1,'off');
        set(axes1,'YMinorGrid','on');

%% Extended FIGURE 2C: surfstat - mixed effects model 

    Y1= [grad1_L_thal_ind_align_sleep; grad1_L_thal_ind_align_awake];
    Y2= [grad2_L_thal_ind_align_sleep; grad2_L_thal_ind_align_awake];

    for for_model = 1
        subj = ([1:29, 1:29])';
        num2str(subj)
        Subj = term(var2fac(subj));
        state = [repmat({'sleep'}, 1,29), repmat({'wake'}, 1,29)]'; 
        State = term(state);
        
        age =  [ data.AgeAsleep; data.AgeAwake]
        Age = term(age);

        fdj =  [ data.FDJ_Sleep; data.FDJ_Awake]
        FDJ = term(fdj);        

        M = 1 + State + random(Subj) + I;
        M_age = 1 + State + Age + random(Subj) + I;
        M_agefdj = 1 + State + Age + FDJ + random(Subj) + I;
        
        C_sleepwake = State.sleep - State.wake; 
        C_wakesleep = State.wake - State.sleep; 
        
        %figure; image(M)
    end
    
    for for_grad1 = 1

        slm_grad1 = SurfStatLinMod(Y1, M);
        slm_grad1_sleepwake = SurfStatT(slm_grad1, C_sleepwake);
        slm_grad1_wakesleep = SurfStatT(slm_grad1, C_wakesleep);

        temp = cii_awake;
        temp.cdata = zeros(50592, 10); 
        temp.cdata(temp.diminfo{1, 1}.models{20, 1}.start: temp.diminfo{1, 1}.models{20, 1}.start + temp.diminfo{1, 1}.models{20, 1}.count -1 ) = slm_grad1_sleepwake.t;
        %ciftisave(temp, 'tval_grad1_sleepwake.cmap.dscalar.nii' );
        
        temp = cii_awake;
        temp.cdata = zeros(50592, 10); 
        temp.cdata(temp.diminfo{1, 1}.models{20, 1}.start: temp.diminfo{1, 1}.models{20, 1}.start + temp.diminfo{1, 1}.models{20, 1}.count -1 ) = slm_grad1_wakesleep.t;
        %ciftisave(temp, 'tval_grad1_wakesleep.cmap.dscalar.nii' );

        %FDR correction
        qval_grad1_sleepwake=SurfStatQ(slm_grad1_sleepwake);
        sig_qval_grad1_sleepwake=sum(qval_grad1_sleepwake.Q<0.05);
        sig_idx_grad1_sleepwake=find(qval_grad1_sleepwake.Q<0.05);
        
        if sig_qval_grad1_sleepwake ~= 0 
            tempq=zeros(1, size(qval_grad1_sleepwake.Q,2));             
            tempq(1, qval_grad1_sleepwake.Q<0.05)=1;

            temp = cii_awake;
            temp.cdata = zeros(50592, 10); 

            temp.cdata(temp.diminfo{1,1}.models{20,1}.start:(temp.diminfo{1,1}.models{20,1}.start-1) + temp.diminfo{1,1}.models{20,1}.count, 1) = tempq';            

            %ciftisave(temp, 'qval_grad1_sleepwake.cmap.dscalar.nii');
            
            grad_values_sigRegion = mean(double(Y1(:,sig_idx_grad1_sleepwake)),2);     
            grad_values_sigRegion_w = [grad_values_sigRegion(1:29), grad_values_sigRegion(30:end)]
            
            figure(); 
            coordLineStyle = 'k.'
            boxplot(grad_values_sigRegion_w, 'Symbol', coordLineStyle); hold on; 
            parallelcoords(grad_values_sigRegion_w, 'Color', 0.7*[1 1 1], 'LineStyle', '-', ... 
                'Marker', '.', 'MarkerSize', 10);
            ylim([-0.011 0.007]);
            %saveas(gcf, 'scatter_qval_grad1_sleepwake.png')  

        end       

        qval_grad1_wakesleep=SurfStatQ(slm_grad1_wakesleep);
        sig_qval_grad1_wakesleep=sum(qval_grad1_wakesleep.Q<0.05);
        sig_idx_grad1_wakesleep=find(qval_grad1_wakesleep.Q<0.05);

        if sig_qval_grad1_wakesleep ~= 0
            tempq=zeros(1, size(qval_grad1_wakesleep.Q,2));             
            tempq(1, qval_grad1_wakesleep.Q<0.05)=1;

            temp = cii_awake;
            temp.cdata = zeros(50592, 10); 

            temp.cdata(temp.diminfo{1,1}.models{20,1}.start:(temp.diminfo{1,1}.models{20,1}.start-1) + temp.diminfo{1,1}.models{20,1}.count, 1) = tempq';            

            %ciftisave(temp, 'qval_grad1_wakesleep.cmap.dscalar.nii');

            grad_values_sigRegion = mean(double(Y1(:,sig_idx_grad1_wakesleep)),2);     
            grad_values_sigRegion_w = [grad_values_sigRegion(1:29), grad_values_sigRegion(30:end)]
            
            figure(); 
            coordLineStyle = 'k.'
            boxplot(grad_values_sigRegion_w, 'Symbol', coordLineStyle); hold on; 
            parallelcoords(grad_values_sigRegion_w, 'Color', 0.7*[1 1 1], 'LineStyle', '-', ... 
                'Marker', '.', 'MarkerSize', 10);
            ylim([-0.011 0.007]);
            %saveas(gcf, 'scatter_qval_grad1_wakesleep.png')  
        end       
        
        
    end
    
    for for_grad2 = 1
        slm_grad2 = SurfStatLinMod(Y2, M);
        slm_grad2_sleepwake = SurfStatT(slm_grad2, C_sleepwake);
        slm_grad2_wakesleep = SurfStatT(slm_grad2, C_wakesleep);

        temp = cii_awake;
        temp.cdata = zeros(50592, 10); 
        temp.cdata(temp.diminfo{1, 1}.models{20, 1}.start: temp.diminfo{1, 1}.models{20, 1}.start + temp.diminfo{1, 1}.models{20, 1}.count -1 ) = slm_grad2_sleepwake.t;
        %ciftisave(temp, 'tval_grad2_sleepwake.cmap.dscalar.nii' );
        
        temp = cii_awake;
        temp.cdata = zeros(50592, 10); 
        temp.cdata(temp.diminfo{1, 1}.models{20, 1}.start: temp.diminfo{1, 1}.models{20, 1}.start + temp.diminfo{1, 1}.models{20, 1}.count -1 ) = slm_grad2_wakesleep.t;
        %ciftisave(temp, 'tval_grad2_wakesleep.cmap.dscalar.nii' );

        %FDR correction
        qval_grad2_sleepwake=SurfStatQ(slm_grad2_sleepwake);
        sig_qval_grad2_sleepwake=sum(qval_grad2_sleepwake.Q<0.05);
        sig_idx_grad2_sleepwake=find(qval_grad2_sleepwake.Q<0.05);

        if sig_qval_grad2_sleepwake ~= 0
            tempq=zeros(1, size(qval_grad2_sleepwake.Q,2));             
            tempq(1, qval_grad2_sleepwake.Q<0.05)=1;

            temp = cii_awake;
            temp.cdata = zeros(50592, 10); 

            temp.cdata(temp.diminfo{1,1}.models{20,1}.start:(temp.diminfo{1,1}.models{20,1}.start-1) + temp.diminfo{1,1}.models{20,1}.count, 1) = tempq';            

            %ciftisave(temp, 'qval_grad2_sleepwake.cmap.dscalar.nii');

            grad_values_sigRegion = mean(double(Y2(:,sig_idx_grad2_sleepwake)),2);     
            grad_values_sigRegion_w = [grad_values_sigRegion(1:29), grad_values_sigRegion(30:end)]
            
            figure(); 
            coordLineStyle = 'k.'
            boxplot(grad_values_sigRegion_w, 'Symbol', coordLineStyle); hold on; 
            parallelcoords(grad_values_sigRegion_w, 'Color', 0.7*[1 1 1], 'LineStyle', '-', ... 
                'Marker', '.', 'MarkerSize', 10);
            ylim([-0.011 0.007]);
            %saveas(gcf, 'scatter_qval_grad2_sleepwake.png')  
            

        end
        
    
        %FDR correction
        qval_grad2_wakesleep=SurfStatQ(slm_grad2_wakesleep);
        sig_qval_grad2_wakesleep=sum(qval_grad2_wakesleep.Q<0.05);
        sig_idx_grad2_wakesleep=find(qval_grad2_wakesleep.Q<0.05);

        if sig_qval_grad2_wakesleep ~= 0
            tempq=zeros(1, size(qval_grad2_wakesleep.Q,2));             
            tempq(1, qval_grad2_wakesleep.Q<0.05)=1;

            temp = cii_awake;
            temp.cdata = zeros(50592, 10); 

            temp.cdata(temp.diminfo{1,1}.models{20,1}.start:(temp.diminfo{1,1}.models{20,1}.start-1) + temp.diminfo{1,1}.models{20,1}.count, 1) = tempq';            

            %ciftisave(temp, 'qval_grad2_wakesleep.cmap.dscalar.nii');

            grad_values_sigRegion = mean(double(Y2(:,sig_idx_grad2_wakesleep)),2);     
            grad_values_sigRegion_w = [grad_values_sigRegion(1:29), grad_values_sigRegion(30:end)]
            
            figure(); 
            coordLineStyle = 'k.'
            boxplot(grad_values_sigRegion_w, 'Symbol', coordLineStyle); hold on; 
            parallelcoords(grad_values_sigRegion_w, 'Color', 0.7*[1 1 1], 'LineStyle', '-', ... 
                'Marker', '.', 'MarkerSize', 10);
            ylim([-0.011 0.007]);
            %saveas(gcf, 'scatter_qval_grad2_wakesleep.png')  
            
        end              
        
    end
    