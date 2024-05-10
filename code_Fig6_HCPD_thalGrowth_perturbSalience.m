
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
% prepare input files -> refer to code_Fig6_prepare_inputs_HCPD_perturbation.m script
% run GNM scripts from https://github.com/StuartJO/GenerativeNetworkModel

%% Figure 6: upload data
load('sourceData_Fig6_HCPD_perturbSalience.mat')

    %% B = an adjacency matrix of the generated network
    
    for for_optimMdl = 1
        adjmat_0_noPerturb = aa_0_noPerturb.OptimMdl{1, 1}.min_maxKS.adjmat{1, 1};
        adjmat_1_Perturb = aa_1_Perturb.OptimMdl{1, 1}.min_maxKS.adjmat{1, 1};
        adjmat_2_Perturb = aa_2_Perturb.OptimMdl{1, 1}.min_maxKS.adjmat{1, 1};
        adjmat_3_Perturb = aa_3_Perturb.OptimMdl{1, 1}.min_maxKS.adjmat{1, 1};
        adjmat_4_Perturb = aa_4_Perturb.OptimMdl{1, 1}.min_maxKS.adjmat{1, 1};
        %adjmat_4_Static = aa_4_Static.OptimMdl{1, 1}.min_maxKS.adjmat{1, 1};
    end
        
    figure; imagesc(adjmat_0_noPerturb);colorbar; caxis([-1 1]); colormap(flipud(rdwhbl))%colormap(hcp_colormap(:,1:3))
    figure; imagesc(adjmat_1_Perturb);colorbar; caxis([-1 1]); colormap(flipud(rdwhbl))%colormap(hcp_colormap(:,1:3))
    figure; imagesc(adjmat_2_Perturb);colorbar; caxis([-1 1]); colormap(flipud(rdwhbl))%colormap(hcp_colormap(:,1:3))
    figure; imagesc(adjmat_3_Perturb);colorbar; caxis([-1 1]); colormap(flipud(rdwhbl))
    figure; imagesc(adjmat_4_Perturb);colorbar; caxis([-1 1]); colormap(flipud(rdwhbl))
    %figure; imagesc(adjmat_4_Static);colorbar; caxis([-1 1]); %colormap(flipud(rdwhbl))%colormap(hcp_colormap(:,1:3))

    %% Figure 6B: segregation index - bar graph 

    %project_detection_community_boview(yeo7_in_sch200,'', schaefer_200_label(1:10242,1)',surfL,1);

    yeo7_in_sch200_vis = find(yeo7_in_sch200 == 1); 
    yeo7_in_sch200_som = find(yeo7_in_sch200 == 2); 
    yeo7_in_sch200_dan = find(yeo7_in_sch200 == 3); 
    yeo7_in_sch200_sal = find(yeo7_in_sch200 == 4); 
    yeo7_in_sch200_lim = find(yeo7_in_sch200 == 5); 
    yeo7_in_sch200_fpn = find(yeo7_in_sch200 == 6); 
    yeo7_in_sch200_dmn = find(yeo7_in_sch200 == 7);

    yeo7_in_sch200_ext = [yeo7_in_sch200_vis; yeo7_in_sch200_som;yeo7_in_sch200_dan];

    % segregation index = abs(within - between)
    for for_btnwtn = 1
        btn_saldan=mean(adjmat_0_noPerturb(yeo7_in_sch200_sal, yeo7_in_sch200_dan),'all');
        btn_salext=mean(adjmat_0_noPerturb(yeo7_in_sch200_sal, yeo7_in_sch200_ext),'all');
        btn_saldmn=mean(adjmat_0_noPerturb(yeo7_in_sch200_sal, yeo7_in_sch200_dmn),'all');
        wtn_sal=mean(adjmat_0_noPerturb(yeo7_in_sch200_sal, yeo7_in_sch200_sal),'all');
        wtn_dan=mean(adjmat_0_noPerturb(yeo7_in_sch200_dan, yeo7_in_sch200_dan),'all');
        wtn_dmn=mean(adjmat_0_noPerturb(yeo7_in_sch200_dmn, yeo7_in_sch200_dmn),'all');
        wtn_ext=mean(adjmat_0_noPerturb(yeo7_in_sch200_ext, yeo7_in_sch200_ext),'all');

%         r_saldan_0 = mean([wtn_sal wtn_dan])-btn_saldan;
%         r_saldmn_0 = mean([wtn_sal wtn_dmn])-btn_saldmn;
%         r_saldansen_0 = mean([wtn_sal wtn_ext])-btn_salext;

        r_saldan_0 = abs(mean([wtn_sal wtn_dan])-btn_saldan);
        r_saldmn_0 = abs(mean([wtn_sal wtn_dmn])-btn_saldmn);
        r_saldansen_0 = abs(mean([wtn_sal wtn_ext])-btn_salext);


        btn_saldan=mean(adjmat_1_Perturb(yeo7_in_sch200_sal, yeo7_in_sch200_dan),'all');
        btn_salext=mean(adjmat_1_Perturb(yeo7_in_sch200_sal, yeo7_in_sch200_ext),'all');
        btn_saldmn=mean(adjmat_1_Perturb(yeo7_in_sch200_sal, yeo7_in_sch200_dmn),'all');
        wtn_sal=mean(adjmat_1_Perturb(yeo7_in_sch200_sal, yeo7_in_sch200_sal),'all');
        wtn_dan=mean(adjmat_1_Perturb(yeo7_in_sch200_dan, yeo7_in_sch200_dan),'all');
        wtn_dmn=mean(adjmat_1_Perturb(yeo7_in_sch200_dmn, yeo7_in_sch200_dmn),'all');
        wtn_ext=mean(adjmat_1_Perturb(yeo7_in_sch200_ext, yeo7_in_sch200_ext),'all');
% 
%         r_saldan_1 = mean([wtn_sal wtn_dan])-btn_saldan;
%         r_saldmn_1 = mean([wtn_sal wtn_dmn])-btn_saldmn;
%         r_saldansen_1 = mean([wtn_sal wtn_ext])-btn_salext;

        r_saldan_1 = abs(mean([wtn_sal wtn_dan])-btn_saldan);
        r_saldmn_1 = abs(mean([wtn_sal wtn_dmn])-btn_saldmn);
        r_saldansen_1 = abs(mean([wtn_sal wtn_ext])-btn_salext);

        btn_saldan=mean(adjmat_2_Perturb(yeo7_in_sch200_sal, yeo7_in_sch200_dan),'all');
        btn_salext=mean(adjmat_2_Perturb(yeo7_in_sch200_sal, yeo7_in_sch200_ext),'all');
        btn_saldmn=mean(adjmat_2_Perturb(yeo7_in_sch200_sal, yeo7_in_sch200_dmn),'all');
        wtn_sal=mean(adjmat_2_Perturb(yeo7_in_sch200_sal, yeo7_in_sch200_sal),'all');
        wtn_dan=mean(adjmat_2_Perturb(yeo7_in_sch200_dan, yeo7_in_sch200_dan),'all');
        wtn_dmn=mean(adjmat_2_Perturb(yeo7_in_sch200_dmn, yeo7_in_sch200_dmn),'all');
        wtn_ext=mean(adjmat_2_Perturb(yeo7_in_sch200_ext, yeo7_in_sch200_ext),'all');
% 
%         r_saldan_2 = mean([wtn_sal wtn_dan])-btn_saldan;
%         r_saldmn_2 = mean([wtn_sal wtn_dmn])-btn_saldmn;
%         r_saldansen_2 = mean([wtn_sal wtn_ext])-btn_salext;

        r_saldan_2 = abs(mean([wtn_sal wtn_dan])-btn_saldan);
        r_saldmn_2 = abs(mean([wtn_sal wtn_dmn])-btn_saldmn);
        r_saldansen_2 = abs(mean([wtn_sal wtn_ext])-btn_salext);


        btn_saldan=mean(adjmat_3_Perturb(yeo7_in_sch200_sal, yeo7_in_sch200_dan),'all');
        btn_salext=mean(adjmat_3_Perturb(yeo7_in_sch200_sal, yeo7_in_sch200_ext),'all');
        btn_saldmn=mean(adjmat_3_Perturb(yeo7_in_sch200_sal, yeo7_in_sch200_dmn),'all');
        wtn_sal=mean(adjmat_3_Perturb(yeo7_in_sch200_sal, yeo7_in_sch200_sal),'all');
        wtn_dan=mean(adjmat_3_Perturb(yeo7_in_sch200_dan, yeo7_in_sch200_dan),'all');
        wtn_dmn=mean(adjmat_3_Perturb(yeo7_in_sch200_dmn, yeo7_in_sch200_dmn),'all');
        wtn_ext=mean(adjmat_3_Perturb(yeo7_in_sch200_ext, yeo7_in_sch200_ext),'all');

%         r_saldan_3 = mean([wtn_sal wtn_dan])-btn_saldan;
%         r_saldmn_3 = mean([wtn_sal wtn_dmn])-btn_saldmn;
%         r_saldansen_3 = mean([wtn_sal wtn_ext])-btn_salext;

        r_saldan_3 = abs(mean([wtn_sal wtn_dan])-btn_saldan);
        r_saldmn_3 = abs(mean([wtn_sal wtn_dmn])-btn_saldmn);
        r_saldansen_3 = abs(mean([wtn_sal wtn_ext])-btn_salext);

        btn_saldan=mean(adjmat_4_Perturb(yeo7_in_sch200_sal, yeo7_in_sch200_dan),'all');
        btn_salext=mean(adjmat_4_Perturb(yeo7_in_sch200_sal, yeo7_in_sch200_ext),'all');
        btn_saldmn=mean(adjmat_4_Perturb(yeo7_in_sch200_sal, yeo7_in_sch200_dmn),'all');
        wtn_sal=mean(adjmat_4_Perturb(yeo7_in_sch200_sal, yeo7_in_sch200_sal),'all');
        wtn_dan=mean(adjmat_4_Perturb(yeo7_in_sch200_dan, yeo7_in_sch200_dan),'all');
        wtn_dmn=mean(adjmat_4_Perturb(yeo7_in_sch200_dmn, yeo7_in_sch200_dmn),'all');
        wtn_ext=mean(adjmat_4_Perturb(yeo7_in_sch200_ext, yeo7_in_sch200_ext),'all');
% 
%         r_saldan_4 = mean([wtn_sal wtn_dan])-btn_saldan;
%         r_saldmn_4 = mean([wtn_sal wtn_dmn])-btn_saldmn;
%         r_saldansen_4 = mean([wtn_sal wtn_ext])-btn_salext;
        
        r_saldan_4 = abs(mean([wtn_sal wtn_dan])-btn_saldan);
        r_saldmn_4 = abs(mean([wtn_sal wtn_dmn])-btn_saldmn);
        r_saldansen_4 = abs(mean([wtn_sal wtn_ext])-btn_salext);

    end

       a2= [r_saldmn_0, r_saldmn_1, r_saldmn_2, r_saldmn_3, r_saldmn_4];   
       a1= [r_saldansen_0, r_saldansen_1, r_saldansen_2, r_saldansen_3, r_saldansen_4];

       %figure; bar( [ a1; a2 ; a3] ); ylim([0 0.5])
       figure; 
       b=bar( [ a1; a2 ], 0.9 ); 
       %ylim([0 0.5]);
        ylim([-0.3 0.5]);
        set(gcf,'color','w'); grid on;
        xticklabels({'Salience - External', 'Salience - Internal' });
        ylabel('Segregation index','FontSize', 12);
    
    
    %% Figure 6B: percentage change 
    %  [ ( old value - new value ) / old value] * 100
    percent_int_1 = [(r_saldmn_1 - r_saldmn_0 )/r_saldmn_0] * 100
    percent_int_2 = [(r_saldmn_0 - r_saldmn_2 )/r_saldmn_0] * 100
    percent_int_3 = [(r_saldmn_0 - r_saldmn_3 )/r_saldmn_0] * 100
    percent_int_4 = [(r_saldmn_0 - r_saldmn_4 )/r_saldmn_0] * 100

    
    percent_ext_1 = [(r_saldansen_0 - r_saldansen_1 )/r_saldansen_0] * 100
    percent_ext_2 = [(r_saldansen_0 - r_saldansen_2 )/r_saldansen_0] * 100
    percent_ext_3 = [(r_saldansen_0 - r_saldansen_3 )/r_saldansen_0] * 100
    percent_ext_4 = [(r_saldansen_0 - r_saldansen_4 )/r_saldansen_0] * 100
    
    
    %% Figure 6C: gradients 
    for for_grad0=1
        conn_matrix=adjmat_0_noPerturb;
        idx_all0 = find(all(conn_matrix==0));
        conn_matrix(idx_all0, :)=[];
        conn_matrix(:,idx_all0)=[];

        gm_0_perturb = GradientMaps('kernel', 'na', 'approach', 'le');
        gm_0_perturb = gm_0_perturb.fit(conn_matrix,'sparsity', 90 );
        %scree_plot(gm_0_perturb.lambda{1});

        per = 1./gm_0_perturb.lambda{1};
        per = per/sum(per)*100
        %scree_plot(per);

        grad_0_perturb = zeros(100,10);
        grad_temp = cell2mat(gm_0_perturb.gradients);
        grad_0_perturb(~logical(all(conn_matrix==0)),:) = grad_temp;

        project_detection_community_boview(grad_0_perturb(:,1)*-1,'', schaefer_200_label(1:10242,1)',surfL,1);
        BoSurfStatColLim([-0.03 0.03]); colormap(flipud(hcp_colormap(:,1:3)));
    end
    for for_grad1=1
        conn_matrix=adjmat_1_Perturb
        idx_all0 = find(all(conn_matrix==0))
        conn_matrix(idx_all0, :)=[]
        conn_matrix(:,idx_all0)=[]
        %conn_matrix(:, idx_all0)=0.000000001

        gm_1_perturb = GradientMaps('kernel', 'na', 'approach', 'le');
        gm_1_perturb = gm_1_perturb.fit(conn_matrix,'sparsity', 90 );
        %scree_plot(gm_1_perturb.lambda{1});

        per = 1./gm_1_perturb.lambda{1};
        per = per/sum(per)*100
        %scree_plot(per);

        grad_1_perturb = zeros(100,10);
        grad_temp = cell2mat(gm_1_perturb.gradients);
        grad_1_perturb(~logical(all(conn_matrix==0)),:) = grad_temp;

        project_detection_community_boview(grad_1_perturb(:,1),'', schaefer_200_label(1:10242,1)',surfL,1);
        BoSurfStatColLim([-0.03 0.03]); colormap(flipud(hcp_colormap(:,1:3)));

    end
    for for_grad2=1
        conn_matrix=adjmat_2_Perturb
        idx_all0 = find(all(conn_matrix==0))
        conn_matrix(idx_all0, :)=[]
        conn_matrix(:,idx_all0)=[]

        gm_2_perturb = GradientMaps('kernel', 'na', 'approach', 'le');
        gm_2_perturb = gm_2_perturb.fit(conn_matrix,'sparsity', 90 );
        %scree_plot(gm_2_perturb.lambda{1});

        per = 1./gm_2_perturb.lambda{1};
        per = per/sum(per)*100
        %scree_plot(per);

        grad_2_perturb = zeros(100,10);
        grad_temp = cell2mat(gm_2_perturb.gradients);
        grad_2_perturb(~logical(all(conn_matrix==0)),:) = grad_temp;

        project_detection_community_boview(grad_2_perturb(:,1)*-1,'', schaefer_200_label(1:10242,1)',surfL,1);
        BoSurfStatColLim([-0.03 0.03]); colormap(flipud(hcp_colormap(:,1:3)));
    end
    for for_grad3=1
        conn_matrix=adjmat_3_Perturb
        idx_all0 = find(all(conn_matrix==0))
        conn_matrix(idx_all0, :)=[]
        conn_matrix(:,idx_all0)=[]
        %conn_matrix(:, idx_all0)=0.000000001

        gm_3_perturb = GradientMaps('kernel', 'na', 'approach', 'le');
        gm_3_perturb = gm_3_perturb.fit(conn_matrix,'sparsity', 90 );
        %scree_plot(gm_3_perturb.lambda{1});

        per = 1./gm_3_perturb.lambda{1};
        per = per/sum(per)*100
        %scree_plot(per);

        grad_3_perturb = zeros(100,10);
        grad_temp = cell2mat(gm_3_perturb.gradients);
        grad_3_perturb(~logical(all(conn_matrix==0)),:) = grad_temp;

        project_detection_community_boview(grad_3_perturb(:,1)*-1,'', schaefer_200_label(1:10242,1)',surfL,1);
        BoSurfStatColLim([-0.03 0.03]); colormap(flipud(hcp_colormap(:,1:3)));


    end
    for for_grad4=1
        conn_matrix=adjmat_4_Perturb
        idx_all0 = find(all(conn_matrix==0))
        conn_matrix(idx_all0, :)=[]
        conn_matrix(:,idx_all0)=[]
        %conn_matrix(:, idx_all0)=0.000000001

        gm_4_perturb = GradientMaps('kernel', 'na', 'approach', 'le');
        gm_4_perturb = gm_4_perturb.fit(conn_matrix,'sparsity', 90 );
        %scree_plot(gm_4_perturb.lambda{1});

        per = 1./gm_4_perturb.lambda{1};
        per = per/sum(per)*100
        %scree_plot(per);

        grad_4_perturb = zeros(100,10);
        grad_temp = cell2mat(gm_4_perturb.gradients);
        grad_4_perturb(~logical(all(conn_matrix==0)),:) = grad_temp;

        project_detection_community_boview(grad_4_perturb(:,1)*-1,'', schaefer_200_label(1:10242,1)',surfL,1);
        BoSurfStatColLim([-0.03 0.03]); colormap(flipud(hcp_colormap(:,1:3)));
    
    end

