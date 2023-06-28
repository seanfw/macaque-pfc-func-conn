%% clean up
close all; clear all; clc;

%% Isolate left hemisphere functional data
% !wb_command -cifti-separate /Users/seanfw/Dropbox/PRIME-DE_Oxford/data/cifti.func_pp_sm3.dtseries.nii COLUMN -metric CORTEX_LEFT /Users/seanfw/Dropbox/PRIME-DE_Oxford/data/oxford_lh_resting.func.gii
% % Isolate right hemisphere functional data
% !wb_command -cifti-separate /Users/seanfw/Dropbox/PRIME-DE_Oxford/data/cifti.func_pp_sm3.dtseries.nii COLUMN -metric CORTEX_RIGHT /Users/seanfw/Dropbox/PRIME-DE_Oxford/data/oxford_rh_resting.func.gii

%% load in the data
func_lh = gifti('/Users/seanfw/Dropbox/PRIME-DE_Oxford/data/oxford_lh_resting.func.gii');
julich_pfc_atlas = gifti('parcellations/prefrontal_2d_10k.label.gii');
julich_full_cortex_atlas = gifti('parcellations/julich_macaque_atlas_2d_filled_10k.label.gii');
lh_inflated_10k = gifti('/Users/seanfw/Dropbox/PRIME-DE_Oxford/cifti_10k_1.5mm/MacaqueYerkes19.L.inflated.10k_fs_LR.surf.gii');
Lyon_atlas_91 = gifti('parcellations/kennedy_atlas_91_10k.label.gii');
areaList_Lyon_list = load('parcellations/areaList_Donahue.mat');
lyon_areas_list = areaList_Lyon_list.areaList_Donahue

% the following is a list of the lobe that each brain region from the Lyon atlas
% belongs to.
% This is needed as we don't have the temporal lobe regions defined
% in the Julich atlas, except for area MT. In the following table, we
% purposely  mislabelled region MT as O (occipital) so that we could
% calculate the functional connectivity with all of the remaining T
% (temporal) lobe areas from the Lyon parcellation.
lyon_lobes = readtable('parcellations/kennedy_areas2lobes.xlsx')


%%
num_regions_julich = max(julich_full_cortex_atlas.cdata);
num_vertices = size(func_lh.cdata,1);
num_timepoints_all = size(func_lh.cdata,2);
num_monkeys = 19;
num_timepoints_per_monkey = num_timepoints_all/num_monkeys;

%% Calculate the principal component timecourses and functional connectivity maps for areas in the Julich atlas

mean_fc_all_regions = nan(num_vertices,num_regions_julich);
pc1_fc_all_regions = nan(num_vertices,num_regions_julich);
pc1_timeseries_all_regions = nan(num_timepoints_all,num_regions_julich);

for current_area = 1:num_regions_julich
    vertices_in_parcel = find(julich_full_cortex_atlas.cdata==current_area);
    vertices_in_parcel_timeseries = func_lh.cdata(vertices_in_parcel,:);
    vertices_in_parcel_timeseries_demeaned = nan(size(vertices_in_parcel_timeseries));
    for current_monkey = 1:num_monkeys
        % demean current monkeys data (seems like it is *almost* demeaned
        % already, so shouldn't make a big difference
        current_monkey_timepoints = (current_monkey-1)*num_timepoints_per_monkey + (1:num_timepoints_per_monkey);
        vertices_in_parcel_timeseries_demeaned(:,current_monkey_timepoints) = vertices_in_parcel_timeseries(:,current_monkey_timepoints) - repmat(mean(vertices_in_parcel_timeseries(:,current_monkey_timepoints),2),1,num_timepoints_per_monkey);

    end
    
    mean_roi_timeseries = mean(func_lh.cdata(vertices_in_parcel,:));
%     mean_fc_all_regions(:,current_area) = corr(func_lh.cdata',mean_roi_timeseries');
   
    [coeff,score,~,~,explained] = pca(vertices_in_parcel_timeseries_demeaned');
    pc1_timeseries_all_regions(:,current_area) = score(:,1);
    % Calculate functional connectivity for PC1 ROI timeseries
    pc1_fc_all_regions(:,current_area) = corr(func_lh.cdata',pc1_timeseries_all_regions(:,current_area));
    fprintf('variance explained by PC1 - region %d - %.02f \n', current_area , explained(1))


end

pc1_fc_all_regions_z = atanh(pc1_fc_all_regions);
% 
% % Save functional connectivity as gifti
% pc1_func_conn_full_cortex_surface = gifti(pc1_fc_all_regions_z);
% save(pc1_func_conn_full_cortex_surface,'pc1_func_conn_full_cortex.func.gifti','Base64Binary');

%% Identify temporal regions from the Lyon atlas

temporal_areas_idx = find(strcmp(lyon_lobes.Lobe, 'T'));

num_regions_Lyon_temporal = length(temporal_areas_idx);


%% Calculate the principal component timecourses and functional connectivity maps for each temporal lobe area of the Lyon atlas

pc1_fc_temporal_regions_Lyon = nan(num_vertices,num_regions_Lyon_temporal);
pc1_timeseries_temporal_regions_Lyon = nan(num_timepoints_all,num_regions_Lyon_temporal);

vertices_in_julich_atlas = find(julich_full_cortex_atlas.cdata);
percent_non_julich_vertices = nan(num_regions_Lyon_temporal,1);


for current_area = 1:num_regions_Lyon_temporal
    
    vertices_in_parcel = find(Lyon_atlas_91.cdata==temporal_areas_idx(current_area));
   
    non_julich_vertices_in_parcel = vertices_in_parcel(~ismember(vertices_in_parcel,vertices_in_julich_atlas));
    
    percent_non_julich_vertices(current_area) = 100*length(non_julich_vertices_in_parcel)/length(vertices_in_parcel);
    vertices_in_parcel_timeseries = func_lh.cdata(non_julich_vertices_in_parcel,:);
    vertices_in_parcel_timeseries_Lyon_demeaned = nan(size(vertices_in_parcel_timeseries));

    for current_monkey = 1:num_monkeys
        % demean current monkeys data (seems like it is *almost* demeaned
        % already, so shouldn't make a big difference
        current_monkey_timepoints = (current_monkey-1)*num_timepoints_per_monkey + (1:num_timepoints_per_monkey);
        vertices_in_parcel_timeseries_Lyon_demeaned(:,current_monkey_timepoints) = vertices_in_parcel_timeseries(:,current_monkey_timepoints) - repmat(mean(vertices_in_parcel_timeseries(:,current_monkey_timepoints),2),1,num_timepoints_per_monkey);

    end
    

    % Calculate functional connectivity for first principal component
    % timeseries of each area
    [coeff,score,~,~,explained] = pca(vertices_in_parcel_timeseries_Lyon_demeaned');
    pc1_timeseries_temporal_regions_Lyon(:,current_area) = score(:,1);
    sprintf('variance explained by PC1 - region %d - %.02f', current_area , explained(1))
    pc1_fc_all_regions_Lyon(:,current_area) = corr(func_lh.cdata',score(:,1));
end

[overlaps_sorted overlap_ind] = sort(percent_non_julich_vertices);

regions_ordered_by_overlap = lyon_areas_list(overlap_ind)
% mean_fc_all_regions_z_Lyon = atanh(mean_fc_all_regions_Lyon);
pc1_fc_all_regions_z_Lyon  = atanh(pc1_fc_all_regions_Lyon);


%% Combine Julich & Lyon (add temporal lobes)

pc1_fc_all_regions_z_Julich_Lyon = [pc1_fc_all_regions_z, pc1_fc_all_regions_z_Lyon];

% % Save functional connectivity as gifti
pc1_func_conn_Julich_Lyon_surface = gifti(pc1_fc_all_regions_z_Julich_Lyon);
% save(pc1_func_conn_full_cortex_surface,'pc1_func_conn_full_cortex.func.gifti','Base64Binary');

%% Calculate functional connectivity matrix

pc1_timeseries_julich_lyon = [pc1_timeseries_all_regions,pc1_timeseries_temporal_regions_Lyon];
[Julich_Lyon_FC,Julich_Lyon_FC_p] = corr(pc1_timeseries_julich_lyon);
Julich_Lyon_FC_z = atanh(Julich_Lyon_FC);

% %% Correct for multiple comparisons 
% 
% % identify unique correlations (remove repeated correlations)
% Julich_FC_unique = triu(Julich_FC,1)
% Julich_FC_unique_ind = find(Julich_FC_unique)
% 
% Julich_FC_unique_vec = Julich_FC(Julich_FC_unique_ind)
% Julich_FC_p_unique_vec = Julich_FC_p(Julich_FC_unique_ind)
% 
% num_tests = length(Julich_FC_p_unique_vec)
% 
% %% False discovery rate correction
% 
% % desired false discovery rate
% q = 0.001  
% % desired method - use more conservative Benjamini & Yekutieli correction
% [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(Julich_FC_p_unique_vec,q,'dep','yes');
% 
% %% Bonferroni correction
% 
% desired_p = 0.05;
% corrected_threshold = desired_p/num_tests;
% % find the significant correlations after Bonferroni correction
% Julich_FC_p_unique_vec_Bonf_sig = Julich_FC_p_unique_vec<corrected_threshold
% 
% % identify the minimum r value
% Julich_FC_unique_vec(Julich_FC_p_unique_vec_Bonf_sig)

%% get julich area names
load julich_areas_list.mat
lyon_areas_list_clean = strrep(lyon_areas_list,'/','_')

julich_lyon_areas_list = [julich_areas_list,lyon_areas_list_clean(temporal_areas_idx)'];
%% Create revised surface plots with colorbrewer2 colormap 

% Colours from Colorbrewer2 - with extra yellow to highlight region of
% interest
eleven_class_RdBu_Ye = [103,0,31;178,24,43;214,96,77;244,165,130;253,219,199;247,247,247;209,229,240;146,197,222;67,147,195;33,102,172;0,0,0;]./255

nine_class_Reds = [0,0,0,;255,245,240;254,224,210;252,187,161;252,146,114;251,106,74;239,59,44;203,24,29;165,15,21;103,0,13]./255;

customcmap = flipud(eleven_class_RdBu_Ye);

for current_region = 1:num_regions_julich
    pc1_func_conn_full_cortex_temp = pc1_func_conn_Julich_Lyon_surface;
    vertices_in_parcel = find(julich_full_cortex_atlas.cdata==current_region);
    % show seed region in black
    pc1_func_conn_full_cortex_temp.cdata(vertices_in_parcel,current_region) = -1;
    % show midline in white
    pc1_func_conn_full_cortex_temp.cdata(isnan(pc1_func_conn_full_cortex_temp.cdata))=0;
    
    myfig = figure('units','normalized','outerposition',[0.5 0.4 0.5 0.6]);
    set(gcf,'color','w');
    colormap(customcmap)
    mysurf = plot(lh_inflated_10k,pc1_func_conn_full_cortex_temp,current_region);
    direction1 = [0 1 0];
    direction2 = [0 0 1];
    direction3 = [1 0 0];
    rotate(mysurf,direction1,60);
    rotate(mysurf,direction2,90);
    rotate(mysurf,direction3,15);
    material([0.5,0.5,0.15])
    h = colorbar;
    set(h,'location','southoutside')
    h.FontSize = 20;
    caxis([-0.25,0.25])
%     ylabel(h, 'functional connectivity (z)')
    saveas(myfig,sprintf('functional_conn_images/func_conn_%s_lateral.png',julich_areas_list{current_region}));

    rotate(mysurf,direction1,180);
    rotate(mysurf,direction3,30);
    light('Position',[0 1 0])



    saveas(myfig,sprintf('functional_conn_images/func_conn_%s_medial.png',julich_areas_list{current_region}));
    close all;

end


for current_region = 1:num_regions_Lyon_temporal
    pc1_func_conn_full_cortex_temp = pc1_func_conn_Julich_Lyon_surface;
    vertices_in_parcel = find(Lyon_atlas_91.cdata==temporal_areas_idx(current_region));
    % show seed region in black
    pc1_func_conn_full_cortex_temp.cdata(vertices_in_parcel,num_regions_julich+current_region) = -1;
    % show midline in white
    pc1_func_conn_full_cortex_temp.cdata(isnan(pc1_func_conn_full_cortex_temp.cdata))=0;
    
    myfig = figure('units','normalized','outerposition',[0.5 0.4 0.5 0.6]);
    set(gcf,'color','w');
    colormap(customcmap)
    mysurf = plot(lh_inflated_10k,pc1_func_conn_full_cortex_temp,num_regions_julich+current_region);
    direction1 = [0 1 0];
    direction2 = [0 0 1];
    direction3 = [1 0 0];
    rotate(mysurf,direction1,60);
    rotate(mysurf,direction2,90);
    rotate(mysurf,direction3,15);
    material([0.5,0.5,0.15])
    h = colorbar;
    set(h,'location','southoutside')
    h.FontSize = 20;
    caxis([-0.25,0.25])
%     ylabel(h, 'functional connectivity (z)')
    saveas(myfig,sprintf('functional_conn_images/func_conn_%s_lateral.png',lyon_areas_list_clean{temporal_areas_idx(current_region)}));

    rotate(mysurf,direction1,180);
    rotate(mysurf,direction3,30);
    light('Position',[0 1 0])



    saveas(myfig,sprintf('functional_conn_images/func_conn_%s_medial.png',lyon_areas_list_clean{temporal_areas_idx(current_region)}));
    close all;

end

%% Save as Excel for Lucija

% Set within-area connectivity Z=1.1 (just past the maximum of the data)
Julich_Lyon_FC_z(Julich_Lyon_FC_z==Inf)=1.1;
filename = 'func_conn_data/connectivity_data_whole_cortex_Julich_Lyon_z.xlsx';
writematrix(Julich_Lyon_FC_z,filename,'Sheet',1,'Range','B2:EI139')
writecell(julich_lyon_areas_list',filename,'Sheet',1,'Range','A2');
writecell(julich_lyon_areas_list,filename,'Sheet',1,'Range','B1');

save func_conn_data/Julich_Lyon_FC_138x138.mat Julich_Lyon_FC_z julich_lyon_areas_list

%% Now create individual monkey/area func conn maps
% 

% 
% pc1_timeseries_indiv_monkeys = nan(num_timepoints_per_monkey,num_regions_julich_lyon,num_monkeys);
% 
% vertices_in_atlas = find(julich_full_cortex_atlas.cdata);
% 
% for current_monkey = 1:num_monkeys
%     
%     current_monkey_timepoints = (1:num_timepoints_per_monkey) + num_timepoints_per_monkey*(current_monkey-1);
%     
%     % 112 Julich areas
%     for current_area = 1:num_regions_julich
%         sprintf('monkey %d - area %d, %s',current_monkey,current_area, julich_areas_list{current_area})
%         vertices_in_parcel = find(julich_full_cortex_atlas.cdata==current_area);
%         vertices_in_parcel_timeseries = func_lh.cdata(vertices_in_parcel,current_monkey_timepoints);
% 
%        [coeff,score,~,~,explained] = pca(vertices_in_parcel_timeseries');
%        pc1_timeseries_indiv_monkeys(:,current_area,current_monkey) = score(:,1);
% 
%     end
%     
%     Julich_Lyon_FC_indiv(:,:,current_monkey) = corr(pc1_timeseries_indiv_monkeys(:,:,current_monkey));
%     
% end
% 
% 
% 
% Julich_FC_indiv_z = atanh(Julich_Lyon_FC_indiv);
% % 
% save func_conn_data/Julich_FC_indiv_z_112x112.mat Julich_FC_indiv_z julich_areas_list
% save func_conn_data/Julich_FC_indiv_r_112x112.mat Julich_FC_indiv julich_areas_list


%% Now create individual monkey/area func conn maps for 
num_regions_julich_lyon = num_regions_julich + num_regions_Lyon_temporal;

Julich_Lyon_FC_indiv = nan(num_regions_julich_lyon,num_regions_julich_lyon,num_monkeys);

pc1_timeseries_indiv_monkeys = nan(num_timepoints_per_monkey,num_regions_julich_lyon,num_monkeys);

% vertices_in_atlas = find(julich_full_cortex_atlas.cdata);

for current_monkey = 1:num_monkeys
    
    current_monkey_timepoints = (1:num_timepoints_per_monkey) + num_timepoints_per_monkey*(current_monkey-1);
    
    % 112 Julich areas
    for current_area = 1:num_regions_julich
        sprintf('monkey %d - area %d, %s',current_monkey,current_area, julich_areas_list{current_area})
        vertices_in_parcel = find(julich_full_cortex_atlas.cdata==current_area);
        vertices_in_parcel_timeseries = func_lh.cdata(vertices_in_parcel,current_monkey_timepoints);

       [coeff,score,~,~,explained] = pca(vertices_in_parcel_timeseries');
       pc1_timeseries_indiv_monkeys(:,current_area,current_monkey) = score(:,1);

    end
    
    % 26 Lyon temporal lobe areas
    for current_area = 1:num_regions_Lyon_temporal
        sprintf('monkey %d - area %d, %s',current_monkey,num_regions_julich + current_area, lyon_areas_list_clean{temporal_areas_idx(current_area)})
        vertices_in_parcel = find(Lyon_atlas_91.cdata==temporal_areas_idx(current_area));
        vertices_in_parcel_timeseries = func_lh.cdata(vertices_in_parcel,current_monkey_timepoints);

       [coeff,score,~,~,explained] = pca(vertices_in_parcel_timeseries');
       pc1_timeseries_indiv_monkeys(:,num_regions_julich+current_area,current_monkey) = score(:,1);

    end
    
    Julich_Lyon_FC_indiv(:,:,current_monkey) = corr(pc1_timeseries_indiv_monkeys(:,:,current_monkey));
    
end



Julich_Lyon_FC_indiv_z = atanh(Julich_Lyon_FC_indiv);
% 


%% Save as Excel

% Set within-area connectivity Z=2 (just past the maximum of the data)
Julich_Lyon_FC_indiv_z(Julich_Lyon_FC_indiv_z==Inf)=2;
filename = 'func_conn_data/connectivity_data_whole_cortex_Julich_Lyon_indiv_z.xlsx';
for current_monkey = 1:num_monkeys
    writematrix(Julich_Lyon_FC_indiv_z(:,:,current_monkey),filename,'Sheet',current_monkey,'Range','B2:EI139')
    writecell(julich_lyon_areas_list',filename,'Sheet',current_monkey,'Range','A2');
    writecell(julich_lyon_areas_list,filename,'Sheet',current_monkey,'Range','B1');
end

% Save as Excel

% Set within-area connectivity r=1 (just past the maximum of the data)
Julich_Lyon_FC_indiv(Julich_Lyon_FC_indiv==Inf)=1;

filename = 'func_conn_data/connectivity_data_whole_cortex_Julich_Lyon_indiv_r.xlsx';
for current_monkey = 1:num_monkeys
    writematrix(Julich_Lyon_FC_indiv(:,:,current_monkey),filename,'Sheet',current_monkey,'Range','B2:EI139')
    writecell(julich_lyon_areas_list',filename,'Sheet',current_monkey,'Range','A2');
    writecell(julich_lyon_areas_list,filename,'Sheet',current_monkey,'Range','B1');
end
% 
