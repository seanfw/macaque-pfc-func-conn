close all; clear all; clc;
%% merge Julich borders
receptor_gradients_dir = '/Users/seanfw/Dropbox/SFW_Wang_projects/Active_projects/Receptor_gradients/walsh-receptor-gradients';

% Find the border files for newly-drawn borders
prefrontal_borders_list = dir('prefrontal_borders_32k/area*');
motor_borders_list = dir('motor_borders_2022_32k/area*');

%  Find the border files for  other areas
parietal_borders_list = dir(fullfile([receptor_gradients_dir '/Julich_borders/Parietal_borders/area*']));
occipital_borders_list = dir(fullfile([receptor_gradients_dir '/Julich_borders/Occipital_borders/area*']));
cingulate_borders_list = dir(fullfile([receptor_gradients_dir '/Julich_borders/Cingulate_borders/area*']));
temporal_borders_list = dir(fullfile([receptor_gradients_dir '/Julich_borders/Temporal_borders/area*']));


prefrontal_borders_list_cell = struct2cell(prefrontal_borders_list);
prefrontal_borders_list_cell = prefrontal_borders_list_cell(1,:);
motor_borders_list_cell = struct2cell(motor_borders_list);
motor_borders_list_cell = motor_borders_list_cell(1,:);
parietal_borders_list_cell = struct2cell(parietal_borders_list);
parietal_borders_list_cell = parietal_borders_list_cell(1,:);
occipital_borders_list_cell = struct2cell(occipital_borders_list);
occipital_borders_list_cell = occipital_borders_list_cell(1,:);
cingulate_borders_list_cell = struct2cell(cingulate_borders_list);
cingulate_borders_list_cell = cingulate_borders_list_cell(1,:);
temporal_borders_list_cell = struct2cell(temporal_borders_list);
temporal_borders_list_cell = temporal_borders_list_cell(1,:);

julich_borders_list = [prefrontal_borders_list_cell, motor_borders_list_cell,parietal_borders_list_cell, occipital_borders_list_cell cingulate_borders_list_cell temporal_borders_list_cell];
julich_areas_list = strrep(julich_borders_list,'.border','');

save julich_areas_list.mat julich_areas_list


julich_pfc_areas_list = strrep(prefrontal_borders_list_cell,'.border','');

save julich_pfc_areas_list.mat julich_pfc_areas_list



% How many border files do we have for each lobe?
num_prefrontal_areas = length(prefrontal_borders_list);
num_motor_areas = length(motor_borders_list);
num_parietal_areas = length(parietal_borders_list);
num_occipital_areas = length(occipital_borders_list);
num_cingulate_areas = length(cingulate_borders_list);
num_temporal_areas = length(temporal_borders_list);

% where are the border files stored?
path_to_prefrontal_borders = 'prefrontal_borders_32k/';
path_to_motor_borders = 'motor_borders_2022_32k/';
% where are the border files stored?
path_to_parietal_borders = fullfile([receptor_gradients_dir '/Julich_borders/Parietal_borders/']);
path_to_occipital_borders = fullfile([receptor_gradients_dir '/Julich_borders/Occipital_borders/']);
path_to_cingulate_borders = fullfile([receptor_gradients_dir '/Julich_borders/Cingulate_borders/']);
path_to_temporal_borders = fullfile([receptor_gradients_dir '/Julich_borders/Temporal_borders/']);

%% Let's merge together all the areas
% Let's put together the prefrontal borders (drawn by Lucija)
% This just creates the prefrontal.border file by putting together the first
% two regions
merge_prefrontal = sprintf('wb_command -border-merge prefrontal_borders_32k/prefrontal.border -border  %s -border  %s ', fullfile([path_to_prefrontal_borders prefrontal_borders_list(1).name]), fullfile([path_to_prefrontal_borders prefrontal_borders_list(2).name]));
[u v]=system(merge_prefrontal);

% This adds on each of the other regions to the prefrontal.border file one at
% a time. It's not the most efficient for the computer, but it's still quite
% quick and it's efficient for the coder (me). 
for current_border = 3:num_prefrontal_areas
    fprintf('adding %s \n',prefrontal_borders_list(current_border).name)
    merge_prefrontal = sprintf('wb_command -border-merge prefrontal_borders_32k/prefrontal.border -border  prefrontal_borders_32k/prefrontal.border -border  %s ', fullfile([path_to_prefrontal_borders prefrontal_borders_list(current_border).name]));
    [u v]=system(merge_prefrontal);

end

% Let's put together the motor borders (drawn by Lucija)
% This just creates the motor.border file by putting together the first
% two regions
merge_motor = sprintf('wb_command -border-merge motor_borders_2022_32k/motor.border -border  %s -border  %s ', fullfile([path_to_motor_borders motor_borders_list(1).name]), fullfile([path_to_motor_borders motor_borders_list(2).name]));
[u v]=system(merge_motor);

% This adds on each of the other regions to the motor.border file one at
% a time. It's not the most efficient for the computer, but it's still quite
% quick and it's efficient for the coder (me). 
for current_border = 3:num_motor_areas
    fprintf('adding %s \n',motor_borders_list(current_border).name)
    merge_motor = sprintf('wb_command -border-merge motor_borders_2022_32k/motor.border -border  motor_borders_2022_32k/motor.border -border  %s ', fullfile([path_to_motor_borders motor_borders_list(current_border).name]));
    [u v]=system(merge_motor);

end

% Now let's do the same for the parietal regions. 
% Let's put together the parietal borders (drawn by Meiqi)
% This just creates the parietal.border file by putting together the first
% two regions
merge_parietal = sprintf('wb_command -border-merge other_borders/parietal.border -border  %s -border  %s ', fullfile([path_to_parietal_borders parietal_borders_list(1).name]), fullfile([path_to_parietal_borders parietal_borders_list(2).name]));
[u v]=system(merge_parietal);

% This adds on each of the other regions to the parietal.border file one at
% a time. It's not the most efficient for the computer, but it's still quite
% quick and it's efficient for the coder (me). 
for current_border = 3:num_parietal_areas;
    fprintf('adding %s \n',parietal_borders_list(current_border).name)
    merge_parietal = sprintf('wb_command -border-merge other_borders/parietal.border -border  other_borders/parietal.border -border  %s ', fullfile([path_to_parietal_borders parietal_borders_list(current_border).name]));
    [u v]=system(merge_parietal);

end

% Now let's put the occipital regions together
% Let's put together the occipital borders (drawn by Lucija and Meiqi)
% This just creates the occipital.border file by putting together the first
% two regions
merge_occipital = sprintf('wb_command -border-merge other_borders/occipital.border -border  %s -border  %s ', fullfile([path_to_occipital_borders occipital_borders_list(1).name]), fullfile([path_to_occipital_borders occipital_borders_list(2).name]));
[u v]=system(merge_occipital);

% This adds on each of the other regions to the occipital.border file one at
% a time. It's not the most efficient for the computer, but it's still quite
% quick and it's efficient for the coder (me). 
for current_border = 3:num_occipital_areas;
    fprintf('adding %s \n',occipital_borders_list(current_border).name)
    merge_occipital = sprintf('wb_command -border-merge other_borders/occipital.border -border  other_borders/occipital.border -border  %s ', fullfile([path_to_occipital_borders occipital_borders_list(current_border).name]));
    [u v]=system(merge_occipital);

end


% Now let's put the cingulate regions together
% Let's put together the cingulate borders (drawn by Nicola)
% This just creates the cingulate.border file by putting together the first
% two regions
merge_cingulate = sprintf('wb_command -border-merge other_borders/cingulate.border -border  %s -border  %s ', fullfile([path_to_cingulate_borders cingulate_borders_list(1).name]), fullfile([path_to_cingulate_borders cingulate_borders_list(2).name]));
[u v]=system(merge_cingulate);

% This adds on each of the other regions to the cingulate.border file one at
% a time. It's not the most efficient for the computer, but it's still quite
% quick and it's efficient for the coder (me). 
for current_border = 3:num_cingulate_areas;
    fprintf('adding %s \n',cingulate_borders_list(current_border).name)
    merge_cingulate = sprintf('wb_command -border-merge other_borders/cingulate.border -border  other_borders/cingulate.border -border  %s ', fullfile([path_to_cingulate_borders cingulate_borders_list(current_border).name]));
    [u v]=system(merge_cingulate);

end

!echo "merging frontal, parietal, occipital, cingulate and MT borders into a 3d file"
!wb_command -border-merge other_borders/frontoparietooccipital.border -border prefrontal_borders_32k/prefrontal.border -border motor_borders_2022_32k/motor.border -border other_borders/parietal.border -border other_borders/occipital.border -border other_borders/cingulate.border -border $RECEPTOR_GRADIENTS_DIR/Julich_borders/Temporal_borders/area_MT.border


%% change Julich borders to ROIs
!echo "filling in borders to create ROIs"
!wb_command -border-to-rois $YERKES_DATA_DIR/MNINonLinear/fsaverage_LR32k/MacaqueYerkes19_v1.2.L.inflated.32k_fs_LR.surf.gii other_borders/frontoparietooccipital.border other_borders/frontoparietooccipital_3d.func.gii

% Let's change the mask to a ROI (so that we only fill in gaps in parts of
% the cortex that we've covered
% Note that this creates a binary mask with 1s in the unfilled regions, and
% zeros in the frontoparietooccipital cortex. Presumably this is because of
% the size of the ROI.
!echo "creating a ROI from the mask border file "
!wb_command -border-to-rois $YERKES_DATA_DIR/MNINonLinear/fsaverage_LR32k/MacaqueYerkes19_v1.2.L.inflated.32k_fs_LR.surf.gii $RECEPTOR_GRADIENTS_DIR/Julich_borders/frontoparietooccipital_mask.border other_borders/non_frontoparietooccipital_mask.func.gii

%% load in 3d ROI file
full_cortex_ROIs = gifti('other_borders/frontoparietooccipital_3d.func.gii');

%% identify and remove overlaps
overlapping_vertices = find(sum(full_cortex_ROIs.cdata,2)>1);
num_overlapping_vertices = length(overlapping_vertices);

for current_vertex = 1:num_overlapping_vertices
    
    overlapping_ROIs = find(full_cortex_ROIs.cdata(overlapping_vertices(current_vertex),:));
    num_overlapping_ROIs = length(overlapping_ROIs);
    max_overlapping_ROIs = max(overlapping_ROIs);
    full_cortex_ROIs.cdata(overlapping_vertices(current_vertex),overlapping_ROIs) = 0;
    full_cortex_ROIs.cdata(overlapping_vertices(current_vertex),max_overlapping_ROIs) = 1;

end

%% merge 3d file into a single 2d surface

num_ROIs = size(full_cortex_ROIs.cdata,2);
num_vertices = size(full_cortex_ROIs.cdata,1);
ROI_idx = 1:num_ROIs;
ROI_idx_mat = repmat(ROI_idx,num_vertices,1);
full_cortex_ROIs_map_3d = full_cortex_ROIs.cdata.*ROI_idx_mat;
full_cortex_ROIs_map_2d = sum(full_cortex_ROIs_map_3d,2);

%% write new gifti file

temp = gifti(fullfile([receptor_gradients_dir '/data/macaque_myelin.func.gii']));
full_cortex_ROIs_2d = temp;
full_cortex_ROIs_2d.cdata = full_cortex_ROIs_map_2d;

save(full_cortex_ROIs_2d,'other_borders/frontoparietooccipital_2d.func.gifti','Base64Binary');
%% plot 
% load in inflated left hemisphere surface
l_inflated = gifti(fullfile([receptor_gradients_dir '/MacaqueYerkes19_v1.2_976nz/MNINonLinear/fsaverage_LR32k/MacaqueYerkes19_v1.2.L.inflated.32k_fs_LR.surf.gii']));

figure
plot(l_inflated,full_cortex_ROIs_2d);


%% Fill in the gaps in the Julich atlas 

% Identify vertices that have a value 0 but lie inside the
% frontoparietoccipital mask.
% This should identify border vertices between the brain sections.

% Find all the vertices without a parcel number - this should be the
% unfilled border areas and unsampled cortex.
border_plus_unsampled = find(full_cortex_ROIs_2d.cdata == 0);

% Let's load in the frontoparietooccipital mask
mask = gifti('other_borders/non_frontoparietooccipital_mask.func.gii');
mask_ind = find(mask.cdata);

% Find all the vertices without a parcel number excluding the unsampled cortex.
julich_border_ind = setdiff(border_plus_unsampled,mask_ind);


julich_borders = temp;
julich_borders.cdata = zeros(num_vertices,1);
julich_borders.cdata(julich_border_ind,1) = 1;


% Visualise borders
close all;
figure
plot(l_inflated,julich_borders);

figure
plot(l_inflated,mask);
% save(julich_borders,'../Julich_borders/julich_borders.func.gii','Base64Binary');
%%

% Let's find the distances for each border vertex to the other vertices

l_midthickness = gifti(fullfile([receptor_gradients_dir '/MacaqueYerkes19_v1.2_976nz/MNINonLinear/fsaverage_LR32k/MacaqueYerkes19_v1.2.L.midthickness.32k_fs_LR.surf.gii']));

num_border_vertices = length(julich_border_ind);

full_cortex_ROIs_2d_filled = full_cortex_ROIs_2d;

for current_vertex = 1:num_border_vertices
   sprintf('Filling in border vertices - %f percent complete', 100*round(current_vertex/num_border_vertices,2))
   % assign each vertex to its closest sample    
   command = [' wb_command -surface-geodesic-distance $YERKES_DATA_DIR/MNINonLinear/fsaverage_LR32k/MacaqueYerkes19_v1.2.L.midthickness.32k_fs_LR.surf.gii ' sprintf('%d',julich_border_ind(current_vertex)-1) ' geo_dist_current_vertex.func.gii'];
   system(command);
   geo_dist_current_vertex = gifti('geo_dist_current_vertex.func.gii');

   [sorted_dist,dist_ind] = sort(geo_dist_current_vertex.cdata);

   % find the closest vertex with a non-zero parcel value
   closest_nonzero_vertex = dist_ind(min(find(full_cortex_ROIs_2d.cdata(dist_ind))));
   % now let's add this border vertex to the closest parcel
   full_cortex_ROIs_2d_filled.cdata(julich_border_ind(current_vertex),1) = full_cortex_ROIs_2d.cdata(closest_nonzero_vertex,1);

end

% save 
save(full_cortex_ROIs_2d_filled,'other_borders/julich_sections_2d_filled.func.gii','Base64Binary');

%% Prepare the label text file

full_cortex_colours = [round(255*rand(num_ROIs,3)),255*ones(num_ROIs,1)];

fileID = fopen('parcellations/full_cortex_label_table.txt','w');

for current_area = 1:num_ROIs
    % create the table with the label name information
    current_label_name = julich_areas_list{current_area};
    current_label_colour = full_cortex_colours(current_area,:);

    fprintf(fileID,'%s \n',sprintf('%s',julich_areas_list{current_area}));
    fprintf(fileID,'%d %d %d %d %d\n',[current_area current_label_colour]);
    
end

fclose(fileID);

%% create the label surface file 
create_label_file = sprintf('wb_command -metric-label-import other_borders/julich_sections_2d_filled.func.gii parcellations/full_cortex_label_table.txt parcellations/julich_macaque_atlas_2d_filled.label.gii');
[u v]=system(create_label_file);


