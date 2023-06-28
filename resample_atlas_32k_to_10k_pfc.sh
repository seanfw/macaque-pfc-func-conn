# echo "downsampling Lucija's atlas from 32k to 10k"
# # Downsample Lucija's map from 32k to 10k to match Ting's connectivity data
# wb_command -label-resample parcellations/prefrontal_2d.label.gii $RECEPTOR_GRADIENTS_DIR/CMI_gradients/Macaque_32k/MacaqueYerkes19.L.sphere.32k_fs_LR.surf.gii $RECEPTOR_GRADIENTS_DIR/CMI_gradients/Macaque_10k/MacaqueYerkes19.L.sphere.10k_fs_LR.surf.gii ADAP_BARY_AREA  parcellations/prefrontal_2d_10k.label.gii -area-surfs $RECEPTOR_GRADIENTS_DIR/CMI_gradients/Macaque_32k/MacaqueYerkes19.L.midthickness.32k_fs_LR.surf.gii $RECEPTOR_GRADIENTS_DIR/CMI_gradients/Macaque_10k/MacaqueYerkes19.L.midthickness.10k_fs_LR.surf.gii

# echo "downsampling Lucija's borders from 32k to 10k"
# # Downsample Lucija's map from 32k to 10k to match Ting's connectivity data
# borders_list=$(ls prefrontal_borders_32k/*.border)
# 
# # echo $borders_list
# 
# 
# for current_border in $borders_list
# do
#     current_border_short=$(echo $current_border | cut -d'/' -f 2)
#     echo "downsampling $current_border_short"
# 
#      wb_command -border-resample $current_border $RECEPTOR_GRADIENTS_DIR/CMI_gradients/Macaque_32k/MacaqueYerkes19.L.sphere.32k_fs_LR.surf.gii $RECEPTOR_GRADIENTS_DIR/CMI_gradients/Macaque_10k/MacaqueYerkes19.L.sphere.10k_fs_LR.surf.gii prefrontal_borders_10k/$current_border_short
#     
#     
# done


echo "downsampling Julich fronto-parieto-occipital atlas from 32k to 10k"
# Downsample Julich map from 32k to 10k to match Ting's connectivity data
wb_command -label-resample parcellations/julich_macaque_atlas_2d_filled.label.gii $RECEPTOR_GRADIENTS_DIR/CMI_gradients/Macaque_32k/MacaqueYerkes19.L.sphere.32k_fs_LR.surf.gii $RECEPTOR_GRADIENTS_DIR/CMI_gradients/Macaque_10k/MacaqueYerkes19.L.sphere.10k_fs_LR.surf.gii ADAP_BARY_AREA  parcellations/julich_macaque_atlas_2d_filled_10k.label.gii -area-surfs $RECEPTOR_GRADIENTS_DIR/CMI_gradients/Macaque_32k/MacaqueYerkes19.L.midthickness.32k_fs_LR.surf.gii $RECEPTOR_GRADIENTS_DIR/CMI_gradients/Macaque_10k/MacaqueYerkes19.L.midthickness.10k_fs_LR.surf.gii
