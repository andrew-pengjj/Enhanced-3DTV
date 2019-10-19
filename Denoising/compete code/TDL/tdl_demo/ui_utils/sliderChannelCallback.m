z = get(schannel, 'Value'); 
z = max(1, min(31, round(z)));
set(schannel, 'Value', z);
set(tchannel, 'String', sprintf('Channel #%3d (slide to change)', z));

set(hi1, 'CData', msi(:, :, z));
set(hi2, 'CData', noisy_msi(:, :, z));
set(hi3, 'CData', clean_msi_bwksvd(:, :, z));
set(hi4, 'CData', clean_msi_bwbm3d(:, :, z));
set(hi5, 'CData', clean_msi_ksvd(:, :, z));
set(hi6, 'CData', clean_msi_nlm3d(:, :, z));
set(hi7, 'CData', clean_msi_bm4d(:, :, z));
set(hi8, 'CData', clean_msi_lrta(:, :, z));
set(hi9, 'CData', clean_msi_parafac(:, :, z));
set(hi10, 'CData', clean_msi_tdl(:, :, z));