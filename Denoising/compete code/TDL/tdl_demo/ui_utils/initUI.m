%% Create UI control objects
hfig = figure;
set(hfig, 'unit', 'normalized', 'position', [0, 0, 1, 1], ...
    'name', 'Tensor Dictionary Learning demo', 'NumberTitle', 'off', ...
    'MenuBar', 'none', 'defaultuicontrolunits', 'normalized');
% Slider for changing channel
tchannel = uicontrol(hfig, 'style', 'text', ...
    'Position', [0.55, 0.8, 0.35, 0.05], ...
    'Horizontal', 'Left', ...
    'String', 'Channel #  1 (slide to change)', ...
    'BackgroundColor', get(hfig, 'color'), ...
    'Fontsize', 12);
schannel = uicontrol(hfig, 'style', 'slider', ...
    'position', [0.55, 0.75, 0.35, 0.05],...
    'max', nbands, 'min', 1, ...
    'sliderstep', [1/(nbands-1), 1/(nbands-1)], ...
    'value', 1, ...
    'callback', 'sliderChannelCallback;');
% Toggle button for zooming
bzoom = uicontrol(hfig, 'style', 'togglebutton', ...
    'position', [0.8, 0.81, 0.05, 0.05],...
    'String', 'Zoom', ...
    'Callback', 'toggleZoomCallback;');
% Toggle button for panning
bpan = uicontrol(hfig, 'style', 'togglebutton', ...
    'position', [0.85, 0.81, 0.05, 0.05],...
    'String', 'Pan', ...
    'Callback', 'togglePanCallback;');

%% Draw the first channel
ha1 = subplot(3, 4, 1);
hi1 = imshow(msi(:, :, 1));
title('Noise-free');

ha2 = subplot(3, 4, 2);
hi2 = imshow(noisy_msi(:, :, 1));
title('Noisy');

ha3 = subplot(3, 4, 5);
hi3 = imshow(clean_msi_bwksvd(:, :, 1));
title('Band-wise KSVD');

ha4 = subplot(3, 4, 6);
hi4 = imshow(clean_msi_bwbm3d(:, :, 1));
title('Band-wise BM3D');

ha5 = subplot(3, 4, 7);
hi5 = imshow(clean_msi_ksvd(:, :, 1));
title('Integral KSVD');

ha6 = subplot(3, 4, 8);
hi6 = imshow(clean_msi_nlm3d(:, :, 1));
title('3D NLM');

ha7 = subplot(3, 4, 9);
hi7 = imshow(clean_msi_bm4d(:, :, 1));
title('BM4D');

ha8 = subplot(3, 4, 10);
hi8 = imshow(clean_msi_lrta(:, :, 1));
title('LRTA');

ha9 = subplot(3, 4, 11);
hi9 = imshow(clean_msi_parafac(:, :, 1));
title('PARAFAC');

ha10 = subplot(3, 4, 12);
hi10 = imshow(clean_msi_tdl(:, :, 1));
title('Tensor DL');

linkaxes([ha1, ha2, ha3, ha4, ha5, ha6, ...
    ha7, ha8, ha9, ha10]);

hzoom = zoom;
hpan = pan;