function visualizeXsect( volumes, fig_title, crop )
    if ~exist('crop','var')
        crop = 1;
    end
    % load constants
    constantsSparseTraj3D
    
    for i=1:length(volumes)
        y = volumes{i};
        if isempty(y)
            continue;
        end
        x = round(size(y)./2);

        p1 = squeeze(y(x(1),:,:));
        p2 = squeeze(y(:,x(2),:))';
        p3 = squeeze(y(:,:,x(3)));

        if crop
            offset = 1;
            p1(:,1:x(3)-offset) = [];
            p1(1:x(2)-offset,:) = [];

            p3(:,1:x(2)-offset) = [];
            p3(1:x(1)-offset,:) = [];

            p2(:,1:x(1)-offset) = [];
            p2(1:x(3)-offset,:) = [];

            x(:) = offset;
        end

        x1 = zeros(size(p1)) + x(1);
        x2 = zeros(size(p2)) + x(2);
        x3 = zeros(size(p3)) + x(3);
        [X1_1,X1_2] = meshgrid(1:size(p1,2),1:size(p1,1));
        [X2_1,X2_2] = meshgrid(1:size(p2,2),1:size(p2,1));
        [X3_1,X3_2] = meshgrid(1:size(p3,2),1:size(p3,1));

        h = figure;
        title(fig_title{i});
        %set(h, 'color', 'white');
        set(h, 'renderer', 'zbuffer');
        hold all
        warp(X1_2,x1,X1_1,p1)
        warp(x2,X2_1,X2_2,p2)
        warp(X3_1,X3_2,x3,p3)
        view([35 30])
        axis off vis3d;
        camproj('persp')
    end

end



