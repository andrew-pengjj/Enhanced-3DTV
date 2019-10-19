function [ S ] = sampling( trajectory, n, rot_deg, line_nbr, line_std, coverage, tol )

    % load constants
    constantsSparseTraj3D

    phantom_size = [n n n];
    S = ones(phantom_size);
    
    if trajectory==RADIAL || trajectory==LIM_ANGLE
        history = zeros(2);
        stop = 0;
        % getting line-mask (radial lines)
        S = zeros(phantom_size);
        scale = n*0.2;
        % scaling parameter to the desired coverage
        while (sum(sum(S(:,:,1)))==0 || abs(sum(sum(S(:,:,1)))/numel(S(:,:,1))-coverage)>tol) && ~stop
            if sum(sum(S(:,:,1)))/numel(S(:,:,1))>coverage
                scale = scale-tol/2;
            else
                scale = scale+tol/2;
            end
            if scale==history(1,1)
                scale = history(1,history(2,:)==min(history(2,:)));
                stop = 1;
            end
            S(:,:,1) = ifftshift( radial( trajectory, n, scale, 1, rot_deg, line_nbr, line_std ) );
            history = circshift(history,[0 1]);
            history(1,end) = scale;
            history(2,end) = abs(sum(sum(S(:,:,1)))/numel(S(:,:,1))-coverage);
        end
        for i2=2:n
            S(:,:,i2) = ifftshift( radial( trajectory, n, scale, i2, rot_deg, line_nbr, line_std ) );
        end
    end

    if trajectory==SPIRAL
        history = zeros(2);
        stop = 0;
        % getting line-mask (spiral lines)
        S = zeros(phantom_size);
        scale = n*6e-3;
        t = linspace(0,10*n,5e5);
        % scaling parameter to the desired coverage
        while (sum(sum(S(:,:,1)))==0 || abs(sum(sum(S(:,:,1)))/numel(S(:,:,1))-coverage)>tol) && ~stop
            if sum(sum(S(:,:,1)))/numel(S(:,:,1))>coverage
                scale = scale+1e-2;
            else
                scale = scale-1e-2;
            end
            if scale==history(1,1)
                scale = history(1,history(2,:)==min(history(2,:)));
                stop = 1;
            end
            S(:,:,1) = ifftshift( spiral( t, n, scale, 1, rot_deg, line_nbr, line_std ));
            history = circshift(history,[0 1]);
            history(1,end) = scale;
            history(2,end) = abs(sum(sum(S(:,:,1)))/numel(S(:,:,1))-coverage);
        end
        for i2=2:n
            S(:,:,i2) = ifftshift(spiral( t, n, scale, i2, rot_deg, line_nbr, line_std ));
        end
    end

    if trajectory==LOG_SPIRAL
        history = zeros(2);
        stop = 0;
        % getting line-mask (log-spiral lines)
        S = zeros(phantom_size);
        scale = 3/n;
        t = linspace(0,10*n,5e5);
        % scaling parameter to the desired coverage
        while (sum(sum(S(:,:,1)))==0 || abs(sum(sum(S(:,:,1)))/numel(S(:,:,1))-coverage)>tol) && ~stop
            if sum(sum(S(:,:,1)))/numel(S(:,:,1))>coverage
                scale = scale+1e-3;
            else
                scale = scale-1e-3;
            end
            if scale==history(1,1)
                scale = history(1,history(2,:)==min(history(2,:)));
                stop = 1;
            end
            S(:,:,1) = ifftshift( log_spiral( t, n, scale, 1, rot_deg, line_nbr, line_std ));
            history = circshift(history,[0 1]);
            history(1,end) = scale;
            history(2,end) = abs(sum(sum(S(:,:,1)))/numel(S(:,:,1))-coverage);
        end
        for i2=2:n
            S(:,:,i2) = ifftshift( log_spiral( t, n, scale, i2, rot_deg, line_nbr, line_std ) );
        end
    end
    
    if trajectory==HELICAL
        scale = 1;
        t = 0:n*1000;
        r = n/2;% 0.01*t;
        x = r.*cos(scale*2*pi*t/10000);
        y = r.*sin(scale*2*pi*t/10000);
        z = t/1000;
        
        temp = zeros([n n n]);
        S = zeros([n n n]);
        for i=1:length(t)
            xr = max(1,min(n,round(x(i)+n/2)));
            yr = max(1,min(n,round(y(i)+n/2)));
            zr = max(1,min(n,round(z(i))));
            
            if temp(yr,xr,zr)==0
                temp(yr,xr,zr) = 1;
                xr = max(1,min(n,round(xr+line_std*randn)));
                yr = max(1,min(n,round(yr+line_std*randn)));
                zr = max(1,min(n,round(zr+line_std*randn)));
                S(yr,xr,zr) = 1;
            end
        end
        for i=1:size(S,3)
            S(:,:,i) = ifftshift(S(:,:,i));
        end
        for i=1:size(S,1)
            for j=1:size(S,2)
                S(i,j,:) = squeeze(ifftshift(S(i,j,:)));
            end
        end
    end
    
    if trajectory==SPHERICAL
        history = zeros(2);
        stop = 0;
        S = zeros(phantom_size);
        scale = 8;
        p0 = [0 0 0];
        p1 = size(S);
        t = 0:1/n:1;
        while (sum(S(:))==0 || abs(mean(S(:))-coverage)>tol) && ~stop
            if mean(S(:))>coverage
                scale = scale+1;
            else
                scale = scale-1;
            end
            if scale==history(1,1)
                scale = history(1,history(2,:)==min(history(2,:)));
                stop = 1;
            end
            S = zeros(phantom_size);
            for sx=0:scale:p1(2)
                for sy=0:scale:p1(1)
                    for sz=0:scale:p1(3)
                        a = p0;
                        b = p1;
                        a(1) = max(1,min(p1(2),a(1)+sx));
                        a(2) = max(1,min(p1(1),a(2)+sy));
                        a(3) = max(1,min(p1(3),a(3)+sz));

                        b(1) = max(1,min(p1(2),b(1)-sx));
                        b(2) = max(1,min(p1(1),b(2)-sy));
                        b(3) = max(1,min(p1(3),b(3)-sz));
                        d = b-a;
                        x = max(1,min(p1(2),round(a(1)+d(1).*t+line_std*randn)));
                        y = max(1,min(p1(1),round(a(2)+d(2).*t+line_std*randn)));
                        z = max(1,min(p1(3),round(a(3)+d(3).*t+line_std*randn)));
                        for i=1:length(t)
                            S(y(i),x(i),z(i)) = 1;
                        end
                    end
                end
            end
            history = circshift(history,[0 1]);
            history(1,end) = scale;
            history(2,end) = abs(mean(S(:))-coverage);
        end
        for i=1:size(S,3)
            S(:,:,i) = ifftshift(S(:,:,i));
        end
        for i=1:size(S,1)
            for j=1:size(S,2)
                S(i,j,:) = squeeze(ifftshift(S(i,j,:)));
            end
        end
    end
    
    if trajectory==RAND_LINES
        S = zeros(phantom_size);
        p = coverage;
        for i2=1:n
            l = rand(1,n);
            S(l<p,:,i2) = 1;
            S(:,:,i2) = ifftshift(S(:,:,i2));
        end
    end
    
    if trajectory==RAND_PTS
        S = zeros(phantom_size);
        p = coverage;
        for i2=1:n
            l = rand(n);
            S(:,:,i2) = double(l<p);
            S(:,:,i2) = ifftshift(S(:,:,i2));
        end
    end
    
    if trajectory==LOW_PASS
        S = zeros(phantom_size);
        p = round((n-round(sqrt(n*n*coverage)))/2);
        for i2=1:n
            S(p:n-p,p:n-p,i2) = 1;
            S(:,:,i2) = ifftshift(S(:,:,i2));
        end
    end
    
    if trajectory==COMPLETE
        S = ones(phantom_size);
    end
    
end

function [ S ] = radial( trajectory, n, scale, i2, rot_deg, line_nbr, line_std ) 
    S = zeros(n);
    L = round(scale); % number of radial lines in the Fourier domain
    for i3=1:line_nbr
        if trajectory==1
            % aperture encompassing all scanning angles (aperture<pi is limited angle)
            aperture  = (pi/180)*180;
            % direction of the scanning beam (middle axis)
            direction = (pi/180)*0+(i2-1)*rot_deg+(i3-1)*2*pi/line_nbr;
        else
            aperture=(pi/180)*90;
            direction=(pi/180)*45+(i2-1)*rot_deg+(i3-1)*2*pi/line_nbr;
        end
        if (pi-aperture)>(aperture/L)
            thc = linspace(-direction-aperture/2, -direction+aperture/2, L);
        else
            thc = linspace(-direction-pi/2, -direction+pi/2-pi/L, L);
        end
        thc = mod(thc,pi);
        for ll = 1:L
            if ((thc(ll) <= pi/4) || (thc(ll) > 3*pi/4))
                yr = round(tan(thc(ll))*(-n/2+1:n/2-1)+n/2+1);
                yr = round(yr + line_std*randn(size(yr)));
                for nn = 1:n-1
                    yr(nn) = max(1,min(yr(nn),n));
                    S(yr(nn),nn+1) = 1;
                end
            else
                xc = round(cot(thc(ll))*(-n/2+1:n/2-1)+n/2+1);
                xc = round(xc + line_std*randn(size(xc)));
                for nn = 1:n-1
                    xc(nn) = max(1,min(xc(nn),n));
                    S(nn+1,xc(nn)) = 1;
                end
            end
        end
    end
end

function [ S ] = spiral( t, n, scale, i2, rot_deg, line_nbr, line_std )
    S = zeros(n);
    for i3=1:line_nbr
        temp = zeros(n);
        x = scale*t .* cos( t+(i2-1)*rot_deg+(i3-1)*2*pi/line_nbr );
        y = scale*t .* sin( t+(i2-1)*rot_deg+(i3-1)*2*pi/line_nbr );
        for i1=1:length(x)
            if round(x(i1)+n/2)<1 || round(y(i1)+n/2)<1 || ...
                    round(x(i1)+n/2)>n || round(y(i1)+n/2)>n
                break
            end
            xr = max(1,min(n,round(x(i1)+n/2)));
            yr = max(1,min(n,round(y(i1)+n/2)));
            if temp(yr,xr)==0
                temp(yr,xr) = 1;
                xr = max(1,min(n,round(x(i1)+n/2+line_std*randn)));
                yr = max(1,min(n,round(y(i1)+n/2+line_std*randn)));
                S(yr,xr) = 1;
            end
        end
    end
end

function [ S ] = log_spiral( t, n, scale, i2, rot_deg, line_nbr, line_std )
    S = zeros(n);
    for i3=1:line_nbr
        temp = zeros(n);
        x = ( exp(scale*t) .* cos(t+(i2-1)*rot_deg+(i3-1)*2*pi/line_nbr) );
        y = ( exp(scale*t) .* sin(t+(i2-1)*rot_deg+(i3-1)*2*pi/line_nbr) );
        % ensure DC is taken
        S(round(n/2),round(n/2)) = 1;
        temp(round(n/2),round(n/2)) = 1;
        for i1=1:length(x)
            if round(x(i1)+n/2)<1 || round(y(i1)+n/2)<1 || ...
                    round(x(i1)+n/2)>n || round(y(i1)+n/2)>n
                break
            end
            xr = max(1,min(n,round(x(i1)+n/2)));
            yr = max(1,min(n,round(y(i1)+n/2)));
            if temp(yr,xr)==0
                temp(yr,xr) = 1;
                xr = max(1,min(n,round(x(i1)+n/2+line_std*randn)));
                yr = max(1,min(n,round(y(i1)+n/2+line_std*randn)));
                S(yr,xr) = 1;
            end
        end
    end
end

