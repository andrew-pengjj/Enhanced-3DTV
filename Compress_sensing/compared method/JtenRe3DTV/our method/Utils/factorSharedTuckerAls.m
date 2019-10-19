function [celT,celU,celCore,out] = factorSharedTuckerAls(celX,rankArr,mode_s,varargin)
%  min \sum_k ||X_k - [D_k; U_1k,U_2k,U_3,U_4k ||^2  
%  s.t. U_1k,U_2k,U_3,U_4k are orthogonal. 
%  where U_3 is shared factor matrix.
%
% ---INPUTs
% celX : cell array and 4-order tensor is contained in each component 
% rankArr: tucker rank for decomposation of each X_k
% mode_s : the shared mode
% varargin:  factorSharedTuckerAls(celX,rankArr,'param'value,...)
%            specifies optional parameters and values. 
%            Valid parameters and their default values are:
%            'tol' - Tolerance on difference in fit {1.0e-4}
%            'maxiters' - Maximum number of iterations {100}
%            'init' - Initial guess [{'random'}|'nvecs'|cell array]
%            'printitn' - Print fit every n iterations {1}
%---OUTPUTs
% celT : the decompostion for each X_k
% celU,celCore: factor matrices and core tensors for each X_k
% out:
%  - iter
%  - time
%  - fitchg
%
%---Syntax:   
%  [celT,celU,celCore,out] = factorSharedTuckerAls(celX,rankArr,mode_s); 
%
%  [celT,celU,celCore,out] = factorSharedTuckerAls(celX,rankArr,mode_s,'tol',1e-4);                                    
%
% wenfei cao and yao wang ( caowenf2006@163.com )
% copyright@xi'an jiaotong university 2014-9-8
% 
%
%% Extracting dim and compute norm
nCls = numel(celX);
nDim = zeros(nCls,1);
for k = 1 : nCls
    nDim(k) = ndims(celX{k});
end
arrNormX = zeros(nCls,1);
for k = 1 : nCls
  arrNormX(k) = norm(tensor(celX{k}));
end
%% Set algorithm parameters from input or by using defaults
params = inputParser;
params.addParamValue('tol',1e-4,@isscalar);
params.addParamValue('maxiters',50,@(x) isscalar(x) & x > 0);
params.addParamValue('initU', 'random', @(x) (iscell(x) || ismember(x,{'random','nvecs','eigs'})));
params.addParamValue('printitn',1,@isscalar);
params.parse(varargin{:});

%% Copy from params object
fitchgtol    = params.Results.tol;
maxiters     = params.Results.maxiters;
initU        = params.Results.initU; % nested cell
printitn     = params.Results.printitn;

%% Error checking 
% Error checking on the rank setting of the shared factor
if ~all(rankArr(:,mode_s) == rankArr(1,mode_s)); 
    error('the rank in third shard factor should be equal!');
end;

% Error checking on maxiters
if maxiters < 0
    error('OPTS.maxiters must be positive');
end
%% Set up and error checking on initial guess for U and shared Us
if iscell(initU) 
    % initU{1}{1-4}    = {U1,U2,U3,U4};
    % initU{2}{1-4}    = {U1,U2,U3,U4}
    % initU{nCls}{1-4} = {U1,U2,U3,U4]  note U3 is shared factor
    if numel(initU) ~= nCls
        error('OPTS.initU does not have %d cells',nCls);
    end
    for k = 1 : nCls
        for n  = 1 : nDim(k)
            if  ~ isequal(size(initU{k}{n}), [size(celX{k},n),rankArr(k,n)])
                   error('Opts.initU{%d}{%d} is the wrong size', k,n); 
            end
        end
   end
   
else
    if strcmp(initU,'random')
        initU = cell(nCls,1);
        for k = 1 : nCls
            for n = 1 : nDim(k)
                 if n ~= mode_s
                       initU{k}{n} = rand(size(celX{k},n),rankArr(k,n));
                 end
            end
        end 
        U_s   = rand(size(celX{k},mode_s),rankArr(1,mode_s));
        for k = 1 : nCls
            initU{k}{mode_s}  = U_s;
        end
        
    elseif strcmp(initU,'nvecs') || strcmp(initU,'eigs') 
        % Compute an orthonormal basis for the dominant
        % Rn-dimensional left singular subspace of
        % X_(n) (1 <= n <= N).
        initU = cell(nCls,1);
        for k = 1 : nCls
           for n = 1 : nDim(k)
               if n ~= mode_s
                   initU{k}{n} = nvecs(tensor(celX{k}),n,rankArr(k,n));
               end
           end
        end
        U_s = nvecs(tensor(celX{1}),mode_s,rankArr(1,mode_s));
        for k = 1 : nCls
            initU{k}{mode_s} = U_s;
        end
    else
        error('The selected initialization method is not supported');
    end
end

%% Set up for iterations - initializing U, Us, and the fit.
U   = initU;
fit = 0;

if printitn > 0
    fprintf('\n Factor Shared Tucker Alternating Least-Squares:\n');
end

%% Main Loop: Iterate until convergence
tic;
rank_s   = rankArr(1,mode_s);
core     = cell(nCls,1);
for iter = 1 : maxiters  
   
    fitold = fit;
    
    %-----------------iterate for each cluster ---------------------%
    for k = 1 : nCls

        clsX = tensor(celX{k});  

        %%-update each factor in each cluster
        for n = 1 : nDim(k)
              if n ~= mode_s
                    Utilde = ttm(clsX,U{k},-n,'t');

                    % Maximize norm(Utilde x_n W') wrt W and
	                % keeping orthonormality of W
                    U{k}{n} = nvecs(Utilde,n,rankArr(k,n));
              end
       end
   end

   %%-undating the shared factor 
   U_mat = sharedFactorUpdating(celX,U,nCls,mode_s,rank_s);
   for jj = 1 : nCls
       U{jj}{mode_s} = U_mat;

   %%-assemble the current approximation for the k-th cluster
       core{jj}      = ttm(tensor(celX{jj}),U{jj},'t');
   end
  
    %-------------------Compute fit ---------------------------------%
    normresidual = zeros(nCls,1);
    for k = 1 : nCls
        normresidual(k) = sqrt( abs(arrNormX(k)^2 - norm(core{k})^2) );
    end
    fit    = mean(1 - (normresidual./arrNormX));
    fitchg = abs(fitold - fit);
        
    if printitn > 0 && mod(iter,printitn) == 0
        fprintf(' Iter %2d: fit = %e fitdelta = %7.1e\n', iter, fit, fitchg);
    end
             
    %------------------Check for convergence ------------------------%
    if (iter > 1) && (fitchg < fitchgtol)
        break;
    end

end

celU    = U;
celCore = core;
celT    = cell(nCls,1);
for k = 1 : nCls
    celT{k} = double(ttensor(core{k}, U{k}));
end
out.iter   = iter;
out.time   = toc;
out.fitchg = fitchg;
out.fit    = fit;
return;

end

function U_mat = sharedFactorUpdating(celX,U,nCls,mode_s,rank_s)
% max_U \sum_{i} ||U*Z_i||^2 s.t. U'U = I;
% This problem is equivalent to the following problem:
% max_U Tr(U(\sum_i Z_i *Z'_i) U') s.t. U'U = I
%
N_s   = size(celX{1},mode_s);
Z     = zeros(N_s,N_s);
for k = 1 : nCls
    tmpU     = ttm(tensor(celX{k}),U{k},-mode_s,'t');
    tmpU_mat = double(tenmat(tmpU,mode_s));
    Z        = Z + tmpU_mat*tmpU_mat';
end

[U, val]  = eigs(Z,rank_s);
U_mat     = U(:,1:rank_s);
end

