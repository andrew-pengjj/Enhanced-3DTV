function ergas = ErrRelGlobAdimSyn(imagery1, imagery2)

%==========================================================================
% Evaluates Erreur Relative Globale Adimensionnelle de Synth¨¨se (ERGAS, 
%   dimensionless global relative error of synthesis)[1] for two MSIs.
%
% Syntax:
%   ergas = ErrRelGlobAdimSyn(imagery1, imagery2)
%
% Input:
%   imagery1 - the reference MSI data array
%   imagery2 - the target MSI data array
% NOTE: MSI data array  is a M*N*K array for imagery with M*N spatial
%	pixels, K bands and DYNAMIC RANGE [0, 255]. If imagery1 and imagery2
%	have different size, the larger one will be truncated to fit the
%	smaller one.
%
% [1] L. Wald, "Data Fusion: Definitions and Architectures¡ªFusion of
%     Images of Different Spatial Resolutions", Les Presses, Ecole des
%     Mines de Paris, Paris, France (2002) 200 pp.
%
% See also StructureSIM, FeatureSIM and SpectAngMapper
%
% by Yi Peng
%==========================================================================

[m, n, k] = size(imagery1);
[mm, nn, kk] = size(imagery2);
m = min(m, mm);
n = min(n, nn);
k = min(k, kk);
imagery1 = imagery1(1:m, 1:n, 1:k);
imagery2 = imagery2(1:m, 1:n, 1:k);

ergas = 0;
for i = 1:k
    ergas = ergas + mse(imagery1(:, :, i) - imagery2(:, :, i))/mean2(imagery1(:, :, i));
end
ergas = 100*sqrt(ergas/k);