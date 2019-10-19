%setup.m
%This function simply compiles the appropriate noiselet .mex files

cd('./utility')
mex('realnoiselet.c')
cd('..')
