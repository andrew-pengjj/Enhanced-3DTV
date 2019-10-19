function idx = tt_sub2ind(siz,subs)
%TT_SUB2IND Converts multidimensional subscripts to linear indices.
%
%   INDS = TT_SUB2IND(SIZ,SUBS) returns the linear indices
%   equivalent to the subscripts in the array SUBS for a tensor of
%   size SIZ.  
%
%   See also TT_IND2SUB, SUB2IND.
%
%MATLAB Tensor Toolbox.
%Copyright 2010, Sandia Corporation. 

% This is the MATLAB Tensor Toolbox by Brett Bader and Tamara Kolda. 
% http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2010) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in tensor_toolbox/LICENSE.txt
% $Id: tt_sub2ind.m,v 1.9 2010/03/19 23:46:31 tgkolda Exp $

mult = [1 cumprod(siz(1:end-1))];
idx = (subs - 1) * mult' + 1;

