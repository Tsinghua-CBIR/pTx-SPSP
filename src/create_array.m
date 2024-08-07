function [b1plusArrayOfMask, posArrayOfMask] = create_array (b1plusMapArray, spatialMask, ...
                                                  fieldOfExcitation, poffset)
% CREATE_ARRAY Create B1+ and position arrays.
%
% Usage: [b1plusArrayOfMask, posArrayOfMask] = create_array (b1plusMapArray,
% spatialMask, fieldOfExcitation, poffset)
%
% Returns
% -------
% b1plusArrayOfMask: Ns-by-Nc array.
% posArrayOfMask: Ns-by-2,3 matrix.
%
% Expects
% -------
% b1plusMapArray: 
% spatialMask: 
% fieldOfExcitation: [fox_y, fox_x, fox_z] (m).
% 
% poffset: [offsety,offsetx,offsetz] in m specifying the offset of FOV with
% respect to grad isocenter. defaults to [0 0 0].
% 
% 
% See also: construct_system_matrix
% 

if nargin < 3
  error('Too few input arguments!')
end
if nargin < 4
  poffset = [0, 0, 0];
end


b1plusArrayOfMask = local_calc_b1plus_array(b1plusMapArray,spatialMask);
posArrayOfMask = local_calc_position_array(spatialMask, fieldOfExcitation, poffset);


% ================
% Local functions
% ================

function b1plusArrayOfMask = local_calc_b1plus_array (b1plusMapArray, spatialMask)
% LOCAL_CALC_B1PLUS_ARRAY Calculate B1+ array.
%
% Usage: b1plusArrayOfMask = calc_b1plus_array (b1plusMapArray, spatialMask)
%
% Returns
% -------
% b1plusArrayOfMask: Ns-by-Nc array.
%
% Expects
% -------
% b1plusMapArray: 
% spatialMask: 
    
%  nDims = length(size(spatialMask));
  nCoils = size(b1plusMapArray,4);      % the last dim is the num of coils.

  for iCoil = 1 : nCoils,
    iB1plusMap = b1plusMapArray(:, :, :, iCoil);
    b1plusArrayOfMask(:,iCoil) = iB1plusMap(spatialMask);
  end

function posArrayOfMask = local_calc_position_array (spatialMask, fieldOfExcitation, ...
                                                  poffset)
% LOCAL_CALC_POSITION_ARRAY Calculate position array of mask from spatial mask and fox.
%
% Usage: posArrayOfMask = calc_position_array (spatialMask, fieldOfExcitation)
%
% Returns
% -------
% posArrayOfMask: Ns-by-2,3 matrix.
%
% Expects
% -------
% spatialMask: logical.
% fieldOfExcitation: [fox_x, fox_y, fox_z] (m).
      
  dim = size(spatialMask);
  nDims = length(dim);
  posVectorCell = cell(1, nDims);

  for iDim = 1 : nDims,
    posVectorCell{iDim} = linspace(-fieldOfExcitation(iDim)./2, ...
                                   fieldOfExcitation(iDim)./2,dim(iDim)) + poffset(iDim);
  end
  
  indexOfMask = find(spatialMask);
  posArrayOfMask = zeros(length(indexOfMask), nDims);

  if (nDims == 2)
    [d1Array, d2Array] = ndgrid(posVectorCell{1}, posVectorCell{2}); %
    posArrayOfMask(:, 1) = d1Array(indexOfMask);
    posArrayOfMask(:, 2) = d2Array(indexOfMask);
    
  else                                    % nDims = 3
    [d1Array, d2Array, d3Array] = ndgrid(posVectorCell{1}, posVectorCell{2}, ...
                                        posVectorCell{3});
    posArrayOfMask(:, 1) = d1Array(indexOfMask);
    posArrayOfMask(:, 2) = d2Array(indexOfMask);
    posArrayOfMask(:, 3) = d3Array(indexOfMask);
  end
