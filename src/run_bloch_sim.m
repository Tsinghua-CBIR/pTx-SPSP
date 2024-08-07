function [mxypat, mzpat] = run_bloch_sim (varargin)

% RUN_BLOCH_SIM Run bloch simulation based on spinor domain rotation
% representation. 
% 
% Usage: [mxypat, mzpat] = run_bloch_sim
% (pulsearr,grad,b1maparr,mask,fox,b0map,z0,b0,dt,poffset,mr0)
%
% Returns
% -------
% mxypat: 2D or 3D complex transverse magnetization (i.e., Mx + iMy) pattern.
% mzpat: 2D or 3D Mz pattern.
%
% Expects
% -------
% pulsearr: nTxchs x ntimepts array
% grad: 2 or 3 x ntimepts array (tesla/m)
% b1maparr: 4D array with the last dim for the num of coils.
% mask: spatial mask, 2D or 3D logical.
% fox: field of excitation [fox_x, fox_y, fox_z] (m)
% b0map: 2D or 3D b0 inhomogeneity map in tesla, and it defaults to [].
% 
% z0: Z value of grad isocenter in m with respect to the middle of fox. ignored
% for 2D simulations. defaults to 0.
% 
% b0: a 1 x ntimepts vector for b0 eddy current, in tesla. defaults to [].
% dt: temporal resolution in sec. defaults to 4 us for varian.
%
% poffset: [offsety,offsetx,offsetz] in mm specifying the offset of FOV with
% respect to grad isocenter. defaults to [0 0 0].
% 
% mr0: 3D or 4D array representing the init magnetization as a function of
% space. The first dim MUST be 3 respectly for Mx, My and Mz. If not specified,
% the init magnetization is assumed to be the equilibrium one, i.e., M0=[0 0 1]'
% for everywhere.
%
%
% See also: perf_bloch_sim create_array run_bloch_simulation run_bloch_sim1d
%
%
% Copyright (C) 2008 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Wed May 14 16:58:08 2008
%
  
[pulsearr,grad,b1maparr,mask,fox,b0map,z0,b0,dt,poffset,Mr0] = parse_inputs(varargin{:});

disp('=> The bloch simulation starts...');

[b1arr,posarr] = create_array(b1maparr,mask,fox,1e-3*poffset);
if size(posarr,2) > 2 % 3D
  posarr(:,3) = posarr(:,3) - z0;
end

B1rt = b1arr * pulsearr;
Bzrt = posarr*grad + b0map(mask)*ones(1,size(grad,2)) + ones(size(b0map(mask)))*b0;

Brt(1,:,:) = real(B1rt);
Brt(2,:,:) = imag(B1rt);
Brt(3,:,:) = Bzrt;

[mxy,mz] = perf_bloch_sim(Brt,Mr0,dt);

mxypat = complex(zeros(size(mask)));
mzpat = zeros(size(mask));
mxypat(mask) = mxy;
mzpat(mask) = mz;

disp('=> The bloch simulation ends.');

% -----
function [pulse,grad,b1map,mask,fox,b0map,z0,b0,dt,poffset,Mr0] = parse_inputs(varargin)
% defaults.
  b0map = [];
  z0 = 0;
  b0 = [];
  dt = 4e-6; % in s
  poffset = [0, 0, 0];
  mr0 = [];
  
  error(nargchk(5,11,nargin,'struct'));
  
  switch nargin
   case 5,
    pulse = varargin{1};
    grad = varargin{2};
    b1map = varargin{3};
    mask = varargin{4};
    fox = varargin{5};
   case 6,
    pulse = varargin{1};
    grad = varargin{2};
    b1map = varargin{3};
    mask = varargin{4};
    fox = varargin{5};
    b0map = varargin{6};
   case 7,
    pulse = varargin{1};
    grad = varargin{2};
    b1map = varargin{3};
    mask = varargin{4};
    fox = varargin{5};
    b0map = varargin{6};
    z0 = varargin{7};
   case 8,
    pulse = varargin{1};
    grad = varargin{2};
    b1map = varargin{3};
    mask = varargin{4};
    fox = varargin{5};
    b0map = varargin{6};
    z0 = varargin{7};
    b0 = varargin{8};
   case 9,
    pulse = varargin{1};
    grad = varargin{2};
    b1map = varargin{3};
    mask = varargin{4};
    fox = varargin{5};
    b0map = varargin{6};
    z0 = varargin{7};
    b0 = varargin{8};
    dt = varargin{9};
   case 10,
    pulse = varargin{1};
    grad = varargin{2};
    b1map = varargin{3};
    mask = varargin{4};
    fox = varargin{5};
    b0map = varargin{6};
    z0 = varargin{7};
    b0 = varargin{8};
    dt = varargin{9};
    poffset = varargin{10};
   case 11,
    pulse = varargin{1};
    grad = varargin{2};
    b1map = varargin{3};
    mask = varargin{4};
    fox = varargin{5};
    b0map = varargin{6};
    z0 = varargin{7};
    b0 = varargin{8};
    dt = varargin{9};
    poffset = varargin{10};
    mr0 = varargin{11};
   otherwise
  end

  if isempty(mask)
      siz=size(b1map);
      mask=true(siz(1:3));
  end
  
  if isempty(b0map)
    b0map = zeros(size(mask));
  end
  
  if isempty(b0)
    b0 = zeros(1,size(grad,2));
  end
  
  if isempty(mr0)
    Mr0 = repmat([0 0 1]',[1 length(mask(mask))]);
  else
    mxx = mr0(1,:);
    myy = mr0(2,:);
    mzz = mr0(3,:);
    mxx=mxx(:);
    myy=myy(:);
    mzz=mzz(:);
    Mr0 = [mxx(mask(:)) myy(mask(:)) mzz(mask(:))].';
  end


  % Check the validity
  glen = size(grad,2);
  plen = size(pulse,2);
  b0len = length(b0);
  maxlen = max([glen,plen,b0len]);

  if glen < maxlen
    grad(:,maxlen) = 0;
  end
  if plen < maxlen
    pulse(:,maxlen) = 0;
  end
  if b0len < maxlen
    b0(:,maxlen) = 0;
  end

