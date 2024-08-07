function grad = design_grad_trapz (kdes, maxamp, maxsr, dt)

% DESIGN_GRAD_TRAPZ Design a trapzoidal or triangle gradient given a k value.
%
% Usage: grad = design_grad_trapz (kdes, maxamp, maxsr, dt)
%
% Returns
% -------
% grad: a vector for a trapzoidal or triangle gradient
%
% Expects
% -------
% kdes: desired k value in rad/m
% maxamp: max grad amp in T/m
% maxsr: max slew rate in T/m/s
% dt: dwell time in sec. defaults to 4 us
%
%
% See also: design_toptgrad1D.m
%
%
% Copyright (C) 2009 CMRR at UMN Author: Xiaoping Wu <xpwu@cmrr.umn.edu>
%Created: Tue Sep 8 14:43:38 2009

if nargin < 4
  dt = 4e-6; % s
end

if nargin < 3
  maxsr = 166; % tesla/m/s
end
          
if nargin < 2
  maxamp = 50e-3; % tesla/m
end

gamma = 2.675e8; 
gsign = sign(kdes);
kdes = abs(kdes);
gsurf= kdes./gamma;
if canUseTriangle(gsurf, maxamp, maxsr)
  grad = createTriangle(gsurf,maxsr,dt);
else
  grad = createTrapzoid(gsurf,maxamp,maxsr,dt); 
end

grad = gsign* grad;
grad= gsurf./(sum(grad)*dt) .* grad;

% -----------------------

function flag = canUseTriangle(gsurf,maxamp,maxsr)
% 
gsurfmax = maxamp.^2 / maxsr;
flag = gsurfmax > gsurf;

% -----------------------
function grad = createTriangle(gsurf, maxsr, dt)
% 
amp = calcAmp(gsurf, maxsr);
grad = [createRampUp(amp,maxsr,dt) createRampDown(amp,maxsr,dt)];

% -----------------------
function grad = createTrapzoid (gsurf, maxamp, maxsr, dt)
% 
  time = (gsurf - maxamp.^2/ maxsr)./ maxamp;
  npts = round(time/dt);
        
  gplat = maxamp* ones(1,npts);
        
  grad = [createRampUp(maxamp,maxsr,dt) gplat createRampDown(maxamp,maxsr,dt)];
      

% -----------------------
function amp = calcAmp (gsurf, maxsr)
amp = (maxsr* gsurf).^0.5;

% -----------------------
function ramp = createRampUp (amp, maxsr, dt)
% 
  time = amp/ maxsr;
  npts = floor(time/dt);
  
  ramp = (0:npts)* dt* maxsr;

function ramp = createRampDown (amp, maxsr, dt)
% 
ramp = fliplr(createRampUp(amp, maxsr, dt));

