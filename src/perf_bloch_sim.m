function [mxy, Mz] = perf_bloch_sim (Brt, Mr0, dt)

% PERF_BLOCH_SIM Perform the bloch simulation, based on spin-domain
% representation, given the overall B field as a function of time and space.
%
% Usage: [mxy, Mz] = perf_bloch_sim (Brt, Mr0, dt)
%
% Returns
% -------
% mxy: np x 1 complex vector containing transverse magnetization, i.e., Mx + iMy.
% Mz: np x 1 vector
%
% Expects
% -------
% Brt: 3 x np x nt array containing the overal B field as a function of space
% and time. in Tesla
% 
% Mr0: 3 x np array representing the init magnetization as a function of space.
% M0 = [Mx0; My0; Mz0].
% 
% dt: temporal resolution in sec.
%
%
% See also: run_bloch_sim run_bloch_sim1d
%
%
% Copyright (C) 2008 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Wed Jul  2 00:47:18 2008
%

np = size(Brt,2);
mxy = complex(zeros(np,1));
Mz = zeros(np,1);

mr0 = [Mr0(1,:)+1i*Mr0(2,:); Mr0(1,:)-1i*Mr0(2,:); Mr0(3,:)];

fprintf('=> %d spatial points detected...\n',np);

parfor idx = 1:np,
  iBt = squeeze(Brt(:,idx,:));
  js = simBlochEq(iBt,dt);
  im = calcM(mr0(:,idx),js);
  
  mxy(idx) = im(1);
  Mz(idx) = im(3);  
end

%% ====================
%%
%%  Local functions
%%
%% ====================

function m = calcM(m0,s)
% calc final magnetization given the init magnetization and final state vector.
% here m = [mxy, mxy*, Mz]', mxy = Mx + iMy.

alfa = s(1);
beta = s(2);

xmat = [conj(alfa)^2           -beta^2       2*conj(alfa)*beta; ...
       -conj(beta)^2            alfa^2       2*alfa*conj(beta); ...
       -conj(alfa)*conj(beta)  -alfa*beta    abs(alfa)^2-abs(beta)^2];

m = xmat * m0;

% -----------------------

function s = simBlochEq(Bt, dt)
% calc the final 2x1 state complex vector given the total B field as a function
% of time (Bt is of 3xnt).

Bt_nz = Bt(:,sum(abs(Bt))>0);
nt = size(Bt_nz,2);
s = [1;0]; % init state is no rotation.

for idx = 1:nt,
  iB = Bt_nz(:,idx);
  iQ = calcRotMat(iB, dt);
  s = iQ * s;
end

% -----------------------

function Q = calcRotMat(B, dt)
% calc the 2x2 unitary spin domain rotation matrix Q given the total B field
% [Bx,By,Bz]'.

gamma = 2.675e8; %

phi = -gamma * dt * norm(B); % flip angle
n = gamma * dt * B / abs(phi); % rotation axis

% Cayley-Klein params
a = cos(phi/2) - 1i* n(3)* sin(phi/2);
b = -1i* (n(1) + 1i*n(2))* sin(phi/2);

Q = [a -conj(b); b conj(a)];
