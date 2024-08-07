function [x, rho, eta, ni] = solve_mlstr (A, b, lambda, tol, x0)

% SOLVE_MLSTR Solve magnitude least squares inverse problem based on tikhonov
% regularization.
%
% Usage: [x, rho, eta, ni] = solve_mlstr (A, b, lambda, tol, x0)
%
% Returns
% -------
% x: solution
% rho: residual norm
% eta: solution norm
% ni: number of iterations
%
% Expects
% -------
% A: system matrix
% 
% b: target vector, only magnitude makes sense
% 
% lambda: regularization param. can be a vector containing a number of
% regularization params.
% 
% tol: tolerance controlling convergence. defaults to 1e-5;
% 
% x0: initial point, defaults to 0.
%
%
% See also: matrix tikhonov solve_lstr solve_mlstr_err solve_mscgls
%
%
% Copyright (C) 2009 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Sat Sep 19 19:23:21 2009
%

if nargin < 4
  tol = 1e-5;
end

A(isnan(A))=0;

if ~isa(A,'matrix')
  A = matrix(A);
end
b = abs(b(:));                    
[U,s,V] = A.getSVD;
% [U,s,V] = svd(A);

if nargin< 5
  x0 = zeros(size(V,1),1);
end
nlambda = length(lambda);
x = complex(zeros(size(V,1),nlambda));
rho = zeros(nlambda,1);
eta = rho;
ni = rho;

% fprintf('=> Computation started for %d lambda(s) ...\n',nlambda);
parfor idx = 1: nlambda,
  [x(:,idx), ni(idx), rho(idx), eta(idx)] = solveMLS (U,s,V,b,lambda(idx),tol,x0);
end
% disp('=> Done...');

%% ====================
%%
%%  Local functions
%%
%% ====================

function [x, ni, rho, eta] = solveMLS (U, s, V, b, lambda, tol, x0)
% solve magnitude least squares tikhonov regularization problem.

A = U*diag(s)*V';

z = exp(1i*angle(A*x0));
x = solveTikhonov(U,s,V, b.*z, lambda);
costNew = calcCost(A, x, b.*z, lambda);
          
costOld = 1e100; % inf
ni = 1;
while ((costOld- costNew)/ costOld > tol),
  ni = ni + 1;
  costOld = costNew;
  
  z = exp(1i*angle(A*x));  
  x = solveTikhonov(U,s,V, b.*z, lambda);
  costNew = calcCost(A, x, b.*z, lambda);
end

[rho,eta] = calcNorm (A, x, b);
% -----------------------

function c = calcCost (A, x, b, lambda)
% Purpose: evaluate cost func
r = [A*x - b; lambda*x];
c= r'*r;
% -----------------------

function [rho, eta] = calcNorm (A, x, b)
% calc residual norm | |A*x| - b | and solution norm | x |
rho = norm(abs(A*x) - b);
eta = norm(x);
% -----------------------

function x = solveTikhonov (U,s,V,b,lambda)
% solve tikhonov regularization problem
beta = U'* b;
zeta = s.* beta;  
x = V* (zeta./ (s.^2 + lambda.^2));
