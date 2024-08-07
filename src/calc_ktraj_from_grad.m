function kTraj = calc_ktraj_from_grad (grad, isTxBool, dt, method)

% CALC_KTRAJ_FROM_GRAD Calculate the k-space trajectory given the gradient
% waveform.
%
% Usage: kTraj = calc_ktraj_from_grad (grad, isTxBool, dt, method)
%
% Example: kTraj = calc_ktraj_from_grad (grad);
%
% Returns
% -------
% kTraj: has the same size as grad.
%
% Expects
% -------
% grad: required and its second dim *MUST* be the time.
% 
% isTxBool: true for excitation k-space trajectory (true by default); false for
% reception k-space trajectory.
% dt: temporal resolution (sec). 4e-6 by default.
% 
% method: a char specifying how to perform the integration. it can be 's' for
% cubic spline interpolation, or 'c' for cumulative sum. defaults to 's'.
% 
% 
% Methods
% -------
% Excitation: kTraj(t) = -gamma * int_t^T grad(s)ds
% Reception:  kTraj(t) =  gamma * int_0^t grad(s)ds
% 
%Copyright (C) 2007 by Xiaoping Wu, Wed Jan  3 00:24:13 2007

if isempty(grad)
    warning('The gradient waveform is empty!')
end

if nargin < 2
  isTxBool = true;
end

if nargin < 3
  dt = 4e-6;                          % sec, default in Varian system.
end
  
if nargin < 4
  method = 's';
end


  gamma = 2.675e8;                    % gyromagnetic ratio for proton.
  
  
  % the k-space trajectory is calculated here.
  % since the 2nd dim of the gradient waveform array is assumed for time, function
  % FLIPLR is used to flip the arrays left-right accordingly.
  

  kTraj = zeros(size(grad));

  switch method
    
   case 's'
    if isTxBool
      for ind= 1: size(grad,1),
        kTraj(ind,:) = fitForKstX(grad(ind,:),dt,gamma);
      end
      
    else                                  % reception
      for ind= 1: size(grad,1),
        kTraj(ind,:) = fitForKstR(grad(ind,:),dt,gamma);
      end
    end
      
   case 'c'
    if isTxBool      
      kTraj = fliplr(-gamma .* dt .* cumsum(fliplr(grad),2));
    else                                  % reception
      kTraj = gamma .* dt .* cumsum(grad,2);
    end
    
   otherwise
   error('Invalid method..')
  end


function ktraj = fitForKstX(grad,dt,gamma)
% fit a 1d grad and calc the corresponding excitation kspace

grad = fliplr(grad);

x_1 = (1:numel(grad));

% --- Create fit

ft_ = fittype('cubicspline');

% Fit this model using new data
cf_ = fit(x_1(:),grad(:),ft_);

ktraj = -gamma.*dt.* integrate(cf_,x_1,x_1(1));
ktraj = ktraj(:).';
ktraj = fliplr(ktraj);

% 
function ktraj = fitForKstR(grad,dt,gamma)
% fit a 1d grad and calc the corresponding reception kspace

x_1 = (1:numel(grad));

% --- Create fit

ft_ = fittype('cubicspline');

% Fit this model using new data
cf_ = fit(x_1(:),grad(:),ft_);

ktraj = gamma.*dt.* integrate(cf_,x_1,x_1(1));
ktraj = ktraj(:).';
