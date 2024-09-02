function [rf,grad,localSAR] = design_pTxSPSP_RF(ID, RF_duration, dt, TR, RFA, type)
% rf_design_function - MATLAB function for designing a non-selective pTx-SPSP RF pulse
%
% Inputs:
%   ID (array) - serial number of calibration, set here so the function
%   could be easily extended to UP training.
%   RF_duration (double) - Duration of the RF pulse, default to 2.88e-3s
%   dt (double) - Dwell time, default to 10e-6s
%   TR (double) - Repetition time, default to 50e-3s
%   RFA - Round flip angle , default ernstAngle at 7T
%   type (string) - Type of RF pulse, default to 'SPINS'
%
% Outputs:
%   rf, grad, localSAR - Resulting RF & grad pulse based on the input parameters, and
%   local SAR according to the vendor-provided 'SARDataUser'
%
% Example:
%   result = rf_design_function(); % Uses default parameters
%   result = rf_design_function(5, 2.88, 10e-6,50e-3, ernstAngle(50e-3), 'KT'); % Uses custom parameters. Pay
%   attention the RF_dutation and dwell time can not be changed in the this demo.

% Set default values
if nargin < 1
    ID = [1];
end
if nargin < 2
    RF_duration = 2.88e-3;
end
if nargin < 3
    dt = 10e-6;
end
if nargin < 4
    TR = 50e-3;
end
if nargin < 5
    RFA = round(ernstAngle(TR));
end
if nargin < 6
    type = 'SPINS';
end

% Display parameter information
fprintf('ID: %d \n', ID);
fprintf('RF Duration: %.2f ms\n', RF_duration*1000);
fprintf('RF Type: %s\n', type);

% Perform different RF designs based on the type
switch type
    case 'SPINS'
        % Call SPINS RF design method
        [rf,grad,localSAR] = design_SPINS_rf(ID,RF_duration,dt, TR, RFA);
    case 'KT'
        % Call Kt-points RF design method
        [rf,grad,localSAR] = design_KT_rf(ID,RF_duration,dt,TR, RFA);
    otherwise
        error('Unknown RF type: %s', type);
end
end