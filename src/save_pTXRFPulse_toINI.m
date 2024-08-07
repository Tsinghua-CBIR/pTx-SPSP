function bSuccess = save_pTXRFPulse_toINI( g_low, b_est, shim_coeff, varargin)
%-------------------------------------------------------------------------------
% Copyright Siemens AG, all rights reserved, confidential
% Project:  NUMARIS/4
% File: \n4_physic_src\pkg\MrServers\MrPhysicsSrv_src\PtxPulseDesign\library_main\in_out\save_pTXRFPulse_toINI.m
% Author:   pfeujodj
% Descrip:  save pulse in pTXExcite-file-format
%
% Variable - Type - Unit
%  
% g_low   - samplesx3, double - [mT/m]
%     gradient trajectory (RO,PE,SS) without oversampling -> g_step
%  
% b_est  - (samples*coils)*npulses, complex - [V]
%  oversampled RF values -> oversample
% 
% shim_coeff -  coils*npulses, complex
%  RF-shimming coefficients
% 
% alpha - double, deg
%  nominal flip angle
%
% g_step - double - [ms]
%  gradient sampling time
%  
% oversample - int
%  RF oversample factor
%  
% rf_form - samplesx1, double
%  basic shape of the RF-pulse (descriptional samples for one channel)
% 
% rf_detail - sectionsx2, int
%  structural description of the RF-pulse sections (consecutive number of samples of type RF-nonoptimized (0), RF-package(1), RF-freeform (2))
% 
%
% History:  090701 - pfeujodj/setsompop - initial version (create_pTXRFPulse.m)
%           091223 - pfeujodj - update for RF shim values
%           100610 - rittdiz7 - integration in pulse design suite
%           100805 - pfeujodj - rename/cleanup (get_pTXRFPulse_fromINI.m / save_pTXRFPulse_toINI.m)
%           140723 - pfeujodj - update to multiple pulses (get_pTXRFPulse_fromINI.m / save_pTXRFPulse_toINI.m)
%-------------------------------------------------------------------------------

FCTNAME = 'save_pTXRFPulse_toINI'; nargVars = 3;
VERSION = 'v2.0, 2014-07-23 (jp)';         % update version info here

% global defines
bSuccess     = 0;       % return value   
iFloatDigits = 5;       % accuracy of written floats	
          
%%  option/narg handling: default options: define whole struct and comment
dopt.RFPULSE_ID         = 0;        % different ID for several pulses used within one sequence 
dopt.RFPULSE_COMMENT    = '-';      % info about creation / properties of pulse
dopt.FACTOROVERSAMPLE   = 1;        % rf dwell time compared to gradient dwell time (fGradRasterTime)
dopt.VERBOSE            = 1;
% RF pulse properties
dopt.INITIALPHASE       = 0.;
dopt.ASYMMETRY          = 0.5;
dopt.PULSENAME          = 'results/pTXRFPulse';     % base output file name
dopt.FAMILY             = 'pTX';
dopt.COMMENT            = 'pTX';
dopt.NOMFLIPANGLE       = -1;              % will be overwritten by the varargin value
dopt.AMPLINT            = 100.;
dopt.ABSINT             = 100.;
dopt.POWERINT           = 100.;
dopt.MINSLICE           = 1.;
dopt.MAXSLICE           = 1.;
dopt.REFGRAD            = 1.;
% Gradient properties
dopt.GRADRASTERTIME     = 10;       % [us] GRADIENT_RASTER_TIME is fix for sequences at 10 us


%% option/narg handling - error: narg > nargVars; error: narg < nargVars
[dopt, narg] 	= handlerOptionsNarg(dopt, varargin, nargin, [0 nargVars]);

%%  option/narg handling: translate options to local parameters
iRFPulseID          = dopt.RFPULSE_ID;
strRFPulseComment   = dopt.RFPULSE_COMMENT;
iFactorOversample   = dopt.FACTOROVERSAMPLE;
bVerbose            = dopt.VERBOSE;
fInitialPhase       = dopt.INITIALPHASE;
fAsymmetry          = dopt.ASYMMETRY;
strPulseName        = dopt.PULSENAME;                   % RF pulse properties 
strFamily           = dopt.FAMILY;          
strComment          = dopt.COMMENT;
fNominalFlipAngle   = dopt.NOMFLIPANGLE;    
fAmplInt            = dopt.AMPLINT;
fAbsInt             = dopt.ABSINT;
fPowerInt           = dopt.POWERINT;
fMinSlice           = dopt.MINSLICE;
fMaxSlice           = dopt.MAXSLICE;
fRefGrad            = dopt.REFGRAD;
fGradRasterTime     = dopt.GRADRASTERTIME;

if (nargin < nargVars)
    %%%help(FCTNAME); demoUsage(FCTNAME); return;   % activate this line for demo
    logprintf(1,'E','wrong number of arguments');
end
%%%%  end: option/narg handling: 

iNUsedChannels = 0;
iNMultiplePulses = 0;
if( isempty(g_low) )
    g_low = zeros(1,3);  % dummy fill to get sequence prep working
end
if( isempty(b_est) )    % can not determine iNUsedChannels from b_est
    if( ~isempty(shim_coeff))
        iNUsedChannels = size(shim_coeff,1);
        iNMultiplePulses = size(shim_coeff,2);
    else
        logprintf(1,'E','%s: can not determine <iNUsedChannels>.', FCTNAME);
    end
    b_est = zeros(iNUsedChannels,iNMultiplePulses);   % dummy fill to get sequence prep working
end
iGLength = size( g_low,1 );
iBLength = iGLength*iFactorOversample;
if( rem(size(b_est,1),iBLength) ~= 0 )
    logprintf(1,'E','%s: wrong data sizes (g_low, b_est).', FCTNAME);
end

g  = [g_low];
iNUsedChannels   = size(b_est,1)/iBLength;
iNMultiplePulses = size(b_est,2);
rf = reshape(b_est,[iBLength iNUsedChannels iNMultiplePulses]);

%% pad with extra zeros to ensure that last gradient point is ZERO (sequence requirement)
if( sum(abs(g(iGLength,:))) ~= 0)
    g  = [ g; zeros(1,size(g,2))];
    rf = [rf; zeros(iFactorOversample, size(rf,2), size(rf,3))];
end

if( ~isempty(shim_coeff) )      % array of complex RF Shim values 
    s_shim_coeff = size(shim_coeff);
    if( s_shim_coeff(1) ~= iNUsedChannels)
        logprintf(1,'E','%s: wrong size of shim_coeff (NUsedChannels=%d)', FCTNAME, iNUsedChannels);
    end
    vfPhaseOffset    = angle(shim_coeff)*180/pi;
    vfAmplitudeScale = abs(shim_coeff);  
else
    vfPhaseOffset    = [];
    vfAmplitudeScale = []; 
end

if( isempty(strRFPulseComment) ) 	% pulse description must be set!
    logprintf(1,'E','%s: set pulse description - option RFPULSE_COMMENT=string', FCTNAME, strRFPulseComment);
end

if( fNominalFlipAngle < 0 )         % design flip angle must be set!
%     logprintf(1,'E','%s: set design flip angle - option NOMFLIPANGLE=value', FCTNAME, fNominalFlipAngle);
end

%% 
iDimRF           = 2;                   % Amplitude/Phase of complex data
iDimGradient     = size(g,2);           % typically 3 (non-zero) gradient channels
iGradientSamples = size(g,1);
iSamples         = size(rf,1);

%% extract phase and amp of rf
abs_rfs     = abs(rf);
angle_rfs   = angle(rf);
index = angle_rfs<0;
angle_rfs(index) = angle_rfs(index)+2*pi;   % phase wrap

%% round g and rf
g         = round(g        *10^iFloatDigits)/10^iFloatDigits;
abs_rfs   = round(abs_rfs  *10^iFloatDigits)/10^iFloatDigits;
angle_rfs = floor(angle_rfs*10^iFloatDigits)/10^iFloatDigits; 

%% find max rf and gradient
fMaxAbsRF = max(max(abs_rfs(:,:,1)));   % maximum over all channels: FIRST Pulse!
fMaxAbsG1 = max(abs(g(:,1)));
fMaxAbsG2 = max(abs(g(:,2)));
fMaxAbsG3 = max(abs(g(:,3)));

%% RFPulseID is contained in filename, typically starting from zero
strFileName = sprintf('%s%d.ini', strPulseName, iRFPulseID);
fid = fopen( strFileName, 'w');

%% header at start
fprintf(fid, '#pTXRFPulse - created by: <%s><%s>', FCTNAME, VERSION);
fprintf(fid, '\n#%s', strRFPulseComment);

% here huge pseudocomment sections will be included 
%   for B0,B1 map, pulse design, SAR calculation ...

%% global parameter section
fprintf(fid,'\n#----------\n');
fprintf(fid,'[pTXPulse]\n');
fprintf(fid,'\n');
fprintf(fid,'NUsedChannels    = %d\n', iNUsedChannels);
fprintf(fid,'NMultiplePulses  = %d\n', iNMultiplePulses);
fprintf(fid,'DimRF            = %d\n', iDimRF);
fprintf(fid,'DimGradient      = %d\n', iDimGradient);
fprintf(fid,'MaxAbsRF         = %g\t\t # scaling for RF amplitude\n', fMaxAbsRF);
fprintf(fid,'InitialPhase     = %g\n', fInitialPhase);
fprintf(fid,'Asymmetry        = %g\n', fAsymmetry);
fprintf(fid,'\n');
fprintf(fid,'PulseName        = %s\t\t #standard RF pulse parameters\n', strPulseName);
fprintf(fid,'Family           = %s\n', strFamily);
fprintf(fid,'Comment          = %s\n', strComment);
fprintf(fid,'NominalFlipAngle = %g\n', fNominalFlipAngle);
fprintf(fid,'Samples          = %d\n', iSamples);
fprintf(fid,'AmplInt          = %g\n', fAmplInt);
fprintf(fid,'AbsInt           = %g\n', fAbsInt);
fprintf(fid,'PowerInt         = %g\n', fPowerInt);
fprintf(fid,'MinSlice         = %g\n', fMinSlice);
fprintf(fid,'MaxSlice         = %g\n', fMaxSlice);
fprintf(fid,'RefGrad          = %g\n', fRefGrad);
fprintf(fid,'\n');


%% gradient section
fprintf(fid,'\n#----------\n');
fprintf(fid,'[Gradient]\t # Gx\t Gy\t Gz\t [mT/m]\n');
fprintf(fid,'\n');
fprintf(fid,'GradientSamples   =  %g\n', iGradientSamples);
fprintf(fid,'MaxAbsGradient[0] =  %g\t %g\t %g\t\t # scaling for G amplitude\n', fMaxAbsG1, fMaxAbsG2, fMaxAbsG3);
fprintf(fid,'GradRasterTime    =  %g   # [us]\n', fGradRasterTime);
fprintf(fid,'\n');
for iG = 1:iGradientSamples
    fprintf(fid,'G[%g]= %g\t %g\t %g\n', iG-1, g(iG,1), g(iG,2), g(iG,3));
end

%% rf section
for iPulse = 1:iNMultiplePulses
    for iCh = 1:iNUsedChannels
        fprintf(fid,'\n#----------\n');
        if( iPulse > 1)
            fprintf(fid,'[pTXPulse_ch%g_%g]\t # Amplitude\t Phase\n', iCh-1, iPulse-1);
        else
            fprintf(fid,'[pTXPulse_ch%g]\t # Amplitude\t Phase\n', iCh-1);
        end
        fprintf(fid,'\n');
        if(~isempty(vfPhaseOffset))
            fprintf(fid,'PhaseOffset       =  %g\n', vfPhaseOffset(iCh,iPulse));
        end
        if(~isempty(vfAmplitudeScale))
            fprintf(fid,'AmplitudeScale    =  %g\n', vfAmplitudeScale(iCh,iPulse));
        end
        for iRF = 1:iSamples
            fprintf(fid,'RF[%g]= %g\t %g\n', iRF-1, abs_rfs(iRF, iCh, iPulse), angle_rfs(iRF, iCh, iPulse));
        end
    end
end
fprintf(fid,'\n');
fprintf(fid,'#EOF\n');
fclose(fid);

if( bVerbose )
%     logprintf(1,'I','writing ''%s'' (%s)', strFileName, VERSION);
end
bSuccess = 1;
return

%--------------------------------------------------------------
function [dopt, narg] = handlerOptionsNarg(defaultOptions, inputOptions, narg, nargLimits)
% 	x = x(dir, filenum, opt())
% 	handler for nargin and function options (optin)
% 	return: 
% 	Functions called: 
% 	Tested for: 
%  	July 2004 -  Josef Pfeuffer

FCTNAME = 'handlerOptionsNarg';

%  option/narg handling 
% --- default options: define whole struct, comment
% defaultOptions

% --- arg handling - error: narg > nargVars; error: narg < nargVars
nargVars		= nargLimits(2);
error(nargchk(0,nargVars+1,narg)) 
if (narg == nargVars+1)
	if prod( size(inputOptions) ) > 1
		logprintf(1,'E','%s: unexpected error<%s>', FCTNAME, size(inputOptions));
	end
	optin	= inputOptions{1};
   	dopt 	= LocalSetopt(defaultOptions, optin);
   	narg 	= narg - 1;    % narg is NOT including options
else
	dopt 	= defaultOptions;
    %logprintf(1,'I','%s: check nargLimits [%d %d]\n', FCTNAME, nargLimits(1), nargLimits(2));
end

return

%------------------------------------------------------------
function optout = LocalSetopt(optall, optin)
%
%   updates default optall struct with input optin
%   optout = setopt(optall, optin)
%
%   effects: warning on unkown option field
%            ignores tagfield
%
%   JP Apr 2000

narg = nargin;
msg = nargchk(2,2,narg);
if ~isempty(msg)
    logprintf(1,'E',msg);
end

tagfield = 'STRUCTNAME';
tagstr = 'optstruct';

optout = optall;
if (isempty(optin))
   return
end
descin = fieldnames(optin);
l_descin = length(descin);

for il=1:l_descin
   if ( strcmp(descin(il), tagfield) &&  ...
	strcmp(getfield(optin,descin{il}), tagstr) )
      % do nothing
      %logprintf(1,'I','field <%s> contains <%s>',descin{il}, ...
      %            getfield(optin,descin{il});
   elseif (  isfield(optall, descin(il)) )
      optout = setfield(optout, descin{il}, getfield(optin,descin{il}));
   else
      logprintf(1,'W','option <%s> ignored',descin{il});
   end
end % for

%--------------------------------------------------------------
function s = opt(varargin)
%
%   builds struct from arguments
%   optstruct = opt('OPTION1',var1,'OPTION2',var2, ...)
%
%   effects: fielddescriptors are converted to UPPERCASE
%
%   JP Apr 2000

tagfield = 'STRUCTNAME';
tagstr = 'optstruct';

len_arg = length(varargin);
s = [];
s = setfield(s,tagfield,tagstr);      % init s - return to call w/o argument
for ilen = 1:2:len_arg-1
   if (~ischar(varargin{ilen}))
      varargin;
      logprintf(1,'E','opt(): descriptor is not a string !');
   end
   fielddescriptor = upper(num2str( varargin{ilen} ));
   s = setfield(s, fielddescriptor, varargin{ilen+1});
end

%------------------------------------------------------------
function demoUsage(FCTNAME)
% %
fprintf('  %%\n  %%\n  %%   %s: entering DEMO mode\n  %%\n  %%\n',FCTNAME)

% demo data for testing
load save_pTXRFPulse_toINI_Data_090627
whos g rf
strComm = '# test pulse for debugging in matlab';

save_pTXRFPulse_toINI( g, rf, [], opt('FACTOROVERSAMPLE',1,'NOMFLIPANGLE',90,'RFPULSE_COMMENT',strComm));

shim_coeff = complex( zeros(1,8));   % 8 channel
aRFShimPhase = 0.111*(1:8);   % phase offset
aRFShimAmp   = 1.11 *(1:8);   % amplitude
shim_coeff(1,:) = aRFShimAmp.*exp(1i*aRFShimPhase/180*pi);
g_low =[];%zeros(1,3);
b_est =[];%zeros(8,1);
save_pTXRFPulse_toINI( g_low, b_est, shim_coeff, opt('NOMFLIPANGLE',90,'RFPULSE_ID',1,'RFPULSE_COMMENT',strComm));

% %load('singleINI.mat');
% load('multipleINI');
% dopt.FACTOROVERSAMPLE = Dpa.oversample;
% dopt.NOMFLIPANGLE     = Dpa.alpha; 
% dopt.ASYMMETRY        = Dpa.pulse_asymmetry;
% dopt.INITIALPHASE     = Dpa.pulse_initialPhase;
% whos g_out b_est shim_coeff
% save_pTXRFPulse_toINI(g_out * 10, b_est, shim_coeff, opt('NOMFLIPANGLE',90,'RFPULSE_ID',0));  % convert g_out [G/cm] -> [mT/m]

%------------------------------------------------------------
