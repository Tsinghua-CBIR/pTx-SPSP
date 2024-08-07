%%% gradPulse.m --- A class for gradient pulses
%% 
%% Filename: gradPulse.m
%% Description: 
%% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
%% Maintainer: 
%% Copyright (C) 2009 CMRR at UMN
%% Created: Thu Sep  3 11:18:56 2009 (CDT)
%% Version: 
%% Last-Updated: Mon Jan 24 14:53:59 2011 (CST)
%%           By: Xiaoping Wu
%%     Update #: 60
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Commentary: 
%% 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Change log:
%% 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Code:

classdef gradPulse < pulse
    properties 
      Type = 'Tx' % for transmit
    end % properties

    properties (Dependent = true, SetAccess = private)
      KspaceTrajectory
      SlewRate
    end
 % -----------------------

    methods 
        function obj = gradPulse (grad,dt,name,type)
        % Constructor
        % Usage: obj = gradPulse (grad,dt,name,type)
        % grad: gradient data in tesla/m. defaults to [].
        % dt: dwell time in s. defaults to 4e-6.
        % name: gradient name, e.g., 'spiral'. defaults to ''.
        % type: can only be 'Tx' for transmit, or 'Rx' for receive. defaults to 'Tx'.
        
      
          if nargin < 3
            name = '';
          end
          
          if nargin < 2
            dt = 4e-6;
          end
          
          if nargin < 1
            grad = [];
          end
          
          args{1}= grad;
          args{2}= dt;
          args{3}= name;
          
          
          obj = obj@pulse(args{:});
          
          if nargin > 3
            obj.Type = type;
          end
          
          
        end
        % -----------------------

        function obj = set.Type(obj,type)
          if ~(strcmpi(type,'Tx')||strcmpi(type,'Rx'))
            error('Type must be Tx or Rx')
          end
        
          obj.Type = type;
        end 
        
        function kspTraj = get.KspaceTrajectory (obj)
          grad= obj.Pulse;
          dt = obj.TimeStep;
          
          mymethod = 'c';                      
          if strcmpi(obj.Type,'Tx')
            kspTraj = calc_ktraj_from_grad(grad,true,dt,mymethod);
          else                                  % reception
            kspTraj = calc_ktraj_from_grad(grad,false,dt,mymethod);
          end
                      
        end
        
        function slewRate = get.SlewRate (obj)
        
        slewRate = diff(obj.Pulse,1,2)/obj.TimeStep;
        slewRate = [zeros(obj.NumOfPulses,1) slewRate];
          
        end

        
        function disp (obj)
          disp@pulse(obj);
          fprintf(1,'Type: %s\n',obj.Type);
        end
        
        function plot (obj, varargin)
          time = obj.createTimeArray;
          
          plot(1e3*time',1e3*obj.Pulse',varargin{:})
          title('Gradient Waveform')
          xlabel('Time (ms)')
          ylabel('Amp (mT/m)')
        end

        function plotSlewRate (obj,varargin)
          time = obj.createTimeArray;
          plot(1e3*time',obj.SlewRate',varargin{:})
          title('Slew Rate')
          xlabel('Time (ms)')
          ylabel('Slew rate (T/m/s)')
        end

        function plotTrajWaveform (obj, varargin)
          time = obj.createTimeArray;
          
          plot(1e3*time', obj.KspaceTrajectory', varargin{:})
          title('Kspace Trajectory Waveform')
          xlabel('Time (ms)')
          ylabel('Amp (rad/m)')
        end        
        
        function plotTraj (obj, varargin)
          
          switch obj.NumOfPulses
           case 1
            disp('Plotting 1D kspace is trivial, please refer to plotTxTraj.')
           case 2
            plot(obj.KspaceTrajectory(1,:),obj.KspaceTrajectory(2,:), ...
                 varargin{:});
           case 3
            plot3(obj.KspaceTrajectory(1,:),obj.KspaceTrajectory(2,:), ...
                  obj.KspaceTrajectory(3,:),varargin{:});
           otherwise
          end
          title('Kspace Trajectory')
            xlabel('Kx (rad/m)')
            ylabel('Ky (rad/m)')
            zlabel('Kz (rad/m)')
            axis square
        end
        
        function write(obj)
        % write gradient shape to a text file for a specified system, e.g. Varian.
          
          for ind = 1: obj.NumOfPulses,
            write_grad_varian(obj.Pulse(ind,:),[obj.Name '_' num2str(ind)]);
          end
                    
        end
                
    end % methods
end % classdef


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% gradPulse.m ends here
