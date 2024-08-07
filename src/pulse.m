%%% pulse.m --- A superclass for pulse stuff, like gradient and RF pulses.
%% 
%% Filename: pulse.m
%% Description: 
%% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
%% Maintainer: 
%% Copyright (C) 2009 CMRR at UMN
%% Created: Thu Sep  3 09:08:39 2009 (CDT)
%% Version: 
%% Last-Updated: Fri Feb  5 12:31:57 2010 (CST)
%%           By: Xiaoping Wu
%%     Update #: 27
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Commentary: 
%% 
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

classdef pulse
    properties 
      Pulse
      TimeStep
      Name = ''
      
    end % properties

    properties (Dependent = true, SetAccess = private)
      NumOfPulses
      NumOfTimePoints
      Duration
    end
 % -----------------------

    methods 
        function obj = pulse (puls,dt,name)
        %
        if nargin > 0
          obj.Pulse = puls;
          obj.TimeStep = dt;
          obj.Name = name;end
          
        end
        % -----------------------

        function npts = get.NumOfTimePoints (obj)
          npts = size(obj.Pulse,2);
        end
        
        function npuls = get.NumOfPulses (obj)
          npuls = size(obj.Pulse,1);
        end
        
        function duration = get.Duration (obj)
          duration = obj.TimeStep * obj.NumOfTimePoints;
        end
        
        
        function disp (obj)
          fprintf(1,'Name: %s\nNumOfPulses:%d\nDuration:%6.2f ms\n',obj.Name, ...
                  obj.NumOfPulses,obj.Duration* 1e3);
        end
        
        function write(obj)
        % write to a text file
          % for ind = 1: obj.NumOfDimensions,
%           % writing ...
%           end          
        disp('-> Not implemented yet...')
        end

        function puls = double(obj)
          puls = obj.Pulse;
        end
    end % methods
    
    methods (Access = protected)
        function timeArray = createTimeArray(obj)          
          timeArray = (1:obj.NumOfTimePoints) * obj.TimeStep;
        end
    end % methods
end % classdef


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pulse.m ends here
