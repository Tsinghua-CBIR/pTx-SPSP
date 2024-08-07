%%% matrix.m --- A class of matrix
%% 
%% Filename: matrix.m
%% Description: 
%% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
%% Maintainer: 
%% Copyright (C) 2009 CMRR at UMN
%% Created: Thu Sep 17 09:15:33 2009 (CDT)
%% Version: 
%% Last-Updated: Thu Jul  2 13:35:51 2015 (CDT)
%%           By: Xiaoping Wu
%%     Update #: 30
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

classdef matrix
    properties 
      Matrix
    end % properties
 
    properties (SetAccess = protected)
      UMatrix
      SingularValueVector
      VMatrix
    end % properties

    methods 
        function obj = matrix (amatrix)
        % Constructor
        % Usage: obj = matrix (amatrix)
        % Matrix: m-by-n matrix data
        
        % Pre init.
        % anything not using output (obj)

          if nargin < 1
            amatrix = [];
          end
          
        % compvalue = classname.staticMethod();

        
        % Object init.
        % call super class, if applicable, before accessing object

        % obj = obj@superClass(args{:});
          obj.Matrix = amatrix;
        
        % Post init.
        % anything including accessing object

        % obj.classMeth
        % obj.Property = compvalue;
        end
       
        function obj = set.Matrix (obj, amatrix)
          obj.Matrix = amatrix;
          [obj.UMatrix,obj.SingularValueVector,obj.VMatrix] = obj.svd;
        end
        
        function disp (obj)
          [m,n] = size(obj.Matrix);
          fprintf(1,'m: %d\nn: %d\n',m,n);
        end
        
        function [U,s,V] = getSVD (obj)
         U = obj.UMatrix;
         s = obj.SingularValueVector;
         V = obj.VMatrix;
        end

    end % methods
    
    methods (Access = protected)
    
        function [U,s,V] = svd (obj)
        % Purpose: singular value decomposition
         if isempty(obj.Matrix)
           U = [];
           s = [];
           V = [];
           return
         end
         
%         disp('-> SVDing ...')
         
         [U,s,V] = csvd(obj.Matrix);
         
%         disp('-> SVD done.')
        end

    end % methods

end % classdef



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% matrix.m ends here
