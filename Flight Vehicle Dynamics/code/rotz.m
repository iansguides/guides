% Author: Ian Harris
% Filename: rotz.m
% Date: 1/27/23
% Version: 2.0.0
% 
% INPUTS
% ------------------------------------------------------------------------- 
% deg: input name - degree change (in deg or rad)
%      format type - double
%      default unit - degrees
%
% type: input name - unit type (either deg or rad allowed)
%       format type - string
%       optional input
%
% 
% OUTPUTS
% -------------------------------------------------------------------------
% Rx: output name - z rotation matrix (unitless)
%     format type - matrix
%     length - 3x3
%
%
% PURPOSE
% -------------------------------------------------------------------------
% The purpose of this function is to form the z rotation matrix for a CCW 
% rotation.


function [Rz] = rotz(ang, type)

    % Check if unit type has been provided or not. If it is not provided, 
    % default to degree.
    if nargin == 1
        type = "deg";
    end
    
    % Check to see if type is equal to either rad or deg.
    if nargin == 2
        if (type ~= "deg") && (type ~= "rad")
            cprintf("red","Error: Unit type must be either ""deg"" or ""rad"".\n");
            Rz = "Error";
            return
        end
    end

    % If the degree type is equal to degree, convert it to radian.
    if type == "deg"
        ang = ang*pi/180;
    end

    % Form the z rotation matrix.
    Rz = [  cos(ang)  sin(ang)  0     
           -sin(ang)  cos(ang)  0 
            0          0        1  ];  
    
end