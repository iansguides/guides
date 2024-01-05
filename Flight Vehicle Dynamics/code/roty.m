% Author: Ian Harris
% Filename: roty.m
% Date: 1/27/23
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
% Ry: output name - y rotation matrix (unitless)
%     format type - matrix
%     length - 3x3
%
%
% PURPOSE
% -------------------------------------------------------------------------
% The purpose of this function is to form the y rotation matrix.


function [Ry] = roty(ang, type)

    % Check if unit type has been provided or not. If it is not provided, 
    % default to degree.
    if nargin == 1
        type = "deg";
    end
    
    % Check to see if type is equal to either rad or deg.
    if nargin == 2
        if (type ~= "deg") && (type ~= "rad")
            cprintf("red","Error: Unit type must be either ""deg"" or ""rad"".\n");
            Rx = "Error";
            return
        end
    end

    % If the degree type is equal to radian, convert it to degree.
    if type == "rad"
        ang = ang*180/pi;
    end

    % Form the x rotation matrix.
    Ry = [  cosd(ang)   0  sind(ang) 
            0           1  0
            -sind(ang)  0  cosd(ang)  ];  
    
end