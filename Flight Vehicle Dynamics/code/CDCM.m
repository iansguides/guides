% Author: Ian Harris
% Filename: CDCM.m
% Date: 02/02/2023
%
% INPUTS
% ------------------------------------------------------------------------- 
% phi: input name - phi (in degrees)
%      format type - double
%      length required - 1  
%
% theta: input name - phi (in degrees)
%        format type - double
%        length required - 1 
%
% psi: input name - phi (in degrees)
%      format type - double
%      length required - 1 
%
% sequence: input name - rotation matrix sequence
%           format type - double array
%           minimum length - 1
%
% 
% OUTPUTS
% -------------------------------------------------------------------------
% DCM: input name - DCM
%      format type - matrix double
%      size - 3x3  
%
%
% REQUIREMENTS
% -------------------------------------------------------------------------
% This function requires that rotx.m, roty.m, and rotz.m must be in the 
% same path to function properly.
%
%
% PURPOSE
% -------------------------------------------------------------------------
% This code allows for quick create of DCMs. The CDCM stands for calculate
% directional cosine matrix.

function [DCM] = CDCM(phi, theta, psi, sequence) 

    % Check if sequence input is defined. If not, then default to [3,2,1].
    if nargin == 3
        sequence = [3,2,1];
    end

    % Check to see if sequence contains any numbers other than 1, 2, or 3.
    if length(sequence) ~= (sum(sequence == 1) + sum(sequence == 2) + sum(sequence == 3))
        cprintf("red","Error: sequence must only contain 1, 2, or 3.\n");
        DCM = "Error";
        return
    end

    % Predefine the DCM matrix.
    DCM = 1;

    for i = 1:length(sequence)
        if sequence(i) == 1
            DCM = DCM * rotx(phi);
        elseif sequence(i) == 2
            DCM = DCM * roty(theta);
        else
            DCM = DCM * rotz(psi);
        end
    end
        
end