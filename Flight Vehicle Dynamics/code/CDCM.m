% Author: Ian Harris
% Filename: CDCM.m
% Date: 01/20/2024
% Version: 3.0.0
%
% INPUTS
% ------------------------------------------------------------------------- 
% rot1: input name - rotation 1
%       format type - double
%       length required - 1  
%
% rot2: input name - rotation 2
%       format type - double
%       length required - 1 
%
% rot3: input name - rotation 3
%       format type - double
%       length required - 1 
%
% sequence: input name - rotation matrix sequence
%           format type - double array
%           length required - 3
% 
% type: input name - unit type (either deg or rad allowed)
%       format type - string
%       optional input
%
% 
% OUTPUTS
% -------------------------------------------------------------------------
% DCM: output name - DCM
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

function [DCM] = CDCM(rot1, rot2, rot3, sequence, type) 

    % Check if sequence input is defined. If not, then default to [3,2,1].
    % If sequence is not defined, then type is not defined either. Default
    % is "deg".
    if nargin == 3
        sequence = [3,2,1];
        type = "deg";
    end

    % Check if type is defined. If not, default is "deg".
    if nargin == 4
        type = "deg";
    end

    % Check to see if sequence contains any numbers other than 1, 2, or 3.
    if length(sequence) ~= (sum(sequence == 1) + sum(sequence == 2) + sum(sequence == 3))
        cprintf("red","Error: sequence must only contain 1, 2, or 3.\n");
        DCM = "Error";
        return
    end

    % Check if length of sequence is equal to 3.
    if length(sequence) ~= 3
        cprintf("red","Error: sequence must be an array of length 3.\n");
        DCM = "Error";
        return
    end

    % Check if lengths of rot1, rot2, and rot3 are equal to 1.
    if (length(rot1) ~= 1) || (length(rot2) ~= 1) || (length(rot3) ~= 1)
        cprintf("red","Error: rot1, rot2, and rot3 must have a length of 1.\n");
        DCM = "Error";
        return
    end

    % Check if type is equal to either "rad" or "deg".
    if (type ~= "deg") && (type ~= "rad")
        cprintf("red","Error: Unit type must be either ""deg"" or ""rad"".\n");
        DCM = "Error";
        return
    end

    % Predefine the DCM matrix.
    DCM = 1;

    % Predefine the rotation array.
    rot_arr = [rot3 rot2 rot1];
    sequence = fliplr(sequence);

    % Calculate the DCM and change the unit type of each rot function to
    % the one specified in CDCM.
    if type == "deg"
        for i = 1:length(sequence)
            if sequence(i) == 1
                DCM = DCM * rotx(rot_arr(i));
            elseif sequence(i) == 2
                DCM = DCM * roty(rot_arr(i));
            else
                DCM = DCM * rotz(rot_arr(i));
            end
        end
    else % type == "rad"
        for i = 1:length(sequence)
            if sequence(i) == 1
                DCM = DCM * rotx(rot_arr(i),"rad");
            elseif sequence(i) == 2
                DCM = DCM * roty(rot_arr(i),"rad");
            else
                DCM = DCM * rotz(rot_arr(i),"rad");
            end
        end
    end
        
end