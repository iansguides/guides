% Author: Ian Harris
% Filename: invDCM.m
% Date: 01/20/2024
% Version: 3.0.0
%
% INPUTS
% ------------------------------------------------------------------------- 
% DCM: input name - DCM (unitless)
%      format type - double matrix
%      length required - 3x3 
%
% sequence: input name - rotation matrix sequence
%           format type - double array
%           length required - 3
% 
% tolerance: input name - DCM error tolerance
%            format type - double
%            optional input
%
% debug: input name - debug mode
%        format type - string
%        optional input
% 
% type: input name - output unit type (either deg or rad allowed)
%       format type - string
%       optional input
%
% 
% OUTPUTS
% -------------------------------------------------------------------------
% angles: output name - rotation angles
%         format type - double matrix 
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
% This code allows for one to backsolve DCMs to find the rotation angles
% used to create a corresponding DCM.

function [angles] = invDCM(DCM, sequence, tolerance, debug, type)
    
    % Check to see if type is defined. If it is, use that output unit type;
    % if not, default to degree. Also if the debug mode is not given,
    % default to "false". Also if tolerance is not defined, default to
    % 0.001.
    if nargin == 2
        tolerance = 0.001;
        debug = "false";
        type = "deg";
    end

    % If the type was not defined, default to degree. If debug was not
    % defined, default to false.
    if nargin == 3
        debug = "false";
        type = "deg";
    end

     % If the type was not defined, default to degree.
    if nargin == 4
        type = "deg";
    end

    % Check to make sure length of sequence is equal to 3.
    if length(sequence) ~= 3
        cprintf("red","Error: sequence must have length of 3.\n");
        angles = "Error";
        return
    end

    % Check to see if sequence contains any numbers other than 1, 2, or 3.
    if length(sequence) ~= (sum(sequence == 1) + sum(sequence == 2) + sum(sequence == 3))
        cprintf("red","Error: sequence must only contain 1, 2, or 3.\n");
        angles = "Error";
        return
    end

    % Check the size of the DCM to make sure it is a 3x3.
    if (size(DCM,1) ~= 3) || (size(DCM,2) ~= 3)
        cprintf("red","Error: DCM must have size of 3x3.\n");
        angles = "Error";
        return
    end

    % Check to see if debug is either "true" or "false".
    if (debug ~= "true") && (debug ~= "false")
        cprintf("red","Error: debug must be either ""true"" or ""false"".\n");
        angles = "Error";
        return
    end

    % Checks to see if type is equal to "deg" or "rad".
    if (type ~= "deg") && (type ~= "rad")
        cprintf("red","Error: type must be either ""rad"" or ""deg"".\n");
        angles = "Error";
        return
    end

    % Check that tolerance is a single value.
    if ~(size(tolerance,1)==1 && size(tolerance,2)==1)
        cprintf("red","Error: tolerance must be a single value, not an array.\n");
        angles = "Error";
        return
    end

    % Check that tolerance is a numeric value.
    if ~isnumeric(tolerance)
        cprintf("red","Error: tolerance must be a numeric value.\n");
        angles = "Error";
        return
    end

    % Define symbolic variables representing rotation angles 1, 2, and 3.
    syms rot1_sym rot2_sym rot3_sym;
    assume(rot1_sym, "real");
    assume(rot2_sym, "real");
    assume(rot3_sym, "real");
    
    % Form array of symbolic rotation variables and predefine the symbolic
    % DCM to map to.
    sym_vars = [rot3_sym, rot2_sym, rot1_sym];
    symDCM = 1;
    sequence = fliplr(sequence);

    % Loop through the sequence and form the symbolic DCM using the
    % corresponding rotations. Pull the rotations from the sym_vars array
    % as the inputs for the rotx(), roty(), and rotz() functions.
    for i = 1:length(sequence)
        if sequence(i) == 1
            symDCM = symDCM * rotx(sym_vars(i),'rad');
        elseif sequence(i) == 2
            symDCM = symDCM * roty(sym_vars(i),'rad');
        elseif sequence(i) == 3
            symDCM = symDCM * rotz(sym_vars(i),'rad');
        end
    end

    % Predefine count matrix. This will contain the number of cosines in
    % the first column and the number of sines in the second column for
    % each index of the symbolic matrix. The number of rows will start out
    % at 9, counting down each row as the MATLAB index counter would for
    % matrix layouts. Also define a 3x3 cell to hold all symbolic variables
    % in each term of the matrix as well as a 3x3 cell to hold the count
    % for each cosine and sine.
    count_mat = [];
    count_cell = cell(3,3);
    sym_var_mat = cell(3,3);
    for i = 1:3
        for j = 1:3
            % Convert the symbolic expression to a string.
            exprString = char(symDCM(i, j));
    
            % Count the number of cosine and sine terms in the string.
            countCos = length(strfind(exprString, 'cos'));
            countSin = length(strfind(exprString, 'sin'));
    
            % Store the counts in the count_mat as well as count_cell.
            count_mat(end+1,1:2) = [countCos, countSin];
            count_cell{i,j} = [countCos, countSin];

            % Add the symbolic variables in each index into the sym_var_mat
            % for later tracking.
            sym_var_mat{i,j} = symvar(symDCM(i,j));
        end
    end

    % Sum the count_mat matrix row-wise to determine which index terms are
    % the most and least complex.
    count_mat_sum = sum(count_mat,2);

    % Reshape the count matrix into a 3x3 to mimic the same position count
    % for each corresponding index.
    count_mat_sum = reshape(count_mat_sum,3,3)';

    % Each 3-rotation DCM will have a single term by itself; we target that
    % here by finding where the sin/cos count is equal to 1. This will be
    % the first equation in our system to solve becauase of its inherit
    % simplicity.
    [i_min,j_min] = ind2sub(size(count_mat_sum), ...
                            find(abs(count_mat_sum)==1));

    % Find the symbolic variable used in the simplest index.
    sym_var_used{1} = sym_var_mat{i_min,j_min};

    % Find all second simplest terms in the DCM. These are the ones where
    % they have a cos and/or sin multiplied together.
    [i_sec,j_sec] = ind2sub(size(count_mat_sum), ...
                            find(abs(count_mat_sum)==2));

    % Find two other indices that include all three rotation variables.
    % There can be overlap, but each index must include a new rotation
    % variable from before. Predefine the counter and an array to hold the
    % index values from the i_sec and j_sec index arrays.
    s = 1;
    st_inds = [];

    % Continue a while loop until all 3 rotation angles are found.
    while length(sym_var_used) < 3
        
        % Get the temporary array to analyze.
        temp_arr = sym_var_mat{i_sec(s),j_sec(s)};
        % First check to see if the simplest term's symbolic variable is in
        % the current, selected index.
        if ismember(sym_var_used{1},temp_arr)
            % Loop through each variable in the array individually; there
            % should be 2 values per each temporary array.
            for r = 1:length(temp_arr)
                % If the sum of the ismember function returns greater than
                % zero, then there is a repeat variable. If equal to zero, then
                % there is a new variable. Save it to the sym_var_used cell
                % array and also save the index for later eq use.
                if sum(ismember(sym_var_used,temp_arr(r))) > 0 
                    continue
                else
                    sym_var_used{end+1} = temp_arr(r);
                    st_inds(end+1) = s;
                    break
                end
            end
        end
        s = s + 1;
    end

    % Get the count of cos and sin in the simplest index.
    trig_count = count_cell{i_min,j_min};
    cos_count = trig_count(1);
    sin_count = trig_count(2);

    % Get the correspding DCM value, symbolic expression, and determine if
    % the symbolic expression has a negative.
    val = DCM(i_min,j_min);
    symInd = symDCM(i_min,j_min);
    neg = length(strfind(char(symInd),'-'));

    % First check if the simplest (single-term) index of the symbolic
    % matrix is sin or cos. Then check if it is negative. If so, add a term
    % of negative 1 when computing the 2nd rot angle. Then,
    % find the alternative angle based off sin or cos for the other
    % possible quadrant.
    if cos_count
        if ~neg
            rot2 = acosd(val);
        else
            rot2 = acosd(-val);
        end
    
        if (rot2 >= 0) && (rot2 <= 90) % in quadrant 1
            rot2_alt = -rot2;
        else % rot2 is in quadrant 2
            rot2_alt = 360-rot2;
        end
    else % sin_count == 1
        if ~neg
            rot2 = asind(val);
        else
            rot2 = asind(-val);
        end
    
        rot2_alt = 180-rot2;
    end

    % Group the two potential angles for the 2nd rotation.
    rot2_arr = [rot2, rot2_alt];

    % Define a temporary array using the first of the two second simplest
    % values.
    temp_arr = sym_var_mat{i_sec(st_inds(1)),j_sec(st_inds(1))};

    % Check to see if the first rotation sym var is present in the temp
    % array. If it is not, switch the index order so rot1 is dealt with
    % first and rot3 is dealt with last.
    if ~ismember(rot1_sym,temp_arr)
        temp_ind = st_inds(1);
        st_inds(1) = st_inds(2);
        st_inds(2) = temp_ind;
    end
     
    % Initialize a matrix to hold the 4 rotation 1 and rotation 3 values.
    rot13_vals = zeros(2,4);

    % Loop through the 2 index values corresponding to i_sec and j_sec to
    % determine the potential angles for each index's new rotation angle.
    for h = 1:length(st_inds)

        % Get the focus symbolic expression, the cos and sin count as well
        % as if the expression is negative or not.
        symInd2 = symDCM(i_sec(st_inds(h)),j_sec(st_inds(h)));
        cos_count = count_cell{i_sec(st_inds(h)),j_sec(st_inds(h))}(1);
        sin_count = count_cell{i_sec(st_inds(h)),j_sec(st_inds(h))}(2);
        neg = length(strfind(char(symInd2),'-'));
    
        % First check if rot2_sym is in a cos term or a sin term. This will
        % alter how the rot1 and rot3 angles are calculated. Then calculate
        % the 2 potential angles for each as well as their alternatives.
        if length(strfind(char(symInd2),'cos(rot2_sym)'))
            if cos_count == 2
                if ~neg
                    rot13_vals(h,1) = acosd(DCM(i_sec(st_inds(h)),j_sec(st_inds(h)))/cosd(rot2));
                    rot13_vals(h,3) = acosd(DCM(i_sec(st_inds(h)),j_sec(st_inds(h)))/cosd(rot2_alt));
                else
                    rot13_vals(h,1) = acosd(DCM(i_sec(st_inds(h)),j_sec(st_inds(h)))/-cosd(rot2));
                    rot13_vals(h,3) = acosd(DCM(i_sec(st_inds(h)),j_sec(st_inds(h)))/-cosd(rot2_alt));
                end
        
                if (rot13_vals(h,1) >= 0) && (rot13_vals(h,1) <= 90) % in quadrant 1
                    rot13_vals(h,2) = -rot13_vals(h,1);
                else % rot2 is in quadrant 2
                    rot13_vals(h,2) = 360-rot13_vals(h,1);
                end
        
                if (rot13_vals(h,3) >= 0) && (rot13_vals(h,3) <= 90) % in quadrant 1
                    rot13_vals(h,4) = -rot13_vals(h,3);
                else % rot2 is in quadrant 2
                    rot13_vals(h,4) = 360-rot13_vals(h,3);
                end
            else % cos count is 1
                if ~neg
                    rot13_vals(h,1) = asind(DCM(i_sec(st_inds(h)),j_sec(st_inds(h)))/cosd(rot2));
                    rot13_vals(h,3) = asind(DCM(i_sec(st_inds(h)),j_sec(st_inds(h)))/cosd(rot2_alt));
                else
                    rot13_vals(h,1) = asind(DCM(i_sec(st_inds(h)),j_sec(st_inds(h)))/-cosd(rot2));
                    rot13_vals(h,3) = asind(DCM(i_sec(st_inds(h)),j_sec(st_inds(h)))/-cosd(rot2_alt));
                end
        
                rot13_vals(h,2) = 180-rot13_vals(h,1);
                rot13_vals(h,4) = 180-rot13_vals(h,3);
            end
        else % index has sin(rot2_sym)
            if sin_count == 2
                if ~neg
                    rot13_vals(h,1) = asind(DCM(i_sec(st_inds(h)),j_sec(st_inds(h)))/sind(rot2));
                    rot13_vals(h,3) = asind(DCM(i_sec(st_inds(h)),j_sec(st_inds(h)))/sind(rot2_alt));
                else
                    rot13_vals(h,1) = asind(DCM(i_sec(st_inds(h)),j_sec(st_inds(h)))/-sind(rot2));
                    rot13_vals(h,3) = asind(DCM(i_sec(st_inds(h)),j_sec(st_inds(h)))/-sind(rot2_alt));
                end
        
                rot13_vals(h,2) = 180-rot13_vals(h,1);
                rot13_vals(h,4) = 180-rot13_vals(h,3);
            else % sin count is 1
                if ~neg
                    rot13_vals(h,1) = acosd(DCM(i_sec(st_inds(h)),j_sec(st_inds(h)))/sind(rot2));
                    rot13_vals(h,3) = acosd(DCM(i_sec(st_inds(h)),j_sec(st_inds(h)))/sind(rot2_alt));
                else
                    rot13_vals(h,1) = acosd(DCM(i_sec(st_inds(h)),j_sec(st_inds(h)))/-sind(rot2));
                    rot13_vals(h,3) = acosd(DCM(i_sec(st_inds(h)),j_sec(st_inds(h)))/-sind(rot2_alt));
                end
        
                if (rot13_vals(h,1) >= 0) && (rot13_vals(h,1) <= 90) % in quadrant 1
                    rot13_vals(h,2) = -rot13_vals(h,1);
                else % rot2 is in quadrant 2
                    rot13_vals(h,2) = 360-rot13_vals(h,1);
                end
        
                if (rot13_vals(h,3) >= 0) && (rot13_vals(h,3) <= 90) % in quadrant 1
                    rot13_vals(h,4) = -rot13_vals(h,3);
                else % rot2 is in quadrant 2
                    rot13_vals(h,4) = 360-rot13_vals(h,3);
                end
            end
        end
    end

    % Extract all of the potential rot1 and rot3 angles.
    rot1_arr = rot13_vals(1,:);
    rot3_arr = rot13_vals(2,:);

    % Predefine an array for dcm error norms as well as angle collections.
    dcm_norms = [];
    angles_collec = {};

    % Next, loop through each potential angle combination, save the angles
    % to angles_collec{}, form the numeric DCM, find its error with the
    % original, and save the norm of the dcm error to the dcm_norms().
    for i = 1:length(rot2_arr)
        for j = 1:length(rot1_arr)
            for k = 1:length(rot3_arr)
                rot_arr = [rot3_arr(j),rot2_arr(i),rot1_arr(k)];
                num_DCM = 1;
                for z = 1:length(sequence)
                    if sequence(z) == 1
                        num_DCM = num_DCM * rotx(rot_arr(z));
                    elseif sequence(z) == 2
                        num_DCM = num_DCM * roty(rot_arr(z));
                    elseif sequence(z) == 3
                        num_DCM = num_DCM * rotz(rot_arr(z));
                    end
                end
                
                DCM_error_norm = norm(DCM - num_DCM);
                dcm_norms(end+1) = DCM_error_norm;
                angles_collec{end+1} = rot_arr;
            end
        end
    end

    % Find all angle sets within a norm tolerance.
    idxs = find(dcm_norms < tolerance);
    
    % Predefine the output angles matrix. Each row will correspond to a set
    % and each column will be a rotation going from rot1 to rot3.
    angles = [];
    for i = 1:length(idxs)
        angles(i,1:3) = angles_collec{idxs(i)};
    end

    angles = fliplr(angles);

    % Check to see what output unit the angles are in, then output that
    % unit. Also use mod() function to get them within the bounds of one
    % positive revolution.
    if type == "rad"
        angles = mod(angles*pi/180,2*pi);
    else
        angles = mod(angles,360);
    end

    if debug == "true"
        % Gather important values.
        min = [i_min,j_min];
        sec1 = [i_sec(st_inds(1)),j_sec(st_inds(1))];
        sec2 = [i_sec(st_inds(2)),j_sec(st_inds(2))];

        % Display values.
        disp("DCM");
        disp(DCM);
        disp("symDCM");
        disp(symDCM);
        disp("min");
        disp(min);
        disp("sec1");
        disp(sec1);
        disp("sec2");
        disp(sec2);
    end

end