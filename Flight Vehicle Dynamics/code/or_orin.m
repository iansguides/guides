% Author: Ian Harris
% Date: 1/19/23
% 
% INPUTS
% ------------------------------------------------------------------------- 
% del: input name - delta (E-W) (in degrees)
%      format type - array
%      length required - 2     
% lam: input name - lambda (N-S) (in degrees)
%      format type - array
%      length required - 2
%
% NOTE: N and E are positive entries and S and W are negative entries
%
% 
% OUTPUTS
% -------------------------------------------------------------------------
% i: output name - inclination (in degrees)
%    format type - double
%    length - 1
% raan: output name - RAAN (in degrees)
%       format type - double
%       length - 1
%
%
% PURPOSE
% -------------------------------------------------------------------------
% The purpose of this function is to take delat and lambda values for
% two object observations and in return, the inclination and RAAN is
% returned for the object in degrees.


function [i, raan] = or_orin(del, lam)

    % Validate user entry. If user entry is not equal to 2 for both arrays,
    % send error message to user with inputed array lengths. If error
    % message is shown, the fn ends.
    if (length(del) ~= 2) || (length(lam) ~= 2)
        cprintf("_red","Error: Array lengths are not equal to 2.\n");
        cprintf("red","    del: length=%i\n", length(del));
        cprintf("red","    lam: length=%i\n", length(lam));
        i = "Error";
        raan = "Error";
        return
    end

    % Predefine velocity matricies and calculate each velocity matrix per
    % data entry group. Should be only two velocities --> one for each
    % observation.
    v_mats = {};
    for a = 1:2
        v_mats{end+1} = [ cosd(lam(a))*cosd(del(a))
                          cosd(lam(a))*sind(del(a))
                          sind(lam(a)) ];
    end

    % Predefine Earth Centric Inertial matrix in i-dir as well as k-dir.
    i_eci = [1 0 0];
    k_eci = [0 0 1];

    % Calculate normal vector to obseration plane using cross() fn.
    n_or = cross(v_mats{1},v_mats{2});

    % Calculate axis of nodes using k_eci and n_or.
    n_an = cross(k_eci, n_or);

    % Compute the inclination of the observation plane.
    i = acosd(dot(n_or,k_eci) / norm(n_or));

    % Compute RAAN of the observation plane.
    raan = atand(n_an(2) / n_an(1));

end