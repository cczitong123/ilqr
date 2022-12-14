function [ result ] = ilqr_single( f, j, dt, N, x0, u0, opt_param )
%ILQR_SINGLE Summary of this function goes here
%   Detailed explanation goes here
    result = ILQRController.ilqr(f, j, dt, N, x0, u0, opt_param);

end

