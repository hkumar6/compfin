function [ x ] = A( psi, k, h, sig )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x = (sig^2)*(psi^2)*(1-psi)^2*k/(h^2);

end

