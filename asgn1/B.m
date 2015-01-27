function [ x ] = B( psi, r, d, k, h )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

x = (r-d)*psi*(1-psi)*k/h;

end

