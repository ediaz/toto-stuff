function [ L ] = LaplacianOperator(h,n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


[Dx,Dz]= DifferenceOperators(h,n);

L = -Dx'*Dx -Dz'*Dz;
end

