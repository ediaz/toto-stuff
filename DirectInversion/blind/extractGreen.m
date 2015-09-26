function [ fgreen ] = extractGreen(green,f,dt,isou)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% this function converts a time-domain green function
% into a frequency domain green function
% for source position isou
% Inputs
% green: green's function cube (t,z,x,xs)
% f    : frequency to extract
% dt   : sampling of time axis
% isou : index of source position
%

[nt,nx,nz,ns] = size(green);


fgreen = complex(zeros([nx,nz]));
for it=[1:nt]
    t = it*dt;
    fgreen = fgreen + exp(-2*pi*f*t)*squeeze(green(it,:,:,isou));
end


end

