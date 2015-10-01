function [ fgreen ] = One_suo_Green(green,f,dt)
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

[nx,nz,nt] = size(green);


fgreen = complex(zeros([nx,nz]));
for it=[1:nt]
    t = it*dt;
    fgreen = fgreen + exp(-2j*pi*f*t)*squeeze(green(:,:,it));
end


end

