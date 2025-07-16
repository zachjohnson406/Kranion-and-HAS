function [ numout ] = was( prompt, numin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
numout = input([prompt,' [',num2str(numin),']  = ']);
if (isempty(numout))
		numout = numin;
end


end