function [micarr, pc] = gen_micarray(firstcoord, d, M)
%This function generates the cartesian coordinates for a linear microphone
%array. It will be subsequently used for the rir_generator.m and
%gen_micsignals.m files to get RIRs. The mic. array will be generated such
%that the x and z coordinates are fixed and we only vary the y-coordinate
%as the microphone spacing.
%
%       micarr = gen_micarray(firstcoord, d, M)
%
%INPUTS:
% firstcoord - first coordinate of the microphone array. This will define
%              the mic. with the rightmost horizontal (x) coordinate, i.e.
%              the right of the mic. array will be this coordinate.
%              Specify as [x ;y ;z] coordinate. Dimension = 3x1. 
%
% d - microphone spacing specified in metres (m)
%
% M - the number of microphones in the array
%
%OUTPUTS:
% micarr - this is the matrix of x,y,z coordinates of the dimensions of the
%          specified microphone array. Dimension = 3 x M. The reference
%          microphone will be the first matrix entry.
%
% pc - centre of the microphone array. Defined in eqn 5.12 (Ivan Tashev
% Sound capture and processing)


micarr = zeros(3,M);
micarr(:,1) = firstcoord;
x1 = firstcoord(1);
y1 = firstcoord(2);
z1 = firstcoord(3);

for m = 2:M
    micarr(:,m) = [x1-((m-1)*d) ;y1; z1]; %adjust the x coordinate accordingly
end

pc = (1/M) * sum(micarr,2);


