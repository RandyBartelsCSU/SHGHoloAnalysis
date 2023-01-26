function [dCsdO, ang] = MieRCS(rad, n, lambda)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates the RCS of a single (stratified) sphere
%
% Lang Wang, 2023, used the code below:
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)
% -------------------------------------------------------------------------
%                              INPUTS
%
% rad            -> radii of the stritified sphere, from center to edge
% n              -> the RI of the sphere and the background
%
% -------------------------------------------------------------------------
%                              OUTPUTS
%
% dCsdO          -> the normalized RCS
% ang            -> angles of the detectors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ns = n(1:end-1);
nm = n(end);
nang = 1800;        % number of far field angles to evaluate
conv = 1;           % convergence factor

[S, C, ang] = calcmie(rad, ns, nm, lambda, nang, ...
    'ConvergenceFactor', conv);

fctr = 2/pi/C.k;
dCsdOp = fctr*squeeze(abs(S(1,1,:).^2));    % parallel
dCsdOn = fctr*squeeze(abs(S(2,2,:).^2));    % perpendicular
dCsdO = 0.5*(dCsdOp + dCsdOn);              % unpolarized

dCsdO=dCsdO./max(dCsdO);

end

