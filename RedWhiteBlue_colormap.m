function [colormap] = RedWhiteBlue_colormap(steps)
% chose a colormap from red to blue with white in the middle.
% The RGB code can be copied and pasted from colorbrewer2
% This one goes from red to blue passing by white
%
% Copyright (C) 2018  Andreu Paredes Gallardo
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
colormap = [
    103,0,31;
    178,24,43;
    214,96,77;
    244,165,130;
    253,219,199;
    247,247,247;
    209,229,240;
    146,197,222;
    67,147,195;
    33,102,172;
    5,48,97];

% Now increase the resolution if needed 
if steps > 11
    step = 10/steps;
    newcolormap(:,1) = interp1(1:11,colormap(:,1),1:step:11);
    newcolormap(:,2) = interp1(1:11,colormap(:,2),1:step:11);
    newcolormap(:,3) = interp1(1:11,colormap(:,3),1:step:11);
else
    newcolormap = colormap;
end
colormap = [];
colormap = flipud(newcolormap/255);