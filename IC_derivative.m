function [dvtv] = IC_derivative(cfg,comp)
% This function computes the average across trials for each component and
% calculates its time derivative. The ratio between the onset and offset
% time windows and the response window (derivative) is calculated - See
% Viola et al 2012 (EEGLAB CIAC plugin).
%
% The function takes the comp structure from fieldtrip and the cfg structure
% as input:
% - cfg.stim_onset: can be a vector of onset times.
% - cfg.stim_duration: druation of the stimulus (number)
% - cfg.response_latency: time window where we expect the neural response
%           (vector, first number is onset re. stimulus onset and second
%           number is offset re. stimulus onset).
% - cfg.offset_window: 1 if we want to calculate both onset and offset
%           ratios. 0 (default) for a single window comprising the full stimulus.
% - cfg.artifact_windlgth: window length for the artifact onset and offset
%           windows. The windows are defined from -cfg.artifact_windlgth to
%           +cfg.artifact_windlgth (stim_duration+cfg.artifact_windlgth if 
%           cfg.offset_window is set to 0).
%
% The function returns the onset and offset ratios as well as the averaged
% component matrix (across trials).
%
% REQUIREMENTS: 
% - "derivative.m" from Scott McKinney
% - "findclosest.m"
%
% Copyright (C) 2018  Andreu Paredes Gallardo
%
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

% Check that cfg is provided
assert(exist('comp', 'var') > 0, 'could not find a comp structure');

if isfield(cfg, 'offset_window') % the default is off
else
    cfg.offset_window = 0;
end

if isfield(cfg, 'artifact_windlgth') % the default is 0.01
else
    cfg.artifact_windlgth = 0.01;
end

% Compute the average for each component across trials
temp = 0;
for trl = 1:length(comp.trial)
    temp =  temp + comp.trial{trl}(:,:);
end
avgcomp = temp/length(comp.trial);

for i = 1:length(comp.label) % loop for each component
    
    % Compute the derivative of the component average
    deriv=derivative(avgcomp(i,:));
    
    % Calculate rms over stimulus and response time window
    for stim = 1:length(cfg.stim_onset)
        response_start = findclosest(comp.time{1},cfg.stim_onset(stim)+cfg.response_latency(1));
        response_end = findclosest(comp.time{1},cfg.stim_onset(stim)+cfg.response_latency(2));
        stim_start = findclosest(comp.time{1},cfg.stim_onset(stim)-cfg.artifact_windlgth);
        stim_end = findclosest(comp.time{1},cfg.stim_onset(stim)+cfg.stim_duration+cfg.artifact_windlgth);
        
        rms_deriv(i,stim)=sqrt(mean(deriv(response_start:response_end).^2)); %%% rms for N1_P2 time window 
        
        if cfg.offset_window == 1 % if the ratio needs to be calculated both in onset and offset
            onset_end = findclosest(comp.time{1},cfg.stim_onset(stim)+cfg.artifact_windlgth);
            offset_start = findclosest(comp.time{1},cfg.stim_onset(stim)+cfg.stim_duration-cfg.artifact_windlgth);
            
            rms_on_deriv(i,stim)=sqrt(mean(deriv(stim_start:onset_end).^2)); %%% rms in stim onset window
            rms_off_deriv(i,stim)=sqrt(mean(deriv(offset_start:stim_end).^2)); %%% rms in stim offset window
            
            onset_ratio(i,stim)=rms_on_deriv(i,stim)/rms_deriv(i,stim); % onset ratio
            offset_ratio(i,stim)=rms_off_deriv(i,stim)/rms_deriv(i,stim); % offset ratio
            
        else % if not, a single window includes the full stimulus (from -cfg.artifact_windlgth to stim_duration+cfg.artifact_windlgth)
            rms_on_deriv(i,stim)=sqrt(mean(deriv(stim_start:stim_end).^2)); %%% rms in stim onset window
            
            onset_ratio(i,stim)=rms_on_deriv(i,stim)/rms_deriv(i,stim); % onset ratio
            offset_ratio = [];
        end

    end
end

% output
dvtv.avgcomp = avgcomp;
dvtv.onset_ratio = onset_ratio;
dvtv.offset_ratio = offset_ratio;

% eof