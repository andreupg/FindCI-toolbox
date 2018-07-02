function [out] = findCI(cfg, comp)
% This function loads the comp structure from FieldTrip ft_componentanalysis
% and returns the components likely to represent the CI artefact. The
% components are selected according to the steps from CIAC (Viola et al 2012).
%
% CONFIGURATION NEEDED:
% - cfg.layout: layout of the topo view
% - cfg.stim_onset: onset times of the stimulus, in seconds (can be a vector)
% - cfg.stim_duration: duration of the stimulus, in seconds (number)
% - cfg.response_latency: vector with the component latency (s) eg.[0.07 0.2]
% - cfg.threshold_ratio: threshold for IC rejection based on the derivative ratio
% - cfg.threshold_corr: threshold for IC rejection based on the correlation
% - cfg.rv: residual variance from the dipole fitting
% - cfg.offset_window: 1 if we want to calculate both onset and offset
%       ratios. 0 (default) for a single window comprising the full stimulus.
% - cfg.artifact_windlgth: window length for the artifact onset and offset
%       windows. The windows are defined from -cfg.artifact_windlgth to
%       +cfg.artifact_windlgth (stim_duration+cfg.artifact_windlgth if 
%       cfg.offset_window is set to 0). Default: 0.01
%
% EXAMPLE
% cfg = [];
% cfg.threshold_ratio = 2.7;
% cfg.threshold_corr = 0.85;
% cfg.rv = rv; % vector containing the residual variance (see
%              % ft_dipolefitting for more info)
% cfg.threshold_rv = 0.1;
% cfg.layout = 'biosemi64.lay';
% cfg.stim_onset = [0:0.5:5];
% cfg.stim_duration = 0.05;
% cfg.response_latency = [0.07 0.2];
% cfg.offset_window = 0;
% cfg.artifact_windlgth = 0.01;
% [CI_comp] = findCI(cfg, comp); % comp is the output from ft_componentanalysis
%
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
assert(exist('comp', 'var') > 0, 'could not find a comp structure');

if isfield(cfg, 'offset_window') % the default is off
else
    cfg.offset_window = 0;
end

if isfield(cfg, 'artifact_windlgth') % the default is 0.01
else
    cfg.artifact_windlgth = 0.01;
end

% Variable alocation and prepare data
% preallocate rejected components
rej_IC = zeros(size(comp.label,1),1);

% cfg for derivative
dvcfg.stim_onset = cfg.stim_onset;
dvcfg.stim_duration = cfg.stim_duration;
dvcfg.response_latency = cfg.response_latency;
dvcfg.offset_window = cfg.offset_window;
dvcfg.artifact_windlgth = cfg.artifact_windlgth;

% IC derivative, ratios, correlation
[dv] = IC_derivative(cfg,comp);

onset_ratio_avg = mean(dv.onset_ratio,2);

if cfg.offset_window == 1 % if the ratio needs to be calculated both in onset and offset
    offset_ratio_avg = mean(dv.offset_ratio,2);
    % Find the IC with the largest ratio
    [largestRatioIC,largestRatioIC_idx] = max([onset_ratio_avg; offset_ratio_avg]);
    if largestRatioIC_idx > length(comp.label)  % detect if max was at offset
        largestRatioIC_idx = largestRatioIC_idx - length(comp.label);
    end
else
   [largestRatioIC,largestRatioIC_idx] = max(onset_ratio_avg); 
end

if cfg.rv(largestRatioIC_idx)<cfg.threshold_rv
    msg = 'The component with the largest ratio has a low residual variance!';
    warning(msg)
end

% correlate ICs with template
template = comp.topo(:,largestRatioIC_idx);
for c = 1:length(comp.label) % for each component
    compcorr(c) = abs(corr(template,comp.topo(:,c)));
end

% Components to be rejected
rej_IC(onset_ratio_avg>cfg.threshold_ratio) = 1; 
if cfg.offset_window == 1 % if the ratio needs to be calculated both in onset and offset
    rej_IC(offset_ratio_avg>cfg.threshold_ratio) = 1;
end
rej_IC(compcorr>cfg.threshold_corr) = 1;
rej_IC(cfg.rv<cfg.threshold_rv) = 0; % do not reject components with low rv

% Output
out.rej_IC = rej_IC; % components to be rejected
out.onset_ratio_avg = onset_ratio_avg; % onset ratio (averaged over all stim presentations)
if cfg.offset_window == 1 % if the ratio needs to be calculated both in onset and offset
    out.offset_ratio_avg = offset_ratio_avg; % offset ratio (averaged over all stim presentations)
end
out.avgcomp = dv.avgcomp; % average of each IC across epochs
out.template = largestRatioIC_idx; % CI template (largest ratio)
out.compcorr = compcorr; % correlation of each component with the template
end
