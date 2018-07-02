function [rej_IC] = CI_comp_GUI(cfg, comp)
% ICA component viewer CI artefact rejection tool
%
% This tool presents a GUI with the time average, first derivative,
% power spectrum, within-trial variance, the topography of the components
% as well as the possibility view the timecourse.
% The components likely to represent the CI artefact are selected according
% to the steps from CIAC (Viola et al 2012). The user can
% toggle components to be kept or to be rejected. Once the user pressess
% "Exit", a list of the components chosen for rejection is returned.
%
% The GUI displays the value for the ratio between the rms of the
% component's derivative onset/offset window and the rms in the component
% latency window. The value is marked in red if this exceeds cfg.threshold_ratio
% or in yellow if it is larger than cfg.warning_ratio. The component is not
% toggled as rejected if the ratio is smaller than the rejection threshold
% (cfg.threshold_ratio).
% The component with the largest ratio is chosen as template. The
% topography of the remaining components is correlated with the template,
% and the correlation value is also displayed in the GUI. The residual 
% variance (rv) between the ICs and the dipole model is used to keep
% components likely to represent neural activity and is displayed in the 
% GUI. Both the correlation and the rv values are color coded as the ratio
% (red-yellow-green).
%
% CONFIGURATION NEEDED:
% cfg.CI_comp: output from findCI.m
% cfg.layout: layout of the topo view
% cfg.stimulus_onset: onset times of the stimulus, in seconds (can be a vector)
% cfg.stimulus_duration: duration of the stimulus, in seconds (number)
% cfg.threshold_ratio: threshold for IC rejection based on the derivative ratio
% cfg.warning_ratio: threshold for visual warning based on the derivative ratio
% cfg.threshold_corr: threshold for IC rejection based on the correlation
% cfg.warning_corr: threshold for visual warning based on the correlation
% cfg.rv: residual variance from the dipole fitting
% cfg.threshold_rv: components with lower rv are not rejected
% cfg.offset_window: 1 if we want to calculate both onset and offset
%           ratios. 0 (default) for a single window comprising the full stimulus.
%
% EXAMPLE
% cfg = [];
% cfg.CI_comp = CI_comp; % CI_comp is the output from find_CI_comp.m
% cfg.threshold_ratio = 2.7;
% cfg.warning_ratio = 1.5;
% cfg.threshold_corr = 0.85;
% cfg.warning_corr = 0.8;
% cfg.rv = rv;
% cfg.threshold_rv = 0.1;
% cfg.layout = 'biosemi64.lay';
% cfg.stim_onset = [0:0.5:5];
% cfg.stim_duration = 0.05;
% [rej_IC] = CI_comp_GUI(cfg, comp); % comp is the output from ft_componentanalysis
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


%% Variable alocation and prepare data
% do the fft on a subset of trials to save time
fft_data = cat(2,comp.trial{1:5:end});

% preallocate variance
comp_var = zeros(size(comp.label,1),length(comp.trial));

% offset ratio?
if isfield(cfg.CI_comp, 'offset_ratio_avg') % the default is off
    show_offset_ratio = 1;
else
    show_offset_ratio = 0;
end

% pass cfg.layout to the variable "lay"
lay = cfg.layout;

cfgtopo = [];
cfgtopo.layout    = lay;     % specify the layout file that should be used for plotting
cfgtopo.comment   = 'no';
cfgtopo.highlight = 'off';
cfgtopo.marker    = 'on';
cfgtopo.style     = 'both';
cfgtopo.zlim = 'maxabs';
if isfield(cfg, 'colormap')
    cfgtopo.colormap  = cfg.colormap;
else
    cfgtopo.colormap = RedWhiteBlue_colormap(100);
end

%% Compute variance across epochs
% Compute variance for epochs
for ep = 1:length(comp.trial) % for each component
    comp_var(:,ep)=var(comp.trial{ep},0,2);
end

%% Rejected components
rej_IC = cfg.CI_comp.rej_IC;

%% Start Loop (GUI)
proceed = 1;
i = 1;
manpos = [0.05 0.05 0.9 0.9]; % figure position
while proceed
f = figure('units','normalized','outerposition', manpos);
    
    % COMPUTE POWER SPECTRUM
    Fs = comp.fsample;
    N = floor(size(fft_data,2));
    FFT = fft(fft_data(i,:));
    FFT = FFT(1:N/2+1);
    POW = (1/(Fs*N)).*abs(FFT).^2;
    POW(2:end-1) = 2*POW(2:end-1);
    
    j = 1;
    k = 1;
    smo = 50;
    stepnum = 10;
    while j<length(POW)-smo
        POW_smo(k)=mean(POW(j:j+smo));
        j = j+stepnum;
        k = k+1;
    end
    
    freq = linspace(0,Fs/2,size(POW_smo,2));
    start_freq = find(freq > 2,1,'first');
    stop_freq  = find(freq < 20,1,'last');
    
    % COMPUTE COMPONENT DERIVATIVE
    deriv=derivative(cfg.CI_comp.avgcomp(i,:));
    
    % PLOT AVERAGE TIME SIGNAL
    t = [0 0.8:0.22:4.76]; % onset times
    subcomp{1}{i} = subplot(2,4,[1 2]);
    for stim = 1:length(t)
        plot([t(stim) t(stim)], [max(cfg.CI_comp.avgcomp(i,:)) min(cfg.CI_comp.avgcomp(i,:))],...
            'color', [0.5 0.5 0.5])
        hold on
        plot([t(stim)+0.05 t(stim)+0.05], [max(cfg.CI_comp.avgcomp(i,:)) min(cfg.CI_comp.avgcomp(i,:))],...
            ':', 'color', [0.5 0.5 0.5])
    end
    plot(comp.time{1},cfg.CI_comp.avgcomp(i,:),'k')
    ylabel('ampl (au)');
    xlabel('Time (s)'); grid on;
    title('Average waveform');
    axis tight;
    
    p = get(subcomp{1}{i}, 'pos');
    p(1) = p(1) - 0.09;
    p(3) = p(3) + 0.1;
    p(4) = p(4) + 0.03;
    set(subcomp{1}{i}, 'pos', p);
    
    % PLOT COMPONENT DERIVATIVE
    subcomp{2}{i} = subplot(2,4,[5 6]);
    for stim = 1:length(t)
        plot([t(stim) t(stim)], [max(deriv(3:end-1)) min(deriv(3:end-1))],...
            'color', [0.5 0.5 0.5])
        hold on
        plot([t(stim)+0.05 t(stim)+0.05], [max(deriv(3:end-1)) min(deriv(3:end-1))],...
            ':', 'color', [0.5 0.5 0.5])
    end
    plot(comp.time{1}(3:end-1),deriv(3:end-1),'k')
    ylabel('ampl (au)');
    xlabel('Time (s)'); grid on;
    title('Average waveform: 1st derivative');
    axis tight;
    
    p = get(subcomp{2}{i}, 'pos');
    p(1) = p(1) - 0.09;
    p(3) = p(3) + 0.1;
    p(4) = p(4) + 0.03;
    set(subcomp{2}{i}, 'pos', p);
    
    % PLOT POWER SPECTRUM
    subcomp{3}{i} = subplot(2,4,3);
    plot(freq(start_freq:stop_freq),log10(POW_smo(start_freq:stop_freq)),'k');
    ylabel('(dB/Hz)');
    set(gca,'TickDir','out','XTick',0:2:20)
    xlabel('Frequency (Hz)'); grid on;
    title('Power spectrum');
    axis tight;
    
    p = get(subcomp{3}{i}, 'pos');
    p(4) = p(4) + 0.03;
    set(subcomp{3}{i}, 'pos', p);
    
    % PLOT VARIANCE OVER TIME
    subcomp{4}{i} = subplot(2,4,[7 8]);
    scatter([1:length(comp.trial)],comp_var(i,:),'k.');
    xlabel('Trial number'); ylabel('Variance'); title('Within-trial variance');
    axis tight; set(gca, 'tickdir', 'out');
    
    p = get(subcomp{4}{i}, 'pos');
    p(4) = p(4) + 0.03;
    set(subcomp{4}{i}, 'pos', p);
    
    % PLOT COMPONENT TOPOGRAPHY
    subcomp{5}{i} = subplot(2,4,4);
    cfgtopo.component = i;       % specify the component(s) that should be plotted
    ft_topoplotIC(cfgtopo, comp);
    
    p = get(subcomp{5}{i}, 'pos');
    p(4) = p(4) + 0.03;
    set(subcomp{5}{i}, 'pos', p);
    
    % SHOW TIMECOURSE OF THIS COMPONENT
    tc = uicontrol('Units','normalized','Position',[0.92 0.8 0.075 0.045],...
        'Style','pushbutton','String','Timecourse','Callback',{@tc_callback, 1});
    
    % REJECT COMPONENT
    if rej_IC(i) == 1
        bkgr_col = 'r';
    else
        bkgr_col = 'g';
    end
    
    pos = [0.92 0.65 0.075 0.045];
    rej = uicontrol('Units','normalized', 'Tag', 'rej1', 'Position',...
        pos(1,:),'Style','pushbutton','String','Keep', ...
        'Backgroundcolor',bkgr_col,'Callback',@rej_callback);
    
    % MOVE TO NEXT AND PREV
    btn_space = 0.03; % space between groups of buttons
    btn_width = 0.05; % button width
    btn_hei = 0.04; % button heigth
    left_x = 0.05; % x location of the first button
    
    % Previous rejected component
    pos = [left_x 0.01 btn_width btn_hei];
    if i > 1
       rejprev = uicontrol('Units','normalized','Position',pos,'Style',...
           'pushbutton','String','Prev rej','Callback',@lastrejplot_callback);
    else
        rejprev = uicontrol('Units','normalized','Position',pos,'Style',...
            'pushbutton','String','');
    end
    
    % Next rejected component
    pos = [left_x+btn_width 0.01 btn_width btn_hei];
    if i < size(comp.label,1)-0
        rejnext = uicontrol('Units','normalized','Position',pos,'Style',...
            'pushbutton','String','Next rej','Callback',@nextrejplot_callback);
    else
        rejnext = uicontrol('Units','normalized','Position',pos,'Style',...
            'pushbutton','String','');
    end
    
    % Previous component
    pos = [left_x+2*btn_width+btn_space 0.01 btn_width btn_hei];
    if i > 1
        prev = uicontrol('Units','normalized','Position',pos,'Style',...
            'pushbutton','String','Prev','Callback',@lastplot_callback);
    else
        prev = uicontrol('Units','normalized','Position',pos,'Style',...
            'pushbutton','String','');
    end
    
    % Next component
    pos = [left_x+3*btn_width+btn_space 0.01 btn_width btn_hei];
    if i < size(comp.label,1)-0
        next = uicontrol('Units','normalized','Position',pos,'Style',...
            'pushbutton','String','Next','Callback',@nextplot_callback);
    else
        next = uicontrol('Units','normalized','Position',pos,'Style',...
            'pushbutton','String','');
    end
    
    % Trial selection box
    trialx = uicontrol('Units','normalized','Position',[0.31 0.01 0.05 0.04],...
        'Style','edit','String',num2str(i));
    gotrialx = uicontrol('Units','normalized','Position',[0.36 0.01 0.05 0.04],...
        'Style','pushbutton','String','go','Callback',@plotx_callback);

    % DISPLAY INFO
    on_ratio=uicontrol('Units','normalized');
    set(on_ratio,'position',[0.92 0.4 0.075 0.045]);
    set(on_ratio,'style','text');
    set(on_ratio,'string',{'Onset ratio',cfg.CI_comp.onset_ratio_avg(i)});
    if cfg.CI_comp.onset_ratio_avg(i) > cfg.threshold_ratio
        set(on_ratio,'Backgroundcolor',[255 99 71]/255)
    elseif cfg.CI_comp.onset_ratio_avg(i) > cfg.warning_ratio
        set(on_ratio,'Backgroundcolor',[255 165 0]/255)
    end
    
    if show_offset_ratio ==1
        off_ratio=uicontrol('Units','normalized');
        set(off_ratio,'position',[0.92 0.32 0.075 0.045]);
        set(off_ratio,'style','text');
        set(off_ratio,'string',{'Offset ratio',cfg.CI_comp.offset_ratio_avg(i)});
        if cfg.CI_comp.offset_ratio_avg(i) > cfg.threshold_ratio
            set(off_ratio,'Backgroundcolor',[255 99 71]/255)
        elseif cfg.CI_comp.offset_ratio_avg(i) > cfg.warning_ratio
            set(off_ratio,'Backgroundcolor',[255 165 0]/255)
        end
    end
    
    
    corrinfo=uicontrol('Units','normalized');
    set(corrinfo,'position',[0.92 0.24 0.075 0.045]);
    set(corrinfo,'style','text');
    set(corrinfo,'string',{'Topo corr',cfg.CI_comp.compcorr(i)});
    if cfg.CI_comp.compcorr(i) > cfg.threshold_corr
        set(corrinfo,'Backgroundcolor',[255 99 71]/255)
    elseif cfg.CI_comp.compcorr(i) > cfg.warning_corr
        set(corrinfo,'Backgroundcolor',[255 165 0]/255)
    end
    
    templinfo=uicontrol('Units','normalized');
    set(templinfo,'position',[0.92 0.16 0.075 0.045]);
    set(templinfo,'style','text');
    set(templinfo,'string',{'Template IC',cfg.CI_comp.template});
    
    if isfield(cfg, 'rv')
        rvinfo=uicontrol('Units','normalized');
        set(rvinfo,'position',[0.92 0.10 0.075 0.045]);
        set(rvinfo,'style','text');
        set(rvinfo,'string',{'RV',[num2str(cfg.rv(i)*100) '%']});
        if cfg.rv(i) < cfg.threshold_rv
            set(rvinfo,'Backgroundcolor',[173 255 47]/255)
        else
            set(rvinfo,'Backgroundcolor',[255 165 0]/255)
        end
    end
    
    if rej_IC(i) == 1
        set(rej,'Backgroundcolor','r','String','Reject')
    end
    
    % SAVE AND EXIT
    ext = uicontrol('Units','normalized','Position',[0.80 0.01 0.075 0.05],...
        'Style','pushbutton','String','Exit','Callback',@exit_callback);
    orient landscape
    uiwait
end

%% callback functions
    function rej_callback(h, evt)
        rej = findobj('Tag','rej1');
        if (rej_IC(i) == 0),
            set(rej,'Backgroundcolor','r','String','Reject'),
            rej_IC(i)=1;
        else
            set(rej,'Backgroundcolor','g','String','Keep'),
            rej_IC(i)=0;
        end
    end

    function tc_callback(h, evt, whichcomp)
        cfgtc = [];
        cfgtc.layout = lay;
        cfgtc.viewmode = 'butterfly';
        cfgtc.channel = [i];
        % cfgtc.ylim = [-5e-13 5e-13];
        ft_databrowser(cfgtc, comp);
    end

    function nextplot_callback(h, evt)
        manpos = get(f,'OuterPosition');
        i = i + 1;
        close;
    end

    function lastplot_callback(h, evt)
        manpos = get(f,'OuterPosition');
        i = i - 1; 
        close;
    end

    function nextrejplot_callback(h, evt)
        rej_idx = find(rej_IC);
        manpos = get(f,'OuterPosition');
        i = rej_idx(find(rej_idx>i,1,'first'));
        close;
    end

    function lastrejplot_callback(h, evt)
        rej_idx = find(rej_IC);
        manpos = get(f,'OuterPosition');
        i = rej_idx(find(rej_idx<i,1,'last'));
        close;
    end

    function plotx_callback(h, evt)
        manpos = get(f,'OuterPosition');
        i = str2num(trialx.String); close;
    end

    function exit_callback(h, evt)
        close(f);
        proceed = 0;
    end

end