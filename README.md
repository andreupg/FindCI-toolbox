# FindCI toolbox

This repository contains a set of matlab functions to help in the process of identifying independent components (ICs) representing the cochlear implant (CI) artifact when using independent component analysis (ICA) and the [Fieldtrip toolbox](http://www.fieldtriptoolbox.org/). The functions are an implementation, using the fieldtrip structure, of the procedure proposed and validated by [Viola et. al. (2012)](http://dx.doi.org/10.1016/j.heares.2011.12.010). An EEGLAB implementation can be found [here](http://www.debener.de/CIAC_tutorial/ciacplugin.html). It is highly recommended to read the description provided by [Viola et. al. (2012)](http://dx.doi.org/10.1016/j.heares.2011.12.010) before using this toolbox.

## Main functions ##

### findCI.m ###
Main function of the toolbox. The function requires a `comp` structure as obtained from [`ft_componentanalysis`](http://www.fieldtriptoolbox.org/reference/ft_componentanalysis) as well as the residual variance (rv) vector obtained from [`ft_dipolefitting`](http://www.fieldtriptoolbox.org/reference/ft_dipolefitting).
#### Inputs ####
* `cfg`
	* `cfg.layout`: layout of the topo view.
	* `cfg.stim_onset`: onset of the stimulus, in seconds (can be a vector if there are multiple stimulus within a trial).
	* `cfg.stim_duration`: duration of the stimulus, in seconds (number).
	* `cfg.response_latency`: vector with the component latency in seconds (i.e. [start end]).
	* `cfg.threshold_ratio`: threshold for IC rejection based on the derivative ratio.
	* `cfg.threshold_corr`: threshold for IC rejection based on the correlation.
	* `cfg.rv`: residual variance from the dipole fitting.
	* `cfg.offset_window`: can be set to 0 (default) or to 1. Set to 0 for a single window comprising the full stimulus (recommended for short stimulus, not overlapping with the response window). 1 if we want to calculate both onset and offset ratios (recommended when the response window overlaps the stimulus, or when the stimulus is relatively long). 
	* `cfg.artifact_windlgth`: window length, in seconds, for the artifact onset and offset windows. The windows are defined from -cfg.artifact_windlgth to +cfg.artifact_windlgth (cfg.stim_duration + cfg.artifact_windlgth if cfg.offset_window is set to 0). Default: 0.01.
* `comp`

#### Outputs ####
* `out.rej_IC`: vector with the ICs tagged for rejection
* `out.onset_ratio_avg`: vector with the ratio between the rms of the component's derivative onset window and the rms in the component's latency window.
* `out.offset_ratio_avg`: vector with the ratio between the rms of the component's derivative offset window and the rms in the component's latency window (if `cfg.offset_window` is set to 0, the offset ratio is not calculated and only the onset ratio, calculated over the total duration of the stimulus, is returned).
* `out.avgcomp`: average of each IC across all epochs.
* `out.template`: IC with the largest rms ratio (CI template).
* `out.compcorr`: correlation coefficients for each IC and the template.

### CI_comp_GUI.m ###
Graphical interface for visualization of the ICs. The interface also displays the time average, first derivative, power spectrum, within-trial variance, the topography of the components and it offers the possibility to plot the timecourse for each IC and trial. The components likely to represent the CI artefact are selected according to the steps proposed by [Viola et. al. (2012)](http://dx.doi.org/10.1016/j.heares.2011.12.010). The user can toggle components to be kept or to be rejected. Once the user pressess `Exit`, a list of the components chosen for rejection is returned.
#### Inputs ####
* `cfg`
	* `cfg.CI_comp`: output from `findCI.m`.
	* `cfg.layout`: layout of the topo view.
	* `cfg.stim_onset`: onset of the stimulus, in seconds (can be a vector if there are multiple stimulus within a trial).
	* `cfg.stim_duration`: duration of the stimulus, in seconds (number).
	* `cfg.threshold_ratio`: threshold for IC rejection based on the derivative ratio.
	* `cfg.warning_ratio`: threshold for visual warning based on the derivative ratio.
	* `cfg.threshold_corr`: threshold for IC rejection based on the correlation.
	* `cfg.warning_corr`: threshold for visual warning based on the correlation.
	* `cfg.rv`: residual variance from the dipole fitting.
	* `cfg.offset_window`: can be set to 0 (default) or to 1. Set to 0 for a single window comprising the full stimulus (recommended for short stimulus, not overlapping with the response window). 1 if we want to calculate both onset and offset ratios (recommended when the response window overlaps the stimulus, or when the stimulus is relatively long).
* `comp`

#### Output ####
The interface returns the vector `rej_IC`, containing the ICs tagged for rejection.

## Requirements ##

* Matlab (tested with Matlab 2016b),
* [Fieldtrip toolbox](http://www.fieldtriptoolbox.org/),

## Issues or questions? ##
If you experience any issues or have any questions, please use the [issue tracking page](https://bitbucket.org/andreupg/ci_artifact-toolbox/issues?status=new&status=open).

## License ##
This toolbox was written by Andreu Paredes Gallardo and belongs to the Technical University of Denmark. The toolbox is distributed under the [GPL v3 licence](https://www.gnu.org/licenses/). A copy of the license is distributed with the toolbox.

The function `derivative.m` from Scott McKinney (c) 2010 is distributed with this toolbox. The function is also available [here](https://www.mathworks.com/matlabcentral/fileexchange/28920-derivative?focused=5169606&tab=function). This function is distributed with its own license, which can be found at the beginning of the file.


