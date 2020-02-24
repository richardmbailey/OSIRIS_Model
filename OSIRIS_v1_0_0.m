
clc;
clear all;

% -------------------------- Choose run settings --------------------------

% All model inputs (parameter values, forcing, etc.) are located in the 
% input files stored in the 'congif_files' directory. 
% prefix            : Prefix for output data and configuration files
% n_LHS=1           : Number of Latin Hyper-cube samples (for uncertainty 
%                     analysis); set to 1 for central parameter values only
%                     When n_LHS>1, first run is always central values only
%                     with no noise
% plot_network      : set true() or false() to show network plots
% plot_timeseries   : set true() or false() to show timeseries plots
% (when n_LHS>1, only data for last LHS sample are plotted)

prefix = 'test';
n_LHS = 2;
plot_network = true();
plot_timeseries = true(); 

% -------------------------------------------------------------------------


%Create configuration files
config_dir = 'config_files';
infile = config_gen(prefix,n_LHS, config_dir, plot_network, plot_timeseries);

% Run model for each LHS iteration
for LHS_i = 1:n_LHS 
    infile_LHS = [prefix '_' int2str(LHS_i) '.mat'];
    OSIRIS_model_main(infile_LHS, LHS_i, n_LHS);
end


