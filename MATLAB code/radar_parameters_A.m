% == Radar parameters A ==
% This file contains all radar parameters for configuration "A".
% To load in a script, simply write the name of this file

num_samples = 96; % number of samples per chirp = ADC samples in mmWave Studio
num_frames = 500; % number of frames = number of chirps
T_f = 40*1e-3; % frame periodicity in mmWave Studio
T_s = T_f/num_samples; % periodicity of chirp samples

t_meas = num_frames*T_f; % length of radar meas in seconds