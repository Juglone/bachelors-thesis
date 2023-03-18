% == Radar parameters B ==
% This file contains all radar parameters for configuration "A".
% To load in a script, simply write the name of this file

%f_start = 77*1e9; % [Hz] start frequency of radar chirp
%idle_time = 100*1e-6; % [s] time between previous and start of next chirp 
%tx_start_time = 0*1e-6; % [s] time from ramp start at which transmitter is turned on
%adc_start_time = 3*1e-6; % [s] time from ramp start when the ADC starts sampling data
num_samples = 96; % number of samples per chirp = ADC samples in mmWave Studio
%slope = 169.993; % slope of chirp

num_frames = 12000;  % number of frames = number of chirps
T_f = 10*1e-3; % frame periodicity in mmWave Studio
T_s = T_f/num_samples; % periodicity of chirp samples

%f_sample = 1/T_f; % sampling frequency = frame frequency = chirp frequency
%f_nyquist = f_sample/2; % nyquist frequency
