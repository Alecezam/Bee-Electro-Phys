
function olfact()
% The main script for the olfactometr
% olfact - starts the package with the new olfactometer harware
% olfact(HardwarePresent) - if HardwarePresent is 0 olfactometer starts in 
%                           simulation mode 
% global variables:
% hw   - hardware and graphics
% fn   - file names
% sess - session parameters
% delete all timers and close all windows

delete(timerfind);
delete(instrfind);
close all;
clear all;
global  fn 

% Set up file names information
cd('C:\Users\basuchap\Desktop\Ephys_05042018');

prog_fold      = pwd;

% folder for protocols
fn.protocol_fold  = fullfile(prog_fold,   'protocol_curr');
fn.config_fold = fullfile('C:/', 'Experiment', 'Behavioral_rigs_config');
addpath(fn.protocol_fold)

% loading config information
config

init_olfact_v2
%init_exper