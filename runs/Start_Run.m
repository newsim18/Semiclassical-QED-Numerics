clear;

warning('off','all');

% PARAMETERS OF THE PROBLEM
mq = (0.001)^(1/2); %mass of field quanta
E0 = 1; % initial electric field strength
w0 = 1; % only for the sauter pulse
a = 0; % initial time
b = 1; % final time
d = 0.01; % k-mode step size
timestep = 0.001; % time step size
q = 1; %charge of field quanta
m = mq;
delta_E0 = 10^(-3); %time perturbation for classical current perturbation expression
tspan = a:timestep:b; % time interval, if you want the Sauter pulse start at t0=-10
Kmin = -100;  % momentum k-value lower cutoff
Kmax = 100;  % momentum k-value upper cutoff
storeY = 1; % toggle storing of .y in file; 0 = Dont Store, 1 = Store
profileType = 0; % 0 = AC, 1 = SP
timeSample = 10; %frequency of saving time: 1=every time loop, 10=every 10th loop, etc.

addpath('../src');
run(fullfile('../src','Solve_LRE.m'));