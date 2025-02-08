
% Backreaction and Linear Response for Spin 1/2 Semiclassical QED
% Ian M. Newsome
% Eric M. Grotzke

format long

% Number of equations required for the time dependent mode functions h1,h2.
% (Kmax-Kmin) is total range of k-modes times 4 columns per mode due to the
% partitioning of Reh1, Imh1, Reh2, Imh2, and this is divided by k-step
% size "d" plus four columns for zero mode k=0, with vector potential and 
% electric field for last two columns
N=4*(Kmax-Kmin)/d+6;

% Initial conditions defined here for h1, h2 modes
Za=zeros(1,N);

for n=0:4:N-3
    K = (n*d)/4+Kmin;
    W = (K^2+mq^2)^(1/2);
    Za(n+1) = sqrt((W-K)/(2*W))*cos(W*a);      % Re(h1)
    Za(n+2) = -sqrt((W-K)/(2*W))*sin(W*a);     % Im(h1)
    Za(n+3) = -sqrt((W+K)/(2*W))*cos(W*a);     % Re(h2)
    Za(n+4) = sqrt((W+K)/(2*W))*sin(W*a);      % Im(h2)
end

% Create folder (if none exists) for saving output data according to profileType
% 0: Asymptotically Constant Background Field Profile
% 1: Sauter Pulse Background Field Profile

if profileType == 1
    dataFolderPath = '../data/SP';
else
    dataFolderPath = '../data/AC';
end

folderName = sprintf('E0=%d_m2=%g', floor(E0), m^2);
mkdir(dataFolderPath,folderName);

% Create output file
filename = sprintf('E0=%g_m2=%g_q=%g_K=%d_trange=%dto%d_tstep=%g_d=%g_e14e14.mat', E0, m^2, q, Kmax, a, b, timestep, d);
fullPath = fullfile(dataFolderPath,folderName,filename);
saveData = matfile(fullPath,'Writable',true);

% Print what values we are running with
fprintf('E0=%g m2=%g q=%g K=%d trange=%dto%d tstep=%g d=%g\n', E0, m^2, q, Kmax, a, b, timestep, d)

%-------------------------------------------------------------------------%

% Solving the semclassical backreaction equation for A, E, h1, h2

% Fix the error tolerance for ODE45 solver
opts = odeset('RelTol',1e-14,'AbsTol',1e-14);

% Set initial mode data for ODE45 solver from initial conditions
y_in=Za;

% storeY declares if full evolving mode data set is stored to file
if storeY
    disp('Saving y to file')
    saveData.y = NaN(floor(length(tspan)/timeSample), N);
    saveData.y(1,:) = y_in;
end

% Full range of modes considered, discretized by k-step size d
k = Kmin:d:Kmax;

% Define quantities of interest
E = zeros(length(tspan),1);
A = zeros(length(tspan),1);
deltaA = zeros(length(tspan),1);
deltaAdot = zeros(length(tspan),1);
deltaJq = zeros(length(tspan),1);
matlength = (length(y_in)-2)/4;

% Define storage for intermediate calculations
Reh1_half = NaN(1,matlength);
Imh1_half = NaN(1,matlength);
Reh2_half = NaN(1,matlength);
Imh2_half = NaN(1,matlength);
Reh1star_half = NaN(1,matlength);
Imh1star_half = NaN(1,matlength);
Reh2star_half = NaN(1,matlength);
Imh2star_half = NaN(1,matlength);
tprimeIntegrand1a_initial = NaN(1,matlength);
tprimeIntegrand1b_initial = NaN(1,matlength);
tprimeIntegrand2a_half = NaN(1,matlength);
tprimeIntegrand2b_half = NaN(1,matlength);
tprimeIntegrand3a_final = NaN(1,matlength);
tprimeIntegrand3b_final = NaN(1,matlength);
Int_tprime1a_initial = NaN(1,matlength);
Int_tprime1b_initial = NaN(1,matlength);
Int_tprime2a_half = NaN(1,matlength);
Int_tprime2b_half = NaN(1,matlength);
Int_tprime3a_final = NaN(1,matlength);
Int_tprime3b_final = NaN(1,matlength);
prev_Int_tprime1a = zeros(1,matlength);
prev_Int_tprime1b = zeros(1,matlength);
tpiece1a_initial = NaN(1,matlength);
tpiece1b_initial = NaN(1,matlength);
tpiece2a_half = NaN(1,matlength);
tpiece2b_half = NaN(1,matlength);
tpiece3a_final = NaN(1,matlength);
tpiece3b_final = NaN(1,matlength);
t_total_1a_initial = NaN(1,matlength);
t_total_1b_initial = NaN(1,matlength);
t_total_2a_half = NaN(1,matlength);
t_total_2b_half = NaN(1,matlength);
t_total_3a_final = NaN(1,matlength);
t_total_3b_final = NaN(1,matlength);
t_total_initial = NaN(1,matlength);
t_total_half = NaN(1,matlength);
t_total_final = NaN(1,matlength);

deltaA(1) = 0;         % IC for vector potential perturbation
deltaAdot(1) = 0;      % IC for vector potential perturbation derivative (- deltaE)

% Defining the classical current perturbation terms
if profileType == 1
    % Sauter Pulse Classical Current
    deltaJ_C = @(t) 2 * q * delta_E0 * w0 * sech(w0 * q * t)^2 * tanh(w0 * q * t);
else
    % Asymptotically Constant Classical Current
    deltaJ_C = @(t) -(q / (1 + q * t)^2) * delta_E0;
end

active_row = 3;         % Time steps to calculate per ODE45 evolution call, minimum required 3
saveRowIndex = 2;       % When saving y data, start at the second row in the file (the first is the initial data)

% Solving backreaction equation
for i=1:active_row-1:length(tspan)-active_row+1

    [t,y] = ode45(@(t,y)SBE(t,y,mq,d,Kmin,N,E0,w0,profileType),tspan(i:i+active_row-1),y_in,opts);

    for c = 0:active_row-2

        fprintf("On run %i\n",i+c)
        
        y_in = y(c+2,:);
    
        if storeY && mod(i+c,timeSample) == 0
            saveData.y(saveRowIndex,:) = y_in;
            saveRowIndex = saveRowIndex + 1;
        end
    
        E(i+c+1,1) = - y_in(N);
        A(i+c+1,1) = y_in(N-1);

    end

    % Here the y-file is partitioned into subarrays for real and imaginary
    % parts of each h1, h2 mode to be used for integration
    length_y = length(y(1,:))-2;
    Reh1 = y(:,1:4:length_y);
    Imh1 = y(:,2:4:length_y);
    Reh2 = y(:,3:4:length_y);
    Imh2 = y(:,4:4:length_y);
    Reh1star = Reh1;
    Imh1star = Imh1;
    Reh2star = Reh2;
    Imh2star = Imh2;

    % Solving linear response equation
    for c = 0:active_row-2
        
        % BEGIN Double integral
    
        % The t'-integral
        % Defining integrands for t'-integrals
        tprimeIntegrand1a_initial(1,:) = (Reh1star(c+1,:) .* Imh2star(c+1,:) + Imh1star(c+1,:) .* Reh2star(c+1,:)) * deltaA(i+c);
        tprimeIntegrand1b_initial(1,:) = (Reh1star(c+1,:) .* Reh2star(c+1,:) - Imh1star(c+1,:) .* Imh2star(c+1,:)) * deltaA(i+c);
            
        Reh1star_half(1,:) = (Reh1star(c+1,:) + Reh1star(c+2,:)) / 2;
        Imh1star_half(1,:) = (Imh1star(c+1,:) + Imh1star(c+2,:)) / 2;
        Reh2star_half(1,:) = (Reh2star(c+1,:) + Reh2star(c+2,:)) / 2;
        Imh2star_half(1,:) = (Imh2star(c+1,:) + Imh2star(c+2,:)) / 2;

        tprimeIntegrand2a_half(1,:) = (Reh1star_half(1,:) .* Imh2star_half(1,:) + Imh1star_half(1,:) .* Reh2star_half(1,:)) * (deltaA(i+c) + (timestep/2) * deltaAdot(i+c));
        tprimeIntegrand2b_half(1,:) = (Reh1star_half(1,:) .* Reh2star_half(1,:) - Imh1star_half(1,:) .* Imh2star_half(1,:)) * (deltaA(i+c) + (timestep/2) * deltaAdot(i+c));
    
        tprimeIntegrand3a_final(1,:) = (Reh1star(c+2,:) .* Imh2star(c+2,:) + Imh1star(c+2,:) .* Reh2star(c+2,:)) * (deltaA(i+c) + timestep * deltaAdot(i+c));
        tprimeIntegrand3b_final(1,:) = (Reh1star(c+2,:) .* Reh2star(c+2,:) - Imh1star(c+2,:) .* Imh2star(c+2,:)) * (deltaA(i+c) + timestep * deltaAdot(i+c));
    
        % Implementing midpoint method for solving t'-integral
        Int_tprime1a_initial(1,:) = prev_Int_tprime1a(1,:);
        Int_tprime1b_initial(1,:) = prev_Int_tprime1b(1,:);

        Int_tprime2a_half(1,:) = (timestep/2) * (tprimeIntegrand1a_initial(1,:) + tprimeIntegrand2a_half(1,:)) / 2 + Int_tprime1a_initial(1,:);
        Int_tprime2b_half(1,:) = (timestep/2) * (tprimeIntegrand1b_initial(1,:) + tprimeIntegrand2b_half(1,:)) / 2 + Int_tprime1b_initial(1,:);

        Int_tprime3a_final(1,:) = (timestep/2) * (tprimeIntegrand2a_half(1,:) + tprimeIntegrand3a_final(1,:)) / 2 + Int_tprime2a_half(1,:);
        Int_tprime3b_final(1,:) = (timestep/2) * (tprimeIntegrand2a_half(1,:) + tprimeIntegrand3b_final(1,:)) / 2 + Int_tprime2b_half(1,:);

        % Storing i'th integral contribution to be used for later integration steps
        % due to the time dependent upper limit of integration
        prev_Int_tprime1a(1,:) = Int_tprime3a_final(1,:);
        prev_Int_tprime1b(1,:) = Int_tprime3b_final(1,:);

        % Creating arrays corresponding to functions of time "t"
        % multiplying t' integral
        tpiece1a_initial(1,:) = Imh1(c+1,:) .* Imh2(c+1,:) - Reh1(c+1,:) .* Reh2(c+1,:);
        tpiece1b_initial(1,:) = Reh1(c+1,:) .* Imh2(c+1,:) + Imh1(c+1,:) .* Reh2(c+1,:);

        Reh1_half(1,:) = (Reh1(c+1,:) + Reh1(c+2,:))./2;
        Imh1_half(1,:) = (Imh1(c+1,:) + Imh1(c+2,:))./2;
        Reh2_half(1,:) = (Reh2(c+1,:) + Reh2(c+2,:))./2;
        Imh2_half(1,:) = (Imh2(c+1,:) + Imh2(c+2,:))./2;
    
        tpiece2a_half(1,:) = Imh1_half(1,:) .* Imh2_half(1,:) - Reh1_half(1,:) .* Reh2_half(1,:);
        tpiece2b_half(1,:) = Reh1_half(1,:) .* Imh2_half(1,:) + Imh1_half(1,:) .* Reh2_half(1,:);
    
        tpiece3a_final(1,:) = Imh1(c+2,:) .* Imh2(c+2,:) - Reh1(c+2,:) .* Reh2(c+2,:);
        tpiece3b_final(1,:) = Reh1(c+2,:) .* Imh2(c+2,:) + Imh1(c+2,:) .* Reh2(c+2,:);

        % Performing multiplication
        t_total_1a_initial(1,:) = tpiece1a_initial(1,:) .* Int_tprime1a_initial(1,:);
        t_total_1b_initial(1,:) = tpiece1b_initial(1,:) .* Int_tprime1b_initial(1,:);
    
        t_total_2a_half(1,:) = tpiece2a_half(1,:) .* Int_tprime2a_half(1,:);
        t_total_2b_half(1,:) = tpiece2b_half(1,:) .* Int_tprime2b_half(1,:);
    
        t_total_3a_final(1,:) = tpiece3a_final(1,:) .* Int_tprime3a_final(1,:);
        t_total_3b_final(1,:) = tpiece3b_final(1,:) .* Int_tprime3b_final(1,:);

        % Adding contributions together
        t_total_initial(1,:) = t_total_1a_initial(1,:) + t_total_1b_initial(1,:);
        t_total_half(1,:) = t_total_2a_half(1,:) + t_total_2b_half(1,:);
        t_total_final(1,:) = t_total_3a_final(1,:) + t_total_3b_final(1,:);
    
        % The k-integral
        doubleInt_initial = simp13(t_total_initial(1,:),d);
        doubleInt_half = simp13(t_total_half(1,:),d);
        doubleInt_final = simp13(t_total_final(1,:),d);

        % END Double integral

        % BEGIN RK4 method to solve system of first order DE's
    
        % d/dt deltaA = chi
        F = @(t,deltaA,deltaAdot) deltaAdot;
    
        % d^2/dt^2 deltaA = d/dt chi
        G_initial = @(t,deltaA,deltaAdot) deltaJ_C(t) - ((q^2)/(pi)) * deltaA - ((4*q^2)/pi) * doubleInt_initial;
        G_half = @(t,deltaA,deltaAdot) deltaJ_C(t) - ((q^2)/(pi)) * deltaA - ((4*q^2)/pi) * doubleInt_half;
        G_final = @(t,deltaA,deltaAdot) deltaJ_C(t) - ((q^2)/(pi)) * deltaA - ((4*q^2)/pi) * doubleInt_final;
    
        k1 = timestep * F(tspan(i+c),deltaA(i+c),deltaAdot(i+c));
        L1 = timestep * G_initial(tspan(i+c),deltaA(i+c),deltaAdot(i+c));
    
        k2 = timestep * F(tspan(i+c)+timestep/2,deltaA(i+c)+k1/2,deltaAdot(i+c)+L1/2);
        L2 = timestep * G_half(tspan(i+c)+timestep/2,deltaA(i+c)+k1/2,deltaAdot(i+c)+L1/2);
    
        k3 = timestep * F(tspan(i+c)+timestep/2,deltaA(i+c)+k2/2,deltaAdot(i+c)+L2/2);
        L3 = timestep * G_half(tspan(i+c)+timestep/2,deltaA(i+c)+k2/2,deltaAdot(i+c)+L2/2);
    
        k4 = timestep * F(tspan(i+c)+timestep,deltaA(i+c)+k3,deltaAdot(i+c)+L3);
        L4 = timestep * G_final(tspan(i+c)+timestep,deltaA(i+c)+k3,deltaAdot(i+c)+L3);
    
        deltaA(i+c+1)    = deltaA(i+c) + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        deltaAdot(i+c+1) = deltaAdot(i+c) + (L1 + 2 * L2 + 2 * L3 + L4) / 6;

        % END RK4

        deltaJq(i+c) = - ((q^2)/(pi)) * deltaA(i+c) - ((4*q^2)/pi) * doubleInt_initial;

    end

end

saveData.A = A;
saveData.E = E;
saveData.deltaA = deltaA;
saveData.deltaE = - deltaAdot;
saveData.deltaJq = deltaJq;
