function Y_sol = SSAdapt_run_delay(simtime,dt,parms)

parms.samenoise = false;

%% Establish Vectors etc
N_neurons = parms.N_neurons;
W = parms.W;
I_in = parms.I_in;
beta = parms.beta;
tau_r = parms.tau_r;
tau_a = parms.tau_a;
Ak = parms.Ak;
A0 = parms.A0;
k = parms.k;
h = parms.h;
n = parms.n;
LR = parms.LR;
noiseamp = parms.noiseamp;
noisefreq = parms.noisefreq;




%initial conditions: no activity. 
%                    -can add different initial conditions if desired.
%r_init = zeros(1,N_neurons);
r_init = rand(1,N_neurons);
%a_init = zeros(1,N_neurons);
a_init = rand(1,N_neurons);

y0 = [r_init,a_init];       %combine initial conditions 
                             %into one long vector for ode45
                                            
%time and solution arrays for ode45
t_tot = simtime/dt;             %number of iterations
tspan = [0:dt:simtime]';         %interval of integration


% %Do these need to be here?....
% r_sol = zeros(t_tot+1,N_neurons);       %solution array for E cells
% a_sol = zeros(t_tot+1,N_neurons);     %solution array for theta  
% T = tspan;                     %column vector of time points
% Y_sol = [r_sol,a_sol];          %solution array for ode45

%% Noise for Input
%noisefilter = [0.5 10];
%Inoise = noiseamp.*wgn(length(tspan),1,0);
switch parms.samenoise
    case true
        numsignals = 1;
    case false
        numsignals = N_neurons;
end 
[ Inoise,noiseT ] = OUNoise(noisefreq,noiseamp,simtime,dt./10,dt,numsignals);
%% System of Equations

    function dy = SSadapt_eqs(t, y,Z)
        %% indices
        r_indexlow = 1;
        r_indexhigh = N_neurons;
        a_indexlow = r_indexhigh + 1;
        a_indexhigh = r_indexhigh + N_neurons;
        
        %% separate the input vector into its e, i, and theta components
        r = y(r_indexlow:r_indexhigh);   
        a = y(a_indexlow:a_indexhigh);
        %% THE DIFF.EQS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Z_correct = [Z(3,1);Z(3,1);Z(1,1);Z(1,1)];

        I_tot = W*r + Z_correct .* LR  - a + I_in + interp1(noiseT,Inoise,t)';

        F_I = zeros(4,1);
        F_I =  k.*(heaviside(I_tot-h).*(I_tot-h)).^n;
        Ainf = zeros(4,1);
        Ainf(1:2) = beta(1)./(1+exp(-Ak(1).*(r(1:2)-A0(1))));
        Ainf(3:4) = beta(2)./(1+exp(-Ak(2).*(r(3:4)-A0(2))));
        dr = (-r + F_I)./tau_r;
        da = (-a + Ainf)./tau_a;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% combine the outputs into one vector
        dy = [dr;da];
    end

%% RUN THE ODE SOLVER

options=odeset('RelTol',1e-12,'AbsTol',1e-12);
Y_sol = dde23(@SSadapt_eqs,[10], y0,tspan);

end

function [ X,T ] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals)
rng('shuffle');
simtimevector = 0:dt:duration;
SimTimeLength = length(simtimevector);
randnums = randn(numsignals,SimTimeLength);
savecounter = 1;
X_t = sigma.*randn(numsignals,1); %Start at random value - gaussian distribted around sigma
clear X
clear T
for tt=1:SimTimeLength
    dX = -theta.*X_t.*dt + sqrt(2.*theta).*sigma.*randnums(:,tt).*sqrt(dt);
    X_t = X_t + dX;
    
    if mod(simtimevector(tt),save_dt)==0
        X(savecounter,:) = X_t; %Pre-allocate
        T(savecounter,:) = simtimevector(tt);
        savecounter = savecounter+1;
    end
end
end