% Neuron producing output to second neuron (one-way)

% Simulation parameters
secs=100;
dt=1e-3;
init=2;

% Common encoder parameters
tau=2.5e-2;
v0=0;
th=1;

% Neuron 1
neurons{1}.tau=tau;                 % Time-constant
neurons{1}.v0=v0;                   % Reset value
neurons{1}.th=th;                   % Threshold value
neurons{1}.mu=1.02;                 % Noise statistics N(mu,sigma^2)
neurons{1}.sigma=0.065;             %  |

% Neuron 2
neurons{2}.tau=tau;
neurons{2}.v0=v0;
neurons{2}.th=th;
neurons{2}.mu=1.015;
neurons{2}.sigma=0.15;

% Neuron 3
neurons{3}.tau=tau;
neurons{3}.v0=v0;
neurons{3}.th=th;
neurons{3}.mu=1.269;
neurons{3}.sigma=0.307;

% Neuron 4
neurons{4}.tau=tau;
neurons{4}.v0=v0;
neurons{4}.th=th;
neurons{4}.mu=0.892;
neurons{4}.sigma=6.20;

% Run simulation
sim=spfm_sim(neurons,secs,dt,init);

% Plot results
spfm_pNeurons(sim);
