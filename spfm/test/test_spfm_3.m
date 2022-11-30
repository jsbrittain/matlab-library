% Example 2 - Halliday 1998

% Simulation parameters
secs=100;
dt=1e-3;
init=2;

% Common encoder parameters
tau=2.5e-2;
v0=0;
th=1;

% Common 1 (10Hz, COV=0.1)
neurons{1}.tau=tau;                 % Time-constant
neurons{1}.v0=v0;                   % Reset value
neurons{1}.th=th;                   % Threshold value
neurons{1}.mu=1.02;                 % Noise statistics N(mu,sigma^2)
neurons{1}.sigma=0.065;             %  |
neurons{1}.synapse{1}.neuron=3;     % Post-synaptic neuron
neurons{1}.synapse{1}.mag=6.63;     % Pulse magnitude
neurons{1}.synapse{1}.width=0.002;  % Pulse width (secs)
neurons{1}.synapse{2}.neuron=4;     % Post-synaptic neuron
neurons{1}.synapse{2}.mag=6.63;     % Pulse magnitude
neurons{1}.synapse{2}.width=0.002;  % Pulse width (secs)

% Common 2 (25Hz, COV=0.1)
neurons{2}.tau=tau;
neurons{2}.v0=v0;
neurons{2}.th=th;
neurons{2}.mu=1.269;
neurons{2}.sigma=0.307;
neurons{2}.synapse{1}.neuron=3;
neurons{2}.synapse{1}.mag=3.98;
neurons{2}.synapse{1}.width=0.002;
neurons{2}.synapse{2}.neuron=4;
neurons{2}.synapse{2}.mag=3.98;
neurons{2}.synapse{2}.width=0.002;

% Neuron 1
neurons{3}.tau=tau;
neurons{3}.v0=v0;
neurons{3}.th=th;
neurons{3}.mu=0.523;
neurons{3}.sigma=6.30;

% Neuron 2
neurons{4}.tau=tau;
neurons{4}.v0=v0;
neurons{4}.th=th;
neurons{4}.mu=0.523;
neurons{4}.sigma=6.30;

% Run simulation
sim=spfm_sim(neurons,secs,dt,init);

% Plot results
spfm_pNeurons(sim);

% Spectral analysis
sptrain3=round(sim.neurons{3}.spikes/dt);
sptrain4=round(sim.neurons{4}.spikes/dt);

% Single pair of spike trains
seg_pwr=10; rate=1/dt; opt_str=' ';
[f1(:,:,1),t1(:,:,1),cl1(1),sp1(:,:,1)]=sp2_ma(sptrain3,sptrain4,rate,seg_pwr,opt_str);
[plf1,plv1]=f_pool(sp1(:,:,1),cl1(1)); cl1(1).what='Single pair';
% Plot single segment analysis
figure; psp2(f1(:,:,1),t1(:,:,1),cl1(1),100,80,40,1,'');

% Pooled spectra
L=20;
for ind=2:L
    % Display progress
    disp([' Pair set ' int2str(ind) ' of ' int2str(L)]);
    % Run simulation
    sim=spfm_sim(neurons,secs,dt,init);
    sptrain3=round(sim.neurons{3}.spikes/dt);
    sptrain4=round(sim.neurons{4}.spikes/dt);
    % Perform single pair analysis
	[f1(:,:,ind),t1(:,:,ind),cl1(ind),sp1(:,:,ind)]=sp2_ma(sptrain3,sptrain4,rate,seg_pwr,opt_str);
    % Add results to pool
	[plf1,plv1]=f_pool(sp1(:,:,ind),cl1(ind),plf1,plv1);
end
% Process pooled spectral coefficients
[f2,t2,cl2]=pl_fout(plf1,plv1); cl1(1).what='Pooled';
% Plot results
figure; psp2(f2,t2,cl2,100,80,40,1,'');
