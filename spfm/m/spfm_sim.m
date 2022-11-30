function sim=spfm_sim(neurons,secs,dt,init);
%function sim=spfm_sim(neurons,secs,dt,init);
%
% Sigma Pulse Frequency Modulation (SPFM) Network Simulator
%
% Input sequence to neurons may consist of independent Gaussian white noise
% and synaptic connections from other neurons.  Synaptic connections are
% modelled as pulse inputs of a specified magnitude and width.
%
% Input parameter
%       neurons     Array structure containing neuron model parameters and
%                   connectivity patterns.
%                       .tau        Time constant
%                       .v0         Reset value
%                       .th         Threshold value
%                       .mu         Mean of Gaussian white noise input
%                       .sigma      Std.dev. of Gaussian white noise input
%                       .synapse{}  Array structure of synaptic connections
%                           .neuron Post-synaptic neuron
%                           .mag    Pulse magnitude
%                           .width  Pulse width
%       secs        Duration of simulation (secs)
%       dt          Sampling interval
%       init        Startup period (secs) to allow divergence of neurons (default: 0secs)
%
% Output parameter
%       sim         Array structure containing output parameters.
%                       .neuron{}   Structure array of neuron output parameters
%                           .spikes Spike train
%                           .v      Membrane potential fluctuations
%                       .config     Structure containing input configuration
%                       .time       Time vector
%                       .runtime    Computation time of the simulation (secs)
%
% Based on the sigma encoder described in
%   D.M. Halliday (1998) "Generation and characterization of correlated
%   spike trains". Computers in Biology and Medicine (28) pp.143-152
%
%function sim=spfm_sim(neurons,secs,dt,init);

% Check input parameters
if (~exist('init'))
    init=0;
end;
try
    iscell(neurons);
    for m=1:length(neurons)
        neurons{m}.tau;
        neurons{m}.v0;
        neurons{m}.th;
        neurons{m}.mu;
        neurons{m}.sigma;
        if (isfield(neurons{m},'synapse'))
            for k=1:length(neurons{m}.synapse)
                neurons{m}.synapse{k}.neuron;
                neurons{m}.synapse{k}.mag;
                neurons{m}.synapse{k}.width;
            end;
        end;
    end;
catch
    error(' Neuron structure not of required form.');
end;

% Display progress
disp('SPFM Network Simulator');
disp('  Initialising');

% Simulation parameters
tic;
N=secs/dt;
time=[0:N]*dt;
M=length(neurons);
trange=[1:(init/dt+N)];
Nprime=length(trange);

% Initialise encoders
for m=1:M
    % Initialise output spike train
    sim.neurons{m}.spikes=[];
    % Initialise output membrane fluctuations
    sim.neurons{m}.v=zeros(Nprime+1,1);
    sim.neurons{m}.v(1)=neurons{m}.v0;
    % Initialise encoder input
    sim.neurons{m}.x=neurons{m}.mu+neurons{m}.sigma*randn(Nprime,1);
    % Ensure synapse structure present
    if (~isfield(neurons{m},'synapse'))
        neurons{m}.synapse={};
    end;
end;

% Run simulation
fprintf('\n Sim: ');
progress=0;
for n=trange       % Step through time
    % Display progress
    if (floor(100*n/Nprime)>progress)
        progress=100*n/Nprime;
        %disp(['    ' int2str(progress) '% complete.']);
        if ( mod(progress,10)==0 )
            fprintf('%g',progress);
        else
            fprintf('.');
        end;
    end;
    % Traverse neurons
    for m=1:M
        % Run encoder
        sim.neurons{m}.v(n+1)=spfm_encoder(sim.neurons{m}.v(n),sim.neurons{m}.x(n),neurons{1}.tau,dt);
        % Check for threshold crossing
        if (sim.neurons{m}.v(n+1)>neurons{m}.th)
            % Mark spike and reset membrane potential
            sim.neurons{m}.spikes=[sim.neurons{m}.spikes; (n*dt)];%+(neurons{m}.th-sim.neurons{m}.v(n))/(sim.neurons{m}.v(n+1)-sim.neurons{m}.v(n))];
            sim.neurons{m}.v(n+1)=sim.neurons{m}.v(n+1)-neurons{m}.th;
            % Recurse synapses adding pulse input to post-synaptic neurons
            for k=1:length(neurons{m}.synapse)
                nindex=neurons{m}.synapse{k}.neuron;
                mag=neurons{m}.synapse{k}.mag;
                width=neurons{m}.synapse{k}.width;
                sim.neurons{nindex}.x(n+[1:min(width/dt,N-n)])=sim.neurons{nindex}.x(n+[1:min(width/dt,N-n)])+mag;
            end;
        end;
    end;
end;

% Remove initialisation period from output
if (init>0)
    trange=(init/dt)+1+[0:N];
    for m=1:M
        sim.neurons{m}.v=sim.neurons{m}.v(trange);
        sim.neurons{m}.x=sim.neurons{m}.x(trange(1:end-1));
        sim.neurons{m}.spikes=sim.neurons{m}.spikes(sim.neurons{m}.spikes>=init)-init+dt;
    end;
end;

% Output common parameters
sim.time=time;
sim.config.secs=secs;
sim.config.dt=dt;
sim.config.neurons=neurons;

% Finish
fprintf('\n  Simulation complete.\n');
sim.runtime=toc;
