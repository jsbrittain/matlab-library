% HMM model testing

%% Load / generate data

dataset=0;
switch (dataset)
    % Simulated sinusoidal HMM process (Penny & Roberts; Eqs.23-24)
    case 0,
        % Define states and transition probabilities
        prior=[1 1 1]/3;
        A=[0.98 0.02 0.00; 0.00 0.98 0.02; 0.02 0.00 0.98];
        B=[0.9 0.05 0.05; 0.1 0.9 0; 0 0.9 0.1];
        % Simulate data
        rate=125; N=200*rate;
        [y,state] = hmm_generate(prior,A,B,N);
    case 1,
        % Random initialisation
        prior=hmm_rownorm(rand(1,3));
        A=hmm_rownorm(rand(3));
        B=hmm_rownorm(rand(3));
        % Simulate data
        rate=125; N=40*rate;
        [y,state] = hmm_generate(prior,A,B,N);
end;

% Assumptions
nstates = 3;

% Estimation parameters
method = 1;

%% Probability of observation given the model

[logp,alphati,betati] = hmm_probobs(y,prior,A,B,method);

%% Viterbi algorithm (estimate underlying states)

[q,qlogp] = hmm_viterbi(y,prior,A,B);
accuracy = length(find(q==state))/length(q);

%% State estimation

switch (0)
    case 0,         % Random initialisation
        prior0=hmm_rownorm(rand(size(prior)));
        A0=hmm_rownorm(rand(size(A)));
        B0=B;%hmm_rownorm(rand(size(B)));
    case 1,         % Uniform initialisation
        prior0=hmm_rownorm(ones(size(prior)));
        A0=hmm_rownorm(ones(size(A)));
        B0=B;%hmm_rownorm(ones(size(B)));
    case 2,         % Correct solution
        prior0 = prior; A0 = A; B0 = B;
end;
[priorhat,Ahat,Bhat,gammati,epsilontij] = hmm_estimatemodel(y,prior0,A0,B0,method);

%% Plot results

time=(1:N)/rate;
subplot(2,1,1);
    h=plotyy(time,state,time,y); set(h,'ylim',[1 nstates]+[-0.1 0.1]);
    set(get(h(1),'children'),'linewidth',3);
    title(['Observation / state (probability log p=' num2str(round(10*logp)/10) ')']);
subplot(2,1,2);
    h=plotyy(time,state,time,q); set(h,'ylim',[1 nstates]+[-0.1 0.1]);
    set(get(h(1),'children'),'linewidth',3);
    set(get(h(2),'children'),'color',[1 0 0]);
    title(['Actual / estimated states (accuracy ' num2str(round(1000*accuracy)/10) '%)']);
