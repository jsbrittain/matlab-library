function spfm_pNeurons(sim);
%function spfm_pNeurons(sim);
%
% Sigma Pulse Frequency Modulation (SPFM) Network Simulator
%
% Plot routine to display individual membrane potential fluctuations and
% spike firing times.
%
% Input parameter
%       sim         SPFM network simulation output
%
%function spfm_pNeurons(sim);

% Plot results
figure;
unitheight=3;
M=length(sim.neurons);
disp(' ');
disp('Node statistics');
N=sim.config.secs/sim.config.dt+1;
for m=1:M
    subplot(M,1,m);
    % Plot membrane potential
    plot(sim.time,sim.neurons{m}.v,'k');
    hold on;
    % Plot marker delineating neuron spike
    if (~isempty(sim.neurons{m}.spikes))
        ypos=min(sim.neurons{m}.v)+1.1*(max(sim.neurons{m}.v)-min(sim.neurons{m}.v));
        plot(sim.neurons{m}.spikes,ypos,'k.');
    end;
    xlim(sim.time([1 end]));
    ymax=min(sim.neurons{m}.v)+1.2*(max(sim.neurons{m}.v)-min(sim.neurons{m}.v));
    ylim([min(sim.neurons{m}.v) ymax]);
    % Calculate statistics
    isi=diff(sim.neurons{m}.spikes);
    mu=1/mean(isi); sd=1/std(isi,1);
    % Display results
    title(['NEURON ' int2str(m) ': ' int2str(length(sim.neurons{m}.spikes)) ' SPIKES (\mu=' int2str(round(mu)) ', \sigma=' num2str(round(sd*10)/10) ', c.o.v=' num2str(round(mu/sd*10)/10) ')']);
    disp(sprintf('  node %2g:    mu=%4.1f, sd=%5.1f, cov=%3.1f',m,mu,sd,mu/sd));
end;
xlabel('TIME (SECS)');
