classdef ANN < handle
    
    % Feedforward neural network implementation
    % Class implementation with error backpropagation and training
    %
    % Based on the python tutorials:
    %   https://machinelearningmastery.com/implement-backpropagation-algorithm-scratch-python/
    %   http://neuralnetworksanddeeplearning.com/chap1.html
    %
    % Useage:
    %   network = ANN();
    %   network.SetArchitecture( [ input_nodes hidden_nodes{can be vector} output_nodes ] );
    %   network.SetData( input_data, output_labels );
    %   network.Initialise();
    %   predicted = network.Predict();          % <-- Untrained performance
    %   network.Train('random',1000);
    %   predicted = network.Predict();          % <-- Trained performance
    %
    % John-Stuart Brittain (University of Birmingham)
    % Created Dec 2021
    
    properties (SetAccess = private)        % Access / GetAccess / SetAccess / SetObservable / Hidden
        layer_nodes
        input_data
        output_labels
        trial_count
        weight_init_method
        activation_fcn
        cost
    
        network_type
        transfer
        transfer_deriv
        weights
        biases
        activation
        output
        errors
        delta
        learning_rate
        max_iterations
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%
        %%% Neural network setup and implementation
        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Constructor
        function this = ANN()
        end
        
        % Set network architecture (layers)
        function SetArchitecture( this, layer_nodes )
            this.layer_nodes = layer_nodes;
        end
        
        % Provide data for training
        function SetData( this, input_data, output_labels )
            if ( isempty(this.layer_nodes) )
                warning('Network architecture not defined - specify with SetArchitecture first.');
                return;
            end
            if ( this.layer_nodes(1) ~= size(input_data,1) )
                warning('Input data must match input layer size (input data is expected to be of size: layer_nodes(1) x trials )');
                return;
            end
            if ( this.layer_nodes(end) ~= size(output_labels,1) )
                warning('Output labels must match output layer size (output labels are expected to be of size: layer_nodes(end) x trials )');
                return;
            end
            if ( size(input_data,2) ~= size(output_labels,2) )
                warning('Input and output trial counts must match (data and labels are expected to be of size: n x trials )');
                return;
            end
            this.input_data = input_data;
            this.output_labels = output_labels;
            this.trial_count = size(input_data,2);
        end
        
        % Specify activation function
        function SetTransferFunction( this, choice )
            if ( isempty(choice) )
                choice = 'Sigmoid';
            end
            this.activation_fcn = choice;
            switch ( choice )
                case 'Perceptron'
                    this.transfer = @(a) a>0;
                case 'Sigmoid'
                    this.transfer = @(a) 1./(1+exp(-a));
                    this.transfer_deriv = @(a) this.transfer(a).*(1-this.transfer(a));	% Derivative of sigmoid
                case 'tanh'
                    this.transfer = @(a) tanh(a);
                    this.transfer_deriv = @(a) 1-tanh(a);
                case 'ReLU'                         % Rectified Linear Unit
                    this.transfer = @(a) ((a>0)+0).*a;
                    this.transfer_deriv = @(a) ((a>0)+0);            % Technically non-differentiable at (a=0), but eh.
                    
                    this.weight_init_method = 'he';
                case 'pReLU'                        % Parametric RLU (shallow slope <0)
                    negslope = 0.01;
                    this.transfer = @(a) ((a>0)+0).*a + ((a<=0)+0).*a*negslope;
                    this.transfer_deriv = @(a) ((a>0)+0) + ((a<=0)+0).*negslope;
                otherwise
                    this.transfer_fcn = 'Undefined';
                    warning(['Unrecognised transfer function specified --- ' choice]);
                    return;
            end
        end
        
        % Initialise network
        function Initialise(this)
            if ( isempty(this.activation_fcn) )
                this.SetTransferFunction([]);
            end
            if ( isempty(this.learning_rate) )
                this.learning_rate = 0.5;
            end
            if ( isempty(this.weight_init_method) )
                this.weight_init_method = 'xavier';
            end
            if ( isempty(this.max_iterations) )
                this.max_iterations = 500;
            end
            this.InitialiseWeights( this.weight_init_method );
        end
        
        % Initialise weights
        function InitialiseWeights(this,method)
            % Populate weights
            switch ( method )
                case 'gaussian'
                    for layer = (2:(length(this.layer_nodes)))
                        this.weights{layer} = randn(this.layer_nodes(layer),this.layer_nodes(layer-1));
                        this.biases{layer} = zeros(this.layer_nodes(layer),1);
                    end
                case 'xavier'               % Xavier Weight Initilisation (used for Sigmoid or TanH activation functions)
                    for layer = (2:(length(this.layer_nodes)))
                        n = this.layer_nodes(layer-1);                              % Number of inputs to layer
                        bound = 1/sqrt(n);                                          % Upper and lower bounds
                        this.weights{layer} = 2*bound*rand(this.layer_nodes(layer),this.layer_nodes(layer-1))-bound;  % Uniform within bounds
                        this.biases{layer} = 0.00*ones(this.layer_nodes(layer),1);	% Can be zero but consider 0.01 for asymetric activation functions, such as ReLU
                    end
                case 'normxavier'           % Normalised Xavier Weight Initilisation
                    for layer = (2:(length(this.layer_nodes)))
                        n = this.layer_nodes(layer-1);                              % Number of inputs to layer
                        m = this.layer_nodes(layer);                                % Number of outputs from the layer
                        bound = sqrt(6)/sqrt(n+m);                                  % Upper and lower bounds
                        this.weights{layer} = 2*bound*rand(this.layer_nodes(layer),this.layer_nodes(layer-1))-bound;  % Uniform within bounds
                        this.biases{layer} = 0.00*ones(this.layer_nodes(layer),1);	% Can be zero but consider 0.01 for asymetric activation functions, such as ReLU
                    end
                case 'he'                   % He Weight Initialisation (used for ReLU nodes)
                    for layer = (2:(length(this.layer_nodes)))
                        n = this.layer_nodes(layer-1);                              % Number of inputs to layer
                        sd = sqrt(2/n);
                        this.weights{layer} = sd*randn(this.layer_nodes(layer),this.layer_nodes(layer-1));
                        this.biases{layer} = 0.01*ones(this.layer_nodes(layer),1);
                    end
            end
        end
        
        % Forward propagate
        function ForwardPropagate(this)
            % Feedforward and calculate Cost
            this.output{1} = this.input_data;
            for layer = (2:length(this.layer_nodes))
                this.activation{layer} = this.weights{layer}*this.output{layer-1} + this.biases{layer};
                this.output{layer} = this.transfer( this.activation{layer} );
            end
            this.cost = mean( sum((this.output_labels - this.output{end}).^2,1), 2 );
        end
        
        % Backpropagate error
        function BackPropagateError(this)
            % Start with output layer
            this.errors{length(this.output)} = this.output{end} - this.output_labels;
            this.delta{length(this.output)} = this.errors{end}.*this.transfer_deriv(this.output{end});
            % Move to hidden layers
            for layer = ((length(this.layer_nodes)-1):-1:2)
                this.errors{layer} = zeros(size(this.output{layer}));
                % Recurse each neuron in the layer
                for neuron = (1:this.layer_nodes(layer))
                    error = zeros(1,size(this.errors{layer},2));
                    for nextlayerneuron = (1:this.layer_nodes(layer+1))
                        error = error + this.weights{layer+1}(nextlayerneuron,neuron)*this.delta{layer+1}(nextlayerneuron,:);
                    end
                    this.errors{layer}(neuron,:) = error;
                end
                this.delta{layer} = this.errors{layer}.*this.transfer_deriv(this.activation{layer});        % NB: This is 'activation' in one tutorial but 'output' in another; 'activation' is intuituvely correct, and it gives better results
            end
        end
        
        % Update weights
        function UpdateWeights(this)
            for layer = (2:length(this.layer_nodes))
                for neuron = (1:this.layer_nodes(layer))
                    % Takes sum over trials : this has been added to facilitate batch updating
                    for lastlayerneuron = (1:size(this.output{layer-1},1))
                        this.weights{layer}(neuron,lastlayerneuron) = this.weights{layer}(neuron,lastlayerneuron) - this.learning_rate * mean( this.delta{layer}(neuron,:) .* this.output{layer-1}(lastlayerneuron,:) );
                    end
                    this.biases{layer}(neuron) = this.biases{layer}(neuron) - this.learning_rate * mean( this.delta{layer}(neuron,:) );
                end
            end
        end
        
        % Predict
        function out = Predict(this)
            this.ForwardPropagate();
            out = this.output{end};
        end
        
        % Prune
        function Prune(this,threshold,make_sparse)
            % Check parameters
            if ( ~exist('make_sparse','var') )
                make_sparse = [];
            end
            if ( isempty(make_sparse) )
                make_sparse = false;
            end
            % Prune weights
            for layer = (2:length(this.layer_nodes))
                this.weights{layer}(abs(this.weights{layer})<threshold) = 0;
            end
            % Make weight matrices sparse (on option)
            if ( make_sparse )
                this.MakeSparse();
            end
        end
        
        % Make Sparse
        function MakeSparse(this)
            for layer = (2:length(this.layer_nodes))
                this.weights{layer} = sparse(this.weights{layer});
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%
        %%% Training routines
        %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Solver settings --- Max Iterations
        function SetMaxIterations( iters )
            this.max_iterations = iters;
        end
    
        % Train
        function Train(this,varargin)
            % Can specify with options
            %  Training method = varargin{1} (see TrainSGD)
            %  Training method options = varargin(2:end)
            this.TrainSGD(varargin{:});
        end
        
        % Train (stochastic gradient descent - actually, implements several versions)
        function TrainSGD(this,varargin)
            % Determine method
            if (~exist('varargin','var'))
                varargin = {};
            end
            if (isempty(varargin))
                varargin = {[]};
            end
            method = varargin{1};
            if (isempty(method))
                method = 'minibatch';
            end
            % Determine initial cost
            this.ForwardPropagate();
            % Iterative training
            this.learning_rate = 10.0;
            pc = Inf; pc_th = 0.0001;
            iter = 0; oldcost = this.cost;
            fprintf('\nIter %g cost = %03f [pc=%03f]',iter,this.cost,pc);
            while ( (abs(pc) > pc_th) || (this.learning_rate>0.0125) )
                iter = iter + 1;
                % Keep track of weights and biases incase of inversion
                w = this.weights;
                b = this.biases;
                
                % Select training method
                switch ( method )
                    case 'batch'
                        this.TrainMiniBatch(size(this.input_data,2));
                    case 'onehit'
                        varargin{2} = 1;        % One-trial at a time
                        this.TrainMiniBatch(varargin{2:end});
                    case 'minibatch'
                        this.TrainMiniBatch(varargin{2:end});
                    case 'random'
                        this.TrainRandomSampling(varargin{2:end});
                    otherwise
                        error(['Unknown training method specified - ' method]);
                end
                pc = (oldcost-this.cost)/oldcost;
                oldcost = this.cost;
                fprintf('\nIter %g cost = %03g lrate=%g [pc=%g%%]',iter,this.cost,this.learning_rate,100*pc);
                
                if ( pc < 0 )
                    % Cost increase - revert and shrink learning rate
                    fprintf('--- Reverting to old cost');
                    this.weights = w;
                    this.biases = b;
                    this.learning_rate = this.learning_rate / 2;
                    this.ForwardPropagate();
                    pc = Inf;
                    oldcost = this.cost;
                    fprintf('\nIter %g cost = %03g lrate=%g [pc=%g%%]',iter,this.cost,this.learning_rate,100*pc);
                else
                    if ( abs(pc) < pc_th )
                        this.learning_rate = this.learning_rate / 2;
                    end
                end
                
                % Single hit learning on onehitch of 1
                if ( strcmp(method,'onehit') || (iter >= this.max_iterations) )
                    break;
                end
            end
            fprintf('\nTraining complete.\n');
        end
        
        % Train in mini-batches
        function TrainMiniBatch(this,batchsize)
            % Initialise batch if not specified
            if (~exist('batchsize','var'))
                batchsize = [];
            end
            if (isempty(batchsize))
                batchsize = floor(size(this.output_labels,2)/10);
            end
            % Backup output labels
            output_labels_all = this.output_labels;
            input_data_all = this.input_data;
            % Process one at a time (can be run in batch so only provide one sample at a time)
            batches = floor(size(this.output_labels,2)/batchsize);
            for batch = (1:batches)
                trials = (1:batchsize)+(batch-1)*batchsize;
                this.SetData( input_data_all(:,trials), output_labels_all(:,trials) );
                this.ForwardPropagate();
                this.BackPropagateError();
                this.UpdateWeights();
            end
            % Restore output labels
            this.SetData( input_data_all, output_labels_all );
            % Get cost over all data
            this.ForwardPropagate();
        end
        
        % Train in mini-batches
        function TrainRandomSampling(this,batchsize)
            % Initialise batch if not specified
            if (~exist('batchsize','var'))
                batchsize = [];
            end
            if (isempty(batchsize))
                batchsize = floor(size(this.output_labels,2)/10);
            end
            % Backup output labels
            output_labels_all = this.output_labels;
            input_data_all = this.input_data;
            % Process random selection
            trials = randi(size(input_data_all,2),[1 batchsize]);
            this.SetData( input_data_all(:,trials), output_labels_all(:,trials) );
            this.ForwardPropagate();
            this.BackPropagateError();
            this.UpdateWeights();
            % Restore output labels
            this.SetData( input_data_all, output_labels_all );
            % Get cost over all data
            this.ForwardPropagate();
        end
        
    end
end
