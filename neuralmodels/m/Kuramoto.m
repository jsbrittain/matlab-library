%%%
%%% Be aware that this class is defined using a handle implementation
%%%  and therefore operates via pass-by-reference (as all good classes
%%%  should!). This means operations are internally persistent, unlike
%%%  pass-by-value classes encountered in Matlab, where updated class
%%%  objects must be returned and updated explicitely in the calling code.
%%%

classdef Kuramoto < handle
    
    properties (SetAccess = public)
        omega               % Frequencies (radians)
        dt                  % Sampling rate of simulation
        theta_init          % Starting phases
        theta               % Current phase of all nodes
        K                   % Connectivity matrix: K(to,from) = c;
        spikeK              % Spiking connectivity matrix
        spike_rate          % Spiking rate ( [] = off )
        spiking             % Vector of instantaneous spikes
        noise_sd            % Noise SD per node
        param               % List of model parametrs (useful for optimisation)
    end
    properties (SetAccess = private, Hidden)
        % Maintain a random stream whose seed can be reset to force
        % reproducilibility; this approach avoids affecting Matlab's global
        % rng() function that would affect, for instance, probabilistic
        % solvers that call this function!
        rndStream = RandStream( 'twister' );
    end
    
    methods
        
        % Constructor
        function obj = Kuramoto(cols,rows)
            if ( nargin ~= 0 )
                if (~exist('rows','var'))
                    rows = 1;
                end;
                obj(cols,rows) = Kuramoto;
            end;
        end
        
        % Set nodeal frequencies (Hz)
        function setFreqs( this, freqs )
            this.omega = freqs*2*pi;
        end
        
        % Get nodeal frequencies (Hz)
        function freqs = getFreqs( this )
            freqs = this.omega/2/pi;
        end
        
        % Set nodal frequencies (radians)
        function set.omega( this, value )
            this.omega = value;
        end
        
        % Set node noise levels
        function set.noise_sd( this, value )
            if (isempty(this.omega))
                error('No nodes present in network; set frequencies first.');
            end;
            if ( isscalar( value ) )
                this.noise_sd = value*ones(1,this.M());
            elseif ( length(value) == this.M() )
                this.noise_sd = value;
            else
                error('Improper input.');
            end;
        end
        
        % Set node connectivity matrix
        function set.K( this, value )
            if (isempty(this.omega))
                error('No nodes present in network; set frequencies first.');
            end;
            if ( isscalar( value ) )
                this.K = value*ones(this.M());
            elseif (all( size(value) == this.M() ) && (ismatrix(value)))
                this.K = value;
            else
                error('Improper input.');
            end;
        end
        
        % Set initial phase of nodes
        function set.theta( this, value )
            if (isempty(this.omega))
                error('No nodes present in network; set frequencies first.');
            end;
            if ( isscalar( value ) )
                this.theta = value*ones(1,this.M());
            elseif ( length(value) == this.M() )
                this.theta = value;
            else
                error('Improper input.');
            end;
        end
        
        % Get number of nodes
        function len = M( this )
            len = length( this.omega );
        end
        
        % Initialise simulation
        function Initialise( this )
            if (isempty( this.theta ))
                this.theta = 0;
            end;
            if (isempty( this.K ))
                this.K = 0;
            end;
            if (isempty( this.spikeK ))
                this.spikeK = zeros(this.M);
            end;
            if (isempty( this.spike_rate ))
                this.spike_rate = zeros(1,this.M);
            end;
            if (isempty( this.noise_sd ))
                this.noise_sd = 0;
            end;
        end
        
        function resetPhases( this )
            this.theta = this.theta_init;
        end
        
        function resetRandomSeed( this, val )
            if (~exist('val','var'))
                val = 0;
            end;
            reset(this.rndStream,val);
        end
        
        % Perform simulation run over specified time period (does not update)
        function [thetaT,spikeT] = run( this, T )
            N = floor( T/this.dt );
            thetaT = nan(N,this.M());
            spikeT = false(N,this.M());
            for n = (0:(N-1))
                this.step();
                thetaT(n+1,:) = this.theta;
                if ( ~isempty(this.spiking) )
                    spikeT(n+1,:) = this.spiking;
                end;
            end
        end
        
        % Simulation step (updates phase vector "theta" in object)
        function step( this )
            % Integration implemented by Euler-Maruyama method,
            %  a simple generalisation of Euler's method for ODEs. The
            %  major difference is the separation of drift and diffusion
            %  elements, scaling the noise process by sqrt(dt) instead of dt.
            % Since SDE integration schemes are not equivalent as dt->0,
            %  it may be necessary to implement further schemes in the
            %  future...
            
            % Spiking probability
            if ( ~isempty( this.spikeK ) )
                % Firing probability envelope
                switch ( 2 )
                    case 1, penv = (2*cos(this.theta))-1;                   % Cosine as 0->1
                    case 2, penv = 2*cos(this.theta); penv(penv>1) = 1;     % Flattened cosine
                end;
                pspike = penv.* ( this.spike_rate*this.dt );    % Poisson process
                this.spiking = ( rand(1,this.M) < pspike );
                
                spike_cont = (double(this.spiking)*this.spikeK).*sin(this.theta);
            else
                spike_cont = 0;
            end;
            
            % Phase update
            this.theta = angle(exp(1i*( ...
                                this.theta + ...
                      this.dt*( this.omega + ...
                                sum(this.K.*sin(repmat(this.theta,this.M,1)-repmat(this.theta',1,this.M)),2).' + ...
                                spike_cont ) + ...
                sqrt(this.dt)*( this.noise_sd.*randn(this.rndStream,1,this.M) ) )));
        end
        
        % Report vectorised parameters
        function x = vecParams( this )
            x = [ this.omega(:); this.K(:); this.noise_sd(:) ];
        end
        
        % Un-vectorise parameters
        function unvecParams( this, x )
            M = sqrt(1+length(x))-1;
            m = 0;
            this.omega    = x( 1:M ).'; m = m + M;
            this.K        = reshape( x( m+(1:M^2) ), M, M ); m = m + M^2;
            this.noise_sd = x( m+(1:M) ).';
        end
        
        % Create new copy of the object (required as Handle class)
        function newobj = copy( this )
            newobj = Kuramoto();
            p = properties(this);
            for i = 1:length(p)
                newobj.(p{i}) = this.(p{i});
            end
        end
        
        % Add parameter to model --- this is useful for optimisation
        function addParameter( this, field, sub, value, lb, ub, name )
            if (~exist('name','var')), name = ['Param_' num2str(length(this.param)+1)]; end;
            if (~exist('ub','var')),   ub   = []; end;
            if (~exist('lb','var')),   lb   = []; end;
            if ( lb == ub )
                % No search space, so not a parameter -- reject!
                warning(['Rejecting parameter: ' name ' as lower_bound = upper_bound = ' num2str(lb)]);
                return;
            end;
            this.param(end+1).name  = name;
            this.param(end).field   = field;
            this.param(end).sub     = sub;
            this.param(end).ind     = cellfun( @(x) sub2ind(size(this.(field)),x{:}), cellfun( @num2cell, num2cell(sub,2), 'uni', 0 ), 'uni', 1 );
            this.param(end).value   = value;
            this.param(end).lb      = lb;
            this.param(end).ub      = ub;
        end
        
        % Get parameter vector
        function [x,bnds] = getParameters( this )
            x = nan(length(this.param),1);
            bnds = nan(length(this.param),2);
            for k = (1:length(x))
                x(k) = this.param(k).value;
                bnds(k,1) = this.param(k).lb;
                bnds(k,2) = this.param(k).ub;
            end;
        end
        
        % Update parameters
        function updateParameters( this, x )
            assert( length(x) == length(this.param), 'Inconsistent parameters list!' );
            for k = (1:length(x))
                this.param(k).value = x(k);
            end;
            unpackParameters( this );
        end
        
        % Unpack parameters
        function unpackParameters( this )
            for k = (1:length(this.param))
                this.(this.param(k).field)(this.param(k).ind) = this.param(k).value;
            end;
        end
        
        % Get parameters
        function param = getParamStruct( this )
            param = this.param;
        end
        
        % Replace parameters
        function replaceParamStruct( this, param )
            this.param = param;
        end
    end
end
