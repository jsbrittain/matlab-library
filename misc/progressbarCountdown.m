classdef progressbarCountdown < handle
    
    properties (SetAccess = public)
        N
    end
    properties (SetAccess = private, Hidden)
        timer
        iterTime
    end
    
    methods
        
        function this = progressbarCountdown( N )
            this.N = N;
            start( this );
        end;
        
        function start( this )
            this.timer = tic;
        end;
        
        function update( this, n )
            % Report expected completion time based on first iteration
            if ( n == 1 )
                % Time for first iteration
                this.iterTime = toc( this.timer );
                [hrs,mins,secs] = this.SecsToHMS( this.iterTime );
                fprintf('  Time to complete first iteration: %02g:%02g:%02g [hh:mm:ss]\n',hrs,mins,secs);
                
                % Time for whole set
                [hrs,mins,secs] = this.SecsToHMS( this.N*this.iterTime );
                fprintf('  Estimated time to complete [N = %g]: %02g:%02g:%02g [hh:mm:ss]\n',this.N,hrs,mins,secs);
            end;
            % Time remaining
            [hrs,mins,secs] = this.SecsToHMS( (this.N-n)*this.iterTime );
            fprintf('  Time remaining: %02g:%02g:%02g [hh:mm:ss]\n',hrs,mins,secs);
        end;
        
        function stop( this )
            %
        end;
        
    end
    
    methods ( Static )
       
        function [hrs,mins,secs] = SecsToHMS( t )
            hrs = floor(t/60/60);
            mins = floor(t/60-60*hrs);
            secs = floor(t-60*60*hrs-60*mins);
        end
        
    end
    
end
