classdef progressbar < handle
    % This is a wrapper class for several progressbar implementations that
    % follow a common template
    
    properties ( SetAccess = private, Hidden )
        pbar
    end
    
    methods
        
        function this = progressbar( N, mode )
            
            % Check for empty call --- done when children created
            if ( ~exist('N','var') )
                return;
            end;
            
            % Determine mode
            default_mode = 'percent';
            if (~exist( 'mode', 'var' ))
                mode = [];
            end;
            if (isempty(mode))
                mode = default_mode;
            end;
            
            % Return instance of relevant progressbar
            switch ( mode )
                case 'percent',   this.pbar = progressbarPercent( N );
                case 'countdown', this.pbar = progressbarCountdown( N );
                otherwise, error('Unknown progressbar mode (%s) specified.',mode);
            end;
            
        end;
        
        function update( this, n )
            this.pbar.update( n );
        end;
        
        function stop( this )
            this.pbar.stop();
        end
        
    end
    
end
