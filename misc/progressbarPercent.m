classdef progressbarPercent < handle
    
    properties (SetAccess = public)
        N
    end
    properties (SetAccess = private, Hidden)
        progress = 0
    end
    
    methods
        
        function this = progressbarPercent( N )
            this.N = N;
            start( this );
        end;
        
        function start( this )
            fprintf('<');
        end;
        
        function update( this, n )
            newprogress = round(100*n/this.N);
            decile = 10*floor(newprogress/10);
            if ( newprogress > this.progress )
                if (( decile > (10*floor(this.progress/10)) ) && (decile~=100) )
                %if ((mod(this.progress,10)==0) && (this.progress~=100))
                    fprintf('%g',decile);
                else
                    fprintf('.');
                end;
                this.progress = newprogress;
            end;
        end
        
        function stop( this )
            fprintf('>\n');
        end
        
    end
    
end
