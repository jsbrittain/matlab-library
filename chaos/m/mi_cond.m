function I = mi_cond( pxyz )
%function I = mi_cond( pxyz )
%
% Conditional Mutual Information
%
% takes histogram
%   p(x,y,z)
% computes
%   I(X;Y|Z)
%
%function I = mi_cond( pxyz )

% Marginal probabilities
pxyz = pxyz/sum(pxyz(:));   % Ensure normalisation
pz = sum(sum(pxyz,1),2);
pxz = sum(pxyz,2);
pyz = sum(pxyz,1);

if ( true )
    % Fast way - but takes up more memory - not usually a problem as
    % histograms are already discretised
    I = sum(sum(sum( pxyz.*log( repmat(pz,[size(pxyz,1) size(pxyz,2) 1]).*pxyz./repmat(pxz,[1 size(pxyz,2) 1])./repmat(pyz,[size(pxyz,1) 1 1])) ,1),2),3);
else
    % Slow way - but works just as well
    pxz = squeeze(pxz);
    pyz = squeeze(pyz);
    I = 0;
    for x = (1:size(pxyz,1))
        for y = (1:size(pxyz,2))
            for z = (1:size(pxyz,3))
                I = I + pxyz(x,y,z)*log(pz(z)*pxyz(x,y,z)/pxz(x,z)/pyz(y,z));
            end;
        end;
    end;
end;
