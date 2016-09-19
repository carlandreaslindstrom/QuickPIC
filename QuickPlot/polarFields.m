function [E, B, F] = polarFields(xs, ys, E, B)
    
    % XY PROJECTION
        
    % test if XY projection exists
    if isfield(E.X, 'XY') && isfield(E.Y, 'XY') && isfield(B.X, 'XY') && isfield(B.Y, 'XY')

        % sinusoids as function of position
        [XX, YY] = meshgrid(xs, ys);
        RR_XY = sqrt(XX.^2 + YY.^2);
        cos_XY = XX./RR_XY;
        sin_XY = YY./RR_XY;

        % E-fields (must be transposed to work)
        E.R.XY = E.X.XY'.*cos_XY + E.Y.XY'.*sin_XY;
        E.TH.XY = - E.X.XY'.*sin_XY + E.Y.XY'.*cos_XY;

        % B-fields
        B.R.XY = B.X.XY'.*cos_XY + B.Y.XY'.*sin_XY;
        B.TH.XY = - B.X.XY'.*sin_XY + B.Y.XY'.*cos_XY;

        % x/y forces
        F.X.XY = E.X.XY' - B.Y.XY';
        F.Y.XY = E.Y.XY' + B.X.XY';

        % radial/azimuthal forces
        F.R.XY = F.X.XY.*cos_XY + F.Y.XY.*sin_XY;
        F.TH.XY = -F.X.XY.*sin_XY + F.Y.XY.*cos_XY;

    end
    
    % XZ PROJECTION
    if isfield(E.X, 'XZ') && isfield(E.Y, 'XZ') && isfield(B.X, 'XZ') && isfield(B.Y, 'XZ')
        
        % E-fields
        E.R.XZ = E.X.XZ;
        E.TH.XZ = E.Y.XZ;

        % B-fields
        B.R.XZ = - B.X.XZ;
        B.TH.XZ = B.Y.XZ;

        % forces
        F.R.XZ = E.R.XZ - B.TH.XZ;
        F.TH.XZ = E.TH.XZ + B.R.XZ;

    end
    
    % YZ PROJECTION
    if isfield(E.X, 'YZ') && isfield(E.Y, 'YZ') && isfield(B.X, 'YZ') && isfield(B.Y, 'YZ')
        
        % E-fields
        E.R.YZ = E.Y.YZ;
        E.TH.YZ = -E.X.YZ;

        % B-fields
        B.R.YZ = -B.Y.YZ;
        B.TH.YZ = -B.X.YZ;

        % forces
        F.R.YZ = E.R.YZ - B.TH.YZ;
        F.TH.YZ = E.TH.YZ + B.R.YZ;
        
    end
end

