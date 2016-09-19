function [] = plotOutput( foldername, step )
    
    % GET CONFIG
    addpath('..');
    outputfolder = CONFIG('outputs');
    
    % GET INPUT PARAMETERS
    rp = rpinputParser(foldername);
    Nbeams = numel(rp.beam)
    Nspecies = numel(rp.plasma)
    
    % constants
    SI_c = 299792458; % [m/s] speed of light
    SI_e = 1.60217662e-19; % [C] electron charge
    SI_me = 9.10938356e-31; % [kg] electron mass
    SI_eps0 = 8.85418782e-12; % [m^-3 kg^-1 s^4 A^2] permittivity of free space
    
    % plasma parameters
    omega_p = sqrt(rp.n0 * SI_e^2/(SI_me*SI_eps0));
    k_p = SI_c/omega_p;
    lambda_p = 2*pi/k_p;
    
    % DEFINE AXES
    xaxis = linspace(-rp.sim.dim.x/2, rp.sim.dim.x/2, 2^rp.sim.ind.x);
    yaxis = linspace(-rp.sim.dim.y/2, rp.sim.dim.y/2, 2^rp.sim.ind.y);
    zaxis = linspace(-rp.beam{1}.offset.z, rp.sim.dim.z-rp.beam{1}.offset.z, 2^rp.sim.ind.z);
    [xlb, ylb, zlb] = deal('x [\mum]', 'y [\mum]', 'z [\mum]');
    
    % INFO TEXT
    info = ['\\bfSIMULATION DETAILS\\rm \n\n' ...
            '\\bfBeam #1\\rm \n' ...
            'Energy = ' num2str(rp.beam{1}.gamma * SI_me * SI_c^2 / (1e9 * SI_e),4) ' GeV \n' ...
            'Gamma = ' num2str(rp.beam{1}.gamma,5) ' \n' ...
            'Number of particles = ' num2str(rp.beam{1}.N,3) ' \n' ...
            'Charge = ' num2str(rp.beam{1}.N * rp.beam{1}.charge * SI_e * 1e9,4) ' nC \n\n' ...
            '\\bfPlasma\\rm \n' ...
            'Density = ' num2str(rp.n0*1e-6,4) ' cm^{-3} \n' ...
            'Wavelength = ' num2str(lambda_p*1e6,4) ' um \n' ...
            'Wavenumber = ' num2str(k_p,6) ' m^{-1} \n' ...
            'Frequency = ' num2str(omega_p*1e-12,4) ' THz'];
        
    
    % EXTRACTING FIELDS AND DENSITIES
    
    % shortcut for extraction
    geth = @(var) h5read([outputfolder '/' foldername '/' var '/' var '_' num2str(step,'%04d') '.h5']);
    
    % extract E-FIELDS projections
    for proj = {'XY','XZ','YZ'}
        % fields
        for field = {'E','B'}
            for dim = {'X','Y','Z'}
                qp.(field{:}).(dim{:}).(proj{:}) = geth(['F' field{:} dim{:} '-' proj{:}]);
            end
        end
        % plasmas
        for i = 1:Nspecies
            Q = geth(['QEP' num2str(i) '-' proj{:}]);
            if i == 1 % first species
                qp.Q.P.(proj{:}) = Q;
            else % otherwise add
                qp.Q.P.(proj{:}) = qp.Q.P.(proj{:}) + Q;
            end
        end
        % beams
        qp.Q.B.(proj{:}) = geth(['QEB-' proj{:}]);
    end
    
    %{
    [qp.E.X.XY, qp.E.Y.XY, qp.E.Z.XY] = deal(geth('FEX-XY')', geth('FEY-XY')', geth('FEZ-XY')');
    [qp.E.X.XZ, qp.E.Y.XZ, qp.E.Z.XZ] = deal(geth('FEX-XZ'), geth('FEY-XZ'), geth('FEZ-XZ'));
    [qp.E.X.YZ, qp.E.Y.YZ, qp.E.Z.YZ] = deal(geth('FEX-YZ'), geth('FEY-YZ'), geth('FEZ-YZ'));
    
    % extract B-FIELDS projections
    [qp.B.X.XY, qp.B.Y.XY, qp.B.Z.XY] = deal(geth('FBX-XY')', geth('FBY-XY')', geth('FBZ-XY')');
    [qp.B.X.XZ, qp.B.Y.XZ, qp.B.Z.XZ] = deal(geth('FBX-XZ'), geth('FBY-XZ'), geth('FBZ-XZ'));
    [qp.B.X.YZ, qp.B.Y.YZ, qp.B.Z.YZ] = deal(geth('FBX-YZ'), geth('FBY-YZ'), geth('FBZ-YZ'));
    
    % extract BEAM density projections
    [qp.Q.B.XY, qp.Q.B.XZ, qp.Q.B.YZ] = deal(geth('QEB-XY'), geth('QEB-XZ'), geth('QEB-YZ'));
    
    % extract PLASMA density projections
    Nspecies = 0;
    while true
        if exist([outputfolder '/' foldername '/QEP' num2str(Nspecies+1) '-XY'],'file') == 7
            Nspecies =  Nspecies + 1;
            continue;
        end
        break;
    end
    [qp.Q.P.XY, qp.Q.P.XZ, qp.Q.P.YZ] = deal(zeros(size(qp.Q.B.XY)), zeros(size(qp.Q.B.XZ)), zeros(size(qp.Q.B.YZ)));
    for i = 1:Nspecies
        qp.Q.P.XY = qp.Q.P.XY + geth(['QEP' num2str(i) '-XY']);
        qp.Q.P.XZ = qp.Q.P.XZ + geth(['QEP' num2str(i) '-XZ']);
        qp.Q.P.YZ = qp.Q.P.YZ + geth(['QEP' num2str(i) '-YZ']);
    end
    %}
    
    % ADD POLAR FIELDS
    [qp.E, qp.B, qp.F] = polarFields2(xaxis, yaxis, qp.E, qp.B); 
    
    % PLOTTING
    set(0,'DefaultFigureWindowStyle','docked');
    
    % set figure
    function newFigure(num, name)    
        fig = figure(num);
        clf(fig);
        set(gcf,'color','w');
        left_color = [0.5,0.5,0.5];
        right_color = [0,0,0];
        set(fig,'defaultAxesColorOrder',[left_color; right_color]);
        set(gcf,'name',name,'numbertitle','off')
    end

    function fieldPlot(row, col, num, field, ttl, xlbl, ylbl, xs, ys, lineout)
        % projection plot
        ax = subplot(row,col,num); 
        switch ttl(1)
            case 'E'
                clbl = 'E-field [GV/m]';
                scaling = SI_me*SI_c*omega_p/(SI_e*1e9); % norm to GV/m
            case 'B'
                clbl = 'B-field [T]';
                scaling = SI_me*omega_p/SI_e; % norm to T
            case 'F'
                clbl = 'Force [GeV/m]';
                scaling = SI_me*SI_c*omega_p/(SI_e*1e9); % norm to GeV/m
            otherwise
                clbl = '';
                scaling = 1;
        end
        field = scaling*field;
        imagesc(xs, ys, field);
        set(gca,'ydir','normal');
        c = colorbar;
        colormap(ax, cmap);
        cex = max(max(max(field)),-min(min(field)));
        caxis([-cex,cex]);
        c.Label.String = clbl;
        xlabel(xlbl);
        ylabel(ylbl);
        title(ttl);
        xlim([min(xs),max(xs)]);
        ylim([min(ys),max(ys)]);
        
        
        % horizontal lineout
        if (nargin > 9) && (lineout==true)
            yout = 0;
            line(get(gca,'xlim'), yout*[1,1],'LineStyle',':','Color',[0.5,0.5,0.5]);
            yyaxis right;
            [~, index] = min(abs(ys-yout));
            plot(xs, field(index,:),'-');
        end
    end

    function particlePlot(row, col, num, density, ttl, xlbl, ylbl, xs, ys, lineout)
        ax = subplot(row,col,num);
        scaling = rp.n0*1e-6;
        imagesc(xs, ys, density*scaling);
        set(gca,'ydir','normal');
        c = colorbar;
        cmin = scaling*min([min(min(qp.Q.B.XY+qp.Q.P.XY)), min(min(qp.Q.B.XZ+qp.Q.P.XZ)), min(min(qp.Q.B.YZ+qp.Q.P.YZ))]);
        cmax = scaling*max([max(max(qp.Q.B.XY+qp.Q.P.XY)), max(max(qp.Q.B.XZ+qp.Q.P.XZ)), max(max(qp.Q.B.YZ+qp.Q.P.YZ))]);
        cmin = min(0,max(-5*scaling,cmin));
        cmax = min(5*scaling, max(cmax,0));
        colormap(ax, cmap(cmin,cmax));
        caxis([cmin,cmax]);
        c.Label.String = 'Density [cm^{-3}]';
        xlabel(xlbl);
        ylabel(ylbl);
        title(ttl);
        
        % horizontal lineout
        if (nargin > 9) && (lineout==true)
            yout = 0;
            line(get(gca,'xlim'), yout*[1,1],'LineStyle',':','Color',[0.5,0.5,0.5]);
            yyaxis right;
            [~, index] = min(abs(ys-yout));
            plot(xs, density(index,:),'-');
        end
    end

    function textPlot(row, col, num, string)
        ax = subplot(row, col, num); 
        cla(ax)
        set(ax, 'visible', 'off');
        text(0, 0.5, sprintf(string)); 
    end

    
    newFigure(1, 'Beam/Plasma');
    particlePlot(2, 2, 1, (qp.Q.B.XZ+qp.Q.P.XZ)', 'Q_{beam} + Q_{plasma}', xlb, zlb, xaxis, zaxis);
    line(get(gca,'xlim'),(rp.dump.slice.z+min(zaxis))*[1,1]);
    particlePlot(2, 2, 3, qp.Q.B.XY+qp.Q.P.XY, 'Q_{beam} + Q_{plasma}', xlb, ylb, xaxis, yaxis);
    line(get(gca,'xlim'),(rp.dump.slice.y+min(yaxis))*[1,1]);
    line((rp.dump.slice.x+min(xaxis))*[1,1],get(gca,'ylim'));
    particlePlot(2, 2, 4, qp.Q.B.YZ+qp.Q.P.YZ, 'Q_{beam} + Q_{plasma}', zlb, ylb, zaxis, yaxis);
    line((rp.dump.slice.z+min(zaxis))*[1,1],get(gca,'ylim'));
    textPlot(2, 2, 2, info);
    
    newFigure(2, 'Fields (XY)'); % XY projection
    [nx, ny] = deal(3,2);
    fieldPlot(nx, ny, 1, qp.E.R.XY, 'E_r', xlb, ylb, xaxis, yaxis);
    fieldPlot(nx, ny, 2, qp.B.R.XY, 'B_r', xlb, ylb, xaxis, yaxis);
    fieldPlot(nx, ny, 3, qp.E.TH.XY, 'E_{\theta}', xlb, ylb, xaxis, yaxis);
    fieldPlot(nx, ny, 4, qp.B.TH.XY, 'B_{\theta}', xlb, ylb, xaxis, yaxis);
    fieldPlot(nx, ny, 5, qp.E.Z.XY, 'E_z', xlb, ylb, xaxis, yaxis);
    fieldPlot(nx, ny, 6, qp.B.Z.XY, 'B_z', xlb, ylb, xaxis, yaxis);
    
    newFigure(3, 'Fields (XZ)'); % XZ projection
    [nx, ny] = deal(3,2);
    fieldPlot(nx, ny, 1, qp.E.X.XZ, 'E_x', zlb, xlb, zaxis, xaxis);
    fieldPlot(nx, ny, 2, qp.B.X.XZ, 'B_x', zlb, xlb, zaxis, xaxis);
    fieldPlot(nx, ny, 3, qp.E.Y.XZ, 'E_y', zlb, xlb, zaxis, xaxis);
    fieldPlot(nx, ny, 4, qp.B.Y.XZ, 'B_y', zlb, xlb, zaxis, xaxis);
    fieldPlot(nx, ny, 5, qp.E.Z.XZ, 'E_z', zlb, xlb, zaxis, xaxis, true);
    fieldPlot(nx, ny, 6, qp.B.Z.XZ, 'B_z', zlb, xlb, zaxis, xaxis);
    
    newFigure(4, 'Fields (YZ)'); % YZ projection
    fieldPlot(nx, ny, 1, qp.E.X.YZ, 'E_x', zlb, ylb, zaxis, yaxis);
    fieldPlot(nx, ny, 2, qp.B.X.YZ, 'B_x', zlb, ylb, zaxis, yaxis);
    fieldPlot(nx, ny, 3, qp.E.Y.YZ, 'E_y', zlb, ylb, zaxis, yaxis);
    fieldPlot(nx, ny, 4, qp.B.Y.YZ, 'B_y', zlb, ylb, zaxis, yaxis);
    fieldPlot(nx, ny, 5, qp.E.Z.YZ, 'E_z', zlb, ylb, zaxis, yaxis, true);
    fieldPlot(nx, ny, 6, qp.B.Z.YZ, 'B_z', zlb, ylb, zaxis, yaxis);
    
    newFigure(5, 'Wake forces'); % YZ projection
    fieldPlot(nx, ny, 1, qp.F.R.XY, 'F_r', xlb, ylb, xaxis, yaxis);
    fieldPlot(nx, ny, 2, qp.F.TH.XY, 'F_{\theta}', xlb, ylb, xaxis, yaxis);
    fieldPlot(nx, ny, 3, qp.F.R.XZ, 'F_r', zlb, xlb, zaxis, xaxis);
    fieldPlot(nx, ny, 4, qp.F.TH.XZ, 'F_{\theta}', zlb, xlb, zaxis, xaxis);
    fieldPlot(nx, ny, 5, qp.F.R.YZ, 'F_r', zlb, ylb, zaxis, yaxis);
    fieldPlot(nx, ny, 6, qp.F.TH.YZ, 'F_{\theta}', zlb, ylb, zaxis, yaxis);
    
    newFigure(6, 'Summary'); % YZ projection
    particlePlot(3, 2, 1, (qp.Q.B.XY+qp.Q.P.XY)', 'Q_{beam} + Q_{plasma}', xlb, ylb, xaxis, yaxis);
    particlePlot(3, 2, 2, qp.Q.B.YZ+qp.Q.P.YZ, 'Q_{beam} + Q_{plasma}', zlb, ylb, zaxis, yaxis);
    fieldPlot(3, 2, 3, qp.F.R.XY, 'F_r', xlb, ylb, xaxis, yaxis);
    fieldPlot(3, 2, 4, qp.F.R.YZ, 'F_r', zlb, ylb, zaxis, yaxis);
    fieldPlot(3, 2, 6, qp.E.Z.YZ, 'E_z', zlb, ylb, zaxis, yaxis, true);
    textPlot(3, 2, 5, info);
    
    figure(18);
    set(gcf,'color','w');
    xdump = (rp.dump.slice.x+min(xaxis));
    ydump = (rp.dump.slice.y+min(yaxis));
    zdump = (rp.dump.slice.z+min(zaxis));
    [~, xind] = min(abs(xaxis-xdump));
    [~, yind] = min(abs(yaxis-ydump));
    [~, zind] = min(abs(zaxis-zdump));
    xslice = xaxis(xind);
    yslice = yaxis(yind);
    zslice = zaxis(zind);
    %q = zeros(numel(xaxis), numel(yaxis), numel(zaxis));
    %q(xind,:,:) = qp.Q.B.YZ+qp.Q.P.YZ;
    %q(:,yind,:) = qp.Q.B.XZ+qp.Q.P.XZ;
    %q(:,:,zind) = qp.Q.B.XY+qp.Q.P.XY;
    %h = slice(xaxis,yaxis,zaxis,q,xslice,yslice,zslice);
    %xlabel('x [\mum]'); ylabel('y [\mum]'); zlabel('z [\mum]');
    
    Q = zeros(numel(zaxis), numel(xaxis), numel(yaxis));
    Q(xind,:,:) = (qp.Q.B.YZ+qp.Q.P.YZ)';
    Q(:,zind,:) = (qp.Q.B.XY+qp.Q.P.XY)';
    Q(:,:,yind) = qp.Q.B.XZ+qp.Q.P.XZ;
    h = slice(zaxis,xaxis,yaxis,Q,zslice,xslice,yslice);
    xlabel('z [\mum]'); ylabel('x [\mum]'); zlabel('y [\mum]');
    set(h, 'EdgeColor','none');
    colormap(cmap(-3,0.5));
    caxis([-3,0.5]);
    colorbar;
    axis tight;
    axis equal;
    
    figure(19);
    set(gcf,'color','w');
    Fr = zeros(numel(zaxis), numel(xaxis), numel(yaxis));
    Fr(xind,:,:) = qp.F.R.YZ';
    Fr(:,zind,:) = qp.F.R.XY';
    Fr(:,:,yind) = qp.F.R.XZ;
    h = slice(zaxis,xaxis,yaxis,Fr,zslice,xslice,yslice);
    xlabel('z [\mum]'); ylabel('x [\mum]'); zlabel('y [\mum]');
    set(h, 'EdgeColor','none');
    colormap(cmap(-0.5,0.5));
    caxis([-0.5,0.5]);
    colorbar;
    axis tight;
    axis equal;
    
    figure(20);
    set(gcf,'color','w');
    Ez = zeros(numel(zaxis), numel(xaxis), numel(yaxis));
    Ez(xind,:,:) = qp.E.Z.YZ';
    Ez(:,zind,:) = qp.E.Z.XY';
    Ez(:,:,yind) = qp.E.Z.XZ;
    h = slice(zaxis,xaxis,yaxis,Ez,zslice,xslice,yslice);
    xlabel('z [\mum]'); ylabel('x [\mum]'); zlabel('y [\mum]');
    set(h, 'EdgeColor','none');
    colormap(cmap(-0.5,0.5));
    caxis([-0.5,0.5]);
    colorbar;
    axis tight;
    axis equal;
    %cmin = min([min(min(qp.Q.B.XY+qp.Q.P.XY)), min(min(qp.Q.B.XZ+qp.Q.P.XZ)), min(min(qp.Q.B.YZ+qp.Q.P.YZ))]);
    %cmax = max([max(max(qp.Q.B.XY+qp.Q.P.XY)), max(max(qp.Q.B.XZ+qp.Q.P.XZ)), max(max(qp.Q.B.YZ+qp.Q.P.YZ))]);
    %cmin = min(0,max(-5,cmin));
    %cmax = min(5, max(cmax,0));
    
    % TODO: make this 3D slice plot a function
    
end

