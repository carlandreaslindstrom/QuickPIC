function [x_slope, y_slope] = qproj( project, sim, location, QS_setting)
    
    % PRODUCES A PROJECTION DOWNSTREAM, ASSUMING THE FACET SPECTROMETER
    
    % location = 'IP2A'/'BPM_3265'/'BPM_3315'/'ELAN'
    % QS_setting = {E_QS, z_ob, z_im}
    
    % define spectrometer params
    z_exit = 1994.38; % plasma exit
    z_IP2A = 1997.05; % [m] z-location of CMOS_ELAN
    z_QS1 = 1999.21; % [m] z-location of QS1 (middle)
    l_QS1 = 1; % [m] length of QS1
    z_BPM_3265 = 1999.92; % [m] z-location of DS BPM
    z_QS2 = 2004.21; % [m] z-location of QS2 (middle)
    l_QS2 = 1; % [m] length of QS2 
    z_BPM_3315 = 2004.92; % [m] z-location of Far DS BPM
    z_ELAN = 2015.22; % [m] z-location of CMOS_ELAN
    
    % transport matrix functions
    Mq = @(k,l) real([cos(sqrt(k)*l), sin(sqrt(k)*l)/sqrt(k), 0, 0; ...
                      -sin(sqrt(k)*l)*sqrt(k), cos(sqrt(k)*l), 0, 0; ...
                      0, 0, cosh(sqrt(k)*l), sinh(sqrt(k)*l)/sqrt(k); ...
                      0, 0, sinh(sqrt(k)*l)*sqrt(k), cosh(sqrt(k)*l)]);
    Md = @(s) [1, s, 0, 0; ...
               0, 1, 0, 0; ...
               0, 0, 1, s; ...
               0, 0, 0, 1];
    
    % shortcut matrices
    M_exit_IP2A = Md(z_IP2A - z_exit); % matrix from exit to QS1
    M_exit_QS1 = Md((z_QS1 - l_QS1/2) - z_exit); % matrix from exit to QS1
    M_QS1_BPM3265 = Md(z_BPM_3265 - (z_QS1 + l_QS1/2)); % matrix from QS1 to DS BPM
    M_QS1_QS2 = Md((z_QS2 - l_QS2/2) - (z_QS1 + l_QS1/2)); % matrix from QS1 to QS2
    M_QS2_BPM3315 = Md(z_BPM_3315 - (z_QS2 + l_QS2/2)); % matrix from QS2 to Far DS BPM
    M_QS2_ELAN = Md(z_ELAN - (z_QS2 + l_QS2/2)); % matrix from QS2 to CMOS_ELAN
      
    % spectrometer setting (quad strengths)
    E0 = 20.35; % [GeV] nominal energy
    if ~exist('QS_setting','var')
        E_QS = E0; % [GeV]
        z_ob = z_exit; % [m]
        z_im = z_ELAN; % [m]
    else
        E_QS = QS_setting{1}; % [GeV] absolute energy set point (focused energy)
        z_ob = QS_setting{2}; % [m] z-location of object plane
        z_im = QS_setting{3}; % [m] z-location of image plane
    end
    disp(['QS energy set point: ' num2str(E_QS, 4) ' GeV (absolute)']);
    disp(['Imaged from z = ' num2str(z_ob, 6) ' m to z = ' num2str(z_im, 6) ' m']);
    
    % find quad setting
    Rij = @(M,i,j) M(i,j);
    R12 = @(M) Rij(M, 1, 2);
    R34 = @(M) Rij(M, 3, 4);
    
    % spectrometer matrix
    M_QS = @(k1,k2) Md(z_im - (z_QS2 + l_QS2/2)) * Mq(k2, l_QS2) * M_QS1_QS2 * Mq(k1, l_QS1) * Md((z_QS1 - l_QS1/2) - z_ob);
    
    % find imaging condition based on settings
    imaging = @(k) R12(M_QS(k(1),k(2)))^2 + R34(M_QS(k(1),k(2)))^2;
    k0 = [0.3, -0.2];
    [k, ~] = fminsearch(imaging, k0, optimset('TolX',1e-6,'TolFun',1e-6));
    k_QS1 = k(1);
    k_QS2 = k(2);
    
    % NOTE: ignoring effect of beryllium/diamond and aluminium windows
    % Possible implementation: scale width of beam (but not absolute
    %                           kick) by the appropriate amount.
    
    % calculate transfer matrix
    switch location
        case 'IP2A'
            M = M_exit_IP2A;
        case 'BPM_3265'
            M = M_QS1_BPM3265 *  Mq(k_QS1, l_QS1) * M_exit_QS1;
        case 'BPM_3315'
            M = M_QS2_BPM3315 * Mq(k_QS2*E_QS/E0, l_QS2) * M_QS1_QS2 * Mq(k_QS1*E_QS/E0, l_QS1) * M_exit_QS1;
        case 'ELAN'
            M = M_QS2_ELAN * Mq(k_QS2*E_QS/E0, l_QS2) * M_QS1_QS2 * Mq(k_QS1*E_QS/E0, l_QS1) * M_exit_QS1;
    end
    
    
    
    % get phase space from simulation (final slice @ z_exit)
    addpath('..');
    outputfolder = CONFIG('outputs');
    simpath = [outputfolder '/' project '/' sim];
    rp = rpinputParser(project, sim);
    final_step = floor(rp.sim.dist.end/rp.sim.dist.step);
    phaseSpace0 = readBeams(simpath, final_step, rp);
    
    % apply fix: the laser is offset, not the beam
    x0 = rp.beam{1}.offset.x - rp.sim.dim.x/2; % [um] initial x offset
    y0 = rp.beam{1}.offset.y - rp.sim.dim.y/2; % [um] initial y offset
    phaseSpace0.X = phaseSpace0.X - x0; % [um] shift x-offset back
    phaseSpace0.Y = phaseSpace0.Y - y0; % [um] shift y-offset back
    dx_laser = -x0; % [um] laser x-offset
    dy_laser = -y0; % [um] laser y-offset
    
    % transfer the phase space
    phaseSpace = phaseSpace0;
    N = numel(phaseSpace.X);
    for i = 1:N
        v0 = [phaseSpace.X(i); phaseSpace.XP(i); phaseSpace.Y(i); phaseSpace.YP(i)];
        v0 = v0*1e-6; % scale from [um] to [m];
        v = M*v0;
        v = v*1e6; % scale from [m] back to [um];
        phaseSpace.X(i) = v(1);
        phaseSpace.XP(i) = v(2);
        phaseSpace.Y(i) = v(3);
        phaseSpace.YP(i) = v(4);
    end
    
    % find centroids
    x_centroid = mean(phaseSpace.X);
    y_centroid = mean(phaseSpace.Y);
    
    % helper functions
    bins = @(data) ceil(sqrt(numel(data)));
    
    % figure
    figure(2);
    set(gcf,'color','w');
    wbgyr = interp1([0 0.25 0.5 0.75 1],[1 1 1; 0 0 1; 0 1 0; 1 1 0; 1 0 0],linspace(0, 1, 256));
    colormap(wbgyr);
    
    % project input on screen
    [xbins0, ybins0] = deal(bins(phaseSpace0.X), bins(phaseSpace0.Y));
    [h0, axs0] = hist3([phaseSpace0.X, phaseSpace0.Y], [xbins0, ybins0]);
    subplot(2,1,1);
    imagesc(axs0{1},axs0{2}, h0');
    hold on;
    plot(dx_laser, dy_laser, 'k+','MarkerSize',20,'LineWidth',2);
    hold off;
    legend(['Laser offset: (' num2str(dx_laser,3) ', ' num2str(dy_laser,3) ') \mum']);
    set(gca,'ydir','normal');
    xlabel('x [\mum]');
    ylabel('y [\mum]');
    title('Initial XY projection (at plasma exit)');
    axis equal;
    colorbar;
    
    % project output on screen
    [xbins, ybins] = deal(bins(phaseSpace.X), bins(phaseSpace.Y));
    [h, axs] = hist3([phaseSpace.X, phaseSpace.Y], [xbins, ybins]);
    subplot(2,1,2);
    imagesc(axs{1}*1e-3,axs{2}*1e-3, h');
    hold on;
    plot(x_centroid*1e-3, y_centroid*1e-3, 'k+','MarkerSize',20,'LineWidth',2);
    hold off;
    legend(['Centroid offset: (' num2str(x_centroid*1e-3,3) ', ' num2str(y_centroid*1e-3,3) ') mm']);
    set(gca,'ydir','normal');
    xlabel('x [mm]');
    ylabel('y [mm]');
    title(['Final XY projection (at ' location ')'], 'Interpreter','None');
    axis equal;
    colorbar;
    
    % beam offset at screen per channel offset
    [x_slope, y_slope] = deal(NaN);
    if abs(dx_laser) > 0  % x
        x_slope = -x_centroid/dx_laser;
        disp(['x-slope = (channel offset/full centroid offset) = ' num2str(x_slope,3) ' (at ' location ')']);
    end
    if abs(dy_laser) > 0 % y
        y_slope = -y_centroid/dy_laser;
        disp(['y-slope = (channel offset/full centroid offset) = ' num2str(y_slope,3) ' (at ' location ')']);
    end
    
end

