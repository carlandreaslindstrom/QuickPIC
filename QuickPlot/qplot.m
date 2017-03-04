function [] = qplot(project, sim, steps, data, ROIs)

    %% PLOTS QuickPIC SIMULATIONS

    % USAGE:
    % qplot(project, sim, steps, {{field, proj}, {...}}, {xlims, ylims, zlims})
    
    % field: 'EX', 'EY', 'EZ', 'ER', 'ETH' for E-FIELDS in (x/y/z/r/theta)
    %         'BX', 'BY', 'BZ', 'BR', 'BTH' for B-FIELDS in (x/y/z/r/theta)
    %         'FX', 'FY', 'FZ', 'FR', 'FTH' for FORCES in (x/y/z/r/theta)
    %         'QB' for BEAM DENSITY (all beams)
    %         'QP' for PLASMA DENSITY (all plasmas)
    %         'Q' for PLASMA + DEAM DENSITY
    % phase_space: 'A_B' where A, B can be (X/Y/Z/XP/YP/E) => underscore gives phase space
    
    % proj: 'XY', 'XZ', 'YZ' for slice plots 
    %       'XYZ' for 3D combo slice plot
    
    % steps: first:last
    %     or {first:last, delay}
    %     or {first:last, delay, loop?}
    
    % ROIs: {xmin:xmax, ymin:ymax, zmin:zmax, xpmin:xpmax, ypmin:ypmax, Emin:Emax}
    %       [] or blank gives default: full ROI
    
    
    %% CONSTANTS
    SI_c = 299792458; % [m/s] speed of light
    SI_e = 1.60217662e-19; % [C] electron charge
    SI_me = 9.10938356e-31; % [kg] electron mass
    SI_eps0 = 8.85418782e-12; % [m^-3 kg^-1 s^4 A^2] permittivity of free space
    
    
    %% GET CONFIG
    addpath('..');
    
    
    %% GET INPUT PARAMETERS
    rp = rpinputParser(project, sim);
    
    % plasma parameters
    omega_p = sqrt(rp.n0 * SI_e^2/(SI_me*SI_eps0));
    
    %% DEFINE AXES
    
    % axes and labels
    [axes, axlbl] = qaxes( project, sim );
    
    % find ROIs (or set to default)
    if ~( exist('ROIs','var') && numel(ROIs)>=1 && ~isempty(ROIs{1}) )
        ROIs{1} = axes.X; % x-ROI
    end
    if ~( exist('ROIs','var') && numel(ROIs)>=2 && ~isempty(ROIs{2}) )
        ROIs{2} = axes.Y; % x-ROI
    end
    if ~( exist('ROIs','var') && numel(ROIs)>=3 && ~isempty(ROIs{3}) )
        ROIs{3} = axes.Z; % x-ROI
    end
    for i = 4:6
        if ~( exist('ROIs','var') && numel(ROIs)>=i && ~isempty(ROIs{i}) )
            ROIs{i} = []; % x', y', E ROIs
        end
    end
    ROI.X = [min(ROIs{1}), max(ROIs{1})];
    ROI.Y = [min(ROIs{2}), max(ROIs{2})];
    ROI.Z = [min(ROIs{3}), max(ROIs{3})];
    ROI.XP = [min(ROIs{4}), max(ROIs{4})];
    ROI.YP = [min(ROIs{5}), max(ROIs{5})];
    ROI.E = [min(ROIs{6}), max(ROIs{6})];
    
    % slice info
    dump.field.X = (rp.dump.slice.field.x + min(axes.X));
    dump.field.Y = (rp.dump.slice.field.y + min(axes.Y));
    dump.field.Z = (rp.dump.slice.field.z + min(axes.Z));
    dump.beam.X = (rp.dump.slice.beam.x + min(axes.X));
    dump.beam.Y = (rp.dump.slice.beam.y + min(axes.Y));
    dump.beam.Z = (rp.dump.slice.beam.z + min(axes.Z));
    dump.plasma.X = (rp.dump.slice.plasma.x + min(axes.X));
    dump.plasma.Y = (rp.dump.slice.plasma.y + min(axes.Y));
    dump.plasma.Z = (rp.dump.slice.plasma.z + min(axes.Z));

    %% PLOT HELPER FUNCTIONS
    
    % color label and scale for a given data type
    function [clbl, scale, cmin, cmax] = colorScale(type, cmin, cmax)
        % color label and scale
        switch type
            case 'E'
                clbl = 'E-field [GV/m]';
                scale = SI_me*SI_c*omega_p/(SI_e*1e9); % norm to GV/m
            case 'B'
                clbl = 'B-field [T]';
                scale = SI_me*omega_p/SI_e; % norm to T
            case 'F'
                clbl = 'Force [GeV/m]';
                scale = SI_me*SI_c*omega_p/(SI_e*1e9); % norm to GeV/m
            case 'G'
                clbl = 'Magnetic field gradient [T/m]';
                scale = SI_me*omega_p/(SI_e); % norm to GeV/m
            case 'Q'
                clbl = 'Slice density [cm^{-3}]';
                scale = rp.n0*1e-6;
        end
         % symmetrize field colormaps
        if ~strcmp(type,'Q')
            [cmin, cmax] = deal(min(cmin,-cmax), max(cmax,-cmin));
        end
    end
    
    % select dump type given a plot label
    function dmp = dumpType(type, comp)
        switch type
            case {'E','B','F','G'} % fields
                dmp = dump.field;
            case 'Q' % charge
                switch comp
                    case 'B' % beams
                        dmp = dump.beam;
                    case 'P' % plasma
                        dmp = dump.plasma;
                    otherwise % fail if plasma and beam slices are different
                        assert(dump.beam.X==dump.plasma.X);
                        assert(dump.beam.Y==dump.plasma.Y);
                        assert(dump.beam.Z==dump.plasma.Z);
                        dmp = dump.beam;
                end
        end
    end

    %% PLOTTING FUNCTIONS
    
    % plot fields 2D
    function ax = slicePlot(row, col, num, slicedData, lbl, s, xdir, ydir, cmin, cmax, lineout, type, comp)
        ax = subplot(row,col,num); 
        
        % scalings and color labels
        if ~isnan(cmin) && ~isnan(cmax)
            [clbl, scale, cmin, cmax] = colorScale(type, cmin, cmax);
        else
            [clbl, scale, cmin, cmax] = colorScale(type, min(min(slicedData)), max(max(slicedData)));
            cmin = cmin*scale; % scale limits
            cmax = cmax*scale;
        end
        
        % show image
        imagesc(axes.(xdir), axes.(ydir), scale*slicedData);
        set(gca,'ydir','normal');
        if strcmp(xdir,'Z')
            set(gca,'xdir','reverse');
        end
        
        % axes
        axis equal;
        xlabel(axlbl.(xdir)); ylabel(axlbl.(ydir));
        xlim(ROI.(xdir)); ylim(ROI.(ydir));
        
        % title
        proj = strrep(strrep('XYZ',xdir,''),ydir,'');
        dmp = dumpType(type, comp);
        form_lbl = [lbl ' \rm(' lower(proj) ' = ' num2str(dmp.(proj),3) ' um, s = ' num2str(s*1e3,3) ' mm)'];
        title(form_lbl);
        
        % colormap
        colormap(ax, cmap(cmin, cmax));
        caxis([cmin, cmax]);
        c = colorbar;
        c.Label.String = clbl;
        
        % horizontal lineout
        if exist('lineout','var') && ~isnan(lineout)
            warning ('off','all'); % ignore axis equal warning
            line(get(gca,'xlim'), lineout*[1,1],'LineStyle',':','Color',left_color,'LineWidth',1.2);
            yyaxis right;
            [~, index] = min(abs(axes.(ydir)-lineout));
            plot(axes.(xdir), scale*slicedData(index,:),'-','LineWidth',1.3);
        end
    end

    % 3D slice plots
    function ax = slice3Dplot(row, col, num, slicedData, lbl, s, cmin, cmax, type, comp)
        ax = subplot(row, col, num);
        
        % scalings and color labels
        if ~isnan(cmin) && ~isnan(cmax)
            [clbl, scale, cmin, cmax] = colorScale(type, cmin, cmax);
        else
            cmin = min([min(min(slicedData.XY)), min(min(slicedData.XZ)), min(min(slicedData.YZ))]);
            cmax = max([max(max(slicedData.XY)), max(max(slicedData.XZ)), max(max(slicedData.YZ))]);
            [clbl, scale, cmin, cmax] = colorScale(type, cmin, cmax);
            cmin = cmin*scale; % scale limits
            cmax = cmax*scale;
        end
        
        % find dump type
        dmp = dumpType(type, comp);
                
        % make slice/projection plots
        [XX, YY] = meshgrid(axes.X, axes.Y); % Z-slice (XY projection)
        ZZ = dmp.Z + zeros(size(slicedData.XY'));
        surf(ZZ, XX, YY, scale*slicedData.XY', 'EdgeColor','none');
        hold on;
        
        YY = dmp.Y + zeros(size(slicedData.XZ)); % Y-slice (XZ projection)
        surf(axes.Z, axes.X, YY, scale*slicedData.XZ, 'EdgeColor','none');
        
        [YY, ZZ] = meshgrid(axes.Y, axes.Z); % X-slice (YZ projection)
        XX = dmp.X + zeros(size(slicedData.YZ));
        surf(ZZ, XX, YY, scale*slicedData.YZ', 'EdgeColor','none');
        hold off;
        
        % axes, limits and labels
        set(gca,'ydir','reverse');
        axis equal; % fix axes
        xlim(ax, ROI.Z); ylim(ax, ROI.X); zlim(ax, ROI.Y); 
        xlabel('z [\mum]'); ylabel('x [\mum]'); zlabel('y [\mum]');
        
        % title
        form_lbl = [lbl ' \rm(s = ' num2str(s*1e3,3) ' mm)'];
        title(form_lbl);
        
        % apply colormap
        colormap(ax, cmap(cmin,cmax));
        caxis([cmin,cmax]);
        c = colorbar;
        c.Label.String = clbl;
        
        % fix up limits and
        rotate3d on;
    end

    % plot phase space
    function ax = phaseSpacePlot(row, col, num, xdir, ydir, s, cmn, cmx)
        ax = subplot(row, col, num);
        
        % determine nature of each axis (for the units shown)
        xspatial = (strcmp(xdir,'X') || strcmp(xdir,'Y') || strcmp(xdir,'Z'));
        xenergy = strcmp(xdir,'E');
        xangle = (strcmp(xdir,'XP') || strcmp(xdir,'YP'));
        yspatial = (strcmp(ydir,'X') || strcmp(ydir,'Y') || strcmp(ydir,'Z'));
        yenergy = strcmp(ydir,'E');
        yangle = (strcmp(ydir,'XP') || strcmp(ydir,'YP'));
        if  xspatial && yspatial
            unitstr = '[\mum^{-2}]';
        elseif xangle && yangle
            unitstr = '[\murad^{-2}]';
        elseif xspatial && yenergy || yspatial && xenergy
            unitstr = '[\mum^{-1} GeV^{-1}]';
        elseif xspatial && yangle || yspatial && xangle
            unitstr = '[\mum^{-1} \murad^{-1}]';
        elseif xenergy && yangle || yenergy && xangle
            unitstr = '[\murad^{-1} GeV^{-1}]';
        end
        
        % make 2D histogram
        xdata = qp.Q.total.phase.(xdir);
        ydata = qp.Q.total.phase.(ydir);
        xbins = ceil(sqrt(numel(xdata)));
        ybins = ceil(sqrt(numel(ydata)));
        [h, axs] = hist3([xdata, ydata], [xbins, ybins]);
        
        % make electrons negative density (only works if all beams are same charge)
        charge = rp.beam{1}.charge;
        h = h*charge;
        
        % fix charge normalization
        resx = (max(axs{1})-min(axs{1}))/numel(axs{1}); % [um/urad/GeV] x-length per histogram pixel
        resy = (max(axs{2})-min(axs{2}))/numel(axs{2}); % [um/urad/GeV] y-length per histogram pixel
        
        Nsamp = rp.dump.phase.beam.sample; % dump fraction
        [Ntotal, Nmacros] = deal(1);
        for j = 1:numel(rp.beam)
            Ntotal = Ntotal + rp.beam{j}.N; % total beam particles
            Nmacros = Nmacros + rp.beam{j}.res.x * rp.beam{j}.res.y * rp.beam{j}.res.z/Nsamp;
        end
        Npermacro = Ntotal/Nmacros; % number of particles represented by each dump particle
        nb = Npermacro/(resx*resy); % [axis-dependent unit] charge surface density per pixel
        h = h*nb; % scale the image
        
        % making our own 2D histogram
        imagesc(axs{1},axs{2}, h');
        set(gca,'ydir','normal');
        if strcmp(xdir,'Z')
            set(gca,'xdir','reverse');
        end
        
        % axis equal if both spatial axes
        if xspatial && yspatial
            axis equal;
        end
        
        % fix ROI
        if numel(ROI.(xdir))>0
            xlim([min(ROI.(xdir)), max(ROI.(xdir))]);
        end
        if numel(ROI.(ydir))>0
            ylim([min(ROI.(ydir)), max(ROI.(ydir))]);
        end
        
        % labels
        xlabel(axlbl.(xdir));
        ylabel(axlbl.(ydir));
        title(['Phase space \rm(s = ' num2str(s*1e3,3) ' mm)'])
        
        % colormap
        if ~isnan(cmn) && ~isnan(cmx)
            cmin = cmn;
            cmax = cmx;
        else
            cmin = min(min(h));
            cmax = max(max(h));
        end
        colormap(ax, cmap(cmin, cmax));
        caxis([cmin,cmax]);
        c = colorbar;
        c.Label.String = ['Projected density ' unitstr];
        
    end

    %% PLOTTING
    
    % default docked
    set(0,'DefaultFigureWindowStyle','docked');
    
    % new figure
    set(gcf,'color','w');
    clf;
    
    % set axis colors
    left_color = [0.5,0.5,0.5];
    left_color = [229 78 14]/256;
    right_color = [0,0,0];
    set(gcf,'defaultAxesColorOrder',[left_color; right_color]);
    
    % determine subplot layout
    Nplots = numel(data);
    figlims = get(gcf, 'Position');
    if figlims(3) < figlims(4) % distribute vertically if aspect ratio < 1
        Nrows = ceil(sqrt(Nplots));
        Ncols = ceil(Nplots/Nrows);
    else
        Ncols = ceil(sqrt(Nplots));
        Nrows = ceil(Nplots/Ncols);
    end
    
    % parse steps (animation possible with cell input)
    cycle = false;
    if isa(steps, 'cell')
        if steps{2}>0
            animationPause = @() pause(steps{2});
            if numel(steps)>2
                cycle = steps{3}; % repeat the animation?
            end
        else
            animationPause = @pause;
        end
        
        steps = steps{1}; %
    elseif isa(steps, 'double')
        animationPause = @drawnow;
    end
    % cycle steps
    onceFlag = true;
    link3Daxisflag = false;
    ax3D = zeros(Nplots,1);
    while onceFlag || cycle
        for step = steps
            % find distance travelled
            s = rp.sim.time.step * step / omega_p * SI_c; % [m]
            
            % get data
            qp = qextract(project, sim, step);
            try 
                qp = qextract(project, sim, step);
            catch
                disp(['Cannot get data for step ' num2str(step) '.']);
                continue; % (skip if non-existant)
            end; 

            % cycle through plots
            for i = 1:numel(data)
                prj_raw = data{i}{2};
                prj = sort(prj_raw);
                type = data{i}{1}(1); % field/force/density: 'E', 'B', 'F' or 'Q'
                comp = data{i}{1}(2:end); % field/force/density component: 'X', 'Y', 'Z', 'R', 'TH', 'P', 'B'
                if strcmp(type,'Q') && isempty(comp) % combine plasma and beams if not specified
                    comp = 'total'; 
                end
                lbl = [type '_{' strrep(lower(comp),'th','\theta') '}']; % format the plot title
                % color limits
                if numel(data{i})>=3
                    if numel(data{i}{3})>1
                        cmn = min(data{i}{3});
                        cmx = max(data{i}{3});
                    else
                        cmx = max(data{i}{3},0);
                        cmn = min(data{i}{3},0);
                    end
                else
                    [cmn, cmx] = deal(NaN);
                end
                
                if numel(prj)==2 % 2D slice
                    if numel(data{i})>=4 % lineouts
                        lineout = data{i}{4};
                    else
                        lineout = NaN;
                    end
                    [xdir, ydir] = deal(prj_raw(1), prj_raw(2)); % directions to plot: x/y/z
                    assert(strcmp(xdir,'X') || strcmp(xdir,'Y') || strcmp(xdir,'Z')); % test x-dir
                    assert(strcmp(ydir,'X') || strcmp(ydir,'Y') || strcmp(ydir,'Z')); % test y-dir
                    f = qp.(type).(comp).(prj);
                    if strcmp(prj,prj_raw); f = f'; end % swap axes if required
                    slicePlot(Nrows, Ncols, i, f, lbl, s, xdir, ydir, cmn, cmx, lineout, type, comp);
                elseif strcmp(sort(prj),'XYZ') 
                    % 3D slice plot
                    ax3D(i) = slice3Dplot(Nrows, Ncols, i, qp.(type).(comp), lbl, s, cmn, cmx, type, comp);
                    link3Daxisflag = true;
                elseif numel(strfind(prj, '_')) % phase space
                    xydirs = strsplit(prj_raw,'_');
                    xdir = xydirs{1};
                    ydir = xydirs{2};
                    phaseSpacePlot(Nrows, Ncols, i, xdir, ydir, s, cmn, cmx);
                end
            end
            
            animationPause();
        end
        onceFlag = false;
    end
    
    % link 3D axes for rotation
    if link3Daxisflag
        hlink = linkprop(ax3D, 'View');
        setappdata(gcf, 'link3D', hlink);
    end

end
