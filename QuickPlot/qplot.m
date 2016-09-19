function [] = qplot(project, sim, data, steps, ROIs)

    % (proj, sim, {{field, projection}, {...}}, {xlims, ylims, zlims})
    
    % field: 'EX', 'EY', 'EZ', 'ER', 'ETH' for E-FIELDS in (x/y/z/r/theta)
    %       'BX', 'BY', 'BZ', 'BR', 'BTH' for B-FIELDS in (x/y/z/r/theta)
    %       'FX', 'FY', 'FZ', 'FR', 'FTH' for FORCES in (x/y/z/r/theta)
    %       'QB' for BEAM DENSITY (total)
    %       'QP' for PLASMA DENSITY (total)
    
    % proj: 'XY', 'XZ', 'YZ' for slice plots 
    %       'XYZ' for 3D combo slice plot
    
    % steps: first:last
    %     or {first:last, delay}
    %     or {first:last, delay, loop?}
    
    % ROIs: {xmin:xmax, ymin:ymax, zmin:zmax}
    %       [] or blank gives default: full ROI
    
    
    %% CONSTANTS
    SI_c = 299792458; % [m/s] speed of light
    SI_e = 1.60217662e-19; % [C] electron charge
    SI_me = 9.10938356e-31; % [kg] electron mass
    SI_eps0 = 8.85418782e-12; % [m^-3 kg^-1 s^4 A^2] permittivity of free space
    
    
    %% GET CONFIG
    addpath('..');
    outputfolder = CONFIG('outputs');
    
    
    %% GET INPUT PARAMETERS
    rp = rpinputParser(project, sim);
    Nspecies = numel(rp.plasma);
    
    % plasma parameters
    omega_p = sqrt(rp.n0 * SI_e^2/(SI_me*SI_eps0));
    k_p = SI_c/omega_p;
    %lambda_p = 2*pi/k_p;
    
    
    %% DEFINE AXES
    
    % labels
    [axlbl.X, axlbl.Y, axlbl.Z] = deal('x [\mum]', 'y [\mum]', 'z [\mum]');
    
    % axes and offsets
    axes.X = linspace(-rp.sim.dim.x/2, rp.sim.dim.x/2, 2^rp.sim.ind.x);
    axes.Y = linspace(-rp.sim.dim.y/2, rp.sim.dim.y/2, 2^rp.sim.ind.y);
    axes.Z = linspace(-rp.beam{1}.offset.z, rp.sim.dim.z-rp.beam{1}.offset.z, 2^rp.sim.ind.z);
    
    % find ROIs (or set to default)
    if ~( exist('ROIs','var') && numel(ROIs)>=1 && ~isempty(ROIs{1}) )
        ROIs{1} = axes.X;
    end
    if ~( exist('ROIs','var') && numel(ROIs)>=2 && ~isempty(ROIs{2}) )
        ROIs{2} = axes.Y;
    end
    if ~( exist('ROIs','var') && numel(ROIs)>=3 && ~isempty(ROIs{3}) )
        ROIs{3} = axes.Z;
    end
    ROI.X = [min(ROIs{1}), max(ROIs{1})];
    ROI.Y = [min(ROIs{2}), max(ROIs{2})];
    ROI.Z = [min(ROIs{3}), max(ROIs{3})];
    
    % slice info
    dump.X = (rp.dump.slice.x + min(axes.X));
    dump.Y = (rp.dump.slice.y + min(axes.Y));
    dump.Z = (rp.dump.slice.z + min(axes.Z));
    
    
    %% EXTRACT FIELDS AND DENSITIES
    
    % function to extract all fields and densities 
    % (can be optimized by only extracting fields needed)    
    function qp = extract(step)
        % TODO: add testing for whether projections exist or not
        
        % extraction helper function
        simpath = [outputfolder '/' project '/' sim];
        geth5 = @(var) h5read([simpath '/' var '/' var '_' num2str(step,'%04d') '.h5']);
    
        for proj = {'XY','XZ','YZ'} % {'XZ'} %
            % fields
            for field = {'E','B'}
                for dim = {'X','Y','Z'}
                    qp.(field{:}).(dim{:}).(proj{:}) = geth5(['F' field{:} dim{:} '-' proj{:}]);
                end
            end

            % plasmas
            for species = 1:Nspecies
                Q = geth5(['QEP' num2str(species) '-' proj{:}]);
                if species == 1 % first species
                    qp.Q.P.(proj{:}) = Q;
                else % otherwise add
                    qp.Q.P.(proj{:}) = qp.Q.P.(proj{:}) + Q;
                end
            end

            % beams
            qp.Q.B.(proj{:}) = geth5(['QEB-' proj{:}]);

            % both (beams + plasma)
            qp.Q.total.(proj{:}) = qp.Q.B.(proj{:}) + qp.Q.P.(proj{:});
        end
        
        % add forces and fields in polar coordinates
        [qp.E, qp.B, qp.F] = polarFields(axes.X, axes.X, qp.E, qp.B);
    end
    

    %% PLOT HELPER FUNCTIONS
    
    % plot fields
    function ax = slicePlot(row, col, num, slicedData, lbl, xdir, ydir, cmin, cmax, lineout)
        ax = subplot(row,col,num); 
        
        % get color scale
        
        
        % scalings and color labels
        if ~isnan(cmin) && ~isnan(cmax)
            [clbl, scale, cmin, cmax] = colorScale(lbl(1), cmin, cmax);
        else
            [clbl, scale, cmin, cmax] = colorScale(lbl(1), min(min(slicedData)), max(max(slicedData)));
            cmin = cmin*scale; % scale limits
            cmax = cmax*scale;
        end
        
        % plot the data
        imagesc(axes.(xdir), axes.(ydir), scale*slicedData);
        set(gca,'ydir','normal');
        if strcmp(xdir,'Z')
            set(gca,'xdir','reverse');
        end
        
        % axes
        axis equal;
        xlabel(axlbl.(xdir)); ylabel(axlbl.(ydir));
        xlim(ROI.(xdir)); ylim(ROI.(ydir));
        title(lbl);
        
        % colormap
        colormap(ax, cmap(cmin, cmax));
        caxis([cmin, cmax]);
        c = colorbar;
        c.Label.String = clbl;
        
        % horizontal lineout
        if exist('lineout','var') && ~isnan(lineout)
            line(get(gca,'xlim'), lineout*[1,1],'LineStyle',':','Color',left_color,'LineWidth',1.2);
            yyaxis right;
            [~, index] = min(abs(axes.(ydir)-lineout));
            plot(axes.(xdir), scale*slicedData(index,:),'-','LineWidth',1.3);
        end
    end

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
            case 'Q'
                clbl = 'Density [cm^{-3}]';
                scale = rp.n0*1e-6;
        end
         % symmetrize field colormaps
        if ~strcmp(type,'Q')
            [cmin, cmax] = deal(min(cmin,-cmax), max(cmax,-cmin));
        end
    end
    
    % 3D slice plots
    function ax = slice3Dplot(row, col, num, slicedData, lbl, cmin, cmax)
        ax = subplot(row, col, num);
        
        % scalings and color labels
        if ~isnan(cmin) && ~isnan(cmax)
            [clbl, scale, cmin, cmax] = colorScale(lbl(1), cmin, cmax);
        else
            cmin = min([min(min(slicedData.XY)), min(min(slicedData.XZ)), min(min(slicedData.YZ))]);
            cmax = max([max(max(slicedData.XY)), max(max(slicedData.XZ)), max(max(slicedData.YZ))]);
            [clbl, scale, cmin, cmax] = colorScale(lbl(1), cmin, cmax);
            cmin = cmin*scale; % scale limits
            cmax = cmax*scale;
        end
        
        % make slice/projection plots
        [XX, YY] = meshgrid(axes.X, axes.Y); % Z-slice (XY projection)
        ZZ = dump.Z + zeros(size(slicedData.XY'));
        surf(ZZ, XX, YY, scale*slicedData.XY', 'EdgeColor','none');
        hold on;
        
        YY = dump.Y + zeros(size(slicedData.XZ)); % Y-slice (XZ projection)
        surf(axes.Z, axes.X, YY, scale*slicedData.XZ, 'EdgeColor','none');
        
        [YY, ZZ] = meshgrid(axes.Y, axes.Z); % X-slice (YZ projection)
        XX = dump.X + zeros(size(slicedData.YZ));
        surf(ZZ, XX, YY, scale*slicedData.YZ', 'EdgeColor','none');
        hold off;
        
        % axes, limits and labels
        set(gca,'ydir','reverse');
        axis equal; % fix axes
        xlim(ax, ROI.Z); ylim(ax, ROI.X); zlim(ax, ROI.Y); 
        xlabel('z [\mum]'); ylabel('x [\mum]'); zlabel('y [\mum]');
        title(lbl); % add title
        
        % apply colormap
        colormap(ax, cmap(cmin,cmax));
        caxis([cmin,cmax]);
        c = colorbar;
        c.Label.String = clbl;
        
        % fix up limits and
        rotate3d on;
    end


    %% PLOTTING
    
    % default docked
    set(0,'DefaultFigureWindowStyle','docked');
    
    % new figure
    set(gcf,'color','w');
    left_color = [0.5,0.5,0.5];
    left_color = [229 78 14]/256;
    
    right_color = [0,0,0];
    set(gcf,'defaultAxesColorOrder',[left_color; right_color]);
    
    % determine subplot layout
    Nplots = numel(data);
    %Nrows = ceil(sqrt(Nplots));
    %Ncols = ceil(Nplots/Nrows);
    Ncols = ceil(sqrt(Nplots));
    Nrows = ceil(Nplots/Ncols);
    
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
            % find time and distance
            t = rp.sim.time.step * step / omega_p *1e12; % [ps] time of s
            
            % get data 
            try 
                qp = extract(step);
            catch
                disp(['Cannot get data for step ' num2str(step) '.']);
                continue; 
            end; % (skip if non-existant)

            % cycle through plots
            for i = 1:numel(data)
                prj = sort(data{i}{2});
                type = data{i}{1}(1); % field/force/density: 'E', 'B', 'F' or 'Q'
                comp = data{i}{1}(2:end); % field/force/density component: 'X', 'Y', 'Z', 'R', 'TH', 'P', 'B'
                if strcmp(type,'Q') && isempty(comp); comp = 'total'; end % combine plasma and beams if not specified
                lbl1 = [type '_{' strrep(lower(comp),'th','\theta') '}']; % format the plot title
                lbl = [lbl1 ' (t = ' num2str(t,'%3.1f') ' ps)']; % add time
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
                    % lineouts
                    if numel(data{i})>=4
                        lineout = data{i}{4};
                    else
                        lineout = NaN;
                    end
                    [xdir, ydir] = deal(data{i}{2}(1), data{i}{2}(2)); % directions to plot: x/y/z
                    f = qp.(type).(comp).(prj);
                    if strcmp(prj,data{i}{2}); f = f'; end % swap axes if required
                    slicePlot(Nrows, Ncols, i, f, lbl, xdir, ydir, cmn, cmx, lineout);
                elseif strcmp(prj,'XYZ') % 3D slice plot
                    ax3D(i) = slice3Dplot(Nrows, Ncols, i, qp.(type).(comp), lbl, cmn, cmx);
                    link3Daxisflag = true;
                end
            end

            animationPause();
        end
        onceFlag = false;
    end
    
    if link3Daxisflag
        hlink = linkprop(ax3D, 'View');
        setappdata(gcf, 'link3D', hlink);
    end
    
    % TODO: allow different slice planes in rp struct (as in rpinput)
    
end
