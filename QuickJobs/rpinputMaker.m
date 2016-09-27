function rp = rpinputMaker( beams, n0, plasmas, z_end, Nsteps, dims, dumps, ress )

    % beams = {{E, N, q, [emx emy sigE/E], {'tri',[sigx sigy sigz]}, [dx dy dz], evolve, quiet, beamres}, {...}}
    %         { ...                    ... {'bi',[sigx sigy],[z, f]}, ...}, {...}}
    %         { ...                    ... {'twiss',[betx bety],[alfx alfy],sigz}, ...}, {...}}
    %         { ...                    ... {'arb',file}, ...}, {...}}
    % NOTE: E has unit [GeV]
    % NOTE: (sigx/y/z, dx/y/z) have units [um]
    % NOTE: (z, f) for 'bi' are column vectors, where z has unit [um]
    % NOTE: (betx/y, alfx/y) have units [m]
    % NOTE: sigE/E = 0 by default
    % NOTE: dx/y = 0 by default
    % NOTE: dz = 3*sigz by default
    % NOTE: q can be given as 'e+'/'e-' or +1/-1;
    % NOTE: evolve defaults to true
    
    % NOTE: n0 has units [cm^-3]
    
    % plasmas = {{'uni', [s,n], res}, {...}}
    %           {{'hc', [r1,r2,n], [s,n], res}, {...}}
    %           {{'cyl', [r1,n], [s,n], res}, {...}}
    %           {{'arb', [r,n], [s,n], res}, {...}}
    % NOTE: n is relative to n0
    % NOTE: n = 1 by default
    % NOTE: longitudinal profile [s,n] is disabled if not listed
    % NOTE: (r, n) for 'arb' are column vectors, where r has unit [um]
    % NOTE: (s, n) are column vectors, where s has unit [um] and n is relative
    % NOTE: (r1, r2) hava units [um]
    
    % dims = [dimx, dimy, dimz]
    
    % ress = [resx, resy, resz]
    % NOTE: resolutions ress = [512, 512, 512] by default
    
    % dumps = {slice, phasespace}
    %       slice = {freq, [x y z]}
    %               {{fieldFreq, [x y z]}, {beamFreq, [x y z]}, {plasmaFreq, [x y z]}}
    %       phasespace = {freq, sample}
    %                    {{beamFreq, beamSample}, {plasmaFreq, plasmaSample}}            
    % NOTE: [freq, fieldFreq, beamFreq, plasmaFreq] = 1 [units of timestep] by default
    % NOTE: [x,y,z] = (center, center, first beam's z-centroid) by default
    % NOTE: [freq, beamFreq, plasmaFreq] = 0 by default
    % NOTE: [sample, beamSample, plasmaSample] = 128 by default (sparse sampling rate)
    % EXAMPLES:
    % rp = rpinputMaker({{20.35,2.5e9,'e+',[10 10],{'tri',[20,20,35]}}},3e16,{{'uni'},{'hc',[200,260,0.5]}},0.0126,10,[700 700 500],{10});
    
    
    %% GUIDE if no arguments
    if nargin == 0
        disp('ARGUMENT GUIDE to the RPINPUT MAKER:');
        disp('rpinputMaker(beams, n0, plasmas, z_end, Nsteps, [dimx, dimy, dimz], dumps, [resx, resy, resz])');
        disp('beams = {{E, N, q, [emNx emNy sigE/E], {''tri'', [sigx sigy sigz]}, [dx dy dz], evolve?, quiet?, beamres}, {...}}');
        disp('                                       {''bi'', [sigx sigy],[z, f]}');
        disp('                                       {''twiss'', [betx bety], [alfx alfy], sigz}');
        disp('                                       {''arb'', file}');
        disp('plasmas = {{''uni'', [s n], res}, {...}}');
        disp('           {''hc'', [r1 r2 n], [s n], res}');
        disp('           {''cyl'', [r1 n], [s n], res}');
        disp('           {''arb'', [r n], [s n], res}');
        disp('dumps = {slice, phasespace}');
        disp('         slice = {freq [x y z]}');
        disp('                 {{fieldFreq [x y z]}, {beamFreq [x y z]}, {plasmaFreq [x y z]}}');
        disp('    phasespace = {{beamFreq beamSample}, {plasmaFreq plasmaSample})');
        return;
    end
    
    
    %% CONFIG
    templatesfolder = [CONFIG('qjobsfolder') '/templates'];
    tempfilesfolder = [CONFIG('qjobsfolder') '/tempfiles'];
    
    % assertions
    assert(numel(beams)>0,'Must have at least one beam.');
    assert(numel(plasmas)>0,'Must have at least one plasma.');
    
    % constants
    SI_c = 299792458; % [m/s] speed of light
    SI_e = 1.60217662e-19; % [C] electron charge
    SI_me = 9.10938356e-31; % [kg] electron mass
    SI_eps0 = 8.85418782e-12; % [m^-3 kg^-1 s^4 A^2] permittivity of free space
    
    
    %% PARSE BEAMS
    fprintf('\nBEAMS\n');
    fprintf('%d beam(s)\n', numel(beams));
    for i = 1:numel(beams)
        beam = beams{i};
        % gamma: convert from energy [GeV]
        rp.beam{i}.gamma = beam{1} * 1e9*SI_e/(SI_me*SI_c^2); 
        % number of particles
        rp.beam{i}.N = beam{2};
        % charge sign
        if ischar(beam{3}) && strcmpi(beam{3},'e+')
            rp.beam{i}.charge = 1;
        elseif ischar(beam{3}) && strcmpi(beam{3},'e-')
            rp.beam{i}.charge = -1;
        else
            rp.beam{i}.charge = beam{3};
        end
        % emittance
        rp.beam{i}.emit = struct('x', beam{4}(1), 'y', beam{4}(2));
        % energy spread
        if numel(beam{4})>2
            rp.beam{i}.espread = beam{4}(3);
        else
            rp.beam{i}.espread = 0; % default zero
        end
        % beam size (x,y,z)
        rp.beam{i}.file = '';
        rp.beam{i}.params = {0, 0, 0, 0, 0};
        switch beam{5}{1}
            case 'tri'
                rp.beam{i}.routine = 1; % routine
                rp.beam{i}.params{2} = beam{5}{2}; % sigmas
                rp.beam{i}.params{3} = [beam{4}, rp.beam{i}.espread]; % emittance and E-spread
            case 'bi'
                rp.beam{i}.params{2} = [beam{5}{2}, size(beam{5}{3},1)]; % sigmas and size of z-profile
                rp.beam{i}.params{3} = [beam{4}, rp.beam{i}.espread, zeros(1,6)]; % emittance, E-spread and centroids
                rp.beam{i}.params{4} = beam{5}{3}(:,2); % z-profile f(z)
                rp.beam{i}.params{5} = beam{5}{3}(:,1); % z-profile zs
                if numel(beam{5}) > 3 && beam{5}{4} % routine
                    rp.beam{i}.routine = 3;
                    rp.beam{i}.rand = 1;
                else
                    rp.beam{i}.routine = 2; % default
                    rp.beam{i}.rand = false;
                end
            case 'arb'
                rp.beam{i}.routine = 4; % routine
                rp.beam{i}.file = beam{5}{2}; % file path
            case 'twiss'
                rp.beam{i}.routine = 5; % routine
                rp.beam{i}.params{2} = [beam{5}{2}, beam{5}{3}, beam{5}{4}]; % twiss params and sigma_z
                rp.beam{i}.params{3} = [beam{4}, rp.beam{i}.espread]; % emittance and E-spread
        end
        % beam offsets (x,y,z) : (0,0,0)
        if numel(beam) > 5
            rp.beam{i}.offset = struct('x', beam{6}(1), 'y', beam{6}(2), 'z', beam{6}(3));
        else
            rp.beam{i}.offset = struct('x', 0, 'y', 0);
            sigmas_offset = 3; % number of sigmas to offset in z
            switch beam{5}{1}
                case 'tri'
                    rp.beam{i}.offset.z = sigmas_offset*beam{5}{2}(3);
                case 'twiss'
                    rp.beam{i}.offset.z = sigmas_offset*beam{5}{4};
                otherwise
                    rp.beam{i}.offset.z = 0;
            end
        end
        % evolve? : true
        if numel(beam) > 6
            rp.beam{i}.evolve = beam{7};
        else
            rp.beam{i}.evolve = true; % default: true
        end
        % quiet? : true
        if numel(beam) > 7
            rp.beam{i}.quiet = beam{8};
        else
            rp.beam{i}.quiet = true; % default: true
        end
        % beam resolution (macroparticles in x,y,z) : (256,256,256)
        if numel(beam) > 8
            rp.beam{i}.res = struct('x', beam{9}(1), 'y', beam{9}(2), 'z', beam{9}(3));
        else
            rp.beam{i}.res = struct('x', 256, 'y', 256, 'z', 256); % default: true
        end
        
        fprintf('> N = %G, gamma = %g, charge = %+d, \n', rp.beam{i}.N, rp.beam{i}.gamma, rp.beam{i}.charge);
        fprintf('  emN = (%g, %g) um, E-spread = %g%%, evolve = %s, quiet = %s, \n', rp.beam{i}.emit.x, rp.beam{i}.emit.y, rp.beam{i}.espread*100, bool2str(rp.beam{i}.evolve), bool2str(rp.beam{i}.quiet));
        switch beam{5}{1}
            case 'tri'
                fprintf('  tri-Gaussian, sigmas = (%g, %g, %g) um, \n', rp.beam{i}.params{2}(1), rp.beam{i}.params{2}(2), rp.beam{i}.params{2}(3));
            case 'bi'
                fprintf('  bi-Gaussian (random = %s), sigmas = (%g, %g) um, %d-step z-profile, \n', bool2str(rp.beam{i}.rand), rp.beam{i}.params{2}(1), rp.beam{i}.params{2}(2), rp.beam{i}.params{2}(3));
            case 'arb'
                fprintf('  arbitrary distribution, given by file: "%s", \n', rp.beam{i}.file);
            case 'twiss'
                fprintf('  Twiss params, beta = (%g, %g) m, beta = (%g, %g) m, sigma_z = %g um, \n', rp.beam{i}.params{2}(1), rp.beam{i}.params{2}(2), rp.beam{i}.params{2}(3),rp.beam{i}.params{2}(4), rp.beam{i}.params{2}(5));
        end
        fprintf('  offset = (%g, %g, %g) um, resolution = (%d, %d, %d)\n\n', rp.beam{i}.offset.x, rp.beam{i}.offset.y, rp.beam{i}.offset.z, rp.beam{i}.res.x, rp.beam{i}.res.y, rp.beam{i}.res.z);
        
    end
    
    %% PARSE PLASMAS
    disp('PLASMA');
    fprintf('Density: %2.2e /cc\n', n0);
    fprintf('%d species\n', numel(plasmas));
    rp.n0 = n0; % plasma density [cm^-3]
    for i = 1:numel(plasmas)
        plasma = plasmas{i};
        
        rp.plasma{i}.r = struct('args', [0,0,0], 'ns', 0, 'rs', 0); % defaults r-args and profile
        rp.plasma{i}.z = struct('enable', false, 'ns', 0, 'ss', 0); % long. profile?
        rp.plasma{i}.res = 2048; % default resolution (per dimension)
        
        switch plasma{1}
            case 'uni'
                rp.plasma{i}.type = 0; % type
                if numel(plasma) > 1
                    rp.plasma{i}.z.enable = true;
                    rp.plasma{i}.z.ss = plasma{2}(:,1)'; % long. profile ss
                    rp.plasma{i}.z.ns = plasma{2}(:,2)'; % long. profile ns
                end
                if numel(plasma) > 2
                    rp.plasma{i}.res = plasma{3}; % resolution
                end
            case 'hc'
                rp.plasma{i}.type = 70; % type
                if numel(plasma{2})>2 % r-args
                    rp.plasma{i}.r.args = [plasma{2}(1), plasma{2}(3), plasma{2}(2)];
                else
                    rp.plasma{i}.r.args = [plasma{2}(1), 1, plasma{2}(2)]; % default relative density : 1
                end
                if numel(plasma) > 2
                    rp.plasma{i}.z.enable = true;
                    rp.plasma{i}.z.ss = plasma{3}(:,1)'; % long. profile ss
                    rp.plasma{i}.z.ns = plasma{3}(:,2)'; % long. profile ns
                end
                if numel(plasma) > 3
                    rp.plasma{i}.res = plasma{4}; % resolution
                end
            case 'cyl'
                rp.plasma{i}.type = 19; % type
                if numel(plasma{2})>1 % r-args
                    rp.plasma{i}.r.args = [plasma{2}(1:2),0];
                else
                    rp.plasma{i}.r.args = [plasma{2}(1),1,0]; % default relative density : 1
                end
                if numel(plasma) > 2
                    rp.plasma{i}.z.enable = true;
                    rp.plasma{i}.z.ss = plasma{3}(:,1)'; % long. profile ss
                    rp.plasma{i}.z.ns = plasma{3}(:,2)'; % long. profile ns
                end
                if numel(plasma) > 3
                    rp.plasma{i}.res = plasma{4}; % resolution
                end
            case 'arb'
                rp.plasma{i}.type = 21; % type
                rp.plasma{i}.r.rs = plasma{2}(:,1)'; % radial profile ss
                rp.plasma{i}.r.ns = plasma{2}(:,2)'; % radial profile ns
                if numel(plasma) > 2
                    rp.plasma{i}.z.enable = true;
                    rp.plasma{i}.z.ss = plasma{3}(:,1)'; % long. profile ss
                    rp.plasma{i}.z.ns = plasma{3}(:,2)'; % long. profile ns
                end
                if numel(plasma) > 3
                    rp.plasma{i}.res = plasma{4}; % resolution
                end
        end
        % print info
        switch plasma{1}
            case 'uni'
                fprintf('> Uniform plasma, ');
            case 'hc'
                fprintf('> Hollow channel, r1 = %g um, r2 = %g um, rel. density = %g, ', rp.plasma{i}.r.args(1), rp.plasma{i}.r.args(3), rp.plasma{i}.r.args(2));
            case 'cyl'
                fprintf('> Cylindrical, r = %g um, rel. density = %g, ', rp.plasma{i}.r.args(1), rp.plasma{i}.r.args(2));
            case 'arb'
                fprintf('> Piece-wise profile, %d-part radial profile, ', numel(rp.plasma{i}.r.rs));
        end
        fprintf('1D-resolution = %d, z-profile = %s', rp.plasma{i}.res, bool2str(rp.plasma{i}.z.enable));
        if rp.plasma{i}.z.enable
            fprintf(' (%d-part profile)', numel(rp.plasma{i}.z.ss));
        end
        disp(' ');
    end
    
    %% PARSE SIMULATION
    fprintf('\nSIMULATION\n');
    omega_p = sqrt(n0*1e6 * SI_e^2/(SI_me*SI_eps0)); % [Hz] plasma frequency
    k_p = omega_p/SI_c; % [m^-1] plasma wavenumber
    
    % time and step
    rp.sim.time.step = round(z_end/Nsteps/SI_c*omega_p); % [1/omega_p]
    rp.sim.time.end = rp.sim.time.step*Nsteps; % [1/omega_p]
    % distance
    rp.sim.dist.end = rp.sim.time.end*SI_c/omega_p; % [m]
    rp.sim.dist.step = rp.sim.time.step*SI_c/omega_p; % [m]
    % dimensions
    rp.sim.dim = struct('x', dims(1), 'y', dims(2), 'z', dims(3));
    % resolution (indices)
    if nargin > 7
        rp.sim.ind = struct('x', round(log2(ress(1))), 'y', round(log2(ress(2))), 'z', round(log2(ress(3))));
    else
        rp.sim.ind = struct('x', 9, 'y', 9, 'z', 9); % default: log2(512)
    end
    % display for information
    fprintf('Dimensions = (%1.f, %1.f, %1.f) um,\n', rp.sim.dim.x, rp.sim.dim.y, rp.sim.dim.z);
    fprintf('Resolution = (2^%d, 2^%d, 2^%d) = (%d, %d, %d),\n', rp.sim.ind.x, rp.sim.ind.y, rp.sim.ind.z, 2^rp.sim.ind.x, 2^rp.sim.ind.y, 2^rp.sim.ind.z);
    fprintf('Step unit (1/k_p) = %2.2f um\n', 1/k_p*1e6); % display step unit
    fprintf('Doing %d steps: \n  0 : %g : %g ps \n  0 : %g : %g mm \n  0 : %g : %g step units\n\n', Nsteps, rp.sim.time.step/omega_p*1e12, rp.sim.time.end/omega_p*1e12, rp.sim.dist.step*1e3, rp.sim.dist.end*1e3, rp.sim.time.step, rp.sim.time.end);
    
    
    %% PARSE DUMPS
    % slice dumps
    slices = struct('freq',1,'x',0,'y',0,'z',rp.beam{1}.offset.z); % default slice
    if exist('dumps','var') && numel(dumps)>=1 && ~isempty(dumps{1})
        if iscell(dumps{1}{1}) % individual slice definitions
            indSlices = {slices, slices, slices}; % [field, beam, plasma]
            for i = 1:numel(indSlices)
                if numel(dumps{1})>=i % does it exist?
                    if numel(dumps{1}{i})>=1 % dump frequency
                        indSlices{i}.freq = dumps{1}{i}{1};
                    end
                    if numel(dumps{1}{i})>=2 % dump slice planes
                        indSlices{i}.x = dumps{1}{i}{2}(1);
                        indSlices{i}.y = dumps{1}{i}{2}(2);
                        indSlices{i}.z = dumps{1}{i}{2}(3);
                    end
                end
            end
            [rp.dump.slice.field, rp.dump.slice.beam, rp.dump.slice.plasma] = deal(indSlices{1}, indSlices{2}, indSlices{3});
        else % common slice definition
            if numel(dumps{1})>=1 % common dump frequency
                slices.freq = dumps{1}{1};
            end
            if numel(dumps{1})>=2 % common dump slice planes
                slices.x = dumps{1}{2}(1);
                slices.y = dumps{1}{2}(2);
                slices.z = dumps{1}{2}(3);
            end
            [rp.dump.slice.field, rp.dump.slice.beam, rp.dump.slice.plasma] = deal(slices);
        end
    else
        [rp.dump.slice.field, rp.dump.slice.beam, rp.dump.slice.plasma] = deal(slices);
    end
    % phase space dumps
    phases = struct('freq',0,'sample',128); % default
    if exist('dumps','var') && numel(dumps)>=2 && ~isempty(dumps{2})
        if iscell(dumps{2}{1}) % individual slice definitions
            indPhases = {phases, phases}; % [field, beam, plasma]
            for i = 1:numel(indPhases)
                if numel(dumps{2})>=i % does it exist?
                    if numel(dumps{2}{i})>=1 % dump frequency
                        indPhases{i}.freq = dumps{2}{i}{1};
                    end
                    if numel(dumps{2}{i})>=2 % dump slice planes
                        indPhases{i}.sample = dumps{2}{i}{2};
                    end
                end
            end
            [rp.dump.phase.beam, rp.dump.phase.plasma] = deal(indPhases{1}, indPhases{2});
        else % common phase space definition
            if numel(dumps{2})>=1 % common dump frequency
                phases.freq = dumps{2}{1};
            end
            if numel(dumps{2})>=2 % common dump sample
                phases.sample = dumps{2}{2};
            end
            [rp.dump.phase.beam, rp.dump.phase.plasma] = deal(phases);
        end
    else
        [rp.dump.phase.beam, rp.dump.phase.plasma] = deal(phases);
    end
	
    % print info
    function n = numDumps(freq) 
        if freq == 0; n = 0; else n = floor(Nsteps/freq); end
    end
    disp('DUMPS');
    fprintf('> Field slices: every %d step(s) = %d dump(s), ', rp.dump.slice.field.freq, numDumps(rp.dump.slice.field.freq));
    fprintf('at (%1.f, %1.f, %1.f) um\n', rp.dump.slice.field.x, rp.dump.slice.field.y, rp.dump.slice.field.z);
    fprintf('> Beam slices: every %d step(s) = %d dump(s), ', rp.dump.slice.beam.freq, numDumps(rp.dump.slice.beam.freq));
    fprintf('at (%1.f, %1.f, %1.f) um\n', rp.dump.slice.beam.x, rp.dump.slice.beam.y, rp.dump.slice.beam.z);
    fprintf('> Plasma slices: every %d step(s) = %d dump(s), ', rp.dump.slice.plasma.freq, numDumps(rp.dump.slice.plasma.freq));
    fprintf('at (%1.f, %1.f, %1.f) um\n', rp.dump.slice.plasma.x, rp.dump.slice.plasma.y, rp.dump.slice.plasma.z);
    fprintf('> Beam phase space: every %d step(s) = %d dump(s), ', rp.dump.phase.beam.freq, numDumps(rp.dump.phase.beam.freq));
    fprintf('1/%1.f of the particles\n', rp.dump.phase.beam.sample);
    fprintf('> Plasma phase space: every %d step(s) = %d dump(s), ', rp.dump.phase.plasma.freq, numDumps(rp.dump.phase.plasma.freq));
    fprintf('1/%1.f of the particles\n', rp.dump.phase.plasma.sample);
    
    
    %% LASER (off)
    rp.laser = false; 
    
    
    %% total number of beams, plasmas and neutrals
    Nbeams = numel(rp.beam);
    Nspecies = numel(rp.plasma);
    Nneutrals = 0;

    
    %% HELPER FUNCTIONS
    
    % helper function for converting booleans to strings
    function str = bool2str(logic)
        if logic; str = 'true'; else str = 'false'; end
    end
    
    % helper function for converting lists to correct strings
    cutstr = @(s) s(1:(end-1));
    list2str = @(l) cutstr(num2str(l,'%G, '));
    num2strG = @(s) num2str(s,'%G');
    


    
    
    %% generate rpinput just like with the jobfile
    % fill in template
    filepath = [tempfilesfolder '/rpinput'];
    fID = fopen(filepath, 'w');
    
    
    %% PART 1: THE START (9 inputs)
    template{1} = fileread([templatesfolder '/rpinput_template_1.txt']); % [9 inputs]
    
    % rarely changed:
    t(2:4) = {rp.sim.dim.x, rp.sim.dim.y, rp.sim.dim.z}; % [um] box size in: {x,y,z}
    t(5:7) = {rp.sim.ind.x, rp.sim.ind.y, rp.sim.ind.z}; % box resolution binary exponent in: 2^{x,y,z}
    
    % never changed:
    t{1} = 1; % number of stages
    t{8} = 'conducting'; % boundary conditions (conducting/periodic)
    
    % auto-generated:
    t{9} = Nbeams; % number of beams
    
    % write Part 1 to file
    t = cellfun(num2strG, t, 'UniformOutput', false);
    
    fprintf(fID, template{1}, t{1}, t{2}, t{3}, t{4}, t{5}, t{6}, t{7}, t{8}, t{9});
    
    
    %% PART 2 (loop): BEAMS (37 inputs)
    template{2} = fileread([templatesfolder '/rpinput_template_2.txt']);
    
    for i = 1:Nbeams
        % often changed:
        t{1} = bool2str(rp.beam{i}.evolve); % beam evolution?
        t(3:5) = {rp.beam{i}.res.x, rp.beam{i}.res.y, rp.beam{i}.res.z}; % beam macro particle resolution in {x, y, z}
        t{6} = rp.beam{i}.charge; % [e] particle charge (1 = positron, -1 = electron)
        t{8} = rp.beam{i}.gamma; % Lorentz factor
        t{9} = rp.beam{i}.N; % number of beam particles
        t{13} = rp.beam{i}.routine; % initialization routine (1-5: see rpinput for guide)
        t{14} = ''; % beam profile file path
        t{15} = bool2str(rp.beam{i}.quiet); % quiet start?
        t(16:18) = {rp.sim.dim.x/2 + rp.beam{i}.offset.x, rp.sim.dim.y/2 + rp.beam{i}.offset.y, rp.beam{i}.offset.z}; % beam centroid from center {x,y} and front {z}
        
        t{20} = list2str(rp.beam{i}.params{2}); % sigma in {x,y,z} (1), sigma in {x,y} and size {z} (2-3), Twiss {ax,bx,ay,by,z} (5)
        t{22} = list2str(rp.beam{i}.params{3}); % emittances, etc.
        t{24} = list2str(rp.beam{i}.params{4}); % z-profile f(z)
        t{26} = list2str(rp.beam{i}.params{5}); % z-profile zs
        
        t{27} = bool2str(false); % use shifter?
        t{29} = list2str(0);
        t{30} = list2str(0);
        t{31} = list2str(0);
        
        t{32} = bool2str(false); % use destroyer?
        t{34} = list2str(0);
        t{35} = list2str(0);
        t{36} = list2str(0);
        
        
        % rarely changed:
        t{7} = 1; % [m_e] particle mass (1 = electron/positron)
        t(10:12) = {0, 0, 0}; % [c] drift velocities (Lorentz betas) in {x, y, z}
        t{37} = bool2str(false); % radiation damping?
        t{2} = 1e7; % minimum number of beam particles
        
        % auto-generated:
        t{19} = numel(rp.beam{i}.params{2}); % size of Init Array 2 
        t{21} = numel(rp.beam{i}.params{3}); % size of Init Array 3
        t{23} = numel(rp.beam{i}.params{4}); % size of Init Array 4
        t{25} = numel(rp.beam{i}.params{5}); % size of Init Array 5
        t{28} = 1;
        t{33} = 1;
        
        
        % write Part 2 to file
        
        t = cellfun(num2strG, t, 'UniformOutput', false);
        fprintf(fID, template{2}, t{1},  t{2},  t{3},  t{4},  t{5},  t{6},  t{7},  t{8},  t{9},  t{10}, ...
                                  t{11}, t{12}, t{13}, t{14}, t{15}, t{16}, t{17}, t{18}, t{19}, t{20}, ...
                                  t{21}, t{22}, t{23}, t{24}, t{25}, t{26}, t{27}, t{28}, t{29}, t{30}, ...
                                  t{31}, t{32}, t{33}, t{34}, t{35}, t{36}, t{37});
        
    end
    
    
    %% PART 3: LASER AND PLASMA PRE-SETUP (4 inputs)
    template{3} = fileread([templatesfolder '/rpinput_template_3.txt']);
    
    % often changed:
    t{1} = bool2str(rp.laser); % laser on?
    t{2} = Nspecies; % number of plasma species
    t{3} = Nneutrals; % number of neutrals
    t{4} = rp.n0; % [cm^-3] plasma density
    
    % write Part 3 to file
    t = cellfun(num2strG, t, 'UniformOutput', false);
    fprintf(fID, template{3}, t{1},  t{2},  t{3},  t{4});
    
    
    %% PART 4 (loop): PLASMA SPECIES (19 inputs)
    template{4} = fileread([templatesfolder '/rpinput_template_4.txt']);
    
    for i = 1:Nspecies
        
        % often changed:
        t{8} = rp.plasma{i}.type; % r-profile type (1 = uniform, 70 = hollow channel)
        t(9:11) = {rp.plasma{i}.r.args(1), rp.plasma{i}.r.args(2), rp.plasma{i}.r.args(3)}; % profile arguments
        
        % rarely change
        t{2} = rp.plasma{i}.res; % number of plasma macro particles per direction in a 2D slice
        t{14} = list2str(rp.plasma{i}.r.ns); % (type 21 only) list of n(r)/n0 (where r=0 is the box middle)
        t{15} = list2str(rp.plasma{i}.r.rs); % (type 21 only) list of corresponding r [um]
        t{16} = bool2str(rp.plasma{i}.z.enable); % use longitudinal profile?
        t{18} = list2str(rp.plasma{i}.z.ns); % list of n(s)/n0
        t{19} = list2str(rp.plasma{i}.z.ss); % list of corresponding s [um]
        
        % almost never change
        t{3} = -1; % [e] plasma particle charge
        t{4} = 1; % [m_e] plasma particle mass
        t(5:6) = {0,0}; % [c] thermal velocity in {x,y}
        t{7} = 1; % ion/electron ratio (1 = quasi-neutral, 0 = electron cloud)
        
        % never change:
        t{1} = -0.08; % load balance TH (mystery: don't change)
        t{12} = 0; % random seed? (mystery: don't change)
        
        % auto-generated:
        t{13} = numel(rp.plasma{i}.r.ns); % (type 21 only) number of radial profile points
        t{17} = numel(rp.plasma{i}.z.ns); % number of long. profile points

        % write Part 4 to file
        t = cellfun(num2strG, t, 'UniformOutput', false);
        fprintf(fID, template{4}, t{1},  t{2},  t{3},  t{4},  t{5},  t{6},  t{7},  t{8},  t{9},  t{10}, ...
                                  t{11}, t{12}, t{13}, t{14}, t{15}, t{16}, t{17}, t{18}, t{19});
    
    end
    
    
    %% PART 5: DUMP PARAMETERS and END (64 inputs)
    template{5} = fileread([templatesfolder '/rpinput_template_5.txt']);
    
    % often changed:
    t{1} = rp.sim.time.end; % [1/omega_p] full simulation time
    t{2} = rp.sim.time.step; % [1/omega_p] time step
    
    t{33} = 0; % full E-field dump freq
    t{34} = rp.dump.slice.field.freq; % E-field slice freq
    t(35:37) = {rp.sim.dim.x/2+rp.dump.slice.field.x, rp.sim.dim.y/2+rp.dump.slice.field.y, rp.dump.slice.field.z}; % E-field slice
    t{38} = 0; % full B-field dump freq
    t{39} = t{34}; % B-field slice freq
    t(40:42) = t(35:37); % B-field slice
    
    t{43} = 0; % full beam density dump freq
    t{44} = rp.dump.slice.beam.freq; % beam density slice freq
    t(45:47) = {rp.sim.dim.x/2+rp.dump.slice.beam.x, rp.sim.dim.y/2+rp.dump.slice.beam.y, rp.dump.slice.beam.z}; % beam density slice
    
    t{50} = 0; % full plasma density dump freq
    t{51} = rp.dump.slice.plasma.freq; % plasma density slice freq
    t(52:54) = {rp.sim.dim.x/2+rp.dump.slice.plasma.x, rp.sim.dim.y/2+rp.dump.slice.plasma.y, rp.dump.slice.plasma.z}; % plasma density slice
    
    t{55} = bool2str(rp.dump.phase.beam.freq>0); % beam phase space dump?
    t{56} = rp.dump.phase.beam.freq; % beam phase space freq
    t{57} = rp.dump.phase.beam.sample; % beam phase space sampling rate
    
    t{58} = bool2str(rp.dump.phase.plasma.freq>0); % plasma phase space dump?
    t{59} = rp.dump.phase.plasma.freq; % plasma phase space freq
    t{60} = rp.dump.phase.plasma.sample; % plasma phase space sampling rate
    
    % rarely changed:
    t(3:12) = {0}; % potential diagnostics
    t(13:17) = {0}; % ponderomotive potential diagnostics
    t(18:22) = {0}; % chi diagnostics
    t(23:32) = {0}; % current diagnostics
    t(61:64) = {bool2str(false), 0, bool2str(false), 0}; % restart file
    t{48} = 0; % dump frequency of beam z-centroid
    t{49} = rp.beam{1}.res.z; % beam centroid diagnostics resolution
    
    % write Part 5 to file
    t = cellfun(num2strG, t, 'UniformOutput', false);
    fprintf(fID, template{5}, t{1},  t{2},  t{3},  t{4},  t{5},  t{6},  t{7},  t{8},  t{9},  t{10}, ...
                              t{11}, t{12}, t{13}, t{14}, t{15}, t{16}, t{17}, t{18}, t{19}, t{20}, ...
                              t{21}, t{22}, t{23}, t{24}, t{25}, t{26}, t{27}, t{28}, t{29}, t{30}, ...
                              t{31}, t{32}, t{33}, t{34}, t{35}, t{36}, t{37}, t{38}, t{39}, t{40}, ...
                              t{41}, t{42}, t{43}, t{44}, t{45}, t{46}, t{47}, t{48}, t{49}, t{50}, ...
                              t{51}, t{52}, t{53}, t{54}, t{55}, t{56}, t{57}, t{58}, t{59}, t{60}, ...
                              t{61}, t{62}, t{63}, t{64});
                          
    %% close file
    fclose(fID);
    
end
