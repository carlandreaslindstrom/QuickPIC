function [ phaseSpace ] = readBeams( simfolder, step, rp, beams )
    
    % constants
    SI_c = 299792458; % [m/s] speed of light
    SI_e = 1.60217662e-19; % [C] electron charge
    SI_me = 9.10938356e-31; % [kg] electron mass
    SI_eps0 = 8.85418782e-12; % [m^-3 kg^-1 s^4 A^2] permittivity of free space
    
    % adds together the beams listed (all by default)
    Nbeams = numel(rp.beam);
    if ~exist('beams','var')
        beams = 1:Nbeams;
    end
    
    k_p = sqrt(rp.n0 * SI_e^2/(SI_me*SI_eps0))/SI_c; % plasma wavenumber
    
    x0 = rp.beam{1}.offset.x - rp.sim.dim.x/2; % initial x offset
    y0 = rp.beam{1}.offset.y - rp.sim.dim.y/2; % initial y offset
    
    try
        phaseSpace = struct('X',[],'Y',[],'Z',[],'XP',[],'YP',[],'E',[], 'Q', []);
        for i = beams
            filename = [simfolder '/RAW-BEAM/' num2str(i,'%02d') '/RAW-BEAM-' num2str(i,'%02d') '_' num2str(step,'%04d') '.h5'];
            
            phaseSpace.X = [phaseSpace.X; x0 + hdf5read(filename, '/x1')/k_p*1e6]; % [um]
            phaseSpace.Y = [phaseSpace.Y; y0 + hdf5read(filename, '/x2')/k_p*1e6]; % [um]
            phaseSpace.Z = [phaseSpace.Z; hdf5read(filename, '/x3')/k_p*1e6]; % [um]
            phaseSpace.XP = [phaseSpace.XP; hdf5read(filename, '/p1')/k_p*1e6]; % [urad]
            phaseSpace.YP = [phaseSpace.YP; hdf5read(filename, '/p2')/k_p*1e6]; % [urad]
            phaseSpace.E = [phaseSpace.E; hdf5read(filename, '/p3') * SI_me*SI_c^2/(1e9*SI_e)]; % [GeV]
            phaseSpace.Q = [phaseSpace.Q; ones(size(phaseSpace.X)) * rp.beam{1}.charge]; % [e] charge
        end
    catch
        phaseSpace = false;
    end

end

