function rp = rpinputParser(project, sim)

    % Convert a job parameter input file (rpinput) to a usable structure (rp)
    % Author: Carl A. Lindstrom (Uni. Oslo, 5.9.2016)
    
    % add CONFIG to the path
    addpath('..');
    
    % read file into a cell array
    outputfolder = CONFIG('outputs');
    filename = [outputfolder '/' project '/' sim '/rpinput'];
    %fileID = fopen(filename,'r');
    %a = textscan(fileID,'%s',10000,'Delimiter','\n');
    
    clean0 = fileread(filename); % read text file into character array
    clean1 = regexprep(clean0, ',[\s]{0,}([A-Z]+)',',\n $1'); % insert appropriate newlines
    %clean2 = regexprep(b, '[\s]{0,}=[\s]{0,}',' = '); % fix ugly equal signs (not important)
    a = textscan(clean1,'%s',10000,'Delimiter','\n');
    carr = a{1};
    
    % getting indices of a certain class
    nth = @(n,list) list(n);
    findiN = @(arr, word, n) nth(n,find(not(cellfun('isempty',strfind(arr, word))),n));
    findi = @(arr, word) findiN(arr, word, 1);
    indexof = @(word,n) findiN(carr, word, n);% nth(n,find(not(cellfun('isempty',strfind(carr, word))),n));
    firstIndexAfter = @(word, line) (line-1) +findi(carr(line:end), word); % @(word, line) (line-1) + find(not(cellfun('isempty',strfind(carr(line:end), word))),1);
    indicesN = @(word, n) (indexof(word,n)+1):(firstIndexAfter('/', indexof(word,n))-1);
    %indices = @(word) indicesN(word,1);
   
    cell2str = @(c) c{:};
    delEmpties = @(list) list(~cellfun(@isempty,list));
    str2bool = @(str) strcmp(str,'true');
    paramsN = @(section, param, n) nth(2,strtrim(strsplit(carr{nth(1,indicesN(section,n))-1+findi(carr(indicesN(section,n)),param)},'=')));
    params2listN = @(section, param, n) cellfun(@str2double, delEmpties(strtrim(strsplit(cell2str(paramsN(section, param, n)),','))));
    params2list = @(section, param) params2listN(section, param, 1); % cellfun(@str2double, delEmpties(strtrim(strsplit(cell2str(params(section, param)),','))));
    params2boolN = @(section, param, n) str2bool(delEmpties(strsplit(cell2str(strsplit(cell2str(paramsN(section, param, n)),',')),'.')));
    params2bool = @(section, param) params2boolN(section, param, 1);
    params2strN = @(section, param, n) cell2str(strrep(delEmpties(strsplit(cell2str(paramsN(section, param, n)),',')), '''', ''));
    %params2str = @(section, param) params2strN(section, param, 1);
    
    % BEAMS
    %Nbeams = params2list('&Num_Beams','NBeams'); % get beam number from rpinput
    Nbeams = 0; % get beam number by counting
    for i = 1:size(carr,1)
        Nbeams = Nbeams + sum(strfind(carr{i}, '&Beam')) - sum(strfind(carr{i}, '&Beam_'));
    end
    for i = 1:Nbeams
        rp.beam{i}.evolve = params2boolN('&Beam', 'BEAM_EVOLUTION', i);
        rp.beam{i}.quiet = params2boolN('&Beam', 'QUIET_START', i);
        rp.beam{i}.gamma = params2listN('&Beam', 'Gamma', i);
        rp.beam{i}.N = params2listN('&Beam', 'Num_Particle', i);
        rp.beam{i}.charge = params2listN('&Beam', 'Charge', i);
        param_arr3 = params2listN('&Beam', 'Parameter_Array(3', i);
        rp.beam{i}.emit.x = param_arr3(1);
        rp.beam{i}.emit.y = param_arr3(2);
        rp.beam{i}.espread = param_arr3(3);
        rp.beam{i}.file = params2strN('&Beam', 'BEAM_PROFILE', i);
        rp.beam{i}.params{2} = params2listN('&Beam', 'Parameter_Array(2', i);
        rp.beam{i}.params{3} = param_arr3;
        rp.beam{i}.params{4} = params2listN('&Beam', 'Parameter_Array(4', i);
        rp.beam{i}.params{5} = params2listN('&Beam', 'Parameter_Array(5', i);
        rp.beam{i}.routine = params2listN('&Beam', 'Init_Routine', i);
        param_arr1 = params2listN('&Beam', 'Parameter_Array(1', i);
        rp.beam{i}.offset.x = param_arr1(1);
        rp.beam{i}.offset.y = param_arr1(2);
        rp.beam{i}.offset.z = param_arr1(3);
        rp.beam{i}.res.x = params2listN('&Beam', 'NPX', i);
        rp.beam{i}.res.y = params2listN('&Beam', 'NPY', i);
        rp.beam{i}.res.z = params2listN('&Beam', 'NPZ', i);
    end
    
    % PLASMAS
    rp.n0 = params2list('&Plasma','Plasma_Density') * 1e6; % convert density to [m^-3]
    % Nspecies = params2list('&Plasma','Nspecies'); % get beam number from the file
    Nspecies = 0; % get beam number by counting
    for i = 1:size(carr,1)
        Nspecies = Nspecies + sum(strfind(carr{i}, '&Species')) - sum(strfind(carr{i}, '&Species_'));
    end
    for i = 1:Nspecies
        rp.plasma{i}.type = params2listN('&Species', 'Profile_type', i);
        rp.plasma{i}.r.args = [params2listN('&Species','argx1',i), params2listN('&Species','argx2',i), params2listN('&Species','argx3',i)];
        rp.plasma{i}.r.ns = params2listN('&Species','Prof_Parameter(1',i);
        rp.plasma{i}.r.rs = params2listN('&Species','Prof_Parameter(2',i);
        rp.plasma{i}.z.enable = params2boolN('&Species','Density_Variation',i);
        rp.plasma{i}.z.ns = params2listN('&Species','Density_Variation_Fs',i);
        rp.plasma{i}.z.ss = params2listN('&Species','Density_Variation_s',i);
        rp.plasma{i}.res = params2listN('&Species','NP2',i);
    end
    
    % useful plasma numbers
    SI_c = 299792458; % [m/s] speed of light
    SI_e = 1.60217662e-19; % [C] electron charge
    SI_me = 9.10938356e-31; % [kg] electron mass
    SI_eps0 = 8.85418782e-12; % [m^-3 kg^-1 s^4 A^2] permittivity of free space
    omega_p = sqrt(rp.n0 * SI_e^2/(SI_me*SI_eps0));
    k_p = SI_c/omega_p;
    
    % SIMULATION BOX
    rp.sim.dim.x = params2list('&Simulation_Sys','Box_X');
    rp.sim.dim.y = params2list('&Simulation_Sys','Box_Y');
    rp.sim.dim.z = params2list('&Simulation_Sys','Box_Z');
    rp.sim.ind.x = params2list('&Simulation_Sys','INDX');
    rp.sim.ind.y = params2list('&Simulation_Sys','INDY');
    rp.sim.ind.z = params2list('&Simulation_Sys','INDZ');
    rp.sim.time.end = params2list('&Simulation_time','TEND');
    rp.sim.time.step = params2list('&Simulation_time','DT');
    rp.sim.dist.end = rp.sim.time.end * k_p;
    rp.sim.dist.step = rp.sim.time.step * k_p;
    
    % DUMP (assuming all slice frequencies are equal, basing on E-field only)
    rp.dump.slice.field.freq = params2list('&Field_Diag','DFESLICE'); % assuming DFBSLICE is the same
    rp.dump.slice.field.x = params2list('&Field_Diag','EX0');
    rp.dump.slice.field.y = params2list('&Field_Diag','EY0');
    rp.dump.slice.field.z = params2list('&Field_Diag','EZ0');
    rp.dump.slice.beam.freq = params2list('&Beam_Diag','DFQEBSLICE');
    rp.dump.slice.beam.x = params2list('&Beam_Diag','QEBX0');
    rp.dump.slice.beam.y = params2list('&Beam_Diag','QEBY0');
    rp.dump.slice.beam.z = params2list('&Beam_Diag','QEBZ0');
    rp.dump.slice.plasma.freq = params2list('&Plasma_Diag','DFQEPSLICE');
    rp.dump.slice.plasma.x = params2list('&Plasma_Diag','QEPX0');
    rp.dump.slice.plasma.y = params2list('&Plasma_Diag','QEPY0');
    rp.dump.slice.plasma.z = params2list('&Plasma_Diag','QEPZ0');
    rp.dump.phase.beam.freq = params2list('&Beam_Phase_Space_Diag','DFPHA_BEAM')*params2bool('&Beam_Phase_Space_Diag','DUMP_PHA_BEAM');
    rp.dump.phase.beam.sample = params2list('&Beam_Phase_Space_Diag','DSAMPLE_BEAM');
    rp.dump.phase.plasma.freq = params2list('&Plasma_Phase_Space_Diag','DFPHA_PLASMA')*params2bool('&Plasma_Phase_Space_Diag','DUMP_PHA_PLASMA');
    rp.dump.phase.plasma.sample = params2list('&Plasma_Phase_Space_Diag','DSAMPLE_PLASMA');
    
    % LASER
    rp.laser = params2bool('&laser_input', 'laser_on');
    
end

