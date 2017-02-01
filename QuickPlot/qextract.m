function [ qp ] = qextract( project, sim, step )
    
    % find output folder
    addpath('..');
    outputfolder = CONFIG('outputs');
    
    % get input parameters
    rp = rpinputParser(project, sim);
    
    % axes and offsets
    axes = qaxes(project, sim);
    
    % extraction helper function
    simpath = [outputfolder '/' project '/' sim];
    geth5 = @(var) h5read([simpath '/' var '/' var '_' num2str(step,'%04d') '.h5']);

    % find available projections
    for proj = {'XY','XZ','YZ'}

        % fields
        for field = {'E','B'}
            for dim = {'X','Y','Z'}
                % skip if no folders exist
                if numel(dir([simpath '/F' field{:} dim{:} '-' proj{:}]))>0
                    qp.(field{:}).(dim{:}).(proj{:}) = geth5(['F' field{:} dim{:} '-' proj{:}]);
                end
            end
        end

        % plasmas
        for species = 1:numel(rp.plasma)
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

    % phase space
    qp.Q.total.phase = readBeams(simpath, step, rp); % all beams

    % add forces and fields in polar coordinates
    [qp.E, qp.B, qp.F] = polarFields(axes.X, axes.X, qp.E, qp.B);

end

