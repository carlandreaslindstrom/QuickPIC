function [] = checkJobs()

    % Check the status of jobs on the cluster
    % Author: Carl A. Lindstrom (Uni. Oslo, 25.08.2016)
    
    % add functions and sshfrommatlab*
    addpath('..');
    addpath(CONFIG('sshfrommatlab'));
    
    % connect to host
    fprintf(['Connecting to cluster: ' CONFIG('host') '\n']);
    chnl = sshfrommatlab(CONFIG('username'), CONFIG('host'), CONFIG('password'));
    
    % check status
    [~, output] = sshfrommatlabissue(chnl, ['qstat -u ' CONFIG('username')]);
    if(numel(strtrim(output))==0)
        disp('No jobs in the queue.');
    else
        disp(' ');
    end
    
    % close connection to host
    sshfrommatlabclose(chnl);
    disp(' ')
    disp('Connection closed.')
    
end

