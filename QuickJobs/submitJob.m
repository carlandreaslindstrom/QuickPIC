function [] = submitJob(project, simname, timelimit, tasks, RAM)
    
    % Submit a job to the cluster
    % Author: Carl A. Lindstrom (Uni. Oslo, 25.08.2016)

    % add functions CONFIG and sshfrommatlab*
    addpath('..');
    addpath(CONFIG('sshfrommatlab'));
    tempfilesfolder = [CONFIG('qjobsfolder') '/tempfiles'];
    
    % defaults
    if ~exist('tasks','var'); tasks = 128; end % # of nodes
    if ~exist('RAM','var'); RAM = 1024; end % [MB]
    
    % connect to host
    disp(['Connecting to cluster: ' CONFIG('host')]);
    chnl = sshfrommatlab(CONFIG('username'), CONFIG('host'), CONFIG('password'));
    
    % get QuickPIC and binary folders
    qpfolder = CONFIG('qpfolder');
    projfolder = [qpfolder '/' project];
    binfolder = CONFIG('binfolder');
    
    
    % make the project folder if it does not exist
    fprintf('Project folder exists already: ');
    existcmd = ['if [ -d ' projfolder ' ]; then echo ''1''; else echo ''0''; fi;'];
    [~, output] = sshfrommatlabissue(chnl, existcmd);
    if ~str2double(output) % project does not exist: make it
        sshfrommatlabissue(chnl, ['mkdir ' projfolder]);
    end
    disp(' ');
    
    % find existing simulations (to avoid overwriting)
    disp('Project folder: ');
    [~, sims] = sshfrommatlabissue(chnl, ['ls ' projfolder]);
    
    disp(' ');disp(' ');
    while ismember(simname, sims)
        newname = [simname '_COPY'];
        disp(['Name "' simname '" busy: changing to "' newname '"']);
        simname = newname;
    end
    simfolder = [projfolder '/' simname];
    
    % make new simulation directory
    disp(['Making a new folder: ' simname]);
    sshfrommatlabissue(chnl, ['mkdir ' simfolder]);
    
    % copy QuickPIC executable
    disp('Copying QuickPIC executable.');
    sshfrommatlabissue(chnl, ['cp ' binfolder '/qpic.e ' simfolder]); 
    
    % make job file (define hours and cores)
    disp('Making the job file (locally).');
    jobfilepath = jobfileMaker(project, simname, timelimit, tasks, RAM);
    rpinputpath = [tempfilesfolder '/rpinput'];
    
    % send job and rpinput file
    disp('Sending job and input files to the cluster.');
    sshpasscmd = [CONFIG('sshpass') ' -p ' CONFIG('password') ' '];
    scpcmd = ['scp ' jobfilepath ' ' rpinputpath ' ' CONFIG('username') '@' CONFIG('host') ':' simfolder];
    system([sshpasscmd scpcmd]);
    
    % submit job
    disp('Submitting the job:');
    enableExpress = true; % disable/enable express queueing
    if enableExpress && timelimit <= 2 && tasks <= 8
        subcmd = ['qsub -l express ' simfolder '/qpic.e.cmd'];
        disp('Running express job.');
    else
        subcmd = ['qsub ' simfolder '/qpic.e.cmd'];
    end
    sshfrommatlabissue(chnl, subcmd);
    disp(' ');
    
    % check status
    sshfrommatlabissue(chnl, ['qstat -u ' CONFIG('username')]);
    disp(' ');
    
    % close connection to host
    sshfrommatlabclose(chnl);
    disp(' ')
    disp('Connection closed.')
    
end

