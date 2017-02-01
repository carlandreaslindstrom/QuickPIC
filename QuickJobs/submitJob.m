function [chnl] = submitJob(project, simname, timelimit, tasks, RAM, chnl, close)
    
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
    if ~exist('chnl','var')
        disp(['Connecting to cluster: ' CONFIG('host')]);
        chnl = sshfrommatlab(CONFIG('username'), CONFIG('host'), CONFIG('password'));
    end
    
    % get QuickPIC and binary folders
    qpfolder = CONFIG('qpfolder');
    projfolder = [qpfolder '/' project];
    binfolder = CONFIG('binfolder');
    
    % make the project folder if it does not exist
    existcmd = ['if [ -d ' projfolder ' ]; then echo ''1''; else echo ''0''; fi;'];
    output = evalc('sshfrommatlabissue(chnl, existcmd);');
    if ~str2double(output) % project does not exist: make it
        disp('Project does not exist yet.');
        sshfrommatlabissue(chnl, ['mkdir ' projfolder]);
    end
    
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
    rp = rpinputParser('','',true);
    if rp.plasma{1}.type==71 % needs a special qpic.e-file
        sshfrommatlabissue(chnl, ['cp ' binfolder '/qpic.e.hca ' simfolder '/qpic.e']); 
    else % just a normal one
        sshfrommatlabissue(chnl, ['cp ' binfolder '/qpic.e ' simfolder]); 
    end
    
    % is express job?
    enableExpress = true; % disable/enable express queueing
    isExpressJob = enableExpress && timelimit <= 2 && tasks <= 8;
    highPri = false; % high priority queue 
    
    % make job file (define hours and cores)
    disp('Making the job file (locally).');
    jobfilepath = jobfileMaker(project, simname, timelimit, tasks, RAM, highPri);
    rpinputpath = [tempfilesfolder '/rpinput'];
    
    % send job and rpinput file
    disp('Sending job and input files to the cluster.');
    sshpasscmd = [CONFIG('sshpass') ' -p ' CONFIG('password') ' '];
    scpcmd = ['scp ' jobfilepath ' ' rpinputpath ' ' CONFIG('username') '@' CONFIG('host') ':' simfolder];
    system([sshpasscmd scpcmd]);
    
    % submit job
    disp('Submitting the job:');
    if isExpressJob
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
    if ~exist('chnl','var') || (exist('close','var') && close)
        sshfrommatlabclose(chnl);
        disp(' ')
        disp('Connection closed.')
    end
    
end

