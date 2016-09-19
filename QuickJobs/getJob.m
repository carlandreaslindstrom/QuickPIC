function [] = getJob(project, simname)
    
    % Download finished job output from the cluster
    % Author: Carl A. Lindstrom (Uni. Oslo, 25.08.2016)
    
    % add functions CONFIG and sshfrommatlab*
    addpath('..');
    addpath(CONFIG('sshfrommatlab'));
    
    % prepare simulation directory on host
    qpfolder = CONFIG('qpfolder');
    simfolder = [qpfolder '/' project '/' simname];
    tarpath = [simfolder '/' simname '.tqz'];
    outputfolder = CONFIG('outputs');
    
    % connect to host
    fprintf(['Connecting to cluster: ' CONFIG('host') '\n\n']);
    chnl = sshfrommatlab(CONFIG('username'), CONFIG('host'), CONFIG('password'));
    
    % compress file to tarball (.tqz) on host
    fprintf('Compressing simulation:');
    existcmd = ['[ -d "' simfolder '" ] && '];
    tarcmd = ['tar -cvzf ' tarpath ' -C ' simfolder ' . --exclude ' simname '.tqz --exclude qpic.e'];
    [~, output] = sshfrommatlabissue(chnl, [existcmd tarcmd]);
    
    % if no
    if size(output,1) == 1 && strcmp(output{1},'') 
        % SIMULATION DOES NOT EXIST
        fprintf(['\nThe requested simulation (' project '/' simname ') does not exist.\n']);
        
    else
        % SIMULATION EXISTS
            
        % display the tarball size
        ducmd = ['du -h ' tarpath];
        sshfrommatlabissue(chnl, ducmd);

        % download the tarball to local machine
        fprintf(['\n\nDownloading file from ' CONFIG('host') '... ']);
        sshpasscmd = [CONFIG('sshpass') ' -p ' CONFIG('password') ' '];
        scpcmd = ['scp ' CONFIG('username') '@' CONFIG('host') ':' tarpath ' ' outputfolder];
        system([sshpasscmd scpcmd]);
        disp('Done');

        % extract the tarball locally and then delete it
        fprintf('Extracting file locally... ');
        simfolder_local = [outputfolder '/' project '/' simname];
        tarpath_local = [outputfolder '/' simname '.tqz'];
        mkdir(simfolder_local);
        system(['tar -xzf ' tarpath_local ' -C ' simfolder_local]);
        disp('Done');
        disp('Deleting the compressed file locally.');
        system(['rm -r ' tarpath_local]);
    
    end
    
    % close connection to host
    sshfrommatlabclose(chnl);
    fprintf('\nConnection closed. \n');
    
end

