function [] = getProject(project)
    
    % Download finished project folder from the cluster
    % Author: Carl A. Lindstrom (Uni. Oslo, 25.08.2016)
    
    % add functions CONFIG and sshfrommatlab*
    addpath('..');
    addpath(CONFIG('sshfrommatlab'));
    
    % prepare simulation directory on host
    qpfolder = CONFIG('qpfolder');
    projfolder = [qpfolder '/' project];
    tarpath = [projfolder '.tqz'];
    outputfolder = CONFIG('outputs');
    
    % connect to host
    fprintf(['Connecting to cluster: ' CONFIG('host') '\n\n']);
    chnl = sshfrommatlab(CONFIG('username'), CONFIG('host'), CONFIG('password'));
    
    % compress file to tarball (.tqz) on host
    fprintf('Compressing simulation:');
    existcmd = ['[ -d "' projfolder '" ] && '];
    tarcmd = ['tar -cvzf ' tarpath ' -C ' projfolder ' . --exclude ' project '.tqz --exclude qpic.e'];
    [~, output] = sshfrommatlabissue(chnl, [existcmd tarcmd]);
    
    % if no
    if size(output,1) == 1 && strcmp(output{1},'') 
        % SIMULATION DOES NOT EXIST
        fprintf(['\nThe requested project (' project ') does not exist.\n']);
        
    else
        % SIMULATION EXISTS
            
        % display the tarball size
        ducmd = ['du -h ' tarpath];
        sshfrommatlabissue(chnl, ducmd);
        
        % download the tarball to local machine
        fprintf(['\n\nDownloading file from ' CONFIG('host') '... ']);
        sshpasscmd = [CONFIG('sshpass') ' -p ' CONFIG('password') ' '];
        projfolder_local = [outputfolder '/' project];
        scpcmd = ['scp ' CONFIG('username') '@' CONFIG('host') ':' tarpath ' ' projfolder_local '.tqz'];
        system([sshpasscmd scpcmd]);
        disp('Done');

        % extract the tarball locally and then delete it
        fprintf('Extracting file locally... ');
        tarpath_local = [outputfolder '/' project '.tqz'];
        mkdir(projfolder_local);
        system(['tar -xzf ' tarpath_local ' -C ' projfolder_local]);
        disp('Done');
        disp('Deleting the compressed file locally.');
        system(['rm -r ' tarpath_local]);
        
        % delete tarball on cluster
        fprintf('Deleting compressed file on cluster... ');
        sshfrommatlabissue(chnl, ['rm -r ' tarpath]);
        disp('Done');
    
    end
    
    % close connection to host
    sshfrommatlabclose(chnl);
    fprintf('\nConnection closed. \n');
    
end

