function val = CONFIG(key)
    
    % IMPORTANT NOTE TO GET STARTED:
    % 1. CHANGE ALL THE BELOW INFO TO FIT YOUR ACCOUNTS AND FOLDERS, ETC.
    % 2. DOWNLOAD AND INSTALL sshpass (MAC: 'sudo port install sshpass')
    %      (this is an add-on used to connect to the cluster throught ssh 
    %       without constantly asking for passwords).
    % 3. CHANGE THE NAME OF THIS FILE FROM 'CONFIG_public.m" TO "CONFIG.m"
    %      (it is used as a function by the framework to get info on your
    %      personal configuration).
    % 4. REMEMBER TO NOT UPLOAD THIS FILE ANYWHERE (like GitHub), AS IT 
    %      CONTAINS PRIVATE INFORMATION.
    
    % User-changable parameters for running QuickPlot and QuickJobs
    % Author: Carl A. Lindstrom (Uni. Oslo, 25.08.2016)

    % CLUSTER LOGIN INFO
    conf.username = 'your_username'; % username on cluster
    conf.password = 'your_password'; % password (don't share!)
    conf.host = 'hoffman2.idre.ucla.edu'; % host server name
    
    % LOCAL INFO
    % Note: install 'sshpass' (which automates ssh passwords) using 
    %  'sudo port install sshpass' on Macs
    %  'sudo apt-get install sshpass' on Linux (TBC)
    conf.mainpath = '/Users/calind/PhD/QuickPIC'; % location of the framework (and this file)
    conf.sshpass = '/opt/local/bin/sshpass'; % location of sshpass (locally)
    conf.outputs = [conf.mainpath '/outputs']; % location of local output folder (put wherever it suits you)
    conf.qplotfolder = [conf.mainpath '/QuickPlot']; % location of local QuickPlot folder (don't change)
    conf.qjobsfolder = [conf.mainpath '/QuickJobs']; % location of local QuickJobs folder (don't change)
    conf.sshfrommatlab = [conf.qjobsfolder '/sshfrommatlab']; % location of sshfrommatlab (don't change: included in the framework)
    
    % REMOTE INFO
    conf.qpfolder = '/u/d2scratch/your_username/your_QP_folder'; % location of QuickPIC folder where outputs are generated
    conf.binfolder = '/u/project/mori/your_username/bin'; % location of the "qpic.e" executable to be copied to the project folder
    
    
    % return the requested config value
    val = conf.(key);
    
end