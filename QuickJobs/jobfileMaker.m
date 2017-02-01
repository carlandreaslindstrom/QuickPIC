function filepath = jobfileMaker(project, simname, timelimit, tasks, totalRAM, highp) % time in hours, RAM in GB
    
    % add path to CONFIG file
    addpath('../');
    templatesfolder = [CONFIG('qjobsfolder') '/templates'];
    tempfilesfolder = [CONFIG('qjobsfolder') '/tempfiles'];
    
    % defaults
    if ~exist('tasks','var'); tasks = 128; end
    if ~exist('totalRAM','var'); totalRAM = 64; end % [GB] default: 0.5 GB per node
    if ~exist('highp','var'); highp = false; end
    
    % extract RAM per node (rounded to nearest MB)
    RAM = round(totalRAM*1024/tasks); % [MB] 
    
    % prepare parameters
    username = CONFIG('username');
    simfolder = [CONFIG('qpfolder') '/' project '/' simname];
    RAM = num2str(RAM);
    hours = num2str(floor(timelimit)); % [hours]
    minutes = num2str(floor(mod(timelimit,1)*60),'%02d'); % [minutes]
    timestamp = datetime('now','Format','eee MMM dd HH:mm:ss yyyy','TimeZone','America/Los_Angeles');
    if tasks <= 8
        sharing = 'shared'; % single-node
    else
        sharing = 'dc*'; % multi-node
    end
    tasks = num2str(tasks);
    
    % get template to fill in
    if highp
        template = fileread([templatesfolder '/jobtemplate_highp.txt']);
    else
        template = fileread([templatesfolder '/jobtemplate.txt']);
    end
    
    % fill in template
    filepath = [tempfilesfolder '/qpic.e.cmd'];
    fID = fopen(filepath, 'w');
    fprintf(fID, template, simfolder, tasks, RAM, hours, minutes, username, timestamp, sharing);
    fclose(fID);
    
end

