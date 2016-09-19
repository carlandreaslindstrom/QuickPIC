function [data, exitflag] = h5read( filename )

    % READS h5 FILES
    % Based on script by Erik Adli
    
    char_slash = findstr(filename, '/');
    char_dot = findstr(filename, '.');
    dataname = filename(char_slash(end):char_dot(end)-1);
    dataname( dataname == '_' ) = '-'; % for some reason QuickPIC has slightly different filename than dataname

    try
        data = hdf5read(filename, dataname);
        exitflag = true;
    catch
        data = [];
        exitflag = false;
    end

end

