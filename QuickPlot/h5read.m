function [data] = h5read( filename )

    % READS h5 FILES
    try
        char_slash = [1 strfind(filename, '/')];
        char_dot = strfind(filename, '.');
        dataname = filename(char_slash(end):char_dot(end)-1);
        dataname( dataname == '_' ) = '-'; % QuickPIC requires different dataname than filename
        data = hdf5read(filename, dataname);
    catch
        data = false;
    end

end

