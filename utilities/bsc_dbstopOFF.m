function bsc_dbstopOFF()
    jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
    retval = ~isempty(jDesktop.getClient('Command Window'));
    if ~retval
     dbclear if error
     fprintf('\n debug mode now turned off');
    else
        fprintf('\n Interactive seession detected, leaving debug mode on \n')
    end
end