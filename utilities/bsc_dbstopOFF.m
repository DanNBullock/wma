function bsc_dbstopOFF()
% Feel free to comment this function out in your code.  It's purpose is to either
% turn debug mode on or off depending on whether or not an interactive session is
% detected.
    jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
    retval = ~isempty(jDesktop.getClient('Command Window'));
    if ~retval
     dbclear if error
     fprintf('\n debug mode now turned off');
    else
        fprintf('\n Interactive seession detected, leaving debug mode on \n')
    end
end
