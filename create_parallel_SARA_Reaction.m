% function create_parallel(parallel_fun, job_list1, parallel_fun_dir, job_name)
%
% This function serves to create parallel Jobs for a given script
% to parallalise it for a given list of jobs.
%
% Therefore, it creates executable sh-scripts, which execute each single
% job, as well as one script which passes all the single job scripts to
% condor (the framework which handles the parallelisation).
%
% This version requsts 600MB free on the host machine. If your job does not
% require free memory (also matlab already needs quite a bit), use
% create_parallel_nomem.m.
%
% IMPORTANT: The function will only run under linux, so all pathes must be
% written in linux style!!!
% Modifying it so that it will
%
% PARAMETERS
%   parallel_fun: name of the function (e.g. 'example_function', NOT
%               'example_script.m', 'example_function()', etc.). This
%               function must take exactly 1 argument as input. For an
%               example, see parallel_example_fun.m;
%
%   job_list1: 1xn vector of elements which will be one after the other
%               passed to the function parallel_fun.
%
%   parallel_fun_dir: path to parallel_fun, will be used to add to the
%               matlab path on the remote machines (the nathans)
%
% OPTIONAL
%   job_name: string to identify your job. Subdirectories with this name will
%       be created insed the specified batch directory and the log
%       directory, so that multiple jobs can be executed at the same time.
%       If no jobname is given, the name of prallel_fun together with the
%       current date and time will be used as job name.
%
% !!! PLEASE GO THROUGH THE TOP OF THIS FUNCTION TO CHECK FURTHER PARAMETERS !!!
%
% Example call (IMPORTANT: last parameter is path to parallel_example_fun.m, NOT to matlab_condor_parallel):
%
%   create_parallel('parallel_example_fun', [1:6, 9], '/analysis/kai/PATH_TO_PARALLEL_EXAMPLE_FUN/') % CHANGE last parameter to the path to your function
%
% TODO: split the logdir variable into windows and linux, in case you want
% to use it under windows.
%
% Kai Goergen, 2009/11/17

% History:
%	Kai Goergen, 2014/11/02
%       Added memory request of 600MB as default
%	Kai Goergen, 2013/06/25
%		Changed data01 -> analysis
%	Kai Goergen, 2013/06/24
%		Adapted so that matlab2013 is now used by default
%   Kai Goergen, 2010/01/19
%       Changed position of project name, to have logs & batchs in one
%       common folder
%       Changed in condor_generate_submit.sh that condor logs will now be
%       written to ~/condor/condor_logs/log*
%   Kai Goergen, 2009/11/17
%       Added job number to send_all_file
%   Kai Goergen, 2009/11/12
%       Added job_name to file single job names
%       Changed naming to avoid confusion:
%           start_all_jobs => send_condor_all_jobs.sh
%           start_jobxx_jobname => run_jobxx_jobname
%       Added quick start version for matlab
%   Kai Goergen, 2009/11/06
%       Added check for parallel function ends with .m or (), and remove it
%   Kai Goergen, 2009/10/16
%       Added job_name as parameter, display log directory at the end now


function create_parallel_SARA_Reaction(parallel_fun, job_list1, job_list2, parallel_fun_dir, job_name)
%% basic input checks
% check if parallel_fun ends with .m or ()
if length(parallel_fun) >= 2 && (strcmp(parallel_fun(end-1:end), '.m') || strcmp(parallel_fun(end-1:end), '()'))
    parallel_fun = parallel_fun(1:end-2);
end
miceopt = {'Alladin','Chief','Esmeralda','Frey'} %options for mice

%% set parameters

% default name for jobs
if ~exist('job_name', 'var')
    job_name = [parallel_fun '_' datestr(now, 'yyyymmddTHHMMSS')];
end

% set directory where the sh-files should be created
% will be created, if it does not exist
project_dir = '/home/ehvbeest/SVM_REACTDEC' % must be the ABSOLUT path

% batch files will be written to job_name/batchdirs
% batch_dir = '~/condor/batch_dir/'  % must be the ABSOLUT path
batch_dir = [project_dir '/' job_name '/batch_dir'] % add jobname

% set log folder
% will be created, if it does not exist
%log_file_dir = '~/condor/logs/'  % use ABSOLUT path (or e.g. starting with ~, which is your home)
log_file_dir = [project_dir '/' job_name '/logs/'] % add jobname

%% location of scripts
% set location of execute_matlab_process.sh
execute_matlab_process_sh = '"$TMPDIR"/BashScripts/execute_Compiledmatlab_process_SARA_Reaction.sh' % must be the ABSOLUTE path
% set location of execute_matlab_process_sh
generate_submit = '"$TMDPIR"/BashScripts/SubmitterOfAll.sh' % must be the ABSOLUT path


%% PROCESSING STARTS FROM HERE (no more parameters to check)
%% create batch & log folder

display('create batch & log folder')
[success, message] = mkdir(batch_dir);
if ~success
    error(['Could not create directory for batch_dir: ' message])
end

if ispc
    error('As the path to the log directory is different under linux and windows, this is a running this script and Windows will not work at the moment. Please run it from a Linux machine (sorry)')
else
    [success, message] = mkdir(log_file_dir);
    if ~success
        error(['Could not create directory for log_file_dir: ' message])
    end
end

%% Check if main batch files already exists

% set initial behaviour, if you want to overwrite sh-files which are there
overwrite_file = 'ask';

display('creating batch files')

% check if main sh-file to start all jobs exists
filename_all = sprintf('send_all_jobs.sh');
fullfilename_all = [batch_dir '/' filename_all];
if exist(fullfilename_all, 'file')
    display(' ')
    display(['File ' fullfilename_all ' already exist.'])
    overwrite_file = input(['Should it be overwritten? [y, n, a (all)]: '], 's');
    if ~(strcmpi(overwrite_file, 'y') || strcmpi(overwrite_file, 'a'))
        error(['File ' filename_all ' already exists and should not be overwritten. Solve problem and start again.'])
    end
    delete(fullfilename_all)
end

%% create the batch files

% The main batch file handles passing the single job batch files to condor
% (using condor_generate_submit.sh)
display(['Creating main batch file: ' fullfilename_all])
fid_commit_all = fopen(fullfilename_all, 'w');

% ensure that the right shell is used !#/bin/sh
fprintf(fid_commit_all, '#!/bin/bash\n');
% add comment that THIS file submits the stuff to condor
fprintf(fid_commit_all, '#\n');
fprintf(fid_commit_all, '#This bash-script submits all jobs to the server, instead of running them locally.\n');
fprintf(fid_commit_all, ['#If you want to submit only some jobs to the server, simply add a "#" in front of \n' ...
    '#the ones you like to ommit and execute the script then.\n']);
fprintf(fid_commit_all, '#\n');

% create all single job batchfiles, and add for each a call in the main
% batch file
for job_ind = 1:length(job_list1)
    for job_ind2 = 1:length(job_list2)
            curr_job = [job_list1(job_ind), job_list2(job_ind2)];            
            %% create batchfile for current job
            % create/overwrite file
            filename = sprintf('run_job%02i_%02i_%s.sh', curr_job(1), curr_job(2),  job_name);
            fullfilename = [batch_dir '/' filename]
            
            display(['Creating Batch file for Job ' num2str(curr_job(1)) '_' num2str(curr_job(2)) ': ' fullfilename])
            if exist(fullfilename, 'file')
                if ~strcmpi(overwrite_file, 'a')
                    display(' ')
                    display(['File ' fullfilename ' already exist.'])
                    overwrite_file = input(['Should it be overwritten? [y, n, a (all)]: '], 's');
                    if ~(strcmpi(overwrite_file, 'y') || strcmpi(overwrite_file, 'a'))
                        error(['File ' filename ' already exists and should not be overwritten. Solve problem and start again.'])
                    end
                end
                delete(fullfilename)
            end
            
            % open single subject file
            fid_single = fopen(fullfilename , 'w');
            
            % ensure that the right shell is used !#/bin/bash
            fprintf(fid_single, '#PBS -S /bin/bash\n');
            % add a comment what this script does
            jobnameline = ['#PBS -N ' job_name '\n']
            try
                fprintf(fid_single, jobnameline);
            catch ME
                display(ME)
            end
            
            fprintf(fid_single, '#PBS -j oe\n');
            fprintf(fid_single, '#PBS -lnodes=1:mem64gb\n');
            fprintf(fid_single, '#PBS -lwalltime=72:00:00\n');
            fprintf(fid_single, '#PBS -o $HOME/logs/\n');
            fprintf(fid_single, '#\n');
            
            fprintf(fid_single,['cp -r $HOME/TMPData/' miceopt{job_ind} '* "$TMPDIR"\n']);
            fprintf(fid_single,'cp -r $HOME/BashScripts* "$TMPDIR"\n');
            fprintf(fid_single,'cp -r $HOME/MVPA_Scripts* "$TMPDIR"\n');
            fprintf(fid_single,'cp -r $HOME/SVM_STIMDEC "$TMPDIR"\n');
            % get the command to start the job
            % this command will be saved in the job script
            
            fprintf(fid_single,'cd "$TMPDIR"\n');
            fprintf(fid_single,['chmod +x ' execute_matlab_process_sh '\n'])
            line = sprintf('%s %s %i %i %s %s', execute_matlab_process_sh, parallel_fun, curr_job(1), curr_job(2),  log_file_dir, parallel_fun_dir);
            fprintf(fid_single, '%s\n', line);
            
            % finally: pass exit status of execute_matlab_process.sh to condor
            fprintf(fid_single, 'exit $?\n');
            
            fclose(fid_single);
            
            %% add the produced single job batch file to the main batch file
            display(['Adding ' fullfilename ' to main batch file.'])
            % add this script to the list off all scripts
            % be careful to use the linux path
            fprintf(fid_commit_all, '# Job %i_%i\n', curr_job(1), curr_job(2));
            line = sprintf('%s %s', 'qsub ', fullfilename);
            fprintf(fid_commit_all, '%s\n\n', line);
    end
end
fclose(fid_commit_all);

%% files written successfully, instruct user how to proceed

% display(' ')
% display('#########################################')
% display(['All Batchfiles successfully created in ' batch_dir])
% display(' ')
% display('To start the processing:')
%
% % check if we are on a nathan
% if ~ispc
%     [s, m] = system('hostname');
%     if length(m) > 6 && ~isempty(strfind(m, 'lisa'))
%         on_cartesius = 1;
%     end
% end
%
% if exist('on_cartesius', 'var') && on_cartesius
%     % on a Cartesius, short instruction
%     setenv('LD_LIBRARY_PATH'); % make condor_q etc. work
%
%     display('1. Your already on Lisa. Good for you! Makes it even easier!')
%     display(' !! However, check that you REALLY want to start this job, please !!')
%     display('2. Copy the following 2 lines in your matlab shell to start the process:')
%     display(' ')
%     % if the remote machine is not linux, change the following line
%     command{1} = ['!chmod a+x ' batch_dir '*'];
%     display(command{1})
%     command{2} = ['!' batch_dir '/' filename_all];
%     display(command{2})
%     display(' ')
%     r = input('Or type "start" to start the lines now (anything else to skip): ', 's');
%     if strcmp(r, 'start')
%         display(command{1})
%         eval(command{1})
%         display(command{2})
%         eval(command{2})
%     end
% else
%     % not on a nathan, long version
%     display('1. Log into one Nathan (e.g. use Putty)')
%     display('2. In a shell, type the two following lines (e.g. copy & paste):')
%     display(' ')
%     % if the remote machine is not linux, change the following line
%     display(['chmod +777 ' batch_dir '/*'])
%     display([batch_dir '/' filename_all])
% end
%
% display(' ')
% display('3. Wait for messages and relax 8-)')
% display(['  Your log-files will be delivered to ' log_file_dir])
% display(' ')
end
