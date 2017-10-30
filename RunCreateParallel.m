%CreateParallel Bash scripts for ori discrimination o nthe S03
%% Preparation
% parallel_fun = 'parallelFun_s03_part1.m'
% job_list1 = [1,2,3];
% job_list2 = [1,2,3];
% job_list3 = [1];
% parallel_fun_dir = '~/MVPA_Scripts/';
% job_name = 'OriDecoderPrep';
% 
% create_parallel_S03(parallel_fun, job_list1, job_list2, job_list3, parallel_fun_dir, job_name)

%% Decoder
% parallel_fun = 'parallelFun_s03_part2.m'
% job_list1 = [1,2,3];
% job_list2 = [1,2,3];
% job_list3 = [1:36];
% parallel_fun_dir = '~/MVPA_Scripts/';
% job_name = 'OriDecoder';
% 
% disp('Running create_parallel_S03')
% create_parallel_S03(parallel_fun, job_list1, job_list2, job_list3, parallel_fun_dir, job_name)

% %% ALL STIMULUS
parallel_fun = 'parallelFun_SARA'
job_list1 = [1,2,3,4]; %MOUSE
job_list2 = [1,2,3,4]; %TW
job_list3 = [1,2,3]; %ReactionOpt
parallel_fun_dir = '"$TMPDIR"/MVPA_Scripts/';
job_name = 'StimDecoder';
% for i_one = job_list1
%     for i_two = job_list2
%         FG_MVPA_SVM_RFE(i_one,i_two)
%     end
% end

disp('Running create_parallel_SARA')
create_parallel_SARA(parallel_fun, job_list1, job_list2, job_list3, parallel_fun_dir, job_name)

%% ALL REACTION
% parallel_fun = 'parallelFun_SARA_Reaction'
% job_list1 = [1,2,3,4]; %MOUSE
% job_list2 = [1,2,3,4]; %TW
% % job_list3 = [1,2,3]; %ReactionOpt
% parallel_fun_dir = '"$TMPDIR"/MVPA_Scripts/';
% job_name = 'ReactionDecoder';
% % for i_one = job_list1
% %     for i_two = job_list2
% %         FG_MVPA_SVM_RFE(i_one,i_two)
% %     end
% % end
% 
% disp('Running create_parallel_SARA')
% create_parallel_SARA_Reaction(parallel_fun, job_list1, job_list2, parallel_fun_dir, job_name)
