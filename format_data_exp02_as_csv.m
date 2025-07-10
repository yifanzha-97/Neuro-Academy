% code that runs data curation
% author: steeve laquitaine
% 
% usage:
%
% setup packages
% addpath(genpath('~/Desktop/projects/project/mgl'));
    
% setup parameters
exp_i = 1;
n_subjects = 5;
write_path = '/Users/zhayf/Matlab/Project/';
exp_name = 'data02_direction1prior';

% run curation
format_files_by_experiment(exp_i, exp_name, write_path, n_subjects);


function format_files_by_experiment(exp_i, exp_name, write_path, n_subjects)
    
    % initialize
    data_by_subject = cell(n_subjects);

    % collect all subject datasets
    for subject_i = 1 : n_subjects
        
        % log
        fprintf(['(format_files_by_experiment) processing subject' num2str(subject_i) '\n'])

        % start timer
        tic

        % get this subject's dataset
        [data, headers] = format_files_by_subject(exp_i, subject_i, write_path);

        % report time
        toc 

        % tape datasets
        data_by_subject{subject_i} = data;
    end

    % stack subject datasets
    data_for_experiment = vertcat(data_by_subject{:});

    % format for writing
    dataset_with_header = cell2table(data_for_experiment,'VariableNames',headers);

    % create path
    mkdir(write_path);
    
    % write dataset
    writetable(dataset_with_header, [write_path exp_name '.csv']);
end

function [data_for_subject, headers] = format_files_by_subject(exp_i, subject_i, write_path)
    
    % setup project path
    proj_path = '/Users/zhayf/Matlab/Project/';

    % setup data path
    data_path = [proj_path];

    % Extract the metadata from file name
    % get the experiments
    experiment_set = {dir([data_path 'data0*']).name};
    exp = experiment_set{exp_i};

    % get the subjects
    subject_set = {dir([data_path exp '/data/sub*']).name};
    subject = str2num(subject_set{subject_i}(4:5));

    % get the .mat file to parse
    subject_path = [data_path exp '/data/' subject_set{subject_i}];
    mat_files = {dir([subject_path '/steeve*.mat']).name};

    % init dataset
    data_by_file = cell(length(mat_files));

    for file_i = 1 : length(mat_files)

        % curate file
        this_file = mat_files{file_i};
        [data, headers] = format_mat_file(exp_i, subject_i, this_file, file_i, data_path, write_path);
        
        % tape formatted file
        data_by_file{file_i} = data;
    end

    % concatenate subject's files
    data_for_subject = vertcat(data_by_file{:});
end

function [dataset, headers] = format_mat_file(exp_i, subject_i, mat_file, file_i, data_path, write_path)

    % % get the experiments
    experiment_set = {dir([data_path 'data0*']).name};
    exp = experiment_set{exp_i};

    % % get the subjects
    subject_set = {dir([data_path exp '/data/sub*']).name};
    subject = str2num(subject_set{subject_i}(4:5));

    % parse the experiment_id (ONLY FOR EXP_I = 1)
    % exp_id = str2num(mat_file(11:12));
    % exp_id_str = mat_file(8:12);

    % parse the session id
    session_id = str2num(mat_file(23:24));
    session_id_str = mat_file(19:24);

    % parse the run id
    try
        run_id = str2num(mat_file(29:30));
    catch except
        keyboard
    end
    run_id_str = mat_file(26:30);

    % parse the prior std condition
    prior_std = str2num(mat_file(36:38));
    prior_std_str = mat_file(32:38);

    % parse the prior mean condition
    prior_mean = str2num(mat_file(44:46));
    prior_mean_str = mat_file(40:46);

    % parse the run date
    % run_date = mat_file(92:97);
    other = mat_file(47:length(mat_file));

    % reconstruct file name
    file_name = [data_path ...
                experiment_set{exp_i} ...
                '/data/' ...
                subject_set{subject_i} ...
                '/steeve' ...
                '_data_' ...
                subject_set{subject_i} ...
                '_' session_id_str ...
                '_' run_id_str ...
                '_' prior_std_str ...
                '_' prior_mean_str ...
                other];

    % get run data
    data = getTaskParameters(file_name);

    % format task variables
    % ---------------------
    n_trials = length(data{2}.randVars.prodcoor);

    % set trial indices
    trial_index = linspace(1, n_trials, n_trials)';

    % format motion direction task variable
    motion_direction = nan(n_trials,1);
    if isfield(data{2}.randVars,'myRandomDir')
        motion_direction = data{2}.randVars.myRandomDir';
    end

    % format motion coherence task variable
    motion_coherence = nan(n_trials,1);
    if isfield(data{2}.randVars,'myRandomCoh')
        motion_coherence = data{2}.randVars.myRandomCoh';
    end

    % format initial response arrow task variable
    response_arrow_start_angle = nan(n_trials,1);
    if isfield(data{2}.randVars,'initAngledeg')
        response_arrow_start_angle = data{2}.randVars.initAngledeg';
    end

    % format trial time task variable
    trial_time = nan(n_trials,1);
    if isfield(data{2},'trialTime')
        trial_time = data{2}.trialTime';
    end


    % format behavioral data
    % ----------------------
    % format estimate data
    estimate = nan(n_trials,2);
    if isfield(data{2}.randVars,'prodcoor')
        for ix = 1:n_trials
            estimate(ix, :) = data{2}.randVars.prodcoor{ix};
        end
    end

    % format reaction time
    reaction_time = nan(n_trials,1);
    if isfield(data{2},'reactionTime')
        reaction_time = data{2}.reactionTime';
    end

    % format raw response time
    raw_response_time = nan(n_trials,1);
    if isfield(data{2},'responseTimeRaw')
        raw_response_time = data{2}.responseTimeRaw';
    end


    % write dataset
    % -------------
    dataset = [
        num2cell(trial_index) ...
        num2cell(trial_time) ...
        num2cell(response_arrow_start_angle) ...
        num2cell(motion_direction) ...
        num2cell(motion_coherence) ...    
        num2cell(estimate) ...
        num2cell(reaction_time) ...
        num2cell(raw_response_time) ...
        num2cell(repmat(prior_std, n_trials,1)) ...  
        num2cell(repmat(prior_mean, n_trials,1)) ...  
        num2cell(repmat(subject, n_trials,1)) ...    
        cellstr(repmat(experiment_set{exp_i}, n_trials,1)) ...
        num2cell(repmat(session_id, n_trials,1)) ...
        num2cell(repmat(run_id, n_trials,1)) ...
        ];
    headers = {
        'trial_index' ...
        'trial_time' ...
        'response_arrow_start_angle' ...
        'motion_direction' ...
        'motion_coherence' ...    
        'estimate_x' ...
        'estimate_y' ...
        'reaction_time' ...
        'raw_response_time' ...
        'prior_std' ...
        'prior_mean' ...
        'subject_id' ...
        'experiment_name' ...
        'session_id' ...
        'run_id' ...
        };
end


