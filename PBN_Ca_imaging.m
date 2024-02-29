
clear all;  % clear all variables
close all;  % close all open graphs
clc   % clear command window

%% data metadata
%%%%FILLIN with the path to where the data is stored
%%%%FILLIN name the file
file_name_1 = ''; %% FILLIN data file (exported from IDPS) pathway
t_point_dff_filename = ''; %% FILLIN output file pathway
zscore_filename = ''; %% %% FILLIN output file pathway
dff_zscore_combined_filename = ''; %% FILLIN output file pathway
single_cell_filename = ''; %% FILLIN output file pathway
%%FILLIN time point in sec
t_point = [];

%% NO need change
%%%% Sorting data 
%%%% Generate variable files with same stimulus
All_cell = readtable(file_name_1, 'ReadVariableNames', true, 'ReadRowNames', true);
status = readcell(file_name_1,'Range','2:2'); %%%% excel use '3:3', csv use '2:2'
cell_num = width(All_cell);
status = status(1,2:end);
T = [];

for i = 1:cell_num
    %% CHOICE 1: accepted cell
     if strcmp(status(1,i),'accepted') 
         T = [T,All_cell(2:end,i)];
     %% CHOICE 2: undecided cell
%         if strcmp(status(1,i),'undecided')  
%           T = [T,All_cell(2:end,i)];
    end
end

% varnames = T.Properties.VariableNames;
rownames = T.Properties.RowNames;
varnames = T.Properties.VariableNames;
cell_num = width(T);

%%%%
%%%% Moving average
%        T = table2array(T);
%        k = 3; %%%%FILLIN windlow length
%        T = movmean(T,k,1);
%        T = array2table(T);
stim_data ={};
t_OL = round(str2double(rownames'),1);
%%%% FILLIN time before targeted time point in INDEX
Pre_t = ; %%%FILIIN in index
Post_t = ; %%%FILLIN in index
sheet_name = [];
num_t = length(t_point);
T_row_mean = [];
end_index = height(T);

for i = 1:num_t
    t_point_n = round(t_point(i),1);
    t = find (t_OL == t_point_n,1);
    if ~isempty(t)
        if t - Pre_t >=0
             Pre = t - Pre_t;
        else
            Pre = 0;
        end
            
        if t + Post_t <= end_index
            Post = t + Post_t;
        else
            Post = end_index;
        end
    
        idx = [Pre:1:Post];
        rownames_t =T(idx,:).Properties.RowNames;
        rownames_t{1,:}=num2str(t);
       
        T_row = array2table(T{idx,:}, 'VariableNames', varnames, 'RowNames', rownames_t);
        row_mean = mean(table2array(T_row),2);
        T_row = addvars(T_row,row_mean,'NewVariableNames', 'Mean');
        stim_data{end+1} = T_row;
    else ('not found' + string(t_point_n))
    end

end

%%%%%Generating 10s before Stim Baseline
baseline ={}; 
for n = 1:length(stim_data)
    baseline_i = stim_data{1,n}(1:Pre_t,:);
    baseline_i = table2array(baseline_i);
    baseline{end+1}= baseline_i;
end
baseline_c ={};

for i = 1 : cell_num
    baseline_c_i=[];
    for n = 1:length(stim_data)
         baseline_c_ii = baseline{1,n}(:,i);
         baseline_c_i  = [baseline_c_i,baseline_c_ii];
    end
    baseline_c{end+1} = mean(baseline_c_i,2);
end

%%dF/F for each type of Stim
T_n_mean_all = [];
T_Z_mean_all = [];
T_n_all ={};
T_Z_all = {};
T_og = T;
for n = 1 : length(stim_data)
    cell_num = width(stim_data{1,n});
    T = stim_data{1,n};
    T_n =[];
    T_Z =[];
    %%%% Calculating Normalization
    Baseline_pre = 1; %% Customizable
    Baseline_post = 90; %% Customizable
    for i = 1:cell_num-1
%            norm = (T{:,i} - mean(T_og{100:11000,i})) / mean(T_og{100:11000,i}); %% CHOICE 1: Whole Trace based baseline for manual ROI data
           norm = (T{:,i}); %% CHOICE 2 : RAW extract for PCA or CNMFe data
        T_n = [T_n , norm]; 
    end

    %%%% Calculating Z score
    for i = 1:cell_num         
        z = (T{:,i} - mean(T{:,i})) / std(T{:,i});
        T_Z = [T_Z , z]; 
    end

    T_n_mean = mean(T_n,2);
    T_n = [T_n,T_n_mean];
    T_Z_mean =mean(T_Z,2);
    T_Z = [T_Z,T_Z_mean];
    T_n_mean_all = [T_n_mean_all,T_n_mean];%% command out with inconsistent length
    T_Z_mean_all =[T_Z_mean_all,T_Z_mean];%% command out with inconsistent length
    varnames_n = T.Properties.VariableNames;
    varnames_z = T.Properties.VariableNames;
    rownames_n = T.Properties.RowNames;
%     varnames_n{end +1} ='F_mean';
    varnames_z{end +1} ='z_mean';
    T_n_all{end+1} = T_n;
    T_Z_all{end+1} = T_Z;
    sheet_name = char(stim_data{1,n}.Properties.RowNames(1));
    T_n = array2table(T_n, 'VariableNames', varnames_n, 'RowNames', rownames_n);
    writetable(T_n, t_point_dff_filename, 'writeVariableNames',1, 'writeRowNames', 1, "Sheet",sheet_name, "WriteMode", "inplace");
    T_Z = array2table(T_Z, 'VariableNames', varnames_z, 'RowNames', rownames_n);
    writetable(T_Z, zscore_filename, 'writeVariableNames',1, 'writeRowNames', 1, "Sheet",sheet_name, "WriteMode", "inplace");
end


% %% command out with inconsistent length
T_n_mean_all = array2table(T_n_mean_all);
writetable(T_n_mean_all, dff_zscore_combined_filename, 'writeVariableNames',1, 'writeRowNames', 1, "Sheet",  'T_dff_mean', "WriteMode", "inplace");
T_Z_mean_all = array2table(T_Z_mean_all);
writetable(T_Z_mean_all, dff_zscore_combined_filename, 'writeVariableNames',1, 'writeRowNames', 1, "Sheet",  'T_Z_mean', "WriteMode", "inplace");

sc_n = [];
sc_Z = [];
for i = 1 : cell_num
    sc_n_i=[];
    sc_Z_i=[];
    for n = 1:length(stim_data)
         sc_n_ii = T_n_all{1,n}(:,i);
         sc_n_i  = [sc_n_i,sc_n_ii];
         sc_Z_ii = T_Z_all{1,n}(:,i);
         sc_Z_i  = [sc_Z_i,sc_Z_ii];
    end
    sc_n = [sc_n,mean(sc_n_i,2)];
    sc_Z = [sc_Z,mean(sc_Z_i,2)];
end 
varnames{end+1} = 'mean';
sc_n_all = array2table(sc_n, "VariableNames",varnames);
writetable(sc_n_all, single_cell_filename, 'writeVariableNames',1,  "Sheet",  'single_cell_dff', "WriteMode", "inplace");
sc_Z_all = array2table(sc_Z, "VariableNames",varnames);
writetable(sc_n_all, single_cell_filename, 'writeVariableNames',1, "Sheet",  'single_cell_Z', "WriteMode", "inplace");

