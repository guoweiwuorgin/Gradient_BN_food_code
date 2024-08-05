%% process pipeline
%% define work dir
clear;clc;
path_dir = 'E:\youyi_fucha\first_test\';
sub_cate = {'CONtrol';'SUB'};
script_dir = 'E:\youyi_fucha\code';
figure_dir = 'E:\youyi_fucha\first_test\figure';
mkdir(figure_dir)
%% Downsample preprocessed images
out_path='E:\youyi_fucha\first_test\gradient_voxel\final_re';
mkdir(out_path);
for group=1:2
    group_all_sub = dir([path_dir sub_cate{group} '\sub-*']);   
    for sub=1:size(group_all_sub)
       subname = group_all_sub(sub).name;
       if ~exist([out_path,'\Resliced_4\CON-',subname],'dir')
           sub_data = [group_all_sub(sub).folder filesep subname];
           if exist([sub_data filesep 'func' filesep subname '_task-rest_space-MNI152NLin2009cAsym_desc-residual_smooth_bold.nii.gz'],'file')
           filename = [sub_data filesep 'func' filesep subname '_task-rest_space-MNI152NLin2009cAsym_desc-residual_smooth_bold.nii.gz'];
           gunzip(filename,[sub_data filesep 'func'])
           end
           filename = [sub_data filesep 'func' filesep subname '_task-rest_space-MNI152NLin2009cAsym_desc-residual_smooth_bold.nii'];
           cd([sub_data filesep 'func'])
           x_reslice([script_dir,'\Reslice_group_mask.nii'],filename,4);
           %mkdir([out_path,'\Resliced\',subname]);
           if group==1
            movefile('r*.nii',[out_path,'\Resliced_4\CON-',subname]);
           else
            movefile('r*.nii',[out_path,'\Resliced_4\PA-',subname]);   
           end
       end
       subname
    end
end
%% Generate connectivity matrix and calculate gradients
cd([out_path,'\Resliced_4']);
list = dir;
list(1:2) = [];
data_dir = 'E:\youyi_fucha\first_test\gradient_voxel\';
for i = 110:length(list)
    tic
    cd([out_path,'\Resliced_4\',list(i).name]);
    filename = dir('*.nii');
    M = x_gen_matrix_voxel([script_dir,'\Reslice_group_mask.nii'],filename.name);
    n = length(M);
    M_spar = M;
    tmp = sort(M);
    tmp = M - repmat(tmp(round(n*0.9),:),n,1);
    M_spar(tmp<0) = 0;
    M_spar = M_spar';
    
    M_cos = 1 - squareform(pdist(M_spar,'cosine'));
    M_normalized = 1 - acos(M_cos)/pi;
    [embedding,result] = x_compute_diffusion_map(M_normalized,0.5,30);
    
    ind_dir = [data_dir,'Gradient_SameLength_4\',list(i).name];
    mkdir(ind_dir);
        filename = [ind_dir,'\gradient.mat'];
        save(filename,'embedding','result');
    toc
end
%% arrange emb
cd([data_dir,'Gradient_SameLength_4']);
list = dir;
list(1:2) = [];
n_sub = length(list);
emb_all = cell(n_sub,1);
res_all = cell(n_sub,1);
for i = 1:n_sub
    cd([data_dir,'Gradient_SameLength_4\',list(i).name]);
    load('gradient.mat');
    emb_all{i} = embedding;
    res_all{i} = result;    
end
% aligned across subjects
[realigned, xfms] = mica_iterativeAlignment(emb_all,100);
realigned = real(realigned);
xfms = cellfun(@real,xfms,'UniformOutput',false);
gradient_emb = cell(30,1);
for i = 1:30
    gradient_emb{i} = squeeze(realigned(:,i,:));
end
%calculate order sequence
seq = zeros(n_sub,30);
for i = 1:n_sub
    tmp = abs(xfms{i});
    [~,I] = sort(tmp,'descend');
    seq_tmp = zeros(5,1);
    seq_tmp(1) = I(1);
    for j = 2:30
        for k = 1:30
            if isempty(find(seq_tmp == I(k,j), 1))
                seq_tmp(j) = I(k,j);
                break;
            end
        end
    end
    seq(i,:) = seq_tmp;
end
%calculate explaination rate
exprate = zeros(n_sub,30);
for i = 1:n_sub
    tmp = res_all{i}.lambdas./sum(res_all{i}.lambdas);
    exprate(i,:) = tmp(seq(i,:));
end


% reorder gradient according to explaination ratio
exprate_mean = mean(exprate);
[~,I] = sort(exprate_mean,'descend');
gradient_emb_reordered = cell(30,1);
for i = 1:30
    gradient_emb_reordered{i} = gradient_emb{I(i)};
end

exprate_reordered = exprate(:,I);
save([out_path '\exprate_reordered.mat'],'exprate_reordered')

% visual check mean gradient
%load('E:\youyi_fucha\exprate_reordered.mat')
hdr_mask = spm_vol([script_dir,'\Reslice_group_mask.nii']);
vol_mask = spm_read_vols(hdr_mask);
hdr = hdr_mask;
hdr.dt(1) = 64;
ind = find(vol_mask);
cd(out_path);
for i = 1:3
    vol = zeros(hdr.dim);
    vol(ind) = mean(gradient_emb_reordered{i},2);
    hdr.fname = ['mean_g_',num2str(i),'.nii'];    
    spm_write_vol(hdr,vol);
end

for i = 1:3
    gradient_emb_reordered{i} = -gradient_emb_reordered{i};
end
save([out_path,'\gradient_emb_correct1_3.mat'],'gradient_emb_reordered');

% generate mean map for both groups
hdr_mask = spm_vol([script_dir,'\Reslice_group_mask.nii']);
vol_mask = spm_read_vols(hdr_mask);
hdr = hdr_mask;
hdr.dt(1) = 64;
ind = find(vol_mask);
cd(out_path);
sub_info = [repmat(1,50,1) ;repmat(2,59,1)];
for i = 1:3
    vol = zeros(hdr.dim);
    vol(ind) = mean(gradient_emb_reordered{i}(:,sub_info(:,1)==2),2);
    hdr.fname = ['mean_g_',num2str(i),'_SUB.nii'];    
    spm_write_vol(hdr,vol);
    vol(ind) = mean(gradient_emb_reordered{i}(:,sub_info(:,1)==1),2);
    hdr.fname = ['mean_g_',num2str(i),'_hc.nii'];   
    spm_write_vol(hdr,vol);
end
mean(exprate_reordered(:,1))
std(exprate_reordered(:,1))
mean(exprate_reordered(sub_info(:,1)==2,1))
std(exprate_reordered(sub_info(:,1)==2,1))
mean(exprate_reordered(sub_info(:,1)==1,1))
std(exprate_reordered(sub_info(:,1)==1,1))

mean(exprate_reordered(:,2))
std(exprate_reordered(:,2))
mean(exprate_reordered(sub_info(:,1)==2,2))
std(exprate_reordered(sub_info(:,1)==2,2))
mean(exprate_reordered(sub_info(:,1)==1,2))
std(exprate_reordered(sub_info(:,1)==1,2))

% draw explaination ratio
close all

close all
cum_exprate_hc = cumsum(mean_exprate_hc);
cum_exprate_mdd = cumsum(mean_exprate_mdd);
ph = plot(1.2:1:30.2,cum_exprate_hc,'.-','color',color_hc,'linewidth',0.5);
hold on
pm = plot(0.8:1:29.8,cum_exprate_mdd,'.-','color',color_mdd,'linewidth',0.5);
hold off
ph.MarkerSize = 10;
pm.MarkerSize = 10;
l = legend;
l.String = {'Controls','SUB'};
l.Location = 'northeast';
xlabel('Component of diffusion embedding');
ylabel('Cumulative explained ratio');
set(gca,'LineWidth',0.5);
set(gca,'FontName','Arial','FontSize',6);
set(gca,'YLim',[0,1],'YTick',0:0.25:1);
set(gca,'XLim',[0,31]);
box off
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 8.9 8]);
print(gcf,[figure_dir,'CompCum.tif'],'-dtiff','-r1000');

% calculate gradient range
emb_range = zeros(n_sub,3);
for i = 1:3
    for j = 1:n_sub
        emb_range(j,i) = max(gradient_emb_reordered{i}(:,j)) - min(gradient_emb_reordered{i}(:,j));
    end
end
save([out_path,'\emb_range.mat'],'emb_range');

% calculate variation
emb_std = zeros(n_sub,3);
for i = 1:3
    for j = 1:n_sub
        emb_std(j,i) = std(gradient_emb_reordered{i}(:,j));
    end
end
save([out_path,'\emb_std.mat'],'emb_std');
