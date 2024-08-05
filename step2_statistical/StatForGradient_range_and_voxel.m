%% process pipeline
%% define work dir
clear;clc;
data_dir = 'D:\youyi_fucha\first_test\gradient_voxel\final_re\';
script_dir = 'D:\youyi_fucha\code\';
figure_dir = 'D:\youyi_fucha\first_test\figure\';
[subdata,b,sub_info] = xlsread('D:\youyi_fucha\first_test\gradient_voxel\Step_01_subject_data\gradient_use.xlsx');
% %% sub-network differences in mean maps
load([data_dir,'gradient_emb_correct1_3.mat']);
for n=2:size(sub_info,1)
    if strcmp(sub_info{n,22},'CON')
        if sub_info{n,3}==1
          sub_infodata(n-1,1) = 1;
        else
          sub_infodata(n-1,1) = 0;
        end  
    else
        sub_infodata(n-1,1) = 2;
    end
end
id_hc = sub_infodata(:,1)==1;
id_sub = sub_infodata(:,1)==2;
for gradient=1:3
    mean_gradient_sub = -mean(gradient_emb_reordered{gradient}(:,id_sub),2)';
    mean_gradient_hc = -mean(gradient_emb_reordered{gradient}(:,id_hc),2)';
    save(['mean_gradient_hc_' num2str(gradient) '.txt'],'mean_gradient_hc','-ascii');
    save(['mean_gradient_sub_' num2str(gradient) '.txt'],'mean_gradient_hc','-ascii');
    close all
    h_pos = histfit(mean_gradient_sub,[],'kernel');
    h_pos(1).FaceColor = [211 211 211]/256;
    h_pos(1).FaceAlpha = 0.6;
    h_pos(2).Color = [211 211 211]/256;
    h_pos(2).LineWidth = 0.5;
    hold on
    h_pos = histfit(mean_gradient_hc,[],'kernel');
    h_pos(1).FaceColor = [238 99 99]/256;
    h_pos(1).FaceAlpha = 0.6;
    h_pos(2).Color = [238 99 99]/256;
    h_pos(2).LineWidth = 0.5;
    box off
    print(gcf,[figure_dir,'Gradient_',num2str(gradient),'.tif'],'-dtiff','-r1000')
end
%%
vol_mod = spm_read_vols(spm_vol([script_dir,'mod_mask.nii']));
vol_mask = spm_read_vols(spm_vol([script_dir,'\Reslice_group_mask.nii']));
ind_mod = vol_mod(vol_mask>0);
for n=1:3
    mean_gradient_sub = mean(gradient_emb_reordered{n}(:,id_sub),2)';
    mean_gradient_hc = mean(gradient_emb_reordered{n}(:,id_hc),2)';
    gradient_score_subnetwork = cell(8,2);
    for i = 1:8
        gradient_score_subnetwork{i,1} = mean_gradient_sub(ind_mod==i);
        gradient_score_subnetwork{i,2} = mean_gradient_hc(ind_mod==i);
    end
    save(['gradient-' num2str(n) '_score_subnetwork.mat'],'gradient_score_subnetwork')
    % test difference
    t = zeros(8,1);
    p = zeros(8,1);
    df = zeros(8,1);
    cohen_d = zeros(8,1);
    for i = 1:8
        [h,p(i),ci,stats] = ttest(gradient_score_subnetwork{i,1},gradient_score_subnetwork{i,2});
        t(i) = stats.tstat;
        df(i) = stats.df;
        cohen_d(i) = mean(gradient_score_subnetwork{i,1}-gradient_score_subnetwork{i,2})/stats.sd;
    end
    [~,~,rnk] = unique(p);
    q_fdr = (p * 8) ./ rnk;
    mod_name = {'VIS','SMN','DAN','VAN','LIM','FPN','DMN','SUB'};
    save([data_dir,'\results\SubNetworkDiffinMeanMap-' num2str(n) '.mat'],'t','p','df','mod_name','q_fdr','cohen_d');
end
%% check skewness and kurtosis
load('exprate_reordered.mat')
load('emb_std.mat')
load('emb_range.mat')
skew_stat = zeros(3,3,2);
for i = 1:3
skew_stat(i,1,1) = skewness(exprate_reordered(id_sub,i));
skew_stat(i,1,2) = skewness(exprate_reordered(id_hc,i));
skew_stat(i,2,1) = skewness(emb_range(id_sub,i));
skew_stat(i,2,2) = skewness(emb_range(id_hc,i));
skew_stat(i,3,1) = skewness(emb_std(id_sub,i));
skew_stat(i,3,2) = skewness(emb_std(id_hc,i));
end

skew_voxel = zeros(3,18933,2);
for i = 1:3
    for j = 1:18933
        skew_voxel(i,j,1) = skewness(gradient_emb_reordered{i}(j,id_sub));
        skew_voxel(i,j,2) = skewness(gradient_emb_reordered{i}(j,id_hc));
    end
end

kurt_stat = zeros(3,3,2);
for i = 1:3
kurt_stat(i,1,1) = kurtosis(exprate_reordered(id_sub,i));
kurt_stat(i,1,2) = kurtosis(exprate_reordered(id_hc,i));
kurt_stat(i,2,1) = kurtosis(emb_range(id_sub,i));
kurt_stat(i,2,2) = kurtosis(emb_range(id_hc,i));
kurt_stat(i,3,1) = kurtosis(emb_std(id_sub,i));
kurt_stat(i,3,2) = kurtosis(emb_std(id_hc,i));
end

kurt_voxel = zeros(3,18933,2);
for i = 1:3
    for j = 1:18933
        kurt_voxel(i,j,1) = kurtosis(gradient_emb_reordered{i}(j,id_sub));
        kurt_voxel(i,j,2) = kurtosis(gradient_emb_reordered{i}(j,id_hc));
    end
end

%% check variance
var_stat = zeros(3,3);
var_p = zeros(3,3);
for i = 1:3
[~,var_p(i,1),~,stat] = vartest2(exprate_reordered(id_sub,i),exprate_reordered(id_hc,i));
var_stat(i,1) = stat.fstat;
[~,var_p(i,2),~,stat] = vartest2(emb_range(id_sub,i),emb_range(id_hc,i));
var_stat(i,2) = stat.fstat;
[~,var_p(i,3),~,stat] = vartest2(emb_std(id_sub,i),emb_std(id_hc,i));
var_stat(i,3) = stat.fstat;
end
%% Between-group difference for global topology
% explained ratio
stat = zeros(3,7);
group = id_hc;
group(group(1:50)==0)=[];
des = [group subdata(logical(id_hc+id_sub),4:5) subdata(logical(id_hc+id_sub),25)];%
n_hc = length(find(id_hc));
n_sub = length(find(id_sub));
sub_index = logical(id_hc+id_sub);
for i = 1:3
    stat(i,1) = mean(exprate_reordered(id_sub,i));
    stat(i,2) = std(exprate_reordered(id_sub,i));
    stat(i,3) = mean(exprate_reordered(id_hc,i));
    stat(i,4) = std(exprate_reordered(id_hc,i));    
    stat_result = regstats(exprate_reordered(sub_index,i),des,'linear',{'tstat','r'});
    stat(i,5) = stat_result.tstat.t(2);
    stat(i,6) = stat_result.tstat.t(2) * sqrt(1/n_hc + 1/n_sub);
    stat(i,7) = stat_result.tstat.pval(2);
    exprate_reordered_all(:,i) = exprate_reordered(sub_index,i);
end
disp(stat)
% gradient range
load([data_dir,'emb_range.mat']);
stat = zeros(3,7);
for i = 1:3
    stat(i,1) = mean(emb_range(id_sub,i));
    stat(i,2) = std(emb_range(id_sub,i));
    stat(i,3) = mean(emb_range(id_hc,i));
    stat(i,4) = std(emb_range(id_hc,i));    
    stat_result = regstats(emb_range(sub_index,i),des,'linear',{'tstat','r'});
    stat(i,5) = stat_result.tstat.t(2);
    stat(i,6) = stat_result.tstat.t(2) * sqrt(1/n_hc + 1/n_sub);
    stat(i,7) = stat_result.tstat.pval(2);
    emb_range_all(:,i) = emb_range(sub_index,i);
end
disp(stat)
% gradient variance
load([data_dir,'emb_std.mat']);
stat = zeros(3,7);
for i = 1:3
    stat(i,1) = mean(emb_std(id_sub,i));
    stat(i,2) = std(emb_std(id_sub,i));
    stat(i,3) = mean(emb_std(id_hc,i));
    stat(i,4) = std(emb_std(id_hc,i));    
    stat_result = regstats(emb_std(sub_index,i),des,'linear',{'tstat','r'});
    stat(i,5) = stat_result.tstat.t(2);
    stat(i,6) = stat_result.tstat.t(2) * sqrt(1/n_hc + 1/n_sub);
    stat(i,7) = stat_result.tstat.pval(2);
end
disp(stat)
%% compute ECC
clear gradient_all
gradient_all.gradient1 = reshape(gradient_emb_reordered{1},18933*109,1);
gradient_all.gradient2 = reshape(gradient_emb_reordered{2},18933*109,1);
gradient_all.gradient3 = reshape(gradient_emb_reordered{3},18933*109,1);
gradient_all.gradient4 = reshape(gradient_emb_reordered{4},18933*109,1);
gradient_all.sub = {};
for n=1:109
    subname_tmp = repmat({['subname-' num2str(n)]},18933,1);
    gradient_all.sub = [gradient_all.sub ;subname_tmp;];
end
gradient_all = struct2table(gradient_all);
Gradient_ECC = compute_eccentricity(gradient_all);
Gradient_ECC = reshape(Gradient_ECC.distance,18933,109);
%Gradient_ECC  = reshape(Gradient_ECC,18933,109);
T_value = zeros(1,n_voxel);
sub_index = logical(id_hc+id_sub);
res_map_ECC = zeros(n_voxel,n_sub);
for j = 1:n_voxel
    stat_result = regstats(Gradient_ECC(j,sub_index)',des,'linear',{'tstat','r'});
    T_value(1,j) = stat_result.tstat.t(2);
    res_map_ECC(j,:) = stat_result.r;
end
dfe = stat_result.tstat.dfe;
Z_value = spm_t2z(T_value,dfe);
% estimate GRF correction, save cluster nii, and draw pic
hdr_mask = spm_vol([script_dir,'Reslice_group_mask.nii']);
vol_mask = spm_read_vols(hdr_mask);
ind = find(vol_mask);
vox = [4 4 4];
voxel_p = 0.001;
cluster_p = 0.05;
tail = 1;
zthrd = norminv(1 - voxel_p/tail);
surfacefile = 'D:\youyi_fucha\code\MDD_ConnectomeGradient-main\MDD_ConnectomeGradient-main\0_GradientCalculation\FSaverage_midthickness_32K.nv';
mkdir([data_dir,'BetweenGroupDiff']);
cd([data_dir,'BetweenGroupDiff']);

sig_index=1;
R_volume = zeros([hdr_mask.dim,n_sub]);
for j = 1:n_sub
    tmp_vol = zeros(hdr_mask.dim);
    tmp_vol(ind) = res_map_ECC(:,j);
    R_volume(:,:,:,j) = tmp_vol;
end

[cluster_size,dlh,fwhm] = x_GRF(R_volume,dfe,vol_mask,vox,voxel_p,cluster_p,tail);
cluster_size
Z_vol = zeros(hdr_mask.dim);
Z_vol(ind) = Z_value(1,:);
hdr_z = hdr_mask;
hdr_z.fname = ['g',num2str(1),'_T2_z_ECC.nii'];
hdr_z.descrip = sprintf('SPM{Z_[1]}{dLh_%f}{FWHMx_%fFWHMy_%fFWHMz_%fmm}',...
    dlh, fwhm(1), fwhm(2), fwhm(3));
hdr_z.dt(1) = 64;
spm_write_vol(hdr_z,Z_vol);

Z_vol(Z_vol < zthrd & Z_vol > -zthrd) = 0;
[L, num] = bwlabeln(Z_vol,26);
n = 0;
for x = 1:num
    theCurrentCluster = L == x;
%         length(find(theCurrentCluster))
    if length(find(theCurrentCluster)) <= cluster_size
        n = n + 1;
        Z_vol(logical(theCurrentCluster)) = 0;
    end
    length(find(theCurrentCluster))
    if length(find(theCurrentCluster)) >= cluster_size
        length(find(theCurrentCluster))
        signifcant_voxel{1,sig_index} = find(theCurrentCluster==1);
        sig_index=sig_index+1;
    end
end
hdr_z.fname = ['g1_T2_z_cluster_ECC.nii'];
spm_write_vol(hdr_z,Z_vol);
niifile = hdr_z.fname;
txtfile = ['g1_T2_z_cluster_ECC_T2_z_cluster.txt'];

hdr_z.fname = ['ECC_T2_z_cluster_mask.nii'];
Z_vol(Z_vol~=0) = 1;
spm_write_vol(hdr_z,Z_vol);

%% Reginal difference
% estimate group difference
[n_voxel,n_sub] = size(gradient_emb_reordered{1});
T_value = zeros(3,n_voxel);
res_map = cell(3,1);
sub_index = logical(id_hc+id_sub);
n_sub = length(group);
for i = 1:3
    res_map{i} = zeros(n_voxel,n_sub);
    for j = 1:n_voxel
        stat_result = regstats(gradient_emb_reordered{i}(j,sub_index)',des,'linear',{'tstat','r'});
        T_value(i,j) = stat_result.tstat.t(2);
        res_map{i}(j,:) = stat_result.r;
    end
end
dfe = stat_result.tstat.dfe;
Z_value = spm_t2z(T_value,dfe);
% estimate GRF correction, save cluster nii, and draw pic
hdr_mask = spm_vol([script_dir,'Reslice_group_mask.nii']);
vol_mask = spm_read_vols(hdr_mask);
ind = find(vol_mask);
vox = [4 4 4];
voxel_p = 0.001;
cluster_p = 0.05;
tail = 1;
zthrd = norminv(1 - voxel_p/tail);
surfacefile = 'D:\youyi_fucha\code\MDD_ConnectomeGradient-main\MDD_ConnectomeGradient-main\0_GradientCalculation\FSaverage_midthickness_32K.nv';
mkdir([data_dir,'BetweenGroupDiff']);
cd([data_dir,'BetweenGroupDiff']);
for i = 1:3
    sig_index=1;
    R_volume = zeros([hdr_mask.dim,n_sub]);
    for j = 1:n_sub
        tmp_vol = zeros(hdr_mask.dim);
        tmp_vol(ind) = res_map{i}(:,j);
        R_volume(:,:,:,j) = tmp_vol;
    end
    [cluster_size,dlh,fwhm] = x_GRF(R_volume,dfe,vol_mask,vox,voxel_p,cluster_p,tail);
    cluster_size
    Z_vol = zeros(hdr_mask.dim);
    Z_vol(ind) = Z_value(i,:);
    hdr_z = hdr_mask;
    hdr_z.fname = ['g',num2str(i),'_T2_z.nii'];
    hdr_z.descrip = sprintf('SPM{Z_[1]}{dLh_%f}{FWHMx_%fFWHMy_%fFWHMz_%fmm}',...
        dlh, fwhm(1), fwhm(2), fwhm(3));
    hdr_z.dt(1) = 64;
    spm_write_vol(hdr_z,Z_vol);
    
    Z_vol(Z_vol < zthrd & Z_vol > -zthrd) = 0;
    [L, num] = bwlabeln(Z_vol,26);
    n = 0;
    for x = 1:num
        theCurrentCluster = L == x;
%         length(find(theCurrentCluster))
        if length(find(theCurrentCluster)) <= cluster_size
            n = n + 1;
            Z_vol(logical(theCurrentCluster)) = 0;
        end
        if length(find(theCurrentCluster)) >= cluster_size
            length(find(theCurrentCluster))
            signifcant_voxel{i,sig_index} = find(theCurrentCluster==1);
            sig_index=sig_index+1;
        end
    end
    hdr_z.fname = ['g',num2str(i),'_T2_z_cluster.nii'];
    spm_write_vol(hdr_z,Z_vol);
    niifile = hdr_z.fname;
    txtfile = ['g',num2str(i),'_T2_z_cluster.txt'];
    
    hdr_z.fname = ['g',num2str(i),'_T2_z_cluster_mask.nii'];
    Z_vol(Z_vol~=0) = 1;
    spm_write_vol(hdr_z,Z_vol);
%     BrainNet_nii2txt(surfacefile,niifile,txtfile);
%     BrainNet_MapCfg('D:\Dropbox\Projects\BrainNet\otherfile\Surface\FSaverage_inflated_32K.nv',...
%         txtfile,[script_dir,'BNVCfg_T2_cluster.mat'],['g',num2str(i),'_T2_z_cluster.tif']);
end


R_volume = zeros([hdr_mask.dim,n_sub]);
for n=1:3
    mkdir(['g',num2str(n),'_all_gradient']);
    cd(['g',num2str(n),'_all_gradient']);
    for j = 1:n_sub
        hdr_z = hdr_mask;
        tmp_vol = zeros(hdr_mask.dim);
        tmp_vol(ind) = gradient_emb_reordered{n}(:,j);
        hdr_z.fname = ['g',num2str(n),'_data-',num2str(j),'_all_gradient.nii'];
        hdr_z.descrip = sprintf('SPM{Z_[1]}{dLh_%f}{FWHMx_%fFWHMy_%fFWHMz_%fmm}',...
        dlh, fwhm(1), fwhm(2), fwhm(3));
        hdr_z.dt(1) = 64;
        spm_write_vol(hdr_z,tmp_vol);
    end
    cd('../');
end


size(R_volume)
all_data = reshape(R_volume,37*45*37,85);


% 
% gradient_2_hipp=mean(all_data(signifcant_voxel,:),1)';
% 
% gradient_2_corr1=mean(all_data(signifcant_voxel{1},:),1)';
% gradient_2_corr2=mean(all_data(signifcant_voxel{2},:),1)';
% gradient_2_corr3=mean(all_data(signifcant_voxel{3},:),1)';
% gradient_2_corr4=mean(all_data(signifcant_voxel{4},:),1)';


for i = 1:3
    txt = Z_value(i,:);
    save(['g',num2str(i),'_T2_z.txt'],'txt','-ascii');
end

%% generate surrogate maps for between-group z maps



%% between-group difference among moduele

% draw figure for modular distribution of z
mod_color = [118 113 113;162 81 172;120 154 192;64 152 50;224 101 254;221.4000  227.7000  180.9000;239 185 67;217 113 125;]/255;
mod_label = {'SUB','VIS','SMN','DAN','VAN','LIB','FPN','DMN'};
for i = 1:8
    close all
    h_pos = histfit(gradient_score_subnetwork{i,1},[],'kernel');
    h_pos(2).Color = mod_color(i,:);
    h_pos(2).LineWidth = 0.5;
    hold on
    yFill = [h_pos(2).YData,zeros(1,length(h_pos(2).YData))];
    xFill = [h_pos(2).XData,fliplr(h_pos(2).XData)];
    fill(xFill,yFill,mod_color(i,:),'FaceAlpha',0.5,'EdgeColor',mod_color(i,:));
    
    h_pos = histfit(gradient_score_subnetwork{i,2},[],'kernel');
    
    
    
   
    print(gcf,[figure_dir,'neg_t2_dist_',mod_label{i},'.tif'],'-dtiff','-r1000')
end
