% This code generates the figures for the AD dCA1-vCA1 paper
%% Split up animals and ROIs, initialize variables, load histology and age data
cd('Y:\Udaysankar Chockanathan\Aim 1\kilosort_results_all'); 
load('kilosort_preload_final.mat'); load('ks_all_data5.mat'); load('VR_reformat.mat'); load('histology_results.mat'); 
dCA1_sites = [1,3,4,6,8,11,14,16,18,20,23,25]; %rows of ks_all_data that have dCA1 recordings
vCA1_sites = [2,5,7,9,10,12,13,15,17,19,21,22,24,26,27]; %rows of ks_all_data that have vCA1 recordings
WT_sites = 16:27;
AD_sites = 1:15;
WT_dCA1_sites = intersect(WT_sites,dCA1_sites);
WT_vCA1_sites = intersect(WT_sites,vCA1_sites); 
AD_dCA1_sites = intersect(AD_sites,dCA1_sites);
AD_vCA1_sites = intersect(AD_sites,vCA1_sites); 
VR_ephys_dCA1_sites = []; VR_ephys_vCA1_sites = []; 
mod_rasters = cell(size(ks_all_data,1),1); 
mod_rasters_nds = mod_rasters;
mod_durat = zeros(size(ks_all_data,1),1); 
mod_durat_nds = mod_durat;
mod_vel = ks_all_data(:,8); 
for an = 1:size(ks_all_data,1)
    mod_rasters{an} = ks_all_data{an,1}(:,3); 
    mod_rasters_nds{an} = ks_all_data{an,1}(:,2); 
    mod_durat(an) = ks_all_data{an,12}; 
    mod_durat_nds(an) = ks_all_data{an,11}; 
end
% Animal ages in days
ks_all_data{1,24} = 463; ks_all_data{2,24} = 463; ks_all_data{3,24} = 463; 
ks_all_data{4,24} = 470; ks_all_data{5,24} = 470; ks_all_data{6,24} = 470; ks_all_data{7,24} = 470; 
ks_all_data{8,24} = 526; ks_all_data{9,24} = 526; 
ks_all_data{10,24} = 561; ks_all_data{11,24} = 561; ks_all_data{12,24} = 561; 
ks_all_data{13,24} = 599; ks_all_data{14,24} = 599; 
ks_all_data{15,24} = 599; 
ks_all_data{16,24} = 456; ks_all_data{17,24} = 456; ks_all_data{18,24} = 456; ks_all_data{19,24} = 456; 
ks_all_data{20,24} = 463; ks_all_data{21,24} = 463;
ks_all_data{22,24} = 558; ks_all_data{23,24} = 558; ks_all_data{24,24} = 558;
ks_all_data{25,24} = 586; ks_all_data{26,24} = 586; ks_all_data{27,24} = 586;
% Format histology data
VR_histo = zeros(size(VR_reformatting,1),5,2);
VR_histo(1,[1,3,5],1) = [dCA1_fraction(1),dHpc_fraction(1),all_fraction(1)]; 
VR_histo(2,[2,4,5],1) = [vCA1_fraction(1),vHpc_fraction(1),all_fraction(1)]; 
VR_histo(3,[1,3,5],1) = [dCA1_fraction(1),dHpc_fraction(1),all_fraction(1)]; 
VR_histo(4,[1,3,5],1) = [dCA1_fraction(2),dHpc_fraction(2),all_fraction(2)];
VR_histo(5,[2,4,5],1) = [vCA1_fraction(2),vHpc_fraction(2),all_fraction(2)];
VR_histo(6,[1,3,5],1) = [dCA1_fraction(2),dHpc_fraction(2),all_fraction(2)]; 
VR_histo(7,[2,4,5],1) = [vCA1_fraction(2),vHpc_fraction(2),all_fraction(2)];
VR_histo(8,[1,3,5],1) = [dCA1_fraction(3),dHpc_fraction(3),all_fraction(3)]; 
VR_histo(9,[2,4,5],1) = [vCA1_fraction(3),vHpc_fraction(3),all_fraction(3)];
VR_histo(10,[2,4,5],1) = [vCA1_fraction(4),vHpc_fraction(4),all_fraction(4)];
VR_histo(11,[1,3,5],1) = [dCA1_fraction(4),dHpc_fraction(4),all_fraction(4)]; 
VR_histo(12,[2,4,5],1) = [vCA1_fraction(4),vHpc_fraction(4),all_fraction(4)];
VR_histo(13,[2,4,5],1) = [vCA1_fraction(5),vHpc_fraction(5),all_fraction(5)];
VR_histo(14,[1,3,5],1) = [dCA1_fraction(5),dHpc_fraction(5),all_fraction(5)]; 
VR_histo(15,[2,4,5],1) = [vCA1_fraction(6),vHpc_fraction(6),all_fraction(6)];
VR_histo(1,[1,3,5],2) = [dCA1_density(1),dHpc_density(1),all_density(1)]*(1.6^2)*(1000)^2; % 1um = 1.6 pixels; *(1000)^2 to convert to plaques per mm^2
VR_histo(2,[2,4,5],2) = [vCA1_density(1),vHpc_density(1),all_density(1)]*(1.6^2)*(1000)^2;  
VR_histo(3,[1,3,5],2) = [dCA1_density(1),dHpc_density(1),all_density(1)]*(1.6^2)*(1000)^2;   
VR_histo(4,[1,3,5],2) = [dCA1_density(2),dHpc_density(2),all_density(2)]*(1.6^2)*(1000)^2;  
VR_histo(5,[2,4,5],2) = [vCA1_density(2),vHpc_density(2),all_density(2)]*(1.6^2)*(1000)^2;  
VR_histo(6,[1,3,5],2) = [dCA1_density(2),dHpc_density(2),all_density(2)]*(1.6^2)*(1000)^2;   
VR_histo(7,[2,4,5],2) = [vCA1_density(2),vHpc_density(2),all_density(2)]*(1.6^2)*(1000)^2;  
VR_histo(8,[1,3,5],2) = [dCA1_density(3),dHpc_density(3),all_density(3)]*(1.6^2)*(1000)^2;   
VR_histo(9,[2,4,5],2) = [vCA1_density(3),vHpc_density(3),all_density(3)]*(1.6^2)*(1000)^2;  
VR_histo(10,[2,4,5],2) = [vCA1_density(4),vHpc_density(4),all_density(4)]*(1.6^2)*(1000)^2;  
VR_histo(11,[1,3,5],2) = [dCA1_density(4),dHpc_density(4),all_density(4)]*(1.6^2)*(1000)^2;   
VR_histo(12,[2,4,5],2) = [vCA1_density(4),vHpc_density(4),all_density(4)]*(1.6^2)*(1000)^2;  
VR_histo(13,[2,4,5],2) = [vCA1_density(5),vHpc_density(5),all_density(5)]*(1.6^2)*(1000)^2;  
VR_histo(14,[1,3,5],2) = [dCA1_density(5),dHpc_density(5),all_density(5)]*(1.6^2)*(1000)^2;   
VR_histo(15,[2,4,5],2) = [vCA1_density(6),vHpc_density(6),all_density(6)]*(1.6^2)*(1000)^2;  
%% Physical location of units
for k = 1:size(ks_all_data,1)
    cd(strcat('Y:\Udaysankar Chockanathan\Aim 1\kilosort_results_all\',ks_all_data{k,5})); 
    load(strcat(ks_all_data{k,5},'_chanMap.mat')); 
    cd('preAutoMerge');     
    %find the com of each template
    templates = double(readNPY('templates.npy')); 
    dV = squeeze(max(templates,[],2)-min(templates,[],2)); 
    template_com = dV*[xcoords,ycoords]./[sum(dV,2),sum(dV,2)]; 
    spike_template_assignments = double(readNPY('spike_templates.npy')); 
    spike_cluster_assignments = double(readNPY('spike_clusters.npy')); 
    cluster_com = zeros(size(ks_all_data{k,1},1),2);
    for c = 1:size(cluster_com,1)
        good_clusters = cell2mat(ks_all_data{k,1}(:,1)); 
        cc = good_clusters(c); 
        templates_in_cluster = spike_template_assignments(spike_cluster_assignments == cc);
        template_counts_in_cluster = histcounts(templates_in_cluster,-.5:1:max(spike_template_assignments)+.5,'Normalization','probability'); 
        if k == 17
            tmp1 = template_counts_in_cluster; tmp1(108) = []; 
            tmp2 = template_com; tmp2(108,:) = []; 
            cluster_com(c,:) = tmp1*tmp2; 
        else
            cluster_com(c,:) = template_counts_in_cluster*template_com;
        end
    end   
    ks_all_data{k,6} = cluster_com; 
end    
%% Generate mean waveforms
wave_window = 123;
half_window = floor(wave_window/2); 
for VR_ephys = 1:size(ks_all_data,1)  
    cd(strcat('Y:\Udaysankar Chockanathan\Aim 1\kilosort_results_all\',ks_all_data{VR_ephys,5},'\preAutoMerge'));  
    fileID = fopen(strcat(ks_all_data{VR_ephys,5},'_binary.dat'));    
    ks_all_data{VR_ephys,25} = zeros(size(ks_all_data{VR_ephys,1},1),wave_window,128); %number of units x window size x channels
    for channel = 1:128
        frewind(fileID); 
        if channel > 1 %this reads in the first (channel-1) values to offset the reading appropriately            
            tmp = fread(fileID,[(channel-1),1],'int16'); 
        end
        skip = 2*127; %specify the skip in bytes (int16 has 2 bytes per value)
        precision = '1*int16'; %specify the precision (so the program will read ONE value and skip 128 values)
        bin_data = fread(fileID,[1,ks_all_data{VR_ephys,11}],precision,skip); %load in the current channel raw data (will not be the same as the corresponding amplifier_data row because it gets remapped in spike_sorting_pipeline1_1
        for unit = 1:size(ks_all_data{VR_ephys,1},1) %for each unit
            spk_times = ks_all_data{VR_ephys,1}{unit,2}; %all the spike times (in 30kHz samples) of the current unit
            spk_times(spk_times <= wave_window) = []; %eliminate spikes that are too close to the beginning
            spk_times(spk_times >= ks_all_data{VR_ephys,11}-wave_window) = []; %eliminate spikes that are too close to the end
            num_spikes = numel(spk_times); %number of spikes 
            tmp_waveform = zeros(num_spikes,wave_window); 
            for s = 1:num_spikes
                tmp_waveform(s,:) = bin_data(spk_times(s)-half_window:spk_times(s)+half_window); %store the current waveform
            end
            ks_all_data{VR_ephys,25}(unit,:,channel) = mean(tmp_waveform,1);      
            unit
            channel
            VR_ephys
        end
    end
    fclose(fileID);
end
%% High pass raw waveforms
for VR_ephys = 1:size(ks_all_data,1)
    ks_all_data{VR_ephys,26} = zeros(size(ks_all_data{VR_ephys,25})); 
    for channel = 1:size(ks_all_data{VR_ephys,25},3)
        for unit = 1:size(ks_all_data{VR_ephys,25},1)
            ks_all_data{VR_ephys,26}(unit,:,channel) = highpass(ks_all_data{VR_ephys,25}(unit,:,channel),100,30000);
        end
    end
end
%% Calculate peak-trough times for each unit
for VR_ephys = 1:size(ks_all_data,1)  
    mag_all_channels = squeeze(abs(max(ks_all_data{VR_ephys,26},[],2)-min(ks_all_data{VR_ephys,26},[],2))); %magnitude of max-min of each unit across all channels
    for unit = 1:size(ks_all_data{VR_ephys,26},1)
        [~,w_idx] = max(mag_all_channels(unit,:));
        woi = squeeze(ks_all_data{VR_ephys,26}(unit,:,w_idx));        
        [~,peaktime] = max(abs(woi)); %time of peak 
        loc_extrema = sort([find(islocalmax(woi)),find(islocalmin(woi))]);    
        troughtime = loc_extrema(find(loc_extrema > peaktime,1));
        peak_trough_time = (troughtime-peaktime)/30; %peak-trough time in miliseconds 
        ks_all_data{VR_ephys,27}(unit) = peak_trough_time; 
        ks_all_data{VR_ephys,28}(unit) = peaktime; 
        VR_ephys
        unit
    end             
end
%% Calculate bursting index for each unit
for VR_ephys = 1:size(ks_all_data,1)
    for c = 1:size(ks_all_data{VR_ephys,1},1) %for each good unit
        spk_binaries = zeros(1,ks_all_data{VR_ephys,12});
        spk_binaries(ks_all_data{VR_ephys,1}{c,3}) = 1;     
        for m = 1:numel(ks_all_data{VR_ephys,1}{c,3})
            curr_st = ks_all_data{VR_ephys,1}{c,3}(m);
            if and(curr_st > 601,curr_st < mod_durat(VR_ephys)-601)
                store_corr(m,:) = xcorr(spk_binaries((curr_st-600):(curr_st+600)),spk_binaries((curr_st-600):(curr_st+600)),300); %ACG at individual spikes
            end
        end        
        ACG = mean(store_corr,1); 
        ks_all_data{VR_ephys,38}(c,:) = ACG;
        ks_all_data{VR_ephys,39}(1,c) = (sum(ACG(304:306))/3)/(sum(ACG(501:601))/101);
        VR_ephys
        c
    end
end
%% Separate interneurons and pyramidal cells
trough_peak_time_threshold = 0.5;
burst_index_threshold = 1.9; 
for VR_ephys = 1:size(ks_all_data,1)
    ks_all_data{VR_ephys,40} = []; ks_all_data{VR_ephys,41} = []; ks_all_data{VR_ephys,42} = [];
    for c = 1:size(ks_all_data{VR_ephys,1},1) %for each unit
        if ks_all_data{VR_ephys,27}(c) > trough_peak_time_threshold
            ks_all_data{VR_ephys,40} = [ks_all_data{VR_ephys,40};ks_all_data{VR_ephys,1}(c,:)]; %put the spike time data for excitatory cells in one column           
            ks_all_data{VR_ephys,42} = [ks_all_data{VR_ephys,42};1]; 
        else
            ks_all_data{VR_ephys,41} = [ks_all_data{VR_ephys,41};ks_all_data{VR_ephys,1}(c,:)]; %put the spike time data for interneurons in another column
            ks_all_data{VR_ephys,42} = [ks_all_data{VR_ephys,42};0]; 
        end
    end
end
%% Plot trough-peak time vs. bursting index
g = figure; set(g,'Position',[0 0 1000 1000]); 
plot_x = cell2mat(ks_all_data(:,27)');
plot_y = cell2mat(ks_all_data(:,39)');
subplot(3,3,1); 
for i = 1:numel(cell2mat(ks_all_data(:,27)'))
    if and(plot_x(i) > 0.5,plot_y(i) > 1.9)
        scatter(plot_x(i)+randn*.005,plot_y(i),.5,[1,165/255,0],'filled'); hold on;
    else
        scatter(plot_x(i)+randn*.005,plot_y(i),.5,'b','filled'); hold on;
    end
end    
ylim([0 5]); xlabel('Trough-peak time (ms)'); ylabel('Burst index'); 
f = figure; set(f,'Position',[0 0 2000 1500]); 
subplot(3,3,2); histogram(cell2mat(ks_all_data(:,27)'),[0:.04:1.6],'Normalization','probability','EdgeColor','none','FaceColor','k','FaceAlpha',1); xlim([0 1.6]); xlabel('Trough-peak time(ms)'); ylabel('Probability'); 
subplot(3,3,3); histogram(cell2mat(ks_all_data(:,39)'),[0:.05:5],'Normalization','probability','EdgeColor','none','FaceColor','k','FaceAlpha',1); xlim([0 5]); xlabel('Burst index'); ylabel('Probability'); 
plot_color = []; plot_x = []; plot_y = [];
for VR_ephys = 1:size(ks_all_data,1)
    if any(VR_ephys == WT_dCA1_sites); plot_color = [plot_color;repmat([0 0 0],numel(ks_all_data{VR_ephys,27}),1)];
    elseif any(VR_ephys == AD_dCA1_sites); plot_color = [plot_color;repmat([1 0 0],numel(ks_all_data{VR_ephys,27}),1)];
    end
    if any(VR_ephys == [WT_dCA1_sites,AD_dCA1_sites])
        plot_x = [plot_x;ks_all_data{VR_ephys,27}'+randn(numel(cell2mat(ks_all_data(VR_ephys,27))),1)*.005];
        plot_y = [plot_y;ks_all_data{VR_ephys,39}']; 
    end
end
neuron_order = randperm(numel(plot_x));
subplot(3,3,4); scatter(plot_x(neuron_order),plot_y(neuron_order),6,plot_color(neuron_order,:),'filled'); set(gca,'YScale','log'); xlim([0 1.6]); ylim([1e-1 7]); xlabel('Trough-peak time (ms)'); ylabel('Burst index');
subplot(3,3,5)
histogram(cell2mat(ks_all_data(WT_dCA1_sites,27)'),[0:.01:1.6],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(cell2mat(ks_all_data(AD_dCA1_sites,27)'),[0:.01:1.6],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
xlabel('Trough-peak time (ms)'); ylabel('Cumulative Probability'); xlim([0 1.6]); ylim([0 1]); legend({'WT dCA1','AD dCA1'},'location','southeast'); legend boxoff;
[p,~] = ranksum(cell2mat(ks_all_data(WT_dCA1_sites,27)'),cell2mat(ks_all_data(AD_dCA1_sites,27)')); text(1,.8,strcat('p = ',num2str(p))); 
subplot(3,3,6)
histogram(cell2mat(ks_all_data(WT_dCA1_sites,39)'),[0:.01:5],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(cell2mat(ks_all_data(AD_dCA1_sites,39)'),[0:.01:5],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
xlabel('Burst index'); ylabel('Cumulative Probability'); xlim([0 5]); ylim([0 1]); legend({'WT dCA1','AD dCA1'},'location','southeast'); legend boxoff;
[p,~] = ranksum(cell2mat(ks_all_data(WT_dCA1_sites,39)'),cell2mat(ks_all_data(AD_dCA1_sites,39)')); text(1,.8,strcat('p = ',num2str(p))); 
plot_color = []; plot_x = []; plot_y = [];
for VR_ephys = 1:size(ks_all_data,1)
    if any(VR_ephys == WT_vCA1_sites); plot_color = [plot_color;repmat([.8 .8 .8],numel(ks_all_data{VR_ephys,27}),1)];
    elseif any(VR_ephys == AD_vCA1_sites); plot_color = [plot_color;repmat([1 .8 .8],numel(ks_all_data{VR_ephys,27}),1)];
    end
    if any(VR_ephys == [WT_vCA1_sites,AD_vCA1_sites])
        plot_x = [plot_x;ks_all_data{VR_ephys,27}'+randn(numel(cell2mat(ks_all_data(VR_ephys,27))),1)*.005]; 
        plot_y = [plot_y;ks_all_data{VR_ephys,39}']; 
    end
end
neuron_order = randperm(numel(plot_x)); 
subplot(3,3,7); scatter(plot_x(neuron_order),plot_y(neuron_order),6,plot_color(neuron_order,:),'filled'); set(gca,'YScale','log'); xlim([0 1.6]); ylim([1e-1 7]); xlabel('Trough-peak time (ms)'); ylabel('Burst index');
subplot(3,3,8)
histogram(cell2mat(ks_all_data(WT_vCA1_sites,27)'),[0:.01:1.6],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(cell2mat(ks_all_data(AD_vCA1_sites,27)'),[0:.01:1.6],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
xlabel('Trough-peak time (ms)'); ylabel('Cumulative Probability'); xlim([0 1.6]); ylim([0 1]); legend({'WT vCA1','AD vCA1'},'location','southeast'); legend boxoff;
[p,~] = ranksum(cell2mat(ks_all_data(WT_vCA1_sites,27)'),cell2mat(ks_all_data(AD_vCA1_sites,27)')); text(1,.8,strcat('p = ',num2str(p))); 
subplot(3,3,9)
histogram(cell2mat(ks_all_data(WT_vCA1_sites,39)'),[0:.01:5],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(cell2mat(ks_all_data(AD_vCA1_sites,39)'),[0:.01:5],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
xlabel('Burst index'); ylabel('Cumulative Probability'); xlim([0 5]); ylim([0 1]); legend({'WT vCA1','AD vCA1'},'location','southeast'); legend boxoff;
[p,~] = ranksum(cell2mat(ks_all_data(WT_vCA1_sites,39)'),cell2mat(ks_all_data(AD_vCA1_sites,39)')); text(1,.8,strcat('p = ',num2str(p))); 
%Example waveforms
h = figure;
VR_ephys = 17; 
mag_all_channels = squeeze(abs(max(ks_all_data{VR_ephys,26},[],2)-min(ks_all_data{VR_ephys,26},[],2))); %magnitude of max-min of each unit across all channels
for unit = [1,57]
    [~,w_idx] = max(mag_all_channels(unit,:));
    woi = squeeze(ks_all_data{VR_ephys,26}(unit,:,w_idx));        
    [~,peaktime] = max(abs(woi)); %time of peak 
    loc_extrema = sort([find(islocalmax(woi)),find(islocalmin(woi))]);    
    troughtime = loc_extrema(find(loc_extrema > peaktime,1));
    peak_trough_time = (troughtime-peaktime)/30; %peak-trough time in miliseconds 
    plot(1:123,woi,'LineWidth',2); hold on;
    text(10,(-1)*unit+10,num2str([peak_trough_time,peaktime,troughtime])); 
end
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 3');
saveas(f,'cell_class','fig'); print('cell_class.pdf','-bestfit','-painters','-dpdf'); close(f);
saveas(g,'cell_class_scatter','fig'); print('cell_class_scatter.pdf','-bestfit','-painters','-dpdf'); close(g);
saveas(h,'example_waveforms','fig'); print('example_waveforms.pdf','-bestfit','-painters','-dpdf'); close(h);
%% Plot mean waveforms, rasters, raw data intervals, and extracted spikes
f = figure; set(f,'Position',[0 0 400 1000]); subplot(14,1,1:5); 
soi = 12; %AD vCA1: 12; %AD dCA1: 1; %WT vCA1: 19; %WT dCA1: 18; %session of interest
raw_limits = 30000*[2254.15,2254.4]; %AD vCA1: 30000*[2254.15,2254.4]; %AD dCA1: 30000*[2290.27,2290.52]; %WT vCA1: [15179700,15187200]; %WT dCA1: [72171300,72178800];
uoi = [103,318]; %AD vCA1: [103,318]; %AD dCA1: [542,255]; %WT vCA1: [31,70]; %WT dCA1: [3,7]; %units of interest
limits2 = [475000 480000]; %WT dCA1: [450000 455000];%WT vCA1: [670000 675000]; %AD dCA1: [475000 480000]; %AD vCA1: [475000 480000]; 
cd(strcat('Y:\Udaysankar Chockanathan\Aim 1\kilosort_results_all\',ks_all_data{soi,5},'\preAutoMerge'));
fileID = fopen(strcat(ks_all_data{soi,5},'_binary.dat')); 
fseek(fileID,2*128*raw_limits(1),'bof'); 
rd_trunc = fread(fileID,[128,diff(raw_limits)],'int16'); 
fclose(fileID); 
rd_coi = rd_trunc([67,93,68,92],:)';%AD vCA1: rd_trunc([67,93,68,92],:)'; %AD dCA1: rd_trunc([50,51,64,63],:)'; %WT vCA1: rd_trunc([10,22,11,21],:)'; %WT dCA1: rd_trunc([5,27,6,26],:)'; %channels on which units of interest are best seen
cd(strcat('Y:\Udaysankar Chockanathan\Aim 1\kilosort_results_all\',ks_all_data{soi,5})); load(strcat(ks_all_data{soi,5},'_chanMap.mat')); cd('preAutoMerge');     
%find the templates and the com of each template
templates = double(readNPY('templates.npy')); 
dV = squeeze(max(templates,[],2)-min(templates,[],2)); 
template_com = dV*[xcoords,ycoords]./[sum(dV,2),sum(dV,2)]; 
spike_template_assignments = double(readNPY('spike_templates.npy')); 
spike_cluster_assignments = double(readNPY('spike_clusters.npy'));     
cluster_com = zeros(size(ks_all_data{soi,1},1),2);
cluster_waveform = zeros(size(templates,3),size(templates,2),size(ks_all_data{soi,1},1)); %number of channels by number of timepoints by number of good units
%Plot the dots of the electrode array
xcoords_mod = xcoords; xcoords_mod(33:end) = xcoords_mod(33:end)-60; xcoords_mod(65:end) = xcoords_mod(65:end)-60; xcoords_mod(97:end) = xcoords_mod(97:end)-60; 
scatter(xcoords_mod,2*ycoords,18,[.5 .5 .5],'filled'); hold on;
ycoords_mod = ycoords; ycoords_mod(33:end) = ycoords_mod(33:end)-1; ycoords_mod(65:end) = ycoords_mod(65:end)-1; ycoords_mod(97:end) = ycoords_mod(97:end)-1; 
wdw_size = size(templates,2); 
d_shank_colors = jet(size(ks_all_data{soi,1},1)); d_shank_colors = d_shank_colors(randperm(size(d_shank_colors,1)),:);
for c = 1:size(ks_all_data{soi,1},1) %for each good unit
    good_clusters = cell2mat(ks_all_data{soi,1}(:,1)); 
    cc = good_clusters(c); 
    templates_in_cluster = spike_template_assignments(spike_cluster_assignments == cc);
    template_counts_in_cluster = histcounts(templates_in_cluster,-.5:1:max(spike_template_assignments)+.5,'Normalization','probability');        
    cluster_com(c,:) = template_counts_in_cluster*template_com;
    for i = 1:size(templates,3)
        cluster_waveform(i,:,c) = template_counts_in_cluster*templates(:,:,i);
    end        
    [~,idx] = mink(sum(([xcoords,ycoords]-cluster_com(c,:)).^2,2),6); 
    for i3 = 1:numel(idx)
        plot(xcoords_mod(idx(i3)):.3:(xcoords_mod(idx(i3))+wdw_size*.3-.3),2*ycoords(idx(i3))+200*cluster_waveform(idx(i3),:,c),'Color',d_shank_colors(c,:),'LineWidth',2);
        hold on; 
    end
    xticks([]); yticks([]); xlim([-25 351]); ylim([-1700 50]); 
end   
%Plot rasters
subplot(14,1,[6:10]); 
for c = 1:size(ks_all_data{soi,1},1)
    tmpA = ks_all_data{soi,1}{c,3};
    in_interval_spikes2 = tmpA(and(tmpA > limits2(1),tmpA < limits2(2))); 
    rasterPlotter(in_interval_spikes2,size(ks_all_data{soi,1},1)-c+1,d_shank_colors(c,:)); hold on; 
    ylim([0 size(ks_all_data{soi,1},1)+1]); yticks([]); yticklabels([]); ylabel('Units'); xlabel('Time (ms)'); 
end
%Plot raw data interval
cd ..; load(strcat(ks_all_data{soi,5},'_rez.mat')); 
raw_colors = [[200/255 0 1];[0 161/255 0];[0 0 1];[1 0 0];[.8 .8 .8]];
wdw_size = 123;
half_wdw_size = 61;
%Plot running behavior
subplot(14,1,11);
plot([0:10:diff(limits2)],ks_all_data{soi,8}(limits2(1)/10:limits2(2)/10),'k','LineWidth',2);
xlim([0 diff(limits2)]); xlabel('Time (ms)'); ylabel('Velocity (cm/s)'); 
%Identify in-interval spikes
in_interval_spikes = cell(numel(uoi),1);
for j = 1:numel(uoi)
    tmp_spikes = rez.st3((rez.st3(:,2) == uoi(j)+1),1); 
    in_interval_spikes{j} = tmp_spikes(and(tmp_spikes >= raw_limits(1),tmp_spikes <= raw_limits(2)));     
end
subplot(14,1,12:13); 
for i = 1:4
    plot((1:size(rd_coi,1))/30000,(i-1)*300+rd_coi(:,i),'Color',[.5 .5 .5]); hold on;    
end
for i = 1:4
    for j = 1:numel(uoi)
        adj_spikes_tmp = in_interval_spikes{j}-raw_limits(1); 
        for k = 1:numel(adj_spikes_tmp)
            ioi = adj_spikes_tmp(k)-half_wdw_size:adj_spikes_tmp(k)+half_wdw_size;
            ioi(ioi <= 0) = []; ioi(ioi > 7500) = [];
            plot(ioi/30000,(i-1)*300+rd_coi(ioi,i),'Color',raw_colors(j,:)); hold on; 
        end
    end
end
xlim([0 size(rd_coi,1)/30000]); ylim([-500 1200]); xlabel('Time (s)'); ylabel('Voltage (\muV)');
%Plot extracted spikes
subplot(14,1,14); 
for j = 1:numel(uoi)
    rasterPlotter(in_interval_spikes{j},j,raw_colors(j,:)); hold on;
end
xlim(raw_limits)
% cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figures 1 and 2'); 
% saveas(f,'mean_waveforms_rasters1','fig'); print('mean_waveforms_rasters1.pdf','-bestfit','-painters','-dpdf'); close(f);
%% CCG
%Make CCG
CCG = cell(numel(uoi)); 
spk_binaries = zeros(2,ks_all_data{soi,12}); 
tmp_spikes_store = cell(numel(uoi),1); 
for j = 1:numel(uoi)
    tmp_spikes = rez.st3((rez.st3(:,2) == uoi(j)+1),1); 
    tmp_spikes = round(tmp_spikes/ks_all_data{soi,4});
    tmp_spikes_store{j} = tmp_spikes; 
    spk_binaries(j,tmp_spikes) = 1;
end
for j = 1:numel(uoi)
    for k = 1:numel(uoi)
        store_corr = zeros(numel(tmp_spikes_store{j}),51);
        for m = 1:numel(tmp_spikes_store{j})      
            curr_st = tmp_spikes_store{j}(m);
            if and(curr_st > 51,curr_st < mod_durat(soi)-51)
                store_corr(m,:) = xcorr(spk_binaries(j,(curr_st-50):(curr_st+50)),spk_binaries(k,(curr_st-50):(curr_st+50)),25);                            
            end
        end
        CCG{j,k} = mean(store_corr,1); 
    end
end
%Plot CCG
g = figure; set(gcf,'Position',[0 0 800 800]); 
for j = 1:numel(uoi)^2
    subplot(2,2,j); 
    if or(j == 1,j == 4)
        bar(-25:-1,CCG{j}(1:25),'EdgeColor','none','FaceColor',raw_colors(sqrt(j),:)); hold on;
        bar(1:25,CCG{j}(27:51),'EdgeColor','none', 'FaceColor',raw_colors(sqrt(j),:)); hold on; 
    else
        bar(-25:25,CCG{j},'EdgeColor','none','FaceColor',[.5 .5 .5]); 
    end
    axis tight; axis square; xlabel('Lag (ms)'); ylabel('Correlation')
end
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figures 1 and 2');
saveas(g,'CCG_4','fig'); print('CCG_4.pdf','-bestfit','-painters','-dpdf'); close(g); 
%% Plot running behavior
f = figure; set(gcf,'Position',[0 0 500 1500]); 
subplot(6,1,1); 
r_limits = [650000 680000]; soi = 2;
plot([0:10:diff(r_limits)],ks_all_data{soi,8}(r_limits(1)/10:r_limits(2)/10),'k','LineWidth',2);
xlim([0 diff(r_limits)]); xlabel('Time (ms)'); ylabel('Velocity (cm/s)'); 
subplot(6,1,2); 
r_limits = [410000 440000]; soi = 17;
plot([0:10:diff(r_limits)],ks_all_data{soi,8}(r_limits(1)/10:r_limits(2)/10),'r','LineWidth',2);
xlim([0 diff(r_limits)]); xlabel('Time (ms)'); ylabel('Velocity (cm/s)'); 
percent_running = zeros(size(ks_all_data,1),1); 
running_velocity = percent_running; 
for VR_ephys = 1:size(ks_all_data,1)
    percent_running(VR_ephys) = sum(ks_all_data{VR_ephys,8} > 0)/numel(ks_all_data{VR_ephys,8});
    running_velocity(VR_ephys) = mean(ks_all_data{VR_ephys,8}(ks_all_data{VR_ephys,8} > 0));
end
subplot(6,1,3:4); 
scatter(ones(numel(WT_sites),1),percent_running(WT_sites),72,'k','filled'); hold on;
scatter(3*ones(numel(AD_sites),1),percent_running(AD_sites),72,'k','filled');
[p,~] = ranksum(percent_running(WT_sites),percent_running(AD_sites)); text(2,.3,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 .55]); xticks([1 3]); xticklabels({'WT','APP/PS1'}); ylabel('Fraction of time spent running'); 
subplot(6,1,5:6); 
scatter(ones(numel(WT_sites),1),running_velocity(WT_sites),72,'k','filled'); hold on;
scatter(3*ones(numel(AD_sites),1),running_velocity(AD_sites),72,'k','filled');
[p,~] = ranksum(running_velocity(WT_sites),running_velocity(AD_sites)); text(2,14,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 16]); xticks([1 3]); xticklabels({'WT','APP/PS1'}); ylabel('Running velocity (cm/s)'); 
% cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figures 1 and 2');
% saveas(f,'Running','fig'); print('Running.pdf','-bestfit','-painters','-dpdf'); close(f); 
%% Identify start and end points of laps
forward_lap_bounds = cell(size(VR_reformatting,1),1); 
backward_lap_bounds = forward_lap_bounds; 
for VR_ephys = 1:size(VR_reformatting,1)     
    % Reformat displacement so it matches frames
    tmp_displacement = ks_all_data{VR_reformatting{VR_ephys,5},7}';
    tmp2_displacement = repmat(tmp_displacement,10,1); tmp2_displacement = tmp2_displacement(:)';
    tmp3_displacement = [tmp2_displacement,tmp2_displacement(end)*ones(1,ks_all_data{VR_reformatting{VR_ephys,5},12}-numel(tmp2_displacement))];
    displacement = tmp3_displacement(VR_reformatting{VR_ephys,6});
    % Identify boundary points between laps
    start_frame = min(VR_reformatting{VR_ephys,4}); 
    end_frame = max(VR_reformatting{VR_ephys,4}); 
    tmp_start_idx = find(diff(VR_reformatting{VR_ephys,4} == start_frame) == -1); %end indices of "islands" of start frames
    tmp_end_idx = find(diff(VR_reformatting{VR_ephys,4} == end_frame) == 1)+1; %start indices of "islands" of end frames
    potential_laps = []; 
    t = 1;
    while t < numel(tmp_start_idx)
        tmpA = tmp_end_idx-tmp_start_idx(t); %number of frames between current start frame island and every end frame island    
        tmpB = tmpA(tmpA > 0); %same as tmpA, but only for future end frame islands
        if isempty(tmpB) == 0 %if you're not on the final time of the start frame (at end of recording)
            tmpC = find(tmpA == min(tmpB)); %index of the first future end frame island (tmpA index is analogous to tmp_end_idx index; it's an index of indices)
            if all(tmp_end_idx(tmpC) < tmp_start_idx((t+1):end)) %if the closest end frame island is closer than all future start frame islands
                if tmp_end_idx(tmpC)-tmp_start_idx(t) ~= 1 %if the lap didn't just involve a backwards "jump" in frame from start to end frame
                    if displacement(tmp_end_idx(tmpC)) > displacement(tmp_start_idx(t))
                        forward_lap_bounds{VR_ephys} = [forward_lap_bounds{VR_ephys};[tmp_start_idx(t),tmp_end_idx(tmpC)]]; 
                    else
                        backward_lap_bounds{VR_ephys} = [backward_lap_bounds{VR_ephys};[tmp_start_idx(t),tmp_end_idx(tmpC)]]; 
                    end                                 
                end
            end
        end
        t = t+1; 
    end    
end
num_laps = cellfun(@numel,forward_lap_bounds)/2;
f = figure; set(gcf,'Position',[0 0 500 1500]); 
subplot(6,1,1:2); 
scatter(ones(numel(WT_sites),1),num_laps(WT_sites),72,'k','filled'); hold on;
scatter(3*ones(numel(AD_sites),1),num_laps(AD_sites),72,'k','filled');
[p,~] = ranksum(num_laps(WT_sites),num_laps(AD_sites)); text(2,14,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([-5 110]); xticks([1 3]); xticklabels({'C57BL/6','APP/PS1'}); ylabel('Number of laps completed'); 
% cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figures 1 and 2');
% saveas(f,'Num_laps','fig'); print('Num_laps.pdf','-bestfit','-painters','-dpdf'); close(f);
%% Calculate consistent firing rate and correlations
converter_original = 2.^(0:99); 
all_bin_sizes = [5,10,25,50,75,100]; 
all_word_lengths = 16; %[4:2:24]; %even though word lengths aren't explicitly necessary for calculating FR/corr, it tells us which small-population sessions to exclude to keep it consistent with entropy/KLD analyses
mod_rasters_FR_corr = cell(size(VR_reformatting,1),1); 
FR_corr_VR_frames = cell(size(VR_reformatting,1),1); 
corr_consistent = cell(size(VR_reformatting,1),numel(converter_original),numel(all_bin_sizes));
corr_consistent_pyramidal = corr_consistent; corr_consistent_interneuron = corr_consistent; 
FR_consistent = corr_consistent; 
for bin_counter = 1:numel(all_bin_sizes)
    FR_corr_bin_size = all_bin_sizes(bin_counter); 
    for word_length_counter = 1:numel(all_word_lengths)
        word_length = all_word_lengths(word_length_counter); 
        for VR_ephys = 1:size(VR_reformatting,1) 
            for i = 1:size(ks_all_data{VR_reformatting{VR_ephys,5},1},1)
                all_st = ks_all_data{VR_reformatting{VR_ephys,5},1}{i,3}; %all spike times (in ms) for a particular unit
                in_interval_st = all_st(and(all_st >= VR_reformatting{VR_ephys,6}(1),all_st <= VR_reformatting{VR_ephys,6}(end))); %only spike times that fall within bounds of current ephys session
                in_interval_st = in_interval_st-VR_reformatting{VR_ephys,6}(1)+1; %correct for the fact that ephys session might not start at first milisecond. 
                mod_rasters_FR_corr{VR_ephys}{i,1} = in_interval_st; 
            end
            tmp1_VR_frames = VR_reformatting{VR_ephys,4}; 
            if mod(numel(tmp1_VR_frames),FR_corr_bin_size) == 0
                excess = 0; 
            else
                excess = numel(tmp1_VR_frames)-FR_corr_bin_size*floor(numel(tmp1_VR_frames)/FR_corr_bin_size);            
            end
            tmp2_VR_frames = tmp1_VR_frames(1:(end-excess)); 
            tmp3_VR_frames = reshape(tmp2_VR_frames,FR_corr_bin_size,numel(tmp2_VR_frames)/FR_corr_bin_size); 
            tmp4_VR_frames = mode(tmp3_VR_frames,1); 
            FR_corr_VR_frames{VR_ephys} = tmp4_VR_frames;    
            if size(mod_rasters_FR_corr{VR_ephys},1) >= word_length
                %Generate binary spike trains from array of spike times for all units
                spike_train_binned_all = zeros(size(mod_rasters_FR_corr{VR_ephys},1),numel(FR_corr_VR_frames{VR_ephys})); %matrix of 0's - #units by # of bins in recording 
                for units = 1:size(mod_rasters_FR_corr{VR_ephys},1) %for each unit 
                    [spike_train_binned_all(units,:),~] = histcounts(mod_rasters_FR_corr{VR_ephys}{units,1},0.5:FR_corr_bin_size:(FR_corr_bin_size*numel(FR_corr_VR_frames{VR_ephys})+.5)); %count the number of spikes that occur in each bin, starting at .5ms
                end   
                spike_train_binned_all = spike_train_binned_all > 0; %this sets all the bins with at least 1 spike to 1 (and leaves the rest as 0's)
                spike_train_binned_all = zscore(spike_train_binned_all,0,2); 
                tmp_corr = corr(spike_train_binned_all'); tmp_corr = tmp_corr-eye(size(tmp_corr,1)).*tmp_corr;
                corr_consistent{VR_ephys,word_length,bin_counter} = squareform(tmp_corr); 
                tmp_pyramidal = tmp_corr(find(ks_all_data{VR_ephys,42}),:); tmp_pyramidal = tmp_pyramidal(:,find(ks_all_data{VR_ephys,42}));
                corr_consistent_pyramidal{VR_ephys,word_length,bin_counter} = squareform(tmp_pyramidal);
                tmp_interneuron = tmp_corr(find(~ks_all_data{VR_ephys,42}),:); tmp_interneuron = tmp_interneuron(:,find(~ks_all_data{VR_ephys,42}));
                corr_consistent_interneuron{VR_ephys,word_length,bin_counter} = squareform(tmp_interneuron);
                FR_consistent{VR_ephys,word_length,bin_counter} = ((sum(spike_train_binned_all,2)/size(spike_train_binned_all,2))*1000/FR_corr_bin_size)';                 
            end           
            FR_corr_bin_size
            VR_ephys
        end
    end
end
%% Plot Consistent Mean Firing rate
FR = cell(size(ks_all_data,1),1); 
FR_pyramidal = FR; FR_interneuron = FR;
for an = 1:size(ks_all_data,1)    
    FR{an} = cellfun(@numel,mod_rasters{an})*1000/double(ks_all_data{an,12});
    if isempty(ks_all_data{an,40}) == 0
        FR_pyramidal{an} = cellfun(@numel,ks_all_data{an,40}(:,3))*1000/double(ks_all_data{an,12});        
    end
    FR_interneuron{an} = cellfun(@numel,ks_all_data{an,41}(:,3))*1000/double(ks_all_data{an,12}); 
end
f = figure; set(f,'Position',[0 0 1500 1500]); 
subplot(3,3,1); 
histogram(cell2mat(FR(WT_dCA1_sites)),[0:.01:75],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(cell2mat(FR(AD_dCA1_sites)),[0:.01:75],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
xlabel('FR (Hz)'); ylabel('Cumulative Probability'); xlim([0 70]); ylim([0 1]); legend({'WT dCA1','AD dCA1'},'location','southeast'); legend boxoff;
[p,~] = ranksum(cell2mat(FR(WT_dCA1_sites)),cell2mat(FR(AD_dCA1_sites))); text(20,.8,strcat('p = ',num2str(p))); 
subplot(3,3,4); 
histogram(cell2mat(FR(WT_vCA1_sites)),[0:.01:75],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(cell2mat(FR(AD_vCA1_sites)),[0:.01:75],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); 
xlabel('FR (Hz)'); ylabel('Cumulative Probability'); xlim([0 70]); ylim([0 1]); legend({'WT vCA1','AD vCA1'},'location','southeast'); legend boxoff;
[p,~] = ranksum(cell2mat(FR(WT_vCA1_sites)),cell2mat(FR(AD_vCA1_sites))); text(20,.6,strcat('p = ',num2str(p))); 
n = [numel(cell2mat(FR(WT_dCA1_sites))),numel(cell2mat(FR(WT_vCA1_sites))),numel(cell2mat(FR(AD_dCA1_sites))),numel(cell2mat(FR(AD_vCA1_sites)))]; text(40,.5,num2str(n))
subplot(3,3,2); 
histogram(cell2mat(FR_pyramidal(WT_dCA1_sites)),[0:.01:75],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(cell2mat(FR_pyramidal(AD_dCA1_sites)),[0:.01:75],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
xlabel('FR (Hz)'); ylabel('Cumulative Probability'); xlim([0 70]); ylim([0 1]); legend({'WT dCA1','AD dCA1'},'location','southeast'); legend boxoff;
[p,~] = ranksum(cell2mat(FR_pyramidal(WT_dCA1_sites)),cell2mat(FR_pyramidal(AD_dCA1_sites))); text(20,.8,strcat('p = ',num2str(p))); 
subplot(3,3,5); 
histogram(cell2mat(FR_pyramidal(WT_vCA1_sites)),[0:.01:75],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(cell2mat(FR_pyramidal(AD_vCA1_sites)),[0:.01:75],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); 
xlabel('FR (Hz)'); ylabel('Cumulative Probability'); xlim([0 70]); ylim([0 1]); legend({'WT vCA1','AD vCA1'},'location','southeast'); legend boxoff;
[p,~] = ranksum(cell2mat(FR_pyramidal(WT_vCA1_sites)),cell2mat(FR_pyramidal(AD_vCA1_sites))); text(20,.6,strcat('p = ',num2str(p))); 
n = [numel(cell2mat(FR_pyramidal(WT_dCA1_sites))),numel(cell2mat(FR_pyramidal(WT_vCA1_sites))),numel(cell2mat(FR_pyramidal(AD_dCA1_sites))),numel(cell2mat(FR_pyramidal(AD_vCA1_sites)))]; text(40,.5,num2str(n))
subplot(3,3,3); 
histogram(cell2mat(FR_interneuron(WT_dCA1_sites)),[0:.01:75],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(cell2mat(FR_interneuron(AD_dCA1_sites)),[0:.01:75],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
xlabel('FR (Hz)'); ylabel('Cumulative Probability'); xlim([0 70]); ylim([0 1]); legend({'WT dCA1','AD dCA1'},'location','southeast'); legend boxoff;
[p,~] = ranksum(cell2mat(FR_interneuron(WT_dCA1_sites)),cell2mat(FR_interneuron(AD_dCA1_sites))); text(20,.8,strcat('p = ',num2str(p))); 
subplot(3,3,6); 
histogram(cell2mat(FR_interneuron(WT_vCA1_sites)),[0:.01:75],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(cell2mat(FR_interneuron(AD_vCA1_sites)),[0:.01:75],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); 
xlabel('FR (Hz)'); ylabel('Cumulative Probability'); xlim([0 70]); ylim([0 1]); legend({'WT vCA1','AD vCA1'},'location','southeast'); legend boxoff;
[p,~] = ranksum(cell2mat(FR_interneuron(WT_vCA1_sites)),cell2mat(FR_interneuron(AD_vCA1_sites))); text(20,.6,strcat('p = ',num2str(p))); 
n = [numel(cell2mat(FR_interneuron(WT_dCA1_sites))),numel(cell2mat(FR_interneuron(WT_vCA1_sites))),numel(cell2mat(FR_interneuron(AD_dCA1_sites))),numel(cell2mat(FR_interneuron(AD_vCA1_sites)))]; text(40,.5,num2str(n))
subplot(3,3,7);
histogram(cell2mat(FR(WT_dCA1_sites)),[0:.01:75],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(cell2mat(FR(WT_vCA1_sites)),[0:.01:75],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); 
xlabel('FR (Hz)'); ylabel('Cumulative Probability'); xlim([0 70]); ylim([0 1]); legend({'WT dCA1','WT vCA1'},'location','southeast'); legend boxoff;
[p,~] = ranksum(cell2mat(FR(WT_dCA1_sites)),cell2mat(FR(WT_vCA1_sites))); text(20,.6,strcat('p = ',num2str(p))); 
subplot(3,3,8);
histogram(cell2mat(FR_pyramidal(WT_dCA1_sites)),[0:.01:75],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(cell2mat(FR_pyramidal(WT_vCA1_sites)),[0:.01:75],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); 
xlabel('FR (Hz)'); ylabel('Cumulative Probability'); xlim([0 70]); ylim([0 1]); legend({'WT dCA1','WT vCA1'},'location','southeast'); legend boxoff;
[p,~] = ranksum(cell2mat(FR_pyramidal(WT_dCA1_sites)),cell2mat(FR_pyramidal(WT_vCA1_sites))); text(20,.6,strcat('p = ',num2str(p))); 
subplot(3,3,9);
histogram(cell2mat(FR_interneuron(WT_dCA1_sites)),[0:.01:75],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(cell2mat(FR_interneuron(WT_vCA1_sites)),[0:.01:75],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); 
xlabel('FR (Hz)'); ylabel('Cumulative Probability'); xlim([0 70]); ylim([0 1]); legend({'WT dCA1','WT vCA1'},'location','southeast'); legend boxoff;
[p,~] = ranksum(cell2mat(FR_interneuron(WT_dCA1_sites)),cell2mat(FR_interneuron(WT_vCA1_sites))); text(20,.6,strcat('p = ',num2str(p))); 
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 4');
saveas(f,'mean_FR','fig'); print('mean_FR.pdf','-bestfit','-painters','-dpdf'); close(f);
%% Compare variance of mean FRs
num_bootstrap_samples = 100000; 
FR_WT_dCA1_consol = cell2mat(FR(WT_dCA1_sites));
FR_AD_dCA1_consol = cell2mat(FR(AD_dCA1_sites));
FR_WT_vCA1_consol = cell2mat(FR(WT_vCA1_sites));
FR_AD_vCA1_consol = cell2mat(FR(AD_vCA1_sites));
WT_dCA1_bs_idx = randi(numel(FR_WT_dCA1_consol),num_bootstrap_samples,numel(FR_WT_dCA1_consol)); 
AD_dCA1_bs_idx = randi(numel(FR_AD_dCA1_consol),num_bootstrap_samples,numel(FR_AD_dCA1_consol)); 
WT_vCA1_bs_idx = randi(numel(FR_WT_vCA1_consol),num_bootstrap_samples,numel(FR_WT_vCA1_consol)); 
AD_vCA1_bs_idx = randi(numel(FR_AD_vCA1_consol),num_bootstrap_samples,numel(FR_AD_vCA1_consol)); 
var_WT_AD_dCA1 = zeros(num_bootstrap_samples,1); 
var_WT_AD_vCA1 = var_WT_AD_dCA1;
var_WT_dCA1_vCA1 = var_WT_AD_dCA1;
for b = 1:num_bootstrap_samples
    WT_dCA1_bs = FR_WT_dCA1_consol(WT_dCA1_bs_idx(b,:));
    AD_dCA1_bs = FR_AD_dCA1_consol(AD_dCA1_bs_idx(b,:)); 
    WT_vCA1_bs = FR_WT_vCA1_consol(WT_vCA1_bs_idx(b,:)); 
    AD_vCA1_bs = FR_AD_vCA1_consol(AD_vCA1_bs_idx(b,:));     
    var_WT_AD_dCA1(b) = var(WT_dCA1_bs)/var(AD_dCA1_bs); 
    var_WT_AD_vCA1(b) = var(WT_vCA1_bs)/var(AD_vCA1_bs); 
    var_WT_dCA1_vCA1(b) = var(WT_dCA1_bs)/var(WT_vCA1_bs); 
end
f = figure; set(f,'Position',[0 300 1500 500]); 
subplot(1,3,1); histogram(var_WT_dCA1_vCA1,[0:.001:3],'normalization','probability','EdgeColor','none','FaceColor',[.8 .8 .8],'FaceAlpha',1); xlabel('<\sigma^2_{WT dCA1}/\sigma^2_{WT vCA1}>'); ylabel('Probability'); 
L_bound = round(prctile(var_WT_dCA1_vCA1,2.5),3); R_bound = round(prctile(var_WT_dCA1_vCA1,97.5),3); 
p = sum(var_WT_dCA1_vCA1 > 1)/num_bootstrap_samples;
text(.5,.001,num2str([prctile(var_WT_dCA1_vCA1,2.5), prctile(var_WT_dCA1_vCA1,97.5), p]));
var_WT_dCA1_vCA1(and(var_WT_dCA1_vCA1 > L_bound,var_WT_dCA1_vCA1 < R_bound)) = 1000;
hold on; histogram(var_WT_dCA1_vCA1,[0:.001:3],'normalization','probability','EdgeColor','none','FaceColor',[0 191/255 1],'FaceAlpha',1); 
hold on; plot([1 1],[0 1],'Color',[1 140/255 0],'LineWidth',1.5); ylim([0 .006]); 
subplot(1,3,2); histogram(var_WT_AD_dCA1,[0:.001:3],'normalization','probability','EdgeColor','none','FaceColor',[.8 .8 .8],'FaceAlpha',1); xlabel('<\sigma^2_{WT dCA1}/\sigma^2_{AD dCA1}>'); ylabel('Probability'); 
L_bound = round(prctile(var_WT_AD_dCA1,2.5),3); R_bound = round(prctile(var_WT_AD_dCA1,97.5),3); 
p = sum(var_WT_AD_dCA1 > 1)/num_bootstrap_samples;
text(.5,.001,num2str([prctile(var_WT_AD_dCA1,2.5), prctile(var_WT_AD_dCA1,97.5), p]));
var_WT_AD_dCA1(and(var_WT_AD_dCA1 > L_bound,var_WT_AD_dCA1 < R_bound)) = 1000;
hold on; histogram(var_WT_AD_dCA1,[0:.001:3],'normalization','probability','EdgeColor','none','FaceColor',[0 191/255 1],'FaceAlpha',1); 
hold on; plot([1 1],[0 1],'Color',[1 140/255 0],'LineWidth',1.5); ylim([0 .006]); 
subplot(1,3,3); histogram(var_WT_AD_vCA1,[0:.001:3],'normalization','probability','EdgeColor','none','FaceColor',[.8 .8 .8],'FaceAlpha',1); xlabel('<\sigma^2_{WT vCA1}/\sigma^2_{AD vCA1}>'); ylabel('Probability'); 
L_bound = round(prctile(var_WT_AD_vCA1,2.5),3); R_bound = round(prctile(var_WT_AD_vCA1,97.5),3); 
p = sum(var_WT_AD_vCA1 > 1)/num_bootstrap_samples;
text(.5,.001,num2str([prctile(var_WT_AD_vCA1,2.5), prctile(var_WT_AD_vCA1,97.5), p]));
var_WT_AD_vCA1(and(var_WT_AD_vCA1 > L_bound,var_WT_AD_vCA1 < R_bound)) = 1000;
hold on; histogram(var_WT_AD_vCA1,[0:.001:3],'normalization','probability','EdgeColor','none','FaceColor',[0 191/255 1],'FaceAlpha',1); 
hold on; plot([1 1],[0 1],'Color',[1 140/255 0],'LineWidth',1.5); ylim([0 .006]); 
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 4'); 
saveas(f,'var_FR','fig'); print('var_FR.pdf','-bestfit','-painters','-dpdf'); close(f); 
%% Compare plaque burden
f = figure; set(gcf,'Position',[0 0 800 800]); 
subplot(3,2,1); scatter(ones(6,1),dHpc_fraction(1:6),72,'b','filled'); hold on; scatter(2*ones(6,1),vHpc_fraction(1:6),72,'r','filled');
[p,~] = signrank(dHpc_fraction,vHpc_fraction); text(1.5,4e-3,strcat('p = ',num2str(p))); xlim([0 3]); ylabel('Fractional plaque area'); xticks([1 2]); xticklabels({'dHpc','vHpc'});
subplot(3,2,2); scatter(ones(6,1),dCA1_fraction(1:6),72,'b','filled'); hold on; scatter(2*ones(6,1),vCA1_fraction(1:6),72,'r','filled');
[p,~] = signrank(dCA1_fraction,vCA1_fraction); text(1.5,.03,strcat('p = ',num2str(p))); xlim([0 3]); ylabel('Fractional plaque area'); xticks([1 2]); xticklabels({'dCA1','vCA1'});
subplot(3,2,3); scatter(ones(6,1),dHpc_density(1:6)*(1.6^2)*(1000)^2,72,'b','filled'); hold on; scatter(2*ones(6,1),vHpc_density(1:6)*(1.6^2)*(1000)^2,72,'r','filled');
[p,~] = signrank(dHpc_density,vHpc_density); text(1.5,1.5e-5,strcat('p = ',num2str(p))); xlim([0 3]); ylabel('Plaque density (per mm^2)'); xticks([1 2]); xticklabels({'dHpc','vHpc'});
subplot(3,2,4); scatter(ones(6,1),dCA1_density(1:6)*(1.6^2)*(1000)^2,72,'b','filled'); hold on; scatter(2*ones(6,1),vCA1_density(1:6)*(1.6^2)*(1000)^2,72,'r','filled');
[p,~] = signrank(dCA1_density,vCA1_density); text(1.5,2.2e-5,strcat('p = ',num2str(p))); xlim([0 3]); ylabel('Plaque density (per mm^2)'); xticks([1 2]); xticklabels({'dCA1','vCA1'});
subplot(3,2,5); scatter(ones(2,1),dCA1_density(7:8)); hold on; scatter(2*ones(2,1),vCA1_density(7:8)); hold on; scatter(3*ones(6,1),dCA1_density(1:6)); hold on; scatter(4*ones(6,1),vCA1_density(1:6)); 
xticks([1,2,3,4]); xticklabels({'dCA1 WT','vCA1 WT','dCA1 AD','vCA1 AD'}); ylabel('Plaque density'); 
for i = 1:6
    subplot(3,2,1); plot([1,2],[dHpc_fraction(i);vHpc_fraction(i)],'k','LineWidth',1); 
    subplot(3,2,2); plot([1,2],[dCA1_fraction(i);vCA1_fraction(i)],'k','LineWidth',1); set(gca,'YScale','log'); ylim([1e-3 1e-1]); 
    subplot(3,2,3); plot([1,2],[dHpc_density(i);vHpc_density(i)]*(1.6^2)*(1000)^2,'k','LineWidth',1); ylim([0 90]); 
    subplot(3,2,4); plot([1,2],[dCA1_density(i);vCA1_density(i)]*(1.6^2)*(1000)^2,'k','LineWidth',1);
end
% cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 3'); 
% saveas(f,'plaque_burden','fig'); print('plaque_burden.pdf','-bestfit','-painters','-dpdf'); close(f);
%% Calculate place field maps during running
spatial_info_running = cell(size(VR_reformatting,1),1); 
spatial_info_running_shuffled = spatial_info_running; 
FR_map_orig_running = spatial_info_running;
FR_map_reordered_by_peak_FR_running = spatial_info_running; 
FR_map_reordered_by_peak_FR_running_idx = spatial_info_running; 
field_size_running = spatial_info_running;
is_place_cell = spatial_info_running; 
is_place_cell3 = spatial_info_running; 
si_mod_vel = spatial_info_running; 
sparsity = spatial_info_running;
prob_in_bin_run_store = spatial_info_running;
in_interval_st_run_store = spatial_info_running;
spatial_bin_size = 20; 
num_frames = 2221; 
spatial_bin_boundaries = 0:spatial_bin_size:num_frames; 
vel_bin_boundaries = 1:1.1:11;
num_shuffles = 100; 
track_length = 186; %in cm
spaces = numel(spatial_bin_boundaries)-1; 
conv_mode = 'square'; 
x = -10:10; sd = 10; gauss_conv = exp(-.5*(x/sd).^2)/(sd*sqrt(2*pi));
for VR_ephys = 1:size(VR_reformatting,1)    
    spatial_info_running{VR_ephys,1} = zeros(size(ks_all_data{VR_reformatting{VR_ephys,5},1},1),1); 
    FR_map_orig_running{VR_ephys,1} = zeros(size(ks_all_data{VR_reformatting{VR_ephys,5},1},1),numel(spatial_bin_boundaries)-5); 
    sparsity{VR_ephys,1} = spatial_info_running{VR_ephys,1};
    %Align running velocity to VR frames
    vel_step = 10; 
    tmp_simv = repmat(mod_vel{VR_reformatting{VR_ephys,5},1},vel_step,1); tmp_simv = tmp_simv(:)';
    if max(VR_reformatting{VR_ephys,6}) > numel(tmp_simv)
        excess = max(VR_reformatting{VR_ephys,6}) - numel(tmp_simv);
        tmp_simv = [tmp_simv,zeros(1,excess)]; 
    end
    spatial_info_mod_vel = tmp_simv(VR_reformatting{VR_ephys,6}) >= 1;   
    si_mod_vel{VR_ephys,1} = tmp_simv(VR_reformatting{VR_ephys,6});
    frames_run_only = VR_reformatting{VR_ephys,4}(spatial_info_mod_vel);
    in_interval_running_ms = sum(spatial_info_mod_vel); %number of ms spent running within the bounds of current ephys sesssion
    for i = 1:size(ks_all_data{VR_reformatting{VR_ephys,5},1},1)
        %Calculate FR in each frame for each unit
        all_st = ks_all_data{VR_reformatting{VR_ephys,5},1}{i,3}; %all spike times (in ms) for a particular unit
        in_interval_st = all_st(and(all_st >= VR_reformatting{VR_ephys,6}(1),all_st <= VR_reformatting{VR_ephys,6}(end))); %only spike times that fall within bounds of current ephys session
        in_interval_st = in_interval_st-VR_reformatting{VR_ephys,6}(1)+1; %correct for the fact that ephys session might not start at first milisecond.                       
        spatial_info_mod_vel2 = VR_reformatting{VR_ephys,10} >= 5; 
        st_stat_or_run = spatial_info_mod_vel2(in_interval_st); %specifies whether each spike in in_interval_st occured during stationary (returns 0) or running (returns 1)
        in_interval_st_run = in_interval_st(st_stat_or_run); %returns only the spikes that occur when running
        in_interval_st_run_store{VR_ephys,1}{i} = in_interval_st_run; 
        [spikes_in_bin_run,~] = histcounts(VR_reformatting{VR_ephys,4}(in_interval_st_run),spatial_bin_boundaries); 
        [time_in_bin_run,~] = histcounts(VR_reformatting{VR_ephys,4}(spatial_info_mod_vel2),spatial_bin_boundaries); 
        FRib_run = 1000*spikes_in_bin_run(:)./time_in_bin_run(:);
        %Smooth and normalize
        smoothed_FRib_run = conv(ones(1,5),FRib_run); %convolution
        smoothed_FRib_run = smoothed_FRib_run(5:(end-4)); %eliminate points of partial overlap
        smoothed_norm_FRib_run = (smoothed_FRib_run-min(smoothed_FRib_run))/(max(smoothed_FRib_run)-min(smoothed_FRib_run)); %normalize between 0 and 1
        smoothed_norm_FRib_run(isnan(smoothed_norm_FRib_run)) = 0; %eliminate nan
        FR_map_orig_running{VR_ephys,1}(i,:) = smoothed_norm_FRib_run;  
        %Calculate spatial information
        prob_in_bin_run = time_in_bin_run./sum(time_in_bin_run); 
        prob_in_bin_run = prob_in_bin_run(1:(end-4)); 
        prob_in_bin_run_store{VR_ephys,1} = prob_in_bin_run; 
        s_lambda_run = nansum(smoothed_norm_FRib_run'.*prob_in_bin_run); 
        s_lambda_i_run = smoothed_norm_FRib_run'; 
        spatial_info_running{VR_ephys,1}(i) = nansum(-prob_in_bin_run.*s_lambda_i_run.*log2(s_lambda_run./s_lambda_i_run));
        s_lambda_squared_run = nansum(((smoothed_norm_FRib_run').^2).*prob_in_bin_run);
        sparsity{VR_ephys,1}(i) = (s_lambda_run.^2)/s_lambda_squared_run;
    end
end
%Generate spatial information of shuffled spike trains during running
for VR_ephys = 1:size(VR_reformatting,1) 
    spatial_info_running_shuffled{VR_ephys,1} = zeros(size(ks_all_data{VR_reformatting{VR_ephys,5},1},1),num_shuffles); 
    spike_pool = cell2mat(in_interval_st_run_store{VR_ephys,1}(:)');
    num_spike_pool = numel(spike_pool);
    for i = 1:size(ks_all_data{VR_reformatting{VR_ephys,5},1},1)
        for t = 1:num_shuffles   
            st_run_shuffled = spike_pool(randi(num_spike_pool,1,numel(in_interval_st_run_store{VR_ephys,1}{i})));
            [spikes_in_bin_run_shuffled,~] = histcounts(VR_reformatting{VR_ephys,4}(st_run_shuffled),spatial_bin_boundaries); 
            FRib_run_shuff = 1000*spikes_in_bin_run_shuffled(:)./time_in_bin_run(:); 
            if strcmp(conv_mode,'gauss') == 1
                smoothed_FRib_run_shuff = conv(gauss_conv,FRib_run_shuff); %convolution
                smoothed_FRib_run_shuff = smoothed_FRib_run_shuff(numel(gauss_conv):(end-numel(gauss_conv)+1)); %eliminate points of partial overlap
            elseif strcmp(conv_mode,'square') == 1
                smoothed_FRib_run_shuff = conv(ones(1,5),FRib_run_shuff); %convolution
                smoothed_FRib_run_shuff = smoothed_FRib_run_shuff(5:(end-4)); 
            end
            smoothed_norm_FRib_run_shuff = (smoothed_FRib_run_shuff-min(smoothed_FRib_run_shuff))/(max(smoothed_FRib_run_shuff)-min(smoothed_FRib_run_shuff));            
            smoothed_norm_FRib_run_shuff(isnan(smoothed_FRib_run_shuff)) = 0;
            lambda_run_shuffled = sum(smoothed_norm_FRib_run_shuff'.*prob_in_bin_run_store{VR_ephys}); 
            lambda_i_run_shuffled = smoothed_norm_FRib_run_shuff'; 
            spatial_info_running_shuffled{VR_ephys,1}(i,t) = nansum(-prob_in_bin_run_store{VR_ephys}.*lambda_i_run_shuffled.*log2(lambda_run_shuffled./lambda_i_run_shuffled));
        end   
    end
    VR_ephys
end    
%% Plot single lap firing patterns and spatial information
marker_distance = 50; %in cm
counter = 0; 
f = figure; set(gcf,'Position',[0 0 3000 1500]); 
ylim_list = [40 80;35 75;60 100;60 100];
g = figure; set(gcf,'Position',[0 0 3000 1500]); 
for VR_ephys = [25,1,22,15] %WT dCA1, AD dCA1, WT vCA1, AD vCA1
%     [~,coi_max] = max(spatial_info_running{VR_ephys,1});
    [~,coi_med] = min(abs(spatial_info_running{VR_ephys,1}-median(spatial_info_running{VR_ephys,1})));        
    counter = counter+1; 
    figure(f); subplot(2,4,counter+4);
    plot(VR_reformatting{VR_ephys,4}*track_length/max(VR_reformatting{VR_ephys,4}),VR_reformatting{VR_ephys,6}/60000,'Color',[.6 .6 .6]); hold on;
    all_st = ks_all_data{VR_reformatting{VR_ephys,5},1}{coi_med,3}; %all spike times (in ms) for a particular unit
    in_interval_st = all_st(and(all_st >= VR_reformatting{VR_ephys,6}(1),all_st <= VR_reformatting{VR_ephys,6}(end))); %only spike times that fall within bounds of current ephys session
    in_interval_st = in_interval_st-VR_reformatting{VR_ephys,6}(1)+1; %correct for the fact that ephys session might not start at first milisecond.                       
    spatial_info_mod_vel2 = VR_reformatting{VR_ephys,10} >= 5; 
    st_stat_or_run = spatial_info_mod_vel2(in_interval_st); %specifies whether each spike in in_interval_st occured during stationary (returns 0) or running (returns 1)
    in_interval_st_run = in_interval_st(st_stat_or_run); %returns only the spikes that occur when running
    scatter(VR_reformatting{VR_ephys,4}(in_interval_st_run)*track_length/max(VR_reformatting{VR_ephys,4}),in_interval_st_run/60000,4,'k','filled');      
    ylabel('Time (min)'); xlabel('Position on track (cm)'); xlim([0 track_length]); title('Median spatial information neuron'); ylim(ylim_list(counter,:)); 
    figure(g); subplot(3,4,counter+4); 
    tmpA = VR_reformatting{VR_ephys,4}(in_interval_st_run)*track_length/max(VR_reformatting{VR_ephys,4});
    tmpA(or(in_interval_st_run/60000 < ylim_list(counter,1),in_interval_st_run/60000 > ylim_list(counter,2))) = [];
    [spikes_in_bin,~] = histcounts(tmpA,[0:10:(10+track_length)],'normalization','count'); 
    tmpB = VR_reformatting{VR_ephys,4}(spatial_info_mod_vel2)*track_length/max(VR_reformatting{VR_ephys,4}); 
    tmpB(or(VR_reformatting{VR_ephys,6}(spatial_info_mod_vel2)/60000 < ylim_list(counter,1),VR_reformatting{VR_ephys,6}(spatial_info_mod_vel2)/60000 > ylim_list(counter,2))) = [];
    [time_in_bin,bin_bounds] = histcounts(tmpB,[0:10:(10+track_length)],'normalization','count');     
    barx = mean([bin_bounds(1:end-1);bin_bounds(2:end)],1);
    norm_FR = spikes_in_bin./time_in_bin; norm_FR = norm_FR/sum(norm_FR); 
    bar(barx,norm_FR,1,'FaceColor','k','EdgeColor','None'); 
    xlabel('Position on track (cm)'); ylabel('Normalized spike probability'); xlim([0 190]); %ylim([0 .01]);
end
figure(g);  
subplot(3,4,[9:10]); 
histogram(cell2mat(spatial_info_running(WT_dCA1_sites)),[0:.001:.44],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(cell2mat(spatial_info_running(AD_dCA1_sites)),[0:.001:.44],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r');
histogram(cell2mat(spatial_info_running_shuffled(WT_dCA1_sites)),[0:.001:.44],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k','LineStyle','--'); hold on;
histogram(cell2mat(spatial_info_running_shuffled(AD_dCA1_sites)),[0:.001:.44],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r','LineStyle','--');
[p,~] = ranksum(cell2mat(spatial_info_running(WT_dCA1_sites)),cell2mat(spatial_info_running(AD_dCA1_sites))); text(.2,.6,strcat('p = ',num2str(p))); 
tmpA = cell2mat(spatial_info_running_shuffled(WT_dCA1_sites)); [p,~] = ranksum(cell2mat(spatial_info_running(WT_dCA1_sites)),tmpA(:)); text(.2,.5,strcat('p = ',num2str(p))); 
tmpA = cell2mat(spatial_info_running_shuffled(AD_dCA1_sites)); [p,~] = ranksum(cell2mat(spatial_info_running(AD_dCA1_sites)),tmpA(:)); text(.2,.4,strcat('p = ',num2str(p))); 
xlabel('Spatial information (bits/s)'); ylabel('Cumulative Probability'); xlim([.05 .37]); ylim([0 1]); 
subplot(3,4,[11:12]); 
histogram(cell2mat(spatial_info_running(WT_vCA1_sites)),[0:.001:.44],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(cell2mat(spatial_info_running(AD_vCA1_sites)),[0:.001:.44],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); 
histogram(cell2mat(spatial_info_running_shuffled(WT_vCA1_sites)),[0:.001:.44],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8],'LineStyle','--'); hold on;
histogram(cell2mat(spatial_info_running_shuffled(AD_vCA1_sites)),[0:.001:.44],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8],'LineStyle','--');
[p,~] = ranksum(cell2mat(spatial_info_running(WT_vCA1_sites)),cell2mat(spatial_info_running(AD_vCA1_sites))); text(.2,.6,strcat('p = ',num2str(p))); 
tmpA = cell2mat(spatial_info_running_shuffled(WT_vCA1_sites)); [p,~] = ranksum(cell2mat(spatial_info_running(WT_vCA1_sites)),tmpA(:)); text(.2,.5,strcat('p = ',num2str(p))); 
tmpA = cell2mat(spatial_info_running_shuffled(AD_vCA1_sites)); [p,~] = ranksum(cell2mat(spatial_info_running(AD_vCA1_sites)),tmpA(:)); text(.2,.4,strcat('p = ',num2str(p))); 
xlabel('Spatial information (bits/s)'); ylabel('Cumulative Probability'); xlim([.05 .4]); ylim([0 1]); 
figure(f); cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 4');
saveas(f,'single_lap','fig'); figure(f); print('single_lap.pdf','-bestfit','-painters','-dpdf'); close(f); 
figure(g); cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 4');
saveas(g,'SI','fig'); print('SI.pdf','-bestfit','-painters','-dpdf'); close(g); 
%% Plot Plaque burden vs. FR
%Combine different recording sessions from same region/animal
%All cells
f = figure; set(gcf,'Position',[0 0 1500 1000]); 
dCA1_FR_consolidated(1) = mean([FR{1};FR{3}]); vCA1_FR_consolidated(1) = mean(FR{2});
dCA1_FR_consolidated(2) = mean([FR{4};FR{6}]); vCA1_FR_consolidated(2) = mean([FR{5};FR{7}]);
dCA1_FR_consolidated(3) = mean(FR{8}); vCA1_FR_consolidated(3) = mean(FR{9});
dCA1_FR_consolidated(4) = mean(FR{11}); vCA1_FR_consolidated(4) = mean([FR{10};FR{12}]); 
dCA1_FR_consolidated(5) = mean(FR{14}); vCA1_FR_consolidated(5) = mean(FR{13}); 
dCA1_FR_consolidated(6) = NaN; vCA1_FR_consolidated(6) = mean(FR{15}); 
model_IFd = fitlm(dCA1_density(1:5)*(1.6^2)*(1000)^2,dCA1_FR_consolidated(1:5)'); 
model_IFv = fitlm(vCA1_density(1:6)*(1.6^2)*(1000)^2,vCA1_FR_consolidated(1:6)'); 
newX = [0:70]'; 
[IFd_newY,IFd_CI] = predict(model_IFd,newX); 
[IFv_newY,IFv_CI] = predict(model_IFv,newX); 
subplot(2,3,1); scatter(dCA1_density(1:5)*(1.6^2)*(1000)^2,dCA1_FR_consolidated(1:5),150,'r','filled'); xlabel('Plaque density'); ylabel('Average FR (Hz)');
[R,P] = corr(dCA1_density(1:5)*(1.6^2)*(1000)^2,dCA1_FR_consolidated(1:5)','Type','Spearman'); 
text(1e-5,7,strcat('p =',num2str(P(1,1)),'r=',num2str(R(1,1)))); xlim([0 70]); ylim([0 15]); hold on;
plot(newX,IFd_newY,'k','LineWidth',2); plot(newX,IFd_CI(:,1),'Color','k'); plot(newX,IFd_CI(:,2),'Color','k'); 
subplot(2,3,4); scatter(vCA1_density(1:6)*(1.6^2)*(1000)^2,vCA1_FR_consolidated(1:6)',150,[1 .8 .8],'filled'); xlabel('Plaque density'); ylabel('Average FR (Hz)');
[R,P] = corr(vCA1_density(1:6)*(1.6^2)*(1000)^2,vCA1_FR_consolidated(1:6)','Type','Spearman'); 
text(1e-5,7,strcat('p =',num2str(P(1,1)),'r=',num2str(R(1,1)))); xlim([0 70]); ylim([0 15]); hold on;
plot(newX,IFv_newY,'k','LineWidth',2); plot(newX,IFv_CI(:,1),'Color','k'); plot(newX,IFv_CI(:,2),'Color','k'); 
%Pyramidal cells only
dCA1_FR_pyramidal_consolidated(1) = mean([FR_pyramidal{1};FR_pyramidal{3}]); vCA1_FR_pyramidal_consolidated(1) = mean(FR_pyramidal{2});
dCA1_FR_pyramidal_consolidated(2) = mean([FR_pyramidal{4};FR_pyramidal{6}]); vCA1_FR_pyramidal_consolidated(2) = mean([FR_pyramidal{5};FR_pyramidal{7}]);
dCA1_FR_pyramidal_consolidated(3) = mean(FR_pyramidal{8}); vCA1_FR_pyramidal_consolidated(3) = mean(FR_pyramidal{9});
dCA1_FR_pyramidal_consolidated(4) = mean(FR_pyramidal{11}); vCA1_FR_pyramidal_consolidated(4) = mean([FR_pyramidal{10};FR_pyramidal{12}]); 
dCA1_FR_pyramidal_consolidated(5) = mean(FR_pyramidal{14}); vCA1_FR_pyramidal_consolidated(5) = mean(FR_pyramidal{13}); 
dCA1_FR_pyramidal_consolidated(6) = NaN; vCA1_FR_pyramidal_consolidated(6) = mean(FR_pyramidal{15}); 
model_IFd = fitlm(dCA1_density(1:5)*(1.6^2)*(1000)^2,dCA1_FR_pyramidal_consolidated(1:5)'); 
model_IFv = fitlm(vCA1_density(1:6)*(1.6^2)*(1000)^2,vCA1_FR_pyramidal_consolidated(1:6)'); 
newX = [0:70]'; 
[IFd_newY,IFd_CI] = predict(model_IFd,newX); 
[IFv_newY,IFv_CI] = predict(model_IFv,newX); 
subplot(2,3,2); scatter(dCA1_density(1:5)*(1.6^2)*(1000)^2,dCA1_FR_pyramidal_consolidated(1:5),150,'r','filled'); xlabel('Plaque density'); ylabel('Average FR (Hz)');
[R,P] = corr(dCA1_density(1:5)*(1.6^2)*(1000)^2,dCA1_FR_pyramidal_consolidated(1:5)','Type','Spearman'); 
text(1e-5,7,strcat('p =',num2str(P(1,1)),'r=',num2str(R(1,1)))); xlim([0 70]); ylim([0 15]); hold on;
plot(newX,IFd_newY,'k','LineWidth',2); plot(newX,IFd_CI(:,1),'Color','k'); plot(newX,IFd_CI(:,2),'Color','k'); 
subplot(2,3,5); scatter(vCA1_density(1:6)*(1.6^2)*(1000)^2,vCA1_FR_pyramidal_consolidated(1:6)',150,[1 .8 .8],'filled'); xlabel('Plaque density'); ylabel('Average FR (Hz)');
[R,P] = corr(vCA1_density(1:6)*(1.6^2)*(1000)^2,vCA1_FR_pyramidal_consolidated(1:6)','Type','Spearman'); 
text(1e-5,7,strcat('p =',num2str(P(1,1)),'r=',num2str(R(1,1)))); xlim([0 70]); ylim([0 15]); hold on;
plot(newX,IFv_newY,'k','LineWidth',2); plot(newX,IFv_CI(:,1),'Color','k'); plot(newX,IFv_CI(:,2),'Color','k'); 
%Interneurons only
dCA1_FR_interneuron_consolidated(1) = mean([FR_interneuron{1};FR_interneuron{3}]); vCA1_FR_interneuron_consolidated(1) = mean(FR_interneuron{2});
dCA1_FR_interneuron_consolidated(2) = mean([FR_interneuron{4};FR_interneuron{6}]); vCA1_FR_interneuron_consolidated(2) = mean([FR_interneuron{5};FR_interneuron{7}]);
dCA1_FR_interneuron_consolidated(3) = mean(FR_interneuron{8}); vCA1_FR_interneuron_consolidated(3) = mean(FR_interneuron{9});
dCA1_FR_interneuron_consolidated(4) = mean(FR_interneuron{11}); vCA1_FR_interneuron_consolidated(4) = mean([FR_interneuron{10};FR_interneuron{12}]); 
dCA1_FR_interneuron_consolidated(5) = mean(FR_interneuron{14}); vCA1_FR_interneuron_consolidated(5) = mean(FR_interneuron{13}); 
dCA1_FR_interneuron_consolidated(6) = NaN; vCA1_FR_interneuron_consolidated(6) = mean(FR_interneuron{15}); 
model_IFd = fitlm(dCA1_density(1:5)*(1.6^2)*(1000)^2,dCA1_FR_interneuron_consolidated(1:5)'); 
model_IFv = fitlm(vCA1_density(1:6)*(1.6^2)*(1000)^2,vCA1_FR_interneuron_consolidated(1:6)'); 
newX = [0:70]'; 
[IFd_newY,IFd_CI] = predict(model_IFd,newX); 
[IFv_newY,IFv_CI] = predict(model_IFv,newX); 
subplot(2,3,3); scatter(dCA1_density(1:5)*(1.6^2)*(1000)^2,dCA1_FR_interneuron_consolidated(1:5),150,'r','filled'); xlabel('Plaque density'); ylabel('Average FR (Hz)');
[R,P] = corr(dCA1_density(1:5)*(1.6^2)*(1000)^2,dCA1_FR_interneuron_consolidated(1:5)','Type','Spearman'); 
text(1e-5,7,strcat('p =',num2str(P(1,1)),'r=',num2str(R(1,1)))); xlim([0 70]); ylim([0 15]); hold on;
plot(newX,IFd_newY,'k','LineWidth',2); plot(newX,IFd_CI(:,1),'Color','k'); plot(newX,IFd_CI(:,2),'Color','k'); 
subplot(2,3,6); scatter(vCA1_density(1:6)*(1.6^2)*(1000)^2,vCA1_FR_interneuron_consolidated(1:6)',150,[1 .8 .8],'filled'); xlabel('Plaque density'); ylabel('Average FR (Hz)');
[R,P] = corr(vCA1_density(1:6)*(1.6^2)*(1000)^2,vCA1_FR_interneuron_consolidated(1:6)','Type','Spearman'); 
text(1e-5,7,strcat('p =',num2str(P(1,1)),'r=',num2str(R(1,1)))); xlim([0 70]); ylim([0 15]); hold on;
plot(newX,IFv_newY,'k','LineWidth',2); plot(newX,IFv_CI(:,1),'Color','k'); plot(newX,IFv_CI(:,2),'Color','k'); 
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 3'); 
saveas(f,'plaque_burden_FR','fig'); print('plaque_burden_FR.pdf','-bestfit','-painters','-dpdf'); close(f);
%% Plot consistent correlations
wl = 16; bin_count = 2;
f = figure; set(f,'Position',[0 0 1500 1500]); 
num_subsamples = 100; 
% All neurons
subsampled_corr = cell(size(corr_consistent,1),1); 
subsampled_corr_allbins = cell(size(corr_consistent,1),numel(all_bin_sizes)); 
for VR_ephys = 1:size(corr_consistent,1)                          
    if numel(corr_consistent{VR_ephys,wl,bin_count}) >= num_subsamples
        rand_idx = randperm(numel(corr_consistent{VR_ephys,wl,bin_count}),num_subsamples);
        subsampled_corr{VR_ephys} = corr_consistent{VR_ephys,wl,bin_count}(rand_idx); 
        for bin = 1:numel(all_bin_sizes)
            subsampled_corr_allbins{VR_ephys,bin} = corr_consistent{VR_ephys,wl,bin}(rand_idx)'; 
        end
    end
end
subplot(3,3,1);
histogram(cell2mat(subsampled_corr(WT_dCA1_sites)),[-.2:.001:.45],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(cell2mat(subsampled_corr(AD_dCA1_sites)),[-.2:.001:.45],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
[p,~] = ranksum(cell2mat(subsampled_corr(WT_dCA1_sites)'),cell2mat(subsampled_corr(AD_dCA1_sites)')); text(.2,.7,strcat('p = ',num2str(p))); 
xlabel('Correlation'); ylabel('Cumulative Probability'); xlim([-.05 .4]); ylim([0 1]); axis square
subplot(3,3,4);
histogram(cell2mat(subsampled_corr(WT_vCA1_sites)),[-.2:.001:.45],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(cell2mat(subsampled_corr(AD_vCA1_sites)),[-.2:.001:.45],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); 
[p,~] = ranksum(cell2mat(subsampled_corr(WT_vCA1_sites)'),cell2mat(subsampled_corr(AD_vCA1_sites)')); text(.2,.6,strcat('p = ',num2str(p))); 
xlabel('Correlation'); ylabel('Cumulative Probability'); xlim([-.05 .2]); ylim([0 1]); axis square
subplot(3,3,7); 
errorbar(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(WT_dCA1_sites,2:end))),std(cell2mat(subsampled_corr_allbins(WT_dCA1_sites,2:end)),0,1)./sqrt(size(cell2mat(subsampled_corr_allbins(WT_dCA1_sites,2:end)),1)),'Color','k','LineWidth',2); hold on; 
scatter(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(WT_dCA1_sites,2:end))),'ko','filled');
errorbar(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(WT_vCA1_sites,2:end))),std(cell2mat(subsampled_corr_allbins(WT_vCA1_sites,2:end)),0,1)./sqrt(size(cell2mat(subsampled_corr_allbins(WT_vCA1_sites,2:end)),1)),'Color',[.8 .8 .8],'LineWidth',2); hold on; 
scatter(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(WT_vCA1_sites,2:end))),'filled','o','MarkerEdgeColor',[.8 .8 .8],'MarkerFaceColor',[.8 .8 .8])
errorbar(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(AD_dCA1_sites,2:end))),std(cell2mat(subsampled_corr_allbins(AD_dCA1_sites,2:end)),0,1)./sqrt(size(cell2mat(subsampled_corr_allbins(AD_dCA1_sites,2:end)),1)),'Color','r','LineWidth',2); hold on; 
scatter(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(AD_dCA1_sites,2:end))),'ro','filled')
errorbar(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(AD_vCA1_sites,2:end))),std(cell2mat(subsampled_corr_allbins(AD_vCA1_sites,2:end)),0,1)./sqrt(size(cell2mat(subsampled_corr_allbins(AD_vCA1_sites,2:end)),1)),'Color',[1 .8 .8],'LineWidth',2); hold on; 
scatter(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(AD_vCA1_sites,2:end))),'filled','o','MarkerEdgeColor',[1 .8 .8],'MarkerFaceColor',[1 .8 .8])
xlim([0 110]); ylim([0 .22]); xlabel('Bin size (ms)'); ylabel('Correlation'); axis square;
ypos = .16; idx = 0;
for i = 2:numel(all_bin_sizes)
    [p,~] = ranksum(cell2mat(subsampled_corr_allbins(WT_dCA1_sites,i)),cell2mat(subsampled_corr_allbins(AD_dCA1_sites,i))); text(110,ypos,strcat('p = ',num2str(p))); ypos = ypos - .01;
    [p,~] = ranksum(cell2mat(subsampled_corr_allbins(WT_vCA1_sites,i)),cell2mat(subsampled_corr_allbins(AD_vCA1_sites,i))); text(110,ypos,strcat('p = ',num2str(p))); ypos = ypos - .01;
end
% Pyramidal neurons
num_subsamples = 100; 
subsampled_corr = cell(size(corr_consistent_pyramidal,1),1); 
subsampled_corr_allbins = cell(size(corr_consistent_pyramidal,1),numel(all_bin_sizes)); 
for VR_ephys = 1:size(corr_consistent_pyramidal,1)                          
    if numel(corr_consistent_pyramidal{VR_ephys,wl,bin_count}) >= num_subsamples
        rand_idx = randperm(numel(corr_consistent_pyramidal{VR_ephys,wl,bin_count}),num_subsamples);
        subsampled_corr{VR_ephys} = corr_consistent_pyramidal{VR_ephys,wl,bin_count}(rand_idx); 
        for bin = 1:numel(all_bin_sizes)
            subsampled_corr_allbins{VR_ephys,bin} = corr_consistent_pyramidal{VR_ephys,wl,bin}(rand_idx)'; 
        end
    end
end
subplot(3,3,2);
histogram(cell2mat(subsampled_corr(WT_dCA1_sites)),[-.2:.001:.45],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(cell2mat(subsampled_corr(AD_dCA1_sites)),[-.2:.001:.45],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
[p,~] = ranksum(cell2mat(subsampled_corr(WT_dCA1_sites)'),cell2mat(subsampled_corr(AD_dCA1_sites)')); text(.2,.7,strcat('p = ',num2str(p))); 
xlabel('Correlation'); ylabel('Cumulative Probability'); xlim([-.05 .4]); ylim([0 1]); axis square
subplot(3,3,5);
histogram(cell2mat(subsampled_corr(WT_vCA1_sites)),[-.2:.001:.45],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(cell2mat(subsampled_corr(AD_vCA1_sites)),[-.2:.001:.45],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); 
[p,~] = ranksum(cell2mat(subsampled_corr(WT_vCA1_sites)'),cell2mat(subsampled_corr(AD_vCA1_sites)')); text(.2,.6,strcat('p = ',num2str(p))); 
xlabel('Correlation'); ylabel('Cumulative Probability'); xlim([-.05 .2]); ylim([0 1]); axis square
subplot(3,3,8); 
errorbar(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(WT_dCA1_sites,2:end))),std(cell2mat(subsampled_corr_allbins(WT_dCA1_sites,2:end)),0,1)./sqrt(size(cell2mat(subsampled_corr_allbins(WT_dCA1_sites,2:end)),1)),'Color','k','LineWidth',2); hold on; 
scatter(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(WT_dCA1_sites,2:end))),'ko','filled');
errorbar(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(WT_vCA1_sites,2:end))),std(cell2mat(subsampled_corr_allbins(WT_vCA1_sites,2:end)),0,1)./sqrt(size(cell2mat(subsampled_corr_allbins(WT_vCA1_sites,2:end)),1)),'Color',[.8 .8 .8],'LineWidth',2); hold on; 
scatter(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(WT_vCA1_sites,2:end))),'filled','o','MarkerEdgeColor',[.8 .8 .8],'MarkerFaceColor',[.8 .8 .8])
errorbar(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(AD_dCA1_sites,2:end))),std(cell2mat(subsampled_corr_allbins(AD_dCA1_sites,2:end)),0,1)./sqrt(size(cell2mat(subsampled_corr_allbins(AD_dCA1_sites,2:end)),1)),'Color','r','LineWidth',2); hold on; 
scatter(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(AD_dCA1_sites,2:end))),'ro','filled')
errorbar(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(AD_vCA1_sites,2:end))),std(cell2mat(subsampled_corr_allbins(AD_vCA1_sites,2:end)),0,1)./sqrt(size(cell2mat(subsampled_corr_allbins(AD_vCA1_sites,2:end)),1)),'Color',[1 .8 .8],'LineWidth',2); hold on; 
scatter(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(AD_vCA1_sites,2:end))),'filled','o','MarkerEdgeColor',[1 .8 .8],'MarkerFaceColor',[1 .8 .8])
xlim([0 110]); ylim([0 .22]); xlabel('Bin size (ms)'); ylabel('Correlation'); axis square;
ypos = .16; idx = 0;
for i = 2:numel(all_bin_sizes)
    [p,~] = ranksum(cell2mat(subsampled_corr_allbins(WT_dCA1_sites,i)),cell2mat(subsampled_corr_allbins(AD_dCA1_sites,i))); text(110,ypos,strcat('p = ',num2str(p))); ypos = ypos - .01;
    [p,~] = ranksum(cell2mat(subsampled_corr_allbins(WT_vCA1_sites,i)),cell2mat(subsampled_corr_allbins(AD_vCA1_sites,i))); text(110,ypos,strcat('p = ',num2str(p))); ypos = ypos - .01;
end
% Interneurons
num_subsamples = 100; 
subsampled_corr = cell(size(corr_consistent_interneuron,1),1); 
subsampled_corr_allbins = cell(size(corr_consistent_interneuron,1),numel(all_bin_sizes)); 
for VR_ephys = 1:size(corr_consistent_interneuron,1)                          
    if numel(corr_consistent_interneuron{VR_ephys,wl,bin_count}) >= num_subsamples
        rand_idx = randperm(numel(corr_consistent_interneuron{VR_ephys,wl,bin_count}),num_subsamples);
        subsampled_corr{VR_ephys} = corr_consistent_interneuron{VR_ephys,wl,bin_count}(rand_idx); 
        for bin = 1:numel(all_bin_sizes)
            subsampled_corr_allbins{VR_ephys,bin} = corr_consistent_interneuron{VR_ephys,wl,bin}(rand_idx)'; 
        end
    end
end
subplot(3,3,3);
histogram(cell2mat(subsampled_corr(WT_dCA1_sites)),[-.2:.001:.45],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(cell2mat(subsampled_corr(AD_dCA1_sites)),[-.2:.001:.45],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
[p,~] = ranksum(cell2mat(subsampled_corr(WT_dCA1_sites)'),cell2mat(subsampled_corr(AD_dCA1_sites)')); text(.2,.7,strcat('p = ',num2str(p))); 
xlabel('Correlation'); ylabel('Cumulative Probability'); xlim([-.05 .4]); ylim([0 1]); axis square
subplot(3,3,6);
histogram(cell2mat(subsampled_corr(WT_vCA1_sites)),[-.2:.001:.45],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(cell2mat(subsampled_corr(AD_vCA1_sites)),[-.2:.001:.45],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); 
[p,~] = ranksum(cell2mat(subsampled_corr(WT_vCA1_sites)'),cell2mat(subsampled_corr(AD_vCA1_sites)')); text(.2,.6,strcat('p = ',num2str(p))); 
xlabel('Correlation'); ylabel('Cumulative Probability'); xlim([-.05 .2]); ylim([0 1]); axis square
subplot(3,3,9); 
errorbar(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(WT_dCA1_sites,2:end))),std(cell2mat(subsampled_corr_allbins(WT_dCA1_sites,2:end)),0,1)./sqrt(size(cell2mat(subsampled_corr_allbins(WT_dCA1_sites,2:end)),1)),'Color','k','LineWidth',2); hold on; 
scatter(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(WT_dCA1_sites,2:end))),'ko','filled');
errorbar(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(WT_vCA1_sites,2:end))),std(cell2mat(subsampled_corr_allbins(WT_vCA1_sites,2:end)),0,1)./sqrt(size(cell2mat(subsampled_corr_allbins(WT_vCA1_sites,2:end)),1)),'Color',[.8 .8 .8],'LineWidth',2); hold on; 
scatter(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(WT_vCA1_sites,2:end))),'filled','o','MarkerEdgeColor',[.8 .8 .8],'MarkerFaceColor',[.8 .8 .8])
errorbar(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(AD_dCA1_sites,2:end))),std(cell2mat(subsampled_corr_allbins(AD_dCA1_sites,2:end)),0,1)./sqrt(size(cell2mat(subsampled_corr_allbins(AD_dCA1_sites,2:end)),1)),'Color','r','LineWidth',2); hold on; 
scatter(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(AD_dCA1_sites,2:end))),'ro','filled')
errorbar(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(AD_vCA1_sites,2:end))),std(cell2mat(subsampled_corr_allbins(AD_vCA1_sites,2:end)),0,1)./sqrt(size(cell2mat(subsampled_corr_allbins(AD_vCA1_sites,2:end)),1)),'Color',[1 .8 .8],'LineWidth',2); hold on; 
scatter(all_bin_sizes(2:end),mean(cell2mat(subsampled_corr_allbins(AD_vCA1_sites,2:end))),'filled','o','MarkerEdgeColor',[1 .8 .8],'MarkerFaceColor',[1 .8 .8])
xlim([0 110]); ylim([0 .22]); xlabel('Bin size (ms)'); ylabel('Correlation'); axis square;
ypos = .16; idx = 0;
for i = 2:numel(all_bin_sizes)
    [p,~] = ranksum(cell2mat(subsampled_corr_allbins(WT_dCA1_sites,i)),cell2mat(subsampled_corr_allbins(AD_dCA1_sites,i))); text(110,ypos,strcat('p = ',num2str(p))); ypos = ypos - .01;
    [p,~] = ranksum(cell2mat(subsampled_corr_allbins(WT_vCA1_sites,i)),cell2mat(subsampled_corr_allbins(AD_vCA1_sites,i))); text(110,ypos,strcat('p = ',num2str(p))); ypos = ypos - .01;
end
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 5'); 
saveas(f,'Consistent correlations','fig'); orient(f,'landscape'); print('Consistent_correlations.pdf','-bestfit','-painters','-dpdf'); close(f); 
%% Figures: Correlation matrices
% CMap = [linspace(232/255,1,41)',linspace(197/255,1,41)',linspace(240/255,1,41)'];
CMap = [linspace(209/255,1,41)',linspace(139/255,1,41)',linspace(225/255,1,41)'];
CMap(end,:) = []; 
CMap = [CMap;linspace(1,50/255,123)',linspace(1,205/255,123)',linspace(1,50/255,123)'];     
cnames = {'WT dCA1','WT vCA1','AD dCA1','AD vCA1'}; n = 0;
%all neurons
for VR_ephys = [25,27,11,15] 
    n = n+1; 
    tmpC = squareform(corr_consistent{VR_ephys,wl,2}); tmpC = tmpC-diag(diag(tmpC));
    if VR_ephys == 17
        tmpC = tmpC(55:end,55:end); 
    end
    f = figure; imagesc(tmpC,[-.1 .3]); colormap(CMap);     
    title(ks_all_data{VR_ephys,5}); colorbar; title(cnames{n}); xlabel('Units'); ylabel('Units'); axis square; 
    cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 5'); 
    saveas(f,strcat(cnames{n},'corrmat'),'fig'); print(strcat(cnames{n},'corrmat.pdf'),'-bestfit','-painters','-dpdf'); close(f);     
end
%pyramidals only
n = 0; 
for VR_ephys = [25,27,11,15] 
    n = n+1; 
    tmpC = squareform(corr_consistent_pyramidal{VR_ephys,wl,2}); tmpC = tmpC-diag(diag(tmpC));
    if VR_ephys == 17
        tmpC = tmpC(55:end,55:end); 
    end
    f = figure; imagesc(tmpC,[-.1 .3]); colormap(CMap);     
    title(ks_all_data{VR_ephys,5}); colorbar; title(cnames{n}); xlabel('Units'); ylabel('Units'); axis square; 
    cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 5'); 
    saveas(f,strcat(cnames{n},'corrmat_pyramidal'),'fig'); print(strcat(cnames{n},'corrmat_pyramidal.pdf'),'-bestfit','-painters','-dpdf'); close(f); 
end
%% Calculate network properties
threshold_levels = 0:.001:.1;
word_length = 16; 
trunc_network_size = 16;
num_subsamples = 100;
range_network_sizes = 3:20;
degree_subsampled = cell(numel(threshold_levels),size(mod_rasters,1)); 
clustering_subsampled = degree_subsampled;
degree_subsampled_pyramidal = degree_subsampled;
clustering_subsampled_pyramidal = degree_subsampled; 
for VR_ephys = 1:size(VR_reformatting,1)
    if size(mod_rasters{VR_ephys},1) >= word_length
        idx = 0;
        for network_threshold = threshold_levels
            idx = idx+1;
            bin_corr_mat = squareform(corr_consistent{VR_ephys,16}); bin_corr_mat(bin_corr_mat < network_threshold) = 0; bin_corr_mat(bin_corr_mat >= network_threshold) = 1;
            bin_corr_mat_pyramidal = squareform(corr_consistent_pyramidal{VR_ephys,16}); bin_corr_mat(bin_corr_mat < network_threshold) = 0; bin_corr_mat(bin_corr_mat >= network_threshold) = 1;
            for subsample = 1:num_subsamples
                preserved_nodes = randperm(size(bin_corr_mat,1),trunc_network_size);
                bin_corr_mat_trunc = bin_corr_mat(preserved_nodes,preserved_nodes); 
                degree_subsampled{idx,VR_ephys}(subsample,1) = mean2(bin_corr_mat_trunc); 
                clustering_subsampled{idx,VR_ephys}(subsample,1) = mean(clustering_coef_bu(bin_corr_mat_trunc)); 
                preserved_nodes_pyramidal = randperm(size(bin_corr_mat_pyramidal,1),trunc_network_size);
                bin_corr_mat_trunc_pyramidal = bin_corr_mat(preserved_nodes_pyramidal,preserved_nodes_pyramidal); 
                degree_subsampled_pyramidal{idx,VR_ephys}(subsample,1) = mean2(bin_corr_mat_trunc_pyramidal); 
                clustering_subsampled_pyramidal{idx,VR_ephys}(subsample,1) = mean(clustering_coef_bu(bin_corr_mat_trunc_pyramidal)); 
            end
        end
    end
    VR_ephys
end
%% Plot network visualizations
network_threshold = .06;
soi_list = [16,22,1,9];
color_list = [0 0 0;.8 .8 .8;1 0 0;1 .8 .8];
trunc_network_size = 16;
for soi = 1:numel(soi_list)
    VR_ephys = soi_list(soi);
    f = figure; set(f,'Position',[0 200 500 500]); 
    vertices = ks_all_data{VR_ephys,6}; 
    vertices_jittered = vertices+50*randn(size(vertices,1),size(vertices,2));
    bin_corr_mat = squareform(corr_consistent{VR_ephys,16}); bin_corr_mat(bin_corr_mat < network_threshold) = NaN; bin_corr_mat(bin_corr_mat >= network_threshold) = 1;
    preserved_nodes = randperm(size(bin_corr_mat,1),trunc_network_size);
    bin_corr_mat_trunc = bin_corr_mat(preserved_nodes,preserved_nodes); 
    for i = 1:size(bin_corr_mat_trunc,1)
        for j = 1:size(bin_corr_mat_trunc,1)
            if isnan(bin_corr_mat_trunc(i,j)) == 0                    
                plot(vertices_jittered([i,j],1),vertices_jittered([i,j],2),'Color',color_list(soi,:)); hold on;
            end
        end
    end
    for i = 1:size(bin_corr_mat_trunc,1)
        scatter(vertices_jittered(i,1),vertices_jittered(i,2),100,color_list(soi,:),'filled'); hold on; 
    end         
    xlim([-200 600]); ylim([-900 150]); xticks([]); yticks([]); 
    plot([-150 -50],[-850 -850],'k','LineWidth',2); 
end
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 5'); 
saveas(gcf,'Network4','fig'); orient(gcf,'landscape'); print('Network4.pdf','-bestfit','-painters','-dpdf'); 
soi_list = [16,27,1,13];
for soi = 1:numel(soi_list)
    VR_ephys = soi_list(soi);
    f = figure; set(f,'Position',[0 200 500 500]); 
    vertices = ks_all_data{VR_ephys,6}; 
    vertices_jittered = vertices+50*randn(size(vertices,1),size(vertices,2));
    bin_corr_mat = squareform(corr_consistent_pyramidal{VR_ephys,16}); bin_corr_mat(bin_corr_mat < network_threshold) = NaN; bin_corr_mat(bin_corr_mat >= network_threshold) = 1;
    preserved_nodes = randperm(size(bin_corr_mat,1),trunc_network_size);
    bin_corr_mat_trunc = bin_corr_mat(preserved_nodes,preserved_nodes); 
    for i = 1:size(bin_corr_mat_trunc,1)
        for j = 1:size(bin_corr_mat_trunc,1)
            if isnan(bin_corr_mat_trunc(i,j)) == 0                    
                plot(vertices_jittered([i,j],1),vertices_jittered([i,j],2),'Color',color_list(soi,:)); hold on;
            end
        end
    end
    for i = 1:size(bin_corr_mat_trunc,1)
        scatter(vertices_jittered(i,1),vertices_jittered(i,2),100,color_list(soi,:),'filled'); hold on; 
    end         
    xlim([-200 600]); ylim([-900 150]); xticks([]); yticks([]); 
    plot([-150 -50],[-850 -850],'k','LineWidth',2); 
end
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 5'); 
saveas(gcf,'Network4_pyramidal','fig'); orient(gcf,'landscape'); print('Network4_pyramidal.pdf','-bestfit','-painters','-dpdf');  
%% Plot network properties
tioi = 6;%11;
f = figure; set(gcf,'Position',[0 -500 500 1500]); 
subplot(3,2,1); 
plot(threshold_levels,nanmean(cell2mat(degree_subsampled(:,AD_dCA1_sites)')),'r','LineWidth',2); hold on;
plot(threshold_levels,nanmean(cell2mat(degree_subsampled(:,WT_dCA1_sites)')),'k','LineWidth',2); hold on;
plot(threshold_levels,nanmean(cell2mat(degree_subsampled(:,AD_vCA1_sites)')),'Color',[1 .8 .8],'LineWidth',2); hold on;
plot(threshold_levels,nanmean(cell2mat(degree_subsampled(:,WT_vCA1_sites)')),'Color',[.8 .8 .8],'LineWidth',2); hold on;
ylabel('Relative Degree'); xlabel('Threshold'); ylim([.5e-2 1]); xlim([1e-2, .1]); set(gca,'YScale','log'); axis square;
subplot(3,2,2); 
plot(threshold_levels,nanmean(cell2mat(clustering_subsampled(:,AD_dCA1_sites)')),'r','LineWidth',2); hold on;
plot(threshold_levels,nanmean(cell2mat(clustering_subsampled(:,WT_dCA1_sites)')),'k','LineWidth',2); hold on;
plot(threshold_levels,nanmean(cell2mat(clustering_subsampled(:,AD_vCA1_sites)')),'Color',[1 .8 .8],'LineWidth',2); hold on;
plot(threshold_levels,nanmean(cell2mat(clustering_subsampled(:,WT_vCA1_sites)')),'Color',[.8 .8 .8],'LineWidth',2); hold on;
ylabel('Clustering Coefficient'); xlabel('Threshold'); ylim([1e-2 1]); xlim([1e-2, .1]); set(gca,'YScale','log'); axis square;
subplot(3,2,3); 
histogram(cell2mat(degree_subsampled(tioi,WT_dCA1_sites)'),[0:.001:1.2],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(cell2mat(degree_subsampled(tioi,AD_dCA1_sites)'),[0:.001:1.2],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
[p,~] = ranksum(cell2mat(degree_subsampled(tioi,WT_dCA1_sites)'),cell2mat(degree_subsampled(tioi,AD_dCA1_sites)')); text(.6,.7,strcat('p = ',num2str(p))); 
xlabel('Relative degree'); ylabel('Cumulative Probability'); xlim([.3 .95]); ylim([0 1]); axis square
subplot(3,2,4); 
histogram(cell2mat(clustering_subsampled(tioi,WT_dCA1_sites)'),[0:.001:1.2],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(cell2mat(clustering_subsampled(tioi,AD_dCA1_sites)'),[0:.001:1.2],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
[p,~] = ranksum(cell2mat(clustering_subsampled(tioi,WT_dCA1_sites)'),cell2mat(clustering_subsampled(tioi,AD_dCA1_sites)')); text(.6,.7,strcat('p = ',num2str(p))); 
xlabel('Clustering coefficient'); ylabel('Cumulative Probability'); xlim([.5 1]); ylim([0 1]); axis square
subplot(3,2,5); 
histogram(cell2mat(degree_subsampled(tioi,WT_vCA1_sites)'),[0:.001:1.2],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(cell2mat(degree_subsampled(tioi,AD_vCA1_sites)'),[0:.001:1.2],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
[p,~] = ranksum(cell2mat(degree_subsampled(tioi,WT_vCA1_sites)'),cell2mat(degree_subsampled(tioi,AD_vCA1_sites)')); text(.6,.7,strcat('p = ',num2str(p))); 
xlabel('Relative degree'); ylabel('Cumulative Probability'); xlim([.05 1]); ylim([0 1]); axis square
subplot(3,2,6); 
histogram(cell2mat(clustering_subsampled(tioi,WT_vCA1_sites)'),[0:.001:1.2],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(cell2mat(clustering_subsampled(tioi,AD_vCA1_sites)'),[0:.001:1.2],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
[p,~] = ranksum(cell2mat(clustering_subsampled(tioi,WT_vCA1_sites)'),cell2mat(clustering_subsampled(tioi,AD_vCA1_sites)')); text(.6,.7,strcat('p = ',num2str(p))); 
xlabel('Clustering coefficient'); ylabel('Cumulative Probability'); xlim([0 1.05]); ylim([0 1]); axis square
g = figure; set(gcf,'Position',[0 -500 500 1500]); 
subplot(3,2,1); 
plot(threshold_levels,nanmean(cell2mat(degree_subsampled_pyramidal(:,AD_dCA1_sites)')),'r','LineWidth',2); hold on;
plot(threshold_levels,nanmean(cell2mat(degree_subsampled_pyramidal(:,WT_dCA1_sites)')),'k','LineWidth',2); hold on;
plot(threshold_levels,nanmean(cell2mat(degree_subsampled_pyramidal(:,AD_vCA1_sites)')),'Color',[1 .8 .8],'LineWidth',2); hold on;
plot(threshold_levels,nanmean(cell2mat(degree_subsampled_pyramidal(:,WT_vCA1_sites)')),'Color',[.8 .8 .8],'LineWidth',2); hold on;
ylabel('Relative Degree'); xlabel('Threshold'); ylim([.5e-2 1]); xlim([1e-2, .1]); set(gca,'YScale','log'); axis square;
subplot(3,2,2); 
plot(threshold_levels,nanmean(cell2mat(clustering_subsampled_pyramidal(:,AD_dCA1_sites)')),'r','LineWidth',2); hold on;
plot(threshold_levels,nanmean(cell2mat(clustering_subsampled_pyramidal(:,WT_dCA1_sites)')),'k','LineWidth',2); hold on;
plot(threshold_levels,nanmean(cell2mat(clustering_subsampled_pyramidal(:,AD_vCA1_sites)')),'Color',[1 .8 .8],'LineWidth',2); hold on;
plot(threshold_levels,nanmean(cell2mat(clustering_subsampled_pyramidal(:,WT_vCA1_sites)')),'Color',[.8 .8 .8],'LineWidth',2); hold on;
ylabel('Clustering Coefficient'); xlabel('Threshold'); ylim([1e-2 1]); xlim([1e-2, .1]); set(gca,'YScale','log'); axis square;
subplot(3,2,3); 
histogram(cell2mat(degree_subsampled_pyramidal(tioi,WT_dCA1_sites)'),[0:.001:1.2],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(cell2mat(degree_subsampled_pyramidal(tioi,AD_dCA1_sites)'),[0:.001:1.2],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
[p,~] = ranksum(cell2mat(degree_subsampled_pyramidal(tioi,WT_dCA1_sites)'),cell2mat(degree_subsampled_pyramidal(tioi,AD_dCA1_sites)')); text(.6,.7,strcat('p = ',num2str(p))); 
xlabel('Relative degree'); ylabel('Cumulative Probability'); xlim([.3 .95]); ylim([0 1]); axis square
subplot(3,2,4); 
histogram(cell2mat(clustering_subsampled_pyramidal(tioi,WT_dCA1_sites)'),[0:.001:1.2],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(cell2mat(clustering_subsampled_pyramidal(tioi,AD_dCA1_sites)'),[0:.001:1.2],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
[p,~] = ranksum(cell2mat(clustering_subsampled_pyramidal(tioi,WT_dCA1_sites)'),cell2mat(clustering_subsampled_pyramidal(tioi,AD_dCA1_sites)')); text(.6,.7,strcat('p = ',num2str(p))); 
xlabel('Clustering coefficient'); ylabel('Cumulative Probability'); xlim([.5 1]); ylim([0 1]); axis square
subplot(3,2,5); 
histogram(cell2mat(degree_subsampled_pyramidal(tioi,WT_vCA1_sites)'),[0:.001:1.2],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(cell2mat(degree_subsampled_pyramidal(tioi,AD_vCA1_sites)'),[0:.001:1.2],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
[p,~] = ranksum(cell2mat(degree_subsampled_pyramidal(tioi,WT_vCA1_sites)'),cell2mat(degree_subsampled_pyramidal(tioi,AD_vCA1_sites)')); text(.6,.7,strcat('p = ',num2str(p))); 
xlabel('Relative degree'); ylabel('Cumulative Probability'); xlim([.05 1]); ylim([0 1]); axis square
subplot(3,2,6); 
histogram(cell2mat(clustering_subsampled_pyramidal(tioi,WT_vCA1_sites)'),[0:.001:1.2],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(cell2mat(clustering_subsampled_pyramidal(tioi,AD_vCA1_sites)'),[0:.001:1.2],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
[p,~] = ranksum(cell2mat(clustering_subsampled_pyramidal(tioi,WT_vCA1_sites)'),cell2mat(clustering_subsampled_pyramidal(tioi,AD_vCA1_sites)')); text(.6,.7,strcat('p = ',num2str(p))); 
xlabel('Clustering coefficient'); ylabel('Cumulative Probability'); xlim([0 1.05]); ylim([0 1]); axis square
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 5');
saveas(f,'Degree_and_clustering','fig'); orient(f,'landscape'); print('Degree_and_clustering.pdf','-bestfit','-painters','-dpdf'); close(f); 
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 5');
saveas(g,'Degree_and_clustering_pyramidal','fig'); orient(g,'landscape'); print('Degree_and_clustering_pyramidal.pdf','-bestfit','-painters','-dpdf'); close(g); 
%% Entropy and max entropy 
ent_trials = 500;
converter_original = 2.^(0:99); 
all_bin_sizes = [10]; entropy_bin_size = all_bin_sizes; 
all_word_lengths = 16; %[4:2:24];
entropy_consistent = cell(size(VR_reformatting,1),numel(converter_original),numel(all_bin_sizes));
entropy_model_Kpairwise_consistent = entropy_consistent; 
entropy_model_pairwise_consistent = entropy_consistent;
entropy_model_ind_consistent = entropy_consistent;
entropy_model_Kpairwise_sampled_consistent = entropy_consistent; 
entropy_model_pairwise_sampled_consistent = entropy_consistent;
entropy_model_ind_sampled_consistent = entropy_consistent;
params_Kpairwise_consistent = entropy_consistent;
params_pairwise_consistent = entropy_consistent;
params_ind_consistent = entropy_consistent;
convergence_check_Kpairwise = entropy_consistent; 
convergence_check_pairwise = entropy_consistent; 
convergence_check_ind = entropy_consistent;
KL_Kpairwise_consistent = entropy_consistent; 
KL_pairwise_consistent = entropy_consistent;
KL_ind_consistent = entropy_consistent;
num_coactive_probs_consistent = entropy_consistent; 
num_coactive_overunder_Kpairwise_consistent = entropy_consistent;
num_coactive_overunder_pairwise_consistent = entropy_consistent;
num_coactive_overunder_ind_consistent = entropy_consistent;
num_coactive_KLD_Kpairwise_consistent = entropy_consistent;
num_coactive_KLD_pairwise_consistent = entropy_consistent;
num_coactive_KLD_ind_consistent = entropy_consistent;
num_coactive_probs_Kpairwise_consistent = entropy_consistent; 
num_coactive_probs_pairwise_consistent = entropy_consistent;
num_coactive_probs_ind_consistent = entropy_consistent; 
triplet_corr_empirical = entropy_consistent;
triplet_corr_Kpairwise = entropy_consistent;
triplet_corr_pairwise = entropy_consistent;
triplet_corr_ind = entropy_consistent;
KL_Kpairwise_sampled = entropy_consistent;
KL_pairwise_sampled = entropy_consistent;
KL_ind_sampled = entropy_consistent;
num_coactive_probs_sampled_Kpairwise = entropy_consistent;
num_coactive_KLD_sampled_Kpairwise = entropy_consistent;
num_coactive_overunder_sampled_Kpairwise = entropy_consistent;
num_coactive_error_Kpairwise = entropy_consistent;
num_coactive_error_per_pattern_Kpairwise = entropy_consistent;
num_coactive_probs_sampled_pairwise = entropy_consistent;
num_coactive_KLD_sampled_pairwise = entropy_consistent;
num_coactive_overunder_sampled_pairwise = entropy_consistent;
num_coactive_error_pairwise = entropy_consistent;
num_coactive_error_per_pattern_pairwise = entropy_consistent;
num_coactive_probs_sampled_ind = entropy_consistent;
num_coactive_KLD_sampled_ind = entropy_consistent;
num_coactive_overunder_sampled_ind = entropy_consistent;
num_coactive_error_ind = entropy_consistent;
num_coactive_error_per_pattern_ind = entropy_consistent;
subsample_composition_consistent = entropy_consistent; subsample_FR_consistent = entropy_consistent; multibins = entropy_consistent;
mod_rasters_entropy = cell(size(VR_reformatting,1),1); 
entropy_VR_frames = cell(size(VR_reformatting,1),1); entropy_VR_vel = entropy_VR_frames; 
for VR_ephys = 1:size(VR_reformatting,1) 
    for i = 1:size(ks_all_data{VR_reformatting{VR_ephys,5},1},1)
        all_st = ks_all_data{VR_reformatting{VR_ephys,5},1}{i,3}; %all spike times (in ms) for a particular unit
        in_interval_st = all_st(and(all_st >= VR_reformatting{VR_ephys,6}(1),all_st <= VR_reformatting{VR_ephys,6}(end))); %only spike times that fall within bounds of current ephys session
        in_interval_st = in_interval_st-VR_reformatting{VR_ephys,6}(1)+1; %correct for the fact that ephys session might not start at first milisecond. 
        mod_rasters_entropy{VR_ephys}{i,1} = in_interval_st; 
    end
    tmp1_VR_frames = VR_reformatting{VR_ephys,4}; tmp1_VR_vel = VR_reformatting{VR_ephys,10};  
    if mod(numel(tmp1_VR_frames),entropy_bin_size) == 0
        excess = 0; 
    else
        excess = numel(tmp1_VR_frames)-entropy_bin_size*floor(numel(tmp1_VR_frames)/entropy_bin_size);            
    end
    tmp2_VR_frames = tmp1_VR_frames(1:(end-excess)); tmp2_VR_vel = tmp1_VR_vel(1:(end-excess)); 
    tmp3_VR_frames = reshape(tmp2_VR_frames,entropy_bin_size,numel(tmp2_VR_frames)/entropy_bin_size); tmp3_VR_vel = reshape(tmp2_VR_vel,entropy_bin_size,numel(tmp2_VR_vel)/entropy_bin_size); 
    tmp4_VR_frames = mode(tmp3_VR_frames,1); tmp4_VR_vel = mode(tmp3_VR_vel,1); 
    entropy_VR_frames{VR_ephys} = tmp4_VR_frames; entropy_VR_vel{VR_ephys} = tmp4_VR_vel;         
end
for bin_counter = 1:numel(all_bin_sizes)
    entropy_bin_size = all_bin_sizes(bin_counter); 
    for word_length_counter = 1:numel(all_word_lengths)
        word_length = all_word_lengths(word_length_counter); 
        converter = converter_original(1:word_length);                       
        word_library = de2bi(0:(2^word_length)-1,word_length)';
        coactive_library = sum(word_library,1); 
        for VR_ephys = 1:size(VR_reformatting,1)             
            if size(mod_rasters_entropy{VR_ephys},1) >= numel(converter) 
                %Initialize matrices for this combination of an,word length, and bin size
                params_pairwise_consistent{VR_ephys,word_length,bin_counter} = zeros(ent_trials,(((word_length^2)-word_length)/2)+word_length); 
                params_ind_consistent{VR_ephys,word_length,bin_counter} = zeros(ent_trials,word_length); 
                KL_pairwise_consistent{VR_ephys,word_length,bin_counter} = zeros(ent_trials,1); 
                KL_ind_consistent{VR_ephys,word_length,bin_counter} = zeros(ent_trials,1); 
                entropy_consistent{VR_ephys,word_length,bin_counter} = zeros(1,ent_trials); 
                subsample_FR_consistent{VR_ephys,word_length,bin_counter} = entropy_consistent{VR_ephys,word_length,bin_counter}; multibins{VR_ephys,word_length,bin_counter} = entropy_consistent{VR_ephys,word_length,bin_counter};                        
                %Generate binary spike trains from array of spike times for all units
                spike_train_binned_all = zeros(size(mod_rasters_entropy{VR_ephys},1),numel(entropy_VR_frames{VR_ephys})); %matrix of 0's - #units by # of bins in recording 
                for units = 1:size(mod_rasters_entropy{VR_ephys},1) %for each unit 
                    [spike_train_binned_all(units,:),~] = histcounts(mod_rasters_entropy{VR_ephys}{units,1},0.5:entropy_bin_size:(entropy_bin_size*numel(entropy_VR_frames{VR_ephys})+.5)); %count the number of spikes that occur in each bin, starting at .5ms
                end   
                [multibins{VR_ephys,word_length,bin_counter},~] = histcounts(spike_train_binned_all(:),-.5:1:max(spike_train_binned_all(:)));
                spike_train_binned_all = spike_train_binned_all > 0; %this sets all the bins with at least 1 spike to 1 (and leaves the rest as 0's)
                unit_subset = zeros(word_length,ent_trials);
                for m = 1:ent_trials
                    unit_subset(:,m) = randperm(size(mod_rasters_entropy{VR_ephys},1),word_length)';
                end
                subsample_composition_massive{VR_ephys,word_length,bin_counter} = unit_subset;                  
                for i12 = 1:size(unit_subset,2) %for each of the subsamples           
                    spike_train_subset_binned_all = spike_train_binned_all(unit_subset(:,i12),:); %create a new binary spike train matrix, but only with the units of the current subsample
                    spike_train_subset_numerized = converter*spike_train_subset_binned_all; %convert each binary word into a unique integer between 0 and 1023
                    % Entropy
                    [N,~] = histcounts(spike_train_subset_numerized,-.5:1:sum(converter)+.5); 
                    P = N/sum(N); entropy_consistent{VR_ephys,word_length,bin_counter}(i12) = nansum(-P.*log2(P)); 
                    for w = 0:word_length
                        num_coactive_probs_consistent{VR_ephys,word_length,bin_counter}(i12,(w+1)) = sum(P(coactive_library == w));
                    end             
                    if word_length <= 19 
                        % Max entropy model: KPairwise
                        model_Kpairwise_consistent = maxent.createModel(word_length,'kising'); 
                        [model_Kpairwise_consistent,convergence_check_Kpairwise{VR_ephys,word_length,bin_counter}(i12)] = maxent.trainModel(model_Kpairwise_consistent,spike_train_subset_binned_all,'silent',true); 
                        params_Kpairwise_consistent{VR_ephys,word_length,bin_counter}(i12,:) = model_Kpairwise_consistent.factors; 
                        word_probs_Kpairwise = exp(maxent.getLogProbability(model_Kpairwise_consistent,word_library)); 
                        KL_Kpairwise_consistent{VR_ephys,word_length,bin_counter}(i12) = KLDiv(P,word_probs_Kpairwise);                         
                        for w = 0:word_length
                            num_coactive_probs_Kpairwise_consistent{VR_ephys,word_length,bin_counter}(i12,(w+1)) = sum(word_probs_Kpairwise(coactive_library == w));
                            num_coactive_KLD_Kpairwise_consistent{VR_ephys,word_length,bin_counter}(i12,:) = KLDiv_UC2(P,word_probs_Kpairwise,coactive_library); 
                            num_coactive_overunder_Kpairwise_consistent{VR_ephys,word_length,bin_counter}(i12,(w+1)) = sum(word_probs_Kpairwise(coactive_library == w) > P(coactive_library == w))/numel(coactive_library); 
                        end                        
                        entropy_model_Kpairwise_consistent{VR_ephys,word_length,bin_counter}(i12) = nansum(-word_probs_Kpairwise.*log2(word_probs_Kpairwise));
                        % Max entropy model: Pairwise
                        model_pairwise_consistent = maxent.createModel(word_length,'ising'); 
                        [model_pairwise_consistent,convergence_check_pairwise{VR_ephys,word_length,bin_counter}(i12)] = maxent.trainModel(model_pairwise_consistent,spike_train_subset_binned_all,'silent',true); 
                        params_pairwise_consistent{VR_ephys,word_length,bin_counter}(i12,:) = model_pairwise_consistent.factors; 
                        word_probs_pairwise = exp(maxent.getLogProbability(model_pairwise_consistent,word_library)); 
                        KL_pairwise_consistent{VR_ephys,word_length,bin_counter}(i12) = KLDiv(P,word_probs_pairwise);                         
                        for w = 0:word_length
                            num_coactive_probs_pairwise_consistent{VR_ephys,word_length,bin_counter}(i12,(w+1)) = sum(word_probs_pairwise(coactive_library == w));
                            num_coactive_KLD_pairwise_consistent{VR_ephys,word_length,bin_counter}(i12,:) = KLDiv_UC2(P,word_probs_pairwise,coactive_library); 
                            num_coactive_overunder_pairwise_consistent{VR_ephys,word_length,bin_counter}(i12,(w+1)) = sum(word_probs_pairwise(coactive_library == w) > P(coactive_library == w))/numel(coactive_library); 
                        end                        
                        entropy_model_pairwise_consistent{VR_ephys,word_length,bin_counter}(i12) = nansum(-word_probs_pairwise.*log2(word_probs_pairwise));                             
                        % Max entropy model: Independent
                        model_ind_consistent = maxent.createModel(word_length,'indep'); 
                        [model_ind_consistent,convergence_check_ind{VR_ephys,word_length,bin_counter}(i12)] = maxent.trainModel(model_ind_consistent,spike_train_subset_binned_all,'silent',true); 
                        params_ind_consistent{VR_ephys,word_length,bin_counter}(i12,:) = model_ind_consistent.factors; 
                        word_probs_ind = exp(maxent.getLogProbability(model_ind_consistent,word_library)); 
                        KL_ind_consistent{VR_ephys,word_length,bin_counter}(i12) = KLDiv(P,word_probs_ind);                         
                        for w = 0:word_length
                            num_coactive_probs_ind_consistent{VR_ephys,word_length,bin_counter}(i12,(w+1)) = sum(word_probs_ind(coactive_library == w));
                            num_coactive_KLD_ind_consistent{VR_ephys,word_length,bin_counter}(i12,:) = KLDiv_UC2(P,word_probs_ind,coactive_library); 
                            num_coactive_overunder_ind_consistent{VR_ephys,word_length,bin_counter}(i12,(w+1)) = sum(word_probs_ind(coactive_library == w) > P(coactive_library == w))/numel(coactive_library); 
                        end                        
                        entropy_model_ind_consistent{VR_ephys,word_length,bin_counter}(i12) = nansum(-word_probs_ind.*log2(word_probs_ind));                             
                        %FR of each subsample
                        subsample_FR_consistent{VR_ephys,word_length,bin_counter}(i12) = sum(sum(spike_train_subset_binned_all));                                            
                        %Triplet correlations
                        spike_train_Kpairwise = double(maxent.generateSamples(model_Kpairwise_consistent,size(spike_train_subset_binned_all,2)));
                        spike_train_pairwise = double(maxent.generateSamples(model_pairwise_consistent,size(spike_train_subset_binned_all,2)));
                        spike_train_ind = double(maxent.generateSamples(model_ind_consistent,size(spike_train_subset_binned_all,2)));                        
%                         CXYZ_empirical = nan(word_length,word_length,word_length); CXYZ_Kpairwise = CXYZ_empirical; CXYZ_pairwise = CXYZ_empirical; CXYZ_ind = CXYZ_empirical;
%                         parfor ix = 1:word_length
%                             for iy = 1:word_length
%                                 for iz = 1:word_length
%                                     if iy > ix
%                                         if iz > iy
%                                             CXYZ_empirical(ix,iy,iz) = triple_corr(spike_train_subset_binned_all(ix,:),spike_train_subset_binned_all(iy,:),spike_train_subset_binned_all(iz,:));
%                                             CXYZ_Kpairwise(ix,iy,iz) = triple_corr(spike_train_Kpairwise(ix,:),spike_train_Kpairwise(iy,:),spike_train_Kpairwise(iz,:));
%                                             CXYZ_pairwise(ix,iy,iz) = triple_corr(spike_train_pairwise(ix,:),spike_train_pairwise(iy,:),spike_train_pairwise(iz,:));
%                                             CXYZ_ind(ix,iy,iz) = triple_corr(spike_train_ind(ix,:),spike_train_ind(iy,:),spike_train_ind(iz,:));
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                         CXYZ_empirical = CXYZ_empirical(:)'; 
%                         CXYZ_Kpairwise = CXYZ_Kpairwise(:)'; 
%                         CXYZ_pairwise = CXYZ_pairwise(:)'; 
%                         CXYZ_ind = CXYZ_ind(:)'; 
%                         triplet_corr_empirical{VR_ephys,word_length,bin_counter}(i12,:) = CXYZ_empirical(~isnan(CXYZ_empirical));
%                         triplet_corr_Kpairwise{VR_ephys,word_length,bin_counter}(i12,:) = CXYZ_Kpairwise(~isnan(CXYZ_Kpairwise));
%                         triplet_corr_pairwise{VR_ephys,word_length,bin_counter}(i12,:) = CXYZ_pairwise(~isnan(CXYZ_pairwise));
%                         triplet_corr_ind{VR_ephys,word_length,bin_counter}(i12,:) = CXYZ_ind(~isnan(CXYZ_ind));
                        spike_train_Kpairwise_numerized = converter*spike_train_Kpairwise; [N,~] = histcounts(spike_train_Kpairwise_numerized,-.5:1:sum(converter)+.5); word_probs_sampled_Kpairwise = N/sum(N);
                        spike_train_pairwise_numerized = converter*spike_train_pairwise; [N,~] = histcounts(spike_train_pairwise_numerized,-.5:1:sum(converter)+.5); word_probs_sampled_pairwise = N/sum(N);
                        spike_train_ind_numerized = converter*spike_train_ind; [N,~] = histcounts(spike_train_ind_numerized,-.5:1:sum(converter)+.5); word_probs_sampled_ind = N/sum(N);
                        entropy_model_Kpairwise_sampled_consistent{VR_ephys,word_length,bin_counter}(i12) = nansum(-word_probs_sampled_Kpairwise.*log2(word_probs_sampled_Kpairwise));
                        entropy_model_pairwise_sampled_consistent{VR_ephys,word_length,bin_counter}(i12) = nansum(-word_probs_sampled_pairwise.*log2(word_probs_sampled_pairwise));
                        entropy_model_ind_sampled_consistent{VR_ephys,word_length,bin_counter}(i12) = nansum(-word_probs_sampled_ind.*log2(word_probs_sampled_ind));                        
                        KL_Kpairwise_sampled{VR_ephys,word_length,bin_counter}(i12) = KLDiv_inf(P,word_probs_sampled_Kpairwise); 
                        KL_pairwise_sampled{VR_ephys,word_length,bin_counter}(i12) = KLDiv_inf(P,word_probs_sampled_pairwise); 
                        KL_ind_sampled{VR_ephys,word_length,bin_counter}(i12) = KLDiv_inf(P,word_probs_sampled_ind); 
                        Pnan = P; Pnan(Pnan == 0) = NaN;
                        for w = 0:word_length
                            coactive_occurrence = sum(~isnan(Pnan(coactive_library == w))); %number of words in current num_coactive batch that actually occured
                            num_coactive_probs_sampled_Kpairwise{VR_ephys,word_length,bin_counter}(i12,(w+1)) = sum(word_probs_sampled_Kpairwise(coactive_library == w)); 
                            num_coactive_KLD_sampled_Kpairwise{VR_ephys,word_length,bin_counter}(i12,:) = KLDiv_UC2(P,word_probs_sampled_Kpairwise,coactive_library);
                            num_coactive_overunder_sampled_Kpairwise{VR_ephys,word_length,bin_counter}(i12,(w+1)) = sum(word_probs_sampled_Kpairwise(coactive_library == w) > P(coactive_library == w))/numel(coactive_library); 
                            num_coactive_error_Kpairwise{VR_ephys,word_length,bin_counter}(i12,(w+1)) = nansum(abs(word_probs_sampled_Kpairwise(coactive_library == w)-Pnan(coactive_library == w))./Pnan(coactive_library == w))/numel(coactive_library);                            
                            num_coactive_error_per_pattern_Kpairwise{VR_ephys,word_length,bin_counter}(i12,(w+1)) = nansum(abs(word_probs_sampled_pairwise(coactive_library == w)-Pnan(coactive_library == w))./Pnan(coactive_library == w))/coactive_occurrence;                            
                            
                            num_coactive_probs_sampled_pairwise{VR_ephys,word_length,bin_counter}(i12,(w+1)) = sum(word_probs_sampled_pairwise(coactive_library == w)); 
                            num_coactive_KLD_sampled_pairwise{VR_ephys,word_length,bin_counter}(i12,:) = KLDiv_UC2(P,word_probs_sampled_pairwise,coactive_library);
                            num_coactive_overunder_sampled_pairwise{VR_ephys,word_length,bin_counter}(i12,(w+1)) = sum(word_probs_sampled_pairwise(coactive_library == w) > P(coactive_library == w))/numel(coactive_library);                           
                            num_coactive_error_pairwise{VR_ephys,word_length,bin_counter}(i12,(w+1)) = nansum(abs(word_probs_sampled_pairwise(coactive_library == w)-Pnan(coactive_library == w))./Pnan(coactive_library == w))/numel(coactive_library);                            
                            num_coactive_error_per_pattern_pairwise{VR_ephys,word_length,bin_counter}(i12,(w+1)) = nansum(abs(word_probs_sampled_pairwise(coactive_library == w)-Pnan(coactive_library == w))./Pnan(coactive_library == w))/coactive_occurrence;                            
                            
                            num_coactive_probs_sampled_ind{VR_ephys,word_length,bin_counter}(i12,(w+1)) = sum(word_probs_sampled_ind(coactive_library == w)); 
                            num_coactive_KLD_sampled_ind{VR_ephys,word_length,bin_counter}(i12,:) = KLDiv_UC2(P,word_probs_sampled_ind,coactive_library);
                            num_coactive_overunder_sampled_ind{VR_ephys,word_length,bin_counter}(i12,(w+1)) = sum(word_probs_sampled_ind(coactive_library == w) > P(coactive_library == w))/numel(coactive_library); 
                            num_coactive_error_ind{VR_ephys,word_length,bin_counter}(i12,(w+1)) = nansum(abs(word_probs_sampled_ind(coactive_library == w)-Pnan(coactive_library == w))./Pnan(coactive_library == w))/numel(coactive_library);                            
                            num_coactive_error_per_pattern_ind{VR_ephys,word_length,bin_counter}(i12,(w+1)) = nansum(abs(word_probs_sampled_ind(coactive_library == w)-Pnan(coactive_library == w))./Pnan(coactive_library == w))/coactive_occurrence;                            
                        end
                    end
                    entropy_bin_size
                    word_length
                    VR_ephys
                    i12
                end                    
            end            
        end
    end
end
%% Entropy and max entropy over multiple word lengths
ent_trials = 500;
converter_original = 2.^(0:99); 
all_bin_sizes = [10]; entropy_bin_size = all_bin_sizes; 
all_word_lengths = 4:2:16;
entropy_consistent2 = cell(size(VR_reformatting,1),numel(converter_original),numel(all_bin_sizes));
KL_Kpairwise_consistent2 = entropy_consistent2; 
KL_pairwise_consistent2 = entropy_consistent2;
KL_ind_consistent2 = entropy_consistent2;
KL_Kpairwise_sampled2 = entropy_consistent2;
KL_pairwise_sampled2 = entropy_consistent2;
KL_ind_sampled2 = entropy_consistent2;
subsample_FR_consistent2 = entropy_consistent2; multibins = entropy_consistent2;
mod_rasters_entropy = cell(size(VR_reformatting,1),1); 
entropy_VR_frames = cell(size(VR_reformatting,1),1); entropy_VR_vel = entropy_VR_frames; 
for VR_ephys = 1:size(VR_reformatting,1) 
    for i = 1:size(ks_all_data{VR_reformatting{VR_ephys,5},1},1)
        all_st = ks_all_data{VR_reformatting{VR_ephys,5},1}{i,3}; %all spike times (in ms) for a particular unit
        in_interval_st = all_st(and(all_st >= VR_reformatting{VR_ephys,6}(1),all_st <= VR_reformatting{VR_ephys,6}(end))); %only spike times that fall within bounds of current ephys session
        in_interval_st = in_interval_st-VR_reformatting{VR_ephys,6}(1)+1; %correct for the fact that ephys session might not start at first milisecond. 
        mod_rasters_entropy{VR_ephys}{i,1} = in_interval_st; 
    end
    tmp1_VR_frames = VR_reformatting{VR_ephys,4}; tmp1_VR_vel = VR_reformatting{VR_ephys,10};  
    if mod(numel(tmp1_VR_frames),entropy_bin_size) == 0
        excess = 0; 
    else
        excess = numel(tmp1_VR_frames)-entropy_bin_size*floor(numel(tmp1_VR_frames)/entropy_bin_size);            
    end
    tmp2_VR_frames = tmp1_VR_frames(1:(end-excess)); tmp2_VR_vel = tmp1_VR_vel(1:(end-excess)); 
    tmp3_VR_frames = reshape(tmp2_VR_frames,entropy_bin_size,numel(tmp2_VR_frames)/entropy_bin_size); tmp3_VR_vel = reshape(tmp2_VR_vel,entropy_bin_size,numel(tmp2_VR_vel)/entropy_bin_size); 
    tmp4_VR_frames = mode(tmp3_VR_frames,1); tmp4_VR_vel = mode(tmp3_VR_vel,1); 
    entropy_VR_frames{VR_ephys} = tmp4_VR_frames; entropy_VR_vel{VR_ephys} = tmp4_VR_vel;         
end
for bin_counter = 1:numel(all_bin_sizes)
    entropy_bin_size = all_bin_sizes(bin_counter); 
    for word_length_counter = 1:numel(all_word_lengths)
        word_length = all_word_lengths(word_length_counter); 
        converter = converter_original(1:word_length);                       
        word_library = de2bi(0:(2^word_length)-1,word_length)';
        coactive_library = sum(word_library,1); 
        for VR_ephys = 1:size(VR_reformatting,1)             
            if size(mod_rasters_entropy{VR_ephys},1) >= numel(converter) 
                %Initialize matrices for this combination of an,word length, and bin size
                KL_pairwise_consistent2{VR_ephys,word_length,bin_counter} = zeros(ent_trials,1); 
                KL_ind_consistent2{VR_ephys,word_length,bin_counter} = zeros(ent_trials,1); 
                entropy_consistent2{VR_ephys,word_length,bin_counter} = zeros(1,ent_trials); 
                subsample_FR_consistent2{VR_ephys,word_length,bin_counter} = entropy_consistent2{VR_ephys,word_length,bin_counter}; multibins{VR_ephys,word_length,bin_counter} = entropy_consistent2{VR_ephys,word_length,bin_counter};                        
                %Generate binary spike trains from array of spike times for all units
                spike_train_binned_all = zeros(size(mod_rasters_entropy{VR_ephys},1),numel(entropy_VR_frames{VR_ephys})); %matrix of 0's - #units by # of bins in recording 
                for units = 1:size(mod_rasters_entropy{VR_ephys},1) %for each unit 
                    [spike_train_binned_all(units,:),~] = histcounts(mod_rasters_entropy{VR_ephys}{units,1},0.5:entropy_bin_size:(entropy_bin_size*numel(entropy_VR_frames{VR_ephys})+.5)); %count the number of spikes that occur in each bin, starting at .5ms
                end   
                [multibins{VR_ephys,word_length,bin_counter},~] = histcounts(spike_train_binned_all(:),-.5:1:max(spike_train_binned_all(:)));
                spike_train_binned_all = spike_train_binned_all > 0; %this sets all the bins with at least 1 spike to 1 (and leaves the rest as 0's)
                unit_subset = zeros(word_length,ent_trials);
                for m = 1:ent_trials
                    unit_subset(:,m) = randperm(size(mod_rasters_entropy{VR_ephys},1),word_length)';
                end
                subsample_composition_massive{VR_ephys,word_length,bin_counter} = unit_subset;                  
                for i12 = 1:size(unit_subset,2) %for each of the subsamples           
                    spike_train_subset_binned_all = spike_train_binned_all(unit_subset(:,i12),:); %create a new binary spike train matrix, but only with the units of the current subsample
%                     spike_train_subset_binned_stationary = spike_train_subset_binned_all(:,entropy_VR_vel{VR_ephys} <= running_threshold);                                                                 
                    spike_train_subset_numerized = converter*spike_train_subset_binned_all; %convert each binary word into a unique integer between 0 and 1023
                    % Entropy
                    [N,~] = histcounts(spike_train_subset_numerized,-.5:1:sum(converter)+.5); 
                    P = N/sum(N); entropy_consistent2{VR_ephys,word_length,bin_counter}(i12) = nansum(-P.*log2(P));        
                    if word_length <= 19 
                        % Max entropy model: KPairwise
                        model_Kpairwise_consistent = maxent.createModel(word_length,'kising'); 
                        [model_Kpairwise_consistent,~] = maxent.trainModel(model_Kpairwise_consistent,spike_train_subset_binned_all,'silent',true); 
                        word_probs_Kpairwise = exp(maxent.getLogProbability(model_Kpairwise_consistent,word_library)); 
                        KL_Kpairwise_consistent2{VR_ephys,word_length,bin_counter}(i12) = KLDiv(P,word_probs_Kpairwise);                         
                        % Max entropy model: Pairwise
                        model_pairwise_consistent = maxent.createModel(word_length,'ising'); 
                        [model_pairwise_consistent,~] = maxent.trainModel(model_pairwise_consistent,spike_train_subset_binned_all,'silent',true); 
                        word_probs_pairwise = exp(maxent.getLogProbability(model_pairwise_consistent,word_library)); 
                        KL_pairwise_consistent2{VR_ephys,word_length,bin_counter}(i12) = KLDiv(P,word_probs_pairwise);                         
                        % Max entropy model: Independent
                        model_ind_consistent = maxent.createModel(word_length,'indep'); 
                        [model_ind_consistent,~] = maxent.trainModel(model_ind_consistent,spike_train_subset_binned_all,'silent',true); 
                        word_probs_ind = exp(maxent.getLogProbability(model_ind_consistent,word_library)); 
                        KL_ind_consistent2{VR_ephys,word_length,bin_counter}(i12) = KLDiv(P,word_probs_ind);                         
                        % FR of each subsample
                        subsample_FR_consistent2{VR_ephys,word_length,bin_counter}(i12) = sum(sum(spike_train_subset_binned_all));                                            
                        spike_train_Kpairwise = double(maxent.generateSamples(model_Kpairwise_consistent,size(spike_train_subset_binned_all,2)));
                        spike_train_pairwise = double(maxent.generateSamples(model_pairwise_consistent,size(spike_train_subset_binned_all,2)));
                        spike_train_ind = double(maxent.generateSamples(model_ind_consistent,size(spike_train_subset_binned_all,2)));                        
                        spike_train_Kpairwise_numerized = converter*spike_train_Kpairwise; [N,~] = histcounts(spike_train_Kpairwise_numerized,-.5:1:sum(converter)+.5); word_probs_sampled_Kpairwise = N/sum(N);
                        spike_train_pairwise_numerized = converter*spike_train_pairwise; [N,~] = histcounts(spike_train_pairwise_numerized,-.5:1:sum(converter)+.5); word_probs_sampled_pairwise = N/sum(N);
                        spike_train_ind_numerized = converter*spike_train_ind; [N,~] = histcounts(spike_train_ind_numerized,-.5:1:sum(converter)+.5); word_probs_sampled_ind = N/sum(N);
                        KL_Kpairwise_sampled2{VR_ephys,word_length,bin_counter}(i12) = KLDiv_inf(P,word_probs_sampled_Kpairwise); 
                        KL_pairwise_sampled2{VR_ephys,word_length,bin_counter}(i12) = KLDiv_inf(P,word_probs_sampled_pairwise); 
                        KL_ind_sampled2{VR_ephys,word_length,bin_counter}(i12) = KLDiv_inf(P,word_probs_sampled_ind);                                                
                    end
                    entropy_bin_size
                    word_length
                    VR_ephys
                    i12
                end                    
            end            
        end
    end
end
%% Figures: Entropy per spike
wl = 16;
entropy_per_spike_consistent = cell(size(VR_reformatting,1),1);
entropy_per_spike_consistent2 = cell(size(VR_reformatting,1),wl);
for VR_ephys = 1:size(VR_reformatting,1)
    entropy_per_spike_consistent{VR_ephys} = 100*entropy_consistent2{VR_ephys,wl}./(subsample_FR_consistent2{VR_ephys,wl}/numel(entropy_VR_frames{VR_ephys}));
    for w = 1:wl
        if isempty(entropy_consistent2{VR_ephys,w}) == 0
            entropy_per_spike_consistent2{VR_ephys,w} = 100*entropy_consistent2{VR_ephys,w}./(subsample_FR_consistent2{VR_ephys,w}/numel(entropy_VR_frames{VR_ephys}));
        end
    end
end
entropy_WT_dCA1_sites = intersect(WT_dCA1_sites,find(~(cellfun(@isempty,entropy_consistent2(:,wl)))));
entropy_WT_vCA1_sites = intersect(WT_vCA1_sites,find(~(cellfun(@isempty,entropy_consistent2(:,wl)))));
entropy_AD_dCA1_sites = intersect(AD_dCA1_sites,find(~(cellfun(@isempty,entropy_consistent2(:,wl)))));
entropy_AD_vCA1_sites = intersect(AD_vCA1_sites,find(~(cellfun(@isempty,entropy_consistent2(:,wl)))));
epsc2_WT_dCA1 = NaN(1,wl);
epsc2_WT_vCA1 = epsc2_WT_dCA1;
epsc2_AD_dCA1 = epsc2_WT_dCA1;
epsc2_AD_vCA1 = epsc2_WT_dCA1;
error_epsc2_WT_dCA1 = NaN(1,wl);
error_epsc2_WT_vCA1 = epsc2_WT_dCA1;
error_epsc2_AD_dCA1 = epsc2_WT_dCA1;
error_epsc2_AD_vCA1 = epsc2_WT_dCA1;
for word_length = 4:2:wl    
    epsc2_WT_dCA1(1,word_length) = nanmean(cell2mat(entropy_per_spike_consistent2(entropy_WT_dCA1_sites,word_length)'));
    epsc2_WT_vCA1(1,word_length) = nanmean(cell2mat(entropy_per_spike_consistent2(entropy_WT_vCA1_sites,word_length)'));
    epsc2_AD_dCA1(1,word_length) = nanmean(cell2mat(entropy_per_spike_consistent2(entropy_AD_dCA1_sites,word_length)'));
    epsc2_AD_vCA1(1,word_length) = nanmean(cell2mat(entropy_per_spike_consistent2(entropy_AD_vCA1_sites,word_length)'));
    error_epsc2_WT_dCA1(1,word_length) = nanstd(cell2mat(entropy_per_spike_consistent2(entropy_WT_dCA1_sites,word_length)'))/sqrt(sum(~isnan(cell2mat(entropy_per_spike_consistent2(entropy_WT_dCA1_sites,word_length)'))));
    error_epsc2_WT_vCA1(1,word_length) = nanstd(cell2mat(entropy_per_spike_consistent2(entropy_WT_vCA1_sites,word_length)'))/sqrt(sum(~isnan(cell2mat(entropy_per_spike_consistent2(entropy_WT_vCA1_sites,word_length)'))));
    error_epsc2_AD_dCA1(1,word_length) = nanstd(cell2mat(entropy_per_spike_consistent2(entropy_AD_dCA1_sites,word_length)'))/sqrt(sum(~isnan(cell2mat(entropy_per_spike_consistent2(entropy_AD_dCA1_sites,word_length)'))));
    error_epsc2_AD_vCA1(1,word_length) = nanstd(cell2mat(entropy_per_spike_consistent2(entropy_AD_vCA1_sites,word_length)'))/sqrt(sum(~isnan(cell2mat(entropy_per_spike_consistent2(entropy_AD_vCA1_sites,word_length)'))));
    [p,~] = ranksum(cell2mat(entropy_per_spike_consistent2(entropy_WT_vCA1_sites,word_length)'),cell2mat(entropy_per_spike_consistent2(entropy_AD_vCA1_sites,word_length)'))
end
f = figure; set(gcf,'Position',[0 500 2000 500]); 
subplot(1,4,1); 
histogram(cell2mat(entropy_per_spike_consistent(AD_dCA1_sites)'),[0:1:900],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
histogram(cell2mat(entropy_per_spike_consistent(WT_dCA1_sites)'),[0:1:900],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
xlabel('Entropy (bits/(s x spike))'); ylabel('Cumulative Probability'); xlim([300 800]); legend({'WT dCA1','WT vCA1','AD dCA1','AD vCA1'},'location','southeast');
[p,~] = ranksum(cell2mat(entropy_per_spike_consistent(WT_dCA1_sites)'),cell2mat(entropy_per_spike_consistent(AD_dCA1_sites)')); text(600,.6,strcat('p = ',num2str(p))); 
subplot(1,4,2); 
histogram(cell2mat(entropy_per_spike_consistent(AD_vCA1_sites)'),[0:1:900],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
histogram(cell2mat(entropy_per_spike_consistent(WT_vCA1_sites)'),[0:1:900],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
xlabel('Entropy (bits/(s x spike)'); ylabel('Cumulative Probability'); xlim([300 800]); legend({'WT dCA1','WT vCA1','AD dCA1','AD vCA1'},'location','southeast');
[p,~] = ranksum(cell2mat(entropy_per_spike_consistent(WT_vCA1_sites)'),cell2mat(entropy_per_spike_consistent(AD_vCA1_sites)')); text(600,.4,strcat('p = ',num2str(p)));
subplot(1,4,3); 
errorbar(6:2:wl,epsc2_AD_dCA1(6:2:wl),error_epsc2_AD_dCA1(6:2:wl),'r','LineWidth',2); hold on;
errorbar(6:2:wl,epsc2_WT_dCA1(6:2:wl),error_epsc2_WT_dCA1(6:2:wl),'k','LineWidth',2); hold on;
xlim([4 18]); ylim([475 600]); xlabel('Pattern length'); ylabel('Entropy (bits/spike)'); 
subplot(1,4,4); 
errorbar(6:2:wl,epsc2_AD_vCA1(6:2:wl),error_epsc2_AD_vCA1(6:2:wl),'Color',[1 .8 .8],'LineWidth',2); hold on;
errorbar(6:2:wl,epsc2_WT_vCA1(6:2:wl),error_epsc2_WT_vCA1(6:2:wl),'Color',[.8 .8 .8],'LineWidth',2); hold on;
xlim([4 18]); ylim([460 575]); xlabel('Pattern length'); ylabel('Entropy (bits/spike)'); 
% cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 6'); 
% saveas(f,'Consistent entropy_per_spike','fig'); orient(f,'landscape'); print('Consistent entropy_per_spike.pdf','-bestfit','-painters','-dpdf'); close(f); 
%% Entropy with similar FR distributions
FR_delta = .1; %maximum allowable difference between the WT and AD FR's for a given subsample
new_subsample_FR = cell(VR_ephys,1); 
new_subsample_entropy = new_subsample_FR;
for i = 1:size(ks_all_data,1)
    new_subsample_FR{i} = (subsample_FR_consistent2{i,wl}/numel(entropy_VR_frames{i}))/(wl*.01); %average FR of a single unit in the ensemble in Hz
    new_subsample_entropy{i} = entropy_consistent2{i,wl}/.01;
end
% For each entropy sample from the WT_dCA1 group, find a corresponding sample with a similar finding rate in the AD_dCA1 group
all_AD_dCA1_FR = cell2mat(new_subsample_FR(AD_dCA1_sites));
all_AD_dCA1_entropy = cell2mat(new_subsample_entropy(AD_dCA1_sites));
matched_FR_WT_dCA1_entropy = []; 
matched_FR_AD_dCA1_entropy = []; 
matched_FR_WT_dCA1_FR = []; 
matched_FR_AD_dCA1_FR = []; 
for i = 1:numel(WT_dCA1_sites)
    if isempty(new_subsample_FR{WT_dCA1_sites(i)}) == 0
        for j = 1:ent_trials
            tmp_FR_WT = new_subsample_FR{WT_dCA1_sites(i)}(j); %for each WT dCA1 sample, take the firing rate
            AD_candidates = find(and(tmp_FR_WT > (all_AD_dCA1_FR-FR_delta),tmp_FR_WT < (all_AD_dCA1_FR+FR_delta))); %all the AD subsamples where the firing rate is within 0.1 of the FR of the current WT subsample
            if isempty(AD_candidates) == 0
                chosen_idx = randperm(numel(AD_candidates),1);             
                matched_FR_AD_dCA1_entropy = [matched_FR_AD_dCA1_entropy;all_AD_dCA1_entropy(AD_candidates(chosen_idx))];
                matched_FR_AD_dCA1_FR = [matched_FR_AD_dCA1_FR;all_AD_dCA1_FR(AD_candidates(chosen_idx))]; 
                matched_FR_WT_dCA1_entropy = [matched_FR_WT_dCA1_entropy;new_subsample_entropy{WT_dCA1_sites(i)}(j)];
                matched_FR_WT_dCA1_FR = [matched_FR_WT_dCA1_FR;tmp_FR_WT]; 
            else
                matched_FR_AD_dCA1_entropy = [matched_FR_AD_dCA1_entropy;NaN];
                matched_FR_AD_dCA1_FR = [matched_FR_AD_dCA1_FR;NaN]; 
                matched_FR_WT_dCA1_entropy = [matched_FR_WT_dCA1_entropy;NaN];
                matched_FR_WT_dCA1_FR = [matched_FR_WT_dCA1_FR;NaN]; 
            end
        end
    end
end
% For each entropy sample from the WT_vCA1 group, find a corresponding sample with a similar finding rate in the AD_vCA1 group
all_WT_vCA1_FR = cell2mat(new_subsample_FR(WT_vCA1_sites));
all_WT_vCA1_entropy = cell2mat(new_subsample_entropy(WT_vCA1_sites));
matched_FR_AD_vCA1_entropy = []; 
matched_FR_WT_vCA1_entropy = []; 
all_AD_vCA1_num_coactive_probs = cell2mat(new_subsample_num_coactive_probs(AD_vCA1_sites)); 
all_AD_vCA1_KLD = cell2mat(new_subsample_KLD(AD_vCA1_sites)); 
matched_FR_AD_vCA1_FR = []; 
matched_FR_WT_vCA1_FR = []; 
for i = 1:numel(AD_vCA1_sites)
    if isempty(new_subsample_FR{AD_vCA1_sites(i)}) == 0
        for j = 1:ent_trials
            tmp_FR_AD = new_subsample_FR{AD_vCA1_sites(i)}(j); %for each AD vCA1 sample, take the firing rate
            WT_candidates = find(and(tmp_FR_AD > (all_WT_vCA1_FR-FR_delta),tmp_FR_AD < (all_WT_vCA1_FR+FR_delta))); %all the WT subsamples where the firing rate is within 0.1 of the FR of the current AD subsample
            if isempty(WT_candidates) == 0
                chosen_idx = randperm(numel(WT_candidates),1);             
                matched_FR_WT_vCA1_entropy = [matched_FR_WT_vCA1_entropy;all_WT_vCA1_entropy(WT_candidates(chosen_idx))];
                matched_FR_WT_vCA1_FR = [matched_FR_WT_vCA1_FR;all_WT_vCA1_FR(WT_candidates(chosen_idx))];        
                matched_FR_AD_vCA1_entropy = [matched_FR_AD_vCA1_entropy;new_subsample_entropy{AD_vCA1_sites(i)}(j)];
                matched_FR_AD_vCA1_FR = [matched_FR_AD_vCA1_FR;tmp_FR_AD]; 
            else
                matched_FR_WT_vCA1_entropy = [matched_FR_WT_vCA1_entropy;NaN];
                matched_FR_WT_vCA1_FR = [matched_FR_WT_vCA1_FR;NaN];       
                matched_FR_AD_vCA1_entropy = [matched_FR_AD_vCA1_entropy;NaN];
                matched_FR_AD_vCA1_FR = [matched_FR_AD_vCA1_FR;NaN]; 
            end
        end
    end
end
%Plot
f = figure; set(gcf,'Position',[0 -500 1000 1000]); 
subplot(2,2,1); 
histogram(matched_FR_WT_dCA1_FR,[0:.01:10],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(matched_FR_AD_dCA1_FR,[0:.01:10],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
xlabel('Mean FR (Hz)'); ylabel('Cumulative Probability'); xlim([0 10]); ylim([0 1]); axis square; 
text(1,.5,strcat(num2str(nanmean(matched_FR_WT_dCA1_FR)),',',num2str(nanstd(matched_FR_WT_dCA1_FR))));
text(1,.4,strcat(num2str(nanmean(matched_FR_AD_dCA1_FR)),',',num2str(nanstd(matched_FR_AD_dCA1_FR))));
[p,~] = signrank(matched_FR_WT_dCA1_FR,matched_FR_AD_dCA1_FR); text(1,.6,strcat('p = ',num2str(p)));
subplot(2,2,2); 
histogram(matched_FR_WT_dCA1_entropy,[0:1:600],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(matched_FR_AD_dCA1_entropy,[0:1:600],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
xlabel('Entropy (bits/s)'); ylabel('Cumulative Probability'); xlim([100 600]); set(gca,'XScale','log'); ylim([0 1]); axis square; 
text(400,.5,strcat(num2str(nanmean(matched_FR_WT_dCA1_entropy)),',',num2str(nanstd(matched_FR_WT_dCA1_entropy))));
text(400,.4,strcat(num2str(nanmean(matched_FR_AD_dCA1_entropy)),',',num2str(nanstd(matched_FR_AD_dCA1_entropy))));
[p,~] = signrank(matched_FR_WT_dCA1_entropy,matched_FR_AD_dCA1_entropy); text(400,.6,strcat('p = ',num2str(p)));
subplot(2,2,3); 
histogram(matched_FR_WT_vCA1_FR,[0:.01:10],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(matched_FR_AD_vCA1_FR,[0:.01:10],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
xlabel('Mean FR (Hz)'); ylabel('Cumulative Probability'); xlim([0 10]); ylim([0 1]); axis square; 
text(1,.5,strcat(num2str(nanmean(matched_FR_WT_vCA1_FR)),',',num2str(nanstd(matched_FR_WT_vCA1_FR))));
text(1,.4,strcat(num2str(nanmean(matched_FR_AD_vCA1_FR)),',',num2str(nanstd(matched_FR_AD_vCA1_FR))));
[p,~] = signrank(matched_FR_WT_vCA1_FR,matched_FR_AD_vCA1_FR); text(1,.6,strcat('p = ',num2str(p)));
subplot(2,2,4); 
histogram(matched_FR_WT_vCA1_entropy,[0:1:600],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(matched_FR_AD_vCA1_entropy,[0:1:600],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
xlabel('Entropy (bits/s)'); ylabel('Cumulative Probability'); xlim([100 600]); set(gca,'XScale','log'); ylim([0 1]); axis square; 
text(400,.2,strcat(num2str(nanmean(matched_FR_WT_vCA1_entropy)),',',num2str(nanstd(matched_FR_WT_vCA1_entropy))));
text(400,.1,strcat(num2str(nanmean(matched_FR_AD_vCA1_entropy)),',',num2str(nanstd(matched_FR_AD_vCA1_entropy)))); 
[p,~] = signrank(matched_FR_WT_vCA1_entropy,matched_FR_AD_vCA1_entropy); text(400,0,strcat('p = ',num2str(p)));
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 6'); 
saveas(f,'matched_FR_entropy','fig'); orient(f,'portrait'); print('matched_FR_entropy.pdf','-bestfit','-painters','-dpdf'); close(f); 
%% Figures: word probabilities
word_length2 = 16; 
converter = converter_original(1:word_length2); 
ent_ex_VR_ephys = [16,17,1,10]; %WT_dCA1, WT_vCA1, AD_dCA1, AD_vCA1 respectively
color_log = [0 0 0; .8 .8 .8; 1 0 0; 1 .8 .8];
entropy_histogram = []; 
word_xlim = []; 
for i = 1:numel(ent_ex_VR_ephys)
    VR_ephys = ent_ex_VR_ephys(i);
    color_oi = color_log(i,:); 
    spike_train_binned_all = zeros(size(mod_rasters{VR_ephys},1),numel(entropy_VR_frames{VR_ephys})); %matrix of 0's - #units by # of bins in recording 
    for units = 1:size(mod_rasters{VR_ephys},1) %for each unit 
        [spike_train_binned_all(units,:),~] = histcounts(mod_rasters{VR_ephys}{units,1},0.5:entropy_bin_size:(entropy_bin_size*numel(entropy_VR_frames{VR_ephys})+.5)); %count the number of spikes that occur in each bin, starting at .5ms
    end   
    [multibins{VR_ephys,word_length2,bin_counter},~] = histcounts(spike_train_binned_all(:),-.5:1:max(spike_train_binned_all(:)));
    spike_train_binned_all = spike_train_binned_all > 0; %this sets all the bins with at least 1 spike to 1 (and leaves the rest as 0's)
    unit_subset = subsample_composition_massive{VR_ephys,word_length2,bin_counter};
    i12 = randi(ent_trials); 
    spike_train_subset_binned_all = spike_train_binned_all(unit_subset(:,i12),:); %create a new binary spike train matrix, but only with the units of the current subsample
    spike_train_subset_numerized = converter*spike_train_subset_binned_all;
    [N,~] = histcounts(spike_train_subset_numerized,-.5:1:sum(converter)+.5);
    P = N/sum(N);
    f = figure; set(f,'Position',[0 600 1500 400]); 
    subplot(1,4,1); bar(sort(P,'descend'),'FaceColor',color_oi); set(gca,'YScale','log'); 
    [~,word_xlim(i)] = find(sort(P,'descend'),1,'last'); hold on;
    xlim([0 6000]); ylim([.5e-6 1]); ylabel('P(pattern)'); xlabel('Pattern'); set(gca,'YScale','log'); axis square; 
    %Visualize an example pattern
    P_nonzero = P; P_nonzero(P_nonzero == 0) = []; 
    [~,idx_median] = min(abs(P-prctile(P_nonzero,50))); Pmed = P(idx_median); idx_median = idx_median-1; median_word_oi = de2bi(idx_median,word_length,2);    
    P(P == 0) = NaN; 
    [Pmin,idx_min] = nanmin(P); idx_min = idx_min-1; min_word_oi = de2bi(idx_min,word_length,2);
    [~,idx_temp] = nanmax(P); P(idx_temp) = NaN; [Pmax,idx_max] = nanmax(P); idx_max = idx_max-1; max_word_oi = de2bi(idx_max,word_length,2); 
    idx_oi = [idx_max,idx_median,idx_min]; 
    P_oi = [Pmax,Pmed,Pmin]; 
    vertices = ks_all_data{VR_ephys,6};
    vertices_jittered = vertices+80*randn(size(vertices,1),size(vertices,2));
    for j = 2:4
        word_oi = de2bi(idx_oi(j-1),word_length,2); 
        subplot(1,4,j);  
        for u = 1:size(ks_all_data{VR_reformatting{VR_ephys,5},1},1)
            idx2 = find(u == unit_subset(:,i12)); %where in the word of interest is the current unit            
            if any(u == unit_subset(:,i12))
                if word_oi(idx2) == 0
                    scatter(vertices_jittered(u,1),vertices_jittered(u,2),72,color_oi); hold on;
                elseif word_oi(idx2) == 1
                    scatter(vertices_jittered(u,1),vertices_jittered(u,2),72,color_oi,'filled'); hold on;
                end
            end
        end         
        text(0,150,num2str(P_oi(j-1))); hold on;
        plot([0,100],[100,100],'k','LineWidth',2); 
        xlim([-200 600]); ylim([-900 150]); xticks([]); yticks([]); axis square; 
    end  
end
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 6'); 
saveas(gcf,'Word_probabilities_4','fig'); orient(gcf,'landscape'); print('Word_probabilities_4.pdf','-bestfit','-painters','-dpdf'); close(gcf); 
%% Figures: Maximum entropy KLD
wl = 16;
f = figure; set(gcf,'Position',[0 700 2000 500]); 
subplot(1,4,1); 
histogram(cell2mat(KL_pairwise_consistent2(AD_dCA1_sites,wl)),[0:.0001:.11],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
histogram(cell2mat(KL_pairwise_consistent2(WT_dCA1_sites,wl)),[0:.0001:.11],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
xlabel('KL Divergence'); ylabel('Cumulative Probability'); xlim([0 .105]); ylim([0 1]); legend({'WT dCA1','WT vCA1','AD dCA1','AD vCA1'},'location','southeast'); axis square;
[p,~] = ranksum(cell2mat(KL_pairwise_consistent2(WT_dCA1_sites,wl)),cell2mat(KL_pairwise_consistent2(AD_dCA1_sites,wl))); text(.05,.6,strcat('p = ',num2str(p))); 
subplot(1,4,2); 
histogram(cell2mat(KL_pairwise_consistent2(AD_vCA1_sites,wl)),[0:.0001:.025],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
histogram(cell2mat(KL_pairwise_consistent2(WT_vCA1_sites,wl)),[0:.0001:.025],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
xlabel('KL Divergence'); ylabel('Cumulative Probability'); xlim([0 .025]); ylim([0 1]); legend({'WT dCA1','WT vCA1','AD dCA1','AD vCA1'},'location','southeast'); axis square;
[p,~] = ranksum(cell2mat(KL_pairwise_consistent2(WT_vCA1_sites,wl)),cell2mat(KL_pairwise_consistent2(AD_vCA1_sites,wl))); text(.01,.4,strcat('p = ',num2str(p))); 
KLwl_WT_dCA1 = NaN(1,wl);
KLwl_WT_vCA1 = KLwl_WT_dCA1;
KLwl_AD_dCA1 = KLwl_WT_dCA1;
KLwl_AD_vCA1 = KLwl_WT_dCA1;
error_KLwl_WT_dCA1 = NaN(1,wl);
error_KLwl_WT_vCA1 = KLwl_WT_dCA1;
error_KLwl_AD_dCA1 = KLwl_WT_dCA1;
error_KLwl_AD_vCA1 = KLwl_WT_dCA1;
for word_length = 4:2:wl    
    KLwl_WT_dCA1(1,word_length) = nanmean(cell2mat(KL_pairwise_consistent2(entropy_WT_dCA1_sites,word_length)));
    KLwl_WT_vCA1(1,word_length) = nanmean(cell2mat(KL_pairwise_consistent2(entropy_WT_vCA1_sites,word_length)));
    KLwl_AD_dCA1(1,word_length) = nanmean(cell2mat(KL_pairwise_consistent2(entropy_AD_dCA1_sites,word_length)));
    KLwl_AD_vCA1(1,word_length) = nanmean(cell2mat(KL_pairwise_consistent2(entropy_AD_vCA1_sites,word_length)));
    error_KLwl_WT_dCA1(1,word_length) = nanstd(cell2mat(KL_pairwise_consistent2(entropy_WT_dCA1_sites,word_length)))./sqrt(numel(cell2mat(KL_pairwise_consistent2(entropy_WT_dCA1_sites,word_length))));
    error_KLwl_WT_vCA1(1,word_length) = nanstd(cell2mat(KL_pairwise_consistent2(entropy_WT_vCA1_sites,word_length)))./sqrt(numel(cell2mat(KL_pairwise_consistent2(entropy_WT_vCA1_sites,word_length))));
    error_KLwl_AD_dCA1(1,word_length) = nanstd(cell2mat(KL_pairwise_consistent2(entropy_AD_dCA1_sites,word_length)))./sqrt(numel(cell2mat(KL_pairwise_consistent2(entropy_AD_dCA1_sites,word_length))));
    error_KLwl_AD_vCA1(1,word_length) = nanstd(cell2mat(KL_pairwise_consistent2(entropy_AD_vCA1_sites,word_length)))./sqrt(numel(cell2mat(KL_pairwise_consistent2(entropy_AD_vCA1_sites,word_length))));    
end
subplot(1,4,3); 
errorbar(6:2:wl,KLwl_WT_dCA1(6:2:wl),error_KLwl_WT_dCA1(6:2:wl),'k','LineWidth',2); hold on;
errorbar(6:2:wl,KLwl_AD_dCA1(6:2:wl),error_KLwl_AD_dCA1(6:2:wl),'r','LineWidth',2); hold on;
xlim([4 18]); ylim([0 .025]); xlabel('Pattern length'); ylabel('KL Divergence'); axis square;
subplot(1,4,4); 
errorbar(6:2:wl,KLwl_WT_vCA1(6:2:wl),error_KLwl_WT_vCA1(6:2:wl),'Color',[.8 .8 .8],'LineWidth',2); hold on;
errorbar(6:2:wl,KLwl_AD_vCA1(6:2:wl),error_KLwl_AD_vCA1(6:2:wl),'Color',[1 .8 .8],'LineWidth',2);
xlim([4 18]); ylim([0 9e-3]); xlabel('Pattern length'); ylabel('KL Divergence'); axis square;
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 7'); 
saveas(f,'KL_divergence','fig'); orient(f,'landscape'); print('KL_divergence.pdf','-bestfit','-painters','-dpdf'); close(f); 
%% Figures: predicted and empirical pattern probabilities
word_length2 = 16; 
converter = converter_original(1:word_length2); 
ent_ex_VR_ephys = AD_vCA1_sites; %[16,17,1,7]; %WT_dCA1, WT_vCA1, AD_dCA1, AD_vCA1 respectively
color_log = [0 0 0; .8 .8 .8; 1 0 0; 1 .8 .8];
entropy_histogram = []; 
word_xlim = []; 
for i = 1:numel(ent_ex_VR_ephys)
    VR_ephys = ent_ex_VR_ephys(i);
    if isempty(subsample_composition_massive{VR_ephys,word_length2}) == 0
        color_oi = [1 .8 .8]; %color_log(i,:); 
        spike_train_binned_all = zeros(size(mod_rasters{VR_ephys},1),numel(entropy_VR_frames{VR_ephys})); %matrix of 0's - #units by # of bins in recording 
        for units = 1:size(mod_rasters{VR_ephys},1) %for each unit 
            [spike_train_binned_all(units,:),~] = histcounts(mod_rasters{VR_ephys}{units,1},0.5:entropy_bin_size:(entropy_bin_size*numel(entropy_VR_frames{VR_ephys})+.5)); %count the number of spikes that occur in each bin, starting at .5ms
        end   
        spike_train_binned_all = spike_train_binned_all > 0; %this sets all the bins with at least 1 spike to 1 (and leaves the rest as 0's)
        unit_subset = subsample_composition_massive{VR_ephys,word_length2,bin_counter};
        i12 = randi(ent_trials); 
        spike_train_subset_binned_all = spike_train_binned_all(unit_subset(:,i12),:); %create a new binary spike train matrix, but only with the units of the current subsample
        spike_train_subset_numerized = converter*spike_train_subset_binned_all;
        [N,~] = histcounts(spike_train_subset_numerized,-.5:1:sum(converter)+.5);
        P = N/sum(N);           
        % Max entropy model: Pairwise
        model_pairwise_consistent = maxent.createModel(word_length,'ising'); 
        [model_pairwise_consistent,~] = maxent.trainModel(model_pairwise_consistent,spike_train_subset_binned_all,'silent',true); 
        word_probs_pairwise = exp(maxent.getLogProbability(model_pairwise_consistent,word_library));                                 
        f = figure; set(f,'Position',[0 600 400 400]); 
        scatter(P,word_probs_pairwise,8,color_oi,'filled','o'); set(gca,'YScale','log','XScale','log'); 
        xlabel('P_{pattern, empirical}'); ylabel('P_{pattern, predicted}'); ylim([1e-12 1]); xlim([1e-6 1]); 
        hold on; plot([1e-12 1],[1e-12 1],'Color','b','LineWidth',.5); axis square;
    end
end
set(gca,'Ylim',[1e-17 1],'XLim',[1e-6 1]); 
% cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 7'); 
% saveas(gcf,'pred_emp_probs_4','fig'); orient(gcf,'landscape'); print('pred_emp_probs_4.pdf','-bestfit','-painters','-dpdf'); %close(gcf); 
%% Figures: predicted and empirical pattern probabilities color coded by num_coactive
word_length2 = 16; 
word_library = de2bi(0:(2^word_length2)-1,word_length2)';
coactive_library = sum(word_library,1); 
converter = converter_original(1:word_length2); 
ent_ex_VR_ephys = [16 16 16 16 16 16 16 16 16 16]; %[16,17,1,7]; %WT_dCA1, WT_vCA1, AD_dCA1, AD_vCA1 respectively
color_log = [0 0 0; .8 .8 .8; 1 0 0; 1 .8 .8];
entropy_histogram = []; 
word_xlim = []; 
coactive_colors = zeros(numel(coactive_library),3); 
% ccl = hsv(17); 
% for w = 1:numel(coactive_library)
%     coactive_colors(w,:) = ccl((coactive_library(w)+1),:);
% end
ccl = hsv(8); 
for w = 1:numel(coactive_library)
    if coactive_library(w) <= size(ccl,1)-1
        coactive_colors(w,:) = ccl((coactive_library(w)+1),:);
    elseif coactive_library(w) > size(ccl,1)-1
        coactive_colors(w,:) = ccl(end,:); 
    end
end
for i = 1:numel(ent_ex_VR_ephys)
    VR_ephys = ent_ex_VR_ephys(i);
    if isempty(subsample_composition_massive{VR_ephys,word_length2}) == 0
        spike_train_binned_all = zeros(size(mod_rasters{VR_ephys},1),numel(entropy_VR_frames{VR_ephys})); %matrix of 0's - #units by # of bins in recording 
        for units = 1:size(mod_rasters{VR_ephys},1) %for each unit 
            [spike_train_binned_all(units,:),~] = histcounts(mod_rasters{VR_ephys}{units,1},0.5:entropy_bin_size:(entropy_bin_size*numel(entropy_VR_frames{VR_ephys})+.5)); %count the number of spikes that occur in each bin, starting at .5ms
        end   
        spike_train_binned_all = spike_train_binned_all > 0; %this sets all the bins with at least 1 spike to 1 (and leaves the rest as 0's)
        unit_subset = subsample_composition_massive{VR_ephys,word_length2,bin_counter};
        i12 = randi(ent_trials); 
        spike_train_subset_binned_all = spike_train_binned_all(unit_subset(:,i12),:); %create a new binary spike train matrix, but only with the units of the current subsample
        spike_train_subset_numerized = converter*spike_train_subset_binned_all;
        [N,~] = histcounts(spike_train_subset_numerized,-.5:1:sum(converter)+.5);
        P = N/sum(N);           
        % Max entropy model: Pairwise
        model_pairwise_consistent = maxent.createModel(word_length,'ising'); 
        [model_pairwise_consistent,~] = maxent.trainModel(model_pairwise_consistent,spike_train_subset_binned_all,'silent',true); 
        word_probs_pairwise = exp(maxent.getLogProbability(model_pairwise_consistent,word_library)); 
        f = figure; set(f,'Position',[0 0 1000 1000]); 
%         scatter(P,word_probs_pairwise,42,coactive_colors,'filled','o'); set(gca,'YScale','log','XScale','log');         
        for j = 0:word_length2
            scatter(P(coactive_library == j),word_probs_pairwise(coactive_library == j),150,coactive_colors(coactive_library == j,:),'filled','o'); hold on;
        end
        set(gca,'YScale','log','XScale','log'); set(gca,'Ylim',[1e-17 1],'XLim',[5e-7 1]);         
        xlabel('P_{pattern, empirical}'); ylabel('P_{pattern, predicted}'); 
        hold on; plot([1e-12 1],[1e-12 1],'Color',[.8 .8 .8],'LineWidth',.5); axis square; 
        for j = 0:size(ccl,1)-1
            scatter(.9,10.^(-.5*(j+2)),150,ccl((j+1),:),'filled','o')
        end
    end
end
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 8'); 
saveas(gcf,'pred_emp_probs_num_coactive4','fig'); orient(gcf,'landscape'); print('pred_emp_probs_num_coactive4.pdf','-bestfit','-painters','-dpdf'); close(gcf); 
%% Figures: word probability distributions (empirical vs. pairwise model resampled predictions)
[norm_coactive,~] = histcounts(coactive_library,-.5:1:wl+1); 
error_num_coactive_error_pairwise_WT_dCA1 = nanstd(cell2mat(num_coactive_error_pairwise(WT_dCA1_sites,wl)))./sqrt(numel(cell2mat(KL_pairwise_consistent2(entropy_WT_dCA1_sites,word_length))));
error_num_coactive_error_pairwise_AD_dCA1 = nanstd(cell2mat(num_coactive_error_pairwise(AD_dCA1_sites,wl)))./sqrt(numel(cell2mat(KL_pairwise_consistent2(entropy_AD_dCA1_sites,word_length))));
error_num_coactive_error_pairwise_WT_vCA1 = nanstd(cell2mat(num_coactive_error_pairwise(WT_vCA1_sites,wl)))./sqrt(numel(cell2mat(KL_pairwise_consistent2(entropy_WT_vCA1_sites,word_length))));
error_num_coactive_error_pairwise_AD_vCA1 = nanstd(cell2mat(num_coactive_error_pairwise(AD_vCA1_sites,wl)))./sqrt(numel(cell2mat(KL_pairwise_consistent2(entropy_AD_vCA1_sites,word_length))));
f = figure; set(gcf,'Position',[0 0 1000 1000]); 
subplot(2,2,1); 
plot(0:wl,mean(cell2mat(num_coactive_probs_consistent(WT_dCA1_sites,wl))),'k','LineWidth',2); hold on;
plot(0:wl,mean(cell2mat(num_coactive_probs_consistent(AD_dCA1_sites,wl))),'r','LineWidth',2); hold on;
plot(0:wl,mean(cell2mat(num_coactive_probs_sampled_pairwise(WT_dCA1_sites,wl))),'k','LineWidth',2,'LineStyle','--'); hold on;
plot(0:wl,mean(cell2mat(num_coactive_probs_sampled_pairwise(AD_dCA1_sites,wl))),'r','LineWidth',2,'LineStyle','--'); hold on;
xlabel('Number of coactive units'); ylabel('Pattern probability'); set(gca,'YScale','log'); xlim([0 17]); axis square; 
legend({'WT dCA1 Empirical','AD dCA1 Empirical','WT dCA1 Predicted','AD dCA1 Predicted'},'location','southwest'); legend boxoff; 
subplot(2,2,2); 
plot(0:wl,mean(cell2mat(num_coactive_probs_consistent(WT_vCA1_sites,wl))),'Color',[.8 .8 .8],'LineWidth',2); hold on;
plot(0:wl,mean(cell2mat(num_coactive_probs_consistent(AD_vCA1_sites,wl))),'Color',[1 .8 .8],'LineWidth',2); 
plot(0:wl,mean(cell2mat(num_coactive_probs_sampled_pairwise(WT_vCA1_sites,wl))),'Color',[.8 .8 .8],'LineWidth',2,'LineStyle','--'); hold on;
plot(0:wl,mean(cell2mat(num_coactive_probs_sampled_pairwise(AD_vCA1_sites,wl))),'Color',[1 .8 .8],'LineWidth',2,'LineStyle','--'); 
xlabel('Number of coactive units'); ylabel('Pattern probability'); set(gca,'YScale','log'); xlim([0 17]); axis square; 
legend({'WT vCA1 Empirical','AD vCA1 Empirical','WT vCA1 Predicted','AD vCA1 Predicted'},'location','southwest'); legend boxoff; 
subplot(2,2,3); 
[x,y] = bargraph2line(-.5:1:(wl-.5),nanmean(cell2mat(num_coactive_error_pairwise(WT_dCA1_sites,wl))));
plot(x,y,'k','LineWidth',2); hold on;
errorbar(0:wl,nanmean(cell2mat(num_coactive_error_pairwise(WT_dCA1_sites,wl))),error_num_coactive_error_pairwise_WT_dCA1,'k','LineWidth',2,'LineStyle','none'); hold on;
[x,y] = bargraph2line(-.5:1:(wl-.5),nanmean(cell2mat(num_coactive_error_pairwise(AD_dCA1_sites,wl))));
plot(x,y,'r','LineWidth',2); hold on;
errorbar(0:wl,nanmean(cell2mat(num_coactive_error_pairwise(AD_dCA1_sites,wl))),error_num_coactive_error_pairwise_AD_dCA1,'r','LineWidth',2,'LineStyle','none'); hold on;
xlabel('Number of coactive units'); ylabel('Total prediction error'); xlim([0 17]); axis square; 
legend({'WT dCA1','AD dCA1'},'location','northeast'); legend boxoff; 
subplot(2,2,4); 
[x,y] = bargraph2line(-.5:1:(wl-.5),nanmean(cell2mat(num_coactive_error_pairwise(WT_vCA1_sites,wl))));
plot(x,y,'Color',[.8 .8 .8],'LineWidth',2); hold on;
errorbar(0:wl,nanmean(cell2mat(num_coactive_error_pairwise(WT_vCA1_sites,wl))),error_num_coactive_error_pairwise_WT_vCA1,'Color',[.8 .8 .8],'LineWidth',2,'LineStyle','none'); hold on;
[x,y] = bargraph2line(-.5:1:(wl-.5),nanmean(cell2mat(num_coactive_error_pairwise(AD_vCA1_sites,wl))));
plot(x,y,'Color',[1 .8 .8],'LineWidth',2); hold on;
errorbar(0:wl,nanmean(cell2mat(num_coactive_error_pairwise(AD_vCA1_sites,wl))),error_num_coactive_error_pairwise_AD_vCA1,'Color',[1 .8 .8],'LineWidth',2,'LineStyle','none'); hold on;
xlabel('Number of coactive units'); ylabel('Total prediction error'); xlim([0 17]); axis square; 
legend({'WT vCA1','AD vCA1'},'location','northeast'); legend boxoff; 
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 8'); 
saveas(f,'num_coactive_prediction_error','fig'); orient(f,'landscape'); print('num_coactive_prediction_error.pdf','-bestfit','-painters','-dpdf'); close(f); 
% error per pattern
error_num_coactive_error_per_pattern_pairwise_WT_dCA1 = nanstd(cell2mat(num_coactive_error_per_pattern_pairwise(WT_dCA1_sites,wl)))./sqrt(numel(cell2mat(KL_pairwise_consistent2(entropy_WT_dCA1_sites,word_length))));
error_num_coactive_error_per_pattern_pairwise_AD_dCA1 = nanstd(cell2mat(num_coactive_error_per_pattern_pairwise(AD_dCA1_sites,wl)))./sqrt(numel(cell2mat(KL_pairwise_consistent2(entropy_AD_dCA1_sites,word_length))));
error_num_coactive_error_per_pattern_pairwise_WT_vCA1 = nanstd(cell2mat(num_coactive_error_per_pattern_pairwise(WT_vCA1_sites,wl)))./sqrt(numel(cell2mat(KL_pairwise_consistent2(entropy_WT_vCA1_sites,word_length))));
error_num_coactive_error_per_pattern_pairwise_AD_vCA1 = nanstd(cell2mat(num_coactive_error_per_pattern_pairwise(AD_vCA1_sites,wl)))./sqrt(numel(cell2mat(KL_pairwise_consistent2(entropy_AD_vCA1_sites,word_length))));
g = figure; set(gcf,'Position',[0 0 1000 500]); 
subplot(1,2,1); 
[x,y] = bargraph2line(-.5:1:(wl-.5),nanmean(cell2mat(num_coactive_error_per_pattern_pairwise(WT_dCA1_sites,wl))));
plot(x,y,'k','LineWidth',2); hold on;
errorbar(0:wl,nanmean(cell2mat(num_coactive_error_per_pattern_pairwise(WT_dCA1_sites,wl))),error_num_coactive_error_per_pattern_pairwise_WT_dCA1,'k','LineWidth',2,'LineStyle','none'); hold on;
[x,y] = bargraph2line(-.5:1:(wl-.5),nanmean(cell2mat(num_coactive_error_per_pattern_pairwise(AD_dCA1_sites,wl))));
plot(x,y,'r','LineWidth',2); hold on;
errorbar(0:wl,nanmean(cell2mat(num_coactive_error_per_pattern_pairwise(AD_dCA1_sites,wl))),error_num_coactive_error_per_pattern_pairwise_AD_dCA1,'r','LineWidth',2,'LineStyle','none'); hold on;
xlabel('Number of coactive units'); ylabel('Prediction error per pattern'); xlim([0 17]); axis square; 
subplot(1,2,2); 
[x,y] = bargraph2line(-.5:1:(wl-.5),nanmean(cell2mat(num_coactive_error_per_pattern_pairwise(WT_vCA1_sites,wl))));
plot(x,y,'Color',[.8 .8 .8],'LineWidth',2); hold on;
errorbar(0:wl,nanmean(cell2mat(num_coactive_error_per_pattern_pairwise(WT_vCA1_sites,wl))),error_num_coactive_error_per_pattern_pairwise_WT_vCA1,'Color',[.8 .8 .8],'LineWidth',2,'LineStyle','none'); hold on;
[x,y] = bargraph2line(-.5:1:(wl-.5),nanmean(cell2mat(num_coactive_error_per_pattern_pairwise(AD_vCA1_sites,wl))));
plot(x,y,'Color',[1 .8 .8],'LineWidth',2); hold on;
errorbar(0:wl,nanmean(cell2mat(num_coactive_error_per_pattern_pairwise(AD_vCA1_sites,wl))),error_num_coactive_error_per_pattern_pairwise_AD_vCA1,'Color',[1 .8 .8],'LineWidth',2,'LineStyle','none'); hold on;
xlabel('Number of coactive units'); ylabel('Prediction error per pattern'); xlim([0 17]); axis square; 
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 8'); 
saveas(g,'num_coactive_prediction_error_per_pattern','fig'); orient(g,'landscape'); print('num_coactive_prediction_error_per_pattern.pdf','-bestfit','-painters','-dpdf'); close(g); 
%% Integrate Corr, Entropy, and KLD as scatter plot and build classifier
wl = 16; 
%Generate subsampled correlations
corrsub = cell(size(VR_reformatting,1),1); 
for VR_ephys = 1:size(VR_reformatting,1)
    if isempty(corr_consistent{VR_ephys,wl,2}) == 0
        corrsub{VR_ephys} = zeros(1,ent_trials); 
        for i12 = 1:ent_trials
            tmpcorr = squareform(corr_consistent{VR_ephys,wl,2});
            tmpcorr2 = tmpcorr(subsample_composition_massive{VR_ephys,wl}(:,i12),subsample_composition_massive{VR_ephys,wl}(:,i12));
            tmp_multiplier = tril(NaN(numel(subsample_composition_massive{VR_ephys,wl}(:,i12)),numel(subsample_composition_massive{VR_ephys,wl}(:,i12)))); 
            tmp_multiplier(~isnan(tmp_multiplier)) = 1; 
            tmpcorr3 = tmp_multiplier.*tmpcorr2;
            corrsub{VR_ephys}(i12) = nanmean(tmpcorr3(:));         
        end
    end
end 
%% Classifier for all sessions
class_WT_dCA1_allsess = [cell2mat(entropy_per_spike_consistent(entropy_WT_dCA1_sites)')',cell2mat(corrsub(entropy_WT_dCA1_sites)')',cell2mat(KL_pairwise_consistent2(entropy_WT_dCA1_sites,wl))]; %class_WT_dCA1_allsess = class_WT_dCA1_allsess(randperm(size(class_WT_dCA1_allsess,1),2000)',:);
class_WT_vCA1_allsess = [cell2mat(entropy_per_spike_consistent(entropy_WT_vCA1_sites)')',cell2mat(corrsub(entropy_WT_vCA1_sites)')',cell2mat(KL_pairwise_consistent2(entropy_WT_vCA1_sites,wl))]; %class_WT_vCA1_allsess = class_WT_vCA1_allsess(randperm(size(class_WT_vCA1_allsess,1),2000)',:);
class_AD_dCA1_allsess = [cell2mat(entropy_per_spike_consistent(entropy_AD_dCA1_sites)')',cell2mat(corrsub(entropy_AD_dCA1_sites)')',cell2mat(KL_pairwise_consistent2(entropy_AD_dCA1_sites,wl))]; %class_AD_dCA1_allsess = class_AD_dCA1_allsess(randperm(size(class_AD_dCA1_allsess,1),2000)',:);
class_AD_vCA1_allsess = [cell2mat(entropy_per_spike_consistent(entropy_AD_vCA1_sites)')',cell2mat(corrsub(entropy_AD_vCA1_sites)')',cell2mat(KL_pairwise_consistent2(entropy_AD_vCA1_sites,wl))]; %class_AD_vCA1_allsess = class_AD_vCA1_allsess(randperm(size(class_AD_vCA1_allsess,1),2000)',:);
num_samples = 100; 
allsess_log_WT_DV = zeros(num_samples,1); 
allsess_log_AD_DV = allsess_log_WT_DV; 
allsess_log_WT_DV_null = allsess_log_WT_DV;
allsess_log_AD_DV_null = allsess_log_WT_DV;
allsess_log_all = allsess_log_WT_DV;
allsess_log_null_all = allsess_log_WT_DV;
for i = 1:num_samples
    train_WT_dCA1_idx = randperm(size(class_WT_dCA1_allsess,1),round(.8*size(class_WT_dCA1_allsess,1)));
    train_AD_dCA1_idx = randperm(size(class_AD_dCA1_allsess,1),round(.8*size(class_AD_dCA1_allsess,1)));
    train_WT_vCA1_idx = randperm(size(class_WT_vCA1_allsess,1),round(.8*size(class_WT_vCA1_allsess,1)));
    train_AD_vCA1_idx = randperm(size(class_AD_vCA1_allsess,1),round(.8*size(class_AD_vCA1_allsess,1)));
    test_WT_dCA1_idx = setdiff(1:size(class_WT_dCA1_allsess,1),train_WT_dCA1_idx); 
    test_AD_dCA1_idx = setdiff(1:size(class_AD_dCA1_allsess,1),train_AD_dCA1_idx); 
    test_WT_vCA1_idx = setdiff(1:size(class_WT_vCA1_allsess,1),train_WT_vCA1_idx); 
    test_AD_vCA1_idx = setdiff(1:size(class_AD_vCA1_allsess,1),train_AD_vCA1_idx); 
    train_WT_dCA1 = [class_WT_dCA1_allsess(train_WT_dCA1_idx,:),ones(round(.8*size(class_WT_dCA1_allsess,1)),1)]; 
    train_AD_dCA1 = [class_AD_dCA1_allsess(train_AD_dCA1_idx,:),2*ones(round(.8*size(class_AD_dCA1_allsess,1)),1)]; 
    train_WT_vCA1 = [class_WT_vCA1_allsess(train_WT_vCA1_idx,:),3*ones(round(.8*size(class_WT_vCA1_allsess,1)),1)]; 
    train_AD_vCA1 = [class_AD_vCA1_allsess(train_AD_vCA1_idx,:),4*ones(round(.8*size(class_AD_vCA1_allsess,1)),1)]; 
    train_WT_dCA1_null = [train_WT_dCA1(:,1:3),2*randi(2,size(train_WT_dCA1,1),1)-1];
    train_AD_dCA1_null = [train_AD_dCA1(:,1:3),2*randi(2,size(train_AD_dCA1,1),1)];
    train_WT_vCA1_null = [train_WT_vCA1(:,1:3),2*randi(2,size(train_WT_vCA1,1),1)-1];
    train_AD_vCA1_null = [train_AD_vCA1(:,1:3),2*randi(2,size(train_AD_vCA1,1),1)];
    train_WT_dCA1_null_all = [train_WT_dCA1(:,1:3),randi(4,size(train_WT_dCA1,1),1)];
    train_AD_dCA1_null_all = [train_AD_dCA1(:,1:3),randi(4,size(train_AD_dCA1,1),1)];
    train_WT_vCA1_null_all = [train_WT_vCA1(:,1:3),randi(4,size(train_WT_vCA1,1),1)];
    train_AD_vCA1_null_all = [train_AD_vCA1(:,1:3),randi(4,size(train_AD_vCA1,1),1)];
    record1 = zeros(size(test_WT_dCA1_idx,1),1); record1_null = record1; record1_all = record1; record1_null_all = record1; 
    record2 = zeros(size(test_WT_vCA1_idx,1),1); record2_null = record2; record2_all = record2; record2_null_all = record2; 
    record3 = zeros(size(test_AD_dCA1_idx,1),1); record3_null = record3; record3_all = record3; record3_null_all = record3;
    record4 = zeros(size(test_AD_vCA1_idx,1),1); record4_null = record4; record4_all = record4; record4_null_all = record4;
    for j = 1:100
        %% WT D vs. V
        %find the nearest three neighbors 
        test1 = class_WT_dCA1_allsess(test_WT_dCA1_idx(j),:);
        distances = sum((test1-[train_WT_dCA1(:,1:3);train_WT_vCA1(:,1:3)]).^2,2);
        distances_all = sum((test1-[train_WT_dCA1(:,1:3);train_AD_dCA1(:,1:3);train_WT_vCA1(:,1:3);train_AD_vCA1(:,1:3)]).^2,2); 
        [~,I] = mink(distances,3);    
        [~,I_all] = mink(distances_all,5); 
        %find their identities
        id_combined = [train_WT_dCA1(:,4);train_WT_vCA1(:,4)];
        id_combined_null = [train_WT_dCA1_null(:,4);train_WT_vCA1_null(:,4)]; 
        id_combined_all = [train_WT_dCA1(:,4);train_AD_dCA1(:,4);train_WT_vCA1(:,4);train_AD_vCA1(:,4)];
        id_combined_null_all = [train_WT_dCA1_null_all(:,4);train_AD_dCA1_null_all(:,4);train_WT_vCA1_null_all(:,4);train_AD_vCA1_null_all(:,4)];
        %assign the mode identity to the test point
        prediction = mode(id_combined(I));
        prediction_null = mode(id_combined_null(I)); 
        prediction_all = mode(id_combined_all(I_all)); 
        prediction_null_all = mode(id_combined_null_all(I_all));
        record1(j) = prediction == 1;
        record1_null(j) = prediction_null == 1; 
        record1_all(j) = prediction_all == 1;
        record1_null_all(j) = prediction_null_all == 1;         
        test2 = class_WT_vCA1_allsess(test_WT_vCA1_idx(j),:);
        distances = sum((test2-[train_WT_dCA1(:,1:3);train_WT_vCA1(:,1:3)]).^2,2);
        distances_all = sum((test2-[train_WT_dCA1(:,1:3);train_AD_dCA1(:,1:3);train_WT_vCA1(:,1:3);train_AD_vCA1(:,1:3)]).^2,2); 
        [~,I] = mink(distances,3);    
        [~,I_all] = mink(distances_all,5);         
        id_combined = [train_WT_dCA1(:,4);train_WT_vCA1(:,4)];
        id_combined_null = [train_WT_dCA1_null(:,4);train_WT_vCA1_null(:,4)]; 
        id_combined_all = [train_WT_dCA1(:,4);train_AD_dCA1(:,4);train_WT_vCA1(:,4);train_AD_vCA1(:,4)];    
        id_combined_null_all = [train_WT_dCA1_null_all(:,4);train_AD_dCA1_null_all(:,4);train_WT_vCA1_null_all(:,4);train_AD_vCA1_null_all(:,4)];        
        prediction = mode(id_combined(I));
        prediction_null = mode(id_combined_null(I)); 
        prediction_all = mode(id_combined_all(I_all)); 
        prediction_null_all = mode(id_combined_null_all(I_all));        
        record2(j) = prediction == 3;
        record2_null(j) = prediction_null == 3; 
        record2_all(j) = prediction_all == 3;
        record2_null_all(j) = prediction_null_all == 3; 
        %% AD D vs. V
        test3 = class_AD_dCA1_allsess(test_AD_dCA1_idx(j),:);
        distances = sum((test3-[train_AD_dCA1(:,1:3);train_AD_vCA1(:,1:3)]).^2,2);
        distances_all = sum((test3-[train_WT_dCA1(:,1:3);train_AD_dCA1(:,1:3);train_WT_vCA1(:,1:3);train_AD_vCA1(:,1:3)]).^2,2); 
        [~,I] = mink(distances,3);   
        [~,I_all] = mink(distances_all,5);         
        id_combined = [train_AD_dCA1(:,4);train_AD_vCA1(:,4)];
        id_combined_null = [train_AD_dCA1_null(:,4);train_AD_vCA1_null(:,4)]; 
        id_combined_all = [train_WT_dCA1(:,4);train_AD_dCA1(:,4);train_WT_vCA1(:,4);train_AD_vCA1(:,4)];      
        id_combined_null_all = [train_WT_dCA1_null_all(:,4);train_AD_dCA1_null_all(:,4);train_WT_vCA1_null_all(:,4);train_AD_vCA1_null_all(:,4)];        
        prediction = mode(id_combined(I));
        prediction_null = mode(id_combined_null(I)); 
        prediction_all = mode(id_combined_all(I_all));  
        prediction_null_all = mode(id_combined_null_all(I_all));        
        record3(j) = prediction == 2;
        record3_null(j) = prediction_null == 2;
        record3_all(j) = prediction_all == 2;
        record3_null_all(j) = prediction_null_all == 2;
        test4 = class_AD_vCA1_allsess(test_AD_vCA1_idx(j),:);
        distances = sum((test4-[train_AD_dCA1(:,1:3);train_AD_vCA1(:,1:3)]).^2,2);
        distances_all = sum((test4-[train_WT_dCA1(:,1:3);train_AD_dCA1(:,1:3);train_WT_vCA1(:,1:3);train_AD_vCA1(:,1:3)]).^2,2); 
        [~,I] = mink(distances,3);    
        [~,I_all] = mink(distances_all,5);         
        id_combined = [train_AD_dCA1(:,4);train_AD_vCA1(:,4)];
        id_combined_null = [train_AD_dCA1_null(:,4);train_AD_vCA1_null(:,4)]; 
        id_combined_all = [train_WT_dCA1(:,4);train_AD_dCA1(:,4);train_WT_vCA1(:,4);train_AD_vCA1(:,4)];        
        id_combined_null_all = [train_WT_dCA1_null_all(:,4);train_AD_dCA1_null_all(:,4);train_WT_vCA1_null_all(:,4);train_AD_vCA1_null_all(:,4)];        
        prediction = mode(id_combined(I));
        prediction_null = mode(id_combined_null(I)); 
        prediction_all = mode(id_combined_all(I_all));  
        prediction_null_all = mode(id_combined_null_all(I_all));
        record4(j) = prediction == 4;
        record4_null(j) = prediction_null == 4;
        record4_all(j) = prediction_all == 4; 
        record4_null_all(j) = prediction_null_all == 4;         
    end
    allsess_log_WT_DV(i) = mean([record1,record2]); 
    allsess_log_AD_DV(i) = mean([record3,record4]); 
    allsess_log_WT_DV_null(i) = mean([record1_null,record2_null]); 
    allsess_log_AD_DV_null(i) = mean([record3_null,record4_null]); 
    allsess_log_all(i) = mean([record1_all,record2_all,record3_all,record4_all]);
    allsess_log_null_all(i) = mean([record1_null_all,record2_null_all,record3_null_all,record4_null_all]); 
end
%% Figure: Scatter plot of corr, entropy, and KLD and classifier results
f = figure; set(f,'Position',[0 0 1500 1000]);
subplot(3,3,1); bins = [0:.002:.12];
histogram(cell2mat(corrsub(entropy_WT_dCA1_sites)')',bins,'DisplayStyle','stairs','normalization','probability','LineWidth',1.5,'EdgeColor','k'); hold on;
histogram(cell2mat(corrsub(entropy_AD_dCA1_sites)')',bins,'DisplayStyle','stairs','normalization','probability','LineWidth',1.5,'EdgeColor','r'); hold on;
histogram(cell2mat(corrsub(entropy_WT_vCA1_sites)')',bins,'DisplayStyle','stairs','normalization','probability','LineWidth',1.5,'EdgeColor',[.8 .8 .8]); hold on;
histogram(cell2mat(corrsub(entropy_AD_vCA1_sites)')',bins,'DisplayStyle','stairs','normalization','probability','LineWidth',1.5,'EdgeColor',[1 .8 .8]); hold on;
xlim([0 .12]); xlabel('Correlation'); ylabel('Probability'); 
subplot(3,3,4); bins = [0:10:800];
histogram(cell2mat(entropy_per_spike_consistent(entropy_WT_dCA1_sites)')',bins,'DisplayStyle','stairs','normalization','probability','LineWidth',1.5,'EdgeColor','k'); hold on;
histogram(cell2mat(entropy_per_spike_consistent(entropy_AD_dCA1_sites)')',bins,'DisplayStyle','stairs','normalization','probability','LineWidth',1.5,'EdgeColor','r'); hold on;
histogram(cell2mat(entropy_per_spike_consistent(entropy_WT_vCA1_sites)')',bins,'DisplayStyle','stairs','normalization','probability','LineWidth',1.5,'EdgeColor',[.8 .8 .8]); hold on;
histogram(cell2mat(entropy_per_spike_consistent(entropy_AD_vCA1_sites)')',bins,'DisplayStyle','stairs','normalization','probability','LineWidth',1.5,'EdgeColor',[1 .8 .8]); hold on;
xlim([200 800]); xlabel('Entropy (bits/spike)'); ylabel('Probability'); 
subplot(3,3,7); bins = [0:.001:.08];
histogram(cell2mat(KL_pairwise_consistent2(entropy_WT_dCA1_sites,wl)),bins,'DisplayStyle','stairs','normalization','probability','LineWidth',1.5,'EdgeColor','k'); hold on;
histogram(cell2mat(KL_pairwise_consistent2(entropy_AD_dCA1_sites,wl)),bins,'DisplayStyle','stairs','normalization','probability','LineWidth',1.5,'EdgeColor','r'); hold on;
histogram(cell2mat(KL_pairwise_consistent2(entropy_WT_vCA1_sites,wl)),bins,'DisplayStyle','stairs','normalization','probability','LineWidth',1.5,'EdgeColor',[.8 .8 .8]); hold on;
histogram(cell2mat(KL_pairwise_consistent2(entropy_AD_vCA1_sites,wl)),bins,'DisplayStyle','stairs','normalization','probability','LineWidth',1.5,'EdgeColor',[1 .8 .8]); hold on;
xlim([0 .08]); xlabel('KLD'); ylabel('Probability'); 
subplot(3,3,[2,3,5,6]); 
scatter3(class_WT_dCA1_allsess(:,1),class_WT_dCA1_allsess(:,2),class_WT_dCA1_allsess(:,3),24,'k','filled','o'); hold on;
scatter3(class_AD_dCA1_allsess(:,1),class_AD_dCA1_allsess(:,2),class_AD_dCA1_allsess(:,3),24,'r','filled','o'); hold on;
scatter3(class_WT_vCA1_allsess(:,1),class_WT_vCA1_allsess(:,2),class_WT_vCA1_allsess(:,3),24,[.8 .8 .8],'filled','o'); hold on;
scatter3(class_AD_vCA1_allsess(:,1),class_AD_vCA1_allsess(:,2),class_AD_vCA1_allsess(:,3),24,[1 .8 .8],'filled','o'); hold on;
xlabel('Entropy (bits/spike)'); zlabel('KLD'); ylabel('Correlation'); set(gca,'zscale','log'); 
xlim([200 800]); zlim([0 .06]); ylim([0 .12]); 
subplot(3,3,[8,9]); 
boxplot([allsess_log_null_all;allsess_log_all;allsess_log_WT_DV_null;allsess_log_WT_DV;allsess_log_AD_DV_null;allsess_log_AD_DV],[zeros(100,1);ones(100,1);3*ones(100,1);4*ones(100,1);6*ones(100,1);7*ones(100,1)],'PlotSTyle','traditional','Colors','ckckck','MedianStyle','target','Symbol','','Labels',{'null' '4-way' 'null' 'WT D vs. V' 'null' 'AD D vs. V'});
ylim([0 .8]); ylabel('Accuracy'); 
% cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 9'); 
% saveas(f,'Integrated','fig'); orient(f,'landscape'); print('Integrated.pdf','-bestfit','-painters','-dpdf'); close(f); 
%% Calculate power in SWR frequency band
tic
gauss_kernel = fspecial('gaussian',[1 32*30],4*30); %make a gaussian kernel with size 32ms and SD 4ms
for VR_ephys = 26:27 %:size(ks_all_data,1)
    cd(strcat('Y:\Udaysankar Chockanathan\Aim 1\kilosort_results_all\',ks_all_data{VR_ephys,5},'\preAutoMerge'));  
    fileID = fopen(strcat(ks_all_data{VR_ephys,5},'_binary.dat'));    
    working_total = zeros(1,ks_all_data{VR_ephys,11}); 
    for channel = 1:128
        frewind(fileID); 
        if channel > 1 %this reads in the first (channel-1) values to offset the reading appropriately            
            tmp = fread(fileID,[(channel-1),1],'int16'); 
        end
        skip = 2*127; %specify the skip in bytes (int16 has 2 bytes per value)
        precision = '1*int16'; %specify the precision (so the program will read ONE value and skip 128 values)
        bin_data = fread(fileID,[1,ks_all_data{VR_ephys,11}],precision,skip); %load in the current channel raw data (will not be the same as the corresponding amplifier_data row because it gets remapped in spike_sorting_pipeline1_1
        to_add = bandpass(bin_data,[150 250],30000).^2;
        working_total = working_total+to_add;
        VR_ephys
        channel
    end
    smoothed_total = conv(gauss_kernel,working_total);
    smoothed_total = smoothed_total(numel(gauss_kernel):end);
    sqrt_smoothed_total = sqrt(smoothed_total);
    ks_all_data{VR_ephys,36} = sqrt_smoothed_total;
    fclose(fileID);
end
toc
%% Identify SWR epochs
SWR_epochs = cell(size(ks_all_data,1),1); 
for VR_ephys = 1:size(ks_all_data,1) %identify SWR start/end times
    is_SWR = ks_all_data{VR_ephys,36} > (mean(ks_all_data{VR_ephys,36})+2*std(ks_all_data{VR_ephys,36})); 
    tt_to_SWR = find(diff(is_SWR) > 0)+1;
    tt_from_SWR = find(diff(is_SWR) < 0)+1;     
    if and(is_SWR(1) == 1,is_SWR(end) == 0) %if it starts, but doesn't end, on a SWR epoch
        tt_from_SWR(1) = [];
    elseif and(is_SWR(1) == 0,is_SWR(end) == 1) %if it ends, but doesn't start, on a SWR epoch
        tt_to_SWR(end) = []; 
    elseif and(is_SWR(1) == 1,is_SWR(end) == 1) % if it starts and ends on a SWR epoch        
        tt_from_SWR(1) = []; %remove the first "from SWR" transition time
        tt_to_SWR(end) = []; %remove the last "to SWR" transition time
    end
    tmp_SWR_epochs = [tt_to_SWR',tt_from_SWR'];
    short_SWR = diff(tmp_SWR_epochs,1,2) < 15*30; %identify those putative SWRs that are less than 15ms long
    tmp_SWR_epochs(short_SWR,:) = []; %eliminate those short SWRs
    SWR_start_times = round(tmp_SWR_epochs(:,1)/30000); %pull together the SWR start times and convert to seconds
    SWR_start_times(SWR_start_times == 0) = []; 
    SWR_vel = mod_vel{VR_ephys}(SWR_start_times); %the velocity at which the mouse was moving at each SWR start time
    tmp_SWR_epochs(SWR_vel ~= 0,:) = []; %eliminate those SWR epochs when the mouse was moving (nonzero velocity)
    SWR_epochs{VR_ephys} = tmp_SWR_epochs;
    VR_ephys
end
%% Discard SWRs that occur too closely to prior SWRs
SWR_unamb_epochs = cell(size(ks_all_data,1),1); 
min_sep = 30000; %minimum separation between two SWRs
for VR_ephys = 1:size(ks_all_data,1)
    SWR_unamb_epochs{VR_ephys}(1,:) = SWR_epochs{VR_ephys}(1,:); 
    new_idx = 1;
    for i = 2:(size(SWR_epochs{VR_ephys},1)-1)
        if SWR_epochs{VR_ephys}(i,1)-SWR_epochs{VR_ephys}((i-1),2) > min_sep
            new_idx = new_idx+1; 
            SWR_unamb_epochs{VR_ephys}(new_idx,:) = SWR_epochs{VR_ephys}(i,:); 
        end
    end    
end
%% Plot representative epochs of SWRs
%plot the SWR start times vs. SWR durations to identify promising 2s intervals
% VR_ephys = [16,1,17,5];
% leveler = repmat([1:size(SWR_epochs{VR_ephys,1},1)]',1,2); leveler = leveler'; leveler = leveler(:); 
% etp = SWR_epochs{VR_ephys,1}'; etp = etp(:)./30000;
% epoch_length = diff(SWR_epochs{VR_ephys,1},1,2); 
% figure; scatter(SWR_epochs{VR_ephys}(:,1)/30000,epoch_length/30000)
%% pull out the raw data for those intervals and band pass
epochs_to_plot = [74.4,78.5,208.8,290];
interval_to_plot = 3*30000;
SWR_to_plot = zeros(4,interval_to_plot); 
raw_to_plot = SWR_to_plot; 
idx = 0; 
for VR_ephys = [16,1,17,5]
    idx = idx+1; 
    cd(strcat('Y:\Udaysankar Chockanathan\Aim 1\kilosort_results_all\',ks_all_data{VR_ephys,5},'\preAutoMerge'));  
    fileID = fopen(strcat(ks_all_data{VR_ephys,5},'_binary.dat'));       
    max_SWR_log = 0; 
    SWR_filter = []; 
    raw_filter = []; 
    for channel = 1:128
        frewind(fileID); 
        fread(fileID,[128*epochs_to_plot(idx)*30000-1,1],'int16'); %this reads in all the values from all channels before the epochs_to_plot
        if channel > 1 %this reads in the first (channel-1) values to offset the reading appropriately            
            tmp = fread(fileID,[(channel-1),1],'int16'); 
        end
        skip = 2*127; %specify the skip in bytes (int16 has 2 bytes per value)
        precision = '1*int16'; %specify the precision (so the program will read ONE value and skip 128 values)
        bin_data = fread(fileID,[1,interval_to_plot],precision,skip); %load in the current channel raw data (will not be the same as the corresponding amplifier_data row because it gets remapped in spike_sorting_pipeline1_1
        tmp_SWR_filter = bandpass(bin_data,[150 250],30000); %filter bin_data at the SWR band
        tmp_raw_filter = lowpass(bin_data,400,30000); %filter bin_data at the 1-400 band
%       identify channel that has the highest SWR power
        if max(tmp_SWR_filter.^2) > max_SWR_log %if the max SWR power of the current channel is larger than any before it, replace the stored values
            max_SWR_log = max(tmp_SWR_filter);
            SWR_filter = tmp_SWR_filter;
            raw_filter = tmp_raw_filter;             
        end        
        VR_ephys
        channel
    end
    SWR_to_plot(idx,:) = SWR_filter; 
    raw_to_plot(idx,:) = tmp_raw_filter; 
end
%% plot the raw and band passed ephys data and highlight the SWR intervals
idx = 0;
for VR_ephys = 16 %[16,1,17,5]
    idx = idx+1; 
    SWR_idx = find(and(SWR_epochs{VR_ephys}(:,1) > epochs_to_plot(idx)*30000,SWR_epochs{VR_ephys}(:,1) < epochs_to_plot(idx)*30000+interval_to_plot)); %find rows of SWR_epochs that are within the plotted interval
    f = figure; set(f,'Position',[0 0 1000 500]); 
    subplot(2,1,1); plot(10/30000:10/30000:3,SWR_to_plot(idx,1:10:end),'-k');
    subplot(2,1,2); plot(10/30000:10/30000:3,raw_to_plot(idx,1:10:end),'-k'); hold on;
    for i = 1:numel(SWR_idx)
        plot((SWR_epochs{VR_ephys}(SWR_idx(i),:)-epochs_to_plot(idx)*30000)/30000,[20,20],'LineWidth',2); hold on;
    end
end
xlim([0 3]); xlabel('Time (s)'); ylabel('Voltage (\muV)'); 
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 10'); 
saveas(gcf,'SWR_schematic_WT_dCA1','fig'); orient(f,'landscape'); print('SWR_schematic_WT_dCA1.pdf','-bestfit','-painters','-dpdf'); close(gcf); 
%% Plot SWR abundance and duration
SWR_duration = zeros(size(ks_all_data,1),1); 
for VR_ephys = 1:size(ks_all_data,1)
    SWR_duration(VR_ephys) = mean(diff(SWR_epochs{VR_ephys},1,2))/30000;
end
SWR_rate = (cellfun(@numel,SWR_epochs)/2)./(mod_durat/1000); % How often do SWRs happen in each group
f = figure; set(f,'Position',[0 0 2000 1000]); 
subplot(2,4,1); 
scatter(ones(numel(WT_dCA1_sites),1),SWR_rate(WT_dCA1_sites),72,'k','filled'); hold on;
scatter(3*ones(numel(WT_vCA1_sites),1),SWR_rate(WT_vCA1_sites),72,[.8 .8 .8],'filled');
[p,~] = ranksum(SWR_rate(WT_dCA1_sites),SWR_rate(WT_vCA1_sites)); text(2,1.5,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 2]); xticks([1 3]); xticklabels({'dCA1','vCA1'}); ylabel('SWR rate (events/s)'); 
subplot(2,4,2); 
scatter(ones(numel(WT_sites),1),SWR_rate(WT_sites),72,'k','filled'); hold on;
scatter(3*ones(numel(AD_sites),1),SWR_rate(AD_sites),72,'k','filled');
[p,~] = ranksum(SWR_rate(WT_sites),SWR_rate(AD_sites)); text(2,1.5,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 2]); xticks([1 3]); xticklabels({'WT','AD'}); ylabel('SWR rate (events/s)'); 
subplot(2,4,3); 
scatter(ones(numel(WT_dCA1_sites),1),SWR_rate(WT_dCA1_sites),72,'k','filled'); hold on;
scatter(3*ones(numel(AD_dCA1_sites),1),SWR_rate(AD_dCA1_sites),72,'r','filled');
[p,~] = ranksum(SWR_rate(WT_dCA1_sites),SWR_rate(AD_dCA1_sites)); text(2,1.5,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 2]); xticks([1 3]); xticklabels({'WT dCA1','AD dCA1'}); ylabel('SWR rate (events/s)'); 
subplot(2,4,4); 
scatter(ones(numel(WT_vCA1_sites),1),SWR_rate(WT_vCA1_sites),72,[.8 .8 .8],'filled'); hold on;
scatter(3*ones(numel(AD_vCA1_sites),1),SWR_rate(AD_vCA1_sites),72,[1 .8 .8],'filled');
[p,~] = ranksum(SWR_rate(WT_vCA1_sites),SWR_rate(AD_vCA1_sites)); text(2,1.5,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 2]); xticks([1 3]); xticklabels({'WT vCA1','AD vCA1'}); ylabel('SWR rate (events/s)'); 
subplot(2,4,5); 
scatter(ones(numel(WT_dCA1_sites),1),SWR_duration(WT_dCA1_sites),72,'k','filled'); hold on;
scatter(3*ones(numel(WT_vCA1_sites),1),SWR_duration(WT_vCA1_sites),72,[.8 .8 .8],'filled');
[p,~] = ranksum(SWR_duration(WT_dCA1_sites),SWR_duration(WT_vCA1_sites)); text(2,.1,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 .2]); xticks([1 3]); xticklabels({'dCA1','vCA1'}); ylabel('SWR duration (s)'); 
subplot(2,4,6); 
scatter(ones(numel(WT_sites),1),SWR_duration(WT_sites),72,'k','filled'); hold on;
scatter(3*ones(numel(AD_sites),1),SWR_duration(AD_sites),72,'k','filled');
[p,~] = ranksum(SWR_duration(WT_sites),SWR_duration(AD_sites)); text(2,.1,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 .2]); xticks([1 3]); xticklabels({'WT','AD'}); ylabel('SWR duration (s)'); 
subplot(2,4,7); 
scatter(ones(numel(WT_dCA1_sites),1),SWR_duration(WT_dCA1_sites),72,'k','filled'); hold on;
scatter(3*ones(numel(AD_dCA1_sites),1),SWR_duration(AD_dCA1_sites),72,'r','filled');
[p,~] = ranksum(SWR_duration(WT_dCA1_sites),SWR_duration(AD_dCA1_sites)); text(2,.1,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 .2]); xticks([1 3]); xticklabels({'WT dCA1','AD dCA1'}); ylabel('SWR duration (s)'); 
subplot(2,4,8); 
scatter(ones(numel(WT_vCA1_sites),1),SWR_duration(WT_vCA1_sites),72,[.8 .8 .8],'filled'); hold on;
scatter(3*ones(numel(AD_vCA1_sites),1),SWR_duration(AD_vCA1_sites),72,[1 .8 .8],'filled');
[p,~] = ranksum(SWR_duration(WT_vCA1_sites),SWR_duration(AD_vCA1_sites)); text(2,.1,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 .2]); xticks([1 3]); xticklabels({'WT vCA1','AD vCA1'}); ylabel('SWR duration (s)'); 
% cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 10'); 
% saveas(f,'SWR_abundance_and_duration','fig'); orient(f,'landscape'); print('SWR_abundance_and_duration.pdf','-bestfit','-painters','-dpdf'); close(f); 
%% Generate PSTH data for SWRs
num_shuffles = 500; %CONSIDER CHANGING TO 5000
mod_index = cell(size(ks_all_data,1),1); 
mod_index_p_value = mod_index; 
mod_directionality = mod_index; 
for VR_ephys = 1:size(ks_all_data,1)
    for u = 1:size(ks_all_data{VR_ephys,1},1) %for each unit
        spikes_for_PSTH = cell(size(SWR_unamb_epochs{VR_ephys},1),1); 
        for s = 1:size(SWR_unamb_epochs{VR_ephys},1) %for each SWR
            soi = round(SWR_unamb_epochs{VR_ephys}(s,1)/30); %current SWR start time (in ms)
            start_PSTH = soi-500; %start of the PSTH window for the current SWR
            end_PSTH = soi+500; %end of the PSTH window for the current SWR
            all_spikes_from_uoi = ks_all_data{VR_ephys,1}{u,3}; %all spikes from current unit            
            tmp_spikes = all_spikes_from_uoi(and(all_spikes_from_uoi > start_PSTH,all_spikes_from_uoi < end_PSTH)); %only the spikes that occured in the PSTH interval
            spikes_for_PSTH{s} = tmp_spikes-soi;
        end
        [N_empirical,~] = histcounts(cell2mat(spikes_for_PSTH'),[-.5:201:200.5]); 
        [N_background,~] = histcounts(cell2mat(spikes_for_PSTH'),[-500.5:500:-.5]); 
        FR_empirical = N_empirical/201;
        FR_background = N_background/500;
        if FR_empirical > FR_background        
            mod_directionality{VR_ephys}(u) = 1;
        elseif FR_empirical < FR_background
            mod_directionality{VR_ephys}(u) = -1; 
        elseif FR_empirical == FR_background
            mod_directionality{VR_ephys}(u) = 0; 
        end
        N_shuffled = zeros(num_shuffles,1); 
        circ_shuffle_spikes = cell(size(SWR_unamb_epochs{VR_ephys},1),1); 
        for shuffle = 1:num_shuffles %for each shuffle
            parfor s = 1:size(SWR_unamb_epochs{VR_ephys},1) %for each SWR
                tmp_shuffle_spikes = spikes_for_PSTH{s}+randi(1000); %add a random number to each spike in the current SWR
                tmp_shuffle_spikes(tmp_shuffle_spikes > 500) = tmp_shuffle_spikes(tmp_shuffle_spikes > 500)-1000; %make the out-of-bounds spikes wrap around
                circ_shuffle_spikes{s} = tmp_shuffle_spikes; %store this shuffled spike train
            end
            [N_shuffled(shuffle),~] = histcounts(cell2mat(circ_shuffle_spikes'),[-.5:201:200.5]);            
        end
        mod_index_null = (N_shuffled-repmat(mean(N_shuffled),num_shuffles,1)).^2; %null distribution of modulation index
        mod_index{VR_ephys}(u) = (N_empirical-mean(N_shuffled)).^2; %empirical modulation index
        mod_index_p_value{VR_ephys}(u) = mean(mod_index{VR_ephys}(u) > mod_index_null); %modulation index p-value    
        VR_ephys
        u
    end    
end
%% Plot single trial PSTH for a representative positively-, negatively-, and non-modulated cell
gauss_kernel = fspecial('gaussian',[1 80],10); %make a gaussian kernel with size 80ms and SD 10ms
VR_ephys = [18,2,16]; 
uoi = [130,76,40]; 
idx = 0;
ttp = 125; %trials to plot
f = figure; set(f,'Position',[0 0 1000 500]); 
for i = 1:3    
    spikes_for_PSTH = cell(size(SWR_unamb_epochs{VR_ephys(i)},1),1); 
    rowtrack = 0;
    subplot(2,3,i);     
    for s = 1:size(SWR_unamb_epochs{VR_ephys(i)},1) %for each SWR
        soi = round(SWR_unamb_epochs{VR_ephys(i)}(s,1)/30); %current SWR start time (in ms)
        start_PSTH = soi-500; %start of the PSTH window for the current SWR
        end_PSTH = soi+500; %end of the PSTH window for the current SWR
        all_spikes_from_uoi = ks_all_data{VR_ephys(i),1}{uoi(i),3}; %all spikes from current unit            
        tmp_spikes = all_spikes_from_uoi(and(all_spikes_from_uoi > start_PSTH,all_spikes_from_uoi < end_PSTH)); %only the spikes that occured in the PSTH interval
        spikes_for_PSTH{s} = tmp_spikes-soi;           
        if and(isempty(spikes_for_PSTH{s}) == 0,rowtrack <= ttp)
            rowtrack = rowtrack+1; 
            rasterPlotter(spikes_for_PSTH{s},rowtrack,'k'); hold on;
        end
    end
    plot([0 0],[0 ttp],'Color',[0 .5 0],'LineWidth',0.5); ylim([0 ttp]); xlabel('Time (ms)'); ylabel('SWR events'); 
    [N,~] = histcounts(cell2mat(spikes_for_PSTH'),[-499.5:1:499.5]); 
    N = 1000*N./size(SWR_unamb_epochs{VR_ephys(i)},1); %convert from spikes per bin to spikes per s
    N = [repmat(N(1),1,100),N,repmat(N(end),1,100)]; %padding either side to minimize edge effects
    N_smoothed = conv(gauss_kernel,N);     
    subplot(2,3,i+3); plot(-499:1:499,N_smoothed(140:1138),'k','LineWidth',2); hold on; xlabel('Time (ms)'); ylabel('FR (Hz)'); 
    plot([0 0],[0 ttp],'Color',[0 .5 0],'LineWidth',0.5); ylim([0 max(N_smoothed+5)]); 
    text(300,6,strcat('p = ',num2str(1-mod_index_p_value{VR_ephys(i)}(uoi(i))))); 
    text(300,7,strcat('mod index = ',num2str(mod_index{VR_ephys(i)}(uoi(i))))); 
end
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 10'); 
saveas(f,'single_trial_SWR_PSTH','fig'); orient(f,'landscape'); print('single_trial_SWR_PSTH.pdf','-bestfit','-painters','-dpdf'); close(f); 
%% Generate smoothed/normalized PSTH for all cells
z_score_stored = cell(size(ks_all_data,1),1); 
for VR_ephys = 1:size(ks_all_data,1)     
    for u = 1:size(ks_all_data{VR_ephys,1},1) %for each unit
        spikes_for_PSTH = cell(size(SWR_unamb_epochs{VR_ephys},1),1); 
        for s = 1:size(SWR_unamb_epochs{VR_ephys},1) %for each SWR
            soi = round(SWR_unamb_epochs{VR_ephys}(s,1)/30); %current SWR start time (in ms)
            start_PSTH = soi-500; %start of the PSTH window for the current SWR
            end_PSTH = soi+500; %end of the PSTH window for the current SWR
            all_spikes_from_uoi = ks_all_data{VR_ephys,1}{u,3}; %all spikes from current unit            
            tmp_spikes = all_spikes_from_uoi(and(all_spikes_from_uoi > start_PSTH,all_spikes_from_uoi < end_PSTH)); %only the spikes that occured in the PSTH interval
            spikes_for_PSTH{s} = tmp_spikes-soi;
        end
        [N,~] = histcounts(cell2mat(spikes_for_PSTH'),[-499.5:1:499.5]); 
        N = 1000*N./size(SWR_unamb_epochs{VR_ephys},1); %convert from spikes per bin to spikes per s
        N = [repmat(N(1),1,100),N,repmat(N(end),1,100)]; %padding either side to minimize edge effects
        N_smoothed = conv(gauss_kernel,N); %smooth
        N_smoothed = N_smoothed(160:1118); %unpad - was 140:1138        
        tmpz = zscore(N_smoothed);
        tmpz = tmpz-mean(tmpz(1:479));               
        z_score_stored{VR_ephys,1}(u,:) = tmpz; %zscored        
    end    
    VR_ephys
end
%% Calculate PSTH kinetics
FWHM = cell(size(ks_all_data,1),1); 
lag = cell(size(ks_all_data,1),1); 
for VR_ephys = 1:size(ks_all_data,1)
    for u = 1:size(ks_all_data{VR_ephys,1},1) %for each unit
        if and(mod_directionality{VR_ephys}(u) >= 0,mod_index_p_value{VR_ephys}(u) > 0.95)
            [tmp_extremum,idx_extremum] = max(z_score_stored{VR_ephys,1}(u,481:end));
            halfmax = tmp_extremum/2; 
            bin_hm = z_score_stored{VR_ephys,1}(u,:) >= halfmax; %binary vector indicating where the function is greater than the half-max   
            left_index = find(diff(bin_hm(1:(idx_extremum+480))) == 1,1,'last'); %last index before extremum where the function crosses half-max
            right_index = find(diff(bin_hm((idx_extremum+480):end)) == -1,1,'first')+idx_extremum+479; %first index after extremum where the function crosses half-max
            lag{VR_ephys}(u) = idx_extremum; 
            if isempty(right_index) == 0
                FWHM{VR_ephys}(u) = right_index-left_index;
            else
                FWHM{VR_ephys}(u) = NaN;
            end
        elseif and(mod_directionality{VR_ephys}(u) < 0,mod_index_p_value{VR_ephys}(u) > 0.95)
            [tmp_extremum,idx_extremum] = min(z_score_stored{VR_ephys,1}(u,481:end));
            halfmax = tmp_extremum/2; 
            bin_hm = z_score_stored{VR_ephys,1}(u,:) <= halfmax; %binary vector indicating where the function is greater than the half-max   
            left_index = find(diff(bin_hm(1:(idx_extremum+480))) == 1,1,'last'); %last index before extremum where the function crosses half-max
            right_index = find(diff(bin_hm((idx_extremum+480):end)) == -1,1,'first')+idx_extremum+479; %first index after extremum where the function crosses half-max          
            lag{VR_ephys}(u) = idx_extremum; 
            if isempty(right_index) == 0
                FWHM{VR_ephys}(u) = right_index-left_index;
            else
                FWHM{VR_ephys}(u) = NaN;
            end
        elseif mod_index_p_value{VR_ephys}(u) <= 0.95
            FWHM{VR_ephys}(u) = NaN;
            lag{VR_ephys}(u) = NaN;
        end
    end
end
%% Collate and plot all the PSTHs belonging to the same strain/group
WT_dCA1_PETH = cell2mat(z_score_stored(WT_dCA1_sites)); WT_dCA1_mod_index = cell2mat(mod_index(WT_dCA1_sites)'); WT_dCA1_mod_index_p_value = cell2mat(mod_index_p_value(WT_dCA1_sites)'); WT_dCA1_mod_directionality = cell2mat(mod_directionality(WT_dCA1_sites)'); WT_dCA1_lag = cell2mat(lag(WT_dCA1_sites)'); WT_dCA1_FWHM = cell2mat(FWHM(WT_dCA1_sites)'); 
AD_dCA1_PETH = cell2mat(z_score_stored(AD_dCA1_sites)); AD_dCA1_mod_index = cell2mat(mod_index(AD_dCA1_sites)'); AD_dCA1_mod_index_p_value = cell2mat(mod_index_p_value(AD_dCA1_sites)'); AD_dCA1_mod_directionality = cell2mat(mod_directionality(AD_dCA1_sites)'); AD_dCA1_lag = cell2mat(lag(AD_dCA1_sites)'); AD_dCA1_FWHM = cell2mat(FWHM(AD_dCA1_sites)'); 
WT_vCA1_PETH = cell2mat(z_score_stored(WT_vCA1_sites)); WT_vCA1_mod_index = cell2mat(mod_index(WT_vCA1_sites)'); WT_vCA1_mod_index_p_value = cell2mat(mod_index_p_value(WT_vCA1_sites)'); WT_vCA1_mod_directionality = cell2mat(mod_directionality(WT_vCA1_sites)'); WT_vCA1_lag = cell2mat(lag(WT_vCA1_sites)'); WT_vCA1_FWHM = cell2mat(FWHM(WT_vCA1_sites)'); 
AD_vCA1_PETH = cell2mat(z_score_stored(AD_vCA1_sites)); AD_vCA1_mod_index = cell2mat(mod_index(AD_vCA1_sites)'); AD_vCA1_mod_index_p_value = cell2mat(mod_index_p_value(AD_vCA1_sites)'); AD_vCA1_mod_directionality = cell2mat(mod_directionality(AD_vCA1_sites)'); AD_vCA1_lag = cell2mat(lag(AD_vCA1_sites)'); AD_vCA1_FWHM = cell2mat(FWHM(AD_vCA1_sites)'); 
% put all the cells that are positively/negatively/unmodulated together 
WT_dCA1_PETH_pos_mod = WT_dCA1_PETH(and(WT_dCA1_mod_directionality > 0,WT_dCA1_mod_index_p_value > 0.95),:); 
WT_dCA1_PETH_pos_mod_index = WT_dCA1_mod_index(and(WT_dCA1_mod_directionality > 0,WT_dCA1_mod_index_p_value > 0.95)); 
WT_dCA1_PETH_neg_mod = WT_dCA1_PETH(and(WT_dCA1_mod_directionality < 0,WT_dCA1_mod_index_p_value > 0.95),:); 
WT_dCA1_PETH_neg_mod_index = WT_dCA1_mod_index(and(WT_dCA1_mod_directionality < 0,WT_dCA1_mod_index_p_value > 0.95)); 
WT_dCA1_PETH_non_mod = WT_dCA1_PETH(WT_dCA1_mod_index_p_value <= 0.95,:); 
WT_dCA1_PETH_non_mod_index = WT_dCA1_mod_index(WT_dCA1_mod_index_p_value <= 0.95); 
AD_dCA1_PETH_pos_mod = AD_dCA1_PETH(and(AD_dCA1_mod_directionality > 0,AD_dCA1_mod_index_p_value > 0.95),:); 
AD_dCA1_PETH_pos_mod_index = AD_dCA1_mod_index(and(AD_dCA1_mod_directionality > 0,AD_dCA1_mod_index_p_value > 0.95)); 
AD_dCA1_PETH_neg_mod = AD_dCA1_PETH(and(AD_dCA1_mod_directionality < 0,AD_dCA1_mod_index_p_value > 0.95),:); 
AD_dCA1_PETH_neg_mod_index = AD_dCA1_mod_index(and(AD_dCA1_mod_directionality < 0,AD_dCA1_mod_index_p_value > 0.95)); 
AD_dCA1_PETH_non_mod = AD_dCA1_PETH(AD_dCA1_mod_index_p_value <= 0.95,:); 
AD_dCA1_PETH_non_mod_index = AD_dCA1_mod_index(AD_dCA1_mod_index_p_value <= 0.95); 
WT_vCA1_PETH_pos_mod = WT_vCA1_PETH(and(WT_vCA1_mod_directionality > 0,WT_vCA1_mod_index_p_value > 0.95),:); 
WT_vCA1_PETH_pos_mod_index = WT_vCA1_mod_index(and(WT_vCA1_mod_directionality > 0,WT_vCA1_mod_index_p_value > 0.95)); 
WT_vCA1_PETH_neg_mod = WT_vCA1_PETH(and(WT_vCA1_mod_directionality < 0,WT_vCA1_mod_index_p_value > 0.95),:); 
WT_vCA1_PETH_neg_mod_index = WT_vCA1_mod_index(and(WT_vCA1_mod_directionality < 0,WT_vCA1_mod_index_p_value > 0.95)); 
WT_vCA1_PETH_non_mod = WT_vCA1_PETH(WT_vCA1_mod_index_p_value <= 0.95,:); 
WT_vCA1_PETH_non_mod_index = WT_vCA1_mod_index(WT_vCA1_mod_index_p_value <= 0.95); 
AD_vCA1_PETH_pos_mod = AD_vCA1_PETH(and(AD_vCA1_mod_directionality > 0,AD_vCA1_mod_index_p_value > 0.95),:); 
AD_vCA1_PETH_pos_mod_index = AD_vCA1_mod_index(and(AD_vCA1_mod_directionality > 0,AD_vCA1_mod_index_p_value > 0.95)); 
AD_vCA1_PETH_neg_mod = AD_vCA1_PETH(and(AD_vCA1_mod_directionality < 0,AD_vCA1_mod_index_p_value > 0.95),:); 
AD_vCA1_PETH_neg_mod_index = AD_vCA1_mod_index(and(AD_vCA1_mod_directionality < 0,AD_vCA1_mod_index_p_value > 0.95)); 
AD_vCA1_PETH_non_mod = AD_vCA1_PETH(AD_vCA1_mod_index_p_value <= 0.95,:); 
AD_vCA1_PETH_non_mod_index = AD_vCA1_mod_index(AD_vCA1_mod_index_p_value <= 0.95); 
% sort the ones that are positively modulated by greatest to least mod index
[~,WT_dCA1_PETH_pos_mod_sorted_idx] = sort(WT_dCA1_PETH_pos_mod_index,'descend'); 
[~,AD_dCA1_PETH_pos_mod_sorted_idx] = sort(AD_dCA1_PETH_pos_mod_index,'descend'); 
[~,WT_vCA1_PETH_pos_mod_sorted_idx] = sort(WT_vCA1_PETH_pos_mod_index,'descend'); 
[~,AD_vCA1_PETH_pos_mod_sorted_idx] = sort(AD_vCA1_PETH_pos_mod_index,'descend'); 
% sort those that are negatively modulated by least to greatest mod_index 
[~,WT_dCA1_PETH_neg_mod_sorted_idx] = sort(WT_dCA1_PETH_neg_mod_index,'ascend'); 
[~,AD_dCA1_PETH_neg_mod_sorted_idx] = sort(AD_dCA1_PETH_neg_mod_index,'ascend'); 
[~,WT_vCA1_PETH_neg_mod_sorted_idx] = sort(WT_vCA1_PETH_neg_mod_index,'ascend'); 
[~,AD_vCA1_PETH_neg_mod_sorted_idx] = sort(AD_vCA1_PETH_neg_mod_index,'ascend'); 
% shuffle around those that are non modulated
WT_dCA1_PETH_non_mod_sorted_idx = randperm(size(WT_dCA1_PETH_non_mod,1));
AD_dCA1_PETH_non_mod_sorted_idx = randperm(size(AD_dCA1_PETH_non_mod,1));
WT_vCA1_PETH_non_mod_sorted_idx = randperm(size(WT_vCA1_PETH_non_mod,1));
AD_vCA1_PETH_non_mod_sorted_idx = randperm(size(AD_vCA1_PETH_non_mod,1));
% plot
numColors = 128; % Define the number of colors for the colormap
% Generate RGB values for blue to white to red
blueToRed = [linspace(0, 1, numColors/2)', linspace(0, 1, numColors/2)', linspace(1, 1, numColors/2)'];
whiteColor = ones(1,3); % Add white color, repeat 4 times to fill the gap
redToBlue = [linspace(1, 1, numColors/2)', linspace(1, 0, numColors/2)', linspace(1, 0, numColors/2)'];
% Concatenate the colors to form the custom colormap
customColormap = [blueToRed; whiteColor; redToBlue];
f = figure; set(f,'Position',[0 0 2000 700]); 
xaxis = -499:499; 
subplot(1,4,1); 
yaxis = 1:(numel(WT_dCA1_PETH_pos_mod_index)+numel(WT_dCA1_PETH_neg_mod_index)+size(WT_dCA1_PETH_non_mod,1)); 
imagesc(xaxis,yaxis,[WT_dCA1_PETH_pos_mod(WT_dCA1_PETH_pos_mod_sorted_idx,:);WT_dCA1_PETH_neg_mod(WT_dCA1_PETH_neg_mod_sorted_idx,:);WT_dCA1_PETH_non_mod(WT_dCA1_PETH_non_mod_sorted_idx,:)],[-6 6]); 
colormap(customColormap); colorbar('southoutside'); hold on; plot([0 0],[0 max(yaxis)+1],'k','LineWidth',1); xlabel('Time (ms)'); ylabel('Units'); yticks([]); title('WT dCA1'); 
subplot(1,4,2); 
yaxis = 1:(numel(AD_dCA1_PETH_pos_mod_index)+numel(AD_dCA1_PETH_neg_mod_index)+size(AD_dCA1_PETH_non_mod,1)); 
imagesc(xaxis,yaxis,[AD_dCA1_PETH_pos_mod(AD_dCA1_PETH_pos_mod_sorted_idx,:);AD_dCA1_PETH_neg_mod(AD_dCA1_PETH_neg_mod_sorted_idx,:);AD_dCA1_PETH_non_mod(AD_dCA1_PETH_non_mod_sorted_idx,:)],[-6 6]); 
colormap(customColormap); colorbar('southoutside'); hold on; plot([0 0],[0 max(yaxis)+1],'k','LineWidth',1); xlabel('Time (ms)'); ylabel('Units'); yticks([]); title('AD dCA1'); 
subplot(1,4,3); 
yaxis = 1:(numel(WT_vCA1_PETH_pos_mod_index)+numel(WT_vCA1_PETH_neg_mod_index)+size(WT_vCA1_PETH_non_mod,1)); 
imagesc(xaxis,yaxis,[WT_vCA1_PETH_pos_mod(WT_vCA1_PETH_pos_mod_sorted_idx,:);WT_vCA1_PETH_neg_mod(WT_vCA1_PETH_neg_mod_sorted_idx,:);WT_vCA1_PETH_non_mod(WT_vCA1_PETH_non_mod_sorted_idx,:)],[-6 6]); 
colormap(customColormap); colorbar('southoutside'); hold on; plot([0 0],[0 max(yaxis)+1],'k','LineWidth',1); xlabel('Time (ms)'); ylabel('Units'); yticks([]); title('WT vCA1'); 
subplot(1,4,4); 
yaxis = 1:(numel(AD_vCA1_PETH_pos_mod_index)+numel(AD_vCA1_PETH_neg_mod_index)+size(AD_vCA1_PETH_non_mod,1)); 
imagesc(xaxis,yaxis,[AD_vCA1_PETH_pos_mod(AD_vCA1_PETH_pos_mod_sorted_idx,:);AD_vCA1_PETH_neg_mod(AD_vCA1_PETH_neg_mod_sorted_idx,:);AD_vCA1_PETH_non_mod(AD_vCA1_PETH_non_mod_sorted_idx,:)],[-6 6]); 
colormap(customColormap); colorbar('southoutside'); hold on; plot([0 0],[0 max(yaxis)+1],'k','LineWidth',1); xlabel('Time (ms)'); ylabel('Units'); yticks([]); title('AD vCA1'); 
% cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 10'); 
% saveas(f,'all_SWR_PSTH','fig'); orient(f,'landscape'); print('all_SWR_PSTH.pdf','-bestfit','-painters','-dpdf'); close(f); 
%% Plot SWR-cell modulation fraction data 
%pie chart of how many cells are positively/negatively/non-modulated
g = figure; set(g,'Position',[0 0 1500 2000]); 
WT_dCA1_pos_mod = 0;
WT_dCA1_neg_mod = 0;
WT_vCA1_pos_mod = 0;
WT_vCA1_neg_mod = 0;
AD_dCA1_pos_mod = 0;
AD_dCA1_neg_mod = 0;
AD_vCA1_pos_mod = 0; 
AD_vCA1_neg_mod = 0; 
WT_dCA1_pos_mod_pyramidal = 0;
WT_dCA1_neg_mod_pyramidal = 0;
WT_vCA1_pos_mod_pyramidal = 0;
WT_vCA1_neg_mod_pyramidal = 0;
AD_dCA1_pos_mod_pyramidal = 0;
AD_dCA1_neg_mod_pyramidal = 0;
AD_vCA1_pos_mod_pyramidal = 0; 
AD_vCA1_neg_mod_pyramidal = 0; 
WT_dCA1_pos_mod_interneuron = 0;
WT_dCA1_neg_mod_interneuron = 0;
WT_vCA1_pos_mod_interneuron = 0;
WT_vCA1_neg_mod_interneuron = 0;
AD_dCA1_pos_mod_interneuron = 0;
AD_dCA1_neg_mod_interneuron = 0;
AD_vCA1_pos_mod_interneuron = 0; 
AD_vCA1_neg_mod_interneuron = 0; 
for VR_ephys = 1:size(ks_all_data,1)
    if any(VR_ephys == WT_dCA1_sites) == 1
        WT_dCA1_pos_mod = WT_dCA1_pos_mod+sum(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == 1));
        WT_dCA1_neg_mod = WT_dCA1_neg_mod+sum(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == -1));
        WT_dCA1_pos_mod_pyramidal = WT_dCA1_pos_mod_pyramidal+sum(all([mod_index_p_value{VR_ephys} > 0.95;mod_directionality{VR_ephys} == 1;logical(ks_all_data{VR_ephys,42}')],1));
        WT_dCA1_neg_mod_pyramidal = WT_dCA1_neg_mod_pyramidal+sum(all([mod_index_p_value{VR_ephys} > 0.95;mod_directionality{VR_ephys} == -1;logical(ks_all_data{VR_ephys,42}')],1));
        WT_dCA1_pos_mod_interneuron = WT_dCA1_pos_mod_interneuron+sum(all([mod_index_p_value{VR_ephys} > 0.95;mod_directionality{VR_ephys} == 1;~logical(ks_all_data{VR_ephys,42}')],1));
        WT_dCA1_neg_mod_interneuron = WT_dCA1_neg_mod_interneuron+sum(all([mod_index_p_value{VR_ephys} > 0.95;mod_directionality{VR_ephys} == -1;~logical(ks_all_data{VR_ephys,42}')],1));
    elseif any(VR_ephys == AD_dCA1_sites) == 1
        AD_dCA1_pos_mod = AD_dCA1_pos_mod+sum(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == 1));
        AD_dCA1_neg_mod = AD_dCA1_neg_mod+sum(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == -1));
        AD_dCA1_pos_mod_pyramidal = AD_dCA1_pos_mod_pyramidal+sum(all([mod_index_p_value{VR_ephys} > 0.95;mod_directionality{VR_ephys} == 1;logical(ks_all_data{VR_ephys,42}')],1));
        AD_dCA1_neg_mod_pyramidal = AD_dCA1_neg_mod_pyramidal+sum(all([mod_index_p_value{VR_ephys} > 0.95;mod_directionality{VR_ephys} == -1;logical(ks_all_data{VR_ephys,42}')],1));
        AD_dCA1_pos_mod_interneuron = AD_dCA1_pos_mod_interneuron+sum(all([mod_index_p_value{VR_ephys} > 0.95;mod_directionality{VR_ephys} == 1;~logical(ks_all_data{VR_ephys,42}')],1));
        AD_dCA1_neg_mod_interneuron = AD_dCA1_neg_mod_interneuron+sum(all([mod_index_p_value{VR_ephys} > 0.95;mod_directionality{VR_ephys} == -1;~logical(ks_all_data{VR_ephys,42}')],1));
    elseif any(VR_ephys == WT_vCA1_sites) == 1
        WT_vCA1_pos_mod = WT_vCA1_pos_mod+sum(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == 1));
        WT_vCA1_neg_mod = WT_vCA1_neg_mod+sum(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == -1));
        WT_vCA1_pos_mod_pyramidal = WT_vCA1_pos_mod_pyramidal+sum(all([mod_index_p_value{VR_ephys} > 0.95;mod_directionality{VR_ephys} == 1;logical(ks_all_data{VR_ephys,42}')],1));
        WT_vCA1_neg_mod_pyramidal = WT_vCA1_neg_mod_pyramidal+sum(all([mod_index_p_value{VR_ephys} > 0.95;mod_directionality{VR_ephys} == -1;logical(ks_all_data{VR_ephys,42}')],1));
        WT_vCA1_pos_mod_interneuron = WT_vCA1_pos_mod_interneuron+sum(all([mod_index_p_value{VR_ephys} > 0.95;mod_directionality{VR_ephys} == 1;~logical(ks_all_data{VR_ephys,42}')],1));
        WT_vCA1_neg_mod_interneuron = WT_vCA1_neg_mod_interneuron+sum(all([mod_index_p_value{VR_ephys} > 0.95;mod_directionality{VR_ephys} == -1;~logical(ks_all_data{VR_ephys,42}')],1));
    elseif any(VR_ephys == AD_vCA1_sites) == 1
        AD_vCA1_pos_mod = AD_vCA1_pos_mod+sum(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == 1));
        AD_vCA1_neg_mod = AD_vCA1_neg_mod+sum(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == -1));
        AD_vCA1_pos_mod_pyramidal = AD_vCA1_pos_mod_pyramidal+sum(all([mod_index_p_value{VR_ephys} > 0.95;mod_directionality{VR_ephys} == 1;logical(ks_all_data{VR_ephys,42}')],1));
        AD_vCA1_neg_mod_pyramidal = AD_vCA1_neg_mod_pyramidal+sum(all([mod_index_p_value{VR_ephys} > 0.95;mod_directionality{VR_ephys} == -1;logical(ks_all_data{VR_ephys,42}')],1));
        AD_vCA1_pos_mod_interneuron = AD_vCA1_pos_mod_interneuron+sum(all([mod_index_p_value{VR_ephys} > 0.95;mod_directionality{VR_ephys} == 1;~logical(ks_all_data{VR_ephys,42}')],1));
        AD_vCA1_neg_mod_interneuron = AD_vCA1_neg_mod_interneuron+sum(all([mod_index_p_value{VR_ephys} > 0.95;mod_directionality{VR_ephys} == -1;~logical(ks_all_data{VR_ephys,42}')],1));
    end
end
WT_dCA1_non_mod = numel(cell2mat(mod_index_p_value(WT_dCA1_sites)'))-WT_dCA1_pos_mod-WT_dCA1_neg_mod; 
AD_dCA1_non_mod = numel(cell2mat(mod_index_p_value(AD_dCA1_sites)'))-AD_dCA1_pos_mod-AD_dCA1_neg_mod; 
WT_vCA1_non_mod = numel(cell2mat(mod_index_p_value(WT_vCA1_sites)'))-WT_vCA1_pos_mod-WT_vCA1_neg_mod; 
AD_vCA1_non_mod = numel(cell2mat(mod_index_p_value(AD_vCA1_sites)'))-AD_vCA1_pos_mod-AD_vCA1_neg_mod; 
WT_dCA1_non_mod_pyramidal = sum(cell2mat(ks_all_data(WT_dCA1_sites,42)))-WT_dCA1_pos_mod_pyramidal-WT_dCA1_neg_mod_pyramidal;
AD_dCA1_non_mod_pyramidal = sum(cell2mat(ks_all_data(AD_dCA1_sites,42)))-AD_dCA1_pos_mod_pyramidal-AD_dCA1_neg_mod_pyramidal;
WT_vCA1_non_mod_pyramidal = sum(cell2mat(ks_all_data(WT_vCA1_sites,42)))-WT_vCA1_pos_mod_pyramidal-WT_vCA1_neg_mod_pyramidal;
AD_vCA1_non_mod_pyramidal = sum(cell2mat(ks_all_data(AD_vCA1_sites,42)))-AD_vCA1_pos_mod_pyramidal-AD_vCA1_neg_mod_pyramidal;
WT_dCA1_non_mod_interneuron = sum(~cell2mat(ks_all_data(WT_dCA1_sites,42)))-WT_dCA1_pos_mod_interneuron-WT_dCA1_neg_mod_interneuron;
AD_dCA1_non_mod_interneuron = sum(~cell2mat(ks_all_data(AD_dCA1_sites,42)))-AD_dCA1_pos_mod_interneuron-AD_dCA1_neg_mod_interneuron;
WT_vCA1_non_mod_interneuron = sum(~cell2mat(ks_all_data(WT_vCA1_sites,42)))-WT_vCA1_pos_mod_interneuron-WT_vCA1_neg_mod_interneuron;
AD_vCA1_non_mod_interneuron = sum(~cell2mat(ks_all_data(AD_vCA1_sites,42)))-AD_vCA1_pos_mod_interneuron-AD_vCA1_neg_mod_interneuron;
rgbmatrix = [1 140/255 0;0 191/255 1;.6 .6 .6];
subplot(4,3,1); 
h = pie([WT_dCA1_pos_mod,WT_dCA1_neg_mod,WT_dCA1_non_mod],{num2str(WT_dCA1_pos_mod),num2str(WT_dCA1_neg_mod),num2str(WT_dCA1_non_mod)}); title('WT dCA1'); 
for k = 1:3; set(h(k*2-1),'FaceColor',rgbmatrix(k,:)); end
subplot(4,3,4); 
h = pie([AD_dCA1_pos_mod,AD_dCA1_neg_mod,AD_dCA1_non_mod],{num2str(AD_dCA1_pos_mod),num2str(AD_dCA1_neg_mod),num2str(AD_dCA1_non_mod)}); title('AD dCA1'); 
for k = 1:3; set(h(k*2-1),'FaceColor',rgbmatrix(k,:)); end
subplot(4,3,7); 
h = pie([WT_vCA1_pos_mod,WT_vCA1_neg_mod,WT_vCA1_non_mod],{num2str(WT_vCA1_pos_mod),num2str(WT_vCA1_neg_mod),num2str(WT_vCA1_non_mod)}); title('WT vCA1'); 
for k = 1:3; set(h(k*2-1),'FaceColor',rgbmatrix(k,:)); end
subplot(4,3,10); 
h = pie([AD_vCA1_pos_mod,AD_vCA1_neg_mod,AD_vCA1_non_mod],{num2str(AD_vCA1_pos_mod),num2str(AD_vCA1_neg_mod),num2str(AD_vCA1_non_mod)}); title('AD vCA1'); 
for k = 1:3; set(h(k*2-1),'FaceColor',rgbmatrix(k,:)); end
subplot(4,3,2); 
h = pie([WT_dCA1_pos_mod_pyramidal,WT_dCA1_neg_mod_pyramidal,WT_dCA1_non_mod_pyramidal],{num2str(WT_dCA1_pos_mod_pyramidal),num2str(WT_dCA1_neg_mod_pyramidal),num2str(WT_dCA1_non_mod_pyramidal)}); title('WT dCA1 pyramidal'); 
for k = 1:3; set(h(k*2-1),'FaceColor',rgbmatrix(k,:)); end
subplot(4,3,5); 
h = pie([AD_dCA1_pos_mod_pyramidal,AD_dCA1_neg_mod_pyramidal,AD_dCA1_non_mod_pyramidal],{num2str(AD_dCA1_pos_mod_pyramidal),num2str(AD_dCA1_neg_mod_pyramidal),num2str(AD_dCA1_non_mod_pyramidal)}); title('AD dCA1 pyramidal'); 
for k = 1:3; set(h(k*2-1),'FaceColor',rgbmatrix(k,:)); end
subplot(4,3,8); 
h = pie([WT_vCA1_pos_mod_pyramidal,WT_vCA1_neg_mod_pyramidal,WT_vCA1_non_mod_pyramidal],{num2str(WT_vCA1_pos_mod_pyramidal),num2str(WT_vCA1_neg_mod_pyramidal),num2str(WT_vCA1_non_mod_pyramidal)}); title('WT vCA1 pyramidal'); 
for k = 1:3; set(h(k*2-1),'FaceColor',rgbmatrix(k,:)); end
subplot(4,3,11); 
h = pie([AD_vCA1_pos_mod_pyramidal,AD_vCA1_neg_mod_pyramidal,AD_vCA1_non_mod_pyramidal],{num2str(AD_vCA1_pos_mod_pyramidal),num2str(AD_vCA1_neg_mod_pyramidal),num2str(AD_vCA1_non_mod_pyramidal)}); title('AD vCA1 pyramidal'); 
for k = 1:3; set(h(k*2-1),'FaceColor',rgbmatrix(k,:)); end
subplot(4,3,3); 
h = pie([WT_dCA1_pos_mod_interneuron,WT_dCA1_neg_mod_interneuron,WT_dCA1_non_mod_interneuron],{num2str(WT_dCA1_pos_mod_interneuron),num2str(WT_dCA1_neg_mod_interneuron),num2str(WT_dCA1_non_mod_interneuron)}); title('WT dCA1 interneuron'); 
for k = 1:3; set(h(k*2-1),'FaceColor',rgbmatrix(k,:)); end
subplot(4,3,6); 
h = pie([AD_dCA1_pos_mod_interneuron,AD_dCA1_neg_mod_interneuron,AD_dCA1_non_mod_interneuron],{num2str(AD_dCA1_pos_mod_interneuron),num2str(AD_dCA1_neg_mod_interneuron),num2str(AD_dCA1_non_mod_interneuron)}); title('AD dCA1 interneuron'); 
for k = 1:3; set(h(k*2-1),'FaceColor',rgbmatrix(k,:)); end
subplot(4,3,9); 
h = pie([WT_vCA1_pos_mod_interneuron,WT_vCA1_neg_mod_interneuron,WT_vCA1_non_mod_interneuron],{num2str(WT_vCA1_pos_mod_interneuron),num2str(WT_vCA1_neg_mod_interneuron),num2str(WT_vCA1_non_mod_interneuron)}); title('WT vCA1 interneuron'); 
for k = 1:3; set(h(k*2-1),'FaceColor',rgbmatrix(k,:)); end
subplot(4,3,12); 
h = pie([AD_vCA1_pos_mod_interneuron,AD_vCA1_neg_mod_interneuron,AD_vCA1_non_mod_interneuron],{num2str(AD_vCA1_pos_mod_interneuron),num2str(AD_vCA1_neg_mod_interneuron),num2str(AD_vCA1_non_mod_interneuron)}); title('AD vCA1 interneuron'); 
for k = 1:3; set(h(k*2-1),'FaceColor',rgbmatrix(k,:)); end
legend({'Positively modulated units','Negatively modulated units','Non-modulated units'},'location','southeast'); 
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 10'); 
saveas(g,'SWR_pie','fig'); orient(g,'landscape'); print('SWR_pie.pdf','-bestfit','-painters','-dpdf'); close(g); 
%% Plot and compare the proportions of positively, negatively, and unmodulated units
fraction_mod = zeros(size(ks_all_data,1),1); 
fraction_pos_mod = fraction_mod;
fraction_neg_mod = fraction_mod; 
fraction_non_mod = fraction_mod; 
for VR_ephys = 1:size(ks_all_data,1)
    fraction_mod(VR_ephys) = mean(mod_index_p_value{VR_ephys} > 0.95); 
    fraction_pos_mod(VR_ephys) = mean(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == 1));
    fraction_neg_mod(VR_ephys) = mean(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == -1)); 
    fraction_non_mod(VR_ephys) = mean(mod_index_p_value{VR_ephys} <= 0.95); 
end
f = figure; set(f,'Position',[0 0 1500 1500]); 
subplot(3,4,1); 
scatter(ones(numel(WT_dCA1_sites),1),fraction_mod(WT_dCA1_sites),72,'k','filled'); hold on;
scatter(3*ones(numel(WT_vCA1_sites),1),fraction_mod(WT_vCA1_sites),72,[.8 .8 .8],'filled');
[p,~] = ranksum(fraction_mod(WT_dCA1_sites),fraction_mod(WT_vCA1_sites)); text(2,1,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 1.2]); xticks([1 3]); xticklabels({'WT dCA1','WT vCA1'}); ylabel('Fraction +/- modulated by SWRs'); axis square
subplot(3,4,2); 
scatter(ones(numel(WT_dCA1_sites),1),fraction_pos_mod(WT_dCA1_sites),72,'k','filled'); hold on;
scatter(3*ones(numel(WT_vCA1_sites),1),fraction_pos_mod(WT_vCA1_sites),72,[.8 .8 .8],'filled');
[p,~] = ranksum(fraction_pos_mod(WT_dCA1_sites),fraction_pos_mod(WT_vCA1_sites)); text(2,1,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 1.2]); xticks([1 3]); xticklabels({'WT dCA1','WT vCA1'}); ylabel('Fraction + modulated by SWRs'); axis square
subplot(3,4,3); 
scatter(ones(numel(WT_dCA1_sites),1),fraction_neg_mod(WT_dCA1_sites),72,'k','filled'); hold on;
scatter(3*ones(numel(WT_vCA1_sites),1),fraction_neg_mod(WT_vCA1_sites),72,[.8 .8 .8],'filled');
[p,~] = ranksum(fraction_neg_mod(WT_dCA1_sites),fraction_neg_mod(WT_vCA1_sites)); text(2,1,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 1.2]); xticks([1 3]); xticklabels({'WT dCA1','WT vCA1'}); ylabel('Fraction - modulated by SWRs'); axis square
subplot(3,4,4); 
scatter(ones(numel(WT_dCA1_sites),1),fraction_non_mod(WT_dCA1_sites),72,'k','filled'); hold on;
scatter(3*ones(numel(WT_vCA1_sites),1),fraction_non_mod(WT_vCA1_sites),72,[.8 .8 .8],'filled');
[p,~] = ranksum(fraction_non_mod(WT_dCA1_sites),fraction_non_mod(WT_vCA1_sites)); text(2,1,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 1.2]); xticks([1 3]); xticklabels({'WT dCA1','WT vCA1'}); ylabel('Fraction non-modulated by SWRs'); axis square
subplot(3,4,5)
scatter(ones(numel(WT_dCA1_sites),1),fraction_mod(WT_dCA1_sites),72,'k','filled'); hold on;
scatter(3*ones(numel(AD_dCA1_sites),1),fraction_mod(AD_dCA1_sites),72,'r','filled');
[p,~] = ranksum(fraction_mod(WT_dCA1_sites),fraction_mod(AD_dCA1_sites)); text(2,1,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 1.2]); xticks([1 3]); xticklabels({'WT','AD'}); ylabel('Fraction +/- modulated by SWRs'); axis square
subplot(3,4,6)
scatter(ones(numel(WT_dCA1_sites),1),fraction_pos_mod(WT_dCA1_sites),72,'k','filled'); hold on;
scatter(3*ones(numel(AD_dCA1_sites),1),fraction_pos_mod(AD_dCA1_sites),72,'r','filled');
[p,~] = ranksum(fraction_pos_mod(WT_dCA1_sites),fraction_pos_mod(AD_dCA1_sites)); text(2,1,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 1.2]); xticks([1 3]); xticklabels({'WT','AD'}); ylabel('Fraction + modulated by SWRs'); axis square
subplot(3,4,7)
scatter(ones(numel(WT_dCA1_sites),1),fraction_neg_mod(WT_dCA1_sites),72,'k','filled'); hold on;
scatter(3*ones(numel(AD_dCA1_sites),1),fraction_neg_mod(AD_dCA1_sites),72,'r','filled');
[p,~] = ranksum(fraction_neg_mod(WT_dCA1_sites),fraction_neg_mod(AD_dCA1_sites)); text(2,1,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 1.2]); xticks([1 3]); xticklabels({'WT','AD'}); ylabel('Fraction - modulated by SWRs'); axis square
subplot(3,4,8)
scatter(ones(numel(WT_dCA1_sites),1),fraction_non_mod(WT_dCA1_sites),72,'k','filled'); hold on;
scatter(3*ones(numel(AD_dCA1_sites),1),fraction_non_mod(AD_dCA1_sites),72,'r','filled');
[p,~] = ranksum(fraction_non_mod(WT_dCA1_sites),fraction_non_mod(AD_dCA1_sites)); text(2,1,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 1.2]); xticks([1 3]); xticklabels({'WT','AD'}); ylabel('Fraction non-modulated by SWRs'); axis square
subplot(3,4,9)
scatter(ones(numel(WT_vCA1_sites),1),fraction_mod(WT_vCA1_sites),72,[.8 .8 .8],'filled'); hold on;
scatter(3*ones(numel(AD_vCA1_sites),1),fraction_mod(AD_vCA1_sites),72,[1 .8 .8],'filled');
[p,~] = ranksum(fraction_mod(WT_vCA1_sites),fraction_mod(AD_vCA1_sites)); text(2,1,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 1.2]); xticks([1 3]); xticklabels({'WT','AD'}); ylabel('Fraction +/- modulated by SWRs'); axis square
subplot(3,4,10)
scatter(ones(numel(WT_vCA1_sites),1),fraction_pos_mod(WT_vCA1_sites),72,[.8 .8 .8],'filled'); hold on;
scatter(3*ones(numel(AD_vCA1_sites),1),fraction_pos_mod(AD_vCA1_sites),72,[1 .8 .8],'filled');
[p,~] = ranksum(fraction_pos_mod(WT_vCA1_sites),fraction_pos_mod(AD_vCA1_sites)); text(2,1,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 1.2]); xticks([1 3]); xticklabels({'WT','AD'}); ylabel('Fraction + modulated by SWRs'); axis square
subplot(3,4,11)
scatter(ones(numel(WT_vCA1_sites),1),fraction_neg_mod(WT_vCA1_sites),72,[.8 .8 .8],'filled'); hold on;
scatter(3*ones(numel(AD_vCA1_sites),1),fraction_neg_mod(AD_vCA1_sites),72,[1 .8 .8],'filled');
[p,~] = ranksum(fraction_neg_mod(WT_vCA1_sites),fraction_neg_mod(AD_vCA1_sites)); text(2,1,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 1.2]); xticks([1 3]); xticklabels({'WT','AD'}); ylabel('Fraction - modulated by SWRs'); axis square
subplot(3,4,12)
scatter(ones(numel(WT_vCA1_sites),1),fraction_non_mod(WT_vCA1_sites),72,[.8 .8 .8],'filled'); hold on;
scatter(3*ones(numel(AD_vCA1_sites),1),fraction_non_mod(AD_vCA1_sites),72,[1 .8 .8],'filled');
[p,~] = ranksum(fraction_non_mod(WT_vCA1_sites),fraction_non_mod(AD_vCA1_sites)); text(2,1,strcat('p = ',num2str(p))); 
xlim([0 4]); ylim([0 1.2]); xticks([1 3]); xticklabels({'WT','AD'}); ylabel('Fraction non-modulated by SWRs'); axis square
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 10'); 
saveas(f,'fraction_SWR_mod','fig'); orient(f,'landscape'); print('fraction_SWR_mod.pdf','-bestfit','-painters','-dpdf'); close(f); 
%% Plot the modulation index
f = figure; set(gcf,'Position',[0 0 1500 1500]); 
subplot(3,3,1); 
histogram(log10(WT_dCA1_mod_index),-4:.1:11,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(log10(AD_dCA1_mod_index),-4:.1:11,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
[p,~] = ranksum(WT_dCA1_mod_index,AD_dCA1_mod_index); text(2,1,strcat('p = ',num2str(p))); 
xlabel('log_{10}Modulation index'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([-4 11]); axis square; title('All units'); 
subplot(3,3,2); 
histogram(log10(WT_vCA1_mod_index),-4:.1:11,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(log10(AD_vCA1_mod_index),-4:.1:11,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
[p,~] = ranksum(WT_vCA1_mod_index,AD_vCA1_mod_index); text(2,1,strcat('p = ',num2str(p))); 
xlabel('log_{10}Modulation index'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([-4 11]); axis square; title('All units'); 
subplot(3,3,4); 
histogram(log10(WT_dCA1_PETH_pos_mod_index),-4:.1:11,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(log10(AD_dCA1_PETH_pos_mod_index),-4:.1:11,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
[p,~] = ranksum(WT_dCA1_PETH_pos_mod_index,AD_dCA1_PETH_pos_mod_index); text(2,1,strcat('p = ',num2str(p))); 
xlabel('log_{10}Modulation index'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([-4 11]); axis square; title('+Mod Units Only'); 
subplot(3,3,5); 
histogram(log10(WT_dCA1_PETH_neg_mod_index),-4:.1:11,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(log10(AD_dCA1_PETH_neg_mod_index),-4:.1:11,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
[p,~] = ranksum(WT_dCA1_PETH_neg_mod_index,AD_dCA1_PETH_neg_mod_index); text(2,1,strcat('p = ',num2str(p))); 
xlabel('log_{10}Modulation index'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([-4 11]); axis square; title('-Mod Units Only'); 
subplot(3,3,6); 
histogram(log10(WT_dCA1_PETH_non_mod_index),-4:.1:11,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(log10(AD_dCA1_PETH_non_mod_index),-4:.1:11,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
[p,~] = ranksum(WT_dCA1_PETH_non_mod_index,AD_dCA1_PETH_non_mod_index); text(2,1,strcat('p = ',num2str(p))); 
xlabel('log_{10}Modulation index'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([-4 11]); axis square; title('Non-Mod Units Only'); 
subplot(3,3,7); 
histogram(log10(WT_vCA1_PETH_pos_mod_index),-4:.1:11,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(log10(AD_vCA1_PETH_pos_mod_index),-4:.1:11,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
[p,~] = ranksum(WT_vCA1_PETH_pos_mod_index,AD_vCA1_PETH_pos_mod_index); text(2,1,strcat('p = ',num2str(p))); 
xlabel('log_{10}Modulation index'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([-4 11]); axis square; title('+Mod Units Only'); 
subplot(3,3,8); 
histogram(log10(WT_vCA1_PETH_neg_mod_index),-4:.1:11,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(log10(AD_vCA1_PETH_neg_mod_index),-4:.1:11,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
[p,~] = ranksum(WT_vCA1_PETH_neg_mod_index,AD_vCA1_PETH_neg_mod_index); text(2,1,strcat('p = ',num2str(p))); 
xlabel('log_{10}Modulation index'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([-4 11]); axis square; title('-Mod Units Only'); 
subplot(3,3,9); 
histogram(log10(WT_vCA1_PETH_non_mod_index),-4:.1:11,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(log10(AD_vCA1_PETH_non_mod_index),-4:.1:11,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
[p,~] = ranksum(WT_vCA1_PETH_non_mod_index,AD_vCA1_PETH_non_mod_index); text(2,1,strcat('p = ',num2str(p))); 
xlabel('log_{10}Modulation index'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([-4 11]); axis square; title('Non-Mod Units Only'); 
% cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 10'); 
% saveas(f,'modulation_index','fig'); orient(f,'landscape'); print('modulation_index.pdf','-bestfit','-painters','-dpdf'); close(f);
%% Plot the lag and FWHM
WT_dCA1_lag_nonan = WT_dCA1_lag(mod_index_p_value{VR_ephys} > 0.95); WT_dCA1_lag_nonan(isnan(WT_dCA1_lag_nonan)) = []; WT_dCA1_FWHM_nonan = WT_dCA1_FWHM(mod_index_p_value{VR_ephys} > 0.95); WT_dCA1_FWHM_nonan(isnan(WT_dCA1_FWHM_nonan)) = []; 
AD_dCA1_lag_nonan = AD_dCA1_lag(mod_index_p_value{VR_ephys} > 0.95); AD_dCA1_lag_nonan(isnan(AD_dCA1_lag_nonan)) = []; AD_dCA1_FWHM_nonan = AD_dCA1_FWHM(mod_index_p_value{VR_ephys} > 0.95); AD_dCA1_FWHM_nonan(isnan(AD_dCA1_FWHM_nonan)) = []; 
WT_vCA1_lag_nonan = WT_vCA1_lag(mod_index_p_value{VR_ephys} > 0.95); WT_vCA1_lag_nonan(isnan(WT_vCA1_lag_nonan)) = []; WT_vCA1_FWHM_nonan = WT_vCA1_FWHM(mod_index_p_value{VR_ephys} > 0.95); WT_vCA1_FWHM_nonan(isnan(WT_vCA1_FWHM_nonan)) = []; 
AD_vCA1_lag_nonan = AD_vCA1_lag(mod_index_p_value{VR_ephys} > 0.95); AD_vCA1_lag_nonan(isnan(AD_vCA1_lag_nonan)) = []; AD_vCA1_FWHM_nonan = AD_vCA1_FWHM(mod_index_p_value{VR_ephys} > 0.95); AD_vCA1_FWHM_nonan(isnan(AD_vCA1_FWHM_nonan)) = []; 
WT_dCA1_lag_pos_mod_nonan = WT_dCA1_lag(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == 1)); WT_dCA1_lag_pos_mod_nonan(isnan(WT_dCA1_lag_pos_mod_nonan)) = [];
AD_dCA1_lag_pos_mod_nonan = AD_dCA1_lag(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == 1)); AD_dCA1_lag_pos_mod_nonan(isnan(AD_dCA1_lag_pos_mod_nonan)) = [];
WT_vCA1_lag_pos_mod_nonan = WT_vCA1_lag(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == 1)); WT_vCA1_lag_pos_mod_nonan(isnan(WT_vCA1_lag_pos_mod_nonan)) = [];
AD_vCA1_lag_pos_mod_nonan = AD_vCA1_lag(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == 1)); AD_vCA1_lag_pos_mod_nonan(isnan(AD_vCA1_lag_pos_mod_nonan)) = [];
WT_dCA1_lag_neg_mod_nonan = WT_dCA1_lag(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == -1)); WT_dCA1_lag_neg_mod_nonan(isnan(WT_dCA1_lag_neg_mod_nonan)) = [];
AD_dCA1_lag_neg_mod_nonan = AD_dCA1_lag(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == -1)); AD_dCA1_lag_neg_mod_nonan(isnan(AD_dCA1_lag_neg_mod_nonan)) = [];
WT_vCA1_lag_neg_mod_nonan = WT_vCA1_lag(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == -1)); WT_vCA1_lag_neg_mod_nonan(isnan(WT_vCA1_lag_neg_mod_nonan)) = [];
AD_vCA1_lag_neg_mod_nonan = AD_vCA1_lag(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == -1)); AD_vCA1_lag_neg_mod_nonan(isnan(AD_vCA1_lag_neg_mod_nonan)) = [];
WT_dCA1_FWHM_pos_mod_nonan = WT_dCA1_FWHM(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == 1)); WT_dCA1_FWHM_pos_mod_nonan(isnan(WT_dCA1_FWHM_pos_mod_nonan)) = [];
AD_dCA1_FWHM_pos_mod_nonan = AD_dCA1_FWHM(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == 1)); AD_dCA1_FWHM_pos_mod_nonan(isnan(AD_dCA1_FWHM_pos_mod_nonan)) = [];
WT_vCA1_FWHM_pos_mod_nonan = WT_vCA1_FWHM(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == 1)); WT_vCA1_FWHM_pos_mod_nonan(isnan(WT_vCA1_FWHM_pos_mod_nonan)) = [];
AD_vCA1_FWHM_pos_mod_nonan = AD_vCA1_FWHM(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == 1)); AD_vCA1_FWHM_pos_mod_nonan(isnan(AD_vCA1_FWHM_pos_mod_nonan)) = [];
WT_dCA1_FWHM_neg_mod_nonan = WT_dCA1_FWHM(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == -1)); WT_dCA1_FWHM_neg_mod_nonan(isnan(WT_dCA1_FWHM_neg_mod_nonan)) = [];
AD_dCA1_FWHM_neg_mod_nonan = AD_dCA1_FWHM(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == -1)); AD_dCA1_FWHM_neg_mod_nonan(isnan(AD_dCA1_FWHM_neg_mod_nonan)) = [];
WT_vCA1_FWHM_neg_mod_nonan = WT_vCA1_FWHM(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == -1)); WT_vCA1_FWHM_neg_mod_nonan(isnan(WT_vCA1_FWHM_neg_mod_nonan)) = [];
AD_vCA1_FWHM_neg_mod_nonan = AD_vCA1_FWHM(and(mod_index_p_value{VR_ephys} > 0.95,mod_directionality{VR_ephys} == -1)); AD_vCA1_FWHM_neg_mod_nonan(isnan(AD_vCA1_FWHM_neg_mod_nonan)) = [];
f = figure; set(gcf,'Position',[0 0 1000 1500]); 
subplot(3,2,1); 
histogram(WT_dCA1_lag_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(AD_dCA1_lag_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
[p,~] = ranksum(WT_dCA1_lag_nonan,AD_dCA1_lag_nonan); text(2,1,strcat('p = ',num2str(p))); text(300,.1,strcat('WT = ',num2str(numel(WT_dCA1_lag_nonan)),';AD = ',num2str(numel(AD_dCA1_lag_nonan)))); 
xlabel('SWR modulation lag (ms)'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([0 480]); axis square; title('All units'); 
subplot(3,2,2); 
histogram(WT_vCA1_lag_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(AD_vCA1_lag_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
[p,~] = ranksum(WT_vCA1_lag_nonan,AD_vCA1_lag_nonan); text(2,1,strcat('p = ',num2str(p))); text(300,.1,strcat('WT = ',num2str(numel(WT_vCA1_lag_nonan)),';AD = ',num2str(numel(AD_vCA1_lag_nonan)))); 
xlabel('SWR modulation lag (ms)'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([0 480]); axis square; title('All units'); 
subplot(3,2,3); 
histogram(WT_dCA1_lag_pos_mod_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(AD_dCA1_lag_pos_mod_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
[p,~] = ranksum(WT_dCA1_lag_pos_mod_nonan,AD_dCA1_lag_pos_mod_nonan); text(2,1,strcat('p = ',num2str(p))); text(300,.1,strcat('WT = ',num2str(numel(WT_dCA1_lag_pos_mod_nonan)),';AD = ',num2str(numel(AD_dCA1_lag_pos_mod_nonan))));  
xlabel('SWR modulation lag (ms)'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([0 480]); axis square; title('Positively modulated units'); 
subplot(3,2,4); 
histogram(WT_vCA1_lag_pos_mod_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(AD_vCA1_lag_pos_mod_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
[p,~] = ranksum(WT_vCA1_lag_pos_mod_nonan,AD_vCA1_lag_pos_mod_nonan); text(2,1,strcat('p = ',num2str(p))); text(300,.1,strcat('WT = ',num2str(numel(WT_vCA1_lag_pos_mod_nonan)),';AD = ',num2str(numel(AD_vCA1_lag_pos_mod_nonan))));  
xlabel('SWR modulation lag (ms)'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([0 480]); axis square; title('Positively modulated units'); 
subplot(3,2,5); 
histogram(WT_dCA1_lag_neg_mod_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(AD_dCA1_lag_neg_mod_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
[p,~] = ranksum(WT_dCA1_lag_neg_mod_nonan,AD_dCA1_lag_neg_mod_nonan); text(2,1,strcat('p = ',num2str(p))); text(300,.1,strcat('WT = ',num2str(numel(WT_dCA1_lag_neg_mod_nonan)),';AD = ',num2str(numel(AD_dCA1_lag_neg_mod_nonan))));  
xlabel('SWR modulation lag (ms)'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([0 480]); axis square; title('Negatively modulated units'); 
subplot(3,2,6); 
histogram(WT_vCA1_lag_neg_mod_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(AD_vCA1_lag_neg_mod_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
[p,~] = ranksum(WT_vCA1_lag_neg_mod_nonan,AD_vCA1_lag_neg_mod_nonan); text(2,1,strcat('p = ',num2str(p))); text(300,.1,strcat('WT = ',num2str(numel(WT_vCA1_lag_neg_mod_nonan)),';AD = ',num2str(numel(AD_vCA1_lag_neg_mod_nonan))));  
xlabel('SWR modulation lag (ms)'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([0 480]); axis square; title('Negatively modulated units'); 
% cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 10'); 
% saveas(f,'SWR_mod_lag','fig'); orient(f,'landscape'); print('SWR_mod_lag.pdf','-bestfit','-painters','-dpdf'); close(f);
g = figure; set(gcf,'Position',[0 0 1000 1500]); 
subplot(3,2,1); 
histogram(WT_dCA1_FWHM_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(AD_dCA1_FWHM_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
[p,~] = ranksum(WT_dCA1_FWHM_nonan,AD_dCA1_FWHM_nonan); text(2,1,strcat('p = ',num2str(p))); text(300,.1,strcat('WT = ',num2str(numel(WT_dCA1_FWHM_nonan)),';AD = ',num2str(numel(AD_dCA1_FWHM_nonan)))); 
xlabel('SWR modulation FWHM (ms)'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([0 480]); axis square; title('All units'); 
subplot(3,2,2); 
histogram(WT_vCA1_FWHM_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(AD_vCA1_FWHM_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
[p,~] = ranksum(WT_vCA1_FWHM_nonan,AD_vCA1_FWHM_nonan); text(2,1,strcat('p = ',num2str(p))); text(300,.1,strcat('WT = ',num2str(numel(WT_vCA1_FWHM_nonan)),';AD = ',num2str(numel(AD_vCA1_FWHM_nonan))));  
xlabel('SWR modulation FWHM (ms)'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([0 480]); axis square; title('All units'); 
subplot(3,2,3); 
histogram(WT_dCA1_FWHM_pos_mod_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(AD_dCA1_FWHM_pos_mod_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
[p,~] = ranksum(WT_dCA1_FWHM_pos_mod_nonan,AD_dCA1_FWHM_pos_mod_nonan); text(2,1,strcat('p = ',num2str(p))); text(300,.1,strcat('WT = ',num2str(numel(WT_dCA1_FWHM_pos_mod_nonan)),';AD = ',num2str(numel(AD_dCA1_FWHM_pos_mod_nonan))));   
xlabel('SWR modulation FWHM (ms)'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([0 480]); axis square; title('Positively modulated units'); 
subplot(3,2,4); 
histogram(WT_vCA1_FWHM_pos_mod_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(AD_vCA1_FWHM_pos_mod_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
[p,~] = ranksum(WT_vCA1_FWHM_pos_mod_nonan,AD_vCA1_FWHM_pos_mod_nonan); text(2,1,strcat('p = ',num2str(p))); text(300,.1,strcat('WT = ',num2str(numel(WT_vCA1_FWHM_pos_mod_nonan)),';AD = ',num2str(numel(AD_vCA1_FWHM_pos_mod_nonan))));   
xlabel('SWR modulation FWHM (ms)'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([0 480]); axis square; title('Positively modulated units'); 
subplot(3,2,5); 
histogram(WT_dCA1_FWHM_neg_mod_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(AD_dCA1_FWHM_neg_mod_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
[p,~] = ranksum(WT_dCA1_FWHM_neg_mod_nonan,AD_dCA1_FWHM_neg_mod_nonan); text(2,1,strcat('p = ',num2str(p))); text(300,.1,strcat('WT = ',num2str(numel(WT_dCA1_FWHM_neg_mod_nonan)),';AD = ',num2str(numel(AD_dCA1_FWHM_neg_mod_nonan))));  
xlabel('SWR modulation FWHM (ms)'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([0 480]); axis square; title('Negatively modulated units'); 
subplot(3,2,6); 
histogram(WT_vCA1_FWHM_neg_mod_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(AD_vCA1_FWHM_neg_mod_nonan,0:1:480,'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
[p,~] = ranksum(WT_vCA1_FWHM_neg_mod_nonan,AD_vCA1_FWHM_neg_mod_nonan); text(2,1,strcat('p = ',num2str(p))); text(300,.1,strcat('WT = ',num2str(numel(WT_vCA1_FWHM_neg_mod_nonan)),';AD = ',num2str(numel(AD_vCA1_FWHM_neg_mod_nonan))));  
xlabel('SWR modulation FWHM (ms)'); ylabel('Cumulative Probability'); ylim([0 1]); xlim([0 480]); axis square; title('Negatively modulated units'); 
% cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 10'); 
% saveas(g,'SWR_mod_FWHM','fig'); orient(g,'landscape'); print('SWR_mod_FWHM.pdf','-bestfit','-painters','-dpdf'); close(g);
%% Plot FR vs. mod index
f = figure; set(gcf,'Position',[0 0 2000 500]); 
subplot(1,4,1); scatter(cell2mat(FR(WT_dCA1_sites)),cell2mat(mod_index(WT_dCA1_sites)'),12,'k','filled'); set(gca,'XScale','log','YScale','log'); axis square; xlim([1e-1 1e2]); xlabel('Mean FR (Hz)'); ylim([1e-4, 1e12]); ylabel('Modulation index'); 
subplot(1,4,2); scatter(cell2mat(FR(AD_dCA1_sites)),cell2mat(mod_index(AD_dCA1_sites)'),12,'r','filled'); set(gca,'XScale','log','YScale','log'); axis square; xlim([1e-1 1e2]); xlabel('Mean FR (Hz)'); ylim([1e-4, 1e12]); ylabel('Modulation index'); 
subplot(1,4,3); scatter(cell2mat(FR(WT_dCA1_sites)),cell2mat(mod_index(WT_dCA1_sites)'),12,[.8 .8 .8],'filled'); set(gca,'XScale','log','YScale','log'); axis square; xlim([1e-1 1e2]); xlabel('Mean FR (Hz)'); ylim([1e-4, 1e12]); ylabel('Modulation index');  
subplot(1,4,4); scatter(cell2mat(FR(AD_dCA1_sites)),cell2mat(mod_index(AD_dCA1_sites)'),12,[1 .8 .8],'filled'); set(gca,'XScale','log','YScale','log'); axis square; xlim([1e-1 1e2]); xlabel('Mean FR (Hz)'); ylim([1e-4, 1e12]); ylabel('Modulation index'); 
%% Find random epochs of same length that occur outside of SWRs
nonSWR_epochs = cell(size(ks_all_data,1),1); %contains start and end times of random epochs that do not contain SWRs anywhere in the interval
for VR_ephys = 1:size(ks_all_data,1)
    nonSWR_epochs{VR_ephys} = zeros(size(SWR_epochs{VR_ephys},1),2);    
    all_nonSWR_epochs = zeros((size(SWR_epochs,1)-1),2);
    for i = 1:(size(SWR_epochs{VR_ephys},1)-1)
        all_nonSWR_epochs(i,:) = [(SWR_epochs{VR_ephys}(i,2)+1),(SWR_epochs{VR_ephys}((i+1),1)-1)]; %start and end times of all epochs between SWR epochs
    end
    duration_nonSWR = diff(all_nonSWR_epochs,1,2)+1; %duration (in 30kHz samples) of each nonSWR epoch
    for i = 1:size(SWR_unamb_epochs{VR_ephys},1) %for each unamb_SWR epoch        
        nonSWR_eoi = all_nonSWR_epochs(duration_nonSWR >= (diff(SWR_unamb_epochs{VR_ephys}(i,:))+1),:); %only the nonSWR epochs that are at least as big as the current SWR epoch    
        randidx = randi(size(nonSWR_eoi,1)); %pick a random interval in nonSWR_eoi
        tmp_nonSWR_start = randi([nonSWR_eoi(randidx,1),(nonSWR_eoi(randidx,2)-diff(SWR_unamb_epochs{VR_ephys}(i,:)))]); 
        tmp_nonSWR_end = tmp_nonSWR_start+diff(SWR_unamb_epochs{VR_ephys}(i,:)); %find an end point that makes the same duration as SWR epoch            
        nonSWR_epochs{VR_ephys}(i,:) = [tmp_nonSWR_start,tmp_nonSWR_end];
    end    
    VR_ephys
end
%% Construct spike trains that occur pre- and post-SWR
preSWR_spikes = cell(size(ks_all_data,1),1); 
postSWR_spikes = preSWR_spikes; 
periSWR_spikes = preSWR_spikes;
prenonSWR_spikes = preSWR_spikes;
postnonSWR_spikes = preSWR_spikes;
perinonSWR_spikes = preSWR_spikes; 
peri_window = 500*30; %number of timepoints around start of SWR to count spikes
for VR_ephys = 1:size(ks_all_data,1)
    preSWR_spikes{VR_ephys,1} = cell(size(ks_all_data{VR_ephys,1},1),1); %will contain all spikes for current unit that occur in the window before the start of an SWR
    postSWR_spikes{VR_ephys,1} = preSWR_spikes{VR_ephys,1}; %will contain all spikes for the current unit that occur in the window after the start of an SWR
    periSWR_spikes{VR_ephys,1} = preSWR_spikes{VR_ephys,1}; %will contain all spikes for the current unit that occur in the window around the start of an SWR
    prenonSWR_spikes{VR_ephys,1} = preSWR_spikes{VR_ephys,1};
    postnonSWR_spikes{VR_ephys,1} = preSWR_spikes{VR_ephys,1}; 
    perinonSWR_spikes{VR_ephys,1} = preSWR_spikes{VR_ephys,1}; 
    for u = 1:size(ks_all_data{VR_ephys,1},1) %for each unit        
        st_delimiter = 0; %will specify the total duration of SWR epochs (for pre and post)
        st_delimiter2 = 0; %will specify the total duration of SWR epochs (for peri)
        for i = 1:size(SWR_unamb_epochs{VR_ephys},1) %for each unamb_SWR epoch
            %construct preSWR spike train
            spikes_in_epoch = ks_all_data{VR_ephys,1}{u,2}(and(ks_all_data{VR_ephys,1}{u,2} < SWR_unamb_epochs{VR_ephys}(i,1),...
                ks_all_data{VR_ephys,1}{u,2} >= SWR_unamb_epochs{VR_ephys}(i,1)-peri_window)); %find all the spikes from this unit that occurred in a pre-SWR epoch
            spikes_in_epoch = spikes_in_epoch-(SWR_unamb_epochs{VR_ephys}(i,1)-peri_window-1)+st_delimiter; %shift the spikes so you can concatenate them into one spike train
            preSWR_spikes{VR_ephys}{u} = [preSWR_spikes{VR_ephys}{u},spikes_in_epoch]; %concatenate the spikes from this SWR_unamb epoch to that of the previous SWR_unamb epochs            
            %construct postSWR spike train
            spikes_in_epoch = ks_all_data{VR_ephys,1}{u,2}(and(ks_all_data{VR_ephys,1}{u,2} >= SWR_unamb_epochs{VR_ephys}(i,1),...
                ks_all_data{VR_ephys,1}{u,2} < SWR_unamb_epochs{VR_ephys}(i,1)+peri_window)); %find all the spikes from this unit that occurred from a post-SWR epoch
            spikes_in_epoch = spikes_in_epoch-(SWR_unamb_epochs{VR_ephys}(i,1)-1)+st_delimiter; %shift the spikes so you can concatenate them into one spike train
            postSWR_spikes{VR_ephys}{u} = [postSWR_spikes{VR_ephys}{u},spikes_in_epoch]; %concatenate the spikes from this nonSWR epoch to that of the previous nonSWR epochs              
            %construct prenonSWR spike train
            spikes_in_epoch = ks_all_data{VR_ephys,1}{u,2}(and(ks_all_data{VR_ephys,1}{u,2} < nonSWR_epochs{VR_ephys}(i,1),...
                ks_all_data{VR_ephys,1}{u,2} >= nonSWR_epochs{VR_ephys}(i,1)-peri_window)); %find all the spikes from this unit that occurred in a pre-SWR epoch
            spikes_in_epoch = spikes_in_epoch-(nonSWR_epochs{VR_ephys}(i,1)-peri_window-1)+st_delimiter; %shift the spikes so you can concatenate them into one spike train
            prenonSWR_spikes{VR_ephys}{u} = [prenonSWR_spikes{VR_ephys}{u},spikes_in_epoch]; %concatenate the spikes from this nonSWR epoch to that of the previous nonSWR epochs            
            %construct postnonSWR spike train
            spikes_in_epoch = ks_all_data{VR_ephys,1}{u,2}(and(ks_all_data{VR_ephys,1}{u,2} >= nonSWR_epochs{VR_ephys}(i,1),...
                ks_all_data{VR_ephys,1}{u,2} < nonSWR_epochs{VR_ephys}(i,1)+peri_window)); %find all the spikes from this unit that occurred from a post-SWR epoch
            spikes_in_epoch = spikes_in_epoch-(nonSWR_epochs{VR_ephys}(i,1)-1)+st_delimiter; %shift the spikes so you can concatenate them into one spike train
            postnonSWR_spikes{VR_ephys}{u} = [postnonSWR_spikes{VR_ephys}{u},spikes_in_epoch]; %concatenate the spikes from this nonSWR epoch to that of the previous nonSWR epochs  
            st_delimiter = st_delimiter+peri_window; %calculate the total duration of the spike train thus far (for pre and post)
            %construct periSWR spike train
            spikes_in_epoch = ks_all_data{VR_ephys,1}{u,2}(and(ks_all_data{VR_ephys,1}{u,2} < SWR_unamb_epochs{VR_ephys}(i,1)+peri_window,...
                ks_all_data{VR_ephys,1}{u,2} >= SWR_unamb_epochs{VR_ephys}(i,1)-peri_window)); %find all the spikes from this unit that occurred in a pre-SWR epoch
            spikes_in_epoch = spikes_in_epoch-(SWR_unamb_epochs{VR_ephys}(i,1)-peri_window-1)+st_delimiter2; %shift the spikes so you can concatenate them into one spike train
            periSWR_spikes{VR_ephys}{u} = [periSWR_spikes{VR_ephys}{u},spikes_in_epoch]; %concatenate the spikes from this SWR_unamb epoch to that of the previous SWR_unamb epochs   
            %construct perinonSWR spike train
            spikes_in_epoch = ks_all_data{VR_ephys,1}{u,2}(and(ks_all_data{VR_ephys,1}{u,2} < nonSWR_epochs{VR_ephys}(i,1)+peri_window,...
                ks_all_data{VR_ephys,1}{u,2} >= nonSWR_epochs{VR_ephys}(i,1)-peri_window)); %find all the spikes from this unit that occurred in a pre-SWR epoch
            spikes_in_epoch = spikes_in_epoch-(nonSWR_epochs{VR_ephys}(i,1)-peri_window-1)+st_delimiter2; %shift the spikes so you can concatenate them into one spike train
            perinonSWR_spikes{VR_ephys}{u} = [perinonSWR_spikes{VR_ephys}{u},spikes_in_epoch]; %concatenate the spikes from this nonSWR epoch to that of the previous nonSWR epochs   
            st_delimiter = st_delimiter+2*peri_window; %calculate the total duration of the spike train thus far (for peri)
        end
    end
    VR_ephys
end
%% Correlations of preSWR vs. postSWR epochs
converter_original = 2.^(0:99);
all_bin_sizes = [5,10,50,100];
all_word_lengths = 16; %[4:2:24]; %even though word lengths aren't explicitly necessary for calculating FR/corr, it tells us which small-population sessions to exclude to keep it consistent with entropy/KLD analyses
x = -40:40; sd = 10; GG = exp(-.5*(x/sd).^2)/(sd*sqrt(2*pi));
mod_rasters_FR_corr_preSWR = cell(size(VR_reformatting,1),1); 
mod_rasters_FR_corr_postSWR = mod_rasters_FR_corr_preSWR;
mod_rasters_FR_corr_prenonSWR = mod_rasters_FR_corr_preSWR;
mod_rasters_FR_corr_postnonSWR = mod_rasters_FR_corr_preSWR;
FR_corr_VR_frames_preSWR = cell(size(VR_reformatting,1),1); 
corr_consistent_preSWR = cell(size(VR_reformatting,1),numel(converter_original),numel(all_bin_sizes));
corr_consistent_postSWR = corr_consistent_preSWR;
corr_consistent_periSWR = corr_consistent_preSWR;
corr_consistent_prenonSWR = corr_consistent_preSWR;
corr_consistent_postnonSWR = corr_consistent_preSWR;
corr_consistent_perinonSWR = corr_consistent_preSWR;
FR_consistent_preSWR = corr_consistent_preSWR;
FR_consistent_postSWR = corr_consistent_preSWR;
FR_consistent_prenonSWR = corr_consistent_preSWR;
FR_consistent_postnonSWR = corr_consistent_preSWR;
for bin_counter = 2
    FR_corr_bin_size = all_bin_sizes(bin_counter); 
    for word_length_counter = 1:numel(all_word_lengths)
        word_length = all_word_lengths(word_length_counter); 
        for VR_ephys = 1:size(VR_reformatting,1) 
            for i = 1:size(ks_all_data{VR_reformatting{VR_ephys,5},1},1)
                all_st = preSWR_spikes{VR_ephys}{i}; %all preSWR spike times (in 30kHz) for a particular unit
                in_interval_st = round(all_st/30); %spike times in ms
                in_interval_st(in_interval_st == 0) = []; %eliminate spikes at time 0
                mod_rasters_FR_corr_preSWR{VR_ephys}{i,1} = in_interval_st; 
                all_st = postSWR_spikes{VR_ephys}{i}; %all postSWR spike times (in 30kHz) for a particular unit
                in_interval_st = round(all_st/30); %spike times in ms
                in_interval_st(in_interval_st == 0) = []; %eliminate spikes at time 0
                mod_rasters_FR_corr_postSWR{VR_ephys}{i,1} = in_interval_st; 
                all_st = prenonSWR_spikes{VR_ephys}{i}; %all prenonSWR spike times (in 30kHz) for a particular unit
                in_interval_st = round(all_st/30); %spike times in ms
                in_interval_st(in_interval_st == 0) = []; %eliminate spikes at time 0
                mod_rasters_FR_corr_prenonSWR{VR_ephys}{i,1} = in_interval_st; 
                all_st = postnonSWR_spikes{VR_ephys}{i}; %all postnonSWR spike times (in 30kHz) for a particular unit
                in_interval_st = round(all_st/30); %spike times in ms
                in_interval_st(in_interval_st == 0) = []; %eliminate spikes at time 0
                mod_rasters_FR_corr_postnonSWR{VR_ephys}{i,1} = in_interval_st; 
            end                                             
            if size(mod_rasters_FR_corr_preSWR{VR_ephys},1) >= word_length
                spike_train_binned_all_preSWR_tmp = zeros(size(mod_rasters_FR_corr_preSWR{VR_ephys},1),(peri_window/30)*size(SWR_unamb_epochs{VR_ephys},1)); 
                spike_train_binned_all_postSWR_tmp = spike_train_binned_all_preSWR_tmp; 
                spike_train_binned_all_preSWR = spike_train_binned_all_preSWR_tmp; 
                spike_train_binned_all_postSWR = spike_train_binned_all_preSWR_tmp; 
                spike_train_binned_all_prenonSWR_tmp = spike_train_binned_all_preSWR_tmp; 
                spike_train_binned_all_postnonSWR_tmp = spike_train_binned_all_preSWR_tmp; 
                spike_train_binned_all_prenonSWR = spike_train_binned_all_preSWR_tmp; 
                spike_train_binned_all_postnonSWR = spike_train_binned_all_preSWR_tmp; 
                for units = 1:size(mod_rasters_FR_corr_preSWR{VR_ephys},1) %for each unit 
                    spike_train_binned_all_preSWR_tmp(units,mod_rasters_FR_corr_preSWR{VR_ephys}{units,1}) = 1;
                    tmp = conv(GG,spike_train_binned_all_preSWR_tmp(units,:)); 
                    spike_train_binned_all_preSWR(units,:) = tmp(1:size(spike_train_binned_all_preSWR(units,:),2)); 
                    spike_train_binned_all_postSWR_tmp(units,mod_rasters_FR_corr_postSWR{VR_ephys}{units,1}) = 1;
                    tmp = conv(GG,spike_train_binned_all_postSWR_tmp(units,:)); 
                    spike_train_binned_all_postSWR(units,:) = tmp(1:size(spike_train_binned_all_postSWR(units,:),2)); 
                    spike_train_binned_all_prenonSWR_tmp(units,mod_rasters_FR_corr_prenonSWR{VR_ephys}{units,1}) = 1;
                    tmp = conv(GG,spike_train_binned_all_prenonSWR_tmp(units,:)); 
                    spike_train_binned_all_prenonSWR(units,:) = tmp(1:size(spike_train_binned_all_prenonSWR(units,:),2)); 
                    spike_train_binned_all_postnonSWR_tmp(units,mod_rasters_FR_corr_postnonSWR{VR_ephys}{units,1}) = 1;
                    tmp = conv(GG,spike_train_binned_all_postnonSWR_tmp(units,:)); 
                    spike_train_binned_all_postnonSWR(units,:) = tmp(1:size(spike_train_binned_all_postnonSWR(units,:),2)); 
                end                
                spike_train_binned_all_periSWR = zscore([spike_train_binned_all_preSWR,spike_train_binned_all_postSWR],0,2);
                tmp_corr = corr(spike_train_binned_all_periSWR'); tmp_corr = tmp_corr-eye(size(tmp_corr,1)).*tmp_corr; tmp_corr(logical(eye(size(tmp_corr,1)))) = 0;
                corr_consistent_periSWR{VR_ephys,word_length,bin_counter} = squareform(tmp_corr);                              
                FR_consistent_preSWR{VR_ephys,word_length,bin_counter} = ((sum(spike_train_binned_all_preSWR,2)/size(spike_train_binned_all_preSWR,2))*1000/FR_corr_bin_size)'; %MAY NEED TO RECALCULATE THIS IF YOU STICK WITH CONVOLUTION as opposed to binning
                spike_train_binned_all_preSWR = zscore(spike_train_binned_all_preSWR,0,2); 
                tmp_corr = corr(spike_train_binned_all_preSWR'); tmp_corr = tmp_corr-eye(size(tmp_corr,1)).*tmp_corr; tmp_corr(logical(eye(size(tmp_corr,1)))) = 0;
                corr_consistent_preSWR{VR_ephys,word_length,bin_counter} = squareform(tmp_corr); 
                FR_consistent_postSWR{VR_ephys,word_length,bin_counter} = ((sum(spike_train_binned_all_postSWR,2)/size(spike_train_binned_all_postSWR,2))*1000/FR_corr_bin_size)'; %MAY NEED TO RECALCULATE THIS IF YOU STICK WITH CONVOLUTION as opposed to binning                
                spike_train_binned_all_postSWR = zscore(spike_train_binned_all_postSWR,0,2); 
                tmp_corr = corr(spike_train_binned_all_postSWR'); tmp_corr = tmp_corr-eye(size(tmp_corr,1)).*tmp_corr; tmp_corr(logical(eye(size(tmp_corr,1)))) = 0;
                corr_consistent_postSWR{VR_ephys,word_length,bin_counter} = squareform(tmp_corr);       
                
                spike_train_binned_all_perinonSWR = zscore([spike_train_binned_all_prenonSWR,spike_train_binned_all_postnonSWR],0,2);
                tmp_corr = corr(spike_train_binned_all_perinonSWR'); tmp_corr = tmp_corr-eye(size(tmp_corr,1)).*tmp_corr; tmp_corr(logical(eye(size(tmp_corr,1)))) = 0;
                corr_consistent_perinonSWR{VR_ephys,word_length,bin_counter} = squareform(tmp_corr);                              
                FR_consistent_prenonSWR{VR_ephys,word_length,bin_counter} = ((sum(spike_train_binned_all_prenonSWR,2)/size(spike_train_binned_all_prenonSWR,2))*1000/FR_corr_bin_size)'; %MAY NEED TO RECALCULATE THIS IF YOU STICK WITH CONVOLUTION as opposed to binning
                spike_train_binned_all_prenonSWR = zscore(spike_train_binned_all_prenonSWR,0,2); 
                tmp_corr = corr(spike_train_binned_all_prenonSWR'); tmp_corr = tmp_corr-eye(size(tmp_corr,1)).*tmp_corr; tmp_corr(logical(eye(size(tmp_corr,1)))) = 0;
                corr_consistent_prenonSWR{VR_ephys,word_length,bin_counter} = squareform(tmp_corr); 
                FR_consistent_postnonSWR{VR_ephys,word_length,bin_counter} = ((sum(spike_train_binned_all_postnonSWR,2)/size(spike_train_binned_all_postnonSWR,2))*1000/FR_corr_bin_size)'; %MAY NEED TO RECALCULATE THIS IF YOU STICK WITH CONVOLUTION as opposed to binning                
                spike_train_binned_all_postnonSWR = zscore(spike_train_binned_all_postnonSWR,0,2); 
                tmp_corr = corr(spike_train_binned_all_postnonSWR'); tmp_corr = tmp_corr-eye(size(tmp_corr,1)).*tmp_corr; tmp_corr(logical(eye(size(tmp_corr,1)))) = 0;
                corr_consistent_postnonSWR{VR_ephys,word_length,bin_counter} = squareform(tmp_corr);   
            end
            FR_corr_bin_size
            VR_ephys
        end
    end
end
%% Plot preSWR vs. postSWR correlations
wl = 16; bin_count = 2;
num_subsamples = 100; 
subsampled_corr_postSWR = cell(size(corr_consistent_postSWR,1),1); 
subsampled_corr_preSWR = subsampled_corr_postSWR;
subsampled_FR_postSWR = subsampled_corr_postSWR;
subsampled_FR_preSWR = subsampled_corr_postSWR; 
subsampled_corr_postnonSWR = subsampled_corr_postSWR;
subsampled_corr_prenonSWR = subsampled_corr_postSWR;
subsampled_FR_postnonSWR = subsampled_corr_postSWR;
subsampled_FR_prenonSWR = subsampled_corr_postSWR; 
for VR_ephys = 1:size(corr_consistent_postSWR,1)                          
    if numel(corr_consistent_postSWR{VR_ephys,wl,bin_count}) >= num_subsamples
        rand_idx = randperm(numel(corr_consistent_postSWR{VR_ephys,wl,bin_count}),num_subsamples);
        subsampled_corr_postSWR{VR_ephys} = corr_consistent_postSWR{VR_ephys,wl,bin_count}(rand_idx); subsampled_corr_postSWR{VR_ephys}(isnan(subsampled_corr_postSWR{VR_ephys})) = []; 
        subsampled_corr_preSWR{VR_ephys} = corr_consistent_preSWR{VR_ephys,wl,bin_count}(rand_idx); 
        rand_idx = randi(numel(FR_consistent_postSWR{VR_ephys,wl,bin_count}),1,num_subsamples); 
        subsampled_FR_postSWR{VR_ephys} = FR_consistent_postSWR{VR_ephys,wl,bin_count}(rand_idx); 
        subsampled_FR_preSWR{VR_ephys} = FR_consistent_preSWR{VR_ephys,wl,bin_count}(rand_idx); 
        rand_idx = randperm(numel(corr_consistent_postnonSWR{VR_ephys,wl,bin_count}),num_subsamples);
        subsampled_corr_postnonSWR{VR_ephys} = corr_consistent_postnonSWR{VR_ephys,wl,bin_count}(rand_idx); subsampled_corr_postnonSWR{VR_ephys}(isnan(subsampled_corr_postnonSWR{VR_ephys})) = []; 
        subsampled_corr_prenonSWR{VR_ephys} = corr_consistent_prenonSWR{VR_ephys,wl,bin_count}(rand_idx); 
        rand_idx = randi(numel(FR_consistent_postnonSWR{VR_ephys,wl,bin_count}),1,num_subsamples); 
        subsampled_FR_postnonSWR{VR_ephys} = FR_consistent_postnonSWR{VR_ephys,wl,bin_count}(rand_idx); 
        subsampled_FR_prenonSWR{VR_ephys} = FR_consistent_prenonSWR{VR_ephys,wl,bin_count}(rand_idx); 
    end
end
g = figure; set(g,'Position',[0 -500 1000 500]); 
subplot(1,2,1); 
histogram(cell2mat(subsampled_corr_postSWR(WT_dCA1_sites)'),[-.2:.001:1],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(cell2mat(subsampled_corr_postSWR(AD_dCA1_sites)'),[-.2:.001:1],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
text(.5,.5,strcat(num2str(mean(cell2mat(subsampled_corr_postSWR(WT_dCA1_sites)'))),',',num2str(std(cell2mat(subsampled_corr_postSWR(WT_dCA1_sites)')))));
text(.5,.4,strcat(num2str(mean(cell2mat(subsampled_corr_postSWR(AD_dCA1_sites)'))),',',num2str(std(cell2mat(subsampled_corr_postSWR(AD_dCA1_sites)')))));
[p,~] = ranksum(cell2mat(subsampled_corr_postSWR(WT_dCA1_sites)'),cell2mat(subsampled_corr_postSWR(AD_dCA1_sites)'));
text(.5,.3,strcat('p = ',num2str(p))); 
ylim([0 1]); xlim([-.1 .8]); 
subplot(1,2,2); 
histogram(cell2mat(subsampled_corr_postSWR(WT_vCA1_sites)'),[-.2:.001:1],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(cell2mat(subsampled_corr_postSWR(AD_vCA1_sites)'),[-.2:.001:1],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
text(.5,.5,strcat(num2str(mean(cell2mat(subsampled_corr_postSWR(WT_vCA1_sites)'))),',',num2str(std(cell2mat(subsampled_corr_postSWR(WT_vCA1_sites)')))));
text(.5,.4,strcat(num2str(mean(cell2mat(subsampled_corr_postSWR(AD_vCA1_sites)'))),',',num2str(std(cell2mat(subsampled_corr_postSWR(AD_vCA1_sites)')))));
[p,~] = ranksum(cell2mat(subsampled_corr_postSWR(WT_vCA1_sites)'),cell2mat(subsampled_corr_postSWR(AD_vCA1_sites)'));
text(.5,.3,strcat('p = ',num2str(p))); 
ylim([0 1]); xlim([-.1 .8]); 
g2 = figure; set(g2,'Position',[0 -500 1000 500]); 
subplot(1,2,1); 
histogram(cell2mat(subsampled_corr_postnonSWR(WT_dCA1_sites)'),[-.2:.001:1],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','k'); hold on;
histogram(cell2mat(subsampled_corr_postnonSWR(AD_dCA1_sites)'),[-.2:.001:1],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor','r'); hold on;
text(.5,.5,strcat(num2str(mean(cell2mat(subsampled_corr_postnonSWR(WT_dCA1_sites)'))),',',num2str(std(cell2mat(subsampled_corr_postnonSWR(WT_dCA1_sites)')))));
text(.5,.4,strcat(num2str(mean(cell2mat(subsampled_corr_postnonSWR(AD_dCA1_sites)'))),',',num2str(std(cell2mat(subsampled_corr_postnonSWR(AD_dCA1_sites)')))));
[p,~] = ranksum(cell2mat(subsampled_corr_postnonSWR(WT_dCA1_sites)'),cell2mat(subsampled_corr_postnonSWR(AD_dCA1_sites)'));
text(.5,.3,strcat('p = ',num2str(p))); 
ylim([0 1]); xlim([-.1 .8]); 
subplot(1,2,2); 
histogram(cell2mat(subsampled_corr_postnonSWR(WT_vCA1_sites)'),[-.2:.001:1],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[.8 .8 .8]); hold on;
histogram(cell2mat(subsampled_corr_postnonSWR(AD_vCA1_sites)'),[-.2:.001:1],'Normalization','cdf','DisplayStyle','stairs','LineWidth',2,'EdgeColor',[1 .8 .8]); hold on;
text(.5,.5,strcat(num2str(mean(cell2mat(subsampled_corr_postnonSWR(WT_vCA1_sites)'))),',',num2str(std(cell2mat(subsampled_corr_postnonSWR(WT_vCA1_sites)')))));
text(.5,.4,strcat(num2str(mean(cell2mat(subsampled_corr_postnonSWR(AD_vCA1_sites)'))),',',num2str(std(cell2mat(subsampled_corr_postnonSWR(AD_vCA1_sites)')))));
[p,~] = ranksum(cell2mat(subsampled_corr_postnonSWR(WT_vCA1_sites)'),cell2mat(subsampled_corr_postnonSWR(AD_vCA1_sites)'));
text(.5,.3,strcat('p = ',num2str(p))); 
ylim([0 1]); xlim([-.1 .8]); 
cd('Y:\Udaysankar Chockanathan\Aim 1\Figures\Paper figures\Figure 10'); 
saveas(g,'SWR_corr','fig'); orient(g,'landscape'); print('SWR_corr.pdf','-bestfit','-painters','-dpdf'); close(g);
saveas(g2,'nonSWR_corr','fig'); orient(g2,'landscape'); print('nonSWR_corr.pdf','-bestfit','-painters','-dpdf'); close(g2);