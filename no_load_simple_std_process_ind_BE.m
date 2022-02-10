% LMO no-load event detection
cd 'G:\Dunn lab OT data'
NLB_tags = {'*/*NL*.mat','*/*NF*.mat','*/*noload*.mat','*/*noforce*.mat'};
output_dir = 'G:\Dunn lab OT data\HT-R13DD_binding_events\no_load\revision\output_dir';

% fnames = [];
% for f=1:numel(NLB_tags)
%     fnames = [fnames; dir(NLB_tags{f})];
% end

for f = 1:numel(fnames)
    addpath(fnames(f).folder)
end
fnames = {fnames.name}';
mkdir(output_dir)
cd(output_dir)
%%
output_dir = pwd;
time_res = 10; %in kHz
thresh = 0.64;
smooth_window = 100; %20 ms
maybe_thresh = 0.8; %flood fill threshold
std_window = 15;

NLB_BE = []; NLB_ids = []; fnames_ids = [];
for f=1:numel(fnames)

    load(fnames{f})
    
    dec_num = 40/time_res;
    F = [decimate(Force.ForceOne.calibratedArray,dec_num) ; ...
    decimate(Force.ForceTwo.calibratedArray,dec_num)]';

    freq_trace = zeros(size(F'));
    for t=1:2
       freq_trace(t,:) = movmean(movstd(F(:,t),std_window).^2,smooth_window)';
       freq_trace(t,:) = freq_trace(t,:)/nanmean(freq_trace(t,:)); 
    end
    
    freq_trace = [freq_trace; mean(freq_trace,1)]; %third column is sum
    
    %isbound initially for each trap if below "thresh" for the single-trap
    %signals or below "thresh_both" for the mean trap signal. 
    isbound = zeros(size(freq_trace));
    isbound(freq_trace < thresh) = 1;
    
    %now processes each binding event to expand its edges until the noise
    %level hits "maybe_thresh"
    for t=1:2
        trace = freq_trace(t,:);
        bes = bwlabel(isbound(t,:));
        for be=1:max(bes)
            %find the "end" of the binding event when the trace first
            %reaches maybe_thresh
            id = find(bes==be,1);
            trace_consider = trace; trace_consider(1:id)=0;
            id_end = find(trace_consider > maybe_thresh,1,'first');
            
            trace_consider = trace; trace_consider(id:end) = 0;
            id_start = find(trace_consider > maybe_thresh,1,'last');
            
            isbound(t,id_start:id_end) = 1;
        end
        
        %if this is the trace for the first trap, get rid of all of the second trap
        %binding events that aren't detected here to make downstream
        %processing faster.
        if t==1
           not_bound = ~isbound(1,:);
           isbound(2,not_bound) = 0; %prelim filter to make next step faster 
        end
    end

   isbound_all = isbound(1,:) & isbound(2,:);
   isbound = isbound_all; clear isbound_all
   
    Fbound1 = F(:,1); Fbound1(~isbound) = NaN; %just for plotting
    Fbound2 = F(:,2); Fbound2(~isbound) = NaN;
    
    BE_lifetime = regionprops(isbound, 'Area');
    BE_lifetime = [BE_lifetime.Area]/time_res;
    
    X1 = 0;
    X2 = size(Force.ForceOne.calibratedArray,2)/40;
    N = size(F,1);
    time_ms = linspace(X1, X2, N);
    
    close
    figure,
    ax1 = subplot(5,1,1);
    plot(time_ms, F(:,1),'.'), hold on
    plot(time_ms, Fbound1,'.')
    title(fnames{f})
    
    ax2 = subplot(5,1,2);
    plot(time_ms, F(:,2),'.'), hold on
    plot(time_ms, Fbound2,'.')
    
    ax3 = subplot(5,1,3);
    plot(time_ms, freq_trace(1,:)), hold on
    freq_trace(1,~isbound) = NaN;
    plot(time_ms, freq_trace(1,:));

    ax4 = subplot(5,1,4);
    plot(time_ms, freq_trace(2,:)), hold on
    freq_trace(2,~isbound) = NaN;
    plot(time_ms, freq_trace(2,:));

    linkaxes([ax1,ax2, ax3, ax4],'x')
    
    ax5 = subplot(5,1,5);
    histogram(BE_lifetime,'BinWidth',10)
    
    savefig([output_dir '\' fnames{f}(1:end-4) '_NLB.fig'])
    
    %characterize files by filename, but manually check and fix
    %mis-classifications
    if contains(fnames{f},'neg')
        NLB_id = 0*ones(size(BE_lifetime));
    elseif contains(fnames{f},'R13-') || contains(fnames{f},'GHRdelta')
        NLB_id = 1*ones(size(BE_lifetime));
    elseif contains(fnames{f},'GHRDD') | contains(fnames{f},'x2')
        NLB_id = 2*ones(size(BE_lifetime));
    elseif contains(fnames{f},'HT-R13DD') | contains(fnames{f},'MR13DD')
        NLB_id = 3*ones(size(BE_lifetime));
    else
        NLB_id = 4*ones(size(BE_lifetime));
    end
    
    NLB_BE = [NLB_BE, BE_lifetime];
    NLB_ids = [NLB_ids, NLB_id];
    fnames_ids = [fnames_ids,f*ones(size(NLB_id))];
    clear BE_lifetime
    clear NLB_id
end
%%
final_results = struct;
final_results.params.method = 'var';
% final_results.params.freq_morlet_modes = [2,6];
final_results.params.fnames = fnames;
final_results.params.thresh = thresh;
final_results.params.time_res_kHz = time_res;
final_results.params.smooth_window_pix = smooth_window;
final_results.params.maybe_thresh = maybe_thresh;
final_results.params.std_window = std_window;

final_results.results.all.NLB_BE = NLB_BE;
final_results.results.all.NLB_ids = NLB_ids;
final_results.results.all.fname_ids = fnames_ids;
final_results.results.all.fname = fnames;

final_results.results.R13DD = NLB_BE(NLB_ids==3);
final_results.results.GHR13DD = NLB_BE(NLB_ids==2);
final_results.results.R13andGHR13 = NLB_BE(NLB_ids==1);
final_results.results.neg = NLB_BE(NLB_ids==0);
