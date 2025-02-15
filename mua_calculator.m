function [hpc, rsc, p_swr, down_start_time,up_start_time,SWR_indices] = mua_calculator(Y_sol)

    % run threshold calculator to identify UP and DOWN states
    [thresh,cross,bihist,diptest,crossup,crossdown] = BimodalThresh(Y_sol.y(1,:)','Schmidt');
    up_start_time = cross.upints(:,1);
    down_start_time = cross.downints(:,1);
    up_state_durations = cross.upints(:,2) - cross.upints(:,1);
    down_state_durations = cross.downints(:,2) - cross.downints(:,1);
    
    mean_mua = [];
    for i = 1:length(up_state_durations)
        % for each up state calculate the mean firing rate in it
        mean_mua = [mean_mua, mean(Y_sol.y(1,up_start_time(i):up_start_time(i)+up_state_durations(i)))];
    end

    % run fixed thresholding for SWR events
    threshold = 2.5;
    SWR_indices = find(Y_sol.y(3,:) > threshold);
    
    %SWR_indices = T(SWR_indices);
    if isempty(SWR_indices)
            result = 'no swrs';
            SWR_d = 0
            iSWR_d = 0
            return;
    end
    SWR_start = find(diff(SWR_indices) > 0);
    
    SWR_beginning = SWR_indices(find(diff(SWR_indices) > 1)+1);
    SWR_beginning = transpose([SWR_indices(1);transpose(SWR_beginning)]);
    
    SWR_ending = SWR_indices(find(diff(SWR_indices) > 1));
    SWR_ending = transpose([transpose(SWR_ending);SWR_indices(end)]);
    
    SWR_timestamps = SWR_indices(SWR_start);
    ups = ones(size(SWR_timestamps))-0.4;
    SWR_durations = SWR_ending-SWR_beginning;
    
    if SWR_beginning(1)<SWR_ending(1)
        iSWR_durations = zeros(length(SWR_beginning),1);
        for i = 1:length(SWR_beginning)-1
            iSWR_durations(i) = SWR_beginning(i+1)-SWR_ending(i);
        end
    elseif SWR_beginning(1)>SWR_ending
        iSWR_durations = SWR_beginning-SWR_ending;
    else 
        disp('no SWRS')
    
    end 
    
    
    SWR_peak = zeros(length(SWR_beginning), 1);  % Initialize SWR peak timestamps
    
    for i = 1:length(SWR_beginning)
        % Extract firing rates within the current SWR interval
        swr_interval_firing_rates = Y_sol.y(3,SWR_beginning(i):SWR_ending(i));
    
        % Find the timestamp of the peak firing rate within the interval
        [~, max_idx] = max(swr_interval_firing_rates);
        SWR_peak(i) = SWR_beginning(i) + max_idx - 1;  % Adjust index to global timestamp
    end

    % find timestamps of pairs of UP and preceding DOWN states

    preceding_down = zeros(size(up_start_time));
    preceding_duration = zeros(size(up_start_time));
    
    for j = 1:numel(up_start_time)
        ups = up_start_time(j);  % Corrected indexing
        preceding_index = find(down_start_time < ups, 1, 'last');
        
        % if found assign the corresponding down start time val
        if ~isempty(preceding_index)
            preceding_down(j) = down_start_time(preceding_index);
            preceding_duration(j) = down_state_durations(preceding_index);
        else 
            preceding_down(j) = NaN;
            preceding_duration(j) = NaN;
        end
    end
    prec_down_start = preceding_down;
    prec_down_duration = preceding_duration;

    % find average MUA and pSWR in each UP/DOWN pair

    num_bins = 15;
    p_swr = zeros(1,num_bins*2);
    rsc = zeros(1,num_bins*2);
    hpc = zeros(1,num_bins*2);
    swr_counts_D = zeros(length(prec_down_start),num_bins);
    mean_ctx_D = zeros(length(prec_down_start),num_bins);
    mean_hpc_D = zeros(length(prec_down_start),num_bins);
    for i = 1:length(prec_down_start)
        if isnan(prec_down_start(i))
            continue;  % Skip to the next iteration of the loop
        end
        % calculate time bins for the current UP state
        time_bins = linspace(prec_down_start(i),prec_down_start(i)+ prec_down_duration(i), num_bins+1);
        %iterate over each time bin
        for j = 1:num_bins
            % find the swrs within the current time bin
            swr_in_bin = SWR_indices >= time_bins(j)& SWR_indices < time_bins(j+1);
            % count the number of swrs in the current time bin
            swr_counts_D(i,j) = sum(swr_in_bin);
    
            % calculate the mean firing rate in the current time bin
            CTX_rates_in_bin = Y_sol.y(1,round(time_bins(j)):round(time_bins(j+1)));
            mean_ctx_D(i,j) = mean(CTX_rates_in_bin);
            HPC_rates_in_bin = Y_sol.y(3,round(time_bins(j)):round(time_bins(j+1)));
            mean_hpc_D(i,j) = mean(HPC_rates_in_bin);
        end
    end
    swr_counts_U = zeros(length(up_start_time),num_bins);
    mean_ctx_U = zeros(length(up_start_time),num_bins);
    mean_hpc_U = zeros(length(up_start_time),num_bins);
    for i = 1:length(up_start_time)
        % calculate time bins for the current UP state
        time_bins = linspace(up_start_time(i),up_start_time(i)+ up_state_durations(i), num_bins+1);
    
        %iterate over each time bin
        for j = 1:num_bins
            % find the swrs within the current time bin
            swr_in_bin = SWR_indices >= time_bins(j)& SWR_indices < time_bins(j+1);
            % count the number of swrs in the current time bin
            swr_counts_U(i,j) = sum(swr_in_bin);
    
            % calculate the mean firing rate in the current time bin
            CTX_rates_in_bin = Y_sol.y(1,round(time_bins(j)):round(time_bins(j+1)));
            mean_ctx_U(i,j) = mean(CTX_rates_in_bin);
            HPC_rates_in_bin = Y_sol.y(3,round(time_bins(j)):round(time_bins(j+1)));
            mean_hpc_U(i,j) = mean(HPC_rates_in_bin);
        end
    end
    
    rsc_mua = [mean(mean_ctx_D(:,1:end),1),mean(mean_ctx_U(:,1:end),1)];
    hpc_mua = [mean(mean_hpc_D(:,1:end),1),mean(mean_hpc_U(:,1:end),1)];
    p_distribution = [sum(swr_counts_D(:,1:end),1)/length(SWR_peak),sum(swr_counts_U(:,1:end),1)/length(SWR_peak)];
    p_swr(:) = p_distribution;
    rsc(:) = rsc_mua;
    hpc(:) = hpc_mua;
end