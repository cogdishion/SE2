%Expected to be called by function bootstrap_SpeakEasy many times, in
%order to generate replicate partitions
%
%This is the fundamental SpeakEasy2 routine that counts up actual labels, compares to
%predicted, and then selects most unexpected label for node identifier.  It has 4 variations on this.
%1) there is the basic mode where popular labels propagate,
%2) At some times it can split up ('bubble') labels into different sets,
%3) it can nurture these bubbles to try to form stable cluster and finally
%4) it can also merge two distinct labels into a single label if they have extensive connections
%
% you might be annoyed that the code looks more complicated than it needs
% to - the need for that comes from the fact we only want to update various subsets
% of nodes, which makes the indexing for counting up labels and assigning them to nodes a pain


function [partitionID secondary_labels_scores  secondary_labels_ID  max_labels_output ]=SpeakEasy2_core(ADJ,IC_store_relevant,main_iter,kin,ktot,options,varargin)


%set up labels and storage of labels for each node
listener_history=zeros(options.timesteps,length(ADJ));  %each column is history represents a single node, and each row is a timestep
listener_history(1,:)=IC_store_relevant;  %load up random initial conditions
score_history=zeros(options.timesteps,length(ADJ));  %each column is history represents a single node, and each row is a timestep

%prep to record bubbling and fusing - just tracking - doesn't affect operation
bubble_history=[1 zeros(1,size(listener_history,1))];  %later on this will hold a 1 when bubbling occured
merge_history= [1 zeros(1,size(listener_history,1))];  %later on this will hold a 1 when fusion occures

% prep storage for solutions
preintervention_state=[];
postintervention_state=[];

%switches that track bubbling
max_labels_this_bubbling_sesion=[]; %always initialize to empty matrix - will determine when to stop bubbling phase
splits_per_phase_history=0;  %initial value 0 that is appended with number of splits
splitting_has_peaked_state=0;  %this always be initialized to zero and it's value indicates how many splits have occurred after the first naturally occuring decrease in the number of labels, indicating any further splits are not necessary, as they've already had max effect
post_peak_split_counter=0;     %should always be initialize to zero since it's a counter
post_peak_split_limit=2;       %split this many times after #communities appears to have peaked- can even set to zero for slight speed improvment, but goal is just o make sure bubbled well. No advantage>2.

%switches that track merging
possible_cut=[]; %intialize empty,  holds a variable set dynamically for "do_fusion" section defining the min module-level overlap sufficient to merge a pair of modules
switch_to_fusion_mode=0;  %always initialize to zero - this toggles behavior between meta clustering or splitting

%if you have enabled overlapping clusters with options.multicommunity>1
if options.multicommunity>1
    secondary_labels_scores=zeros(options.multicommunity,length(ADJ),options.target_partitions +options.discard_transient );
    secondary_labels_ID=    zeros(options.multicommunity,length(ADJ),options.target_partitions +options.discard_transient ); %todel
    subop_tracker=1;
    
else  %disjoint clustering
    secondary_labels_ID=[];
    secondary_labels_scores=[];
end

%useful for figures  - could ultimately be removed as long as you also delete refs to these
number_of_labels=[];    %for each time point counts #unique lables
median_cluster_size=[]; %for each time point counts medians size of all clusters
mean_cluster_size=[];
max_labels_output=[];  %for figure
adjustedrand_all_time=zeros(1,options.timesteps);  %only used diagnostically if an official solution variable 'official' is loaded
nmi_all_time=zeros(1,options.timesteps);

%useful for figures - stats on evolving solutions
number_of_labels=[];    %for each time point counts #unique lables
median_cluster_size=[]; %for each time point counts medians size of all clusters



i=1; %counter for main loop, which continues until generating the target_paritions # of solutions (aka clusters aka partitions)
while size(postintervention_state,1) < (options.target_partitions +options.discard_transient )%&&post_peak_split_counter<=  %main label updating loop
    tic;
    
    i=i+1;
    listener_history(i,:)=listener_history(i-1,:);  %most labesl in listener_history(i,:) will be updated
    nodes_to_update=1:length(ADJ);  %these next four lines select a unique set of nodes that will be updated at this timestsep
    to_ignore_frac=.1;   %  this is the randomly selected fraction of nodes that are not updated .5 to .1 all work about as well on LFR, but the more you ignore, the longer it takes to complete
    nodes_to_ignore=randperm(length(ADJ),ceil(to_ignore_frac*length(ADJ)));
    nodes_to_update(nodes_to_ignore)=[];
    
    
    if size(listener_history,1)-i<300  %we don't know total run time, so keep adding preallocation chunks
        listener_history=cat(1,listener_history,zeros(1000,size(listener_history,2)));
        bubble_history=cat(2, bubble_history, zeros(1,1000));
        merge_history= cat(2, merge_history,  zeros(1,1000));
    end
    
    
    if i>0
        labels_unq=unique(listener_history(i,:));
        number_of_labels(end+1)=length(labels_unq);
        %evaluting stability of labels, so need some time to pass - specific value not very important, just a few timesteps, can even be zero
        % cluster_stability_measure(i)=(number_of_labels(i-1)-number_of_labels(i))/number_of_labels(i); %stability comparison between points as close as i and i-2 works
    end
    
    %determine which of four things to do at each timestep
    if i<=20   %just do typical updating at first few steps... value could be 5+ doesn't really matter
        do_bubbling=0;
        do_fusion=0;
        do_typical=1;
        do_nurture=0;
        record_mode(i)=5;  %just a tracking variable to see what is done at what timestep
        
    else   %after some initial period of time, select one of four options
        do_bubbling=0;
        do_fusion=0;
        do_typical=0;
        do_nurture=0;
        record_mode(i)=4;
        
        if switch_to_fusion_mode==0  %i.e. do bubbling or nurture because you're NOT merging/fusing
            
            %don't bubble too soon after merging       and  wait a bit since any prior bubbling
            if i-max(find(merge_history==1))>2   &&   i-max(find(bubble_history==1))>=15   %bubbling
                do_bubbling=1;
                do_fusion=0;
                do_typical=0;
                do_nurture=0;
                record_mode(i)=1;
                
            elseif    i-max(find(merge_history==1))>=2   &&  i-max(find(bubble_history==1))<=4    %nurture-style updating
                do_bubbling=0;
                do_fusion=0;
                do_typical=0;
                do_nurture=1;
                record_mode(i)=2;
            end
            
        elseif    switch_to_fusion_mode==1    &&  i-max(find(merge_history==1))>=2   &&  i-max(find(bubble_history==1))>3% %fusion - cluster-level merging
            do_bubbling=0;
            do_fusion=1;
            do_typical=0;
            do_nurture=0;
            record_mode(i)=3;
        end
        
        
        if  sum([do_bubbling do_fusion do_nurture])==0;   %if not merging or bubbling, just do typical
            do_bubbling=0;
            do_fusion=0;
            do_typical=1;
            do_nurture=0;
            record_mode(i)=4;
        end
    end
    
    
    
    
    
    
    
    
    
    
    
    
    if do_typical==1
        
        %'actual_counts' contains have signed number of total receipts of each label, so if you have 5 positive incoming connections from label #2, and 4 negative inputs for label#2, actual_counts of label #2 will ==1
        [actual_counts active_labels]=actual_label_counts(ADJ,listener_history(i-1,:),options); %actual_counts is sparse for sparse input
        expected_counts=expected_label_counts(actual_counts,kin,0);%seems was set to 1 for lfr testing %set param ==1 for fullADJ increases speed... pretty much makes sense   %last parameter ==0  means expected ONLY for labels received by a given node (faster than providing expected values for all labels, including ones a node never hears)
        act_less_exp=(actual_counts-expected_counts);
        [best_act_less_exp_score idx_best_label]=max(act_less_exp,[],1);
        listener_history(i,nodes_to_update)= active_labels(idx_best_label(nodes_to_update));
        score_history(i,nodes_to_update)=best_act_less_exp_score(nodes_to_update);
        
        [active_labels act_less_exp]=adjust_for_ignored_labels(listener_history(i,:),active_labels,act_less_exp);  %this line matters a lot
        %Explaination of why the line above and the "adjust_for_ignored_labels" function matters a lot:
        %we generate updated labels for all nodes, but we only apply some
        %of them, selected by 'to_update_pseudo'.
        %Some labels may carry over from the pervious timestep but never appear in
        %active_labels -i.e. it will have more unique elements than listener_history(i,:), although
        %listener_history will have some labels that don't appear in
        %"active_labels", because nodes with those labels were not updated      
        
    end  %end do_typical
    
    
    
    
    
    
    
    
    
    
    if do_nurture==1 % will only update unstable nodes - typically done near times of bubbling
        %for future code update: probably clearer and faster to just update all
        %nodes, but then rever to old label assignments for stable nodes...
        %here we're trying to save ime by not calculating updates for them,
        %but they are few in number, and started doing this when thought
        %maybe we'd only be updating 10% of nodes, but now we're doing 90%
        %so doesn't make as much sense
        
        %only update some labels
        percentile_cut_point=.9;  % if this value is .9 means we DO NOT update the top 10% most secure nodes i.e. those with clear label selection
        percentile_marker=sort(best_act_less_exp_score);
        percentile_marker=percentile_marker(ceil(length(percentile_marker)*percentile_cut_point));
        low_fit_nodes_all= find(best_act_less_exp_score<percentile_marker ) ;
        
        if length(low_fit_nodes_all)==0  %sometimes , in small clusters or subclusters, you can get ties in "best_act_less_exp_score" such that all or no nodes end up in "low_fit_nodes_all"
            
            %in this case we just arbitrarily assign some nodes to be low or high fit.  Ratio doesn't seem to matter.  Could also just update all nodes.
            %disp('all nodes end uphigh selected')
%            low_fit_nodes_all=randsample(1:length(ADJ),ceil(length(ADJ)/2),1); %original, but calls randperm from stats toolbox
            low_fit_nodes_all=randperm(length(ADJ));
            low_fit_nodes_all=low_fit_nodes_all(1:ceil(length(ADJ)/2));
            high_fit_nodes_all=1:length(ADJ); %faster this way than setdiff
            high_fit_nodes_all(low_fit_nodes_all)=[];
            
        else
            high_fit_nodes_all=1:length(ADJ);
            high_fit_nodes_all(low_fit_nodes_all)=[];  %when all nodes are high fit (i.e. no low_fit nodes), active labels comes back empty
            %since high_fit_nodes are relatively rare in future update might be slightly faster to select them at first
        end
        
        %reiterating the above - sometimes you get many similar values in best_act_less_exp_score - so when you
        %take the nth percentile value as opposed to the nth value by idx
        %you can either put everything in low_fit_nodes_all (using <) or
        %everything if <=.  If you put nothing in, then
        %high_fit==length(adj) and all nodes get removed in
        %actual_label_counts and you get and empty return
        
        kin_limited_alt=kin-sum(ADJ(high_fit_nodes_all,:));
        [actual_counts active_labels]=actual_label_counts(ADJ,listener_history(i-1,:),options,high_fit_nodes_all); %actual_counts is sparse for sparse input
        expected_counts=expected_label_counts(actual_counts,kin_limited_alt,0);%seems was set to 1 for lfr testing %set param ==1 for fullADJ increases speed... pretty much makes sense   %last parameter ==0  means expected ONLY for labels received by a given node (faster than providing expected values for all labels, including ones a node never hears)
        act_less_exp=(actual_counts-expected_counts);
        [best_act_less_exp_score idx_best_label]=max(act_less_exp,[],1);     
        listener_history(i,nodes_to_update)= active_labels(idx_best_label(nodes_to_update));
        score_history(i,nodes_to_update)=best_act_less_exp_score(nodes_to_update);
        
        listener_history(i,high_fit_nodes_all)=listener_history(i-1,high_fit_nodes_all);   %we DO NOT take updated labels for high_fit_nodes_all aka nodes with strong label selection    
        [active_labels act_less_exp]=adjust_for_ignored_labels(listener_history(i,:),active_labels,act_less_exp);  %this line matters a lot fore reasons stated in paragraph in do_typical       
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    if do_bubbling==1
        low_fit_nodes_all=[];
        
        i=i+1; %increment "i" because we're artificially (via bubbling) creating next set of labels
        listener_history(i,:)=listener_history(i-1,:);  %most of these will be updated, but ome will not, so just have them continue from prior timestep
        bubble_history(i)=1; %used in determining what mode to go into
        
        splits_per_phase_history(end+1)=splits_per_phase_history(end)+1;  %not used - just for visualization - record how many times we bubble during each bubbling phase
        number_of_labels(end+1)=length(labels_unq);
        preintervention_state(end+1,:)=listener_history(i,:);
        
        [label_unq_loc todel]=    splitlist(listener_history(i,:));
        median_cluster_size(end+1)=median(cellfun(@length,label_unq_loc));  % NOT just diagnostic - used to determine how much to split big clusters
        prior_biggest_label=max(listener_history(i,:));
        
        percentile_marker_bubble=sort(best_act_less_exp_score);  %we get best_act_less_exp_score from when it is computed in do_nuture stage
        percentile_cut_point_bubble=.9;  %if value is ==.9 we're going to select nodes in the lowes 90% of maxval_values and split them into new clusters
        percentile_marker_bubble=percentile_marker_bubble(ceil(length(percentile_marker_bubble)*percentile_cut_point_bubble));
        
        %for each clusters we're going to take all the low fit nodes and assign them (in random set) to novel labels
        for k=1:length(labels_unq)
            if length(label_unq_loc{k})>4  %don't split small clusters - tried also filtering based on cluster density and didn't help
                
                low_fit_nodes=label_unq_loc{k}(  find(best_act_less_exp_score(label_unq_loc{k})  <=percentile_marker_bubble  ));
                low_fit_nodes_all=[low_fit_nodes_all; low_fit_nodes'];
                
                if length(low_fit_nodes)>0
                    
                    split_labels_to_insert=randi(prior_biggest_label+[1   max(2,  min([10 ceil([ length(label_unq_loc{k})/median_cluster_size(end)]) ])   )] , length(low_fit_nodes)  , 1 );
                    listener_history(i,low_fit_nodes)=split_labels_to_insert;
                    prior_biggest_label=max(split_labels_to_insert); % want to insert new labels bigger than this
                end
            end
        end
        
        max_labels_this_bubbling_sesion(end+1)=length(unique(listener_history(i,:)));  %working with the #labels a few time steps after bubbling might be more logical, but this acts siilarly and is easier to read
        
        %if statements below determine if we should stop bubling
        if splitting_has_peaked_state==0
            %this above if statement might not ben eeded
            if length(max_labels_this_bubbling_sesion)>2  &&   ~isempty(find(max_labels_this_bubbling_sesion(1:end-1)>.9*max_labels_this_bubbling_sesion(end)))   %.8 is just pragmatic - could do 1.0 for even-handedness, but not necesary
                %rationale for above line: check if we are NOT producing more subclusters than previous time point, becuase if we are NOT, then we DON'T need to bubble any more
                splitting_has_peaked_state=1;
                max_labels_output(end+1)=max(max_labels_this_bubbling_sesion);   %just for diagnostic experiments
            end
        end
        
        if splitting_has_peaked_state==1
            post_peak_split_counter=post_peak_split_counter+1;
            if post_peak_split_counter>=post_peak_split_limit
                switch_to_fusion_mode=1;  %means stop bubling annd start fusing
                post_peak_split_counter=0;%reset for next round of bubbling
            end
        end
        
        possible_cut_this_round=[];
        
    end  %end bubbling
    
    
    
    
    
    
    
    
    
    
    if do_fusion==1
        
        i=i+1;% increment i because we're going to record down merged labels as the next row in listener_history, without updating ID's in  the usual way
        number_of_labels(end+1)=length(labels_unq);
        merge_history(i)=1;  % used in decisions about which mode to go into
        listener_history(i,:)=listener_history(i-1,:);  %most of these will be updated, but this fills in non-updated nodes
        splits_per_phase_history(end)=0;   %reset split counter for next time in do_bubblinig
        post_peak_split_counter=0; %reset to zero but don't use in this section - used in bubbling
        splitting_has_peaked_state=0;
        max_labels_this_bubbling_sesion=[];
        
        [actual_counts active_labels]=actual_label_counts(ADJ,listener_history(i,:),options); %actual_counts is sparse for sparse input
        expected_counts=expected_label_counts(actual_counts,kin,0); %set to zero experimentally for large network   %last parameter ==1  means compute expected for all labels, not just labels received, which (may?) matter in do_fusion phase
        actual_by_group_compare=[actual_label_counts(actual_counts',listener_history(i,:),options)]';
        expected_by_group_compare=expected_label_counts(actual_by_group_compare,sum(actual_by_group_compare),0);%note last param=0 is the low ram version
        
        act_less_exp_group=actual_by_group_compare-expected_by_group_compare;
        act_less_exp_group_diag=diag(act_less_exp_group);
        act_less_exp_group(1:length(act_less_exp_group)+1:numel(act_less_exp_group))=-Inf;    %diagonal contains within label connections, which are not relevant to merging
        [merge_target_val merge_target_idx]=max(act_less_exp_group,[],1);  %get the best non-self community
        %the "keyvals" matrix below is holding stats that determine which clusters will be merged
        keyvals=sortrows([merge_target_val(:) merge_target_idx(:) [1:length(merge_target_val)]' act_less_exp_group_diag(merge_target_idx(:))  act_less_exp_group_diag],-1  );  %ranks labels  with most overlap
        %col2 of keyvals is the idx of best merge partner
        %col3 of keyvals is just a list of sequential indices used to references all labels - basically the idx of the listening labels
        
        
        
        [label_idx_in_cells  listener_history(i,:) active_labels] =splitlist(listener_history(i,:));   %bucket_labels
        label_number_at_time=cellfun(@length, label_idx_in_cells);
        label_number_at_time_for_certain_nodes_time0=label_number_at_time(keyvals(:,2));
        label_number_at_time_for_certain_nodes_timeplus=label_number_at_time(keyvals(:,3));
        keyvals(:,1) =[keyvals(:,1)'./(label_number_at_time_for_certain_nodes_time0+label_number_at_time_for_certain_nodes_timeplus)]';
        possible_cut(end+1)=mean(keyvals(find(keyvals(:,1)>0),1));
        
        
        keyvals=sortrows(keyvals,-1);
        to_merge=find(keyvals(:,1)>possible_cut(end));
        merge_pairs=[];
        been_merged=zeros(1,size(keyvals,1)); %this will track labels that have been merged - we only merge a label once every time we enter fusion mode
        keyval_med=median(keyvals(:,1));
        
        
        if length(to_merge)==0  || (length(possible_cut)~=1 && .5*possible_cut(end-1)<=possible_cut(end)) %   ||         to_merge_length_record*1.0<length(to_merge);  %this will be true when there are no cluster overlaps greater than expected, which will happen after we've previously merged for a while
            %  to_merge_length_record=length(ADJ);
            switch_to_fusion_mode=0;  %stop merging/fusing
            postintervention_state(end+1,:)=listener_history(i,:);  %record the
            
            
            if options.multicommunity>1   %store some suboptimal labels
                [temp_label_scores  temp_label_ID] =sort(act_less_exp,1);
                possible_multicom_labels=min([options.multicommunity size(temp_label_scores,1)])-1;
                secondary_labels_scores(end-possible_multicom_labels :end,:,subop_tracker)=      temp_label_scores(end-possible_multicom_labels :end,:);
                secondary_labels_ID    (end-possible_multicom_labels :end,:,subop_tracker)=active_labels(temp_label_ID(end-possible_multicom_labels :end,:));
                subop_tracker=subop_tracker+1;
            end
            
        else
            % to_merge_length_record=length(to_merge);
            cond1=find(min([label_number_at_time_for_certain_nodes_time0(:) label_number_at_time_for_certain_nodes_timeplus(:)],2)>=2); %idea is to only merge labels with larger number of members
            cond2=find(keyvals(:,1)>keyval_med); %we DO NOT process pairs of labels with >connectivity - find it helpful to do some updates instead of merging absolutely everything
            
            indices_to_check_for_potential_merging=intersect(cond1,cond2);  %will see if the labels at these indices have ben previously merged - only try to merge once per label
            
            for k=1:length(indices_to_check_for_potential_merging)
                if (been_merged(keyvals(indices_to_check_for_potential_merging(k),2:3)))==[ 0 0]
                    
                    merge_pairs(end+1,:)=keyvals(indices_to_check_for_potential_merging(k),2:3);
                    listener_history(i,label_idx_in_cells{active_labels(merge_pairs(end,1))})    =active_labels(merge_pairs(end,2));  %merge labels under the label of the first ID of the pir
                    
                    been_merged(keyvals(indices_to_check_for_potential_merging(k),[2,3]))=1;
                end
            end
            
            
            
            if sum(been_merged)==0  %this will be true when there is nothing permissible to merge
                
                switch_to_fusion_mode=0; %set this to zero exit this fusion mode
                postintervention_state(end+1,:)=listener_history(i,:);  %record the final state
                
                
                if options.multicommunity>1   %store some suboptimal labels
                    [ temp_label_scores  temp_label_ID] =sort(act_less_exp,1);
                    
                    possible_multicom_labels=min([options.multicommunity size(temp_label_scores,1)])-1;
                    secondary_labels_scores(end-possible_multicom_labels :end,:,subop_tracker)=      temp_label_scores(end-possible_multicom_labels :end,:);
                    secondary_labels_ID    (end-possible_multicom_labels :end,:,subop_tracker)=active_labels(temp_label_ID(end-possible_multicom_labels :end,:));
                    
                    subop_tracker=subop_tracker+1;
                end
            end
        end
        
    end  %conditions for fusion
    
    
    timer(i)=toc;
    
end %end main while loop
used_timesteps=i;


%trim off initial solutions
postintervention_state(1:options.discard_transient,:)=[];
max_labels_output(1:options.discard_transient)=[];
if options.multicommunity>1
    secondary_labels_scores(:,:,1:options.discard_transient)=[];
    secondary_labels_ID(:,:,1:options.discard_transient)=[];
end

partitionID=zeros(size(postintervention_state)); %for all clusters, get their locations in the partitions and their numeric identifier
if options.multicommunity>1
    for j=1:size(postintervention_state,1)
        [todel partitionID(j,:) todel2 secondary_labels_ID(:,:,j) ] =splitlist(postintervention_state(j,:),secondary_labels_ID(:,:,j));
    end
    
else
    for j=1:size(postintervention_state,1)
        [todel, partitionID(j,:)  ] =splitlist(postintervention_state(j,:));
    end
end




max_labels_output=mean(max_labels_output);
