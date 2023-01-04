%This function is called by SpeakEasy2.m to generate partitions for consensus
%clustering.  First gets initial conditions for a set of runs (seeding.m)
%then calls the core label updating function(SpeakEasycore2_core.m) a number of times
%finally passes results of all runs to "select_representative_partitions.m"
%to pick consensus solution

function [partition_codes partition_codes_overlapping cell_partition cell_partition_overlapping  multicom_nodes_all  median_of_all_ARI est_node_confidence]=SpeakEasy2_bootstrap(ADJ, main_iter, options) %#ok<NCOMMA>

kin=full(sum(ADJ));
kout=full(sum(ADJ,2));
ktot=(kin+kout');
if main_iter==1
    if options.multicommunity>1
        disp('attempting overlapping clustering')
    end
    IC_store=seeding(ADJ,1,main_iter, kin, options);  %repeated label
    disp('completed generating initial labels')
else
    IC_store=seeding(ADJ,1,main_iter, kin, options);  %one label per node... coould switch this to repeated labels, but less needed with smaller clusters
end


if main_iter==1 & options.max_threads>0
    disp('starting level 1 clustering; independent runs might not be displayed in order - that is okay')
end


parfor(i=1:options.independent_runs,options.max_threads)  %this is the main loop, feeding in initial conditions to start each independent run

    if options.seed_set_by_user==1
    rng(i,'twister')
    end
    
   % if options.subcluster==1
    if main_iter==1
        disp(' ')
        disp(['starting independent run #' num2str(i) ' of ' num2str(options.independent_runs)])
    end
    
    [partitionID secondary_labels_scores secondary_labels_ID max_labels_output ]=SpeakEasy2_core(ADJ,IC_store(i,:),main_iter,kin,ktot,options);
    partitionID_store{i,1}=partitionID; %parfor isn't happy unless this is filled seprately like this as opposed to direct funtion output
    secondary_labels_scores_in_cells{i,1}=secondary_labels_scores;
    secondary_labels_ID_in_cells{i,1}=secondary_labels_ID;
end

secondary_labels_ID=cat(3,secondary_labels_ID_in_cells{:});
secondary_labels_scores=cat(3,secondary_labels_scores_in_cells{:});
partitionID_store=cell2mat(partitionID_store);



if options.multicommunity~=1 | length(ADJ)<20001
    
    subset_nodes_for_NMI=length(ADJ);
    
else
    subset_nodes_for_NMI=10000;
    
end
[partition_codes partition_codes_overlapping cell_partition cell_partition_overlapping multicom_nodes_all median_of_all_ARI]=select_representative_partition(ADJ, partitionID_store', main_iter, secondary_labels_scores, secondary_labels_ID, subset_nodes_for_NMI, options);


if options.node_confidence==1
est_node_confidence=node_confidence(partitionID_store,partition_codes_overlapping);
else 
est_node_confidence=[];
end



if main_iter==1
    disp(['generated ' num2str(size(partitionID_store,1)) ' partitions at level ' num2str(main_iter)])
end

