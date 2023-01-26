% is list of nodeID's for nodes that show up in more than one community
%
%Inputs:
%"partitions" has columns which are alternate clustering results
%"partitionID" are alternate clustering results with clusters stored in cells, which make adjusted rand index faster to compute
%"accept_multi" higher values are more stringent criterion for overlapping clusters, 1==disjoint clusters
%secondary_labels_scores

function [nodes_and_partition_identifiers_disjoint nodes_and_partition_identifiers_overlapping cell_partition cell_partition_overlapping multicom_nodes_all median_of_all_ARI]=select_representative_partition_old(ADJ,partitions, main_iter, secondary_labels_scores, secondary_labels_ID, subset_nodes_for_NMI,options)

%%first task is to compare all partitions and find representative one
adjustedrand_pairwise=zeros(size(partitions,2));  %holds all possible adjusted rand index comparisons of partitions
adjustedrand_pairwise_subset=zeros(size(partitions,2));  %holds all possible adjusted rand index comparisons of partitions

nmi_pairwise=zeros(size(partitions,2));

core_info = evalc('feature(''numcores'')');
usable_cores=max([1 str2num(core_info(end-4:end-1))]);



if length(ADJ)>subset_nodes_for_NMI
    partitions_store=partitions;  %need full info later
%    partitions=partitions(randsample(length(ADJ),subset_nodes_for_NMI,1),:);   %     original version needs stats package
partition_picker_idx=randperm(length(ADJ));
partitions=partitions_store(partition_picker_idx(1:subset_nodes_for_NMI),:);

end


for i=1:size(partitions,2)  %calculate similarity of partitions
    
    
    if size(partitions,2)>100
        if mod(i,ceil(size(partitions,2)/10))==0
            disp(['ARI computation ' num2str(round(100*i/size(partitions,2))) ' % complete'])
        end
    end

    parfor (j=1:size(partitions,2),options.max_threads)  %logic is this is pretty fast, rarely have more than 100 partitions, and setting up a new cluster with a bunch of workers (usable_cores variable) rarely worth this initiation time, but if you already have one up via max_threads, fine to use that      
  
        if j<i %metric is symmetric, so save time by computing only half
  %          nmi_pairwise(i,j)=similarity(partitions(:,i),partitions(:,j),'type','nmi')  ;     
            nmi_pairwise(i,j)=discrete_nmi(partitions(:,i),partitions(:,j))  ;     

        end
    end
    
    
    
end
nmi_pairwise=nmi_pairwise+nmi_pairwise';  %axis==0 since we're going to take max

median_of_all_ARI=median(nmi_pairwise(:));
mean_of_all_ARI=mean(nmi_pairwise(:));
if main_iter==1
    disp(['median of all NMIs ' num2str(median_of_all_ARI)])
    disp(['mean of all NMIs ' num2str(mean_of_all_ARI)])
end


[most_similar_var most_similar_idx]=max(sum(nmi_pairwise));


if length(ADJ)>subset_nodes_for_NMI

partitions=partitions_store;
end

winning_partition=partitions(:,most_similar_idx);
winning_partition_memberIDs_unq=unique(partitions(:,most_similar_idx));




%%strategy below is to use info from all partitions(not just winning one) to
%find nodes in multiple communities.  However, community defs are set by
%the winning partition.  To reiterate, the selection of which nodes are
%multicommunity is separate from the question of what communities of which they are members.

if options.multicommunity>1
    %dimensions of secondary_labels_scores are [options.multicommunity x #nodesInADJ x options.target_partitions]
    
    if options.graphics==1
        figure
        imagesc (secondary_labels_scores(:,:,1))
        xlabel('nodes')
        title('some randomly selected score')
 
        figure
        imagesc (sum(secondary_labels_scores,3))
        title('raw overlaps, so influenced by #partitions')
    end
    
    secondary_labels_scores=sum(secondary_labels_scores,3);  %the final row of this has the scores for the top (discrete) communities (not overlapping)
    secondary_labels_scores=secondary_labels_scores./ repmat( secondary_labels_scores(end,:),size(secondary_labels_scores,1),1 );  %we compare suboptimal scores normalized to the top score for each node
    
    if options.graphics==1
        figure
        plot(sort( secondary_labels_scores(end-1,:) ))
        ylabel('sorted ratio of 2nd best to best')
    end
    
    secondary_labels_scores=    secondary_labels_scores(end-options.multicommunity+1:end-1,:);  %trim  off last row which holds the discrete solutions
    secondary_labels_ID_winning=secondary_labels_ID    (end-options.multicommunity+1:end-1,:,most_similar_idx);
    
    if options.graphics==1
        figure
        imagesc(secondary_labels_scores)
        title('secondary_labels_scores')
    end
    
    cutoff_for_multicom_status=2/(options.multicommunity+1);  %justified by concept that labels should be split roughly n-ways
    [row_not_used multicom_nodes_all]=(find(secondary_labels_scores>=cutoff_for_multicom_status));
    multicom_nodes_all=multicom_nodes_all(:);  %multicom_nodes_all are raw indices
    multicom_lin_idx=                 (find(secondary_labels_scores>=cutoff_for_multicom_status));
    multicom_nodes_labels=secondary_labels_ID_winning(multicom_lin_idx);  %"multicom_nodes_labels" has the labels for each node in "multicom_nodes_all"
    multicom_nodes_labels=multicom_nodes_labels(:); %force column because this can switch orientation if "secondary_labels_ID_winning" is  single row whichi then cuases problems
    
    close all
    
    cell_partition={}; %each cell in this will contain the indices of nodes in a given cluster
    for i=1:length(winning_partition_memberIDs_unq)
        
        cell_partition_overlapping{i,1}=        cat(1,find(winning_partition==winning_partition_memberIDs_unq(i)), multicom_nodes_all(find(multicom_nodes_labels==winning_partition_memberIDs_unq(i))) );  %seems like same content as partitionID(most_similar_idx) but reordered       
        cell_partition{i,1}=find(winning_partition==winning_partition_memberIDs_unq(i));
    end
    
    
else  %for discrete clusters just load indices of each label into cells - they will be ordered for visualization later
    cell_partition={}; %get the indices of nodes which end up in the same cluster
    
  %  save selectstuff.mat
    winning_partition_sorted=sortrows([[1:length(winning_partition)]' winning_partition],2);
    cell_partition = accumarray(winning_partition_sorted(:,2),winning_partition_sorted(:,1),[],@(v){v}); % split data according to matched bins.
%     for i=1:length(winning_partition_memberIDs_unq)
%         cell_partition{i,1}=find(winning_partition==winning_partition_memberIDs_unq(i));  %seems like same content as partitionID(most_similar_idx) but reordered
%     end
%     
%     selectC{1}
%     cell_partition{1}
%     selectC{2}
%     cell_partition{2}
% 
%     selectC{3}
%         cell_partition{3}

    
    cell_partition_overlapping=cell_partition;  %we still produce this output for convenience in the discrete case but it's NOT overlapping
    multicom_nodes_all=[];
end

if options.multicommunity>1
    if main_iter==1
        disp(['overlapping vs discrete length (at lev 1): ' num2str(sum(cellfun(@length,cell_partition_overlapping))) ' vs ' num2str(sum(cellfun(@length,cell_partition)))])
    end
end

%%from here on it's just arranging the output
[trash, idx_large_partition] = sort(cellfun('size', cell_partition, 1), 'descend');
cell_partition=cell_partition(idx_large_partition); % reorder partitions by size, with no change to contents

[trash, idx_large_partition] = sort(cellfun('size', cell_partition_overlapping, 1), 'descend');
cell_partition_overlapping=cell_partition_overlapping(idx_large_partition); % reorder partitions by size,

cluster_density=cell(length(cell_partition),1); %sort order of nodes within each cluster for display purposes
for i=1:length(cell_partition)
    cluster_density{i}=mean(ADJ(cell_partition{i},cell_partition{i}));
    [tmp,ind]=sort(cluster_density{i}, 'descend');
    cell_partition{i}=cell_partition{i}(ind);
end


cluster_density_overlapping=cell(length(cell_partition),1); %sort order of nodes within each cluster for display purposes
if length(ADJ)>100000 && options.max_threads~=0  %not worth it to bootup cluster on small data
    parfor ( i=1:length(cell_partition_overlapping),options.max_threads)
        cluster_density_overlapping{i}=mean(ADJ(cell_partition_overlapping{i},cell_partition_overlapping{i}));
        [tmp,ind]=sort(cluster_density_overlapping{i}, 'descend');
        cell_partition_overlapping{i}=cell_partition_overlapping{i}(ind);
    end
    
else
    for i=1:length(cell_partition_overlapping)
        cluster_density_overlapping{i}=mean(ADJ(cell_partition_overlapping{i},cell_partition_overlapping{i}));
        [tmp,ind]=sort(cluster_density_overlapping{i}, 'descend');
        cell_partition_overlapping{i}=cell_partition_overlapping{i}(ind);
    end
end


partition_marker_sorted_hard=[];
partition_marker_sorted_overlapping=[];
for i=1:length(cell_partition)
    partition_marker_sorted_hard(end+1:end+length(cell_partition{i}))=i;
    partition_marker_sorted_overlapping(end+1:end+length(cell_partition_overlapping{i}))=i;
end
nodes_and_partition_identifiers_disjoint=sortrows([vertcat(cell_partition{:}) partition_marker_sorted_hard']);

nodes_and_partition_identifiers_overlapping=sortrows([vertcat(cell_partition_overlapping{:}) partition_marker_sorted_overlapping']);
