%Goal is to estimate the reliability of each nodes placement in a given
%cluster by using the replicate partitions.
%
%Inputs: expects partitions to be input as rows
%"winning partition is an actual partition, NOT just a ref to where it is in "partitions"
%"winning partition" has two columns (nodeID and label)
%
%Output: a column of cells, one for each cluster, range [0 1] ,with 1 being max stability/confidence
function node_confidence_cell=node_confidence(partitions, winning_parition)

node_confidence=zeros(size(winning_parition)); %will hold same values as "node_confidence_cell'
for i=2:size(partitions,1)%ensure labels in each partition do not repeat
    
    partitions(i,:)=partitions(i,:)+max(partitions(i-1,:));
    
end


winning_parition_unq=unique(winning_parition(:,2));
for i=1:length(winning_parition_unq) %for each cluster
    
    current_node_idx=winning_parition(find(winning_parition(:,2)==winning_parition_unq(i)),1);
       
    labels_in_cluster_all_partitions=partitions(:,current_node_idx);  %all labels of nodes in a given cluster, from all replicate partitions
    
    xcoord=repmat([1:size(labels_in_cluster_all_partitions,2)],size(labels_in_cluster_all_partitions,1),1);
    %i.e. if there are four nodes in a cluster and we do 3 replicate paritions "xcoord" is:
    %[1 2 3 4]
    %[1 2 3 4]
    %[1 2 3 4]
    sorted_labels=sortrows([labels_in_cluster_all_partitions(:) xcoord(:)]);
    endpoints=[find(sorted_labels(1:end-1,1)-sorted_labels(2:end,1)~=0) ; length(sorted_labels)];
    startpoints=[1; endpoints(1:end-1)+1];
    counts_of_labels_in_cluster_all_partitions=endpoints-startpoints+1; 
    
    %the measure we're using for reliablity is how often a given node in a
    %given cluster had the same label as most other nodes in that cluster
    if length(current_node_idx)>1
        cooccurence_per_unique_label=(counts_of_labels_in_cluster_all_partitions-1) ./(repmat(size(labels_in_cluster_all_partitions,2), length(counts_of_labels_in_cluster_all_partitions), 1)-1); %-1 due to node always being it's own label

    else
        
        cooccurence_per_unique_label=zeros(size(counts_of_labels_in_cluster_all_partitions)); %for node in it's own community - you get nans from divid zero in cooccurence_per_unique_label otherwise
    end
    
    hold_repeated_cooc_vals=zeros(size(sorted_labels,1),1);
    for j=1:length(counts_of_labels_in_cluster_all_partitions)
        hold_repeated_cooc_vals(startpoints(j):endpoints(j)) =cooccurence_per_unique_label(j);
        
    end
    
    
    raw_node_confidence=[accumarray(xcoord(:),hold_repeated_cooc_vals,[],@sum)/size(partitions,1)];
    node_confidence(current_node_idx)=raw_node_confidence; %put back in original order
    node_confidence_cell{i,1}=raw_node_confidence;
end


