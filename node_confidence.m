%partitions shold be as rows
%winning parititinos is an actual partition, NOT just a ref to where it is in "partitions"
%
%output is a column, length(ADJ) [0 1] with 1 being max stability/confidence
function node_confidence_cell=node_confidence(partitions, winning_parition)

node_confidence=zeros(size(winning_parition));
%imagesc(partitions)
%length(find(partitions(1:end-1,:)-partitions(2:end,:))~=0)
%pause
for i=2:size(partitions,1)%ensure labels in each partition do not overlap
    
    partitions(i,:)=partitions(i,:)+max(partitions(i-1,:));
    
end


winning_parition_unq=unique(winning_parition(:,2));

for i=1:length(winning_parition_unq)
    
    current_node_idx=winning_parition(find(winning_parition(:,2)==winning_parition_unq(i)),1);
    
    
    labels_in_cluster_all_partitions=partitions(:,current_node_idx);
    
    xcoord=repmat([1:size(labels_in_cluster_all_partitions,2)],size(labels_in_cluster_all_partitions,1),1);
    
    sorted_labels=sortrows([labels_in_cluster_all_partitions(:) xcoord(:)]);
    
    
    max(current_node_idx);
    endpoints=[find(sorted_labels(1:end-1,1)-sorted_labels(2:end,1)~=0) ; length(sorted_labels)];
    startpoints=[1; endpoints(1:end-1)+1];
    labels_themselves=sorted_labels(endpoints,1);
    counts_of_labels=endpoints-startpoints+1;
    
    
    if length(current_node_idx)>1
        cooc_val_per_unique_label=(counts_of_labels-1) ./(repmat(size(labels_in_cluster_all_partitions,2), length(counts_of_labels), 1)-1);
    else
        
        cooc_val_per_unique_label=zeros(size(counts_of_labels)); %for node in it's own community - you ge nans from divid zero in cooc_val_per_unique_label otherwise
    end
    %[cooc_val_per_unique_label labels_themselves]
    
    hold_repeated_cooc_vals=zeros(size(sorted_labels,1),1);
    for j=1:length(counts_of_labels)
        hold_repeated_cooc_vals(startpoints(j):endpoints(j)) =cooc_val_per_unique_label(j);
        
    end
    
    % [hold_repeated_cooc_vals xcoord(:)]
    
    raw_node_confidence=[accumarray(xcoord(:),hold_repeated_cooc_vals,[],@sum)/size(partitions,1)];
    node_confidence(current_node_idx)=raw_node_confidence;
    node_confidence_cell{i,1}=raw_node_confidence;
end

