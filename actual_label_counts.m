%returns the number of each label that each cell "heard" over the past
%several time-steps (specified by options.nback)
%
%another way of thinking about this is we're summing select rows of the
%ADJ - all those in a particular label
%
%return matrix format is labels as rows, nodes as columns, and values as
%the number of labels input to a given node from a given label.  The labels
%of each row (the actual numerals) are provided in "labels_unq"
%
%to count up how many labels of each type that each node "hears",
%add up rows of ADJthat have same label identifier


function [labelID_x_nodes labels_unq comm_size]=actual_label_counts(ADJ,current_listener_history,options,varargin)


%section is useful to understand what is going on here from a simpler and
%much slower perspective
% diagnose=0;
% if diagnose==1
%  simple_version=zeros(length(unique(current_listener_history)),length(ADJ));
%
%  labels_unq=unique(current_listener_history);
%
%  for i=1:length(ADJ)
%
%      label_idx=find(labels_unq==current_listener_history(i))
%
%          simple_version(label_idx,i)=1;
%
%  end
%  figure
%  imagesc(simple_version)
%  title('simple version')
%  ylabel('unique labels')
%  pause
% end

if length(varargin)==1
    
    ignore_nodes=varargin{1};
    new_largest_label=max(current_listener_history)+1;
    current_listener_history(ignore_nodes)=new_largest_label;
    
end

current_listener_history=current_listener_history(:);
[discard indices]=sort(current_listener_history);%sort faster than sortrows
temp=[current_listener_history (1:length(current_listener_history))'];

sorted_labels=temp(indices,:);
transition_spots=find ( [sorted_labels(1:end-1,1)-sorted_labels(2:end,1) ] ~=0);
transitions=[1; 1+transition_spots]; %first element of subsequent set of labels
labels_unq=sorted_labels(transitions); %faster than finding labels in full set, since we know where new labels will be

comm_size=histc(current_listener_history,labels_unq);

future_markers=zeros(size(sorted_labels,1),1);
future_markers(transitions)=1;
future_markers=cumsum(future_markers);
node_identifiers=zeros(size(sorted_labels,1),1);
node_identifiers(sorted_labels(:,2))=future_markers;


%idea is to consider the labels from different time-steps to be different
%(even though they are not) then add the row of temp created with these pseudo-different labels, to add up the rows that are infact related to the %same core label
%the obvious way to do this is:
%temp=zeros(length(labels_unq), length(node_identifiers));
%temp=bsxfun(@eq, sparse(1:length(labels_unq)).', node_identifiers');
%but a faster way is:
running_sum=sparse(node_identifiers, 1:length(node_identifiers), ones(length(node_identifiers),1), length(labels_unq), length(node_identifiers));


if exist('ignore_nodes')
   running_sum(end,:)=[];
   labels_unq(end)=[];
   comm_size(end)=[];
   
end
% %we can do that becaue node identifiers are numbers sequentially starting at 1
% %in each row of "nodes_by_labels_all_times" we tick off positions (using a 1) where that label occurs in the full list of labels

    if options.memory_efficient==1       %weird this is working faster even on full matrices
  
        labelID_x_nodes=running_sum*ADJ;     
    else%
        labelID_x_nodes=full(running_sum)*ADJ;
    end

