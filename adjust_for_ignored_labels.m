%This function is needed because we don't update all labels at every 
%time step (to improve convergence), we might not update members of
%some labels at all, which will lead to "active_labels" having more
%unique labels.
%
%That is just due to nature of selective/random updating, and we need to account for it.
%
%concretely, this deletes some elements (corresponding to select labels) of "active_labels" and some rows
%(corresponding to labels) of "actual_minus_expected"

function  [active_labels actual_minus_expected]=adjust_for_ignored_labels(listener_history,active_labels,actual_minus_expected)

sampled_labels= unique(listener_history);
if length(sampled_labels)~=length(active_labels)
    
    % length(sampled_labels)  %this will be smaller than length(active_labels) because it randomly didn't sample various labels
    excluded_label=setdiff(active_labels,sampled_labels);   

    if ~isempty(excluded_label)
        not_used_idx=find(ismember(active_labels,excluded_label)==1)   ;
        active_labels(not_used_idx)=[];
        actual_minus_expected(not_used_idx,:)=[];  
    end
end
