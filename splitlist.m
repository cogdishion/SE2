%returns the sets of indices for each unique value in a matrix (first output).
%also can take a list of integers and relabels them so they start with 1 and are
%sequential (2nd output). 
%also you might want to realign labels in secondary_labels_ID in sync with
%'original_list' which is what varargin handles

function [list_in_cells list_renumbered cell_contents_value varargout]=splitlist(original_list, varargin)


if min(size(original_list))~=1
    error('only works when first input is a row or column')
end

original_list=(original_list(:))';

if length(varargin)==1
    another_matrix_to_adjust=varargin{1};
    
    [row_w_no_values col_w_no_values]=find(another_matrix_to_adjust==0);  %because sometimes we coulnd't generate an n-best secondary solution, so there are unused rows in 'another_matrix_to_adjust' from SpeakEasy2_core.m filled with zeros
    another_matrix_to_adjust(1:max(row_w_no_values),:)=[];
    
    
    
    another_matrix_to_adjust_unique=unique(another_matrix_to_adjust(:));
    if length(find(ismember(another_matrix_to_adjust_unique,original_list)==1))~=length(another_matrix_to_adjust_unique)
        error('serious problem in multicom where you have a secondary label for a node not among the primary labels')
    end
    
    original_and_another=cat(1, original_list,another_matrix_to_adjust);
    
    
    
    [source_labels label_indices]=sort(original_and_another(:)');
    shift_in_label_list=find([source_labels(2:end)-source_labels(1:end-1)  ]~=0);
    shift_in_label_list=[0 shift_in_label_list length(source_labels)];
    
    combined_list_renumbered=zeros(size(another_matrix_to_adjust,1)+1,size(another_matrix_to_adjust,2));
    for sink_labels=1:length(shift_in_label_list)-1
        %list_in_cells{sink_labels}= label_indices(shift_in_label_list(sink_labels)+1:shift_in_label_list(sink_labels+1))  ;
        combined_list_renumbered(label_indices(shift_in_label_list(sink_labels)+1:shift_in_label_list(sink_labels+1)))=sink_labels;
    end
    
    combined_list_renumbered=reshape(combined_list_renumbered,size(another_matrix_to_adjust,1)+1,size(another_matrix_to_adjust,2));
    another_matrix_to_adjust_renumbered=combined_list_renumbered(2:end,:);
    varargout{1}=cat(1, zeros(max(row_w_no_values),size( another_matrix_to_adjust_renumbered ,2) )  ,another_matrix_to_adjust_renumbered);  %stick any original zeros back on
end


[source_labels label_indices]=sort(original_list);
shift_in_label_list=find([source_labels(2:end)-source_labels(1:end-1)  ]~=0);
shift_in_label_list=[0 shift_in_label_list length(source_labels)];

list_renumbered=zeros(size(original_list));
list_in_cells=cell(1,length(shift_in_label_list)-1);
length(shift_in_label_list);
for sink_labels=1:length(shift_in_label_list)-1

    
    list_in_cells{sink_labels}=label_indices(shift_in_label_list(sink_labels)+1:shift_in_label_list(sink_labels+1))  ;
    list_renumbered(list_in_cells{sink_labels})=sink_labels;

end

cell_contents_value=1:length(shift_in_label_list)-1;
