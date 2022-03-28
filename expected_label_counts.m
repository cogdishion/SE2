%input format for "actual_counts" is labels as rows and nodes as columns,
%and values as the number of a given label received by a given node
%
%returns the number of labels each node should expect to hear at random for
%each node in an ADJ, given the current set of labels
%


function expected=expected_label_counts(actual_counts,kin,do_full_version)

counts_normkout=sum(actual_counts,2);  %',2)' needed for case of size-1 clusters %counts_normkout is the total number of transmitted labels (for each label); sum(actual_counts) will be sum(ADJ); thus we're counting labels, which are influenced by k -  hence the name
counts_normkout_norm1=full(counts_normkout./sum(counts_normkout));  %proportions of various labels normalized to 1
%counts_normkout_norm1 incorporates not just proportions of labels, but also k, because we care about how frequently a label is transmitted (and hence received or expected to be received)

%default
if do_full_version==0  %if matrix is very sparse, or too large to store a full ADJ, we only care about generating expected counts of labels if a node actually receives some of that label
 
    [x y]=find(actual_counts);   %x refers to a lable and y refers to a node
    
    
   
    
    % of note for future debugging:
    %     x  %x for some reason turns into a row when it's all 1's, which then turns [full(counts_normkout_norm1(x(:)))] into a row, when it's usually a column which then screws up the dot product
    % [full(counts_normkout_norm1(x(:)))]
    % [kin(y)']
    
  %  if options.nback==1  %save multiplying by 1
        expected=sparse(x,y,(counts_normkout_norm1(x(:)))'.* kin(y),              size(actual_counts,1),size(actual_counts,2));  %same format as lxn new
        % counts_normkout_norm1(x(:)))'   and     kin(y) are two long columns
        
       % the_full_version=([full(counts_normkout_norm1)]*full(kin));  % this is the full version, even though currently embedded in the sparse routine
%        sum(sum(expected))
%sum(sum(the_full_version))
%'exp'
%pause
%      
%     else
%         expected=sparse(x,y,options.nback*((counts_normkout_norm1(x(:)))'.* kin(y)),              size(actual_counts,1),size(actual_counts,2));  %same format as lxn new
%     end
    
    %the two lines below are easier to understand and perform the same operation as the above, but are slightly slower
    %scaled_kin=nback*([full(counts_normkout_norm1(x))]'.* kin(y)); %scales normalized counts by total input (some nodes have more inputs and thus you would expect more of all labels)
    %expected=sparse(x,y,scaled_kin, size(actual_counts,1),size(actual_counts,2));  %same format as lxn
    
    
    
    
else  %not so memory efficient, but can be faster
 
     %   if options.nback==1  %save multiplying by 1
           %     [x y]=find(actual_counts);   %x refers to a label and y refers to a node
            %        expectedtest=sparse(x,y,(counts_normkout_norm1(x(:)))'.* kin(y),              size(actual_counts,1),size(actual_counts,2));  %same format as lxn new


          expected=([full(counts_normkout_norm1)]*full(kin));  % this is the full version, even though currently embedded in the sparse routine
% 
%         else
%           expected=options.nback*([full(counts_normkout_norm1)]*full(kin));
% 
%      
%         end
        

        
end
        
            
%             'actual       expected'
%             [ sum(sum(actual_counts)) sum(sum(expected))]  %the reason this version of sum(sum(expected)) isn't equal to sum(sum(actual)) is that we only calculate some expected values - those where we actually have an incoming "actual" value - it actually is equal to other version of 'expected' caluclated by this function at those points
%             'expectedv2'
%             sum(sum(expectedv2))
%             'in expect routine'
%             
%             
%             sum(expected(find(actual_counts==0)))
%             'see no vals'
%             
%             
%             sum(expectedv2(find(actual_counts==0)))+sum(sum(expected))
%             'quote missing negatives observed here'
%             pause

    
    
    
   





























