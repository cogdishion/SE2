%By checking the network density and weighted/unweighted format, we can
%potentially save time on later operations.
%For larger matrices we only sample their properties as it would be slow to
%check every link.

%Outputs:
%is_ADJ_weighted==1 for case when entries are not entirely 0 or 1
%force_efficient==1 if matrix is larger than max_ADJ_size or if it is smaller than max_ADJ_size and low density
function [ADJ is_ADJ_weighted is_ADJ_symmetric]=ADJcharacteristics(ADJ, ADJsum)

smallest_val=min(min(ADJ));
max_all_ADJ=max(max(ADJ));
 if length(find(isnan(ADJ))~=0)
        error('you have NaNs in your data, please fix')
    end
    
    
    if size(ADJ,1)~= size(ADJ,2)
        error('your ADJ is not square, please fix')        
    end
    
    if max_all_ADJ>1 || min(min(ADJ))<-1
        error('your connection strength is outside [-1 1] please fix or normalize to that range')
    end
    
     noinputs=find(sum(abs(ADJ))==0);  %because if a node has no inputs it can't receive labels, we set it to receive it's own label, i.e. remain isolated
    if ~isempty(noinputs)
        ADJ((noinputs-1)*length(ADJ)+noinputs)=1;
        disp('note, you had a few nodes with no inputs, to which we assigned self inputs')
    end
   

fraction_to_be_full=.1;  %matrices less dense than this are cohnverted to sparse

if size(ADJ)<1000  %to conserve memory, test a fraction of links unless ADJ is small
    sample_of_links=ADJ(:);
    
else
    sample_of_links=ADJ(ceil(numel(ADJ)*rand(1,100000)));  %get random ADJ values
end

ADJ_density=length(find(sample_of_links~=0))/length(sample_of_links);
disp(['approximate edge density is ' num2str(ADJ_density)])

if length(find(sample_of_links~=0))/length(sample_of_links)>fraction_to_be_full  %dense ADJ
    
    if length(find(sample_of_links>.9999))+length(find(abs(sample_of_links)<.0001))~=length(sample_of_links);  %in case there are rounding errors
        disp('input type treated as weighted full')   %ADJ with links of only +1 and -1 will still be considered weighted
        is_ADJ_weighted=1;
        
    else
        disp('input type treated as unweighted full')
        is_ADJ_weighted=0;
        
        
    end
    
else    %non-dense ADJ
    if ~issparse(ADJ)
        ADJ=sparse(ADJ);
    end
    
    %in this case we can afford to test all links
    [todel1 todel1 edge_val]=find(ADJ);
    if length(find(edge_val>.9999))~=length(edge_val);  %in case there are numeric issues
        disp('weighted low density')
        is_ADJ_weighted=1;
        
    else
        disp('unweighted low density')
        is_ADJ_weighted=0;
    end
end



if length(ADJ)<1000
   if max(max(ADJ-ADJ'))<.00001;
        disp('ADJ is symmetric')
        is_ADJ_symmetric=1;

   else
       disp('ADJ not symmetric and considered directed')
       is_ADJ_symmetric=0;
   end
   
   
else  %on a larger matrix, samples some points.... one worry is that on ultra sparse matrices,this will say it is symmetric by just finding a bunch of zeros
   
   % randrow=ceil(length(ADJ)*rand(1,5000));
   % randcol=ceil(length(ADJ)*rand(1,5000));
   % ADJ_subset=ADJ((randrow-1)*length(ADJ)+randcol)-ADJ((randcol-1)*length(ADJ)+randrow);
   % if mean(abs(ADJ_subset))<.0001

  randrowcol=ceil(rand(1,1000)*length(ADJ));
    ADJ_subset=ADJ(randrowcol,randrowcol);

   if max(max(ADJ_subset-ADJ_subset'))<.00001;
   
   disp('ADJ is probably symmetric')
        is_ADJ_symmetric=1;
    else
        disp('ADJ not symmetric and considered directed - be sure you know what this entails')
        is_ADJ_symmetric=0;
        
        
    end
end

 
    
if is_ADJ_weighted==1

    ADJ=ADJ./max([abs(min(min(ADJ))) max_all_ADJ]);  %make largest abs val ==1
   
    if length(ADJ)<1000

    sample_edges=ADJ(find(ADJ));
    else
    sample_edges=ADJ_subset(find(ADJ_subset));
    end

    if skewns(sample_edges)>=2  
           disp('high skew to edge weight distribution; reweighting main diag')
   
    mean_link_weight_by_node=(sum(ADJ)-diag(ADJ)')./(sum(sign(ADJ))-sign(diag(ADJ)'));
    ADJ(1:length(ADJ)+1:numel(ADJ))=mean_link_weight_by_node;
        
       if smallest_val>=0
               disp('adding very small offset to all edges')

    offset=mean(mean_link_weight_by_node);

            if issparse(ADJ)
                all_edges=(find(ADJ));
                ADJ(all_edges)=ADJ(all_edges)+offset;

            else
                    ADJ=(1-offset)*ADJ+offset*sign(ADJ);

            end
        end    
            
    else
            ADJ(1:length(ADJ)+1:numel(ADJ))=1;

    end

  
else %i.e. unweighted
        ADJ(1:length(ADJ)+1:numel(ADJ))=1;

    
end


end



function output= skewns(x)

output=(sum((x-mean(x)).^3)./length(x)) ./ (var(x,1).^1.5);
end

% 
%     if length(ADJ)<max_ADJ_size        %if matrix is small enough that we could use full multiplicationn
%         
%         if ADJ_density>.03             %do not change this value; this value was selected by observing the two multiplication strategies it switches between require equal time at .03 density
%             
%             memory_efficient=0;         %ADJ dense enough that full mult is faster than pairwise (and small enough it will fit in memory)
%         else
%             if ~issparse(ADJ)
%             ADJ=sparse(ADJ);
%             end
%             memory_efficient=1;    %'ultrasparse and small ADJ, so pairwise multi is faster'
%         end
%         
%     else
%         if ADJ_density<=.03  %you still might want ADJ sparse if it's huge and RAM usage is an issue, but hoping the user handles that
%             memory_efficient=1;
%             
%         else
%             memory_efficient=0; %if you have a huge matrix is still may be wise to put it as sparse even if density >.03.  .03 is the point at which it becomes faster to compute with a full, but you may need to convert to sparse simply for memory reasons
%             
%         end
%     end

  %memory_efficient=1;

