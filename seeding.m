%setup initial labels for SE2. There can be repeated labels (IC_type==1) or
%one label per node (IC_type==2)
function IC_store=seeding(ADJ, IC_type, main_iter, ADJsum,options);

if IC_type==1
        reps_needed=ceil(length(ADJ)/options.target_clusters);
        repeating_labels=repmat(1:options.target_clusters,1,reps_needed);
     
        IC_store=zeros(options.independent_runs,length(ADJ));
        %for i=1:options.independent_runs*options.nback
        
        has_connection=find((ADJsum)-(diag(ADJ))'~=0);
%        has_connection=find(ADJsum>1);
%whos
%options.target_clusters
        for i=1:options.independent_runs
            IC_store(i,has_connection)=repeating_labels(randperm(length(has_connection)));
        %  length(unique(IC_store(i,:)))
            if length(has_connection)<length(ADJ)
            biggest_label=max(IC_store(i,:));
            no_connection=setdiff(1:size(IC_store,2),has_connection);
            IC_store(i,no_connection)=(1:length(no_connection))+biggest_label+1;
            

            end
          %  length(unique(IC_store(i,:)))
          %  whos
          %  pause
        end
        
        if main_iter==1
        disp(['produced about ' num2str(length(unique(IC_store(1,:)))) ' seed labels, while goal was ' num2str(options.target_clusters) ])      
        
        if options.verbose==1
           disp('finished IC gen') 
        end
        end
     
elseif IC_type==2
    IC_store=randperm(length(ADJ));
end

