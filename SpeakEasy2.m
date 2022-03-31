%%This is the suggested fuction through which to interact with SpeakEasy2
%clustering.  You can run SpeakEasy without setting any parameters, or
%control exactly how it functions, using the parameters below, which
%control how the algorithm is applied, but not how it operates.
%
%Description of outputs:
% varargout{1} - "partition_tags{x}" is a two column matrix with node ID's in the first column and numeric label ID's in the 2nd, {x} refers to subcuster level
% varargout{2} - "partition_cells{x}" is a sequence of cells, each containing a list of nodes in a given cluster.  The clusters are ordered in size.
% varargout{3} - "convenient_order" lists nodes in conveneint order for visualizing clusters  - i.e.as a quick check that clusters are plausible, try running: imagesc(ADJ(convenient_node_ordering{1},convenient_node_ordering{1})
% varargout{4} - list of multi-community nodes

%Example uses
% >> SpeakEasy2(rand(30)); % just to verify it runs
% >> load ADJdemo
% >> [partition  partition_cell_version  convenient_node_order]=SpeakEasy2(ADJ);  %output results to workspace
% >> SpeakEasy2(ADJ,'subcluster',2)  %clustering primary clusters again
% >> SpeakEasy2(ADJ,'maxthreads',3)  %three independent threads - requires parallel toolbox

function [varargout]=SpeakEasy2(ADJ,varargin); %#ok<NCOMMA>

%%Settings.  In this section change how SE2 is applied, but not its fundamental behavior.
%The most likely inputs to be adjusted are at the top.
options = inputParser;
options.CaseSensitive = false;

%settings related to how many times SpeekEasy is applied
addOptional(options,'filename','SpeakEasy2_results')
addOptional(options,'independent_runs',10);   %number of completely independent initil conditions (IC's)
addOptional(options,'subcluster',1);           %if you want to sub-cluster your primary clusters, make this 2+
addOptional(options,'multicommunity',1)  %rename after testing... should be equal to max number of communities per node (so will be 2 or greater if you want overlapping output
addOptional(options,'target_partitions',5);
addOptional(options,'target_clusters',max([min([10 length(ADJ)]) round(length(ADJ)/100)]));  %if you enable bubbling and run long enough this doesn't matter

%core and probably don't need to ever change, as these are generally arbitrary, fast and accurate
addOptional(options,'minclust',5);         %min size for sub-clustering (i.e. if cluster is already this small, do NOT further subcluster) just matke this more than 3 as suits your purose, but
addOptional(options,'timesteps',100);     %no longer matters as it did in original SpeakEasy, as total time is adaptive now
addOptional(options,'subcluster_split',2)  %when bubbling, split clusters into this many parts - can be several more, but doesn't really matter - probably remove
addOptional(options,'contiguous_labels',0)  % set to 1 forces all communities to be connected - sounds like something you want, but generall not really worth it/needed
addOptional(options,'discard_transient',3)  %disregard the first few solutions that are not at equilibrium

%these merely affect how we do the clustering, but not the result
addOptional(options,'max_threads',0);        %you need the parallel package for this - will create copies of the ADJ and run in parallel, so be sure you have enough memory
addOptional(options,'memory_efficient',1); %setting to zero may improve sped on full matrices that fit in memory
addOptional(options,'random_seed',rng(randi(10000,1)));   %for repro
addOptional(options,'autoshutdown',1); %shutsdown parpool unless you set to 0 - which yiou might want to do if doing a bunch of runs on mulitple networks

%for extra output
addOptional(options,'verbose',0);
addOptional(options,'graphics',0);
addOptional(options,'node_confidence',0);

%addOptional(options,'do_memmap',0);   %only works for non-sparse inputs  -
%would have to rewrite ADJ variable because can't pass memmap obj by reference - clunky
%% Do some error checks on inputs and parameters
%parse(options,varargin{:});

% if options.Results.do_memmap==1
%     ADJ_mapped=do_memmap(ADJ);
%     ADJ_mapped
%     ADJ_mapped+1
%     ADJ=ADJ_mapped.data.data;
%     whos
%     class(ADJ)
%     pause
% end

%%  Call SpeakEasy to generate primary clusters (and subsequently if options.layers>1)
parse(options,varargin{:});
%rng(options.Results.random_seed);
%rng(randi(10000,1))
for main_iter=1:options.Results.subcluster   %main loop over clustering / subclustering
    
    if main_iter==1  %only need to check ADJ characteristics once
        
        [ADJ is_ADJ_weighted is_ADJ_symmetric]=ADJ_characteristics(ADJ); %some memory properties optimized on characteristics of ADJ
        addOptional(options,'is_ADJ_weighted',is_ADJ_weighted);
        addOptional(options,'is_ADJ_symmetric',is_ADJ_symmetric);
        parse(options,varargin{:});
        options=options.Results;
        
        if options.independent_runs*options.target_partitions>100
            disp('you probably only need max value of 100 paritions, so you may be able to reduce option.independnent runs or options.timesteps')
        end
        
        if options.minclust<3
            error('you set minclust < 3 - this is logically odd and may cause problems, so please increase');         %min size for sub-clustering(if already small, dont need to subcluster)
        end
        
        if options.timesteps<10
            error('you set timesteps < 10 - this is not recommended - go for 20+');
        end
        
        if options.max_threads==1
            error('set parallel equal to desired number of threads, not just 1, which indicates a single thread');
        end
        
        if options.multicommunity==1
            disp('Alert - you set the multicom option equal to 1 - this will NOT enable overlapping community detection.')
            disp('If you want to do that select an integer greater than 1, equal to the # community a node may belong to')
        end
        
    end
    
    
    
    if main_iter==1
        disp('calling main routine at level 1')  %in identically named outputs, last one so multicom will overwrite discree, if it is different
        [partition_tags{main_iter} partition_tags{main_iter} partition_cells{main_iter} partition_cells{main_iter} multicom_nodes_all{main_iter}  median_of_all_NMI{main_iter} confidence_score_temp{main_iter}]=SpeakEasy2_bootstrap(ADJ,main_iter,options);
        %partition_tags" comes out sorted by nodeID, whereas when we created it in sub_iter from partition_cells, it is unsorted, then we sort it
        
        if options.node_confidence==1
            
            for m=1:length(partition_cells{main_iter})   %partition_cells is collection of NodeID's, so max 1000
                %three columns - col1:NodeID (not IDX, except for discrete output!) col2:confidence   col3:cluster label
                nodeID_and_confidence_score_temp{m,1}= [partition_cells{main_iter}{m} confidence_score_temp{main_iter}{m}  repmat(m,length(partition_cells{main_iter}{m}),1)  ];
                %you get specific confidence value for each node in each
                %cluster, so there will be some NodeID's with multiple different confidence scores if multicommunity is enabled
                %practically speaking the fourth column of nodeID_and_confidence_score_temp has to be calculated in separate step
                
            end
            
            confidence_score{main_iter}=(vertcat(nodeID_and_confidence_score_temp{:}));
            
            clear nodeID_and_confidence_score_temp;
            clear confidence_score_temp
        end
        
        
    else  %main_iter>1 subclustering and sub-subclustering etc
        disp(['doing subclustering at level ' num2str(main_iter)])
        
        %could make next line parfor to process subclusters in  parallel, but you never know how large those are, so can be problematic as you could get big replicated subsets of ADJ
        for sub_iter=1:length(partition_cells{main_iter-1})  %go through each of the previously determined clusters, hence{main_iter-1}
            
            current_nodes=partition_cells{main_iter-1}{sub_iter};  %these are the top-level Node ID's of the set of nodes we will subcluster
            
            if length(current_nodes)>options.minclust   %subcluster big clusters; used to have mean cluster density criterion to avoid subclustering a really dense cluster... reasons to do each, you could add back here if desired
                
                %"partition_tags_temp" isn't used directly - just "partition_cells_temp" - which later defines "partition_tags"
                [partition_tags_temp{sub_iter,1} partition_tags_temp{sub_iter,1} partition_cells_temp{sub_iter,1} partition_cells_temp{sub_iter,1} todel1 todel2 confidence_score_temp{sub_iter,1} ]=SpeakEasy2_bootstrap(ADJ(current_nodes,current_nodes),main_iter,options);
                partition_tags_temp{sub_iter,1}(:,1)=current_nodes(partition_tags_temp{sub_iter,1}(:,1));                       %the subclusters you've found are at various locations in the full ADJ - need to adjust indexes to reflect this
                %seems "partition_cells_temp" should be updated to "current_nodes" too, which is done in the for-k loop or else condition
                
                
                %"sub_iter" can also be thought of as idx of cluster we're working on
                for k=1:length(partition_cells_temp{sub_iter,1})  %update to main node ID's
                    if options.node_confidence==1
                        nodeID_and_confidence_score_temp{sub_iter,1}{k,1}= [current_nodes(partition_cells_temp{sub_iter,1}{k}) confidence_score_temp{sub_iter}{k}   ];
                    end
                    partition_cells_temp{sub_iter,1}{k}=    current_nodes(partition_cells_temp{sub_iter,1}{k});  %"k" has the main clusters of the sub cluster, so really subsubclusters since we're working on a select part of the ADJ to begin with
                    
                end
                
                
            else %if module is too small for subclustering
                partition_cells_temp{sub_iter,1}{1}=current_nodes;
                if options.node_confidence==1
                    nodeID_and_confidence_score_temp{sub_iter,1}{1}=[current_nodes(:) zeros(length(current_nodes),1)];
                    %would like to know what clusterID each confidence score will be a member
                    %of, since you might have overlapping clusters, meaning the same nodeID will have mulitle confidence scores, but final clusterID's
                    %aren't assigned until we're out of this loop
                end
            end
            
            if options.node_confidence==1
                tracker_confid=(vertcat(nodeID_and_confidence_score_temp{:}));
            end
            tracker_partition_cells=vertcat(partition_cells_temp{:});
            
        end
        
        
        %now we've made it past the sub_iter loop, we're going to combine the results (subclusters) so they're more convenient
        partition_cells{main_iter}=vertcat(partition_cells_temp{:});
        partition_cells_temp=[];
        
        if options.node_confidence==1
            confidence_score{main_iter}=(vertcat(nodeID_and_confidence_score_temp{:}));
            confidence_score_temp=[];
        end
        
        for m=1:length(partition_cells{main_iter}) %indices stored in the mth cell get assigned label==m
            nodeID_and_new_cluster_label{m,1}=[partition_cells{main_iter}{m} repmat(m,length(partition_cells{main_iter}{m}),1) ]; %two cols - nodeID and assignment
            if options.node_confidence==1
                nodeID_and_new_confidence{m,1}  = [confidence_score{main_iter}{m} repmat(m,length(partition_cells{main_iter}{m}),1) ]     ;
                %nodeID's will be in sync:          nodeID_and_new_confidence{m,1}  = [partition_cells{main_iter}{m} confidence_score{main_iter}{m} repmat(m,length(partition_cells{main_iter}{m}),1) ]
            end
        end
        if options.node_confidence==1
            confidence_score{main_iter}=(vertcat(nodeID_and_new_confidence{:})); %was below
        end
        
        partition_tags{main_iter}=vertcat(nodeID_and_new_cluster_label{:});
        
        
    end %for if/else main_iter processing
    
    
    
    convenient_order{main_iter}=vertcat(partition_cells{main_iter}{:});
    partition_tags{main_iter}=sortrows(partition_tags{main_iter});
    
end


if options.max_threads>0 && options.autoshutdown==1
    delete(gcp('nocreate'))
end

if options.node_confidence==1
    save([options.filename '_node_confidence.mat'], 'confidence_score')
end


if nargout==0
    
    save([options.filename '.mat'], 'partition_tags', 'partition_cells', 'convenient_order')
else
    varargout{1}=partition_tags ;
    varargout{2}=partition_cells ;
    varargout{3}=convenient_order ;
end

if options.multicommunity>1
    
    
    if nargout>3
        varargout{4}=multicom_nodes_all;
    else
        save([options.filename '_multicom_nodes_all.mat'], 'multicom_nodes_all')
    end
    
end
