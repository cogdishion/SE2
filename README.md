Looking to cluster data or networks reliably and accurately?  Because data and networks vary so widely, we designed a method to handle a range of data without parameter tweaking. <br />
<br />
Instructions for clustering networks with SpeakEasy2: Champagne. 
Tested under 2018b.  Parallel toolbox not required, but helpful for large datasets.

To cluster data, 
Download SpeakEasy2: Champagne (SE2), navigate to the SpeakEasy2 folder in matlab (or put it on your matlab path) and type: <br />
load ADJdemo %at the matlab prompt to load an example network.  <br />
Then run SpeakEasy2 by typing… <br />
SpeakEasy2(ADJ) <br />
You’ve clustered some data!  You’ll see results have been saved in the folder as SpeakEasy2_results.mat.  If you’d like to verify that SpeakEasy2 is producing plausible results, copy this into the matlab terminal:

load SpeakEasy2_results <br />
figure('Position',[0 0 1000 400]) <br />
subplot(1,2,1); imagesc(ADJ) <br />
title('original input ADJ') <br />
xlabel('nodes') <br />
ylabel('nodes') <br />
subplot(1,2,2); imagesc(ADJcolorful(convenient_order{1},convenient_order{1})) <br />
title('ADJ reorganized by SE2 clusters') <br />
xlabel('reordered nodes') <br />
ylabel('reordered nodes') <br />
colormap(cmap); <br />

![example_SE2_results](https://user-images.githubusercontent.com/46224527/160920014-849a3173-c757-484d-9c93-93ac46d0c5f2.jpg)

On the left you can see the input adjacency matrix for the network input to SE2, while the image on the right shows the same data arranged according to the clusters found by SE2.  You can see that nodes have been organized into groups (along the diagonal).  These have been color-coded with the ground truth communities, and you can see in most cases the original communities are received (nodes in a group tend to be the same color).

If you want to apply SE2 to some of your own data, just substitute your adjacency matrix for the ADJ from ADJdemo.mat.

If you are interested in further customizing the application of SpeakEasy2, continue reading, but the basic function above works for many data.



Runtime - 
If you are running a multi-core machine, have the matlab parallel toolbox installed, and wish to explicitly parallelize SE2 for faster execution, use the max_threads parameter.
  For instance, >>SpeakEasy2(ADJdemo,’max_threads’,50) will utilize up to 50 threads.  Performance gains from multiple cores/threads are substantial.  However, a copy of the ADJ is sent to each worker, so for very large networks, you may have to limit the number of workers to fit in memory.


Accuracy - 
The more independent partitions that are generated by SE2, the more stable and accurate the final partition solution will be.  In practice, 10 estimated partitions are generally sufficient, but if you wish to generate more partitions to slightly improve accuracy, you can use the independent_runs parameter to increase this.
  For instance, >>SpeakEasy2(ADJdemo,’independent_runs’, 20) will double the default number of estimated partitions.


Overlapping clusters - 
By default SE2 returns disjoint communities, but you can enable varying levels of overlapping output with the “multicommunity” parameter.
  For instance, >>SpeakEasy2(ADJdemo,’multicommunity’, 4) allows nodes to be members in up to 4 communities.  The usual saved file will indicate nodes with multiple communities (it’s length will be greater than the number of nodes), and if you just want a list of the multi-community nodes, that is saved as well.


Subclustering - 
Occasionally your data will have multiple scales, or large batch effects that might look like clusters, and you need to get underneath them for meaningful results.  (Often the case for bulk RNAseq, but generally not single cell RNAseq)
You can apply SpeakEasy2 to each of the top-tier clusters with subcluster.  This may be done multiple times if you wish (i.e. sub-sub clustering).
   For instance, >>SpeakEasy2(ADJdemo,’subcluster’, 3) produces three levels of clustering - the usual top-level clusters, as well as sub-clusters and sub-sub clusters.  You can access lower levels of clusters as follows:<br />
  >> load SpeakEasy2_results <br />
partition_tags{1}; %these are the top-level communities <br />
partition_tags{2}; %sub clusters <br />
partition_tags{3}; %sub-sub clusters <br />

You may not wish to split communities below a certain size, which is facilitated by minclust.  For instance >>SpeakEasy2(ADJdemo,’subcluster’, 3,’minclust’,20) will not subcluster primary clusters with less than 20 nodes.


I/O -
If you wish to access the partition without loading a saved file, you can run a command with 1-3 outputs:

load ADJdemo <br />
[node_tags  node_groups convenient_order]=SpeakEasy2(ADJ);

In this case node_tags{1} will contain two columns - the first is the node ID and the 2nd is an arbitrary numeric cluster assignment.  For instance (hypothetically) it could contain: <br />
1 1 <br />
2 1 <br />
3 2 <br />
4 3 <br />
5 3 <br />

-the first output matrix ('node_tags') above denotes a 5-node network, with three clusters.

The second output ('node_groups') lists the nodes associated with the three clusters (it is a cell of length 3).

The third output ('convenient_order') orders nodes for convenient visualization (i.e. imagesc(ADJ(convenient_order{1},convenient_order{1}) is visually interpretable).


Background reading -

This original SpeakEasy algorithm is [described and still available.](https://www.cs.rpi.edu/~szymansk/SpeakEasy/index.html)  and results [published](https://www.nature.com/articles/srep16361).

