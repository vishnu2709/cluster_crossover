# cluster_crossover
Program to generate a new cluster from two given clusters

Current Status:- Works for all the current test files (Ni4O1, Ni4O4, Ni4O5)
Potential Bugs/Weaknesses :- The crossover uses a tolerance value to get a rough cluster. The rough cluster is then used to find the true cluster. However, this takes the assumption that the top most atom of the cluster will be a part of the rough cluster. While this is indeed a reasonable assumption, it is possible that it may not be so. Due to lack of any such test files which show such behvaiour, we are not very worried about this issue, but it is something that we will have to work on in the near future.
