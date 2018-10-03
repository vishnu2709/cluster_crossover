# cluster_crossover
Program to generate a new cluster from two given clusters

Current Status:- Works for all the current test files (Ni4O1, Ni4O4, Ni4O5) <br>
Potential Bugs/Weaknesses :- The crossover uses a tolerance value to get a rough cluster. The rough cluster is then used to find the true cluster. However, this takes the assumption that the top most atom of the cluster will be a part of the rough cluster. While this is indeed a reasonable assumption, there is a theoretical possibility that the top most atom may not be a part of the rough cluster. Due to lack of any such test files which show such behaviour, we are not very worried about this issue, but it is something that we will have to work on in the near future. <br> <br>

One possible fix is to make the tolerance check stricter. This will ensure that the entire cluster (along with some substrate atoms) will be part of the rough cluster. However, even this is not a guaranteed solution
