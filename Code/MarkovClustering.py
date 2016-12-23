__author__ = 'Pranshu'
import numpy as np

def main():
    # Input parameters
    datasets = ['attweb_net.txt', 'physics_collaboration_net.txt', 'yeast_undirected_metabolic.txt']
    loop_var = True
    while loop_var:
        dataset = input("Enter dataset complete file name ")
        if dataset not in datasets:
            print("Incorrect dataset specified ")
        else:
            break

    while loop_var:
        exp = input("Enter Expansion parameter (natural numbers > 1 as input) ")
        if exp.isdigit() and int(exp) > 1:
            exp = int(exp)
            break
        else:
            print("Incorrect Expansion parameter specified.")

    while loop_var:
        inf = input("Enter Inflation parameter > 1 ")
        if float(inf) > 1:
            inf = float(inf)
            break
        else:
            print("Incorrect Inflation parameter specified.")

    #Load the dataset
    if dataset == "physics_collaboration_net.txt":
        edges = np.loadtxt(dataset, dtype=str)
    else:
        edges = np.loadtxt(dataset, dtype=int)

    # Assign index to each graph node
    node_index = {}
    idx = -1
    for edge in edges:
        if (edge[0] not in node_index):
            idx += 1
            node_index[edge[0]] = idx
        if (edge[1] not in node_index):
            idx += 1
            node_index[edge[1]] = idx
    total_nodes = node_index.get(max(node_index, key=node_index.get)) + 1
    print("No of Nodes",total_nodes)

    # Create Adjacency Matrix
    AdjacencyMatrix = np.zeros((total_nodes, total_nodes), dtype=int)
    for edge in edges:
        AdjacencyMatrix[node_index.get(edge[0])][node_index.get(edge[1])] = 1
        AdjacencyMatrix[node_index.get(edge[1])][node_index.get(edge[0])] = 1

    # Add Self Loop
    for i in range(0, total_nodes):
        AdjacencyMatrix[i][i] = 1

    #Normalize the data
    AdjacencyMatrix = AdjacencyMatrix/np.sum(AdjacencyMatrix, axis=0)
    # print("Next",AdjacencyMatrix)

    for itr in range(100):
        # Store copy of old data
        prev = AdjacencyMatrix
        # Expand
        AdjacencyMatrix = np.linalg.matrix_power(AdjacencyMatrix, exp)
        # Inflate
        AdjacencyMatrix = np.power(AdjacencyMatrix, inf)
        # Normalize
        AdjacencyMatrix = AdjacencyMatrix/np.sum(AdjacencyMatrix, axis=0)
        # Check Convergence
        if(np.array_equal(AdjacencyMatrix, prev)):
            # print("MCL converged in iteration ", itr+1)
            break
        # Prune the data
        for index, val in np.ndenumerate(AdjacencyMatrix):
            if (val <= 0.0000001):
                AdjacencyMatrix[index[0]][index[1]] = 0

    # Find Clusters
    clusters = set([])
    cluster_map = {}
    for row in AdjacencyMatrix:
        myList = np.where(row > 0)[0]
        if myList.any():
            clusters.add(tuple(myList))
    # print(clusters)
    print("Number of Clusters formed", len(clusters))
    for idx, clusters in enumerate(clusters):
        for node in clusters:
            cluster_map[node+1] = idx+1
    # print(cluster_map)

    # Write in clu file
    fileObject = open("file_M"+str(exp)+"E"+str(inf)+".clu", 'w')
    fileObject.write("*Partition PartitionName")
    fileObject.write("\n")
    fileObject.write("*Vertices " + str(int(total_nodes)))
    fileObject.write("\n")

    for key in cluster_map:
        fileObject.write(str(cluster_map.get(key)))
        fileObject.write("\n")

main()