import xxhash
from ete3 import Tree
from Neighbour_Joining import NJ


def read_fasta(input_f):

    seq = {}
    with open(input_f, 'r') as in_f:
        for line in in_f:
            if line.startswith('>'):
                contigs = line.replace('>', '').strip()  # Obtain names for each contigs without '>' and '\n'
                seq[contigs] = ''  # Initialization of values for the keys of dictionary
            else:
                seq[contigs] += line.replace('\n', '').strip()  # Obtain the sequence corresponding to the contigs

    return seq  # Return a dict whose keys are names of contigs and values are sequences


def sliding_window(seq, k):  # Obtain a list contains all k-mers

    windows = {}  # Initialisation of a dictionary, whose keys are contigs and whose values are k-mers

    for contigs in seq:
        L = len(seq[contigs])
        windows[contigs] = []
        for i in range(L - k + 1):
            windows[contigs].append(seq[contigs][i:i + k])

    windowsList = list(windows.values())  # Turn the dictionary into a list
    windowsList_combine = []

    for i in range(len(windowsList)):  # Combine all subsequences together, important for the contigs
        for item in windowsList[i]:
            windowsList_combine.append(item)

    return windowsList_combine


def calculate_Jaccard_distance(seqA, seqB, k):  # Calculate the distance based on the strings(k_mers)

    intersection = 0
    A_kmers = sliding_window(seqA, k)  # Obtain the list of k-mers of seqA
    B_kmers = sliding_window(seqB, k)

    for i in range(len(A_kmers)):
        for j in range(len(B_kmers)):
            if A_kmers[i] == B_kmers[j]:
                intersection += 1  # Get the number of intersections of the two k-mers

    union = len(A_kmers) + len(B_kmers) - intersection
    # If x+intersection=len(A) and y+intersection=len(B), then union=x+y+intersection=len(A)+len(B)-intersection
    distance = 1 - intersection / union  # The Jaccard distance

    return distance


def calculate_Jaccard_distance_hash(seqA, seqB, L):  # Calculate the distance based on the numbers(sketch)

    intersection = 0

    for i in range(L):  # L is the sketch size
        for j in range(L):
            if seqA[i] == seqB[j]:
                intersection += 1

    union = 2 * L - intersection
    distance = 1 - intersection / union

    return distance


def hash_function(SeqList):  # Get the list of k-mers encoded with a hash function

    hashList = []

    for i in range(len(SeqList)):
        x = xxhash.xxh32(SeqList[i]).intdigest()
        hashList.append(x)

    return hashList


def sketch(seq, sketch_size, filename, k):

    SeqList = sliding_window(seq, k)  # Get the list of all k-mers from sequences in the dictionary
    hashList = hash_function(SeqList)  # Get the hash values corresponding to the sequence
    hashList.sort()  # Sort hashList from smallest to largest
    hashList1 = hashList[:sketch_size]  # Obtain the specified number of hash values

    with open(filename, 'w') as out_f:
        for line in hashList1:
            out_f.write(str(line) + '\n')  # Write the specified number of hash values to the file


def load_sketch(filename):
    # Read the contents of a sketch and return a list of hash values

    arr_seqs = []

    with open(filename, 'r') as in_f:
        for line in in_f:
            arr_seqs.append(line.strip())  # Put hash values into the list

    return arr_seqs


if __name__ == '__main__':
    input_reference1_f = 'R6.fa'
    input_reference2_f = 'TIGR4.fa'
    input_draft1_f = '14412_3#82.contigs_velvet.fa'
    input_draft2_f = '14412_3#84.contigs_velvet.fa'
    seq1 = read_fasta(input_reference1_f)  # Read the fasta file and get a dictionary containing the sequence
    seq2 = read_fasta(input_reference2_f)
    seq3 = read_fasta(input_draft1_f)
    seq4 = read_fasta(input_draft2_f)
    k = 14  # k of k-mers
    L = 1000  # sketch size
    sketch(seq1, L, 'sketch1.txt', k)  # Write the specified number of k-mers to the file as hash values
    sketch(seq2, L, 'sketch2.txt', k)
    sketch(seq3, L, 'sketch3.txt', k)
    sketch(seq4, L, 'sketch4.txt', k)
    A = load_sketch('sketch1.txt')
    B = load_sketch('sketch2.txt')
    C = load_sketch('sketch3.txt')
    D = load_sketch('sketch4.txt')
    dAB = calculate_Jaccard_distance_hash(A, B, L)
    dAC = calculate_Jaccard_distance_hash(A, C, L)
    dAD = calculate_Jaccard_distance_hash(A, D, L)
    dBC = calculate_Jaccard_distance_hash(B, C, L)
    dBD = calculate_Jaccard_distance_hash(B, D, L)
    dCD = calculate_Jaccard_distance_hash(C, D, L)
    print(dAB, dAC, dAD, dBC, dBD, dCD)

    M_labels = ['A', 'B', 'C', 'D']
    M = [
        [0, dAB, dAC, dAD],
        [dAB, 0, dBC, dBD],
        [dAC, dBC, 0, dCD],
        [dAD, dBD, dCD, 0],
    ]
    newick_NJ = '{};'.format(NJ(M, M_labels))  # Convert the tuple obtained by the NJ function into newick format
    print(Tree(newick_NJ))

""" Comparing isolates:
    1.
    L = 10: DISTANCE = [1.0 1.0 1.0 1.0 1.0 0.5714285714285714]
    L = 100: DISTANCE = [0.9417989417989419 0.9949748743718593 0.9949748743718593
                         0.98989898989899 0.98989898989899 0.5915492957746479]
    L = 500: DISTANCE = [0.9572471324296141 0.98989898989899 0.9878542510121457
                         0.9847715736040609 0.9878542510121457 0.5549132947976878]
    L = 1000: DISTANCE = [0.9434759640781828 0.9832231825114387 0.9816700610997964 
                          0.9868287740628167 0.9868287740628167 0.5228951255539143]
    L = 1500: DISTANCE = [0.5556090515166106 0.9851150202976996 0.9847715736040609 
                          0.9868287740628167 0.9864864864864865 0.5322896281800391]
    L = 2000: DISTANCE = [0.6741796486576068 0.9858012170385395 0.985544002028912 
                          0.9883662114314619 0.9878542510121457 0.5013113525665043]
    L = 5000: DISTANCE = [0.8467304809133894 0.9876493217250456 0.9881614894262876 
                          0.9870340356564019 0.9876493217250456 0.5090204264201581]
                          
MinHash distances are generally relatively approximate compared to full distances, but in some cases MinHash distances 
lack precision and appear to fluctuate widely as the sketch size changes. This is not the case with full distances, 
which are more accurate, but difficult to achieve for large genomic data. Therefore, I believe that full distances 
should be used for small-scale data, and that the sketch size should be increased to obtain more accurate MinHash 
distances for large-scale data.

With the change of sketch size, most MinHash distance don't change significantly, but there are still distances which 
change a lot, for example, dAB changes from 0.56 to 1, so we can say that sketch size can influence the result of 
MinHah distance.


    2.
    for L =1000       
      /-A 
   /-|
  |  |   /-C
--|   \-|
  |      \-D
  |
   \-B
   
    for L = 10
      /-A
   /-|
  |   \-B
--|
  |   /-C
   \-|
      \-D
      
   for L = 100
         /-A
      /-|
   /-|   \-B
  |  |
--|   \-C
  |
   \-D
  
   for L = 500
         /-A
      /-|
   /-|   \-B
  |  |
--|   \-C
  |
   \-D
     
   for L = 2000
         /-A
      /-|
   /-|   \-B
  |  |
--|   \-D
  |
   \-C
   
   for L = 5000
      /-A
   /-|
  |   \-B
--|
  |   /-C
   \-|
      \-D
      
If 'sketch size' = 1000, their relationship in newick format is '((A,(C,D)),B);', This means that C and D are the 
closest relatives, A and (C,D) are relatively closely related and B and (A,C,D) are relatively distantly related.

When we change the sketch size, their MinHash distance changes as well. We found that as the distance changes the
NJ tree also changes, and in some cases even when the distance does not change much the form of the NJ tree changes 
significantly, for example the newick for L = 500 is '(((A,B),C),D);' while the newick for L = 1000 where the distance 
does not change significantly is '((A,(C,D)),B);'. This shows that the NJ tree is sensitive to distance.
"""