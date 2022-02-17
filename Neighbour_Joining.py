import numpy as np
import warnings

warnings.filterwarnings('ignore')
# The program generates an alert that does not affect the result, so it can be filtered out directly


def uij(table):  # Calculate ui,j for future calculation of Q-criterion

    u = np.zeros(len(table))
    for i in range(len(table)):
        u[i] = np.sum(table[i]) / (len(table) - 2)

    return u


def Q_criterion(table):  # Calculate Q-criterion for each pair of sequence

    q = np.zeros((len(table), len(table)))
    u = uij(table)
    for i in range(len(table)):
        for j in range(len(table)):
            q[i][j] = table[i][j] - (u[i] + u[j])
            q[i][i] = 0

    return q  # Obtain a Q-criterion matrix


def minimum_index(table):  # Obtain the index corresponding to the smallest value in the table

    min = float("inf")
    x, y = 0, 0

    for i in range(len(table)):
        for j in range(len(table[i])):
            if table[i][j] < min:
                min = table[i][j]
                x, y = i, j

    if y < x:
        x, y = y, x  # Sort the obtained indexes in order

    return x, y


def new_matrix(table):  # Generate the new distance matrix after merging the rows and columns

    q = Q_criterion(table)
    x, y = minimum_index(q)
    # Obtain the index of the minimum Q-criterion
    # Their corresponding rows and columns are also the rows and columns we need to merge
    matrix = np.zeros((len(table), len(table)))  # Initialisation of a zero matrix (new distance matrix)

    for i in range(len(matrix)):
        for j in range(len(matrix)):
            matrix[i][x] = (table[i][x] + table[i][y] - table[y][x]) / 2  # Calculation of new distance
            matrix[x][j] = (table[x][j] + table[y][j] - table[x][y]) / 2
            if i != x and j != x:
                matrix[i][j] = table[i][j]  # Only the distances in the merged rows and columns have changed

            matrix[i][i] = 0
    mat_a = np.delete(matrix, y, axis=0)
    mat_b = np.delete(mat_a, y, axis=1)  # Delete the row and column corresponding to index y

    return mat_b


def NJ(table, labels):

    while len(labels) > 1:  # Loop to get the new distance matrix and the corresponding x,y of the Q-criterion matrix
        q = Q_criterion(table)
        x, y = minimum_index(q)
        table = new_matrix(table)
        labels[x] = "(" + labels[x] + "," + labels[y] + ")"
        # Merge the labels corresponding to indexes x and y into the labels corresponding to x
        del labels[y]  # Delete the labels corresponding to index y

    return labels[0]  # Return to a tuple
