import spkmeansmodule
import numpy as np
import sys

np.random.seed(0)

def main(args):
    #args = [k, goal ,filename]
    if len(args) != 3 and len(args) != 2:
        invalid_input_error()
    try:
        k = int(args[0])
    except ValueError:
        invalid_input_error()
    goal = args[1] #check if goal is one of the options... and act accordingly
    input_path = args[2]
    vectors = read_from_file(input_path)
    N = len(vectors)
    d = len(vectors[0])
    if goal == "spk":
        spk(k, N, d, vectors)
        exit(0)
    elif goal == "wam":
        spkmeansmodule.cal_by_algo(0, vectors, N, d)
        exit(0)
    elif goal == "ddg":
        spkmeansmodule.cal_by_algo(1, vectors, N, d)
        exit(0)
    elif goal == "lnorm":
        spkmeansmodule.cal_by_algo(2, vectors, N, d)
        exit(0)
    elif goal == "jacobi":
        spkmeansmodule.cal_by_algo(3, vectors, N, d)
        exit(0)
    invalid_input_error()

def invalid_input_error():
    print("Invalid Input!")
    exit(1)


def general_error():
    print("An Error Has Occurred")
    exit(1)

def spk(k, N, d, input_vectors):

    # calculate the T matrix out of the input vectors (steps 1-5 of the spk algorithm).
    # if the input k is 0, then we determine k (step 3), otherwise the k is the input k.
    t_mat = spkmeansmodule.cal_spk_vectors(input_vectors, N, d, k)
    if k == 0:
        if len(t_mat) == 0:
            general_error()
        k = int(t_mat.pop())
    else:
        t_mat.pop()

    K_means(t_mat, k)

def K_means(vectors, k):
    vectors = np.array(vectors)
    n, d = vectors.shape
    centroids = np.zeros((k, d))

    # change centorids
    res_indices = k_means_pp(vectors, centroids, k)
    print(*res_indices, sep=",")

    # convert centroids and vectors to lists, for the c module
    centroids_lst = centroids.tolist()
    vectors_lst = vectors.tolist()

    # do the kmeans algorithm with the chosen first k centroids
    spkmeansmodule.fit(vectors_lst, centroids_lst, n, d, k, 300, 0)
    # print in c - print the k centorids as rows

def k_means_pp(vectors, centroids, k):
    # vectors - a matrix N x d, of N vectors, each one d coordinates. from the input files.
    # centroids - a matrix k x d, initialized with zeros, to store the k centroids.
    # the function changes "centroids" - to contain the k centroids calculated by the kmeans++ algorithm,
    # out of the N vectors.
    # the function returns a list of indices of the k vectors, that was chosen to be the centroids.

    N = len(vectors)
    d = len(vectors[0])
    i = 0
    indices = np.arange(N)  # N vectors: 0,1, ... , N-1

    # initialize k-size lists to store the chosen vectors indices
    chosen_indices = [None] * k

    # select first centroid randomly from the N vectors
    select_random_vector(vectors, centroids, indices, chosen_indices, i)

    # initialize N-size lists to store D values  and probabilty values for each vecto
    D_lst = [None] * N
    P_lst = [None] * N

    while i < k - 1:
        for j in range(N):
            vector = np.copy(vectors[j])
            D_vector = calculate_D(vector, centroids, i)
            D_lst[j] = D_vector

        D_sum = sum(D_lst)

        for j in range(N):
            P_lst[j] = D_lst[j] / D_sum

        i += 1

        # choose the vector randomly, vector in index [i] has P_lst[i] probablity to be chosen
        select_random_vector(vectors, centroids, indices, chosen_indices, i, probablity=P_lst)

    return chosen_indices


def calculate_D(vector, centroids, i):
    min_value = -1
    for j in range(i + 1):
        norm = cal_norm_squared(vector, centroids[j])
        if (min_value == -1 or norm < min_value):
            min_value = norm
    return min_value


def cal_norm_squared(v1, v2):
    args_sum = 0
    for i in range(len(v1)):
        args_sum += ((v1[i] - v2[i]) ** 2)
    return args_sum

def select_random_vector(vectors, centroids, indices, chosen_indices, i, probablity = None):
    # select centroid randomly from the N vectors
    if probablity == None: # each vector has the same probablity to be chosen
        rand_index = np.random.choice(indices)
    if probablity != None: # vector in index [i] has probablity[i] probablity to be chosen
        rand_index = np.random.choice(indices, p = probablity)
    centroids[i] = np.copy(vectors[rand_index])
    chosen_indices[i] = rand_index

def read_from_file(path):
    centroids = []
    try:
        lines_in_file = open(path, 'r').readlines()
    except FileNotFoundError:
        invalid_input_error()
    except:
        general_error()
    for line in lines_in_file:
        if not line:
            continue
        a = line.split(',')
        values = [float(strval) for strval in a]
        centroids.append(values)
    return centroids

def cal_norm(v1, v2):
    args_sum = 0
    for i in range(len(v1)):
        args_sum += (abs(v1[i] - v2[i]) ** 2)
    return np.sqrt(args_sum)


if __name__ == "__main__":
    main(sys.argv[1:])