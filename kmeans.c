#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int N = 0, K = 0, iterations = 0;
int vector_size = 0;
double epsilon = 0.001;

typedef double *Data; /*Data is an alias to a pointer to a double*/

struct linked_list {
    Data d;
    struct linked_list *next;
};

typedef struct linked_list Element; /* Element is an alias to a linked_list */
typedef Element *Link;              /* Link is defined as a pointer to an Element */

Data calc_vector(const char *input_line, int line_size, int vector_size) {
    Data curr_vector = (Data) calloc(vector_size, sizeof(double));
    if (curr_vector == NULL)
    {
        printf("\nAn Error Has Occurred\n");
        exit(0);
    }
    int vector_entrance_index = 0, start = 0, i, j, entrance_length;
    char *vector_entrance = (char *) calloc(1, sizeof(char));
    for (i = 0; i < line_size; i++) {
        if (input_line[i] == ',' || input_line[i] == '\n') {
            entrance_length = i - start;
            vector_entrance = (char *) realloc(vector_entrance, (entrance_length + 1) * sizeof(char)); /* Add +1 for the null terminator */
            if (vector_entrance == NULL)
            {
                printf("\nAn Error Has Occurred\n");
                exit(0);
            }
            for (j = 0; j < entrance_length; j++) { /* copy the current vector entrance */
                vector_entrance[j] = input_line[start + j];
            }
            vector_entrance[entrance_length] = '\0'; /* add the null terminator to the end of the vector entrance */
            curr_vector[vector_entrance_index] = atof(vector_entrance); /* set the entrance of the vector */
            start = i + 1;
            vector_entrance_index++;
        }
    }
    free(vector_entrance);
    return curr_vector;
}

int calc_vector_size(const char *input_line, int line_size) {
    int size = 1, i;
    for (i = 0; i < line_size; i++) {
        if (input_line[i] == ',') {
            size++;
        }
    }
    return size;
}

Link initialize_linked_list(char *input_line, int vector_size, int line_size) {
    Link head = (Link) calloc(1, sizeof(Element));
    if (head == NULL)
    {
        printf("\nAn Error Has Occurred\n");
        exit(0);
    }
    head->d = calc_vector(input_line, line_size, vector_size);
    head->next = NULL;
    return head;
}

void link_next_line(Link curr_element, char *input_line, int vector_size, int line_size) {
    Link next_element = (Link) calloc(1, sizeof(Element));
    if (next_element == NULL)
    {
        printf("\nAn Error Has Occurred\n");
        exit(0);
    }
    next_element->d = calc_vector(input_line, line_size, vector_size);
    next_element->next = NULL;
    curr_element->next = next_element;
}

Data *list_to_matrix(Link head, int n, int vector_size) { /* n = the total number of vectors in the given input */
    Data p;
    Data *a;
    int i, j;
    p = (Data) calloc(n * vector_size, sizeof(double)); /* initialize a one dimensional array of size n * vector_size */
    a = (Data *) calloc(n, sizeof(Data));               /* initialize an array that represents the rows of the matrix */
    if (a == NULL || p == NULL) {
        printf("\nAn Error Has Occurred\n");
        exit(0);
    }
    for (i = 0; i < n; i++) {
        a[i] = p + i * vector_size; /* a[i] = the memory address of the first entrance of the i'th vector */
        for (j = 0; j < vector_size; j++) {
            a[i][j] = head->d[j];
        }
        head = head->next;
    }
    return a;
}

void free_linked_list(Link head) {
    Link tmp;
    while (head != NULL) {
        tmp = head;
        head = head->next;
        free(tmp->d);
        free(tmp);
    }
}

Data *handle_stream(const char *file_name) {
    FILE *file;
    char *input_line = NULL;
    size_t line_size = 0;
    int read;
    int counter = 0;
    Link head = NULL, curr_element = NULL;
    file = fopen(file_name, "r");
    if (file == NULL) {
        printf("\nAn Error Has Occurred\n");
        exit(0);
    }
    while ((read = getline(&input_line, &line_size, file)) != -1) {
        if (counter == 0) {
            vector_size = calc_vector_size(input_line, read);
            head = initialize_linked_list(input_line, vector_size, read);
            curr_element = head;
        } else {
            link_next_line(curr_element, input_line, vector_size, read);
            curr_element = curr_element->next;
        }
        counter++;
        line_size = 0; /* Reset line_size to 0 for the next iteration */
    }
    free(input_line);
    fclose(file);
    N = counter;
    Data *matrix = list_to_matrix(head, N, vector_size);
    free_linked_list(head);
    return matrix;
}

typedef struct cluster {
    Data centroid;
    Data members;
} cluster;

double distance(Data v1, Data v2) {
    double sum = 0;
    int i;
    for (i = 0; i < vector_size; i++) {
        sum += pow(v1[i] - v2[i], 2);
    }
    return sqrt(sum);
}

int find_cluster(cluster **clusters, Data vector) {
    int index = -1, i;
    double curr_distance;
    double min_distance = pow(2, sizeof(double) * 8 - 2) * (1 - pow(2, -52)); /* set min_distance to the max possible double value */
    Data curr_centroid;
    for (i = 0; i < K; i++) {
        curr_centroid = clusters[i]->centroid;
        curr_distance = distance(curr_centroid, vector);
        if (curr_distance < min_distance) {
            min_distance = curr_distance;
            index = i;
        }
    }
    return index;
}

void add_to_cluster(cluster **clusters, int *cluster_members_counter_copy, int index, Data vector) {
    int i;
    Data vector_copy = (Data) malloc(vector_size * sizeof(double)); /* allocate memory for vector_copy */
    if (vector_copy == NULL) {
        printf("\nAn Error Has Occurred\n");
        exit(0);
    }
    for (i = 0; i < vector_size; i++) {
        vector_copy[i] = vector[i]; /* copy vector values to vector_copy */
        clusters[index]->members[(cluster_members_counter_copy[index] - 1) * vector_size + i] = vector_copy[i];
    }
    cluster_members_counter_copy[index]--;
}

    /* decrease the cluster's members count. Why? because essentially we use cluster_members_counter twice:
    1. to calculate how much space to allocate for each cluster members
    2. to keep track of where to paste each vector in members array
    (we paste the vectors in reversed order - from the last position to the
    first position) */

void set_cluster_members_counters(cluster **clusters, Data *matrix, int *cluster_members_counter, int* cluster_members_counter_copy) {
    int i, index;
    for (i = 0; i < N; i++) /* go over each vector in the matrix */
    {
        index = find_cluster(clusters, matrix[i]);
        cluster_members_counter[index]++;
        cluster_members_counter_copy[index]++;
    }
}

void allocate_memory_for_cluster_members(cluster **clusters, int *cluster_members_counter) {
    int i;
    for (i = 0; i < K; i++) {
        clusters[i]->members = (Data) calloc(vector_size * cluster_members_counter[i], sizeof(double));
    }
}

void add_vectors(Data v1, Data v2) {
    int i;
    for (i = 0; i < vector_size; i++) {
        v1[i] = v1[i] + v2[i];
    }
}

Data calc_centroid(cluster **clusters, const int *cluster_members_counter, int index) {
    Data sum = (Data) calloc(vector_size, sizeof(double)); /* initialize 'sum' as the zero vector */
    int num_of_members = cluster_members_counter[index];
    int i;
    if (num_of_members == 0) /* if the cluster is empty then its centroid is not changed */
    {
        return clusters[index]->centroid;
    }
    for (i = 0; i < num_of_members; i++) /* sum all the members in the cluster */
    {
        add_vectors(sum, clusters[index]->members + i * vector_size);
    }
    for (i = 0; i < vector_size; i++) /* divide each entrance of the vector by the number of members in the cluster */
    {
        sum[i] = sum[i] / num_of_members;
    }
    return sum;
}

int main(int argc, char *argv[]) {
    Data *matrix;
    K = atoi(argv[1]);
    if (argc == 4) /* max iterations parameter is given */
    {
        iterations = atoi(argv[2]);
        matrix = handle_stream(argv[3]); /* each row in matrix is a vector given in input */
    } else if (argc == 3) /* max iterations parameter not given */
    {
        iterations = 200;
        matrix = handle_stream(argv[2]); /* each row in matrix is a vector given in input */
    } else /* not enough arguments given */
    {
        printf("Invalid command line arguments.\n");
        return 1;
    }
    if (K <= 1 || K >= N) {
        printf("Invalid number of clusters!\n");
    }
    int valid = 0, i;

    /* allocate memory for 'clusters' to be an array of K pointers to objects of type cluster */
    cluster **clusters = (cluster **) calloc(K, sizeof(cluster *));
    for (i = 0; i < K; i++) {
        clusters[i] = (cluster *) calloc(1, sizeof(cluster)); /* allocate memory for each cluster object in the clusters array */
        clusters[i]->centroid = (double*) calloc(vector_size, sizeof(double)); /* allocate memory for each centroid */
        int j;
        for (j = 0; j < vector_size; j++) {
            clusters[i]->centroid[j] = matrix[i][j]; /* copy the first k vectors from matrix to be the centroids */
        }
    }
    for (i = 0; i < K; i++) { /* allocate an initial memory for every cluster members array */
        clusters[i]->members = (Data) calloc(1, sizeof(double));
        if (clusters[i]->members == NULL)
        {
            printf("\nAn Error Has Occurred\n");
            exit(0);
        }
    }
    while (iterations > 0 && !valid) { /* line 67 on python file */
        valid = 1;
        int i, index;
        int *cluster_members_counter = (int *) calloc(K, sizeof(int)); /* allocate memory for the c_m_c array */
        int *cluster_members_counter_copy = (int *) calloc(K, sizeof(int)); /* allocate memory for the c_m_c_c array */
        if (cluster_members_counter == NULL || cluster_members_counter_copy == NULL)
        {
            printf("\nAn Error Has Occurred\n");
            exit(0);
        }
        set_cluster_members_counters(clusters, matrix, cluster_members_counter, cluster_members_counter_copy);
        allocate_memory_for_cluster_members(clusters, cluster_members_counter);
        for (i = 0; i < N; i++) /* for each vector in the given in input */
        {
            index = find_cluster(clusters, matrix[i]); /* find the vector's cluster */
            add_to_cluster(clusters, cluster_members_counter_copy, index, matrix[i]);
        }
        free(cluster_members_counter_copy);
        for (i = 0; i < K; i++) {
            Data new_centroid = calc_centroid(clusters, cluster_members_counter, i);
            if (valid) {
                double delta = distance(new_centroid, clusters[i]->centroid);
                if (delta > epsilon) {
                    valid = 0;
                }
            }
            free(clusters[i]->centroid);
            clusters[i]->centroid = new_centroid;
            free(clusters[i]->members);
        }
        free(cluster_members_counter);
        iterations--;
    }
    int j;
    for (i = 0; i < K; i++) /* for every cluster */
    {
        for (j = 0; j < vector_size; j++) {
            printf("%.4f\t", clusters[i]->centroid[j]);
        }
        printf("\n");
    }
    for (i = 0; i < K; i++) {
        free(clusters[i]);
    }
    free(clusters);
    free(matrix[0]);
    free(matrix);
    return 0;
}
