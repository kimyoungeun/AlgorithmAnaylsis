#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <time.h> //youngeun

#define INFINITY 10000

struct adj_list
{
	int v_index;
	int w;
	struct adj_list* next;
};

const int max_num_lists = 512; //  We expect no more than 512 lists
struct adj_list* all_lists;
int num_all_lists = 0; // initialization

enum v_color {
	WHITE, // 0
	GRAY, // 1
	BLACK // 2
};

struct DFS_vertex
{
	int dis; // discovery time
	int fin; // finish time
	int color;
	int pi; // index of parent vertex
};

int getVeticesNames(char** v_names, char* filename);
void getAdjMat(int** Adj_mat, int num_v, char* filename);
void printMat(int** Adj_mat, int num_v);
void getAdjArray(struct adj_list** Adj_array, int num_v, int** Adj_mat);
void printAdjArray(struct adj_list** Adj_array, int num_v, char* v_names);
int DFS_VISIT(struct adj_list** Adj_array, int num_v, struct DFS_vertex* DFS_vertices, int u, int time);
void DFS(struct adj_list** Adj_array, int num_v, struct DFS_vertex* DFS_vertices, int source_index, int* indexOrder);
void printDFS(struct DFS_vertex* DFS_vertices, int num_v, char* v_names);
void getTransAdjMat(int** Adj_mat, int** Adj_mat_trans, int num_v);
void swapItems(int* Array, int i, int j);
int partition(int* Array, int p, int r, int* originalIndex);
void quickSort(int* Array, int p, int r, int* originalIndex);
int findDecendants(struct DFS_vertex* DFS_vertices, int num_v, int rootIndex, int*decendants);
int findSSCs(struct DFS_vertex* DFS_vertices, int num_v, int**SSCs);
void printSSCs(int** SSCs, int num_SSCs, char* v_names);

//youngeun
void heapSort(int* Array, int p, int r, int* originalIndex);
void heapify(int* Array, int size, int mid, int* originalIndex);
void buildMaxHeap(int* Array, int size, int* originalIndex);

void countingSort(int* Array, int p, int r, int* originalIndex);

void dijkstra(int** new_mat, int** Adj_mat, int num_v);

bool bellmanford(int** new_mat, int** Adj_mat, int num_v);

void floydwarshall(int** new_mat, int** Adj_mat, int num_v);

void new_printMat(int** new_mat, int num_v, char* v_names);
//youngeun

int main()
{
	// input data
	char filename[] = "HW4.dat";

	// getting the number and names of vertices
	char* v_names = NULL; // assume each vertex has a single-charactor name
	int num_v = getVeticesNames(&v_names, filename); // num_v: number of all vertices in the given graph

	// getting adjacency matrix
	int** Adj_mat = malloc(sizeof(int*)*num_v);
	for (int i = 0; i<num_v; i++) Adj_mat[i] = malloc(sizeof(int)*num_v);
	getAdjMat(Adj_mat, num_v, filename);
	printf("Adjacency Matrix\n");
	printMat(Adj_mat, num_v);
    
	printf("\n");

	// getting the array of adjacency list from the adjacency matrix
	struct adj_list** Adj_array;
	Adj_array = malloc(sizeof(struct adj_list*)*num_v);
	all_lists = malloc(sizeof(struct adj_list)*max_num_lists);
	getAdjArray(Adj_array, num_v, Adj_mat);
	printAdjArray(Adj_array, num_v, v_names);

	printf("\n");

	// DFS
	struct DFS_vertex* DFS_vertices;
	DFS_vertices = malloc(sizeof(struct DFS_vertex)*num_v);
	int source_index = 0;
	printf("With %c as the source vertex,\n", v_names[source_index]);
	DFS(Adj_array, num_v, DFS_vertices, source_index, NULL);
	printDFS(DFS_vertices, num_v, v_names);

	printf("\n");
    
	// Transpose of adjacency matrix
	int** Adj_mat_trans = malloc(sizeof(int*)*num_v);
	for (int i = 0; i<num_v; i++) Adj_mat_trans[i] = malloc(sizeof(int)*num_v);
	getTransAdjMat(Adj_mat, Adj_mat_trans, num_v);
	printf("Transpose of Adjacency Matrix\n");
	printMat(Adj_mat_trans, num_v);

	printf("\n");

	// getting the array of adjacency list from the transpose of the adjacency matrix
	struct adj_list** Adj_array_trans;
	Adj_array_trans = malloc(sizeof(struct adj_list*)*num_v);
	getAdjArray(Adj_array_trans, num_v, Adj_mat_trans);
	printAdjArray(Adj_array_trans, num_v, v_names);

	printf("\n");

	// sorting the indices of vertices according to the finishi time
	int* finishTimes = malloc(sizeof(int)*num_v);
	int* originalIndex = malloc(sizeof(int)*num_v);
	for (int i = 0; i<num_v; i++)
	{
		finishTimes[i] = DFS_vertices[i].fin;
		originalIndex[i] = i;
	}
	//quickSort(finishTimes, 0, num_v - 1, originalIndex);
    heapSort(finishTimes, 0, num_v - 1, originalIndex); //youngeun
    //countingSort(finishTimes, 0, num_v - 1, originalIndex); //youngeun
	// try heapSort() instead of quickSort
	// try countingSort() instead of quickSort
	for (int i = 0; i<(int)((double)num_v / 2.0); i++)
	{
		swapItems(originalIndex, i, num_v - i - 1); // to get the reversed order
	}

	// Second DFS
	struct DFS_vertex* DFS_vertices_second;
	DFS_vertices_second = malloc(sizeof(struct DFS_vertex)*num_v);
	source_index = originalIndex[0];
	printf("With %c as the source vertex,\n", v_names[source_index]);
	DFS(Adj_array_trans, num_v, DFS_vertices_second, source_index, originalIndex);
	printDFS(DFS_vertices_second, num_v, v_names);

	printf("\n");

	// Find SSCs
	int** SSCs = malloc(sizeof(int*)*num_v);
	for (int i = 0; i<num_v; i++) SSCs[i] = malloc(sizeof(int)*num_v);
	int num_SSCs = findSSCs(DFS_vertices_second, num_v, SSCs);
	printSSCs(SSCs, num_SSCs, v_names);
    
    //youngeun
    clock_t start, end;
    printf("\n<Dijkstra algorithm>\n");
    int** new_mat = malloc(sizeof(int*)*num_v);
    for (int i = 0; i<num_v; i++) new_mat[i] = malloc(sizeof(int)*num_v);
    
    start = clock();
    dijkstra(new_mat ,Adj_mat, num_v);
    end= clock();
    new_printMat(new_mat, num_v, v_names);
    printf("Dijkstra algorithm: %.4lf 초\n", (end-start)/(double)1000);
    
    //youngeun
    printf("\n<Bellman-Ford algorithm>\n");
    new_mat = malloc(sizeof(int*)*num_v);
    for (int i = 0; i<num_v; i++) new_mat[i] = malloc(sizeof(int)*num_v);
    
    start = clock();
    bellmanford(new_mat, Adj_mat, num_v);
    end = clock();
    new_printMat(new_mat, num_v, v_names);
    printf("Bellman-Ford algorithm: %.4lf 초\n", (end-start)/(double)1000);
    
    //youngeun
    printf("\n<Floyd-Warshall algorithm>\n");
    new_mat = malloc(sizeof(int*)*num_v);
    for (int i = 0; i<num_v; i++) new_mat[i] = malloc(sizeof(int)*num_v);
    
    start = clock();
    floydwarshall(new_mat, Adj_mat, num_v);
    end = clock();
    new_printMat(new_mat, num_v, v_names);
    printf("Floyd-Warshall algorithm: %.4lf 초\n", (end-start)/(double)1000);

	// freeing allocated memories
	for (int i = 0; i<num_v; i++)
	{
		free(Adj_mat[i]);
		free(Adj_mat_trans[i]);
		free(SSCs[i]);
	}
	free(Adj_mat);
	free(all_lists);
	free(Adj_array);
	free(DFS_vertices);
	free(Adj_mat_trans);
	free(Adj_array_trans);
	free(finishTimes);
	free(originalIndex);
	free(DFS_vertices_second);
	free(SSCs);

	return 0;
}

int getVeticesNames(char** v_names, char* filename)
{
	FILE* fp;
	fp = fopen(filename, "r");

	int num_v = 0;

	char buffer[256];
	fgets(buffer, 256, fp); // first line
	for (int i = 0; i<strlen(buffer) - 1; i++)
	{
		if (buffer[i] != '\t')
		{
			num_v++;
			*v_names = realloc(*v_names, sizeof(char)*num_v);
			(*v_names)[num_v - 1] = buffer[i];
		}
	}

	fclose(fp);

	return num_v;
}

void getAdjMat(int** Adj_mat, int num_v, char* filename)
{
	FILE* fp;
	fp = fopen(filename, "r");

	char buffer[256];
	fgets(buffer, 256, fp); // skip the first line

	for (int i = 0; i<num_v; i++)
	{
		for (int j = -1; j<num_v; j++)
		{
			fscanf(fp, "%s", buffer);
			if ((j != -1) && (strcmp(buffer, "INF") != 0))
				Adj_mat[i][j] = atoi(buffer);
			else if (strcmp(buffer, "INF") == 0)
				Adj_mat[i][j] = INFINITY; // just setting the maximum value of integer
		}
	}

	fclose(fp);
	return;
}

void printMat(int** Adj_mat, int num_v)
{
	for (int i = 0; i<num_v; i++)
	{
		for (int j = 0; j<num_v; j++)
		{
			if (Adj_mat[i][j] == INFINITY)
				printf("INF\t");
			else
				printf("%d\t", Adj_mat[i][j]);
		}
		printf("\n");
	}
}

void new_printMat(int** new_mat, int num_v, char* v_names) //youngeun
{
    for(int i =0; i< num_v; i++)
        printf("\t%c", v_names[i]);
    printf("\n");
    
    for(int i =0; i< num_v; i++){
        printf("%c\t", v_names[i]);
        for (int j = 0; j<num_v; j++)
        {
            if (new_mat[i][j] >= 9997)
                printf("INF\t");
            else
                printf("%d\t", new_mat[i][j]);
        }
        printf("\n");
    }
}

void getAdjArray(struct adj_list** Adj_array, int num_v, int** Adj_mat)
{
	for (int i = 0; i<num_v; i++)
	{
		num_all_lists++;
		if (num_all_lists > max_num_lists) exit(1);
		all_lists[num_all_lists - 1].v_index = i;
		all_lists[num_all_lists - 1].w = 0;
		all_lists[num_all_lists - 1].next = NULL;
		Adj_array[i] = &all_lists[num_all_lists - 1];

		for (int j = 0; j<num_v; j++)
		{
			if (i != j)
			{
				if (Adj_mat[i][j] != INFINITY)
				{
					num_all_lists++;
					if (num_all_lists > max_num_lists) exit(1);
					all_lists[num_all_lists - 1].v_index = j;
					all_lists[num_all_lists - 1].w = Adj_mat[i][j];
					all_lists[num_all_lists - 1].next = NULL;
					all_lists[num_all_lists - 2].next = &all_lists[num_all_lists - 1];
				}
			}
		}
	}

	return;
}

void printAdjArray(struct adj_list** Adj_array, int num_v, char* v_names)
{
	// print the array of adjacency list (Adj_array)
	printf("The array of adjacency list\n");

	for (int i = 0; i<num_v; i++)
	{
		printf("%c: ", v_names[i]);
		struct adj_list* cur_list;
		cur_list = Adj_array[i]->next;
		while (cur_list != NULL)
		{
			printf("%c, %d; ", v_names[cur_list->v_index], cur_list->w);
			cur_list = cur_list->next;
		}
		printf("\n");
	}
}

int DFS_VISIT(struct adj_list** Adj_array, int num_v, struct DFS_vertex* DFS_vertices, int u, int time)
{
	// change color to gray
    DFS_vertices[u].color = GRAY; //youngeun
    
	time++;
	// update discovery time
    DFS_vertices[u].dis = time; //youngeun
    
	struct adj_list* cur_list = Adj_array[u]->next;
	while (cur_list != NULL)
	{
		if (DFS_vertices[cur_list->v_index].color == WHITE)
		{
			DFS_vertices[cur_list->v_index].pi = u;
			time = DFS_VISIT(Adj_array, num_v, DFS_vertices, cur_list->v_index, time);
		}
		cur_list = cur_list->next;
	}
	// change color to black
    DFS_vertices[u].color= BLACK; //youngeun
    
	time++;
	// update finish time
    DFS_vertices[u].fin = time; //youngeun

	return time;
}

void DFS(struct adj_list** Adj_array, int num_v, struct DFS_vertex* DFS_vertices, int source_index, int* indexOrder)
{
	for (int u = 0; u<num_v; u++)
	{
		//each vertex's color has to be initialized to WHITE;
        DFS_vertices[u].color = WHITE; //youngeun
		DFS_vertices[u].pi = -1; //NULL;
	}
	int time = 0;
    
	// do DFS_VISIT for the source vertex first
	time = DFS_VISIT(Adj_array, num_v, DFS_vertices, source_index, time);

	for (int i = 0; i<num_v; i++)
	{
		if (indexOrder == NULL)
		{
			//if i-th vertex' color is white
            if(DFS_vertices[i].color == WHITE) //youngeun
				time = DFS_VISIT(Adj_array, num_v, DFS_vertices, i, time);
		}
		else
		{
			//if (indexOrder[i])-th vertex' color is white
            if(DFS_vertices[indexOrder[i]].color == WHITE) //youngeun
				time = DFS_VISIT(Adj_array, num_v, DFS_vertices, indexOrder[i], time);
		}

	}
}

void printDFS(struct DFS_vertex* DFS_vertices, int num_v, char* v_names)
{
	printf("DFS result\n");
	for (int u = 0; u<num_v; u++)
	{
		if (DFS_vertices[u].pi == -1)
			printf("%c: d=%d, f=%d, pi=%d(root)\n", v_names[u], DFS_vertices[u].dis, DFS_vertices[u].fin, DFS_vertices[u].pi);
		else
			printf("%c: d=%d, f=%d, pi=%d(%c)\n", v_names[u], DFS_vertices[u].dis, DFS_vertices[u].fin, DFS_vertices[u].pi, v_names[DFS_vertices[u].pi]);
	}
}

void getTransAdjMat(int** Adj_mat, int** Adj_mat_trans, int num_v)
{
	for (int i = 0; i<num_v; i++)
	{
		for (int j = 0; j<num_v; j++)
		{
			Adj_mat_trans[i][j] = Adj_mat[j][i];
		}
	}

	return;
}

void swapItems(int* Array, int i, int j)
{
	int temp = Array[i];
	Array[i] = Array[j];
	Array[j] = temp;
}

int partition(int* Array, int p, int r, int* originalIndex)
{
	int x = Array[r];
	int i = p - 1;
	for (int j = p; j<r; j++)
	{
		if (Array[j] <= x)
		{
			i++;
			swapItems(Array, i, j);
			swapItems(originalIndex, i, j);
		}
	}
	swapItems(Array, i + 1, r);
	swapItems(originalIndex, i + 1, r);

	return i + 1;
}

void quickSort(int* Array, int p, int r, int* originalIndex)
{
	if (p < r)
	{
		int q = partition(Array, p, r, originalIndex);
		quickSort(Array, p, q - 1, originalIndex);
		quickSort(Array, q + 1, r, originalIndex);
	}
}

int findDecendants(struct DFS_vertex* DFS_vertices, int num_v, int rootIndex, int*decendants)
{
	int num_decendants = 0;

	for (int u = 0; u<num_v; u++)
	{
		if ((DFS_vertices[u].dis > DFS_vertices[rootIndex].dis) && (DFS_vertices[u].fin < DFS_vertices[rootIndex].fin))
		{
			num_decendants++;
			decendants[num_decendants - 1] = u;
		}
	}

	return num_decendants;
}

int findSSCs(struct DFS_vertex* DFS_vertices, int num_v, int**SSCs)
{
	// assuming SSCs has (num_v x num_v) dimension (maximum dimension)
	// use findDecendants function
	// how to indicate there is no further element? (hint: look at printSSCs() function)
    
    int num_SSCs = 0;
    
    //youngeun
    for(int u =0 ;u < num_v ; u++){
        if(DFS_vertices[u].pi == -1){
            int* decendants = malloc(sizeof(int)*num_v);
            int num_decendants = findDecendants(DFS_vertices, num_v, u, decendants);
            SSCs[num_SSCs][0] =u;
            
            int j =0;
            for(j =0; j<num_decendants ; j++)
                SSCs[num_SSCs][j+1]= decendants[j];
            
            SSCs[num_SSCs][j+1]= -1;
            
            //printf("SSCs 갯수: %d\n", num_SSCs);
            //printf("root의 index: %d\n", u);
            //printf("decendats의 갯수는?: %d\n", num_decendants);
            num_SSCs++;
        }
    }
    //youngeun
    
    //int* decendants = malloc(sizeof(int)*num_v);
    //int num_decendants = findDecendants(DFS_vertices, num_v, u, decendants);

	return num_SSCs;
}

void printSSCs(int** SSCs, int num_SSCs, char* v_names)
{
	printf("SSCs\n");
	for (int i = 0; i<num_SSCs; i++)
	{
		int index = 0;
		printf("SSC%d: ", i + 1);
		while (SSCs[i][index] >= 0) // if it meats negative value, it means there is no further element.
		{
			printf("%c, ", v_names[SSCs[i][index]]);
			index++;
		}
		printf("\n");
	}
}

void heapSort(int* Array, int p, int r, int* originalIndex){
    printf("Implement Heap-sort function\n\n");
    int size = r-p+1;
    buildMaxHeap(Array, size, originalIndex);
    
    for(int i= size; i >= 1 ; i--){
        swapItems(Array, 0, size-1);
        swapItems(originalIndex, 0, size-1);
        size--;
        heapify(Array, size, 0, originalIndex);
    }
    
    /*printf("sorting 한 후 finish time: ");
    for(int i=0; i< 7; i++){
        printf("%d ", Array[i]);
    }
    printf("\n\n");*/
}

void buildMaxHeap(int* Array, int size, int* originalIndex){
    for(int mid = size/2-1; mid>=0; mid--){
        heapify(Array, size, mid, originalIndex);
    }
}

void heapify(int* Array, int size, int mid, int* originalIndex){
    int parent_node = mid;
    int left_node = parent_node*2+1;
    int right_node = parent_node*2+2;
    int largest_node = parent_node;
    int temp;
    
    if(left_node < size && Array[left_node] > Array[largest_node]){
        largest_node = left_node;
    }
    
    if(right_node < size && Array[right_node] > Array[largest_node]){
        largest_node = right_node;
    }
    if(parent_node != largest_node){
        swapItems(Array, largest_node, parent_node);
        swapItems(originalIndex, largest_node, parent_node);
        heapify(Array, size, largest_node, originalIndex);
    }
}

void countingSort(int* Array, int p, int r, int* originalIndex){
    printf("Implement Counting-sort function\n\n");
    int size = r-p+1;
    int B[size];
    int B2[size];
    
    int k =0;
    for(int i=0; i< size; i++){
        if(Array[i]> k){
            k = Array[i];
        }
    }
   
    int C[k+1];
    for(int i=0; i<= k; i++)
        C[i] = 0;

    for(int j=0; j<size; j++){
        C[Array[j]]++;
    }
    
    for(int i=1; i<=k; i++){
        C[i] = C[i] + C[i-1];
    }
    
    for(int j=size-1; j >=0; j--){
        B[C[Array[j]]-1] = Array[j];
        B2[C[Array[j]]-1] = originalIndex[j];
        C[Array[j]]--;
    }

    for(int j=0; j<size; j++){
        Array[j]= B[j];
        originalIndex[j] = B2[j];
    }
    
    /*printf("sorting 한 후 Array의 finish time: ");
    for(int i=0; i< size; i++){
        printf("%d ", Array[i]);
    }
    printf("\n\n");
    
    printf("sorting 한 후 origanalInex의 finish time: ");
    for(int i=0; i< size; i++){
        printf("%d ", originalIndex[i]);
    }
    printf("\n\n");*/
}

//dijkstra
void dijkstra(int** new_mat ,int** Adj_mat, int num_v){
    int cost[num_v][num_v], visited[num_v], nextnode;
    
    for(int k=0;k<num_v;k++){
        for(int i=0;i<num_v;i++)
            for(int j=0;j<num_v;j++)
                if(Adj_mat[i][j]==0)
                    cost[i][j]=0;
                else
                    cost[i][j]=Adj_mat[i][j];
        
        for(int i=0;i<num_v;i++){
            new_mat[k][i]=cost[k][i];
            visited[i]=0;
        }
        
        new_mat[k][k]=0;
        visited[k]=1;
        int count=1;
        
        while(count<num_v-1){
            int mindistance=INFINITY;
            for(int i=0;i<num_v;i++)
                if(new_mat[k][i] < mindistance && !visited[i]){
                    mindistance=new_mat[k][i];
                    nextnode=i;
                }
            
            visited[nextnode]=1;
            for(int i=0;i<num_v;i++)
                if(!visited[i])
                    if(mindistance+cost[nextnode][i]<new_mat[k][i]){
                        new_mat[k][i]=mindistance+cost[nextnode][i];
                    }
            count++;
        }
    }
}

//bellmanford
bool bellmanford(int** new_mat, int** Adj_mat, int num_v){
    for(int first = 0; first < num_v; first++){
        for(int second =0; second < num_v; second++){ //initialize
            if(second == first)  new_mat[first][second]= 0;
            else new_mat[first][second]= INFINITY;
        }
        
        for(int i = 1; i < num_v; i++){
            for(int j = 0; j < num_v; j++){
                for(int k = 0; k < num_v; k++){ //relaxation
                    if(Adj_mat[j][k] != INFINITY && new_mat[first][j] != INFINITY){
                        if(new_mat[first][k] > new_mat[first][j] + Adj_mat[j][k])
                            new_mat[first][k] = new_mat[first][j] + Adj_mat[j][k];
                    }
                }
            }
        }
        
        for(int j= 0; j < num_v; j++){
            for(int k = 0; k < num_v; k++){ //relaxation check
                if(Adj_mat[j][k] != INFINITY && new_mat[first][j] != INFINITY){
                    if(new_mat[first][k] > new_mat[first][j] + Adj_mat[j][k])
                        return false; //negative-weight cycle 존재하는 경우
                }
            }
        }
        
    }
    return true; //no such cycle인 경우
    
    /*for(int i=0; i< num_v; i++){
        for(int j=0; j< num_v; j++)
            printf("%12d ", new_mat[i][j]);
        printf("\n");
    }
    printMat(new_mat, num_v);*/
}

//floydwarshall
void floydwarshall(int** new_mat, int** Adj_mat, int num_v){
    for(int i=0; i< num_v; i++){ //initialize
        for(int j=0; j< num_v; j++){
            if(i == j) new_mat[i][j]=0;
            else new_mat[i][j] = Adj_mat[i][j];
        }
    }
    
    for(int k=0; k <num_v; k++){ //dynamic programming
        for(int i=0; i< num_v; i++){
            for(int j=0; j< num_v; j++){
                    if(new_mat[i][k] + new_mat[k][j] < new_mat[i][j] && new_mat[i][k] != INFINITY && new_mat[k][j] != INFINITY)
                        new_mat[i][j]= new_mat[i][k] + new_mat[k][j];
            }
        }
    }
    
    /*for(int i=0; i< num_v; i++){
        for(int j=0; j< num_v; j++)
            printf("%12d ", new_mat[i][j]);
        printf("\n");
    }
    printMat(new_mat, num_v);*/
}
