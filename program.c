

#include <stdio.h>
#include <stdlib.h>

#define FILENAME "DirectedConnectedGraph.txt" // <- ***NAME OF FILE***

// The file containing the edges of the graph should either be named DirectedGraph.txt, or the symbolic constant above should be
// changed to the name of the file being used. 

// This program runs in O(V+E) time, which scales linearly. Please see the detailed description above each function explaining why. (as well as in-line comments)
// Basically, any nested for/while loops only iterate proportionally to O(V+E) or O(V+V), which equals O(V). DFS is performed with a counter
// to stop once every vertex has been explored, & the out degree sequence is found using three arrays in addition to the adjacency list
// array & the array of Vertex pointers. The first array is initialiized with the vertices adjacent to each SCC. Then, using this array,
// the 2nd array is initalized by adding the current SCC to the linked list associated with each adjacent SCC, while incrementing a counter
// for the current SCC to find the out degree. The actual linked lists created in this process are not used any further since the graph is
// a directed graph, & this process adds to each SCC the SCCs that are pointing to it rather than vice versa, but the counters still function
// correctly. The next step is merely implemented to output the vertices contained in each SCC when printing the out-degrees. This is simply done 
// by creating another array of Node* & iterating once through the array of Vertex*, whilst adding each vertex to the SCC it is contained in.
// Finally, the first of the 3 arrays created is reused in order to order the SCCs based on their out-degrees by having the index of the array
// represent the out-degree and the value of each node added to the array be a SCC #.  

// Note that when implementing DFS counters are used for the time as well as counting the number of vertices explored as well
// as incrementing for the strongly-connect-components, and each of these are somewhat hidden inside of a Vertex struct at index
// 0 (vertexArray[0] or *vertexArray) of the array of vertices; this was done in order to lesson the number of parameters used for
// depth-first-search. Also note that for each dummy node created in the arrays of adjacency lists, a counter is hidden in front
// of the dummy Node * struct (method taught in 2050) to be used to keep track of the number of adjacent vertices of each vertex. 
// This wasn't an incredibly necessary addition, but I felt that it made things easier in DFS.

// struct 1: used to represent the graph
typedef struct node{
	int index;
	struct node *next;
} Node;

// struct 2: used to store information for each vertex in the graph; is not connected at all to Node* but the structs are used alongside each other
typedef struct vertex{
	int color;   //0 = white, 1 = black the first time DFS is performed, vice versa the 2nd time DFS is performed
	int discovery;
	int finish;
	int SCC;
} Vertex;


// Function Prototypes:
void output(char *fileName);
Node **readFile(char *fileName, int vertices, Node ***altList);
Node **createList(int vertices);
int addAdjacent(Node** list, int vertex, int adjacent);
Node* createNode(int vertex);
void freeList(Node **list, int vertices);
void freeListSCC(Node** array, int numOfSCC);
Vertex **initVertexArray(int vertices);
int outputVertices(char *fileName);
void DFS(Vertex** vertexArray, Node** adjArray, int n, int *finishingOrder);
void SCC(Vertex** vertexArray, Node** adjArray, int n);
void outDegreeSequence(Node** list, Vertex** vertexArray, int vertices);


// In order to create maximum data encapsulation, the user just has one function to call which will output the product
// FILENAME is the symbollic constant defined for the name of the file
int main(void) {
	output(FILENAME);
}

// Function output() takes the name of the file containing edges of a graph, & outputs the result
/* 
	Parameters:
		char *fileName is the name of a text file containing the edges of a directed graph
*/
// EFFICIENCY: should be O(V+E) as readFile() takes O(V+E) time, iterating through the file takes O(E) time, freeing takes
// O(V+E) time, & depth-first-search takes O(V+E) time both times it is performed
// so total time is dependant soley on the number of edges and vertices, making time linear
void output(char *fileName){
	puts("");
	
	// O(V)
	int vertices = outputVertices(fileName); // creates an integer variable & stores number of vertices in it by scanning the file
	if (!vertices){return;}
	
	// O(V)
	Vertex ** vertexArray = initVertexArray(vertices); // initVertexArray returns an array of Vertex* using the number of vertices
	if (vertexArray == NULL){return;}				   // NOTE: vertexArray is indexed starting at 1 instead of 0 so the 0 index can
													   // hold a counter, a counter for time, & a counter for the SCCs
	// O(V+E)
	Node** altList; // altList is a 2nd list used as the transpose of the graph
	Node** list = readFile(fileName, vertices, &altList); // readFile initializes both list (the adjacency list for the graph)  
	if (list == NULL || altList == NULL){		  // & altList (the adjacency list for the transpose of the graph)
		puts("\n\n\n  Program author: Noah Free, nsfq94\n\n");
		return;
	}

	// O(V)
	int* finishingOrder = malloc(sizeof(int)*(2*vertices)); // finishingOrder is an array used to find finishing order of nodes in DFS
	for(int i = 0; i < 2*vertices; i++){					// the size of finishingOrder is equal to the time it will take for DFS (2*vertices)
		finishingOrder[i] = -1; // each index in finishingOrder is initialized to -1
	}
	
	// O(V+E) since the for loop ends when DFS has been called for each vertex
	for (int i = 1; ((*vertexArray)->finish != vertices) && (vertexArray[i]->color != 1); i++){ 
		DFS(vertexArray, list, i, finishingOrder); // DFS called on each vertex; since DFS is recursive,DFS will only be called once for
	}											   // each vertex; this is done by using the counter (*vertexArray)->finish, incremented for
												   // each call of DFS
	
	// the first index of vertexArray (vertexArray[0]) is used as a counter, a counter for time, and a counter for the SCCs
	(*vertexArray)->finish = 0;	// counter reset to 0 before running DFS on the transpose using SCC() function
	(*vertexArray)->SCC = 0;
	
	// O(V+E) for the same reason as above
	for (int i = (2*vertices - 1); (*vertexArray)->finish != vertices; i--){
		if (finishingOrder[i] != -1){
			if (vertexArray[finishingOrder[i]]->color == 1){ // now color 1 is white & color 0 is black rather than vice versa as before
				SCC(vertexArray, altList, finishingOrder[i]); // SCC, which is a modified DFS function, is recursive & will be called once for each vertex
				(*vertexArray)->SCC += 1; // when SCC returns, the SCC counter increments since SCC has already been called for a strongly-connect-component
			}
		}
	}
	
	// O(V+E)
	outDegreeSequence(list, vertexArray, vertices); // outDegreeSequence takes the list, the vertexArray, & the # of vertices & prints the outDegreeSequence for SCCs
	
	free(finishingOrder); // only 1 call to malloc() was used for finishingOrder, so only 1 call to free() is required
	// O(V)
	for (int i = 1; i <= vertices; i++){ // for loop is used to free memory allocated for vertexArray; print statement can be uncommented to print info for each vertex
		// printf("\n  Vertex: %d\n  Discovery: %d\n  Finish: %d\n  Strongly-Connected-Component: %d\n", i, vertexArray[i]->discovery, vertexArray[i]->finish, vertexArray[i]->SCC);
		free(vertexArray[i]);
	}
	free((*vertexArray)); // the only reason the [0] index is freed after the rest is if the print statement in the for loop above is uncommented
	free(vertexArray); 
	
	// O(V+E)
	freeList(list, vertices); // function freeList is used to free all memory allocated for list, & then called again below for altList (the transpose)
	// O(V+E)
	freeList(altList, vertices);
	list = NULL;
	
	puts("\n\n\n  Program author: Noah Free, nsfq94\n\n");
}

// Function createList() takes the number of vertices and initializes a list of dummy nodes for the list by calling
// the function createNode(), which runs in O(1) time; the created list is returns unless malloc() fails, in which
// case NULL will be returned & an error message printed
/*
	Parameters:
		int vertices is the number of vertices in a graph
*/
// EFFICIENCY: O(V) since the number of elements initialized in the array depends on the number of vertices
Node** createList(int vertices){
	Node** list = malloc(sizeof(Node *)*vertices);
	if (list == NULL){
		puts("\n  Error: Unable to allocate memory.\n");
		return NULL;
	}
	else {
		// O(V) time
		for (int i = 0; i < vertices; i++){
			list[i] = createNode(i+1); // createList() calls createNode() for each vertex in the malloced list
			if (list[i] == NULL){
				for (int j = 0; j < i; j++){
					free(((int *)list[j] - 1)); // memory is freed if malloc every fails; cast to (int *) because an integer counter is hidden
				}								// in front of each dummy node using a method taught in 2050
				free(list);
				puts("\n  Error: Unable to allocate memory.\n");
				return NULL;
			}
		}
	}
	return list;
}


// Function readFile() reads the inputted file & creates a list of vertices & their adjacent vertices which is returned to the calling program
// The transpose of the first list/graph is also created and set equal to altList
// If malloc ever fails during the function all memory previously allocated is freed & NULL is returned.
/* 
	Parameters:
		char *fileName is the name of a text file that contains pairs of vertices separated by a space that represents an edge on a graph
		int vertices is the number of vertices in the graph, also the largest number in the input file
		Node ***altList is a pointer to an array of adjacency lists that is used to return a 2nd list
*/
// EFFICIENCY: O(V+E) since the function first calls createList() twice, which each take O(V) time; it then iterates through the file, 
// which takes time dependant on number of edges, while adding each edge to both lists/arrays in O(1) time
Node** readFile(char *fileName, int vertices, Node ***altList){
	FILE* fPtr = fopen(fileName, "r");
	if (!fPtr){
		puts("\n  Error: File is not valid.\n");
		return NULL;
	}
	
	Node** list = createList(vertices); // function createList is called to create an empty list for the graph
	if (list == NULL){
		fclose(fPtr);
		return NULL;
	}
	Node** tempList = createList(vertices); // function createList is called to create an empty list for the transpose
	if (tempList == NULL){
		free(list);
		fclose(fPtr);
		return NULL;
	}
	
	int v1, v2;
	while (!(feof(fPtr))){ // while loop iterates through the inputted file
		fscanf(fPtr, "%d %d\n", &v1, &v2);
		// if addAdjacent() returns 0 malloc() has failed
		if (!(addAdjacent(list, v1, v2))){ // addAdjacent is called in the if statement, & it returns 1 if successfull
			puts("\n  Error: Unable to allocate memory.\n");
			freeList(list, vertices); // function freeList is used to free memory if malloc fails even though lists are not complete
			freeList(tempList, vertices);
			fclose(fPtr);
			*altList = NULL;
			return NULL;
		}
		if (!(addAdjacent(tempList, v2, v1))){ // in the same while loop addAdjacent is called for the transpose, switching v1 & v2
			puts("\n  Error: Unable to allocate memory.\n");
			freeList(tempList, vertices);
			freeList(list, vertices);
			fclose(fPtr);
			*altList = NULL;
			return NULL;
		}
	}
	
	fclose(fPtr);
	*altList = tempList; // if malloc never fails, the address of altList is set equal to the address of tempList, & list is returned
	return list;
}


// Function addAdjacent() adds a vertex to the adjacency list of another vertex
// 1 is returned if the function is successful, 0 is returned if malloc() fails
/* 
	Parameters:
		Node** list is an array of linked lists for which the adjacent vertex needs to be added
		int vertex is the vertex that the adjacent vertex is adjacent to
		int adjacent is the adjacent vertex which is added to the adjacency list
*/
// EFFICIENCY: O(1) since the function time is not dependant on number of vertices or edges; it merely inserts an adjacent vertex at the inputted vertex/index
int addAdjacent(Node** list, int vertex, int adjacent){
	Node *new = malloc(sizeof(Node));
	if (new == NULL) return 0;
	else {
		(*(((int *)(list[vertex - 1])) - 1))++; // the counter hidden in front of each dummy node in the inputted list is incremented by 1
		new->next = list[(vertex - 1)]->next; // the created node is inserted into the corresponding linked list 
		list[(vertex - 1)]->next = new;
		new->index = adjacent; 
	}
	return 1;
}


// Function createNode() creates a dummy node that contains the inputted vertex & a counter for the out degree 
// The created node is returned unless malloc() fails, in which case NULL will be returned
/*
	Parameters:
		int vertex is the vertex that the node should hold
*/	
// EFFICIENCY: O(1) since this function is independant of the size of V or E; it merely takes an input & returns an output in constant time
// NOTE: this function uses a technique taught in CS2050 where an int is hidden in front of another data type, in this case a Node*
Node* createNode(int vertex){
	int *temp = malloc(sizeof(Node) + sizeof(int)); // an int* is malloced instead of a Node* in order to hide an integer counter in front of the Node
	if (temp != NULL){
		*temp = 0; // counter is set equal to 0
		Node *node = (Node *)(temp + 1); // the dummy node is set equal to the int* incremented by the size of an int
		node->index = vertex; // inputted vertex is put in the dummy node
		node->next = NULL;
		return node;
	}
	else{
		return NULL;
	}
}


// Function freeList() frees all memory alloquated for the inputted list containing adjacency lists
/*
	Parameters:
		Node** list is a linked list containing vertices & their adjacency lists
		int vertices is the number of vertices in the array
*/
// EFFICIENCY: O(V+E) where V = vertices & E = edges, since the function iterates through the array & frees the address
// of every adjacent element for each vertex, so time is proportional to the number of vertices plus the number of edges
void freeList(Node** list, int vertices){
	Node *tempNode, *freeNode; // 2 nodes are used in order to free all the memory allocated for the inputted list
	for (int i = 0; i < vertices; i++){ // is called once for each vertex in the list
		tempNode = list[i]->next; // tempNode is set equal to first node in linked list after the dummy node
		if (tempNode != NULL){ // if tempNode is equal to NULL the adjacency list for this vertex is empty
			while (tempNode->next != NULL){ // otherwise, free each node until tempNode->next equals NULL
				freeNode = tempNode;
				tempNode = tempNode->next;
				free(freeNode);
			}
			free(tempNode); // then free tempNode
		}
		free((((int *)(list[i])) - 1)); // because an integer was stored in front of each dummy node, the dummy nodes have to be freed in a special way
	}
	free(list); // finally, the memory allocated for the list itself is freed
}

// Function freeListSCC() frees all memory alloquated for the inputted SCC list
/*
	Parameters:
		Node** array is a linked list containing several linked lists
		int numOfSCC is the number of strongly connected components in the inputted array, also number of linked lists the array contains
*/
// EFFICIENCY: O(V) or O(V+E) since the function iterates through the array & frees all allocated memory, which is reliant on the number 
// of SCCs plus either the number of vertices or the number of edges, & how many SCCs there are depends on the number of vertices
void freeListSCC(Node** array, int numOfSCC){
	Node *tempNode, *freeNode;
	for (int i = 0; i < numOfSCC; i++){
		tempNode = array[i];
		while (tempNode != NULL){ // this function works the same as freeList but the dummy nodes do not hide an integer in front of them, so
			freeNode = tempNode;  // freeing the memory is a little more simple
			tempNode = tempNode->next;
			free(freeNode);
		}
	}
	free(array);
}

// Function initVertexArray() takes in the number of vertices and outputs an array of struct pointers, or NULL is malloc() fails
/*
	Parameters: 
		int vertices is the number of vertices in a graph
*/
// EFFICIENCY: O(V) time since it allocates memory for an array & then allocates memory for each 
// vertex in the graph
Vertex **initVertexArray(int vertices){
	Vertex ** vertexArray = malloc(sizeof(Vertex *)*(vertices + 1)); // memory for an extra space is allocated bebcause the first index is used for counters
	if (vertexArray != NULL){
		for (int i = 0; i <= vertices; i++){
		vertexArray[i] = malloc(sizeof(Vertex)); // for each vertex, allocate memory & then initialize the struct values 
			if (vertexArray[i] != NULL){
				vertexArray[i]->color = 0;
				vertexArray[i]->discovery = 0;
				vertexArray[i]->finish = 0;
				vertexArray[i]->SCC = -1;
			}
			else {
				for (int j = 1; j < i; j++){ // all memory is freed if malloc every fails
					free(vertexArray[j]);
				}
				puts("\n  Error: Unable to allocate memory.\n");
				return NULL;
			}
		}	
	}
	else {
		puts("\n  Error: Unable to allocate memory.\n");
		return NULL;
	}
		
	return vertexArray;
}

// Function outputVertices() takes the name of a file & outputs the largest vertex in the file, which is equal to the number of vertices
/*
	Parameters:
		char *fileName is the name of a file that contains the edges of a graph in "x y\n" format
*/
// EFFICIENCY: this function takes O(E) time since it scans the file and compares each vertex of each edge to the current largest vertex
int outputVertices(char *fileName){
	int vertices = 0; // vertices is used to find the largest value in the file
	int x, y; // x & y are used to retrieve values from file & then compare to vertices
	
	FILE* fPtr = fopen(fileName, "r");
	if (fPtr != NULL){
		// while loop iterates through the file, scanning for the largest element (O(E) time)
		while (!feof(fPtr)){ 
			fscanf(fPtr, "%d %d\n", &x, &y);
			if (vertices < x){ // vertices is set equal to x if x is larger
				vertices = x;
			}
			if (vertices < y){ // vertices is set equal to y if y is larger
				vertices = y;
			}
		}
		if (!vertices){ // if the number of vertices is still equal to 0, then the file is empty
			puts("\n  File empty: Input file does not contain any edges.\n");
		}
	}
	else{
		puts("\n  Error: Unable to open file.\n");
		return 0;
	}
	fclose(fPtr);
	return vertices;
}

// Function DFS() runs Depth-First-Search on the graph represented by the inputted arrays, starting with int n-1
/*
	Parameters:
		Vertex** vertexArray is an array of Vertex pointers (Vertex *)
		Node** adjArray is an array of adjacency lists for each vertex
		int n is the starting vertex for DFS
		int *finishingOrder is an integer array containing the list of vertices in order based on finishing time
*/
// Depth-First-Search is performed recursively, searching each unsearched adjacent vertex & then backtracking, all
// while labeling discovery & finishing times
void DFS(Vertex** vertexArray, Node** adjArray, int n, int *finishingOrder){
	Node *temp = adjArray[n-1]->next; // a temp node is set equal to the first node in n's adjacency list
	Vertex *vertex = vertexArray[n]; // vertexArray[n] corresponds to adjArray[n-1]
	(*vertexArray)->discovery += 1; // time is incremented by 1
	vertex->color = 1; // the color of the vertex being searched is changed from white to black (0 to 1)
	vertex->discovery = (*vertexArray)->discovery; // the discovery time of the vertex being searched is set equal to the time counter
	
	for (int i = 0; i < *(((int *)(adjArray[n-1]))-1); i++){ // for loop iterates for each vertex in a vertex's adjacency list, using the stored counter in front of the cummy node
		if (vertexArray[temp->index]->color == 0){ // DFS is only called on the vertex if its color is white/0 (it hasn't yet been searched
			DFS(vertexArray, adjArray, temp->index, finishingOrder);
		}
		temp = temp->next;
	}
	(*vertexArray)->discovery += 1; // time counter incremented by 1
	(*vertexArray)->finish += 1; // vertex counter is incremented by 1
	vertex->finish = (*vertexArray)->discovery; // the finishing time of the vertex is set equal to the time counter
	finishingOrder[(*vertexArray)->discovery - 1] = n; // the vertex's number is placed in the finishingOrder array based on its finishing time
}

// Function SCC() performs Depth-First-Search but does so on the transpose of the graph while
// finding each vertex's strongly-connect-component
/*
	Parameters:
		Vertex** vertexArray is an array of Vertex pointers (Vertex *)
		Node** adjArray is an array of adjacency lists for each vertex
		int n is the starting vertex for DFS
*/
// Depth-First-Search is performed recursively, searching each unsearched adjacent vertex & then backtracking, all
// while labeling each vertex's strongly-connect-component
void SCC(Vertex** vertexArray, Node** adjArray, int n){
	Node *temp = adjArray[n-1]->next; // a temp node is set equal to the first node in n's adjacency list
	Vertex *vertex = vertexArray[n]; // vertexArray[n] corresponds to adjArray[n-1]
	vertex->color = 0; // since the same vertexArray is used as before, the vertices' colors are represented by 1 = white/unsearched, 0 = black/searched
	
	for (int i = 0; i < *(((int *)(adjArray[n-1]))-1); i++){ // for loop iterates for each vertex in a vertex's adjacency list, using the stored counter in front of the cummy node
		if (vertexArray[temp->index]->color == 1){ // SCC is only called on the vertex if its color is white/1 (it hasn't yet been searched)
			SCC(vertexArray, adjArray, temp->index);
		}
		temp = temp->next;
	}
	(*vertexArray)->finish += 1; // the vertex counter is incremented by 1
	vertex->SCC = (*vertexArray)->SCC; // the vertex's SCC is set equal to the SCC counter
}

// Function outDegreeSequence() computes the out-degrees for the strongly connected components and prints them
// In this function, 3 different arrays are created and used in order to find the out degree of each SCC as well
// as print the vertices contained in each SCC; the first array is then reused to order the SCCs based on out-degrees
/*
	Parameters:
		Node **list is an array of linked lists representing adjacency lists for vertices in a graph
		Vertex **vertexArray is an array of Vertex pointers, each representing a vertex & containing that vertex's SCC
		int vertices is the number of vertices in the graph
*/
// EFFICIENCY: O(V+E) because a similar method to Machine Problem One is used in order to avoid quadratic time
// Note that this function frees ALL memory that is allocated in the function, both if malloc() fails and if it doesn't fail
void outDegreeSequence(Node** list, Vertex** vertexArray, int vertices){
	int boolean = 0; // this boolean variable is changed to 1 if memory allocation is every unsuccessful during this function
	Node *new, *temp; // new & temp are used throughout the function for traversing/freeing the linked lists
	Node **tempArray = malloc(sizeof(Node *)*((*vertexArray)->SCC)); // tempArray is the first of the arrays created in this function;
	if (tempArray == NULL) boolean = 1;								 // it is used temporarily to hold the adjacent vertices of each vertex 
																     // in each strongly-connect-component
																	 
	// O(V), proportional to number of SCCs which is proportional to number of vertices
	for (int i = 0; (i < ((*vertexArray)->SCC)) && !boolean; i++){ // tempArray is initialized with dummy nodes based on number of SCCs
		tempArray[i] = malloc(sizeof(Node));
		if (tempArray[i] != NULL){
			tempArray[i]->index = 0;
			tempArray[i]->next = NULL;
		}
		else {
			for (int j = 0; j < i; j++){ // all memory is freed if malloc fails, & boolean value is set to true, numbered 2 to show place of error
				free(tempArray[j]);
			}
			free(tempArray);
			boolean = 2;
		}
	}
	
	// O(V+E) since the array of adjency lists is iterated through completely, & it contains the vertices and their 
	// adjacent vertices, which depends on the number of edges in the graph
	for (int i = 0; (i < vertices) && !boolean; i++){ // for loop iterates for every vertex in the graph
		temp = list[i]->next;
		for (int j = 0; (j < *(((int *)(list[i]))-1)) && !boolean; j++){ // using the hidden counter in front of the dummy nodes in list,
																		 // the for loop iterates for every adjacent vertex of the current vertex 
			if (vertexArray[temp->index]->SCC != vertexArray[i+1]->SCC){ // if the SCC of the adjacent element does not equal the SCC of the current vertex
				new = malloc(sizeof(Node));								 // then the adjacent element is added to tempArray at the SCC's index
				if (new != NULL){
					new->index = temp->index;
					new->next = tempArray[vertexArray[i+1]->SCC]->next;
					tempArray[vertexArray[i+1]->SCC]->next = new;
				}
				else {
					freeListSCC(tempArray, (*vertexArray)->SCC); // function freeListSCC is called if malloc fails
					boolean = 3; // & boolean is set to true
				}
			}
			temp = temp->next; // temp is moved to the next pointer
		}
	}
	
	Node **array = malloc(sizeof(Node *)*((*vertexArray)->SCC)); // now a 2nd SCC array is allocated in order to find out-degrees for the SCCs
	if (array != NULL &&  boolean) free(array); // if boolean is true, then malloc has already failed somewhere, so array needs to be freed
	if (array == NULL && !boolean){
		freeListSCC(tempArray, (*vertexArray)->SCC); // if boolean is false & array equals NULL, then tempArray needs to be freed
		boolean = 4;
	}
	
	// O(V), proportional to number of SCCs which is proportional to number of vertices
	for (int i = 0; (i < ((*vertexArray)->SCC)) && !boolean; i++){ // this loop allocates memory for array[] based on the number of SCCs
		array[i] = malloc(sizeof(Node)); // dummy nodes are created 
		if (array[i] != NULL){
			array[i]->index = 0;
			array[i]->next = NULL;
		}
		else { // if malloc() fails, all memory allocated in this function is freed & boolean is set to 1
			for (int j = 0; j < i; j++){ 
				free(array[j]);
			}
			free(array);
			freeListSCC(tempArray, (*vertexArray)->SCC);
			boolean = 5;
		}
	}
	
	// O(V+E) since tempArray[] is iterated through completely, comparing the SCC of each vertex in tempArray[] to the current
	// SCC and then adding the current SCC onto the linked list for the vertex's SCC if it is not already contained in it
	for (int i = 0; (i < ((*vertexArray)->SCC)) && !boolean; i++){ // iterates once for each SCC
		temp = tempArray[i]->next; // temp is set equal to first adjacent element of tempArray[i]
		while ((temp != NULL) && !boolean){ // iterates until temp reaches the end of the linked list in each array index
			// basically, this if statement checks if the linked list for the SCC of the current adjacent element, & if it's not empty then
			// compares the SCC of the last added element in the SCC list of the adjacent element in tempArray to i, which is the current SCC
			if (array[vertexArray[temp->index]->SCC]->next == NULL || array[vertexArray[temp->index]->SCC]->next->index != i){
				new = malloc(sizeof(Node));
				if (new != NULL){
					new->next = array[vertexArray[temp->index]->SCC]->next; // then, the current SCC (i) is added to the linked list of the adjacent
					array[vertexArray[temp->index]->SCC]->next = new; 		// elements's SCC
					new->index = i;
					array[i]->index += 1; // the index of i, the current SCC, is incremented; this serves as a counter for the out-degree of the SCC
				}
				else {
					freeListSCC(array, (*vertexArray)->SCC); // all memory is freed if malloc fails
					freeListSCC(tempArray, (*vertexArray)->SCC);
					boolean = 6;
				}
			}
			temp = temp->next; // the temp pointer is moved to the next node
		}
	}
	
	Node **finalArray = malloc(sizeof(Node *)*((*vertexArray)->SCC)); // finalArray is used in order to print the vertices contained in each SCC
	if (finalArray != NULL && boolean) free(finalArray); // if boolean is true, then malloc has already failed somewhere, so array needs to be freed								  // (the out-degree of each SCC has already been found in the last loop)
	if (finalArray == NULL && !boolean){ // if finalArray is NULL, then ALL previously allocated memory in the function needs to be freed
		freeListSCC(tempArray, (*vertexArray)->SCC);
		freeListSCC(array, (*vertexArray)->SCC);
		boolean = 7;
	}
	
	// O(V), proportional to number of SCCs which is proportional to number of vertices
	for (int i = 0; (i < ((*vertexArray)->SCC)) && !boolean; i++){ 
		finalArray[i] = malloc(sizeof(Node)); // a dummy node is created for the head of each linked list in the array
		if (finalArray[i] != NULL){
			finalArray[i]->index = 0;
			finalArray[i]->next = NULL;
		}
		else {
			// once again, all memory is freed if malloc fails
			for (int j = 0; j < i; j++){
				free(finalArray[j]);
			}
			free(finalArray);
			freeListSCC(tempArray, (*vertexArray)->SCC);
			freeListSCC(array, (*vertexArray)->SCC);
			boolean = 8;
		}
	}
	
	// O(V) since each vertex in the graph is added to the linked list associated with the vertex's SCC
	for (int i = vertices; (i >= 1) && !boolean; i--){ // this loop adds each vertex in the graph to the linked list associated with the vertex's
		temp = malloc(sizeof(Node));				   // SCC in the newly created finalArray, starting with the highest numbered vertex
		if (temp != NULL){
			temp->index = i; // the index of the newly created node is set equal to i (a vertex)
			temp->next = finalArray[vertexArray[i]->SCC]->next; // the node is then added into the linked list associate with the vertex's SCC,
			finalArray[vertexArray[i]->SCC]->next = temp;		// which corresponds to finalArray[vertexArray[i]->SCC]->next
		}
		else {
			freeListSCC(finalArray, (*vertexArray)->SCC); // ALL memory is freed if malloc every fails
			freeListSCC(array, (*vertexArray)->SCC);
			freeListSCC(tempArray, (*vertexArray)->SCC);
			boolean = 9;
		}
	}
	
	// O(V+E) since this for loop is iterating through tempArray, which contains all vertices & their adjacent vertices, to free
	// the vertices in order to reuse the array
	for (int i = 0; (i < (*vertexArray)->SCC) && !boolean; i++){ // iterates based on the number of SCCs
		temp = tempArray[i]->next;
		while (temp != NULL){ // loop is used to free every node in the list besides the dummy nodes (since the list is about to be reused)
			new = temp; // variables new & temp are used to traverse each list and free each node
			temp = temp->next;
			free(new);
		}
		tempArray[i]->next = NULL;
	}
	
	// O(V) since this for loop runs proportionally to the number of SCC which is proportional to the number of vertices
	// This for loop is used to order the SCCs based on their out-degrees
	for (int i = 0; (i < (*vertexArray)->SCC) && !boolean; i++){ // loop iterates for every SCC 
		temp = malloc(sizeof(Node));
		if (temp != NULL){
			temp->index = i;
			// the method used here is similar to that described in the solution to Machine Assignment 1; in order to order the SCC's out-degrees,
			// the list of SCCs containing their out-degrees is iterated through, and each SCC number is placed in a node and added to the linked
			// list that corresponds to that node's out degree, so an SCC with out degree 0 will be added to tempArray at index 0, an SCC with out
			// degree 1 will be added to tempArray at index 1, and so on. If more than 1 SCC have same out degree, they will be placed next to each
			// other in the same linked list
			temp->next = tempArray[array[i]->index]->next;
			tempArray[array[i]->index]->next = temp;
		}
		else {
			freeListSCC(finalArray, (*vertexArray)->SCC); // ALL memory allocated in this function is freed if malloc fails
			freeListSCC(array, (*vertexArray)->SCC);
			freeListSCC(tempArray, (*vertexArray)->SCC);
			boolean = 10;
		}
	}
	
	// O(V) since for each SCC in the graph (which is dependant on the number of vertices) it's out degree is printed
	// along with the vertices contained in it, which technically takes O(V+V) of O(2V) time, but that scales to O(V)
	for (int i = 0; (i < ((*vertexArray)->SCC)) && !boolean; i++){ // loop iterates based on the number of SCCs
		new = tempArray[i]->next;
		while (new != NULL){ // this while loop goes through tempArray & prints the SCC for each out degree
			temp = finalArray[new->index]->next;
			printf("\n  Strongly Connected Component %d Out-Degree:  %d\n  (SCC %d represents vertices ", (new->index + 1), array[new->index]->index, (new->index + 1));
			while (temp != NULL){ // this while loop goes through finalArray & prints the vertices contained in each SCC
				if (temp->next != NULL){
					printf("%d, ", temp->index);
				}
				else {
					printf("%d)\n", temp->index);
				}
				temp = temp->next; // temp is then reassigned to the next node (or NULL) in finalArray
			}
			new = new->next; // new is then reassigned to the next node (or NULL) in tempArray
		}
	}
	
	if (!boolean){ // if boolean is equal to 1-10, then memory allocated failed somewhere in the function, but all memory has 
				   // already been freed in that case; otherwise, boolean wouuld still equal 1, & tempArray, array, & vertexArray
				   // are done being used, so freeListSCC() is called on each of them to free memory allocated
		// O(V+E)
		freeListSCC(tempArray, ((*vertexArray)->SCC));
		// O(V)
		freeListSCC(array, ((*vertexArray)->SCC));
		// O(V)
		freeListSCC(finalArray, ((*vertexArray)->SCC));
	}
	else {
		printf("\n  Error: Unable to allocate memory for SCC list. (%d)\n", boolean); // a message is printed if memory allocation failed
	}
}