/*Implementation of
A Simple Linear 3-Edge-Connected Components Algorithm with reduction
(Standard Version)

Author:
Yung H. Tsin
School of Computer Science
University of Windsor

Publication:
Theory of Computing Systems

Publication Date:
2005

Implemented in C by:
Nima Norouzi
School of Computer Sience
University of Windsor
2006
*/

/*This program compute all 3edge connected components of a biconnected graph.
  You need to make an executable program file after compiling this program by
  using gcc compiler or another C compiler. A desired biconnected graph must
	be ready in a text file having the number of vertices at the first line
	and the adjacency lists of vertices in the following; each edge is determined
	by ">", and each list starts at a new	line. The name of the file has to be
	passed to the executable file as the first argument; e.g "./3edge graph.txt"
	Content example for an acceptable text file(i.e. graph.txt):
5
1>2>5>2
2>5>4>4>3>3>1>1
3>2>2
4>2>2
5>2>1
*/

/*Use this command "ulimit -s" to set an appropriate stack size for
  your Linux in order to NOT encounter stack overflows during recursive
	calls at run times.
*/

/*Use command "time" before any execution command to get an accurate
  running time; e.g. "time ./3edge graph.txt"
*/

/*To display connected components, remove block comments(*/ /*) from lines
  containing "//PRINT".
*/

/*To have the algorithm say only the input graph is 3-edge-connected or NOT,
  remove block comments(*/ /*) from lines containing "//YesOrNo".
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <inttypes.h>

#include "commonC.h"

typedef int int32_t; //ensures we get big enough ints.

#define max_num 2147483647
/*Since we used "int" variables to store vertices, the lergest positive number
can be represented by a 32bit "int" data type variable is 2147483647.
*/
#define max_length 10
/*the number 2147483647 has 10 digits.
*/

struct adjacent_with_u_in_G {
   /*to define a linked list containing the identifiers
   of all vertices adjacent to vertex w in order to represent edges
   that incident with vertex w.*/
   int u;
	 /*to represent an edge (w,u) in the list of vertex w;
     note that this edge is repeated in the list of vertex u as well.
   */
   struct adjacent_with_u_in_G *more;
};
typedef struct adjacent_with_u_in_G* adjacentG;

adjacentG *LG, *LB, *LBend;
/* LG is a pointer to the pointer adjacentG. later on when we find out what is
   the vertex size of our graph, by using malloc we make LG to be
   a list of adjacentG pointers in order to construct our AJACANCY LISTs. At
   this stage LG could be considered as a list with only one element; in
   another word LG is a 1*1 array, and the only element it has is the
   adjacentG pointer. when we make LG a list of adjacentG pointers, we consider
   size of the list as the number of vertices. And while we construct
   adjacency lists we use each pointer of the list for adjacency of a vertex,
   and consider it as head of the link list(adjacency list).
*/
adjacentG edge, edge2, temp;
/*In order to handle edges of AJACANCY LISTs.
*/
int count, Pu, compNum=0, degree, bedge, tmp2;
/*bedge is used while we analysis lists of back edges to determine degree of u
*/
int *pre, *lowpt, *nd, *next_on_path, *next_sigma_element;
char *visited, *outgoing_tree_edge;

int parent, child;
//******************************************************************************

void abrt(char ch[50]) {
	printf("\n%s", ch);
//   getchar();
   exit(1);
}

void absorb_path(int x0, int xi, int end) {
   int xi_1 = x0;
   if (xi_1 != xi && xi_1 != end) {
      while (xi_1 != xi) {

   xi_1 = next_sigma_element[x0];
   next_sigma_element[x0] = next_sigma_element[xi];
   next_sigma_element[xi] = xi_1;
   /*using the variable xi_1 (temporalily) to swip next_sigma_element[x0] with next_sigma_element[xi]
   */

         /*going to append the entire LB[u](here LB[xi]) to LB[w](here LB[x0]).
				 */
         if (LB[x0]==NULL) {
            LB[x0] = LB[xi];
            LBend[x0] = LBend[xi];
         }
         else {
            LBend[x0]->more = LB[xi];
            LBend[x0] = LBend[xi];
         }
				 /*end of appending LB[u] to LB[w]
				 */
         xi_1 = xi;
         if (xi != end)
            xi = next_on_path[xi];
      }
   }
}

struct List *list;
struct List *list2;

struct Frame {
	int w;
	int v;
	int u;
	adjacentG edge;
	int start;
};

void addToStack(int w,int v, int u, adjacentG edge, int start, struct List *stack) {
	struct Frame *frame;
	frame = malloc(sizeof(struct Frame));
	frame->w = w;
	frame->v = v;
	frame->u = u;
	frame->edge = edge;
	frame->start = start;
	listAppend(stack, frame);
}

void three_edge_connectP(int w,int v, struct Frame *frame, struct List *stack);

void three_edge_connect(int w,int v) {
	struct Frame *frame;
	struct List *stack;

	stack = constructEmptyList(0, NULL);
	addToStack(w, v, 0, NULL, 0, stack);
	while(stack->length > 0) {
		frame = stack->list[--stack->length];
		three_edge_connectP(frame->w, frame->v, frame, stack);
		free(frame);
	}
}

void three_edge_connectP(int w,int v, struct Frame *frame, struct List *stack) {

   int u;
   adjacentG edge;

   if(frame->start == 1) {
	   u = frame->u;
	   edge = frame->edge;
	   goto next; //naughty hack!
   }

   nd[w] = 1;
   visited[w] = 'Y';

   next_sigma_element[w] = w;
   /*to indicate elements of the set sigma(w) and print out the component of w
   */
   next_on_path[w] = w;
   /*To represent W-path and U-path; next_on_path[w] is the next vertex that can
     be absorbed on W_path. IF next_on_path[w] == w, that means W_path is NULL.
	   In another word next_on_path[w] is a child of w on current W_path in each
		 new resulting graph.
   */

   pre[w]=count;
   lowpt[w]=pre[w];
   count = count + 1;

   edge = LG[w];
   while (edge != NULL) {
	 /*for every edge e=(w,u)=(w,edge->u) of LG[w] do the followings.
	 */
      u = edge->u;

//(Please ignore this!)      deg[w] = deg[w] + 1; // it is obvious that after facing each new edge, deg[w] must be increased


//1
      if (visited[u]=='N') {
    	  addToStack(frame->w, frame->v, u, edge, 1, stack);
    	  addToStack(u, w, 0, NULL, 0, stack);
    	  return;
    	 next:
         //three_edge_connect(u,w);
         nd[w] = nd[w] + nd[u];
         degree=0;
         bedge=0;
         if (next_on_path[u] == u) {
	 /*if U-path is null
	 */
            while (bedge<=1 && LB[u]!=NULL) {
	    /*Scan the list of back edges until more than one back edge is found or the list gets exhausted.
            */
               if (pre[u] > pre[LB[u]->u]) {
               /*The current back edge is not a self loop.
                 Note that the current edge is the head of list.
                 If this is the first one that we encounter then we memorize it
                 using pointer temp, and set the next edge as the head of list
                 in order to temporally eliminate the back edge from the list;
                 however we add it again to the list after the while loop.
                 If this is the second one that we encounter then bedge would
                 become equal to 2, and the while loop ends having the current
                 back edge (the second one) as the head of list.
               */
                  bedge = bedge + 1;
                  if (bedge==1) {
                     temp = LB[u];
                     LB[u] = LB[u]->more;
                  }
               }
               else {
               /*The current back edge is a self loop or an outdated one.
                 We free up the allocated memory, and set the next edge on the
                 list as the head of list in order to eliminate the self-loop
                 from the list.
               */
                  edge2 = LB[u];
                  LB[u] = LB[u]->more;
                  free(edge2);
               }
            }
            /*lB[U] is always the head of list, and at this time either it is
              pointing to the second back edge we encountered or it is pointing
              to NULL because in the previous while loop the list got exhausted
              and we encountered at most one back edge; note that if we have seen
              the first back edge at the end of list, still LB[u] is pointing
              to NULL because we always eliminate the first back edge from the
              list temporally. Hence we must now add the first back edge to the
              list.
              Note that if the list got exhausted and we have seen a
              self-loop at the end of the list, therefore deg[u] is at most 2 and we
              have eliminated the last edge (self-loop) from the list and
              LBend[u] is now pointing
              to that eliminated self-loop but we don't update LBend[u] (set
              it to NULL or to the possible back-edge we found) because after this
              we might use LB[u] only one more time and we append the possible back-edge
              we found to LB[w] in the following steps and since LB[u] has at most
              one element we can do the appending without using LBend[u].
            */
            if (bedge!=0) { //if we encountered one or two back-edges, the first has been
                            //eliminated. so we add it again.
               temp->more = LB[u];
               LB[u] = temp;
            }
            if (bedge<=1)
            /*since the u-path is null, if no back-edge or exactly one back-edge was found in
              the list of back edges, degree of u is 2; it means degree of u is at most 2.
            */
               degree = 2; //degree of u is at most 2.
         }
         else {
         /*if u-path is not null
         */
            while (bedge==0 && LB[u]!=NULL) {
            /*Scan the list of back edges until no back edge is found or the
              list gets exhausted.
            */
               if (pre[u] > pre[LB[u]->u])
               /*the current back edge is not a self loop.
               */
                  bedge = bedge + 1;
               else {
               /*The current back edge is a self loop or an outdated one.
               */
                  edge2 = LB[u];
                  LB[u] = LB[u]->more;
                  free(edge2);
               }
            }
            if (bedge==0)
               /*since the U-path is not null, if no back edge was found in the
                 list of back edges, degree of u is 2.
               */
               degree = 2;
         }
//1.1
         if (degree==2) { //if degree of u is at most 2.
/*//YesOrNo
   printf("It's a NO instance!");
	 exit(1);
*///YesOrNo
//            if (LB[u]==NULL) {//U_path is not NULL
            if (next_on_path[u] != u) {
               Pu = next_on_path[u];
               /*equivalent to (Pu=Pu-u) what proposed in the paper.
                 that means the next vertex can be absorbed on U-path is
                 the one after u, and u is not absorbed at all because it
                 has just been spitted out.
               */
            }
            else {//U_path is NULL
            /*Since next_on_path[u] == u, that means U-path is NULL, and there
              might be only u to be absorbed by w; but u has just been spitted
              out (deg(u)<=2) and became an isolated vertex; hence u can not be
              absorbed any more. To tell vertex w that it can not absorb u we do
              this: Pu = w. Note that although U_path is NULL, but u itself
              could be available to be absorbed on W_path by w. That is why when
              u is not also available to be absorbed, we set Pu to w(Pu = w;).
              In this way if we later set next_on_path[w] to Pu (1.4) that means
              W-path is NULL as well.
            */
               Pu = w;
            /*If at this point LB[u] is NULL that means U_path is not
              NULL and there is not any outgoing back edge of u. So not only u
	      has not absorbed any vertex, but also nothing has been appended
              so far to LB[u] because the last descendent of u(let say y) has a
              degree bigger than 2, and is waiting for an incoming back edge
              (x,y) of x(x is an ancestor of u) to be explored by DFS so that y
              is added to sigma(x). Also any other descendent of u (if exists)
              has a degree of 2, and its LB is NULL same as u.
              Note that when a vertex x absorbs a vertex y, LB[y] can not be
              NULL, and LB[y] is appended to LB[x]. Sometimes y is a child of x,
              and degree of y is 2; so x does not absorb y, however if LB[y] is
              not NULL, LB[y] is appended to LB[x] when DSF backup from y to x.
              if LB[u] is not NULL, then we are going to append the entire LB[u]
              to LB[w]. (At this time LB[u] has only one element)
            */
              if (LB[u] != NULL) { //if there is an outgoing back-edge of u;
                                   //we know there might be at most one.
                if (LB[w]==NULL)
                   LBend[w]=LB[u];
                LB[u]->more = LB[w];
                LB[w] = LB[u];
                /*Note that we can not bring the if statement after
                  "LB[w] = LB[u];" because in that case LB[w] is assigned to
                  something which in not NULL, hence the if statement never
                  get satisfied.
                */
              }
            }
            compNum = compNum + 1;
//PRINT
            list2 = constructEmptyList(0, (void (*)(void *))destructInt);
        	listAppend(list, list2);
        	listAppend(list2, constructInt(u));
            logDebug( "\nNew component found: %d",u);
//PRINT
            tmp2 = next_sigma_element[u];
            while (tmp2 != u) {
//PRINT
               logDebug( ",%d",tmp2);
               listAppend(list2, constructInt(tmp2));
//PRINT
               tmp2 = next_sigma_element[tmp2];
            }
         }//end of if (degree==2)
         else
            Pu = u;
            /*since deg(u) is not 2 then the next vertex can be absorbed on
              U-path is u.
            */
//1.2
         if (lowpt[w]<=lowpt[u])
//1.3
            absorb_path(w, Pu, 0);//(w+Pu)
         else {
          	lowpt[w]=lowpt[u];
//1.4
            absorb_path(w, next_on_path[w], 0);//(Pw)
            next_on_path[w] = Pu;
         }
      }//end of if u is NOT visited
      else {
			/*when u is visited
			*/
         if (u == v && outgoing_tree_edge[w]=='1') {
            outgoing_tree_edge[w] = '0';
						/*Once we encounter to an (w,u) edge, such that u is the
              parent of w, we consider this edge as the incoming tree edge from
							parent of w or outgoing_tree_edge of w.
              From now, no more outgoing_tree-edge (w,u) must be encountered!
							It means any other (if exists) (w,u) edge, which u=v=parent(w),
							is a parallel edge and also an outgoing back-edge!
						*/
         }
//1.5.0
       	 else if (pre[w]>pre[u]) {
				 /*if (w,u) is an outgoing back-edge of w
				 */
				    /*going to append the outgoing back_edge to LB[w]
						*/
            if (( edge2 = (adjacentG)malloc(sizeof(struct adjacent_with_u_in_G)) ) == NULL)
               abrt("Not enough memory to allocate bufferLB");
   		      edge2->u = u;
      	    edge2->more = LB[w];
            if (LB[w]==NULL)
               LBend[w]=edge2;
         	  LB[w] = edge2;
						/*set the new node as the head of the linked list of LB[w]
						*/
						/*end of appending the outgoing back_edge to LB[w].
						*/

				    if (pre[u]<lowpt[w]) {
//1.5
	             absorb_path(w, next_on_path[w], 0);//Pw
               next_on_path[w] = w;
               lowpt[w] = pre[u];
            }
         }
//1.6.0
         else if (next_on_path[w]!=w) {
            /*When pre[w]<pre[u], it means (w,u) is an incoming back-edge of w
              , however we first check to make sure Pw is not Null (by
              next_on_path[w]!=w) because if it is Null there would be nothing
              to be absorbed by w at this time.
            */
            parent = w;
            child = next_on_path[w];
            while ((parent!=child) && (pre[child]<=pre[u]) && (pre[u]<=pre[child] + nd[child] - 1) ) {
               /*while parent_path in not NULL and child is an ancestor of u
               */
               parent = child;
               child = next_on_path[child];
            }
//1.6
            absorb_path(w, next_on_path[w], parent);//Pw[w..u]
	       /*starting from child of w absorb everything on Pw[w..u]
                 until x(here variable parent) is absorbed. x lies on Pw[w..u]!
               */
            if (parent == next_on_path[parent])
               /*if X_path is NULL then Pw gets NULL as well.
               */
               next_on_path[w] = w;
            else
               /*if X_path is not NULL then Pw is set to X_path.
               */
               next_on_path[w] = next_on_path[parent];
         }
      }
      edge = edge->more;
   }
}//end of three-edge-connect procedure
//******************************************************************************

int main(int argc, char **argv) {

	list = constructEmptyList(0, (void (*)(void *))destructList);


   if (argc < 3) {
      printf("No input!\n");
      abort();
   }

   //set the logging
   if(strcmp(argv[1], "INFO") == 0) {
   		setLogLevel(LOGGING_INFO);
   	}
   	if(strcmp(argv[1], "DEBUG") == 0) {
   		setLogLevel(LOGGING_DEBUG);
   	}

   FILE *in;
   char in_filename[5000];
   int Vnum;
   int edgeNum = 0; /*initilizing the number of edges in G*/
   int r, n, ch, v, indx;
   int next_list = 1, filename_length = strlen(argv[2]);
   char ch2[max_length];
   double tsum;
   clock_t first, end;
	 logInfo("\nComputing 3edge connected components using\nDr. Tsin's algorithm(The one with reduction)...\n");
   first=clock();  //save CPU clock to variable first

//   if ((in_filename = (char*) malloc(filename_length * sizeof(char))) == NULL)
//      abrt("Not enough memory to allocate buffer1");
   for (indx=0; indx < filename_length; indx++)
      in_filename[indx] = argv[2][indx];
   in_filename[filename_length] = '\0';
   if ((in = fopen(in_filename, "rt")) == NULL)
      abrt("Cannot open input file.");
	 indx = 0;
   while ( (ch = fgetc(in)) != 10) {
      ch2[indx] = (char)ch;
      indx = indx + 1;
   }
   ch2[indx] = '\0';
   Vnum = atoi(ch2) + 1;

//*********************************Memory allocation

   if (( LG = (adjacentG*)malloc(Vnum * sizeof(struct adjacent_with_u_in_G *)) ) == NULL)
         abrt("Not enough memory to allocate buffer");

   if (( LB = (adjacentG*)malloc(Vnum * sizeof(struct adjacent_with_u_in_G *)) ) == NULL)
		abrt("Not enough memory to allocate buffer");

   if (( LBend = (adjacentG*)malloc(Vnum * sizeof(struct adjacent_with_u_in_G *)) ) == NULL)
         abrt("Not enough memory to allocate buffer");

   if (( lowpt = (int *)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( pre = (int *)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( nd = (int *)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( next_on_path = (int *)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( next_sigma_element = (int *)malloc(Vnum * sizeof(int)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( visited = (char *)malloc(Vnum * sizeof(char)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   if (( outgoing_tree_edge = (char *)malloc(Vnum * sizeof(char)) ) == NULL)
      abrt("Not enough memory to allocate buffer");

   for (indx = 0; indx < Vnum; indx++) {
      LG[indx] = NULL;
      LB[indx] = NULL;
      LBend[indx] = NULL;
      visited[indx] = 'N';
      outgoing_tree_edge[indx] = '1';
	 }
	 indx = 0;
   while ( (ch = fgetc(in)) != EOF) {
      if (indx == 10 && (ch != 62 || ch != 10))
			/*10 is a carriage return and 62 is a '>', after we reach to the maximum
			length for a number we expect new line or '>'.
			*/
   		   abrt("Your input file has an error!");
      else if (ch == 62) {
         ch2[indx] = '\0';
         indx = 0;
			   n = atoi(ch2);
         if (next_list) {
				      v = n;
							/*determine the next list
							*/
            next_list = 0;
         }
         else {
            if (( edge = (adjacentG)malloc(sizeof(struct adjacent_with_u_in_G)) ) == NULL)
               abrt("Not enough memory to allocate buffer23");
            edge->more = NULL;
   	        edge->u = n;
      	    edge->more = LG[v];
         	  LG[v] = edge;
            edgeNum = edgeNum + 1;
         }
		  }
      else if (ch == 10) {
         ch2[indx] = '\0';
         indx = 0;
				 if (!next_list) {
   			    n = atoi(ch2);
            if (( edge = (adjacentG)malloc(sizeof(struct adjacent_with_u_in_G)) ) == NULL)
               abrt("Not enough memory to allocate buffer24");
            edge->more = NULL;
            edge->u = n;
            edge->more = LG[v];
            LG[v] = edge;
            edgeNum = edgeNum + 1;
            next_list = 1;
         }
      }
      else {
         ch2[indx] = (char)ch;
         indx = indx + 1;
      }
   }
	 if (!next_list) {
      ch2[indx] = '\0';
      n = atoi(ch2);
      if (( edge = (adjacentG)malloc(sizeof(struct adjacent_with_u_in_G)) ) == NULL)
         abrt("Not enough memory to allocate buffer25");
      edge->more = NULL;
      edge->u = n;
      edge->more = LG[v];
      LG[v] = edge;
      edgeNum = edgeNum + 1;
   }
   edgeNum = edgeNum / 2;
   fclose(in);

	 logInfo( "\nComplexity of the given graph:\n|V| + |E| = %d + %d = %d\n",Vnum-1,edgeNum,Vnum+edgeNum-1);

   count = 1;
//   r = 1;
   for (r=1;r<Vnum;r++) {
      if (visited[r]=='N') {

/*//YesOrNo
         if (r>1) {
            printf("It's a NO instance!");
            exit(1);
         }
*///YesOrNo

         three_edge_connect(r,0);
         compNum++;
//PRINT
         list2 = constructEmptyList(0, (void (*)(void *))destructInt);
         listAppend(list, list2);
         listAppend(list2, constructInt(r));
         logDebug( "\nNew component found: %d",r);
//PRINT
         tmp2 = next_sigma_element[r];
         while (tmp2 != r) {
//PRINT
            logDebug( ",%d",tmp2);
            listAppend(list2, constructInt(tmp2));
//PRINT
            tmp2 = next_sigma_element[tmp2];
         }
      }
   }

   logDebug("Found components\n");

   int32_t i, j;
   FILE *fileHandle;

   fileHandle = fopen(argv[3], "w");
   fprintf(fileHandle, "" INT_STRING " ", list->length);
   for(i=0; i<list->length; i++) {
   	   list2 = list->list[i];
   	   fprintf(fileHandle, "" INT_STRING " ", list2->length);
   	   for(j=0; j<list2->length; j++) {
   		   fprintf(fileHandle, "" INT_STRING " ", ((int32_t *)list2->list[j])[0]);
   	   }
   }
   fprintf(fileHandle, "\n");
   fclose(fileHandle);
   destructList(list);


/*//YesOrNo
   printf("It's a YES instance!");
   end=clock(); //save again CPU clock to variable end
	 tsum = (end-first)/CLOCKS_PER_SEC;//CLK_TCK;  //compute total elapsed time
   printf("\nElapsed Time: %f",tsum);
	 exit(1);
*///YesOrNo

   end=clock();         //save again CPU clock to variable end
   tsum = (end-first)/CLOCKS_PER_SEC;  //compute total elapsed time
   logInfo( "\nElapsed Time: %f",tsum);
   logInfo( "\nConnected Components: %d\n",compNum);

   return 0;
}
