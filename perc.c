// CITS 3402 Project
// AUTH : Nathan Graves 21256779 & Matthew Hall TODO
// Compiled with gcc-4.9 -fopen -o {filename} perc.c stack.c 
//
//
// USAGE : ./perc n p {OPTIONS} :
// 		p is probabilty double from 0-1
// 		OPTIONS:
// 		-o t  : USE OpenMp  where t is the number of threads required
// 		-c t  : Compare OpenMP Para with regular t is Number of threads
// 		-p    : print Seeded Matrix and largest cluster Representation ( n<200 ) (Not defined output with openmp)
// 		-v    : verticle or horizontal percolation allowed ( without its vertical AND Horizontal ) 
// 		TODO
// 		-b    :Bond Site perculation IDEA
//				X S X S X S
// 				S N S N S N     N = Node 
// 				X S X S X S	S = Site Percolation Site
// 				S N S N S N     X = VOID
// 				X S X S X S
// 				S N S N S N
//      TODO: KNOWN BUGS:
//      	./perc n 1 drops some elements with large value of n  10000+ for eg.  -Uncomment lines 205 - 210 to Help Diagnose 
//		./perc n 1 -o (or -c) malloc errors (accessing already freed memory) 
//
//	ISSUES   : parallel code seems to always be slower regardless of number of threads cache hit/miss problems in way array is currentlty split
//		   ie: Split in rows so a 10*10 matrix with 4 threads thread 1 will do up to  mat[2][10] thread 2 [4][10] thread3 [6][10] and thread 4 gets the biggest chunk till mat[10][10] 
//		   becasue of the leftover after division. 
//			
//		  This is not including that the peices havent been put back together yet still thinking about implementation of dfs across threads becasue of the wrap around nature
//		  there is some stuff on SO for dfs with parralell regions but use of different data structures needed IE graph or tree.
//
//
//
//	Effiency : Recursive Vs Iterative DFS ?
//		   Minimizing cache miss somehow? 
//		   derefrenicng vs Lookups  -- ALL Relativaly Minimal in comparaison
//		
//	Paperwork: 
//		  Graphs on speed up etc. 
//		  see project Description.
//
//	More: ? 
//
//	DUE DATE 8th of October .
//
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "stack.h"
#include <string.h>
#include <omp.h>


float percProb(){
    return (double)rand() / (double)((unsigned)RAND_MAX + 1);
}

void  printLargestCluster(int **mat,int n ,int m, int ldx){ 
    for(int i =0; i<n; i++){
        for(int j =0 ; j < m; j++){
  
            if(mat[i][j]==ldx ) printf(GRN  "%c" RESET, 'x');
            else  printf(RED  "%c" RESET  , 'o');
           
        }
        printf("\n");
    }
}
void printMatrix(int **mat, int n){

    for(int i =0; i<n; i++){
        for(int j =0 ; j < n; j++){
          
            if(mat[i][j]==1) printf(GRN  "%d" RESET, mat[i][j] );
            else  printf(RED  "%d" RESET  ,**mat );
            
        }
        printf("\n");
    }
}

//Site Seeding 
void SeedMatrix(int **mat, int n, float p){
    
	//Seed RNG
      srand(time(NULL));
 
      //Seeding Matrix
    for(int i =0; i<n;  i++){
	    for(int j =0;j < n ; j++){
    		    mat[i][j] =  (int)( percProb()< p);
	    }
    }
}

int findCluster( int size, int l,  int **mat, int print, int percCond, int tid){   

   	//Set all sites to unvisited

     int **visited = malloc(l * sizeof(int *));
     int *t = malloc(l * size  * sizeof(int));
        for(int i = 0; i < size; i++){
	 visited[i] = t +( i * size);
	}

    for(int i= 0; i < l; i++){
	    for(int j =0 ; j < size; j++){
		    visited[i][j] = 0; 
	     }
    }

    //Variables for tracking cluster Stats	
    int clusIdx = 1;
    int clusHeight =0;
    int clusWidth = 0;
    int elementsInCluster = 0; 
    int largestCluster = 0; 
    int lcidx =0;	
    int percolates = 0;

    //Loop through every element
    for(int i =0; i< l; i++){
    	for(int j=0 ;j < size ;j++){    
	//Move along becasue this site is Part of another cluster        
	if(visited[i][j] == 1) continue; 
	//Mark as Visited       
	visited[i][j] = 1;;    
	// If This site has a bond	
	if(mat[i][j]  == 1){
		//Time This Cluster
		//clock_t t = clock(); 
		//Initiate a new stack -- Space optimiation by subtracing the largest cluster As garunteed not to be larger than that;
		 Stack *s;
		 s= malloc(sizeof(*s)*l);
		 stackInit(s, l*size - largestCluster);

		//printf("-------------------Start Of New Cluster ID= %d-----------------\n", clusIdx);
		push(s, i * l  + j );		
	        int rowsOccupied[l];
	        int colsOccupied[size];

		for(int m=0; m<l; m++){
		  rowsOccupied[m] = 0;
	        }

		for(int m=0; m < size; m++){
		  colsOccupied[m] = 0;
		}

		rowsOccupied[i] = 1;
		colsOccupied[j] = 1;
       			         
         
  // if(tid!=0)   printf("%d bp5\n",tid);    
		 // Search for connected neighbours	
		 while(!isEmpty(s)) {
			int cur = pop(s);
			elementsInCluster++; 

		        //printf("Thread %d - Popping Top Element off the stack CUR=%d : mat[%d][%d] 1d equiv is %d\n",tid, cur, cur/l, cur%size , cur);
			//Convert to row and col	
			int row = cur/size;
			int col = cur%size;
			
			int above =  (row-1 + l) % l ;
			int below =  (row+1 + l) % l ;
			int left   =  (col-1 + size) % size ;
			int right =  (col+1 + size) % size ;
			//To Help Visualize individual clusters ( Comment out for large clusters)
			if(print == 1)
		//	mat[row][col] = clusIdx; 
			//Search Above,Below,Left and Right For Bond Sites

			//printf("Thread %d, row = %d col = %d above = %d below =%d\n",tid ,row ,col, above ,below);
			
			if( mat[above][col] == 1 && visited[above][col] ==0) {
			//	printf("Thread %d Pushing Above Element  mat[%d][%d] (1d = ele %d) to the stack\n",tid ,above,col, above*size+col );
				push(s, above * size  + col );
			      	visited[above][col] = 1;
				rowsOccupied[above] = 1;       		
			}
			if( mat[below][col] == 1 && visited[below][col] ==0) {
			//	printf("Thread %d Pushing Below Element  mat[%d][%d] (1d = ele %d) to the stack \n",tid,below,col, below*size+col );
				push(s, below * size  + col );
			       	visited[below][col] = 1;
				rowsOccupied[below]= 1;	
			}
			if( mat[row][left] == 1 && visited[row][left] ==0) {
			//	printf("Thread %d Pushing Left Element mat[%d][%d] (1d = ele %d) to the stack \n",tid ,row,left, row*size+left );
				push(s,  row * size  + left );
			      	visited[row][left] = 1;
				colsOccupied[left] = 1;
			}
			if( mat[row][right] == 1 && visited[row][right] ==0) {
		  	//	printf("Thread %d Pushing Right Element  mat[%d][%d] (1d = ele %d) to the stack \n",tid, row,right, row*size+right ); 
		       		push(s, row * size  + right );
			       	visited[row][right] = 1;
				colsOccupied[right] = 1;
			}
			//printf("Cur is[%d][%d] Above is [%d][%d] -Below is [%d][%d]- Left is [%d][%d] - Right is [%d][%d]\n", row, col, above, col, below, col, row,left,row,right  );	
		 }	

	 //freeStack
	 stackDestroy(s);	 
	 free(s);
         
	//MEMORY LEAK HELP - UNCOMMENT FOLLOWING LINES AND RUN ON LARGE ARRAY 10,000+ seems to miss a few elements for some reason
	// 	 for(int i = 0 ; i < size; i++){
	//	 for(int j =0; j< size; j++){
	//	if(visited[i][j] == 0) printf("Not Visited[%d][%d] = %d\n", i,j,visited[i][j] );
        //	}}	
	
       	//End Of cnnected neigbur search
	 for(int count=0; count<size; count++){
		 clusWidth += colsOccupied[count];
	  }
	 for(int count=0; count<l; count++){
		 clusHeight += rowsOccupied[count];
	 }

	 //DOES IT PERCOLATE
         if(percCond ==1 ){ 
	  if(clusWidth==size || clusHeight ==size)
		percolates = 1; //Use this line for horizontal or vertical percolation	
	 }
	 else if(clusWidth==size && clusHeight==size) percolates = 1; //Use this line for two-way percolation

	 //printf("----Custer %d Has %d Elements, It is %d elements wide and %d elements high\n", clusIdx, elementsInCluster , clusWidth, clusHeight);	
         
	 //largest Cluster Check
	 if(elementsInCluster > largestCluster){
		 lcidx = clusIdx;
		 largestCluster = elementsInCluster;              
	 }
	 //Reset Counters
         elementsInCluster=0;
	 clusIdx++;
	 clusWidth =0;
	 clusHeight =0;
	}
    //End Of Matrix    
	}
     }
    free(t);    
    free(visited);
    printf("Largest Cluster = %d\n", largestCluster);
 
    //Visualize Largest Cluster
    if(print == 1){
    printLargestCluster(mat, l, size , lcidx);
    }
  
    return percolates;
}

int runParalell(int n, int  **mat,  int percCond, int printMat,  int numThreads ){

     //OPEN MP START
     int nthreads, tid;   
     clock_t time  =clock();
     omp_set_num_threads(numThreads);  
     printf("numThreads = %d\n", numThreads);
     #pragma omp parallel private(nthreads, tid)
     {
     int arrPartSize = n/numThreads;	
     int leftOvers= n - (numThreads* arrPartSize);     
     tid = omp_get_thread_num();
    
     int end;
     int start; 
     int peiceSize;

     for(int i = 0; i < numThreads ; i++){
	  start = (arrPartSize *tid);
	  end   = start + arrPartSize;
        if(tid == numThreads-1) end += leftOvers;
     }
     printf("Thread %d Ready for work on my part of matrix which is from element %d to element %d\n", tid, start, end);
     peiceSize = end -start;
     clock_t t =clock();
     //Find Clusters 
     int perc  = findCluster(n ,peiceSize,   mat+start , printMat, percCond, tid); 
     t = clock() - t;
     double time_taken = ((double)t/CLOCKS_PER_SEC); // in seconds
     printf("Time taken in finding Clusters is %f\n", time_taken);
     if(tid == 0) {
	//printf("Master Thread WAITS FOR CHILDREN - Check Percolation Here  %d\n",omp_get_num_threads());
     }
     if( perc == 1) 
	printf("My part Percolates Part %d\n", tid) ;    
     else printf("My part Does Not Percolate %d\n",tid);	

     }
     //All Threads Rejoined Master Here
     time = clock() -time;
     double totTime = ((double)time/CLOCKS_PER_SEC); //Seconds
     printf(RED "Time Taken For  all peices is %f\n" RESET , totTime); 	
     return 0;
}

int runNormal(int n, int **mat, int printMat, int percCond){
    
     clock_t t =clock();
      int perc = findCluster(n , n,  mat , printMat, percCond, 0); 
     t = clock() - t;
     double  time_taken = ((double)t/CLOCKS_PER_SEC); // in seconds
     printf(RED "Time taken in series is %f\n"RESET , time_taken);
    
     if( perc == 1) 
	printf("Matrix Percolates\n") ;    
     else printf("Matrix Does Not Percolate\n");	
 
    return 0;     
}


int main(int argc , char* argv[]){

    if(argc < 3){
    printf(BLU "Usage: ./perc {MatrixSize} {SeedingProbability} {Optional Flags}\n-p: Prints to console a visual representation of matrix- ONLY USE FOR SMALL MATRIX\n-v: Matrix only has to percolate vertically or horizontially -DEFAULT IS BOTH\n-b: Bond Percolation - DEFAULT IS SITE PERCOLATION\n-o Run with openMp\n-c Compare with regular\n"RESET );
    exit(EXIT_FAILURE);
    }

    //Command Line Options
    int n = atoi(argv[1]);
    float prob =atof(argv[2]);
    int printMat=0;
    int bPerc =0;
    int percCond =0;
    int ompflag =0; 
    int compareFlag =0;
    int idx =3;
    int numThreads = 4; //TODO Test with more 
    while( idx < argc){
	if(strncmp(argv[idx] ,"-p",2)==0 && n <=200) printMat =1;
	if(strncmp(argv[idx], "-v",2)==0) percCond =1;
	if(strncmp(argv[idx], "-b",2)==0) bPerc =1; 
	if(strncmp(argv[idx], "-o",2)==0) {
		ompflag =1;
		numThreads = atoi(argv[idx+1]) ;  
	}
	if(strncmp(argv[idx], "-c",2)==0){
	       	compareFlag = 1;
		numThreads = atoi(argv[idx+1]) ;
	}	
    	idx++;
    }

     printf("nthreads= %d\n", numThreads);  
    //Matrix Memory Allocation; 
     int **mat = malloc(n* sizeof(int *));
     int *temp = malloc(n * n  * sizeof(int));
     for(int i = 0; i < n; i++){
	 mat[i] = temp + ( i * n);
     }

     //Seed The Matrix	
     clock_t t;
     t = clock();
     SeedMatrix(mat,n, prob); 
     t = clock() - t;
     double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
     printf("Time taken in Seeding Matrix is %f\n", time_taken);
      
     //Copy of matrix for comparison  
     int **mat2 = malloc(n* sizeof(int *));
     int *temp2 = malloc(n * n  * sizeof(int));
     for(int i = 0; i < n; i++){
	 mat2[i] = temp2 + ( i * n);
     }

     for(int i = 0 ; i < n ; i++){
	for(int j =0 ; j < n ; j++){
	mat2[i][j] = mat[i][j];
	}
     }
     //Print Matrix if requested
     if(printMat == 1){
     printMatrix(mat,n); 
     }  
   
     //FIND CLUSTERS
     if(ompflag ==0){
     runNormal(n,mat,printMat,percCond); 
     }
     if(ompflag ==1 || compareFlag== 1 ){ 
     runParalell(n,mat2,printMat,percCond,numThreads);
     }
     //Free Memory 
     free(temp);
     free(mat);
     return 0; 
}
