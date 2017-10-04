#include <stdlib.h>
#include <stdio.h>
#include <string.h>  /* for strlen() */
#include "stack.h"

 int main(int argc, char* argv[])
{



      if(argc != 2){
	      printf(RED "Usage: ./prog sizeofstack\n" RESET);
	}	      
      int n = atoi(argv[1]);
      
      
      Stack s; 
      stackInit(&s , n);

      for(int i = 0; i < n; i++) {
      push(&s, i);
      printf(BLU "pushing %d to stack\n" RESET,i);
		      }


      for(int i = 0; i < n; i++) {
      printf(RED"popping %d to stack\n"RESET,pop(&s));
		      }
     return 0;
    
} 

