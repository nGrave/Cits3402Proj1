#include <stdlib.h>
#include <stdio.h>
#include <string.h>  /* for strlen() */
#include "stack.h"

 int main(int argc, char* argv[])
{

      int n = atoi(argv[1]);
      Stack s; 
      stackInit(&s , n);

      for(int i = 0; i < n; i++) {
      push(&s, i);
      printf("\rpushing %d to stack",i);
		      }


      for(int i = 0; i < n; i++) {
      printf("\r popping %d to stack",pop(&s));
		      }
     return 0;
    
} 

