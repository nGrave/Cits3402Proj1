
#ifndef _STACK_H
#define _STACK_H

//Max Size Of Stack, Need No Bigger than n*m size matrix



//Colour formatting for visualising clusters (small Matricies only)
#define EMPTY -1
#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RESET "\x1B[0m"


// Stack And Opertations
typedef struct {
    int *c;
    int top;
    int maxSize;
 } Stack;
int SizeOfStack(Stack *s);

void stackDestroy(Stack* s);

void stackInit(Stack* s, int maxSize);
	
void push(Stack* s  , int num);

int pop(Stack* s);

int isEmpty(Stack* s);

int isFull(Stack* s);

#endif  /* not defined _STACK_H */
