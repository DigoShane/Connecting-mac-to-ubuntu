#include <stdio.h>
#include <math.h>
#include "try.h"

int main()
{
	int X = 20;
	int y = 2;
	int Z,A;
	Z = add(X,y);
	A = multiply(X,y);
	printf("%d + %d=%d\n",X,y,Z);
	printf("%d X %d=%d\n",X,y,A); 
}
