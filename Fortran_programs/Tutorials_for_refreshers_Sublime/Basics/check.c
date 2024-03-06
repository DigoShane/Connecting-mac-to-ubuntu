#include <stdio.h>

int main() 
{

   char str[100];
   double i;

   printf( "Enter a value :");
   scanf("%s %lf", str, &i);

   i =55;

   printf( "\nYou entered: %s %lf", str, i);

   return 0;
}