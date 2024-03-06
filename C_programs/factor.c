#include <stdio.h>
#include <math.h>
#include "myfncs.h" /*STILL HAVE TO WRITE THIS FUNCTION*/

int main(){
int option,inte;
char duh;
char c='Y';
/*FILE *f=fopen("factory.txt", "r");
fclose(f);*/
double i,j,ans; /*variables*/
/*double facts[5];*/
/* Defining a 1x5 array for the factors*/
do {/* This do ends at line 60*/
printf ("1. Multiply two numbers\n");
printf ("2. Add two numbers\n");
printf ("3. Subtract two numbers\n");
printf ("4. no. of factors of number\n");
printf ("5. Totient function of a number\n");
printf ("6. Find the factors of a number\n");
printf ("please the no. corresponding to your choice");
scanf ("%d",&option);

if (option == 1){
printf ("Please enter the numbers you want to multiply\n");
scanf ("%lf%lf\n", &i, &j); /* %lf is for double*/
ans=multiply(i,j); /* include this function in myfncs.h*/
printf("%lf x %lf=%lf\n",i,j,ans);
}
else if (option ==2){
printf ("Please enter the numbers you want to add\n");
scanf ("%lf%lf\n", &i, &j);
/*& refers to the adress of the variable i, so i was stored in the memory at some adress location whihc is given by &i*/
/*The no. that the user inputs, we want stored as the value of i, so the only way of doing it is to store it at that adress*/
/*Thus any operation we want done on i, the program will do the operation on the value stored at that adress location.*/
ans=add(i,j); /* inclucde this function in myfncs.h*/
printf("%lf + %lf=%lf\n",i,j,ans);
}
else if (option ==3){
printf ("Please enter the numbers you want to subtract\n");
scanf ("%lf%lf\n", &i, &j);
ans=substract(i,j); /* include this fn in "myfncs.h"*/
printf("%lf - %lf=%lf\n",i,j,ans);
}

else if (option ==5){
printf ("Please enter the number whose euler's totient you want evaluated\n");
scanf ("%d\n", &inte);
ans=Toit(inte);
printf("Toitient of %d is %lf\n",inte,ans);
}

else if (option ==4){
printf ("Please enter the number whose sum of factors you want evaluated\n");
scanf ("%d\n", &inte);
ans=factorsmof(inte); /*input this function in "myfncs.h"*/
printf("Sum of factors of %d is %lf\n",inte,ans);
}

else if (option ==6){
printf ("Please enter the integer whose factors you want stored in a text file factory \n");
scanf ("%d\n", &inte);
factors(inte);
FILE *f=fopen ("factory.txt", "r");
char input[512];
int line =0;
while( fgets( input, 512, f) ){ 
line++;
printf("Line:%d -> %s", line, input);
};
printf("\n\ndidnt think i would go thought the effort of priting it on screen did ya\n");

}

else {
printf ("Please input a proper option\n");
}

printf ("Do you want to continue, ['y','Y']\n");
scanf(" %c", &c);
}while(c=='y' || c=='Y');/*The corresponding starting point is at line 12*/

return 0;
}
