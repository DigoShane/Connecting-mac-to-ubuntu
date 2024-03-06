/*If I initialize int *p;, int a=5; and then set p=&a; Then my pointer 'p' (thats how pointer are initialized) stores the adress of a. */
/*thus printf("%d",*p) and printf("%d",p).. the first one will return the value at the adress pointed by p i.e.,5.. while the second one will
return the adress of a. */
/*Similarly we cna use pointers to set the value of variables. Consider the following
int *p; (pointer)   int a;
*p=12 (read as value at adress pointed by p is 12, thus now a=12. "Dereferencing" )
Also instead of writing 
int *p;int a=20;p=&a;
\equiv int *p=&a;
*/
#include<stdio.h>
int main(){
int a=10;
int b=13;
int *p;
p=&a;
printf("The adress of the pointer p is %p*\n",p);
printf ("The adress of a is: %p\n", &a);
/*The above two won't get executed*/
printf("The above is the adress of a");
printf("a= %d\n",a);
printf("the no. that the pointer p is referencing has value %d\n",*p);
*p=12;
printf("the no. that the pointer p is referencing has value %d\n",*p);
printf("We used the pointer for a to store vlaue 12 at a.");
printf("a is now = %d",a);
*p=b;
printf("*p=b does not change the adress that p points to, the command just says that whatever value is there in b is to be now stored in the adress that p points to.");
/*printf("Adress of p is %d\n",p);*/
printf("value at p is %d\n",*p);
printf("Lets do some pointer arithmetic");
/*printf("%d\n",p);
printf("%d\n",p+1);*/
}
