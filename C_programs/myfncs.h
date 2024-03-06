


int add (double x, double y){
return(x+y);/* Simple enough */
}

int multiply (double x, double y){
return(x*y);
}

int substract (double x, double y){
return(x-y);
}

int Toit (int inte){
/* We want to find the Euler Toitinet*/
/*so its the nos. less than inte which are co prime to it*/
/*since #coprime+#prime=inte*/
/* we calculate #prime and use the above formula */
int i;
int count=0;
for (i=1; i<inte; i=i+1 ){
if ( inte % i ==0){ /*To check if 'i' divides 'inte'*/
count=count+1;
}}
return(inte-count);
}

int factorsmof (int inte){
/*We just calculate the sum of the factors*/
/*here 'count' is being used as the prototypical sum or S in S=S+#*/
int i;
int count=0;
for (i=1; i<inte; i=i+1 ){
if ( inte % i ==0){
count=count+i;
}}
return(count);
}

/*This function has 2 inputs, integer and filename(array of 100 character)*/

int factors( int x ){
FILE *pToFile2=fopen("factory.txt", "a");
int i;
for (i=1; i<x; i++){
if ( x % i ==0){
if (pToFile2 != NULL){
fprintf(pToFile2, "%d\n", i);
} else {
printf ("File Not Found!\n");
}
}
}
}
