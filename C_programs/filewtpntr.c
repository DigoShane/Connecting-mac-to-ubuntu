#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include "myfncs.h"

int main(){
FILE * fPointer;
FILE *pToFile=fopen("doc_1.txt", "r");
FILE *pToFile1=fopen("doc_1.txt", "a");
/*We need a pointer to read/write to the file, where do we read/write to?
Thus for that we ned the file pointer.
r is to read from a file 
w is to write to a file
a is add more crap
*/

fPointer = fopen("bacon.txt", "w"); /*This helps us open a file and add crap to it, the w is as defined above*/
fprintf(fPointer, "Bacon is good\n"); /*prints to a file. Note that now the pointer points to the beginning of the file, 1st line*/
fclose(fPointer);/*Good practise to close the file*/

/*If we were to change "I love cheese" to "bacon is good" then the next time we ran this, the progrma would essentially
delete the previous file bacon.txt, create a new file bacon.txt and write to it "bacon is good". This is becuase we are
using the "w" command, if we wanted to add anything, we would have had to use the "a" command whcich stands for append.
*/

/*Now we want to read from the file*/
/*Due to the structure of the C program, we can't just use fopen with "r", we must read text form the file character by character*/
int line =0;/*to count the no. of lines*/
char input[512];/*we want to store the contents of the file in this character array with 512 slots*/
/*syntax: fgets(char *str, int n, FILE *stream)
str -> This is the pointer to an array of chars where the string read is stored.n -> This is the macimum no. of characters to be read (including the final null character). Usually, the length of the array passed as str is used.
stream -> This is the pointer to a FILE object that identifies the stream where the characters are read from.

If we consider scanf vs fgets, then scanf does have a limitation on the no. of characters, it takes in form the user. Thus if the user enters more than the allocated memory, we will end up writing past the allocated memeory to where the meory of other program gets 'modfied' which might cause the systme to crash or even worse BSOD in windows.

fgets is also known as the blocking fn. When the compiler reaches the fgets line
fgets waits for the enter key to be pressed and all the characters upto and
including the enter key is stored in the specified file. The function fgets
returns a non-zero value. 

The while loop syntax is such that it will keep running till the argument of the
while function is a zero value (which stands for False in Boolean). Thus as fgets gives us a nonzero value, we can use it as argument for the while function
*/

while( fgets( input, 512, pToFile) ){ /*ends at 48-ish*/
line++;
printf("Line:%d -> %s", line, input);
}/*starts at 45-ish*/
/*fPointer = fopen("bacon.txt", "r");fopen opens the file and "r" is specifiying that we want to read it*/
fclose(pToFile);
printf("\n\nend of Program\n");

/*Now lets append to a file in some more detail*/
int num=40;
char single = 'C';
char myWord[4] = "Cat";
 
if (pToFile1 != NULL){
fprintf(pToFile1, "%d %c %s %s\n", num, single, myWord, myWord);
} else {
printf ("File Not Found!\n");
}


return 0;

}
