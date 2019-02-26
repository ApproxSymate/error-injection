#include <stdio.h>
#include "bitflip.h"
void main(){
double a = 10.1234567891;
printf("a %.10f masked_a %.10f\n",a, bit_mask(a,16)); 
}
