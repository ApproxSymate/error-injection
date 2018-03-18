#include "bitflip.h"
//#define TEST
#define LIMIT 16

int generate_bit_position()
{
        int num;
        srand(time(NULL));
        num = rand();
        return (num % LIMIT);
}
double bitflip_double (double num)
{
        double temp_num = num;
        long a = *((long*)&temp_num);

        long temp_result = bitflip_long (a);
  //      printf( "\n %ld %ld \n", a , temp_result);
        return *((double*)&temp_result);
}

float bitflip_float (float num)
{
        float temp_num = num;
        long a = *((long*)&temp_num);

        long temp_result = bitflip_int (a);
        //sprintf( "\n %ld %ld \n", a , temp_result);
        return *((float*)&temp_result);
}

int bitflip_int (int num)
{
        int pos = generate_bit_position();
  //      printf(" %d ", pos);
        return num ^ (1<<pos);
}

long bitflip_long (long num)
{
        int pos = generate_bit_position();
        //printf(" %d ", pos);
        return num ^ (1<<pos);
}


int bitflip_int_at_pos (int num, int pos)
{
        //int pos = generate_bit_position();
  //      printf(" %d ", pos);
        return num ^ (1<<pos);
}

#ifdef TEST
int main()
{
        double num = 0.8;
        printf("%lf\n",bitflip_float(num));
}
#endif