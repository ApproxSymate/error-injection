#include "bitflip.h"
//#define TEST


int ipow(int base, int exp)
{
    int result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }

    return result;
}
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

int limitprecision_int (int num, int precision) // saturated max value for signed int
{
		int limit = ipow(2,precision-1);
		if(num > limit-1)
			return limit -1;
		if (num < (0-limit))
			return 0 - limit;
		return num;
}

#ifdef TEST
int main()
{
        double num = 0.8;
        printf("%lf\n",bitflip_float(num));
        int num2 = 15;
        printf ("%d %d \n", limitprecision_int(num2,4), limitprecision_int(num2,5));
}
#endif
