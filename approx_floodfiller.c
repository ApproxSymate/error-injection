/*
 * FloodFiller.c
 *
 *  Created on: Oct 26, 2017
 *      Author: himeshi
 *
 * Portions Copyright 2017 National University of Singapore
 *
 * The C version of the imagefill Java code of the EnerJ benchmark.
 * https://sampa.cs.washington.edu/new/research/approximation/enerj.html
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bitflip.h"
#include <math.h>
//#include "klee/klee.h"

#define WIDTH 32
#define HEIGHT 32
#define MAX_STACK_SIZE 500
int *xstack;
int *ystack;
int stackSize;
int image[WIDTH][HEIGHT];

int image_bitflip[WIDTH][HEIGHT];


int targetColor;
int maxStackSize;

int fill(int x, int y);

int fill_bitflip(int x, int y);

int getPix(int x, int y);
void setPix(int x, int y, int c);
int popx();
int popy();
void fillLine(int x1, int x2, int y);
void push(int x, int y);



int main() {
  int i, j;
  char c;

  FILE *fp = fopen("floodfiller_input.txt", "r");

  for (i = 0; i < HEIGHT; ++i) {
    for (j = 0; j < WIDTH; ++j) {
      c = getc(fp);
      if (c != '\n' && c != '\r'){
			image[i][j] = c - '0';
			image_bitflip[i][j] = c - '0';
		}
    }
  }

  //initialize
  targetColor = 2;
  stackSize = 0;
  maxStackSize = MAX_STACK_SIZE;
  xstack = (int *) malloc(sizeof(int) * maxStackSize);
  ystack = (int *) malloc(sizeof(int) * maxStackSize);

  fill(0, 0);



  targetColor = 2;
  stackSize = 0;
  maxStackSize = MAX_STACK_SIZE;
  xstack = (int *) malloc(sizeof(int) * maxStackSize);
  ystack = (int *) malloc(sizeof(int) * maxStackSize);

  fill_bitflip(0, 0);

	double diff_bitflip = 0.0;

  for (i = 0; i < HEIGHT; ++i) {
    for (j = 0; j < WIDTH; ++j) {
			//printf("%d, ",image[i][j] );
			//printf("%d, ",image_bitflip[i][j] );
			//printf("%d \n",image[i][j]-image_bitflip[i][j]);
			diff_bitflip += abs(image[i][j]-image_bitflip[i][j]);
    }
  }
  printf("%lf, ", diff_bitflip/(WIDTH*HEIGHT));

  //printf("done\n");
}

int getPix(int x, int y) {
  if (x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT)
    return -1;
  else
    return image[x][y];
}

void setPix(int x, int y, int c) {
  if (x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT)
    return;
  else
    image[x][y] = c;
}


int getPix_bitflip(int x, int y) {
  if (x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT)
    return -1;
  else
    return bitflip_int(image_bitflip[x][y]);
}


void setPix_bitflip(int x, int y, int c) {
  if (x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT)
    return;
  else
    image_bitflip[x][y] = bitflip_int(c);
}


int popx() {
  if (stackSize == 0)
    return -1;
  else
    return xstack[stackSize - 1];
}

int popy() {
  int value = ystack[stackSize - 1];
  stackSize--;
  return value;
}

void fillLine(int x1, int x2, int y) {
  int x;
  if (x1 > x2) {
    int t = x1;
    x1 = x2;
    x2 = t;
  }
  for (x = x1; x <= x2; x++)
    setPix(x, y, targetColor);
}

void fillLine_bitflip(int x1, int x2, int y) {
  int x;
  if (x1 > x2) {
    int t = x1;
    x1 = x2;
    x2 = t;
  }
  for (x = x1; x <= x2; x++)
#ifdef AGGRESSIVE
    setPix_bitflip(x, bitflip_int(y), targetColor);
#else
    setPix_bitflip(x, (y), targetColor);
#endif
}


void push(int x, int y) {
  stackSize++;
  if (stackSize == MAX_STACK_SIZE) {
    int *newXStack = (int *) malloc(sizeof(int) * maxStackSize * 2);
    int *newYStack = (int *) malloc(sizeof(int) * maxStackSize * 2);
    memcpy(newXStack, xstack, sizeof(int) * maxStackSize);
    memcpy(newYStack, ystack, sizeof(int) * maxStackSize);
    xstack = newXStack;
    ystack = newYStack;
    maxStackSize *= 2;
  }
  xstack[stackSize - 1] = x;
  ystack[stackSize - 1] = y;
}

int fill(int x, int y) {
  int width = WIDTH;
  int height = HEIGHT;
  int color = getPix(x, y);
  fillLine(x, x, y);
  int newColor = getPix(x, y);
  setPix(x, y, color);

  if (color == newColor)
    return -1;

  stackSize = 0;
  push(x, y);

  while (1) {
    x = popx();
    if (x == -1)
      return 1;
    y = popy();

    if (getPix(x, y)!=color) continue;
    int x1 = x;
    int x2 = x;

    while (getPix(x1, y)==color && x1>=0) x1--; // find start of scan-line
    x1++;

    while (getPix(x2, y)==color && x2<width) x2++;  // find end of scan-line
    x2--;

    fillLine(x1,x2,y); // fill scan-line

    int inScanLine = 0;
    int i;
    for (i = x1; i <= x2; i++) { // find scan-lines above this one
      if (!inScanLine && y > 0 && getPix(i, y - 1) == color) {
	push(i, y - 1);
	inScanLine = 1;
      } else if (inScanLine && y > 0 && getPix(i, y - 1) != color)
	inScanLine = 0;
    }

    inScanLine = 0;
    for (i = x1; i <= x2; i++) { // find scan-lines below this one
      if (!inScanLine && y < height - 1 && getPix(i, y + 1) == color) {
	push(i, y + 1);
	inScanLine = 1;
      } else if (inScanLine && y < height - 1
		 && getPix(i, y + 1) != color)
	inScanLine = 0;
    }
  }
}


int fill_bitflip(int x, int y) {
  int width = WIDTH;
  int height = HEIGHT;
#ifdef AGGRESSIV
  int color = getPix_bitflip(x, y);
#else
  int color = bitflip_int(getPix_bitflip(x, y));
#endif
  fillLine_bitflip(x, x, y);
#ifdef AGGRESSIVE
  int newColor = (getPix_bitflip(x, y));
#else
  int newColor = bitflip_int(getPix_bitflip(x, y)); //negate bitflip effect
#endif

  setPix_bitflip(x, y, color);

  if (color == newColor)
    return -1;

  stackSize = 0;
  push(bitflip_int(x), bitflip_int(y));

  while (1) {
    x = popx();
    if (x == -1)
      return 1;
    y = popy();

    if (getPix_bitflip(x, y)!=bitflip_int(color)) continue;
    int x1 = x;
    int x2 = x;

    while (getPix_bitflip(x1, y)==bitflip_int(color) && x1>=0) x1--; // find start of scan-line
    x1++;

    while (getPix_bitflip(x2, y)==bitflip_int(color) && x2<width) x2++;  // find end of scan-line
    x2--;

    fillLine_bitflip(x1,x2,y); // fill scan-line

    int inScanLine = 0;
    int i;
    for (i = x1; i <= x2; i++) { // find scan-lines above this one
      if (!inScanLine && y > 0 && getPix_bitflip(i, y - 1) == bitflip_int(color)) {
	push(bitflip_int(i), bitflip_int(y - 1));
	inScanLine = 1;
} else if (inScanLine && y > 0 && getPix_bitflip(i, y - 1) != bitflip_int(color))
	inScanLine = 0;
    }

    inScanLine = 0;
    for (i = x1; i <= x2; i++) { // find scan-lines below this one
      if (!inScanLine && y < height - 1 && getPix_bitflip(i, y + 1) == bitflip_int(color)) {
	push(bitflip_int(i), bitflip_int(y + 1));
	inScanLine = 1;
      } else if (inScanLine && y < height - 1
		 && getPix_bitflip(i, y + 1) != bitflip_int(color))
	inScanLine = 0;
    }
  }
}
