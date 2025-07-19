#include<float.h>
#include<math.h>
#include<omp.h>
#include<stdbool.h>
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>

double linearize(int c)
{
	return (c < 11) ? c / 3294.6 : pow((c + 14.025) / 269.025, 2.4);
}

double relativeLuminance[256][256][256];

void memoize()
{
	// Coefficients for relative luminance formula
	// G - 0	R - 1	B - 2
	double coefficients[3] = { 0.7152, 0.2126, 0.0722 };

	// Results of linearize(x) in the interval [0, 255]
	// multiplied by the coefficients above
	double values[3][256];
	//memoize R, G, and B
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 11; j++)
			values[i][j] = coefficients[i] * j / 3294.6;
		for (int j = 11; j < 256; j++)
			values[i][j] = coefficients[i] * pow((j + 14.025) / 269.025, 2.4);
	}
	
	// Memoize relative luminance
	for (int i = 0; i < 256; i++)
	for (int j = 0; j < 256; j++)
	for (int k = 0; k < 256; k++)
		relativeLuminance[i][j][k] = values[0][i] + values[1][j] + values[2][k] + 0.05;
}

int getChannelValue(char channel[])
{
	int result;
	
	do {
		printf("%s: ", channel);
		scanf("%x", &result);
	} while(0 > result || result > 255);
	
	return result;
}

double getRatio()
{
	double result;
	do {
		printf("Ratio (1 to 21 inclusive): ");
		scanf("%lf", &result);
	} while(1 > result || result > 21);
	
	return result;
}

void swap(char *x, char *y)
{
    int temp = *x;
    *x = *y;
    *y = temp;
}

void update(double* a, double* b, double c)
{
	if(*a > c) *a = c; // lower bound
	if(*b < c) *b = c; // upper bound
}

int min(int a, int b)
{
	return a > b ? b : a; 
}

void findPair()
{
	double 	ratio = getRatio(), 
			minimum = DBL_MAX,
			current, rl1, rl2;
	int 	G1, R1, B1, G2, R2, B2;

	// rl1 / rl2 ≥ 1
	// rl1 ≥ rl2
	// relativeLuminance[x0][x1][x2] ≥ relativeLuminance[x3][x4][x5]
	// relativeLuminance[x0][x1][x2] - relativeLuminance[x3][x4][x5] ≥ 0

	// #pragma omp parallel for
	// for(int x0 = 0; x0 < 256; x0++)
	// for(int x1 = 0; x1 < 256; x1++)
	// for(int x2 = 1; x2 < 256; x2++)
	// {
	// 	rl1 = relativeLuminance[x0][x1][x2];
	// 	for(int x3 = x0  ; x3 < 256; x3++)
	// 	for(int x4 = x1  ; x4 < 256; x4++)
	// 	for(int x5 = x2-1; x5 < 256; x5++)
	// 	{
	// 		current = fabs(rl1 - ratio * relativeLuminance[x3][x4][x5]);

	// 		if (minimum > current)
	// 		{
	// 			minimum = current;
	// 			G1 = x0;
	// 			R1 = x1;
	// 			B1 = x2;
	// 			G2 = x3;
	// 			R2 = x4;
	// 			B2 = x5;
	// 		}
	// 	}
	// }

	double	p_ratio = pow(ratio, 1/2.4), r_rl2;
	int 	x0_min, x0_max,
			x1_min, x1_max,
			x2_min, x2_max,
			rl1_compo_max = (int)ceil(255 / p_ratio)+1;
		
	#pragma omp parallel for
	for(int x3 = 0; x3 < rl1_compo_max; x3++)
	{
		x0_min = min((int)floor(p_ratio * x3), 255);
		x0_max = min((int)ceil (ratio   * x3), 255);
		for(int x4 = 0; x4 < rl1_compo_max; x4++)
		{
			x1_min = min((int)floor(p_ratio * x4), 255);
			x1_max = min((int)ceil (ratio   * x4), 255);
			for(int x5 = 0; x5 < rl1_compo_max; x5++)
			{
				x2_min = min((int)floor(p_ratio * x5), 255);
				x2_max = min((int)ceil (ratio   * x5), 255);
				r_rl2 = ratio * relativeLuminance[x3][x4][x5];
				for(int x0 = x0_min; x0 < x0_max; x0++)
				for(int x1 = x1_min; x1 < x1_max; x1++)
				for(int x2 = x2_min; x2 < x2_max; x2++)
				{
					current = fabs(relativeLuminance[x0][x1][x2] - r_rl2);

					if (minimum > current)
					{
						minimum = current;
						G1 = x0;
						R1 = x1;
						B1 = x2;
						G2 = x3;
						R2 = x4;
						B2 = x5;
					}
				}
			}
		}
	}
	
	printf("#%02X %02X %02X - #%02X %02X %02X\n", R1, G1, B1, R2, G2, B2);
}

void findPartner()
{
	double 	rl1, ratio = getRatio(), 
			minimum = DBL_MAX, current;
	int 	R1 = getChannelValue("R1"),
			G1 = getChannelValue("G1"),
			B1 = getChannelValue("B1"),
			/*\                                 /*\
			|*|---------------------------------|*|
			|*| x[0][0] - index 0  x[1][0] - G2 |*|
			|*| x[0][1] - index 1  x[1][1] - R2 |*|
			|*| x[0][2] - index 2  x[1][2] - B2 |*|
			|*|---------------------------------|*|
			\*/                               /*\*/
			x[2][3];
	
	rl1 = relativeLuminance[R1][G1][B1];
	for(x[0][0] = 0; x[0][0] < 256; x[0][0]++) 
	for(x[0][1] = 0; x[0][1] < 256; x[0][1]++) 
	for(x[0][2] = 0; x[0][2] < 256; x[0][2]++)
	{	
		current = fabs(rl1 - ratio * relativeLuminance[x[0][0]][x[0][1]][x[0][2]]);
		
		if (minimum > current)
		{
			minimum = current;
			for(int i = 0; i < 3; i++)
				x[1][i] = x[0][i];
			//printf("\r#%02X%02X%02X - #%02X%02X%02X", R1, G1, B1, R2, G2, B2);
			//fflush(stdout);
		}
	}
	
	printf("#%02X%02X%02X - #%02X%02X%02X\n", R1, G1, B1, x[1][0], x[1][1], x[1][2]);
}

void calculate()
{
	double 	rl1, rl2, minimum = DBL_MAX,
			target_ratio = getRatio();
	int 	R1 = getChannelValue("R1"),
			G1 = getChannelValue("G1"),
			B1 = getChannelValue("B1"),
			R2 = getChannelValue("R2"),
			G2 = getChannelValue("G2"),
			B2 = getChannelValue("B2");
	
	rl1 = relativeLuminance[G1][R1][B1];
	rl2 = relativeLuminance[G2][R2][B2];
	
	printf("#%02X%02X%02X\n", R1, G1, B1);
	
	printf("linearize(R1) = %.20f\n", linearize(R1));
	printf("linearize(G1) = %.20f\n", linearize(G1));
	printf("linearize(B1) = %.20f\n", linearize(B1));
	
	printf("Relative Luminance: %.20f\n", rl1);
	
	printf("#%02X%02X%02X\n", R2, G2, B2);
	
	printf("linearize(R2) = %.20f\n", linearize(R2));
	printf("linearize(G2) = %.20f\n", linearize(G2));
	printf("linearize(B2) = %.20f\n", linearize(B2));
	
	printf("Relative Luminance: %.20f\n", rl2);
	
	printf("Target Contrast Ratio: %.20f\n", target_ratio);
	
	printf("Calculated Contrast Ratio: %.20f\n", (rl1 > rl2) ? (rl1 / rl2) : (rl2 / rl1));
	
	printf("objective: %.20f\n", rl1 - target_ratio * rl2);
}

void findMidway()
{
	double 	rl1, rl2, ratio, 
			minimum = DBL_MAX, current;
	char 	R1 = getChannelValue("R1"),
			G1 = getChannelValue("G1"),
			B1 = getChannelValue("B1"),
			R2 = getChannelValue("R2"),
			G2 = getChannelValue("G2"),
			B2 = getChannelValue("B2"),
			/*\                                 /*\
			|*|---------------------------------|*|
			|*| x[0][0] - index 0  x[1][0] - G3 |*|
			|*| x[0][1] - index 1  x[1][1] - R3 |*|
			|*| x[0][2] - index 2  x[1][2] - B3 |*|
			|*|---------------------------------|*|
			\*/                               /*\*/
			x[2][3];
		
	if(relativeLuminance[G1][R1][B1] < relativeLuminance[G2][R2][B2])
	{
		swap(&G1, &G2);
		swap(&R1, &R2);
		swap(&B1, &B2);
	}
	
	rl1 = (relativeLuminance[G1][R1][B1] - relativeLuminance[G2][R2][B2]) / 2;	
	for (x[0][0] = G1; x[0][0] < G1+1; x[0][0]++) 
	for (x[0][1] = R1; x[0][1] < R1+1; x[0][1]++) 
	for (x[0][2] = B1; x[0][2] < B1+1; x[0][2]++)
	{
		current = fabs(rl1 - relativeLuminance[x[0][0]][x[0][1]][x[0][2]]);
		
		if (minimum > current)
		{
			minimum = current;
			for(int i = 0; i < 3; i++)
				x[1][i] = x[0][i];
		}
	}
	
	printf("\n#%02X%02X%02X\n", x[1][0], x[1][1], x[1][2]);
}


void getGradient(char G1, char R1, char B1, char G2, char R2, char B2)
{
	if(G1 == G2 && R1 == R2 && B1 == B2)
		return;
		
	double 	rl1, rl2, ratio, 
			minimum = DBL_MAX, current;
			/*\                                 /*\
			|*|---------------------------------|*|
			|*| x[0][0] - index 0  x[1][0] - G3 |*|
			|*| x[0][1] - index 1  x[1][1] - R3 |*|
			|*| x[0][2] - index 2  x[1][2] - B3 |*|
			|*|---------------------------------|*|
			\*/                               /*\*/
	char	x[2][3];
	
	if(relativeLuminance[G1][R1][B1] < relativeLuminance[G2][R2][B2])
	{
		swap(&G1, &G2);
		swap(&R1, &R2);
		swap(&B1, &B2);
	}
	
	rl1 = (relativeLuminance[G1][R1][B1] - relativeLuminance[G2][R2][B2]) / 2;	
	for (x[0][0] = G1; x[0][0] < G1+1; x[0][0]++) 
	for (x[0][1] = R1; x[0][1] < R1+1; x[0][1]++) 
	for (x[0][2] = B1; x[0][2] < B1+1; x[0][2]++)
	{
		current = fabs(rl1 - relativeLuminance[x[0][0]][x[0][1]][x[0][2]]);
		
		if (minimum > current)
		{
			minimum = current;
			for(int i = 0; i < 3; i++)
				x[1][i] = x[0][i];
		}
	}
	
	printf("\n#%02X%02X%02X\n", x[1][0], x[1][1], x[1][2]);
	
	
	if(G1 == x[1][0] && R1 == x[1][1] && B1 == x[1][2])
		return;
	if(G2 == x[1][0] && R2 == x[1][1] && B2 == x[1][2])
		return;
		
	getGradient(G1, R1, B1, x[1][0], x[1][1], x[1][2]);
	getGradient(x[1][0], x[1][1], x[1][2], G2, R2, B2);
}

int main(void)
{
	int choice = 9;
	while(choice != 0)
	{
		system("cls");
		printf("0: Exit\n");
		printf("1: Find the 2 colors whose contrast ratio      \n");
		printf("   is as close to the given one as possible    \n");
		printf("2: Find the color that will make the contrast  \n");
		printf("   ratio with the given color as close to      \n");
		printf("   the given contrast ratio as possible        \n");
		printf("3: Calculate the relative luminances and       \n");
		printf("   contrast ratio of the 2 given colors        \n");
		printf("4: Find the 3rd color whose relative           \n");
		printf("   luminance is in between that of             \n");
		printf("   the 2 given colors                          \n");
		printf("5: Get colors for a smooth gradient            \n");
		printf("   between the 2 given colors                  \n");
		printf("Enter your choice: ");
		scanf("%d", &choice);
		
		// Memoize relative luminances to reduce
		// the contrast ratio formula into
		// look-up, subtraction, and divide
		memoize();
		
		switch(choice)
		{
			case 0:
				break;
			case 1:
				findPair();
				break;
			case 2:
				findPartner();
				break;
			case 3:
				calculate();
				break;
			case 4:
				findMidway();
				break;
			case 5:
				char 	R1 = getChannelValue("R1"),
						G1 = getChannelValue("G1"),
						B1 = getChannelValue("B1"),
						R2 = getChannelValue("R2"),
						G2 = getChannelValue("G2"),
						B2 = getChannelValue("B2");
				printf("%d, %d, %d, ", R1, G1, B1);
				getGradient(G1, R1, B1, G2, R2, B2);
				printf("%d, %d, %d, ", R2, G2, B2);
				break;
		}
		system("pause");
	}
}
