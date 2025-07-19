#include<float.h>
#include<omp.h>
#include<stdbool.h>
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#include<cmath>
#include<map>
#include<numeric>
#include<vector>
#include<string>

std::vector<std::vector<std::vector<double> > > relative_luminance(
	256, std::vector<std::vector<double> >(
		256, std::vector<double>(
			256,0.0
		)
	)
);

std::map<double, std::tuple<uint8_t, uint8_t, uint8_t>> RelativeLuminances;

double linearize(int c)
{
	return (c < 11) ? c / 3294.6 : std::pow((c + 14.025) / 269.025, 2.4);
}

void memoize()
{
    std::vector<double> values(256,0.0);
    std::vector<double> coefficients = { 0.7152, 0.2126, 0.0722, 1.000 };
    std::vector<double> linearized_components(4,0.0);
	//memoize R, G, and B
	for (int i = 0; i < 11; i++)
		values[i] = i / 3294.6;
	for (int i = 11; i < 256; i++)
		values[i] = pow((i + 14.025) / 269.025, 2.4);
	
	// Memoize relative luminance
	for (int i = 0; i < 256; i++)
    for (int j = 0; j < 256; j++)
    for (int k = 0; k < 256; k++)
    {
        linearized_components = {values[i], values[j], values[k], 0.05};
        relative_luminance[i][j][k] =    std::inner_product(
             coefficients.begin()
            ,coefficients.end()
            ,linearized_components.begin()
            ,0
        );
		RelativeLuminances[relative_luminance[i][j][k]] = std::make_tuple((uint8_t)i,(uint8_t)j,(uint8_t)k);
    }
}

int getChannelValue(std::string channel)
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

struct Color_Pair{
    uint8_t G1; 
    uint8_t R1; 
    uint8_t B1;
    uint8_t G2; 
    uint8_t R2; 
    uint8_t B2;
};

void findPair()
{
	// rl1 / rl2 ≥ ratio
	// rl1 ≥ ratio * rl2
	// relative_luminance[x0][x1][x2] ≥ ratio * relative_luminance[x3][x4][x5]
	// relative_luminance[x0][x1][x2] - ratio * relative_luminance[x3][x4][x5] ≥ 0
	double ratio = getRatio(), minimum = DBL_MAX, current, obj;
    std::map<double, std::tuple<uint8_t, uint8_t, uint8_t>>::iterator
    rl1 = RelativeLuminances.lower_bound(ratio * 0.05),
    rl2 = RelativeLuminances.begin(),
    rl2lb = RelativeLuminances.lower_bound((1.05/ratio)),
    color1, color2;
	Color_Pair result;

    rl2lb++;
    for(; rl2 != rl2lb; rl2++)
    {
        obj = ratio * (rl2->first);
        rl1 = RelativeLuminances.lower_bound(obj);
        current = fabs(rl1->first - obj);
        if(minimum > current)
        {
            minimum = current;
            color1 = rl1;
            color2 = rl2;
        }
        rl1--;
        current = fabs(rl1->first - obj);
        if(minimum > current)
        {
            minimum = current;
            color1 = rl1;
            color2 = rl2;
        }
    }

	result.G1 = std::get<0>(rl1->second);
	result.R1 = std::get<1>(rl1->second);
	result.B1 = std::get<2>(rl1->second);
	result.G2 = std::get<0>(rl2->second);
	result.R2 = std::get<1>(rl2->second);
	result.B2 = std::get<2>(rl2->second);
	
	printf("#%02X%02X%02X - #%02X%02X%02X\n", result.R1, result.G1, result.B1, result.R2, result.G2, result.B2);
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
	
	rl1 = relative_luminance[R1][G1][B1];
	for(x[0][0] = 0; x[0][0] < 256; x[0][0]++) 
	for(x[0][1] = 0; x[0][1] < 256; x[0][1]++) 
	for(x[0][2] = 0; x[0][2] < 256; x[0][2]++)
	{	
		current = fabs(rl1 - ratio * relative_luminance[x[0][0]][x[0][1]][x[0][2]]);
		
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
	
	rl1 = relative_luminance[G1][R1][B1];
	rl2 = relative_luminance[G2][R2][B2];
	
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
		
	if(relative_luminance[G1][R1][B1] < relative_luminance[G2][R2][B2])
	{
		swap(&G1, &G2);
		swap(&R1, &R2);
		swap(&B1, &B2);
	}
	
	rl1 = (relative_luminance[G1][R1][B1] - relative_luminance[G2][R2][B2]) / 2;	
	for (x[0][0] = G1; x[0][0] < G1+1; x[0][0]++) 
	for (x[0][1] = R1; x[0][1] < R1+1; x[0][1]++) 
	for (x[0][2] = B1; x[0][2] < B1+1; x[0][2]++)
	{
		current = fabs(rl1 - relative_luminance[x[0][0]][x[0][1]][x[0][2]]);
		
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
	
	if(relative_luminance[G1][R1][B1] < relative_luminance[G2][R2][B2])
	{
		swap(&G1, &G2);
		swap(&R1, &R2);
		swap(&B1, &B2);
	}
	
	rl1 = (relative_luminance[G1][R1][B1] - relative_luminance[G2][R2][B2]) / 2;	
	for (x[0][0] = G1; x[0][0] < G1+1; x[0][0]++) 
	for (x[0][1] = R1; x[0][1] < R1+1; x[0][1]++) 
	for (x[0][2] = B1; x[0][2] < B1+1; x[0][2]++)
	{
		current = fabs(rl1 - relative_luminance[x[0][0]][x[0][1]][x[0][2]]);
		
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
