#include<cfloat>
#include<cmath>
#include<cstdint>
#include<cstdio>
#include<ctime>

#include<map>
#include<numeric>
#include<vector>
#include<string>

std::map<std::tuple<uint8_t, uint8_t, uint8_t>, double> RelativeLuminances_1; 
std::map<double, std::tuple<uint8_t, uint8_t, uint8_t>> RelativeLuminances_2;

double 
linearize(
	int c
)
{
	return (c < 11) ? c / 3294.6 : std::pow((c + 14.025) / 269.025, 2.4);
}

void 
memoize()
{
	std::vector<double> 
	values(256,0.0),
	coefficients = { 0.7152, 0.2126, 0.0722, 1.000 }, 
	linearized(4,0.0);
	
	std::tuple<uint8_t,uint8_t,uint8_t> 
	color;
	
	double rl;
	
	// Linearize
	for (int i = 0; i < 11; i++)
		values[i] = i / 3294.6;
	for (int i = 11; i < 256; i++)
		values[i] = pow((i + 14.025) / 269.025, 2.4);
	
	// Memoize relative luminance
	for (int i = 0; i < 256; i++)
	for (int j = 0; j < 256; j++)
	for (int k = 0; k < 256; k++)
	{
		linearized = {values[i], values[j], values[k], 0.05};
		color = std::make_tuple((uint8_t)i,(uint8_t)j,(uint8_t)k);
		rl = std::inner_product(
			 coefficients.begin()
			,coefficients.end()
			,linearized.begin()
			,0.0
		);
		RelativeLuminances_1[color] = rl;
		RelativeLuminances_2[rl] = color;
	}
}

uint8_t 
getChannelValue(
	std::string channel
)
{
	uint8_t result;
	
	do {
		printf("%s: ", channel);
		scanf("%x", &result);
	} while(0 > result || result > 255);
	
	return result;
}

struct Color_Pair{
	std::tuple<uint8_t, uint8_t, uint8_t> 
	first, 
	second;
};

Color_Pair 
getColors()
{
	uint8_t 	
	R1 = getChannelValue("R1"),
	G1 = getChannelValue("G1"),
	B1 = getChannelValue("B1"),
	R2 = getChannelValue("R2"),
	G2 = getChannelValue("G2"),
	B2 = getChannelValue("B2");
	fflush(stdin);

	Color_Pair colors;
	
	colors.first  = std::make_tuple(G1, R1, B1), 
	colors.second = std::make_tuple(G2, R2, B2);
	return colors;
}

std::tuple<uint8_t, uint8_t, uint8_t> 
getColor()
{
	uint8_t 	
	R1 = getChannelValue("R1"),
	G1 = getChannelValue("G1"),
	B1 = getChannelValue("B1");
	fflush(stdin);

	return std::make_tuple(G1, R1, B1);
}

void 
printColors(	 
	std::tuple<uint8_t,uint8_t,uint8_t> color1,
	std::tuple<uint8_t,uint8_t,uint8_t> color2
)
{
	printf(
		"#%02X%02X%02X - #%02X%02X%02X\n", 
		std::get<1>(color1), 
		std::get<0>(color1), 
		std::get<2>(color1),
		std::get<1>(color2), 
		std::get<0>(color2), 
		std::get<2>(color2)
	);
	fflush(stdout);
}

void 
printColor(	 
	std::tuple<uint8_t,uint8_t,uint8_t> color
)
{
	printf(
		"#%02X%02X%02X\n", 
		std::get<1>(color), 
		std::get<0>(color), 
		std::get<2>(color)
	);
	fflush(stdout);
}

double 
getRatio()
{
	double result;
	do {
		printf("Ratio (1 to 21 inclusive): ");
		scanf("%lf", &result);
	} while(1 > result || result > 21);
	
	return result;
}

void 
swap(
	std::tuple<uint8_t,uint8_t,uint8_t> *x, 
	std::tuple<uint8_t,uint8_t,uint8_t> *y
)
{
	std::tuple<uint8_t,uint8_t,uint8_t> temp = *x;
	*x = *y;
	*y = temp;
}

void 
findPair()
{
	// ratio ≥ rl1 / rl2
	// ratio * rl2 ≥ rl1 ≥ rl2
	// |ratio * rl2 - rl1| ≥ 0
	double 
	ratio = getRatio(), 
	minimum = DBL_MAX, 
	current, obj;
	
	std::map<double, std::tuple<uint8_t, uint8_t, uint8_t>>::iterator
	rl1 = RelativeLuminances_2.lower_bound(ratio * 0.05),
	rl2 = RelativeLuminances_2.begin(),
	rl2_upper = RelativeLuminances_2.lower_bound((1.05/ratio)),
	color1, 
	color2;
	
	// Checking til 1.05/ratio is enough 
	// because 1.05 is the maximum contrast ratio
	rl2_upper++;
	for(; rl2 != rl2_upper; rl2++)
	{
		// Ideally, only the color that has a relative luminance of ratio * rl2 should be checked
		// but ratio * rl2 may not be in RelativeLuminances_2,
		// so the colors that have relative luminances just above and below ratio * rl2 are checked.
		obj = ratio * (rl2->first);
		rl1 = RelativeLuminances_2.lower_bound(obj);
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
	
	printColors(
		rl1->second,
		rl2->second
	);
}

void 
findPartner(
	std::tuple<uint8_t, uint8_t, uint8_t> color1
)
{
	std::tuple<uint8_t,uint8_t,uint8_t> 
	color2;
	
	double 
	ratio = getRatio(),
	rl1 = RelativeLuminances_1[color1],
	rl2 = RelativeLuminances_2.lower_bound(rl1/ratio)->first,
	minimum = DBL_MAX,
	current;

	color2 = (fabs(rl2 * ratio - rl1) < fabs((rl2-1) * ratio - rl1)) ? 
			 RelativeLuminances_2[rl2] : 
			 RelativeLuminances_2[rl2-1];
	
	printColors(
		color1,
		color2
	);
}

void calculate(	 
	 std::tuple<uint8_t,uint8_t,uint8_t> color1
	,std::tuple<uint8_t,uint8_t,uint8_t> color2
)
{
	uint8_t 	
	R1 = std::get<1>(color1),
	G1 = std::get<0>(color1),
	B1 = std::get<2>(color1),
	R2 = std::get<1>(color2),
	G2 = std::get<0>(color2),
	B2 = std::get<2>(color2);

	double 	
	rl1 = RelativeLuminances_1[color1], 
	rl2 = RelativeLuminances_1[color2], 
	minimum = DBL_MAX,
	target_ratio = getRatio();
	
	printColor(color1);
	
	printf("linearize(R1) = %.20f\n", linearize(R1));
	printf("linearize(G1) = %.20f\n", linearize(G1));
	printf("linearize(B1) = %.20f\n", linearize(B1));
	
	printf("Relative Luminance: %.20f\n", rl1);
	
	printColor(color2);
	
	printf("linearize(R2) = %.20f\n", linearize(R2));
	printf("linearize(G2) = %.20f\n", linearize(G2));
	printf("linearize(B2) = %.20f\n", linearize(B2));
	
	printf("Relative Luminance: %.20f\n", rl2);
	
	printf("Target Contrast Ratio: %.20f\n", target_ratio);
	
	printf("Calculated Contrast Ratio: %.20f\n", (rl1 > rl2) ? (rl1 / rl2) : (rl2 / rl1));
	
	printf("objective: %.20f\n", fabs(rl1 - target_ratio * rl2));
}

std::tuple<uint8_t,uint8_t,uint8_t> 
findMidway(	 
	std::tuple<uint8_t,uint8_t,uint8_t> color1,
	std::tuple<uint8_t,uint8_t,uint8_t> color2
)
{
	if(RelativeLuminances_1[color1] < RelativeLuminances_1[color2])
		swap(
			&color1, 
			&color2
		);

	std::tuple<uint8_t,uint8_t,uint8_t>
	holder = RelativeLuminances_2
			 .lower_bound((RelativeLuminances_1[color1] - RelativeLuminances_1[color2]) / 2)
			 ->second;
	
	printColor(holder);
	return holder;
}

void 
getGradient(
	 std::tuple<uint8_t,uint8_t,uint8_t> color1
	,std::tuple<uint8_t,uint8_t,uint8_t> color2
)
{
	if(color1 == color2)
		return;

	std::tuple<uint8_t, uint8_t, uint8_t>
	holder = findMidway(color1, color2);
	
	if(color1 == holder || color2 == holder)
		return;
		
	getGradient(color1,holder);
	getGradient(holder, color2);
}

int 
main(
	void
)
{
	int 
	choice = 9;
	while(choice != 0)
	{
		system("cls");
		printf("0: Exit\n");
		printf("1: Find the 2 colors whose contrast ratio	  \n");
		printf("   is as close to the given one as possible   \n");
		printf("2: Find the color that will make the contrast \n");
		printf("   ratio with the given color as close to     \n");
		printf("   the given contrast ratio as possible       \n");
		printf("3: Calculate the relative luminances and      \n");
		printf("   contrast ratio of the 2 given colors	      \n");
		printf("4: Find the 3rd color whose relative          \n");
		printf("   luminance is in between that of            \n");
		printf("   the 2 given colors						  \n");
		printf("5: Get colors for a smooth gradient           \n");
		printf("   between the 2 given colors                 \n");
		printf("Enter your choice: ");
		scanf("%d", &choice);
		
		// Memoize relative luminances into 2 red-and-black trees to:
		// 1. simplify the objective formula, and
		// 2. reduce the amount of look-ups.
		memoize();
		
		switch(choice)
		{
			case 0:
				break;
			case 1:
				findPair();
				break;
			case 2:
				std::tuple<uint8_t, uint8_t, uint8_t>
				color = getColor();
				findPartner(color);
				break;
			case 3:
				Color_Pair colors = getColors();
				calculate(colors.first, colors.second);
				break;
			case 4:
				Color_Pair colors = getColors();
				findMidway(colors.first, colors.second);
				break;
			case 5:
				Color_Pair colors = getColors();
				printColor(colors.first);
				getGradient(colors.first, colors.second);
				printColor(colors.second);
				break;
		}
		system("pause");
	}
}
