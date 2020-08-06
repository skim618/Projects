/* ENGGEN131 Project - C Project - 2019 */
/* The Warehouse */

/* 
Name: Sooyong Kim
UPI: skim618
ID: 573360750
 */

#include "project.h"

/*
The TimeWorked function calculates the length of a worker’s shift, by computing 
the number of seconds that have elapsed between them clocking in and clocking out.
The function will take four inputs (Clock in - minutes and seconds and Clock out - minutes and secodns).
The function will compute and return how much time has elapsed, in seconds, between those
two times.
Author: Sooyong Kim
*/
int TimeWorked(int minuteA, int secondA, int minuteB, int secondB)
{
	int totalSecondsA, totalSecondsB, secondsWorked;
	totalSecondsA = secondA + (minuteA * 60);
	totalSecondsB = secondB + (minuteB * 60);
	if (totalSecondsA > totalSecondsB) {
		secondsWorked = totalSecondsA - totalSecondsB;
	} else if (totalSecondsA < totalSecondsB) {
		secondsWorked = totalSecondsB - totalSecondsA;
	} else {
		secondsWorked = 0;
	}
	return secondsWorked;	
}

/*
The WarehouseAddress function takes an input of the maximum limit of a prime number. 
The function calculates the prime number closest but not exact to the limit.
The function uses a series of loops with conditional statements starting from the limit
downwards. The function will return the closest prime number.
Author: Sooyong Kim
*/
int WarehouseAddress(int maximum)
{
	int i, prime, count, number;
	number = maximum - 1;
	prime = 0;
	
	while (number > 0) {
		count = 0;
		for (i = 1; i <= number; i++) {
			if ((number % i) == 0) {
				count = count + 1;
			}
		}
		if (count == 2) {
			prime = number;
			break;
		}
	number--;	
	}	
	return prime;
}

/*
The Advertise function primarily works like the sliding LED displays people see
in shops in cities. This function Takes the first letter and brings it to 
the end, while shifting all the other letters one position to the left.
This function achives this by a use of a while loop and and storing the 
value of a specific position by the value to its left. A 'temp' (temporary)
variable is used to store first value initially which will be then placed at the end
at the very end of the loop.
Author: Sooyong Kim
*/
void Advertise(char *words)
{
	int i;
	i = 0;
	char temp;
	temp = 'i';
	temp = words[0];
	while (words[i] != '\0') {
		if (words[i+1] != '\0') {
			words[i] = words[i+1];
		} else {
			words[i] = temp;
		}
		i++;
	}
}

/*
The WinningBid function takes an array of values and its length as an input,
and will return the smallest DISTINCT number in that array. This function 
computes this by sorting the array from smallest to largest using bubble sort.
The function will then analyse every element and return the first distinct value 
in the array using another nested for loop.
Author: Sooyong Kim
*/
int WinningBid(int *values, int length)
{
	int i = 0;
	int j = 0;
	for (j = 0; j < length; j++) {
		for (i = 0; i+1 < length; i++) {
			int temp;
		    if (values[i] > values[i+1]) {
			temp = values[i];
			values[i] = values[i+1];
			values[i+1] = temp;
			}
		}
	}
	
	int result, count;
	result = -1;
	for (i = 0; i < length; i++) {
		count = 0;
		for (j = 0; j < length; j++) {
			if ((values [i] == values[j])&&(i != j)) {
				count = 1;
				break;
			}
		}
		if (!count) {
			result = values[i];
			break;
		}
	}
	return result;
}

/*
The BoxDesign function returns a 1D character array called 'design' that resembles a box 
structure of dimensions 'width' and 'height'. The Function returns the design of the box with
the middle of the box marked with an X.
The function computes this by first designing the box in a 2D character array called 'Space' 
then converting the 2D array into the 1D 'design' array.
Author: Sooyong Kim
*/
void BoxDesign(char *design, int width, int height)
{
	int maxLength, evenWidth, oddWidth, evenHeight, oddHeight;
	char Space[MAX_ARRAY_SIZE][MAX_ARRAY_SIZE];
	/* Obtain MaxLength of the 1D Design Array, and categorise the width and height
	   into even or odd
	*/
	evenWidth = 0;
	evenHeight = 0;
	oddHeight = 0;
	oddWidth = 0;
	maxLength = (width + 1) * height;
	if ((width % 2) == 0) {
		evenWidth = width;
	} else {
		oddWidth = width;
	}
	if ((height % 2) == 0) {
		evenHeight = height;
	} else {
		oddHeight = height;
	}

	// Form the box structure of the design
	int i, j;
	for (i = 0; i < height; i++) {
		for (j = 0; j <= width; j++) {
			if (j == width) {
				Space[i][j] = 10;
			} else {
			Space[i][j] = 42;
			}
		}
	}
	// Add the blank spaces only if width and height is above 2
	if ((width > 2) && (height > 2)) {
		for (i = 1; i < height - 1; i++) {
			for (j = 1; j < width - 1; j++) {
			Space[i][j] = 32;
			}
		}
		// According to the dimensions of the design, add the middle X
		if ((width == oddWidth) && (height == oddHeight)) {
		Space[height/2][width/2] = 88;
		} else if ((width == oddWidth) && (height == evenHeight)) {
		Space[height/2][width/2] = 88;
		Space[(height/2)-1][width/2] = 88;
		} else if ((width == evenWidth) && (height == oddHeight)) {
		Space[height/2][width/2] = 88;
		Space[height/2][(width/2)-1] = 88;
		} else if ((width == evenWidth) && (height == evenHeight)) {
		Space[height/2][width/2] = 88;
		Space[height/2][(width/2)-1] = 88;
		Space[(height/2)-1][width/2] = 88;
		Space[(height/2)-1][(width/2)-1] = 88;
		}
	} 
	// Convert the 2D Array Space into the original 1D char array design.
	int count = 0;

	for (i = 0; i < height; i++) {
		for (j = 0; j <= width; j++) {
			design[count] = Space[i][j];
			count++;
		}
	}
}

/*
The WorkerRoute function takes in a 2D array of values as an input. 
The 1 in the array will represent the worker, and the 2 will represent
the destination of the worker. This function computes the fastest route for the worker,
by moving horizontally and veritcally only. The function will pinpoint the locations
of numbers 1 and 2, and move the worker horizontally until it aligns it self 
with the horizontal (column) coordinate of the destination, and will move 
vertically until it reaches the destination. The function uses a series of while loops.
Author: Sooyong Kim
*/
void WorkerRoute(int warehouse[10][10])
{
	int i, j, row1, col1, row2, col2;
	row1 = 0;
	col1 = 0;
	row2 = 0;
	col2 = 0; 
	for (i = 0; i < 10; i++) {
		for (j = 0; j < 10; j++) {
			if (warehouse[i][j] == 1) {
				row1 = i;
				col1 = j;
			} else if (warehouse[i][j] == 2) {
				row2 = i;
				col2 = j;
			}
		}
	}
	if (col1 < col2) {
		while (col1 < col2) {
			warehouse[row1][col1+1] = 3;
			col1++;
		} 
	} else {
		while (col1 > col2) {
			warehouse[row1][col1-1] = 3;
			col1--;
		}
	}
	if (row1 < row2) {
		while (row1 < row2) {
			warehouse[row1+1][col1] = 3;
			row1++;
		} 
	} else {
		while (row1 > row2) {
			warehouse[row1-1][col1] = 3;
			row1--;
		}
	}
	warehouse[row2][col2] = 2;
}

/*
The MakeMove function takes in a 2D value of arrays which represnts a warehouse map, and 
either of 4 string inputs, w(up), d(left), s(down), and a(right). In the map there is a worker,
defined as number 5, and will move accordingly to the command string inputted. 
There are many different constraints for this movement, and the objective of this function is to get the 
worker moving the boxes(defined as the number 3) in to its targets(defined as 2). The worker may not move 
in front of a wall, or any obstacles that block its path unless if its a box. 
The funciton uses many different conditional statements to account for all the different possibilities 
of movement and restrictions there are.
Author: Sooyong Kim
*/
int MakeMove(int warehouse[10][10], char move)
{
	// Determine the position of the worker or a worker on a target and record it.
	int i, j; 
	int rowWorker = 0;
	int colWorker = 0;
	for (i = 0; i < 10; i++) {
		for (j = 0; j < 10; j++) {
			if ((warehouse[i][j] == 5) || (warehouse[i][j] == 6)) {
				rowWorker = i;
				colWorker = j;
			}
		}
	}
	//Initialise temp variable to store the position of the worker to swap with the upcoming spaces.
	int temp;
	temp = warehouse[rowWorker][colWorker]; 
	//MOVING RIGHT!!!
	if ((move == 'd') && (warehouse[rowWorker][colWorker+1] != 1)) {
		if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker][colWorker+1] == 3) && (warehouse[rowWorker][colWorker+2] == 0)) {
			warehouse[rowWorker][colWorker+2] = warehouse[rowWorker][colWorker+1];
			warehouse[rowWorker][colWorker] = 0;
			warehouse[rowWorker][colWorker+1] = temp;
		} else if ((warehouse[rowWorker][colWorker+1] == 3) && (warehouse[rowWorker][colWorker+2] == 1)) {
		return 0;
		} else if ((warehouse[rowWorker][colWorker+1] == 3) && (warehouse[rowWorker][colWorker+2] == 3)) {
		return 0;
		} else if ((warehouse[rowWorker][colWorker+1] == 3) && (warehouse[rowWorker][colWorker+2] == 4)) {
		return 0;
		} else if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker][colWorker+1] == 3) && (warehouse[rowWorker][colWorker+2] == 2)) {
			warehouse[rowWorker][colWorker+2] = 4;
			warehouse[rowWorker][colWorker] = 0;
			warehouse[rowWorker][colWorker+1] = temp;
		} else if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker][colWorker+1] == 2)) {
			warehouse[rowWorker][colWorker+1] = 6;
			warehouse[rowWorker][colWorker] = 0;
		} else if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker][colWorker+1] == 4) && (warehouse[rowWorker][colWorker+2] == 0)) {
			warehouse[rowWorker][colWorker+2] = 3;
			warehouse[rowWorker][colWorker+1] = 6;
			warehouse[rowWorker][colWorker] = 0;
		} else if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker][colWorker+1] == 4) && (warehouse[rowWorker][colWorker+2] == 2)) {
			warehouse[rowWorker][colWorker+2] = 4;
			warehouse[rowWorker][colWorker+1] = 6;
			warehouse[rowWorker][colWorker] = 0;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker][colWorker+1] == 4) && (warehouse[rowWorker][colWorker+2] == 2)) {
			warehouse[rowWorker][colWorker+2] = 4;
			warehouse[rowWorker][colWorker+1] = 6;
			warehouse[rowWorker][colWorker] = 2;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker][colWorker+1] == 4) && (warehouse[rowWorker][colWorker+2] == 0)) {
			warehouse[rowWorker][colWorker+2] = 3;
			warehouse[rowWorker][colWorker+1] = 6;
			warehouse[rowWorker][colWorker] = 2;
		} else if ((warehouse[rowWorker][colWorker+1] == 4) && (warehouse[rowWorker][colWorker+2] == 1)) {
			return 0;
		} else if ((warehouse[rowWorker][colWorker+1] == 4) && (warehouse[rowWorker][colWorker+2] == 3)) {
			return 0;
		} else if ((warehouse[rowWorker][colWorker+1] == 4) && (warehouse[rowWorker][colWorker+2] == 4)) {
			return 0;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker][colWorker+1] == 3) && (warehouse[rowWorker][colWorker+2] == 0)) {
			warehouse[rowWorker][colWorker+2] = 3;
			warehouse[rowWorker][colWorker+1] = 5;
			warehouse[rowWorker][colWorker] = 2;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker][colWorker+1] == 3) && (warehouse[rowWorker][colWorker+2] == 2)) {
			warehouse[rowWorker][colWorker+2] = 4;
			warehouse[rowWorker][colWorker+1] = 5;
			warehouse[rowWorker][colWorker] = 2;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker][colWorker+1] == 0)) {
			warehouse[rowWorker][colWorker+1] = 5;
			warehouse[rowWorker][colWorker] = 2;
		} else {
		warehouse[rowWorker][colWorker] = warehouse[rowWorker][colWorker+1];
		warehouse[rowWorker][colWorker+1] = temp;
		}
	//MOVING UP!!!
	} else if ((move == 'w') && (warehouse[rowWorker-1][colWorker] != 1)) {
		if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker-1][colWorker] == 3) && (warehouse[rowWorker-2][colWorker] == 0)) {
			warehouse[rowWorker-2][colWorker] = warehouse[rowWorker-1][colWorker];
			warehouse[rowWorker][colWorker] = 0;
			warehouse[rowWorker-1][colWorker] = temp;
		} else if ((warehouse[rowWorker-1][colWorker] == 3) && (warehouse[rowWorker-2][colWorker] == 1)) {
		return 0;
		} else if ((warehouse[rowWorker-1][colWorker] == 3) && (warehouse[rowWorker-2][colWorker] == 3)) {
		return 0;
		} else if ((warehouse[rowWorker-1][colWorker] == 3) && (warehouse[rowWorker-2][colWorker] == 4)) {
		return 0;
		} else if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker-1][colWorker] == 3) && (warehouse[rowWorker-2][colWorker] == 2)) {
			warehouse[rowWorker-2][colWorker] = 4;
			warehouse[rowWorker][colWorker] = 0;
			warehouse[rowWorker-1][colWorker] = temp;
		} else if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker-1][colWorker] == 2)) {
			warehouse[rowWorker-1][colWorker] = 6;
			warehouse[rowWorker][colWorker] = 0;
		} else if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker-1][colWorker] == 4) && (warehouse[rowWorker-2][colWorker] == 0)) {
			warehouse[rowWorker-2][colWorker] = 3;
			warehouse[rowWorker-1][colWorker] = 6;
			warehouse[rowWorker][colWorker] = 0;
		} else if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker-1][colWorker] == 4) && (warehouse[rowWorker-2][colWorker] == 2)) {
			warehouse[rowWorker-2][colWorker] = 4;
			warehouse[rowWorker-1][colWorker] = 6;
			warehouse[rowWorker][colWorker] = 0;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker-1][colWorker] == 4) && (warehouse[rowWorker-2][colWorker] == 2)) {
			warehouse[rowWorker-2][colWorker] = 4;
			warehouse[rowWorker-1][colWorker] = 6;
			warehouse[rowWorker][colWorker] = 2;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker-1][colWorker] == 4) && (warehouse[rowWorker-2][colWorker] == 0)) {
			warehouse[rowWorker-2][colWorker] = 3;
			warehouse[rowWorker-1][colWorker] = 6;
			warehouse[rowWorker][colWorker] = 2;
		} else if ((warehouse[rowWorker-1][colWorker] == 4) && (warehouse[rowWorker-2][colWorker] == 1)) {
			return 0;
		} else if ((warehouse[rowWorker-1][colWorker] == 4) && (warehouse[rowWorker-2][colWorker] == 3)) {
			return 0;
		} else if ((warehouse[rowWorker-1][colWorker] == 4) && (warehouse[rowWorker-2][colWorker] == 4)) {
			return 0;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker-1][colWorker] == 3) && (warehouse[rowWorker-2][colWorker] == 0)) {
			warehouse[rowWorker-2][colWorker] = 3;
			warehouse[rowWorker-1][colWorker] = 5;
			warehouse[rowWorker][colWorker] = 2;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker-1][colWorker] == 3) && (warehouse[rowWorker-2][colWorker] == 2)) {
			warehouse[rowWorker-2][colWorker] = 4;
			warehouse[rowWorker-1][colWorker] = 5;
			warehouse[rowWorker][colWorker] = 2;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker-1][colWorker] == 0)) {
			warehouse[rowWorker-1][colWorker] = 5;
			warehouse[rowWorker][colWorker] = 2;
		} else {
		warehouse[rowWorker][colWorker] = warehouse[rowWorker-1][colWorker];
		warehouse[rowWorker-1][colWorker] = temp;
		}
	// MOVING LEFT!!
	} else if ((move == 'a') && (warehouse[rowWorker][colWorker-1] != 1)) {
		if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker][colWorker-1] == 3) && (warehouse[rowWorker][colWorker-2] == 0)) {
			warehouse[rowWorker][colWorker-2] = warehouse[rowWorker][colWorker-1];
			warehouse[rowWorker][colWorker] = 0;
			warehouse[rowWorker][colWorker-1] = temp;
		} else if ((warehouse[rowWorker][colWorker-1] == 3) && (warehouse[rowWorker][colWorker-2] == 1)) {
		return 0;
		} else if ((warehouse[rowWorker][colWorker-1] == 3) && (warehouse[rowWorker][colWorker-2] == 3)) {
		return 0;
		} else if ((warehouse[rowWorker][colWorker-1] == 3) && (warehouse[rowWorker][colWorker-2] == 4)) {
		return 0;
		} else if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker][colWorker-1] == 3) && (warehouse[rowWorker][colWorker-2] == 2)) {
			warehouse[rowWorker][colWorker-2] = 4;
			warehouse[rowWorker][colWorker] = 0;
			warehouse[rowWorker][colWorker-1] = temp;
		} else if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker][colWorker-1] == 2)) {
			warehouse[rowWorker][colWorker-1] = 6;
			warehouse[rowWorker][colWorker] = 0;
		} else if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker][colWorker-1] == 4) && (warehouse[rowWorker][colWorker-2] == 0)) {
			warehouse[rowWorker][colWorker-2] = 3;
			warehouse[rowWorker][colWorker-1] = 6;
			warehouse[rowWorker][colWorker] = 0;
		} else if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker][colWorker-1] == 4) && (warehouse[rowWorker][colWorker-2] == 2)) {
			warehouse[rowWorker][colWorker-2] = 4;
			warehouse[rowWorker][colWorker-1] = 6;
			warehouse[rowWorker][colWorker] = 0;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker][colWorker-1] == 4) && (warehouse[rowWorker][colWorker-2] == 2)) {
			warehouse[rowWorker][colWorker-2] = 4;
			warehouse[rowWorker][colWorker-1] = 6;
			warehouse[rowWorker][colWorker] = 2;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker][colWorker-1] == 4) && (warehouse[rowWorker][colWorker-2] == 0)) {
			warehouse[rowWorker][colWorker-2] = 3;
			warehouse[rowWorker][colWorker-1] = 6;
			warehouse[rowWorker][colWorker] = 2;
		} else if ((warehouse[rowWorker][colWorker-1] == 4) && (warehouse[rowWorker][colWorker-2] == 1)) {
			return 0;
		} else if ((warehouse[rowWorker][colWorker-1] == 4) && (warehouse[rowWorker][colWorker-2] == 3)) {
			return 0;
		} else if ((warehouse[rowWorker][colWorker-1] == 4) && (warehouse[rowWorker][colWorker-2] == 4)) {
			return 0;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker][colWorker-1] == 3) && (warehouse[rowWorker][colWorker-2] == 0)) {
			warehouse[rowWorker][colWorker-2] = 3;
			warehouse[rowWorker][colWorker-1] = 5;
			warehouse[rowWorker][colWorker] = 2;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker][colWorker-1] == 3) && (warehouse[rowWorker][colWorker-2] == 2)) {
			warehouse[rowWorker][colWorker-2] = 4;
			warehouse[rowWorker][colWorker-1] = 5;
			warehouse[rowWorker][colWorker] = 2;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker][colWorker-1] == 0)) {
			warehouse[rowWorker][colWorker-1] = 5;
			warehouse[rowWorker][colWorker] = 2;
		} else {
		warehouse[rowWorker][colWorker] = warehouse[rowWorker][colWorker-1];
		warehouse[rowWorker][colWorker-1] = temp;
		}
	//MOVING DOWN!!
	} else if ((move == 's') && (warehouse[rowWorker+1][colWorker] != 1)) {
		if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker+1][colWorker] == 3) && (warehouse[rowWorker+2][colWorker] == 0)) {
			warehouse[rowWorker+2][colWorker] = warehouse[rowWorker+1][colWorker];
			warehouse[rowWorker][colWorker] = 0;
			warehouse[rowWorker+1][colWorker] = temp;
		} else if ((warehouse[rowWorker+1][colWorker] == 3) && (warehouse[rowWorker+2][colWorker] == 1)) {
		return 0;
		} else if ((warehouse[rowWorker+1][colWorker] == 3) && (warehouse[rowWorker+2][colWorker] == 3)) {
		return 0;
		} else if ((warehouse[rowWorker+1][colWorker] == 3) && (warehouse[rowWorker+2][colWorker] == 4)) {
		return 0;
		} else if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker+1][colWorker] == 3) && (warehouse[rowWorker+2][colWorker] == 2)) {
			warehouse[rowWorker+2][colWorker] = 4;
			warehouse[rowWorker][colWorker] = 0;
			warehouse[rowWorker+1][colWorker] = temp;
		} else if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker+1][colWorker] == 2)) {
			warehouse[rowWorker+1][colWorker] = 6;
			warehouse[rowWorker][colWorker] = 0;
		} else if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker+1][colWorker] == 4) && (warehouse[rowWorker+2][colWorker] == 0)) {
			warehouse[rowWorker+2][colWorker] = 3;
			warehouse[rowWorker+1][colWorker] = 6;
			warehouse[rowWorker][colWorker] = 0;
		} else if ((warehouse[rowWorker][colWorker] == 5) && (warehouse[rowWorker+1][colWorker] == 4) && (warehouse[rowWorker+2][colWorker] == 2)) {
			warehouse[rowWorker+2][colWorker] = 4;
			warehouse[rowWorker+1][colWorker] = 6;
			warehouse[rowWorker][colWorker] = 0;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker+1][colWorker] == 4) && (warehouse[rowWorker+2][colWorker] == 2)) {
			warehouse[rowWorker+2][colWorker] = 4;
			warehouse[rowWorker+1][colWorker] = 6;
			warehouse[rowWorker][colWorker] = 2;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker+1][colWorker] == 4) && (warehouse[rowWorker+2][colWorker] == 0)) {
			warehouse[rowWorker+2][colWorker] = 3;
			warehouse[rowWorker+1][colWorker] = 6;
			warehouse[rowWorker][colWorker] = 2;
		} else if ((warehouse[rowWorker+1][colWorker] == 4) && (warehouse[rowWorker+2][colWorker] == 1)) {
			return 0;
		} else if ((warehouse[rowWorker+1][colWorker] == 4) && (warehouse[rowWorker+2][colWorker] == 3)) {
			return 0;
		} else if ((warehouse[rowWorker+1][colWorker] == 4) && (warehouse[rowWorker+2][colWorker] == 4)) {
			return 0;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker+1][colWorker] == 3) && (warehouse[rowWorker+2][colWorker] == 0)) {
			warehouse[rowWorker+2][colWorker] = 3;
			warehouse[rowWorker+1][colWorker] = 5;
			warehouse[rowWorker][colWorker] = 2;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker+1][colWorker] == 3) && (warehouse[rowWorker+2][colWorker] == 2)) {
			warehouse[rowWorker+2][colWorker] = 4;
			warehouse[rowWorker+1][colWorker] = 5;
			warehouse[rowWorker][colWorker] = 2;
		} else if ((warehouse[rowWorker][colWorker] == 6) && (warehouse[rowWorker+1][colWorker] == 0)) {
			warehouse[rowWorker+1][colWorker] = 5;
			warehouse[rowWorker][colWorker] = 2;
		} else {
		warehouse[rowWorker][colWorker] = warehouse[rowWorker+1][colWorker];
		warehouse[rowWorker+1][colWorker] = temp;
		}
	} 
	// Count how many objectives there are left, and return 1 ONLY if 
	// there are no more targets and the worker is on at least one target
	int countObjctves = 0;
	int count4 = 0;
	int count6 = 0;
	for (i = 0; i < 10; i++) {
		for (j = 0; j < 10; j++) {
			if (warehouse[i][j] == 2) {
				countObjctves++;
			} else if (warehouse[i][j] == 4) {
				count4++;
			} else if (warehouse[i][j] == 6) {
				count6++;
			}
		}
	}
	if ((countObjctves == 0) && (count4 == 0) ){
		return 0;
	} else if ((countObjctves == 0) && (count6 == 1)) {
		return 1;
	} else return 0;
}

