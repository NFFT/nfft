#include "foo.h"
#include<stdio.h>

int foo_special();

int main() {
	printf("%d %d\n", foo(), foo_special());
	return 0;
}
