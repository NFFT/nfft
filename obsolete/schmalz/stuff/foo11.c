int foo_private();

int foo() {
	return foo_private();
}

int foo_special() {
	return 0;
}
