CC = g++
CFLAGS = -I./ -std=c++17 -march=native -O3 -Wall -Wextra -Wpedantic -Weffc++
CFLAGS_COV = -I./ -std=c++17 --coverage -g -O0 -Wall -Wextra -Wpedantic -Weffc++

OBJS=tests/test.o tests/catch.o
FIXP_HEADER=FixedPoint.h

%.o: %.cc
	$(CC) $(CFLAGS) -c $^ -o $@

.PHONY: run_test run_test_with_coverage clean

run_test: tests/catch_test.out
	@tests/catch_test.out

run_test_with_coverage: $(FIXP_HEADER)
	$(CC) $(CFLAGS_COV) tests/test.cc tests/catch.cc -o tests/catch_test_cov.out
	tests/catch_test_cov.out
	@gcovr -r . --html --html-details -o coverage.html
	@echo "Generated coverage files, entry 'coverage.html'."

tests/catch_test.out: $(FIXP_HEADER) $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o tests/catch_test.out

tests/test.o: $(FIXP_HEADER) tests/test.cc
	$(CC) $(CFLAGS) -c tests/test.cc -o tests/test.o

clean:
	-@rm -v tests/*.out
	-@rm -v tests/*.o
	-@rm -v tests/*.gcda
	-@rm -v *.html
	-@rm -v *.gcno
	-@rm -v *.gcda
