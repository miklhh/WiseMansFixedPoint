CC = g++
CFLAGS = -I./ -std=c++17 -O3 -march=native -Wall -Wextra -Wpedantic -Weffc++
CFLAGS_COV = -I./ -std=c++17 --coverage -g -O0 -Wall -Wextra -Wpedantic -Weffc++
FIXP_HEADER=FixedPoint.h

.PHONY: run_unit_test run_unit_test_with_coverage run_profile clean

%.o: %.cc
	$(CC) $(CFLAGS) -c $^ -o $@

run_unit_test: tests/catch_test.out
	@tests/catch_test.out

run_unit_test_with_coverage: $(FIXP_HEADER)
	$(CC) $(CFLAGS_COV) tests/test.cc tests/catch.cc -o tests/catch_test_cov.out
	tests/catch_test_cov.out
	@gcovr -r . --html --html-details -o coverage.html
	@echo "Generated coverage files, entry 'coverage.html'."

run_profile: tests/CoordinateDescend_AD_ML_CPPMEX_Q_Profile.out
	@tests/CoordinateDescend_AD_ML_CPPMEX_Q_Profile.out

OBJS=tests/test.o tests/catch.o
tests/catch_test.out: $(FIXP_HEADER) $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@

tests/CoordinateDescend_AD_ML_CPPMEX_Q_Profile.out: tests/CoordinateDescend_AD_ML_CPPMEX_Q_Profile.cc $(FIXP_HEADER)
	$(CC) $(CFLAGS) $< -o $@

tests/test.o: $(FIXP_HEADER) tests/test.cc
	$(CC) $(CFLAGS) -c tests/test.cc -o tests/test.o

clean:
	@rm -fv tests/*.out
	@rm -fv tests/*.o
	@rm -fv tests/*.gcda
	@rm -fv *.html
	@rm -fv *.gcno
	@rm -fv *.gcda
