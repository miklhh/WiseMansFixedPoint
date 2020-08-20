CC = g++
CFLAGS = -I./ -std=c++14 -O3 -Wall -Wextra -Wpedantic -Weffc++

OBJS=tests/test.o tests/catch.o
HEADER=FixedPoint.h

%.o: %.cc
	$(CC) $(CFLAGS) -c $^ -o $@

.PHONY: run_test clean

run_test: $(SRC) tests/catch_test.out
	@tests/catch_test.out

tests/catch_test.out: $(HEADER) $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o tests/catch_test.out

tests/test.o: $(HEADER) tests/test.cc
	$(CC) $(CFLAGS) -c tests/test.cc -o tests/test.o

clean:
	-@rm -v tests/catch.o
	-@rm -v tests/catch_test.out
	-@rm -v tests/test.o
