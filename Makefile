CXX = clang++ -std=c++14

INCLUDE = -I.

MAIN = test_rng

all: test_rng

test_rng: test_rng.C
	$(CXX) $(INCLUDE) $@.C -o $@

clean:
	$(RM) *~ test_rng
