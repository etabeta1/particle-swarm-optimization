CXX=g++

# Contains all the sources to run tests
TEST_FOLDER=tests

# Contains the main source code of the library
SWARM_FOLDER=src

# Contains the main application source code (e.g. reproducible simulations, demos, etc.)
APP_FOLDER=app

# Contains compiled binaries
BINARY_FOLDER=bin

# Contains generated documentation
DOCS_FOLDER=docs

COMMON_FLAGS=-I$(SWARM_FOLDER) -fopenmp -Wall -Wextra -std=c++23 -march=native

TEST_FLAGS=$(COMMON_FLAGS) -g -O0 -fsanitize=address,undefined
DIST_FLAGS=$(COMMON_FLAGS) -O2

all: test docs dist

test:
	mkdir -p $(BINARY_FOLDER)
	for test_file in $(TEST_FOLDER)/*.cpp; do \
		$(CXX) -o $(BINARY_FOLDER)/$$(basename $$test_file .cpp) $$test_file $(TEST_FLAGS); \
		$(BINARY_FOLDER)/$$(basename $$test_file .cpp); \
	done

clean:
	rm -rf $(BINARY_FOLDER) $(DOCS_FOLDER)/html $(DOCS_FOLDER)/latex

dist:
	mkdir -p $(BINARY_FOLDER)
	for app_file in $(APP_FOLDER)/*.cpp; do \
		$(CXX) -o $(BINARY_FOLDER)/$$(basename $$app_file .cpp) $$app_file $(DIST_FLAGS); \
	done

docs:
	doxygen Doxyfile
