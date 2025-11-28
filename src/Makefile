TARGET=main
SRCFILE=main.cpp
CXXFLAGS=-std=c++23

all: $(TARGET)

$(TARGET): $(SRCFILE)
	g++ -o $(TARGET) $(SRCFILE) $(CXXFLAGS)

run: $(TARGET)
	./$(TARGET)
