TARGET=main
SRCFILE=main.cpp

all: $(TARGET)

$(TARGET): $(SRCFILE)
	g++ -o $(TARGET) $(SRCFILE) -fopenmp

run: $(TARGET)
	./$(TARGET)