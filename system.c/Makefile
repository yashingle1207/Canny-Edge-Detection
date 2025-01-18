SYSTEMC = /opt/pkg/systemc-2.3.1

INCLUDE = -I. -I$(SYSTEMC)/include
LIBRARY = -L$(SYSTEMC)/lib-linux64 -Wl,-rpath,$(SYSTEMC)/lib-linux64
CFLAGS = $(INCLUDE) $(LIBRARY) -lsystemc -O2 -std=c++03


CC = g++
RM = rm -f

TARGET	= hw3
OFILES = main.o stim.o mult.o mon.o

all: $(TARGET)

$(TARGET): $(OFILES)
	$(CC) $^ $(CFLAGS) -o $@

main.o: main.cpp top.h stim.h mult.h mon.h
	$(CC) -c main.cpp $(CFLAGS) -o $@

stim.o: stim.cpp stim.h
	$(CC) -c stim.cpp $(CFLAGS) -o $@

mon.o: mon.cpp mon.h 
	$(CC) -c mon.cpp $(CFLAGS) -o $@

mult.o: mult.cpp mult.h
	$(CC) -c mult.cpp $(CFLAGS) -o $@

clean:
	$(RM) *.o $(TARGET)

test: $(TARGET)
	./$(TARGET)

