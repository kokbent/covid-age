CFLAGS=--ansi --pedantic -O2 -std=c++11 
#CFLAGS=--ansi --pedantic -g 
INCLUDE= -I../../src/
LDFLAGS=  ../../src/*.o

model: chicago_yr1.cpp
	g++ $(CFLAGS) chicago_yr1.cpp $(INCLUDE) $(LDFLAGS) -o model

