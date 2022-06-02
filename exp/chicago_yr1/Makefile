CFLAGS=--ansi --pedantic -O2 -std=c++11 
#CFLAGS=--ansi --pedantic -g 
INCLUDE= -I../../src/
LDFLAGS=  ../../src/*.o

model: ../../src/Event_Driven_NUCOVID.h ../../src/Utility.h ../../src/Time_Series.h chicago_yr1.cpp
	g++ $(CFLAGS) chicago_yr1.cpp $(INCLUDE) $(LDFLAGS) -o model

