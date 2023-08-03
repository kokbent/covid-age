CFLAGS=--ansi --pedantic -O2 -std=c++11 
#CFLAGS=--ansi --pedantic -g 
INCLUDE= -I../../src/
LDFLAGS=  ../../src/*.o

model: ../../src/Event_Driven_NUCOVID.h ../../src/Utility.h ../../src/Time_Series.h chicago_yr1.cpp
	g++ $(CFLAGS) chicago_yr1.cpp $(INCLUDE) $(LDFLAGS) -o model

model_mod: ../../src/Event_Driven_NUCOVID.h ../../src/Utility.h ../../src/Time_Series.h chicago_yr1_mod.cpp
	g++ $(CFLAGS) chicago_yr1_mod.cpp $(INCLUDE) $(LDFLAGS) -o model_mod

model_cereal: ../../src/NUCOVID_cereal.h ../../src/Utility.h ../../src/Time_Series.h chicago_yr1_cereal.cpp
	g++ $(CFLAGS) chicago_yr1_cereal.cpp $(INCLUDE) $(LDFLAGS) -o model_cereal
