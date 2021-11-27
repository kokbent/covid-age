#include "Time_Series.h"
#include "Event_Driven_NUCOVID.h"

int main() { 

    int N = 2e6;
    double Ki    = 0.14; // S -> E transition rate
    double Kasymp= 0.4/4.0; // E -> I (asymp) rate
    double Kpres = 0.6/4.0; // E -> I (symp) rate
    double Kmild = 0.92/3.0;
    double Kseve = 0.08/3.0;
    double Khosp = 1.0/4.5;
    double Kcrit = 1.0/5.0;
    double Kdeath = 1.0/5.0;
    vector<double> Krec = {1.0/9.0, 1.0/9.0, 1.0/5.0, 1.0/9.0, 1.0/2.5}; // I -> R transition rate
    double Pcrit = 0.55;
    double Pdeath = 0.2;
    vector<double> Pdetect = {0.2/6, 0.2/6, 0.2, 0.5}; // I -> R transition rate
    //vector<double> Pdetect = {0, 0, 0, 0}; // I -> R transition rate
    double immunity_duration = 10000;

    for(int i=0; i<1; i++ ) {
        Event_Driven_NUCOVID sim(N, Ki, Kasymp, Kpres, Kmild, Kseve, Khosp, Kcrit,
                                 Kdeath, Krec, Pcrit, Pdeath, Pdetect, immunity_duration);
        sim.rng.seed(time(0)); // this simulator has its own RNG which must be seeded as well
        sim.rand_infect(10);
        sim.run_simulation(365, i);
        //cout << sim.current_epidemic_size() << endl;
    }

    return 0;
}
