#ifndef ED_NUCOVID_H
#define ED_NUCOVID_H

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <queue>
#include "Utility.h"
#include <climits>

using namespace std;

typedef enum {
    PRE, ASY, SYMM, SYMS, HOS, CRI, HPC, DEA, RECA, RECM, RECH, RECC, CON, IMM
} eventType;

class Event {
    public:
        double time;
        eventType type;
//        Node* node;
        Event(const Event& o) {  time=o.time; type=o.type; }
        Event(double t, eventType e) {time=t; type=e; }
        Event& operator=(const Event& o) { time=o.time; type=o.type; return *this; }
};

class compTime {
    public:
        bool operator() (const Event* lhs, const Event* rhs) const {
            return (lhs->time>rhs->time);
        }

        bool operator() (const Event& lhs, const Event& rhs) const {
            return (lhs.time>rhs.time);
        }
};


class Event_Driven_NUCOVID {
    public:
        // S -> E -> I -> R -> S
        typedef enum {
            SUSCEPTIBLE, 
            EXPOSED, 
            ASYMPTOMATIC, 
            PRESYMPTOMATIC, 
            SYMPTOMATIC_MILD,
            SYMPTOMATIC_SEVERE,
            HOSPITALIZED,
            HOSPITALIZED_CRIT,
            CRITICAL,
            DEATH,
            RESISTANT, 
            STATE_SIZE // STATE_SIZE must be last
        } stateType;

        // constructor
        Event_Driven_NUCOVID ( int n, vector<double> ki, double ka, double kp, double km, double ks, 
                               double kh, double kc, double kd, vector<vector<double>> kr,
                               vector<double> pc, vector<double> pd, vector<vector<double>> pdets,
                               double fia, double fid) {

            N = n;
            Ki = ki;
            Kasym = ka;
            Kpres = kp; 
            Kmild = km;
            Ksevere = ks;
            Khosp = kh;
            Kcrit = kc;
            Kdeath = kd;
            Krec = kr;
            Pcrit = pc;
            Pdeath = pd;
            Pdetect = pdets;
            cumu_symptomatic = 0;

            frac_infectiousness_As = fia;
            frac_infectiousness_det = fid;
            reset();
        }

        int N;                      // population size
        vector<double> Ki;                  // param for exponential exposed duration
        double Kasym;               // param for exponential time to infectious (asymptomatic)
        double Kpres;                // param for exponential time to infectious (presymptomatic)
        double Kmild;               // param for exponential time to mild from presymp
        double Ksevere;             // param for exponential time to severe from presymp 
        double Khosp;
        double Kcrit;
        double Kdeath;
        vector<vector<double>> Krec;                // param for exponential time to recovery
        vector<double> Pcrit;
        vector<double> Pdeath;
        vector<vector<double>> Pdetect;
        double frac_infectiousness_As;
        double frac_infectiousness_det;
        size_t cumu_symptomatic;

        priority_queue<Event, vector<Event>, compTime > EventQ; // event queue
        vector<int> state_counts;   // S, E, I, R counts
        double Now;                 // Current "time" in simulation

        mt19937 rng;              // RNG

        void run_simulation(double duration, int serial) {
            double start_time = Now;
            int day = (int) Now;
            cout << setprecision(3) << fixed;
            cout << "serial\ttime\tKi\tS\tE\tAP\tSYM\tHOS\tCRIT\tDEA\tR\tcumu_sym" << endl; 
            while (next_event() and Now < start_time + duration) {
                if ((int) Now > day) {
                    cout << serial << "\t" << (int) Now << "\t"  << Ki[day] << "\t" 
                                          << state_counts[SUSCEPTIBLE] << "\t" 
                                          << state_counts[EXPOSED] << "\t" 
                                          << state_counts[ASYMPTOMATIC] + state_counts[PRESYMPTOMATIC]<< "\t" 
                                          << state_counts[SYMPTOMATIC_MILD] + state_counts[SYMPTOMATIC_SEVERE] << "\t" 
                                          << state_counts[HOSPITALIZED] + state_counts[HOSPITALIZED_CRIT]<< "\t"
                                          << state_counts[CRITICAL] << "\t"
                                          << state_counts[DEATH] << "\t"
                                          << state_counts[RESISTANT] << "\t"
                                          << cumu_symptomatic << endl;
                    day = (int) Now;
                }

                continue;
            }
        }

        int current_epidemic_size() {
            return state_counts[EXPOSED] + state_counts[ASYMPTOMATIC] + state_counts[PRESYMPTOMATIC] + 
                   state_counts[SYMPTOMATIC_MILD] + state_counts[SYMPTOMATIC_SEVERE] +
                   state_counts[HOSPITALIZED] + state_counts[CRITICAL];
        }

        void reset() {
            Now = 0.0;
           
            //vector<Node*> nodes = network->get_nodes();
            //for (unsigned int i = 0; i < network->size(); i++) nodes[i]->set_state(SUSCEPTIBLE);

            state_counts.clear();
            state_counts.resize(STATE_SIZE, 0);
            state_counts[SUSCEPTIBLE] = N;
            
            EventQ = priority_queue<Event, vector<Event>, compTime > ();
        }

        void rand_infect(int k) {   // randomly infect k people
            for (unsigned int i = 0; i < k; i++) {
                infect();
            }
            return;
        }
        
        double get_detection_modifier(double p, double r) {return p * r + (1.0 - p);}
        double get_Ki(int day) { return day < Ki.size() ? Ki[day] : Ki[Ki.size() - 1]; }
        double get_Pdet(int day, int ind) { return day < Pdetect.size() ? Pdetect[day][ind] : Pdetect[Pdetect.size() - 1][ind]; }
        double get_Pcrit(int day) { return day < Pcrit.size() ? Pcrit[day] : Pcrit[Pcrit.size() - 1]; }
        double get_Pdeath(int day) { return day < Pdeath.size() ? Pdeath[day] : Pdeath[Pdeath.size() - 1]; }
        double get_Krec(int day, int ind) { return day < Krec.size() ? Krec[day][ind] : Krec[Krec.size() - 1][ind]; }

        void infect() {
            assert(state_counts[SUSCEPTIBLE] > 0);
            state_counts[SUSCEPTIBLE]--;  // decrement susceptible groupjj
            state_counts[EXPOSED]++;      // increment exposed group

            // time to become infectious (duel between presymp and asymp)
            double Tpres = rand_exp(Kpres, &rng) + Now;
            double Tasym = rand_exp(Kasym, &rng) + Now;
            double Ti;
            double Tsym;
            double Tr;
            double Th = 0;
            double Tcr = -1;
            double Thc = 0;
            double Td = 0;
            vector<double> Times;
            vector<double> Ki_modifier;

            char symp_path;

            // Start differentiating the two paths
            if (Tpres < Tasym) {
                // SYMPTOMATIC PATH
                // time to become infectious
                Ti = Tpres;
                add_event(Tpres, PRE);
                Times.push_back(Ti);
                Ki_modifier.push_back(get_detection_modifier( get_Pdet((int) Ti, 1), frac_infectiousness_det ));

                double Tmild = rand_exp(Kmild, &rng) + Tpres;
                double Tsevere = rand_exp(Ksevere, &rng) + Tpres;

                if (Tmild < Tsevere) {
                    // Mild SYMPTOMATIC PATH
                    Tsym = Tmild;
                    add_event(Tsym, SYMM);
                    Times.push_back(Tsym);
                    Ki_modifier.push_back(get_detection_modifier( 1 - (1 - get_Pdet((int) Ti, 1)) * (1 - get_Pdet((int) Tsym, 2)), frac_infectiousness_det ));

                    Tr = rand_exp(get_Krec((int) Tsym, 1), &rng) + Tsym;
                    add_event(Tr, RECM);
                } else {
                    // Severe SYMPTOMATIC PATH
                    Tsym = Tsevere;
                    add_event(Tsym, SYMS);
                    Times.push_back(Tsym);
                    Ki_modifier.push_back(get_detection_modifier( 1 - (1 - get_Pdet((int) Ti, 1)) * (1 - get_Pdet((int) Tsym, 3)), frac_infectiousness_det ));

                    Th = rand_exp(Khosp, &rng) + Tsym;
                    add_event(Th, HOS);
                    Times.push_back(Th);
                    Ki_modifier.push_back(get_detection_modifier( 1 - (1 - get_Pdet((int) Ti, 1)) * (1 - get_Pdet((int) Th, 4)) , 0 ));

                    if (rand_uniform(0, 1, &rng) > get_Pcrit((int) Th)) {
                        // Hospitalized and recovered
                        Tr = rand_exp(get_Krec((int) Th, 2), &rng) + Th;
                        add_event(Tr, RECH);
                    } else {
                        // Hospitalized and become critical
                        Tcr = rand_exp(Kcrit, &rng) + Th;
                        add_event(Tcr, CRI);

                        if (rand_uniform(0, 1, &rng) > get_Pdeath((int) Tcr)) {
                            // Critical and recovered
                            Thc = rand_exp(get_Krec((int) Tcr, 3), &rng) + Tcr;
                            add_event(Thc, HPC);

                            Tr = rand_exp(get_Krec((int) Thc, 4), &rng) + Thc;
                            add_event(Tr, RECC);
                        } else {
                            // Critical and die
                            Tr = rand_exp(Kdeath, &rng) + Tcr;
                            add_event(Tr, DEA);
                        } 
                    }

                }

            } else {
                // ASYMPTOMATIC PATH
                // time to become infectious
                Ti = Tasym;
                add_event(Tasym, ASY);
                Times.push_back(Ti);
                Ki_modifier.push_back(get_detection_modifier( get_Pdet((int) Ti, 0), frac_infectiousness_det ) * frac_infectiousness_As);

                // time to recovery
                Tr = rand_exp(get_Krec((int) Ti, 0), &rng) + Ti;
                add_event(Tr, RECA);
            }
            
            // time to next contact
            int bin = 0;
            double Tc = rand_exp(get_Ki((int) Ti) * Ki_modifier[bin], &rng) + Ti;
            while ( Tc < Tr ) {     // does contact occur before recovery?
                add_event(Tc, CON); // potential transmission event
                while (bin < Times.size() - 1 and Times[bin+1] < Tc) {bin++;} // update bin if necessary
                Tc += rand_exp(get_Ki((int) Tc) * Ki_modifier[bin], &rng);
            }

            // time to become susceptible again (not used for now)
            //double Ts = Tr + immunity_duration; 
            //add_event(Ts, IMM);
            //return;
        }

        int next_event() {
            if ( EventQ.empty() ) return 0;
            Event event = EventQ.top(); // get the element
            EventQ.pop();               // remove from Q

            Now = event.time;           // advance time
            switch(event.type) {
                case PRE:
                    state_counts[EXPOSED]--;      // decrement exposed class
                    state_counts[PRESYMPTOMATIC]++;   // increment symptomatic class
                    break;
                case ASY:
                    state_counts[EXPOSED]--;      // decrement exposed class
                    state_counts[ASYMPTOMATIC]++;   // increment asymptomatic class
                    break;
                case SYMM:
                    state_counts[PRESYMPTOMATIC]--;      // decrement exposed class
                    state_counts[SYMPTOMATIC_MILD]++;   // increment symptomatic class
                    cumu_symptomatic++;
                    break;
                case SYMS:
                    state_counts[PRESYMPTOMATIC]--;      // decrement exposed class
                    state_counts[SYMPTOMATIC_SEVERE]++;   // increment symptomatic class
                    cumu_symptomatic++;
                    break;
                case HOS:
                    state_counts[SYMPTOMATIC_SEVERE]--;      // decrement exposed class
                    state_counts[HOSPITALIZED]++;   // increment symptomatic class
                    break;
                case CRI:
                    state_counts[HOSPITALIZED]--;      // decrement exposed class
                    state_counts[CRITICAL]++;   // increment symptomatic class
                    break;
                case HPC:
                    state_counts[CRITICAL]--;      // decrement exposed class
                    state_counts[HOSPITALIZED_CRIT]++;   // increment symptomatic class
                    break;
                case DEA:
                    state_counts[CRITICAL]--;      // decrement exposed class
                    state_counts[DEATH]++;   // increment symptomatic class
                    break;
                case RECA:
                    state_counts[ASYMPTOMATIC]--;      // decrement asymptomatic class
                    state_counts[RESISTANT]++;   // increment recovery class
                    break;
                case RECM:
                    state_counts[SYMPTOMATIC_MILD]--;      // decrement symptomatic class
                    state_counts[RESISTANT]++;   // increment recovery class
                    break;
                case RECH:
                    state_counts[HOSPITALIZED]--;      // decrement symptomatic class
                    state_counts[RESISTANT]++;   // increment recovery class
                    break;
                case RECC:
                    state_counts[HOSPITALIZED_CRIT]--;      // decrement symptomatic class
                    state_counts[RESISTANT]++;   // increment recovery class
                    break;
                case IMM:
                    state_counts[RESISTANT]--;      // decrement recovery class
                    state_counts[SUSCEPTIBLE]++;   // increment susceptible class
                    break;
                case CON:
                    {
                        const int rand_contact = rand_uniform_int(0, N, &rng);
                        if (rand_contact < state_counts[SUSCEPTIBLE]) infect();
                    }
                    break;
                default:    
                    cerr << "Unknown event type encountered in simulator: " << event.type << "\nQuitting.\n";
            }
            return 1;
        }

        void add_event( double time, eventType type) {
            EventQ.push( Event(time,type) );
            return;
        }

};
#endif
