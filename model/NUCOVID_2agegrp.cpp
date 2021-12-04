#include "Time_Series.h"
#include "Event_Driven_NUCOVID.h"

vector<vector<double>> transpose2dVector( vector<vector<double>> vec2d ) {
    vector<vector<double>> tvec2d;

    for (size_t i = 1; i < vec2d.size(); i++) {
        if (not vec2d[i].size() == vec2d[1].size()) {
            cerr << "Vector of vectors not of same size." << endl;
            exit(1);
        }
    }

    for (size_t i = 0; i < vec2d[1].size(); i++) {
        vector<double> v;
        for (size_t j = 0; j < vec2d.size(); j++) { v.push_back(vec2d[j][i]); }    
        tvec2d.push_back(v);
    }

    return tvec2d;
}

vector<Node*> initialize_2nodes() {
    vector<Node*> nodes;
    int N = 1546204;
    vector<double> Ki; // S -> E transition rate

    double ini_Ki = 1.1;
    vector<TimeSeriesAnchorPoint> Ki_ap = {
        {0, 1.0     },
        {28, 0.6263 },
        {33, 0.3526 },
        {37, 0.09   },
        {68, 0.07   },
        {98, 0.07   },
        {129, 0.11  },
        {163, 0.11  },
        {194, 0.16  },
        {217, 0.19  },
        {237, 0.19  },
        {272, 0.115 },
        {311, 0.117 },
        {342, 0.1156},
        {368, 0.1223},
        {400, 0.1223}
    };
    for (size_t i = 0; i < Ki_ap.size(); i++) { Ki_ap[i].value = Ki_ap[i].value * ini_Ki; }
    Ki = stepwiseTimeSeries(Ki_ap);

    double Kasymp= 0.42/3.677037; // E -> I (asymp) rate
    double Kpres = 0.58/3.677037; // E -> I (symp) rate
    double Kmild = 0.97/3.409656;
    double Kseve = 0.03/3.409656;
    double Khosp = 1.0/4.076704;
    double Kcrit = 1.0/5.592791;
    double Kdeath = 1.0/5.459323;

    vector<double> Krec_asym(400, 1.0/9.0);
    vector<double> Krec_symm(400, 1.0/9.0);
    vector<double> Krec_crit(400, 1.0/9.671261);
    vector<double> Krec_hpc (400, 1.0/2.194643);
    vector<TimeSeriesAnchorPoint> Krh_ap = {
        {0, 1/5.78538},
        {270, 1/4.7},
        {400, 1/4.7}
    };
    vector<double> Krec_hosp = stepwiseTimeSeries(Krh_ap);
    vector<vector<double>> Krec{Krec_asym, Krec_symm, Krec_hosp, Krec_crit, Krec_hpc };
    Krec = transpose2dVector(Krec);
    
    vector<double> Pcrit(400, 0.25);
    vector<double> Pdeath(400, 0.1);

    vector<double> Pdet_asym;
    vector<double> Pdet_pres;
    vector<double> Pdet_symm;
    vector<TimeSeriesAnchorPoint> Pdsm_ap = {
        {0, 0.000630373},
        {48, 0.03520381},
        {78, 0.08413338},
        {109, 0.1474772},
        {139, 0.1528583},
        {170, 0.1064565},
        {201, 0.1514133},
        {231, 0.1608695},
        {262, 0.4160216},
        {400, 0.4160216}
    };
    Pdet_symm = stepwiseTimeSeries(Pdsm_ap);

    vector<double> Pdet_syms;
    vector<TimeSeriesAnchorPoint> Pdss_ap = {
        {0, 0.009849325},
        {31, 0.1456243 },
        {48, 0.5841068 },
        {78, 0.6877389 },
        {109, 0.9820229},
        {139, 0.5239712},
        {170, 0.5520378},
        {201, 0.7033732},
        {231, 0.881767 },
        {400, 0.881767 }
    };
    Pdet_syms = stepwiseTimeSeries(Pdss_ap);

    for (size_t i = 0; i < Pdet_symm.size(); i++) {
        Pdet_asym.push_back(Pdet_symm[i]/6);
        Pdet_pres.push_back(Pdet_symm[i]/6);
    }

    vector<vector<double>> Pdetect{Pdet_asym, Pdet_pres, Pdet_symm, Pdet_syms};
    Pdetect = transpose2dVector(Pdetect);

    double frac_infectiousness_As = 0.8;
    double frac_infectiousness_det = 0.15;

    Node* n1 = new  Node(0, N, Ki, Kasymp, Kpres, Kmild, Kseve, Khosp, Kcrit,
                         Kdeath, Krec, Pcrit, Pdeath, Pdetect,
                         frac_infectiousness_As, frac_infectiousness_det);
    nodes.push_back(n1);
    
    N = 1198476;
    Kasymp = 0.23 / 3.677037;
    Kpres = 0.77/3.677037; // E -> I (symp) rate
    Kmild = 0.804/3.409656;
    Kseve = 0.196/3.409656;

    for(size_t i = 0; i < Pdeath.size(); i++) Pdeath[i] = 0.5;
    
    ini_Ki = 0.80;
    vector<TimeSeriesAnchorPoint> Ki_ap2 = {
        {0, 1.0     },
        {28, 0.6263 },
        {33, 0.3526 },
        {37, 0.09   },
        {68, 0.07   },
        {98, 0.07   },
        {129, 0.11  },
        {163, 0.11  },
        {194, 0.11  },
        {217, 0.11  },
        {237, 0.14 },
        {272, 0.115 },
        {311, 0.117 },
        {342, 0.1156},
        {368, 0.1223},
        {400, 0.1223}
    };
    for (size_t i = 0; i < Ki_ap2.size(); i++) { Ki_ap2[i].value = Ki_ap2[i].value * ini_Ki; }
    Ki = stepwiseTimeSeries(Ki_ap2);

    Node* n2 = new  Node(1, N, Ki, Kasymp, Kpres, Kmild, Kseve, Khosp, Kcrit,
                         Kdeath, Krec, Pcrit, Pdeath, Pdetect,
                         frac_infectiousness_As, frac_infectiousness_det);
    nodes.push_back(n2);

    return(nodes);
}

void runsim (int serial) {
    cout << "Running Sim " << to_string(serial) << endl;

    vector<string> out_buffer;
    vector<Node*> nodes = initialize_2nodes();
    vector<vector<double>> infection_matrix;

    for(int i = 0; i < 2; i++) {
        vector<double> inf_prob;
        for (int j = 0; j < 2; j++) {
            double r = i == j ? 0.9 : 0.1;
            inf_prob.push_back(r);
        }
        infection_matrix.push_back(inf_prob);
    }

    Event_Driven_NUCOVID sim(nodes, infection_matrix);
    sim.rng.seed(time(0));
    sim.Now = 9;
    sim.rand_infect(3, nodes[0]);
    sim.rand_infect(7, nodes[1]);
    out_buffer = sim.run_simulation(371, false);
    string out_fname = "../out/daily_output." + to_string(serial);
    write_buffer(out_buffer, out_fname, true);

    return;
}

int main(int argc, char* argv[]) { 
    int serial = stoi(argv[1]);

    run1sim(serial);

    return 0;
}
