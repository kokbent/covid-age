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

int main() { 

    int N = 2456274;

    vector<double> Ki; // S -> E transition rate
    double ini_Ki = 1.0522295;
    vector<TimeSeriesAnchorPoint> Ki_ap = {
        {0, 1.0     },
        {28, 0.6263 },
        {33, 0.3526 },
        {37, 0.09   },
        {68, 0.07   },
        {98, 0.07   },
        {129, 0.11  },
        {163, 0.11  },
        {194, 0.11  },
        {217, 0.13  },
        {237, 0.198 },
        {272, 0.115 },
        {311, 0.117 },
        {342, 0.1156},
        {368, 0.1223},
        {400, 0.1223}
    };
    for (size_t i = 0; i < Ki_ap.size(); i++) { Ki_ap[i].value = Ki_ap[i].value * ini_Ki; }
    Ki = stepwiseTimeSeries(Ki_ap);

    double Kasymp= 0.41/3.677037; // E -> I (asymp) rate
    double Kpres = 0.59/3.677037; // E -> I (symp) rate
    double Kmild = 0.92/3.409656;
    double Kseve = 0.08/3.409656;
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
    
    //vector<double> Krec = {1.0/9.0, 1.0/9.0, 1.0/5.0, 1.0/9.0, 1.0/2.5}; // I -> R transition rate
    vector<TimeSeriesAnchorPoint> Pcrit_ap = {
        {0, 0.4947052},
        {48, 0.4074691  },
        {62, 0.3531274  },
        {78, 0.2697153  },
        {109, 0.147316  },
        {139, 0.2186977 },
        {170, 0.187974  },
        {201, 0.11607   },
        {231, 0.1658555 },
        {262, 0.09531188},
        {311, 0.07186472},
        {400, 0.07186472}
    };
    vector<double> Pcrit = stepwiseTimeSeries(Pcrit_ap);

    vector<TimeSeriesAnchorPoint> Pdead_ap = {
        {0, 0.2033216},
        {78,  0.2832796 },
        {109, 0.181298  },
        {139, 0.09841214},
        {170, 0.06885639},
        {201, 0.1266077 },
        {231, 0.1636148 },
        {262, 0.1933856 },
        {292, 0.1559855 },
        {323, 0.02952469},
        {400, 0.2033216}
    };
    vector<double> Pdeath = stepwiseTimeSeries(Pdead_ap);

    for (size_t i = 0; i < Pcrit.size(); i++) {
        Pcrit[i] = Pcrit[i] + Pdeath[i];
        Pdeath[i] = Pdeath[i] / Pcrit[i];
    }
    
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

    //vector<double> Pdetect = {0.2/6, 0.2/6, 0.2, 0.5}; // I -> R transition rate
    //vector<double> Pdetect = {0, 0, 0, 0}; // I -> R transition rate
    double frac_infectiousness_As = 0.8;
    double frac_infectiousness_det = 0.15;
    double immunity_duration = 10000;

    for(int i=0; i<1; i++ ) {
        Event_Driven_NUCOVID sim(N, Ki, Kasymp, Kpres, Kmild, Kseve, Khosp, Kcrit,
                                 Kdeath, Krec, Pcrit, Pdeath, Pdetect, 
                                 frac_infectiousness_As, frac_infectiousness_det);
        sim.rng.seed(time(0)); // this simulator has its own RNG which must be seeded as well
        sim.Now = 9;
        sim.rand_infect(10);
        sim.run_simulation(365, i);
        //cout << sim.current_epidemic_size() << endl;
    }

    return 0;
}
