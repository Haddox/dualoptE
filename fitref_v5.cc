#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <map>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <sstream>
#include <unordered_map>

using namespace std;

/////////////////////////
//  AAs
int
aa2idx (char AA) {
	static std::map<char,int> aamap;
	if (aamap.empty()) {
		aamap['A']=0;  aamap['C']=1;  aamap['D']=2;  aamap['E']=3;  aamap['F']=4;
		aamap['G']=5;  aamap['H']=6;  aamap['I']=7;  aamap['K']=8;  aamap['L']=9;
		aamap['M']=10; aamap['N']=11; aamap['P']=12; aamap['Q']=13; aamap['R']=14;
		aamap['S']=15; aamap['T']=16; aamap['V']=17; aamap['W']=18; aamap['Y']=19;
	}
	return aamap[AA];
}

char
idx2aa (int ii) {
	static std::map<int,char> aamap;
	if (aamap.empty()) {
		aamap[0]='A';  aamap[1]='C';  aamap[2]='D';  aamap[3]='E';  aamap[4]='F';
		aamap[5]='G';  aamap[6]='H';  aamap[7]='I';  aamap[8]='K';  aamap[9]='L';
		aamap[10]='M'; aamap[11]='N'; aamap[12]='P'; aamap[13]='Q'; aamap[14]='R';
		aamap[15]='S'; aamap[16]='T'; aamap[17]='V'; aamap[18]='W'; aamap[19]='Y';
	}
	return aamap[ii];
}

// based on layer design:
// 		<Action selector_logic="surface AND loop" aas="DEGHKNPQRST"/>
bool
ispolar (int ii) {
	static std::map<int,bool> ispolar;
	if (ispolar.empty()) {
		ispolar[0]=0;  ispolar[1]=0;  ispolar[2]=1;  ispolar[3]=1;  ispolar[4]=0;
		ispolar[5]=1;  ispolar[6]=1;  ispolar[7]=0;  ispolar[8]=1;  ispolar[9]=0;
		ispolar[10]=0; ispolar[11]=1; ispolar[12]=1; ispolar[13]=1; ispolar[14]=1;
		ispolar[15]=1; ispolar[16]=1; ispolar[17]=0; ispolar[18]=0; ispolar[19]=0;
	}
	return ispolar[ii];
}


char
three2one (std::string three) {
	static std::map<std::string,char> aamap;
	if (aamap.empty()) {
		aamap["ALA"]='A';  aamap["CYS"]='C';  aamap["ASP"]='D';  aamap["GLU"]='E';  aamap["PHE"]='F';
		aamap["GLY"]='G';  aamap["HIS"]='H';  aamap["ILE"]='I';  aamap["LYS"]='K';  aamap["LEU"]='L';
		aamap["MET"]='M'; aamap["ASN"]='N'; aamap["PRO"]='P'; aamap["GLN"]='Q'; aamap["ARG"]='R';
		aamap["SER"]='S'; aamap["THR"]='T'; aamap["VAL"]='V'; aamap["TRP"]='W'; aamap["TYR"]='Y';
	}
	return aamap[three];
}


// aa frequencies
double
get_freq( int ii ) {
	static std::map<int,double> aamap;
	if (aamap.empty()) {
		aamap[0]=7.4;  aamap[1]=3.3;  aamap[2]=5.9;  aamap[3]=5.8;  aamap[4]=4.0;
		aamap[5]=7.4;  aamap[6]=2.9;  aamap[7]=3.8;  aamap[8]=7.2;  aamap[9]=7.6;
		aamap[10]=1.8; aamap[11]=4.4; aamap[12]=5.0; aamap[13]=3.7; aamap[14]=4.2;
		aamap[15]=8.1; aamap[16]=6.2; aamap[17]=6.8; aamap[18]=1.3; aamap[19]=3.3;
	}
	return aamap[ii];
}

// aa frequencies
double
get_sum( std::vector< double > const &ref ) {
    double sum = 0.0;
    for(std::vector<double>::const_iterator it = ref.begin(); it != ref.end(); ++it)
        sum += *it;
	return sum;
}

enum TOPOS {
    TOPO_3H = 0,
    TOPO_4H,
    TOPO_EHEE,
    TOPO_FERR
};


/////////////////////////
//  correlations
void 
rank_vector(
		std::vector<double> const &vec,
		std::vector<double> &rank
) {
	int N = vec.size();
	std::vector<std::pair<double,int>> vsort(N);
	rank.resize(N);

	for (int i=0; i<N; ++i) {
		vsort[i] = std::make_pair( vec[i],i );
	}

    std::sort(vsort.begin(), vsort.end(), 
		[] (std::pair<double, int> const &a, std::pair<double, int> const &b) { return a.first < b.first; });

	// fix ties
	int currrank = 1;
    for (int i = 1; i <= N; i++) {
		if ( i==N || vsort[i].first != vsort[i-1].first ) {
			double tierank = (currrank + i-1.0)/2.0;
			for (int j=currrank; j<i; ++j) rank[ vsort[j].second ] = tierank;
			currrank=i;
		}
    }
}

double 
pearson(std::vector<double> const &E, std::vector<double> const &S) {
	// correl
	double Esum=0, Esum2=0, Ssum=0, Ssum2=0, ESsum=0;

	// slow (but numerically stable) version
	int N = E.size();
	for (int i=0; i<N; ++i) {
		Esum += E[i];
		Ssum += S[i];
	}
	Esum /= N;
	Ssum /= N;

	for (int i=0; i<N; ++i) {
		Esum2 += (E[i]-Esum)*(E[i]-Esum);
		Ssum2 += (S[i]-Ssum)*(S[i]-Ssum);
	}
	Esum2 = std::sqrt(Esum2 / N);
	Ssum2 = std::sqrt(Ssum2 / N);

	for (int i=0; i<N; ++i) {
		ESsum += ((S[i]-Ssum)/Ssum2)*((E[i]-Esum)/Esum2);
	}
	ESsum /= N;

	return ESsum;
}

double 
spearman(std::vector<double> v1, std::vector<double> v2) {
	int n = v1.size();
	std::vector<double> R1;
	std::vector<double> R2;
    std::vector<double> d;

	rank_vector(v1,R1);
	rank_vector(v2,R2);

	return pearson(R1,R2);
}



/////////////////////////
//  NelderMead: simplex optimization
class NelderMead {
private:
	double ALPHA, BETA, GAMMA;

	std::vector< double > wts_best;
	double score_best;

public:
	NelderMead() {
		ALPHA=1;
		BETA=0.5;
		GAMMA=2.0;

		score_best = 1e6;
		wts_best.resize(20,0.0);
	}

	void
	checkpoint(
		std::vector< double > const &wts_i,
		double score_i
	) {
		if (score_i<score_best) {
			score_best = score_i;
			wts_best = wts_i;
		}
	}

	void
	recover_best(
		std::vector< double > &wts_i,
		double &score_i
	) {
		score_i = score_best;
		wts_i = wts_best;
	}

	// Helper function - find the lowest, 2nd highest and highest position
	void
	FindMarkers (
		std::vector<double> const &Y,
		int &idx_lo,
		int &idx_nhi,
		int &idx_hi
	) {
		int N = Y.size()-1;

		idx_lo=0;
		if (Y[0]>Y[1]) {
			idx_hi=0; idx_nhi=1;
		} else {
			idx_hi=1; idx_nhi=0;
		}

		for(int i=0; i<=N; ++i) {
			if (Y[i] < Y[idx_lo]) idx_lo= i;
			if (Y[i] > Y[idx_hi]) {
				idx_nhi = idx_hi; idx_hi = i;
			} else if (Y[i] > Y[idx_nhi] && idx_hi != i) {
				idx_nhi = i;
			}
		}
	}

	void CalcCentroid (
		std::vector< std::vector<double> > const&P,
		int idx_hi,
		std::vector<double> &C
	) {
		int N = P.size()-1;

		C.resize(N);
		for (int j=0; j<N; ++j) {
			C[j]=0.0;
			for(int i=0; i<=N; ++i) {
				if (i!=idx_hi) C[j] += P[i][j];
			}
			C[j] /= N;
		}
	}

	void
	AdjustCentroid (
		std::vector<double> &C,
		std::vector< std::vector<double> > const&P,
		int idx_hi_o,
		int idx_hi
	) {
		int N = P.size()-1;

		if (idx_hi_o != idx_hi) {
			for (int j=0; j<N; ++j) {
				C[j] += (P[idx_hi_o][j] - P[idx_hi][j]) / N;
			}
		}
	}

	void
	CalcReflection(
		std::vector<double> const &p1,
		std::vector<double> const &p2,
		double scale,
		std::vector<double> &p_refl
	) {
		int N = p1.size();

		p_refl.resize(N);
		for (int j=0; j<N; ++j) {
			p_refl[j] = p1[j] + scale*(p1[j]-p2[j]);
		}
	}

	void
	run(
		std::vector<double> const& guesses,
		std::vector<double> const& scales,
		double (*evaluatorFn)(std::vector<double> const&),
		double tol = 5e-4,
		int itmax = 2000,
		bool verbose = false
	) {
		int N = guesses.size();
		std::vector< std::vector<double> > vertices(N+1);
		std::vector<double> scores(N+1,0);

		// construct vertices
		vertices[0] = guesses;
		for (int i=0; i<N; ++i) {
			vertices[i+1] = guesses;
			vertices[i+1][i] += scales[i];
		}

		// evaluate vertices
		for (int i=0; i<=N; ++i) {
			scores[i] = (*evaluatorFn)(vertices[i]);
			checkpoint( vertices[i], scores[i] );
		}

		// run
		bool recalc=true;
		std::vector<double> C(N,0); // centroid
		std::vector<double> PR(N,0); // reflection point
		double YR = 0; // reflection point value
		int ihi_o = 0;

		int iter = 1;

		while (iter++<itmax) {
			int ilo, inhi, ihi;
			FindMarkers(scores, ilo, inhi, ihi);

			// Stopping conditions
			double rtol = 2*abs(scores[ihi]-scores[ilo])/(abs(scores[ihi])+abs(scores[ilo])+1e-12);
			if (rtol<tol) break;

			if (recalc) {
				CalcCentroid(vertices, ihi, C);
			} else {
				AdjustCentroid(C, vertices, ihi_o, ihi);
			}

			recalc = false;

			// Determine the reflection point, evaluate its value
			CalcReflection(C, vertices[ihi], ALPHA, PR);
			YR = (*evaluatorFn)(PR);
			checkpoint( PR, YR );

			if (YR < scores[ilo]) {
				std::vector<double> PE(N,0);
				CalcReflection(C, PR, -GAMMA, PE);
				double YE = (*evaluatorFn)(PE);
				checkpoint( PE, YE );
				if (YE < scores[ilo]) {
					vertices[ihi] = PE; scores[ihi] = YE;
				} else {
					vertices[ihi] = PR; scores[ihi] = YR;
				}
			} else if (YR >= scores[inhi]) {
				if (YR < scores[ihi] ) {
					vertices[ihi] = PR; scores[ihi] = YR;
				}

				std::vector<double> PC(N,0);
				CalcReflection(C, vertices[ihi], -BETA, PC);
				double YC=(*evaluatorFn)(PC);
				checkpoint( PC, YC );

				if (YC < scores[ihi]) {
					vertices[ihi] = PC; scores[ihi] = YC;
				} else {
					for(int i=0; i<=N; ++i) {
						if (i!=ilo) {
							std::vector<double> vi_old = vertices[i];
							CalcReflection(vertices[ilo], vi_old, -BETA, vertices[i]);
							scores[i] = (*evaluatorFn)(vertices[i]);
							checkpoint( vertices[i], scores[i] );
						}
					}
					recalc = true;
				}
			} else {
				vertices[ihi] = PR; scores[ihi] = YR;
			}

			ihi_o = ihi;
		}
	}
};

/////////////////////////
//  PointstabilityData: one data point for 1-at-a-time recov
class PointstabilityData {
public:
	PointstabilityData() {
		native_ = 0;
		data_.resize(20,0.0);
	}

	double &
	operator[](int i) {
		return data_[i];
	}

	double const &
	operator[](int i) const {
		return data_[i];
	}

	void native( int i) { native_=i; }
	int native( ) { return native_; }

	void sasa( double i ) { sasa_=i; }
	double sasa( ) { return sasa_; }

	void normalize() {
		double sum=0.0;
		for (int i=0; i<20; ++i) sum += data_[i];
		for (int i=0; i<20; ++i) data_[i] /= sum;
	}

private:
	int native_;
	std::vector<double> data_;
    double sasa_;
};

/////////////////////////
//  StabilityData: one data point for design stability
class StabilityData {
public:
	StabilityData() {
		stability_ = base_energy_ = 0.0;
	}

	StabilityData(std::string tag_in, std::string seq_in, double stab_in, double ener_in, int hyd_in) {
		tag_ = tag_in;
		seq_ = seq_in;
		stability_ = stab_in;
		base_energy_ = ener_in;
        hydrophobicity_ = hyd_in;
	}

	int
	length() const {
		return seq_.length();
	}

	char const &
	operator[](int i) const {
		return seq_[i];
	}

	void stability( double x ) { stability_=x; }
	double stability( ) { return stability_; }

	void tag( std::string tag_in ) { tag_ = tag_in; }
	std::string tag( ) { return tag_; }

	void energy( double x ) { base_energy_=x; }
	double energy( ) { return base_energy_; }

	void hydrophobicity( int x ) { hydrophobicity_=x; }
	double hydrophobicity( ) { return hydrophobicity_; }

private:
	std::string seq_,tag_;
	double stability_, base_energy_;
    int hydrophobicity_;
};

/////////////////////////
//  DDGData: one data point for design stability
class DDGData {
public:
	DDGData() {
		Ecalc_ = Eexp_ = 0.0;
		aafrom_ = aato_ = 'X';
	}

	DDGData(char aafrom, char aato, double Ecalc, double Eexp) {
		aafrom_ = aafrom; aato_ = aato;
		Ecalc_ = Ecalc; Eexp_ = Eexp;
	}

	char aafrom() { return aafrom_; }
	char aato() { return aato_; }

	double ddg_calc() { return Ecalc_; }
	double ddg_exp() { return Eexp_; }
	void ddg_calc(double x) { Ecalc_=x; }
	void ddg_exp(double x) { Eexp_=x; }

private:
	char aafrom_, aato_;
	double Ecalc_, Eexp_;
};


////////////////////////////////////
//  global data
std::vector <PointstabilityData> ENERGIES, COUNTS;
std::vector <std::vector <StabilityData>> STABILITIES;
std::vector <DDGData> DDGS;
double W_RECOV=1.0, W_FIXBB=1.0, W_STABILITY=1.0, W_DDG=1.0, W_KL1=0.5, W_KL2=0.5;
int MAXITER=2000;
int NSTRUCT=0;
double KT=0.5;
double SASA_CUT=1e20;
////////////////////////////////////

//
void
evaluate_natives( 
	std::vector<double> const &values,
	double &score_recov,
	double &seq_recov,
	double &score_diverge,
    int verbose
) {
	std::vector<double> native_distr_bur(20,0.0);
	std::vector<double> decoy_distr_bur(20,0.0);
	std::vector<double> native_distr_surf(20,0.0);
	std::vector<double> decoy_distr_surf(20,0.0);
    int Nsurf=0, Nbur=0;

	double kT = KT;

	int N = ENERGIES.size();

    score_recov = seq_recov = score_diverge = 0.0;

	for (int i=0; i<N; ++i) {
		double e_min=1e6, prof_wt=0.0, prof_max=0.0, probsum=0.0;
		int aa_min = -1;
        bool is_surface = (ENERGIES[i].sasa() > SASA_CUT);

		for (int j=0; j<20; ++j) {
            bool is_surface_nonpolar = is_surface && !ispolar(j);
            if (is_surface_nonpolar) continue;

			double e_j = ENERGIES[i][j]	+ values[j];
			if (e_min > e_j) {
				e_min = e_j;
				aa_min = j;
			}
		}

		for (int j=0; j<20; ++j) {
            bool is_surface_nonpolar = is_surface && !ispolar(j);

			double e_j = ENERGIES[i][j]	+ values[j];
			double p_j = std::exp( - (e_j - e_min)/kT );
            if (is_surface_nonpolar) p_j = 0.0;
			probsum += p_j;
			prof_wt += COUNTS[i][j] * p_j;
			prof_max += COUNTS[i][j] * COUNTS[i][j];
		}
  		score_recov -= prof_wt/probsum;

        if (!is_surface) {
    		native_distr_bur[ COUNTS[i].native() ]++;
    		decoy_distr_bur[ aa_min ]++;
            Nbur++;
        } else {
       		native_distr_surf[ COUNTS[i].native() ]++;
       		decoy_distr_surf[ aa_min ]++;
            Nsurf++;
        }

        if (aa_min == COUNTS[i].native()) {
            seq_recov += 1.0;
        }
	}

	for (int j=0; j<20; ++j) {
		double p_i_bur = (native_distr_bur[j] + native_distr_surf[j]) / N;
		double q_i_bur = (decoy_distr_bur[j] + decoy_distr_surf[j]) / N;
		if (p_i_bur == 0) { p_i_bur = 0.5/N; }
		if (q_i_bur == 0) { q_i_bur = 0.5/N; }
		score_diverge += p_i_bur * log( p_i_bur/q_i_bur );
	}
	score_recov /= N;
	seq_recov /= N;
}

//
void 
evaluate_stability( 
	std::vector<double> const &values,
	std::vector<double> &splits,
	int verbose
) {
	int Nset = STABILITIES.size();
    splits.clear();
    splits.resize(Nset, 0.0);

    //std::vector<int> hyd_bins = {794, 2529, 2664, 2769, 2860, 2946, 3033, 3134, 3255, 3454, 4430};
    std::vector<int> hyd_bins = {794, 4430};

    for (int i=0; i<Nset; ++i) {
        int N = STABILITIES[i].size();
        if (verbose >= 2) { std::cout << "stab set " << i << ": "; }
        int nvalidbins = 0;
        for (int bin=0; bin<hyd_bins.size()-1; ++bin) {
            int nstable=0;
            int hydmin = hyd_bins[bin];
            int hydmax = hyd_bins[bin+1];

            std::vector< double > Es;
            std::vector< double > Ss;
            double Ej,Sj;

            for (int j=0; j<N; ++j) {
                if ( STABILITIES[i][j].hydrophobicity() >= hydmin &&  STABILITIES[i][j].hydrophobicity() < hydmax) {
                    Ej = STABILITIES[i][j].energy();
                    for (int k=0; k<STABILITIES[i][j].length(); ++k) Ej  += values[aa2idx(STABILITIES[i][j][k])];
                    Ej /= STABILITIES[i][j].length(); // energy per residue
                    Sj = STABILITIES[i][j].stability();
                    if (Sj > 1.0) nstable++;
                    Es.push_back( Ej );
                    Ss.push_back( Sj );
                }
            }
            double split_i = Ss.size() * pearson( Es, Ss );

            if (nstable < 30) {
                split_i = 0.0;       
            } else {
                nvalidbins += Ss.size();
            }
            if (verbose >= 2) { std::cout << split_i << " (" << nstable << "/" << Ss.size() << ")"; }
            splits[i] += split_i;
        }
        splits[i] /= nvalidbins;
        if (verbose >= 2) { std::cout << " = " << splits[i] << std::endl; }
    }
}

void 
evaluate_ddgs( 
	std::vector<double> const &values,
	double &pearson0,
	double &spearman0,
	int verbose
) {
	int N = DDGS.size();

	std::vector< double > Eexp(N,0);
	std::vector< double > Ecalc(N,0);
	for (int i=0; i<N; ++i) {
		Eexp[i] = DDGS[i].ddg_exp();
		Ecalc[i] = DDGS[i].ddg_calc() + values[aa2idx( DDGS[i].aato() )] - values[aa2idx( DDGS[i].aafrom() )];

        if (verbose>=3) std::cout << "ddg " << Eexp[i] << " " << Ecalc[i] << std::endl;
	}

	pearson0 = -pearson( Eexp, Ecalc );
	spearman0 = -spearman( Eexp, Ecalc );
    if (verbose>=3) std::cout << "ddg correl " << -pearson0 << " " << -spearman0 << std::endl;
}

void 
evaluate_fixbb( 
	std::vector<double> const &values,
	int iter,
	double &score_fixbbrecov,
	double &score_diverge,
	int verbose
) {
	//1 write weights file beta_nov16_iter.wts
	{
		std::stringstream ss1;
		ss1 << "opte_iter_" << iter << ".wts";
		std::string weightfile = ss1.str();
		std::ofstream ofwts( weightfile.c_str(), ofstream::out);

		ofwts << "METHOD_WEIGHTS ref ";
		for (int i=0; i<20; ++i) { ofwts << values[i] << " "; }
		ofwts << "\n";
	}

	//2 call run.sh
	{
		std::stringstream ss2;
		ss2 << "./runALL_fastbb.sh " << iter;
		std::string cmd = ss2.str();
		system(cmd.c_str());
	}

	//3 read scores_iter
	{
		double score=0.0;
		std::stringstream ss3;
		ss3 << "scores_" << iter;
		std::string scorefile = ss3.str();
		std::ifstream ifener (scorefile.c_str(), ifstream::in);
		ifener >> score >> score_fixbbrecov >> score_diverge;
		score_fixbbrecov *= -1;
	}
}


// point evaluator
double
evaluator( std::vector<double> const &values, int verbose ) {
	static int ITER=1;
	double score_recov=0.0, seq_recov=0.0, score_fixbbrecov=0.0, score_diverge1=0.0, score_diverge2=0.0;
	std::vector<double> score_stability(4,0.0);
	double score_ddg1=0.0, score_ddg2=0.0;

    // build reference weights!
    //  biased set accounts for distribution diffs
    //  weight is fitted
	std::vector<double> refset(20), refset_bias(20);

	double sum = get_sum(values) - values[0];
    double wt_freq = values[0];
	for (int i=0; i<20; ++i) {
		refset[i] = values[i+1] - sum/20;
		refset_bias[i] = refset[i] - sum/20 + wt_freq*log(get_freq(i));
    }

	if (W_RECOV > 0) 
        evaluate_natives(refset_bias, score_recov, seq_recov, score_diverge1, verbose);
	if (W_FIXBB > 0) 
        evaluate_fixbb(refset_bias, ITER, score_fixbbrecov, score_diverge2, verbose );
	if (W_STABILITY > 0) 
        evaluate_stability(refset, score_stability, verbose );
	if (W_DDG > 0) 
        evaluate_ddgs(refset, score_ddg1, score_ddg2, verbose );

    double score_stabtot = 0.25*(
        score_stability[TOPO_3H] 
        + score_stability[TOPO_4H] 
        + score_stability[TOPO_EHEE] 
        + score_stability[TOPO_FERR] 
    );

    double score_seqrecovtot =
        W_RECOV*score_recov + 
        W_KL1*score_diverge1;

	double score = 
        W_RECOV*score_recov + 
        W_STABILITY*score_stabtot +
        W_FIXBB*score_fixbbrecov + 
        W_DDG*score_ddg1 + 
        W_KL1*score_diverge1 + 
        W_KL2*score_diverge2;

	if (verbose>=1) {
		std::cout << NSTRUCT << " " << score << " " << score_seqrecovtot << " " << score_stabtot << std::endl;
		std::cout << "  fixbb: " << score_fixbbrecov  << " div:" << score_diverge2 << std::endl;
        std::cout << "  ssm:" << score_recov << " (" << 100*seq_recov << "%recov) " << " div:" << score_diverge1 << std::endl;
        std::cout << "  stab: " << score_stability[0] << " " << score_stability[1] << " " << score_stability[2] << " " << score_stability[3] << std::endl;
	}

	ITER++;
	return (score);
}

// point evaluator
double
evaluator( std::vector<double> const &values ) {
	return evaluator( values, 0 );
}

void 
parse_natives(std::vector<std::string> tags) {
	for (int arg=0;arg<tags.size(); ++arg) {
		std::string tag(tags[arg]);

		// temp storage for 1 struct worth of data
		std::vector <PointstabilityData> energies_i, counts_i;
		std::vector <int> resids;

		std::string scorefile = tag+".ENERGIES";
		ifstream ifener (scorefile.c_str(), ifstream::in);
		std::string buf, line;

		while ( getline (ifener,line) ) {
			std::stringstream ss(line);
			std::vector<string> tokens; // Create vector to hold our words
			while (ss >> buf) tokens.push_back(buf);

			if (tokens.size() < 23) continue;

			PointstabilityData energies_ij;
			energies_ij.native( atoi(tokens[1].c_str()) );
			for (int i=0; i<20; ++i) {
				energies_ij[i] = atof(tokens[i+2].c_str());
			}
            energies_ij.sasa( atof(tokens[22].c_str()) );
			energies_i.push_back( energies_ij );
		}

		std::string countfile = tag+".COUNTS";
		ifstream ifcount (countfile.c_str(), ifstream::in);

		while ( getline (ifcount,line) ) {
			stringstream ss(line);
			std::vector<string> tokens; // Create vector to hold our words
			while (ss >> buf) tokens.push_back(buf);

			if (tokens.size() < 22) continue;

			PointstabilityData counts_ij;
			counts_ij.native( aa2idx(tokens[1][0]) );
			for (int i=0; i<20; ++i) {
				counts_ij[i] = atof(tokens[i+2].c_str());
			}
			counts_ij.normalize();
			counts_i.push_back( counts_ij );
		}

		for (int i=0; i<energies_i.size(); ++i) {
			if (counts_i.size() > energies_i[i].native()-1) {
				ENERGIES.push_back( energies_i[i] );
				COUNTS.push_back( counts_i[energies_i[i].native()-1] );
                NSTRUCT++;
			}
		}
	}
}

std::vector<std::string> 
splitstring (const std::string &s, char delim) {
    std::vector<std::string> result;
    std::stringstream ss (s);
    std::string item;
    while (getline (ss, item, delim)) {
        result.push_back (item);
    }
    return result;
}

void 
parse_stability(std::vector<std::string> tags) {
	std::map <std::string, int> tag2set;
	std::map <std::string, int> tag2idx;

	if (tags.size() < 2) exit(1);

    std::vector <std::vector <StabilityData>> allSTABILITIES(4);

    // sequence, topo, hydropathy, %hydro, buried_npsa, tag, stability
    ifstream ifref (tags[0].c_str(), ifstream::in);
    std::string buf, line;

    int counter = 0;
    while ( getline (ifref,line) ) {
        stringstream ss(line);
        std::vector<string> tokens; // Create vector to hold our words
        while (ss >> buf) tokens.push_back(buf);
        if (tokens.size() != 7) continue;

        int bin=0;
        if (tokens[1] == "3h") bin = TOPO_3H;
        if (tokens[1] == "4h") bin = TOPO_4H;
        if (tokens[1] == "EHEE") bin = TOPO_EHEE;
        if (tokens[1] == "ferr") bin = TOPO_FERR;

        double stabscore = std::atof(tokens[6].c_str());
        // squeeze input stabilities to [0,2.25]
        if (stabscore < 0) stabscore = 0;
        if (stabscore > 2.25) stabscore = 2.25;
        int hydro_score = std::atoi(tokens[4].c_str());
        std::string tag = tokens[5];
        allSTABILITIES[bin].push_back( StabilityData( tokens[3], tokens[0], stabscore, 0.0, hydro_score ));
        if (tag2set.find( tag ) != tag2set.end()) { 
            std::cerr << "Duplicate token " << tag << std::endl;
            exit(1);
        }
        tag2set [ tag ] = bin;
        tag2idx [ tag ] = allSTABILITIES[bin].size()-1;
    }

	// files 2-N are the scorefiles
    STABILITIES.resize(4);
	for (int arg=1; arg<tags.size(); ++arg) {
		ifstream ifref (tags[arg].c_str(), ifstream::in);
		std::string buf, line;
		getline (ifref,line); // skip 1st line
		int ref_field = -1;
		while ( getline (ifref,line) ) {
			stringstream ss(line);
			std::vector<string> tokens; // Create vector to hold our words
			while (ss >> buf) tokens.push_back(buf);
			if (tokens.size() < 3) continue;
			if (tokens[ tokens.size()-1] == "description") {
				for (int j=0; j<tokens.size(); ++j) {
					if (tokens[j] == "ref") ref_field = j;
				}
				continue;
			}

            std::string fulltag = tokens[ tokens.size()-1];
            std::vector<std::string> bits = splitstring( fulltag,'_' );
            std::string tag = bits[0]+"_"+bits[1]+"_"+bits[2]+"_"+bits[3];

			if (tag2idx.find( tag ) == tag2idx.end()) {
				std::cerr << "ERROR! parse_stability: " << tag << " not found" << std::endl;
				exit(1);
			}
			int set = tag2set[ tag ];
			int idx = tag2idx[ tag ];
			if (ref_field>=0) {
				allSTABILITIES[set][idx].energy( std::atof(tokens[1].c_str()) - std::atof(tokens[ref_field].c_str()) );
			} else {
				allSTABILITIES[set][idx].energy( std::atof(tokens[1].c_str()) );
			}
            STABILITIES[set].push_back( allSTABILITIES[set][idx] );
            NSTRUCT++;
		}
	}

    int stabcount=0, ecount=0;
    for (int i=0; i<4; ++i) {
        stabcount += allSTABILITIES[i].size();
        ecount += STABILITIES[i].size();
    }
    //if (stabcount != ecount) {
    //    std::cout << "WARNING: Scores not found for all stabilities (expected " << stabcount << " found " << ecount << ")" << std::endl;
    //}

}

void 
parse_ddgs(std::vector<std::string> tags) {
	std::map <std::string, int> tag2idx;

	if (tags.size() < 2) exit(1);
    std::vector <DDGData> allDDGS;

	// 1st file is the reference
	{
		//id,aafrom,pos,aato,ddg
		ifstream ifref (tags[0].c_str(), ifstream::in);
		std::string buf, line;
		while ( getline (ifref,line) ) {
			stringstream ss(line);
			std::vector<string> tokens; // Create vector to hold our words
			while (ss >> buf) tokens.push_back(buf);
			if (tokens.size() != 5) continue;

			std::string tag = tokens[0]+"_"+tokens[1]+tokens[2]+tokens[3];
			allDDGS.push_back( DDGData( tokens[1][0], tokens[3][0], 0.0, std::atof(tokens[4].c_str()) ));
			tag2idx [ tag ] = allDDGS.size()-1;
		}
	}

	// files 2-N are the scorefiles
	for (int arg=1;arg<tags.size(); ++arg) {
		//name,seq,stabilityscore
		ifstream ifref (tags[arg].c_str(), ifstream::in);

		int start = tags[arg].rfind("/")+1;
		int stop = tags[arg].find(".ddg", start+1);
		std::string recordtag = tags[arg].substr( start, stop-start );

		std::string buf, line;
		int ref_field = -1;
		double score_wt=0.0,score_mut=0.0;
		int count_wt=0, count_mut=0;
		while ( getline (ifref,line) ) {
			stringstream ss(line);
			std::vector<string> tokens; // Create vector to hold our words
			while (ss >> buf) tokens.push_back(buf);

			if (tokens[2]=="WT:") {
				score_wt += std::atof(tokens[3].c_str());
				for (int i=4; i<tokens.size()-1;i+=2) { if (tokens[i] == "ref:") score_wt -= std::atof(tokens[i+1].c_str()); }
				count_wt++;
			} else {
				score_mut += std::atof(tokens[3].c_str());
				for (int i=4; i<tokens.size()-1;i+=2) { if (tokens[i] == "ref:") score_mut -= std::atof(tokens[i+1].c_str()); }
				count_mut++;
			}
		}
		if (tag2idx.find( recordtag ) == tag2idx.end()) {
			std::cerr << "ERROR! ddg tag " << recordtag << " not found." << std::endl;
			exit(1);
		}
		int idx = tag2idx[ recordtag ];
        if (count_mut != count_wt) {
            std::cerr << "error in ddg file " << tags[arg] << std::endl;
        }
		allDDGS[idx].ddg_calc( score_mut/count_mut - score_wt/count_wt );
        DDGS.push_back( allDDGS[idx] );
        NSTRUCT++;
	}
    //if (allDDGS.size() != DDGS.size()) {
    //    std::cout << "WARNING: Scores not found for all DDGs (expected " << allDDGS.size() << " found " << DDGS.size() << ")" << std::endl;
    //}
}

//
int main(int ARGV, char **ARGC) {
	if (ARGV < 4) {
		std::cerr << "usage: " << ARGC[0] << " <native_tags> <stability_tags> <ddg_tags>"
            << 	" wt_1atatime_recov wt_stability wt_ddg wt_fixbb_recov wt_kldiv MAXITER\n";

		return 1;
	}

    // beta16 refs, no correction
	double start[] = { 
        0.0, 
        // beta16
        1.83940,3.61960,-2.37160,-2.73480,1.0402,0.83697,-0.45461,0.73287,-1.51070,0.18072,0.60916,-0.93687,-2.41190,-0.18838,-1.2888,-0.77834,-1.08740,1.93420,1.69060,1.2797
        // beta16.nostab
        //2.3386,3.2718,-2.2837,-2.5358,1.4028,1.2108,0.134426,1.0317,-1.6738,0.729516,1.2334,-0.873554,-5.1227,-1.0644,-1.281,-1.1772,-1.425,2.085,3.035,0.964136
        // beta16.nostab-refit
        //2.21484,3.4346,-2.39781,-2.68179,1.31859,1.23958,0.0616106,0.913459,-1.82268,0.606044,1.16192,-0.983199,-5.34561,-1.19243,-1.41861,-1.31616,-1.63995,1.97455,3.1978,0.863196
        // new
        //2.7777,3.5846,-2.088,-2.6565,1.291,0.45199,0.30925,0.98694,-1.6719,0.40764,0.67518,-0.63745,-4.1627,-1.0534,-1.1989,-1.0738,-1.5912,2.0715,2.6058,0.9722
    };
	double scale[] = { 
        0.0,  // set to nonzero to allow correction
        5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 
        5, 5, 5, 5, 5, 
        5, 5, 5, 5, 5 };

    // uncomment to randomize start
	//srand (time(NULL));
	//for (int i=0; i<20; ++i) {
	//	start[i] = -2 + 4*(static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
	//}

	std::vector<double> startV(start, start + sizeof(start) / sizeof(start[0]) );
    std::vector<double> workingV = startV;
	std::vector<double> scaleV(scale, scale + sizeof(scale) / sizeof(scale[0]) );

	std::vector<std::string> nativetags, stabilitytags, ddgtags;
	std::string line;

    // read wts
	if (ARGV>4) { W_RECOV=std::atof( ARGC[4] ); }
	if (ARGV>5) { W_STABILITY=std::atof( ARGC[5] ); }
	if (ARGV>6) { W_DDG=std::atof( ARGC[6] ); }
	if (ARGV>7) { W_FIXBB=std::atof( ARGC[7] ); }
	if (ARGV>8) { W_KL1=W_KL2=std::atof( ARGC[8] ); }
	if (ARGV>9) { MAXITER=std::atof( ARGC[9] ); }
	if (ARGV>10) { SASA_CUT=std::atof( ARGC[10] ); }

	// 1 parse native recovery data
	if (W_RECOV != 0) {
		ifstream ifs_nat (ARGC[1], ifstream::in);
		while ( getline (ifs_nat,line) ) nativetags.push_back( line );
		parse_natives ( nativetags );
	}

	// 2 parse stability data
	if (W_STABILITY != 0) {
		ifstream ifs_stability (ARGC[2], ifstream::in);
		while ( getline (ifs_stability,line) ) stabilitytags.push_back( line );
		parse_stability ( stabilitytags );
	}

	// 3 parse ddg data
	if (W_DDG != 0) {
		ifstream ifs_ddgs (ARGC[3], ifstream::in);
		while ( getline (ifs_ddgs,line) ) ddgtags.push_back( line );
		parse_ddgs ( ddgtags );
	}


	// write start score
    //evaluator(startV, 1);

	int NCYCLES = 1;
	std::vector< double > final_ref(20,0);
	double final_score=0;
	for (int i=1; i<=NCYCLES; ++i) {
		NelderMead optimizer;
		optimizer.run( workingV, scaleV, &evaluator, 5e-4, MAXITER );
		optimizer.recover_best(final_ref,final_score);
		workingV = final_ref;
	}

	// write final score
	evaluator(final_ref, 1);

	double sum = get_sum(final_ref) - final_ref[0];
    //std::cout << std::endl << "REF:" << std::endl;
	//for (int i=0; i<20; ++i) {
    //    std::cout << idx2aa(i) << ": " << startV[i+1] << " -> " << final_ref[i+1] - sum/20 << std::endl;
    //}

	std::cout << std::endl <<  "METHOD_WEIGHTS ref " <<std::setprecision(5);
	for (int i=0; i<20; ++i) {
		std::cout << final_ref[i+1] - sum/20 << " ";
	}
	std::cout << "\n";

    double wt_freq = final_ref[0];
    if (wt_freq != 0) {
        std::cout << std::endl <<  "[biased " << wt_freq << "] METHOD_WEIGHTS ref " <<std::setprecision(5);
        for (int i=0; i<20; ++i) {
            std::cout << final_ref[i+1] - sum/20  + wt_freq*log(get_freq(i)) << " ";
        }
        std::cout << "\n";
    }
}

