// testing merge
#include<iostream>
#include<fenv.h>
#include<fstream>
#include<cmath>
#include<string>
#include "using.h"
#include "Awk.h"
#include "Array.h"
#include "Matrix.h"
#include "GetOpt.h"
#include "util.h"
#include "stats.h"
#include "Dict.h"

bool dump_flag(false);
bool verbose(false);

struct HawkesPoint{
  string mark;
  double time;
};

ostream& operator<<(ostream& os, const HawkesPoint& h){
  os << format("mark: %s time: %f\n",h.mark.c_str(),time);
  return os;
}

        
class SigmaRho{
  Matrix<double>& omega1; // inputs
  Array<HawkesPoint>& data;
  double lambda;
  
  double& rho; // initial values and outputs. Set at construction
  double& _sigma;
  int N;
  
  double t(int i){return data[i].time;}

  double sigma(double rho){
    if(rho == 0) {
      double sum = 0;
      for(int j = 1;j < N;j++){
        sum += 1 - t(j);
        // if(sum < 0){
        //   cout << format("negative sum at t(%d) = %f\n",j,t(j));
        // }
      }
      cout << format("sigma(0): %f, sum: %f\n",sum,(N-lambda)/sum);
      return (N-lambda)/sum;
    } // r > 0
    double sum = 0; // r > 0 here
    for(int j = 1;j < N;j++)sum += 1-exp(-rho*(1-t(j)));
    return (N-lambda)*rho/sum;
  }

  double dsigma(double rho, double _sigma){
    double S_0 = 0;
    double dS_0 = 0;
    for(int j = 1;j < N;j++){
      double t_Nj = t(N) - t(j);
      double e_Nj = exp(-rho*t_Nj);
      S_0 += 1-e_Nj;
      dS_0 += t_Nj*e_Nj;
    }
    return ((N-lambda)-_sigma*dS_0)/S_0;
  }
      
  double dQ(double rho){ // compute dQ_drho 
    if(rho == 0){
      cerr << "dQ: rho = 0. Bailing out\n";
      exit(1);
    }
    double _sigma = sigma(rho);

    double _dsigma_drho = dsigma(rho,_sigma);
    double dQ_drho = 0;
    for(int i = 2;i <= N;i++){
      for(int j = 1;j < i;j++){
        double t_ij = t(i)-t(j);
        double e_ij = exp(-rho*t_ij);
        double k_hat = _sigma/rho*(1-e_ij);
        dQ_drho += (omega1(i,j)/(2*k_hat))*(_dsigma_drho - k_hat + _sigma*t_ij*e_ij);
      }
    }
    dQ_drho /= rho;
    return dQ_drho; 
  }

  double real_Q(double r){ // testing only
    double sum = 0;
    double k;
    double s = sigma(r);
    for(int i = 2;i <= N;i++){
      for(int j = 1;j < i;j++){
        if(r == 0) k = s*(t(i)-t(j));
        else k = (s/r)*(1-exp(-r*(t(i)-t(j))));
        sum += omega1(i,j)*(k*log(k) - k - gammln(k));
      }
    }
    return sum;
  }
  
  double Q_comp(double r){ // approximate version of real_Q
    double sum = 0;
    double s = sigma(r);
    for(int i = 2;i < N;i++){
      for(int j = 1;j < i;j++){
        if(r == 0) sum += omega1(i,j)*log(s*(t(i)-t(j)));
        else sum += omega1(i,j)*log((s/r)*(1-exp(-r*(t(i)-t(j)))));
      }
    }
    return sum/2;
  }
  
public:
  SigmaRho(int n, Matrix<double>& om, Array<HawkesPoint>& d,
           double l, double& r, double& s) :
    omega1(om),data(d),lambda(l),rho(r),_sigma(s),N(n) {}

  void solve(int max_iters, double eps){
    int niters = 0;
    double f_min = dQ(.001);
    double f_max = dQ(1);
    cout << format("f(.001) = %f, f(1) = %f\n",f_min, f_max);
    double r_min = 0;
    double r_max = 1;
    double new_r, new_f;

    //step1: find an r s.t. f_min & f_max have opp. sign
    output("rho_test.plt");
    exit(0);
    
    while(niters++ < max_iters){
      if(f_min*f_max <= 0) break;
      cout << format("dQ(%f) = %f\n",r_max,dQ(r_max));
      r_max += .1;
      f_max = dQ(r_max);
    }
    if(f_min*f_max > 0){
      cerr << "failed to find r_max after "<<niters<<" iterations."<<endl;
      output("rho_test.plt");
      exit(1);
    }
    else cout << "rmax = "<<r_max<<", f_max = "<<f_max<<endl;
        
    niters = 0;
    while(niters++ < max_iters && fabs(r_min-r_max) > eps){
      new_r = (r_min+r_max)/2;
      new_f = dQ(new_r);
      //      cout << format("sr_iteration %d: r: %f f: %f\n",niters,new_r,new_f);
      if(f_min*new_f > 0) {
        f_min = new_f;
        r_min = new_r;
      }
      else{
        f_max = new_f;
        r_max = new_r;
      }
    }
    if(niters >= max_iters){
      cerr << "SigmaRho.solve did not converge after "<<niters<<
        " iterations with error = "<<fabs(r_min-r_max)<<endl;
      exit(1);
    }
    rho = new_r;
    _sigma = sigma(rho);
  }

  void output(char* file, double max_rho = 10){
    ofstream out(file);
    if(!out.good()){
      cerr << "Can't open "<<file<<endl;
      exit(1);
    }
    double r_temp;
    out << "omega1:\n"<<omega1;
    out << format("   rho\t   dQ_drho  Q(rho)    real_Q\n");
    for(r_temp = .1; r_temp < max_rho;r_temp += .1){
      out << format("%8.4f %8.4f %8.4f %8.4f\n", r_temp, dQ(r_temp), Q_comp(r_temp), real_Q(r_temp));
    }
    out.close();
  }

};


struct Mark{
  string name;
  double rho;
  double sigma;
  int nkids;
  Array<int> kids;

  Mark(void) : name(""),rho(0),sigma(0),nkids(0),kids(10) {}
};

ostream& operator<<(ostream& os, const Mark& m){
  os << format("%s: rho: %f sigma:%f children: ",m.name.c_str(),m.rho,m.sigma);
  for(int i = 0;i <m.nkids;i++)os << m.kids[i] <<", ";
  os <<endl;
  return os;
}
 
int main(int argc, char** argv){

  string data_file = "hawkes_test_data.txt";
  string model_input_file = "";
  string model_output_file = "model.out";
  string file_dir = "";
  int ndata = 100;
  double sigma_0 = 1;
  double rho_0 = log(2)/.1; // nominal half-life = .1
  double lambda_0 = 50; // nominal no. of base-process events in [0,1].
  int max_iters = 100;
  int simulation = 1;
  int sr_iters = 50;
  double sr_err = 1.0e-8;
  
  GetOpt cl(argc,argv); // parse command line
  cl.get("data_file",data_file); cout << "data file: "<<data_file<<endl;
  cl.get("model_output_file",model_output_file); cout << "model output to "<<model_output_file<<endl;
  cl.get("file_dir",file_dir); cout << "file directory is "<<file_dir<<endl;
  cl.get("ndata",ndata); cout << "ndata: "<<ndata<<endl;
  cl.get("sigma_0",sigma_0); cout << "sigma_0: "<<sigma_0<<endl;
  cl.get("rho_0",rho_0); cout << "rho_0: "<<rho_0<<endl;
  cl.get("lambda_0",lambda_0); cout << "lambda_0: "<<lambda_0<<endl;
  cl.get("max_iters",max_iters); cout << "max iterations: "<<max_iters<<endl;
  cl.get("simulation",simulation); cout << "simulation mode: "<<simulation<<endl;
  cl.get("sr_iters",sr_iters);cout<<"sr_iters: "<<sr_iters<<endl;
  cl.get("sr_err",sr_err); cout <<"sr_err: "<<sr_err<<endl;
  
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  if(file_dir !=""){
    data_file = file_dir+"/"+data_file;
    model_input_file = file_dir+"/"+model_input_file;
    model_output_file = file_dir+"/"+model_output_file;
  }
  std::ofstream model_output_stream(model_output_file.c_str());
  if(model_output_file != "" && !model_output_stream.good()){
    cerr << "can't open "<<model_output_file<<endl;
    exit(1);
  }
  Array<HawkesPoint> data(10000);
  Dict<string,Mark> marks;
  
  
  Awk awk;
  if(!awk.open(data_file.c_str())){
    cerr << "Can't open "<<data_file<<endl;
    exit(1);
  }

  int i = 0;
  while((awk.nf = awk.next()) != -1 && (ndata == 0?  true : i <= ndata)){ // exit on EOF
    if(awk.nf != 2) continue;
    data[i].mark = string(awk[1]);
    Mark& m = marks[data[i].mark];
    if(m.name == ""){
      m.name = data[i].mark;
      m.sigma = sigma_0;
      m.rho = rho_0;
    }
    m.kids[m.nkids++] = i;
    data[i].time = atof(awk[2]);
    i++;
  }
  int N = i-1; // t[0] is a reference time, not an event!  
  for(int i = 0;i <= N;i++) data[i].time = (data[i].time-data[0].time)/data[N].time; // normalize times to [0,1].


  /* The base process starts at 0 and generates the first arrival
   * at t[1], as well as others throughout the time interval [0,1]   
   * the total number of arrivals is N.  Each arrival (except t[N]) generates
   * a child process, which *may* generate one or more events, so there are
   * N-1 child processes.   
   */
  

  Array<double> omega(N); // the overall probability of state j
  Matrix<double> omega1(N+1,N); // posterior prob of state j at time i
  Matrix<double> k_hat_0(N+1,N); // prior "
  
  //  Matrix<double> k_hat(N,N); // posterior expected # events in state j up to time i
  // ColVector<double> A_k(2);
  // Matrix<double> AtA(2,2); // normal matrix
  // ColVector<double> Atb(2);
  // ColVector<double> delta(2);

  // set parameters to nominal
  double rho = rho_0;
  double sigma = sigma_0;
  double lambda = lambda_0;
  double sum = N*(N+1);
  for(int j = 0;j < N;j++) omega.fill((N - j)/sum);
  
  // begin EM iteration
  SigmaRho sr(N,omega1,data,lambda,rho,sigma);
  int niters = 0;
  while(niters <= max_iters){
    // compute posterior probability matrix omega1
    k_hat_0.fill(0.0);
    for(int i = 0;i < marks.nkeys();i++){ // set all the marks data
      marks(i).sigma = sigma;
      marks(i).rho = rho;
    }
    omega1(0,0) /*= k_hat(0,0) = k_hat_0(0,0)*/ = 1.0;
    double log_likelihood = 0.0;
    double score = 0.0;
    for(int i = 1;i <= N;i++){
      double t_i = data[i].time;
      double t_ij;
      double row_sum = 0.0;
      int mode = 0; // max_j omega1(i,j)
      for(int j = 0; j < i;j++){
        double sigma = marks[data[i].mark].sigma;
        double rho = marks[data[i].mark].rho;
        if(j == 0){
          k_hat_0(i,0) = lambda*t_i;
          t_ij = t_i;
        }
        else{
          t_ij = t_i - data[j].time;
          k_hat_0(i,j) = sigma*(1-exp(-rho*t_ij))/rho;
        }
        double& kk = k_hat_0(i,j);
        omega1(i,j) = omega[j]*exp(kk*(log(kk)-1) - gammln(kk)-log(t_ij));
        assert(omega1(i,j) > 0);
        if(omega1(i,j) > omega1(i,mode)) mode = j;
        row_sum += omega1(i,j);
      }
      log_likelihood += log(row_sum);
      for(int j = 0;j < N;j++) omega1(i,j) /= row_sum;
      int truth = simulation? atoi(data[i].mark.c_str()) : -1;
      //      cout << format("i = %d: mode: process %d (%f) truth: process %d\n",i,mode,omega1(i,mode),truth);
      if(simulation)score += log(i*omega1(i,truth));
    }
    cout << "log likelihood at iteration "<<niters<<": "<<log_likelihood<<endl;
    cout << format("scoring rate: %f bits/obs\n",score/((N-1)*log(2)));
    //    sr.output("rho_test.plt");
    //    exit(0);
    if(niters < max_iters){
      // now we re-estimate the parameters
      cout << "\nbegin iteration "<<niters<<endl;
      for(int j = 0;j < N;j++){
        omega[j] = 0.0;
        for(i = j+1;i <= N;i++) omega[j] += omega1(i,j);
        omega[j] /= N;
      }
      lambda = omega[0]/data[N].time;//Note: t[N] = 1. This is for clarity
      cout << "update for lambda: "<<lambda<<endl;
      sr.solve(sr_iters,sr_err); // updates sigma and rho
      cout << format("at iteration %d, lambda = %f, sigma = %f, rho = %f\n",
                     niters,lambda,sigma,rho);
    } // end re-estimation
    niters++;
  } // on to the next EM iteration
  // output results here
  cout<<format("lambda: %f, sigma: %f, rho: %f\n",lambda,sigma,rho); 
}  
