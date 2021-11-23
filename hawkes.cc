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

// struct SigRho{
//   double sigma;
//   double rho;
// };

// struct Parameters{
//   double lambda;
//   Array<double> omega;
//   Dict<string,SigRho> decay_params;

//   Parameters(int block_size = 100):omega(block_size){}
// };
// ostream& operator<<(ostream& os, const Parameters& p){
//   os <<format("\nlambda: %f decay_paramers:\n",p.lambda);
//   for(int i= 0;i < p.decay_params.nkeys();i++){
//     os << format("mark: %s sigma: %f rho: %f\n",decay_params.key(i).c_str(),
//                  decay_params(i).sigma,decay_params(i).rho);
//   }
//   return os;
// }

struct HawkesPoint{
  string mark;
  double time;
};
ostream& operator<<(ostream& os, const HawkesPoint& h){
  os << format("mark: %s time: %f\n",h.mark.c_str(),time);
  return os;
}

struct FullMark{
  string name;
  double rho;
  double sigma;
  int nkids;
  Array<int> kids;

  FullMark(void) : name(""),rho(0),sigma(0),nkids(0),kids(10) {}
};
ostream& operator<<(ostream& os, const FullMark& m){
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
  double sigma_0 = 2.0;
  double rho_0 = 1;
  double lambda_0 = .1;
  int max_iters = 5;
  int simulation = 1;

  GetOpt cl(argc,argv); // parse command line
  cl.get("data_file",data_file); cout << "data file: "<<data_file<<endl;
  cl.get("model_output_file",model_output_file); cout << "model output to "<<model_output_file<<endl;
  cl.get("file_dir",file_dir); cout << "file directory is "<<file_dir<<endl;
  cl.get("ndata",ndata); cout << "ndata: "<<ndata<<endl;
  cl.get("sigma_0",sigma_0); cout << "sigma_0: "<<sigma_0<<endl;
  cl.get("rho_0",rho_0); cout << "rho_0: "<<rho_0<<endl;
  cl.get("lambda_0",lambda_0); cout << "lambda_0: "<<lambda_0<<endl;
  cl.get("max_iters",max_iters); cout << "max iterations: "<<max_iters<<endl;
  cl.get("simulation",simulation); cout << "simulation: "<<simulation<<endl;
  
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
  Dict<string,FullMark> full_marks;
  
  
  Awk awk;
  if(!awk.open(data_file.c_str())){
    cerr << "Can't open "<<data_file<<endl;
    exit(1);
  }
  int i = 0;
  while((awk.nf = awk.next()) != -1 && (ndata == 0?  true : i < ndata)){ // exit on EOF
    if(awk.nf != 2) continue;
    data[i].mark = string(awk[1]);
    if(data[i].mark == "0") data[i].mark = "base";
    else data[i].mark = "1";
    FullMark& m = full_marks[data[i].mark];
    if(m.name == ""){
      m.name = data[i].mark;
      m.sigma = sigma_0;
      m.rho = rho_0;
    }
    m.kids[m.nkids++] = i;
    data[i].time = atof(awk[2]);
    i++;
  }
  ndata = i;
  Array<double> omega(ndata);
  Matrix<double> omega1(ndata,ndata); // posterior prob of state j at time i
  Matrix<double> k_hat(ndata,ndata); // posterior expected # events in state j up to time i
  Matrix<double> k_hat_0(ndata,ndata); // prior "
  ColVector<double> A_k(2);
  Matrix<double> AtA(2,2); // normal matrix
  ColVector<double> Atb(2);
  ColVector<double> delta(2);
  double lambda = lambda_0;
  double sum = ndata*(ndata+1);
  for(int i = 0;i < ndata;i++) omega.fill((ndata - i)/sum);
  
  // begin EM iteration
  int niters = 0;
  while(niters <= max_iters){
    // cout << "begin omega computation for iteration "<<niters<<endl;
    // compute posterior probability matrix omega1
    k_hat_0.fill(0.0);
    omega1(0,0) = k_hat(0,0) = k_hat_0(0,0) = 1.0;
    double log_likelihood = 0.0;
    double score = 0.0;
    for(int i = 1;i < ndata;i++){
      double t_i = data[i].time;
      double t_ij;
      double row_sum = 0.0;
      int mode = 0;
      for(int j = 0; j < i;j++){
        double& sigma = full_marks[data[i].mark].sigma;
        double& rho = full_marks[data[i].mark].rho;
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
        if(omega1(i,j) > omega1(i,mode)) mode = j; // find the largest entry in row i
        row_sum += omega1(i,j);
      }
      
      log_likelihood += log(row_sum);
      for(int j = 0;j < ndata;j++) omega1(i,j) /= row_sum;
      int truth;
      if(simulation){
        truth = data[i].mark == "base"? 0 : atoi(data[i].mark.c_str()); // convert the string to an integer
                                                                        //(only for simulation data)
        if(i > 1)score += log(i*omega1(i,truth));
      }
    }
    cout << "log likelihood at iteration "<<niters<<": "<<log_likelihood<<endl;
    if(simulation) cout << format("scoring rate at iteration %d: %f bits/obs\n",niters,score/((ndata-2)*log(2)));

    if(niters < max_iters){
      // now we re-estimate the parameters
      cout << "begin re-estimation for iteration "<<niters<<endl;
      k_hat.fill(0.0);
      for(int j = 0;j < ndata;j++){
        omega[j] = 0.0;
        for(i = j+1;i < ndata;i++){
          omega[j] += omega1(i,j);
          k_hat(i,j) = k_hat(i-1,j) + omega1(i,j);
        }
        omega[j] /= ndata;
      }
      lambda = 0.0;
      for(i = 0;i < ndata;i++) lambda += omega1(i,0)/data[i].time;
      lambda /= ndata;
      cout << "update for lambda: "<<lambda<<endl;

      // now re-estimate rho and sigma using non-linear least squares
      for(int m = 0;m < full_marks.nkeys();m++){
        FullMark& mark = full_marks(m); 
        if(mark.nkids < 2 || mark.name == "base") continue; // this only applies to children
        cout << "processing FullMark "<<mark.name<<": "<<mark;
        double& sigma = mark.sigma;
        double& rho = mark.rho;
        cout << "re-estimation is using sigma = "<<sigma<<endl;
        double rms_error = 0;
        int neq = 0;
        AtA.fill(0);
        Atb.fill(0);
        for(int k = 0; k < mark.nkids;k++){ // get the equations for this mark
          int j = mark.kids[k]; // next child for this mark
          for(int i = j+1;i < ndata;i++){ // loop over all subsequent events
            double t_ij = data[i].time - data[j].time;
            A_k[0] = k_hat_0(i,j)/sigma;
            A_k[1] = t_ij*(sigma/rho - k_hat_0(i,j)) - k_hat_0(i,j)/rho;
            double b = k_hat(i,j) - k_hat_0(i,j); // residual
            //            cout << k_hat(i,j)<<" "<<k_hat_0(i,j)<<" "<<b<<endl;
            AtA += A_k*A_k.Tr();
            Atb[0] += A_k[0]*b; Atb[1] += A_k[1]*b;
            rms_error += b*b;
            neq++;
          }
        }
        Matrix<double> AtAinv = inv(AtA);
        delta = AtAinv*Atb;
        cout << "AtAinv:\n"<<AtAinv<<"Atb:\n"<<Atb<<endl;
        cout << "raw delta:\n"<<delta<<endl;
        rms_error = sqrt(rms_error/neq);
        double f = .1*sigma/fabs(delta[0]); // bound changes to .1 of originals
        if(f < 1) delta[0] *= f;
        sigma += delta[0];
        f = .1*rho/fabs(delta[1]);
        if(f < 1) delta[1] *= f;
        rho += delta[1];
        cout << "bounded delta:\n"<<delta;
        cout << format("update for mark %s at iteration %d: sigma: %f rho: %f rms_error: %f\n",
                       mark.name.c_str(),niters,mark.sigma,mark.rho,rms_error);
      } // on to the next mark
    } // end re-estimation
    niters++;
    //    if(niters > 2) exit(0);
  } // on to the next EM iteration
  // output results here
}  
