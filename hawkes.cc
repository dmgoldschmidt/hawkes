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
  os << format("mark: %s time: %f\n",mark.c_str(),time);
  return os;
}

struct Mark{
  string name;
  double rho;
  double sigma;
  int nkids;
  Array<int> kids;

  Mark(void) : name(""),rho(0),sigma(0),nkids(0),kids(10) {}
};
ostream& operator<<(ostream& os, const Mark& m){
  os << format("%s: rho: %f sigma:%f children: ",name.c_str(),rho,sigma);
  for(int i = 0;i <nkids;i++)os << kids[i],", ";
  os <<endl;
  return os;
}

  
 
int main(int argc, char** argv){

  string data_file = "marked_events.txt";
  string model_input_file = "";
  string model_output_file = "model.out";
  string file_dir = "";
  int ndata = 100;
  double sigma_0 = 2.0;
  double rho_0 = 1;
  double lambda_0 = .1;
  int max_iters = 10;

  GetOpt cl(argc,argv); // parse command line
  cl.get("data_file",data_file); cout << "data file: "<<data_file<<endl;
  cl.get("model_output_file",model_output_file); cout << "model output to "<<model_output_file<<endl;
  cl.get("file_dir",file_dir); cout << "file directory is "<<file_dir<<endl;
  cl.get("ndata",ndata);
  cl.get("sigma_0",sigma_0);
  cl.get("rho_0",rho_0);
  cl.get("lambda_0",lambda_0);
  cl.get("max_iters",max_iters);
  
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
  if(!awk.open(data_file)){
    cerr << "Can't open "<<data_file<<endl;
    exit(1);
  }
  int i = 0;
  while((nf = awk.next()) != -1 && (ndata == 0?  true : i < ndata)){ // exit on EOF
    if(nf != 2) continue;
    data[i].mark = string(awk[1]);
    Mark m = marks[data[i].mark];
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
  lambda = lambda_0;
  omega.fill(1.0/ndata);

  // begin EM iteration
  niters = 0;
  while(niters <= max_iters){
    // compute posterior probability matrix omega1
    k_hat.fill(0.0);
    k_hat_0.file(0.0);
    omega1(0,0) = k_hat(0,0) = k_hat_0(0,0) = 1.0;
    double log_likelihood = 0.0;
    for(i = 1;i < ndata;i++){
      double t_i = data[i].time;
      double t_ij;
      row_sum = 0.0;
      for(j = 0; j < i;j++){
        double sigma = marks[data[i].mark].sigma;
        double rho = marks[data[i].mark].rho;
        if(j == 0){
          k_hat_0(i,0) = lambda*t_i;
          t_ij = t_i;
        }
        else{
          t_ij = (j == 0 ? t_i : t_i - data[j].time);
          k_hat_0(i,j) = sigma*(1-exp(-rho*t_ij))/rho;
        }
        omega1(i,j) = omega[j]*sqrt(k_hat_0(i,j))/t_ij;
        row_sum += omega1(i,j);
      }
      log_likelihood += log(row_sum);
      for(j = 0;j < ndata;j++) omega(i,j) /= row_sum;
    }
    cout << "log likelihood at iteration "<<niters<<": "<<log_likelihood<<endl;
    if(niters < max_iters){
      // now we re-estimate the parameters
      cout << "begin iteration "<<niters<<endl;
      for(j = 0;j < ndata;j++){
        omega[j] = 0.0;
        for(i = j;i < ndata;i++) omega[j] += omega1(i,j);
        omega[j] \= ndata;
      }
      lambda = 0.0;
      for(i = 0;i < ndata;i++) lambda += omega1(i,0)/data[i].time;
      lambda /= ndata;

      // now re-estimate rho and sigma using non-linear least squares
      for(int m = 0;m < marks.nkeys();m++){
        Mark mark = marks(m);
        if(mark.nkids < 2) continue;
        double sigma& = mark.sigma;
        double rho& = mark.rho;
        double rms_error = 0;
        AtA.fill(0);
        Atb.fill(0);
        for(int j = 0;j < mark.nkids;j++){ // get the equations for this mark
          int j = mark.kids[m]; // next child for this mark
          for(int i = j+1;i < ndata;i++){ // loop over all subsequent events
            double t_ij = data[i].time - data[j].time;
            A_k[0] = k_hat_0(i,j)/sigma;
            A_k[1] = t_ij*(sigma/rho - k_hat_0(i,j)) - k_hat_0(i,j)/rho;
            double b = k_hat(i,j) - k_hat_0(i,j); // residual
            AtA += A_k*A_k.Tr();
            Atb += b*A_k;
            rms_error += b*b;
          }
        }
        delta = AtA.inv()*Atb;
        double f = .1*sigma/fabs(delta[0]); // bound changes to .1 of originals
        if(f < 1) delta[0] *= f;
        sigma += delta[0];
        f = .1*rho/fabs(delta[1]);
        if(f < 1) delta[0] *= f;
        rho += delta[1];
      } // on to the next mark
    } // end re-estimation
  } // on to the next EM iteration
  // output results here
}  
