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
#include "FlexHeap.h"



struct Process{
  int no; // index in the data table
  int parent_no; // index of the parent in the data table
  int generation; // parent genration + 1
  double start_time;
};

ostream& operator<<(ostream& os, const Process& p){
  os << format("process no. %d: parent: %d generation: %d start: %f\n",p.no,p.parent,p.generation,p.start_time);
  return os;
}

struct HawkesPoint{
  int mark; // processes[mark] is the generating process for this event
  double time;
  bool operator< (const HawkesPoint& h){ return time < h.time;}
};

ostream& operator<<(ostream& os, const HawkesPoint& h){
  os << format("mark: %s time: %f\n",h.mark.c_str(),time);
  return os;
}
  
 
int main(int argc, char** argv){

  string model_input_file = "";
  string output_file = "hawkes_test_data.txt";
  string file_dir = "";
  int ndata = 100;
  double sigma_0 = 1
  double half_life = .1;
  double lambda_0 = 1;
  double total_time = 20;

  GetOpt cl(argc,argv); // parse command line
  cl.get("data_file",data_file); cout << "data file: "<<data_file<<endl;
  cl.get("model_output_file",model_output_file); cout << "model output to "<<model_output_file<<endl;
  cl.get("file_dir",file_dir); cout << "file directory is "<<file_dir<<endl;
  cl.get("ndata",ndata); cout << "ndata: "<<ndata<<endl;
  cl.get("sigma_0",sigma_0); cout << "sigma_0: "<<sigma_0<<endl;
  cl.get("half_life",rho_0); cout << "half_life: "<<half_life<<endl;
  cl.get("lambda_0",lambda_0); cout << "lambda_0: "<<lambda_0<<endl;
  cl.get("max_iters",max_iters); cout << "max iterations: "<<max_iters<<endl;

  rho = log(2)/half_life; // decay rate of child processes
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  if(file_dir !=""){
    data_file = file_dir+"/"+data_file;
    model_input_file = file_dir+"/"+model_input_file;
    output_file = file_dir+"/"+output_file;
  }
  std::ofstream output(output_file.c_str());
  if(output_file != "" && !output.good()){
    cerr << "can't open "<< output_file<<endl;
    exit(1);
  }
  
  FlexHeap<HawkesPoint> data(10000);
  Array<Process> processes(10000);
  int next_process = 0;
  int next_point = 0;
  Process fake(0,-1,0,0.0); // fake parent for base process
  data[0] = HawkesPoint(0,0.0);
  processes[0] = fake;
  LFSR128 rng(seed);
  
  while(next_point < data.len()){
    Process& parent = processes[data[next_point].mark]; // the child becomes a parent
    double t0 = data[next_point].time;
    Process child(next_process,parent.parent_no,parent.generation+1,t0);
    processes[next_process] = child;
    double t = t0;
    while(true){
      lambda = sigma*exp(-rho*(t-t0));
      t -= log(1-rng.uniform())/lambda; // delta_t is exponentially distributed with rate lambda
      if(t > total_time) break;
      data[next_point++] = HawkesPoint(next_process,t);
    }
    next_process++;
  } // OK, on to create the next child process
  data.sort();
  
  for(int i = 0;i < data.nitems();i++){
    output<<data[i].mark<<" "<<data[i].time<<endl;
  }
}  
