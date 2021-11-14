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
  int no; // index in the processes table
  int parent_no; // index of the parent in the processes table
  int generation; // parent genration + 1
  double start_time;
  Process(int n=0, int p=-1, int g=1, double t=0.0): no(n),parent_no(p),generation(g),start_time(t) {}
  Process& operator=(const Process& p){
    no = p.no;
    parent_no = p.parent_no;
    generation = p.generation;
    start_time = p.start_time;
    return *this;
  }
};

ostream& operator<<(ostream& os, const Process& p){
  os << format("process no. %d: parent_no: %d generation: %d start: %f\n",p.no,p.parent_no,p.generation,p.start_time);
  return os;
}

struct HawkesPoint{
  int mark; // processes[mark] is the generating process for this event
  double time;
  bool operator< (const HawkesPoint& h) const { return time < h.time;}
  HawkesPoint(int m = 0, double t = 0.0) : mark(m), time(t) {}
  HawkesPoint& operator=(const HawkesPoint& h){
    mark = h.mark;
    time = h.time;
    return *this;
  }
};

ostream& operator<<(ostream& os, const HawkesPoint& h){
  os << format("mark: %d time: %f\n",h.mark,h.time);
  return os;
}
  
 
int main(int argc, char** argv){

  string model_input_file = "";
  string output_file = "hawkes_test_data.txt";
  string file_dir = "";
  int ndata = 100;
  double sigma_0 = 1;
  double half_life = .1;
  double lambda_0 = 1;
  double total_time = 20;
  int seed = 12345;
  
  GetOpt cl(argc,argv); // parse command line
  cl.get("output_file",output_file); cout << "model output to "<<output_file<<endl;
  cl.get("file_dir",file_dir); cout << "file directory is "<<file_dir<<endl;
  cl.get("ndata",ndata); cout << "ndata: "<<ndata<<endl;
  cl.get("sigma_0",sigma_0); cout << "sigma_0: "<<sigma_0<<endl;
  cl.get("half_life",half_life); cout << "half_life: "<<half_life<<endl;
  cl.get("lambda_0",lambda_0); cout << "lambda_0: "<<lambda_0<<endl;
  cl.get("seed",seed); cout << "seed: "<<seed<<endl;
  
  double rho = log(2)/half_life; // decay rate of child processes
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  if(file_dir !=""){
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
  //  int next_point = 0;
  Process fake; // fake parent for base process
  data.add(HawkesPoint());
  processes[0] = fake;
  LFSR128 rng(seed);
  
  while(next_process < data.nitems()){
    Process& parent = processes[data[next_process].mark]; // the child becomes a parent
    double t0 = data[next_process].time;
    cout << "new parent: "<<parent<<" time: "<<t0<<endl;
    Process child(next_process,parent.no,parent.generation+1,t0);
    processes[next_process] = child;
    double t = t0;
    while(true){
      double sigma = sigma_0;
      double lambda = sigma*exp(-rho*(t-t0));
      
      t -= log(1-rng.uniform())/lambda; // delta_t is exponentially distributed with rate lambda
      if(t > total_time) break;
      data.add(HawkesPoint(next_process,t));
    }
    next_process++;
  } // OK, create the next child process
  data.sort();
  
  for(int i = 0;i < data.nitems();i++){
    output<<data[i].mark<<" "<<data[i].time<<endl;
  }
  output.close();
}  
