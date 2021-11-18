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
  int process_no; // index in the processes table
  int parent_no; // index of the parent in the processes table
  int generation; // parent genration + 1
  int seq_no; // sequence no. in the current child process
  double parent_start;
  Process(int n=0, int p=-1, int g=0, int s = 0, double t=0.0):
    process_no(n),parent_no(p),generation(g),seq_no(s), parent_start(t) {}
  Process& operator=(const Process& p){
    process_no = p.process_no;
    parent_no = p.parent_no;
    generation = p.generation;
    seq_no = p.seq_no;
    parent_start = p.parent_start;
    return *this;
  }
};

ostream& operator<<(ostream& os, const Process& p){
  os << format("process no. %d: parent_no: %d generation: %d sequence_no: %d start: %f",
               p.process_no,p.parent_no,p.generation,p.seq_no,p.parent_start);
  return os;
}

struct HawkesPoint{
  Process mark; // mark is the generating process for this event
  double time;
  bool operator< (const HawkesPoint& h) const { return time < h.time;}
  HawkesPoint(Process m = 0, double t = 0.0) : mark(m), time(t) {}
  HawkesPoint& operator=(const HawkesPoint& h){
    mark = h.mark;
    time = h.time;
    return *this;
  }
};

ostream& operator<<(ostream& os, const HawkesPoint& h){
  os << format("mark: %s time: %f\n",h.mark,h.time);
  return os;
}
  
 
int main(int argc, char** argv){

  string model_input_file = "";
  string output_file = "hawkes_test_data.txt";
  string file_dir = "";
  int ndata = 100;
  double sigma_0 = 2;
  double half_life = 2;
  double lambda_0 = 1;
  double max_time = 5;
  int seed = 12345;
  
  GetOpt cl(argc,argv); // parse command line
  cl.get("output_file",output_file); cout << "model output to "<<output_file<<endl;
  cl.get("file_dir",file_dir); cout << "file directory is "<<file_dir<<endl;
  cl.get("ndata",ndata); cout << "ndata: "<<ndata<<endl;
  cl.get("sigma_0",sigma_0); cout << "sigma_0: "<<sigma_0<<endl;
  cl.get("half_life",half_life); cout << "half_life: "<<half_life<<endl;
  cl.get("lambda_0",lambda_0); cout << "lambda_0: "<<lambda_0<<endl;
  cl.get("max_time",max_time); cout << "max_time: "<<max_time<<endl;
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
  //  Array<Process> processes(10000);
  int nparents = 0; // no. of parents processed so far
  int npoints = 0; // no. of data points stored so far (store next data point at data[npoints]
  Process fake; // fake parent for base process
  data[npoints++] = HawkesPoint();
  //  processes[0] = fake;
  LFSR128 rng(seed);
  
  while(nparents < npoints){
    Process& parent = data[nparents].mark; // the child becomes a parent
    double t0 = data[nparents].time;
    Process& this_process = data[npoints].mark;
    if(nparents > 0) this_process = Process(nparents,parent.process_no,parent.generation+1,0,t0);
    else this_process = fake;
    cout << "new process: "<<this_process <<" time: "<<t0<<endl;
    //    processes[nparents] = this_process;
    double t = t0;
    while(true){
      double sigma = sigma_0;
      double lambda = this_process.process_no == 0? lambda_0 : sigma*exp(-rho*(t-t0))/pow(2,this_process.generation);
      t -= log(1-rng.uniform())/lambda; // delta_t is exponentially distributed with rate lambda
      if(t > max_time) break;
      this_process.seq_no++;
      data[npoints++] = HawkesPoint(this_process,t);
    }
    cout << format("process %d has %d points\n",this_process.process_no,this_process.seq_no);
    nparents++;
  } // OK, create the next child process
  data.sort(npoints);
  cout << "sorted (but not re-indexed) data:\n";
  for(int i = 0;i < npoints;i++) cout << data[i].mark<<" "<<data[i].time<<endl;
  Array<int> new_process_no(10000);
  new_process_no.fill(-1);
  int new_i;
  for(int i = 0;i < npoints-1;i++){ // re-number processes after sorting
    Process m = data[i+1].mark;
    if(new_process_no[m.process_no] == -1) new_process_no[m.process_no] = i; // start of new process
    int pno = (m.parent_no >= 0 ? m.parent_no : 0);
    output<<new_process_no[pno]<<" "<<data[i+1].time<<endl;
  }
  output.close();
  for(int i = 0;i < new_process_no.len();i++){
    cout << format("%d -> %d\n",i,new_process_no[i]);
  }
}  
