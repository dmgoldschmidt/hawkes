INC = /home/david/code/include
SRC = /home/david/code/src
CFLAGS = -g -O0 -Wall -I/home/david/code/include -std=c++14 -fno-stack-protector
CFLAGS1 = -03 -I/home/david/code/include -std=c++14
VPATH =  /home/david/code/include
CC = g++
#CC = clang++ -DFILAMENT_REQUIRES_CXXABI=true
 
util.o: util.h util.cc
	$(CC) $(CFLAGS) -c -o util.o util.cc -lz
GetOpt.o: GetOpt.h GetOpt.cc
	$(CC) $(CFLAGS) -c -o GetOpt.o GetOpt.cc
Awk.o:	Awk.h Awk.cc util.o 
	$(CC) $(CFLAGS) -c -o Awk.o Awk.cc -lz 
test_gz: test_gzstream.cc gzstream.o Awk.o util.o GetOpt.o
	$(CC) $(CFLAGS) -o test_gz test_gzstream.cc gzstream.o Awk.o util.o GetOpt.o -lz
read_lanl: read_lanl.cc gzstream.o Awk.o util.o GetOpt.o matrix.o GammaFunc.o Array.h KeyIndex.h Histogram.h GammaFunc.h
	$(CC) $(CFLAGS) -o read_lanl read_lanl.cc gzstream.o Awk.o util.o GetOpt.o matrix.o GammaFunc.o -lz
matrix.o: matrix.cc Matrix.h Array.h
	$(CC) $(CFLAGS) -c matrix.cc
digamma.o: digamma.cc
	$(CC) $(CFLAGS) -c digamma.cc
GammaFunc.o: GammaFunc.cc GammaFunc.h
	$(CC) $(CFLAGS) -c GammaFunc.cc
fit_gamma: fit_gamma.cc digamma.o util.o GetOpt.o Awk.o matrix.o GammaFunc.o Matrix.h Array.h  GetOpt.h util.h GammaFunc.h
	$(CC) $(CFLAGS) -o fit_gamma fit_gamma.cc digamma.o util.o GetOpt.o Awk.o matrix.o gzstream.o GammaFunc.o -lz
hashstats: hashstats.cc Awk.o util.o GetOpt.o KeyIndex.h Array.h
	$(CC) $(CFLAGS) -o hashstats hashstats.cc Awk.o util.o GetOpt.o gzstream.o -lz 
EM.o: EM.cc EM.h Array.h Matrix.h matrix.cc matrix.o
	$(CC) $(CFLAGS) -c EM.cc 
Gaussian.o: Gaussian.cc Gaussian.h Matrix.h matrix.cc matrix.o
	$(CC) $(CFLAGS) -c Gaussian.cc
testEM: testEM.cc EM.o GetOpt.o util.o Array.h Matrix.h matrix.cc Awk.o Gaussian.o
	$(CC) $(CFLAGS) -o testEM testEM.cc EM.o GetOpt.o util.o matrix.o Awk.o Gaussian.o -lz -lm 
Cstrings.o: Cstrings.cc Cstrings.h KeyIndex.h util.o
	$(CC) $(CFLAGS) -c Cstrings.cc
Student_t.o: Student_t.cc
	$(CC) $(CFLAGS) -c Student_t.cc 
test_Cstrings: test_Cstrings.cc Cstrings.o
	$(CC) $(CFLAGS) -o test_Cstrings test_Cstrings.cc Cstrings.o util.o
netflow_stats: netflow_stats.cc Cstrings.o Awk.o Array.h Matrix.h Heap.h GetOpt.o util.o Student_t.o
	$(CC) $(CFLAGS) -o netflow_stats netflow_stats.cc Cstrings.o GetOpt.o Awk.o util.o Student_t.o \
-lz 
busiest_stats: busiest_stats.cc Cstrings.o Awk.o Array.h Heap.h GetOpt.o util.o matrix.o Student_t.o
	$(CC) $(CFLAGS) -o busiest_stats busiest_stats.cc Cstrings.o GetOpt.o Awk.o util.o matrix.o Student_t.o -lz
get_connections: get_connections.cc Array.h Heap.h KeyIndex.h Cstrings.o Awk.o GetOpt.o util.o
	$(CC) $(CFLAGS) -o get_connections get_connections.cc Cstrings.o GetOpt.o Awk.o util.o -lz
get_quantiles: get_quantiles.cc Array.h Matrix.h Heap.h KeyIndex.h ParseNetflowRecord.h Cstrings.o Awk.o GetOpt.o util.o matrix.o
	$(CC) $(CFLAGS) -o get_quantiles get_quantiles.cc Cstrings.o GetOpt.o Awk.o util.o matrix.o -lz
k-means: k-means.cc Array.h Awk.o GetOpt.o util.o
	$(CC) $(CFLAGS) -o k-means k-means.cc GetOpt.o Awk.o util.o -lz
test_MatrixWelford: test_MatrixWelford.cc matrix.o util.o Awk.o
	g++ $(CFLAGS) -o test_MatrixWelford test_MatrixWelford.cc matrix.o util.o Awk.o -lz
compare_gammas: compare_gammas.cc Awk.o Matrix.h GetOpt.o util.o Awk.o matrix.o
	g++ $(CFLAGS) -o compare_gammas compare_gammas.cc GetOpt.o util.o matrix.o Awk.o -lz -lm
PEcluster: PEcluster.cc Array.h Matrix.h Heap.h Awk.o GetOpt.o util.o matrix.o EM.o Gaussian.o
	$(CC) $(CFLAGS) -fopenmp -o PEcluster PEcluster.cc GetOpt.o Awk.o util.o matrix.o EM.o Gaussian.o -lz
test_Matrix: test_Matrix.cc Matrix.h util.o matrix.o
	$(CC) $(CFLAGS) -o test_Matrix test_Matrix.cc util.o matrix.o
rare_webips: rare_webips.cc Matrix.h Heap.h Awk.o util.o GetOpt.o
	$(CC) $(CFLAGS) -o rare_webips rare_webips.cc GetOpt.o util.o Awk.o\
 -lz
assign_cluster: assign_cluster.cc GetOpt.o util.o Array.h Heap.h Matrix.h Interpolate.h ParseNetflowRecord.h Gaussian.o Awk.o matrix.o
	$(CC) $(CFLAGS) -o assign_cluster assign_cluster.cc GetOpt.o util.o Awk.o Gaussian.o matrix.o -lz
score_for_malware: score_for_malware.cc GetOpt.o util.o Array.h Matrix.h Heap.h Interpolate.h Gaussian.o Awk.o matrix.o
	$(CC) $(CFLAGS) -o score_for_malware score_for_malware.cc GetOpt.o util.o Awk.o Gaussian.o matrix.o -lz
lookup_webip: lookup_webip.cc GetOpt.o util.o Heap.h Array.h BinarySearch.h Awk.o
	$(CC) $(CFLAGS) -o lookup_webip lookup_webip.cc GetOpt.o util.o Awk.o -lz
NewCstrings.o: NewCstrings.cc NewCstrings.h Array.h KeyIndex.h util.h
	$(CC) $(CFLAGS) -c NewCstrings.cc
get_wsa: get_wsa.cc NewCstrings.o util.o GetOpt.o Awk.o util.o
	$(CC) $(CFLAGS) -o get_wsa get_wsa.cc NewCstrings.o GetOpt.o util.o Awk.o -lz
get_fts_flowsets: get_fts_flowsets.cc ParseNetflowRecord.h stats.h Matrix.h Array.h FlexHeap.h Awk.o GetOpt.o util.o NewCstrings.o 
	$(CC) $(CFLAGS) -o get_fts_flowsets get_fts_flowsets.cc Awk.o GetOpt.o util.o NewCstrings.o -lz
wsa_lookup: wsa_lookup.cc Array.h BinarySearch.h GetOpt.o util.o Awk.o
	$(CC) $(CFLAGS) -o wsa_lookup wsa_lookup.cc GetOpt.o util.o Awk.o -lz
features2cdf: features2cdf.cc Array.h Matrix.h FlexHeap.h Awk.o GetOpt.o util.o matrix.o 
	$(CC) $(CFLAGS) -o features2cdf features2cdf.cc GetOpt.o Awk.o util.o matrix.o -lz
score_fts_flowsets: score_fts_flowsets.cc GetOpt.o util.o Array.h Matrix.h Gaussian.o Awk.o matrix.o
	$(CC) $(CFLAGS) -o score_fts_flowsets score_fts_flowsets.cc GetOpt.o util.o Awk.o Gaussian.o matrix.o -lz
ref_histogram: ref_histogram.cc util.o Awk.o Array.h GetOpt.o
	$(CC) $(CFLAGS) -o ref_histogram ref_histogram.cc GetOpt.o util.o Awk.o -lz
hawkes_model: hawkes.cc matrix.o util.o Awk.o GetOpt.o stats.h Dict.h Array.h Matrix.h
	$(CC) $(CFLAGS) -o hawkes_model hawkes.cc matrix.o util.o Awk.o GetOpt.o libcrc/lib/libcrc.a -lz
simulate: simulate.cc matrix.o util.o Awk.o GetOpt.o stats.h Array.h Matrix.h FlexHeap.h
	$(CC) $(CFLAGS) -o simulate simulate.cc matrix.o util.o Awk.o GetOpt.o -lz
test_plot: test_plot.cc
	$(CC) $(CFLAGS) -I/usr/include/python3.8 -o test_plot test_plot.cc -lpython3.8



