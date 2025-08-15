#define _GNU_SOURCE


#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <pthread.h>
#include <sched.h>
#include <unistd.h>


#include "filter.h"
#include "signal.h"
#include "timing.h"


#define MAXWIDTH 40
#define THRESHOLD 2.0
#define ALIENS_LOW  50000.0
#define ALIENS_HIGH 150000.0

volatile int shared_band = 0;
typedef struct {
  // no start or end bands?
    int thread_id;
    //int start_band;
    //int end_band;
    int filter_order;
    int num_bands;
    double bandwidth;
    double Fs;
    signal* sig;
    double* band_power;
    //int* shared_band;
}
thread_data_t;


int num_threads;
int num_bands;
pthread_t* tid;
thread_data_t* thread_data;


void usage() {
  printf("usage: band_scan text|bin|mmap signal_file Fs filter_order num_bands num_threads num_processors\n");
}


double avg_power(double* data, int num) {


  double ss = 0;
  for (int i = 0; i < num; i++) {
    ss += data[i] * data[i];
  }


  return ss / num;
}


double max_of(double* data, int num) {


  double m = data[0];
  for (int i = 1; i < num; i++) {
    if (data[i] > m) {
      m = data[i];
    }
  }
  return m;
}


double avg_of(double* data, int num) {


  double s = 0;
  for (int i = 0; i < num; i++) {
    s += data[i];
  }
  return s / num;
}



void remove_dc(double* data, int num) {

  int i, j;
  double dc = avg_of(data,num);

  printf("Removing DC component of %lf\n",dc);


  for (i=0; i < num; i+=4){
    data[i] -= dc;
    data[i+1] -= dc;
    data[i+2] -= dc;
    data[i+3] -= dc;
  }
  for (j = num-(num%4); j<num; j++){
    data[j] -= dc;
  }
}



void* worker(void* arg) {
  thread_data_t* data = (thread_data_t*)arg;
  //long myid     = (long)arg;
  //int blocksize =  num_bands/ num_threads;
  // Bind thread to a specific core
  cpu_set_t set;
  CPU_ZERO(&set);
  CPU_SET(data->thread_id % sysconf(_SC_NPROCESSORS_ONLN), &set);
  
  int band;

  while (1){
    band = __sync_fetch_and_add(&shared_band, 1);
    
    if (band >= data->num_bands)
      break;
  // for (band = data->start_band; band < data->end_band; band++) {
    double filter_coeffs[data->filter_order + 1];
  
    double low_band = band * data->bandwidth + 0.0001;
    double high_band = (band + 1) * data->bandwidth - 0.0001;

    
    generate_band_pass(data->Fs,
                        low_band,
                        high_band,
                        data->filter_order,
                        filter_coeffs);


    hamming_window(data->filter_order, filter_coeffs);


    convolve_and_compute_power(data->sig->num_samples,
                                data->sig->data,
                                data->filter_order,
                                filter_coeffs,
                                &data->band_power[band]);
  }


  pthread_exit(NULL);
}



int analyze_signal(signal* sig, int filter_order, int num_bands, double* lb, double* ub) {
  double Fc = (sig->Fs) / 2;
  double bandwidth = Fc / num_bands;


  remove_dc(sig->data, sig->num_samples);


  double signal_power = avg_power(sig->data, sig->num_samples);


  printf("signal average power:     %lf\n", signal_power);


  resources rstart;
  get_resources(&rstart, THIS_PROCESS);
  double start = get_seconds();
  unsigned long long tstart = get_cycle_count();


  double* band_power = aligned_alloc(64, num_bands * sizeof(double));//[num_bands]; //double band_power[num_bands]; // fix to dynamic memory allocation?
  pthread_t threads[num_threads];
  thread_data_t thread_data[num_threads];

  //volatile int shared_band = 0;
  // Divide bands across threads

  //int bands_per_thread = num_bands / num_threads;
  for (int i = 0; i < num_threads; i++) {
    thread_data[i].thread_id = i;
    thread_data[i].num_bands = num_bands; 
    //thread_data[i].shared_band = &shared_band;
    //thread_data[i].start_band = i * bands_per_thread;
    //thread_data[i].end_band = (i == num_threads - 1) ? num_bands : (i + 1) * bands_per_thread;
    thread_data[i].Fs = sig->Fs;
    thread_data[i].bandwidth = bandwidth;
    thread_data[i].filter_order = filter_order;
    thread_data[i].sig = sig;
    thread_data[i].band_power = band_power;


    pthread_create(&threads[i], NULL, worker, &thread_data[i]);
  }


  // Wait for threads to finish
  for (int i = 0; i < num_threads; i++) {
    pthread_join(threads[i], NULL);
  }


  unsigned long long tend = get_cycle_count();
  double end = get_seconds();


  resources rend;
  get_resources(&rend, THIS_PROCESS);


  resources rdiff;
  get_resources_diff(&rstart, &rend, &rdiff);


  double max_band_power = max_of(band_power, num_bands);
  double avg_band_power = avg_of(band_power, num_bands);

  int wow = 0;
  *lb = -1;
  *ub = -1;

  
  for (int band = 0; band < num_bands; band++) {
    double band_low = band * bandwidth + 0.0001;
    double band_high = (band + 1) * bandwidth - 0.0001;


    printf("%5d %20lf to %20lf Hz: %20lf ",
            band, band_low, band_high, band_power[band]);


    for (int i = 0; i < MAXWIDTH * (band_power[band] / max_band_power); i++) {   // this is the issue
      printf("*");

    
    }


    if ((band_low >= ALIENS_LOW && band_low <= ALIENS_HIGH) ||
        (band_high >= ALIENS_LOW && band_high <= ALIENS_HIGH)) {


      // Band of interest
      if (band_power[band] > THRESHOLD * avg_band_power) {
        printf("(WOW)");
        wow = 1;
        if (*lb < 0) {
          *lb = band_low;
        }
        *ub = band_high;
      }
     
      else {
        printf("(meh)");
      }
    }
   
    else {
      printf("(meh)");
    }


    printf("\n");
  }


  printf("Resource usages:\n"
          "User time        %lf seconds\n"
          "System time      %lf seconds\n"
          "Page faults      %ld\n"
          "Page swaps       %ld\n"
          "Blocks of I/O    %ld\n"
          "Signals caught   %ld\n"
          "Context switches %ld\n",
          rdiff.usertime,
          rdiff.systime,
          rdiff.pagefaults,
          rdiff.pageswaps,
          rdiff.ioblocks,
          rdiff.sigs,
          rdiff.contextswitches);


  printf("Analysis took %llu cycles (%lf seconds) by cycle count, timing overhead=%llu cycles\n"
          "Note that cycle count only makes sense if the thread stayed on one core\n",
          tend - tstart, cycles_to_seconds(tend - tstart), timing_overhead());
  printf("Analysis took %lf seconds by basic timing\n", end - start);


  return wow;
}


int main(int argc, char* argv[]) {
  if (argc != 8) {
    usage();
    return -1;
  }


  char sig_type = toupper(argv[1][0]);
  char* sig_file = argv[2];
  double Fs = atof(argv[3]);
  int filter_order = atoi(argv[4]);
  num_bands = atoi(argv[5]);
  num_threads = atoi(argv[6]);


  printf("type:     %s\n"
        "file:     %s\n"
        "Fs:       %lf Hz\n"
        "order:    %d\n"
        "bands:    %d\n",
        sig_type == 'T' ? "Text" : (sig_type == 'B' ? "Binary" : (sig_type == 'M' ? "Mapped Binary" : "UNKNOWN TYPE")),
        sig_file,
        Fs,
        filter_order,
        num_bands);


printf("Load or map file\n");
  signal* sig;
  switch (sig_type) {
    case 'T': sig = load_text_format_signal(sig_file); break;
    case 'B': sig = load_binary_format_signal(sig_file); break;
    case 'M': sig = map_binary_format_signal(sig_file); break;
    default:
      printf("Unknown signal type\n");
      return -1;
  }


  if (!sig) {
    printf("Unable to load or map file\n");
    return -1;
  }

  sig->Fs = Fs;

  tid = malloc(sizeof(pthread_t) * num_threads);
  thread_data = malloc(sizeof(thread_data_t) * num_threads);

  double lb, ub;
  if (analyze_signal(sig, filter_order, num_bands, &lb, &ub)) {
    printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n", lb, ub, (ub + lb) / 2.0);
  }
 
  else {
    printf("No aliens detected.\n");
  }

  free_signal(sig);
  free(tid);
  free(thread_data);

  return 0;
}