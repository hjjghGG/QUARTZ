#pragma once

#include <unordered_map>
#include <map>
#include <random>
#include <chrono>
#include <cstdint>
#include <ostream>
#include <iostream>


//#define RS_IN_CELL 1
//#define ORG_COLLAPSE 1
#define KCU_CELL_SIZE 5
#define KCU_EPOCH_SIZE 5000
#define KCU_epsilon  0.015
#define KCU_epsilon2 0.3


typedef struct KCU_item KCUITEM;

struct RandomGenerator {
    std::mt19937 engine;
    std::uniform_real_distribution<double> uniform_dist;

    RandomGenerator() : uniform_dist(0.0, 1.0) {
        unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        engine.seed(seed);
    }

    double uniform() {
        return uniform_dist(engine);
    }

    long long uniform_int(long long min, long long max) {
        std::uniform_int_distribution<long long> dist(min, max);
        return dist(engine);
    }
};

// random
static RandomGenerator g_rand;

class sampleEntry
{
public:
    float latency;
    unsigned short epoch;
};

struct KCU_cell
{
  unsigned int item;
  unsigned int epoch;//last update
  int vect_cnt;//cnt

#ifdef RS_IN_CELL
    std::vector<sampleEntry> samples;
    unsigned int nb_processed;
#endif
};


struct KCU_item
{
  struct  KCU_cell cells[KCU_CELL_SIZE];//cells
  std::unordered_map<unsigned int, uint8_t> index_map;
  KCU_cell max_cell;
  double alpha;
  struct dds_psketch* sketch;//UDDSketch
    unsigned int l_x;


};


typedef struct KCU_type{

  int k;//size
  int max_buckets; //max_bucket
  long long a,b;
#ifdef KCU_SIZE
  KCUITEM items[KCU_SIZE];
  KCUITEM* hashtable[KCU_TBLSIZE];
#else
  KCUITEM * items;
  KCUITEM **hashtable;

#endif
} KCU_type;

extern  KCU_type * KCU_Init(float fPhi, double alpha,unsigned long samples_num, int sketch_bucket_num);
extern void KCU_Destroy(KCU_type *);
extern void KCU_UpdateLatency(KCU_type *, unsigned int, double);
extern void KCU_UpdateLatencywGSample(KCU_type * kcu, unsigned int newitem, double latency);
extern double KCU_QuantileQuery(KCU_type *kcu,unsigned int newitem, double rank) ;
#ifdef RS_IN_CELL
extern dds_psketch *dds_psketch_tmp(KCU_type *kcu, unsigned int newitem);
#endif
extern double KCU_QuantileQuerywGsamples(dds_psketch *tmpsketch, double rank);
extern std::vector<int> KCU_ComputeMem(KCU_type * kcu);
extern void KCU_cells_update(KCU_item *item, unsigned int new_item, unsigned int new_epoch, double latency);
extern void KCU_collapse(KCU_type * kcu);
extern std::map<uint32_t, uint32_t> KCU_Output(KCU_type * kcu);
extern std::pair<int, std::map<uint32_t, uint32_t>> KCU_Output2(KCU_type * kcu);
extern std::map<uint32_t, uint32_t> KCU_Output3(KCU_type * kcu);
