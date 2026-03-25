#include <cstdlib>
#include <cstdio>
#include "KC_pu.h"
#include <climits>
extern "C" {
#include "dds_mapstore.h"
#include "dds_psketch.h"
}


int ts = 0;
int collapse_times = 0;
double W;
double theta_;
double skip;
unsigned long long next_sample_index = 0;
unsigned s_num;

#define FILTERING 1
//#define FAST_SAMPLE_PLUS 1

// Global
unsigned m_skips[1 << 16];
unsigned short skip_idx;
unsigned skip_remaining;
double p = 0.001;

double w_skips[1 << 16];
unsigned short w_skip_idx;
unsigned w_skip_remaining;


double logrnd_skips[1 << 16];
unsigned short logrnd_skip_idx;
unsigned logrnd_skip_remaining;


unsigned rand_skips[1 << 16];
unsigned short rand_skip_idx;
unsigned rand_skip_remaining;


// To run before any insert (c'tor)
void init_skips(double p) {
    // pick packets w.p. p
    for (int i = 0; i < (1 << 16); ++i) {
        m_skips[i] = log(g_rand.uniform()) / log(1 - p) + 1;
    }
    skip_remaining = m_skips[0];
    skip_idx = 0;
}
long hash31(int64_t a, int64_t b, int64_t x)
{

    int64_t result;
    long lresult;
    result=(a * x) + b;
    result = ((result >> 31) + result) & 2147483647;
    lresult=(long) result;

    return(lresult);
}

void init_fastSampling(int sample_size) {
    double r;
    r = g_rand.uniform() ;
    W =  exp(log(r)/ (double) sample_size);
    next_sample_index = 0;
    skip = 0.0;

    for (int i = 0; i < (1 << 16); ++i) {
        r = g_rand.uniform();

        w_skips[i] = exp(log(r) / (double) sample_size);
        r = g_rand.uniform();

        logrnd_skips[i] = log(r);

        rand_skips[i] =  g_rand.uniform_int(0, sample_size - 1);
    }

    w_skip_idx = 0;
    logrnd_skip_idx = 0;
    rand_skip_idx = 0;
}
//link each hash bucket to a corresponding KCUITEM
void KCU_LinkHashtableToItems(KCU_type *kcu) {
    if (!kcu || !kcu->hashtable || !kcu->items) {
        return;
    }

    for (int i = 0; i < kcu->k; ++i) {
        // 1 to 1
        kcu->hashtable[i] = &kcu->items[i];
    }
}
#ifdef ORG_COLLAPSE
    int buckets_per_sketch = 100;
#endif

KCU_type *KCU_Init(float fPhi, double alpha,unsigned long samples_num, int sketch_bucket_num) {
    int i;
    int k = 1 + (int) 1.0 / fPhi;
    theta_ = fPhi;
    KCU_type *result = (KCU_type *) calloc(1, sizeof(KCU_type));
    s_num = samples_num;
    result->a = (long long) 698124007;
    result->b = (long long) 5125833;
    if (k < 1) k = 1;
    result->k = k;
    result->max_buckets = sketch_bucket_num;
#ifdef ORG_COLLAPSE
    buckets_per_sketch = sketch_bucket_num/k;
#endif


#ifdef FAST_SAMPLE_PLUS
    init_fastSampling(samples_num);
#endif

    //k hash bucket
    result->hashtable = (KCUITEM **) calloc(result->k, sizeof(KCUITEM *));
    result->items = (KCUITEM *) calloc(k, sizeof(KCUITEM));
    KCU_LinkHashtableToItems(result);

    for (i = 0; i < k; i++) {

        result->items[i].l_x=0;
        result->items[i].index_map = std::unordered_map<unsigned int, uint8_t>();
        result->items[i].index_map.clear();
        result->items[i].max_cell.item = UINT_MAX;
        result->items[i].max_cell.epoch = 0;
        result->items[i].max_cell.vect_cnt = 0;
        int j;
        for (j = 0; j < KCU_CELL_SIZE; j++) {
            result->items[i].cells[j].item = UINT_MAX;
            result->items[i].cells[j].epoch = 0;
            result->items[i].cells[j].vect_cnt = 0;


#ifdef RS_IN_CELL
            result->items[i].cells[j].samples = std::vector<sampleEntry>();
            result->items[i].cells[j].samples.clear();
            result->items[i].cells[j].nb_processed=0;
#endif

        }
        result->items[i].alpha = alpha;
        result->items[i].sketch = dds_psketch_init(COLLAPSINGALL_MAPSTORE, alpha, 0.0,INT_MAX);

    }
    return (result);
}


void dds_sketch_update_mask(dds_psketch* sketch , double value , long count) {
    dds_sketch_update((dds_sketch*)sketch, value, count);
#ifdef ORG_COLLAPSE
    if ( dds_store_get_num_buckets(sketch->store)>buckets_per_sketch+1) {

        collapse_times++;
        dds_mapstore_collapse((struct dds_mapstore *) sketch->store);

    }
#endif

}
size_t KCU_ComputeMemoryWithoutDDSketch(KCU_type *kcu) {
    // struct
    size_t total_memory = 0;
    total_memory += sizeof(KCU_type);
    total_memory += sizeof(KCU_item) * kcu->k;
    total_memory += sizeof(KCU_item*) * kcu->k;
    total_memory += (sizeof(unsigned int)+sizeof(uint8_t))*KCU_CELL_SIZE*kcu->k;

    return total_memory;
}
std::vector<int> KCU_ComputeMem(KCU_type *kcu) {
    int mem = 0;
    //bucket
    int bucket_memory_usage = kcu->max_buckets * 12; //long+int

    int udd_memory_usage = kcu->k * (56+56+32);//psketch+mapsotre+mapping

    int KCU_memory =(int) KCU_ComputeMemoryWithoutDDSketch(kcu);

    int RS_memory = 0;
#ifdef RS_IN_CELL
    //RS_IN_CELL
    RS_memory += (sizeof(std::vector<sampleEntry>)+ sizeof(unsigned))* kcu->k;//cell attribute
    RS_memory += (int) (pow(0.01, -1)) * kcu->k*sizeof(sampleEntry); //samplesize

#endif
    mem = bucket_memory_usage + udd_memory_usage + KCU_memory+ RS_memory;
    std::cout << "entries taking in account " << kcu->k << ", sketch mem:[Mbyte] " << mem / (double) (
        pow(10.0, 6.0)) << ","<<
            bucket_memory_usage <<","<<
            udd_memory_usage <<","<<
            KCU_memory <<","<<
            RS_memory <<","<<

            std::endl;
    std::vector<int> v;
    v.emplace_back(static_cast<int>(mem));
    v.emplace_back(bucket_memory_usage);
    v.emplace_back(udd_memory_usage);
    v.emplace_back(KCU_memory);
    v.emplace_back(RS_memory);
    return v;
}

#ifdef RS_IN_CELL
//https://en.wikipedia.org/wiki/Reservoir_sampling#An_optimal_algorithm
void KCU_addToSample(KCU_cell *kcu, double latency) {
    ++kcu->nb_processed;
    sampleEntry se;
    se.latency = latency;
    se.epoch = ts/KCU_EPOCH_SIZE +1;
    if (kcu->samples.empty()) {
        kcu->nb_processed = 0;
    }

    if (s_num >= kcu->nb_processed) {
        kcu->samples.push_back(se);
    } else {
        unsigned long j =g_rand.uniform_int(0, kcu->nb_processed - 1);
        if (j < s_num) {
            kcu->samples[j] = se;
        }
    }
}



void KCU_addToSampleFAST(KCU_cell *kcu, double latency) {
    ++kcu->nb_processed;
    sampleEntry se;
    se.latency = latency;
    se.epoch = ts/KCU_EPOCH_SIZE +1;

    if (kcu->samples.empty()) {
        kcu->nb_processed = 0;
    }
    if (s_num >= kcu->nb_processed) {
        kcu->samples.push_back(se);
        if (s_num == kcu->samples.size()) {
            skip = floor(logrnd_skips[logrnd_skip_idx++ %(1<<16)] / log(1.0-W));
            next_sample_index = kcu->nb_processed +skip +1;
        }
        return;
    }

    if (kcu->nb_processed>=next_sample_index) {

        int j = rand_skips[rand_skip_idx++ %(1<<16)];
        kcu->samples[j] = se;

        W *= w_skips[w_skip_idx++ %(1<<16)];
        skip = floor(logrnd_skips[logrnd_skip_idx++ %(1<<16)] / log(1.0-W));
        next_sample_index = kcu->nb_processed +skip +1;
    }

}
#endif



int sort_cell_reverse(KCU_item *item,int current_pos) {
    while (current_pos < KCU_CELL_SIZE) {
        //swap
        item->index_map[item->cells[current_pos].item] = current_pos - 1;
        std::swap(item->cells[current_pos], item->cells[current_pos - 1]);
        current_pos++;
    }
    return current_pos;
}

int sort_cell(int current_pos, KCU_item *item) {
    while (current_pos > 0 && item->cells[current_pos].vect_cnt > item->cells[current_pos - 1].vect_cnt) {
        //swap
        item->index_map[item->cells[current_pos].item] = current_pos - 1;
        item->index_map[item->cells[current_pos - 1].item] = current_pos;
        std::swap(item->cells[current_pos], item->cells[current_pos - 1]);
        current_pos--;
    }
    return current_pos;
}

void KCU_cells_update(KCU_item *item, unsigned int new_item, unsigned int new_epoch, double latency) {
    auto it = item->index_map.find(new_item);
    // if existed, increment
    if (it != item->index_map.end()) {
        int item_index = it->second;

        if (item->cells[item_index].epoch < new_epoch) {
            item->cells[item_index].item = new_item;
            item->cells[item_index].vect_cnt++;
            item->cells[item_index].epoch = new_epoch;

        }
#ifdef  RS_IN_CELL
    // anyway, sample
        KCU_addToSample(&item->cells[item_index],latency);
#endif
        // move
        int current_pos = item_index;
        current_pos = sort_cell(current_pos, item);

        // process max_cell
        if (current_pos == 0 && item->cells[0].vect_cnt > std::max(3,1 + (int) (theta_ * new_epoch))) {
            // enter max_cell
            if (item->max_cell.item == UINT_MAX) {
                item->max_cell = item->cells[0];
                item->cells[0] = KCU_cell{UINT_MAX, 0, 0};
#ifdef RS_IN_CELL
                //initial for safe
                item->cells[0].samples = std::vector<sampleEntry>();
                item->cells[0].nb_processed=0;
#endif

                item->index_map.erase(item->cells[0].item);

                item->index_map[item->cells[current_pos].item] = KCU_CELL_SIZE - 1;
                sort_cell_reverse(item,current_pos+1);

                //update lx and enter sketch
                item->l_x = ts;
                dds_sketch_update_mask(item->sketch, latency, 1);

            } else {
                //replace
                if (item->max_cell.vect_cnt < item->cells[0].vect_cnt) {

                    double thresh = KCU_epsilon2 * (new_epoch - item->max_cell.epoch) /
                                    (item->max_cell.vect_cnt+1);
                    if (g_rand.uniform()< thresh) {

                        item->index_map.erase(item->cells[0].item);
                        item->index_map[item->max_cell.item] = 0;

                        std::swap(item->max_cell, item->cells[0]);
                        dds_psketch_destroy((dds_psketch *)item->sketch);

                        item->sketch = dds_psketch_init(COLLAPSINGALL_MAPSTORE, item->alpha, 0.0,INT_MAX);


                        dds_sketch_update_mask(item->sketch, latency, 1);

                        //update lx
                        item->l_x = ts;
                    }
                }
            }
        }
    } else {
        int last_index = KCU_CELL_SIZE - 1;
        // insert into empty cell
        if (item->cells[last_index].item == UINT_MAX) {
            item->index_map[new_item] = last_index;
            item->cells[last_index].item = new_item;
            item->cells[last_index].vect_cnt = 1;
            item->cells[last_index].epoch = new_epoch;
#ifdef  RS_IN_CELL
            // insert into empty, don't need to initial
            KCU_addToSample(&item->cells[last_index],latency);
#endif

            sort_cell(last_index, item);

        } else {
            //full
            int rand_index = g_rand.uniform_int(1, KCU_CELL_SIZE - 1);
            if ((item->cells[rand_index].epoch>=new_epoch))return;
            double thresh =1-exp(2.0*KCU_epsilon
                *(-1.0*(new_epoch-item->cells[rand_index].epoch ))/item->cells[rand_index].vect_cnt);
            if (g_rand.uniform()<thresh) {
                //decrement and replace
              //  if (--item->cells[rand_index].vect_cnt <= 0)
            //    {

                    item->index_map.erase(item->cells[rand_index].item);
                    item->cells[rand_index].item = new_item;
                    item->cells[rand_index].vect_cnt = 1;
                    item->cells[rand_index].epoch = new_epoch;

#ifdef  RS_IN_CELL
                // should sample if replace, but need to reset
                item->cells[rand_index].samples.clear();
                item->cells[rand_index].nb_processed = 0;
                KCU_addToSample(&item->cells[rand_index],latency);
#endif

                    //move to last
                    item->index_map[new_item] = KCU_CELL_SIZE - 1;
                    sort_cell_reverse(item,rand_index+1);
           //     }
            }
        }
    }
}


void KCU_UpdateLatencywGSample(KCU_type *kcu, unsigned int newitem, double latency) {
    if (--skip_remaining == 0) {
        KCU_UpdateLatency(kcu, newitem, latency);
        skip_remaining = m_skips[++skip_idx];
    }
}


void KCU_collapse(KCU_type *kcu) {
    int cnt[kcu->k];
    int num_buckets = 0;

    for (int i = 0; i < kcu->k; i++) {
        if (kcu->items[i].sketch == nullptr) {
            cnt[i] = INT_MAX;
            continue;
        }

        int b = dds_store_get_num_buckets(kcu->items[i].sketch->store);
        num_buckets += b;

        if (b <= kcu->max_buckets/kcu->k) {
            cnt[i] = INT_MAX;
            continue;
        }

        cnt[i] = -b;
    }

    int minIndex = 0;


    while (num_buckets > kcu->max_buckets) {

        // find target
        minIndex = 0;
        bool found_valid = false;
        for (int i = 0; i < kcu->k; i++) {
            if (cnt[i] != INT_MAX) {
                if (!found_valid || cnt[i] < cnt[minIndex]) {
                    minIndex = i;
                    found_valid = true;
                }
            }
        }

        if (!found_valid) {
            break;
        }


        cnt[minIndex] = INT_MAX;
        int before_buckets = dds_store_get_num_buckets(kcu->items[minIndex].sketch->store);
        num_buckets -= before_buckets;
        dds_mapstore_collapse((struct dds_mapstore *) kcu->items[minIndex].sketch->store);
        int after_buckets = dds_store_get_num_buckets(kcu->items[minIndex].sketch->store);
        ++collapse_times;
        num_buckets += after_buckets;
        struct dds_bucket_id_mapping *mapping = dds_sketch_get_id_mapping(
            (struct dds_sketch *) kcu->items[minIndex].sketch);
        if (mapping == nullptr) {
            continue;
        }

        double gamma_value = mapping->gamma_value;
        double new_gamma = pow(gamma_value, 2.0);
        dds_bucket_id_mapping_update(mapping, new_gamma);

    }

}

uint16_t make_fingerprint(int64_t a, int64_t b,int value) {
    int64_t result;

    result = (a * (int64_t)(value)) + b;
    result = ((result >> 8) + result) & 0xFFFF;

    return (uint16_t)(result);
}

void KCU_UpdateLatency(KCU_type *kcu, unsigned int newitem, double latency) {
    int h;
    KCUITEM *il;
    ++ts;
#ifndef ORG_COLLAPSE
    if (ts % KCU_EPOCH_SIZE == 0) {
        KCU_collapse(kcu);
    }
#endif

    h = hash31(kcu->a, kcu->b, newitem) % kcu->k;
    il = kcu->hashtable[h];
    unsigned epoch = ts / KCU_EPOCH_SIZE + 1;

    bool should_skip = false;
    //first check max_cell
    if (il->max_cell.item == newitem) {
        if (il->max_cell.epoch < epoch) {
            il->max_cell.vect_cnt++;
            il->max_cell.epoch = epoch;
        }
        //whenever insert

        dds_sketch_update_mask(il->sketch, latency, 1);

#ifdef  RS_IN_CELL
        KCU_addToSample(&il->max_cell,latency);
#endif
    } else {

        KCU_cells_update(il, newitem, epoch, latency);

    }


}


double KCU_QuantileQuery(KCU_type *kcu,unsigned int newitem, double rank) {
    int h;
    KCUITEM *il;

    h = hash31(kcu->a, kcu->b, newitem) % kcu->k;
    il = kcu->hashtable[h];
    if (newitem != il->max_cell.item) // item is not monitored (not in hashtable)
        return -1;

    bool is_accurate = true;

    double val = dds_sketch_get_quantile((dds_sketch *) il->sketch, rank, &is_accurate);
    return val;
}

struct dds_psketch *dds_psketch_copy(const struct dds_psketch *source) {
    if (!source) {
        return NULL;
    }
    // parameter
    enum store_type store_type = source->store_type;
    double alpha = source->id_map->alpha;
    double min_addressable_value = source->min_addressable_value;
    int max_store_size = ((struct dds_mapstore *) source->store)->max_num_buckets;

    // create copy
    struct dds_psketch *copy = dds_psketch_init(store_type, alpha, min_addressable_value, max_store_size);
    if (!copy) {
        return NULL;
    }
    // copy zero_bucket
    copy->zero_bucket = source->zero_bucket;

    // store data
    struct dds_mapstore *source_store = (struct dds_mapstore *) source->store;
    struct dds_mapstore *copy_store = (struct dds_mapstore *) copy->store;

    // all buckets
    dict_itor *itor = dict_itor_new(source_store->buckets);
    for (dict_itor_first(itor); dict_itor_valid(itor); dict_itor_next(itor)) {
        int bid = *(int *) dict_itor_key(itor);
        long count = **(long **) dict_itor_datum(itor);

        // store in new bucket
        dds_store_insert(copy->store, bid, count);
    }
    dict_itor_free(itor);

    // for safe
    copy_store->total_counts = source_store->total_counts;
    copy_store->num_collapses = source_store->num_collapses;
    copy_store->is_collapsed = source_store->is_collapsed;

    return copy;
}

#ifdef RS_IN_CELL
dds_psketch *dds_psketch_tmp(KCU_type *kcu, unsigned int newitem) {
    int h;
    KCUITEM *il;
    std::vector<double> sample_buff;

    h = hash31(kcu->a, kcu->b, newitem) % kcu->k;
    il = kcu->hashtable[h];
    if (newitem != il->max_cell.item) // item is not in hashtable
        return nullptr;

    dds_psketch *tmpsketch = dds_psketch_copy(il->sketch);

    for (auto const &value: il->max_cell.samples) {
        if (value.epoch < il->l_x/KCU_EPOCH_SIZE +1) {
            //std::cout << "To sample buffer:! flow: " << value.flow <<" latency: " <<value.latency << " ts: " << value.ts << std::endl;
            sample_buff.push_back(value.latency);
            dds_sketch_update_mask((dds_sketch *) tmpsketch, value.latency,
                              (long) (il->max_cell.nb_processed/ il->max_cell.samples.size()));
        }
    }
    return tmpsketch;
}
#endif


double KCU_QuantileQuerywGsamples(dds_psketch *tmpsketch, double rank) {

    bool is_accurate = true;
    // with cell sample
    double val = dds_sketch_get_quantile((dds_sketch *) tmpsketch, rank, &is_accurate);
    return val;
}

std::map<uint32_t, uint32_t> KCU_Output(KCU_type *kcu) {
    return KCU_Output2(kcu).second;
}
std::pair<int, std::map<uint32_t, uint32_t>> KCU_Output2(KCU_type * kcu) {
    int buckets = 0;
    std::map<uint32_t, uint32_t> res;

    for (int i = 0; i < kcu->k; ++i)
        if (!dds_sketch_is_empty((dds_sketch *) kcu->items[i].sketch)) {
            res.insert(std::pair<uint32_t, uint32_t>(kcu->items[i].max_cell.item, kcu->items[i].max_cell.vect_cnt));
           // collision_times += dds_sketch_get_num_collapses((dds_sketch *) kcu->items[i].sketch);
            buckets +=dds_store_get_num_buckets(kcu->items[i].sketch->store);
        }
    return std::make_pair(buckets, res);
}
std::map<uint32_t, uint32_t> KCU_Output3(KCU_type * kcu) {
    std::map<uint32_t, uint32_t> res;

    for (int i = 0; i < kcu->k; ++i){

      //  for (KCU_cell c :kcu->items[i].cells) {
            res.insert(std::pair<uint32_t, uint32_t>(
                kcu->items[i].cells[0].item,
                kcu->items[i].cells[0].vect_cnt));
    //    }
        if (!dds_sketch_is_empty((dds_sketch *) kcu->items[i].sketch)) {
            res.insert(std::pair<uint32_t, uint32_t>
                (kcu->items[i].max_cell.item, kcu->items[i].max_cell.vect_cnt));

        }
    }
    return res;
}

void KCU_Destroy(KCU_type *kcu) {
    for (int i = 0; i < kcu->k; i++) {
        delete(kcu->items[i].sketch); // should done
    }
    free(kcu->items);
    free(kcu->hashtable);
    free(kcu);
}
